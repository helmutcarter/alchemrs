use alchemrs::{
    analysis::{BlockEstimate, ConvergencePoint},
    prep::EquilibrationResult,
    DeltaFMatrix, DhdlSeries, FreeEnergyEstimate, OverlapMatrix, StatePoint, SwitchingTrajectory,
    UNkMatrix,
};
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::error::to_py_err;

#[pyclass(name = "StatePoint", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyStatePoint {
    pub(crate) inner: StatePoint,
}

#[pymethods]
impl PyStatePoint {
    #[new]
    fn new(lambdas: Vec<f64>, temperature_k: f64) -> PyResult<Self> {
        StatePoint::new(lambdas, temperature_k)
            .map(|inner| Self { inner })
            .map_err(to_py_err)
    }

    #[getter]
    fn lambdas(&self) -> Vec<f64> {
        self.inner.lambdas().to_vec()
    }

    #[getter]
    fn temperature_k(&self) -> f64 {
        self.inner.temperature_k()
    }
}

#[pyclass(name = "DhdlSeries", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyDhdlSeries {
    pub(crate) inner: DhdlSeries,
}

#[pymethods]
impl PyDhdlSeries {
    #[new]
    fn new(
        state: Py<PyStatePoint>,
        time_ps: &Bound<'_, PyAny>,
        values: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let py = time_ps.py();
        let state = state.borrow(py).inner.clone();
        let time_ps = extract_vec1(time_ps)?;
        let values = extract_vec1(values)?;
        DhdlSeries::new(state, time_ps, values)
            .map(|inner| Self { inner })
            .map_err(to_py_err)
    }

    #[getter]
    fn state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.state().clone(),
        }
    }

    #[getter]
    fn time_ps<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.time_ps().to_vec().into_pyarray(py)
    }

    #[getter]
    fn values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.values().to_vec().into_pyarray(py)
    }
}

#[pyclass(name = "UNkMatrix", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyUNkMatrix {
    pub(crate) inner: UNkMatrix,
}

#[pymethods]
impl PyUNkMatrix {
    #[new]
    #[pyo3(signature = (data, time_ps, evaluated_states, sampled_state=None, lambda_labels=None))]
    fn new(
        py: Python<'_>,
        data: &Bound<'_, PyAny>,
        time_ps: &Bound<'_, PyAny>,
        evaluated_states: Vec<Py<PyStatePoint>>,
        sampled_state: Option<Py<PyStatePoint>>,
        lambda_labels: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let matrix = extract_matrix(data)?;
        let n_samples = matrix.len();
        let n_states = matrix.first().map(|row| row.len()).unwrap_or(0);
        if matrix.iter().any(|row| row.len() != n_states) {
            return Err(PyValueError::new_err(
                "data must be rectangular with shape (n_samples, n_states)",
            ));
        }
        let flat = matrix.into_iter().flatten().collect::<Vec<_>>();
        let time_ps = extract_vec1(time_ps)?;
        let evaluated_states = evaluated_states
            .into_iter()
            .map(|state| state.borrow(py).inner.clone())
            .collect::<Vec<_>>();
        let sampled_state = sampled_state.map(|state| state.borrow(py).inner.clone());

        UNkMatrix::new_with_labels(
            n_samples,
            n_states,
            flat,
            time_ps,
            sampled_state,
            evaluated_states,
            lambda_labels,
        )
        .map(|inner| Self { inner })
        .map_err(to_py_err)
    }

    #[getter]
    fn data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        matrix_to_pyarray2(
            py,
            self.inner.n_samples(),
            self.inner.n_states(),
            self.inner.data(),
        )
    }

    #[getter]
    fn time_ps<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.time_ps().to_vec().into_pyarray(py)
    }

    #[getter]
    fn n_samples(&self) -> usize {
        self.inner.n_samples()
    }

    #[getter]
    fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    #[getter]
    fn sampled_state(&self) -> Option<PyStatePoint> {
        self.inner
            .sampled_state()
            .cloned()
            .map(|inner| PyStatePoint { inner })
    }

    #[getter]
    fn evaluated_states(&self) -> Vec<PyStatePoint> {
        self.inner
            .evaluated_states()
            .iter()
            .cloned()
            .map(|inner| PyStatePoint { inner })
            .collect()
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }
}

#[pyclass(name = "SwitchingTrajectory", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PySwitchingTrajectory {
    pub(crate) inner: SwitchingTrajectory,
}

#[pymethods]
impl PySwitchingTrajectory {
    #[new]
    #[pyo3(signature = (initial_state, final_state, reduced_work, lambda_path=None, dvdl_path=None, rms_dvdl_path=None))]
    fn new(
        py: Python<'_>,
        initial_state: Py<PyStatePoint>,
        final_state: Py<PyStatePoint>,
        reduced_work: f64,
        lambda_path: Option<&Bound<'_, PyAny>>,
        dvdl_path: Option<&Bound<'_, PyAny>>,
        rms_dvdl_path: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let initial_state = initial_state.borrow(py).inner.clone();
        let final_state = final_state.borrow(py).inner.clone();
        let lambda_path = lambda_path
            .map(extract_vec1)
            .transpose()?
            .unwrap_or_default();
        let dvdl_path = dvdl_path.map(extract_vec1).transpose()?.unwrap_or_default();
        let rms_dvdl_path = rms_dvdl_path
            .map(extract_vec1)
            .transpose()?
            .unwrap_or_default();

        SwitchingTrajectory::new_with_profile_and_rms(
            initial_state,
            final_state,
            reduced_work,
            lambda_path,
            dvdl_path,
            rms_dvdl_path,
        )
        .map(|inner| Self { inner })
        .map_err(to_py_err)
    }

    #[getter]
    fn initial_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.initial_state().clone(),
        }
    }

    #[getter]
    fn final_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.final_state().clone(),
        }
    }

    #[getter]
    fn reduced_work(&self) -> f64 {
        self.inner.reduced_work()
    }

    #[getter]
    fn lambda_path<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.lambda_path().to_vec().into_pyarray(py)
    }

    #[getter]
    fn dvdl_path<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.dvdl_path().to_vec().into_pyarray(py)
    }

    #[getter]
    fn rms_dvdl_path<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.rms_dvdl_path().to_vec().into_pyarray(py)
    }
}

#[pyclass(name = "FreeEnergyEstimate", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyFreeEnergyEstimate {
    pub(crate) inner: FreeEnergyEstimate,
}

#[pymethods]
impl PyFreeEnergyEstimate {
    #[getter]
    fn delta_f(&self) -> f64 {
        self.inner.delta_f()
    }

    #[getter]
    fn uncertainty(&self) -> Option<f64> {
        self.inner.uncertainty()
    }

    #[getter(from_state)]
    fn source_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.from_state().clone(),
        }
    }

    #[getter]
    fn to_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.to_state().clone(),
        }
    }
}

#[pyclass(name = "DeltaFMatrix", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyDeltaFMatrix {
    pub(crate) inner: DeltaFMatrix,
}

#[pymethods]
impl PyDeltaFMatrix {
    #[getter]
    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        matrix_to_pyarray2(
            py,
            self.inner.n_states(),
            self.inner.n_states(),
            self.inner.values(),
        )
    }

    #[getter]
    fn uncertainties<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray2<f64>>>> {
        self.inner
            .uncertainties()
            .map(|values| {
                matrix_to_pyarray2(py, self.inner.n_states(), self.inner.n_states(), values)
            })
            .transpose()
    }

    #[getter]
    fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    #[getter]
    fn states(&self) -> Vec<PyStatePoint> {
        self.inner
            .states()
            .iter()
            .cloned()
            .map(|inner| PyStatePoint { inner })
            .collect()
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }
}

#[pyclass(name = "OverlapMatrix", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyOverlapMatrix {
    pub(crate) inner: OverlapMatrix,
}

#[pymethods]
impl PyOverlapMatrix {
    #[getter]
    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        matrix_to_pyarray2(
            py,
            self.inner.n_states(),
            self.inner.n_states(),
            self.inner.values(),
        )
    }

    #[getter]
    fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    #[getter]
    fn states(&self) -> Vec<PyStatePoint> {
        self.inner
            .states()
            .iter()
            .cloned()
            .map(|inner| PyStatePoint { inner })
            .collect()
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }
}

#[pyclass(name = "EquilibrationResult", module = "alchemrs._alchemrs")]
#[derive(Clone, Copy)]
pub struct PyEquilibrationResult {
    pub(crate) inner: EquilibrationResult,
}

#[pymethods]
impl PyEquilibrationResult {
    #[getter]
    fn t0(&self) -> usize {
        self.inner.t0
    }

    #[getter]
    fn g(&self) -> f64 {
        self.inner.g
    }

    #[getter]
    fn neff_max(&self) -> f64 {
        self.inner.neff_max
    }
}

#[pyclass(name = "ConvergencePoint", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyConvergencePoint {
    pub(crate) inner: ConvergencePoint,
}

#[pymethods]
impl PyConvergencePoint {
    #[getter]
    fn n_windows(&self) -> usize {
        self.inner.n_windows()
    }

    #[getter]
    fn delta_f(&self) -> f64 {
        self.inner.delta_f()
    }

    #[getter]
    fn uncertainty(&self) -> Option<f64> {
        self.inner.uncertainty()
    }

    #[getter(from_state)]
    fn source_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.from_state().clone(),
        }
    }

    #[getter]
    fn to_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.to_state().clone(),
        }
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }
}

#[pyclass(name = "BlockEstimate", module = "alchemrs._alchemrs")]
#[derive(Clone)]
pub struct PyBlockEstimate {
    pub(crate) inner: BlockEstimate,
}

#[pymethods]
impl PyBlockEstimate {
    #[getter]
    fn block_index(&self) -> usize {
        self.inner.block_index()
    }

    #[getter]
    fn n_blocks(&self) -> usize {
        self.inner.n_blocks()
    }

    #[getter]
    fn delta_f(&self) -> f64 {
        self.inner.delta_f()
    }

    #[getter]
    fn uncertainty(&self) -> Option<f64> {
        self.inner.uncertainty()
    }

    #[getter(from_state)]
    fn source_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.from_state().clone(),
        }
    }

    #[getter]
    fn to_state(&self) -> PyStatePoint {
        PyStatePoint {
            inner: self.inner.to_state().clone(),
        }
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStatePoint>()?;
    module.add_class::<PyDhdlSeries>()?;
    module.add_class::<PyUNkMatrix>()?;
    module.add_class::<PySwitchingTrajectory>()?;
    module.add_class::<PyFreeEnergyEstimate>()?;
    module.add_class::<PyDeltaFMatrix>()?;
    module.add_class::<PyOverlapMatrix>()?;
    module.add_class::<PyEquilibrationResult>()?;
    module.add_class::<PyConvergencePoint>()?;
    module.add_class::<PyBlockEstimate>()?;
    Ok(())
}

fn extract_vec1(value: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    if let Ok(array) = value.extract::<PyReadonlyArray1<'_, f64>>() {
        return Ok(array.as_slice()?.to_vec());
    }
    value
        .extract::<Vec<f64>>()
        .map_err(|_| PyTypeError::new_err("expected a 1D float sequence or numpy.ndarray"))
}

fn extract_matrix(value: &Bound<'_, PyAny>) -> PyResult<Vec<Vec<f64>>> {
    if let Ok(array) = value.extract::<PyReadonlyArray2<'_, f64>>() {
        let view = array.as_array();
        return Ok(view.outer_iter().map(|row| row.to_vec()).collect());
    }
    value
        .extract::<Vec<Vec<f64>>>()
        .map_err(|_| PyTypeError::new_err("expected a 2D float sequence or numpy.ndarray"))
}

pub(crate) fn matrix_to_pyarray2<'py>(
    py: Python<'py>,
    n_rows: usize,
    n_cols: usize,
    values: &[f64],
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let array = Array2::from_shape_vec((n_rows, n_cols), values.to_vec())
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
    Ok(array.into_pyarray(py))
}

pub(crate) fn wrap_dhdl_series(inner: DhdlSeries) -> PyDhdlSeries {
    PyDhdlSeries { inner }
}

pub(crate) fn wrap_unk_matrix(inner: UNkMatrix) -> PyUNkMatrix {
    PyUNkMatrix { inner }
}

pub(crate) fn wrap_switching_trajectory(inner: SwitchingTrajectory) -> PySwitchingTrajectory {
    PySwitchingTrajectory { inner }
}

pub(crate) fn wrap_equilibration_result(inner: EquilibrationResult) -> PyEquilibrationResult {
    PyEquilibrationResult { inner }
}

pub(crate) fn wrap_free_energy_estimate(inner: FreeEnergyEstimate) -> PyFreeEnergyEstimate {
    PyFreeEnergyEstimate { inner }
}

pub(crate) fn wrap_delta_f_matrix(inner: DeltaFMatrix) -> PyDeltaFMatrix {
    PyDeltaFMatrix { inner }
}

pub(crate) fn wrap_overlap_matrix(inner: OverlapMatrix) -> PyOverlapMatrix {
    PyOverlapMatrix { inner }
}

pub(crate) fn wrap_block_estimate(inner: BlockEstimate) -> PyBlockEstimate {
    PyBlockEstimate { inner }
}

pub(crate) fn wrap_convergence_point(inner: ConvergencePoint) -> PyConvergencePoint {
    PyConvergencePoint { inner }
}

pub(crate) fn extract_vec1_any(value: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    extract_vec1(value)
}
