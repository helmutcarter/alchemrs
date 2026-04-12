use alchemrs::{
    BarEstimator, BarFit, BarMethod, BarOptions, IexpEstimator, IexpFit, IexpOptions,
    IntegrationMethod, MbarEstimator, MbarFit, MbarOptions, MbarSolver, NesEstimator, NesFit,
    NesOptions, TiEstimator, TiFit, TiOptions,
};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::data::{
    matrix_to_pyarray2, wrap_block_estimate, wrap_delta_f_matrix, wrap_free_energy_estimate,
    wrap_overlap_matrix, PyBlockEstimate, PyDhdlSeries, PyFreeEnergyEstimate, PyOverlapMatrix,
    PySwitchingTrajectory, PyUNkMatrix,
};
use crate::error::to_py_err;

#[pyclass(name = "TI", module = "alchemrs._alchemrs.estimators")]
#[derive(Clone)]
pub struct PyTiEstimator {
    inner: TiEstimator,
}

#[pymethods]
impl PyTiEstimator {
    #[new]
    #[pyo3(signature = (method="trapezoidal", parallel=false))]
    fn new(method: &str, parallel: bool) -> PyResult<Self> {
        Ok(Self {
            inner: TiEstimator::new(TiOptions {
                method: parse_integration_method(method)?,
                parallel,
            }),
        })
    }

    fn fit(&self, py: Python<'_>, series: Vec<Py<PyDhdlSeries>>) -> PyResult<PyTiFit> {
        let series = clone_dhdl_series(py, series);
        self.inner
            .fit(&series)
            .map(|inner| PyTiFit { inner })
            .map_err(to_py_err)
    }

    fn estimate(
        &self,
        py: Python<'_>,
        series: Vec<Py<PyDhdlSeries>>,
    ) -> PyResult<PyFreeEnergyEstimate> {
        let series = clone_dhdl_series(py, series);
        self.inner
            .estimate(&series)
            .map(wrap_free_energy_estimate)
            .map_err(to_py_err)
    }

    fn block_average(
        &self,
        py: Python<'_>,
        series: Vec<Py<PyDhdlSeries>>,
        n_blocks: usize,
    ) -> PyResult<Vec<PyBlockEstimate>> {
        let series = clone_dhdl_series(py, series);
        self.inner
            .block_average(&series, n_blocks)
            .map(|values| values.into_iter().map(wrap_block_estimate).collect())
            .map_err(to_py_err)
    }
}

#[pyclass(name = "TiFit", module = "alchemrs._alchemrs.estimators")]
pub struct PyTiFit {
    inner: TiFit,
}

#[pymethods]
impl PyTiFit {
    #[getter]
    fn method(&self) -> &'static str {
        self.inner.method().as_str()
    }

    #[getter]
    fn delta_f(&self) -> f64 {
        self.inner.delta_f()
    }

    #[getter]
    fn uncertainty(&self) -> Option<f64> {
        self.inner.uncertainty()
    }

    fn result(&self) -> PyResult<PyFreeEnergyEstimate> {
        self.inner
            .result()
            .map(wrap_free_energy_estimate)
            .map_err(to_py_err)
    }
}

#[pyclass(name = "BAR", module = "alchemrs._alchemrs.estimators")]
#[derive(Clone)]
pub struct PyBarEstimator {
    inner: BarEstimator,
}

#[pymethods]
impl PyBarEstimator {
    #[new]
    #[pyo3(signature = (method="false-position", maximum_iterations=10_000, relative_tolerance=1e-7, parallel=false))]
    fn new(
        method: &str,
        maximum_iterations: usize,
        relative_tolerance: f64,
        parallel: bool,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: BarEstimator::new(BarOptions {
                maximum_iterations,
                relative_tolerance,
                method: parse_bar_method(method)?,
                parallel,
            }),
        })
    }

    fn fit(&self, py: Python<'_>, windows: Vec<Py<PyUNkMatrix>>) -> PyResult<PyBarFit> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .fit(&windows)
            .map(|inner| PyBarFit { inner })
            .map_err(to_py_err)
    }

    fn estimate(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
    ) -> PyResult<crate::data::PyDeltaFMatrix> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .estimate(&windows)
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn block_average(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
        n_blocks: usize,
    ) -> PyResult<Vec<PyBlockEstimate>> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .block_average(&windows, n_blocks)
            .map(|values| values.into_iter().map(wrap_block_estimate).collect())
            .map_err(to_py_err)
    }
}

#[pyclass(name = "BarFit", module = "alchemrs._alchemrs.estimators")]
pub struct PyBarFit {
    inner: BarFit,
}

#[pymethods]
impl PyBarFit {
    #[getter]
    fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }

    fn result(&self) -> PyResult<crate::data::PyDeltaFMatrix> {
        self.inner
            .result()
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }
}

#[pyclass(name = "MBAR", module = "alchemrs._alchemrs.estimators")]
#[derive(Clone)]
pub struct PyMbarEstimator {
    inner: MbarEstimator,
}

#[pymethods]
impl PyMbarEstimator {
    #[new]
    #[pyo3(signature = (max_iterations=10_000, tolerance=1e-7, initial_f_k=None, solver="fixed-point", parallel=false))]
    fn new(
        max_iterations: usize,
        tolerance: f64,
        initial_f_k: Option<Vec<f64>>,
        solver: &str,
        parallel: bool,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: MbarEstimator::new(MbarOptions {
                max_iterations,
                tolerance,
                initial_f_k,
                parallel,
                solver: parse_mbar_solver(solver)?,
            }),
        })
    }

    fn fit(&self, py: Python<'_>, windows: Vec<Py<PyUNkMatrix>>) -> PyResult<PyMbarFit> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .fit(&windows)
            .map(|inner| PyMbarFit { inner })
            .map_err(to_py_err)
    }

    fn estimate(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
    ) -> PyResult<crate::data::PyDeltaFMatrix> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .estimate(&windows)
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn estimate_with_uncertainty(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
    ) -> PyResult<crate::data::PyDeltaFMatrix> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .estimate_with_uncertainty(&windows)
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn block_average(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
        n_blocks: usize,
    ) -> PyResult<Vec<PyBlockEstimate>> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .block_average(&windows, n_blocks)
            .map(|values| values.into_iter().map(wrap_block_estimate).collect())
            .map_err(to_py_err)
    }
}

#[pyclass(name = "MbarFit", module = "alchemrs._alchemrs.estimators")]
pub struct PyMbarFit {
    inner: MbarFit,
}

#[pymethods]
impl PyMbarFit {
    #[getter]
    fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }

    #[getter]
    fn free_energies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.free_energies().to_vec().into_pyarray(py)
    }

    #[getter]
    fn state_counts<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.state_counts().to_vec().into_pyarray(py)
    }

    #[getter]
    fn log_weights<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let n_states = self.inner.n_states();
        let n_samples = self.inner.log_weights().len() / n_states;
        matrix_to_pyarray2(py, n_samples, n_states, self.inner.log_weights())
    }

    fn result(&self) -> PyResult<crate::data::PyDeltaFMatrix> {
        self.inner
            .result()
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn result_with_uncertainty(&self) -> PyResult<crate::data::PyDeltaFMatrix> {
        self.inner
            .result_with_uncertainty()
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn overlap_matrix(&self) -> PyResult<PyOverlapMatrix> {
        self.inner
            .overlap_matrix()
            .map(wrap_overlap_matrix)
            .map_err(to_py_err)
    }

    fn overlap_eigenvalues<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        self.inner
            .overlap_eigenvalues()
            .map(|values| values.into_pyarray(py))
            .map_err(to_py_err)
    }

    fn overlap_scalar(&self) -> PyResult<f64> {
        self.inner.overlap_scalar().map_err(to_py_err)
    }
}

#[pyclass(name = "IEXP", module = "alchemrs._alchemrs.estimators")]
#[derive(Clone)]
pub struct PyIexpEstimator {
    inner: IexpEstimator,
}

#[pymethods]
impl PyIexpEstimator {
    #[new]
    #[pyo3(signature = (parallel=false))]
    fn new(parallel: bool) -> Self {
        Self {
            inner: IexpEstimator::new(IexpOptions { parallel }),
        }
    }

    fn fit(&self, py: Python<'_>, windows: Vec<Py<PyUNkMatrix>>) -> PyResult<PyIexpFit> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .fit(&windows)
            .map(|inner| PyIexpFit { inner })
            .map_err(to_py_err)
    }

    fn estimate(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
    ) -> PyResult<crate::data::PyDeltaFMatrix> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .estimate(&windows)
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn estimate_with_uncertainty(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
    ) -> PyResult<crate::data::PyDeltaFMatrix> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .estimate_with_uncertainty(&windows)
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn block_average(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
        n_blocks: usize,
    ) -> PyResult<Vec<PyBlockEstimate>> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .block_average(&windows, n_blocks)
            .map(|values| values.into_iter().map(wrap_block_estimate).collect())
            .map_err(to_py_err)
    }

    fn reverse_block_average(
        &self,
        py: Python<'_>,
        windows: Vec<Py<PyUNkMatrix>>,
        n_blocks: usize,
    ) -> PyResult<Vec<PyBlockEstimate>> {
        let windows = clone_unk_windows(py, windows);
        self.inner
            .reverse_block_average(&windows, n_blocks)
            .map(|values| values.into_iter().map(wrap_block_estimate).collect())
            .map_err(to_py_err)
    }
}

#[pyclass(name = "IexpFit", module = "alchemrs._alchemrs.estimators")]
pub struct PyIexpFit {
    inner: IexpFit,
}

#[pymethods]
impl PyIexpFit {
    #[getter]
    fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }

    fn result(&self) -> PyResult<crate::data::PyDeltaFMatrix> {
        self.inner
            .result()
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn result_with_uncertainty(&self) -> PyResult<crate::data::PyDeltaFMatrix> {
        self.inner
            .result_with_uncertainty()
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }
}

#[pyclass(name = "NES", module = "alchemrs._alchemrs.estimators")]
#[derive(Clone)]
pub struct PyNesEstimator {
    inner: NesEstimator,
}

#[pymethods]
impl PyNesEstimator {
    #[new]
    #[pyo3(signature = (n_bootstrap=0, seed=0))]
    fn new(n_bootstrap: usize, seed: u64) -> Self {
        Self {
            inner: NesEstimator::new(NesOptions { n_bootstrap, seed }),
        }
    }

    fn fit(
        &self,
        py: Python<'_>,
        trajectories: Vec<Py<PySwitchingTrajectory>>,
    ) -> PyResult<PyNesFit> {
        let trajectories = clone_trajectories(py, trajectories);
        self.inner
            .fit(&trajectories)
            .map(|inner| PyNesFit { inner })
            .map_err(to_py_err)
    }

    fn estimate(
        &self,
        py: Python<'_>,
        trajectories: Vec<Py<PySwitchingTrajectory>>,
    ) -> PyResult<PyFreeEnergyEstimate> {
        let trajectories = clone_trajectories(py, trajectories);
        self.inner
            .estimate(&trajectories)
            .map(wrap_free_energy_estimate)
            .map_err(to_py_err)
    }
}

#[pyclass(name = "NesFit", module = "alchemrs._alchemrs.estimators")]
pub struct PyNesFit {
    inner: NesFit,
}

#[pymethods]
impl PyNesFit {
    fn result(&self) -> PyResult<PyFreeEnergyEstimate> {
        self.inner
            .result()
            .map(wrap_free_energy_estimate)
            .map_err(to_py_err)
    }
}

pub fn register(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(py, "estimators")?;
    module.add_class::<PyTiEstimator>()?;
    module.add_class::<PyTiFit>()?;
    module.add_class::<PyBarEstimator>()?;
    module.add_class::<PyBarFit>()?;
    module.add_class::<PyMbarEstimator>()?;
    module.add_class::<PyMbarFit>()?;
    module.add_class::<PyIexpEstimator>()?;
    module.add_class::<PyIexpFit>()?;
    module.add_class::<PyNesEstimator>()?;
    module.add_class::<PyNesFit>()?;
    parent.add_submodule(&module)?;
    Ok(())
}

fn clone_dhdl_series(py: Python<'_>, values: Vec<Py<PyDhdlSeries>>) -> Vec<alchemrs::DhdlSeries> {
    values
        .into_iter()
        .map(|value| value.borrow(py).inner.clone())
        .collect()
}

fn clone_unk_windows(py: Python<'_>, values: Vec<Py<PyUNkMatrix>>) -> Vec<alchemrs::UNkMatrix> {
    values
        .into_iter()
        .map(|value| value.borrow(py).inner.clone())
        .collect()
}

fn clone_trajectories(
    py: Python<'_>,
    values: Vec<Py<PySwitchingTrajectory>>,
) -> Vec<alchemrs::SwitchingTrajectory> {
    values
        .into_iter()
        .map(|value| value.borrow(py).inner.clone())
        .collect()
}

pub(crate) fn parse_integration_method(method: &str) -> PyResult<IntegrationMethod> {
    match method {
        "trapezoidal" => Ok(IntegrationMethod::Trapezoidal),
        "simpson" => Ok(IntegrationMethod::Simpson),
        "cubic-spline" => Ok(IntegrationMethod::CubicSpline),
        "pchip" => Ok(IntegrationMethod::Pchip),
        "akima" => Ok(IntegrationMethod::Akima),
        "gaussian-quadrature" => Ok(IntegrationMethod::GaussianQuadrature),
        _ => Err(PyValueError::new_err(
            "method must be one of: 'trapezoidal', 'simpson', 'cubic-spline', 'pchip', 'akima', 'gaussian-quadrature'",
        )),
    }
}

fn parse_bar_method(method: &str) -> PyResult<BarMethod> {
    match method {
        "false-position" => Ok(BarMethod::FalsePosition),
        "self-consistent-iteration" => Ok(BarMethod::SelfConsistentIteration),
        "bisection" => Ok(BarMethod::Bisection),
        _ => Err(PyValueError::new_err(
            "method must be one of: 'false-position', 'self-consistent-iteration', 'bisection'",
        )),
    }
}

fn parse_mbar_solver(solver: &str) -> PyResult<MbarSolver> {
    match solver {
        "fixed-point" => Ok(MbarSolver::FixedPoint),
        "lbfgs" => Ok(MbarSolver::Lbfgs),
        _ => Err(PyValueError::new_err(
            "solver must be one of: 'fixed-point', 'lbfgs'",
        )),
    }
}
