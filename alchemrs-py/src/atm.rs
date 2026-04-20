use alchemrs::{
    AtmBindingEstimate, AtmDirection, AtmEstimator, AtmFit, AtmLogQMatrix, AtmOptions, AtmSample,
    AtmSampleSet, AtmSchedule, AtmState, AtmUncertaintyMethod, UwhamOptions,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::data::{
    extract_vec1_any, wrap_delta_f_matrix, wrap_free_energy_estimate, PyDeltaFMatrix,
    PyFreeEnergyEstimate, PyStatePoint,
};
use crate::error::to_py_err;

#[pyclass(name = "AtmState", module = "alchemrs._alchemrs.atm")]
#[derive(Clone)]
pub struct PyAtmState {
    pub(crate) inner: AtmState,
}

struct PyAtmStateInit {
    state_id: usize,
    direction: String,
    lambda1: f64,
    lambda2: f64,
    alpha: f64,
    u0: f64,
    w0: f64,
    temperature_k: f64,
    lambda3: Option<f64>,
    u1: Option<f64>,
}

#[pymethods]
impl PyAtmState {
    #[new]
    #[pyo3(signature = (state_id, direction, lambda1, lambda2, alpha, u0, w0, temperature_k, lambda3=None, u1=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        state_id: usize,
        direction: String,
        lambda1: f64,
        lambda2: f64,
        alpha: f64,
        u0: f64,
        w0: f64,
        temperature_k: f64,
        lambda3: Option<f64>,
        u1: Option<f64>,
    ) -> PyResult<Self> {
        let init = PyAtmStateInit {
            state_id,
            direction,
            lambda1,
            lambda2,
            alpha,
            u0,
            w0,
            temperature_k,
            lambda3,
            u1,
        };
        AtmState::new(
            init.state_id,
            parse_direction(&init.direction)?,
            init.lambda1,
            init.lambda2,
            init.lambda3,
            init.alpha,
            init.u0,
            init.u1,
            init.w0,
            init.temperature_k,
        )
        .map(|inner| Self { inner })
        .map_err(to_py_err)
    }

    #[getter]
    fn state_id(&self) -> usize {
        self.inner.state_id()
    }

    #[getter]
    fn direction(&self) -> &'static str {
        match self.inner.direction() {
            AtmDirection::Forward => "forward",
            AtmDirection::Reverse => "reverse",
        }
    }

    #[getter]
    fn lambda1(&self) -> f64 {
        self.inner.lambda1()
    }

    #[getter]
    fn lambda2(&self) -> f64 {
        self.inner.lambda2()
    }

    #[getter]
    fn lambda3(&self) -> Option<f64> {
        self.inner.lambda3()
    }

    #[getter]
    fn alpha(&self) -> f64 {
        self.inner.alpha()
    }

    #[getter]
    fn u0(&self) -> f64 {
        self.inner.u0()
    }

    #[getter]
    fn u1(&self) -> Option<f64> {
        self.inner.u1()
    }

    #[getter]
    fn w0(&self) -> f64 {
        self.inner.w0()
    }

    #[getter]
    fn temperature_k(&self) -> f64 {
        self.inner.temperature_k()
    }
}

#[pyclass(name = "AtmSchedule", module = "alchemrs._alchemrs.atm")]
#[derive(Clone)]
pub struct PyAtmSchedule {
    pub(crate) inner: AtmSchedule,
}

#[pymethods]
impl PyAtmSchedule {
    #[new]
    fn new(py: Python<'_>, states: Vec<Py<PyAtmState>>) -> PyResult<Self> {
        let states = states
            .into_iter()
            .map(|state| state.borrow(py).inner.clone())
            .collect::<Vec<_>>();
        AtmSchedule::new(states)
            .map(|inner| Self { inner })
            .map_err(to_py_err)
    }

    #[getter]
    fn states(&self) -> Vec<PyAtmState> {
        self.inner
            .states()
            .iter()
            .cloned()
            .map(|inner| PyAtmState { inner })
            .collect()
    }

    #[getter]
    fn direction(&self) -> &'static str {
        match self.inner.direction() {
            AtmDirection::Forward => "forward",
            AtmDirection::Reverse => "reverse",
        }
    }

    #[getter]
    fn temperature_k(&self) -> f64 {
        self.inner.temperature_k()
    }

    #[getter]
    fn uses_multi_softplus(&self) -> bool {
        self.inner.uses_multi_softplus()
    }
}

#[pyclass(name = "AtmSample", module = "alchemrs._alchemrs.atm")]
#[derive(Clone)]
pub struct PyAtmSample {
    pub(crate) inner: AtmSample,
}

#[pymethods]
impl PyAtmSample {
    #[new]
    fn new(
        state_id: usize,
        potential_energy_kcal_per_mol: f64,
        perturbation_energy_kcal_per_mol: f64,
    ) -> PyResult<Self> {
        AtmSample::new(
            state_id,
            potential_energy_kcal_per_mol,
            perturbation_energy_kcal_per_mol,
        )
        .map(|inner| Self { inner })
        .map_err(to_py_err)
    }

    #[getter]
    fn state_id(&self) -> usize {
        self.inner.state_id()
    }

    #[getter]
    fn potential_energy_kcal_per_mol(&self) -> f64 {
        self.inner.potential_energy_kcal_per_mol()
    }

    #[getter]
    fn perturbation_energy_kcal_per_mol(&self) -> f64 {
        self.inner.perturbation_energy_kcal_per_mol()
    }
}

#[pyclass(name = "AtmSampleSet", module = "alchemrs._alchemrs.atm")]
#[derive(Clone)]
pub struct PyAtmSampleSet {
    pub(crate) inner: AtmSampleSet,
}

#[pymethods]
impl PyAtmSampleSet {
    #[new]
    fn new(
        py: Python<'_>,
        schedule: Py<PyAtmSchedule>,
        samples: Vec<Py<PyAtmSample>>,
    ) -> PyResult<Self> {
        let schedule = schedule.borrow(py).inner.clone();
        let samples = samples
            .into_iter()
            .map(|sample| sample.borrow(py).inner.clone())
            .collect::<Vec<_>>();
        AtmSampleSet::new(schedule, samples)
            .map(|inner| Self { inner })
            .map_err(to_py_err)
    }

    #[getter]
    fn schedule(&self) -> PyAtmSchedule {
        PyAtmSchedule {
            inner: self.inner.schedule().clone(),
        }
    }

    #[getter]
    fn samples(&self) -> Vec<PyAtmSample> {
        self.inner
            .samples()
            .iter()
            .cloned()
            .map(|inner| PyAtmSample { inner })
            .collect()
    }

    fn to_log_q_matrix(&self) -> PyResult<PyAtmLogQMatrix> {
        self.inner
            .to_log_q_matrix()
            .map(wrap_atm_log_q_matrix)
            .map_err(to_py_err)
    }
}

#[pyclass(name = "AtmLogQMatrix", module = "alchemrs._alchemrs.atm")]
#[derive(Clone)]
pub struct PyAtmLogQMatrix {
    pub(crate) inner: AtmLogQMatrix,
}

#[pymethods]
impl PyAtmLogQMatrix {
    #[new]
    #[pyo3(signature = (n_observations, n_states, log_q, states, sampled_counts, lambda_labels=None))]
    fn new(
        py: Python<'_>,
        n_observations: usize,
        n_states: usize,
        log_q: Vec<f64>,
        states: Vec<Py<PyStatePoint>>,
        sampled_counts: Vec<usize>,
        lambda_labels: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let states = states
            .into_iter()
            .map(|state| state.borrow(py).inner.clone())
            .collect::<Vec<_>>();
        AtmLogQMatrix::new_with_labels(
            n_observations,
            n_states,
            log_q,
            states,
            sampled_counts,
            lambda_labels,
        )
        .map(|inner| Self { inner })
        .map_err(to_py_err)
    }

    #[getter]
    fn n_observations(&self) -> usize {
        self.inner.n_observations()
    }

    #[getter]
    fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    #[getter]
    fn log_q(&self) -> Vec<f64> {
        self.inner.log_q().to_vec()
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
    fn sampled_counts(&self) -> Vec<usize> {
        self.inner.sampled_counts().to_vec()
    }

    #[getter]
    fn lambda_labels(&self) -> Option<Vec<String>> {
        self.inner.lambda_labels().map(|labels| labels.to_vec())
    }
}

#[pyclass(name = "ATM", module = "alchemrs._alchemrs.atm")]
#[derive(Clone)]
pub struct PyAtmEstimator {
    inner: AtmEstimator,
}

#[pymethods]
impl PyAtmEstimator {
    #[new]
    #[pyo3(signature = (max_iterations=10_000, tolerance=1e-7, parallel=false, uncertainty="analytical", n_bootstrap=256, seed=0))]
    fn new(
        max_iterations: usize,
        tolerance: f64,
        parallel: bool,
        uncertainty: &str,
        n_bootstrap: usize,
        seed: u64,
    ) -> PyResult<Self> {
        let uncertainty = parse_uncertainty_method(uncertainty, n_bootstrap, seed)?;
        Ok(Self {
            inner: AtmEstimator::new(AtmOptions {
                uwham: UwhamOptions {
                    max_iterations,
                    tolerance,
                    parallel,
                },
                uncertainty,
            }),
        })
    }

    fn fit(&self, py: Python<'_>, data: Py<PyAtmLogQMatrix>) -> PyResult<PyAtmFit> {
        self.inner
            .fit(&data.borrow(py).inner)
            .map(|inner| PyAtmFit { inner })
            .map_err(to_py_err)
    }

    fn fit_leg(&self, py: Python<'_>, samples: Py<PyAtmSampleSet>) -> PyResult<PyAtmFit> {
        self.inner
            .fit_leg(&samples.borrow(py).inner)
            .map(|inner| PyAtmFit { inner })
            .map_err(to_py_err)
    }

    fn estimate(&self, py: Python<'_>, data: Py<PyAtmLogQMatrix>) -> PyResult<PyDeltaFMatrix> {
        self.inner
            .estimate(&data.borrow(py).inner)
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn estimate_leg(
        &self,
        py: Python<'_>,
        samples: Py<PyAtmSampleSet>,
    ) -> PyResult<PyFreeEnergyEstimate> {
        self.inner
            .estimate_leg(&samples.borrow(py).inner)
            .map(wrap_free_energy_estimate)
            .map_err(to_py_err)
    }

    fn estimate_binding(
        &self,
        py: Python<'_>,
        leg1: Py<PyAtmSampleSet>,
        leg2: Py<PyAtmSampleSet>,
    ) -> PyResult<PyAtmBindingEstimate> {
        self.inner
            .estimate_binding(&leg1.borrow(py).inner, &leg2.borrow(py).inner)
            .map(wrap_atm_binding_estimate)
            .map_err(to_py_err)
    }

    fn estimate_rbfe(
        &self,
        py: Python<'_>,
        leg1: Py<PyAtmSampleSet>,
        leg2: Py<PyAtmSampleSet>,
    ) -> PyResult<PyAtmBindingEstimate> {
        self.inner
            .estimate_rbfe(&leg1.borrow(py).inner, &leg2.borrow(py).inner)
            .map(wrap_atm_binding_estimate)
            .map_err(to_py_err)
    }

    fn estimate_abfe(
        &self,
        py: Python<'_>,
        leg1: Py<PyAtmSampleSet>,
        leg2: Py<PyAtmSampleSet>,
    ) -> PyResult<PyAtmBindingEstimate> {
        self.inner
            .estimate_abfe(&leg1.borrow(py).inner, &leg2.borrow(py).inner)
            .map(wrap_atm_binding_estimate)
            .map_err(to_py_err)
    }
}

#[pyclass(name = "AtmFit", module = "alchemrs._alchemrs.atm")]
pub struct PyAtmFit {
    inner: AtmFit,
}

#[pymethods]
impl PyAtmFit {
    #[getter]
    fn n_observations(&self) -> usize {
        self.inner.n_observations()
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

    #[getter]
    fn free_energies(&self) -> Vec<f64> {
        self.inner.free_energies().to_vec()
    }

    #[getter]
    fn weights(&self) -> Vec<f64> {
        self.inner.weights().to_vec()
    }

    #[getter]
    fn leg_uncertainty(&self) -> Option<f64> {
        self.inner.leg_uncertainty()
    }

    fn result(&self) -> PyResult<PyDeltaFMatrix> {
        self.inner
            .result()
            .map(wrap_delta_f_matrix)
            .map_err(to_py_err)
    }

    fn leg_result(&self) -> PyResult<PyFreeEnergyEstimate> {
        self.inner
            .leg_result()
            .map(wrap_free_energy_estimate)
            .map_err(to_py_err)
    }
}

#[pyclass(name = "AtmBindingEstimate", module = "alchemrs._alchemrs.atm")]
#[derive(Clone)]
pub struct PyAtmBindingEstimate {
    inner: AtmBindingEstimate,
}

#[pymethods]
impl PyAtmBindingEstimate {
    #[getter]
    fn delta_f(&self) -> f64 {
        self.inner.delta_f()
    }

    #[getter]
    fn uncertainty(&self) -> Option<f64> {
        self.inner.uncertainty()
    }

    #[getter]
    fn leg1(&self) -> PyFreeEnergyEstimate {
        wrap_free_energy_estimate(self.inner.leg1().clone())
    }

    #[getter]
    fn leg2(&self) -> PyFreeEnergyEstimate {
        wrap_free_energy_estimate(self.inner.leg2().clone())
    }
}

#[pyfunction]
#[pyo3(signature = (state_ids, direction, lambda1, lambda2, alpha, u0, w0, temperature_k, lambda3=None, u1=None))]
#[allow(clippy::too_many_arguments)]
fn schedule_from_arrays(
    state_ids: Vec<usize>,
    direction: String,
    lambda1: Vec<f64>,
    lambda2: Vec<f64>,
    alpha: Vec<f64>,
    u0: Vec<f64>,
    w0: Vec<f64>,
    temperature_k: Vec<f64>,
    lambda3: Option<Vec<f64>>,
    u1: Option<Vec<f64>>,
) -> PyResult<PyAtmSchedule> {
    let n_states = state_ids.len();
    for (label, len) in [
        ("lambda1", lambda1.len()),
        ("lambda2", lambda2.len()),
        ("alpha", alpha.len()),
        ("u0", u0.len()),
        ("w0", w0.len()),
        ("temperature_k", temperature_k.len()),
    ] {
        if len != n_states {
            return Err(PyValueError::new_err(format!(
                "{label} must have length {n_states}, found {len}"
            )));
        }
    }
    if let Some(values) = lambda3.as_ref() {
        if values.len() != n_states {
            return Err(PyValueError::new_err(format!(
                "lambda3 must have length {n_states}, found {}",
                values.len()
            )));
        }
    }
    if let Some(values) = u1.as_ref() {
        if values.len() != n_states {
            return Err(PyValueError::new_err(format!(
                "u1 must have length {n_states}, found {}",
                values.len()
            )));
        }
    }

    let direction = parse_direction(&direction)?;
    let states = (0..n_states)
        .map(|idx| {
            AtmState::new(
                state_ids[idx],
                direction,
                lambda1[idx],
                lambda2[idx],
                lambda3.as_ref().map(|values| values[idx]),
                alpha[idx],
                u0[idx],
                u1.as_ref().map(|values| values[idx]),
                w0[idx],
                temperature_k[idx],
            )
        })
        .collect::<alchemrs::Result<Vec<_>>>()
        .map_err(to_py_err)?;

    AtmSchedule::new(states)
        .map(|inner| PyAtmSchedule { inner })
        .map_err(to_py_err)
}

#[pyfunction]
fn sample_set_from_arrays(
    py: Python<'_>,
    schedule: Py<PyAtmSchedule>,
    state_ids: Vec<usize>,
    potential_energies_kcal_per_mol: &Bound<'_, PyAny>,
    perturbation_energies_kcal_per_mol: &Bound<'_, PyAny>,
) -> PyResult<PyAtmSampleSet> {
    let schedule = schedule.borrow(py).inner.clone();
    let potential = extract_vec1_any(potential_energies_kcal_per_mol)?;
    let perturbation = extract_vec1_any(perturbation_energies_kcal_per_mol)?;
    let n_samples = state_ids.len();
    if potential.len() != n_samples || perturbation.len() != n_samples {
        return Err(PyValueError::new_err(format!(
            "state_ids, potential_energies_kcal_per_mol, and perturbation_energies_kcal_per_mol must all have the same length; found {n_samples}, {}, {}",
            potential.len(),
            perturbation.len()
        )));
    }

    let samples = state_ids
        .into_iter()
        .zip(potential)
        .zip(perturbation)
        .map(|((state_id, pot_e), pert_e)| AtmSample::new(state_id, pot_e, pert_e))
        .collect::<alchemrs::Result<Vec<_>>>()
        .map_err(to_py_err)?;

    AtmSampleSet::new(schedule, samples)
        .map(|inner| PyAtmSampleSet { inner })
        .map_err(to_py_err)
}

pub fn register(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(py, "atm")?;
    module.add_class::<PyAtmState>()?;
    module.add_class::<PyAtmSchedule>()?;
    module.add_class::<PyAtmSample>()?;
    module.add_class::<PyAtmSampleSet>()?;
    module.add_class::<PyAtmLogQMatrix>()?;
    module.add_class::<PyAtmEstimator>()?;
    module.add_class::<PyAtmFit>()?;
    module.add_class::<PyAtmBindingEstimate>()?;
    module.add_function(wrap_pyfunction!(schedule_from_arrays, &module)?)?;
    module.add_function(wrap_pyfunction!(sample_set_from_arrays, &module)?)?;
    parent.add_submodule(&module)?;
    Ok(())
}

pub(crate) fn wrap_atm_log_q_matrix(inner: AtmLogQMatrix) -> PyAtmLogQMatrix {
    PyAtmLogQMatrix { inner }
}

fn wrap_atm_binding_estimate(inner: AtmBindingEstimate) -> PyAtmBindingEstimate {
    PyAtmBindingEstimate { inner }
}

fn parse_direction(direction: &str) -> PyResult<AtmDirection> {
    match direction {
        "forward" | "leg1" => Ok(AtmDirection::Forward),
        "reverse" | "leg2" => Ok(AtmDirection::Reverse),
        _ => Err(PyValueError::new_err(
            "direction must be one of: 'forward', 'reverse', 'leg1', 'leg2'",
        )),
    }
}

fn parse_uncertainty_method(
    uncertainty: &str,
    n_bootstrap: usize,
    seed: u64,
) -> PyResult<AtmUncertaintyMethod> {
    match uncertainty {
        "analytical" => Ok(AtmUncertaintyMethod::Analytical),
        "bootstrap" => {
            if n_bootstrap == 0 {
                Err(PyValueError::new_err(
                    "n_bootstrap must be positive when uncertainty='bootstrap'",
                ))
            } else {
                Ok(AtmUncertaintyMethod::Bootstrap { n_bootstrap, seed })
            }
        }
        "none" => Ok(AtmUncertaintyMethod::None),
        _ => Err(PyValueError::new_err(
            "uncertainty must be one of: 'analytical', 'bootstrap', 'none'",
        )),
    }
}
