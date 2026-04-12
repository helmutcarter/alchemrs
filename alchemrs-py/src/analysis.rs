use alchemrs::{analysis, BarOptions, IexpOptions, MbarOptions, MbarSolver, NesOptions, TiOptions};
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::data::{
    wrap_convergence_point, wrap_overlap_matrix, PyConvergencePoint, PyDhdlSeries, PyOverlapMatrix,
    PySwitchingTrajectory, PyUNkMatrix,
};
use crate::error::to_py_err;
use crate::estimators::parse_integration_method;

#[pyfunction]
#[pyo3(signature = (windows, max_iterations=10_000, tolerance=1e-7, initial_f_k=None, solver="fixed-point", parallel=false))]
fn overlap_matrix(
    py: Python<'_>,
    windows: Vec<Py<PyUNkMatrix>>,
    max_iterations: usize,
    tolerance: f64,
    initial_f_k: Option<Vec<f64>>,
    solver: &str,
    parallel: bool,
) -> PyResult<PyOverlapMatrix> {
    let windows = clone_unk_windows(py, windows);
    analysis::overlap_matrix(
        &windows,
        Some(MbarOptions {
            max_iterations,
            tolerance,
            initial_f_k,
            parallel,
            solver: parse_mbar_solver(solver)?,
        }),
    )
    .map(wrap_overlap_matrix)
    .map_err(to_py_err)
}

#[pyfunction]
fn overlap_eigenvalues<'py>(
    py: Python<'py>,
    overlap: &PyOverlapMatrix,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    analysis::overlap_eigenvalues(&overlap.inner)
        .map(|values| values.into_pyarray(py))
        .map_err(to_py_err)
}

#[pyfunction]
fn overlap_scalar(overlap: &PyOverlapMatrix) -> PyResult<f64> {
    analysis::overlap_scalar(&overlap.inner).map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (series, method="trapezoidal", parallel=false))]
fn ti_convergence(
    py: Python<'_>,
    series: Vec<Py<PyDhdlSeries>>,
    method: &str,
    parallel: bool,
) -> PyResult<Vec<PyConvergencePoint>> {
    let series = clone_dhdl_series(py, series);
    analysis::ti_convergence(
        &series,
        Some(TiOptions {
            method: parse_integration_method(method)?,
            parallel,
        }),
    )
    .map(|values| values.into_iter().map(wrap_convergence_point).collect())
    .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (windows, method="false-position", maximum_iterations=10_000, relative_tolerance=1e-7, parallel=false))]
fn bar_convergence(
    py: Python<'_>,
    windows: Vec<Py<PyUNkMatrix>>,
    method: &str,
    maximum_iterations: usize,
    relative_tolerance: f64,
    parallel: bool,
) -> PyResult<Vec<PyConvergencePoint>> {
    let windows = clone_unk_windows(py, windows);
    analysis::bar_convergence(
        &windows,
        Some(BarOptions {
            maximum_iterations,
            relative_tolerance,
            method: parse_bar_method(method)?,
            parallel,
        }),
    )
    .map(|values| values.into_iter().map(wrap_convergence_point).collect())
    .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (windows, max_iterations=10_000, tolerance=1e-7, initial_f_k=None, solver="fixed-point", parallel=false))]
fn mbar_convergence(
    py: Python<'_>,
    windows: Vec<Py<PyUNkMatrix>>,
    max_iterations: usize,
    tolerance: f64,
    initial_f_k: Option<Vec<f64>>,
    solver: &str,
    parallel: bool,
) -> PyResult<Vec<PyConvergencePoint>> {
    let windows = clone_unk_windows(py, windows);
    analysis::mbar_convergence(
        &windows,
        Some(MbarOptions {
            max_iterations,
            tolerance,
            initial_f_k,
            parallel,
            solver: parse_mbar_solver(solver)?,
        }),
    )
    .map(|values| values.into_iter().map(wrap_convergence_point).collect())
    .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (windows, parallel=false))]
fn exp_convergence(
    py: Python<'_>,
    windows: Vec<Py<PyUNkMatrix>>,
    parallel: bool,
) -> PyResult<Vec<PyConvergencePoint>> {
    let windows = clone_unk_windows(py, windows);
    analysis::exp_convergence(&windows, Some(IexpOptions { parallel }))
        .map(|values| values.into_iter().map(wrap_convergence_point).collect())
        .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (windows, parallel=false))]
fn dexp_convergence(
    py: Python<'_>,
    windows: Vec<Py<PyUNkMatrix>>,
    parallel: bool,
) -> PyResult<Vec<PyConvergencePoint>> {
    let windows = clone_unk_windows(py, windows);
    analysis::dexp_convergence(&windows, Some(IexpOptions { parallel }))
        .map(|values| values.into_iter().map(wrap_convergence_point).collect())
        .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (trajectories, n_bootstrap=0, seed=0))]
fn nes_convergence(
    py: Python<'_>,
    trajectories: Vec<Py<PySwitchingTrajectory>>,
    n_bootstrap: usize,
    seed: u64,
) -> PyResult<Vec<PyConvergencePoint>> {
    let trajectories = clone_trajectories(py, trajectories);
    analysis::nes_convergence(&trajectories, Some(NesOptions { n_bootstrap, seed }))
        .map(|values| values.into_iter().map(wrap_convergence_point).collect())
        .map_err(to_py_err)
}

pub fn register(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(py, "analysis")?;
    module.add_function(wrap_pyfunction!(overlap_matrix, &module)?)?;
    module.add_function(wrap_pyfunction!(overlap_eigenvalues, &module)?)?;
    module.add_function(wrap_pyfunction!(overlap_scalar, &module)?)?;
    module.add_function(wrap_pyfunction!(ti_convergence, &module)?)?;
    module.add_function(wrap_pyfunction!(bar_convergence, &module)?)?;
    module.add_function(wrap_pyfunction!(mbar_convergence, &module)?)?;
    module.add_function(wrap_pyfunction!(exp_convergence, &module)?)?;
    module.add_function(wrap_pyfunction!(dexp_convergence, &module)?)?;
    module.add_function(wrap_pyfunction!(nes_convergence, &module)?)?;
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

fn parse_bar_method(method: &str) -> PyResult<alchemrs::BarMethod> {
    match method {
        "false-position" => Ok(alchemrs::BarMethod::FalsePosition),
        "self-consistent-iteration" => Ok(alchemrs::BarMethod::SelfConsistentIteration),
        "bisection" => Ok(alchemrs::BarMethod::Bisection),
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
