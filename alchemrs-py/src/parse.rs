use std::path::PathBuf;

use alchemrs::parse;
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::data::{wrap_dhdl_series, wrap_switching_trajectory, wrap_unk_matrix};
use crate::error::to_py_err;

#[pyfunction]
fn infer_temperature(path: PathBuf) -> PyResult<f64> {
    parse::infer_temperature(path).map_err(to_py_err)
}

#[pyfunction]
fn extract_dhdl(path: PathBuf, temperature_k: f64) -> PyResult<crate::data::PyDhdlSeries> {
    parse::extract_dhdl(path, temperature_k)
        .map(wrap_dhdl_series)
        .map_err(to_py_err)
}

#[pyfunction]
fn extract_u_nk(path: PathBuf, temperature_k: f64) -> PyResult<crate::data::PyUNkMatrix> {
    parse::extract_u_nk(path, temperature_k)
        .map(wrap_unk_matrix)
        .map_err(to_py_err)
}

#[pyfunction]
fn extract_u_nk_with_potential<'py>(
    py: Python<'py>,
    path: PathBuf,
    temperature_k: f64,
) -> PyResult<(crate::data::PyUNkMatrix, Bound<'py, PyArray1<f64>>)> {
    parse::extract_u_nk_with_potential(path, temperature_k)
        .map(|(u_nk, potential)| (wrap_unk_matrix(u_nk), potential.into_pyarray(py)))
        .map_err(to_py_err)
}

#[pyfunction]
fn extract_nes_trajectory(
    path: PathBuf,
    temperature_k: f64,
) -> PyResult<crate::data::PySwitchingTrajectory> {
    parse::extract_nes_trajectory(path, temperature_k)
        .map(wrap_switching_trajectory)
        .map_err(to_py_err)
}

pub fn register(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(py, "parse")?;
    module.add_function(wrap_pyfunction!(infer_temperature, &module)?)?;
    module.add_function(wrap_pyfunction!(extract_dhdl, &module)?)?;
    module.add_function(wrap_pyfunction!(extract_u_nk, &module)?)?;
    module.add_function(wrap_pyfunction!(extract_u_nk_with_potential, &module)?)?;
    module.add_function(wrap_pyfunction!(extract_nes_trajectory, &module)?)?;
    parent.add_submodule(&module)?;
    Ok(())
}
