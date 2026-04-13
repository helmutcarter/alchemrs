mod analysis;
mod atm;
mod data;
mod error;
mod estimators;
mod parse;
mod prep;

use pyo3::prelude::*;

#[pymodule]
fn _alchemrs(py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    error::register(py, module)?;
    data::register(module)?;
    atm::register(py, module)?;
    analysis::register(py, module)?;
    estimators::register(py, module)?;
    parse::register(py, module)?;
    prep::register(py, module)?;
    Ok(())
}
