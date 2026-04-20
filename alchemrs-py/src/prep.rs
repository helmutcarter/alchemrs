use alchemrs::prep::{self, DecorrelationOptions, UNkSeriesMethod};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::data::{
    extract_vec1_any, wrap_dhdl_series, wrap_equilibration_result, wrap_unk_matrix, PyDhdlSeries,
    PyEquilibrationResult, PyUNkMatrix,
};
use crate::error::to_py_err;

#[pyclass(name = "DecorrelationOptions", module = "alchemrs._alchemrs.prep")]
#[derive(Clone)]
pub struct PyDecorrelationOptions {
    pub(crate) inner: DecorrelationOptions,
}

#[pymethods]
impl PyDecorrelationOptions {
    #[new]
    #[pyo3(signature = (
        drop_duplicates=true,
        sort=true,
        conservative=true,
        remove_burnin=false,
        fast=false,
        nskip=1,
        lower=None,
        upper=None,
        step=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        drop_duplicates: bool,
        sort: bool,
        conservative: bool,
        remove_burnin: bool,
        fast: bool,
        nskip: usize,
        lower: Option<f64>,
        upper: Option<f64>,
        step: Option<usize>,
    ) -> Self {
        Self {
            inner: DecorrelationOptions {
                drop_duplicates,
                sort,
                conservative,
                remove_burnin,
                fast,
                nskip,
                lower,
                upper,
                step,
            },
        }
    }

    #[getter]
    fn drop_duplicates(&self) -> bool {
        self.inner.drop_duplicates
    }

    #[getter]
    fn sort(&self) -> bool {
        self.inner.sort
    }

    #[getter]
    fn conservative(&self) -> bool {
        self.inner.conservative
    }

    #[getter]
    fn remove_burnin(&self) -> bool {
        self.inner.remove_burnin
    }

    #[getter]
    fn fast(&self) -> bool {
        self.inner.fast
    }

    #[getter]
    fn nskip(&self) -> usize {
        self.inner.nskip
    }

    #[getter]
    fn lower(&self) -> Option<f64> {
        self.inner.lower
    }

    #[getter]
    fn upper(&self) -> Option<f64> {
        self.inner.upper
    }

    #[getter]
    fn step(&self) -> Option<usize> {
        self.inner.step
    }
}

#[pyfunction]
#[pyo3(signature = (series, options=None))]
fn decorrelate_dhdl(
    series: &PyDhdlSeries,
    options: Option<PyRef<'_, PyDecorrelationOptions>>,
) -> PyResult<PyDhdlSeries> {
    let options = options_or_default(options);
    prep::decorrelate_dhdl(&series.inner, &options)
        .map(wrap_dhdl_series)
        .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (series, options=None))]
fn detect_equilibration_dhdl(
    series: &PyDhdlSeries,
    options: Option<PyRef<'_, PyDecorrelationOptions>>,
) -> PyResult<PyEquilibrationResult> {
    let options = options_or_default(options);
    prep::detect_equilibration_dhdl(&series.inner, &options)
        .map(wrap_equilibration_result)
        .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (u_nk, method="de", options=None))]
fn decorrelate_u_nk(
    u_nk: &PyUNkMatrix,
    method: &str,
    options: Option<PyRef<'_, PyDecorrelationOptions>>,
) -> PyResult<PyUNkMatrix> {
    let options = options_or_default(options);
    prep::decorrelate_u_nk(&u_nk.inner, parse_method(method)?, &options)
        .map(wrap_unk_matrix)
        .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (u_nk, observable, options=None))]
fn decorrelate_u_nk_with_observable(
    u_nk: &PyUNkMatrix,
    observable: &Bound<'_, PyAny>,
    options: Option<PyRef<'_, PyDecorrelationOptions>>,
) -> PyResult<PyUNkMatrix> {
    let observable = extract_vec1_any(observable)?;
    let options = options_or_default(options);
    prep::decorrelate_u_nk_with_observable(&u_nk.inner, &observable, &options)
        .map(wrap_unk_matrix)
        .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (u_nk, method="de", options=None))]
fn detect_equilibration_u_nk(
    u_nk: &PyUNkMatrix,
    method: &str,
    options: Option<PyRef<'_, PyDecorrelationOptions>>,
) -> PyResult<PyEquilibrationResult> {
    let options = options_or_default(options);
    prep::detect_equilibration_u_nk(&u_nk.inner, parse_method(method)?, &options)
        .map(wrap_equilibration_result)
        .map_err(to_py_err)
}

#[pyfunction]
#[pyo3(signature = (observable, options=None))]
fn detect_equilibration_observable(
    observable: &Bound<'_, PyAny>,
    options: Option<PyRef<'_, PyDecorrelationOptions>>,
) -> PyResult<PyEquilibrationResult> {
    let observable = extract_vec1_any(observable)?;
    let options = options.map(|value| value.inner.clone()).unwrap_or_default();
    prep::detect_equilibration_observable(&observable, options.fast, options.nskip)
        .map(wrap_equilibration_result)
        .map_err(to_py_err)
}

pub fn register(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(py, "prep")?;
    module.add_class::<PyDecorrelationOptions>()?;
    module.add_function(wrap_pyfunction!(decorrelate_dhdl, &module)?)?;
    module.add_function(wrap_pyfunction!(detect_equilibration_dhdl, &module)?)?;
    module.add_function(wrap_pyfunction!(decorrelate_u_nk, &module)?)?;
    module.add_function(wrap_pyfunction!(decorrelate_u_nk_with_observable, &module)?)?;
    module.add_function(wrap_pyfunction!(detect_equilibration_u_nk, &module)?)?;
    module.add_function(wrap_pyfunction!(detect_equilibration_observable, &module)?)?;
    parent.add_submodule(&module)?;
    Ok(())
}

fn options_or_default(options: Option<PyRef<'_, PyDecorrelationOptions>>) -> DecorrelationOptions {
    options.map(|value| value.inner.clone()).unwrap_or_default()
}

fn parse_method(method: &str) -> PyResult<UNkSeriesMethod> {
    match method {
        "de" => Ok(UNkSeriesMethod::DE),
        "all" => Ok(UNkSeriesMethod::All),
        _ => Err(PyValueError::new_err("method must be one of: 'de', 'all'")),
    }
}
