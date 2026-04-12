use alchemrs::CoreError;
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(alchemrs, AlchemrsError, PyException);
create_exception!(alchemrs, InvalidShapeError, AlchemrsError);
create_exception!(alchemrs, InvalidStateError, AlchemrsError);
create_exception!(alchemrs, InvalidTimeOrderError, AlchemrsError);
create_exception!(alchemrs, NonFiniteValueError, AlchemrsError);
create_exception!(alchemrs, ParseError, AlchemrsError);
create_exception!(alchemrs, ConvergenceError, AlchemrsError);
create_exception!(alchemrs, RequiresOneDimensionalLambdaError, AlchemrsError);
create_exception!(alchemrs, UnsupportedInputError, AlchemrsError);

pub fn register(py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("AlchemrsError", py.get_type::<AlchemrsError>())?;
    module.add("InvalidShapeError", py.get_type::<InvalidShapeError>())?;
    module.add("InvalidStateError", py.get_type::<InvalidStateError>())?;
    module.add(
        "InvalidTimeOrderError",
        py.get_type::<InvalidTimeOrderError>(),
    )?;
    module.add("NonFiniteValueError", py.get_type::<NonFiniteValueError>())?;
    module.add("ParseError", py.get_type::<ParseError>())?;
    module.add("ConvergenceError", py.get_type::<ConvergenceError>())?;
    module.add(
        "RequiresOneDimensionalLambdaError",
        py.get_type::<RequiresOneDimensionalLambdaError>(),
    )?;
    module.add(
        "UnsupportedInputError",
        py.get_type::<UnsupportedInputError>(),
    )?;
    Ok(())
}

pub fn to_py_err(err: CoreError) -> PyErr {
    let message = err.to_string();
    match err {
        CoreError::InvalidShape { .. } => InvalidShapeError::new_err(message),
        CoreError::InvalidState(_) => InvalidStateError::new_err(message),
        CoreError::InvalidTimeOrder(_) => InvalidTimeOrderError::new_err(message),
        CoreError::NonFiniteValue(_) => NonFiniteValueError::new_err(message),
        CoreError::Parse(_) => ParseError::new_err(message),
        CoreError::ConvergenceFailure => ConvergenceError::new_err(message),
        CoreError::RequiresOneDimensionalLambda { .. } => {
            RequiresOneDimensionalLambdaError::new_err(message)
        }
        CoreError::Unsupported(_) => UnsupportedInputError::new_err(message),
    }
}
