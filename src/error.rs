use thiserror::Error;

#[derive(Error, Debug)]
pub enum CoreError {
    #[error("invalid shape: expected {expected}, found {found}")]
    InvalidShape { expected: usize, found: usize },
    #[error("invalid state metadata: {0}")]
    InvalidState(String),
    #[error("invalid time ordering: {0}")]
    InvalidTimeOrder(String),
    #[error("non-finite numerical value encountered: {0}")]
    NonFiniteValue(String),
    #[error("parser error: {0}")]
    Parse(String),
    #[error("estimator failed to converge")]
    ConvergenceFailure,
    #[error("unsupported input: {0}")]
    Unsupported(String),
}

pub type Result<T> = std::result::Result<T, CoreError>;

pub(crate) fn ensure_finite(label: &str, values: &[f64]) -> Result<()> {
    for (idx, value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(CoreError::NonFiniteValue(format!(
                "{label}[{idx}] must be finite"
            )));
        }
    }
    Ok(())
}

pub(crate) fn ensure_finite_or_positive_infinity(label: &str, values: &[f64]) -> Result<()> {
    for (idx, value) in values.iter().enumerate() {
        if value.is_nan() || *value == f64::NEG_INFINITY {
            return Err(CoreError::NonFiniteValue(format!(
                "{label}[{idx}] must be finite or +inf"
            )));
        }
    }
    Ok(())
}

pub(crate) fn ensure_finite_or_nan(label: &str, values: &[f64]) -> Result<()> {
    for (idx, value) in values.iter().enumerate() {
        if value.is_infinite() {
            return Err(CoreError::NonFiniteValue(format!(
                "{label}[{idx}] must be finite or NaN"
            )));
        }
    }
    Ok(())
}
