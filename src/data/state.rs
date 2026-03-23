use crate::error::{CoreError, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct StatePoint {
    lambdas: Vec<f64>,
    temperature_k: f64,
}

impl StatePoint {
    pub fn new(lambdas: Vec<f64>, temperature_k: f64) -> Result<Self> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(CoreError::InvalidState(
                "temperature must be finite and positive".to_string(),
            ));
        }
        for (idx, value) in lambdas.iter().enumerate() {
            if !value.is_finite() {
                return Err(CoreError::InvalidState(format!(
                    "lambda[{idx}] must be finite"
                )));
            }
        }
        Ok(Self {
            lambdas,
            temperature_k,
        })
    }

    pub fn lambdas(&self) -> &[f64] {
        &self.lambdas
    }

    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
}

#[cfg(test)]
mod tests {
    use super::StatePoint;
    use crate::error::CoreError;

    #[test]
    fn state_point_rejects_invalid_temperature() {
        let err = StatePoint::new(vec![0.0], 0.0).unwrap_err();
        assert!(matches!(err, CoreError::InvalidState(_)));
    }

    #[test]
    fn state_point_rejects_nonfinite_lambda() {
        let err = StatePoint::new(vec![f64::NAN], 300.0).unwrap_err();
        assert!(matches!(err, CoreError::InvalidState(_)));
    }
}
