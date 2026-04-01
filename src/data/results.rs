use crate::data::StatePoint;
use crate::error::{ensure_finite, ensure_finite_or_nan, CoreError, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct FreeEnergyEstimate {
    delta_f: f64,
    uncertainty: Option<f64>,
    from_state: StatePoint,
    to_state: StatePoint,
}

impl FreeEnergyEstimate {
    pub fn new(
        delta_f: f64,
        uncertainty: Option<f64>,
        from_state: StatePoint,
        to_state: StatePoint,
    ) -> Result<Self> {
        if !delta_f.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "delta_f must be finite".to_string(),
            ));
        }
        if let Some(value) = uncertainty {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "uncertainty must be finite".to_string(),
                ));
            }
        }
        Ok(Self {
            delta_f,
            uncertainty,
            from_state,
            to_state,
        })
    }

    pub fn delta_f(&self) -> f64 {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DeltaFMatrix {
    values: Vec<f64>,
    uncertainties: Option<Vec<f64>>,
    n_states: usize,
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
}

impl DeltaFMatrix {
    pub fn new(
        values: Vec<f64>,
        uncertainties: Option<Vec<f64>>,
        n_states: usize,
        states: Vec<StatePoint>,
    ) -> Result<Self> {
        Self::new_with_labels(values, uncertainties, n_states, states, None)
    }

    pub fn new_with_labels(
        values: Vec<f64>,
        uncertainties: Option<Vec<f64>>,
        n_states: usize,
        states: Vec<StatePoint>,
        lambda_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        let expected = n_states
            .checked_mul(n_states)
            .ok_or(CoreError::InvalidShape {
                expected: n_states,
                found: n_states,
            })?;
        if values.len() != expected {
            return Err(CoreError::InvalidShape {
                expected,
                found: values.len(),
            });
        }
        if let Some(ref unc) = uncertainties {
            if unc.len() != expected {
                return Err(CoreError::InvalidShape {
                    expected,
                    found: unc.len(),
                });
            }
            ensure_finite_or_nan("delta_f_uncertainty", unc)?;
        }
        if states.len() != n_states {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: states.len(),
            });
        }
        if let Some(labels) = lambda_labels.as_ref() {
            let expected = states
                .first()
                .map(|state| state.lambdas().len())
                .unwrap_or(0);
            if labels.len() != expected {
                return Err(CoreError::InvalidShape {
                    expected,
                    found: labels.len(),
                });
            }
        }
        ensure_finite("delta_f", &values)?;
        Ok(Self {
            values,
            uncertainties,
            n_states,
            states,
            lambda_labels,
        })
    }

    pub fn n_states(&self) -> usize {
        self.n_states
    }

    pub fn values(&self) -> &[f64] {
        &self.values
    }

    pub fn uncertainties(&self) -> Option<&[f64]> {
        self.uncertainties.as_deref()
    }

    pub fn states(&self) -> &[StatePoint] {
        &self.states
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct OverlapMatrix {
    values: Vec<f64>,
    n_states: usize,
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
}

impl OverlapMatrix {
    pub fn new(values: Vec<f64>, n_states: usize, states: Vec<StatePoint>) -> Result<Self> {
        Self::new_with_labels(values, n_states, states, None)
    }

    pub fn new_with_labels(
        values: Vec<f64>,
        n_states: usize,
        states: Vec<StatePoint>,
        lambda_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        let expected = n_states
            .checked_mul(n_states)
            .ok_or(CoreError::InvalidShape {
                expected: n_states,
                found: n_states,
            })?;
        if values.len() != expected {
            return Err(CoreError::InvalidShape {
                expected,
                found: values.len(),
            });
        }
        if states.len() != n_states {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: states.len(),
            });
        }
        if let Some(labels) = lambda_labels.as_ref() {
            let expected = states
                .first()
                .map(|state| state.lambdas().len())
                .unwrap_or(0);
            if labels.len() != expected {
                return Err(CoreError::InvalidShape {
                    expected,
                    found: labels.len(),
                });
            }
        }
        ensure_finite("overlap_matrix", &values)?;
        Ok(Self {
            values,
            n_states,
            states,
            lambda_labels,
        })
    }

    pub fn n_states(&self) -> usize {
        self.n_states
    }

    pub fn values(&self) -> &[f64] {
        &self.values
    }

    pub fn states(&self) -> &[StatePoint] {
        &self.states
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::{DeltaFMatrix, FreeEnergyEstimate, OverlapMatrix};
    use crate::data::StatePoint;
    use crate::error::CoreError;

    #[test]
    fn delta_f_matrix_checks_lengths() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let err = DeltaFMatrix::new(vec![0.0, 1.0, 2.0], None, 2, vec![state]).unwrap_err();
        assert!(matches!(err, CoreError::InvalidShape { .. }));
    }

    #[test]
    fn free_energy_estimate_rejects_nonfinite() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let err = FreeEnergyEstimate::new(f64::INFINITY, None, state.clone(), state).unwrap_err();
        assert!(matches!(err, CoreError::NonFiniteValue(_)));
    }

    #[test]
    fn delta_f_matrix_preserves_lambda_labels() {
        let state0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
        let matrix = DeltaFMatrix::new_with_labels(
            vec![0.0, 1.0, -1.0, 0.0],
            None,
            2,
            vec![state0, state1],
            Some(vec!["coul-lambda".to_string(), "vdw-lambda".to_string()]),
        )
        .unwrap();
        assert_eq!(
            matrix.lambda_labels().unwrap(),
            &["coul-lambda", "vdw-lambda"]
        );
    }

    #[test]
    fn overlap_matrix_preserves_lambda_labels() {
        let state0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
        let matrix = OverlapMatrix::new_with_labels(
            vec![1.0, 0.0, 0.0, 1.0],
            2,
            vec![state0, state1],
            Some(vec!["coul-lambda".to_string(), "vdw-lambda".to_string()]),
        )
        .unwrap();
        assert_eq!(
            matrix.lambda_labels().unwrap(),
            &["coul-lambda", "vdw-lambda"]
        );
    }
}
