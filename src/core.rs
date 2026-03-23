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

#[derive(Debug, Clone, PartialEq)]
pub struct DhdlSeries {
    state: StatePoint,
    time_ps: Vec<f64>,
    values: Vec<f64>,
}

impl DhdlSeries {
    pub fn new(state: StatePoint, time_ps: Vec<f64>, values: Vec<f64>) -> Result<Self> {
        if time_ps.len() != values.len() {
            return Err(CoreError::InvalidShape {
                expected: time_ps.len(),
                found: values.len(),
            });
        }
        ensure_finite("time", &time_ps)?;
        ensure_finite("values", &values)?;
        for idx in 1..time_ps.len() {
            if time_ps[idx] < time_ps[idx - 1] {
                return Err(CoreError::InvalidTimeOrder(format!(
                    "time[{idx}] < time[{prev}]",
                    prev = idx - 1
                )));
            }
        }
        Ok(Self {
            state,
            time_ps,
            values,
        })
    }

    pub fn state(&self) -> &StatePoint {
        &self.state
    }

    pub fn time_ps(&self) -> &[f64] {
        &self.time_ps
    }

    pub fn values(&self) -> &[f64] {
        &self.values
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct UNkMatrix {
    n_samples: usize,
    n_states: usize,
    data: Vec<f64>,
    time_ps: Vec<f64>,
    sampled_state: Option<StatePoint>,
    evaluated_states: Vec<StatePoint>,
}

impl UNkMatrix {
    pub fn new(
        n_samples: usize,
        n_states: usize,
        data: Vec<f64>,
        time_ps: Vec<f64>,
        sampled_state: Option<StatePoint>,
        evaluated_states: Vec<StatePoint>,
    ) -> Result<Self> {
        let expected = n_samples
            .checked_mul(n_states)
            .ok_or(CoreError::InvalidShape {
                expected: n_samples,
                found: n_states,
            })?;
        if data.len() != expected {
            return Err(CoreError::InvalidShape {
                expected,
                found: data.len(),
            });
        }
        if evaluated_states.len() != n_states {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: evaluated_states.len(),
            });
        }
        if time_ps.len() != n_samples {
            return Err(CoreError::InvalidShape {
                expected: n_samples,
                found: time_ps.len(),
            });
        }
        ensure_finite_or_positive_infinity("u_nk", &data)?;
        ensure_finite("time", &time_ps)?;
        for idx in 1..time_ps.len() {
            if time_ps[idx] < time_ps[idx - 1] {
                return Err(CoreError::InvalidTimeOrder(format!(
                    "time[{idx}] < time[{prev}]",
                    prev = idx - 1
                )));
            }
        }
        Ok(Self {
            n_samples,
            n_states,
            data,
            time_ps,
            sampled_state,
            evaluated_states,
        })
    }

    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    pub fn n_states(&self) -> usize {
        self.n_states
    }

    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn time_ps(&self) -> &[f64] {
        &self.time_ps
    }

    pub fn sampled_state(&self) -> Option<&StatePoint> {
        self.sampled_state.as_ref()
    }

    pub fn evaluated_states(&self) -> &[StatePoint] {
        &self.evaluated_states
    }
}

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
}

impl DeltaFMatrix {
    pub fn new(
        values: Vec<f64>,
        uncertainties: Option<Vec<f64>>,
        n_states: usize,
        states: Vec<StatePoint>,
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
        ensure_finite("delta_f", &values)?;
        Ok(Self {
            values,
            uncertainties,
            n_states,
            states,
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
}

#[derive(Debug, Clone, PartialEq)]
pub struct OverlapMatrix {
    values: Vec<f64>,
    n_states: usize,
    states: Vec<StatePoint>,
}

impl OverlapMatrix {
    pub fn new(values: Vec<f64>, n_states: usize, states: Vec<StatePoint>) -> Result<Self> {
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
        ensure_finite("overlap_matrix", &values)?;
        Ok(Self {
            values,
            n_states,
            states,
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
}

fn ensure_finite(label: &str, values: &[f64]) -> Result<()> {
    for (idx, value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(CoreError::NonFiniteValue(format!(
                "{label}[{idx}] must be finite"
            )));
        }
    }
    Ok(())
}

fn ensure_finite_or_positive_infinity(label: &str, values: &[f64]) -> Result<()> {
    for (idx, value) in values.iter().enumerate() {
        if value.is_nan() || *value == f64::NEG_INFINITY {
            return Err(CoreError::NonFiniteValue(format!(
                "{label}[{idx}] must be finite or +inf"
            )));
        }
    }
    Ok(())
}

fn ensure_finite_or_nan(label: &str, values: &[f64]) -> Result<()> {
    for (idx, value) in values.iter().enumerate() {
        if value.is_infinite() {
            return Err(CoreError::NonFiniteValue(format!(
                "{label}[{idx}] must be finite or NaN"
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn dhdl_series_checks_lengths() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let err = DhdlSeries::new(state, vec![0.0], vec![1.0, 2.0]).unwrap_err();
        assert!(matches!(err, CoreError::InvalidShape { .. }));
    }

    #[test]
    fn dhdl_series_checks_time_order() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let err = DhdlSeries::new(state, vec![1.0, 0.5], vec![1.0, 2.0]).unwrap_err();
        assert!(matches!(err, CoreError::InvalidTimeOrder(_)));
    }

    #[test]
    fn unk_matrix_checks_shape() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let err = UNkMatrix::new(2, 2, vec![1.0, 2.0, 3.0], vec![0.0, 1.0], None, vec![state])
            .unwrap_err();
        assert!(matches!(err, CoreError::InvalidShape { .. }));
    }

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
    fn unk_matrix_allows_positive_infinity() {
        let state0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let matrix = UNkMatrix::new(
            2,
            2,
            vec![0.0, f64::INFINITY, 0.0, 1.0],
            vec![0.0, 1.0],
            Some(state0.clone()),
            vec![state0, state1],
        )
        .unwrap();
        assert_eq!(matrix.data()[1], f64::INFINITY);
    }

    #[test]
    fn unk_matrix_rejects_negative_infinity() {
        let state0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let err = UNkMatrix::new(
            1,
            2,
            vec![0.0, f64::NEG_INFINITY],
            vec![0.0],
            Some(state0.clone()),
            vec![state0, state1],
        )
        .unwrap_err();
        assert!(matches!(err, CoreError::NonFiniteValue(_)));
    }
}
