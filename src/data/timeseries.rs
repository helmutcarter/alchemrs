use crate::data::StatePoint;
use crate::error::{ensure_finite, ensure_finite_or_positive_infinity, CoreError, Result};

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
    lambda_labels: Option<Vec<String>>,
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
        Self::new_with_labels(
            n_samples,
            n_states,
            data,
            time_ps,
            sampled_state,
            evaluated_states,
            None,
        )
    }

    pub fn new_with_labels(
        n_samples: usize,
        n_states: usize,
        data: Vec<f64>,
        time_ps: Vec<f64>,
        sampled_state: Option<StatePoint>,
        evaluated_states: Vec<StatePoint>,
        lambda_labels: Option<Vec<String>>,
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
        if let Some(labels) = lambda_labels.as_ref() {
            let expected = sampled_state
                .as_ref()
                .map(|state| state.lambdas().len())
                .or_else(|| evaluated_states.first().map(|state| state.lambdas().len()))
                .unwrap_or(0);
            if labels.len() != expected {
                return Err(CoreError::InvalidShape {
                    expected,
                    found: labels.len(),
                });
            }
        }
        Ok(Self {
            n_samples,
            n_states,
            data,
            time_ps,
            sampled_state,
            evaluated_states,
            lambda_labels,
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

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::{DhdlSeries, UNkMatrix};
    use crate::data::StatePoint;
    use crate::error::CoreError;

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
