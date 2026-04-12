use crate::data::StatePoint;
use crate::error::{ensure_finite_or_negative_infinity, CoreError, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct AtmLogQMatrix {
    n_observations: usize,
    n_states: usize,
    log_q: Vec<f64>,
    states: Vec<StatePoint>,
    sampled_counts: Vec<usize>,
    lambda_labels: Option<Vec<String>>,
}

impl AtmLogQMatrix {
    pub fn new(
        n_observations: usize,
        n_states: usize,
        log_q: Vec<f64>,
        states: Vec<StatePoint>,
        sampled_counts: Vec<usize>,
    ) -> Result<Self> {
        Self::new_with_labels(
            n_observations,
            n_states,
            log_q,
            states,
            sampled_counts,
            None,
        )
    }

    pub fn new_with_labels(
        n_observations: usize,
        n_states: usize,
        log_q: Vec<f64>,
        states: Vec<StatePoint>,
        sampled_counts: Vec<usize>,
        lambda_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        if n_observations == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        if n_states == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }

        let expected = n_observations
            .checked_mul(n_states)
            .ok_or(CoreError::InvalidShape {
                expected: n_observations,
                found: n_states,
            })?;
        if log_q.len() != expected {
            return Err(CoreError::InvalidShape {
                expected,
                found: log_q.len(),
            });
        }
        if states.len() != n_states {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: states.len(),
            });
        }
        if sampled_counts.len() != n_states {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: sampled_counts.len(),
            });
        }
        let observed = sampled_counts.iter().sum::<usize>();
        if observed != n_observations {
            return Err(CoreError::InvalidShape {
                expected: n_observations,
                found: observed,
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

        ensure_finite_or_negative_infinity("log_q", &log_q)?;

        Ok(Self {
            n_observations,
            n_states,
            log_q,
            states,
            sampled_counts,
            lambda_labels,
        })
    }

    pub fn n_observations(&self) -> usize {
        self.n_observations
    }

    pub fn n_states(&self) -> usize {
        self.n_states
    }

    pub fn log_q(&self) -> &[f64] {
        &self.log_q
    }

    pub fn states(&self) -> &[StatePoint] {
        &self.states
    }

    pub fn sampled_counts(&self) -> &[usize] {
        &self.sampled_counts
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::AtmLogQMatrix;
    use crate::data::StatePoint;
    use crate::error::CoreError;

    fn two_states() -> Vec<StatePoint> {
        vec![
            StatePoint::new(vec![0.0], 300.0).unwrap(),
            StatePoint::new(vec![1.0], 300.0).unwrap(),
        ]
    }

    #[test]
    fn atm_log_q_matrix_checks_shape() {
        let err =
            AtmLogQMatrix::new(2, 2, vec![0.0, -1.0, -2.0], two_states(), vec![1, 1]).unwrap_err();
        assert!(matches!(err, CoreError::InvalidShape { .. }));
    }

    #[test]
    fn atm_log_q_matrix_requires_counts_to_sum_to_observations() {
        let err = AtmLogQMatrix::new(2, 2, vec![0.0, -1.0, -2.0, -3.0], two_states(), vec![2, 2])
            .unwrap_err();
        assert!(matches!(err, CoreError::InvalidShape { .. }));
    }

    #[test]
    fn atm_log_q_matrix_rejects_positive_infinity_and_nan() {
        let err = AtmLogQMatrix::new(
            2,
            2,
            vec![0.0, f64::INFINITY, -2.0, -3.0],
            two_states(),
            vec![1, 1],
        )
        .unwrap_err();
        assert!(matches!(err, CoreError::NonFiniteValue(_)));

        let err = AtmLogQMatrix::new(
            2,
            2,
            vec![0.0, f64::NAN, -2.0, -3.0],
            two_states(),
            vec![1, 1],
        )
        .unwrap_err();
        assert!(matches!(err, CoreError::NonFiniteValue(_)));
    }

    #[test]
    fn atm_log_q_matrix_allows_negative_infinity_and_preserves_labels() {
        let matrix = AtmLogQMatrix::new_with_labels(
            2,
            2,
            vec![0.0, f64::NEG_INFINITY, -2.0, -3.0],
            two_states(),
            vec![1, 1],
            Some(vec!["lambda".to_string()]),
        )
        .unwrap();

        assert_eq!(matrix.n_observations(), 2);
        assert_eq!(matrix.n_states(), 2);
        assert_eq!(matrix.sampled_counts(), &[1, 1]);
        assert_eq!(matrix.lambda_labels().unwrap(), &["lambda"]);
    }
}
