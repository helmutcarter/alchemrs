use std::collections::HashMap;

use crate::data::StatePoint;
use crate::error::{ensure_finite, ensure_finite_or_negative_infinity, CoreError, Result};

const ATM_TEMPERATURE_TOLERANCE: f64 = 1.0e-9;
const KB_KCAL_PER_MOL_K: f64 = 0.001_986_209;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtmDirection {
    Forward,
    Reverse,
}

impl AtmDirection {
    pub fn sign(self) -> i8 {
        match self {
            Self::Forward => 1,
            Self::Reverse => -1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtmState {
    state_id: usize,
    direction: AtmDirection,
    lambda1: f64,
    lambda2: f64,
    lambda3: Option<f64>,
    alpha: f64,
    u0: f64,
    u1: Option<f64>,
    w0: f64,
    temperature_k: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtmSchedule {
    states: Vec<AtmState>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtmSample {
    state_id: usize,
    potential_energy_kcal_per_mol: f64,
    perturbation_energy_kcal_per_mol: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtmSampleSet {
    schedule: AtmSchedule,
    samples: Vec<AtmSample>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtmLogQMatrix {
    n_observations: usize,
    n_states: usize,
    log_q: Vec<f64>,
    states: Vec<StatePoint>,
    sampled_counts: Vec<usize>,
    lambda_labels: Option<Vec<String>>,
}

impl AtmState {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        state_id: usize,
        direction: AtmDirection,
        lambda1: f64,
        lambda2: f64,
        lambda3: Option<f64>,
        alpha: f64,
        u0: f64,
        u1: Option<f64>,
        w0: f64,
        temperature_k: f64,
    ) -> Result<Self> {
        let has_multi = lambda3.is_some() || u1.is_some();
        if has_multi && (lambda3.is_none() || u1.is_none()) {
            return Err(CoreError::InvalidState(
                "lambda3 and u1 must either both be present or both be absent".to_string(),
            ));
        }

        ensure_finite(
            "atm_state",
            &[
                lambda1,
                lambda2,
                alpha,
                u0,
                w0,
                temperature_k,
                lambda3.unwrap_or(0.0),
                u1.unwrap_or(0.0),
            ],
        )?;

        if temperature_k <= 0.0 {
            return Err(CoreError::InvalidState(
                "ATM state temperature must be finite and positive".to_string(),
            ));
        }

        Ok(Self {
            state_id,
            direction,
            lambda1,
            lambda2,
            lambda3,
            alpha,
            u0,
            u1,
            w0,
            temperature_k,
        })
    }

    pub fn state_id(&self) -> usize {
        self.state_id
    }

    pub fn direction(&self) -> AtmDirection {
        self.direction
    }

    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    pub fn lambda1(&self) -> f64 {
        self.lambda1
    }

    pub fn lambda2(&self) -> f64 {
        self.lambda2
    }

    pub fn lambda3(&self) -> Option<f64> {
        self.lambda3
    }

    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    pub fn u0(&self) -> f64 {
        self.u0
    }

    pub fn u1(&self) -> Option<f64> {
        self.u1
    }

    pub fn w0(&self) -> f64 {
        self.w0
    }

    fn is_multi(&self) -> bool {
        self.lambda3.is_some()
    }

    fn parameter_vector(&self) -> Vec<f64> {
        let mut parameters = vec![self.lambda1, self.lambda2];
        if let Some(lambda3) = self.lambda3 {
            parameters.push(lambda3);
        }
        parameters.push(self.alpha);
        parameters.push(self.u0);
        if let Some(u1) = self.u1 {
            parameters.push(u1);
        }
        parameters.push(self.w0);
        parameters
    }

    fn to_state_point(&self) -> Result<StatePoint> {
        StatePoint::new(self.parameter_vector(), self.temperature_k)
    }

    fn bias_kcal_per_mol(&self, perturbation_energy_kcal_per_mol: f64) -> f64 {
        if let (Some(lambda3), Some(u1)) = (self.lambda3, self.u1) {
            let f3 = lambda3 * perturbation_energy_kcal_per_mol + self.w0;
            if self.alpha > 0.0 {
                let f2 = self.lambda2 * perturbation_energy_kcal_per_mol
                    + self.w0
                    + (lambda3 - self.lambda2) * u1;
                let f1 = self.lambda1 * perturbation_energy_kcal_per_mol
                    + self.w0
                    + (lambda3 - self.lambda2) * u1
                    + (self.lambda2 - self.lambda1) * self.u0;
                let average =
                    ((self.alpha * f1).exp() + (self.alpha * f2).exp() + (self.alpha * f3).exp())
                        / 3.0;
                average.ln() / self.alpha
            } else {
                f3
            }
        } else {
            let softplus = if self.alpha > 0.0 {
                let ee = 1.0 + (-self.alpha * (perturbation_energy_kcal_per_mol - self.u0)).exp();
                (self.lambda2 - self.lambda1) * ee.ln() / self.alpha
            } else {
                0.0
            };
            softplus + self.lambda2 * perturbation_energy_kcal_per_mol + self.w0
        }
    }
}

impl AtmSchedule {
    pub fn new(mut states: Vec<AtmState>) -> Result<Self> {
        if states.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }

        states.sort_by_key(AtmState::state_id);

        let direction = states[0].direction();
        let temperature_k = states[0].temperature_k();
        let multi = states[0].is_multi();
        let mut previous_state_id = None;
        for state in &states {
            if Some(state.state_id()) == previous_state_id {
                return Err(CoreError::InvalidState(
                    "ATM schedule state ids must be unique".to_string(),
                ));
            }
            previous_state_id = Some(state.state_id());

            if state.direction() != direction {
                return Err(CoreError::InvalidState(
                    "ATM schedule must contain exactly one leg direction".to_string(),
                ));
            }
            if (state.temperature_k() - temperature_k).abs() > ATM_TEMPERATURE_TOLERANCE {
                return Err(CoreError::InvalidState(
                    "ATM schedule must use a single temperature".to_string(),
                ));
            }
            if state.is_multi() != multi {
                return Err(CoreError::InvalidState(
                    "ATM schedule must use either standard or multi-softplus states consistently"
                        .to_string(),
                ));
            }
        }

        Ok(Self { states })
    }

    pub fn states(&self) -> &[AtmState] {
        &self.states
    }

    pub fn direction(&self) -> AtmDirection {
        self.states[0].direction()
    }

    pub fn temperature_k(&self) -> f64 {
        self.states[0].temperature_k()
    }

    pub fn uses_multi_softplus(&self) -> bool {
        self.states[0].is_multi()
    }

    fn ordered_states(&self) -> Vec<&AtmState> {
        let mut ordered = self.states.iter().collect::<Vec<_>>();
        if self.direction() == AtmDirection::Reverse {
            ordered.reverse();
        }
        ordered
    }

    fn state_lookup(&self) -> HashMap<usize, &AtmState> {
        self.states
            .iter()
            .map(|state| (state.state_id(), state))
            .collect()
    }

    fn lambda_labels(&self) -> Vec<String> {
        let mut labels = vec!["lambda1".to_string(), "lambda2".to_string()];
        if self.uses_multi_softplus() {
            labels.push("lambda3".to_string());
        }
        labels.push("alpha".to_string());
        labels.push("u0".to_string());
        if self.uses_multi_softplus() {
            labels.push("u1".to_string());
        }
        labels.push("w0".to_string());
        labels
    }
}

impl AtmSample {
    pub fn new(
        state_id: usize,
        potential_energy_kcal_per_mol: f64,
        perturbation_energy_kcal_per_mol: f64,
    ) -> Result<Self> {
        ensure_finite(
            "atm_sample",
            &[
                potential_energy_kcal_per_mol,
                perturbation_energy_kcal_per_mol,
            ],
        )?;

        Ok(Self {
            state_id,
            potential_energy_kcal_per_mol,
            perturbation_energy_kcal_per_mol,
        })
    }

    pub fn state_id(&self) -> usize {
        self.state_id
    }

    pub fn potential_energy_kcal_per_mol(&self) -> f64 {
        self.potential_energy_kcal_per_mol
    }

    pub fn perturbation_energy_kcal_per_mol(&self) -> f64 {
        self.perturbation_energy_kcal_per_mol
    }
}

impl AtmSampleSet {
    pub fn new(schedule: AtmSchedule, samples: Vec<AtmSample>) -> Result<Self> {
        if samples.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }

        let state_lookup = schedule.state_lookup();
        for sample in &samples {
            if !state_lookup.contains_key(&sample.state_id()) {
                return Err(CoreError::InvalidState(format!(
                    "ATM sample references unknown state_id {}",
                    sample.state_id()
                )));
            }
        }

        Ok(Self { schedule, samples })
    }

    pub fn schedule(&self) -> &AtmSchedule {
        &self.schedule
    }

    pub fn samples(&self) -> &[AtmSample] {
        &self.samples
    }

    pub fn to_log_q_matrix(&self) -> Result<AtmLogQMatrix> {
        let ordered_states = self.schedule.ordered_states();
        let ordered_state_points = ordered_states
            .iter()
            .map(|state| state.to_state_point())
            .collect::<Result<Vec<_>>>()?;
        let ordered_ids = ordered_states
            .iter()
            .map(|state| state.state_id())
            .collect::<Vec<_>>();
        let index_by_state_id = ordered_ids
            .iter()
            .enumerate()
            .map(|(idx, &state_id)| (state_id, idx))
            .collect::<HashMap<_, _>>();
        let state_lookup = self.schedule.state_lookup();
        let beta = 1.0 / (KB_KCAL_PER_MOL_K * self.schedule.temperature_k());
        let mut sampled_counts = vec![0_usize; ordered_states.len()];
        let mut log_q = Vec::with_capacity(self.samples.len() * ordered_states.len());

        for sample in &self.samples {
            let sampled_index = *index_by_state_id.get(&sample.state_id()).ok_or_else(|| {
                CoreError::InvalidState(format!(
                    "ATM sample references unknown state_id {}",
                    sample.state_id()
                ))
            })?;
            sampled_counts[sampled_index] += 1;

            let sampled_state = state_lookup
                .get(&sample.state_id())
                .expect("validated ATM sample state_id");
            let e0 = sample.potential_energy_kcal_per_mol()
                - sampled_state.bias_kcal_per_mol(sample.perturbation_energy_kcal_per_mol());

            for state in &ordered_states {
                let negative_reduced_potential = -beta
                    * (e0 + state.bias_kcal_per_mol(sample.perturbation_energy_kcal_per_mol()));
                log_q.push(negative_reduced_potential);
            }
        }

        AtmLogQMatrix::new_with_labels(
            self.samples.len(),
            ordered_states.len(),
            log_q,
            ordered_state_points,
            sampled_counts,
            Some(self.schedule.lambda_labels()),
        )
    }
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
    use super::{
        AtmDirection, AtmLogQMatrix, AtmSample, AtmSampleSet, AtmSchedule, AtmState,
        KB_KCAL_PER_MOL_K,
    };
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

    #[test]
    fn atm_schedule_rejects_mixed_directions() {
        let states = vec![
            AtmState::new(
                0,
                AtmDirection::Forward,
                0.0,
                0.0,
                None,
                0.0,
                0.0,
                None,
                0.0,
                300.0,
            )
            .unwrap(),
            AtmState::new(
                1,
                AtmDirection::Reverse,
                1.0,
                1.0,
                None,
                0.0,
                0.0,
                None,
                0.0,
                300.0,
            )
            .unwrap(),
        ];

        let err = AtmSchedule::new(states).unwrap_err();
        assert!(matches!(err, CoreError::InvalidState(_)));
    }

    #[test]
    fn atm_sample_set_builds_log_q_matrix_for_forward_leg() {
        let schedule = AtmSchedule::new(vec![
            AtmState::new(
                0,
                AtmDirection::Forward,
                0.0,
                0.0,
                None,
                0.0,
                0.0,
                None,
                0.0,
                300.0,
            )
            .unwrap(),
            AtmState::new(
                1,
                AtmDirection::Forward,
                1.0,
                1.0,
                None,
                0.0,
                0.0,
                None,
                0.0,
                300.0,
            )
            .unwrap(),
        ])
        .unwrap();
        let samples = vec![
            AtmSample::new(0, 10.0, 2.0).unwrap(),
            AtmSample::new(1, 12.0, 2.0).unwrap(),
        ];
        let matrix = AtmSampleSet::new(schedule, samples)
            .unwrap()
            .to_log_q_matrix()
            .unwrap();

        let beta = 1.0 / (KB_KCAL_PER_MOL_K * 300.0);
        let expected = [-beta * 10.0, -beta * 12.0, -beta * 10.0, -beta * 12.0];
        for (lhs, rhs) in matrix.log_q().iter().zip(expected.iter()) {
            assert!((lhs - rhs).abs() < 1.0e-12);
        }
        assert_eq!(matrix.sampled_counts(), &[1, 1]);
        assert_eq!(
            matrix.lambda_labels().unwrap(),
            &["lambda1", "lambda2", "alpha", "u0", "w0"]
        );
    }

    #[test]
    fn atm_sample_set_reverses_state_order_for_reverse_leg() {
        let schedule = AtmSchedule::new(vec![
            AtmState::new(
                0,
                AtmDirection::Reverse,
                0.0,
                0.0,
                None,
                0.0,
                0.0,
                None,
                0.0,
                300.0,
            )
            .unwrap(),
            AtmState::new(
                1,
                AtmDirection::Reverse,
                1.0,
                1.0,
                None,
                0.0,
                0.0,
                None,
                0.0,
                300.0,
            )
            .unwrap(),
        ])
        .unwrap();
        let samples = vec![
            AtmSample::new(0, 10.0, 2.0).unwrap(),
            AtmSample::new(1, 12.0, 2.0).unwrap(),
        ];
        let matrix = AtmSampleSet::new(schedule, samples)
            .unwrap()
            .to_log_q_matrix()
            .unwrap();

        assert_eq!(matrix.sampled_counts(), &[1, 1]);
        assert_eq!(matrix.states()[0].lambdas()[0], 1.0);
        assert_eq!(matrix.states()[1].lambdas()[0], 0.0);
    }
}
