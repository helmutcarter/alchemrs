use crate::data::{find_state_index_exact, StatePoint};
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
pub struct SwitchingTrajectory {
    initial_state: StatePoint,
    final_state: StatePoint,
    reduced_work: f64,
    lambda_path: Vec<f64>,
    dvdl_path: Vec<f64>,
    rms_dvdl_path: Vec<f64>,
}

impl SwitchingTrajectory {
    pub fn new(
        initial_state: StatePoint,
        final_state: StatePoint,
        reduced_work: f64,
    ) -> Result<Self> {
        Self::new_with_profile(
            initial_state,
            final_state,
            reduced_work,
            Vec::new(),
            Vec::new(),
        )
    }

    pub fn new_with_profile(
        initial_state: StatePoint,
        final_state: StatePoint,
        reduced_work: f64,
        lambda_path: Vec<f64>,
        dvdl_path: Vec<f64>,
    ) -> Result<Self> {
        Self::new_with_profile_and_rms(
            initial_state,
            final_state,
            reduced_work,
            lambda_path,
            dvdl_path,
            Vec::new(),
        )
    }

    pub fn new_with_profile_and_rms(
        initial_state: StatePoint,
        final_state: StatePoint,
        reduced_work: f64,
        lambda_path: Vec<f64>,
        dvdl_path: Vec<f64>,
        rms_dvdl_path: Vec<f64>,
    ) -> Result<Self> {
        ensure_finite("reduced_work", &[reduced_work])?;
        if initial_state.temperature_k() != final_state.temperature_k() {
            return Err(CoreError::InvalidState(
                "switching trajectory states must have the same temperature".to_string(),
            ));
        }
        if lambda_path.len() != dvdl_path.len() {
            return Err(CoreError::InvalidShape {
                expected: lambda_path.len(),
                found: dvdl_path.len(),
            });
        }
        ensure_finite("lambda_path", &lambda_path)?;
        ensure_finite("dvdl_path", &dvdl_path)?;
        if !rms_dvdl_path.is_empty() && rms_dvdl_path.len() != lambda_path.len() {
            return Err(CoreError::InvalidShape {
                expected: lambda_path.len(),
                found: rms_dvdl_path.len(),
            });
        }
        ensure_finite("rms_dvdl_path", &rms_dvdl_path)?;
        Ok(Self {
            initial_state,
            final_state,
            reduced_work,
            lambda_path,
            dvdl_path,
            rms_dvdl_path,
        })
    }

    pub fn initial_state(&self) -> &StatePoint {
        &self.initial_state
    }

    pub fn final_state(&self) -> &StatePoint {
        &self.final_state
    }

    pub fn reduced_work(&self) -> f64 {
        self.reduced_work
    }

    pub fn lambda_path(&self) -> &[f64] {
        &self.lambda_path
    }

    pub fn dvdl_path(&self) -> &[f64] {
        &self.dvdl_path
    }

    pub fn rms_dvdl_path(&self) -> &[f64] {
        &self.rms_dvdl_path
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NesMbarSample {
    step_index: usize,
    time_ps: f64,
    lambda_protocol: f64,
    reduced_work: f64,
    reduced_energy_protocol: f64,
    reduced_energies_states: Vec<f64>,
}

impl NesMbarSample {
    pub fn new(
        step_index: usize,
        time_ps: f64,
        lambda_protocol: f64,
        reduced_work: f64,
        reduced_energy_protocol: f64,
        reduced_energies_states: Vec<f64>,
    ) -> Result<Self> {
        if step_index == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        ensure_finite("time", &[time_ps])?;
        ensure_finite("lambda_protocol", &[lambda_protocol])?;
        ensure_finite("reduced_work", &[reduced_work])?;
        ensure_finite("reduced_energy_protocol", &[reduced_energy_protocol])?;
        ensure_finite_or_positive_infinity("reduced_energies_states", &reduced_energies_states)?;
        if reduced_energies_states.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        Ok(Self {
            step_index,
            time_ps,
            lambda_protocol,
            reduced_work,
            reduced_energy_protocol,
            reduced_energies_states,
        })
    }

    pub fn step_index(&self) -> usize {
        self.step_index
    }

    pub fn time_ps(&self) -> f64 {
        self.time_ps
    }

    pub fn lambda_protocol(&self) -> f64 {
        self.lambda_protocol
    }

    pub fn reduced_work(&self) -> f64 {
        self.reduced_work
    }

    pub fn reduced_energy_protocol(&self) -> f64 {
        self.reduced_energy_protocol
    }

    pub fn reduced_energies_states(&self) -> &[f64] {
        &self.reduced_energies_states
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NesMbarTrajectory {
    initial_state: StatePoint,
    final_state: StatePoint,
    target_states: Vec<StatePoint>,
    samples: Vec<NesMbarSample>,
}

impl NesMbarTrajectory {
    pub fn new(
        initial_state: StatePoint,
        final_state: StatePoint,
        target_states: Vec<StatePoint>,
        samples: Vec<NesMbarSample>,
    ) -> Result<Self> {
        if initial_state.temperature_k() != final_state.temperature_k() {
            return Err(CoreError::InvalidState(
                "NES MBAR trajectory states must have the same temperature".to_string(),
            ));
        }
        if target_states.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        if samples.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        for state in &target_states {
            if (state.temperature_k() - initial_state.temperature_k()).abs() > 0.0 {
                return Err(CoreError::InvalidState(
                    "NES MBAR target states must have the same temperature".to_string(),
                ));
            }
        }
        let expected_n_states = target_states.len();
        let mut previous_time = f64::NEG_INFINITY;
        for sample in &samples {
            if sample.reduced_energies_states().len() != expected_n_states {
                return Err(CoreError::InvalidShape {
                    expected: expected_n_states,
                    found: sample.reduced_energies_states().len(),
                });
            }
            if sample.time_ps() < previous_time {
                return Err(CoreError::InvalidTimeOrder(
                    "NES MBAR sample times must be nondecreasing".to_string(),
                ));
            }
            previous_time = sample.time_ps();
        }

        Ok(Self {
            initial_state,
            final_state,
            target_states,
            samples,
        })
    }

    pub fn initial_state(&self) -> &StatePoint {
        &self.initial_state
    }

    pub fn final_state(&self) -> &StatePoint {
        &self.final_state
    }

    pub fn target_states(&self) -> &[StatePoint] {
        &self.target_states
    }

    pub fn samples(&self) -> &[NesMbarSample] {
        &self.samples
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
        if n_states == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
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
        if let Some(sampled) = sampled_state.as_ref() {
            if find_state_index_exact(&evaluated_states, sampled).is_err() {
                return Err(CoreError::InvalidState(
                    "sampled_state not found in evaluated_states".to_string(),
                ));
            }
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
    fn unk_matrix_rejects_empty_state_grid() {
        let err = UNkMatrix::new(0, 0, vec![], vec![], None, vec![]).unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 1,
                found: 0
            }
        ));
    }

    #[test]
    fn unk_matrix_rejects_sampled_state_missing_from_grid_exactly() {
        let sampled = StatePoint::new(vec![0.1 + 0.2], 300.0).unwrap();
        let evaluated = StatePoint::new(vec![0.3], 300.0).unwrap();
        let err =
            UNkMatrix::new(1, 1, vec![0.0], vec![0.0], Some(sampled), vec![evaluated]).unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidState(message)
                if message == "sampled_state not found in evaluated_states"
        ));
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
