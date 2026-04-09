use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use crate::data::{find_state_index_exact, DeltaFMatrix, NesMbarTrajectory, StatePoint};
use crate::error::{CoreError, Result};

#[derive(Debug, Clone)]
pub struct NesMbarOptions {
    pub n_bootstrap: usize,
    pub seed: u64,
    pub sample_stride: usize,
}

impl Default for NesMbarOptions {
    fn default() -> Self {
        Self {
            n_bootstrap: 0,
            seed: 0,
            sample_stride: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct NesMbarEstimator {
    pub options: NesMbarOptions,
}

#[derive(Debug, Clone)]
pub struct NesMbarFit {
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
    relative_free_energies: Vec<f64>,
    uncertainties: Option<Vec<f64>>,
}

impl Default for NesMbarEstimator {
    fn default() -> Self {
        Self::new(NesMbarOptions::default())
    }
}

impl NesMbarEstimator {
    pub fn new(options: NesMbarOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, trajectories: &[NesMbarTrajectory]) -> Result<NesMbarFit> {
        if trajectories.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        if self.options.sample_stride == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }

        let states = estimator_states(trajectories)?;
        let relative_free_energies =
            estimate_relative_free_energies(trajectories, &states, self.options.sample_stride)?;
        let uncertainties = if self.options.n_bootstrap == 0 {
            None
        } else {
            Some(bootstrap_uncertainties(
                trajectories,
                &states,
                self.options.sample_stride,
                self.options.n_bootstrap,
                self.options.seed,
            )?)
        };

        Ok(NesMbarFit {
            states,
            lambda_labels: None,
            relative_free_energies,
            uncertainties,
        })
    }

    pub fn estimate(&self, trajectories: &[NesMbarTrajectory]) -> Result<DeltaFMatrix> {
        self.fit(trajectories)?.result()
    }
}

impl NesMbarFit {
    pub fn relative_free_energies(&self) -> &[f64] {
        &self.relative_free_energies
    }

    pub fn uncertainties(&self) -> Option<&[f64]> {
        self.uncertainties.as_deref()
    }

    pub fn states(&self) -> &[StatePoint] {
        &self.states
    }

    pub fn result(&self) -> Result<DeltaFMatrix> {
        let n_states = self.states.len();
        let mut values = vec![0.0; n_states * n_states];
        let mut uncertainties = self
            .uncertainties
            .as_ref()
            .map(|_| vec![0.0; n_states * n_states]);
        for i in 0..n_states {
            for j in 0..n_states {
                let idx = i * n_states + j;
                values[idx] = self.relative_free_energies[j] - self.relative_free_energies[i];
                if let (Some(sigmas), Some(matrix_unc)) =
                    (self.uncertainties.as_ref(), uncertainties.as_mut())
                {
                    matrix_unc[idx] = (sigmas[i] * sigmas[i] + sigmas[j] * sigmas[j]).sqrt();
                }
            }
        }
        DeltaFMatrix::new_with_labels(
            values,
            uncertainties,
            n_states,
            self.states.clone(),
            self.lambda_labels.clone(),
        )
    }
}

fn estimator_states(trajectories: &[NesMbarTrajectory]) -> Result<Vec<StatePoint>> {
    let first = &trajectories[0];
    let mut states = vec![first.initial_state().clone()];
    for state in first.target_states() {
        if find_state_index_exact(&states, state).is_err() {
            states.push(state.clone());
        }
    }
    for trajectory in trajectories.iter().skip(1) {
        if trajectory.initial_state() != first.initial_state() {
            return Err(CoreError::InvalidState(
                "all NES MBAR trajectories must share the same initial state".to_string(),
            ));
        }
        if trajectory.target_states() != first.target_states() {
            return Err(CoreError::InvalidState(
                "all NES MBAR trajectories must share the same target state grid".to_string(),
            ));
        }
    }
    Ok(states)
}

fn estimate_relative_free_energies(
    trajectories: &[NesMbarTrajectory],
    states: &[StatePoint],
    sample_stride: usize,
) -> Result<Vec<f64>> {
    let initial_idx = 0usize;
    let mut relative = vec![0.0; states.len()];
    let target_states = trajectories[0].target_states();
    for (state_idx, state) in states.iter().enumerate() {
        if state_idx == initial_idx {
            continue;
        }
        let target_idx = find_state_index_exact(target_states, state).map_err(|_| {
            CoreError::InvalidState(
                "NES MBAR target state not found in trajectory target grid".to_string(),
            )
        })?;
        let mut max_log_weight = f64::NEG_INFINITY;
        let mut n_terms = 0usize;
        for trajectory in trajectories {
            for sample in trajectory.samples().iter().step_by(sample_stride) {
                let log_weight = -sample.reduced_work()
                    - (sample.reduced_energies_states()[target_idx]
                        - sample.reduced_energy_protocol());
                if log_weight.is_finite() {
                    max_log_weight = max_log_weight.max(log_weight);
                    n_terms += 1;
                }
            }
        }
        if n_terms == 0 {
            return Err(CoreError::NonFiniteValue(
                "NES MBAR state estimate had no finite nonequilibrium weights".to_string(),
            ));
        }
        let mut sum = 0.0;
        for trajectory in trajectories {
            for sample in trajectory.samples().iter().step_by(sample_stride) {
                let log_weight = -sample.reduced_work()
                    - (sample.reduced_energies_states()[target_idx]
                        - sample.reduced_energy_protocol());
                if log_weight.is_finite() {
                    sum += (log_weight - max_log_weight).exp();
                }
            }
        }
        let mean_log = max_log_weight + sum.ln() - (n_terms as f64).ln();
        relative[state_idx] = -mean_log;
    }
    Ok(relative)
}

fn bootstrap_uncertainties(
    trajectories: &[NesMbarTrajectory],
    states: &[StatePoint],
    sample_stride: usize,
    n_bootstrap: usize,
    seed: u64,
) -> Result<Vec<f64>> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut replicas = vec![Vec::with_capacity(n_bootstrap); states.len()];
    let mut resample = Vec::with_capacity(trajectories.len());
    for _ in 0..n_bootstrap {
        resample.clear();
        for _ in 0..trajectories.len() {
            let idx = rng.gen_range(0..trajectories.len());
            resample.push(trajectories[idx].clone());
        }
        let rel = estimate_relative_free_energies(&resample, states, sample_stride)?;
        for (idx, value) in rel.into_iter().enumerate() {
            replicas[idx].push(value);
        }
    }
    let mut uncertainties = Vec::with_capacity(states.len());
    for samples in replicas {
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let variance = samples
            .iter()
            .map(|value| {
                let diff = *value - mean;
                diff * diff
            })
            .sum::<f64>()
            / samples.len() as f64;
        uncertainties.push(variance.sqrt());
    }
    Ok(uncertainties)
}
