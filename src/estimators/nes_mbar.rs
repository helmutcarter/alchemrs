use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use crate::data::{find_state_index_exact, DeltaFMatrix, NesMbarTrajectory, StatePoint};
use crate::error::{CoreError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NesMbarWeighting {
    #[default]
    Uniform,
    SliceEss,
    SliceEssFiltered,
    InverseVariance,
}

#[derive(Debug, Clone)]
pub struct NesMbarOptions {
    pub n_bootstrap: usize,
    pub seed: u64,
    pub sample_stride: usize,
    pub weighting: NesMbarWeighting,
    pub max_slice_trajectory_fraction: Option<f64>,
    pub min_slice_ess_fraction: Option<f64>,
}

impl Default for NesMbarOptions {
    fn default() -> Self {
        Self {
            n_bootstrap: 0,
            seed: 0,
            sample_stride: 1,
            weighting: NesMbarWeighting::Uniform,
            max_slice_trajectory_fraction: None,
            min_slice_ess_fraction: None,
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

#[derive(Debug, Clone, PartialEq)]
pub struct NesMbarContribution {
    pub index: usize,
    pub fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NesMbarStateDiagnostics {
    pub state: StatePoint,
    pub effective_sample_size: f64,
    pub max_sample_fraction: f64,
    pub top_slices: Vec<NesMbarContribution>,
    pub top_trajectories: Vec<NesMbarContribution>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NesMbarDiagnostics {
    pub states: Vec<NesMbarStateDiagnostics>,
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
            estimate_relative_free_energies(trajectories, &states, &self.options)?;
        let uncertainties = if self.options.n_bootstrap == 0 {
            None
        } else {
            Some(bootstrap_uncertainties(
                trajectories,
                &states,
                &self.options,
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

    pub fn diagnostics(&self, trajectories: &[NesMbarTrajectory]) -> Result<NesMbarDiagnostics> {
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
        let target_states = trajectories[0].target_states();
        let mut diagnostics = Vec::with_capacity(states.len());
        diagnostics.push(NesMbarStateDiagnostics {
            state: states[0].clone(),
            effective_sample_size: f64::INFINITY,
            max_sample_fraction: 0.0,
            top_slices: Vec::new(),
            top_trajectories: Vec::new(),
        });
        for state in states.iter().skip(1) {
            let target_idx = find_state_index_exact(target_states, state).map_err(|_| {
                CoreError::InvalidState(
                    "NES MBAR target state not found in trajectory target grid".to_string(),
                )
            })?;
            diagnostics.push(state_diagnostics(
                trajectories,
                state.clone(),
                self.options.sample_stride,
                target_idx,
            )?);
        }
        Ok(NesMbarDiagnostics {
            states: diagnostics,
        })
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
    options: &NesMbarOptions,
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
        let slice_estimates =
            estimate_slice_weights(trajectories, options.sample_stride, target_idx)?;
        if slice_estimates.is_empty() {
            return Err(CoreError::NonFiniteValue(
                "NES MBAR state estimate had no finite nonequilibrium weights".to_string(),
            ));
        }
        let combined_log_mean = combine_slice_estimates(
            &slice_estimates,
            options.weighting,
            trajectories.len(),
            options.max_slice_trajectory_fraction,
            options.min_slice_ess_fraction,
        )?;
        relative[state_idx] = -combined_log_mean;
    }
    Ok(relative)
}

#[derive(Debug, Clone, Copy)]
struct SliceEstimate {
    log_mean: f64,
    log_variance_of_mean: Option<f64>,
    trajectory_ess: f64,
    max_trajectory_fraction: f64,
}

fn estimate_slice_weights(
    trajectories: &[NesMbarTrajectory],
    sample_stride: usize,
    target_idx: usize,
) -> Result<Vec<SliceEstimate>> {
    let mut slice_logs: Vec<Vec<f64>> = Vec::new();
    let mut slice_by_trajectory: Vec<Vec<Vec<f64>>> = Vec::new();
    for (trajectory_idx, trajectory) in trajectories.iter().enumerate() {
        for (ordinal, sample) in trajectory
            .samples()
            .iter()
            .step_by(sample_stride)
            .enumerate()
        {
            let log_weight = -sample.reduced_work()
                - (sample.reduced_energies_states()[target_idx] - sample.reduced_energy_protocol());
            if !log_weight.is_finite() {
                continue;
            }
            if slice_logs.len() <= ordinal {
                slice_logs.resize_with(ordinal + 1, Vec::new);
                slice_by_trajectory
                    .resize_with(ordinal + 1, || vec![Vec::new(); trajectories.len()]);
            }
            slice_logs[ordinal].push(log_weight);
            slice_by_trajectory[ordinal][trajectory_idx].push(log_weight);
        }
    }

    let estimates = slice_logs
        .into_iter()
        .zip(slice_by_trajectory.into_iter())
        .filter_map(|logs| {
            let (logs, logs_by_trajectory) = logs;
            if logs.is_empty() {
                return None;
            }
            let max_log = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let n = logs.len() as f64;
            let scaled = logs
                .iter()
                .map(|log| (*log - max_log).exp())
                .collect::<Vec<_>>();
            let mean_scaled = scaled.iter().sum::<f64>() / n;
            let log_mean = max_log + mean_scaled.ln();
            let log_variance_of_mean = if logs.len() > 1 {
                let mean = mean_scaled;
                let m2 = scaled
                    .iter()
                    .map(|value| {
                        let diff = *value - mean;
                        diff * diff
                    })
                    .sum::<f64>();
                let variance_of_mean_scaled = (m2 / (n - 1.0)) / n;
                (variance_of_mean_scaled > 0.0)
                    .then_some(2.0 * max_log + variance_of_mean_scaled.ln())
            } else {
                None
            };
            let sum_scaled = scaled.iter().sum::<f64>();
            let sum_scaled_sq = scaled.iter().map(|value| value * value).sum::<f64>();
            let trajectory_ess = if sum_scaled_sq > 0.0 {
                (sum_scaled * sum_scaled) / sum_scaled_sq
            } else {
                0.0
            };
            let mut trajectory_totals = Vec::new();
            for trajectory_logs in logs_by_trajectory {
                if trajectory_logs.is_empty() {
                    continue;
                }
                let total = trajectory_logs
                    .iter()
                    .map(|log| (*log - max_log).exp())
                    .sum::<f64>();
                trajectory_totals.push(total);
            }
            let total_weight = trajectory_totals.iter().sum::<f64>();
            let max_trajectory_fraction = if total_weight > 0.0 {
                trajectory_totals
                    .iter()
                    .map(|value| *value / total_weight)
                    .fold(0.0, f64::max)
            } else {
                1.0
            };
            Some(SliceEstimate {
                log_mean,
                log_variance_of_mean,
                trajectory_ess,
                max_trajectory_fraction,
            })
        })
        .collect::<Vec<_>>();
    Ok(estimates)
}

fn combine_slice_estimates(
    slice_estimates: &[SliceEstimate],
    weighting: NesMbarWeighting,
    n_trajectories: usize,
    max_slice_trajectory_fraction: Option<f64>,
    min_slice_ess_fraction: Option<f64>,
) -> Result<f64> {
    let filtered = filter_slice_estimates(
        slice_estimates,
        n_trajectories,
        max_slice_trajectory_fraction,
        min_slice_ess_fraction,
    );
    let slice_estimates = if filtered.is_empty() {
        slice_estimates
    } else {
        filtered.as_slice()
    };
    match weighting {
        NesMbarWeighting::Uniform => {
            let log_mean = log_mean_exp(
                &slice_estimates
                    .iter()
                    .map(|estimate| estimate.log_mean)
                    .collect::<Vec<_>>(),
            )?;
            validate_combined_log_mean(log_mean)
        }
        NesMbarWeighting::SliceEss => {
            let mut numerator_terms = Vec::with_capacity(slice_estimates.len());
            let mut denominator_terms = Vec::with_capacity(slice_estimates.len());
            for estimate in slice_estimates {
                let ess = estimate.trajectory_ess.max(1.0);
                let log_ess = ess.ln();
                numerator_terms.push(log_ess + estimate.log_mean);
                denominator_terms.push(log_ess);
            }
            let log_numerator = log_sum_exp(&numerator_terms)?;
            let log_denominator = log_sum_exp(&denominator_terms)?;
            validate_combined_log_mean(log_numerator - log_denominator)
        }
        NesMbarWeighting::SliceEssFiltered => {
            let mut numerator_terms = Vec::with_capacity(slice_estimates.len());
            let mut denominator_terms = Vec::with_capacity(slice_estimates.len());
            for estimate in slice_estimates {
                let ess_fraction = (estimate.trajectory_ess / n_trajectories as f64).max(1e-12);
                let log_weight = ess_fraction.ln();
                numerator_terms.push(log_weight + estimate.log_mean);
                denominator_terms.push(log_weight);
            }
            let log_numerator = log_sum_exp(&numerator_terms)?;
            let log_denominator = log_sum_exp(&denominator_terms)?;
            validate_combined_log_mean(log_numerator - log_denominator)
        }
        NesMbarWeighting::InverseVariance => {
            let min_log_variance = slice_estimates
                .iter()
                .filter_map(|estimate| estimate.log_variance_of_mean)
                .filter(|var| var.is_finite())
                .fold(f64::INFINITY, f64::min);
            if !min_log_variance.is_finite() {
                let log_mean = log_mean_exp(
                    &slice_estimates
                        .iter()
                        .map(|estimate| estimate.log_mean)
                        .collect::<Vec<_>>(),
                )?;
                return validate_combined_log_mean(log_mean);
            }
            let log_variance_floor = min_log_variance + 1e-6f64.ln();
            let mut numerator_terms = Vec::with_capacity(slice_estimates.len());
            let mut denominator_terms = Vec::with_capacity(slice_estimates.len());
            for estimate in slice_estimates {
                let log_variance = estimate
                    .log_variance_of_mean
                    .filter(|var| var.is_finite())
                    .unwrap_or(min_log_variance)
                    .max(log_variance_floor);
                let log_precision = -log_variance;
                numerator_terms.push(log_precision + estimate.log_mean);
                denominator_terms.push(log_precision);
            }
            let log_numerator = log_sum_exp(&numerator_terms)?;
            let log_denominator = log_sum_exp(&denominator_terms)?;
            validate_combined_log_mean(log_numerator - log_denominator)
        }
    }
}

fn filter_slice_estimates<'a>(
    slice_estimates: &'a [SliceEstimate],
    n_trajectories: usize,
    max_slice_trajectory_fraction: Option<f64>,
    min_slice_ess_fraction: Option<f64>,
) -> Vec<SliceEstimate> {
    slice_estimates
        .iter()
        .copied()
        .filter(|estimate| {
            let trajectory_ok = max_slice_trajectory_fraction
                .map(|cutoff| estimate.max_trajectory_fraction <= cutoff)
                .unwrap_or(true);
            let ess_ok = min_slice_ess_fraction
                .map(|cutoff| (estimate.trajectory_ess / n_trajectories as f64) >= cutoff)
                .unwrap_or(true);
            trajectory_ok && ess_ok
        })
        .collect()
}

fn validate_combined_log_mean(log_mean: f64) -> Result<f64> {
    if !log_mean.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "NES MBAR combined nonequilibrium weight became non-finite".to_string(),
        ));
    }
    Ok(log_mean)
}

fn log_sum_exp(log_values: &[f64]) -> Result<f64> {
    if log_values.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let max_log = log_values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if !max_log.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "all NES MBAR log-weights were non-finite".to_string(),
        ));
    }
    let sum = log_values
        .iter()
        .filter(|value| value.is_finite())
        .map(|value| (*value - max_log).exp())
        .sum::<f64>();
    Ok(max_log + sum.ln())
}

fn log_mean_exp(log_values: &[f64]) -> Result<f64> {
    Ok(log_sum_exp(log_values)? - (log_values.len() as f64).ln())
}

fn bootstrap_uncertainties(
    trajectories: &[NesMbarTrajectory],
    states: &[StatePoint],
    options: &NesMbarOptions,
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
        let rel = estimate_relative_free_energies(&resample, states, options)?;
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

fn state_diagnostics(
    trajectories: &[NesMbarTrajectory],
    state: StatePoint,
    sample_stride: usize,
    target_idx: usize,
) -> Result<NesMbarStateDiagnostics> {
    let mut log_weights = Vec::new();
    let mut slice_totals = Vec::<f64>::new();
    let mut trajectory_totals = vec![0.0; trajectories.len()];

    for (trajectory_idx, trajectory) in trajectories.iter().enumerate() {
        for (slice_idx, sample) in trajectory
            .samples()
            .iter()
            .step_by(sample_stride)
            .enumerate()
        {
            let log_weight = -sample.reduced_work()
                - (sample.reduced_energies_states()[target_idx] - sample.reduced_energy_protocol());
            if !log_weight.is_finite() {
                continue;
            }
            log_weights.push((trajectory_idx, slice_idx, log_weight));
        }
    }

    if log_weights.is_empty() {
        return Err(CoreError::NonFiniteValue(
            "NES MBAR diagnostics had no finite nonequilibrium weights".to_string(),
        ));
    }

    let max_log = log_weights
        .iter()
        .map(|(_, _, log_weight)| *log_weight)
        .fold(f64::NEG_INFINITY, f64::max);
    let mut total = 0.0;
    for (trajectory_idx, slice_idx, log_weight) in log_weights {
        let weight = (log_weight - max_log).exp();
        total += weight;
        if slice_totals.len() <= slice_idx {
            slice_totals.resize(slice_idx + 1, 0.0);
        }
        slice_totals[slice_idx] += weight;
        trajectory_totals[trajectory_idx] += weight;
    }

    let mut sample_weight_squares = 0.0;
    let mut per_sample_weights = Vec::new();
    for trajectory in trajectories {
        for sample in trajectory.samples().iter().step_by(sample_stride) {
            let log_weight = -sample.reduced_work()
                - (sample.reduced_energies_states()[target_idx] - sample.reduced_energy_protocol());
            if log_weight.is_finite() {
                per_sample_weights.push((log_weight - max_log).exp() / total);
            }
        }
    }
    for weight in &per_sample_weights {
        sample_weight_squares += weight * weight;
    }
    let max_sample_fraction = per_sample_weights.iter().copied().fold(0.0, f64::max);
    let effective_sample_size = if sample_weight_squares > 0.0 {
        1.0 / sample_weight_squares
    } else {
        0.0
    };

    Ok(NesMbarStateDiagnostics {
        state,
        effective_sample_size,
        max_sample_fraction,
        top_slices: top_contributions(&slice_totals, total, 5),
        top_trajectories: top_contributions(&trajectory_totals, total, 5),
    })
}

fn top_contributions(values: &[f64], total: f64, limit: usize) -> Vec<NesMbarContribution> {
    let mut contributions = values
        .iter()
        .enumerate()
        .filter_map(|(index, value)| {
            (*value > 0.0).then_some(NesMbarContribution {
                index,
                fraction: *value / total,
            })
        })
        .collect::<Vec<_>>();
    contributions.sort_by(|left, right| right.fraction.total_cmp(&left.fraction));
    contributions.truncate(limit);
    contributions
}
