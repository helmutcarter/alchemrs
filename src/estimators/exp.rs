use crate::analysis::{self, BlockEstimate};
use crate::data::{find_state_index_exact, DeltaFMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};

use super::common::{ensure_consistent_lambda_labels, sample_variance, work_values};

#[derive(Debug, Clone, Default)]
pub struct IexpOptions {
    pub parallel: bool,
}

#[derive(Debug, Clone, Default)]
pub struct IexpEstimator {
    pub options: IexpOptions,
}

#[derive(Debug, Clone)]
pub struct IexpFit {
    values: Vec<f64>,
    uncertainties: Vec<f64>,
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
}

impl IexpEstimator {
    pub fn new(options: IexpOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<IexpFit> {
        if windows.len() < 2 {
            return Err(CoreError::InvalidShape {
                expected: 2,
                found: windows.len(),
            });
        }

        let states = sorted_sampled_states(windows)?;
        let lambda_labels = ensure_consistent_lambda_labels(windows)?;
        let mut window_by_index: Vec<Option<&UNkMatrix>> = vec![None; states.len()];

        for window in windows {
            let sampled = window.sampled_state().ok_or_else(|| {
                CoreError::InvalidState("sampled_state required for EXP".to_string())
            })?;
            let idx = find_state_index_exact(&states, sampled)?;
            if window_by_index[idx].is_some() {
                return Err(CoreError::InvalidState(
                    "multiple windows for same sampled_state".to_string(),
                ));
            }
            window_by_index[idx] = Some(window);
        }

        let pair_results: Vec<Result<(f64, f64, f64, f64)>> = if self.options.parallel {
            use rayon::prelude::*;
            (0..(states.len() - 1))
                .into_par_iter()
                .map(|idx| {
                    let win_f = window_by_index[idx].ok_or_else(|| {
                        CoreError::InvalidState(format!("missing window for state {idx}"))
                    })?;
                    let win_r = window_by_index[idx + 1].ok_or_else(|| {
                        CoreError::InvalidState(format!("missing window for state {}", idx + 1))
                    })?;
                    let (idx_f0, idx_f1) = adjacent_pair_indices(win_f, &states, idx)?;
                    let (idx_r0, idx_r1) = adjacent_pair_indices(win_r, &states, idx)?;
                    let work_f = work_values(win_f, idx_f0, idx_f1)?;
                    let work_r = work_values(win_r, idx_r1, idx_r0)?;
                    Ok((
                        exp_delta_f(&work_f)?,
                        exp_uncertainty(&work_f)?,
                        exp_delta_f(&work_r)?,
                        exp_uncertainty(&work_r)?,
                    ))
                })
                .collect::<Vec<_>>()
        } else {
            let mut out = Vec::with_capacity(states.len() - 1);
            for idx in 0..(states.len() - 1) {
                let win_f = window_by_index[idx].ok_or_else(|| {
                    CoreError::InvalidState(format!("missing window for state {idx}"))
                })?;
                let win_r = window_by_index[idx + 1].ok_or_else(|| {
                    CoreError::InvalidState(format!("missing window for state {}", idx + 1))
                })?;
                let (idx_f0, idx_f1) = adjacent_pair_indices(win_f, &states, idx)?;
                let (idx_r0, idx_r1) = adjacent_pair_indices(win_r, &states, idx)?;
                let work_f = work_values(win_f, idx_f0, idx_f1)?;
                let work_r = work_values(win_r, idx_r1, idx_r0)?;
                out.push(Ok((
                    exp_delta_f(&work_f)?,
                    exp_uncertainty(&work_f)?,
                    exp_delta_f(&work_r)?,
                    exp_uncertainty(&work_r)?,
                )));
            }
            out
        };

        let mut forward_delta = Vec::with_capacity(states.len() - 1);
        let mut forward_sigma = Vec::with_capacity(states.len() - 1);
        let mut reverse_delta = Vec::with_capacity(states.len() - 1);
        let mut reverse_sigma = Vec::with_capacity(states.len() - 1);
        for result in pair_results {
            let (df_f, sigma_f, df_r, sigma_r) = result?;
            forward_delta.push(df_f);
            forward_sigma.push(sigma_f);
            reverse_delta.push(df_r);
            reverse_sigma.push(sigma_r);
        }

        let n_states = states.len();
        let mut values = vec![0.0; n_states * n_states];
        let mut uncertainties = vec![0.0; n_states * n_states];

        for span in 0..forward_delta.len() {
            for start in 0..(forward_delta.len() - span) {
                let end = start + span + 1;
                let forward = forward_delta[start..end].iter().sum::<f64>();
                let reverse = reverse_delta[start..end].iter().sum::<f64>();
                let forward_var = forward_sigma[start..end]
                    .iter()
                    .map(|sigma| sigma * sigma)
                    .sum::<f64>();
                let reverse_var = reverse_sigma[start..end]
                    .iter()
                    .map(|sigma| sigma * sigma)
                    .sum::<f64>();

                values[start * n_states + end] = forward;
                uncertainties[start * n_states + end] = forward_var.sqrt();
                values[end * n_states + start] = reverse;
                uncertainties[end * n_states + start] = reverse_var.sqrt();
            }
        }

        Ok(IexpFit {
            values,
            uncertainties,
            states,
            lambda_labels,
        })
    }

    pub fn estimate(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        self.fit(windows)?.result()
    }

    pub fn estimate_with_uncertainty(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        self.fit(windows)?.result_with_uncertainty()
    }

    pub fn block_average(
        &self,
        windows: &[UNkMatrix],
        n_blocks: usize,
    ) -> Result<Vec<BlockEstimate>> {
        analysis::exp_block_average(windows, n_blocks, Some(self.options.clone()))
    }

    pub fn reverse_block_average(
        &self,
        windows: &[UNkMatrix],
        n_blocks: usize,
    ) -> Result<Vec<BlockEstimate>> {
        analysis::dexp_block_average(windows, n_blocks, Some(self.options.clone()))
    }
}

impl IexpFit {
    pub fn n_states(&self) -> usize {
        self.states.len()
    }

    pub fn states(&self) -> &[StatePoint] {
        &self.states
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }

    pub fn result(&self) -> Result<DeltaFMatrix> {
        DeltaFMatrix::new_with_labels(
            self.values.clone(),
            None,
            self.n_states(),
            self.states.clone(),
            self.lambda_labels.clone(),
        )
    }

    pub fn result_with_uncertainty(&self) -> Result<DeltaFMatrix> {
        DeltaFMatrix::new_with_labels(
            self.values.clone(),
            Some(self.uncertainties.clone()),
            self.n_states(),
            self.states.clone(),
            self.lambda_labels.clone(),
        )
    }
}

fn sorted_sampled_states(windows: &[UNkMatrix]) -> Result<Vec<StatePoint>> {
    let first = windows[0]
        .sampled_state()
        .ok_or_else(|| CoreError::InvalidState("sampled_state required for EXP".to_string()))?;
    let n_lambda = first.lambdas().len();
    let temperature = first.temperature_k();

    let mut states = Vec::with_capacity(windows.len());
    for window in windows {
        let sampled = window
            .sampled_state()
            .ok_or_else(|| CoreError::InvalidState("sampled_state required for EXP".to_string()))?;
        if sampled.lambdas().len() != n_lambda {
            return Err(CoreError::InvalidShape {
                expected: n_lambda,
                found: sampled.lambdas().len(),
            });
        }
        if sampled.temperature_k() != temperature {
            return Err(CoreError::InvalidState(
                "sampled_state temperatures differ between windows".to_string(),
            ));
        }
        if find_state_index_exact(&states, sampled).is_ok() {
            return Err(CoreError::InvalidState(
                "multiple windows for same sampled_state".to_string(),
            ));
        }
        states.push(sampled.clone());
    }

    states.sort_by(compare_state_points);
    Ok(states)
}

fn adjacent_pair_indices(
    window: &UNkMatrix,
    states: &[StatePoint],
    idx: usize,
) -> Result<(usize, usize)> {
    let from_idx =
        find_state_index_exact(window.evaluated_states(), &states[idx]).map_err(|_| {
            CoreError::InvalidState(format!(
                "window for state {idx} is missing adjacent state {} in evaluated_states",
                idx + 1
            ))
        })?;
    let to_idx =
        find_state_index_exact(window.evaluated_states(), &states[idx + 1]).map_err(|_| {
            CoreError::InvalidState(format!(
                "window for state {idx} is missing adjacent state {} in evaluated_states",
                idx + 1
            ))
        })?;
    Ok((from_idx, to_idx))
}

fn compare_state_points(left: &StatePoint, right: &StatePoint) -> std::cmp::Ordering {
    for (lhs, rhs) in left.lambdas().iter().zip(right.lambdas().iter()) {
        let ordering = lhs.total_cmp(rhs);
        if ordering != std::cmp::Ordering::Equal {
            return ordering;
        }
    }
    left.lambdas().len().cmp(&right.lambdas().len())
}

fn logsumexp(values: &[f64]) -> f64 {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    let sum = values.iter().map(|v| (v - max).exp()).sum::<f64>();
    max + sum.ln()
}

fn exp_delta_f(work: &[f64]) -> Result<f64> {
    if work
        .iter()
        .any(|value| *value == f64::NEG_INFINITY || value.is_nan())
    {
        return Err(CoreError::NonFiniteValue(
            "EXP work values must be finite or +inf".to_string(),
        ));
    }
    let t = work.len() as f64;
    let neg: Vec<f64> = work.iter().map(|w| -w).collect();
    let delta_f = -(logsumexp(&neg) - t.ln());
    if !delta_f.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "EXP delta_f became non-finite; no finite Boltzmann weight remained".to_string(),
        ));
    }
    Ok(delta_f)
}

fn exp_uncertainty(work: &[f64]) -> Result<f64> {
    if work.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    if work
        .iter()
        .any(|value| *value == f64::NEG_INFINITY || value.is_nan())
    {
        return Err(CoreError::NonFiniteValue(
            "EXP work values must be finite or +inf".to_string(),
        ));
    }
    let t = work.len() as f64;
    let mut max_arg = f64::NEG_INFINITY;
    for value in work {
        let arg = -value;
        if arg > max_arg {
            max_arg = arg;
        }
    }
    let mut x = Vec::with_capacity(work.len());
    for value in work {
        x.push((-value - max_arg).exp());
    }
    let mean = x.iter().sum::<f64>() / t;
    if mean == 0.0 || !mean.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "EXP uncertainty became non-finite; no finite Boltzmann weight remained".to_string(),
        ));
    }
    let Some(var) = sample_variance(&x) else {
        return Ok(f64::NAN);
    };
    let std = var.sqrt();
    let dx = std / (t.sqrt());
    let uncertainty = dx / mean;
    if !uncertainty.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "EXP uncertainty became non-finite".to_string(),
        ));
    }
    Ok(uncertainty)
}
