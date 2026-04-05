use crate::analysis::{self, BlockEstimate};
use crate::data::{find_state_index_exact, DeltaFMatrix, UNkMatrix};
use crate::error::{CoreError, Result};

use super::common::{
    ensure_consistent_lambda_labels, ensure_consistent_states, work_values,
    ExpRow,
};

#[derive(Debug, Clone)]
pub struct ExpOptions {
    pub compute_uncertainty: bool,
    pub parallel: bool,
}

impl Default for ExpOptions {
    fn default() -> Self {
        Self {
            compute_uncertainty: true,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ExpEstimator {
    pub options: ExpOptions,
}

impl ExpEstimator {
    pub fn new(options: ExpOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        if windows.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        let states = ensure_consistent_states(windows)?;
        let n_states = states.len();
        let mut window_map: Vec<Option<&UNkMatrix>> = vec![None; n_states];

        for window in windows {
            let sampled = window.sampled_state().ok_or_else(|| {
                CoreError::InvalidState("sampled_state required for EXP".to_string())
            })?;
            let idx = find_state_index_exact(&states, sampled)?;
            if window_map[idx].is_some() {
                return Err(CoreError::InvalidState(
                    "multiple windows for same sampled_state".to_string(),
                ));
            }
            window_map[idx] = Some(window);
        }

        for (idx, entry) in window_map.iter().enumerate() {
            if entry.is_none() {
                return Err(CoreError::InvalidState(format!(
                    "missing window for state {idx}"
                )));
            }
        }

        let mut values = vec![0.0; n_states * n_states];
        let mut uncertainties = if self.options.compute_uncertainty {
            Some(vec![0.0; n_states * n_states])
        } else {
            None
        };

        if self.options.parallel {
            use rayon::prelude::*;
            let rows: Vec<Result<ExpRow>> = (0..n_states)
                .into_par_iter()
                .map(|i| {
                    let window = window_map[i].expect("window present");
                    let mut row = Vec::with_capacity(n_states);
                    let mut row_unc = Vec::with_capacity(n_states);
                    for j in 0..n_states {
                        let work = work_values(window, i, j)?;
                        let delta = exp_delta_f(&work)?;
                        row.push(delta);
                        if self.options.compute_uncertainty {
                            row_unc.push(exp_uncertainty(&work)?);
                        }
                    }
                    Ok((i, row, row_unc))
                })
                .collect::<Vec<_>>();
            for entry in rows {
                let (i, row, row_unc) = entry?;
                for (j, value) in row.into_iter().enumerate() {
                    values[i * n_states + j] = value;
                }
                if let Some(ref mut unc) = uncertainties {
                    for (j, value) in row_unc.into_iter().enumerate() {
                        unc[i * n_states + j] = value;
                    }
                }
            }
        } else {
            for i in 0..n_states {
                let window = window_map[i].expect("window present");
                for j in 0..n_states {
                    let work = work_values(window, i, j)?;
                    let delta = exp_delta_f(&work)?;
                    values[i * n_states + j] = delta;
                    if let Some(ref mut unc) = uncertainties {
                        unc[i * n_states + j] = exp_uncertainty(&work)?;
                    }
                }
            }
        }

        let lambda_labels = ensure_consistent_lambda_labels(windows)?;
        DeltaFMatrix::new_with_labels(values, uncertainties, n_states, states, lambda_labels)
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
    let mut var = 0.0;
    for value in &x {
        let diff = value - mean;
        var += diff * diff;
    }
    var /= t;
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
