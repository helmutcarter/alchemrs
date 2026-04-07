use crate::analysis::{self, BlockEstimate};
use crate::data::{find_state_index_exact, DeltaFMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};

use super::common::{ensure_consistent_lambda_labels, work_values, PairEstimate};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarMethod {
    FalsePosition,
    SelfConsistentIteration,
    Bisection,
}

#[derive(Debug, Clone)]
pub struct BarOptions {
    pub maximum_iterations: usize,
    pub relative_tolerance: f64,
    pub method: BarMethod,
    pub parallel: bool,
}

impl Default for BarOptions {
    fn default() -> Self {
        Self {
            maximum_iterations: 10_000,
            relative_tolerance: 1.0e-7,
            method: BarMethod::FalsePosition,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct BarEstimator {
    pub options: BarOptions,
}

#[derive(Debug, Clone)]
pub struct BarFit {
    values: Vec<f64>,
    uncertainties: Vec<f64>,
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
}

impl BarEstimator {
    pub fn new(options: BarOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<BarFit> {
        if windows.len() < 2 {
            return Err(CoreError::InvalidShape {
                expected: 2,
                found: windows.len(),
            });
        }
        let eval_states = sorted_sampled_states(windows)?;
        let lambda_labels = ensure_consistent_lambda_labels(windows)?;
        let mut window_by_index: Vec<Option<&UNkMatrix>> = vec![None; eval_states.len()];
        for window in windows {
            let sampled = window.sampled_state().ok_or_else(|| {
                CoreError::InvalidState("sampled_state required for BAR".to_string())
            })?;
            let idx = find_state_index_exact(&eval_states, sampled)?;
            if window_by_index[idx].is_some() {
                return Err(CoreError::InvalidState(
                    "multiple windows for same sampled_state".to_string(),
                ));
            }
            window_by_index[idx] = Some(window);
        }

        let pair_results: Vec<Result<PairEstimate>> = if self.options.parallel {
            use rayon::prelude::*;
            (0..(eval_states.len() - 1))
                .into_par_iter()
                .map(|idx| {
                    let win_f = window_by_index[idx].ok_or_else(|| {
                        CoreError::InvalidState(format!("missing window for state {idx}"))
                    })?;
                    let win_r = window_by_index[idx + 1].ok_or_else(|| {
                        CoreError::InvalidState(format!("missing window for state {}", idx + 1))
                    })?;
                    let (idx_f0, idx_f1) = adjacent_pair_indices(win_f, &eval_states, idx)?;
                    let (idx_r0, idx_r1) = adjacent_pair_indices(win_r, &eval_states, idx)?;
                    let w_f = work_values(win_f, idx_f0, idx_f1)?;
                    let w_r = work_values(win_r, idx_r0, idx_r1)?
                        .into_iter()
                        .map(|w| -w)
                        .collect::<Vec<_>>();
                    let (df, ddf) = bar_estimate(
                        &w_f,
                        &w_r,
                        self.options.method,
                        self.options.maximum_iterations,
                        self.options.relative_tolerance,
                    )?;
                    Ok((df, ddf))
                })
                .collect::<Vec<_>>()
        } else {
            let mut out = Vec::new();
            for idx in 0..(eval_states.len() - 1) {
                let win_f = window_by_index[idx].ok_or_else(|| {
                    CoreError::InvalidState(format!("missing window for state {idx}"))
                })?;
                let win_r = window_by_index[idx + 1].ok_or_else(|| {
                    CoreError::InvalidState(format!("missing window for state {}", idx + 1))
                })?;
                let (idx_f0, idx_f1) = adjacent_pair_indices(win_f, &eval_states, idx)?;
                let (idx_r0, idx_r1) = adjacent_pair_indices(win_r, &eval_states, idx)?;
                let w_f = work_values(win_f, idx_f0, idx_f1)?;
                let w_r = work_values(win_r, idx_r0, idx_r1)?
                    .into_iter()
                    .map(|w| -w)
                    .collect::<Vec<_>>();
                let (df, ddf) = bar_estimate(
                    &w_f,
                    &w_r,
                    self.options.method,
                    self.options.maximum_iterations,
                    self.options.relative_tolerance,
                )?;
                out.push(Ok((df, ddf)));
            }
            out
        };

        let mut deltas = Vec::with_capacity(pair_results.len());
        let mut d_deltas = Vec::with_capacity(pair_results.len());
        for result in pair_results {
            let (delta, d_delta) = result?;
            deltas.push(delta);
            d_deltas.push(d_delta);
        }

        let n_states = eval_states.len();
        let mut adelta = vec![0.0; n_states * n_states];
        let mut ad_delta = vec![f64::NAN; n_states * n_states];

        for j in 0..deltas.len() {
            let mut out = Vec::new();
            let mut dout = Vec::new();
            for i in 0..(deltas.len() - j) {
                out.push(deltas[i..=i + j].iter().sum::<f64>());
                if j == 0 {
                    dout.push(d_deltas[i..=i + j].iter().sum::<f64>());
                } else {
                    dout.push(f64::NAN);
                }
            }
            for (i, value) in out.into_iter().enumerate() {
                let row = i;
                let col = i + j + 1;
                adelta[row * n_states + col] = value;
            }
            for (i, value) in dout.into_iter().enumerate() {
                let row = i;
                let col = i + j + 1;
                ad_delta[row * n_states + col] = value;
            }
        }

        for i in 0..n_states {
            for j in (i + 1)..n_states {
                let val = adelta[i * n_states + j];
                adelta[j * n_states + i] = -val;
                let unc = ad_delta[i * n_states + j];
                ad_delta[j * n_states + i] = unc;
            }
        }

        Ok(BarFit {
            values: adelta,
            uncertainties: ad_delta,
            states: eval_states,
            lambda_labels,
        })
    }

    pub fn estimate(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        self.fit(windows)?.result()
    }

    pub fn block_average(
        &self,
        windows: &[UNkMatrix],
        n_blocks: usize,
    ) -> Result<Vec<BlockEstimate>> {
        analysis::bar_block_average(windows, n_blocks, Some(self.options.clone()))
    }
}

fn sorted_sampled_states(windows: &[UNkMatrix]) -> Result<Vec<StatePoint>> {
    let first = windows[0]
        .sampled_state()
        .ok_or_else(|| CoreError::InvalidState("sampled_state required for BAR".to_string()))?;
    let n_lambda = first.lambdas().len();
    let temperature = first.temperature_k();

    let mut states = Vec::with_capacity(windows.len());
    for window in windows {
        let sampled = window
            .sampled_state()
            .ok_or_else(|| CoreError::InvalidState("sampled_state required for BAR".to_string()))?;
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

impl BarFit {
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
            Some(self.uncertainties.clone()),
            self.n_states(),
            self.states.clone(),
            self.lambda_labels.clone(),
        )
    }
}

fn exp_delta(w_f: &[f64]) -> Result<f64> {
    let log_sum = logsumexp(&w_f.iter().map(|w| -w).collect::<Vec<_>>());
    let t = w_f.len() as f64;
    let delta = -(log_sum - t.ln());
    if !delta.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "BAR initial estimate became non-finite".to_string(),
        ));
    }
    Ok(delta)
}

fn bar_zero(w_f: &[f64], w_r: &[f64], delta_f: f64) -> f64 {
    let t_f = w_f.len() as f64;
    let t_r = w_r.len() as f64;
    let m = (t_f / t_r).ln();

    let log_f_f: Vec<f64> = w_f.iter().map(|w| neg_log1pexp(m + w - delta_f)).collect();
    let log_numer = logsumexp(&log_f_f);

    let log_f_r: Vec<f64> = w_r
        .iter()
        .map(|w| neg_log1pexp(-(m - w - delta_f)))
        .collect();
    let log_denom = logsumexp(&log_f_r);

    log_numer - log_denom
}

fn neg_log1pexp(value: f64) -> f64 {
    if value == f64::INFINITY {
        f64::NEG_INFINITY
    } else if value == f64::NEG_INFINITY {
        0.0
    } else if value > 0.0 {
        -value - (-value).exp().ln_1p()
    } else {
        -value.exp().ln_1p()
    }
}

pub(crate) fn bar_estimate(
    w_f: &[f64],
    w_r: &[f64],
    method: BarMethod,
    maximum_iterations: usize,
    relative_tolerance: f64,
) -> Result<(f64, f64)> {
    if w_f.is_empty() || w_r.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    if w_f.iter().any(|v| v.is_nan()) || w_r.iter().any(|v| v.is_nan()) {
        return Err(CoreError::NonFiniteValue(
            "BAR work values must not contain NaN".to_string(),
        ));
    }
    if !w_f.iter().any(|v| v.is_finite()) || !w_r.iter().any(|v| v.is_finite()) {
        return Err(CoreError::NonFiniteValue(
            "BAR requires at least one finite work value in each direction".to_string(),
        ));
    }
    let mut delta_f = 0.0;
    let initial_upper = exp_delta(w_f)?;
    let initial_lower = -exp_delta(w_r)?;
    let (mut lower, mut upper, mut f_lower, mut f_upper) =
        bracket_bar_root(w_f, w_r, initial_lower, initial_upper, maximum_iterations)?;

    if !f_upper.is_finite() || !f_lower.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "BAR overlap function became non-finite while bracketing the solution".to_string(),
        ));
    }

    let mut converged = false;
    for _ in 0..=maximum_iterations {
        let old = delta_f;
        match method {
            BarMethod::FalsePosition => {
                delta_f = if upper == 0.0 && lower == 0.0 {
                    0.0
                } else {
                    upper - f_upper * (upper - lower) / (f_upper - f_lower)
                };
            }
            BarMethod::Bisection => {
                delta_f = (upper + lower) / 2.0;
            }
            BarMethod::SelfConsistentIteration => {
                delta_f = -bar_zero(w_f, w_r, delta_f) + delta_f;
            }
        }
        let f_new = bar_zero(w_f, w_r, delta_f);
        if f_new.abs() < relative_tolerance {
            converged = true;
            break;
        }
        if delta_f != 0.0 {
            let relative_change = ((delta_f - old) / delta_f).abs();
            if relative_change < relative_tolerance {
                converged = true;
                break;
            }
        }

        if matches!(method, BarMethod::FalsePosition | BarMethod::Bisection) {
            if f_upper * f_new < 0.0 {
                lower = delta_f;
                f_lower = f_new;
            } else if f_lower * f_new <= 0.0 {
                upper = delta_f;
                f_upper = f_new;
            } else {
                return Err(CoreError::ConvergenceFailure);
            }
        }
    }

    if !converged {
        return Err(CoreError::ConvergenceFailure);
    }

    let d_delta_f = bar_uncertainty(w_f, w_r, delta_f)?;
    Ok((delta_f, d_delta_f))
}

fn bracket_bar_root(
    w_f: &[f64],
    w_r: &[f64],
    initial_lower: f64,
    initial_upper: f64,
    maximum_iterations: usize,
) -> Result<(f64, f64, f64, f64)> {
    let mut lower = initial_lower.min(initial_upper);
    let mut upper = initial_lower.max(initial_upper);
    let mut f_lower = bar_zero(w_f, w_r, lower);
    let mut f_upper = bar_zero(w_f, w_r, upper);

    if !f_lower.is_finite() || !f_upper.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "BAR overlap function became non-finite while bracketing the solution".to_string(),
        ));
    }
    if f_lower == 0.0 || f_upper == 0.0 || f_lower * f_upper < 0.0 {
        return Ok((lower, upper, f_lower, f_upper));
    }

    let mut expansion = (upper - lower).abs().max(0.1);
    for _ in 0..maximum_iterations {
        lower -= expansion;
        upper += expansion;
        f_lower = bar_zero(w_f, w_r, lower);
        f_upper = bar_zero(w_f, w_r, upper);
        if !f_lower.is_finite() || !f_upper.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "BAR overlap function became non-finite while bracketing the solution".to_string(),
            ));
        }
        if f_lower == 0.0 || f_upper == 0.0 || f_lower * f_upper < 0.0 {
            return Ok((lower, upper, f_lower, f_upper));
        }
        expansion *= 2.0;
    }

    Err(CoreError::ConvergenceFailure)
}

fn bar_uncertainty(w_f: &[f64], w_r: &[f64], delta_f: f64) -> Result<f64> {
    if w_f.iter().any(|v| !v.is_finite()) || w_r.iter().any(|v| !v.is_finite()) {
        return Ok(f64::NAN);
    }
    let t_f = w_f.len() as f64;
    let t_r = w_r.len() as f64;
    let m = (t_f / t_r).ln();
    let c = m - delta_f;

    let exp_arg_f: Vec<f64> = w_f.iter().map(|w| w + c).collect();
    let max_arg_f = exp_arg_f.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let log_f_f: Vec<f64> = exp_arg_f
        .iter()
        .map(|arg| -(((-max_arg_f).exp() + (arg - max_arg_f).exp()).ln()))
        .collect();
    let af_f = (logsumexp(&log_f_f) - max_arg_f).exp() / t_f;
    let af_f2 = (logsumexp(&log_f_f.iter().map(|v| 2.0 * v).collect::<Vec<_>>()) - 2.0 * max_arg_f)
        .exp()
        / t_f;

    let exp_arg_r: Vec<f64> = w_r.iter().map(|w| w - c).collect();
    let max_arg_r = exp_arg_r.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let log_f_r: Vec<f64> = exp_arg_r
        .iter()
        .map(|arg| -(((-max_arg_r).exp() + (arg - max_arg_r).exp()).ln()))
        .collect();
    let af_r = (logsumexp(&log_f_r) - max_arg_r).exp() / t_r;
    let af_r2 = (logsumexp(&log_f_r.iter().map(|v| 2.0 * v).collect::<Vec<_>>()) - 2.0 * max_arg_r)
        .exp()
        / t_r;

    let nrat = (t_f + t_r) / (t_f * t_r);
    let variance = (af_f2 / (af_f * af_f)) / t_f + (af_r2 / (af_r * af_r)) / t_r - nrat;
    Ok(variance.max(0.0).sqrt())
}

fn logsumexp(values: &[f64]) -> f64 {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    let sum = values.iter().map(|v| (v - max).exp()).sum::<f64>();
    max + sum.ln()
}
