use crate::analysis::{self, BlockEstimate};
use crate::data::{find_state_index_exact, DeltaFMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};

use super::common::{
    checked_nonnegative_sqrt_or_nan, ensure_consistent_lambda_labels, sample_covariance,
    sample_variance, work_values,
};

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

#[derive(Debug, Clone)]
struct BarEdgeEstimate {
    delta_f: f64,
    variance: f64,
    sigma: f64,
    derivative: f64,
    forward_weights: Vec<f64>,
    reverse_weights: Vec<f64>,
}

#[derive(Debug, Clone)]
struct BarUncertaintyEstimate {
    variance: f64,
    sigma: f64,
    derivative: f64,
    forward_weights: Vec<f64>,
    reverse_weights: Vec<f64>,
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

        let pair_results: Vec<Result<BarEdgeEstimate>> = if self.options.parallel {
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
                    bar_edge_estimate(
                        &w_f,
                        &w_r,
                        self.options.method,
                        self.options.maximum_iterations,
                        self.options.relative_tolerance,
                    )
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
                out.push(bar_edge_estimate(
                    &w_f,
                    &w_r,
                    self.options.method,
                    self.options.maximum_iterations,
                    self.options.relative_tolerance,
                ));
            }
            out
        };

        let mut edge_estimates = Vec::with_capacity(pair_results.len());
        for result in pair_results {
            edge_estimates.push(result?);
        }
        let deltas = edge_estimates
            .iter()
            .map(|estimate| estimate.delta_f)
            .collect::<Vec<_>>();
        let mut adjacent_covariances = Vec::with_capacity(edge_estimates.len().saturating_sub(1));
        for pair in edge_estimates.windows(2) {
            adjacent_covariances.push(adjacent_bar_covariance(&pair[0], &pair[1])?);
        }

        let n_states = eval_states.len();
        let mut adelta = vec![0.0; n_states * n_states];
        let mut ad_delta = vec![0.0; n_states * n_states];

        for j in 0..deltas.len() {
            let mut out = Vec::new();
            let mut dout = Vec::new();
            for i in 0..(deltas.len() - j) {
                out.push(deltas[i..=i + j].iter().sum::<f64>());
                if j == 0 {
                    dout.push(edge_estimates[i].sigma);
                } else {
                    dout.push(path_uncertainty(
                        &edge_estimates,
                        &adjacent_covariances,
                        i,
                        i + j,
                    ));
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

#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn bar_estimate(
    w_f: &[f64],
    w_r: &[f64],
    method: BarMethod,
    maximum_iterations: usize,
    relative_tolerance: f64,
) -> Result<(f64, f64)> {
    let estimate = bar_edge_estimate(w_f, w_r, method, maximum_iterations, relative_tolerance)?;
    Ok((estimate.delta_f, estimate.sigma))
}

fn bar_edge_estimate(
    w_f: &[f64],
    w_r: &[f64],
    method: BarMethod,
    maximum_iterations: usize,
    relative_tolerance: f64,
) -> Result<BarEdgeEstimate> {
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

    let uncertainty = bar_uncertainty(w_f, w_r, delta_f)?;
    Ok(BarEdgeEstimate {
        delta_f,
        variance: uncertainty.variance,
        sigma: uncertainty.sigma,
        derivative: uncertainty.derivative,
        forward_weights: uncertainty.forward_weights,
        reverse_weights: uncertainty.reverse_weights,
    })
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

fn bar_uncertainty(w_f: &[f64], w_r: &[f64], delta_f: f64) -> Result<BarUncertaintyEstimate> {
    let t_f = w_f.len() as f64;
    let t_r = w_r.len() as f64;
    let m = (t_f / t_r).ln();

    let mut forward_weights = Vec::with_capacity(w_f.len());
    let mut reverse_weights = Vec::with_capacity(w_r.len());
    let mut derivative = 0.0;

    for &work in w_f {
        let weight = fermi_weight(m + work - delta_f);
        derivative += weight * (1.0 - weight) / t_f;
        forward_weights.push(weight);
    }
    for &work in w_r {
        let weight = fermi_weight(-m + work + delta_f);
        derivative += weight * (1.0 - weight) / t_r;
        reverse_weights.push(weight);
    }

    let variance = match (
        sample_variance(&forward_weights),
        sample_variance(&reverse_weights),
    ) {
        (Some(var_f), Some(var_r)) if derivative.is_finite() && derivative > 0.0 => {
            (var_f / t_f + var_r / t_r) / (derivative * derivative)
        }
        _ => f64::NAN,
    };
    let sigma = if variance.is_finite() {
        checked_nonnegative_sqrt_or_nan(variance, variance.abs(), "BAR adjacent uncertainty")
    } else {
        f64::NAN
    };
    Ok(BarUncertaintyEstimate {
        variance,
        sigma,
        derivative,
        forward_weights,
        reverse_weights,
    })
}

fn adjacent_bar_covariance(left: &BarEdgeEstimate, right: &BarEdgeEstimate) -> Result<f64> {
    if !left.derivative.is_finite()
        || !right.derivative.is_finite()
        || left.derivative <= 0.0
        || right.derivative <= 0.0
    {
        return Ok(f64::NAN);
    }

    let shared_covariance = match sample_covariance(&left.reverse_weights, &right.forward_weights)?
    {
        Some(value) => value,
        None => return Ok(f64::NAN),
    };
    let n_shared = left.reverse_weights.len() as f64;
    Ok(-shared_covariance / (n_shared * left.derivative * right.derivative))
}

fn path_uncertainty(
    edge_estimates: &[BarEdgeEstimate],
    adjacent_covariances: &[f64],
    start_edge: usize,
    end_edge: usize,
) -> f64 {
    let mut variance = 0.0;
    let mut scale = 0.0;
    for edge in &edge_estimates[start_edge..=end_edge] {
        if !edge.variance.is_finite() {
            return f64::NAN;
        }
        variance += edge.variance;
        scale += edge.variance.abs();
    }
    for covariance in &adjacent_covariances[start_edge..end_edge] {
        if !covariance.is_finite() {
            return f64::NAN;
        }
        variance += 2.0 * covariance;
        scale += 2.0 * covariance.abs();
    }
    checked_nonnegative_sqrt_or_nan(variance, scale, "BAR cumulative uncertainty")
}

fn fermi_weight(argument: f64) -> f64 {
    if argument >= 0.0 {
        let exp_neg = (-argument).exp();
        exp_neg / (1.0 + exp_neg)
    } else {
        1.0 / (1.0 + argument.exp())
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

#[cfg(test)]
mod tests {
    use super::{
        adjacent_bar_covariance, bar_edge_estimate, path_uncertainty, BarEdgeEstimate, BarMethod,
        BarOptions,
    };

    fn edge(w_f: &[f64], w_r: &[f64]) -> BarEdgeEstimate {
        bar_edge_estimate(
            w_f,
            w_r,
            BarMethod::FalsePosition,
            BarOptions::default().maximum_iterations,
            BarOptions::default().relative_tolerance,
        )
        .expect("BAR edge estimate")
    }

    #[test]
    fn cumulative_uncertainty_includes_neighbor_covariance() {
        let left = edge(&[0.0, 0.4, 0.8, 1.2], &[-0.8, -0.2, 0.4, 1.0]);
        let right = edge(&[-0.6, 0.0, 0.6, 1.2], &[0.2, 0.5, 0.8, 1.1]);

        let covariance = adjacent_bar_covariance(&left, &right).expect("adjacent covariance");
        assert!(covariance.is_finite());
        assert!(covariance.abs() > 1.0e-6);

        let sigma = path_uncertainty(&[left.clone(), right.clone()], &[covariance], 0, 1);
        let quadrature = (left.sigma * left.sigma + right.sigma * right.sigma).sqrt();
        let expected = (left.variance + right.variance + 2.0 * covariance)
            .max(0.0)
            .sqrt();

        assert!((sigma - expected).abs() < 1.0e-12);
        assert!((sigma - quadrature).abs() > 1.0e-6);
    }
}
