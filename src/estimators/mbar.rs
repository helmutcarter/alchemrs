use crate::analysis::{self, BlockEstimate};
use crate::data::{find_state_index_exact, DeltaFMatrix, OverlapMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};
use rayon::prelude::*;
use std::sync::OnceLock;

use super::common::{ensure_consistent_lambda_labels, ensure_consistent_states, CombinedWindows};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MbarSolver {
    #[default]
    FixedPoint,
    Lbfgs,
}

#[derive(Debug, Clone)]
pub struct MbarOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub initial_f_k: Option<Vec<f64>>,
    pub parallel: bool,
    pub solver: MbarSolver,
}

impl Default for MbarOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
            tolerance: 1.0e-7,
            initial_f_k: None,
            parallel: false,
            solver: MbarSolver::FixedPoint,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MbarEstimator {
    pub options: MbarOptions,
}

#[derive(Debug, Clone)]
struct PreparedMbar {
    u_kn: Vec<Vec<f64>>,
    // Keep a sample-major view for the default serial hot loops so denominator and
    // log-weight passes can walk contiguous rows without changing the parallel path.
    u_nk: Vec<f64>,
    n_k: Vec<f64>,
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
}

#[derive(Debug)]
pub struct MbarFit {
    prepared: PreparedMbar,
    f_k: Vec<f64>,
    parallel: bool,
    log_weights: OnceLock<Vec<f64>>,
}

impl MbarEstimator {
    pub fn new(options: MbarOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<MbarFit> {
        if windows.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        let prepared = PreparedMbar::from_windows(windows)?;
        let mut f_k = initial_f_k(&self.options, prepared.states.len())?;

        match self.options.solver {
            MbarSolver::FixedPoint => mbar_solve_fixed_point(
                &prepared.u_kn,
                &prepared.u_nk,
                &prepared.n_k,
                &mut f_k,
                self.options.tolerance,
                self.options.max_iterations,
                self.options.parallel,
            )?,
            MbarSolver::Lbfgs => mbar_solve_lbfgs(
                &prepared.u_nk,
                &prepared.n_k,
                &mut f_k,
                self.options.tolerance,
                self.options.max_iterations,
                self.options.parallel,
            )?,
        }

        Ok(MbarFit::new(prepared, f_k, self.options.parallel))
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
        analysis::mbar_block_average(windows, n_blocks, Some(self.options.clone()))
    }
}

impl MbarFit {
    fn new(prepared: PreparedMbar, f_k: Vec<f64>, parallel: bool) -> Self {
        Self {
            prepared,
            f_k,
            parallel,
            log_weights: OnceLock::new(),
        }
    }

    pub fn n_states(&self) -> usize {
        self.prepared.states.len()
    }

    pub fn states(&self) -> &[StatePoint] {
        &self.prepared.states
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.prepared.lambda_labels.as_deref()
    }

    pub fn free_energies(&self) -> &[f64] {
        &self.f_k
    }

    pub fn state_counts(&self) -> &[f64] {
        &self.prepared.n_k
    }

    pub fn log_weights(&self) -> &[f64] {
        self.log_weights.get_or_init(|| {
            mbar_log_weights(
                &self.prepared.u_kn,
                &self.prepared.u_nk,
                &self.prepared.n_k,
                &self.f_k,
                self.parallel,
            )
        })
    }

    pub fn result(&self) -> Result<DeltaFMatrix> {
        self.delta_f_matrix_inner(None)
    }

    pub fn result_with_uncertainty(&self) -> Result<DeltaFMatrix> {
        let uncertainties = mbar_uncertainty(
            &self.prepared.u_kn,
            &self.prepared.u_nk,
            &self.prepared.n_k,
            &self.f_k,
            self.parallel,
        )?;
        self.delta_f_matrix_inner(Some(uncertainties))
    }

    pub fn overlap_matrix(&self) -> Result<OverlapMatrix> {
        let n_states = self.prepared.states.len();
        let log_w = self.log_weights();
        if log_w.len() % n_states != 0 {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: log_w.len(),
            });
        }
        let n_samples = log_w.len() / n_states;
        let mut wtw = vec![0.0; n_states * n_states];
        let mut weights = vec![0.0; n_states];
        for n in 0..n_samples {
            let base = n * n_states;
            for i in 0..n_states {
                weights[i] = log_w[base + i].exp();
            }
            for i in 0..n_states {
                let wi = weights[i];
                let row_offset = i * n_states;
                for (j, &wj) in weights.iter().enumerate().skip(i) {
                    let idx = row_offset + j;
                    wtw[idx] = wi.mul_add(wj, wtw[idx]);
                }
            }
        }

        let mut values = vec![0.0; n_states * n_states];
        for i in 0..n_states {
            for j in 0..n_states {
                let upper_idx = if i <= j {
                    i * n_states + j
                } else {
                    j * n_states + i
                };
                values[i * n_states + j] = wtw[upper_idx] * self.prepared.n_k[j];
            }
        }

        OverlapMatrix::new_with_labels(
            values,
            n_states,
            self.prepared.states.clone(),
            self.prepared.lambda_labels.clone(),
        )
    }

    pub fn overlap_eigenvalues(&self) -> Result<Vec<f64>> {
        let overlap = self.overlap_matrix()?;
        analysis::overlap_eigenvalues(&overlap)
    }

    pub fn overlap_scalar(&self) -> Result<f64> {
        let overlap = self.overlap_matrix()?;
        analysis::overlap_scalar(&overlap)
    }

    fn delta_f_matrix_inner(&self, uncertainties: Option<Vec<f64>>) -> Result<DeltaFMatrix> {
        let n_states = self.prepared.states.len();
        let mut values = vec![0.0; n_states * n_states];
        for i in 0..n_states {
            for j in 0..n_states {
                values[i * n_states + j] = self.f_k[j] - self.f_k[i];
            }
        }
        DeltaFMatrix::new_with_labels(
            values,
            uncertainties,
            n_states,
            self.prepared.states.clone(),
            self.prepared.lambda_labels.clone(),
        )
    }
}

impl PreparedMbar {
    fn from_windows(windows: &[UNkMatrix]) -> Result<Self> {
        let (u_kn, n_k, states, lambda_labels) = combine_windows(windows)?;
        let u_nk = sample_major_from_state_major(&u_kn);
        Ok(Self {
            u_kn,
            u_nk,
            n_k,
            states,
            lambda_labels,
        })
    }
}

fn initial_f_k(options: &MbarOptions, n_states: usize) -> Result<Vec<f64>> {
    if let Some(ref initial) = options.initial_f_k {
        if initial.len() != n_states {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: initial.len(),
            });
        }
        Ok(initial.clone())
    } else {
        Ok(vec![0.0; n_states])
    }
}

fn combine_windows(windows: &[UNkMatrix]) -> Result<CombinedWindows> {
    let states = ensure_consistent_states(windows)?;
    let lambda_labels = ensure_consistent_lambda_labels(windows)?;
    let n_states = states.len();
    let mut n_k = vec![0.0; n_states];
    let mut offsets = Vec::with_capacity(windows.len());
    let mut total_samples = 0usize;

    for window in windows {
        let sampled = window.sampled_state().ok_or_else(|| {
            CoreError::InvalidState("sampled_state required for MBAR".to_string())
        })?;
        let idx = find_state_index_exact(&states, sampled)?;
        n_k[idx] += window.n_samples() as f64;
        offsets.push(total_samples);
        total_samples += window.n_samples();
    }

    let mut u_kn = vec![vec![0.0; total_samples]; n_states];
    for (win_idx, window) in windows.iter().enumerate() {
        let offset = offsets[win_idx];
        let data = window.data();
        for sample_idx in 0..window.n_samples() {
            let row_offset = sample_idx * n_states;
            let out_idx = offset + sample_idx;
            for k in 0..n_states {
                u_kn[k][out_idx] = data[row_offset + k];
            }
        }
    }

    validate_mbar_input(&u_kn)?;

    Ok((u_kn, n_k, states, lambda_labels))
}

fn mbar_solve_fixed_point(
    u_kn: &[Vec<f64>],
    u_nk: &[f64],
    n_k: &[f64],
    f_k: &mut [f64],
    tolerance: f64,
    max_iterations: usize,
    parallel: bool,
) -> Result<()> {
    let n_states = u_kn.len();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let n_samples = u_kn[0].len();
    let ln_n_k: Vec<f64> = n_k
        .iter()
        .map(|n| if *n > 0.0 { n.ln() } else { f64::NEG_INFINITY })
        .collect();

    let mut log_denominator = vec![0.0; n_samples];
    let mut state_offsets = vec![0.0; n_states];
    let mut f_new = vec![0.0; n_states];
    for _ in 0..max_iterations {
        // Precompute ln(n_k) + f_k once per iteration so the sample loop only pays the
        // u_kn access and the log-sum-exp update.
        fill_state_offsets(&ln_n_k, f_k, &mut state_offsets);
        if parallel {
            fill_log_denominator(u_kn, &state_offsets, &mut log_denominator, true)?;
            fill_free_energies(u_kn, &log_denominator, &mut f_new, true)?;
        } else {
            // The serial path uses the sample-major buffer to keep per-sample scans contiguous.
            fill_log_denominator_serial(u_nk, n_states, &state_offsets, &mut log_denominator)?;
            fill_free_energies_serial(u_nk, n_states, &log_denominator, &mut f_new)?;
        }
        let shift = f_new[0];
        let mut max_delta = 0.0;
        for (new_value, old_value) in f_new.iter_mut().zip(f_k.iter()) {
            *new_value -= shift;
            let delta = (*new_value - *old_value).abs();
            if delta > max_delta {
                max_delta = delta;
            }
        }
        f_k.copy_from_slice(&f_new);
        if max_delta < tolerance {
            return Ok(());
        }
    }

    Err(CoreError::ConvergenceFailure)
}

#[derive(Debug)]
struct LbfgsHistoryEntry {
    s: Vec<f64>,
    y: Vec<f64>,
    rho: f64,
}

fn mbar_solve_lbfgs(
    u_nk: &[f64],
    n_k: &[f64],
    f_k: &mut [f64],
    tolerance: f64,
    max_iterations: usize,
    parallel: bool,
) -> Result<()> {
    const LBFGS_HISTORY_LIMIT: usize = 10;
    const LBFGS_C1: f64 = 1.0e-4;
    const LBFGS_MIN_STEP: f64 = 1.0e-20;
    const LBFGS_MAX_LINE_SEARCH_STEPS: usize = 40;
    const LBFGS_CURVATURE_EPS: f64 = 1.0e-12;

    let n_states = f_k.len();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    if n_states == 1 {
        f_k[0] = 0.0;
        return Ok(());
    }

    let ln_n_k = positive_log_counts(n_k)?;
    let mut reduced_b = reduced_b_from_free_energies(&ln_n_k, f_k);
    let mut full_b = vec![0.0; n_states];
    fill_full_b(&reduced_b, &mut full_b);
    let (mut objective, full_grad) =
        mbar_objective_gradient(u_nk, n_states, n_k, &full_b, parallel);
    let mut gradient = full_grad[1..].to_vec();
    let mut history: Vec<LbfgsHistoryEntry> = Vec::new();

    for _ in 0..max_iterations {
        if max_abs(&gradient) < tolerance {
            write_free_energies_from_reduced_b(&reduced_b, &ln_n_k, f_k);
            return Ok(());
        }

        let mut direction = lbfgs_direction(&gradient, &history);
        let mut directional_derivative = dot(&gradient, &direction);
        if !directional_derivative.is_finite() || directional_derivative >= 0.0 {
            direction = gradient.iter().map(|value| -*value).collect();
            directional_derivative = -dot(&gradient, &gradient);
        }

        let mut step = 1.0;
        let mut accepted = None;
        for _ in 0..LBFGS_MAX_LINE_SEARCH_STEPS {
            let trial_reduced = add_scaled(&reduced_b, &direction, step);
            let mut trial_full = vec![0.0; n_states];
            fill_full_b(&trial_reduced, &mut trial_full);
            let (trial_objective, trial_full_grad) =
                mbar_objective_gradient(u_nk, n_states, n_k, &trial_full, parallel);
            if trial_objective <= objective + LBFGS_C1 * step * directional_derivative {
                accepted = Some((trial_reduced, trial_objective, trial_full_grad));
                break;
            }
            step *= 0.5;
            if step < LBFGS_MIN_STEP {
                break;
            }
        }

        let Some((trial_reduced, trial_objective, trial_full_grad)) = accepted else {
            break;
        };
        let trial_gradient = trial_full_grad[1..].to_vec();
        let s = subtract(&trial_reduced, &reduced_b);
        let y = subtract(&trial_gradient, &gradient);
        let ys = dot(&y, &s);

        reduced_b = trial_reduced;
        objective = trial_objective;
        gradient = trial_gradient;

        if max_abs(&s) < tolerance && max_abs(&gradient) < tolerance {
            write_free_energies_from_reduced_b(&reduced_b, &ln_n_k, f_k);
            return Ok(());
        }

        if ys > LBFGS_CURVATURE_EPS {
            history.push(LbfgsHistoryEntry {
                s,
                y,
                rho: 1.0 / ys,
            });
            if history.len() > LBFGS_HISTORY_LIMIT {
                history.remove(0);
            }
        }
    }

    if max_abs(&gradient) < tolerance {
        write_free_energies_from_reduced_b(&reduced_b, &ln_n_k, f_k);
        Ok(())
    } else {
        Err(CoreError::ConvergenceFailure)
    }
}

fn mbar_uncertainty(
    u_kn: &[Vec<f64>],
    u_nk: &[f64],
    n_k: &[f64],
    f_k: &[f64],
    parallel: bool,
) -> Result<Vec<f64>> {
    let n_states = u_kn.len();
    let theta = mbar_theta(u_kn, u_nk, n_k, f_k, parallel)?;

    let mut uncertainties = vec![0.0; n_states * n_states];
    for i in 0..n_states {
        for j in 0..n_states {
            uncertainties[i * n_states + j] = pair_uncertainty_from_theta(&theta, i, j);
        }
    }
    Ok(uncertainties)
}

fn mbar_theta(
    u_kn: &[Vec<f64>],
    u_nk: &[f64],
    n_k: &[f64],
    f_k: &[f64],
    parallel: bool,
) -> Result<nalgebra::DMatrix<f64>> {
    use nalgebra::{DMatrix, SymmetricEigen};

    let n_states = u_kn.len();
    let n_samples = u_kn[0].len();

    let log_w = mbar_log_weights(u_kn, u_nk, n_k, f_k, parallel);
    let mut w_data = vec![0.0; n_samples * n_states];
    if parallel {
        w_data
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, slot)| *slot = log_w[idx].exp());
    } else {
        for (idx, slot) in w_data.iter_mut().enumerate() {
            *slot = log_w[idx].exp();
        }
    }
    let w = DMatrix::from_row_slice(n_samples, n_states, &w_data);

    let wtw = &w.transpose() * &w;
    let eigen = SymmetricEigen::new(wtw);
    let mut s2 = eigen.eigenvalues;
    for value in s2.iter_mut() {
        if *value < 0.0 {
            *value = 0.0;
        }
    }
    let sigma = DMatrix::from_diagonal(&s2.map(|v| v.sqrt()));
    let v = eigen.eigenvectors;

    let mut ndiag = DMatrix::zeros(n_states, n_states);
    for i in 0..n_states {
        ndiag[(i, i)] = n_k[i];
    }

    let identity = DMatrix::identity(n_states, n_states);
    let a = &identity - &sigma * &v.transpose() * &ndiag * &v * &sigma;
    let a_inv = pseudoinverse(&a, 1.0e-10)?;
    Ok(&v * &sigma * a_inv * &sigma * v.transpose())
}

fn positive_log_counts(n_k: &[f64]) -> Result<Vec<f64>> {
    let mut values = Vec::with_capacity(n_k.len());
    for (state_idx, &count) in n_k.iter().enumerate() {
        if count <= 0.0 || !count.is_finite() {
            return Err(CoreError::InvalidState(format!(
                "MBAR solver requires a positive finite sample count for state {state_idx}"
            )));
        }
        values.push(count.ln());
    }
    Ok(values)
}

fn reduced_b_from_free_energies(ln_n_k: &[f64], f_k: &[f64]) -> Vec<f64> {
    let mut full_b = Vec::with_capacity(f_k.len());
    for (&ln_n, &free_energy) in ln_n_k.iter().zip(f_k.iter()) {
        full_b.push(-(ln_n + free_energy));
    }
    let shift = full_b[0];
    for value in &mut full_b {
        *value -= shift;
    }
    full_b[1..].to_vec()
}

fn write_free_energies_from_reduced_b(reduced_b: &[f64], ln_n_k: &[f64], out: &mut [f64]) {
    out[0] = -ln_n_k[0];
    for (state_idx, slot) in out.iter_mut().enumerate().skip(1) {
        *slot = -reduced_b[state_idx - 1] - ln_n_k[state_idx];
    }
    let shift = out[0];
    for value in out.iter_mut() {
        *value -= shift;
    }
}

fn fill_full_b(reduced_b: &[f64], full_b: &mut [f64]) {
    full_b[0] = 0.0;
    full_b[1..].copy_from_slice(reduced_b);
}

fn mbar_objective_gradient(
    u_nk: &[f64],
    n_states: usize,
    n_k: &[f64],
    full_b: &[f64],
    parallel: bool,
) -> (f64, Vec<f64>) {
    let n_samples = u_nk.len() / n_states;
    let (sum_log_denom, mut gradient) = if parallel {
        (0..n_samples)
            .into_par_iter()
            .fold(
                || (0.0, vec![0.0; n_states]),
                |(objective_sum, mut gradient), sample_idx| {
                    let row = sample_major_row(u_nk, n_states, sample_idx);
                    let log_denom = row_logsumexp_with_bias(row, full_b);
                    for (state_idx, value) in row.iter().enumerate() {
                        gradient[state_idx] -= (-*value - full_b[state_idx] - log_denom).exp();
                    }
                    (objective_sum + log_denom, gradient)
                },
            )
            .reduce(
                || (0.0, vec![0.0; n_states]),
                |(left_obj, mut left_grad), (right_obj, right_grad)| {
                    for (lhs, rhs) in left_grad.iter_mut().zip(right_grad.iter()) {
                        *lhs += *rhs;
                    }
                    (left_obj + right_obj, left_grad)
                },
            )
    } else {
        let mut objective_sum = 0.0;
        let mut gradient = vec![0.0; n_states];
        for sample_idx in 0..n_samples {
            let row = sample_major_row(u_nk, n_states, sample_idx);
            let log_denom = row_logsumexp_with_bias(row, full_b);
            objective_sum += log_denom;
            for (state_idx, value) in row.iter().enumerate() {
                gradient[state_idx] -= (-*value - full_b[state_idx] - log_denom).exp();
            }
        }
        (objective_sum, gradient)
    };

    let inv_n = 1.0 / n_samples as f64;
    for (state_idx, value) in gradient.iter_mut().enumerate() {
        *value = (*value + n_k[state_idx]) * inv_n;
    }
    let objective = inv_n * (sum_log_denom + dot(full_b, n_k));
    (objective, gradient)
}

fn row_logsumexp_with_bias(row: &[f64], full_b: &[f64]) -> f64 {
    let mut max_arg = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for (&value, &bias) in row.iter().zip(full_b.iter()) {
        accumulate_logsumexp(&mut max_arg, &mut sum, -value - bias);
    }
    max_arg + sum.ln()
}

fn lbfgs_direction(gradient: &[f64], history: &[LbfgsHistoryEntry]) -> Vec<f64> {
    if history.is_empty() {
        return gradient.iter().map(|value| -*value).collect();
    }

    let mut q = gradient.to_vec();
    let mut alpha = vec![0.0; history.len()];
    for (idx, entry) in history.iter().enumerate().rev() {
        alpha[idx] = entry.rho * dot(&entry.s, &q);
        axpy(&mut q, -alpha[idx], &entry.y);
    }

    let gamma = if let Some(last) = history.last() {
        let yy = dot(&last.y, &last.y);
        if yy > 0.0 {
            dot(&last.s, &last.y) / yy
        } else {
            1.0
        }
    } else {
        1.0
    };
    scale(&mut q, gamma);

    for (idx, entry) in history.iter().enumerate() {
        let beta = entry.rho * dot(&entry.y, &q);
        axpy(&mut q, alpha[idx] - beta, &entry.s);
    }

    scale(&mut q, -1.0);
    q
}

fn dot(lhs: &[f64], rhs: &[f64]) -> f64 {
    lhs.iter().zip(rhs.iter()).map(|(l, r)| l * r).sum()
}

fn axpy(dst: &mut [f64], alpha: f64, src: &[f64]) {
    for (lhs, rhs) in dst.iter_mut().zip(src.iter()) {
        *lhs += alpha * rhs;
    }
}

fn scale(values: &mut [f64], factor: f64) {
    for value in values.iter_mut() {
        *value *= factor;
    }
}

fn add_scaled(lhs: &[f64], rhs: &[f64], scale: f64) -> Vec<f64> {
    lhs.iter()
        .zip(rhs.iter())
        .map(|(l, r)| l + scale * r)
        .collect()
}

fn subtract(lhs: &[f64], rhs: &[f64]) -> Vec<f64> {
    lhs.iter().zip(rhs.iter()).map(|(l, r)| l - r).collect()
}

fn max_abs(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0, |max_value, value| max_value.max(value.abs()))
}

fn pair_uncertainty_from_theta(
    theta: &nalgebra::DMatrix<f64>,
    from_idx: usize,
    to_idx: usize,
) -> f64 {
    let val =
        theta[(from_idx, from_idx)] + theta[(to_idx, to_idx)] - 2.0 * theta[(from_idx, to_idx)];
    if val < 0.0 && val > -1e-10 {
        0.0
    } else if val < 0.0 {
        val.abs().sqrt()
    } else {
        val.sqrt()
    }
}

fn mbar_log_weights(
    u_kn: &[Vec<f64>],
    u_nk: &[f64],
    n_k: &[f64],
    f_k: &[f64],
    parallel: bool,
) -> Vec<f64> {
    let n_states = u_kn.len();
    let ln_n_k: Vec<f64> = n_k
        .iter()
        .map(|n| if *n > 0.0 { n.ln() } else { f64::NEG_INFINITY })
        .collect();
    let state_offsets = state_offsets(&ln_n_k, f_k);

    let mut log_w = vec![0.0; u_kn[0].len() * n_states];
    if parallel {
        log_w
            .par_chunks_mut(n_states)
            .enumerate()
            .for_each(|(n, row)| fill_log_weights_row(u_kn, &state_offsets, f_k, n, row));
    } else {
        for (n, row) in log_w.chunks_mut(n_states).enumerate() {
            fill_log_weights_row_serial(u_nk, n_states, &state_offsets, f_k, n, row);
        }
    }
    log_w
}

fn fill_log_denominator(
    u_kn: &[Vec<f64>],
    state_offsets: &[f64],
    log_denominator: &mut [f64],
    parallel: bool,
) -> Result<()> {
    if parallel {
        log_denominator
            .par_iter_mut()
            .enumerate()
            .try_for_each(|(n, slot)| -> Result<()> {
                *slot = checked_log_denominator_for_sample(u_kn, state_offsets, n)?;
                Ok(())
            })?;
    } else {
        for (n, slot) in log_denominator.iter_mut().enumerate() {
            *slot = checked_log_denominator_for_sample(u_kn, state_offsets, n)?;
        }
    }
    Ok(())
}

fn fill_log_denominator_serial(
    u_nk: &[f64],
    n_states: usize,
    state_offsets: &[f64],
    log_denominator: &mut [f64],
) -> Result<()> {
    for (n, slot) in log_denominator.iter_mut().enumerate() {
        *slot = checked_log_denominator_for_sample_serial(u_nk, n_states, state_offsets, n)?;
    }
    Ok(())
}

fn fill_free_energies(
    u_kn: &[Vec<f64>],
    log_denominator: &[f64],
    f_new: &mut [f64],
    parallel: bool,
) -> Result<()> {
    if parallel {
        f_new
            .par_iter_mut()
            .enumerate()
            .try_for_each(|(k, slot)| -> Result<()> {
                *slot = checked_free_energy_for_state(u_kn, log_denominator, k)?;
                Ok(())
            })?;
    } else {
        for (k, slot) in f_new.iter_mut().enumerate() {
            *slot = checked_free_energy_for_state(u_kn, log_denominator, k)?;
        }
    }
    Ok(())
}

fn fill_free_energies_serial(
    u_nk: &[f64],
    n_states: usize,
    log_denominator: &[f64],
    f_new: &mut [f64],
) -> Result<()> {
    let mut max_args = vec![f64::NEG_INFINITY; n_states];
    let mut sums = vec![0.0; n_states];
    for (n, denom) in log_denominator.iter().enumerate() {
        let row = sample_major_row(u_nk, n_states, n);
        for (k, value) in row.iter().enumerate() {
            let val = -*denom - *value;
            accumulate_logsumexp(&mut max_args[k], &mut sums[k], val);
        }
    }
    for (state_idx, slot) in f_new.iter_mut().enumerate() {
        *slot = if max_args[state_idx].is_finite() {
            -(max_args[state_idx] + sums[state_idx].ln())
        } else {
            f64::NAN
        };
        if !slot.is_finite() {
            return Err(CoreError::NonFiniteValue(format!(
                "MBAR free energy became non-finite for state {state_idx}"
            )));
        }
    }
    Ok(())
}

fn checked_log_denominator_for_sample(
    u_kn: &[Vec<f64>],
    state_offsets: &[f64],
    sample_idx: usize,
) -> Result<f64> {
    let value = log_denominator_for_sample(u_kn, state_offsets, sample_idx);
    if !value.is_finite() {
        return Err(CoreError::NonFiniteValue(format!(
            "MBAR denominator became non-finite at sample {sample_idx}"
        )));
    }
    Ok(value)
}

fn checked_log_denominator_for_sample_serial(
    u_nk: &[f64],
    n_states: usize,
    state_offsets: &[f64],
    sample_idx: usize,
) -> Result<f64> {
    let value = log_denominator_for_sample_serial(u_nk, n_states, state_offsets, sample_idx);
    if !value.is_finite() {
        return Err(CoreError::NonFiniteValue(format!(
            "MBAR denominator became non-finite at sample {sample_idx}"
        )));
    }
    Ok(value)
}

fn log_denominator_for_sample(u_kn: &[Vec<f64>], state_offsets: &[f64], sample_idx: usize) -> f64 {
    let mut max_arg = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for (k, state_offset) in state_offsets.iter().enumerate() {
        if !state_offset.is_finite() {
            continue;
        }
        let val = state_offset - u_kn[k][sample_idx];
        accumulate_logsumexp(&mut max_arg, &mut sum, val);
    }
    max_arg + sum.ln()
}

fn log_denominator_for_sample_serial(
    u_nk: &[f64],
    n_states: usize,
    state_offsets: &[f64],
    sample_idx: usize,
) -> f64 {
    let mut max_arg = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for (state_offset, value) in state_offsets
        .iter()
        .zip(sample_major_row(u_nk, n_states, sample_idx).iter())
    {
        if !state_offset.is_finite() {
            continue;
        }
        accumulate_logsumexp(&mut max_arg, &mut sum, state_offset - *value);
    }
    max_arg + sum.ln()
}

fn checked_free_energy_for_state(
    u_kn: &[Vec<f64>],
    log_denominator: &[f64],
    state_idx: usize,
) -> Result<f64> {
    let value = free_energy_for_state(u_kn, log_denominator, state_idx);
    if !value.is_finite() {
        return Err(CoreError::NonFiniteValue(format!(
            "MBAR free energy became non-finite for state {state_idx}"
        )));
    }
    Ok(value)
}

fn free_energy_for_state(u_kn: &[Vec<f64>], log_denominator: &[f64], state_idx: usize) -> f64 {
    let mut max_arg = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for (n, denom) in log_denominator.iter().enumerate() {
        let val = -*denom - u_kn[state_idx][n];
        accumulate_logsumexp(&mut max_arg, &mut sum, val);
    }
    -(max_arg + sum.ln())
}

fn fill_log_weights_row(
    u_kn: &[Vec<f64>],
    state_offsets: &[f64],
    f_k: &[f64],
    sample_idx: usize,
    row: &mut [f64],
) {
    let log_denom = log_denominator_for_sample(u_kn, state_offsets, sample_idx);
    for k in 0..u_kn.len() {
        row[k] = f_k[k] - u_kn[k][sample_idx] - log_denom;
    }
}

fn fill_log_weights_row_serial(
    u_nk: &[f64],
    n_states: usize,
    state_offsets: &[f64],
    f_k: &[f64],
    sample_idx: usize,
    row: &mut [f64],
) {
    let log_denom = log_denominator_for_sample_serial(u_nk, n_states, state_offsets, sample_idx);
    for (k, value) in sample_major_row(u_nk, n_states, sample_idx)
        .iter()
        .enumerate()
    {
        row[k] = f_k[k] - *value - log_denom;
    }
}

fn sample_major_row(u_nk: &[f64], n_states: usize, sample_idx: usize) -> &[f64] {
    let start = sample_idx * n_states;
    &u_nk[start..start + n_states]
}

fn sample_major_from_state_major(u_kn: &[Vec<f64>]) -> Vec<f64> {
    if u_kn.is_empty() {
        return Vec::new();
    }
    let n_states = u_kn.len();
    let n_samples = u_kn[0].len();
    let mut u_nk = vec![0.0; n_samples * n_states];
    for sample_idx in 0..n_samples {
        let row = &mut u_nk[sample_idx * n_states..(sample_idx + 1) * n_states];
        for state_idx in 0..n_states {
            row[state_idx] = u_kn[state_idx][sample_idx];
        }
    }
    u_nk
}

fn state_offsets(ln_n_k: &[f64], f_k: &[f64]) -> Vec<f64> {
    let mut offsets = vec![0.0; ln_n_k.len()];
    fill_state_offsets(ln_n_k, f_k, &mut offsets);
    offsets
}

fn fill_state_offsets(ln_n_k: &[f64], f_k: &[f64], offsets: &mut [f64]) {
    for ((slot, ln_n), f) in offsets.iter_mut().zip(ln_n_k.iter()).zip(f_k.iter()) {
        *slot = if ln_n.is_finite() {
            ln_n + f
        } else {
            f64::NEG_INFINITY
        };
    }
}

fn accumulate_logsumexp(max_arg: &mut f64, sum: &mut f64, value: f64) {
    // Stable one-pass log-sum-exp update. This avoids a separate max scan and sum scan in
    // the MBAR hot loops while preserving the usual numerical stability.
    if !max_arg.is_finite() {
        *max_arg = value;
        *sum = 1.0;
    } else if value <= *max_arg {
        *sum += (value - *max_arg).exp();
    } else {
        *sum = (*sum) * (*max_arg - value).exp() + 1.0;
        *max_arg = value;
    }
}

fn pseudoinverse(matrix: &nalgebra::DMatrix<f64>, tol: f64) -> Result<nalgebra::DMatrix<f64>> {
    let svd = matrix.clone().svd(true, true);
    let u = svd.u.ok_or(CoreError::ConvergenceFailure)?;
    let v_t = svd.v_t.ok_or(CoreError::ConvergenceFailure)?;
    let mut sigma_inv = nalgebra::DMatrix::zeros(matrix.nrows(), matrix.ncols());
    for i in 0..svd.singular_values.len() {
        let value = svd.singular_values[i];
        if value > tol {
            sigma_inv[(i, i)] = 1.0 / value;
        }
    }
    Ok(v_t.transpose() * sigma_inv * u.transpose())
}

fn validate_mbar_input(u_kn: &[Vec<f64>]) -> Result<()> {
    if u_kn.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let n_samples = u_kn[0].len();
    for (state_idx, row) in u_kn.iter().enumerate() {
        if row.len() != n_samples {
            return Err(CoreError::InvalidShape {
                expected: n_samples,
                found: row.len(),
            });
        }
        let has_finite = row.iter().any(|value| value.is_finite());
        if !has_finite {
            return Err(CoreError::NonFiniteValue(format!(
                "MBAR state {state_idx} has no finite reduced potentials"
            )));
        }
        for (sample_idx, value) in row.iter().enumerate() {
            if value.is_nan() || *value == f64::NEG_INFINITY {
                return Err(CoreError::NonFiniteValue(format!(
                    "MBAR reduced potential for state {state_idx} sample {sample_idx} must be finite or +inf"
                )));
            }
        }
    }
    for sample_idx in 0..n_samples {
        if !u_kn.iter().any(|row| row[sample_idx].is_finite()) {
            return Err(CoreError::NonFiniteValue(format!(
                "MBAR sample {sample_idx} has no finite reduced potentials"
            )));
        }
    }
    Ok(())
}
