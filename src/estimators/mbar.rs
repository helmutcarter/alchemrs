use crate::analysis::{self, BlockEstimate};
use crate::data::{find_state_index_exact, DeltaFMatrix, OverlapMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};
use rayon::prelude::*;
use std::sync::OnceLock;

use super::common::{ensure_consistent_lambda_labels, ensure_consistent_states, CombinedWindows};

#[derive(Debug, Clone)]
pub struct MbarOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub initial_f_k: Option<Vec<f64>>,
    pub parallel: bool,
}

impl Default for MbarOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
            tolerance: 1.0e-7,
            initial_f_k: None,
            parallel: false,
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

        mbar_solve(
            &prepared.u_kn,
            &prepared.n_k,
            &mut f_k,
            self.options.tolerance,
            self.options.max_iterations,
            self.options.parallel,
        )?;

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
                for j in i..n_states {
                    let idx = row_offset + j;
                    wtw[idx] = wi.mul_add(weights[j], wtw[idx]);
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
        Ok(Self {
            u_kn,
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

fn mbar_solve(
    u_kn: &[Vec<f64>],
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
    for _ in 0..max_iterations {
        // Precompute ln(n_k) + f_k once per iteration so the sample loop only pays the
        // u_kn access and the log-sum-exp update.
        let state_offsets = state_offsets(&ln_n_k, f_k);
        fill_log_denominator(u_kn, &state_offsets, &mut log_denominator, parallel)?;

        let mut f_new = vec![0.0; n_states];
        fill_free_energies(u_kn, &log_denominator, &mut f_new, parallel)?;
        let shift = f_new[0];
        for value in f_new.iter_mut().take(n_states) {
            *value -= shift;
        }

        let mut max_delta = 0.0;
        for k in 0..n_states {
            let delta = (f_new[k] - f_k[k]).abs();
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

fn mbar_uncertainty(
    u_kn: &[Vec<f64>],
    n_k: &[f64],
    f_k: &[f64],
    parallel: bool,
) -> Result<Vec<f64>> {
    use nalgebra::{DMatrix, SymmetricEigen};

    let n_states = u_kn.len();
    let n_samples = u_kn[0].len();

    let log_w = mbar_log_weights(u_kn, n_k, f_k, parallel);
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
    let theta = &v * &sigma * a_inv * &sigma * v.transpose();

    let mut uncertainties = vec![0.0; n_states * n_states];
    for i in 0..n_states {
        for j in 0..n_states {
            let val = theta[(i, i)] + theta[(j, j)] - 2.0 * theta[(i, j)];
            uncertainties[i * n_states + j] = if val < 0.0 && val > -1e-10 {
                0.0
            } else if val < 0.0 {
                val.abs().sqrt()
            } else {
                val.sqrt()
            };
        }
    }
    Ok(uncertainties)
}

fn mbar_log_weights(u_kn: &[Vec<f64>], n_k: &[f64], f_k: &[f64], parallel: bool) -> Vec<f64> {
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
            fill_log_weights_row(u_kn, &state_offsets, f_k, n, row);
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

fn state_offsets(ln_n_k: &[f64], f_k: &[f64]) -> Vec<f64> {
    ln_n_k
        .iter()
        .zip(f_k.iter())
        .map(|(ln_n, f)| {
            if ln_n.is_finite() {
                ln_n + f
            } else {
                f64::NEG_INFINITY
            }
        })
        .collect()
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
