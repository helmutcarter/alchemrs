use crate::analysis::{self, BlockEstimate};
use crate::data::{DeltaFMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};

use super::common::{ensure_consistent_lambda_labels, ensure_consistent_states, find_state_index, CombinedWindows};

#[derive(Debug, Clone)]
pub struct MbarOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub initial_f_k: Option<Vec<f64>>,
    pub compute_uncertainty: bool,
    pub parallel: bool,
}

impl Default for MbarOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
            tolerance: 1.0e-7,
            initial_f_k: None,
            compute_uncertainty: true,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MbarEstimator {
    pub options: MbarOptions,
}

impl MbarEstimator {
    pub fn new(options: MbarOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        if windows.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        let (u_kn, n_k, states, lambda_labels) = combine_windows(windows)?;
        let mut f_k = initial_f_k(&self.options, states.len())?;

        mbar_solve(
            &u_kn,
            &n_k,
            &mut f_k,
            self.options.tolerance,
            self.options.max_iterations,
        )?;

        let n_states = states.len();
        let mut values = vec![0.0; n_states * n_states];
        for i in 0..n_states {
            for j in 0..n_states {
                values[i * n_states + j] = f_k[j] - f_k[i];
            }
        }
        let uncertainties = if self.options.compute_uncertainty {
            Some(mbar_uncertainty(&u_kn, &n_k, &f_k)?)
        } else {
            None
        };
        DeltaFMatrix::new_with_labels(values, uncertainties, n_states, states, lambda_labels)
    }

    pub fn block_average(
        &self,
        windows: &[UNkMatrix],
        n_blocks: usize,
    ) -> Result<Vec<BlockEstimate>> {
        analysis::mbar_block_average(windows, n_blocks, Some(self.options.clone()))
    }
}

pub fn mbar_log_weights_from_windows(
    windows: &[UNkMatrix],
    options: &MbarOptions,
) -> Result<(Vec<f64>, Vec<f64>, Vec<StatePoint>)> {
    if windows.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let (u_kn, n_k, states, _lambda_labels) = combine_windows(windows)?;
    let mut f_k = initial_f_k(options, states.len())?;
    mbar_solve(
        &u_kn,
        &n_k,
        &mut f_k,
        options.tolerance,
        options.max_iterations,
    )?;
    let log_weights = mbar_log_weights(&u_kn, &n_k, &f_k)?;
    Ok((log_weights, n_k, states))
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
        let idx = find_state_index(&states, sampled)?;
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
        for n in 0..n_samples {
            let mut max_arg = f64::NEG_INFINITY;
            for k in 0..n_states {
                if n_k[k] == 0.0 {
                    continue;
                }
                let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
                if val > max_arg {
                    max_arg = val;
                }
            }
            let mut sum = 0.0;
            for k in 0..n_states {
                if n_k[k] == 0.0 {
                    continue;
                }
                let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
                sum += (val - max_arg).exp();
            }
            log_denominator[n] = max_arg + sum.ln();
            if !log_denominator[n].is_finite() {
                return Err(CoreError::NonFiniteValue(format!(
                    "MBAR denominator became non-finite at sample {n}"
                )));
            }
        }

        let mut f_new = vec![0.0; n_states];
        for k in 0..n_states {
            let mut max_arg = f64::NEG_INFINITY;
            for n in 0..n_samples {
                let val = -log_denominator[n] - u_kn[k][n];
                if val > max_arg {
                    max_arg = val;
                }
            }
            let mut sum = 0.0;
            for n in 0..n_samples {
                let val = -log_denominator[n] - u_kn[k][n];
                sum += (val - max_arg).exp();
            }
            f_new[k] = -(max_arg + sum.ln());
            if !f_new[k].is_finite() {
                return Err(CoreError::NonFiniteValue(format!(
                    "MBAR free energy became non-finite for state {k}"
                )));
            }
        }
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

fn mbar_uncertainty(u_kn: &[Vec<f64>], n_k: &[f64], f_k: &[f64]) -> Result<Vec<f64>> {
    use nalgebra::{DMatrix, SymmetricEigen};

    let n_states = u_kn.len();
    let n_samples = u_kn[0].len();

    let log_w = mbar_log_weights(u_kn, n_k, f_k)?;
    let mut w_data = Vec::with_capacity(n_samples * n_states);
    for n in 0..n_samples {
        for k in 0..n_states {
            w_data.push(log_w[n * n_states + k].exp());
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

fn mbar_log_weights(u_kn: &[Vec<f64>], n_k: &[f64], f_k: &[f64]) -> Result<Vec<f64>> {
    let n_states = u_kn.len();
    let n_samples = u_kn[0].len();
    let ln_n_k: Vec<f64> = n_k
        .iter()
        .map(|n| if *n > 0.0 { n.ln() } else { f64::NEG_INFINITY })
        .collect();

    let mut log_w = vec![0.0; n_samples * n_states];
    for n in 0..n_samples {
        let mut max_arg = f64::NEG_INFINITY;
        for k in 0..n_states {
            if n_k[k] == 0.0 {
                continue;
            }
            let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
            if val > max_arg {
                max_arg = val;
            }
        }
        let mut sum = 0.0;
        for k in 0..n_states {
            if n_k[k] == 0.0 {
                continue;
            }
            let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
            sum += (val - max_arg).exp();
        }
        let log_denom = max_arg + sum.ln();
        for k in 0..n_states {
            log_w[n * n_states + k] = f_k[k] - u_kn[k][n] - log_denom;
        }
    }
    Ok(log_w)
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
