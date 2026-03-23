use crate::data::{OverlapMatrix, UNkMatrix};
use crate::error::{CoreError, Result};
use crate::estimators::{mbar_log_weights_from_windows, MbarOptions};
use nalgebra::{DMatrix, Schur};

pub fn overlap_matrix(
    windows: &[UNkMatrix],
    options: Option<MbarOptions>,
) -> Result<OverlapMatrix> {
    let options = options.unwrap_or_default();
    let (log_w, n_k, states) = mbar_log_weights_from_windows(windows, &options)?;
    let n_states = states.len();
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
            for j in 0..n_states {
                wtw[i * n_states + j] += wi * weights[j];
            }
        }
    }

    let mut values = vec![0.0; n_states * n_states];
    for i in 0..n_states {
        for j in 0..n_states {
            values[i * n_states + j] = wtw[i * n_states + j] * n_k[j];
        }
    }

    OverlapMatrix::new(values, n_states, states)
}

pub fn overlap_eigenvalues(overlap: &OverlapMatrix) -> Result<Vec<f64>> {
    let n_states = overlap.n_states();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let matrix = DMatrix::from_row_slice(n_states, n_states, overlap.values());
    let eigenvalues = Schur::new(matrix).complex_eigenvalues();
    let mut values = Vec::with_capacity(n_states);
    for eig in eigenvalues.iter() {
        if eig.im.abs() > 1e-8 {
            return Err(CoreError::InvalidState(
                "overlap eigenvalue has significant imaginary component".to_string(),
            ));
        }
        values.push(eig.re);
    }
    values.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    Ok(values)
}

pub fn overlap_scalar(overlap: &OverlapMatrix) -> Result<f64> {
    if overlap.n_states() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: overlap.n_states(),
        });
    }
    let eigenvalues = overlap_eigenvalues(overlap)?;
    Ok(1.0 - eigenvalues[1])
}
