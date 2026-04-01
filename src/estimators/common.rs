use crate::data::{StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};

pub(crate) type CombinedWindows = (
    Vec<Vec<f64>>,
    Vec<f64>,
    Vec<StatePoint>,
    Option<Vec<String>>,
);
pub(crate) type PairEstimate = (f64, f64);
pub(crate) type ExpRow = (usize, Vec<f64>, Vec<f64>);

pub(crate) fn ensure_consistent_states(windows: &[UNkMatrix]) -> Result<Vec<StatePoint>> {
    let first = windows[0].evaluated_states();
    let first_temperature = first[0].temperature_k();
    for state in first {
        if (state.temperature_k() - first_temperature).abs() > 1e-6 {
            return Err(CoreError::InvalidState(
                "evaluated_states temperatures differ within a window".to_string(),
            ));
        }
    }
    for window in windows.iter().skip(1) {
        let states = window.evaluated_states();
        if states.len() != first.len() {
            return Err(CoreError::InvalidShape {
                expected: first.len(),
                found: states.len(),
            });
        }
        for (a, b) in first.iter().zip(states.iter()) {
            if a.lambdas().len() != b.lambdas().len() {
                return Err(CoreError::InvalidShape {
                    expected: a.lambdas().len(),
                    found: b.lambdas().len(),
                });
            }
            if (a.temperature_k() - b.temperature_k()).abs() > 1e-6 {
                return Err(CoreError::InvalidState(
                    "evaluated_states temperatures differ between windows".to_string(),
                ));
            }
            if !states_match(a, b) {
                return Err(CoreError::InvalidState(
                    "evaluated_states differ between windows".to_string(),
                ));
            }
        }
    }
    for window in windows {
        let sampled = window.sampled_state().ok_or_else(|| {
            CoreError::InvalidState("sampled_state required for estimator".to_string())
        })?;
        if sampled.lambdas().len() != first[0].lambdas().len() {
            return Err(CoreError::InvalidShape {
                expected: first[0].lambdas().len(),
                found: sampled.lambdas().len(),
            });
        }
        if (sampled.temperature_k() - first_temperature).abs() > 1e-6 {
            return Err(CoreError::InvalidState(
                "sampled_state temperature differs from evaluated_states".to_string(),
            ));
        }
        if find_state_index(first, sampled).is_err() {
            return Err(CoreError::InvalidState(
                "sampled_state not found in evaluated_states".to_string(),
            ));
        }
    }
    Ok(first.to_vec())
}

pub(crate) fn ensure_consistent_lambda_labels(
    windows: &[UNkMatrix],
) -> Result<Option<Vec<String>>> {
    let first = windows[0].lambda_labels().map(|labels| labels.to_vec());
    for window in windows.iter().skip(1) {
        match (&first, window.lambda_labels()) {
            (None, None) => {}
            (Some(expected), Some(found)) if expected == found => {}
            _ => {
                return Err(CoreError::InvalidState(
                    "lambda_labels differ between windows".to_string(),
                ))
            }
        }
    }
    Ok(first)
}

pub(crate) fn states_match(a: &StatePoint, b: &StatePoint) -> bool {
    a.lambdas().len() == b.lambdas().len()
        && a.lambdas()
            .iter()
            .zip(b.lambdas().iter())
            .all(|(left, right)| (*left - *right).abs() < 1e-6)
}

pub(crate) fn find_state_index(states: &[StatePoint], target: &StatePoint) -> Result<usize> {
    for (idx, state) in states.iter().enumerate() {
        if states_match(state, target)
            && (state.temperature_k() - target.temperature_k()).abs() < 1e-6
        {
            return Ok(idx);
        }
    }
    Err(CoreError::InvalidState(
        "sampled_state not found in evaluated_states".to_string(),
    ))
}

pub(crate) fn work_values(window: &UNkMatrix, idx_a: usize, idx_b: usize) -> Result<Vec<f64>> {
    let n_states = window.n_states();
    let n_samples = window.n_samples();
    if idx_b >= n_states {
        return Err(CoreError::InvalidShape {
            expected: n_states,
            found: idx_b + 1,
        });
    }
    let data = window.data();
    let mut out = Vec::with_capacity(n_samples);
    for sample_idx in 0..n_samples {
        let offset = sample_idx * n_states;
        let a = data[offset + idx_a];
        let b = data[offset + idx_b];
        if a.is_nan() || b.is_nan() {
            return Err(CoreError::NonFiniteValue(format!(
                "u_nk work input for states ({idx_a}, {idx_b}) at sample {sample_idx} must not be NaN"
            )));
        }
        let delta = b - a;
        if delta.is_nan() {
            return Err(CoreError::NonFiniteValue(format!(
                "work value for states ({idx_a}, {idx_b}) at sample {sample_idx} must not be NaN"
            )));
        }
        out.push(delta);
    }
    Ok(out)
}
