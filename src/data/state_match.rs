use crate::data::StatePoint;
use crate::error::{CoreError, Result};

pub(crate) fn state_points_match_exact(left: &StatePoint, right: &StatePoint) -> bool {
    left.temperature_k() == right.temperature_k()
        && left.lambdas().len() == right.lambdas().len()
        && left
            .lambdas()
            .iter()
            .zip(right.lambdas().iter())
            .all(|(lhs, rhs)| lhs == rhs)
}

pub(crate) fn find_state_index_exact(states: &[StatePoint], target: &StatePoint) -> Result<usize> {
    for (idx, state) in states.iter().enumerate() {
        if state_points_match_exact(state, target) {
            return Ok(idx);
        }
    }
    Err(CoreError::InvalidState(
        "sampled_state not found in evaluated_states".to_string(),
    ))
}

pub(crate) fn find_scalar_lambda_state_index_exact(
    states: &[StatePoint],
    lambda: f64,
) -> Result<usize> {
    for (idx, state) in states.iter().enumerate() {
        if state.lambdas().len() == 1 && state.lambdas()[0] == lambda {
            return Ok(idx);
        }
    }
    Err(CoreError::InvalidState(
        "sampled_state not found in evaluated_states".to_string(),
    ))
}
