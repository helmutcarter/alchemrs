use crate::data::{find_state_index_exact, AtmLogQMatrix, DeltaFMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};
use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use std::fs;
use std::path::Path;
use std::sync::OnceLock;

use super::common::{ensure_consistent_lambda_labels, ensure_consistent_states};

#[derive(Debug, Clone)]
pub struct UwhamOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub parallel: bool,
}

impl Default for UwhamOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
            tolerance: 1.0e-7,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct UwhamEstimator {
    pub options: UwhamOptions,
}

#[derive(Debug, Clone)]
pub(crate) struct UwhamLogQInput {
    log_q: Vec<f64>,
    n_observations: usize,
    n_states: usize,
    size: Vec<f64>,
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
}

#[derive(Debug, Clone)]
struct ReducedUwhamInput {
    log_q_relative: Vec<f64>,
    sampled_indices: Vec<usize>,
    base_sampled_index: usize,
    rho: Vec<f64>,
    sampled_log_q_relative: Vec<f64>,
}

#[derive(Debug)]
pub struct UwhamFit {
    input: UwhamLogQInput,
    free_energies: Vec<f64>,
    weights: OnceLock<Vec<f64>>,
    parallel: bool,
}

#[derive(Debug)]
struct UwhamObjective {
    value: f64,
    gradient: DVector<f64>,
    hessian: DMatrix<f64>,
}

impl UwhamEstimator {
    pub fn new(options: UwhamOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<UwhamFit> {
        if windows.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }

        let input = UwhamLogQInput::from_windows(windows)?;
        fit_log_q_input(input, &self.options)
    }

    pub fn estimate(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        self.fit(windows)?.result()
    }
}

pub(crate) fn fit_log_q_input(input: UwhamLogQInput, options: &UwhamOptions) -> Result<UwhamFit> {
    let free_energies = solve_uwham(&input, options)?;
    Ok(UwhamFit {
        input,
        free_energies,
        weights: OnceLock::new(),
        parallel: options.parallel,
    })
}

impl UwhamFit {
    pub fn n_observations(&self) -> usize {
        self.input.n_observations
    }

    pub fn n_states(&self) -> usize {
        self.input.n_states
    }

    pub fn states(&self) -> &[StatePoint] {
        &self.input.states
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.input.lambda_labels.as_deref()
    }

    pub fn free_energies(&self) -> &[f64] {
        &self.free_energies
    }

    pub fn weights(&self) -> &[f64] {
        self.weights.get_or_init(|| {
            let reduced = ReducedUwhamInput::from_pooled(&self.input)
                .expect("pooled UWHAM input should stay valid after fit");
            uwham_weights(&self.input, &reduced, &self.free_energies, self.parallel)
                .expect("UWHAM weights should be derivable from a converged fit")
        })
    }

    pub fn result(&self) -> Result<DeltaFMatrix> {
        delta_f_matrix_from_ze(
            &self.free_energies,
            self.input.states.clone(),
            self.input.lambda_labels.clone(),
        )
    }

    pub fn write_reference_inputs(&self, output_dir: impl AsRef<Path>) -> Result<()> {
        let output_dir = output_dir.as_ref();
        fs::create_dir_all(output_dir).map_err(|err| {
            CoreError::Unsupported(format!(
                "failed to create UWHAM reference directory {}: {err}",
                output_dir.display()
            ))
        })?;

        write_csv_matrix(
            output_dir.join("logQ.csv"),
            self.input.log_q.clone(),
            self.input.n_observations,
            self.input.n_states,
        )?;
        write_csv_vector(
            output_dir.join("size.csv"),
            self.input.size.iter().copied(),
            "size",
        )?;
        write_csv_matrix(
            output_dir.join("delta_f.csv"),
            {
                let result = self.result()?;
                result.values().to_vec()
            },
            self.input.n_states,
            self.input.n_states,
        )?;
        write_csv_vector(
            output_dir.join("ze.csv"),
            self.free_energies.iter().copied(),
            "ze",
        )?;
        write_metadata(output_dir, &self.input)?;
        Ok(())
    }
}

impl UwhamLogQInput {
    fn from_windows(windows: &[UNkMatrix]) -> Result<Self> {
        let states = ensure_consistent_states(windows)?;
        let lambda_labels = ensure_consistent_lambda_labels(windows)?;
        let n_states = states.len();
        let n_observations: usize = windows.iter().map(UNkMatrix::n_samples).sum();
        let mut size = vec![0.0; n_states];
        let mut log_q = Vec::with_capacity(n_observations * n_states);

        for window in windows {
            let sampled = window.sampled_state().ok_or_else(|| {
                CoreError::InvalidState("sampled_state required for UWHAM".to_string())
            })?;
            let sampled_idx = find_state_index_exact(&states, sampled)?;
            size[sampled_idx] += window.n_samples() as f64;

            for &value in window.data() {
                log_q.push(if value.is_infinite() {
                    f64::NEG_INFINITY
                } else {
                    -value
                });
            }
        }

        Ok(Self {
            log_q,
            n_observations,
            n_states,
            size,
            states,
            lambda_labels,
        })
    }

    pub(crate) fn from_atm(matrix: &AtmLogQMatrix) -> Self {
        Self {
            log_q: matrix.log_q().to_vec(),
            n_observations: matrix.n_observations(),
            n_states: matrix.n_states(),
            size: matrix
                .sampled_counts()
                .iter()
                .map(|&count| count as f64)
                .collect(),
            states: matrix.states().to_vec(),
            lambda_labels: matrix.lambda_labels().map(|labels| labels.to_vec()),
        }
    }
}

impl ReducedUwhamInput {
    fn from_pooled(input: &UwhamLogQInput) -> Result<Self> {
        let sampled_indices = input
            .size
            .iter()
            .enumerate()
            .filter_map(|(idx, &count)| (count > 0.0).then_some(idx))
            .collect::<Vec<_>>();
        if sampled_indices.is_empty() {
            return Err(CoreError::InvalidState(
                "UWHAM requires at least one sampled state".to_string(),
            ));
        }

        let base_index = sampled_indices[0];
        let base_sampled_index = 0;
        let mut log_q_relative = vec![0.0; input.log_q.len()];
        for observation_idx in 0..input.n_observations {
            let row = observation_row(&input.log_q, input.n_states, observation_idx);
            let base_value = row[base_index];
            let dest = observation_idx * input.n_states;
            if base_value == f64::NEG_INFINITY {
                return Err(CoreError::NonFiniteValue(format!(
                    "UWHAM baseline state {base_index} has -inf logQ at observation {observation_idx}"
                )));
            }
            for state_idx in 0..input.n_states {
                log_q_relative[dest + state_idx] = row[state_idx] - base_value;
            }
        }

        let mut rho = Vec::with_capacity(sampled_indices.len());
        let mut sampled_log_q_relative =
            Vec::with_capacity(input.n_observations * sampled_indices.len());
        for &state_idx in &sampled_indices {
            rho.push(input.size[state_idx] / input.n_observations as f64);
        }
        for observation_idx in 0..input.n_observations {
            let row = observation_row(&log_q_relative, input.n_states, observation_idx);
            for &state_idx in &sampled_indices {
                sampled_log_q_relative.push(row[state_idx]);
            }
        }

        Ok(Self {
            log_q_relative,
            sampled_indices,
            base_sampled_index,
            rho,
            sampled_log_q_relative,
        })
    }
}

pub(crate) fn solve_uwham(input: &UwhamLogQInput, options: &UwhamOptions) -> Result<Vec<f64>> {
    let reduced = ReducedUwhamInput::from_pooled(input)?;
    let sampled_count = reduced.sampled_indices.len();
    let variable_count = sampled_count.saturating_sub(1);
    let mut ze_nonbase = DVector::zeros(variable_count);

    if options.max_iterations == 0 {
        return Err(CoreError::ConvergenceFailure);
    }

    for _ in 0..options.max_iterations {
        let objective = evaluate_objective(input, &reduced, &ze_nonbase, options.parallel)?;
        let grad_norm = objective
            .gradient
            .iter()
            .fold(0.0_f64, |acc, value| acc.max(value.abs()));
        if grad_norm < options.tolerance {
            return finalize_free_energies(input, &reduced, &ze_nonbase, options.parallel);
        }

        let fallback_direction = -&objective.gradient;
        let mut directions = Vec::with_capacity(2);
        if let Some(direction) = newton_direction(&objective) {
            let step_norm = direction
                .iter()
                .fold(0.0_f64, |acc, value| acc.max(value.abs()));
            if step_norm < options.tolerance {
                return finalize_free_energies(input, &reduced, &ze_nonbase, options.parallel);
            }
            directions.push(direction);
        }
        directions.push(fallback_direction);

        let mut accepted = false;
        for candidate_direction in directions {
            let directional_derivative = objective.gradient.dot(&candidate_direction);
            if !directional_derivative.is_finite() || directional_derivative >= 0.0 {
                continue;
            }

            let mut step_scale = 1.0;
            while step_scale >= 1.0e-8 {
                let candidate = &ze_nonbase + step_scale * &candidate_direction;
                let candidate_objective =
                    evaluate_objective(input, &reduced, &candidate, options.parallel)?;
                if candidate_objective.value
                    <= objective.value + 1.0e-4 * step_scale * directional_derivative
                {
                    ze_nonbase = candidate;
                    accepted = true;
                    break;
                }
                step_scale *= 0.5;
            }

            if accepted {
                break;
            }
        }

        if !accepted {
            return Err(CoreError::ConvergenceFailure);
        }
    }

    Err(CoreError::ConvergenceFailure)
}

fn finalize_free_energies(
    input: &UwhamLogQInput,
    reduced: &ReducedUwhamInput,
    ze_nonbase: &DVector<f64>,
    parallel: bool,
) -> Result<Vec<f64>> {
    let sampled_ze = insert_base_state(ze_nonbase, reduced.base_sampled_index);
    let log_qsum = log_qsum_per_observation(input, reduced, &sampled_ze, parallel)?;
    let mut ze = vec![0.0; input.n_states];
    for (sampled_idx, &state_idx) in reduced.sampled_indices.iter().enumerate() {
        ze[state_idx] = sampled_ze[sampled_idx];
    }

    let mut z_mean = vec![0.0; input.n_states];
    for (observation_idx, &log_denom) in log_qsum.iter().enumerate() {
        let row = observation_row(&reduced.log_q_relative, input.n_states, observation_idx);
        for state_idx in 0..input.n_states {
            let contrib = if row[state_idx] == f64::NEG_INFINITY {
                0.0
            } else {
                (row[state_idx] - log_denom).exp()
            };
            z_mean[state_idx] += contrib;
        }
    }
    for value in &mut z_mean {
        *value /= input.n_observations as f64;
    }

    for state_idx in 0..input.n_states {
        if input.size[state_idx] > 0.0 {
            continue;
        }
        let mean = z_mean[state_idx];
        if !mean.is_finite() || mean <= 0.0 {
            return Err(CoreError::NonFiniteValue(format!(
                "UWHAM produced a non-positive unsampled-state mean weight for state {state_idx}"
            )));
        }
        ze[state_idx] = mean.ln();
    }

    Ok(ze)
}

fn evaluate_objective(
    input: &UwhamLogQInput,
    reduced: &ReducedUwhamInput,
    ze_nonbase: &DVector<f64>,
    parallel: bool,
) -> Result<UwhamObjective> {
    let sampled_ze = insert_base_state(ze_nonbase, reduced.base_sampled_index);
    let log_qsum = log_qsum_per_observation(input, reduced, &sampled_ze, parallel)?;
    let sampled_count = reduced.sampled_indices.len();
    let variable_count = sampled_count.saturating_sub(1);
    let mut value = 0.0;
    let mut col_sums = vec![0.0; variable_count];
    let mut gram = vec![0.0; variable_count * variable_count];

    for (observation_idx, &log_denom) in log_qsum.iter().enumerate() {
        value += log_denom;
        let row = observation_row(
            &reduced.sampled_log_q_relative,
            sampled_count,
            observation_idx,
        );
        let mut weights = vec![0.0; variable_count];
        let mut variable_idx = 0;
        for sampled_idx in 0..sampled_count {
            if sampled_idx == reduced.base_sampled_index {
                continue;
            }
            let contribution = if row[sampled_idx] == f64::NEG_INFINITY {
                0.0
            } else {
                (reduced.rho[sampled_idx].ln() + row[sampled_idx]
                    - sampled_ze[sampled_idx]
                    - log_denom)
                    .exp()
            };
            weights[variable_idx] = contribution;
            col_sums[variable_idx] += contribution;
            variable_idx += 1;
        }
        for i in 0..variable_count {
            for j in 0..variable_count {
                gram[i * variable_count + j] += weights[i] * weights[j];
            }
        }
    }

    value /= input.n_observations as f64;
    for (sampled_idx, &sampled_ze_value) in sampled_ze.iter().enumerate().take(sampled_count) {
        value += sampled_ze_value * reduced.rho[sampled_idx];
    }

    let mut gradient = DVector::zeros(variable_count);
    let mut variable_idx = 0;
    for sampled_idx in 0..sampled_count {
        if sampled_idx == reduced.base_sampled_index {
            continue;
        }
        gradient[variable_idx] =
            -col_sums[variable_idx] / input.n_observations as f64 + reduced.rho[sampled_idx];
        variable_idx += 1;
    }

    let mut hessian = DMatrix::zeros(variable_count, variable_count);
    for i in 0..variable_count {
        for j in 0..variable_count {
            hessian[(i, j)] = -gram[i * variable_count + j] / input.n_observations as f64;
        }
        hessian[(i, i)] += col_sums[i] / input.n_observations as f64;
    }

    if !value.is_finite() || gradient.iter().any(|value| !value.is_finite()) {
        return Err(CoreError::NonFiniteValue(
            "UWHAM objective became non-finite".to_string(),
        ));
    }

    Ok(UwhamObjective {
        value,
        gradient,
        hessian,
    })
}

fn newton_direction(objective: &UwhamObjective) -> Option<DVector<f64>> {
    if objective.gradient.is_empty() {
        return Some(DVector::zeros(0));
    }

    let identity = DMatrix::identity(objective.gradient.len(), objective.gradient.len());
    let mut damping = 0.0;
    for _ in 0..8 {
        let hessian: DMatrix<f64> = &objective.hessian + damping * &identity;
        if let Some(step) = hessian.lu().solve(&(-&objective.gradient)) {
            return Some(step);
        }
        damping = if damping == 0.0 {
            1.0e-8
        } else {
            damping * 10.0
        };
    }

    None
}

fn log_qsum_per_observation(
    input: &UwhamLogQInput,
    reduced: &ReducedUwhamInput,
    sampled_ze: &[f64],
    parallel: bool,
) -> Result<Vec<f64>> {
    let sampled_count = reduced.sampled_indices.len();
    let mut qsum = vec![0.0; input.n_observations];
    if parallel {
        qsum.par_iter_mut()
            .enumerate()
            .try_for_each(|(observation_idx, slot)| -> Result<()> {
                *slot = log_qsum_for_row(reduced, sampled_ze, sampled_count, observation_idx)?;
                Ok(())
            })?;
    } else {
        for (observation_idx, slot) in qsum.iter_mut().enumerate() {
            *slot = log_qsum_for_row(reduced, sampled_ze, sampled_count, observation_idx)?;
        }
    }
    Ok(qsum)
}

fn log_qsum_for_row(
    reduced: &ReducedUwhamInput,
    sampled_ze: &[f64],
    sampled_count: usize,
    observation_idx: usize,
) -> Result<f64> {
    let row = observation_row(
        &reduced.sampled_log_q_relative,
        sampled_count,
        observation_idx,
    );
    let mut max_term = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for (sampled_idx, &sampled_ze_value) in sampled_ze.iter().enumerate().take(sampled_count) {
        if row[sampled_idx] == f64::NEG_INFINITY {
            continue;
        }
        accumulate_logsumexp(
            &mut max_term,
            &mut sum,
            reduced.rho[sampled_idx].ln() + row[sampled_idx] - sampled_ze_value,
        );
    }
    if max_term == f64::NEG_INFINITY || !sum.is_finite() {
        return Err(CoreError::NonFiniteValue(format!(
            "UWHAM denominator became non-finite at observation {observation_idx}"
        )));
    }
    Ok(max_term + sum.ln())
}

fn insert_base_state(ze_nonbase: &DVector<f64>, base_index: usize) -> Vec<f64> {
    let mut ze = Vec::with_capacity(ze_nonbase.len() + 1);
    let mut source_idx = 0;
    for sampled_idx in 0..=ze_nonbase.len() {
        if sampled_idx == base_index {
            ze.push(0.0);
        } else {
            ze.push(ze_nonbase[source_idx]);
            source_idx += 1;
        }
    }
    ze
}

fn uwham_weights(
    input: &UwhamLogQInput,
    reduced: &ReducedUwhamInput,
    ze: &[f64],
    parallel: bool,
) -> Result<Vec<f64>> {
    let sampled_ze = reduced
        .sampled_indices
        .iter()
        .map(|&state_idx| ze[state_idx])
        .collect::<Vec<_>>();
    let log_qsum = log_qsum_per_observation(input, reduced, &sampled_ze, parallel)?;
    let mut weights = vec![0.0; input.n_observations * input.n_states];
    if parallel {
        weights
            .par_chunks_mut(input.n_states)
            .enumerate()
            .try_for_each(|(observation_idx, row)| -> Result<()> {
                fill_weights_row(reduced, ze, &log_qsum, observation_idx, row);
                Ok(())
            })?;
    } else {
        for (observation_idx, row) in weights.chunks_mut(input.n_states).enumerate() {
            fill_weights_row(reduced, ze, &log_qsum, observation_idx, row);
        }
    }
    Ok(weights)
}

fn fill_weights_row(
    reduced: &ReducedUwhamInput,
    ze: &[f64],
    log_qsum: &[f64],
    observation_idx: usize,
    row: &mut [f64],
) {
    let log_denom = log_qsum[observation_idx];
    let relative_row = observation_row(&reduced.log_q_relative, row.len(), observation_idx);
    for state_idx in 0..row.len() {
        row[state_idx] = if relative_row[state_idx] == f64::NEG_INFINITY {
            0.0
        } else {
            (relative_row[state_idx] - ze[state_idx] - log_denom).exp()
        };
    }
}

fn accumulate_logsumexp(max_term: &mut f64, sum: &mut f64, value: f64) {
    if value == f64::NEG_INFINITY {
        return;
    }
    if *max_term == f64::NEG_INFINITY {
        *max_term = value;
        *sum = 1.0;
        return;
    }
    if value <= *max_term {
        *sum += (value - *max_term).exp();
    } else {
        *sum = *sum * (*max_term - value).exp() + 1.0;
        *max_term = value;
    }
}

fn delta_f_matrix_from_ze(
    ze: &[f64],
    states: Vec<StatePoint>,
    lambda_labels: Option<Vec<String>>,
) -> Result<DeltaFMatrix> {
    let n_states = states.len();
    let mut values = vec![0.0; n_states * n_states];
    for i in 0..n_states {
        for j in 0..n_states {
            values[i * n_states + j] = ze[j] - ze[i];
        }
    }
    DeltaFMatrix::new_with_labels(values, None, n_states, states, lambda_labels)
}

fn observation_row(values: &[f64], n_states: usize, observation_idx: usize) -> &[f64] {
    let start = observation_idx * n_states;
    &values[start..start + n_states]
}

fn write_csv_matrix(
    path: impl AsRef<Path>,
    values: Vec<f64>,
    n_rows: usize,
    n_cols: usize,
) -> Result<()> {
    let expected = n_rows.checked_mul(n_cols).ok_or(CoreError::InvalidShape {
        expected: n_rows,
        found: n_cols,
    })?;
    if values.len() != expected {
        return Err(CoreError::InvalidShape {
            expected,
            found: values.len(),
        });
    }
    let mut content = String::new();
    for row in 0..n_rows {
        for col in 0..n_cols {
            if col > 0 {
                content.push(',');
            }
            content.push_str(&format!("{}", values[row * n_cols + col]));
        }
        content.push('\n');
    }
    fs::write(path.as_ref(), content).map_err(|err| {
        CoreError::Unsupported(format!(
            "failed to write UWHAM reference CSV {}: {err}",
            path.as_ref().display()
        ))
    })?;
    Ok(())
}

fn write_csv_vector(
    path: impl AsRef<Path>,
    values: impl IntoIterator<Item = f64>,
    header: &str,
) -> Result<()> {
    let mut content = String::new();
    content.push_str(header);
    content.push('\n');
    for value in values {
        content.push_str(&format!("{}\n", value));
    }
    fs::write(path.as_ref(), content).map_err(|err| {
        CoreError::Unsupported(format!(
            "failed to write UWHAM reference CSV {}: {err}",
            path.as_ref().display()
        ))
    })?;
    Ok(())
}

fn write_metadata(output_dir: &Path, input: &UwhamLogQInput) -> Result<()> {
    let mut content = String::new();
    content.push_str("# UWHAM reference input export\n");
    content.push_str("# Generated from alchemrs fixture-backed windows.\n");
    content.push_str("# Files:\n");
    content.push_str("# - logQ.csv: pooled log unnormalized densities (R UWHAM input)\n");
    content.push_str("# - size.csv: sample counts per sampled state\n");
    content.push_str("# - delta_f.csv: alchemrs reference pairwise free energies\n");
    content.push_str("# - ze.csv: alchemrs reference log normalizing constants\n\n");
    content.push_str("states\n");
    for state in &input.states {
        let lambdas = state
            .lambdas()
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",");
        content.push_str(&format!("{lambdas},{}\n", state.temperature_k()));
    }
    if let Some(labels) = input.lambda_labels.as_ref() {
        content.push_str("\nlambda_labels\n");
        for label in labels {
            content.push_str(label);
            content.push('\n');
        }
    }
    fs::write(output_dir.join("README.txt"), content).map_err(|err| {
        CoreError::Unsupported(format!(
            "failed to write UWHAM reference metadata {}: {err}",
            output_dir.display()
        ))
    })?;
    Ok(())
}
