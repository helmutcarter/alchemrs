use crate::data::{
    DeltaFMatrix, DhdlSeries, FreeEnergyEstimate, OverlapMatrix, StatePoint, UNkMatrix,
};
use crate::error::{CoreError, Result};
use crate::estimators::{
    mbar_log_weights_from_windows, BarEstimator, BarOptions, ExpEstimator, ExpOptions,
    MbarEstimator, MbarOptions, TiEstimator, TiOptions,
};
use nalgebra::{DMatrix, Schur};
use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq)]
pub struct ConvergencePoint {
    n_windows: usize,
    delta_f: f64,
    uncertainty: Option<f64>,
    from_state: StatePoint,
    to_state: StatePoint,
    lambda_labels: Option<Vec<String>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlockEstimate {
    block_index: usize,
    n_blocks: usize,
    delta_f: f64,
    uncertainty: Option<f64>,
    from_state: StatePoint,
    to_state: StatePoint,
    lambda_labels: Option<Vec<String>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdvisorEstimator {
    Mbar,
    Bar,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScheduleAdvisorOptions {
    pub estimator: AdvisorEstimator,
    pub overlap_min: f64,
    pub block_cv_min: f64,
    pub n_blocks: usize,
    pub suggest_midpoints: bool,
}

impl Default for ScheduleAdvisorOptions {
    fn default() -> Self {
        Self {
            estimator: AdvisorEstimator::Mbar,
            overlap_min: 0.03,
            block_cv_min: 0.15,
            n_blocks: 4,
            suggest_midpoints: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeSeverity {
    Healthy,
    Monitor,
    AddSampling,
    AddWindow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SuggestionKind {
    NoChange,
    ExtendSampling,
    InsertWindow,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdjacentEdgeDiagnostic {
    edge_index: usize,
    from_state: StatePoint,
    to_state: StatePoint,
    lambda_labels: Option<Vec<String>>,
    lambda_delta: Vec<f64>,
    overlap_forward: f64,
    overlap_reverse: f64,
    overlap_min: f64,
    delta_f: f64,
    uncertainty: Option<f64>,
    block_mean: Option<f64>,
    block_stddev: Option<f64>,
    block_cv: Option<f64>,
    severity: EdgeSeverity,
}

impl AdjacentEdgeDiagnostic {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        edge_index: usize,
        from_state: StatePoint,
        to_state: StatePoint,
        lambda_labels: Option<Vec<String>>,
        overlap_forward: f64,
        overlap_reverse: f64,
        delta_f: f64,
        uncertainty: Option<f64>,
        block_mean: Option<f64>,
        block_stddev: Option<f64>,
        block_cv: Option<f64>,
        severity: EdgeSeverity,
    ) -> Result<Self> {
        if !overlap_forward.is_finite() || !overlap_reverse.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "adjacent overlap values must be finite".to_string(),
            ));
        }
        if !delta_f.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "adjacent delta_f must be finite".to_string(),
            ));
        }
        if let Some(value) = uncertainty {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "adjacent uncertainty must be finite".to_string(),
                ));
            }
        }
        if let Some(value) = block_mean {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "block_mean must be finite".to_string(),
                ));
            }
        }
        if let Some(value) = block_stddev {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "block_stddev must be finite".to_string(),
                ));
            }
        }
        if let Some(value) = block_cv {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "block_cv must be finite".to_string(),
                ));
            }
        }
        let lambda_delta = to_state
            .lambdas()
            .iter()
            .zip(from_state.lambdas().iter())
            .map(|(to, from)| to - from)
            .collect();
        Ok(Self {
            edge_index,
            from_state,
            to_state,
            lambda_labels,
            lambda_delta,
            overlap_forward,
            overlap_reverse,
            overlap_min: overlap_forward.min(overlap_reverse),
            delta_f,
            uncertainty,
            block_mean,
            block_stddev,
            block_cv,
            severity,
        })
    }

    pub fn edge_index(&self) -> usize {
        self.edge_index
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }

    pub fn lambda_delta(&self) -> &[f64] {
        &self.lambda_delta
    }

    pub fn overlap_forward(&self) -> f64 {
        self.overlap_forward
    }

    pub fn overlap_reverse(&self) -> f64 {
        self.overlap_reverse
    }

    pub fn overlap_min(&self) -> f64 {
        self.overlap_min
    }

    pub fn delta_f(&self) -> f64 {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn block_mean(&self) -> Option<f64> {
        self.block_mean
    }

    pub fn block_stddev(&self) -> Option<f64> {
        self.block_stddev
    }

    pub fn block_cv(&self) -> Option<f64> {
        self.block_cv
    }

    pub fn severity(&self) -> EdgeSeverity {
        self.severity
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ScheduleSuggestion {
    kind: SuggestionKind,
    edge_index: usize,
    from_state: StatePoint,
    to_state: StatePoint,
    proposed_state: Option<StatePoint>,
    reason: String,
}

impl ScheduleSuggestion {
    pub(crate) fn new(
        kind: SuggestionKind,
        edge_index: usize,
        from_state: StatePoint,
        to_state: StatePoint,
        proposed_state: Option<StatePoint>,
        reason: String,
    ) -> Self {
        Self {
            kind,
            edge_index,
            from_state,
            to_state,
            proposed_state,
            reason,
        }
    }

    pub fn kind(&self) -> SuggestionKind {
        self.kind
    }

    pub fn edge_index(&self) -> usize {
        self.edge_index
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }

    pub fn proposed_state(&self) -> Option<&StatePoint> {
        self.proposed_state.as_ref()
    }

    pub fn reason(&self) -> &str {
        &self.reason
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ScheduleAdvice {
    edges: Vec<AdjacentEdgeDiagnostic>,
    suggestions: Vec<ScheduleSuggestion>,
}

impl ScheduleAdvice {
    pub fn new(edges: Vec<AdjacentEdgeDiagnostic>, suggestions: Vec<ScheduleSuggestion>) -> Self {
        Self { edges, suggestions }
    }

    pub fn edges(&self) -> &[AdjacentEdgeDiagnostic] {
        &self.edges
    }

    pub fn suggestions(&self) -> &[ScheduleSuggestion] {
        &self.suggestions
    }
}

impl BlockEstimate {
    pub fn new(
        block_index: usize,
        n_blocks: usize,
        delta_f: f64,
        uncertainty: Option<f64>,
        from_state: StatePoint,
        to_state: StatePoint,
        lambda_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        if n_blocks == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        if block_index >= n_blocks {
            return Err(CoreError::InvalidShape {
                expected: n_blocks,
                found: block_index + 1,
            });
        }
        if !delta_f.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "delta_f must be finite".to_string(),
            ));
        }
        if let Some(value) = uncertainty {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "uncertainty must be finite".to_string(),
                ));
            }
        }
        Ok(Self {
            block_index,
            n_blocks,
            delta_f,
            uncertainty,
            from_state,
            to_state,
            lambda_labels,
        })
    }

    pub fn block_index(&self) -> usize {
        self.block_index
    }

    pub fn n_blocks(&self) -> usize {
        self.n_blocks
    }

    pub fn delta_f(&self) -> f64 {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

impl ConvergencePoint {
    pub fn new(
        n_windows: usize,
        delta_f: f64,
        uncertainty: Option<f64>,
        from_state: StatePoint,
        to_state: StatePoint,
        lambda_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        if n_windows == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        if !delta_f.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "delta_f must be finite".to_string(),
            ));
        }
        if let Some(value) = uncertainty {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "uncertainty must be finite".to_string(),
                ));
            }
        }
        Ok(Self {
            n_windows,
            delta_f,
            uncertainty,
            from_state,
            to_state,
            lambda_labels,
        })
    }

    pub fn n_windows(&self) -> usize {
        self.n_windows
    }

    pub fn delta_f(&self) -> f64 {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

pub fn overlap_matrix(
    windows: &[UNkMatrix],
    options: Option<MbarOptions>,
) -> Result<OverlapMatrix> {
    let options = options.unwrap_or_default();
    let lambda_labels = windows
        .first()
        .and_then(|window| window.lambda_labels().map(|labels| labels.to_vec()));
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

    OverlapMatrix::new_with_labels(values, n_states, states, lambda_labels)
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

pub fn ti_convergence(
    series: &[DhdlSeries],
    options: Option<TiOptions>,
) -> Result<Vec<ConvergencePoint>> {
    if series.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: series.len(),
        });
    }
    let estimator = TiEstimator::new(options.unwrap_or_default());
    let mut points = Vec::with_capacity(series.len() - 1);
    for end in 2..=series.len() {
        let result = estimator.fit(&series[..end])?;
        points.push(convergence_point_from_scalar(end, &result)?);
    }
    Ok(points)
}

pub fn bar_convergence(
    windows: &[UNkMatrix],
    options: Option<BarOptions>,
) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| BarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub fn mbar_convergence(
    windows: &[UNkMatrix],
    options: Option<MbarOptions>,
) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| MbarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        1,
    )
}

pub fn exp_convergence(
    windows: &[UNkMatrix],
    options: Option<ExpOptions>,
) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub fn dexp_convergence(
    windows: &[UNkMatrix],
    options: Option<ExpOptions>,
) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        true,
        2,
    )
}

pub fn advise_lambda_schedule(
    windows: &[UNkMatrix],
    options: Option<ScheduleAdvisorOptions>,
) -> Result<ScheduleAdvice> {
    if windows.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: windows.len(),
        });
    }

    let options = validate_schedule_advisor_options(options.unwrap_or_default())?;
    let ordered = sort_windows_by_sampled_state(windows)?;
    let trimmed = trim_windows_to_sampled_states(&ordered)?;
    let overlap = overlap_matrix(&trimmed, Some(MbarOptions::default()))?;
    let n_states = overlap.n_states();
    let labels = overlap.lambda_labels().map(|values| values.to_vec());

    let mut edges = Vec::with_capacity(n_states - 1);
    let mut suggestions = Vec::new();
    for edge_index in 0..(n_states - 1) {
        let pair_windows = trim_windows_to_sampled_states(&[
            ordered[edge_index].clone(),
            ordered[edge_index + 1].clone(),
        ])?;
        let pair_result = fit_pair(&pair_windows, options.estimator)?;
        let block_stats = pair_block_stats(&pair_windows, options)?;
        let overlap_forward = overlap.values()[edge_index * n_states + edge_index + 1];
        let overlap_reverse = overlap.values()[(edge_index + 1) * n_states + edge_index];
        let severity = classify_edge(
            overlap_forward.min(overlap_reverse),
            block_stats.map(|(_, _, cv)| cv),
            &options,
        );
        let diagnostic = AdjacentEdgeDiagnostic::new(
            edge_index,
            pair_result.states()[0].clone(),
            pair_result.states()[1].clone(),
            labels.clone(),
            overlap_forward,
            overlap_reverse,
            pair_result.values()[1],
            pair_result.uncertainties().map(|values| values[1]),
            block_stats.map(|(mean, _, _)| mean),
            block_stats.map(|(_, stddev, _)| stddev),
            block_stats.map(|(_, _, cv)| cv),
            severity,
        )?;

        if let Some(suggestion) = schedule_suggestion(&diagnostic, &options)? {
            suggestions.push(suggestion);
        }
        edges.push(diagnostic);
    }

    Ok(ScheduleAdvice::new(edges, suggestions))
}

pub(crate) fn ti_block_average(
    series: &[DhdlSeries],
    n_blocks: usize,
    options: Option<TiOptions>,
) -> Result<Vec<BlockEstimate>> {
    if series.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: series.len(),
        });
    }
    if n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let estimator = TiEstimator::new(options.unwrap_or_default());
    let blocked = series
        .iter()
        .map(|item| split_dhdl_series(item, n_blocks))
        .collect::<Result<Vec<_>>>()?;

    let mut points = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let block = blocked
            .iter()
            .map(|chunks| chunks[block_index].clone())
            .collect::<Vec<_>>();
        let result = estimator.fit(&block)?;
        points.push(block_estimate_from_scalar(block_index, n_blocks, &result)?);
    }
    Ok(points)
}

pub(crate) fn mbar_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<MbarOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| MbarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        1,
    )
}

pub(crate) fn bar_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<BarOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| BarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub(crate) fn exp_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<ExpOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub(crate) fn dexp_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<ExpOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        true,
        2,
    )
}

fn convergence_from_matrix_windows<F>(
    windows: &[UNkMatrix],
    fit: F,
    reverse: bool,
    minimum_windows: usize,
) -> Result<Vec<ConvergencePoint>>
where
    F: Fn(&[UNkMatrix]) -> Result<DeltaFMatrix>,
{
    if windows.len() < minimum_windows {
        return Err(CoreError::InvalidShape {
            expected: minimum_windows,
            found: windows.len(),
        });
    }
    let mut points = Vec::with_capacity(windows.len() - minimum_windows + 1);
    for count in minimum_windows..=windows.len() {
        let subset = if reverse {
            &windows[windows.len() - count..]
        } else {
            &windows[..count]
        };
        let trimmed = trim_windows_to_sampled_states(subset)?;
        let result = fit(&trimmed)?;
        points.push(convergence_point_from_matrix(count, &result, reverse)?);
    }
    Ok(points)
}

fn block_average_from_windows<F>(
    windows: &[UNkMatrix],
    n_blocks: usize,
    fit: F,
    reverse: bool,
    minimum_windows: usize,
) -> Result<Vec<BlockEstimate>>
where
    F: Fn(&[UNkMatrix]) -> Result<DeltaFMatrix>,
{
    if windows.len() < minimum_windows {
        return Err(CoreError::InvalidShape {
            expected: minimum_windows,
            found: windows.len(),
        });
    }
    if n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let blocked = windows
        .iter()
        .map(|window| split_u_nk_window(window, n_blocks))
        .collect::<Result<Vec<_>>>()?;

    let mut points = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let subset = blocked
            .iter()
            .map(|chunks| chunks[block_index].clone())
            .collect::<Vec<_>>();
        let trimmed = trim_windows_to_sampled_states(&subset)?;
        let result = fit(&trimmed)?;
        points.push(block_estimate_from_matrix(
            block_index,
            n_blocks,
            &result,
            reverse,
        )?);
    }
    Ok(points)
}

fn trim_windows_to_sampled_states(windows: &[UNkMatrix]) -> Result<Vec<UNkMatrix>> {
    let states = windows
        .iter()
        .map(|window| {
            window.sampled_state().cloned().ok_or_else(|| {
                CoreError::InvalidState(
                    "sampled_state required for convergence analysis".to_string(),
                )
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut trimmed = Vec::with_capacity(windows.len());
    for window in windows {
        let indices = states
            .iter()
            .map(|state| {
                window
                    .evaluated_states()
                    .iter()
                    .position(|candidate| candidate == state)
                    .ok_or_else(|| {
                        CoreError::InvalidState(
                            "sampled_state not found in evaluated_states".to_string(),
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;

        let mut data = Vec::with_capacity(window.n_samples() * states.len());
        for sample_idx in 0..window.n_samples() {
            let row_offset = sample_idx * window.n_states();
            for &idx in &indices {
                data.push(window.data()[row_offset + idx]);
            }
        }

        trimmed.push(UNkMatrix::new_with_labels(
            window.n_samples(),
            states.len(),
            data,
            window.time_ps().to_vec(),
            window.sampled_state().cloned(),
            states.clone(),
            window.lambda_labels().map(|labels| labels.to_vec()),
        )?);
    }

    Ok(trimmed)
}

fn split_dhdl_series(series: &DhdlSeries, n_blocks: usize) -> Result<Vec<DhdlSeries>> {
    let block_len = block_length(series.values().len(), n_blocks)?;
    let mut out = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let start = block_index * block_len;
        let end = start + block_len;
        out.push(DhdlSeries::new(
            series.state().clone(),
            series.time_ps()[start..end].to_vec(),
            series.values()[start..end].to_vec(),
        )?);
    }
    Ok(out)
}

fn split_u_nk_window(window: &UNkMatrix, n_blocks: usize) -> Result<Vec<UNkMatrix>> {
    let block_len = block_length(window.n_samples(), n_blocks)?;
    let row_width = window.n_states();
    let mut out = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let start = block_index * block_len;
        let end = start + block_len;
        let mut data = Vec::with_capacity(block_len * row_width);
        for sample_idx in start..end {
            let row_start = sample_idx * row_width;
            data.extend_from_slice(&window.data()[row_start..row_start + row_width]);
        }
        out.push(UNkMatrix::new_with_labels(
            block_len,
            row_width,
            data,
            window.time_ps()[start..end].to_vec(),
            window.sampled_state().cloned(),
            window.evaluated_states().to_vec(),
            window.lambda_labels().map(|labels| labels.to_vec()),
        )?);
    }
    Ok(out)
}

fn block_length(n_samples: usize, n_blocks: usize) -> Result<usize> {
    if n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let block_len = n_samples / n_blocks;
    if block_len == 0 {
        return Err(CoreError::InvalidShape {
            expected: n_blocks,
            found: n_samples,
        });
    }
    Ok(block_len)
}

fn validate_schedule_advisor_options(
    options: ScheduleAdvisorOptions,
) -> Result<ScheduleAdvisorOptions> {
    if !options.overlap_min.is_finite() || options.overlap_min < 0.0 {
        return Err(CoreError::InvalidState(
            "schedule advisor overlap_min must be finite and non-negative".to_string(),
        ));
    }
    if !options.block_cv_min.is_finite() || options.block_cv_min < 0.0 {
        return Err(CoreError::InvalidState(
            "schedule advisor block_cv_min must be finite and non-negative".to_string(),
        ));
    }
    if options.n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    Ok(options)
}

fn sort_windows_by_sampled_state(windows: &[UNkMatrix]) -> Result<Vec<UNkMatrix>> {
    for window in windows {
        if window.sampled_state().is_none() {
            return Err(CoreError::InvalidState(
                "sampled_state required for schedule advisor".to_string(),
            ));
        }
    }
    let mut ordered = windows.to_vec();
    ordered.sort_by(|left, right| {
        compare_state_points(
            left.sampled_state().expect("sampled_state checked above"),
            right.sampled_state().expect("sampled_state checked above"),
        )
    });
    for idx in 1..ordered.len() {
        if ordered[idx - 1].sampled_state() == ordered[idx].sampled_state() {
            return Err(CoreError::InvalidState(
                "multiple windows for same sampled_state".to_string(),
            ));
        }
    }
    Ok(ordered)
}

fn compare_state_points(left: &StatePoint, right: &StatePoint) -> Ordering {
    let len_order = left.lambdas().len().cmp(&right.lambdas().len());
    if len_order != Ordering::Equal {
        return len_order;
    }
    for (lhs, rhs) in left.lambdas().iter().zip(right.lambdas().iter()) {
        match lhs.partial_cmp(rhs).unwrap_or(Ordering::Equal) {
            Ordering::Equal => continue,
            other => return other,
        }
    }
    Ordering::Equal
}

fn fit_pair(windows: &[UNkMatrix], estimator: AdvisorEstimator) -> Result<DeltaFMatrix> {
    match estimator {
        AdvisorEstimator::Mbar => MbarEstimator::default().fit(windows),
        AdvisorEstimator::Bar => BarEstimator::default().fit(windows),
    }
}

fn pair_block_stats(
    windows: &[UNkMatrix],
    options: ScheduleAdvisorOptions,
) -> Result<Option<(f64, f64, f64)>> {
    let blocks = match options.estimator {
        AdvisorEstimator::Mbar => MbarEstimator::default().block_average(windows, options.n_blocks),
        AdvisorEstimator::Bar => BarEstimator::default().block_average(windows, options.n_blocks),
    };
    let blocks = match blocks {
        Ok(blocks) => blocks,
        Err(CoreError::InvalidShape { .. }) => return Ok(None),
        Err(error) => return Err(error),
    };
    if blocks.is_empty() {
        return Ok(None);
    }
    let values = blocks
        .iter()
        .map(BlockEstimate::delta_f)
        .collect::<Vec<_>>();
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / values.len() as f64;
    let stddev = variance.sqrt();
    let cv = if mean.abs() > 1.0e-12 {
        stddev / mean.abs()
    } else {
        stddev
    };
    Ok(Some((mean, stddev, cv)))
}

fn classify_edge(
    overlap_min: f64,
    block_cv: Option<f64>,
    options: &ScheduleAdvisorOptions,
) -> EdgeSeverity {
    if overlap_min < options.overlap_min {
        EdgeSeverity::AddWindow
    } else if block_cv.is_some_and(|value| value >= options.block_cv_min) {
        EdgeSeverity::AddSampling
    } else if (options.overlap_min > 0.0 && overlap_min < options.overlap_min * 1.5)
        || block_cv.is_some_and(|value| value >= options.block_cv_min * 0.75)
    {
        EdgeSeverity::Monitor
    } else {
        EdgeSeverity::Healthy
    }
}

fn schedule_suggestion(
    diagnostic: &AdjacentEdgeDiagnostic,
    options: &ScheduleAdvisorOptions,
) -> Result<Option<ScheduleSuggestion>> {
    let suggestion = match diagnostic.severity() {
        EdgeSeverity::AddWindow => {
            let midpoint = if options.suggest_midpoints {
                Some(midpoint_state(
                    diagnostic.from_state(),
                    diagnostic.to_state(),
                )?)
            } else {
                None
            };
            Some(ScheduleSuggestion::new(
                SuggestionKind::InsertWindow,
                diagnostic.edge_index(),
                diagnostic.from_state().clone(),
                diagnostic.to_state().clone(),
                midpoint,
                format!(
                    "adjacent overlap {:.6} is below threshold {:.6}",
                    diagnostic.overlap_min(),
                    options.overlap_min
                ),
            ))
        }
        EdgeSeverity::AddSampling => Some(ScheduleSuggestion::new(
            SuggestionKind::ExtendSampling,
            diagnostic.edge_index(),
            diagnostic.from_state().clone(),
            diagnostic.to_state().clone(),
            None,
            format!(
                "block coefficient of variation {:.6} exceeds threshold {:.6}",
                diagnostic.block_cv().unwrap_or_default(),
                options.block_cv_min
            ),
        )),
        EdgeSeverity::Healthy | EdgeSeverity::Monitor => None,
    };
    Ok(suggestion)
}

fn midpoint_state(from_state: &StatePoint, to_state: &StatePoint) -> Result<StatePoint> {
    let lambdas = from_state
        .lambdas()
        .iter()
        .zip(to_state.lambdas().iter())
        .map(|(from, to)| 0.5 * (from + to))
        .collect::<Vec<_>>();
    StatePoint::new(lambdas, from_state.temperature_k())
}

fn convergence_point_from_scalar(
    n_windows: usize,
    result: &FreeEnergyEstimate,
) -> Result<ConvergencePoint> {
    ConvergencePoint::new(
        n_windows,
        result.delta_f(),
        result.uncertainty(),
        result.from_state().clone(),
        result.to_state().clone(),
        None,
    )
}

fn block_estimate_from_scalar(
    block_index: usize,
    n_blocks: usize,
    result: &FreeEnergyEstimate,
) -> Result<BlockEstimate> {
    BlockEstimate::new(
        block_index,
        n_blocks,
        result.delta_f(),
        result.uncertainty(),
        result.from_state().clone(),
        result.to_state().clone(),
        None,
    )
}

fn convergence_point_from_matrix(
    n_windows: usize,
    result: &DeltaFMatrix,
    reverse: bool,
) -> Result<ConvergencePoint> {
    let n_states = result.n_states();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let index = if reverse {
        (n_states - 1) * n_states
    } else {
        n_states - 1
    };
    let (from_state, to_state) = if reverse {
        (
            result.states().last().unwrap().clone(),
            result.states().first().unwrap().clone(),
        )
    } else {
        (
            result.states().first().unwrap().clone(),
            result.states().last().unwrap().clone(),
        )
    };
    ConvergencePoint::new(
        n_windows,
        result.values()[index],
        result.uncertainties().map(|values| values[index]),
        from_state,
        to_state,
        result.lambda_labels().map(|labels| labels.to_vec()),
    )
}

fn block_estimate_from_matrix(
    block_index: usize,
    n_blocks: usize,
    result: &DeltaFMatrix,
    reverse: bool,
) -> Result<BlockEstimate> {
    let n_states = result.n_states();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let index = if reverse {
        (n_states - 1) * n_states
    } else {
        n_states - 1
    };
    BlockEstimate::new(
        block_index,
        n_blocks,
        result.values()[index],
        result.uncertainties().map(|values| values[index]),
        if reverse {
            result.states().last().unwrap().clone()
        } else {
            result.states().first().unwrap().clone()
        },
        if reverse {
            result.states().first().unwrap().clone()
        } else {
            result.states().last().unwrap().clone()
        },
        result.lambda_labels().map(|labels| labels.to_vec()),
    )
}

#[cfg(test)]
mod tests {
    use super::{advise_lambda_schedule, AdvisorEstimator, EdgeSeverity, ScheduleAdvisorOptions};
    use crate::data::{StatePoint, UNkMatrix};

    fn build_window(
        sampled_lambda: f64,
        rows: &[[f64; 3]],
        evaluated_states: &[StatePoint],
    ) -> UNkMatrix {
        let n_samples = rows.len();
        let data = rows
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let time = (0..n_samples).map(|idx| idx as f64).collect::<Vec<_>>();
        UNkMatrix::new(
            n_samples,
            3,
            data,
            time,
            Some(StatePoint::new(vec![sampled_lambda], 298.0).unwrap()),
            evaluated_states.to_vec(),
        )
        .unwrap()
    }

    #[test]
    fn advise_lambda_schedule_sorts_windows_and_suggests_midpoints() {
        let states = vec![
            StatePoint::new(vec![0.0], 298.0).unwrap(),
            StatePoint::new(vec![0.5], 298.0).unwrap(),
            StatePoint::new(vec![1.0], 298.0).unwrap(),
        ];
        let windows = vec![
            build_window(
                1.0,
                &[
                    [2.8, 0.2, 0.0],
                    [3.0, 0.3, 0.0],
                    [2.9, 0.2, 0.0],
                    [3.1, 0.4, 0.0],
                ],
                &states,
            ),
            build_window(
                0.0,
                &[
                    [0.0, 0.2, 2.8],
                    [0.0, 0.3, 3.0],
                    [0.0, 0.2, 2.9],
                    [0.0, 0.4, 3.1],
                ],
                &states,
            ),
            build_window(
                0.5,
                &[
                    [0.2, 0.0, 0.2],
                    [0.3, 0.0, 0.3],
                    [0.2, 0.0, 0.2],
                    [0.4, 0.0, 0.4],
                ],
                &states,
            ),
        ];

        let advice = advise_lambda_schedule(
            &windows,
            Some(ScheduleAdvisorOptions {
                estimator: AdvisorEstimator::Mbar,
                overlap_min: 0.99,
                block_cv_min: 1.0e30,
                n_blocks: 4,
                suggest_midpoints: true,
            }),
        )
        .unwrap();

        assert_eq!(advice.edges().len(), 2);
        assert_eq!(advice.edges()[0].from_state().lambdas(), &[0.0]);
        assert_eq!(advice.edges()[0].to_state().lambdas(), &[0.5]);
        assert_eq!(advice.edges()[0].severity(), EdgeSeverity::AddWindow);
        assert_eq!(advice.suggestions().len(), 2);
        assert_eq!(
            advice.suggestions()[0].proposed_state().unwrap().lambdas(),
            &[0.25]
        );
    }

    #[test]
    fn advise_lambda_schedule_can_return_healthy_edges() {
        let states = vec![
            StatePoint::new(vec![0.0], 298.0).unwrap(),
            StatePoint::new(vec![0.5], 298.0).unwrap(),
            StatePoint::new(vec![1.0], 298.0).unwrap(),
        ];
        let windows = vec![
            build_window(
                0.0,
                &[
                    [0.0, 0.2, 2.8],
                    [0.0, 0.3, 3.0],
                    [0.0, 0.2, 2.9],
                    [0.0, 0.4, 3.1],
                ],
                &states,
            ),
            build_window(
                0.5,
                &[
                    [0.2, 0.0, 0.2],
                    [0.3, 0.0, 0.3],
                    [0.2, 0.0, 0.2],
                    [0.4, 0.0, 0.4],
                ],
                &states,
            ),
            build_window(
                1.0,
                &[
                    [2.8, 0.2, 0.0],
                    [3.0, 0.3, 0.0],
                    [2.9, 0.2, 0.0],
                    [3.1, 0.4, 0.0],
                ],
                &states,
            ),
        ];

        let advice = advise_lambda_schedule(
            &windows,
            Some(ScheduleAdvisorOptions {
                estimator: AdvisorEstimator::Bar,
                overlap_min: 0.0,
                block_cv_min: 1.0e30,
                n_blocks: 4,
                suggest_midpoints: false,
            }),
        )
        .unwrap();

        assert_eq!(advice.edges().len(), 2);
        assert!(advice.suggestions().is_empty());
        assert!(advice
            .edges()
            .iter()
            .all(|edge| edge.severity() == EdgeSeverity::Healthy));
    }
}
