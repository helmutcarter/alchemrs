use crate::data::{
    find_state_index_exact, state_points_match_exact, DeltaFMatrix, DhdlSeries, FreeEnergyEstimate,
    OverlapMatrix, StatePoint, UNkMatrix,
};
use crate::error::{CoreError, Result};
use crate::estimators::{
    BarEstimator, BarOptions, ExpEstimator, ExpOptions, IntegrationMethod, MbarEstimator,
    MbarOptions, TiEstimator, TiOptions,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProposalStrategy {
    Midpoint,
    FocusedSplit,
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
    relative_overlap: Option<f64>,
    relative_uncertainty: Option<f64>,
    dominant_components: Vec<String>,
    priority_score: f64,
    severity: EdgeSeverity,
}

#[derive(Debug, Clone)]
struct EdgeMetrics {
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
        relative_overlap: Option<f64>,
        relative_uncertainty: Option<f64>,
        dominant_components: Vec<String>,
        priority_score: f64,
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
        if let Some(value) = relative_overlap {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "relative_overlap must be finite".to_string(),
                ));
            }
        }
        if let Some(value) = relative_uncertainty {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "relative_uncertainty must be finite".to_string(),
                ));
            }
        }
        if !priority_score.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "priority_score must be finite".to_string(),
            ));
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
            relative_overlap,
            relative_uncertainty,
            dominant_components,
            priority_score,
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

    pub fn relative_overlap(&self) -> Option<f64> {
        self.relative_overlap
    }

    pub fn relative_uncertainty(&self) -> Option<f64> {
        self.relative_uncertainty
    }

    pub fn dominant_components(&self) -> &[String] {
        &self.dominant_components
    }

    pub fn priority_score(&self) -> f64 {
        self.priority_score
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
    focus_components: Vec<String>,
    proposal_strategy: Option<ProposalStrategy>,
    priority_score: f64,
    reason: String,
}

impl ScheduleSuggestion {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        kind: SuggestionKind,
        edge_index: usize,
        from_state: StatePoint,
        to_state: StatePoint,
        proposed_state: Option<StatePoint>,
        focus_components: Vec<String>,
        proposal_strategy: Option<ProposalStrategy>,
        priority_score: f64,
        reason: String,
    ) -> Self {
        Self {
            kind,
            edge_index,
            from_state,
            to_state,
            proposed_state,
            focus_components,
            proposal_strategy,
            priority_score,
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

    pub fn focus_components(&self) -> &[String] {
        &self.focus_components
    }

    pub fn proposal_strategy(&self) -> Option<ProposalStrategy> {
        self.proposal_strategy
    }

    pub fn priority_score(&self) -> f64 {
        self.priority_score
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TiScheduleAdvisorOptions {
    pub n_blocks: usize,
    pub block_cv_min: f64,
    pub curvature_z_min: f64,
    pub slope_z_min: f64,
    pub interval_uncertainty_z_min: f64,
    pub suggest_midpoints: bool,
}

impl Default for TiScheduleAdvisorOptions {
    fn default() -> Self {
        Self {
            n_blocks: 4,
            block_cv_min: 0.15,
            curvature_z_min: 1.5,
            slope_z_min: 1.5,
            interval_uncertainty_z_min: 1.5,
            suggest_midpoints: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TiEdgeSeverity {
    Healthy,
    Monitor,
    AddSampling,
    AddWindow,
    AddWindowAndSampling,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TiSuggestionKind {
    NoChange,
    ExtendSampling,
    InsertWindow,
    InsertWindowAndExtendSampling,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TiWindowDiagnostic {
    window_index: usize,
    state: StatePoint,
    kept_samples: usize,
    mean_dhdl: f64,
    sem_dhdl: Option<f64>,
    block_mean: Option<f64>,
    block_stddev: Option<f64>,
    block_cv: Option<f64>,
    split_half_dhdl_delta: Option<f64>,
}

impl TiWindowDiagnostic {
    pub fn window_index(&self) -> usize {
        self.window_index
    }

    pub fn state(&self) -> &StatePoint {
        &self.state
    }

    pub fn kept_samples(&self) -> usize {
        self.kept_samples
    }

    pub fn mean_dhdl(&self) -> f64 {
        self.mean_dhdl
    }

    pub fn sem_dhdl(&self) -> Option<f64> {
        self.sem_dhdl
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

    pub fn split_half_dhdl_delta(&self) -> Option<f64> {
        self.split_half_dhdl_delta
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TiIntervalDiagnostic {
    interval_index: usize,
    from_state: StatePoint,
    to_state: StatePoint,
    delta_lambda: f64,
    left_mean_dhdl: f64,
    right_mean_dhdl: f64,
    trapezoid_contribution: f64,
    slope: f64,
    abs_slope: f64,
    curvature: Option<f64>,
    interval_uncertainty: Option<f64>,
    severity: TiEdgeSeverity,
    priority_score: f64,
}

impl TiIntervalDiagnostic {
    pub fn interval_index(&self) -> usize {
        self.interval_index
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }

    pub fn delta_lambda(&self) -> f64 {
        self.delta_lambda
    }

    pub fn left_mean_dhdl(&self) -> f64 {
        self.left_mean_dhdl
    }

    pub fn right_mean_dhdl(&self) -> f64 {
        self.right_mean_dhdl
    }

    pub fn trapezoid_contribution(&self) -> f64 {
        self.trapezoid_contribution
    }

    pub fn slope(&self) -> f64 {
        self.slope
    }

    pub fn abs_slope(&self) -> f64 {
        self.abs_slope
    }

    pub fn curvature(&self) -> Option<f64> {
        self.curvature
    }

    pub fn interval_uncertainty(&self) -> Option<f64> {
        self.interval_uncertainty
    }

    pub fn severity(&self) -> TiEdgeSeverity {
        self.severity
    }

    pub fn priority_score(&self) -> f64 {
        self.priority_score
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TiScheduleSuggestion {
    kind: TiSuggestionKind,
    interval_index: usize,
    from_state: StatePoint,
    to_state: StatePoint,
    proposed_state: Option<StatePoint>,
    priority_score: f64,
    reason: String,
}

impl TiScheduleSuggestion {
    pub fn kind(&self) -> TiSuggestionKind {
        self.kind
    }

    pub fn interval_index(&self) -> usize {
        self.interval_index
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

    pub fn priority_score(&self) -> f64 {
        self.priority_score
    }

    pub fn reason(&self) -> &str {
        &self.reason
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TiScheduleAdvice {
    windows: Vec<TiWindowDiagnostic>,
    intervals: Vec<TiIntervalDiagnostic>,
    suggestions: Vec<TiScheduleSuggestion>,
}

impl TiScheduleAdvice {
    pub fn new(
        windows: Vec<TiWindowDiagnostic>,
        intervals: Vec<TiIntervalDiagnostic>,
        suggestions: Vec<TiScheduleSuggestion>,
    ) -> Self {
        Self {
            windows,
            intervals,
            suggestions,
        }
    }

    pub fn windows(&self) -> &[TiWindowDiagnostic] {
        &self.windows
    }

    pub fn intervals(&self) -> &[TiIntervalDiagnostic] {
        &self.intervals
    }

    pub fn suggestions(&self) -> &[TiScheduleSuggestion] {
        &self.suggestions
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TiMethodRecommendationOptions {
    pub schedule: TiScheduleAdvisorOptions,
    pub disagreement_sigma_factor: f64,
    pub absolute_delta_tolerance: f64,
}

impl Default for TiMethodRecommendationOptions {
    fn default() -> Self {
        Self {
            schedule: TiScheduleAdvisorOptions::default(),
            disagreement_sigma_factor: 1.5,
            absolute_delta_tolerance: 1.0e-2,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TiMethodAssessment {
    method: IntegrationMethod,
    supported: bool,
    delta_f: Option<f64>,
    uncertainty: Option<f64>,
    rationale: String,
}

impl TiMethodAssessment {
    pub fn method(&self) -> IntegrationMethod {
        self.method
    }

    pub fn supported(&self) -> bool {
        self.supported
    }

    pub fn delta_f(&self) -> Option<f64> {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn rationale(&self) -> &str {
        &self.rationale
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TiMethodRecommendation {
    recommended_method: IntegrationMethod,
    reason: String,
    assessments: Vec<TiMethodAssessment>,
    cross_method_spread: Option<f64>,
}

impl TiMethodRecommendation {
    pub fn recommended_method(&self) -> IntegrationMethod {
        self.recommended_method
    }

    pub fn reason(&self) -> &str {
        &self.reason
    }

    pub fn assessments(&self) -> &[TiMethodAssessment] {
        &self.assessments
    }

    pub fn cross_method_spread(&self) -> Option<f64> {
        self.cross_method_spread
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
    MbarEstimator::new(options.unwrap_or_default())
        .fit(windows)?
        .overlap_matrix()
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
    let options = options.unwrap_or_default();
    if matches!(
        options.method,
        crate::estimators::IntegrationMethod::GaussianQuadrature
    ) {
        // Prefix subsets of an n-point Gauss-Legendre schedule are not themselves valid
        // quadrature schedules, so TI convergence does not have a sound interpretation here.
        return Err(CoreError::Unsupported(
            "TI convergence is unsupported for Gaussian quadrature because prefix windows are not valid quadrature schedules".to_string(),
        ));
    }

    let start = match options.method {
        crate::estimators::IntegrationMethod::Trapezoidal => 2,
        crate::estimators::IntegrationMethod::Simpson => 3,
        crate::estimators::IntegrationMethod::CubicSpline => 2,
        crate::estimators::IntegrationMethod::Pchip => 2,
        crate::estimators::IntegrationMethod::Akima => 2,
        crate::estimators::IntegrationMethod::GaussianQuadrature => unreachable!(),
    };
    if series.len() < start {
        return Err(CoreError::InvalidShape {
            expected: start,
            found: series.len(),
        });
    }

    let estimator = TiEstimator::new(options);
    let mut points = Vec::with_capacity(series.len() - start + 1);
    for end in start..=series.len() {
        if matches!(
            estimator.options.method,
            crate::estimators::IntegrationMethod::Simpson
        ) && (end & 1) == 0
        {
            continue;
        }
        let result = estimator.estimate(&series[..end])?;
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
        |subset| BarEstimator::new(options.clone().unwrap_or_default()).estimate(subset),
        false,
        2,
    )
}

pub fn mbar_convergence(
    windows: &[UNkMatrix],
    options: Option<MbarOptions>,
) -> Result<Vec<ConvergencePoint>> {
    if windows.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let trimmed = trim_windows_to_sampled_states(windows)?;
    let total_windows = trimmed.len();
    let base_options = options.unwrap_or_default();
    let mut points = Vec::with_capacity(total_windows);
    let mut initial_f_k = base_options.initial_f_k.clone();

    for count in 1..=total_windows {
        let subset = select_trimmed_windows_range(&trimmed[..count], 0, count)?;
        let mut fit_options = base_options.clone();
        // Each prefix adds one neighboring state, so the previous solution is a good
        // initial guess for the next MBAR solve and usually reduces iteration count.
        fit_options.initial_f_k = adapt_mbar_initial_f_k(initial_f_k.as_deref(), count);

        let fit = MbarEstimator::new(fit_options).fit(&subset)?;
        let result = fit.result_with_uncertainty()?;
        points.push(convergence_point_from_matrix(count, &result, false)?);
        initial_f_k = Some(fit.free_energies().to_vec());
    }

    Ok(points)
}

pub fn exp_convergence(
    windows: &[UNkMatrix],
    options: Option<ExpOptions>,
) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| {
            ExpEstimator::new(options.clone().unwrap_or_default()).estimate_with_uncertainty(subset)
        },
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
        |subset| {
            ExpEstimator::new(options.clone().unwrap_or_default()).estimate_with_uncertainty(subset)
        },
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
    let metrics = match options.estimator {
        AdvisorEstimator::Mbar => schedule_edge_metrics_mbar(&ordered, options)?,
        AdvisorEstimator::Bar => schedule_edge_metrics_bar(&ordered, options)?,
    };

    let overlap_mins = metrics
        .iter()
        .map(|edge| edge.overlap_forward.min(edge.overlap_reverse))
        .collect::<Vec<_>>();
    let uncertainties = metrics
        .iter()
        .map(|edge| edge.uncertainty)
        .collect::<Vec<_>>();

    let mut edges = Vec::with_capacity(metrics.len());
    let mut suggestions = Vec::new();
    for metric in metrics {
        let relative_overlap = neighbor_relative_value(&overlap_mins, metric.edge_index);
        let relative_uncertainty = neighbor_relative_option(&uncertainties, metric.edge_index);
        let priority_score = edge_priority_score(
            metric.overlap_forward.min(metric.overlap_reverse),
            metric.block_cv,
            metric.uncertainty,
            relative_overlap,
            relative_uncertainty,
            &options,
        );
        let dominant_components = dominant_component_names(
            &metric.from_state,
            &metric.to_state,
            metric.lambda_labels.as_deref(),
        );
        let severity = classify_edge(
            metric.overlap_forward.min(metric.overlap_reverse),
            metric.block_cv,
            relative_overlap,
            relative_uncertainty,
            &options,
        );
        let diagnostic = AdjacentEdgeDiagnostic::new(
            metric.edge_index,
            metric.from_state,
            metric.to_state,
            metric.lambda_labels,
            metric.overlap_forward,
            metric.overlap_reverse,
            metric.delta_f,
            metric.uncertainty,
            metric.block_mean,
            metric.block_stddev,
            metric.block_cv,
            relative_overlap,
            relative_uncertainty,
            dominant_components,
            priority_score,
            severity,
        )?;

        if let Some(suggestion) = schedule_suggestion(&diagnostic, &options)? {
            suggestions.push(suggestion);
        }
        edges.push(diagnostic);
    }

    suggestions.sort_by(|left, right| {
        right
            .priority_score()
            .partial_cmp(&left.priority_score())
            .unwrap_or(Ordering::Equal)
    });
    Ok(ScheduleAdvice::new(edges, suggestions))
}

fn schedule_edge_metrics_mbar(
    windows: &[UNkMatrix],
    options: ScheduleAdvisorOptions,
) -> Result<Vec<EdgeMetrics>> {
    let trimmed = trim_windows_to_sampled_states(windows)?;
    let overlap = overlap_matrix(&trimmed, Some(MbarOptions::default()))?;
    let n_states = overlap.n_states();
    let labels = overlap.lambda_labels().map(|values| values.to_vec());

    let mut metrics = Vec::with_capacity(n_states - 1);
    for edge_index in 0..(n_states - 1) {
        let pair_windows =
            select_trimmed_windows_range(&trimmed[edge_index..edge_index + 2], edge_index, 2)?;
        let pair_result = fit_pair(&pair_windows, options.estimator)?;
        let block_stats = pair_block_stats(&pair_windows, options)?;
        let overlap_forward = overlap.values()[edge_index * n_states + edge_index + 1];
        let overlap_reverse = overlap.values()[(edge_index + 1) * n_states + edge_index];
        metrics.push(EdgeMetrics {
            edge_index,
            from_state: pair_result.states()[0].clone(),
            to_state: pair_result.states()[1].clone(),
            lambda_labels: labels.clone(),
            overlap_forward,
            overlap_reverse,
            delta_f: pair_result.values()[1],
            uncertainty: pair_result.uncertainties().map(|values| values[1]),
            block_mean: block_stats.map(|(mean, _, _)| mean),
            block_stddev: block_stats.map(|(_, stddev, _)| stddev),
            block_cv: block_stats.map(|(_, _, cv)| cv),
        });
    }

    Ok(metrics)
}

fn schedule_edge_metrics_bar(
    windows: &[UNkMatrix],
    options: ScheduleAdvisorOptions,
) -> Result<Vec<EdgeMetrics>> {
    let mut metrics = Vec::with_capacity(windows.len() - 1);
    for edge_index in 0..(windows.len() - 1) {
        let pair_states = windows[edge_index..edge_index + 2]
            .iter()
            .map(|window| {
                window.sampled_state().cloned().ok_or_else(|| {
                    CoreError::InvalidState(
                        "sampled_state required for schedule advisor".to_string(),
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let pair_windows =
            select_windows_for_states(&windows[edge_index..edge_index + 2], &pair_states)?;
        let pair_result = fit_pair(&pair_windows, options.estimator)?;
        let pair_overlap = overlap_matrix(&pair_windows, Some(MbarOptions::default()))?;
        let block_stats = pair_block_stats(&pair_windows, options)?;
        metrics.push(EdgeMetrics {
            edge_index,
            from_state: pair_result.states()[0].clone(),
            to_state: pair_result.states()[1].clone(),
            lambda_labels: pair_result.lambda_labels().map(|values| values.to_vec()),
            overlap_forward: pair_overlap.values()[1],
            overlap_reverse: pair_overlap.values()[2],
            delta_f: pair_result.values()[1],
            uncertainty: pair_result.uncertainties().map(|values| values[1]),
            block_mean: block_stats.map(|(mean, _, _)| mean),
            block_stddev: block_stats.map(|(_, stddev, _)| stddev),
            block_cv: block_stats.map(|(_, _, cv)| cv),
        });
    }

    Ok(metrics)
}

pub fn advise_ti_schedule(
    series: &[DhdlSeries],
    options: Option<TiScheduleAdvisorOptions>,
) -> Result<TiScheduleAdvice> {
    if series.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: series.len(),
        });
    }

    let options = validate_ti_schedule_advisor_options(options.unwrap_or_default())?;
    let ordered = sort_dhdl_series_by_state(series)?;

    let mut window_means = Vec::with_capacity(ordered.len());
    let mut window_sems = Vec::with_capacity(ordered.len());
    let mut block_cvs = Vec::with_capacity(ordered.len());
    let mut split_half_dhdl_deltas = Vec::with_capacity(ordered.len());
    let mut windows = Vec::with_capacity(ordered.len());

    for (window_index, item) in ordered.iter().enumerate() {
        let mean = mean_dhdl_values(item.values())?;
        let sem = sem_dhdl_values(item.values());
        let block_stats = dhdl_block_stats(item, options.n_blocks)?;
        let split_half_dhdl_delta = split_half_dhdl_delta(item.values());
        windows.push(TiWindowDiagnostic {
            window_index,
            state: item.state().clone(),
            kept_samples: item.values().len(),
            mean_dhdl: mean,
            sem_dhdl: sem,
            block_mean: block_stats.map(|(block_mean, _, _)| block_mean),
            block_stddev: block_stats.map(|(_, stddev, _)| stddev),
            block_cv: block_stats.map(|(_, _, cv)| cv),
            split_half_dhdl_delta,
        });
        window_means.push(mean);
        window_sems.push(sem);
        block_cvs.push(block_stats.map(|(_, _, cv)| cv));
        split_half_dhdl_deltas.push(split_half_dhdl_delta);
    }

    let lambdas = ordered
        .iter()
        .map(|item| extract_scalar_lambda(item.state(), "TI schedule advisor"))
        .collect::<Result<Vec<_>>>()?;

    let slopes = (0..(ordered.len() - 1))
        .map(|index| {
            let delta_lambda = lambdas[index + 1] - lambdas[index];
            if delta_lambda <= 0.0 {
                return Err(CoreError::InvalidState(
                    "TI schedule advisor requires strictly increasing lambda values".to_string(),
                ));
            }
            Ok((window_means[index + 1] - window_means[index]) / delta_lambda)
        })
        .collect::<Result<Vec<_>>>()?;

    let slope_magnitudes = slopes
        .iter()
        .map(|value: &f64| value.abs())
        .collect::<Vec<_>>();
    let slope_mean = mean_optionless(&slope_magnitudes);
    let slope_std = stddev_optionless(&slope_magnitudes, slope_mean);

    let point_curvatures = ti_point_curvatures(&lambdas, &window_means)?;
    let curvature_abs = point_curvatures
        .iter()
        .flatten()
        .map(|value| value.abs())
        .collect::<Vec<_>>();
    let curvature_mean = mean_optionless(&curvature_abs);
    let curvature_std = stddev_optionless(&curvature_abs, curvature_mean);

    let interval_uncertainties = (0..(ordered.len() - 1))
        .map(|index| {
            ti_interval_uncertainty(
                lambdas[index + 1] - lambdas[index],
                window_sems[index],
                window_sems[index + 1],
            )
        })
        .collect::<Vec<_>>();
    let interval_uncertainty_values = interval_uncertainties
        .iter()
        .flatten()
        .copied()
        .collect::<Vec<_>>();
    let uncertainty_mean = mean_optionless(&interval_uncertainty_values);
    let uncertainty_std = stddev_optionless(&interval_uncertainty_values, uncertainty_mean);

    let mut intervals = Vec::with_capacity(ordered.len() - 1);
    let mut suggestions = Vec::new();
    for interval_index in 0..(ordered.len() - 1) {
        let delta_lambda = lambdas[interval_index + 1] - lambdas[interval_index];
        let slope = slopes[interval_index];
        let abs_slope = slope.abs();
        let curvature = match (
            point_curvatures[interval_index],
            point_curvatures[interval_index + 1],
        ) {
            (Some(left), Some(right)) => Some(left.abs().max(right.abs())),
            (Some(value), None) | (None, Some(value)) => Some(value.abs()),
            (None, None) => None,
        };
        let slope_z = z_score(abs_slope, slope_mean, slope_std);
        let curvature_z = curvature.and_then(|value| z_score(value, curvature_mean, curvature_std));
        let interval_uncertainty = interval_uncertainties[interval_index];
        let uncertainty_z = interval_uncertainty
            .and_then(|value| z_score(value, uncertainty_mean, uncertainty_std));
        let left_block_cv = block_cvs[interval_index];
        let right_block_cv = block_cvs[interval_index + 1];
        let left_split_half_drift = split_half_dhdl_deltas[interval_index];
        let right_split_half_drift = split_half_dhdl_deltas[interval_index + 1];
        let sampling_issue = left_block_cv.is_some_and(|value| value >= options.block_cv_min)
            || right_block_cv.is_some_and(|value| value >= options.block_cv_min)
            || uncertainty_z.is_some_and(|value| value >= options.interval_uncertainty_z_min)
            || left_split_half_drift
                .is_some_and(|value| value > window_sems[interval_index].unwrap_or(0.0) * 2.0)
            || right_split_half_drift
                .is_some_and(|value| value > window_sems[interval_index + 1].unwrap_or(0.0) * 2.0);
        let structure_issue = slope_z.is_some_and(|value| value >= options.slope_z_min)
            || curvature_z.is_some_and(|value| value >= options.curvature_z_min);
        let severity = match (structure_issue, sampling_issue) {
            (true, true) => TiEdgeSeverity::AddWindowAndSampling,
            (true, false) => TiEdgeSeverity::AddWindow,
            (false, true) => TiEdgeSeverity::AddSampling,
            (false, false)
                if slope_z.is_some_and(|value| value >= options.slope_z_min * 0.75)
                    || curvature_z.is_some_and(|value| value >= options.curvature_z_min * 0.75)
                    || uncertainty_z.is_some_and(|value| {
                        value >= options.interval_uncertainty_z_min * 0.75
                    }) =>
            {
                TiEdgeSeverity::Monitor
            }
            _ => TiEdgeSeverity::Healthy,
        };
        let priority_score = ti_priority_score(
            slope_z,
            curvature_z,
            uncertainty_z,
            left_block_cv,
            right_block_cv,
            left_split_half_drift,
            right_split_half_drift,
            window_sems[interval_index],
            window_sems[interval_index + 1],
            options.block_cv_min,
        );
        let diagnostic = TiIntervalDiagnostic {
            interval_index,
            from_state: ordered[interval_index].state().clone(),
            to_state: ordered[interval_index + 1].state().clone(),
            delta_lambda,
            left_mean_dhdl: window_means[interval_index],
            right_mean_dhdl: window_means[interval_index + 1],
            trapezoid_contribution: 0.5
                * (window_means[interval_index] + window_means[interval_index + 1])
                * delta_lambda,
            slope,
            abs_slope,
            curvature,
            interval_uncertainty,
            severity,
            priority_score,
        };
        if let Some(suggestion) = ti_schedule_suggestion(&diagnostic, &options)? {
            suggestions.push(suggestion);
        }
        intervals.push(diagnostic);
    }

    suggestions.sort_by(|left, right| {
        right
            .priority_score()
            .partial_cmp(&left.priority_score())
            .unwrap_or(Ordering::Equal)
    });
    Ok(TiScheduleAdvice::new(windows, intervals, suggestions))
}

pub fn recommend_ti_method(
    series: &[DhdlSeries],
    options: Option<TiMethodRecommendationOptions>,
) -> Result<TiMethodRecommendation> {
    let options = validate_ti_method_recommendation_options(options.unwrap_or_default())?;
    let advice = advise_ti_schedule(series, Some(options.schedule))?;

    let lambdas = advice
        .windows()
        .iter()
        .map(|window| extract_scalar_lambda(window.state(), "TI method recommendation"))
        .collect::<Result<Vec<_>>>()?;
    let means = advice
        .windows()
        .iter()
        .map(|window| window.mean_dhdl())
        .collect::<Vec<_>>();

    let uniform_spacing = has_uniform_spacing(&lambdas);
    let odd_window_count = (lambdas.len() & 1) == 1;
    let monotone_means = is_monotone_series(&means);
    let has_structure_issue = advice.suggestions().iter().any(|suggestion| {
        matches!(
            suggestion.kind(),
            TiSuggestionKind::InsertWindow | TiSuggestionKind::InsertWindowAndExtendSampling
        )
    });
    let has_sampling_issue = advice.suggestions().iter().any(|suggestion| {
        matches!(
            suggestion.kind(),
            TiSuggestionKind::ExtendSampling | TiSuggestionKind::InsertWindowAndExtendSampling
        )
    });

    let mut assessments = Vec::with_capacity(TI_METHOD_CANDIDATES.len());
    for method in TI_METHOD_CANDIDATES {
        let estimator = TiEstimator::new(TiOptions {
            method,
            parallel: false,
        });
        match estimator.fit(series) {
            Ok(fit) => assessments.push(TiMethodAssessment {
                method,
                supported: true,
                delta_f: Some(fit.delta_f()),
                uncertainty: fit.uncertainty(),
                rationale: supported_ti_method_rationale(
                    method,
                    uniform_spacing,
                    odd_window_count,
                    monotone_means,
                ),
            }),
            Err(CoreError::Unsupported(message)) => assessments.push(TiMethodAssessment {
                method,
                supported: false,
                delta_f: None,
                uncertainty: None,
                rationale: message,
            }),
            Err(error @ CoreError::InvalidShape { .. })
                if method != IntegrationMethod::Trapezoidal =>
            {
                assessments.push(TiMethodAssessment {
                    method,
                    supported: false,
                    delta_f: None,
                    uncertainty: None,
                    rationale: error.to_string(),
                })
            }
            Err(error) if method == IntegrationMethod::Trapezoidal => return Err(error),
            Err(error) => assessments.push(TiMethodAssessment {
                method,
                supported: false,
                delta_f: None,
                uncertainty: None,
                rationale: error.to_string(),
            }),
        }
    }

    let supported_deltas = assessments
        .iter()
        .filter_map(|assessment| assessment.delta_f())
        .collect::<Vec<_>>();
    let cross_method_spread = if supported_deltas.is_empty() {
        None
    } else {
        let min = supported_deltas
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let max = supported_deltas
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        Some(max - min)
    };
    let disagreement_limit = options.absolute_delta_tolerance.max(
        options.disagreement_sigma_factor
            * preferred_uncertainty_reference(&assessments).unwrap_or(0.0),
    );
    let methods_disagree = cross_method_spread.is_some_and(|spread| spread > disagreement_limit);
    let severe_methods_disagree =
        cross_method_spread.is_some_and(|spread| spread > disagreement_limit * 3.0);

    let (recommended_method, reason) = choose_ti_method(
        &assessments,
        uniform_spacing,
        odd_window_count,
        monotone_means,
        has_structure_issue,
        has_sampling_issue,
        methods_disagree,
        severe_methods_disagree,
    );

    Ok(TiMethodRecommendation {
        recommended_method,
        reason,
        assessments,
        cross_method_spread,
    })
}

const TI_METHOD_CANDIDATES: [IntegrationMethod; 6] = [
    IntegrationMethod::Trapezoidal,
    IntegrationMethod::Simpson,
    IntegrationMethod::CubicSpline,
    IntegrationMethod::Pchip,
    IntegrationMethod::Akima,
    IntegrationMethod::GaussianQuadrature,
];

fn choose_ti_method(
    assessments: &[TiMethodAssessment],
    uniform_spacing: bool,
    odd_window_count: bool,
    monotone_means: bool,
    has_structure_issue: bool,
    has_sampling_issue: bool,
    methods_disagree: bool,
    severe_methods_disagree: bool,
) -> (IntegrationMethod, String) {
    if supports_ti_method(assessments, IntegrationMethod::GaussianQuadrature) {
        return (
            IntegrationMethod::GaussianQuadrature,
            "sampled lambdas match a supported Gauss-Legendre rule on [0, 1]".to_string(),
        );
    }

    if has_structure_issue {
        return (
            IntegrationMethod::Trapezoidal,
            "TI schedule diagnostics flagged under-resolved intervals, so trapezoidal is the safest current choice".to_string(),
        );
    }

    if uniform_spacing
        && odd_window_count
        && !has_sampling_issue
        && !methods_disagree
        && supports_ti_method(assessments, IntegrationMethod::Simpson)
    {
        return (
            IntegrationMethod::Simpson,
            "uniform odd-count lambda grid with healthy TI diagnostics and low cross-method disagreement".to_string(),
        );
    }

    if methods_disagree {
        if severe_methods_disagree {
            return (
                IntegrationMethod::Trapezoidal,
                "supported methods disagree strongly enough that trapezoidal is the conservative fallback".to_string(),
            );
        }
        if !has_sampling_issue
            && monotone_means
            && supports_ti_method(assessments, IntegrationMethod::Pchip)
        {
            return (
                IntegrationMethod::Pchip,
                "supported methods disagree materially, and monotone dH/dlambda means favor shape-preserving PCHIP".to_string(),
            );
        }
        if !has_sampling_issue && supports_ti_method(assessments, IntegrationMethod::Akima) {
            return (
                IntegrationMethod::Akima,
                "supported methods disagree materially on a nonuniform schedule, so Akima is the most local smooth interpolant available".to_string(),
            );
        }
        return (
            IntegrationMethod::Trapezoidal,
            "supported methods disagree materially and/or the TI diagnostics look unstable, so trapezoidal is the conservative fallback".to_string(),
        );
    }

    if !has_sampling_issue
        && monotone_means
        && supports_ti_method(assessments, IntegrationMethod::Pchip)
    {
        return (
            IntegrationMethod::Pchip,
            "healthy nonuniform schedule with monotone dH/dlambda means favors shape-preserving PCHIP".to_string(),
        );
    }

    if !has_sampling_issue && supports_ti_method(assessments, IntegrationMethod::CubicSpline) {
        return (
            IntegrationMethod::CubicSpline,
            "healthy nonuniform schedule with low cross-method disagreement favors a smooth cubic-spline integral".to_string(),
        );
    }

    if !has_sampling_issue && supports_ti_method(assessments, IntegrationMethod::Akima) {
        return (
            IntegrationMethod::Akima,
            "healthy nonuniform schedule can use Akima as a local smooth interpolant".to_string(),
        );
    }

    (
        IntegrationMethod::Trapezoidal,
        "trapezoidal remains valid on arbitrary TI schedules and is the safest fallback here"
            .to_string(),
    )
}

fn supports_ti_method(assessments: &[TiMethodAssessment], method: IntegrationMethod) -> bool {
    assessments
        .iter()
        .find(|assessment| assessment.method() == method)
        .is_some_and(|assessment| assessment.supported())
}

fn supported_ti_method_rationale(
    method: IntegrationMethod,
    uniform_spacing: bool,
    odd_window_count: bool,
    monotone_means: bool,
) -> String {
    match method {
        IntegrationMethod::Trapezoidal => {
            "supported on arbitrary strictly increasing lambda schedules".to_string()
        }
        IntegrationMethod::Simpson => {
            if uniform_spacing && odd_window_count {
                "supported on this uniform odd-count lambda grid".to_string()
            } else {
                "supported when the lambda grid is uniform and has an odd number of windows"
                    .to_string()
            }
        }
        IntegrationMethod::CubicSpline => {
            "supported as a smooth global cubic interpolant with propagated uncertainty".to_string()
        }
        IntegrationMethod::Pchip => {
            if monotone_means {
                "supported and shape-preserving for this monotone dH/dlambda trend".to_string()
            } else {
                "supported as a shape-preserving cubic Hermite interpolant".to_string()
            }
        }
        IntegrationMethod::Akima => {
            "supported as a local cubic interpolant for nonuniform schedules".to_string()
        }
        IntegrationMethod::GaussianQuadrature => {
            "supported because the sampled lambdas match a Gauss-Legendre quadrature rule"
                .to_string()
        }
    }
}

fn preferred_uncertainty_reference(assessments: &[TiMethodAssessment]) -> Option<f64> {
    assessments
        .iter()
        .find(|assessment| assessment.method() == IntegrationMethod::Trapezoidal)
        .and_then(|assessment| assessment.uncertainty())
        .or_else(|| {
            assessments
                .iter()
                .filter_map(|assessment| assessment.uncertainty())
                .find(|value| value.is_finite())
        })
}

fn has_uniform_spacing(lambdas: &[f64]) -> bool {
    if lambdas.len() < 3 {
        return false;
    }
    let h = (lambdas[lambdas.len() - 1] - lambdas[0]) / ((lambdas.len() - 1) as f64);
    if h == 0.0 {
        return false;
    }
    let tol = h.abs() * 1.0e-8;
    for (index, lambda) in lambdas.iter().enumerate().skip(1) {
        let expected = lambdas[0] + (index as f64) * h;
        if (lambda - expected).abs() > tol {
            return false;
        }
    }
    true
}

fn is_monotone_series(values: &[f64]) -> bool {
    values.windows(2).all(|pair| pair[1] >= pair[0])
        || values.windows(2).all(|pair| pair[1] <= pair[0])
}

pub(crate) fn ti_block_average(
    series: &[DhdlSeries],
    n_blocks: usize,
    options: Option<TiOptions>,
) -> Result<Vec<BlockEstimate>> {
    let options = options.unwrap_or_default();
    let minimum_windows = match options.method {
        crate::estimators::IntegrationMethod::Trapezoidal => 2,
        crate::estimators::IntegrationMethod::Simpson => 3,
        crate::estimators::IntegrationMethod::CubicSpline => 2,
        crate::estimators::IntegrationMethod::Pchip => 2,
        crate::estimators::IntegrationMethod::Akima => 2,
        crate::estimators::IntegrationMethod::GaussianQuadrature => 1,
    };
    if series.len() < minimum_windows {
        return Err(CoreError::InvalidShape {
            expected: minimum_windows,
            found: series.len(),
        });
    }
    if n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let estimator = TiEstimator::new(options);
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
        let result = estimator.estimate(&block)?;
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
        |subset| {
            MbarEstimator::new(options.clone().unwrap_or_default())
                .estimate_with_uncertainty(subset)
        },
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
        |subset| BarEstimator::new(options.clone().unwrap_or_default()).estimate(subset),
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
        |subset| {
            ExpEstimator::new(options.clone().unwrap_or_default()).estimate_with_uncertainty(subset)
        },
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
        |subset| {
            ExpEstimator::new(options.clone().unwrap_or_default()).estimate_with_uncertainty(subset)
        },
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
    let trimmed = trim_windows_to_sampled_states(windows)?;
    let total_windows = trimmed.len();
    let mut points = Vec::with_capacity(windows.len() - minimum_windows + 1);
    for count in minimum_windows..=total_windows {
        let start = if reverse { total_windows - count } else { 0 };
        let subset = select_trimmed_windows_range(&trimmed[start..start + count], start, count)?;
        let result = fit(&subset)?;
        points.push(convergence_point_from_matrix(count, &result, reverse)?);
    }
    Ok(points)
}

fn adapt_mbar_initial_f_k(source: Option<&[f64]>, n_states: usize) -> Option<Vec<f64>> {
    let source = source?;
    if source.is_empty() {
        return Some(vec![0.0; n_states]);
    }
    let mut initial = source.to_vec();
    if initial.len() > n_states {
        initial.truncate(n_states);
    } else if initial.len() < n_states {
        let next = match initial.as_slice() {
            [.., prev, last] => last + (last - prev),
            [last] => *last,
            [] => 0.0,
        };
        initial.resize(n_states, next);
    }
    Some(initial)
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

    let trimmed = trim_windows_to_sampled_states(windows)?;
    let blocked = trimmed
        .iter()
        .map(|window| split_u_nk_window(window, n_blocks))
        .collect::<Result<Vec<_>>>()?;

    let mut points = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let subset = blocked
            .iter()
            .map(|chunks| chunks[block_index].clone())
            .collect::<Vec<_>>();
        let result = fit(&subset)?;
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
    select_windows_for_states(windows, &states)
}

fn select_windows_for_states(
    windows: &[UNkMatrix],
    states: &[StatePoint],
) -> Result<Vec<UNkMatrix>> {
    let mut selected = Vec::with_capacity(windows.len());
    for window in windows {
        let indices = states
            .iter()
            .map(|state| find_state_index_exact(window.evaluated_states(), state))
            .collect::<Result<Vec<_>>>()?;
        if indices.len() == window.n_states() && indices.iter().copied().eq(0..window.n_states()) {
            selected.push(window.clone());
        } else {
            selected.push(select_window_state_indices(window, &indices, states)?);
        }
    }

    Ok(selected)
}

fn select_trimmed_windows_range(
    windows: &[UNkMatrix],
    state_start: usize,
    state_count: usize,
) -> Result<Vec<UNkMatrix>> {
    if windows.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    if state_count == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let states = windows[0]
        .evaluated_states()
        .get(state_start..state_start + state_count)
        .ok_or(CoreError::InvalidShape {
            expected: state_start + state_count,
            found: windows[0].n_states(),
        })?
        .to_vec();

    windows
        .iter()
        .map(|window| select_window_state_range(window, state_start, state_count, &states))
        .collect()
}

fn select_window_state_range(
    window: &UNkMatrix,
    state_start: usize,
    state_count: usize,
    states: &[StatePoint],
) -> Result<UNkMatrix> {
    let mut data = Vec::with_capacity(window.n_samples() * state_count);
    for sample_idx in 0..window.n_samples() {
        let row_offset = sample_idx * window.n_states();
        let range_start = row_offset + state_start;
        let range_end = range_start + state_count;
        data.extend_from_slice(&window.data()[range_start..range_end]);
    }

    UNkMatrix::new_with_labels(
        window.n_samples(),
        state_count,
        data,
        window.time_ps().to_vec(),
        window.sampled_state().cloned(),
        states.to_vec(),
        window.lambda_labels().map(|labels| labels.to_vec()),
    )
}

fn select_window_state_indices(
    window: &UNkMatrix,
    indices: &[usize],
    states: &[StatePoint],
) -> Result<UNkMatrix> {
    let mut data = Vec::with_capacity(window.n_samples() * indices.len());
    if let Some((state_start, state_count)) = contiguous_state_range(indices) {
        for sample_idx in 0..window.n_samples() {
            let row_offset = sample_idx * window.n_states();
            let range_start = row_offset + state_start;
            let range_end = range_start + state_count;
            data.extend_from_slice(&window.data()[range_start..range_end]);
        }
    } else {
        for sample_idx in 0..window.n_samples() {
            let row_offset = sample_idx * window.n_states();
            for &idx in indices {
                data.push(window.data()[row_offset + idx]);
            }
        }
    }

    UNkMatrix::new_with_labels(
        window.n_samples(),
        indices.len(),
        data,
        window.time_ps().to_vec(),
        window.sampled_state().cloned(),
        states.to_vec(),
        window.lambda_labels().map(|labels| labels.to_vec()),
    )
}

fn contiguous_state_range(indices: &[usize]) -> Option<(usize, usize)> {
    let start = *indices.first()?;
    if indices
        .iter()
        .copied()
        .enumerate()
        .all(|(offset, idx)| idx == start + offset)
    {
        Some((start, indices.len()))
    } else {
        None
    }
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

fn validate_ti_schedule_advisor_options(
    options: TiScheduleAdvisorOptions,
) -> Result<TiScheduleAdvisorOptions> {
    if options.n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    for (name, value) in [
        ("block_cv_min", options.block_cv_min),
        ("curvature_z_min", options.curvature_z_min),
        ("slope_z_min", options.slope_z_min),
        (
            "interval_uncertainty_z_min",
            options.interval_uncertainty_z_min,
        ),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(CoreError::InvalidState(format!(
                "TI schedule advisor {name} must be finite and non-negative"
            )));
        }
    }
    Ok(options)
}

fn validate_ti_method_recommendation_options(
    options: TiMethodRecommendationOptions,
) -> Result<TiMethodRecommendationOptions> {
    let schedule = validate_ti_schedule_advisor_options(options.schedule)?;
    for (name, value) in [
        (
            "disagreement_sigma_factor",
            options.disagreement_sigma_factor,
        ),
        ("absolute_delta_tolerance", options.absolute_delta_tolerance),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(CoreError::InvalidState(format!(
                "TI method recommendation {name} must be finite and non-negative"
            )));
        }
    }
    Ok(TiMethodRecommendationOptions {
        schedule,
        ..options
    })
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
        if state_points_match_exact(
            ordered[idx - 1]
                .sampled_state()
                .expect("sampled_state checked above"),
            ordered[idx]
                .sampled_state()
                .expect("sampled_state checked above"),
        ) {
            return Err(CoreError::InvalidState(
                "multiple windows for same sampled_state".to_string(),
            ));
        }
    }
    Ok(ordered)
}

fn sort_dhdl_series_by_state(series: &[DhdlSeries]) -> Result<Vec<DhdlSeries>> {
    let mut ordered = series.to_vec();
    ordered.sort_by(|left, right| compare_state_points(left.state(), right.state()));
    for idx in 1..ordered.len() {
        if state_points_match_exact(ordered[idx - 1].state(), ordered[idx].state()) {
            return Err(CoreError::InvalidState(
                "multiple DhdlSeries for same lambda state".to_string(),
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
    left.temperature_k()
        .partial_cmp(&right.temperature_k())
        .unwrap_or(Ordering::Equal)
}

fn extract_scalar_lambda(state: &StatePoint, operation: &'static str) -> Result<f64> {
    if state.lambdas().len() != 1 {
        return Err(CoreError::RequiresOneDimensionalLambda { operation });
    }
    Ok(state.lambdas()[0])
}

fn fit_pair(windows: &[UNkMatrix], estimator: AdvisorEstimator) -> Result<DeltaFMatrix> {
    match estimator {
        AdvisorEstimator::Mbar => MbarEstimator::default().estimate_with_uncertainty(windows),
        AdvisorEstimator::Bar => BarEstimator::default().estimate(windows),
    }
}

fn mean_dhdl_values(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    Ok(values.iter().sum::<f64>() / values.len() as f64)
}

fn sem_dhdl_values(values: &[f64]) -> Option<f64> {
    if values.len() < 2 {
        return None;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / (values.len() - 1) as f64;
    Some((variance / values.len() as f64).sqrt())
}

fn dhdl_block_stats(series: &DhdlSeries, n_blocks: usize) -> Result<Option<(f64, f64, f64)>> {
    let blocks = match split_dhdl_series(series, n_blocks) {
        Ok(blocks) => blocks,
        Err(CoreError::InvalidShape { .. }) => return Ok(None),
        Err(error) => return Err(error),
    };
    if blocks.is_empty() {
        return Ok(None);
    }
    let values = blocks
        .iter()
        .map(|block| mean_dhdl_values(block.values()))
        .collect::<Result<Vec<_>>>()?;
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

fn split_half_dhdl_delta(values: &[f64]) -> Option<f64> {
    if values.len() < 4 {
        return None;
    }
    let half = values.len() / 2;
    if half == 0 || half == values.len() {
        return None;
    }
    let forward = values[..half].iter().sum::<f64>() / half as f64;
    let reverse = values[values.len() - half..].iter().sum::<f64>() / half as f64;
    Some((forward - reverse).abs())
}

fn ti_point_curvatures(lambdas: &[f64], means: &[f64]) -> Result<Vec<Option<f64>>> {
    if lambdas.len() != means.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: means.len(),
        });
    }
    let mut out = vec![None; lambdas.len()];
    if lambdas.len() < 3 {
        return Ok(out);
    }
    for index in 1..(lambdas.len() - 1) {
        let dx_left = lambdas[index] - lambdas[index - 1];
        let dx_right = lambdas[index + 1] - lambdas[index];
        if dx_left <= 0.0 || dx_right <= 0.0 {
            return Err(CoreError::InvalidState(
                "TI schedule advisor requires strictly increasing lambda values".to_string(),
            ));
        }
        let slope_left = (means[index] - means[index - 1]) / dx_left;
        let slope_right = (means[index + 1] - means[index]) / dx_right;
        let center_dx = 0.5 * (dx_left + dx_right);
        out[index] = Some((slope_right - slope_left) / center_dx);
    }
    Ok(out)
}

fn ti_interval_uncertainty(
    delta_lambda: f64,
    left_sem: Option<f64>,
    right_sem: Option<f64>,
) -> Option<f64> {
    match (left_sem, right_sem) {
        (Some(left), Some(right)) => Some(0.5 * delta_lambda * left.hypot(right)),
        _ => None,
    }
}

fn mean_optionless(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        Some(values.iter().sum::<f64>() / values.len() as f64)
    }
}

fn stddev_optionless(values: &[f64], mean: Option<f64>) -> Option<f64> {
    let mean = mean?;
    if values.len() < 2 {
        return None;
    }
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / values.len() as f64;
    let stddev = variance.sqrt();
    if stddev <= 1.0e-12 {
        None
    } else {
        Some(stddev)
    }
}

fn z_score(value: f64, mean: Option<f64>, stddev: Option<f64>) -> Option<f64> {
    let mean = mean?;
    let stddev = stddev?;
    Some(((value - mean) / stddev).max(0.0))
}

fn ti_priority_score(
    slope_z: Option<f64>,
    curvature_z: Option<f64>,
    uncertainty_z: Option<f64>,
    left_block_cv: Option<f64>,
    right_block_cv: Option<f64>,
    left_split_half_drift: Option<f64>,
    right_split_half_drift: Option<f64>,
    left_sem: Option<f64>,
    right_sem: Option<f64>,
    block_cv_min: f64,
) -> f64 {
    let mut score = 0.0;
    score += slope_z.unwrap_or(0.0);
    score += curvature_z.unwrap_or(0.0) * 2.0;
    score += uncertainty_z.unwrap_or(0.0) * 1.5;
    for value in [left_block_cv, right_block_cv].into_iter().flatten() {
        score += (value / block_cv_min.max(1.0e-12) - 1.0).max(0.0) * 1.5;
    }
    for (delta, sem) in [
        (left_split_half_drift, left_sem),
        (right_split_half_drift, right_sem),
    ] {
        if let (Some(delta), Some(sem)) = (delta, sem) {
            let threshold = (2.0 * sem).max(1.0e-12);
            score += (delta / threshold - 1.0).max(0.0) * 1.5;
        }
    }
    score.max(0.0)
}

fn ti_schedule_suggestion(
    diagnostic: &TiIntervalDiagnostic,
    options: &TiScheduleAdvisorOptions,
) -> Result<Option<TiScheduleSuggestion>> {
    let suggestion = match diagnostic.severity() {
        TiEdgeSeverity::AddWindow => Some(TiScheduleSuggestion {
            kind: TiSuggestionKind::InsertWindow,
            interval_index: diagnostic.interval_index(),
            from_state: diagnostic.from_state().clone(),
            to_state: diagnostic.to_state().clone(),
            proposed_state: if options.suggest_midpoints {
                Some(midpoint_state(
                    diagnostic.from_state(),
                    diagnostic.to_state(),
                )?)
            } else {
                None
            },
            priority_score: diagnostic.priority_score(),
            reason: format!(
                "interval slope/curvature indicates under-resolved TI integrand between {} and {}",
                diagnostic.from_state().lambdas()[0],
                diagnostic.to_state().lambdas()[0]
            ),
        }),
        TiEdgeSeverity::AddSampling => Some(TiScheduleSuggestion {
            kind: TiSuggestionKind::ExtendSampling,
            interval_index: diagnostic.interval_index(),
            from_state: diagnostic.from_state().clone(),
            to_state: diagnostic.to_state().clone(),
            proposed_state: None,
            priority_score: diagnostic.priority_score(),
            reason: "interval uncertainty or block instability suggests more TI sampling"
                .to_string(),
        }),
        TiEdgeSeverity::AddWindowAndSampling => Some(TiScheduleSuggestion {
            kind: TiSuggestionKind::InsertWindowAndExtendSampling,
            interval_index: diagnostic.interval_index(),
            from_state: diagnostic.from_state().clone(),
            to_state: diagnostic.to_state().clone(),
            proposed_state: if options.suggest_midpoints {
                Some(midpoint_state(
                    diagnostic.from_state(),
                    diagnostic.to_state(),
                )?)
            } else {
                None
            },
            priority_score: diagnostic.priority_score(),
            reason: "TI interval is both structurally under-resolved and statistically unstable"
                .to_string(),
        }),
        TiEdgeSeverity::Healthy | TiEdgeSeverity::Monitor => None,
    };
    Ok(suggestion)
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

fn neighbor_relative_value(values: &[f64], index: usize) -> Option<f64> {
    let neighbors = neighbor_values(values, index);
    if neighbors.is_empty() {
        return None;
    }
    let mean = neighbors.iter().sum::<f64>() / neighbors.len() as f64;
    if mean.abs() <= 1.0e-12 {
        None
    } else {
        Some(values[index] / mean)
    }
}

fn neighbor_relative_option(values: &[Option<f64>], index: usize) -> Option<f64> {
    let current = values[index]?;
    let neighbors = neighbor_optional_values(values, index);
    if neighbors.is_empty() {
        return None;
    }
    let mean = neighbors.iter().sum::<f64>() / neighbors.len() as f64;
    if mean.abs() <= 1.0e-12 {
        None
    } else {
        Some(current / mean)
    }
}

fn neighbor_values(values: &[f64], index: usize) -> Vec<f64> {
    let mut neighbors = Vec::with_capacity(2);
    if index > 0 {
        neighbors.push(values[index - 1]);
    }
    if index + 1 < values.len() {
        neighbors.push(values[index + 1]);
    }
    neighbors
}

fn neighbor_optional_values(values: &[Option<f64>], index: usize) -> Vec<f64> {
    let mut neighbors = Vec::with_capacity(2);
    if index > 0 {
        if let Some(value) = values[index - 1] {
            neighbors.push(value);
        }
    }
    if index + 1 < values.len() {
        if let Some(value) = values[index + 1] {
            neighbors.push(value);
        }
    }
    neighbors
}

fn dominant_component_names(
    from_state: &StatePoint,
    to_state: &StatePoint,
    lambda_labels: Option<&[String]>,
) -> Vec<String> {
    dominant_component_indices(from_state, to_state)
        .into_iter()
        .map(|index| {
            lambda_labels
                .and_then(|labels| labels.get(index).cloned())
                .unwrap_or_else(|| format!("lambda[{index}]"))
        })
        .collect()
}

fn dominant_component_indices(from_state: &StatePoint, to_state: &StatePoint) -> Vec<usize> {
    let deltas = from_state
        .lambdas()
        .iter()
        .zip(to_state.lambdas().iter())
        .map(|(from, to)| (to - from).abs())
        .collect::<Vec<_>>();
    let max_delta = deltas.iter().copied().fold(0.0, f64::max);
    if max_delta <= 1.0e-12 {
        return Vec::new();
    }

    deltas
        .iter()
        .enumerate()
        .filter_map(|(index, delta)| {
            if *delta >= max_delta * 0.8 {
                Some(index)
            } else {
                None
            }
        })
        .collect()
}

fn edge_priority_score(
    overlap_min: f64,
    block_cv: Option<f64>,
    uncertainty: Option<f64>,
    relative_overlap: Option<f64>,
    relative_uncertainty: Option<f64>,
    options: &ScheduleAdvisorOptions,
) -> f64 {
    let mut score = 0.0;
    if options.overlap_min > 0.0 {
        score += ((options.overlap_min - overlap_min) / options.overlap_min).max(0.0) * 3.0;
    }
    if let Some(value) = relative_overlap {
        score += (1.0 - value).max(0.0) * 2.0;
    }
    if let Some(value) = block_cv {
        score += (value / options.block_cv_min.max(1.0e-12) - 1.0).max(0.0);
    }
    if let Some(value) = relative_uncertainty {
        score += (value - 1.0).max(0.0) * 0.5;
    } else if uncertainty.is_some() {
        score += 0.1;
    }
    score.max(0.0)
}

fn classify_edge(
    overlap_min: f64,
    block_cv: Option<f64>,
    relative_overlap: Option<f64>,
    relative_uncertainty: Option<f64>,
    options: &ScheduleAdvisorOptions,
) -> EdgeSeverity {
    let localized_overlap_gap = relative_overlap.is_some_and(|value| value < 0.85);
    if overlap_min < options.overlap_min && (localized_overlap_gap || relative_overlap.is_none()) {
        EdgeSeverity::AddWindow
    } else if block_cv.is_some_and(|value| value >= options.block_cv_min)
        || relative_uncertainty.is_some_and(|value| value >= 1.5)
    {
        EdgeSeverity::AddSampling
    } else if (options.overlap_min > 0.0 && overlap_min < options.overlap_min * 1.5)
        || relative_overlap.is_some_and(|value| value < 0.95)
        || block_cv.is_some_and(|value| value >= options.block_cv_min * 0.75)
        || relative_uncertainty.is_some_and(|value| value >= 1.2)
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
            let (proposal_state, proposal_strategy) = if options.suggest_midpoints {
                let (state, strategy) = proposed_insert_state(diagnostic)?;
                (Some(state), Some(strategy))
            } else {
                (None, None)
            };
            Some(ScheduleSuggestion::new(
                SuggestionKind::InsertWindow,
                diagnostic.edge_index(),
                diagnostic.from_state().clone(),
                diagnostic.to_state().clone(),
                proposal_state,
                diagnostic.dominant_components().to_vec(),
                proposal_strategy,
                diagnostic.priority_score(),
                format!(
                    "{} overlap {} is below threshold {}",
                    focus_component_prefix(diagnostic),
                    format_reason_number(diagnostic.overlap_min()),
                    format_reason_number(options.overlap_min)
                ),
            ))
        }
        EdgeSeverity::AddSampling => Some(ScheduleSuggestion::new(
            SuggestionKind::ExtendSampling,
            diagnostic.edge_index(),
            diagnostic.from_state().clone(),
            diagnostic.to_state().clone(),
            None,
            diagnostic.dominant_components().to_vec(),
            None,
            diagnostic.priority_score(),
            sampling_reason(diagnostic, options),
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

fn proposed_insert_state(
    diagnostic: &AdjacentEdgeDiagnostic,
) -> Result<(StatePoint, ProposalStrategy)> {
    let dominant = dominant_component_indices(diagnostic.from_state(), diagnostic.to_state());
    let changed = changed_component_indices(diagnostic.from_state(), diagnostic.to_state());
    if dominant.len() == 1 && changed.len() > 1 {
        let focus_index = dominant[0];
        let lambdas = diagnostic
            .from_state()
            .lambdas()
            .iter()
            .zip(diagnostic.to_state().lambdas().iter())
            .enumerate()
            .map(|(index, (from, to))| {
                if index == focus_index {
                    0.5 * (from + to)
                } else {
                    *from
                }
            })
            .collect::<Vec<_>>();
        Ok((
            StatePoint::new(lambdas, diagnostic.from_state().temperature_k())?,
            ProposalStrategy::FocusedSplit,
        ))
    } else {
        Ok((
            midpoint_state(diagnostic.from_state(), diagnostic.to_state())?,
            ProposalStrategy::Midpoint,
        ))
    }
}

fn focus_component_prefix(diagnostic: &AdjacentEdgeDiagnostic) -> String {
    if diagnostic.dominant_components().is_empty() {
        "adjacent edge".to_string()
    } else {
        format!(
            "components [{}]:",
            diagnostic.dominant_components().join(", ")
        )
    }
}

fn format_reason_number(value: f64) -> String {
    format!("{value:.3}")
}

fn sampling_reason(
    diagnostic: &AdjacentEdgeDiagnostic,
    options: &ScheduleAdvisorOptions,
) -> String {
    let block_cv_trigger = diagnostic
        .block_cv()
        .filter(|value| *value >= options.block_cv_min);
    let relative_uncertainty_trigger = diagnostic
        .relative_uncertainty()
        .filter(|value| *value >= 1.5);

    match (block_cv_trigger, relative_uncertainty_trigger) {
        (Some(block_cv), Some(relative_uncertainty)) => format!(
            "{} block coefficient of variation {} exceeds threshold {} and relative uncertainty {} exceeds attention threshold {}",
            focus_component_prefix(diagnostic),
            format_reason_number(block_cv),
            format_reason_number(options.block_cv_min),
            format_reason_number(relative_uncertainty),
            format_reason_number(1.5)
        ),
        (Some(block_cv), None) => format!(
            "{} block coefficient of variation {} exceeds threshold {}",
            focus_component_prefix(diagnostic),
            format_reason_number(block_cv),
            format_reason_number(options.block_cv_min)
        ),
        (None, Some(relative_uncertainty)) => format!(
            "{} relative uncertainty {} exceeds attention threshold {}",
            focus_component_prefix(diagnostic),
            format_reason_number(relative_uncertainty),
            format_reason_number(1.5)
        ),
        (None, None) => {
            if let Some(relative_uncertainty) = diagnostic.relative_uncertainty() {
                format!(
                    "{} relative uncertainty {} suggests additional sampling",
                    focus_component_prefix(diagnostic),
                    format_reason_number(relative_uncertainty)
                )
            } else if let Some(block_cv) = diagnostic.block_cv() {
                format!(
                    "{} block coefficient of variation {} suggests additional sampling",
                    focus_component_prefix(diagnostic),
                    format_reason_number(block_cv)
                )
            } else {
                format!(
                    "{} sampling stability metrics suggest additional sampling",
                    focus_component_prefix(diagnostic)
                )
            }
        }
    }
}

fn changed_component_indices(from_state: &StatePoint, to_state: &StatePoint) -> Vec<usize> {
    from_state
        .lambdas()
        .iter()
        .zip(to_state.lambdas().iter())
        .enumerate()
        .filter_map(|(index, (from, to))| {
            if (to - from).abs() > 1.0e-12 {
                Some(index)
            } else {
                None
            }
        })
        .collect()
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
    use super::{
        advise_lambda_schedule, advise_ti_schedule, recommend_ti_method, schedule_suggestion,
        AdjacentEdgeDiagnostic, AdvisorEstimator, EdgeSeverity, ProposalStrategy,
        ScheduleAdvisorOptions, TiEdgeSeverity, TiScheduleAdvisorOptions,
    };
    use crate::data::{StatePoint, UNkMatrix};
    use crate::{DhdlSeries, IntegrationMethod};

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

    fn build_window_with_labels(
        sampled_lambdas: Vec<f64>,
        rows: &[[f64; 3]],
        evaluated_states: &[StatePoint],
        labels: Vec<String>,
    ) -> UNkMatrix {
        let n_samples = rows.len();
        let data = rows
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let time = (0..n_samples).map(|idx| idx as f64).collect::<Vec<_>>();
        UNkMatrix::new_with_labels(
            n_samples,
            3,
            data,
            time,
            Some(StatePoint::new(sampled_lambdas, 298.0).unwrap()),
            evaluated_states.to_vec(),
            Some(labels),
        )
        .unwrap()
    }

    fn build_variable_window(
        sampled_lambda: f64,
        rows: &[Vec<f64>],
        evaluated_lambdas: &[f64],
    ) -> UNkMatrix {
        let n_samples = rows.len();
        let n_states = evaluated_lambdas.len();
        let data = rows
            .iter()
            .flat_map(|row| {
                assert_eq!(row.len(), n_states);
                row.iter().copied()
            })
            .collect::<Vec<_>>();
        let time = (0..n_samples).map(|idx| idx as f64).collect::<Vec<_>>();
        let evaluated_states = evaluated_lambdas
            .iter()
            .map(|&lambda| StatePoint::new(vec![lambda], 298.0).unwrap())
            .collect::<Vec<_>>();
        UNkMatrix::new(
            n_samples,
            n_states,
            data,
            time,
            Some(StatePoint::new(vec![sampled_lambda], 298.0).unwrap()),
            evaluated_states,
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
                    [3.0, 0.1, 0.0],
                    [3.1, 0.1, 0.0],
                    [2.9, 0.2, 0.0],
                    [3.0, 0.1, 0.0],
                ],
                &states,
            ),
            build_window(
                0.0,
                &[
                    [0.0, 2.0, 3.0],
                    [0.0, 2.1, 3.1],
                    [0.0, 1.9, 2.9],
                    [0.0, 2.0, 3.0],
                ],
                &states,
            ),
            build_window(
                0.5,
                &[
                    [2.0, 0.0, 0.1],
                    [2.1, 0.0, 0.1],
                    [1.9, 0.0, 0.2],
                    [2.0, 0.0, 0.1],
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
        assert_eq!(advice.suggestions().len(), 1);
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
                    [0.0, 0.1, 0.2],
                    [0.0, 0.1, 0.2],
                    [0.0, 0.1, 0.2],
                    [0.0, 0.1, 0.2],
                ],
                &states,
            ),
            build_window(
                0.5,
                &[
                    [0.1, 0.0, 0.1],
                    [0.1, 0.0, 0.1],
                    [0.1, 0.0, 0.1],
                    [0.1, 0.0, 0.1],
                ],
                &states,
            ),
            build_window(
                1.0,
                &[
                    [0.2, 0.1, 0.0],
                    [0.2, 0.1, 0.0],
                    [0.2, 0.1, 0.0],
                    [0.2, 0.1, 0.0],
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

    #[test]
    fn advise_lambda_schedule_bar_supports_local_neighbor_grids() {
        let windows = vec![
            build_variable_window(
                0.0,
                &[
                    vec![0.0, 0.2],
                    vec![0.0, 0.1],
                    vec![0.0, 0.2],
                    vec![0.0, 0.1],
                ],
                &[0.0, 0.5],
            ),
            build_variable_window(
                0.5,
                &[
                    vec![0.2, 0.0, 0.2],
                    vec![0.1, 0.0, 0.1],
                    vec![0.2, 0.0, 0.2],
                    vec![0.1, 0.0, 0.1],
                ],
                &[0.0, 0.5, 1.0],
            ),
            build_variable_window(
                1.0,
                &[
                    vec![0.2, 0.0],
                    vec![0.1, 0.0],
                    vec![0.2, 0.0],
                    vec![0.1, 0.0],
                ],
                &[0.5, 1.0],
            ),
        ];

        let advice = advise_lambda_schedule(
            &windows,
            Some(ScheduleAdvisorOptions {
                estimator: AdvisorEstimator::Bar,
                overlap_min: 0.0,
                block_cv_min: 1.0e30,
                n_blocks: 2,
                suggest_midpoints: false,
            }),
        )
        .unwrap();

        assert_eq!(advice.edges().len(), 2);
        assert!(advice.suggestions().is_empty());
        assert_eq!(advice.edges()[0].from_state().lambdas(), &[0.0]);
        assert_eq!(advice.edges()[0].to_state().lambdas(), &[0.5]);
        assert_eq!(advice.edges()[1].from_state().lambdas(), &[0.5]);
        assert_eq!(advice.edges()[1].to_state().lambdas(), &[1.0]);
        assert!(advice
            .edges()
            .iter()
            .all(|edge| edge.overlap_forward().is_finite() && edge.overlap_reverse().is_finite()));
    }

    #[test]
    fn advise_lambda_schedule_reports_dominant_multidimensional_components() {
        let states = vec![
            StatePoint::new(vec![0.0, 0.0], 298.0).unwrap(),
            StatePoint::new(vec![0.8, 0.1], 298.0).unwrap(),
            StatePoint::new(vec![1.0, 1.0], 298.0).unwrap(),
        ];
        let labels = vec!["coul-lambda".to_string(), "vdw-lambda".to_string()];
        let windows = vec![
            build_window_with_labels(
                vec![0.0, 0.0],
                &[
                    [0.0, 2.0, 3.0],
                    [0.0, 2.1, 3.1],
                    [0.0, 1.9, 2.9],
                    [0.0, 2.0, 3.0],
                ],
                &states,
                labels.clone(),
            ),
            build_window_with_labels(
                vec![0.8, 0.1],
                &[
                    [2.0, 0.0, 0.1],
                    [2.1, 0.0, 0.1],
                    [1.9, 0.0, 0.2],
                    [2.0, 0.0, 0.1],
                ],
                &states,
                labels.clone(),
            ),
            build_window_with_labels(
                vec![1.0, 1.0],
                &[
                    [3.0, 0.1, 0.0],
                    [3.1, 0.1, 0.0],
                    [2.9, 0.2, 0.0],
                    [3.0, 0.1, 0.0],
                ],
                &states,
                labels.clone(),
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

        assert_eq!(
            advice.edges()[0].dominant_components(),
            &["coul-lambda".to_string()]
        );
        let suggestion = advice
            .suggestions()
            .iter()
            .find(|suggestion| suggestion.edge_index() == 0)
            .expect("edge 0 suggestion");
        assert_eq!(suggestion.focus_components(), &["coul-lambda".to_string()]);
        assert_eq!(
            suggestion.proposal_strategy(),
            Some(ProposalStrategy::FocusedSplit)
        );
        assert_eq!(suggestion.proposed_state().unwrap().lambdas(), &[0.4, 0.0]);
    }

    #[test]
    fn advise_ti_schedule_marks_linear_series_healthy() {
        let s0 = StatePoint::new(vec![0.0], 298.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 298.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 298.0).unwrap();
        let series = vec![
            DhdlSeries::new(s0, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
            DhdlSeries::new(s1, vec![0.0, 1.0, 2.0, 3.0], vec![1.0, 1.0, 1.0, 1.0]).unwrap(),
            DhdlSeries::new(s2, vec![0.0, 1.0, 2.0, 3.0], vec![2.0, 2.0, 2.0, 2.0]).unwrap(),
        ];

        let advice = advise_ti_schedule(
            &series,
            Some(TiScheduleAdvisorOptions {
                n_blocks: 2,
                block_cv_min: 1.0e30,
                curvature_z_min: 10.0,
                slope_z_min: 10.0,
                interval_uncertainty_z_min: 10.0,
                suggest_midpoints: true,
            }),
        )
        .unwrap();

        assert_eq!(advice.intervals().len(), 2);
        assert!(advice.suggestions().is_empty());
        assert!(advice
            .intervals()
            .iter()
            .all(|interval| interval.severity() == TiEdgeSeverity::Healthy));
    }

    #[test]
    fn advise_ti_schedule_can_suggest_new_window() {
        let s0 = StatePoint::new(vec![0.0], 298.0).unwrap();
        let s1 = StatePoint::new(vec![0.33], 298.0).unwrap();
        let s2 = StatePoint::new(vec![0.66], 298.0).unwrap();
        let s3 = StatePoint::new(vec![1.0], 298.0).unwrap();
        let series = vec![
            DhdlSeries::new(s0, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
            DhdlSeries::new(s1, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
            DhdlSeries::new(s2, vec![0.0, 1.0, 2.0, 3.0], vec![10.0, 10.0, 10.0, 10.0]).unwrap(),
            DhdlSeries::new(s3, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
        ];

        let advice = advise_ti_schedule(
            &series,
            Some(TiScheduleAdvisorOptions {
                n_blocks: 2,
                block_cv_min: 1.0e30,
                curvature_z_min: 0.1,
                slope_z_min: 0.1,
                interval_uncertainty_z_min: 10.0,
                suggest_midpoints: true,
            }),
        )
        .unwrap();

        assert!(!advice.suggestions().is_empty());
        assert!(advice
            .suggestions()
            .iter()
            .any(|suggestion| suggestion.proposed_state().is_some()));
    }

    #[test]
    fn advise_ti_schedule_split_half_drift_contributes_to_priority() {
        let s0 = StatePoint::new(vec![0.0], 298.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 298.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 298.0).unwrap();
        let series = vec![
            DhdlSeries::new(s0, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
            DhdlSeries::new(s1, vec![0.0, 1.0, 2.0, 3.0], vec![-5.0, -5.0, 5.0, 5.0]).unwrap(),
            DhdlSeries::new(s2, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
        ];

        let advice = advise_ti_schedule(
            &series,
            Some(TiScheduleAdvisorOptions {
                n_blocks: 2,
                block_cv_min: 1.0e30,
                curvature_z_min: 10.0,
                slope_z_min: 10.0,
                interval_uncertainty_z_min: 10.0,
                suggest_midpoints: true,
            }),
        )
        .unwrap();

        assert_eq!(advice.intervals().len(), 2);
        assert!(advice
            .intervals()
            .iter()
            .all(|interval| interval.severity() == TiEdgeSeverity::AddSampling));
        assert!(advice
            .intervals()
            .iter()
            .all(|interval| interval.priority_score() > 0.0));
    }

    #[test]
    fn recommend_ti_method_prefers_gaussian_quadrature_for_exact_nodes() {
        let l0 = 0.21132486540518713;
        let l1 = 0.7886751345948129;
        let series = vec![
            DhdlSeries::new(
                StatePoint::new(vec![l0], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![1.0, 1.0, 1.0],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![l1], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![2.0, 2.0, 2.0],
            )
            .unwrap(),
        ];

        let recommendation = recommend_ti_method(&series, None).unwrap();
        assert_eq!(
            recommendation.recommended_method(),
            IntegrationMethod::GaussianQuadrature
        );
        assert!(recommendation.reason().contains("Gauss-Legendre"));
    }

    #[test]
    fn recommend_ti_method_prefers_simpson_on_healthy_uniform_grid() {
        let series = vec![
            DhdlSeries::new(
                StatePoint::new(vec![0.0], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![0.0, 0.0, 0.0],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![0.5], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![1.0, 1.0, 1.0],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![1.0], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![2.0, 2.0, 2.0],
            )
            .unwrap(),
        ];

        let recommendation = recommend_ti_method(&series, None).unwrap();
        assert_eq!(
            recommendation.recommended_method(),
            IntegrationMethod::Simpson
        );
    }

    #[test]
    fn recommend_ti_method_prefers_pchip_for_healthy_monotone_nonuniform_grid() {
        let series = vec![
            DhdlSeries::new(
                StatePoint::new(vec![0.0], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![0.3, 0.3, 0.3],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![0.2], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![0.5, 0.5, 0.5],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![0.55], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![0.8, 0.8, 0.8],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![1.0], 298.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![1.1, 1.1, 1.1],
            )
            .unwrap(),
        ];

        let recommendation = recommend_ti_method(&series, None).unwrap();
        assert_eq!(
            recommendation.recommended_method(),
            IntegrationMethod::Pchip
        );
        let simpson = recommendation
            .assessments()
            .iter()
            .find(|assessment| assessment.method() == IntegrationMethod::Simpson)
            .expect("simpson assessment");
        assert!(!simpson.supported());
    }

    #[test]
    fn recommend_ti_method_falls_back_to_trapezoidal_for_under_resolved_schedule() {
        let s0 = StatePoint::new(vec![0.0], 298.0).unwrap();
        let s1 = StatePoint::new(vec![0.33], 298.0).unwrap();
        let s2 = StatePoint::new(vec![0.66], 298.0).unwrap();
        let s3 = StatePoint::new(vec![1.0], 298.0).unwrap();
        let series = vec![
            DhdlSeries::new(s0, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
            DhdlSeries::new(s1, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
            DhdlSeries::new(s2, vec![0.0, 1.0, 2.0, 3.0], vec![10.0, 10.0, 10.0, 10.0]).unwrap(),
            DhdlSeries::new(s3, vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
        ];

        let recommendation = recommend_ti_method(&series, None).unwrap();
        assert_eq!(
            recommendation.recommended_method(),
            IntegrationMethod::Trapezoidal
        );
    }

    #[test]
    fn add_sampling_reason_uses_relative_uncertainty_when_block_cv_is_below_threshold() {
        let from = StatePoint::new(vec![0.3], 298.0).unwrap();
        let to = StatePoint::new(vec![0.4], 298.0).unwrap();
        let diagnostic = AdjacentEdgeDiagnostic::new(
            5,
            from.clone(),
            to.clone(),
            None,
            0.058,
            0.058,
            -6.056,
            Some(0.636),
            None,
            None,
            Some(0.053_021),
            Some(0.667),
            Some(1.811),
            vec!["lambda[0]".to_string()],
            1.072,
            EdgeSeverity::AddSampling,
        )
        .unwrap();

        let suggestion = schedule_suggestion(
            &diagnostic,
            &ScheduleAdvisorOptions {
                estimator: AdvisorEstimator::Mbar,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 4,
                suggest_midpoints: true,
            },
        )
        .unwrap()
        .unwrap();

        assert_eq!(
            suggestion.reason(),
            "components [lambda[0]]: relative uncertainty 1.811 exceeds attention threshold 1.500"
        );
    }
}
