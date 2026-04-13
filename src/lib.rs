//! `alchemrs` is a CLI-first tool for alchemical free energy analysis with a
//! robust Rust API underneath.
//!
//! The library crate provides the scientific engine behind the `alchemrs`
//! command-line interface:
//!
//! - AMBER and GROMACS parsers for `dH/dlambda` and `u_nk` data
//! - preprocessing utilities for burn-in removal, equilibration detection, and decorrelation
//! - estimators including TI, BAR, MBAR, IEXP, DEXP, and NES
//! - overlap diagnostics built on top of MBAR weights
//! - optional SVG plotting helpers behind the `plotting` feature
//!
//! Most users should start with the CLI. For embedding, scripting, custom
//! workflows, and future bindings, the Rust API remains available through the
//! top-level re-exports and the explicit [`parse`], [`prep`], [`estimators`],
//! [`analysis`], [`data`], and [`error`] modules.
//!
//! # Example
//!
//! ```no_run
//! use alchemrs::{
//!     decorrelate_u_nk_with_observable, extract_u_nk_with_potential, DecorrelationOptions,
//!     MbarEstimator, MbarOptions,
//! };
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let (u_nk, epot) = extract_u_nk_with_potential("prod.out", 300.0)?;
//! let u_nk = decorrelate_u_nk_with_observable(&u_nk, &epot, &DecorrelationOptions::default())?;
//! let fit = MbarEstimator::new(MbarOptions::default()).fit(&[u_nk])?;
//! let result = fit.result_with_uncertainty()?;
//!
//! println!("delta_f = {}", result.values()[result.n_states() - 1]);
//! if let Some(labels) = fit.lambda_labels() {
//!     println!("lambda components = {:?}", labels);
//! }
//! # Ok(())
//! # }
//! ```
//!
pub mod analysis;
pub mod data;
pub mod error;
pub mod estimators;
pub mod parse;
#[cfg(feature = "plotting")]
pub mod plot;
pub mod prep;

pub use analysis::{
    advise_lambda_schedule, advise_nes, advise_ti_schedule, bar_convergence, dexp_convergence,
    exp_convergence, mbar_convergence, nes_convergence, overlap_eigenvalues, overlap_matrix,
    overlap_scalar, recommend_ti_method, ti_convergence, AdjacentEdgeDiagnostic, AdvisorEstimator,
    BlockEstimate, ConvergencePoint, EdgeSeverity, NesAdvice, NesAdvisorOptions, NesCurvaturePoint,
    NesProfilePoint, NesSuggestionKind, ProposalStrategy, ScheduleAdvice, ScheduleAdvisorOptions,
    ScheduleSuggestion, SuggestionKind, TiEdgeSeverity, TiIntervalDiagnostic, TiMethodAssessment,
    TiMethodRecommendation, TiMethodRecommendationOptions, TiScheduleAdvice,
    TiScheduleAdvisorOptions, TiScheduleSuggestion, TiSuggestionKind, TiWindowDiagnostic,
};
pub use data::{
    DeltaFMatrix, DhdlSeries, FreeEnergyEstimate, NesMbarSample, NesMbarTrajectory, OverlapMatrix,
    StatePoint, SwitchingTrajectory, UNkMatrix,
};
pub use error::{CoreError, Result};
pub use estimators::{BarEstimator, BarFit, BarMethod, BarOptions};
pub use estimators::{
    IexpEstimator, IexpFit, IexpOptions, MbarEstimator, MbarFit, MbarOptions, MbarSolver,
    NesEstimator, NesFit, NesMbarContribution, NesMbarDiagnostics, NesMbarEstimator, NesMbarFit,
    NesMbarOptions, NesMbarStateDiagnostics, NesMbarWeighting, NesOptions,
};
pub use estimators::{IntegrationMethod, TiEstimator, TiFit, TiOptions};
pub use parse::{
    extract_dhdl, extract_nes_mbar_trajectory, extract_nes_trajectory, extract_u_nk,
    extract_u_nk_with_potential,
};
#[cfg(feature = "plotting")]
pub use plot::{
    render_block_average_svg, render_convergence_svg, render_delta_f_state_svg,
    render_overlap_matrix_svg, render_ti_dhdl_svg, BlockAveragePlotOptions, ConvergencePlotOptions,
    DeltaFStatePlotOptions, OverlapPlotOptions, TiDhdlPlotOptions,
};
pub use prep::{
    decorrelate_dhdl, decorrelate_u_nk, decorrelate_u_nk_with_observable,
    detect_equilibration_dhdl, detect_equilibration_observable, detect_equilibration_u_nk,
    DecorrelationOptions, EquilibrationResult, UNkSeriesMethod,
};
