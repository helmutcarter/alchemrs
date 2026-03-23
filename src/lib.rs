pub mod analysis;
pub mod core;
pub mod estimators;
pub mod parse;
pub mod prep;

pub use analysis::{overlap_eigenvalues, overlap_matrix, overlap_scalar};
pub use core::{
    CoreError, DeltaFMatrix, DhdlSeries, FreeEnergyEstimate, OverlapMatrix, Result, StatePoint,
    UNkMatrix,
};
pub use estimators::{BarEstimator, BarMethod, BarOptions, BarUncertainty};
pub use estimators::{ExpEstimator, ExpOptions, MbarEstimator, MbarOptions};
pub use estimators::{IntegrationMethod, TiEstimator, TiOptions};
pub use parse::amber::{extract_dhdl, extract_u_nk, extract_u_nk_with_potential};
pub use prep::{
    decorrelate_dhdl, decorrelate_u_nk, decorrelate_u_nk_with_observable,
    detect_equilibration_dhdl, detect_equilibration_observable, detect_equilibration_u_nk,
    DecorrelationOptions, EquilibrationResult, UNkSeriesMethod,
};
