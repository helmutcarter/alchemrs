//! `alchemrs` is a Rust-first toolkit for alchemical free energy analysis.
//!
//! The crate provides:
//!
//! - AMBER and GROMACS parsers for `dH/dlambda` and `u_nk` data
//! - preprocessing utilities for burn-in removal, equilibration detection, and decorrelation
//! - estimators including TI, BAR, MBAR, EXP, and DEXP
//! - overlap diagnostics built on top of MBAR weights
//!
//! Most users can work directly with the top-level re-exports, while more explicit module
//! namespaces remain available through [`parse`], [`prep`], [`estimators`], [`analysis`],
//! [`data`], and [`error`].
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
//! let result = MbarEstimator::new(MbarOptions::default()).fit(&[u_nk])?;
//!
//! println!("delta_f = {}", result.values()[result.n_states() - 1]);
//! # Ok(())
//! # }
//! ```
//!
pub mod analysis;
pub mod data;
pub mod error;
pub mod estimators;
pub mod parse;
pub mod prep;

pub use analysis::{overlap_eigenvalues, overlap_matrix, overlap_scalar};
pub use data::{
    DeltaFMatrix, DhdlSeries, FreeEnergyEstimate, OverlapMatrix, StatePoint, UNkMatrix,
};
pub use error::{CoreError, Result};
pub use estimators::{BarEstimator, BarMethod, BarOptions, BarUncertainty};
pub use estimators::{ExpEstimator, ExpOptions, MbarEstimator, MbarOptions};
pub use estimators::{IntegrationMethod, TiEstimator, TiOptions};
pub use parse::{extract_dhdl, extract_u_nk, extract_u_nk_with_potential};
pub use prep::{
    decorrelate_dhdl, decorrelate_u_nk, decorrelate_u_nk_with_observable,
    detect_equilibration_dhdl, detect_equilibration_observable, detect_equilibration_u_nk,
    DecorrelationOptions, EquilibrationResult, UNkSeriesMethod,
};
