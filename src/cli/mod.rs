use std::path::PathBuf;

use alchemrs::{BarMethod, IntegrationMethod, MbarSolver};
use clap::{ArgAction, Parser, Subcommand, ValueEnum};

pub mod commands;
pub mod input;
pub mod output;
pub mod overlap;

#[derive(Debug, Parser)]
#[command(name = "alchemrs")]
#[command(about = "Alchemical free energy analysis CLI")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    AdviseSchedule {
        /// Simulation output files (AMBER `.out` or GROMACS `dhdl.xvg`)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K; inferred from the input files when omitted
        #[arg(long)]
        temperature: Option<f64>,
        /// Estimator to use for adjacent edge estimates
        #[arg(long, value_enum, default_value_t = AdvisorEstimatorArg::Mbar)]
        estimator: AdvisorEstimatorArg,
        /// Apply decorrelation to each window using the selected u_nk observable
        #[arg(long)]
        decorrelate: bool,
        /// Skip this many initial samples before any analysis
        #[arg(long = "remove-burnin", default_value_t = 0)]
        remove_burnin: usize,
        /// Automatically detect equilibration and remove burn-in
        #[arg(long)]
        auto_equilibrate: bool,
        /// Use fast statistical inefficiency estimate
        #[arg(long)]
        fast: bool,
        /// Use conservative subsampling
        #[arg(
            long,
            action = ArgAction::Set,
            default_value_t = true,
            default_missing_value = "true",
            require_equals = true,
            num_args = 0..=1
        )]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
        /// Observable to use for u_nk auto-equilibration and decorrelation
        #[arg(long = "u-nk-observable", value_enum, default_value_t = UNkObservable::De)]
        u_nk_observable: UNkObservable,
        /// Force advise-schedule to treat the inputs as u_nk or dH/dlambda data
        #[arg(long = "input-kind", value_enum, default_value_t = AdviseInputKind::Auto)]
        input_kind: AdviseInputKind,
        /// Minimum adjacent overlap before suggesting a new window
        #[arg(long = "overlap-min", default_value_t = 0.03)]
        overlap_min: f64,
        /// Minimum block coefficient of variation before suggesting more sampling
        #[arg(long = "block-cv-min", default_value_t = 0.15)]
        block_cv_min: f64,
        /// Number of blocks for adjacent-edge block averaging
        #[arg(long = "n-blocks", default_value_t = 4)]
        n_blocks: usize,
        /// Disable midpoint proposals for insert-window suggestions
        #[arg(long = "no-midpoints")]
        no_midpoints: bool,
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        output_format: OutputFormat,
        /// Write output to a file instead of stdout
        #[arg(long)]
        output: Option<PathBuf>,
        /// Write an HTML advisor report to a file
        #[arg(long)]
        report: Option<PathBuf>,
    },
    Ti {
        /// Simulation output files (AMBER `.out` or GROMACS `dhdl.xvg`)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K; inferred from the input files when omitted
        #[arg(long)]
        temperature: Option<f64>,
        /// Integration method
        #[arg(long, value_enum, default_value_t = TiMethod::Trapezoidal)]
        method: TiMethod,
        /// Output units
        #[arg(long, value_enum, default_value_t = OutputUnits::KT)]
        output_units: OutputUnits,
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        output_format: OutputFormat,
        /// Write output to a file instead of stdout
        #[arg(long)]
        output: Option<PathBuf>,
        /// Enable parallel processing
        #[arg(long)]
        parallel: bool,
        /// Apply decorrelation to each window using the dH/dlambda series
        #[arg(long)]
        decorrelate: bool,
        /// Skip this many initial samples before any analysis
        #[arg(long = "remove-burnin", default_value_t = 0)]
        remove_burnin: usize,
        /// Automatically detect equilibration and remove burn-in
        #[arg(long)]
        auto_equilibrate: bool,
        /// Use fast statistical inefficiency estimate
        #[arg(long)]
        fast: bool,
        /// Use conservative subsampling
        #[arg(
            long,
            action = ArgAction::Set,
            default_value_t = true,
            default_missing_value = "true",
            require_equals = true,
            num_args = 0..=1
        )]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
        /// Accepted only to provide a clearer error; TI uses dH/dlambda, not u_nk
        #[arg(long = "u-nk-observable", value_enum, hide = true)]
        u_nk_observable: Option<UNkObservable>,
    },
    Bar {
        /// Simulation output files (AMBER `.out` or GROMACS `dhdl.xvg`)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K; inferred from the input files when omitted
        #[arg(long)]
        temperature: Option<f64>,
        /// BAR method
        #[arg(long, value_enum, default_value_t = BarMethodArg::FalsePosition)]
        method: BarMethodArg,
        /// Output units
        #[arg(long, value_enum, default_value_t = OutputUnits::KT)]
        output_units: OutputUnits,
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        output_format: OutputFormat,
        /// Write output to a file instead of stdout
        #[arg(long)]
        output: Option<PathBuf>,
        /// Include overlap scalar and eigenvalues in the output
        #[arg(long = "overlap-summary")]
        overlap_summary: bool,
        /// Enable parallel processing
        #[arg(long)]
        parallel: bool,
        /// Apply decorrelation to each window using the selected u_nk observable
        #[arg(long)]
        decorrelate: bool,
        /// Skip this many initial samples before any analysis
        #[arg(long = "remove-burnin", default_value_t = 0)]
        remove_burnin: usize,
        /// Automatically detect equilibration and remove burn-in
        #[arg(long)]
        auto_equilibrate: bool,
        /// Use fast statistical inefficiency estimate
        #[arg(long)]
        fast: bool,
        /// Use conservative subsampling
        #[arg(
            long,
            action = ArgAction::Set,
            default_value_t = true,
            default_missing_value = "true",
            require_equals = true,
            num_args = 0..=1
        )]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
        /// Observable to use for u_nk auto-equilibration and decorrelation
        #[arg(long = "u-nk-observable", value_enum, default_value_t = UNkObservable::De)]
        u_nk_observable: UNkObservable,
    },
    Exp {
        /// Simulation output files (AMBER `.out` or GROMACS `dhdl.xvg`)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K; inferred from the input files when omitted
        #[arg(long)]
        temperature: Option<f64>,
        /// Apply decorrelation to each window using the selected u_nk observable
        #[arg(long)]
        decorrelate: bool,
        /// Skip this many initial samples before any analysis
        #[arg(long = "remove-burnin", default_value_t = 0)]
        remove_burnin: usize,
        /// Automatically detect equilibration and remove burn-in
        #[arg(long)]
        auto_equilibrate: bool,
        /// Use fast statistical inefficiency estimate
        #[arg(long)]
        fast: bool,
        /// Use conservative subsampling
        #[arg(
            long,
            action = ArgAction::Set,
            default_value_t = true,
            default_missing_value = "true",
            require_equals = true,
            num_args = 0..=1
        )]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
        /// Observable to use for u_nk auto-equilibration and decorrelation
        #[arg(long = "u-nk-observable", value_enum, default_value_t = UNkObservable::De)]
        u_nk_observable: UNkObservable,
        /// Disable uncertainty estimation
        #[arg(long)]
        no_uncertainty: bool,
        /// Output units
        #[arg(long, value_enum, default_value_t = OutputUnits::KT)]
        output_units: OutputUnits,
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        output_format: OutputFormat,
        /// Write output to a file instead of stdout
        #[arg(long)]
        output: Option<PathBuf>,
        /// Include overlap scalar and eigenvalues in the output
        #[arg(long = "overlap-summary")]
        overlap_summary: bool,
        /// Enable parallel processing
        #[arg(long)]
        parallel: bool,
    },
    Dexp {
        /// Simulation output files (AMBER `.out` or GROMACS `dhdl.xvg`)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K; inferred from the input files when omitted
        #[arg(long)]
        temperature: Option<f64>,
        /// Apply decorrelation to each window using the selected u_nk observable
        #[arg(long)]
        decorrelate: bool,
        /// Skip this many initial samples before any analysis
        #[arg(long = "remove-burnin", default_value_t = 0)]
        remove_burnin: usize,
        /// Automatically detect equilibration and remove burn-in
        #[arg(long)]
        auto_equilibrate: bool,
        /// Use fast statistical inefficiency estimate
        #[arg(long)]
        fast: bool,
        /// Use conservative subsampling
        #[arg(
            long,
            action = ArgAction::Set,
            default_value_t = true,
            default_missing_value = "true",
            require_equals = true,
            num_args = 0..=1
        )]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
        /// Observable to use for u_nk auto-equilibration and decorrelation
        #[arg(long = "u-nk-observable", value_enum, default_value_t = UNkObservable::De)]
        u_nk_observable: UNkObservable,
        /// Disable uncertainty estimation
        #[arg(long)]
        no_uncertainty: bool,
        /// Output units
        #[arg(long, value_enum, default_value_t = OutputUnits::KT)]
        output_units: OutputUnits,
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        output_format: OutputFormat,
        /// Write output to a file instead of stdout
        #[arg(long)]
        output: Option<PathBuf>,
        /// Include overlap scalar and eigenvalues in the output
        #[arg(long = "overlap-summary")]
        overlap_summary: bool,
        /// Enable parallel processing
        #[arg(long)]
        parallel: bool,
    },
    Mbar {
        /// Simulation output files (AMBER `.out` or GROMACS `dhdl.xvg`)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K; inferred from the input files when omitted
        #[arg(long)]
        temperature: Option<f64>,
        /// Apply decorrelation to each window using the selected u_nk observable
        #[arg(long)]
        decorrelate: bool,
        /// Skip this many initial samples before any analysis
        #[arg(long = "remove-burnin", default_value_t = 0)]
        remove_burnin: usize,
        /// Automatically detect equilibration and remove burn-in
        #[arg(long)]
        auto_equilibrate: bool,
        /// Use fast statistical inefficiency estimate
        #[arg(long)]
        fast: bool,
        /// Use conservative subsampling
        #[arg(
            long,
            action = ArgAction::Set,
            default_value_t = true,
            default_missing_value = "true",
            require_equals = true,
            num_args = 0..=1
        )]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
        /// Observable to use for u_nk auto-equilibration and decorrelation
        #[arg(long = "u-nk-observable", value_enum, default_value_t = UNkObservable::De)]
        u_nk_observable: UNkObservable,
        /// Maximum MBAR iterations
        #[arg(long, default_value_t = 10_000)]
        max_iterations: usize,
        /// MBAR relative tolerance
        #[arg(long, default_value_t = 1.0e-7)]
        tolerance: f64,
        /// MBAR solver backend
        #[arg(long, value_enum, default_value_t = MbarCliSolver::FixedPoint)]
        solver: MbarCliSolver,
        /// Disable uncertainty estimation
        #[arg(long)]
        no_uncertainty: bool,
        /// Output units
        #[arg(long, value_enum, default_value_t = OutputUnits::KT)]
        output_units: OutputUnits,
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        output_format: OutputFormat,
        /// Write output to a file instead of stdout
        #[arg(long)]
        output: Option<PathBuf>,
        /// Include overlap scalar and eigenvalues in the output
        #[arg(long = "overlap-summary")]
        overlap_summary: bool,
        /// Enable parallel processing
        #[arg(long)]
        parallel: bool,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum TiMethod {
    Auto,
    Trapezoidal,
    Simpson,
    CubicSpline,
    Pchip,
    Akima,
    GaussianQuadrature,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum AdviseInputKind {
    Auto,
    UNk,
    Dhdl,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MbarCliSolver {
    FixedPoint,
    Lbfgs,
}

impl MbarCliSolver {
    pub fn solver(self) -> MbarSolver {
        match self {
            Self::FixedPoint => MbarSolver::FixedPoint,
            Self::Lbfgs => MbarSolver::Lbfgs,
        }
    }
}

impl TiMethod {
    pub fn explicit_method(self) -> Option<IntegrationMethod> {
        match self {
            TiMethod::Auto => None,
            TiMethod::Trapezoidal => Some(IntegrationMethod::Trapezoidal),
            TiMethod::Simpson => Some(IntegrationMethod::Simpson),
            TiMethod::CubicSpline => Some(IntegrationMethod::CubicSpline),
            TiMethod::Pchip => Some(IntegrationMethod::Pchip),
            TiMethod::Akima => Some(IntegrationMethod::Akima),
            TiMethod::GaussianQuadrature => Some(IntegrationMethod::GaussianQuadrature),
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum OutputUnits {
    KT,
    Kcal,
    Kj,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum OutputFormat {
    Text,
    Json,
    Csv,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AdvisorEstimatorArg {
    Mbar,
    Bar,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum UNkObservable {
    Epot,
    De,
    All,
}

impl UNkObservable {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Epot => "epot",
            Self::De => "de",
            Self::All => "all",
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum BarMethodArg {
    FalsePosition,
    SelfConsistentIteration,
    Bisection,
}

impl From<BarMethodArg> for BarMethod {
    fn from(method: BarMethodArg) -> Self {
        match method {
            BarMethodArg::FalsePosition => BarMethod::FalsePosition,
            BarMethodArg::SelfConsistentIteration => BarMethod::SelfConsistentIteration,
            BarMethodArg::Bisection => BarMethod::Bisection,
        }
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::{Cli, Command};

    #[test]
    fn commands_default_to_conservative_subsampling() {
        for args in [
            ["alchemrs", "ti", "window.out"].as_slice(),
            ["alchemrs", "bar", "window.out"].as_slice(),
            ["alchemrs", "exp", "window.out"].as_slice(),
            ["alchemrs", "dexp", "window.out"].as_slice(),
            ["alchemrs", "mbar", "window.out"].as_slice(),
            ["alchemrs", "advise-schedule", "window.out"].as_slice(),
        ] {
            let cli = Cli::parse_from(args);
            assert!(command_conservative(cli.command));
        }
    }

    #[test]
    fn commands_accept_explicit_false_for_conservative_subsampling() {
        for args in [
            ["alchemrs", "ti", "--conservative=false", "window.out"].as_slice(),
            ["alchemrs", "bar", "--conservative=false", "window.out"].as_slice(),
            ["alchemrs", "exp", "--conservative=false", "window.out"].as_slice(),
            ["alchemrs", "dexp", "--conservative=false", "window.out"].as_slice(),
            ["alchemrs", "mbar", "--conservative=false", "window.out"].as_slice(),
            [
                "alchemrs",
                "advise-schedule",
                "--conservative=false",
                "window.out",
            ]
            .as_slice(),
        ] {
            let cli = Cli::parse_from(args);
            assert!(!command_conservative(cli.command));
        }
    }

    #[test]
    fn ti_command_accepts_gaussian_quadrature_method() {
        let cli = Cli::parse_from([
            "alchemrs",
            "ti",
            "--method",
            "gaussian-quadrature",
            "window.out",
        ]);
        match cli.command {
            Command::Ti { method, .. } => {
                assert!(matches!(method, super::TiMethod::GaussianQuadrature));
            }
            _ => panic!("expected ti command"),
        }
    }

    #[test]
    fn ti_command_accepts_auto_method() {
        let cli = Cli::parse_from(["alchemrs", "ti", "--method", "auto", "window.out"]);
        match cli.command {
            Command::Ti { method, .. } => {
                assert!(matches!(method, super::TiMethod::Auto));
            }
            _ => panic!("expected ti command"),
        }
    }

    #[test]
    fn ti_command_accepts_cubic_spline_method() {
        let cli = Cli::parse_from(["alchemrs", "ti", "--method", "cubic-spline", "window.out"]);
        match cli.command {
            Command::Ti { method, .. } => {
                assert!(matches!(method, super::TiMethod::CubicSpline));
            }
            _ => panic!("expected ti command"),
        }
    }

    #[test]
    fn ti_command_accepts_pchip_method() {
        let cli = Cli::parse_from(["alchemrs", "ti", "--method", "pchip", "window.out"]);
        match cli.command {
            Command::Ti { method, .. } => {
                assert!(matches!(method, super::TiMethod::Pchip));
            }
            _ => panic!("expected ti command"),
        }
    }

    #[test]
    fn ti_command_accepts_akima_method() {
        let cli = Cli::parse_from(["alchemrs", "ti", "--method", "akima", "window.out"]);
        match cli.command {
            Command::Ti { method, .. } => {
                assert!(matches!(method, super::TiMethod::Akima));
            }
            _ => panic!("expected ti command"),
        }
    }

    fn command_conservative(command: Command) -> bool {
        match command {
            Command::AdviseSchedule { conservative, .. }
            | Command::Ti { conservative, .. }
            | Command::Bar { conservative, .. }
            | Command::Exp { conservative, .. }
            | Command::Dexp { conservative, .. }
            | Command::Mbar { conservative, .. } => conservative,
        }
    }
}
