use std::path::PathBuf;

use alchemrs_estimators::{BarMethod, IntegrationMethod};
use clap::{Parser, Subcommand, ValueEnum};

#[derive(Debug, Parser)]
#[command(name = "alchemrs-cli")]
#[command(about = "Alchemical free energy analysis CLI")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    Ti {
        /// AMBER output files (one per lambda window)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K
        #[arg(long, default_value_t = 300.0)]
        temperature: f64,
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
        #[arg(long, default_value_t = true)]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
    },
    Bar {
        /// AMBER output files (one per lambda window)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K
        #[arg(long, default_value_t = 300.0)]
        temperature: f64,
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
        /// Apply decorrelation to each window using EPtot from the AMBER output
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
        #[arg(long, default_value_t = true)]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
    },
    Exp {
        /// AMBER output files (one per lambda window)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K
        #[arg(long, default_value_t = 300.0)]
        temperature: f64,
        /// Apply decorrelation to each window using EPtot from the AMBER output
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
        #[arg(long, default_value_t = true)]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
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
        /// AMBER output files (one per lambda window)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K
        #[arg(long, default_value_t = 300.0)]
        temperature: f64,
        /// Apply decorrelation to each window using EPtot from the AMBER output
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
        #[arg(long, default_value_t = true)]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
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
        /// AMBER output files (one per lambda window)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K
        #[arg(long, default_value_t = 300.0)]
        temperature: f64,
        /// Apply decorrelation to each window using EPtot from the AMBER output
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
        #[arg(long, default_value_t = true)]
        conservative: bool,
        /// Subsample stride for equilibration detection
        #[arg(long, default_value_t = 1)]
        nskip: usize,
        /// Maximum MBAR iterations
        #[arg(long, default_value_t = 10_000)]
        max_iterations: usize,
        /// MBAR relative tolerance
        #[arg(long, default_value_t = 1.0e-7)]
        tolerance: f64,
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
    Trapezoidal,
    Simpson,
}

impl From<TiMethod> for IntegrationMethod {
    fn from(method: TiMethod) -> Self {
        match method {
            TiMethod::Trapezoidal => IntegrationMethod::Trapezoidal,
            TiMethod::Simpson => IntegrationMethod::Simpson,
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
