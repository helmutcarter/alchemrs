use std::path::PathBuf;

use alchemrs_estimators::{
    BarEstimator, BarMethod, BarOptions, ExpEstimator, ExpOptions, IntegrationMethod, MbarEstimator,
    MbarOptions, TiEstimator, TiOptions,
};
use alchemrs_parse::amber::{extract_dhdl, extract_u_nk};
use alchemrs_prep::{decorrelate_dhdl, decorrelate_u_nk, DecorrelationOptions, UNkSeriesMethod};
use clap::{Parser, Subcommand, ValueEnum};

#[derive(Debug, Parser)]
#[command(name = "alchemrs-cli")]
#[command(about = "Alchemical free energy analysis CLI")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
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
        /// Apply decorrelation to each window
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
        /// Apply decorrelation to each window
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
        /// Apply decorrelation to each window
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
    },
    Dexp {
        /// AMBER output files (one per lambda window)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K
        #[arg(long, default_value_t = 300.0)]
        temperature: f64,
        /// Apply decorrelation to each window
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
    },
    Mbar {
        /// AMBER output files (one per lambda window)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Temperature in K
        #[arg(long, default_value_t = 300.0)]
        temperature: f64,
        /// Apply decorrelation to each window
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
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum TiMethod {
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
enum BarMethodArg {
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let load_windows = |inputs: Vec<PathBuf>,
                        temperature: f64,
                        decorrelate: bool,
                        remove_burnin: usize,
                        auto_equilibrate: bool,
                        fast: bool,
                        conservative: bool,
                        nskip: usize|
     -> Result<Vec<alchemrs_core::UNkMatrix>, Box<dyn std::error::Error>> {
        let mut windows = Vec::with_capacity(inputs.len());
        for path in inputs {
            let mut u_nk = extract_u_nk(path, temperature)?;
            if remove_burnin > 0 {
                let n_samples = u_nk.n_samples();
                if remove_burnin >= n_samples {
                    return Err("remove-burnin exceeds series length".into());
                }
                let n_states = u_nk.n_states();
                let data = u_nk.data()[(remove_burnin * n_states)..].to_vec();
                let time = u_nk.time_ps()[remove_burnin..].to_vec();
                u_nk = alchemrs_core::UNkMatrix::new(
                    n_samples - remove_burnin,
                    n_states,
                    data,
                    time,
                    u_nk.sampled_state().cloned(),
                    u_nk.evaluated_states().to_vec(),
                )?;
            }
            if decorrelate {
                let options = DecorrelationOptions {
                    remove_burnin: auto_equilibrate,
                    fast,
                    conservative,
                    nskip,
                    ..DecorrelationOptions::default()
                };
                u_nk = decorrelate_u_nk(&u_nk, UNkSeriesMethod::DE, &options)?;
            }
            windows.push(u_nk);
        }
        Ok(windows)
    };

    match cli.command {
        Command::Ti {
            inputs,
            temperature,
            method,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
        } => {
            let mut series = Vec::with_capacity(inputs.len());
            for path in inputs {
                let mut dhdl = extract_dhdl(path, temperature)?;
                if remove_burnin > 0 {
                    let values = dhdl.values();
                    let time = dhdl.time_ps();
                    if remove_burnin >= values.len() {
                        return Err("remove-burnin exceeds series length".into());
                    }
                    let new_time = time[remove_burnin..].to_vec();
                    let new_values = values[remove_burnin..].to_vec();
                    dhdl = alchemrs_core::DhdlSeries::new(
                        dhdl.state().clone(),
                        new_time,
                        new_values,
                    )?;
                }
                if decorrelate {
                    let options = DecorrelationOptions {
                        remove_burnin: auto_equilibrate,
                        fast,
                        conservative,
                        nskip,
                        ..DecorrelationOptions::default()
                    };
                    dhdl = decorrelate_dhdl(&dhdl, &options)?;
                }
                series.push(dhdl);
            }

            let estimator = TiEstimator::new(TiOptions {
                method: method.into(),
                parallel: false,
            });
            let result = estimator.fit(&series)?;
            let from_lambda = result.from_state().lambdas()[0];
            let to_lambda = result.to_state().lambdas()[0];

            println!("delta_f: {}", result.delta_f());
            match result.uncertainty() {
                Some(value) => println!("uncertainty: {}", value),
                None => println!("uncertainty: none"),
            }
            println!("from_lambda: {}", from_lambda);
            println!("to_lambda: {}", to_lambda);
        }
        Command::Bar {
            inputs,
            temperature,
            method,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
        } => {
            let windows = load_windows(
                inputs,
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            )?;

            let estimator = BarEstimator::new(BarOptions {
                method: method.into(),
                ..BarOptions::default()
            });
            let result = estimator.fit(&windows)?;
            let n = result.n_states();
            let delta = result.values()[0 * n + (n - 1)];
            let sigma = result
                .uncertainties()
                .map(|u| u[0 * n + (n - 1)]);
            let from_lambda = result.states().first().unwrap().lambdas()[0];
            let to_lambda = result.states().last().unwrap().lambdas()[0];

            println!("delta_f: {}", delta);
            match sigma {
                Some(value) => println!("uncertainty: {}", value),
                None => println!("uncertainty: none"),
            }
            println!("from_lambda: {}", from_lambda);
            println!("to_lambda: {}", to_lambda);
        }
        Command::Exp {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            no_uncertainty,
        } => {
            let windows = load_windows(
                inputs,
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            )?;

            let estimator = ExpEstimator::new(ExpOptions {
                compute_uncertainty: !no_uncertainty,
                ..ExpOptions::default()
            });
            let result = estimator.fit(&windows)?;
            let n = result.n_states();
            let delta = result.values()[0 * n + (n - 1)];
            let sigma = result
                .uncertainties()
                .map(|u| u[0 * n + (n - 1)]);
            let from_lambda = result.states().first().unwrap().lambdas()[0];
            let to_lambda = result.states().last().unwrap().lambdas()[0];

            println!("delta_f: {}", delta);
            match sigma {
                Some(value) => println!("uncertainty: {}", value),
                None => println!("uncertainty: none"),
            }
            println!("from_lambda: {}", from_lambda);
            println!("to_lambda: {}", to_lambda);
        }
        Command::Dexp {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            no_uncertainty,
        } => {
            let windows = load_windows(
                inputs,
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            )?;

            let estimator = ExpEstimator::new(ExpOptions {
                compute_uncertainty: !no_uncertainty,
                ..ExpOptions::default()
            });
            let result = estimator.fit(&windows)?;
            let n = result.n_states();
            let delta = result.values()[(n - 1) * n + 0];
            let sigma = result
                .uncertainties()
                .map(|u| u[(n - 1) * n + 0]);
            let from_lambda = result.states().last().unwrap().lambdas()[0];
            let to_lambda = result.states().first().unwrap().lambdas()[0];

            println!("delta_f: {}", delta);
            match sigma {
                Some(value) => println!("uncertainty: {}", value),
                None => println!("uncertainty: none"),
            }
            println!("from_lambda: {}", from_lambda);
            println!("to_lambda: {}", to_lambda);
        }
        Command::Mbar {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            max_iterations,
            tolerance,
            no_uncertainty,
        } => {
            let windows = load_windows(
                inputs,
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            )?;

            let estimator = MbarEstimator::new(MbarOptions {
                max_iterations,
                tolerance,
                compute_uncertainty: !no_uncertainty,
                ..MbarOptions::default()
            });
            let result = estimator.fit(&windows)?;
            let n = result.n_states();
            let delta = result.values()[0 * n + (n - 1)];
            let sigma = result
                .uncertainties()
                .map(|u| u[0 * n + (n - 1)]);
            let from_lambda = result.states().first().unwrap().lambdas()[0];
            let to_lambda = result.states().last().unwrap().lambdas()[0];

            println!("delta_f: {}", delta);
            match sigma {
                Some(value) => println!("uncertainty: {}", value),
                None => println!("uncertainty: none"),
            }
            println!("from_lambda: {}", from_lambda);
            println!("to_lambda: {}", to_lambda);
        }
    }

    Ok(())
}
