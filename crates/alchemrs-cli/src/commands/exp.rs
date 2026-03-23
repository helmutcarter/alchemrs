use std::path::PathBuf;

use alchemrs::{ExpEstimator, ExpOptions, MbarOptions};

use crate::cli::{OutputFormat, OutputUnits};
use crate::input::{load_windows, AnalysisInputOptions};
use crate::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::overlap::summarize_overlap;
use crate::CliResult;

pub struct ExpRunOptions {
    pub no_uncertainty: bool,
    pub output_units: OutputUnits,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
    pub overlap_summary: bool,
    pub parallel: bool,
}

pub fn run_forward(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: ExpRunOptions,
) -> CliResult<()> {
    run(inputs, input_options, run_options, false)
}

pub fn run_reverse(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: ExpRunOptions,
) -> CliResult<()> {
    run(inputs, input_options, run_options, true)
}

fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: ExpRunOptions,
    reverse: bool,
) -> CliResult<()> {
    let loaded = load_windows(inputs, input_options)?;
    let windows = loaded.windows;
    let overlap = if run_options.overlap_summary {
        Some(summarize_overlap(
            &windows,
            Some(MbarOptions {
                parallel: run_options.parallel,
                ..MbarOptions::default()
            }),
        )?)
    } else {
        None
    };
    let estimator = ExpEstimator::new(ExpOptions {
        compute_uncertainty: !run_options.no_uncertainty,
        parallel: run_options.parallel,
    });
    let result = estimator.fit(&windows)?;
    let n_states = result.n_states();
    let delta_index = if reverse {
        (n_states - 1) * n_states
    } else {
        n_states - 1
    };
    let (from_lambda, to_lambda) = if reverse {
        (
            result.states().last().unwrap().lambdas()[0],
            result.states().first().unwrap().lambdas()[0],
        )
    } else {
        (
            result.states().first().unwrap().lambdas()[0],
            result.states().last().unwrap().lambdas()[0],
        )
    };

    print_scalar_result(
        &ScalarResult {
            delta: result.values()[delta_index],
            sigma: result.uncertainties().map(|u| u[delta_index]),
            from_lambda,
            to_lambda,
            units: run_options.output_units,
            temperature: input_options.temperature,
            overlap,
            provenance: OutputProvenance {
                estimator: if reverse { "dexp" } else { "exp" },
                decorrelate: input_options.decorrelate,
                remove_burnin: input_options.remove_burnin,
                auto_equilibrate: input_options.auto_equilibrate,
                fast: input_options.effective_fast(),
                conservative: input_options.effective_conservative(),
                nskip: input_options.nskip,
                u_nk_observable: input_options.u_nk_observable_name(),
            },
            sample_counts: loaded.sample_counts,
        },
        run_options.output_format,
        run_options.output_path.as_deref(),
    )?;
    Ok(())
}
