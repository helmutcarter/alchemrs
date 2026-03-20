use std::path::PathBuf;

use alchemrs_estimators::{ExpEstimator, ExpOptions};

use crate::cli::{OutputFormat, OutputUnits};
use crate::input::{load_windows, AnalysisInputOptions};
use crate::output::{print_scalar_result, ScalarResult};
use crate::CliResult;

pub fn run_forward(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    no_uncertainty: bool,
    output_units: OutputUnits,
    output_format: OutputFormat,
    parallel: bool,
) -> CliResult<()> {
    run(
        inputs,
        input_options,
        no_uncertainty,
        output_units,
        output_format,
        parallel,
        false,
    )
}

pub fn run_reverse(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    no_uncertainty: bool,
    output_units: OutputUnits,
    output_format: OutputFormat,
    parallel: bool,
) -> CliResult<()> {
    run(
        inputs,
        input_options,
        no_uncertainty,
        output_units,
        output_format,
        parallel,
        true,
    )
}

fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    no_uncertainty: bool,
    output_units: OutputUnits,
    output_format: OutputFormat,
    parallel: bool,
    reverse: bool,
) -> CliResult<()> {
    let windows = load_windows(inputs, input_options)?;
    let estimator = ExpEstimator::new(ExpOptions {
        compute_uncertainty: !no_uncertainty,
        parallel,
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
            units: output_units,
            temperature: input_options.temperature,
        },
        output_format,
    );
    Ok(())
}
