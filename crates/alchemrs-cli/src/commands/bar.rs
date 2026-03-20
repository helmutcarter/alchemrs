use std::path::PathBuf;

use alchemrs_estimators::{BarEstimator, BarOptions};

use crate::cli::{BarMethodArg, OutputUnits};
use crate::input::{load_windows, AnalysisInputOptions};
use crate::output::print_scalar_result;
use crate::CliResult;

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    method: BarMethodArg,
    output_units: OutputUnits,
    parallel: bool,
) -> CliResult<()> {
    let windows = load_windows(inputs, input_options)?;
    let estimator = BarEstimator::new(BarOptions {
        method: method.into(),
        parallel,
        ..BarOptions::default()
    });
    let result = estimator.fit(&windows)?;
    let delta_index = result.n_states() - 1;

    print_scalar_result(
        result.values()[delta_index],
        result.uncertainties().map(|u| u[delta_index]),
        result.states().first().unwrap().lambdas()[0],
        result.states().last().unwrap().lambdas()[0],
        output_units,
        input_options.temperature,
    );
    Ok(())
}
