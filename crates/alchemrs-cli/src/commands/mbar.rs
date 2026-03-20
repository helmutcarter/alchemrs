use std::path::PathBuf;

use alchemrs_estimators::{MbarEstimator, MbarOptions};

use crate::cli::OutputUnits;
use crate::input::{load_windows, AnalysisInputOptions};
use crate::output::print_scalar_result;
use crate::CliResult;

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    max_iterations: usize,
    tolerance: f64,
    no_uncertainty: bool,
    output_units: OutputUnits,
    parallel: bool,
) -> CliResult<()> {
    let windows = load_windows(inputs, input_options)?;
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations,
        tolerance,
        compute_uncertainty: !no_uncertainty,
        parallel,
        ..MbarOptions::default()
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
