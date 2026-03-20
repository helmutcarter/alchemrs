use std::path::PathBuf;

use alchemrs_estimators::{TiEstimator, TiOptions};

use crate::cli::{OutputUnits, TiMethod};
use crate::input::{load_dhdl_series, AnalysisInputOptions};
use crate::output::print_scalar_result;
use crate::CliResult;

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    method: TiMethod,
    output_units: OutputUnits,
    parallel: bool,
) -> CliResult<()> {
    let series = load_dhdl_series(inputs, input_options)?;
    let estimator = TiEstimator::new(TiOptions {
        method: method.into(),
        parallel,
    });
    let result = estimator.fit(&series)?;

    print_scalar_result(
        result.delta_f(),
        result.uncertainty(),
        result.from_state().lambdas()[0],
        result.to_state().lambdas()[0],
        output_units,
        input_options.temperature,
    );
    Ok(())
}
