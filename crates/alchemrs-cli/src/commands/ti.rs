use std::path::PathBuf;

use alchemrs_estimators::{TiEstimator, TiOptions};

use crate::cli::{OutputFormat, OutputUnits, TiMethod};
use crate::input::{load_dhdl_series, AnalysisInputOptions};
use crate::output::{print_scalar_result, ScalarResult};
use crate::CliResult;

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    method: TiMethod,
    output_units: OutputUnits,
    output_format: OutputFormat,
    output_path: Option<PathBuf>,
    parallel: bool,
) -> CliResult<()> {
    let series = load_dhdl_series(inputs, input_options)?;
    let estimator = TiEstimator::new(TiOptions {
        method: method.into(),
        parallel,
    });
    let result = estimator.fit(&series)?;

    print_scalar_result(
        &ScalarResult {
            delta: result.delta_f(),
            sigma: result.uncertainty(),
            from_lambda: result.from_state().lambdas()[0],
            to_lambda: result.to_state().lambdas()[0],
            units: output_units,
            temperature: input_options.temperature,
            overlap: None,
        },
        output_format,
        output_path.as_deref(),
    )?;
    Ok(())
}
