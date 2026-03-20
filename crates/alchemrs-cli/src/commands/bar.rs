use std::path::PathBuf;

use alchemrs_estimators::{BarEstimator, BarOptions, MbarOptions};

use crate::cli::{BarMethodArg, OutputFormat, OutputUnits};
use crate::input::{load_windows, AnalysisInputOptions};
use crate::output::{print_scalar_result, ScalarResult};
use crate::overlap::summarize_overlap;
use crate::CliResult;

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    method: BarMethodArg,
    output_units: OutputUnits,
    output_format: OutputFormat,
    overlap_summary: bool,
    parallel: bool,
) -> CliResult<()> {
    let windows = load_windows(inputs, input_options)?;
    let overlap = if overlap_summary {
        Some(summarize_overlap(
            &windows,
            Some(MbarOptions {
                parallel,
                ..MbarOptions::default()
            }),
        )?)
    } else {
        None
    };
    let estimator = BarEstimator::new(BarOptions {
        method: method.into(),
        parallel,
        ..BarOptions::default()
    });
    let result = estimator.fit(&windows)?;
    let delta_index = result.n_states() - 1;

    print_scalar_result(
        &ScalarResult {
            delta: result.values()[delta_index],
            sigma: result.uncertainties().map(|u| u[delta_index]),
            from_lambda: result.states().first().unwrap().lambdas()[0],
            to_lambda: result.states().last().unwrap().lambdas()[0],
            units: output_units,
            temperature: input_options.temperature,
            overlap,
        },
        output_format,
    );
    Ok(())
}
