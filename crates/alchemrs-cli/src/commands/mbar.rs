use std::path::PathBuf;

use alchemrs_estimators::{MbarEstimator, MbarOptions};

use crate::cli::{OutputFormat, OutputUnits};
use crate::input::{load_windows, AnalysisInputOptions};
use crate::output::{print_scalar_result, ScalarResult};
use crate::overlap::summarize_overlap;
use crate::CliResult;

pub struct MbarRunOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub no_uncertainty: bool,
    pub output_units: OutputUnits,
    pub output_format: OutputFormat,
    pub overlap_summary: bool,
    pub parallel: bool,
}

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: MbarRunOptions,
) -> CliResult<()> {
    let windows = load_windows(inputs, input_options)?;
    let overlap = if run_options.overlap_summary {
        Some(summarize_overlap(
            &windows,
            Some(MbarOptions {
                max_iterations: run_options.max_iterations,
                tolerance: run_options.tolerance,
                compute_uncertainty: false,
                parallel: run_options.parallel,
                ..MbarOptions::default()
            }),
        )?)
    } else {
        None
    };
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: run_options.max_iterations,
        tolerance: run_options.tolerance,
        compute_uncertainty: !run_options.no_uncertainty,
        parallel: run_options.parallel,
        ..MbarOptions::default()
    });
    let result = estimator.fit(&windows)?;
    let delta_index = result.n_states() - 1;

    print_scalar_result(
        &ScalarResult {
            delta: result.values()[delta_index],
            sigma: result.uncertainties().map(|u| u[delta_index]),
            from_lambda: result.states().first().unwrap().lambdas()[0],
            to_lambda: result.states().last().unwrap().lambdas()[0],
            units: run_options.output_units,
            temperature: input_options.temperature,
            overlap,
        },
        run_options.output_format,
    );
    Ok(())
}
