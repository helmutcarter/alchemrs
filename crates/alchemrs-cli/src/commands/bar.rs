use std::path::PathBuf;

use alchemrs_estimators::{BarEstimator, BarOptions, MbarOptions};

use crate::cli::{BarMethodArg, OutputFormat, OutputUnits};
use crate::input::{load_windows, AnalysisInputOptions};
use crate::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::overlap::summarize_overlap;
use crate::CliResult;

pub struct BarRunOptions {
    pub method: BarMethodArg,
    pub output_units: OutputUnits,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
    pub overlap_summary: bool,
    pub parallel: bool,
}

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: BarRunOptions,
) -> CliResult<()> {
    let windows = load_windows(inputs, input_options)?;
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
    let estimator = BarEstimator::new(BarOptions {
        method: run_options.method.into(),
        parallel: run_options.parallel,
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
            units: run_options.output_units,
            temperature: input_options.temperature,
            overlap,
            provenance: OutputProvenance {
                estimator: "bar",
                decorrelate: input_options.decorrelate,
                remove_burnin: input_options.remove_burnin,
                auto_equilibrate: input_options.auto_equilibrate,
                fast: input_options.fast,
                conservative: input_options.conservative,
                nskip: input_options.nskip,
            },
        },
        run_options.output_format,
        run_options.output_path.as_deref(),
    )?;
    Ok(())
}
