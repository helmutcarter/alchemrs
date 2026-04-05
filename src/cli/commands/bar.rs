use std::path::PathBuf;

use alchemrs::{BarEstimator, BarOptions, MbarOptions};

use crate::cli::input::{load_windows, AnalysisInputOptions};
use crate::cli::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::cli::overlap::summarize_overlap;
use crate::cli::{BarMethodArg, OutputFormat, OutputUnits};
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
    let estimator = BarEstimator::new(BarOptions {
        method: run_options.method.into(),
        parallel: run_options.parallel,
        ..BarOptions::default()
    });
    let fit = estimator.fit(&windows)?;
    let result = fit.result()?;
    let delta_index = result.n_states() - 1;

    print_scalar_result(
        &ScalarResult {
            delta: result.values()[delta_index],
            sigma: result.uncertainties().map(|u| u[delta_index]),
            from_state: result.states().first().unwrap().lambdas().to_vec(),
            to_state: result.states().last().unwrap().lambdas().to_vec(),
            units: run_options.output_units,
            temperature: input_options.temperature,
            overlap,
            provenance: OutputProvenance {
                estimator: "bar",
                decorrelate: input_options.decorrelate,
                remove_burnin: input_options.remove_burnin,
                auto_equilibrate: input_options.auto_equilibrate,
                fast: input_options.effective_fast(),
                conservative: input_options.effective_conservative(),
                nskip: input_options.nskip,
                u_nk_observable: input_options.u_nk_observable_name(),
                lambda_components: windows
                    .first()
                    .and_then(|window| window.lambda_labels().map(|labels| labels.to_vec())),
            },
            sample_counts: loaded.sample_counts,
        },
        run_options.output_format,
        run_options.output_path.as_deref(),
    )?;
    Ok(())
}
