use std::path::PathBuf;

use alchemrs::{MbarEstimator, MbarOptions};

use crate::cli::input::{load_windows, AnalysisInputOptions};
use crate::cli::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::cli::overlap::summarize_overlap;
use crate::cli::{OutputFormat, OutputUnits};
use crate::CliResult;

pub struct MbarRunOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub no_uncertainty: bool,
    pub output_units: OutputUnits,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
    pub overlap_summary: bool,
    pub parallel: bool,
}

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: MbarRunOptions,
) -> CliResult<()> {
    let loaded = load_windows(inputs, input_options)?;
    let windows = loaded.windows;
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
            from_state: result.states().first().unwrap().lambdas().to_vec(),
            to_state: result.states().last().unwrap().lambdas().to_vec(),
            units: run_options.output_units,
            temperature: input_options.temperature,
            overlap,
            provenance: OutputProvenance {
                estimator: "mbar",
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
