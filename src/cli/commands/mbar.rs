use std::path::PathBuf;

use alchemrs::{MbarEstimator, MbarOptions, MbarSolver};

use crate::cli::input::{load_windows, AnalysisInputOptions};
use crate::cli::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::cli::overlap::summarize_overlap_matrix;
use crate::cli::{OutputFormat, OutputUnits};
use crate::CliResult;

pub struct MbarRunOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub solver: MbarSolver,
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
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: run_options.max_iterations,
        tolerance: run_options.tolerance,
        parallel: run_options.parallel,
        solver: run_options.solver,
        ..MbarOptions::default()
    });
    let fit = estimator.fit(&windows)?;
    let result = if run_options.no_uncertainty {
        fit.result()?
    } else {
        fit.result_with_uncertainty()?
    };
    let overlap = if run_options.overlap_summary {
        Some(summarize_overlap_matrix(&fit.overlap_matrix()?)?)
    } else {
        None
    };
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
                ti_method: None,
                ti_method_reason: None,
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
