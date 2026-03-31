use std::path::PathBuf;

use alchemrs::{TiEstimator, TiOptions};

use crate::cli::input::{load_dhdl_series, AnalysisInputOptions};
use crate::cli::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::cli::{OutputFormat, OutputUnits, TiMethod};
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
    let loaded = load_dhdl_series(inputs, input_options)?;
    let series = loaded.series;
    let estimator = TiEstimator::new(TiOptions {
        method: method.into(),
        parallel,
    });
    let result = estimator.fit(&series)?;

    print_scalar_result(
        &ScalarResult {
            delta: result.delta_f(),
            sigma: result.uncertainty(),
            from_state: result.from_state().lambdas().to_vec(),
            to_state: result.to_state().lambdas().to_vec(),
            units: output_units,
            temperature: input_options.temperature,
            overlap: None,
            provenance: OutputProvenance {
                estimator: "ti",
                decorrelate: input_options.decorrelate,
                remove_burnin: input_options.remove_burnin,
                auto_equilibrate: input_options.auto_equilibrate,
                fast: input_options.effective_fast(),
                conservative: input_options.effective_conservative(),
                nskip: input_options.nskip,
                u_nk_observable: input_options.u_nk_observable_name(),
                lambda_components: None,
            },
            sample_counts: loaded.sample_counts,
        },
        output_format,
        output_path.as_deref(),
    )?;
    Ok(())
}
