use std::path::PathBuf;

use alchemrs::{recommend_ti_method, IntegrationMethod, TiEstimator, TiOptions};

use crate::cli::input::{load_dhdl_series, AnalysisInputOptions};
use crate::cli::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::cli::{OutputFormat, OutputUnits, TiMethod};
use crate::CliResult;

pub struct TiRunOptions {
    pub method: TiMethod,
    pub output_units: OutputUnits,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
    pub parallel: bool,
}

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: TiRunOptions,
) -> CliResult<()> {
    let loaded = load_dhdl_series(inputs, input_options)?;
    let series = loaded.series;
    let (resolved_method, method_reason) = match run_options.method {
        TiMethod::Auto => {
            let recommendation = recommend_ti_method(&series, None)?;
            (
                recommendation.recommended_method(),
                Some(recommendation.reason().to_string()),
            )
        }
        explicit => (
            explicit
                .explicit_method()
                .expect("non-auto TI method must resolve"),
            None,
        ),
    };
    let estimator = TiEstimator::new(TiOptions {
        method: resolved_method,
        parallel: run_options.parallel,
    });
    let fit = estimator.fit(&series)?;
    let result = fit.result()?;

    let scalar = ScalarResult {
        delta: result.delta_f(),
        sigma: result.uncertainty(),
        from_state: result.from_state().lambdas().to_vec(),
        to_state: result.to_state().lambdas().to_vec(),
        units: run_options.output_units,
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
            ti_method: Some(render_ti_method(fit.method())),
            ti_method_reason: method_reason,
            lambda_components: None,
        },
        sample_counts: loaded.sample_counts,
    };

    print_scalar_result(
        &scalar,
        run_options.output_format,
        run_options.output_path.as_deref(),
    )?;
    Ok(())
}

fn render_ti_method(method: IntegrationMethod) -> &'static str {
    method.as_str()
}
