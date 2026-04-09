use std::path::PathBuf;

use alchemrs::{extract_nes_trajectory, NesEstimator, NesOptions};

use crate::cli::input::AnalysisSampleCounts;
use crate::cli::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::cli::{OutputFormat, OutputUnits};
use crate::CliResult;

pub struct NesRunOptions {
    pub n_bootstrap: usize,
    pub seed: u64,
    pub no_uncertainty: bool,
    pub output_units: OutputUnits,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
}

pub fn run(inputs: Vec<PathBuf>, temperature: f64, run_options: NesRunOptions) -> CliResult<()> {
    let mut trajectories = Vec::with_capacity(inputs.len());
    for path in inputs {
        trajectories.push(extract_nes_trajectory(path, temperature)?);
    }

    let estimator = NesEstimator::new(NesOptions {
        n_bootstrap: if run_options.no_uncertainty {
            0
        } else {
            run_options.n_bootstrap
        },
        seed: run_options.seed,
    });
    let result = estimator.estimate(&trajectories)?;

    print_scalar_result(
        &ScalarResult {
            delta: result.delta_f(),
            sigma: if run_options.no_uncertainty {
                None
            } else {
                result.uncertainty()
            },
            from_state: result.from_state().lambdas().to_vec(),
            to_state: result.to_state().lambdas().to_vec(),
            units: run_options.output_units,
            temperature,
            overlap: None,
            provenance: OutputProvenance {
                estimator: "nes",
                decorrelate: false,
                remove_burnin: 0,
                auto_equilibrate: false,
                fast: false,
                conservative: true,
                nskip: 1,
                u_nk_observable: None,
                ti_method: None,
                ti_method_reason: None,
                lambda_components: None,
            },
            sample_counts: AnalysisSampleCounts {
                windows: trajectories.len(),
                samples_in: trajectories.len(),
                samples_after_burnin: trajectories.len(),
                samples_kept: trajectories.len(),
            },
        },
        run_options.output_format,
        run_options.output_path.as_deref(),
    )?;

    Ok(())
}
