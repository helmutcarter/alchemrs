use std::path::PathBuf;

use alchemrs::{extract_nes_mbar_trajectory, NesMbarEstimator, NesMbarOptions};

use crate::cli::input::AnalysisSampleCounts;
use crate::cli::output::{print_scalar_result, OutputProvenance, ScalarResult};
use crate::cli::{OutputFormat, OutputUnits};
use crate::CliResult;

pub struct NesMbarRunOptions {
    pub n_bootstrap: usize,
    pub seed: u64,
    pub sample_stride: usize,
    pub output_units: OutputUnits,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
}

pub fn run(
    inputs: Vec<PathBuf>,
    temperature: f64,
    run_options: NesMbarRunOptions,
) -> CliResult<()> {
    let mut trajectories = Vec::with_capacity(inputs.len());
    let mut total_samples_in = 0usize;
    for path in inputs {
        let trajectory = extract_nes_mbar_trajectory(path, temperature)?;
        total_samples_in += trajectory.samples().len();
        trajectories.push(trajectory);
    }

    let estimator = NesMbarEstimator::new(NesMbarOptions {
        n_bootstrap: run_options.n_bootstrap,
        seed: run_options.seed,
        sample_stride: run_options.sample_stride,
    });
    let result = estimator.estimate(&trajectories)?;
    let n_states = result.n_states();
    let delta_index = n_states - 1;
    let samples_kept = trajectories
        .iter()
        .map(|trajectory| {
            trajectory
                .samples()
                .iter()
                .step_by(run_options.sample_stride)
                .count()
        })
        .sum();

    print_scalar_result(
        &ScalarResult {
            delta: result.values()[delta_index],
            sigma: result.uncertainties().map(|u| u[delta_index]),
            from_state: result.states().first().unwrap().lambdas().to_vec(),
            to_state: result.states().last().unwrap().lambdas().to_vec(),
            units: run_options.output_units,
            temperature,
            overlap: None,
            provenance: OutputProvenance {
                estimator: "nes-mbar",
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
                samples_in: total_samples_in,
                samples_after_burnin: total_samples_in,
                samples_kept,
            },
        },
        run_options.output_format,
        run_options.output_path.as_deref(),
    )?;

    Ok(())
}
