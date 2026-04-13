use std::env;
use std::path::PathBuf;

use alchemrs::{
    extract_nes_mbar_trajectory, extract_nes_trajectory, NesEstimator, NesMbarEstimator,
    NesMbarOptions, NesMbarSample, NesMbarTrajectory, NesMbarWeighting,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (temperature, paths) = parse_args()?;
    if paths.is_empty() {
        return Err("no input files provided".into());
    }

    let switching: Vec<_> = paths
        .iter()
        .map(|path| extract_nes_trajectory(path, temperature))
        .collect::<Result<_, _>>()?;
    let nes_anchor = NesEstimator::default().fit(&switching)?.result()?;

    let trajectories: Vec<_> = paths
        .iter()
        .map(|path| extract_nes_mbar_trajectory(path, temperature))
        .collect::<Result<_, _>>()?;

    let mut mean_work_diff = 0.0f64;
    let mut max_work_diff = 0.0f64;
    for (nes, nes_mbar) in switching.iter().zip(trajectories.iter()) {
        let diff = (nes.reduced_work() - nes_mbar.samples().last().unwrap().reduced_work()).abs();
        mean_work_diff += diff;
        max_work_diff = max_work_diff.max(diff);
    }
    mean_work_diff /= trajectories.len() as f64;

    println!("anchor_nes_delta_f_kT,{:.12}", nes_anchor.delta_f());
    println!("mean_terminal_work_diff_kT,{mean_work_diff:.12}");
    println!("max_terminal_work_diff_kT,{max_work_diff:.12}");
    println!("variant,parameter,delta_f_kT,delta_vs_nes_kT,samples_kept");

    for stride in [1usize, 2, 5, 10, 20, 50, 100] {
        let fit = NesMbarEstimator::new(NesMbarOptions {
            n_bootstrap: 0,
            seed: 0,
            sample_stride: stride,
            weighting: NesMbarWeighting::Uniform,
            max_slice_trajectory_fraction: None,
            min_slice_ess_fraction: None,
        })
        .estimate(&trajectories)?;
        let delta = endpoint_delta(&fit);
        let samples_kept: usize = trajectories
            .iter()
            .map(|trajectory| trajectory.samples().iter().step_by(stride).count())
            .sum();
        println!(
            "stride,{stride},{delta:.12},{:.12},{samples_kept}",
            delta - nes_anchor.delta_f()
        );
    }

    for stride in [1usize, 10, 20, 50, 100] {
        let fit = NesMbarEstimator::new(NesMbarOptions {
            n_bootstrap: 0,
            seed: 0,
            sample_stride: stride,
            weighting: NesMbarWeighting::SliceEss,
            max_slice_trajectory_fraction: None,
            min_slice_ess_fraction: None,
        })
        .estimate(&trajectories)?;
        let delta = endpoint_delta(&fit);
        let samples_kept: usize = trajectories
            .iter()
            .map(|trajectory| trajectory.samples().iter().step_by(stride).count())
            .sum();
        println!(
            "slice_ess,{stride},{delta:.12},{:.12},{samples_kept}",
            delta - nes_anchor.delta_f()
        );
    }

    for stride in [1usize, 10, 20, 50] {
        let fit = NesMbarEstimator::new(NesMbarOptions {
            n_bootstrap: 0,
            seed: 0,
            sample_stride: stride,
            weighting: NesMbarWeighting::SliceEssFiltered,
            max_slice_trajectory_fraction: Some(0.25),
            min_slice_ess_fraction: Some(0.05),
        })
        .estimate(&trajectories)?;
        let delta = endpoint_delta(&fit);
        let samples_kept: usize = trajectories
            .iter()
            .map(|trajectory| trajectory.samples().iter().step_by(stride).count())
            .sum();
        println!(
            "slice_ess_filtered,{stride},{delta:.12},{:.12},{samples_kept}",
            delta - nes_anchor.delta_f()
        );
    }

    for block in [5usize, 10, 20, 50] {
        let end_trajectories: Vec<_> = trajectories
            .iter()
            .map(|trajectory| block_end_trajectory(trajectory, block))
            .collect::<Result<_, _>>()?;
        let fit = NesMbarEstimator::default().estimate(&end_trajectories)?;
        let delta = endpoint_delta(&fit);
        let samples_kept: usize = end_trajectories.iter().map(|t| t.samples().len()).sum();
        println!(
            "block_end,{block},{delta:.12},{:.12},{samples_kept}",
            delta - nes_anchor.delta_f()
        );

        let avg_trajectories: Vec<_> = trajectories
            .iter()
            .map(|trajectory| block_average_trajectory(trajectory, block))
            .collect::<Result<_, _>>()?;
        let fit = NesMbarEstimator::default().estimate(&avg_trajectories)?;
        let delta = endpoint_delta(&fit);
        let samples_kept: usize = avg_trajectories.iter().map(|t| t.samples().len()).sum();
        println!(
            "block_avg,{block},{delta:.12},{:.12},{samples_kept}",
            delta - nes_anchor.delta_f()
        );
    }

    Ok(())
}

fn parse_args() -> Result<(f64, Vec<PathBuf>), Box<dyn std::error::Error>> {
    let mut args = env::args_os().skip(1);
    let mut temperature = 300.0;
    let mut paths = Vec::new();
    while let Some(arg) = args.next() {
        if arg == "--temperature" {
            let value = args
                .next()
                .ok_or("missing value for --temperature")?
                .into_string()
                .map_err(|_| "temperature must be valid UTF-8")?;
            temperature = value.parse::<f64>()?;
            continue;
        }
        paths.push(PathBuf::from(arg));
    }
    Ok((temperature, paths))
}

fn endpoint_delta(matrix: &alchemrs::DeltaFMatrix) -> f64 {
    matrix.values()[matrix.n_states() - 1]
}

fn block_end_trajectory(
    trajectory: &NesMbarTrajectory,
    block_size: usize,
) -> Result<NesMbarTrajectory, Box<dyn std::error::Error>> {
    if block_size == 0 {
        return Err("block_size must be positive".into());
    }
    let samples = trajectory
        .samples()
        .chunks(block_size)
        .map(|chunk| chunk.last().unwrap().clone())
        .collect::<Vec<_>>();
    Ok(NesMbarTrajectory::new(
        trajectory.initial_state().clone(),
        trajectory.final_state().clone(),
        trajectory.target_states().to_vec(),
        samples,
    )?)
}

fn block_average_trajectory(
    trajectory: &NesMbarTrajectory,
    block_size: usize,
) -> Result<NesMbarTrajectory, Box<dyn std::error::Error>> {
    if block_size == 0 {
        return Err("block_size must be positive".into());
    }
    let n_targets = trajectory.target_states().len();
    let mut samples = Vec::new();
    for chunk in trajectory.samples().chunks(block_size) {
        let mut lambda = 0.0;
        let mut energy_protocol = 0.0;
        let mut energies_states = vec![0.0; n_targets];
        for sample in chunk {
            lambda += sample.lambda_protocol();
            energy_protocol += sample.reduced_energy_protocol();
            for (acc, value) in energies_states
                .iter_mut()
                .zip(sample.reduced_energies_states().iter())
            {
                *acc += *value;
            }
        }
        let denom = chunk.len() as f64;
        lambda /= denom;
        energy_protocol /= denom;
        for value in &mut energies_states {
            *value /= denom;
        }
        let last = chunk.last().unwrap();
        samples.push(NesMbarSample::new(
            last.step_index(),
            last.time_ps(),
            lambda,
            last.reduced_work(),
            energy_protocol,
            energies_states,
        )?);
    }
    Ok(NesMbarTrajectory::new(
        trajectory.initial_state().clone(),
        trajectory.final_state().clone(),
        trajectory.target_states().to_vec(),
        samples,
    )?)
}
