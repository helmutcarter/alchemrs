use std::env;
use std::path::PathBuf;

use alchemrs::{extract_nes_mbar_trajectory, NesMbarEstimator};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (temperature, paths) = parse_args()?;
    if paths.is_empty() {
        return Err("no input files provided".into());
    }

    let trajectories = paths
        .iter()
        .map(|path| extract_nes_mbar_trajectory(path, temperature))
        .collect::<Result<Vec<_>, _>>()?;

    let diagnostics = NesMbarEstimator::default().diagnostics(&trajectories)?;
    for state in diagnostics.states.iter().skip(1) {
        println!("state_lambda,{:.6}", state.state.lambdas()[0]);
        println!("effective_sample_size,{:.6}", state.effective_sample_size);
        println!("max_sample_fraction,{:.6}", state.max_sample_fraction);
        println!("top_slices");
        for contribution in &state.top_slices {
            println!("{},{}", contribution.index, contribution.fraction);
        }
        println!("top_trajectories");
        for contribution in &state.top_trajectories {
            println!("{},{}", contribution.index, contribution.fraction);
        }
        println!();
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
