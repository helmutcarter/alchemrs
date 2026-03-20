use std::env;

use alchemrs::{
    decorrelate_u_nk_with_observable, extract_u_nk_with_potential, overlap_matrix, overlap_scalar,
    DecorrelationOptions, MbarEstimator, MbarOptions, UNkMatrix,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let temp_arg = args
        .next()
        .ok_or("usage: amber_mbar <temperature_k> <amber.out> [amber.out ...]")?;
    let temperature_k: f64 = temp_arg.parse()?;

    let mut windows = Vec::new();
    for path in args {
        let (u_nk, epot) = extract_u_nk_with_potential(&path, temperature_k)?;
        let decorrelated =
            decorrelate_u_nk_with_observable(&u_nk, &epot, &DecorrelationOptions::default())?;
        windows.push(decorrelated);
    }
    if windows.len() < 2 {
        return Err("need at least two AMBER output files for MBAR".into());
    }

    run_mbar(&windows)
}

fn run_mbar(windows: &[UNkMatrix]) -> Result<(), Box<dyn std::error::Error>> {
    let estimator = MbarEstimator::new(MbarOptions::default());
    let result = estimator.fit(windows)?;
    let delta_index = result.n_states() - 1;
    let delta_f = result.values()[delta_index];
    let uncertainty = result.uncertainties().map(|values| values[delta_index]);
    let overlap = overlap_matrix(windows, Some(MbarOptions::default()))?;
    let overlap = overlap_scalar(&overlap)?;

    match uncertainty {
        Some(sigma) => println!(
            "MBAR dG = {:.6} +/- {:.6} kT (from lambda={:.4} to lambda={:.4}, overlap={:.6})",
            delta_f,
            sigma,
            result.states().first().unwrap().lambdas()[0],
            result.states().last().unwrap().lambdas()[0],
            overlap
        ),
        None => println!(
            "MBAR dG = {:.6} kT (from lambda={:.4} to lambda={:.4}, overlap={:.6})",
            delta_f,
            result.states().first().unwrap().lambdas()[0],
            result.states().last().unwrap().lambdas()[0],
            overlap
        ),
    }

    Ok(())
}
