use std::env;

use alchemrs::estimators::{TiEstimator, TiOptions};
use alchemrs::parse::amber::extract_dhdl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let temp_arg = args
        .next()
        .ok_or("usage: amber_ti <temperature_k> <amber.out> [amber.out ...]")?;
    let temperature_k: f64 = temp_arg.parse()?;

    let mut series = Vec::new();
    for path in args {
        series.push(extract_dhdl(&path, temperature_k)?);
    }
    if series.len() < 2 {
        return Err("need at least two AMBER output files for TI".into());
    }

    let estimator = TiEstimator::new(TiOptions::default());
    let result = estimator.fit(&series)?;
    println!(
        "TI dG = {:.6} (from lambda={:.4} to lambda={:.4})",
        result.delta_f(),
        result.from_state().lambdas()[0],
        result.to_state().lambdas()[0]
    );

    Ok(())
}
