# Library Guide

The top-level `alchemrs` crate re-exports the main types and functions needed for direct Rust use.

## Common imports

```rust
use alchemrs::{
    DecorrelationOptions, MbarEstimator, MbarOptions, TiEstimator, TiOptions,
    decorrelate_dhdl, decorrelate_u_nk_with_observable, extract_dhdl,
    extract_u_nk_with_potential, overlap_matrix, overlap_scalar,
};
```

## TI example

```rust
use alchemrs::estimators::{TiEstimator, TiOptions};
use alchemrs::parse::amber::extract_dhdl;

fn run_ti(paths: &[String], temperature_k: f64) -> Result<(), Box<dyn std::error::Error>> {
    let mut series = Vec::new();
    for path in paths {
        series.push(extract_dhdl(path, temperature_k)?);
    }

    let estimator = TiEstimator::new(TiOptions::default());
    let result = estimator.fit(&series)?;
    println!("TI dG = {:.6}", result.delta_f());
    Ok(())
}
```

## MBAR example with `EPtot` decorrelation

```rust
use alchemrs::{
    DecorrelationOptions, MbarEstimator, MbarOptions, UNkMatrix,
    decorrelate_u_nk_with_observable, extract_u_nk_with_potential,
    overlap_matrix, overlap_scalar,
};

fn load_windows(
    paths: &[String],
    temperature_k: f64,
) -> Result<Vec<UNkMatrix>, Box<dyn std::error::Error>> {
    let mut windows = Vec::new();
    for path in paths {
        let (u_nk, epot) = extract_u_nk_with_potential(path, temperature_k)?;
        let decorrelated =
            decorrelate_u_nk_with_observable(&u_nk, &epot, &DecorrelationOptions::default())?;
        windows.push(decorrelated);
    }
    Ok(windows)
}

fn run_mbar(windows: &[UNkMatrix]) -> Result<(), Box<dyn std::error::Error>> {
    let estimator = MbarEstimator::new(MbarOptions::default());
    let result = estimator.fit(windows)?;
    let delta_index = result.n_states() - 1;
    println!("MBAR dG = {:.6}", result.values()[delta_index]);
    if let Some(labels) = result.lambda_labels() {
        println!("lambda components = {:?}", labels);
    }

    let overlap = overlap_matrix(windows, Some(MbarOptions::default()))?;
    println!("overlap = {:.6}", overlap_scalar(&overlap)?);
    if let Some(labels) = overlap.lambda_labels() {
        println!("overlap components = {:?}", labels);
    }
    Ok(())
}
```

For GROMACS multidimensional schedules, parser-derived component names are available from the parsed windows too:

```rust
if let Some(labels) = windows[0].lambda_labels() {
    println!("window components = {:?}", labels);
}
```

## Choosing the right preprocessing entry point

- Use `decorrelate_dhdl` for TI inputs.
- Use `decorrelate_u_nk(..., UNkSeriesMethod::DE, ...)` when you want an internal `u_nk`-derived scalar series.
- Use `decorrelate_u_nk_with_observable` when you already have a finite external observable such as `EPtot`.

## Top-level re-exports

The root crate re-exports:

- `analysis`
- `data`
- `error`
- `estimators`
- `parse`
- `prep`

You can use either:

- the direct re-exported functions and types
- the crate namespaces for more explicit imports

depending on your style preference.
