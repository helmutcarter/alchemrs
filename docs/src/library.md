# Library Guide

The top-level `alchemrs` crate re-exports the main types and functions needed for direct Rust use.

Use the library when the CLI is not enough:

- embedding `alchemrs` inside another Rust application
- building custom preprocessing or estimation pipelines
- preparing for later Python or other bindings on top of the same core implementation

## Common imports

```rust
use alchemrs::{
    DecorrelationOptions, MbarEstimator, MbarOptions, TiEstimator, TiOptions,
    decorrelate_dhdl, decorrelate_u_nk_with_observable, extract_dhdl,
    extract_u_nk_with_potential, mbar_convergence,
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
    let fit = estimator.fit(&series)?;
    let result = fit.result()?;
    println!("TI dG = {:.6}", result.delta_f());
    Ok(())
}
```

## MBAR example with `EPtot` decorrelation

```rust
use alchemrs::{
    DecorrelationOptions, MbarEstimator, MbarOptions, UNkMatrix,
    decorrelate_u_nk_with_observable, extract_u_nk_with_potential,
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
    let fit = estimator.fit(windows)?;
    let result = fit.result_with_uncertainty()?;
    let delta_index = result.n_states() - 1;
    println!("MBAR dG = {:.6}", result.values()[delta_index]);
    if let Some(labels) = fit.lambda_labels() {
        println!("lambda components = {:?}", labels);
    }

    let overlap = fit.overlap_matrix()?;
    println!("overlap = {:.6}", fit.overlap_scalar()?);
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

The same `analysis` module also provides plot-ready convergence series for TI, BAR, MBAR, IEXP,
and DEXP:

```rust
let points = mbar_convergence(windows, Some(MbarOptions::default()))?;
for point in points {
    println!(
        "windows={} delta_f={:.6} from={:?} to={:?}",
        point.n_windows(),
        point.delta_f(),
        point.from_state().lambdas(),
        point.to_state().lambdas()
    );
}
```

If you enable the optional `plotting` feature, those convergence points can be rendered directly
to SVG:

```rust
# #[cfg(feature = "plotting")]
# {
use alchemrs::{ConvergencePlotOptions, render_convergence_svg};

let svg = render_convergence_svg(
    &points,
    Some(ConvergencePlotOptions {
        title: "MBAR Convergence".to_string(),
        ..ConvergencePlotOptions::default()
    }),
)?;
assert!(svg.contains("<svg"));
# }
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
