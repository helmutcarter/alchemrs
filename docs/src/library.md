# Library Guide

The top-level `alchemrs` crate re-exports the main types and functions needed for direct Rust use.

Use the library when the CLI is not enough:

- embedding `alchemrs` inside another Rust application
- building custom preprocessing or estimation pipelines
- sharing the same scientific core used by the existing Python bindings

## Common imports

```rust
use alchemrs::{
    AtmDirection, AtmEstimator, AtmSample, AtmSampleSet, AtmSchedule, AtmState,
    DecorrelationOptions, MbarEstimator, MbarOptions, NesEstimator, NesOptions, TiEstimator,
    TiOptions, UwhamEstimator,
    decorrelate_dhdl, decorrelate_u_nk_with_observable, extract_dhdl,
    extract_nes_trajectory, extract_u_nk_with_potential, mbar_convergence, nes_convergence,
};
```

## TI example

```rust
use alchemrs::{TiEstimator, TiOptions, extract_dhdl};

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

## UWHAM example

```rust
use alchemrs::{UwhamEstimator, extract_u_nk};

fn run_uwham(paths: &[String], temperature_k: f64) -> Result<(), Box<dyn std::error::Error>> {
    let mut windows = Vec::new();
    for path in paths {
        windows.push(extract_u_nk(path, temperature_k)?);
    }

    let fit = UwhamEstimator::default().fit(&windows)?;
    let result = fit.result_with_uncertainty()?;
    println!("UWHAM dG = {:.6}", result.values()[result.n_states() - 1]);
    Ok(())
}
```

`UwhamFit::write_reference_inputs(...)` is also available when you want to
export `logQ.csv`, `size.csv`, and the fitted reference outputs for comparison
against external UWHAM implementations.

## MBAR example with `EPtot` decorrelation

```rust
use alchemrs::{
    DecorrelationOptions, MbarEstimator, MbarOptions, MbarSolver, UNkMatrix,
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

`MbarOptions::default()` uses `MbarSolver::Lbfgs`. It automatically falls back to
fixed-point when the evaluated state grid contains zero sampled-count states. To
force fixed-point explicitly:

```rust
let estimator = MbarEstimator::new(MbarOptions {
    solver: MbarSolver::FixedPoint,
    ..MbarOptions::default()
});
```

## ATM leg example

```rust
use alchemrs::{
    AtmDirection, AtmEstimator, AtmSample, AtmSampleSet, AtmSchedule, AtmState,
};

fn run_atm_leg() -> Result<(), Box<dyn std::error::Error>> {
    let leg = AtmSampleSet::new(
        AtmSchedule::new(vec![
            AtmState::new(0, AtmDirection::Forward, 0.0, 0.0, None, 0.0, 0.0, None, 0.0, 300.0)?,
            AtmState::new(1, AtmDirection::Forward, 1.0, 1.0, None, 0.0, 0.0, None, 0.0, 300.0)?,
        ])?,
        vec![
            AtmSample::new(0, 10.0, 2.0)?,
            AtmSample::new(1, 12.0, 2.0)?,
        ],
    )?;

    let estimate = AtmEstimator::default().estimate_leg(&leg)?;
    println!("ATM leg dG = {:.6}", estimate.delta_f());
    println!("ATM leg sigma = {:?}", estimate.uncertainty());
    Ok(())
}
```

For paired binding legs, use `estimate_binding(...)`, `estimate_rbfe(...)`, or
`estimate_abfe(...)` on two `AtmSampleSet` values with opposite directions.

For GROMACS multidimensional schedules, parser-derived component names are available from the parsed windows too:

```rust
if let Some(labels) = windows[0].lambda_labels() {
    println!("window components = {:?}", labels);
}
```

## OpenMM `u_kln` conversion example

If you already have reduced potentials from an OpenMM workflow as a
`u_kln[k][l][n]` tensor, you do not need a dedicated parser. Convert each
sampled-state slice into a `UNkMatrix` window and pass the resulting windows to
the existing estimators.

The runnable example in `examples/openmm_u_kln_mbar.rs` demonstrates this
conversion end to end.

If you prefer Python, see [Python and OpenMM](python.md) for the bindings API,
OpenMM helper functions, and the runnable Python toy-system examples.

The same `analysis` module also provides plot-ready convergence series for TI, BAR, MBAR, IEXP,
DEXP, and NES:

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

## NES example

```rust
use alchemrs::{extract_nes_trajectory, NesEstimator, NesOptions};

fn run_nes(paths: &[String], temperature_k: f64) -> Result<(), Box<dyn std::error::Error>> {
    let mut trajectories = Vec::new();
    for path in paths {
        trajectories.push(extract_nes_trajectory(path, temperature_k)?);
    }

    let estimator = NesEstimator::new(NesOptions::default());
    let result = estimator.estimate(&trajectories)?;
    println!("NES dG = {:.6}", result.delta_f());
    println!("NES sigma = {:?}", result.uncertainty());
    Ok(())
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
