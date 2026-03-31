# alchemrs

`alchemrs` is a Rust-first toolkit for alchemical free energy analysis. The `alchemrs` package contains both the main library crate and the `alchemrs` command-line binary. The library provides modules for parsing AMBER outputs and GROMACS `dhdl.xvg` outputs, preprocessing time series (equilibration trimming and decorrelation), running common estimators (TI, BAR, MBAR, EXP/DEXP), and computing diagnostics like overlap analysis. Fixtures and tests compare results against established reference implementations (`alchemlyb`) to ensure scientific correctness.

Native SVG plotting is available as an optional `plotting` feature, currently covering
convergence traces and overlap heatmaps.

## Library API

The top-level `alchemrs` crate re-exports the common parse, prep, estimator, and analysis entry points:

```rust
use alchemrs::{
    decorrelate_u_nk_with_observable, extract_u_nk_with_potential, DecorrelationOptions,
    MbarEstimator, MbarOptions,
};

let (u_nk, epot) = extract_u_nk_with_potential("prod.out", 300.0)?;
let u_nk = decorrelate_u_nk_with_observable(&u_nk, &epot, &DecorrelationOptions::default())?;
let result = MbarEstimator::new(MbarOptions::default()).fit(&[u_nk])?;
if let Some(labels) = result.lambda_labels() {
    println!("lambda components = {:?}", labels);
}
```

The repo also includes runnable top-level examples:

- `cargo run --example amber_ti -- 300 path/to/lambda0.out path/to/lambda1.out`
- `cargo run --example amber_mbar -- 300 path/to/lambda0.out path/to/lambda1.out path/to/lambda2.out`

Optional plotting helpers can be enabled with:

```bash
cargo test --features plotting
```

## Documentation

Documentation is online at https://helmutcarter.github.io/alchemrs/. The source files are in [`docs/`](docs/) and can be viewed locally:

```bash
cargo install mdbook
mdbook serve docs
```

## CLI

The `alchemrs` binary provides a command-line workflow.

### Build
It can either be invoked through cargo:
```bash
cargo run --release [arguments]
```
or compiled using cargo and called from the binary:
```bash
cargo build --release
./target/release/alchemrs [arguments]
```

### Common flags

- `--remove-burnin <N>`: skip the first `N` samples in each window before analysis.
- `--auto-equilibrate`: use `pymbar`-style equilibration detection as described [here](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00784).
- `--decorrelate`: subsample to reduce correlation before estimating.
  - `ti` uses the parsed `dH/dlambda` series.
  - `bar`, `exp`, `dexp`, and `mbar` use the observable selected by `--u-nk-observable <de|all|epot>`.
- `--u-nk-observable <de|all|epot>`: choose the scalar observable used for `u_nk` auto-equilibration and decorrelation on `bar`, `exp`, `dexp`, and `mbar` runs. The default is `de`.
- `--output-units <kt|kcal|kj>`: output energy units (default `kt`).
- `--output-format <text|json|csv>`: output format for estimator results (default `text`).
- `--output <PATH>`: write the formatted result to a file instead of stdout.
- `--overlap-summary`: include overlap scalar and overlap eigenvalues for BAR/EXP/DEXP/MBAR runs.

### Structured output

JSON output is useful for shell pipelines and downstream tools:

```bash
alchemrs mbar \
  --temperature 300 \
  --decorrelate \
  --overlap-summary \
  --output-format json \
  --output results.json \
  /path/to/*/prod.out
```

For `u_nk`-based estimators (`bar`, `exp`, `dexp`, `mbar`), the observable selected by `--u-nk-observable` is used for both `--auto-equilibrate` and `--decorrelate`, and any retained indices are then applied back to the parsed `u_nk` samples. The default is `de`; `epot` uses an engine-provided potential-energy observable.

Example JSON output:

```json
{
  "delta_f": -113.58,
  "uncertainty": 1.17,
  "from_lambda": 0.0,
  "to_lambda": 1.0,
  "units": "kT",
  "overlap": {
    "scalar": 0.0183,
    "eigenvalues": [1.0, 0.9817]
  },
  "provenance": {
    "estimator": "mbar",
    "temperature_k": 300.0,
    "decorrelate": true,
    "remove_burnin": 0,
    "auto_equilibrate": false,
    "fast": false,
    "conservative": true,
    "nskip": 1,
    "u_nk_observable": "de",
    "lambda_components": null,
    "windows": 15,
    "samples_in": 300,
    "samples_after_burnin": 300,
    "samples_kept": 126
  }
}
```

For multidimensional GROMACS schedules, the parsed `UNkMatrix`, estimator `DeltaFMatrix`, and `OverlapMatrix` all preserve parser-derived lambda component names when available through `.lambda_labels()`.

CSV output is useful for quick ingestion into spreadsheets or tabular tools:

```bash
alchemrs bar \
  --temperature 300 \
  --output-format csv \
  /path/to/*/prod.out
```

CSV columns include estimator parameters after the result fields:

```text
delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues,estimator,temperature_k,decorrelate,remove_burnin,auto_equilibrate,fast,conservative,nskip,u_nk_observable,lambda_components,windows,samples_in,samples_after_burnin,samples_kept
```

### TI (trapezoidal)

```bash
alchemrs ti \
  --temperature 300 \
  --method trapezoidal \
  --remove-burnin 125 \
  /path/to/*/prod.out
```

### BAR

```bash
alchemrs bar \
  --temperature 300 \
  --method false-position \
  --decorrelate \
  /path/to/*/prod.out
```

Note: BAR uncertainties are only computed for adjacent windows; non-adjacent state pairs (for example, `0->1`) are reported as `NaN`, matching `alchemlyb`.

### MBAR

```bash
alchemrs mbar \
  --temperature 300 \
  --decorrelate \
  --max-iterations 10000 \
  --tolerance 1e-7 \
  /path/to/*/prod.out
```

### EXP / DEXP

```bash
alchemrs exp \
  --temperature 300 \
  /path/to/*/prod.out

alchemrs dexp \
  --temperature 300 \
  /path/to/*/prod.out
```

EXP reports FEP results in the forward direction, DEXP reports FEP results in the reverse direction.

## Performance 
Initial results demonstrate 6x performance improvement over `alchemlyb`. 
TODO: Table of performance comparisons, include vs. `alchemlyb`, `pymbar`, `gmx bar` 

## License Information
This project is licensed under either of
  - MIT license ([LICENSE-MIT](LICENSE-MIT))
  - Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))

Copyright © 2026 Helmut Carter
