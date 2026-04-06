# alchemrs

`alchemrs` is a CLI-first tool for alchemical free energy analysis. The package ships an `alchemrs` command-line binary for the main workflow and a Rust library crate for embedding, custom pipelines, and future bindings. The scientific core covers parsing AMBER and GROMACS outputs, preprocessing time series, running TI/BAR/MBAR/EXP/DEXP estimators, and computing diagnostics such as overlap analysis and schedule advice. Fixtures and tests compare results against established reference implementations (`alchemlyb`) to keep the numerical behavior grounded.

Native SVG plotting is available as an optional `plotting` feature.

## CLI

The `alchemrs` binary is the primary entry point.

Commands:

- `advise-schedule`
- `ti`
- `bar`
- `mbar`
- `exp`
- `dexp`

### Build
Pre-built binaries for versioned releases can be found at https://github.com/helmutcarter/alchemrs/releases 
crates.io support will come later in development.

If you want the latest update, or to compile it for your machine:

First, clone this repo
```bash
git clone git@github.com:helmutcarter/alchemrs.git
cd alchemrs
```

then you can either automatically compile and run the binary through `cargo`:
```bash
cargo run --release -- [command] [arguments]
```

or compile using cargo and then call the binary:

```bash
cargo build --release
./target/release/alchemrs [command] [arguments]
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


### Schedule advisor

The CLI also includes a lambda-schedule advisor. For `u_nk` workflows it analyzes adjacent edges and reports whether the current schedule looks healthy, should be monitored, needs more sampling, or likely needs an inserted window. For TI workflows, pass `--input-kind dhdl` to switch the same command over to `dH/dlambda`-based spacing diagnostics.

```bash
alchemrs advise-schedule \
  --temperature 300 \
  --decorrelate \
  --u-nk-observable de \
  --report schedule-report.html \
  --output-format json \
  /path/to/*/prod.out
```

Full CLI usage is documented in [`docs/src/cli.md`](docs/src/cli.md).

## Rust API

The top-level `alchemrs` crate re-exports the common parse, prep, estimator, and analysis entry points used by the CLI and available for direct Rust integration:

```rust
use alchemrs::{
    decorrelate_u_nk_with_observable, extract_u_nk_with_potential, DecorrelationOptions,
    MbarEstimator, MbarOptions,
};

let (u_nk, epot) = extract_u_nk_with_potential("prod.out", 300.0)?;
let u_nk = decorrelate_u_nk_with_observable(&u_nk, &epot, &DecorrelationOptions::default())?;
let fit = MbarEstimator::new(MbarOptions::default()).fit(&[u_nk])?;
let result = fit.result_with_uncertainty()?;
if let Some(labels) = fit.lambda_labels() {
    println!("lambda components = {:?}", labels);
}
```

Use the Rust API when you need embedding, custom orchestration, or tighter control than the CLI exposes. The repo also includes runnable top-level examples:

- `cargo run --example amber_ti -- 300 path/to/lambda0.out path/to/lambda1.out`
- `cargo run --example amber_mbar -- 300 path/to/lambda0.out path/to/lambda1.out path/to/lambda2.out`

Optional plotting helpers can be enabled with:

```bash
cargo build --features plotting --release
```

## CLI Output Example

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

### Advisor details

The CLI also includes a lambda-schedule advisor. For `u_nk` workflows it analyzes adjacent edges and reports whether the current schedule looks healthy, should be monitored, needs more sampling, or likely needs an inserted window. For TI workflows, pass `--input-kind dhdl` to switch the same command over to `dH/dlambda`-based spacing diagnostics.

```bash
alchemrs advise-schedule \
  --temperature 300 \
  --decorrelate \
  --u-nk-observable de \
  --report schedule-report.html \
  --output-format json \
  /path/to/*/prod.out
```

Useful advisor-specific flags:

- `--estimator <mbar|bar>`
- `--input-kind <auto|u-nk|dhdl>`
- `--overlap-min <VALUE>`
- `--block-cv-min <VALUE>`
- `--n-blocks <N>`
- `--no-midpoints`
- `--report <PATH>`

For TI spacing diagnostics:

```bash
alchemrs advise-schedule \
  --temperature 300 \
  --input-kind dhdl \
  --decorrelate \
  --report ti-schedule-report.html \
  --output-format json \
  /path/to/*/prod.out
```

The JSON output includes `sample_counts`, `provenance`, and `suggestions`, plus either `edges` for `u_nk` mode or `windows` and `intervals` for TI mode. `u_nk` edge diagnostics include neighbor-relative metrics, dominant changing components, and a priority score. TI interval diagnostics include `mean_dhdl`, `slope`, `curvature`, `interval_uncertainty`, and block-stability summaries. When `--report` is provided, the CLI also writes a standalone HTML report; the `u_nk` report includes the full SVG lambda-axis and multidimensional component breakdown, while the TI report also includes an integration-method shape gallery showing the applicable TI interpolants on the current lambda grid.

CSV output is useful for quick ingestion into spreadsheets or tabular tools:

```bash
alchemrs bar \
  --temperature 300 \
  --output-format csv \
  /path/to/*/prod.out
```

CSV columns include estimator parameters after the result fields:

```text
delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues,estimator,temperature_k,decorrelate,remove_burnin,auto_equilibrate,fast,conservative,nskip,u_nk_observable,ti_method,ti_method_reason,lambda_components,windows,samples_in,samples_after_burnin,samples_kept
```

### TI

```bash
alchemrs ti \
  --temperature 300 \
  --method trapezoidal \
  --remove-burnin 125 \
  /path/to/*/prod.out
```

`ti` supports `--method <trapezoidal|simpson|cubic-spline|pchip|akima|gaussian-quadrature>`.
For more on picking an integration method, see the docs.

To let the CLI choose a method automatically and record why it was selected:

```bash
alchemrs ti \
  --temperature 300 \
  --method auto \
  --output-format json \
  /path/to/*/prod.out
```

TI outputs now include `ti_method` and, for `--method auto`, `ti_method_reason` in provenance.

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

## Documentation

Documentation is online at https://helmutcarter.github.io/alchemrs/. The source files are in [`docs/`](docs/) and can be viewed locally:

```bash
cargo install mdbook
mdbook serve docs
```

## Performance

Initial internal benchmarks on the bundled fixtures show substantial speedups over
`alchemlyb` for several workflows. Ongoing profiling notes and experiment results live in
[`optimizations.md`](optimizations.md).

## License Information

This project is licensed under either of

- MIT license ([LICENSE-MIT](LICENSE-MIT))
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))

Copyright (c) 2026 Helmut Carter
