# alchemrs

`alchemrs` is a Rust-first toolkit for alchemical free energy analysis. It provides a layered workspace of crates for parsing AMBER outputs, preprocessing time series (equilibration trimming and decorrelation), and running common estimators (TI, BAR, MBAR, EXP/DEXP) plus diagnostics like overlap analysis. A CLI (`alchemrs-cli`) offers a simple command-line workflow, and the top-level `alchemrs` crate re-exports a clean public API for library use. Fixtures and tests compare results against established reference implementations (`alchemlyb`) to ensure scientific correctness.

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
```

## CLI

`alchemrs-cli` provides a command-line workflow for AMBER output files.

### Build

```bash
cargo build -p alchemrs-cli --release
```

### Common flags

- `--remove-burnin <N>`: skip the first `N` samples in each window before analysis.
- `--auto-equilibrate`: use `pymbar`-style equilibration detection as described [here](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00784).
- `--decorrelate`: subsample to reduce correlation before estimating.
  - `ti` uses the parsed `dH/dlambda` series.
  - `bar`, `exp`, `dexp`, and `mbar` use `EPtot` parsed from the AMBER output when decorrelation is requested.
- `--output-units <kt|kcal|kj>`: output energy units (default `kt`).
- `--output-format <text|json|csv>`: output format for estimator results (default `text`).
- `--output <PATH>`: write the formatted result to a file instead of stdout.
- `--overlap-summary`: include overlap scalar and overlap eigenvalues for BAR/EXP/DEXP/MBAR runs.

### Structured output

JSON output is useful for shell pipelines and downstream tools:

```bash
cargo run -p alchemrs-cli --release -- mbar \
  --temperature 300 \
  --decorrelate \
  --overlap-summary \
  --output-format json \
  --output results.json \
  /path/to/*/prod.out
```

For `u_nk`-based estimators (`bar`, `exp`, `dexp`, `mbar`), `--decorrelate` uses `EPtot` as the finite decorrelation observable and then applies those retained indices back to the parsed `u_nk` samples.

Example JSON payload:

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
    "nskip": 1
  }
}
```

CSV output is useful for quick ingestion into spreadsheets or tabular tools:

```bash
cargo run -p alchemrs-cli --release -- bar \
  --temperature 300 \
  --output-format csv \
  /path/to/*/prod.out
```

CSV columns now include estimator provenance after the result fields:

```text
delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues,estimator,temperature_k,decorrelate,remove_burnin,auto_equilibrate,fast,conservative,nskip
```

### TI (trapezoidal)

```bash
cargo run -p alchemrs-cli --release -- ti \
  --temperature 300 \
  --method trapezoidal \
  --remove-burnin 125 \
  /path/to/*/prod.out
```

### BAR

```bash
cargo run -p alchemrs-cli --release -- bar \
  --temperature 300 \
  --method false-position \
  --decorrelate \
  /path/to/*/prod.out
```

Note: BAR uncertainties are only computed for adjacent windows; non-adjacent state pairs (for example, `0->1`) are reported as `NaN`, matching alchemlyb.

### MBAR

```bash
cargo run -p alchemrs-cli --release -- mbar \
  --temperature 300 \
  --decorrelate \
  --max-iterations 10000 \
  --tolerance 1e-7 \
  /path/to/*/prod.out
```

### EXP / DEXP

```bash
cargo run -p alchemrs-cli --release -- exp \
  --temperature 300 \
  /path/to/*/prod.out

cargo run -p alchemrs-cli --release -- dexp \
  --temperature 300 \
  /path/to/*/prod.out
```

EXP reports FEP results in the forward direction, DEXP reports FEP results in the reverse direction.

## Performance 
Initial results demonstrate 6x performance improvement over `alchemlyb`. 
TODO: Table of performance comparisons, include vs. `alchemlyb`, `pymbar`, `gmx bar` 
