# alchemrs

`alchemrs` is a Rust-first toolkit for alchemical free energy analysis. It provides a layered workspace of crates for parsing AMBER outputs, preprocessing time series (equilibration trimming and decorrelation), and running common estimators (TI, BAR, MBAR, EXP/DEXP) plus diagnostics like overlap analysis. A CLI (`alchemrs-cli`) offers a simple command-line workflow, and the top-level `alchemrs` crate re-exports a clean public API for library use. Fixtures and tests compare results against established reference implementations (`alchemlyb`) to ensure scientific correctness.

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
- `--output-units <kt|kcal|kj>`: output energy units (default `kt`).
- `--output-format <text|json|csv>`: output format for estimator results (default `text`).
- `--overlap-summary`: include overlap scalar and overlap eigenvalues for BAR/EXP/DEXP/MBAR runs.

### Structured output

JSON output is useful for shell pipelines and downstream tools:

```bash
cargo run -p alchemrs-cli --release -- mbar \
  --temperature 300 \
  --decorrelate \
  --overlap-summary \
  --output-format json \
  /path/to/*/prod.out
```

CSV output is useful for quick ingestion into spreadsheets or tabular tools:

```bash
cargo run -p alchemrs-cli --release -- bar \
  --temperature 300 \
  --output-format csv \
  /path/to/*/prod.out
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

Note: BAR uncertainty is only computed for adjacent windows; non-adjacent entries (e.g. 0→1) are reported as `NaN`, matching alchemlyb.

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
