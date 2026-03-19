# alchemrs

Rust library for alchemical free energy analysis.

## CLI

`alchemrs-cli` provides a plain-text workflow for AMBER output files.

### Build

```bash
cargo build -p alchemrs-cli --release
```

### Common flags

- `--remove-burnin <N>`: skip the first `N` samples in each window before analysis.
- `--auto-equilibrate`: use pymbar-style equilibration detection (overrides `--fast`/`--conservative`).
- `--decorrelate`: subsample to reduce correlation before estimating.
- `--output-units <kt|kcal|kj>`: output energy units (default `kt`).

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
