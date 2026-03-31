# CLI Guide

The CLI is the `alchemrs` binary in the same package as the `alchemrs` library crate.

Commands:

- `ti`
- `bar`
- `mbar`
- `exp`
- `dexp`

## Build

```bash
cargo build --release
```

## Common workflow

Typical CLI usage is:

1. pass one simulation output per lambda window
2. optionally trim burn-in
3. optionally run auto-equilibration
4. optionally decorrelate
5. fit the requested estimator
6. emit text, JSON, or CSV

## Shared flags

- `--temperature <K>`
- `--remove-burnin <N>`
- `--auto-equilibrate`
- `--decorrelate`
- `--fast`
- `--conservative[=true|false]`
- `--nskip <N>`
- `--output-units <kt|kcal|kj>`
- `--output-format <text|json|csv>`
- `--output <PATH>`
- `--parallel`

For `bar`, `mbar`, `exp`, and `dexp`:

- `--u-nk-observable <de|all|epot>`

For overlap-aware commands:

- `--overlap-summary`

## Observable selection

For `u_nk`-based estimators:

- `de` is the default and matches the `alchemlyb`-style adjacent-state `dE` observable
- `all` sums the full `u_nk` row
- `epot` uses an engine-provided potential-energy observable (`EPtot` for AMBER, `Potential Energy` for GROMACS `dhdl.xvg`)

Use `epot` when:

- you want the CLI's external-observable path
- the `u_nk` matrix contains positive infinity values that make `de` invalid
- you need a preprocessing observable that does not rely on one-dimensional adjacent-state `DE` semantics

For multidimensional `u_nk` schedules:

- `bar`, `mbar`, `exp`, and `dexp` accept them
- CLI output renders `from_lambda` / `to_lambda` as JSON arrays or bracketed text/CSV values when the state has multiple lambda components
- parser-derived lambda component labels are logged in text, JSON, and CSV provenance when available
- `de` remains one-dimensional, so use `all` or `epot` when preprocessing multidimensional schedules

## TI-specific behavior

TI uses `dH/dlambda`, not `u_nk`.

For GROMACS files with multidimensional lambda schedules, the CLI currently rejects TI input because the parser cannot collapse multiple `dH/dlambda` components into a single scalar series safely.

The CLI accepts `--u-nk-observable` on `ti` only to provide a more helpful error message. If supplied, the command fails with a domain-specific explanation instead of a generic unknown-flag parse error.

## Command examples

### TI

```bash
cargo run --release -- ti \
  --temperature 300 \
  --method trapezoidal \
  --decorrelate \
  path/to/*/prod.out
```

### BAR

```bash
cargo run --release -- bar \
  --temperature 300 \
  --decorrelate \
  --u-nk-observable de \
  --overlap-summary \
  path/to/*/prod.out
```

### MBAR with EPtot fallback

```bash
cargo run --release -- mbar \
  --temperature 300 \
  --auto-equilibrate \
  --decorrelate \
  --u-nk-observable epot \
  --output-format json \
  path/to/*/prod.out
```

### EXP / DEXP

```bash
cargo run --release -- exp --temperature 300 path/to/*/prod.out

cargo run --release -- dexp --temperature 300 path/to/*/prod.out
```

## Effective settings

When `--auto-equilibrate` is enabled, the CLI reports the effective preprocessing policy in provenance:

- `fast = true`
- `conservative = false`

That override is deliberate and mirrors the `alchemlyb` / `pymbar` equilibration-detection workflow.
