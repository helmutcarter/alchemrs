# CLI Guide

The CLI crate is `alchemrs-cli`.

Commands:

- `ti`
- `bar`
- `mbar`
- `exp`
- `dexp`

## Build

```bash
cargo build -p alchemrs-cli --release
```

## Common workflow

Typical CLI usage is:

1. pass one AMBER output per lambda window
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
- `--conservative`
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
- `epot` uses `EPtot` parsed from the AMBER output

Use `epot` when:

- you want the CLI’s external-observable path
- the `u_nk` matrix contains positive infinity values that make `de` invalid

## TI-specific behavior

TI uses `dH/dlambda`, not `u_nk`.

The CLI accepts `--u-nk-observable` on `ti` only to provide a more helpful error message. If supplied, the command fails with a domain-specific explanation instead of a generic unknown-flag parse error.

## Command examples

### TI

```bash
cargo run -p alchemrs-cli --release -- ti \
  --temperature 300 \
  --method trapezoidal \
  --decorrelate \
  path/to/*/prod.out
```

### BAR

```bash
cargo run -p alchemrs-cli --release -- bar \
  --temperature 300 \
  --decorrelate \
  --u-nk-observable de \
  --overlap-summary \
  path/to/*/prod.out
```

### MBAR with EPtot fallback

```bash
cargo run -p alchemrs-cli --release -- mbar \
  --temperature 300 \
  --auto-equilibrate \
  --decorrelate \
  --u-nk-observable epot \
  --output-format json \
  path/to/*/prod.out
```

### EXP / DEXP

```bash
cargo run -p alchemrs-cli --release -- exp --temperature 300 path/to/*/prod.out

cargo run -p alchemrs-cli --release -- dexp --temperature 300 path/to/*/prod.out
```

## Effective settings

When `--auto-equilibrate` is enabled, the CLI reports the effective preprocessing policy in provenance:

- `fast = true`
- `conservative = false`

That override is deliberate and mirrors the `alchemlyb` / `pymbar` equilibration-detection workflow.
