# Repository Layout

The repository is centered on a single Cargo package, `alchemrs`.

That package contains:

- the `alchemrs` library crate
- the `alchemrs` command-line binary

## `alchemrs` library layout

The `alchemrs` crate is organized into modules rather than separate library crates:

- `data`: domain types and scientific data-model structures
- `error`: shared crate-level error handling
- `parse`: engine-specific parsers for AMBER outputs and GROMACS `dhdl.xvg` files
- `prep`: time-series preparation such as duplicate cleanup, sorting, equilibration detection, and decorrelation
- `estimators`: TI, BAR, MBAR, IEXP, DEXP, and NES implementations
- `analysis`: overlap diagnostics, schedule advisors, and plot-ready convergence series
- `plot`: optional SVG rendering helpers behind the `plotting` feature

The crate root also re-exports the most common types and functions directly for a flatter API.

## CLI binary

The `alchemrs` binary wraps parsing, preprocessing, estimation, and output formatting into a command-line workflow. The scientific logic lives in the library; the CLI is a thin consumer of that API.

## Layering

The intended dependency direction inside the library is:

```text
parse -> data + error
prep -> data + error
estimators -> data + error
analysis -> data + error + estimators
```

The binary depends on the library through the same package source tree.
