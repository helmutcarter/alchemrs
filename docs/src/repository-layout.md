# Repository Layout

The repository is organized as a Cargo workspace with two Rust packages:

- `alchemrs`: the core library crate and CLI binary
- `alchemrs-py`: the Python binding crate built as a `cdylib`

The surrounding top-level directories then hold tests, examples, fixtures,
documentation, and Python-side packaging/test code.

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

## `alchemrs-py` package

`alchemrs-py` is a separate Rust package in the same workspace. It provides the
native Python extension layer and depends on the core `alchemrs` crate, but the
dependency does not go the other way.

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

The CLI depends on the core library crate, and `alchemrs-py` depends on the
same core library crate through the workspace.
