# Repository Layout

The repository is centered on a single Cargo package, `alchemrs`.

That package contains:

- the `alchemrs` library crate
- the `alchemrs` command-line binary

## `alchemrs` library layout

The `alchemrs` crate is organized into modules rather than separate library crates:

- `core`: domain types and shared error model
- `parse`: engine-specific parsers; today this means AMBER parsing
- `prep`: time-series preparation such as duplicate cleanup, sorting, equilibration detection, and decorrelation
- `estimators`: TI, BAR, MBAR, EXP, and DEXP implementations
- `analysis`: overlap diagnostics built on top of MBAR log weights

The crate root also re-exports the most common types and functions directly for a flatter API.

## CLI binary

The `alchemrs` binary wraps parsing, preprocessing, estimation, and output formatting into a command-line workflow. The scientific logic lives in the library; the CLI is a thin consumer of that API.

## Layering

The intended dependency direction inside the library is:

```text
parse -> core
prep -> core
estimators -> core
analysis -> core + estimators
```

The binary depends on the library through the same package source tree.
