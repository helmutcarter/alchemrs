# Workspace Layout

The workspace is organized into a top-level crate plus several focused internal crates.

## Top-level crate

The root `alchemrs` crate re-exports the most commonly used parse, preprocessing, estimator, and analysis entry points. For many users, this is the main library surface.

It re-exports:

- core data types such as `StatePoint`, `DhdlSeries`, and `UNkMatrix`
- parse functions such as `extract_dhdl`, `extract_u_nk`, and `extract_u_nk_with_potential`
- preprocessing functions such as `decorrelate_dhdl`, `decorrelate_u_nk`, and `detect_equilibration_u_nk`
- estimators such as `TiEstimator`, `BarEstimator`, `MbarEstimator`, and `ExpEstimator`
- overlap helpers such as `overlap_matrix`, `overlap_eigenvalues`, and `overlap_scalar`

## Internal crates

### `alchemrs-core`

Defines the canonical domain types and shared error model.

### `alchemrs-parse`

Parses engine-specific output into the core types. The implemented parser today is the AMBER parser.

### `alchemrs-prep`

Handles time-series preparation:

- duplicate-time cleanup
- sorting
- time slicing
- equilibration detection
- decorrelation / subsampling

### `alchemrs-estimators`

Implements the free energy estimators:

- TI
- BAR
- MBAR
- EXP
- DEXP

### `alchemrs-analysis`

Provides overlap diagnostics built on top of MBAR log weights.

### `alchemrs-cli`

Wraps the parser, preprocessing, estimators, and output formatting into a command-line workflow.

## Layering

The intended dependency direction is:

```text
parse -> core
prep -> core
estimators -> core
analysis -> core + estimators
cli -> parse + prep + estimators + analysis
top-level crate -> re-exports
```

This separation matters because it keeps parsing, preprocessing, estimation, and presentation concerns distinct. The CLI is not where the scientific logic lives; it composes the library crates.
