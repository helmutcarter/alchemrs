# Compatibility and Limitations

This section describes the current scope of the codebase and where it intentionally matches external tools.

## Current compatibility targets

The repository includes fixtures and tests intended to compare behavior against:

- `alchemlyb`
- `pymbar`
- the R `UWHAM` package

Important matches include:

- `u_nk` `DE` decorrelation semantics
- `pymbar`-style statistical inefficiency estimation
- `pymbar`-style automated equilibration detection
- UWHAM free-energy estimates and pooled-input exports validated against committed R reference CSVs
- AMBER NES Jarzynski estimates and preprocessing semantics validated against external reference scripts

Important intentional difference:

- `alchemrs` reports equilibration `Neff_max` using the actual suffix length `N - t0`, so this scalar can differ by 1 from `pymbar` / `alchemlyb`

## Current parser scope

The implemented parser surface supports AMBER outputs, including AMBER nonequilibrium switching outputs for NES, and GROMACS `dhdl.xvg` outputs.

AMBER still has the broader and more battle-tested workflow coverage, but the documented and tested parsing workflow now includes both engines. Multidimensional GROMACS schedules are supported for `u_nk`-based analysis; multidimensional TI is still unsupported.

ATM and OpenMM workflows currently enter through the typed data model and Python helper APIs rather than through a dedicated file parser.

## One-dimensional lambda assumption

The data-model types can store more general state metadata, and `UNkMatrix`-based estimators now support multidimensional lambda states when windows share a consistent evaluated-state grid.

The remaining one-dimensional restriction is TI: `DhdlSeries` is still scalar, so multidimensional `dH/dlambda` workflows are not yet supported.

NES trajectories are also currently one-dimensional in the advisor path because the retained `lambda_path` / `dvdl_path` diagnostics assume one switching coordinate.

## Non-finite `u_nk`

Positive infinity values are allowed in `UNkMatrix`, but not every preprocessing path can use them safely.

Practical consequence:

- `de` and `all` observable choices can fail during decorrelation if the derived scalar series is non-finite
- `epot` can remain usable if the external observable is finite

## Documentation scope

This book documents the current implementation, not a frozen public-stability promise. The project is still early enough that:

- APIs may still evolve
- defaults may still be refined
- more formal hosted documentation infrastructure may still change

## What is not covered here

This repository does not currently provide:

- a large general-purpose plotting/reporting surface beyond the current SVG helpers and advisor HTML reports
- production-grade parser coverage for every major MD engine
- a dedicated CLI surface for ATM or direct UWHAM analysis
- a polished published Python distribution workflow beyond the current repo-local `maturin` setup
- versioned public API guarantees

Those may come later, but they are not part of the implemented surface documented in this book.
