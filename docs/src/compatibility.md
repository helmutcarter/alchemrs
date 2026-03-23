# Compatibility and Limitations

This section describes the current scope of the codebase and where it intentionally matches external tools.

## Current compatibility targets

The repository includes fixtures and tests intended to compare behavior against:

- `alchemlyb`
- `pymbar`

Important matches include:

- `u_nk` `DE` decorrelation semantics
- `pymbar`-style statistical inefficiency estimation
- `pymbar`-style automated equilibration detection
- BAR reporting `NaN` uncertainties for non-adjacent pairs

Important intentional difference:

- `alchemrs` reports equilibration `Neff_max` using the actual suffix length `N - t0`, so this scalar can differ by 1 from `pymbar` / `alchemlyb`

## AMBER-first scope

The implemented parser surface today is AMBER-focused.

The architecture leaves room for other engines later, but the documented and tested workflow in this repository currently centers on AMBER outputs.

## One-dimensional lambda assumption

The data-model types can store more general state metadata, but the estimator layer currently assumes one-dimensional lambda states.

If you need multidimensional alchemical states, the current estimators are not yet the right abstraction surface.

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

- plotting in the Rust library
- production-grade multi-engine parser support
- Python bindings in the current package
- versioned public API guarantees

Those may come later, but they are not part of the implemented surface documented in this book.
