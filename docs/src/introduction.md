# Introduction

`alchemrs` is a Rust-first toolkit for alchemical free energy analysis.

The project is organized as one main library crate plus a separate CLI package. Inside the library, the code is split into focused modules for:

- parsing AMBER output into typed Rust data structures
- preprocessing time series by trimming and decorrelation
- estimating free energies with TI, BAR, MBAR, EXP, and DEXP
- computing overlap diagnostics
- exposing the workflow through a CLI and a clean top-level Rust API

The project is intentionally Rust-native. The `alchemrs` library crate is the canonical implementation surface, and a CLI is included for convenience. Python bindings will come later.

At the moment, the practical workflow is:

1. Parse AMBER outputs into `DhdlSeries` or `UNkMatrix`.
2. Trim and decorrelate those data if needed.
3. Apply an estimator.
4. Inspect uncertainties and overlap diagnostics.

This can be done from the CLI or the library API.

This book documents the current codebase rather than a future roadmap.
