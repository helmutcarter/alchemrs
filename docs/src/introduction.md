# Introduction

`alchemrs` is a Rust-first toolkit for alchemical free energy analysis.

The workspace is split into focused crates for:

- parsing AMBER output into typed Rust data structures
- preprocessing time series by trimming and decorrelation
- estimating free energies with TI, BAR, MBAR, EXP, and DEXP
- computing overlap diagnostics
- exposing the workflow through a CLI and a clean top-level Rust API

The project is intentionally Rust-native. The goal is not to wrap a Python implementation; the Rust crates are the canonical implementation surface.

At the moment, the practical workflow is:

1. Parse AMBER outputs into `DhdlSeries` or `UNkMatrix`.
2. Trim and decorrelate those data if needed.
3. Fit an estimator.
4. Inspect uncertainties and overlap diagnostics.
5. Use the same functionality from the CLI or the library API.

This book documents the current codebase rather than a future roadmap. It focuses on how the library and CLI behave today, including the scientifically important details that are easy to miss if you only read the README.
