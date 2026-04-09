# Introduction

`alchemrs` is a CLI-first tool for alchemical free energy analysis.

The project is organized as one package containing a library crate and a single CLI binary. Inside the library, the code is split into focused modules for:

- parsing supported AMBER and GROMACS outputs into typed Rust data structures
- preprocessing time series by trimming and decorrelation
- estimating free energies with TI, BAR, MBAR, IEXP, DEXP, and NES
- computing overlap diagnostics
- computing schedule advice for `u_nk`, TI, and NES workflows
- exposing the workflow through a CLI and a clean top-level Rust API

The intended user journey is:

1. Start with the `alchemrs` CLI for the standard workflows.
2. Drop to the Rust API when you need embedding, scripting, tighter control, or integration into a larger application.
3. Reuse the same library surface later from bindings rather than maintaining a separate implementation.

At the moment, the practical workflow is:

1. Parse supported engine outputs into `DhdlSeries`, `UNkMatrix`, or `SwitchingTrajectory`.
2. Trim and decorrelate those data if needed.
3. Apply an estimator.
4. Inspect uncertainties and overlap diagnostics.

This can be done from the CLI or the library API, but the CLI is the default front door.

This book documents the current codebase rather than a future roadmap.
