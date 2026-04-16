# Introduction

`alchemrs` is a CLI-first tool for alchemical free energy analysis.

The project is organized as a Cargo workspace with a core Rust package, a CLI
binary, and a Python extension package. Inside the library, the code is split
into focused modules for:

- parsing supported AMBER and GROMACS outputs into typed Rust data structures
- preprocessing time series by trimming and decorrelation
- estimating free energies with TI, BAR, MBAR, IEXP, DEXP, NES, UWHAM, and ATM
- computing overlap diagnostics
- computing schedule advice for `u_nk`, TI, and NES workflows
- exposing the workflow through a CLI, a clean top-level Rust API, and Python bindings

The intended user journey is:

1. Start with the `alchemrs` CLI for the standard workflows.
2. Drop to the Rust API when you need embedding, scripting, tighter control, or integration into a larger application.
3. Use the existing Python package when you want the same core implementation from Python or OpenMM-driven workflows.

At the moment, the practical workflow is:

1. Parse supported engine outputs into `DhdlSeries`, `UNkMatrix`, or `SwitchingTrajectory`, or construct `AtmSampleSet` / `AtmLogQMatrix` directly for ATM workflows.
2. Trim and decorrelate those data if needed.
3. Apply an estimator.
4. Inspect uncertainties, overlap diagnostics, or schedule advice.

This can be done from the CLI or the library API, but the CLI is the default front door.

This book documents the current codebase rather than a future roadmap.
