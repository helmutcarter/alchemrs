# UWHAM Fixtures

This directory is reserved for external-reference fixtures used to validate the
UWHAM estimator against the R `UWHAM` package.

The intended refresh workflow is:

1. Generate a Rust-side reference bundle with `cargo run --example export_uwham_reference -- <output-dir>`.
2. Feed the generated `logQ.csv` and `size.csv` into the R `UWHAM` package using your local helper script.
3. Commit the resulting R-side `ze.csv` and `delta_f.csv` back into this directory.

The repository does not depend on R at test time. The CSVs checked into this
directory should be small, deterministic, and generated from the bundled AMBER
`acetamide_tiny` fixture set.

The current committed R reference uses a reduced five-window subset of that schedule
(`0.2`, `0.3`, `0.4`, `0.5`, and `0.6`). The evaluated-state grid is still the
full parser-derived lambda grid, but only those five windows contribute sampled
observations. This keeps the UWHAM reference test materially cheaper while
preserving a nontrivial solve.
