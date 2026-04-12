# UWHAM Fixtures

This directory is reserved for external-reference fixtures used to validate the
UWHAM estimator against the R `UWHAM` package.

The intended workflow is:

1. Generate pooled `logQ.csv` and `size.csv` with `cargo run --example export_uwham_reference -- <output-dir>`.
2. Run `Rscript scripts/uwham_reference.R <output-dir>` on those CSVs.
3. Commit the resulting `ze.csv` and `delta_f.csv` back into this directory.

The repository does not depend on R at test time. The CSVs checked into this
directory should be small, deterministic, and generated from the bundled AMBER
`acetamide_tiny` fixture set.
