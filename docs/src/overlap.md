# Overlap Diagnostics

Overlap diagnostics live in the `alchemrs::analysis` module.

Available functions:

- `overlap_matrix`
- `overlap_eigenvalues`
- `overlap_scalar`
- `ti_convergence`
- `bar_convergence`
- `mbar_convergence`
- `exp_convergence`
- `dexp_convergence`

## `overlap_matrix`

This function:

1. computes MBAR log weights from the supplied windows
2. constructs the overlap matrix from those weights and state sample counts
3. returns an `OverlapMatrix`

When the input windows carry named lambda dimensions, the returned matrix preserves them via `overlap.lambda_labels()`.

Input:

- slice of `UNkMatrix`
- optional `MbarOptions`

## `overlap_eigenvalues`

Computes the eigenvalues of the overlap matrix and returns them sorted in descending order.

The implementation rejects eigenvalues with significant imaginary components.

## `overlap_scalar`

Computes:

```text
1 - second_largest_eigenvalue
```

This is the scalar overlap diagnostic currently surfaced by the CLI.

## CLI usage

`bar`, `mbar`, `exp`, and `dexp` can include overlap diagnostics with:

```bash
--overlap-summary
```

The CLI then reports:

- overlap scalar
- overlap eigenvalues

in text, JSON, or CSV output.

## Convergence series

The same module also exposes cumulative convergence helpers for plot-ready library output.

Each function returns a `Vec<ConvergencePoint>`:

- `ti_convergence`
- `bar_convergence`
- `mbar_convergence`
- `exp_convergence`
- `dexp_convergence`

Each point includes:

- the number of windows included in that cumulative estimate
- the estimated free energy
- optional uncertainty
- the endpoint states for that estimate
- optional `lambda_labels()` when they are available from the parsed windows
