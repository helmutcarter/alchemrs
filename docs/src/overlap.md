# Overlap Diagnostics

Overlap diagnostics live in `alchemrs-analysis`.

Available functions:

- `overlap_matrix`
- `overlap_eigenvalues`
- `overlap_scalar`

## `overlap_matrix`

This function:

1. computes MBAR log weights from the supplied windows
2. constructs the overlap matrix from those weights and state sample counts
3. returns an `OverlapMatrix`

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
