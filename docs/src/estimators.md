# Estimators

Estimator implementations live in the `alchemrs::estimators` module.

## Shared assumptions

For `UNkMatrix`-based estimators, the current implementation assumes:

- consistent lambda-vector dimensionality across windows
- consistent evaluated-state grids across windows
- `sampled_state` is set and present in `evaluated_states`
- state temperatures are consistent across all windows
- `lambda_labels`, when present, are consistent across windows

If those invariants are violated, the estimators fail early.

## TI

Types:

- `TiEstimator`
- `TiFit`
- `TiOptions`
- `IntegrationMethod::{Trapezoidal, Simpson, GaussianQuadrature}`

Input:

- slice of `DhdlSeries`

Behavior:

- sorts windows by lambda
- computes the mean `dH/dlambda` for each window
- integrates across lambda using trapezoidal, Simpson, or Gaussian quadrature
- `fit(...)` returns `TiFit`
- `result()` materializes the `FreeEnergyEstimate`

Current uncertainty behavior:

- trapezoidal mode reports an uncertainty estimate
- Gaussian quadrature reports an uncertainty estimate from the quadrature weights and per-window SEMs
- Simpson mode currently returns no uncertainty

Gaussian quadrature behavior:

- only supports one-dimensional lambda schedules
- only supports Gauss-Legendre schedules with `1..=16` sampled windows
- validates that the sampled lambdas match the supported quadrature nodes
- reports the free energy over the full `[0, 1]` interval, so `result().from_state()` is `lambda=0` and `result().to_state()` is `lambda=1` even though the sampled windows are interior quadrature nodes

## BAR

Types:

- `BarEstimator`
- `BarFit`
- `BarOptions`
- `BarMethod::{FalsePosition, SelfConsistentIteration, Bisection}`

Input:

- slice of `UNkMatrix`

Behavior:

- expects windows mapped to sampled states on one evaluated-state grid
- computes work values between adjacent states
- `fit(...)` returns `BarFit`
- `result()` materializes the full pairwise `DeltaFMatrix`
- preserves `lambda_labels()` from the input windows when available

Important limitation:

- uncertainty is only computed for adjacent windows
- non-adjacent uncertainties are reported as `NaN`

This matches the current CLI documentation and the reference-comparison tests in the repository.

## MBAR

Types:

- `MbarEstimator`
- `MbarOptions`

Important options:

- `max_iterations`
- `tolerance`
- `initial_f_k`
- `parallel`

Input:

- slice of `UNkMatrix`

Behavior:

- combines windows into a shared reduced-energy representation
- solves for free energies iteratively
- returns an `MbarFit` that can derive a `DeltaFMatrix`, overlap diagnostics, and MBAR log weights without resolving the same window set
- preserves `lambda_labels()` from the input windows when available

Typical usage:

- `fit(...).result()`
- `fit(...).result_with_uncertainty()`
- `fit(...).overlap_matrix()`
- `fit(...).overlap_scalar()`

The analysis layer uses MBAR-derived log weights to compute overlap diagnostics.

## EXP and DEXP

Types:

- `ExpEstimator`
- `ExpFit`
- `ExpOptions`

Input:

- slice of `UNkMatrix`

Behavior:

- computes exponential averaging from each sampled-state window to all evaluated states
- `fit(...)` returns `ExpFit`
- `result()` materializes the full pairwise `DeltaFMatrix`
- `result_with_uncertainty()` includes uncertainty estimates
- preserves `lambda_labels()` from the input windows when available

CLI direction conventions:

- `exp` reports the forward direction
- `dexp` reports the reverse direction

`DEXP` is not a separate estimator type in the library; it is a CLI convention over the same estimator results.

## Parallel execution

Some estimators support a `parallel` option and use Rayon-backed parallel loops internally.

Parallelism affects performance, not the conceptual API.

## Common API shape

The estimator interfaces follow one shared pattern:

- `fit(...) -> Fit`
- `estimate(...) -> primary result`
- `Fit::result() -> primary result`

Estimators with an explicit uncertainty split also expose `estimate_with_uncertainty()` and `Fit::result_with_uncertainty()`.
