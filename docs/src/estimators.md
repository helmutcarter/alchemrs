# Estimators

Estimator implementations live in `alchemrs-estimators`.

## Shared assumptions

For `UNkMatrix`-based estimators, the current implementation assumes:

- one-dimensional lambda states
- consistent evaluated-state grids across windows
- `sampled_state` is set and present in `evaluated_states`
- state temperatures are consistent across all windows

If those invariants are violated, the estimators fail early.

## TI

Types:

- `TiEstimator`
- `TiOptions`
- `IntegrationMethod::{Trapezoidal, Simpson}`

Input:

- slice of `DhdlSeries`

Behavior:

- sorts windows by lambda
- computes the mean `dH/dlambda` for each window
- integrates across lambda using trapezoidal or Simpson integration

Current uncertainty behavior:

- trapezoidal mode reports an uncertainty estimate
- Simpson mode currently returns no uncertainty

## BAR

Types:

- `BarEstimator`
- `BarOptions`
- `BarMethod::{FalsePosition, SelfConsistentIteration, Bisection}`

Input:

- slice of `UNkMatrix`

Behavior:

- expects windows mapped to sampled states on one evaluated-state grid
- computes work values between adjacent states
- fills a full pairwise `DeltaFMatrix`

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
- `compute_uncertainty`
- `parallel`

Input:

- slice of `UNkMatrix`

Behavior:

- combines windows into a shared reduced-energy representation
- solves for free energies iteratively
- returns a full `DeltaFMatrix`

The analysis crate uses MBAR-derived log weights to compute overlap diagnostics.

## EXP and DEXP

Types:

- `ExpEstimator`
- `ExpOptions`

Input:

- slice of `UNkMatrix`

Behavior:

- computes exponential averaging from each sampled-state window to all evaluated states
- returns a full pairwise `DeltaFMatrix`

CLI direction conventions:

- `exp` reports the forward direction
- `dexp` reports the reverse direction

`DEXP` is not a separate estimator type in the library; it is a CLI convention over the same estimator results.

## Parallel execution

Some estimators support a `parallel` option and use Rayon-backed parallel loops internally.

Parallelism affects performance, not the conceptual API.
