# Estimators

Estimator implementations live in the `alchemrs::estimators` module.

The CLI currently exposes TI, BAR, MBAR, IEXP, DEXP, and NES. The library and
Python bindings additionally expose UWHAM and ATM.

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
- `IntegrationMethod::{Trapezoidal, Simpson, CubicSpline, Pchip, Akima, GaussianQuadrature}`
- `recommend_ti_method(...)`
- `TiMethodRecommendation`
- `TiMethodAssessment`
- `TiMethodRecommendationOptions`

Input:

- slice of `DhdlSeries`

Behavior:

- sorts windows by lambda
- computes the mean `dH/dlambda` for each window
- integrates across lambda using trapezoidal, Simpson, a natural cubic spline, PCHIP, Akima interpolation, or Gaussian quadrature
- `fit(...)` returns `TiFit`
- `result()` materializes the `FreeEnergyEstimate`

Current uncertainty behavior:

- trapezoidal mode reports an uncertainty estimate
- Simpson mode reports an uncertainty estimate by propagating the Simpson-rule weights through the per-window SEMs
- cubic-spline mode reports an uncertainty estimate by propagating the implied spline-integration weights through the per-window SEMs
- PCHIP and Akima report an uncertainty estimate by numerically differentiating the integrated result with respect to the window means and propagating the per-window SEMs through that Jacobian
- Gaussian quadrature reports an uncertainty estimate from the quadrature weights and per-window SEMs

Interpolation references:

- `Pchip` follows the piecewise-cubic Hermite shape-preserving literature of Fritsch and Carlson, *Monotone Piecewise Cubic Interpolation*, SIAM Journal on Numerical Analysis 17(2), 1980, DOI: <https://doi.org/10.1137/0717021>, and Fritsch and Butland, *A Method for Constructing Local Monotone Piecewise Cubic Interpolants*, SIAM Journal on Scientific and Statistical Computing 5(2), 1984, DOI: <https://doi.org/10.1137/0905021>
- `Akima` follows Akima, *A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures*, Journal of the ACM 17(4), 1970, DOI: <https://doi.org/10.1145/321607.321609>

Gaussian quadrature behavior:

- only supports one-dimensional lambda schedules
- only supports Gauss-Legendre schedules with `1..=16` sampled windows
- validates that the sampled lambdas match the supported quadrature nodes
- reports the free energy over the full `[0, 1]` interval, so `result().from_state()` is `lambda=0` and `result().to_state()` is `lambda=1` even though the sampled windows are interior quadrature nodes

Method recommendation:

- `recommend_ti_method(...)` lives in `alchemrs::analysis`
- evaluates all currently supported TI methods on the same preprocessed `DhdlSeries`
- combines method eligibility, TI schedule diagnostics, monotonicity, and cross-method spread
- returns the recommended `IntegrationMethod` plus per-method assessments and a human-readable rationale

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

Uncertainty behavior:

- adjacent BAR uncertainties come from the BAR implicit-root equation via a delta-method linearization of the Fermi weights
- cumulative non-adjacent uncertainties include the neighboring-edge covariance induced by the shared intermediate window rather than assuming adjacent BAR edges are independent

## MBAR

Types:

- `MbarEstimator`
- `MbarOptions`

Important options:

- `max_iterations`
- `tolerance`
- `initial_f_k`
- `parallel`
- `solver`

Default solver:

- `MbarSolver::Lbfgs`
- falls back to fixed-point automatically when the evaluated grid contains states with zero sampled counts, because the LBFGS objective requires positive sampled counts
- use `MbarSolver::FixedPoint` when reproducing older fixed-point results or comparing solver behavior

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

## UWHAM

Types:

- `UwhamEstimator`
- `UwhamFit`
- `UwhamOptions`

Important options:

- `max_iterations`
- `tolerance`
- `parallel`

Input:

- slice of `UNkMatrix`

Behavior:

- pools all windows into log unnormalized densities (`log_q`)
- accepts positive infinity in `u_nk` by mapping it to negative infinity in `log_q`
- solves the UWHAM equations on the pooled observations
- `fit(...)` returns `UwhamFit`
- `result()` materializes a `DeltaFMatrix`
- `result_with_uncertainty()` uses the analytical UWHAM covariance
- `weights()`, `covariance()`, and `variances()` expose fitted internals without rerunning the solve
- `write_reference_inputs(...)` exports `logQ.csv`, `size.csv`, `delta_f.csv`, `ze.csv`, and metadata for cross-checking against external UWHAM implementations

## ATM

Types:

- `AtmEstimator`
- `AtmFit`
- `AtmOptions`
- `AtmUncertaintyMethod::{Analytical, Bootstrap { .. }, None}`
- `AtmBindingEstimate`

Input:

- `fit(...)`: `AtmLogQMatrix`
- `fit_leg(...)`, `estimate_leg(...)`, `estimate_binding(...)`: `AtmSampleSet`

Behavior:

- converts ATM schedule/sample data into pooled `log_q`
- delegates the core statistical solve to the UWHAM implementation
- `result()` on `AtmFit` returns a `DeltaFMatrix` across all ATM schedule states
- `leg_result()` and `estimate_leg(...)` collapse that fit to one endpoint-to-endpoint `FreeEnergyEstimate`
- `estimate_binding(...)`, `estimate_rbfe(...)`, and `estimate_abfe(...)` combine two opposite-direction legs and propagate leg uncertainties in quadrature
- preserves ATM-derived lambda component labels when available

Uncertainty behavior:

- `Analytical` uses UWHAM covariance on the fitted leg
- `Bootstrap { n_bootstrap, seed }` is only available when starting from `AtmSampleSet`
- `None` suppresses uncertainty reporting

## IEXP and DEXP

Several different free energy methods utilize perturbation to estimate free energy differences, but when people say Free Energy Perturbation (FEP) they are typically referring to an exponential averaging scheme based on perturbing potential energies from one lambda state to another. This perturbation can take place in either the forward or the reverse direction. In accordance with Klimovich, Shirts, and Mobley's 2016 [Guidelines for the analysis of free energy calculations](https://pmc.ncbi.nlm.nih.gov/articles/PMC4420631/), we define `iexp` as insertion exponential averaging, typically proceeding in the direction of decreasing entropy, and `dexp` as deletion exponential averaging, typically proceeding in the direction of increasing entropy.

Types:

- `IexpEstimator`
- `IexpFit`
- `IexpOptions`

Input:

- slice of `UNkMatrix`

Behavior:

- computes pairwise adjacent exponential averaging along the lambda schedule
- `fit(...)` returns `IexpFit`
- `result()` materializes the full pairwise `DeltaFMatrix`
- `result_with_uncertainty()` includes delta-method uncertainty estimates from the exponential weights `exp(-W)` using an unbiased finite-sample variance
- preserves `lambda_labels()` from the input windows when available

CLI direction conventions:

- `iexp` reports the forward direction
- `dexp` reports the reverse direction

`DEXP` is not a separate estimator type in the library; it is a CLI convention over the same estimator results.

## NES

Types:

- `NesEstimator`
- `NesFit`
- `NesOptions`

Input:

- slice of `SwitchingTrajectory`

Behavior:

- expects all trajectories to share the same initial and final states
- applies the Jarzynski equality to the reduced switching work values
- `fit(...)` returns `NesFit`
- `result()` materializes a scalar `FreeEnergyEstimate`

Uncertainty behavior:

- analytic uncertainty is the default and is computed from the trajectory-level Jarzynski weights `exp(-W)` with an unbiased finite-sample variance estimate
- setting `NesOptions { n_bootstrap: N, .. }` with `N > 0` switches to bootstrap uncertainty across trajectories

Parsing / workflow expectations:

- the intended CLI input is one AMBER nonequilibrium switching trajectory per file
- each trajectory contributes one reduced work value
- when present, the retained `lambda_path` and `dvdl_path` support NES advisor diagnostics, but the estimator itself only needs the reduced work

## Parallel execution

Some estimators support a `parallel` option and use Rayon-backed parallel loops internally to speed up calculations.


## Common API shape

The estimator interfaces follow one shared pattern:

- `fit(...) -> Fit`
- `estimate(...) -> primary result`
- `Fit::result() -> primary result`

Estimators with an explicit uncertainty split also expose `estimate_with_uncertainty()` and `Fit::result_with_uncertainty()`.

ATM extends that shape with:

- `fit_leg(...)`
- `estimate_leg(...)`
- `estimate_binding(...)`
- `estimate_rbfe(...)`
- `estimate_abfe(...)`
