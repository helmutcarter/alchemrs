# Data Model

The `alchemrs::data` module defines the typed structures that the rest of the library operates on.

## `StatePoint`

`StatePoint` stores:

- one or more lambda values
- temperature in Kelvin

Validation rules:

- temperature must be finite and positive
- every lambda value must be finite

The type itself can store multiple lambda dimensions.

Current estimator support is split:

- `UNkMatrix`-based library estimators support multidimensional lambda states when windows share the same evaluated-state grid.
- `DhdlSeries`-based TI remains one-dimensional.
- The CLI supports multidimensional `u_nk` estimator inputs and renders multidimensional endpoints as lambda vectors.

## `DhdlSeries`

`DhdlSeries` represents a time-indexed `dH/dlambda` series for one state.

It stores:

- the associated `StatePoint`
- `time_ps`
- scalar values

Validation rules:

- `time_ps.len() == values.len()`
- all times must be finite
- all values must be finite
- time must be monotonic nondecreasing

This is the primary input type for TI.

## `SwitchingTrajectory`

`SwitchingTrajectory` represents one nonequilibrium switching trajectory.

It stores:

- the initial `StatePoint`
- the final `StatePoint`
- one reduced work value for the full trajectory
- an optional `lambda_path`
- an optional `dvdl_path` in reduced units

Validation rules:

- the initial and final states must share the same temperature
- `lambda_path.len() == dvdl_path.len()`
- all stored values must be finite

This is the primary input type for NES. The reduced work is sufficient for the Jarzynski estimator itself, while `lambda_path` and `dvdl_path` support the NES advisor's profile and curvature diagnostics.

## `UNkMatrix`

`UNkMatrix` stores reduced energies evaluated across states.

Conceptually:

- rows are samples
- columns are evaluated states
- `sampled_state` identifies the state the trajectory was sampled from
- `evaluated_states` identifies the column states
- `lambda_labels`, when present, name the lambda-vector components in order

Validation rules:

- `data.len() == n_samples * n_states`
- `evaluated_states.len() == n_states`
- `time_ps.len() == n_samples`
- time must be monotonic nondecreasing
- `u_nk` values may be finite or positive infinity, but not `NaN` or negative infinity

That last point matters for preprocessing:

- some `u_nk`-derived observables fail on positive infinity
- observable-based preprocessing using an external potential-energy series such as `EPtot` can still work if the external observable is finite

## `AtmState` and `AtmSchedule`

ATM workflows use dedicated schedule types in addition to generic `StatePoint`.

`AtmState` stores:

- `state_id`
- leg direction (`Forward` or `Reverse`)
- standard-softplus parameters `lambda1`, `lambda2`, `alpha`, `u0`, `w0`
- optional multi-softplus parameters `lambda3`, `u1`
- temperature in Kelvin

Validation rules:

- all numeric parameters must be finite
- temperature must be finite and positive
- `lambda3` and `u1` must either both be present or both be absent

`AtmSchedule` stores one ATM leg as an ordered set of `AtmState` values.

Validation rules:

- at least one state
- unique `state_id`
- one leg direction per schedule
- one temperature per schedule
- either all states use standard softplus or all use multi-softplus

Each ATM state can also be converted into a `StatePoint`. The resulting state
vector is ordered as `lambda1`, `lambda2`, optional `lambda3`, `alpha`, `u0`,
optional `u1`, `w0`, and those names are retained as `lambda_labels` when the
ATM data is converted into matrix results.

## `AtmSample` and `AtmSampleSet`

`AtmSample` represents one ATM observation.

It stores:

- `state_id`
- `potential_energy_kcal_per_mol`
- `perturbation_energy_kcal_per_mol`

Validation rules:

- both stored energies must be finite

`AtmSampleSet` combines:

- one `AtmSchedule`
- one or more `AtmSample` values assigned to that schedule

Validation rules:

- at least one sample
- every sample `state_id` must exist in the schedule

This is the most convenient library and Python entry point for ATM leg analysis,
because it can be converted into a pooled `AtmLogQMatrix` automatically.

## `AtmLogQMatrix`

`AtmLogQMatrix` stores pooled log unnormalized densities for ATM/UWHAM-style
analysis.

It stores:

- `n_observations`
- `n_states`
- dense `log_q` data in observation-major order
- the ordered state list
- `sampled_counts` per state
- `lambda_labels`, when present

Validation rules:

- `log_q.len() == n_observations * n_states`
- `states.len() == n_states`
- `sampled_counts.len() == n_states`
- `sum(sampled_counts) == n_observations`
- `log_q` values may be finite or negative infinity, but not `NaN` or positive infinity

## `FreeEnergyEstimate`

Represents a scalar free energy result with:

- `delta_f`
- optional uncertainty
- `from_state`
- `to_state`

This is used by TI, NES, and ATM leg results.

## `DeltaFMatrix`

Represents pairwise free energy differences between all states.

It stores:

- a dense `n_states x n_states` matrix of values
- optional uncertainties
- the ordered state list
- `lambda_labels`, when present, naming the lambda-vector components in state order

This is used by BAR, MBAR, IEXP, DEXP, UWHAM, and ATM matrix results.

## `OverlapMatrix`

Represents the overlap matrix between states and is used by the diagnostics layer.

It stores:

- a dense `n_states x n_states` matrix of overlap values
- the ordered state list
- `lambda_labels`, when present, naming the lambda-vector components in state order

## Error model

The crate uses the shared `CoreError` type from `alchemrs::error` for:

- shape mismatches
- invalid state metadata
- invalid time ordering
- non-finite numerical values
- parse failures converted into the shared error model
- convergence failures
- unsupported inputs

This keeps error handling consistent across parsing, preprocessing, estimation, and analysis.
