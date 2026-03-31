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

## `UNkMatrix`

`UNkMatrix` stores reduced energies evaluated across states.

Conceptually:

- rows are samples
- columns are evaluated states
- `sampled_state` identifies the state the trajectory was sampled from
- `evaluated_states` identifies the column states

Validation rules:

- `data.len() == n_samples * n_states`
- `evaluated_states.len() == n_states`
- `time_ps.len() == n_samples`
- time must be monotonic nondecreasing
- `u_nk` values may be finite or positive infinity, but not `NaN` or negative infinity

That last point matters for preprocessing:

- some `u_nk`-derived observables fail on positive infinity
- observable-based preprocessing using `EPtot` can still work if the external observable is finite

## `FreeEnergyEstimate`

Represents a scalar free energy result with:

- `delta_f`
- optional uncertainty
- `from_state`
- `to_state`

This is used by TI.

## `DeltaFMatrix`

Represents pairwise free energy differences between all states.

It stores:

- a dense `n_states x n_states` matrix of values
- optional uncertainties
- the ordered state list

This is used by BAR, MBAR, EXP, and DEXP.

## `OverlapMatrix`

Represents the overlap matrix between states and is used by the diagnostics layer.

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
