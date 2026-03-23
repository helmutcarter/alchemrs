# Preprocessing

Preprocessing lives in `alchemrs-prep` and is where most of the scientifically important workflow details live.

## Core concepts

The prep crate supports:

- duplicate-time cleanup
- sorting by time
- optional time slicing
- equilibration detection
- decorrelation / subsampling

The main option struct is `DecorrelationOptions`.

Default values are:

- `drop_duplicates = true`
- `sort = true`
- `conservative = true`
- `remove_burnin = false`
- `fast = false`
- `nskip = 1`
- `lower = None`
- `upper = None`
- `step = None`

## Time cleanup

Before decorrelation or equilibration detection, the prep crate can:

- drop duplicate time values
- sort by time

This matters because the timeseries logic assumes one contiguous series ordered in time.

For `UNkMatrix`, cleanup preserves row/state alignment by selecting or reordering entire rows, not individual values.

## Time slicing

`lower`, `upper`, and `step` allow simple time-domain slicing before decorrelation.

These options are applied to:

- the series itself for `DhdlSeries`
- both the selected scalar observable and the aligned `u_nk` rows for `UNkMatrix`

## Statistical inefficiency and `g`

Decorrelation is based on an estimate of statistical inefficiency `g`.

The implementation follows the same basic `pymbar.timeseries` logic:

- center the series
- estimate autocorrelation contributions at increasing lag
- stop once the correlation function becomes non-positive after a minimum lag
- clamp `g` to at least `1`

## `fast`

`fast` changes how `g` is estimated.

- `fast = false` evaluates lag contributions one lag at a time
- `fast = true` increases the lag increment as it goes, which is cheaper but less accurate

This affects:

- plain decorrelation
- equilibration detection

because equilibration detection repeatedly estimates `g` on suffixes of the series.

## `conservative`

`conservative` changes how the code turns `g` into retained sample indices.

- `conservative = true` uses a uniform stride of `ceil(g)`
- `conservative = false` uses rounded multiples of fractional `g`

The conservative mode keeps fewer samples and is intentionally cautious.

## Equilibration detection

Equilibration detection returns:

- `t0`: the chosen start index for equilibrated data
- `g`: statistical inefficiency for the retained suffix
- `neff_max`: the estimated number of effectively uncorrelated samples

The heuristic scans possible start positions and chooses the one that maximizes:

```text
Neff = (N - t0 + 1) / g
```

This is the same style of automated equilibration detection used by `pymbar`. For further reading on the topic, see [Chodera 2016](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00784).

## `remove_burnin`

In the prep crate, `remove_burnin = true` means:

- run equilibration detection
- discard all data before `t0`
- then decorrelate the remaining suffix using the `g` estimated for that suffix

This is different from the CLI flag `--remove-burnin <N>`, which is a fixed-count trim done before any automated detection.

## `u_nk` observables

`u_nk` decorrelation needs a scalar time series. The prep crate supports two native derived observables and one external-observable path.

### `UNkSeriesMethod::DE`

For each sample row:

1. identify the sampled-state column
2. use the next evaluated-state column if it exists
3. otherwise use the previous column for the last state
4. compute:

```text
DE_t = u_nk[t, other] - u_nk[t, sampled]
```

This matches the `alchemlyb` `dE` convention.

Important detail:

- “adjacent” means adjacent in evaluated-state order, not nearest by numeric lambda distance

### `UNkSeriesMethod::All`

For each sample row, sum all evaluated-state reduced energies in that row.

This is a generic matrix-derived scalar, but it is usually less targeted than `DE`.

### `decorrelate_u_nk_with_observable`

This path accepts an external scalar observable, most commonly `EPtot`.

The observable determines which sample indices are retained, and those retained indices are then applied back to the full `u_nk` rows.

This is useful when:

- the matrix contains infinite energy values that make `DE` invalid

## Non-finite `u_nk` behavior

`decorrelate_u_nk` rejects non-finite derived scalar series.

In practice:

- `DE` or `All` can fail if the derived scalar is not finite
- `decorrelate_u_nk_with_observable` can still succeed if the supplied observable itself is finite

That is why the CLI exposes `epot` as an observable choice instead of only supporting `de`.

## CLI preprocessing order

For CLI commands, preprocessing order is:

1. fixed-count `--remove-burnin <N>`
2. `--auto-equilibrate`
3. `--decorrelate`

TI uses `dH/dlambda` as the scalar series.

`BAR`, `MBAR`, `EXP`, and `DEXP` use the observable chosen by `--u-nk-observable <de|all|epot>`.

## CLI auto-equilibrate overrides

When `--auto-equilibrate` is enabled in the CLI, the effective preprocessing flags become:

- `fast = true`
- `conservative = false`

even if the user supplied different `--fast` or `--conservative` values.

This is intentional and mirrors the `alchemlyb` / `pymbar` automated-equilibration workflow:

- `detect_equilibration(..., fast=True)`
- then `subsample_correlated_data(..., conservative=False)`

Those effective values are still recorded in CLI provenance so the output remains auditable.
