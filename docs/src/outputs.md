# Outputs and Provenance

The CLI supports three output formats:

- text
- JSON
- CSV

## Units

Output units can be selected with:

```bash
--output-units <kt|kcal|kj>
```

The CLI converts the reported free energy and uncertainty from reduced units into the selected output unit using the requested temperature.

## Text output

Text output is intended for direct reading in the terminal.

It includes:

- `delta_f`
- uncertainty
- lambda endpoints, which are scalars for one-dimensional schedules and vectors for multidimensional schedules
- sample counts
- `u_nk_observable` when relevant
- overlap summary when requested

## JSON output

JSON output is suitable for shell pipelines and downstream tooling.

The payload contains:

- the scalar result fields
- optional overlap summary
- a provenance object

`from_lambda` and `to_lambda` are encoded as:

- numbers for one-dimensional schedules
- arrays for multidimensional schedules

## CSV output

CSV output appends provenance fields after the result columns so the row remains self-describing.

For multidimensional schedules, the endpoint columns are emitted as quoted bracketed vectors such as `"[0;0;0.8]"`.

## Provenance

In this project, “provenance” means the analysis settings and sample-count metadata that explain how a reported result was produced.

Current provenance fields include:

- estimator name
- temperature
- whether decorrelation ran
- fixed-count burn-in trimming
- whether auto-equilibration ran
- the effective `fast` value
- the effective `conservative` value
- `nskip`
- `u_nk_observable` when applicable
- number of windows
- number of samples before preprocessing
- number of samples after burn-in trimming / auto-equilibration
- number of retained samples after decorrelation

This matters because preprocessing choices directly change the data used by the estimators.

## Example JSON payload

```json
{
  "delta_f": -113.58,
  "uncertainty": 1.17,
  "from_lambda": [0.0, 0.0, 0.8],
  "to_lambda": [0.0, 0.0, 0.9],
  "units": "kT",
  "overlap": {
    "scalar": 0.0183,
    "eigenvalues": [1.0, 0.9817]
  },
  "provenance": {
    "estimator": "mbar",
    "temperature_k": 300.0,
    "decorrelate": true,
    "remove_burnin": 0,
    "auto_equilibrate": false,
    "fast": false,
    "conservative": true,
    "nskip": 1,
    "u_nk_observable": "de",
    "windows": 15,
    "samples_in": 300,
    "samples_after_burnin": 300,
    "samples_kept": 126
  }
}
```
