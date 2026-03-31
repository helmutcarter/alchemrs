# Parsing Outputs

The parser implementation supports AMBER outputs in `alchemrs::parse::amber` and GROMACS `dhdl.xvg` outputs in `alchemrs::parse::gromacs`.

The top-level parser entry points `alchemrs::extract_dhdl`, `alchemrs::extract_u_nk`, and `alchemrs::extract_u_nk_with_potential` auto-detect between those supported formats.

## Available entry points

- `extract_dhdl(path, temperature_k)`
- `extract_u_nk(path, temperature_k)`
- `extract_u_nk_with_potential(path, temperature_k)`

## `extract_dhdl`

Returned value:

- `DhdlSeries`

AMBER parser behavior:

- extracts `temp0`, `clambda`, `dt`, `ntpr`, the `begin time` coordinate, and `DV/DL` values from `NSTEP` blocks
- converts gradients to reduced units by multiplying by `beta = 1 / (k_B T)`

GROMACS parser behavior:

- reads `dhdl.xvg` headers to identify the simulation temperature and lambda
- extracts the `dH/dlambda` series from the legend tagged with `dH` when the file contains exactly one derivative component
- keeps the x-axis values as the returned `time_ps` values
- rejects multidimensional schedules with multiple `dH/dlambda` components because `DhdlSeries` is still scalar

Common failure modes:

- missing temperature or lambda metadata
- missing `dH/dlambda` samples
- temperature mismatch between the file and the requested temperature

## `extract_u_nk`

Returned value:

- `UNkMatrix`

AMBER parser behavior:

- reads MBAR output blocks and constructs a `UNkMatrix`
- requires a valid MBAR lambda list in the file header
- requires `clambda` to be present in the evaluated-state grid
- keeps sample rows with positive infinity values
- fails if no finite MBAR samples remain after filtering unusable rows

GROMACS parser behavior:

- reads `dhdl.xvg` Delta H columns and maps them into reduced-energy rows, including multidimensional lambda tuples from legends like `to (coul, vdw)`
- inserts the sampled lambda state explicitly with zero reduced energy
- sorts the evaluated-state grid by lambda before building the matrix
- requires at least one foreign-lambda Delta H column

## `extract_u_nk_with_potential`

Returned value:

- `(UNkMatrix, Vec<f64>)`

AMBER parser behavior:

- does everything `extract_u_nk` does and additionally extracts `EPtot` from `NSTEP` blocks
- requires the number of `EPtot` samples to match the retained MBAR sample count exactly

GROMACS parser behavior:

- does everything `extract_u_nk` does and additionally extracts the first energy-like legend labeled `Potential Energy` or `Total Energy` from `dhdl.xvg`
- requires that energy observable to be present when this entry point is used

This function exists because an external energy observable is often useful for:

- auto-equilibration
- decorrelation

especially when the `u_nk` matrix contains positive infinity values that make `DE`-based preprocessing invalid.

## Temperature validation

The parser checks the temperature in the AMBER file against the temperature passed by the caller and fails if they differ by more than a small tolerance.

This is deliberate: parser output is temperature-dependent because reduced energies depend on `beta`.

## Current scope

The public parser surface in this repository supports AMBER outputs and GROMACS `dhdl.xvg` outputs. The architecture still leaves room for additional engines later.
