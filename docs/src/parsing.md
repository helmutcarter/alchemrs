# Parsing AMBER Outputs

The current parser implementation is AMBER-specific and lives in `alchemrs::parse::amber`.

## Available entry points

- `extract_dhdl(path, temperature_k)`
- `extract_u_nk(path, temperature_k)`
- `extract_u_nk_with_potential(path, temperature_k)`

## `extract_dhdl`

This parser extracts:

- `temp0`
- `clambda`
- `dt`
- `ntpr`
- the `begin time` coordinate
- `DV/DL` values from `NSTEP` blocks

The extracted gradients are converted to reduced units by multiplying by `beta = 1 / (k_B T)`.

Returned value:

- `DhdlSeries`

Common failure modes:

- missing metadata fields such as `temp0`, `clambda`, `dt`, or `ntpr`
- no `DV/DL` samples found
- temperature mismatch between the file and the requested temperature

## `extract_u_nk`

This parser reads AMBER MBAR output blocks and constructs a `UNkMatrix`.

Key behavior:

- the parser requires a valid MBAR lambda list in the file header
- `clambda` must be present in the evaluated-state grid
- the returned matrix keeps sample rows with positive infinity values
- rows that are entirely unusable can still cause parsing to fail if no finite MBAR samples remain

Returned value:

- `UNkMatrix`

## `extract_u_nk_with_potential`

This parser does everything `extract_u_nk` does, and additionally extracts `EPtot` from `NSTEP` blocks.

Returned value:

- `(UNkMatrix, Vec<f64>)`

The parser requires the number of `EPtot` samples to match the number of retained MBAR samples exactly.

This function exists because `EPtot` is often useful as a finite external observable for:

- auto-equilibration
- decorrelation

especially when the `u_nk` matrix contains positive infinity values that make `DE`-based preprocessing invalid.

## Temperature validation

The parser checks the temperature in the AMBER file against the temperature passed by the caller and fails if they differ by more than a small tolerance.

This is deliberate: parser output is temperature-dependent because reduced energies depend on `beta`.

## Current scope

The public parser surface in this repository is AMBER-only. The architecture leaves room for additional engines later, but the implemented and tested parser today is the AMBER parser.
