# Python and OpenMM

The repository now includes Python bindings and pure-OpenMM example workflows in
addition to the Rust CLI and library APIs.

This page focuses on local development usage from the repository checkout. A
more polished installation workflow can be added later when the Python package
build/release flow is finalized.

## Current layout

Relevant paths:

- `python/alchemrs`: Python package wrapper around the native extension
- `python/tests`: Python-side test coverage
- `python/examples/amber_fixture_analysis.py`: bundled AMBER fixture analysis through the Python API
- `python/examples/openmm_u_kln_mbar.py`: pure-OpenMM equilibrium MBAR toy example
- `python/examples/openmm_nes.py`: pure-OpenMM nonequilibrium switching toy example

## Local usage

The local package is not installed into your Python environment automatically.
For now, run the examples and tests from the repository root with `python/` on
`PYTHONPATH`.

PowerShell:

```powershell
$env:PYTHONPATH=".\python"
python .\python\examples\amber_fixture_analysis.py
python .\python\examples\openmm_u_kln_mbar.py
python .\python\examples\openmm_nes.py
python -m pytest .\python\tests -q
```

## Bundled AMBER fixture workflow

Use `python/examples/amber_fixture_analysis.py` for the simplest end-to-end
Python example against real repository data. It:

- parses the bundled `fixtures/amber/acetamide_tiny` dataset
- runs TI directly from `dH/dlambda`
- runs MBAR from `u_nk` plus `EPtot` decorrelation
- prints `dG` values in `kcal/mol`

## Python API shape

The bindings mirror the Rust crate structure at a high level:

- `alchemrs.parse`
- `alchemrs.prep`
- `alchemrs.estimators`
- `alchemrs.analysis`
- `alchemrs.openmm`

Top-level estimator aliases are also available:

- `alchemrs.MBAR`
- `alchemrs.BAR`
- `alchemrs.TI`
- `alchemrs.IEXP`
- `alchemrs.NES`

## OpenMM support model

`alchemrs` does not currently ship a dedicated OpenMM parser. Instead, the
Python side exposes conversion helpers and runnable examples that turn OpenMM
quantities into the validated `alchemrs` data model.

Current reusable helper module:

- `alchemrs.openmm.kt_in_kcal_per_mol(temperature_k)`
- `alchemrs.openmm.windows_from_u_kln(...)`
- `alchemrs.openmm.switching_trajectories_from_work_values(...)`

This keeps the core design simple:

- equilibrium reduced potentials map to `UNkMatrix`
- nonequilibrium work values map to `SwitchingTrajectory`

Once those objects exist, the standard `prep`, `estimators`, and `analysis`
APIs apply unchanged.

## Equilibrium OpenMM workflow

Use `python/examples/openmm_u_kln_mbar.py` when you want a concrete reference
for:

1. sampling equilibrium configurations at multiple lambda windows
2. evaluating all samples on all lambda states to form `u_kln`
3. converting `u_kln` into `UNkMatrix` windows
4. running `alchemrs.MBAR()`

The example prints `dG` in `kcal/mol` and reports the MBAR overlap scalar.

## Nonequilibrium OpenMM workflow

Use `python/examples/openmm_nes.py` when you want a concrete reference for:

1. running forward-only switching trajectories
2. running separate forward and reverse switching campaigns
3. converting switching work into `SwitchingTrajectory`
4. running `alchemrs.NES()` and `alchemrs.analysis.nes_convergence(...)`

The current example uses a toy harmonic oscillator so it can run with plain
`openmm` from `pip install openmm`.

## Scope and limitations

Current OpenMM support in this repository is example-driven, not parser-driven.

That means:

- no dedicated `src/parse/openmm.rs`
- no claim of drop-in support for arbitrary OpenMM output files
- examples focus on converting known in-memory quantities into `alchemrs`

This is intentional. OpenMM is a simulation engine, not a single canonical
free-energy file format, so the stable integration point is the validated
`alchemrs` data model rather than ad hoc file parsing.
