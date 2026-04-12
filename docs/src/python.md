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
- `python/examples/openmm_atm.py`: pure-OpenMM ATM-style toy example

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
python .\python\examples\openmm_atm.py
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
- `alchemrs.atm`
- `alchemrs.openmm`

Top-level estimator aliases are also available:

- `alchemrs.ATM`
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

For ATM-style workflows, the Python bindings expose a dedicated analysis model:

- `alchemrs.atm.AtmState`
- `alchemrs.atm.AtmSchedule`
- `alchemrs.atm.AtmSample`
- `alchemrs.atm.AtmSampleSet`
- `alchemrs.atm.ATM`
- `alchemrs.atm.schedule_from_arrays(...)`
- `alchemrs.atm.sample_set_from_arrays(...)`

This keeps the core design simple:

- equilibrium reduced potentials map to `UNkMatrix`
- nonequilibrium work values map to `SwitchingTrajectory`
- ATM schedule metadata plus per-sample `state_id`, `potE`, and `pertE` map to `AtmSampleSet`

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

## ATM-style OpenMM workflow

Use `python/examples/openmm_atm.py` when you want a concrete reference for:

1. defining an ATM leg schedule from arrays
2. sampling separate thermodynamic states in OpenMM
3. recording the per-sample quantities used by AToM-style analysis:
   - `state_id`
   - `potE`
   - `pertE`
4. converting those samples into `alchemrs.atm.AtmSampleSet`
5. running `alchemrs.ATM().estimate_leg(...)` and `estimate_rbfe(...)`

The example is intentionally analysis-focused. It does not parse a standard ATM
output format because OpenMM-based ATM workflows do not currently share one.
Instead, it documents the stable integration boundary for custom workflows:

- schedule metadata
- sampled state ids
- potential energies
- perturbation energies

That is the Python-level equivalent of the AToM analysis data that feeds its
UWHAM step.

## Scope and limitations

Current OpenMM support in this repository is example-driven, not parser-driven.

That means:

- no dedicated `src/parse/openmm.rs`
- no claim of drop-in support for arbitrary OpenMM output files
- examples focus on converting known in-memory quantities into `alchemrs`

This is intentional. OpenMM is a simulation engine, not a single canonical
free-energy file format, so the stable integration point is the validated
`alchemrs` data model rather than ad hoc file parsing.
