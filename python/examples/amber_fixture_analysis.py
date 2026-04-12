"""Process bundled AMBER example data through the Python bindings.

Run from the repository root with:

    $env:PYTHONPATH=".\python"
    python .\python\examples\amber_fixture_analysis.py

This example uses the `fixtures/amber/acetamide_tiny` dataset that is already
checked into the repository. It demonstrates:

- parsing AMBER `dH/dlambda` windows for TI
- parsing AMBER `u_nk` windows with `EPtot` for MBAR
- preprocessing through the Python bindings
- reporting `dG` values in `kcal/mol`
"""

from __future__ import annotations

from pathlib import Path

import alchemrs as ar


REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = REPO_ROOT / "fixtures" / "amber" / "acetamide_tiny"
TEMPERATURE_K = 300.0


def kt_in_kcal_per_mol(temperature_k: float) -> float:
    return ar.openmm.kt_in_kcal_per_mol(temperature_k)


def amber_windows() -> list[Path]:
    return sorted(
        FIXTURE_ROOT.glob("*/acetamide.prod.out"),
        key=lambda path: float(path.parent.name),
    )


def run_ti(paths: list[Path]) -> None:
    series = [ar.parse.extract_dhdl(path, TEMPERATURE_K) for path in paths]
    result = ar.TI().estimate(series)
    conversion = kt_in_kcal_per_mol(TEMPERATURE_K)

    print("TI")
    print(f"  windows: {len(series)}")
    print(f"  dG: {result.delta_f * conversion:.6f} kcal/mol")
    print(f"  sigma: {result.uncertainty * conversion:.6f} kcal/mol")
    print(f"  from lambda: {result.from_state.lambdas[0]:.2f}")
    print(f"  to lambda: {result.to_state.lambdas[0]:.2f}")


def run_mbar(paths: list[Path]) -> None:
    windows = []
    for path in paths:
        u_nk, epot = ar.parse.extract_u_nk_with_potential(path, TEMPERATURE_K)
        windows.append(ar.prep.decorrelate_u_nk_with_observable(u_nk, epot))

    fit = ar.MBAR().fit(windows)
    result = fit.result_with_uncertainty()
    conversion = kt_in_kcal_per_mol(TEMPERATURE_K)

    print("MBAR")
    print(f"  windows: {len(windows)}")
    print(
        "  dG: "
        f"{result.values[0, -1] * conversion:.6f} +/- "
        f"{result.uncertainties[0, -1] * conversion:.6f} kcal/mol"
    )
    print(f"  overlap scalar: {fit.overlap_scalar():.6f}")
    print(f"  from lambda: {result.states[0].lambdas[0]:.2f}")
    print(f"  to lambda: {result.states[-1].lambdas[0]:.2f}")


def main() -> None:
    paths = amber_windows()
    run_ti(paths)
    print()
    run_mbar(paths)


if __name__ == "__main__":
    main()
