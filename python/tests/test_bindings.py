from __future__ import annotations

import importlib.util
from pathlib import Path
import random
import sys

import numpy as np
import pytest

import alchemrs as ar


FIXTURES = Path(__file__).resolve().parents[2] / "fixtures" / "amber" / "acetamide_tiny"
EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
TEMPERATURE_K = 300.0


def _amber_windows() -> list[Path]:
    return sorted(FIXTURES.glob("*/acetamide.prod.out"), key=lambda path: float(path.parent.name))


def test_parse_and_detect_equilibration_against_fixture_reference() -> None:
    series = ar.parse.extract_dhdl(FIXTURES / "0.1" / "acetamide.prod.out", TEMPERATURE_K)

    assert series.state.lambdas == pytest.approx([0.1])
    assert series.state.temperature_k == pytest.approx(TEMPERATURE_K)
    assert series.time_ps.ndim == 1
    assert series.values.ndim == 1
    assert len(series.time_ps) == len(series.values)

    eq = ar.prep.detect_equilibration_dhdl(series)
    assert eq.t0 == 0
    assert eq.g == pytest.approx(1.0)
    assert eq.neff_max == pytest.approx(20.0)


def test_ti_matches_fixture_reference() -> None:
    series = [ar.parse.extract_dhdl(path, TEMPERATURE_K) for path in _amber_windows()]

    result = ar.TI().estimate(series)

    assert result.delta_f == pytest.approx(13.50045227, abs=1e-6)
    assert result.uncertainty == pytest.approx(0.81741622, abs=1e-6)
    assert result.from_state.lambdas == pytest.approx([0.0])
    assert result.to_state.lambdas == pytest.approx([1.0])


def test_mbar_and_overlap_match_fixture_references() -> None:
    windows = [ar.parse.extract_u_nk(path, TEMPERATURE_K) for path in _amber_windows()]

    fit = ar.MBAR().fit(windows)
    result = fit.result_with_uncertainty()

    assert result.values[0, -1] == pytest.approx(-113.58992316, abs=1e-5)
    assert result.uncertainties[0, -1] == pytest.approx(1.16998889, abs=1e-6)
    assert fit.overlap_scalar() == pytest.approx(0.01834742779436016, abs=1e-7)


def test_analysis_convergence_and_overlap_module_exports() -> None:
    windows = [ar.parse.extract_u_nk(path, TEMPERATURE_K) for path in _amber_windows()]

    overlap = ar.analysis.overlap_matrix(windows)
    eigenvalues = ar.analysis.overlap_eigenvalues(overlap)
    points = ar.analysis.mbar_convergence(windows)

    assert overlap.values.shape[0] == overlap.values.shape[1] == len(windows)
    assert len(eigenvalues) == len(windows)
    assert len(points) == len(windows)
    assert points[-1].delta_f == pytest.approx(-113.58992316, abs=1e-5)


def test_openmm_helper_windows_from_u_kln_and_switching_trajectories() -> None:
    u_kln = np.array(
        [
            [[0.0, 0.1], [0.5, 0.6]],
            [[0.4, 0.3], [0.0, 0.1]],
        ],
        dtype=float,
    )
    windows = ar.openmm.windows_from_u_kln(u_kln, [0.0, 1.0], 300.0)

    assert len(windows) == 2
    assert windows[0].n_samples == 2
    assert windows[0].n_states == 2
    assert np.allclose(windows[0].data, np.array([[0.0, 0.5], [0.1, 0.6]]))

    trajectories = ar.openmm.switching_trajectories_from_work_values(
        [1.0, 2.0],
        initial_lambdas=[0.0],
        final_lambdas=[1.0],
        temperature_k=300.0,
        lambda_paths=[[0.0, 0.5], [0.0, 0.5]],
        dvdl_paths=[[1.0, 1.5], [2.0, 2.5]],
    )
    assert len(trajectories) == 2
    assert trajectories[0].reduced_work == pytest.approx(1.0)
    assert len(trajectories[0].lambda_path) == 2
    assert len(trajectories[0].dvdl_path) == 2


def test_amber_fixture_python_example_matches_expected_ranges() -> None:
    module = load_example_module("amber_fixture_analysis")
    paths = module.amber_windows()

    assert len(paths) == 15

    series = [ar.parse.extract_dhdl(path, module.TEMPERATURE_K) for path in paths]
    ti_result = ar.TI().estimate(series)
    ti_kcal = ti_result.delta_f * module.kt_in_kcal_per_mol(module.TEMPERATURE_K)

    windows = []
    for path in paths:
        u_nk, epot = ar.parse.extract_u_nk_with_potential(path, module.TEMPERATURE_K)
        windows.append(ar.prep.decorrelate_u_nk_with_observable(u_nk, epot))

    mbar_fit = ar.MBAR().fit(windows)
    mbar_result = mbar_fit.result_with_uncertainty()
    mbar_kcal = mbar_result.values[0, -1] * module.kt_in_kcal_per_mol(module.TEMPERATURE_K)

    assert ti_kcal == pytest.approx(8.048447, abs=1e-6)
    assert mbar_kcal == pytest.approx(-75.004329, abs=1e-6)
    assert mbar_fit.overlap_scalar() == pytest.approx(0.005514, abs=1e-6)


@pytest.mark.skipif(
    importlib.util.find_spec("openmm") is None,
    reason="openmm is required for the runnable OpenMM examples",
)
def test_openmm_equilibrium_example_conversion_and_mbar() -> None:
    module = load_example_module("openmm_u_kln_mbar")
    random.seed(0)
    np.random.seed(0)

    config = module.EquilibriumConfig(
        equilibration_steps=50,
        steps_per_sample=5,
        n_samples_per_window=20,
        lambdas=(0.0, 0.5, 1.0),
    )
    temperature_k = config.temperature.value_in_unit(module.unit.kelvin)

    u_kln = module.collect_u_kln(config)
    windows = ar.openmm.windows_from_u_kln(u_kln, config.lambdas, temperature_k)
    fit = ar.MBAR().fit(windows)
    result = fit.result_with_uncertainty()

    assert u_kln.shape == (3, 3, 20)
    assert len(windows) == 3
    assert all(window.n_samples == 20 for window in windows)
    assert np.isfinite(result.values[0, -1])
    assert np.isfinite(result.uncertainties[0, -1])
    assert np.isfinite(fit.overlap_scalar())


@pytest.mark.skipif(
    importlib.util.find_spec("openmm") is None,
    reason="openmm is required for the runnable OpenMM examples",
)
def test_openmm_nes_example_forward_and_reverse_paths() -> None:
    module = load_example_module("openmm_nes")
    random.seed(0)

    config = module.SwitchingConfig(
        equilibration_steps=50,
        n_steps_per_switch=20,
        n_switches=8,
    )

    forward = module.run_switches(config, lambda_start=0.0, lambda_end=1.0)
    reverse = module.run_switches(config, lambda_start=1.0, lambda_end=0.0)

    forward_result = ar.NES().estimate(forward)
    reverse_result = ar.NES().estimate(reverse)

    assert len(forward) == 8
    assert len(reverse) == 8
    assert all(len(traj.lambda_path) == 20 for traj in forward)
    assert all(len(traj.dvdl_path) == 20 for traj in forward)
    assert np.isfinite(forward_result.delta_f)
    assert np.isfinite(forward_result.uncertainty)
    assert np.isfinite(reverse_result.delta_f)
    assert np.isfinite(reverse_result.uncertainty)


def load_example_module(name: str):
    path = EXAMPLES / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"python_examples_{name}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module
