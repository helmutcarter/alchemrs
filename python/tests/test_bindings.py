from __future__ import annotations

from pathlib import Path

import pytest

import alchemrs as ar


FIXTURES = Path(__file__).resolve().parents[2] / "fixtures" / "amber" / "acetamide_tiny"
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
