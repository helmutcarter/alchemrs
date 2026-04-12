from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

import numpy as np

from ._alchemrs import StatePoint, SwitchingTrajectory, UNkMatrix, atm

if TYPE_CHECKING:
    from . import AtmSampleSet, AtmSchedule


def kt_in_kcal_per_mol(temperature_k: float) -> float:
    return 0.00198720425864083 * temperature_k


def windows_from_u_kln(
    u_kln,
    lambdas: Sequence[float] | Sequence[Sequence[float]],
    temperature_k: float,
    *,
    time_ps=None,
    lambda_labels: Sequence[str] | None = None,
) -> list[UNkMatrix]:
    array = np.asarray(u_kln, dtype=float)
    if array.ndim != 3:
        raise ValueError("u_kln must have shape (sampled_state, evaluated_state, sample)")

    states = _state_grid(lambdas, temperature_k)
    n_sampled, n_evaluated, n_samples = array.shape
    if n_sampled != len(states) or n_evaluated != len(states):
        raise ValueError("u_kln state dimensions must match the lambda grid")

    if time_ps is None:
        time_values = np.arange(n_samples, dtype=float)
    else:
        time_values = np.asarray(time_ps, dtype=float)
        if time_values.ndim != 1:
            raise ValueError("time_ps must be one-dimensional")
        if len(time_values) != n_samples:
            raise ValueError("time_ps length must match the sample dimension of u_kln")

    windows: list[UNkMatrix] = []
    for sampled_idx, sampled_state in enumerate(states):
        windows.append(
            UNkMatrix(
                data=array[sampled_idx].T.copy(),
                time_ps=time_values.copy(),
                sampled_state=sampled_state,
                evaluated_states=states,
                lambda_labels=list(lambda_labels) if lambda_labels is not None else None,
            )
        )

    return windows


def switching_trajectories_from_work_values(
    reduced_work_values: Sequence[float],
    *,
    initial_lambdas: Sequence[float],
    final_lambdas: Sequence[float],
    temperature_k: float,
    lambda_paths: Sequence[Sequence[float]] | None = None,
    dvdl_paths: Sequence[Sequence[float]] | None = None,
    rms_dvdl_paths: Sequence[Sequence[float]] | None = None,
) -> list[SwitchingTrajectory]:
    n = len(reduced_work_values)
    lambda_paths = _optional_paths(lambda_paths, n, "lambda_paths")
    dvdl_paths = _optional_paths(dvdl_paths, n, "dvdl_paths")
    rms_dvdl_paths = _optional_paths(rms_dvdl_paths, n, "rms_dvdl_paths")

    initial_state = StatePoint(list(initial_lambdas), temperature_k)
    final_state = StatePoint(list(final_lambdas), temperature_k)

    trajectories: list[SwitchingTrajectory] = []
    for idx, reduced_work in enumerate(reduced_work_values):
        kwargs = {
            "initial_state": initial_state,
            "final_state": final_state,
            "reduced_work": float(reduced_work),
        }
        if lambda_paths is not None:
            kwargs["lambda_path"] = list(lambda_paths[idx])
        if dvdl_paths is not None:
            kwargs["dvdl_path"] = list(dvdl_paths[idx])
        if rms_dvdl_paths is not None:
            kwargs["rms_dvdl_path"] = list(rms_dvdl_paths[idx])
        trajectories.append(SwitchingTrajectory(**kwargs))

    return trajectories


def atm_schedule_from_lambdas(
    lambdas: Sequence[float],
    *,
    direction: str,
    temperature_k: float,
    alpha: Sequence[float] | float = 0.0,
    u0: Sequence[float] | float = 0.0,
    w0: Sequence[float] | float = 0.0,
    lambda3: Sequence[float] | float | None = None,
    u1: Sequence[float] | float | None = None,
) -> AtmSchedule:
    n_states = len(lambdas)
    if n_states == 0:
        raise ValueError("lambdas must contain at least one state")

    return atm.schedule_from_arrays(
        state_ids=list(range(n_states)),
        direction=direction,
        lambda1=[float(value) for value in lambdas],
        lambda2=[float(value) for value in lambdas],
        alpha=_broadcast_state_values(alpha, n_states, "alpha"),
        u0=_broadcast_state_values(u0, n_states, "u0"),
        w0=_broadcast_state_values(w0, n_states, "w0"),
        temperature_k=[float(temperature_k)] * n_states,
        lambda3=(
            _broadcast_state_values(lambda3, n_states, "lambda3")
            if lambda3 is not None
            else None
        ),
        u1=(
            _broadcast_state_values(u1, n_states, "u1")
            if u1 is not None
            else None
        ),
    )


def atm_sample_set_from_arrays(
    *,
    lambdas: Sequence[float],
    direction: str,
    temperature_k: float,
    state_ids: Sequence[int],
    potential_energies_kcal_per_mol,
    perturbation_energies_kcal_per_mol,
    alpha: Sequence[float] | float = 0.0,
    u0: Sequence[float] | float = 0.0,
    w0: Sequence[float] | float = 0.0,
    lambda3: Sequence[float] | float | None = None,
    u1: Sequence[float] | float | None = None,
) -> AtmSampleSet:
    schedule = atm_schedule_from_lambdas(
        lambdas,
        direction=direction,
        temperature_k=temperature_k,
        alpha=alpha,
        u0=u0,
        w0=w0,
        lambda3=lambda3,
        u1=u1,
    )
    return atm.sample_set_from_arrays(
        schedule,
        state_ids=list(state_ids),
        potential_energies_kcal_per_mol=potential_energies_kcal_per_mol,
        perturbation_energies_kcal_per_mol=perturbation_energies_kcal_per_mol,
    )


def _state_grid(
    lambdas: Sequence[float] | Sequence[Sequence[float]],
    temperature_k: float,
) -> list[StatePoint]:
    if not lambdas:
        raise ValueError("lambdas must contain at least one state")

    first = lambdas[0]
    if isinstance(first, (list, tuple, np.ndarray)):
        return [StatePoint(list(state), temperature_k) for state in lambdas]
    return [StatePoint([float(value)], temperature_k) for value in lambdas]


def _optional_paths(paths, expected_len: int, label: str):
    if paths is None:
        return None
    if len(paths) != expected_len:
        raise ValueError(f"{label} length must match reduced_work_values")
    return paths


def _broadcast_state_values(values, n_states: int, label: str) -> list[float]:
    if np.isscalar(values):
        return [float(values)] * n_states

    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{label} must be a scalar or one-dimensional sequence")
    if len(array) != n_states:
        raise ValueError(f"{label} length must match lambdas")
    return array.tolist()
