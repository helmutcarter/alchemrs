from __future__ import annotations

from typing import Sequence

import numpy as np

from . import StatePoint, SwitchingTrajectory, UNkMatrix


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
