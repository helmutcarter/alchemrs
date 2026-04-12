"""Pure-OpenMM equilibrium MBAR example for alchemrs.

This file avoids `openmmtools` so it can run with just:

    pip install openmm

It uses the same one-particle harmonic oscillator toy system as the NES
example, but samples equilibrium configurations at several lambda windows,
evaluates the reduced potential of every configuration at every lambda state,
and converts the resulting `u_kln[k][l][n]` tensor into the `UNkMatrix`
windows expected by `alchemrs`.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
import random

import alchemrs as ar
import numpy as np

from openmm import Context, CustomExternalForce, LangevinMiddleIntegrator, Platform, System, Vec3
from openmm import unit


K_B_KCAL_PER_MOL_K = 0.00198720425864083


@dataclass
class EquilibriumConfig:
    temperature: unit.Quantity = 300 * unit.kelvin
    friction: unit.Quantity = 10.0 / unit.picoseconds
    timestep: unit.Quantity = 2.0 * unit.femtoseconds
    equilibration_steps: int = 500
    steps_per_sample: int = 20
    n_samples_per_window: int = 200
    spring_constant: unit.Quantity = (
        10.0 * unit.kilocalories_per_mole / unit.angstroms**2
    )
    mass: unit.Quantity = 39.9 * unit.amu
    x_start: unit.Quantity = 0.0 * unit.angstroms
    x_end: unit.Quantity = 5.0 * unit.angstroms
    lambdas: tuple[float, ...] = (0.0, 0.5, 1.0)
    platform_name: str | None = None


def make_harmonic_system(config: EquilibriumConfig) -> System:
    system = System()
    system.addParticle(config.mass)

    force = CustomExternalForce("0.5 * k * (x - x0)^2")
    force.addGlobalParameter(
        "k", config.spring_constant.value_in_unit(unit.kilocalories_per_mole / unit.nanometers**2)
    )
    force.addGlobalParameter("x0", config.x_start.value_in_unit(unit.nanometers))
    force.addParticle(0, [])
    system.addForce(force)

    return system


def make_context(config: EquilibriumConfig) -> Context:
    system = make_harmonic_system(config)
    integrator = LangevinMiddleIntegrator(
        config.temperature, config.friction, config.timestep
    )
    platform = (
        Platform.getPlatformByName(config.platform_name)
        if config.platform_name is not None
        else Platform.getPlatformByName("Reference")
    )
    context = Context(system, integrator, platform)

    sigma = sqrt(
        (K_B_KCAL_PER_MOL_K * config.temperature.value_in_unit(unit.kelvin))
        / config.spring_constant.value_in_unit(unit.kilocalories_per_mole / unit.angstroms**2)
    )
    dx = random.gauss(0.0, sigma) * unit.angstroms
    context.setPositions([Vec3((config.x_start + dx).value_in_unit(unit.nanometers), 0.0, 0.0)])
    context.setVelocitiesToTemperature(config.temperature)
    return context


def lambda_to_x0_nm(config: EquilibriumConfig, lam: float) -> float:
    x_start_nm = config.x_start.value_in_unit(unit.nanometers)
    x_end_nm = config.x_end.value_in_unit(unit.nanometers)
    return x_start_nm + lam * (x_end_nm - x_start_nm)


def reduced_potential(context: Context, temperature_k: float) -> float:
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    return energy / (K_B_KCAL_PER_MOL_K * temperature_k)


def kt_in_kcal_per_mol(temperature_k: float) -> float:
    return K_B_KCAL_PER_MOL_K * temperature_k


def collect_u_kln(config: EquilibriumConfig) -> np.ndarray:
    temperature_k = config.temperature.value_in_unit(unit.kelvin)
    n_windows = len(config.lambdas)
    u_kln = np.zeros((n_windows, n_windows, config.n_samples_per_window), dtype=float)

    for sampled_idx, sampled_lambda in enumerate(config.lambdas):
        context = make_context(config)
        integrator = context.getIntegrator()
        context.setParameter("x0", lambda_to_x0_nm(config, sampled_lambda))
        integrator.step(config.equilibration_steps)

        for sample_idx in range(config.n_samples_per_window):
            integrator.step(config.steps_per_sample)
            for evaluated_idx, evaluated_lambda in enumerate(config.lambdas):
                context.setParameter("x0", lambda_to_x0_nm(config, evaluated_lambda))
                u_kln[sampled_idx, evaluated_idx, sample_idx] = reduced_potential(
                    context, temperature_k
                )
            context.setParameter("x0", lambda_to_x0_nm(config, sampled_lambda))

    return u_kln


def windows_from_u_kln(
    u_kln: np.ndarray,
    lambdas: tuple[float, ...],
    temperature_k: float,
) -> list[ar.UNkMatrix]:
    if u_kln.ndim != 3:
        raise ValueError("u_kln must have shape (sampled_state, evaluated_state, sample)")
    if u_kln.shape[0] != len(lambdas) or u_kln.shape[1] != len(lambdas):
        raise ValueError("u_kln state dimensions must match the lambda grid")

    evaluated_states = [
        ar.StatePoint([lam], temperature_k)
        for lam in lambdas
    ]

    windows: list[ar.UNkMatrix] = []
    for sampled_idx, sampled_lambda in enumerate(lambdas):
        data = u_kln[sampled_idx].T.copy()
        time_ps = np.arange(data.shape[0], dtype=float)
        windows.append(
            ar.UNkMatrix(
                data=data,
                time_ps=time_ps,
                sampled_state=ar.StatePoint([sampled_lambda], temperature_k),
                evaluated_states=evaluated_states,
            )
        )

    return windows


def main() -> None:
    random.seed(0)
    np.random.seed(0)

    config = EquilibriumConfig()
    temperature_k = config.temperature.value_in_unit(unit.kelvin)
    conversion = kt_in_kcal_per_mol(temperature_k)

    u_kln = collect_u_kln(config)
    windows = windows_from_u_kln(u_kln, config.lambdas, temperature_k)

    fit = ar.MBAR().fit(windows)
    result = fit.result_with_uncertainty()

    print("OpenMM u_kln -> alchemrs MBAR")
    print(
        "  dG: "
        f"{result.values[0, -1] * conversion:.6f} +/- "
        f"{result.uncertainties[0, -1] * conversion:.6f} kcal/mol"
    )
    print(f"  overlap scalar: {fit.overlap_scalar():.6f}")
    print(f"  windows: {len(windows)}")
    print(f"  samples per window: {config.n_samples_per_window}")


if __name__ == "__main__":
    main()
