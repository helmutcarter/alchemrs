"""Pure-OpenMM ATM-style example for alchemrs.

This file avoids `openmmtools` so it can run with just:

    pip install openmm

It uses a one-particle harmonic toy system. The base potential is a harmonic
well centered at the origin, and the perturbation energy is defined as the
energy difference to a second displaced harmonic well. For alpha=0 and
lambda1=lambda2=lambda, the ATM bias reduces to `lambda * pertE`, so the
sampled Hamiltonian is a simple interpolation between the two wells.

The important part is the data boundary: for each sampled state we record
`state_id`, `potE`, and `pertE`, then hand those samples to `alchemrs.atm`.
That matches the analysis-side shape used by AToM without assuming any
particular engine output format.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
import random

import alchemrs as ar

from openmm import Context, CustomExternalForce, LangevinMiddleIntegrator, Platform, System, Vec3
from openmm import unit


K_B_KCAL_PER_MOL_K = 0.00198720425864083


@dataclass
class AtmConfig:
    temperature: unit.Quantity = 300 * unit.kelvin
    friction: unit.Quantity = 10.0 / unit.picoseconds
    timestep: unit.Quantity = 2.0 * unit.femtoseconds
    equilibration_steps: int = 200
    steps_per_sample: int = 10
    n_samples_per_state: int = 100
    spring_constant: unit.Quantity = (
        10.0 * unit.kilocalories_per_mole / unit.angstroms**2
    )
    mass: unit.Quantity = 39.9 * unit.amu
    x_start: unit.Quantity = 0.0 * unit.angstroms
    x_end: unit.Quantity = 5.0 * unit.angstroms
    lambdas: tuple[float, ...] = (0.0, 0.5, 1.0)
    platform_name: str | None = None


def make_atm_system(config: AtmConfig) -> System:
    system = System()
    system.addParticle(config.mass)

    force = CustomExternalForce(
        "0.5 * k * ((1 - lam) * x^2 + lam * (x - d)^2)"
    )
    force.addGlobalParameter(
        "k",
        config.spring_constant.value_in_unit(unit.kilocalories_per_mole / unit.nanometers**2),
    )
    force.addGlobalParameter("lam", 0.0)
    force.addGlobalParameter("d", config.x_end.value_in_unit(unit.nanometers))
    force.addParticle(0, [])
    system.addForce(force)

    return system


def make_context(config: AtmConfig) -> Context:
    system = make_atm_system(config)
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


def position_nm(context: Context) -> float:
    state = context.getState(getPositions=True)
    return state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)[0, 0]


def potential_energy_kcal_per_mol(context: Context) -> float:
    state = context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)


def perturbation_energy_kcal_per_mol(config: AtmConfig, x_nm: float, displaced_nm: float) -> float:
    k = config.spring_constant.value_in_unit(
        unit.kilocalories_per_mole / unit.nanometers**2
    )
    base = 0.5 * k * x_nm * x_nm
    displaced = 0.5 * k * (x_nm - displaced_nm) * (x_nm - displaced_nm)
    return displaced - base


def make_schedule(config: AtmConfig, direction: str):
    return ar.atm.schedule_from_arrays(
        state_ids=list(range(len(config.lambdas))),
        direction=direction,
        lambda1=list(config.lambdas),
        lambda2=list(config.lambdas),
        alpha=[0.0] * len(config.lambdas),
        u0=[0.0] * len(config.lambdas),
        w0=[0.0] * len(config.lambdas),
        temperature_k=[config.temperature.value_in_unit(unit.kelvin)] * len(config.lambdas),
    )


def collect_leg_samples(config: AtmConfig, *, direction: str, displacement: unit.Quantity):
    schedule = make_schedule(config, direction)
    displaced_nm = displacement.value_in_unit(unit.nanometers)
    state_ids: list[int] = []
    potential_energies: list[float] = []
    perturbation_energies: list[float] = []

    for state_id, lam in enumerate(config.lambdas):
        context = make_context(config)
        integrator = context.getIntegrator()
        context.setParameter("lam", float(lam))
        context.setParameter("d", displaced_nm)
        integrator.step(config.equilibration_steps)

        for _ in range(config.n_samples_per_state):
            integrator.step(config.steps_per_sample)
            x_nm = position_nm(context)
            state_ids.append(state_id)
            potential_energies.append(potential_energy_kcal_per_mol(context))
            perturbation_energies.append(
                perturbation_energy_kcal_per_mol(config, x_nm, displaced_nm)
            )

    return ar.atm.sample_set_from_arrays(
        schedule,
        state_ids=state_ids,
        potential_energies_kcal_per_mol=potential_energies,
        perturbation_energies_kcal_per_mol=perturbation_energies,
    )


def kt_in_kcal_per_mol(temperature_k: float) -> float:
    return ar.openmm.kt_in_kcal_per_mol(temperature_k)


def main() -> None:
    random.seed(0)

    config = AtmConfig()
    temperature_k = config.temperature.value_in_unit(unit.kelvin)
    conversion = kt_in_kcal_per_mol(temperature_k)

    leg1 = collect_leg_samples(
        config,
        direction="forward",
        displacement=5.0 * unit.angstroms,
    )
    leg2 = collect_leg_samples(
        config,
        direction="reverse",
        displacement=3.0 * unit.angstroms,
    )

    estimator = ar.ATM()
    leg1_result = estimator.estimate_leg(leg1)
    leg2_result = estimator.estimate_leg(leg2)
    binding = estimator.estimate_rbfe(leg1, leg2)

    print("OpenMM ATM-style analysis")
    print(f"  leg1 dG: {leg1_result.delta_f * conversion:.6f} kcal/mol")
    print(f"  leg2 dG: {leg2_result.delta_f * conversion:.6f} kcal/mol")
    print(f"  binding dG: {binding.delta_f * conversion:.6f} kcal/mol")
    print(f"  samples/state: {config.n_samples_per_state}")
    print(f"  states/leg: {len(config.lambdas)}")


if __name__ == "__main__":
    main()
