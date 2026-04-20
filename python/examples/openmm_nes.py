"""Pure-OpenMM NES examples for alchemrs.

This file avoids `openmmtools` so it can run with just:

    pip install openmm

It uses a minimal one-particle harmonic oscillator where the center of the
potential is moved along a lambda protocol through a global context parameter.
For each switching run, the example accumulates the nonequilibrium protocol
work explicitly from the instantaneous reduced-energy change caused by each
lambda update, then converts that into `alchemrs.SwitchingTrajectory`.

This is a toy model intended to demonstrate the data conversion into alchemrs,
not a production alchemical workflow.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import sqrt
import random

import alchemrs as ar

from openmm import Context, CustomExternalForce, LangevinMiddleIntegrator, Platform, System, Vec3
from openmm import unit


K_B_KCAL_PER_MOL_K = 0.00198720425864083


@dataclass
class SwitchingConfig:
    temperature: unit.Quantity = field(default_factory=lambda: 300 * unit.kelvin)
    friction: unit.Quantity = field(default_factory=lambda: 10.0 / unit.picoseconds)
    timestep: unit.Quantity = field(default_factory=lambda: 2.0 * unit.femtoseconds)
    equilibration_steps: int = 200
    n_steps_per_switch: int = 100
    n_switches: int = 50
    spring_constant: unit.Quantity = field(
        default_factory=lambda: 10.0 * unit.kilocalories_per_mole / unit.angstroms**2
    )
    mass: unit.Quantity = field(default_factory=lambda: 39.9 * unit.amu)
    x_start: unit.Quantity = field(default_factory=lambda: 0.0 * unit.angstroms)
    x_end: unit.Quantity = field(default_factory=lambda: 5.0 * unit.angstroms)
    platform_name: str | None = None


def make_harmonic_system(config: SwitchingConfig) -> System:
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


def make_context(config: SwitchingConfig) -> Context:
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
    x0 = config.x_start.value_in_unit(unit.nanometers)
    dx = random.gauss(0.0, sigma) * unit.angstroms
    context.setPositions([Vec3((config.x_start + dx).value_in_unit(unit.nanometers), 0.0, 0.0)])
    context.setVelocitiesToTemperature(config.temperature)
    return context


def reduced_potential(context: Context, temperature_k: float) -> float:
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    return energy / (K_B_KCAL_PER_MOL_K * temperature_k)


def kt_in_kcal_per_mol(temperature_k: float) -> float:
    return K_B_KCAL_PER_MOL_K * temperature_k


def reduced_dvdl(
    context: Context,
    x0_nm: float,
    delta_x_nm: float,
    spring_constant_kcal_per_mol_nm2: float,
    temperature_k: float,
) -> float:
    state = context.getState(getPositions=True)
    x_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)[0, 0]
    dvdl = -spring_constant_kcal_per_mol_nm2 * (x_nm - x0_nm) * delta_x_nm
    return dvdl / (K_B_KCAL_PER_MOL_K * temperature_k)


def run_switch(
    config: SwitchingConfig,
    lambda_start: float,
    lambda_end: float,
) -> ar.SwitchingTrajectory:
    temperature_k = config.temperature.value_in_unit(unit.kelvin)
    x_start_nm = config.x_start.value_in_unit(unit.nanometers)
    x_end_nm = config.x_end.value_in_unit(unit.nanometers)
    delta_x_nm = x_end_nm - x_start_nm
    spring_constant = config.spring_constant.value_in_unit(
        unit.kilocalories_per_mole / unit.nanometers**2
    )
    context = make_context(config)
    integrator = context.getIntegrator()

    context.setParameter("x0", x_start_nm + lambda_start * delta_x_nm)
    integrator.step(config.equilibration_steps)

    reduced_work = 0.0
    lambda_path: list[float] = []
    dvdl_path: list[float] = []

    for step in range(config.n_steps_per_switch):
        lam_old = lambda_start + step * (lambda_end - lambda_start) / config.n_steps_per_switch
        lam_new = lambda_start + (step + 1) * (lambda_end - lambda_start) / config.n_steps_per_switch

        x0_old = x_start_nm + lam_old * delta_x_nm
        x0_new = x_start_nm + lam_new * delta_x_nm

        context.setParameter("x0", x0_old)
        u_old = reduced_potential(context, temperature_k)

        context.setParameter("x0", x0_new)
        u_new = reduced_potential(context, temperature_k)
        reduced_work += u_new - u_old

        lambda_path.append(lam_old)
        dvdl_path.append(
            reduced_dvdl(
                context,
                x0_nm=x0_new,
                delta_x_nm=delta_x_nm,
                spring_constant_kcal_per_mol_nm2=spring_constant,
                temperature_k=temperature_k,
            )
        )
        integrator.step(1)

    return ar.SwitchingTrajectory(
        initial_state=ar.StatePoint([lambda_start], temperature_k),
        final_state=ar.StatePoint([lambda_end], temperature_k),
        reduced_work=reduced_work,
        lambda_path=lambda_path,
        dvdl_path=dvdl_path,
    )


def run_switches(
    config: SwitchingConfig,
    lambda_start: float,
    lambda_end: float,
) -> list[ar.SwitchingTrajectory]:
    return [
        run_switch(config, lambda_start=lambda_start, lambda_end=lambda_end)
        for _ in range(config.n_switches)
    ]


def forward_only_example(config: SwitchingConfig) -> None:
    forward = run_switches(config, lambda_start=0.0, lambda_end=1.0)

    result = ar.NES().estimate(forward)
    convergence = ar.analysis.nes_convergence(forward)
    conversion = kt_in_kcal_per_mol(config.temperature.value_in_unit(unit.kelvin))

    print("Forward-only NES")
    print(f"  dG: {result.delta_f * conversion:.6f} kcal/mol")
    print(f"  sigma: {result.uncertainty * conversion:.6f} kcal/mol")
    print("  trajectories:", len(forward))
    print("  convergence points:", len(convergence))


def forward_and_reverse_example(config: SwitchingConfig) -> None:
    forward = run_switches(config, lambda_start=0.0, lambda_end=1.0)
    reverse = run_switches(config, lambda_start=1.0, lambda_end=0.0)

    forward_result = ar.NES().estimate(forward)
    reverse_result = ar.NES().estimate(reverse)
    conversion = kt_in_kcal_per_mol(config.temperature.value_in_unit(unit.kelvin))

    print("Forward and reverse NES")
    print(
        "  forward dG: "
        f"{forward_result.delta_f * conversion:.6f} +/- "
        f"{forward_result.uncertainty * conversion:.6f} kcal/mol"
    )
    print(
        "  reverse dG: "
        f"{reverse_result.delta_f * conversion:.6f} +/- "
        f"{reverse_result.uncertainty * conversion:.6f} kcal/mol"
    )
    print("  forward trajectories:", len(forward))
    print("  reverse trajectories:", len(reverse))


def main() -> None:
    random.seed(0)
    config = SwitchingConfig()
    forward_only_example(config)
    print()
    forward_and_reverse_example(config)


if __name__ == "__main__":
    main()
