from __future__ import annotations

import alchemrs as ar


TEMPERATURE_K = 300.0


def leg_schedule(direction: str):
    return ar.atm.schedule_from_arrays(
        state_ids=[0, 1],
        direction=direction,
        lambda1=[0.0, 1.0],
        lambda2=[0.0, 1.0],
        alpha=[0.0, 0.0],
        u0=[0.0, 0.0],
        w0=[0.0, 0.0],
        temperature_k=[TEMPERATURE_K, TEMPERATURE_K],
    )


def leg_samples(direction: str):
    schedule = leg_schedule(direction)
    return ar.atm.sample_set_from_arrays(
        schedule,
        state_ids=[0, 1, 0, 1],
        potential_energies_kcal_per_mol=[10.0, 12.0, 11.0, 13.0],
        perturbation_energies_kcal_per_mol=[2.0, 2.0, 2.5, 2.5],
    )


def main() -> None:
    estimator = ar.ATM()

    leg1 = leg_samples("forward")
    leg2 = leg_samples("reverse")

    leg1_result = estimator.estimate_leg(leg1)
    leg2_result = estimator.estimate_leg(leg2)
    binding = estimator.estimate_rbfe(leg1, leg2)

    print("ATM leg analysis")
    print(f"  leg1 dG: {leg1_result.delta_f:.6f} kT")
    print(f"  leg2 dG: {leg2_result.delta_f:.6f} kT")
    print(f"  binding dG: {binding.delta_f:.6f} kT")


if __name__ == "__main__":
    main()
