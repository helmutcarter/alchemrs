use crate::cli::OutputUnits;

const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;
const K_B_KJ_PER_MOL_K: f64 = 0.00831446261815324;

pub fn print_scalar_result(
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: OutputUnits,
    temperature: f64,
) {
    println!(
        "delta_f: {} {}",
        convert_value(delta, units, temperature),
        format_units(units)
    );
    match sigma {
        Some(value) => println!(
            "uncertainty: {} {}",
            convert_value(value, units, temperature),
            format_units(units)
        ),
        None => println!("uncertainty: none"),
    }
    println!("from_lambda: {}", from_lambda);
    println!("to_lambda: {}", to_lambda);
}

fn convert_value(value: f64, units: OutputUnits, temperature: f64) -> f64 {
    match units {
        OutputUnits::KT => value,
        OutputUnits::Kcal => value * (K_B_KCAL_PER_MOL_K * temperature),
        OutputUnits::Kj => value * (K_B_KJ_PER_MOL_K * temperature),
    }
}

fn format_units(units: OutputUnits) -> &'static str {
    match units {
        OutputUnits::KT => "kT",
        OutputUnits::Kcal => "kcal/mol",
        OutputUnits::Kj => "kJ/mol",
    }
}
