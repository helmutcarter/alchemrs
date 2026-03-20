use crate::cli::{OutputFormat, OutputUnits};

const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;
const K_B_KJ_PER_MOL_K: f64 = 0.00831446261815324;

pub struct ScalarResult {
    pub delta: f64,
    pub sigma: Option<f64>,
    pub from_lambda: f64,
    pub to_lambda: f64,
    pub units: OutputUnits,
    pub temperature: f64,
}

pub fn print_scalar_result(result: &ScalarResult, format: OutputFormat) {
    print!("{}", render_scalar_result(result, format));
}

fn render_scalar_result(result: &ScalarResult, format: OutputFormat) -> String {
    let delta = convert_value(result.delta, result.units, result.temperature);
    let sigma = result
        .sigma
        .map(|value| convert_value(value, result.units, result.temperature));
    let units = format_units(result.units);

    match format {
        OutputFormat::Text => {
            render_text(delta, sigma, result.from_lambda, result.to_lambda, units)
        }
        OutputFormat::Json => {
            render_json(delta, sigma, result.from_lambda, result.to_lambda, units)
        }
        OutputFormat::Csv => render_csv(delta, sigma, result.from_lambda, result.to_lambda, units),
    }
}

fn render_text(
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: &'static str,
) -> String {
    let sigma = sigma
        .map(|value| format!("{value} {units}"))
        .unwrap_or_else(|| "none".to_string());
    format!(
        "delta_f: {delta} {units}\nuncertainty: {sigma}\nfrom_lambda: {from_lambda}\nto_lambda: {to_lambda}\n"
    )
}

fn render_json(
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: &'static str,
) -> String {
    let sigma = sigma
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string());
    format!(
        "{{\"delta_f\":{delta},\"uncertainty\":{sigma},\"from_lambda\":{from_lambda},\"to_lambda\":{to_lambda},\"units\":\"{units}\"}}\n"
    )
}

fn render_csv(
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: &'static str,
) -> String {
    let sigma = sigma.map(|value| value.to_string()).unwrap_or_default();
    format!("delta_f,uncertainty,from_lambda,to_lambda,units\n{delta},{sigma},{from_lambda},{to_lambda},{units}\n")
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

#[cfg(test)]
mod tests {
    use super::{render_scalar_result, ScalarResult};
    use crate::cli::{OutputFormat, OutputUnits};

    #[test]
    fn renders_scalar_result_as_json() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: Some(0.25),
                from_lambda: 0.0,
                to_lambda: 1.0,
                units: OutputUnits::KT,
                temperature: 300.0,
            },
            OutputFormat::Json,
        );

        assert_eq!(
            output,
            "{\"delta_f\":1.5,\"uncertainty\":0.25,\"from_lambda\":0,\"to_lambda\":1,\"units\":\"kT\"}\n"
        );
    }

    #[test]
    fn renders_scalar_result_as_csv_without_uncertainty() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: None,
                from_lambda: 0.0,
                to_lambda: 1.0,
                units: OutputUnits::KT,
                temperature: 300.0,
            },
            OutputFormat::Csv,
        );

        assert_eq!(
            output,
            "delta_f,uncertainty,from_lambda,to_lambda,units\n1.5,,0,1,kT\n"
        );
    }
}
