use std::fs;
use std::path::Path;

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
    pub overlap: Option<OverlapSummary>,
}

pub struct OverlapSummary {
    pub scalar: f64,
    pub eigenvalues: Vec<f64>,
}

pub fn print_scalar_result(
    result: &ScalarResult,
    format: OutputFormat,
    output_path: Option<&Path>,
) -> std::io::Result<()> {
    let rendered = render_scalar_result(result, format);
    if let Some(path) = output_path {
        fs::write(path, rendered)
    } else {
        print!("{rendered}");
        Ok(())
    }
}

fn render_scalar_result(result: &ScalarResult, format: OutputFormat) -> String {
    let delta = convert_value(result.delta, result.units, result.temperature);
    let sigma = result
        .sigma
        .map(|value| convert_value(value, result.units, result.temperature));
    let units = format_units(result.units);

    match format {
        OutputFormat::Text => render_text(
            delta,
            sigma,
            result.from_lambda,
            result.to_lambda,
            units,
            result.overlap.as_ref(),
        ),
        OutputFormat::Json => render_json(
            delta,
            sigma,
            result.from_lambda,
            result.to_lambda,
            units,
            result.overlap.as_ref(),
        ),
        OutputFormat::Csv => render_csv(
            delta,
            sigma,
            result.from_lambda,
            result.to_lambda,
            units,
            result.overlap.as_ref(),
        ),
    }
}

fn render_text(
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: &'static str,
    overlap: Option<&OverlapSummary>,
) -> String {
    let sigma = sigma
        .map(|value| format!("{value} {units}"))
        .unwrap_or_else(|| "none".to_string());
    let mut output = format!(
        "delta_f: {delta} {units}\nuncertainty: {sigma}\nfrom_lambda: {from_lambda}\nto_lambda: {to_lambda}\n"
    );
    if let Some(overlap) = overlap {
        output.push_str(&format!("overlap_scalar: {}\n", overlap.scalar));
        output.push_str("overlap_eigenvalues: ");
        output.push_str(
            &overlap
                .eigenvalues
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(", "),
        );
        output.push('\n');
    }
    output
}

fn render_json(
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: &'static str,
    overlap: Option<&OverlapSummary>,
) -> String {
    let sigma = sigma
        .filter(|value| value.is_finite())
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string());
    let overlap = overlap
        .map(|overlap| {
            let eigenvalues = overlap
                .eigenvalues
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(",");
            format!(
                "{{\"scalar\":{},\"eigenvalues\":[{}]}}",
                overlap.scalar, eigenvalues
            )
        })
        .unwrap_or_else(|| "null".to_string());
    format!(
        "{{\"delta_f\":{delta},\"uncertainty\":{sigma},\"from_lambda\":{from_lambda},\"to_lambda\":{to_lambda},\"units\":\"{units}\",\"overlap\":{overlap}}}\n"
    )
}

fn render_csv(
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: &'static str,
    overlap: Option<&OverlapSummary>,
) -> String {
    let sigma = sigma.map(|value| value.to_string()).unwrap_or_default();
    let overlap_scalar = overlap
        .map(|summary| summary.scalar.to_string())
        .unwrap_or_default();
    let overlap_eigenvalues = overlap
        .map(|summary| {
            summary
                .eigenvalues
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(";")
        })
        .unwrap_or_default();
    format!(
        "delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues\n{delta},{sigma},{from_lambda},{to_lambda},{units},{overlap_scalar},{overlap_eigenvalues}\n"
    )
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
    use super::{render_scalar_result, OverlapSummary, ScalarResult};
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
                overlap: None,
            },
            OutputFormat::Json,
        );

        assert_eq!(
            output,
            "{\"delta_f\":1.5,\"uncertainty\":0.25,\"from_lambda\":0,\"to_lambda\":1,\"units\":\"kT\",\"overlap\":null}\n"
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
                overlap: None,
            },
            OutputFormat::Csv,
        );

        assert_eq!(
            output,
            "delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues\n1.5,,0,1,kT,,\n"
        );
    }

    #[test]
    fn renders_scalar_result_as_text_with_overlap() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: Some(0.25),
                from_lambda: 0.0,
                to_lambda: 1.0,
                units: OutputUnits::KT,
                temperature: 300.0,
                overlap: Some(OverlapSummary {
                    scalar: 0.125,
                    eigenvalues: vec![1.0, 0.875, 0.5],
                }),
            },
            OutputFormat::Text,
        );

        assert_eq!(
            output,
            "delta_f: 1.5 kT\nuncertainty: 0.25 kT\nfrom_lambda: 0\nto_lambda: 1\noverlap_scalar: 0.125\noverlap_eigenvalues: 1, 0.875, 0.5\n"
        );
    }
}
