use std::fs;
use std::path::Path;

use crate::cli::input::AnalysisSampleCounts;
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
    pub provenance: OutputProvenance,
    pub sample_counts: AnalysisSampleCounts,
}

pub struct OverlapSummary {
    pub scalar: f64,
    pub eigenvalues: Vec<f64>,
}

pub struct OutputProvenance {
    pub estimator: &'static str,
    pub decorrelate: bool,
    pub remove_burnin: usize,
    pub auto_equilibrate: bool,
    pub fast: bool,
    pub conservative: bool,
    pub nskip: usize,
    pub u_nk_observable: Option<&'static str>,
}

struct RenderedScalarResult<'a> {
    delta: f64,
    sigma: Option<f64>,
    from_lambda: f64,
    to_lambda: f64,
    units: &'static str,
    temperature: f64,
    overlap: Option<&'a OverlapSummary>,
    provenance: &'a OutputProvenance,
    sample_counts: AnalysisSampleCounts,
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
    let rendered = RenderedScalarResult {
        delta: convert_value(result.delta, result.units, result.temperature),
        sigma: result
            .sigma
            .map(|value| convert_value(value, result.units, result.temperature)),
        from_lambda: result.from_lambda,
        to_lambda: result.to_lambda,
        units: format_units(result.units),
        temperature: result.temperature,
        overlap: result.overlap.as_ref(),
        provenance: &result.provenance,
        sample_counts: result.sample_counts,
    };

    match format {
        OutputFormat::Text => render_text(&rendered),
        OutputFormat::Json => render_json(&rendered),
        OutputFormat::Csv => render_csv(&rendered),
    }
}

fn render_text(result: &RenderedScalarResult<'_>) -> String {
    let sigma = result
        .sigma
        .map(|value| format!("{value} {}", result.units))
        .unwrap_or_else(|| "none".to_string());
    let mut output = format!(
        "delta_f: {} {}\nuncertainty: {sigma}\nfrom_lambda: {}\nto_lambda: {}\nwindows: {}\nsamples_in: {}\nsamples_after_burnin: {}\nsamples_kept: {}\n",
        result.delta, result.units, result.from_lambda, result.to_lambda
        ,result.sample_counts.windows
        ,result.sample_counts.samples_in
        ,result.sample_counts.samples_after_burnin
        ,result.sample_counts.samples_kept
    );
    if let Some(u_nk_observable) = result.provenance.u_nk_observable {
        output.push_str(&format!("u_nk_observable: {u_nk_observable}\n"));
    }
    if let Some(overlap) = result.overlap {
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

fn render_json(result: &RenderedScalarResult<'_>) -> String {
    let sigma = result
        .sigma
        .filter(|value| value.is_finite())
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string());
    let overlap = result
        .overlap
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
    let provenance = format!(
        "{{\"estimator\":\"{}\",\"temperature_k\":{},\"decorrelate\":{},\"remove_burnin\":{},\"auto_equilibrate\":{},\"fast\":{},\"conservative\":{},\"nskip\":{},\"u_nk_observable\":{},\"windows\":{},\"samples_in\":{},\"samples_after_burnin\":{},\"samples_kept\":{}}}",
        result.provenance.estimator,
        result.temperature,
        result.provenance.decorrelate,
        result.provenance.remove_burnin,
        result.provenance.auto_equilibrate,
        result.provenance.fast,
        result.provenance.conservative,
        result.provenance.nskip,
        result
            .provenance
            .u_nk_observable
            .map(|value| format!("\"{value}\""))
            .unwrap_or_else(|| "null".to_string()),
        result.sample_counts.windows,
        result.sample_counts.samples_in,
        result.sample_counts.samples_after_burnin,
        result.sample_counts.samples_kept
    );
    format!(
        "{{\"delta_f\":{},\"uncertainty\":{sigma},\"from_lambda\":{},\"to_lambda\":{},\"units\":\"{}\",\"overlap\":{overlap},\"provenance\":{provenance}}}\n",
        result.delta, result.from_lambda, result.to_lambda, result.units
    )
}

fn render_csv(result: &RenderedScalarResult<'_>) -> String {
    let sigma = result
        .sigma
        .map(|value| value.to_string())
        .unwrap_or_default();
    let overlap_scalar = result
        .overlap
        .map(|summary| summary.scalar.to_string())
        .unwrap_or_default();
    let overlap_eigenvalues = result
        .overlap
        .map(|summary| {
            summary
                .eigenvalues
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(";")
        })
        .unwrap_or_default();
    let u_nk_observable = result.provenance.u_nk_observable.unwrap_or_default();
    format!(
        "delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues,estimator,temperature_k,decorrelate,remove_burnin,auto_equilibrate,fast,conservative,nskip,u_nk_observable,windows,samples_in,samples_after_burnin,samples_kept\n{delta},{sigma},{from_lambda},{to_lambda},{units},{overlap_scalar},{overlap_eigenvalues},{},{temperature},{},{},{},{},{},{},{},{},{},{},{}\n",
        result.provenance.estimator,
        result.provenance.decorrelate,
        result.provenance.remove_burnin,
        result.provenance.auto_equilibrate,
        result.provenance.fast,
        result.provenance.conservative,
        result.provenance.nskip,
        u_nk_observable,
        result.sample_counts.windows,
        result.sample_counts.samples_in,
        result.sample_counts.samples_after_burnin,
        result.sample_counts.samples_kept,
        delta = result.delta,
        from_lambda = result.from_lambda,
        to_lambda = result.to_lambda,
        units = result.units,
        temperature = result.temperature
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
    use super::{render_scalar_result, OutputProvenance, OverlapSummary, ScalarResult};
    use crate::cli::input::AnalysisSampleCounts;
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
                provenance: OutputProvenance {
                    estimator: "mbar",
                    decorrelate: true,
                    remove_burnin: 5,
                    auto_equilibrate: false,
                    fast: false,
                    conservative: true,
                    nskip: 1,
                    u_nk_observable: Some("de"),
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 15,
                    samples_in: 300,
                    samples_after_burnin: 225,
                    samples_kept: 120,
                },
            },
            OutputFormat::Json,
        );

        assert_eq!(
            output,
            "{\"delta_f\":1.5,\"uncertainty\":0.25,\"from_lambda\":0,\"to_lambda\":1,\"units\":\"kT\",\"overlap\":null,\"provenance\":{\"estimator\":\"mbar\",\"temperature_k\":300,\"decorrelate\":true,\"remove_burnin\":5,\"auto_equilibrate\":false,\"fast\":false,\"conservative\":true,\"nskip\":1,\"u_nk_observable\":\"de\",\"windows\":15,\"samples_in\":300,\"samples_after_burnin\":225,\"samples_kept\":120}}\n"
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
                provenance: OutputProvenance {
                    estimator: "ti",
                    decorrelate: false,
                    remove_burnin: 0,
                    auto_equilibrate: false,
                    fast: false,
                    conservative: true,
                    nskip: 1,
                    u_nk_observable: None,
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 2,
                    samples_in: 40,
                    samples_after_burnin: 40,
                    samples_kept: 20,
                },
            },
            OutputFormat::Csv,
        );

        assert_eq!(
            output,
            "delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues,estimator,temperature_k,decorrelate,remove_burnin,auto_equilibrate,fast,conservative,nskip,u_nk_observable,windows,samples_in,samples_after_burnin,samples_kept\n1.5,,0,1,kT,,,ti,300,false,0,false,false,true,1,,2,40,40,20\n"
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
                provenance: OutputProvenance {
                    estimator: "bar",
                    decorrelate: true,
                    remove_burnin: 0,
                    auto_equilibrate: false,
                    fast: false,
                    conservative: true,
                    nskip: 1,
                    u_nk_observable: Some("de"),
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 15,
                    samples_in: 300,
                    samples_after_burnin: 300,
                    samples_kept: 120,
                },
            },
            OutputFormat::Text,
        );

        assert_eq!(
            output,
            "delta_f: 1.5 kT\nuncertainty: 0.25 kT\nfrom_lambda: 0\nto_lambda: 1\nwindows: 15\nsamples_in: 300\nsamples_after_burnin: 300\nsamples_kept: 120\nu_nk_observable: de\noverlap_scalar: 0.125\noverlap_eigenvalues: 1, 0.875, 0.5\n"
        );
    }
}
