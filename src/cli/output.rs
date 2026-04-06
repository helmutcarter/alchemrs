use std::fs;
use std::io;
use std::path::Path;

use serde_json::{json, Map, Value};

use crate::cli::input::AnalysisSampleCounts;
use crate::cli::{OutputFormat, OutputUnits};

const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;
const K_B_KJ_PER_MOL_K: f64 = 0.00831446261815324;

pub struct ScalarResult {
    pub delta: f64,
    pub sigma: Option<f64>,
    pub from_state: Vec<f64>,
    pub to_state: Vec<f64>,
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
    pub ti_method: Option<&'static str>,
    pub ti_method_reason: Option<String>,
    pub lambda_components: Option<Vec<String>>,
}

struct RenderedScalarResult<'a> {
    delta: f64,
    sigma: Option<f64>,
    from_state: Vec<f64>,
    to_state: Vec<f64>,
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
) -> io::Result<()> {
    let rendered = render_scalar_result(result, format).map_err(io::Error::other)?;
    if let Some(path) = output_path {
        fs::write(path, rendered)
    } else {
        print!("{rendered}");
        Ok(())
    }
}

fn render_scalar_result(result: &ScalarResult, format: OutputFormat) -> Result<String, String> {
    let rendered = RenderedScalarResult {
        delta: convert_value(result.delta, result.units, result.temperature),
        sigma: result
            .sigma
            .map(|value| convert_value(value, result.units, result.temperature)),
        from_state: result.from_state.clone(),
        to_state: result.to_state.clone(),
        units: format_units(result.units),
        temperature: result.temperature,
        overlap: result.overlap.as_ref(),
        provenance: &result.provenance,
        sample_counts: result.sample_counts,
    };

    match format {
        OutputFormat::Text => Ok(render_text(&rendered)),
        OutputFormat::Json => render_json(&rendered),
        OutputFormat::Csv => render_csv(&rendered),
    }
}

fn render_text(result: &RenderedScalarResult<'_>) -> String {
    let sigma = result
        .sigma
        .map(|value| format!("{value} {}", result.units))
        .unwrap_or_else(|| "none".to_string());
    let from_state = format_state_text(&result.from_state);
    let to_state = format_state_text(&result.to_state);
    let mut output = format!(
        "delta_f: {} {}\nuncertainty: {sigma}\nfrom_lambda: {}\nto_lambda: {}\nwindows: {}\nsamples_in: {}\nsamples_after_burnin: {}\nsamples_kept: {}\n",
        result.delta,
        result.units,
        from_state,
        to_state,
        result.sample_counts.windows,
        result.sample_counts.samples_in,
        result.sample_counts.samples_after_burnin,
        result.sample_counts.samples_kept
    );
    if let Some(u_nk_observable) = result.provenance.u_nk_observable {
        output.push_str(&format!("u_nk_observable: {u_nk_observable}\n"));
    }
    if let Some(ti_method) = result.provenance.ti_method {
        output.push_str(&format!("ti_method: {ti_method}\n"));
    }
    if let Some(ti_method_reason) = result.provenance.ti_method_reason.as_ref() {
        output.push_str(&format!("ti_method_reason: {ti_method_reason}\n"));
    }
    if let Some(lambda_components) = result.provenance.lambda_components.as_ref() {
        output.push_str(&format!(
            "lambda_components: {}\n",
            lambda_components.join(", ")
        ));
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

fn render_json(result: &RenderedScalarResult<'_>) -> Result<String, String> {
    let mut provenance = Map::new();
    provenance.insert("estimator".to_string(), json!(result.provenance.estimator));
    provenance.insert("temperature_k".to_string(), json!(result.temperature));
    provenance.insert(
        "decorrelate".to_string(),
        json!(result.provenance.decorrelate),
    );
    provenance.insert(
        "remove_burnin".to_string(),
        json!(result.provenance.remove_burnin),
    );
    provenance.insert(
        "auto_equilibrate".to_string(),
        json!(result.provenance.auto_equilibrate),
    );
    provenance.insert("fast".to_string(), json!(result.provenance.fast));
    provenance.insert(
        "conservative".to_string(),
        json!(result.provenance.conservative),
    );
    provenance.insert("nskip".to_string(), json!(result.provenance.nskip));
    provenance.insert(
        "u_nk_observable".to_string(),
        option_str_value(result.provenance.u_nk_observable),
    );
    provenance.insert(
        "ti_method".to_string(),
        option_str_value(result.provenance.ti_method),
    );
    provenance.insert(
        "ti_method_reason".to_string(),
        result
            .provenance
            .ti_method_reason
            .as_ref()
            .map(|value| json!(value))
            .unwrap_or(Value::Null),
    );
    provenance.insert(
        "lambda_components".to_string(),
        result
            .provenance
            .lambda_components
            .as_ref()
            .map(|labels| json!(labels))
            .unwrap_or(Value::Null),
    );
    provenance.insert("windows".to_string(), json!(result.sample_counts.windows));
    provenance.insert(
        "samples_in".to_string(),
        json!(result.sample_counts.samples_in),
    );
    provenance.insert(
        "samples_after_burnin".to_string(),
        json!(result.sample_counts.samples_after_burnin),
    );
    provenance.insert(
        "samples_kept".to_string(),
        json!(result.sample_counts.samples_kept),
    );

    let mut payload = Map::new();
    payload.insert("delta_f".to_string(), json!(result.delta));
    payload.insert(
        "uncertainty".to_string(),
        result
            .sigma
            .filter(|value| value.is_finite())
            .map(|value| json!(value))
            .unwrap_or(Value::Null),
    );
    payload.insert("from_lambda".to_string(), json_state(&result.from_state));
    payload.insert("to_lambda".to_string(), json_state(&result.to_state));
    payload.insert("units".to_string(), json!(result.units));
    payload.insert(
        "overlap".to_string(),
        result
            .overlap
            .map(|overlap| {
                json!({
                    "scalar": overlap.scalar,
                    "eigenvalues": overlap.eigenvalues,
                })
            })
            .unwrap_or(Value::Null),
    );
    payload.insert("provenance".to_string(), Value::Object(provenance));

    let mut output =
        serde_json::to_string(&Value::Object(payload)).map_err(|err| err.to_string())?;
    output.push('\n');
    Ok(output)
}

fn render_csv(result: &RenderedScalarResult<'_>) -> Result<String, String> {
    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .from_writer(Vec::new());
    writer
        .write_record([
            "delta_f",
            "uncertainty",
            "from_lambda",
            "to_lambda",
            "units",
            "overlap_scalar",
            "overlap_eigenvalues",
            "estimator",
            "temperature_k",
            "decorrelate",
            "remove_burnin",
            "auto_equilibrate",
            "fast",
            "conservative",
            "nskip",
            "u_nk_observable",
            "ti_method",
            "ti_method_reason",
            "lambda_components",
            "windows",
            "samples_in",
            "samples_after_burnin",
            "samples_kept",
        ])
        .map_err(|err| err.to_string())?;
    let row = vec![
        result.delta.to_string(),
        option_string(result.sigma),
        format_state_csv(&result.from_state),
        format_state_csv(&result.to_state),
        result.units.to_owned(),
        option_string(result.overlap.map(|summary| summary.scalar)),
        result
            .overlap
            .map(format_overlap_eigenvalues)
            .unwrap_or_default(),
        result.provenance.estimator.to_string(),
        result.temperature.to_string(),
        result.provenance.decorrelate.to_string(),
        result.provenance.remove_burnin.to_string(),
        result.provenance.auto_equilibrate.to_string(),
        result.provenance.fast.to_string(),
        result.provenance.conservative.to_string(),
        result.provenance.nskip.to_string(),
        result
            .provenance
            .u_nk_observable
            .unwrap_or_default()
            .to_string(),
        result
            .provenance
            .ti_method
            .unwrap_or_default()
            .to_string(),
        result.provenance.ti_method_reason.clone().unwrap_or_default(),
        result
            .provenance
            .lambda_components
            .as_ref()
            .map(|labels| format!("[{}]", labels.join(";")))
            .unwrap_or_default(),
        result.sample_counts.windows.to_string(),
        result.sample_counts.samples_in.to_string(),
        result.sample_counts.samples_after_burnin.to_string(),
        result.sample_counts.samples_kept.to_string(),
    ];
    writer.write_record(row).map_err(|err| err.to_string())?;
    let bytes = writer.into_inner().map_err(|err| err.to_string())?;
    String::from_utf8(bytes).map_err(|err| err.to_string())
}

fn json_state(state: &[f64]) -> Value {
    if state.len() == 1 {
        json!(state[0])
    } else {
        json!(state)
    }
}

fn format_state_text(state: &[f64]) -> String {
    if state.len() == 1 {
        state[0].to_string()
    } else {
        format!(
            "[{}]",
            state
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

fn format_state_csv(state: &[f64]) -> String {
    if state.len() == 1 {
        state[0].to_string()
    } else {
        format!(
            "[{}]",
            state
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(";")
        )
    }
}

fn option_string(value: Option<f64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn option_str_value(value: Option<&str>) -> Value {
    value.map(|value| json!(value)).unwrap_or(Value::Null)
}

fn format_overlap_eigenvalues(summary: &OverlapSummary) -> String {
    summary
        .eigenvalues
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(";")
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
    use serde_json::Value;

    use super::{render_scalar_result, OutputProvenance, OverlapSummary, ScalarResult};
    use crate::cli::input::AnalysisSampleCounts;
    use crate::cli::{OutputFormat, OutputUnits};

    #[test]
    fn renders_scalar_result_as_json() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: Some(0.25),
                from_state: vec![0.0],
                to_state: vec![1.0],
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
                    ti_method: None,
                    ti_method_reason: None,
                    lambda_components: None,
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 15,
                    samples_in: 300,
                    samples_after_burnin: 225,
                    samples_kept: 120,
                },
            },
            OutputFormat::Json,
        )
        .unwrap();

        let value: Value = serde_json::from_str(&output).unwrap();
        assert_eq!(value["delta_f"].as_f64(), Some(1.5));
        assert_eq!(value["uncertainty"].as_f64(), Some(0.25));
        assert_eq!(value["from_lambda"].as_f64(), Some(0.0));
        assert_eq!(value["to_lambda"].as_f64(), Some(1.0));
        assert_eq!(value["units"].as_str(), Some("kT"));
        assert_eq!(value["overlap"], Value::Null);
        assert_eq!(value["provenance"]["estimator"].as_str(), Some("mbar"));
        assert_eq!(value["provenance"]["temperature_k"].as_f64(), Some(300.0));
        assert_eq!(value["provenance"]["decorrelate"].as_bool(), Some(true));
        assert_eq!(value["provenance"]["remove_burnin"].as_u64(), Some(5));
        assert_eq!(
            value["provenance"]["auto_equilibrate"].as_bool(),
            Some(false)
        );
        assert_eq!(value["provenance"]["fast"].as_bool(), Some(false));
        assert_eq!(value["provenance"]["conservative"].as_bool(), Some(true));
        assert_eq!(value["provenance"]["nskip"].as_u64(), Some(1));
        assert_eq!(value["provenance"]["u_nk_observable"].as_str(), Some("de"));
        assert_eq!(value["provenance"]["ti_method"], Value::Null);
        assert_eq!(value["provenance"]["ti_method_reason"], Value::Null);
        assert_eq!(value["provenance"]["lambda_components"], Value::Null);
        assert_eq!(value["provenance"]["windows"].as_u64(), Some(15));
        assert_eq!(value["provenance"]["samples_in"].as_u64(), Some(300));
        assert_eq!(
            value["provenance"]["samples_after_burnin"].as_u64(),
            Some(225)
        );
        assert_eq!(value["provenance"]["samples_kept"].as_u64(), Some(120));
    }

    #[test]
    fn renders_scalar_result_as_csv_without_uncertainty() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: None,
                from_state: vec![0.0],
                to_state: vec![1.0],
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
                    ti_method: Some("trapezoidal"),
                    ti_method_reason: None,
                    lambda_components: None,
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 2,
                    samples_in: 40,
                    samples_after_burnin: 40,
                    samples_kept: 20,
                },
            },
            OutputFormat::Csv,
        )
        .unwrap();

        assert_eq!(
            output,
            "delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues,estimator,temperature_k,decorrelate,remove_burnin,auto_equilibrate,fast,conservative,nskip,u_nk_observable,ti_method,ti_method_reason,lambda_components,windows,samples_in,samples_after_burnin,samples_kept\n1.5,,0,1,kT,,,ti,300,false,0,false,false,true,1,,trapezoidal,,,2,40,40,20\n"
        );
    }

    #[test]
    fn renders_scalar_result_as_text_with_overlap() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: Some(0.25),
                from_state: vec![0.0],
                to_state: vec![1.0],
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
                    ti_method: None,
                    ti_method_reason: None,
                    lambda_components: None,
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 15,
                    samples_in: 300,
                    samples_after_burnin: 300,
                    samples_kept: 120,
                },
            },
            OutputFormat::Text,
        )
        .unwrap();

        assert_eq!(
            output,
            "delta_f: 1.5 kT\nuncertainty: 0.25 kT\nfrom_lambda: 0\nto_lambda: 1\nwindows: 15\nsamples_in: 300\nsamples_after_burnin: 300\nsamples_kept: 120\nu_nk_observable: de\noverlap_scalar: 0.125\noverlap_eigenvalues: 1, 0.875, 0.5\n"
        );
    }

    #[test]
    fn renders_multidimensional_scalar_result_as_json() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: Some(0.25),
                from_state: vec![0.0, 0.0, 0.8],
                to_state: vec![0.0, 0.0, 0.9],
                units: OutputUnits::KT,
                temperature: 300.0,
                overlap: None,
                provenance: OutputProvenance {
                    estimator: "mbar",
                    decorrelate: false,
                    remove_burnin: 0,
                    auto_equilibrate: false,
                    fast: false,
                    conservative: true,
                    nskip: 1,
                    u_nk_observable: Some("all"),
                    ti_method: None,
                    ti_method_reason: None,
                    lambda_components: Some(vec![
                        "mass-lambda".to_string(),
                        "coul-lambda".to_string(),
                        "vdw-lambda".to_string(),
                    ]),
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 3,
                    samples_in: 30,
                    samples_after_burnin: 30,
                    samples_kept: 30,
                },
            },
            OutputFormat::Json,
        )
        .unwrap();

        let value: Value = serde_json::from_str(&output).unwrap();
        assert_eq!(value["from_lambda"], serde_json::json!([0.0, 0.0, 0.8]));
        assert_eq!(value["to_lambda"], serde_json::json!([0.0, 0.0, 0.9]));
        assert_eq!(value["provenance"]["u_nk_observable"].as_str(), Some("all"));
        assert_eq!(value["provenance"]["ti_method"], Value::Null);
        assert_eq!(value["provenance"]["ti_method_reason"], Value::Null);
        assert_eq!(
            value["provenance"]["lambda_components"],
            serde_json::json!(["mass-lambda", "coul-lambda", "vdw-lambda"])
        );
    }

    #[test]
    fn renders_multidimensional_scalar_result_as_text() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: None,
                from_state: vec![0.0, 0.0, 0.8],
                to_state: vec![0.0, 0.0, 0.9],
                units: OutputUnits::KT,
                temperature: 300.0,
                overlap: None,
                provenance: OutputProvenance {
                    estimator: "mbar",
                    decorrelate: false,
                    remove_burnin: 0,
                    auto_equilibrate: false,
                    fast: false,
                    conservative: true,
                    nskip: 1,
                    u_nk_observable: Some("all"),
                    ti_method: None,
                    ti_method_reason: None,
                    lambda_components: Some(vec![
                        "mass-lambda".to_string(),
                        "coul-lambda".to_string(),
                        "vdw-lambda".to_string(),
                    ]),
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 3,
                    samples_in: 30,
                    samples_after_burnin: 30,
                    samples_kept: 30,
                },
            },
            OutputFormat::Text,
        )
        .unwrap();

        assert_eq!(
            output,
            "delta_f: 1.5 kT\nuncertainty: none\nfrom_lambda: [0, 0, 0.8]\nto_lambda: [0, 0, 0.9]\nwindows: 3\nsamples_in: 30\nsamples_after_burnin: 30\nsamples_kept: 30\nu_nk_observable: all\nlambda_components: mass-lambda, coul-lambda, vdw-lambda\n"
        );
    }

    #[test]
    fn renders_multidimensional_scalar_result_as_csv() {
        let output = render_scalar_result(
            &ScalarResult {
                delta: 1.5,
                sigma: None,
                from_state: vec![0.0, 0.0, 0.8],
                to_state: vec![0.0, 0.0, 0.9],
                units: OutputUnits::KT,
                temperature: 300.0,
                overlap: None,
                provenance: OutputProvenance {
                    estimator: "mbar",
                    decorrelate: false,
                    remove_burnin: 0,
                    auto_equilibrate: false,
                    fast: false,
                    conservative: true,
                    nskip: 1,
                    u_nk_observable: Some("all"),
                    ti_method: None,
                    ti_method_reason: None,
                    lambda_components: Some(vec![
                        "mass-lambda".to_string(),
                        "coul-lambda".to_string(),
                        "vdw-lambda".to_string(),
                    ]),
                },
                sample_counts: AnalysisSampleCounts {
                    windows: 3,
                    samples_in: 30,
                    samples_after_burnin: 30,
                    samples_kept: 30,
                },
            },
            OutputFormat::Csv,
        )
        .unwrap();

        assert_eq!(
            output,
            "delta_f,uncertainty,from_lambda,to_lambda,units,overlap_scalar,overlap_eigenvalues,estimator,temperature_k,decorrelate,remove_burnin,auto_equilibrate,fast,conservative,nskip,u_nk_observable,ti_method,ti_method_reason,lambda_components,windows,samples_in,samples_after_burnin,samples_kept\n1.5,,[0;0;0.8],[0;0;0.9],kT,,,mbar,300,false,0,false,false,true,1,all,,,[mass-lambda;coul-lambda;vdw-lambda],3,30,30,30\n"
        );
    }
}
