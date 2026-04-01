use std::fs;
use std::io;
use std::path::PathBuf;

use alchemrs::{
    advise_lambda_schedule, AdvisorEstimator, EdgeSeverity, ScheduleAdvice, ScheduleAdvisorOptions,
    SuggestionKind,
};
use serde_json::{json, Map, Value};

use crate::cli::input::{load_windows, AnalysisInputOptions, AnalysisSampleCounts};
use crate::cli::{AdvisorEstimatorArg, OutputFormat};
use crate::CliResult;

pub struct AdviseRunOptions {
    pub estimator: AdvisorEstimatorArg,
    pub overlap_min: f64,
    pub block_cv_min: f64,
    pub n_blocks: usize,
    pub suggest_midpoints: bool,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
}

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: AdviseRunOptions,
) -> CliResult<()> {
    let loaded = load_windows(inputs, input_options)?;
    let advice = advise_lambda_schedule(
        &loaded.windows,
        Some(ScheduleAdvisorOptions {
            estimator: run_options.estimator.into(),
            overlap_min: run_options.overlap_min,
            block_cv_min: run_options.block_cv_min,
            n_blocks: run_options.n_blocks,
            suggest_midpoints: run_options.suggest_midpoints,
        }),
    )?;
    let lambda_components = loaded
        .windows
        .first()
        .and_then(|window| window.lambda_labels().map(|labels| labels.to_vec()));
    let rendered = render_advice(
        &advice,
        loaded.sample_counts,
        &input_options,
        &run_options,
        lambda_components,
    )
    .map_err(io::Error::other)?;

    if let Some(path) = run_options.output_path.as_deref() {
        fs::write(path, rendered)?;
    } else {
        print!("{rendered}");
    }
    Ok(())
}

impl From<AdvisorEstimatorArg> for AdvisorEstimator {
    fn from(value: AdvisorEstimatorArg) -> Self {
        match value {
            AdvisorEstimatorArg::Mbar => AdvisorEstimator::Mbar,
            AdvisorEstimatorArg::Bar => AdvisorEstimator::Bar,
        }
    }
}

fn render_advice(
    advice: &ScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
    lambda_components: Option<Vec<String>>,
) -> Result<String, String> {
    match run_options.output_format {
        OutputFormat::Text => Ok(render_text(
            advice,
            sample_counts,
            input_options,
            run_options,
            lambda_components,
        )),
        OutputFormat::Json => render_json(
            advice,
            sample_counts,
            input_options,
            run_options,
            lambda_components,
        ),
        OutputFormat::Csv => render_csv(advice, run_options),
    }
}

fn render_text(
    advice: &ScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
    lambda_components: Option<Vec<String>>,
) -> String {
    let mut output = format!(
        "advisor_estimator: {}\nwindows: {}\nsamples_in: {}\nsamples_after_burnin: {}\nsamples_kept: {}\noverlap_min_threshold: {}\nblock_cv_min_threshold: {}\nn_blocks: {}\n",
        estimator_name(run_options.estimator),
        sample_counts.windows,
        sample_counts.samples_in,
        sample_counts.samples_after_burnin,
        sample_counts.samples_kept,
        run_options.overlap_min,
        run_options.block_cv_min,
        run_options.n_blocks,
    );
    if let Some(u_nk_observable) = input_options.u_nk_observable_name() {
        output.push_str(&format!("u_nk_observable: {u_nk_observable}\n"));
    }
    if let Some(labels) = lambda_components.as_ref() {
        output.push_str(&format!("lambda_components: {}\n", labels.join(", ")));
    }
    output.push_str("edges:\n");
    for edge in advice.edges() {
        output.push_str(&format!(
            "  - edge {}: {} -> {}, overlap_min={}, delta_f={} kT, uncertainty={}, block_cv={}, severity={}\n",
            edge.edge_index(),
            format_state_text(edge.from_state().lambdas()),
            format_state_text(edge.to_state().lambdas()),
            edge.overlap_min(),
            edge.delta_f(),
            option_string(edge.uncertainty()),
            option_string(edge.block_cv()),
            severity_name(edge.severity()),
        ));
    }
    if advice.suggestions().is_empty() {
        output.push_str("suggestions: none\n");
    } else {
        output.push_str("suggestions:\n");
        for suggestion in advice.suggestions() {
            output.push_str(&format!(
                "  - edge {}: {} ({})",
                suggestion.edge_index(),
                suggestion.reason(),
                suggestion_name(suggestion.kind()),
            ));
            if let Some(state) = suggestion.proposed_state() {
                output.push_str(&format!(
                    ", proposed_lambda={}",
                    format_state_text(state.lambdas())
                ));
            }
            output.push('\n');
        }
    }
    output
}

fn render_json(
    advice: &ScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
    lambda_components: Option<Vec<String>>,
) -> Result<String, String> {
    let edges = advice
        .edges()
        .iter()
        .map(|edge| {
            json!({
                "edge_index": edge.edge_index(),
                "from_lambda": json_state(edge.from_state().lambdas()),
                "to_lambda": json_state(edge.to_state().lambdas()),
                "lambda_delta": edge.lambda_delta(),
                "overlap_forward": edge.overlap_forward(),
                "overlap_reverse": edge.overlap_reverse(),
                "overlap_min": edge.overlap_min(),
                "delta_f": edge.delta_f(),
                "uncertainty": edge.uncertainty(),
                "block_mean": edge.block_mean(),
                "block_stddev": edge.block_stddev(),
                "block_cv": edge.block_cv(),
                "severity": severity_name(edge.severity()),
            })
        })
        .collect::<Vec<_>>();

    let suggestions = advice
        .suggestions()
        .iter()
        .map(|suggestion| {
            json!({
                "kind": suggestion_name(suggestion.kind()),
                "edge_index": suggestion.edge_index(),
                "from_lambda": json_state(suggestion.from_state().lambdas()),
                "to_lambda": json_state(suggestion.to_state().lambdas()),
                "proposed_lambda": suggestion.proposed_state().map(|state| json_state(state.lambdas())).unwrap_or(Value::Null),
                "reason": suggestion.reason(),
            })
        })
        .collect::<Vec<_>>();

    let mut provenance = Map::new();
    provenance.insert(
        "advisor_estimator".to_string(),
        json!(estimator_name(run_options.estimator)),
    );
    provenance.insert("decorrelate".to_string(), json!(input_options.decorrelate));
    provenance.insert(
        "remove_burnin".to_string(),
        json!(input_options.remove_burnin),
    );
    provenance.insert(
        "auto_equilibrate".to_string(),
        json!(input_options.auto_equilibrate),
    );
    provenance.insert("fast".to_string(), json!(input_options.effective_fast()));
    provenance.insert(
        "conservative".to_string(),
        json!(input_options.effective_conservative()),
    );
    provenance.insert("nskip".to_string(), json!(input_options.nskip));
    provenance.insert(
        "u_nk_observable".to_string(),
        input_options
            .u_nk_observable_name()
            .map(|name| json!(name))
            .unwrap_or(Value::Null),
    );
    provenance.insert("overlap_min".to_string(), json!(run_options.overlap_min));
    provenance.insert("block_cv_min".to_string(), json!(run_options.block_cv_min));
    provenance.insert("n_blocks".to_string(), json!(run_options.n_blocks));
    provenance.insert(
        "lambda_components".to_string(),
        lambda_components
            .as_ref()
            .map(|labels| json!(labels))
            .unwrap_or(Value::Null),
    );

    let payload = json!({
        "sample_counts": {
            "windows": sample_counts.windows,
            "samples_in": sample_counts.samples_in,
            "samples_after_burnin": sample_counts.samples_after_burnin,
            "samples_kept": sample_counts.samples_kept,
        },
        "provenance": provenance,
        "edges": edges,
        "suggestions": suggestions,
    });
    let mut output = serde_json::to_string(&payload).map_err(|err| err.to_string())?;
    output.push('\n');
    Ok(output)
}

fn render_csv(advice: &ScheduleAdvice, run_options: &AdviseRunOptions) -> Result<String, String> {
    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .from_writer(Vec::new());
    writer
        .write_record([
            "edge_index",
            "from_lambda",
            "to_lambda",
            "lambda_delta",
            "overlap_forward",
            "overlap_reverse",
            "overlap_min",
            "delta_f",
            "uncertainty",
            "block_mean",
            "block_stddev",
            "block_cv",
            "severity",
            "suggestion_kind",
            "suggested_lambda",
            "suggestion_reason",
            "advisor_estimator",
        ])
        .map_err(|err| err.to_string())?;

    for edge in advice.edges() {
        let suggestion = advice
            .suggestions()
            .iter()
            .find(|suggestion| suggestion.edge_index() == edge.edge_index());
        writer
            .write_record([
                edge.edge_index().to_string(),
                format_state_csv(edge.from_state().lambdas()),
                format_state_csv(edge.to_state().lambdas()),
                format_state_csv(edge.lambda_delta()),
                edge.overlap_forward().to_string(),
                edge.overlap_reverse().to_string(),
                edge.overlap_min().to_string(),
                edge.delta_f().to_string(),
                option_string(edge.uncertainty()),
                option_string(edge.block_mean()),
                option_string(edge.block_stddev()),
                option_string(edge.block_cv()),
                severity_name(edge.severity()).to_string(),
                suggestion
                    .map(|suggestion| suggestion_name(suggestion.kind()).to_string())
                    .unwrap_or_else(|| suggestion_name(SuggestionKind::NoChange).to_string()),
                suggestion
                    .and_then(|suggestion| suggestion.proposed_state())
                    .map(|state| format_state_csv(state.lambdas()))
                    .unwrap_or_default(),
                suggestion
                    .map(|suggestion| suggestion.reason().to_string())
                    .unwrap_or_default(),
                estimator_name(run_options.estimator).to_string(),
            ])
            .map_err(|err| err.to_string())?;
    }

    let bytes = writer.into_inner().map_err(|err| err.to_string())?;
    String::from_utf8(bytes).map_err(|err| err.to_string())
}

fn estimator_name(estimator: AdvisorEstimatorArg) -> &'static str {
    match estimator {
        AdvisorEstimatorArg::Mbar => "mbar",
        AdvisorEstimatorArg::Bar => "bar",
    }
}

fn severity_name(severity: EdgeSeverity) -> &'static str {
    match severity {
        EdgeSeverity::Healthy => "healthy",
        EdgeSeverity::Monitor => "monitor",
        EdgeSeverity::AddSampling => "add_sampling",
        EdgeSeverity::AddWindow => "add_window",
    }
}

fn suggestion_name(kind: SuggestionKind) -> &'static str {
    match kind {
        SuggestionKind::NoChange => "no_change",
        SuggestionKind::ExtendSampling => "extend_sampling",
        SuggestionKind::InsertWindow => "insert_window",
    }
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

#[cfg(test)]
mod tests {
    use alchemrs::{
        advise_lambda_schedule, extract_u_nk, AdvisorEstimator, ScheduleAdvisorOptions,
    };
    use serde_json::Value;

    use super::{render_advice, AdviseRunOptions};
    use crate::cli::input::{AnalysisInputOptions, AnalysisSampleCounts};
    use crate::cli::{AdvisorEstimatorArg, OutputFormat, UNkObservable};

    fn sample_advice() -> alchemrs::ScheduleAdvice {
        let base = env!("CARGO_MANIFEST_DIR");
        let windows = vec![
            extract_u_nk(format!("{base}/fixtures/gromacs/lambda_0.xvg"), 298.0).unwrap(),
            extract_u_nk(format!("{base}/fixtures/gromacs/lambda_1.xvg"), 298.0).unwrap(),
        ];
        advise_lambda_schedule(
            &windows,
            Some(ScheduleAdvisorOptions {
                estimator: AdvisorEstimator::Mbar,
                overlap_min: 0.99,
                block_cv_min: 1.0e30,
                n_blocks: 4,
                suggest_midpoints: true,
            }),
        )
        .unwrap()
    }

    #[test]
    fn renders_schedule_advice_as_json() {
        let output = render_advice(
            &sample_advice(),
            AnalysisSampleCounts {
                windows: 2,
                samples_in: 40,
                samples_after_burnin: 35,
                samples_kept: 20,
            },
            &AnalysisInputOptions {
                temperature: 300.0,
                decorrelate: true,
                remove_burnin: 5,
                auto_equilibrate: false,
                fast: false,
                conservative: true,
                nskip: 1,
                u_nk_observable: Some(UNkObservable::De),
            },
            &AdviseRunOptions {
                estimator: AdvisorEstimatorArg::Mbar,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 4,
                suggest_midpoints: true,
                output_format: OutputFormat::Json,
                output_path: None,
            },
            None,
        )
        .unwrap();

        let value: Value = serde_json::from_str(&output).unwrap();
        assert_eq!(value["edges"][0]["severity"], "add_window");
        assert_eq!(value["suggestions"][0]["kind"], "insert_window");
        assert_eq!(
            value["suggestions"][0]["proposed_lambda"],
            serde_json::json!([0.0, 0.0, 0.025, 0.0, 0.0])
        );
    }
}
