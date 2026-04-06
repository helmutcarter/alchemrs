use std::fs;
use std::io;
use std::path::PathBuf;

use alchemrs::{
    advise_lambda_schedule, advise_ti_schedule, AdvisorEstimator, EdgeSeverity, IntegrationMethod,
    ScheduleAdvice, ScheduleAdvisorOptions, SuggestionKind, TiEdgeSeverity, TiIntervalDiagnostic,
    TiScheduleAdvice, TiScheduleAdvisorOptions, TiScheduleSuggestion, TiSuggestionKind,
    TiWindowDiagnostic,
};
use alchemrs::estimators::sample_ti_curve;
use serde_json::{json, Map, Value};

use crate::cli::input::{
    load_dhdl_series, load_windows, AnalysisInputOptions, AnalysisSampleCounts,
};
use crate::cli::{AdviseInputKind, AdvisorEstimatorArg, OutputFormat};
use crate::CliResult;

pub struct AdviseRunOptions {
    pub estimator: AdvisorEstimatorArg,
    pub input_kind: AdviseInputKind,
    pub overlap_min: f64,
    pub block_cv_min: f64,
    pub n_blocks: usize,
    pub suggest_midpoints: bool,
    pub output_format: OutputFormat,
    pub output_path: Option<PathBuf>,
    pub report_path: Option<PathBuf>,
}

enum AdviceKind {
    UNk {
        advice: ScheduleAdvice,
        lambda_components: Option<Vec<String>>,
    },
    Ti {
        advice: TiScheduleAdvice,
    },
}

pub fn run(
    inputs: Vec<PathBuf>,
    input_options: AnalysisInputOptions,
    run_options: AdviseRunOptions,
) -> CliResult<()> {
    let (sample_counts, advice_kind) = match run_options.input_kind {
        AdviseInputKind::UNk => {
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
            let lambda_components = normalize_lambda_components(
                loaded
                    .windows
                    .first()
                    .and_then(|window| window.lambda_labels().map(|labels| labels.to_vec())),
            );
            (
                loaded.sample_counts,
                AdviceKind::UNk {
                    advice,
                    lambda_components,
                },
            )
        }
        AdviseInputKind::Dhdl => {
            let loaded = load_dhdl_series(inputs, input_options)?;
            let advice = advise_ti_schedule(
                &loaded.series,
                Some(TiScheduleAdvisorOptions {
                    n_blocks: run_options.n_blocks,
                    block_cv_min: run_options.block_cv_min,
                    suggest_midpoints: run_options.suggest_midpoints,
                    ..TiScheduleAdvisorOptions::default()
                }),
            )?;
            (loaded.sample_counts, AdviceKind::Ti { advice })
        }
        AdviseInputKind::Auto => {
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
            let lambda_components = normalize_lambda_components(
                loaded
                    .windows
                    .first()
                    .and_then(|window| window.lambda_labels().map(|labels| labels.to_vec())),
            );
            (
                loaded.sample_counts,
                AdviceKind::UNk {
                    advice,
                    lambda_components,
                },
            )
        }
    };
    let rendered =
        render_schedule_advice(&advice_kind, sample_counts, &input_options, &run_options)
            .map_err(io::Error::other)?;

    if let Some(path) = run_options.report_path.as_deref() {
        let report =
            render_schedule_html_report(&advice_kind, sample_counts, &input_options, &run_options)
                .map_err(io::Error::other)?;
        fs::write(path, report)?;
    }

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

fn render_schedule_advice(
    advice: &AdviceKind,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
) -> Result<String, String> {
    match advice {
        AdviceKind::UNk {
            advice,
            lambda_components,
        } => render_advice(
            advice,
            sample_counts,
            input_options,
            run_options,
            lambda_components.clone(),
        ),
        AdviceKind::Ti { advice } => {
            render_ti_advice(advice, sample_counts, input_options, run_options)
        }
    }
}

fn render_schedule_html_report(
    advice: &AdviceKind,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
) -> Result<String, String> {
    match advice {
        AdviceKind::UNk {
            advice,
            lambda_components,
        } => render_html_report(
            advice,
            sample_counts,
            input_options,
            run_options,
            lambda_components.clone(),
        ),
        AdviceKind::Ti { advice } => {
            render_ti_html_report(advice, sample_counts, input_options, run_options)
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

fn render_ti_advice(
    advice: &TiScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
) -> Result<String, String> {
    match run_options.output_format {
        OutputFormat::Text => Ok(render_ti_text(
            advice,
            sample_counts,
            input_options,
            run_options,
        )),
        OutputFormat::Json => render_ti_json(advice, sample_counts, input_options, run_options),
        OutputFormat::Csv => render_ti_csv(advice, run_options),
    }
}

fn render_ti_text(
    advice: &TiScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    _input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
) -> String {
    let mut output = format!(
        "advisor_mode: ti\nwindows: {}\nsamples_in: {}\nsamples_after_burnin: {}\nsamples_kept: {}\nblock_cv_min_threshold: {}\nn_blocks: {}\nintervals:\n",
        sample_counts.windows,
        sample_counts.samples_in,
        sample_counts.samples_after_burnin,
        sample_counts.samples_kept,
        run_options.block_cv_min,
        run_options.n_blocks,
    );
    for interval in advice.intervals() {
        output.push_str(&format!(
            "  - interval {}: {} -> {}, delta_lambda={}, slope={}, curvature={}, interval_uncertainty={}, severity={}, priority_score={}\n",
            interval.interval_index(),
            format_state_text(interval.from_state().lambdas()),
            format_state_text(interval.to_state().lambdas()),
            interval.delta_lambda(),
            interval.slope(),
            option_string(interval.curvature()),
            option_string(interval.interval_uncertainty()),
            ti_severity_name(interval.severity()),
            interval.priority_score(),
        ));
    }
    if advice.suggestions().is_empty() {
        output.push_str("suggestions: none\n");
    } else {
        output.push_str("suggestions:\n");
        for suggestion in advice.suggestions() {
            output.push_str(&format!(
                "  - interval {}: {} ({})",
                suggestion.interval_index(),
                suggestion.reason(),
                ti_suggestion_name(suggestion.kind()),
            ));
            if let Some(state) = suggestion.proposed_state() {
                output.push_str(&format!(
                    ", proposed_lambda={}",
                    format_state_text(state.lambdas())
                ));
            }
            output.push_str(&format!(", priority_score={}", suggestion.priority_score()));
            output.push('\n');
        }
    }
    output
}

fn render_ti_json(
    advice: &TiScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
) -> Result<String, String> {
    let windows = advice
        .windows()
        .iter()
        .map(|window| {
            json!({
                "window_index": window.window_index(),
                "lambda": json_state(window.state().lambdas()),
                "mean_dhdl": window.mean_dhdl(),
                "sem_dhdl": window.sem_dhdl(),
                "block_mean": window.block_mean(),
                "block_stddev": window.block_stddev(),
                "block_cv": window.block_cv(),
                "forward_reverse_delta": window.forward_reverse_delta(),
            })
        })
        .collect::<Vec<_>>();
    let intervals = advice
        .intervals()
        .iter()
        .map(|interval| {
            json!({
                "interval_index": interval.interval_index(),
                "from_lambda": json_state(interval.from_state().lambdas()),
                "to_lambda": json_state(interval.to_state().lambdas()),
                "delta_lambda": interval.delta_lambda(),
                "left_mean_dhdl": interval.left_mean_dhdl(),
                "right_mean_dhdl": interval.right_mean_dhdl(),
                "trapezoid_contribution": interval.trapezoid_contribution(),
                "slope": interval.slope(),
                "abs_slope": interval.abs_slope(),
                "curvature": interval.curvature(),
                "interval_uncertainty": interval.interval_uncertainty(),
                "severity": ti_severity_name(interval.severity()),
                "priority_score": interval.priority_score(),
            })
        })
        .collect::<Vec<_>>();
    let suggestions = advice
        .suggestions()
        .iter()
        .map(|suggestion| {
            json!({
                "kind": ti_suggestion_name(suggestion.kind()),
                "interval_index": suggestion.interval_index(),
                "from_lambda": json_state(suggestion.from_state().lambdas()),
                "to_lambda": json_state(suggestion.to_state().lambdas()),
                "proposed_lambda": suggestion.proposed_state().map(|state| json_state(state.lambdas())).unwrap_or(Value::Null),
                "priority_score": suggestion.priority_score(),
                "reason": suggestion.reason(),
            })
        })
        .collect::<Vec<_>>();
    let payload = json!({
        "advisor_mode": "ti",
        "sample_counts": {
            "windows": sample_counts.windows,
            "samples_in": sample_counts.samples_in,
            "samples_after_burnin": sample_counts.samples_after_burnin,
            "samples_kept": sample_counts.samples_kept,
        },
        "provenance": {
            "decorrelate": input_options.decorrelate,
            "remove_burnin": input_options.remove_burnin,
            "auto_equilibrate": input_options.auto_equilibrate,
            "fast": input_options.effective_fast(),
            "conservative": input_options.effective_conservative(),
            "nskip": input_options.nskip,
            "input_kind": advise_input_kind_name(run_options.input_kind),
            "n_blocks": run_options.n_blocks,
            "block_cv_min": run_options.block_cv_min,
        },
        "windows": windows,
        "intervals": intervals,
        "suggestions": suggestions,
    });
    let mut output = serde_json::to_string(&payload).map_err(|err| err.to_string())?;
    output.push('\n');
    Ok(output)
}

fn render_ti_csv(
    advice: &TiScheduleAdvice,
    _run_options: &AdviseRunOptions,
) -> Result<String, String> {
    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .from_writer(Vec::new());
    writer
        .write_record([
            "interval_index",
            "from_lambda",
            "to_lambda",
            "delta_lambda",
            "left_mean_dhdl",
            "right_mean_dhdl",
            "trapezoid_contribution",
            "slope",
            "abs_slope",
            "curvature",
            "interval_uncertainty",
            "severity",
            "priority_score",
            "suggestion_kind",
            "suggested_lambda",
            "suggestion_reason",
        ])
        .map_err(|err| err.to_string())?;
    for interval in advice.intervals() {
        let suggestion = advice
            .suggestions()
            .iter()
            .find(|suggestion| suggestion.interval_index() == interval.interval_index());
        writer
            .write_record([
                interval.interval_index().to_string(),
                format_state_csv(interval.from_state().lambdas()),
                format_state_csv(interval.to_state().lambdas()),
                interval.delta_lambda().to_string(),
                interval.left_mean_dhdl().to_string(),
                interval.right_mean_dhdl().to_string(),
                interval.trapezoid_contribution().to_string(),
                interval.slope().to_string(),
                interval.abs_slope().to_string(),
                option_string(interval.curvature()),
                option_string(interval.interval_uncertainty()),
                ti_severity_name(interval.severity()).to_string(),
                interval.priority_score().to_string(),
                suggestion
                    .map(|suggestion| ti_suggestion_name(suggestion.kind()).to_string())
                    .unwrap_or_else(|| ti_suggestion_name(TiSuggestionKind::NoChange).to_string()),
                suggestion
                    .and_then(|suggestion| suggestion.proposed_state())
                    .map(|state| format_state_csv(state.lambdas()))
                    .unwrap_or_default(),
                suggestion
                    .map(|suggestion| suggestion.reason().to_string())
                    .unwrap_or_default(),
            ])
            .map_err(|err| err.to_string())?;
    }
    let bytes = writer.into_inner().map_err(|err| err.to_string())?;
    String::from_utf8(bytes).map_err(|err| err.to_string())
}

fn render_text(
    advice: &ScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
    lambda_components: Option<Vec<String>>,
) -> String {
    let lambda_components = normalize_lambda_components(lambda_components);
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
            "  - edge {}: {} -> {}, overlap_min={}, relative_overlap={}, delta_f={} kT, uncertainty={}, relative_uncertainty={}, block_cv={}, severity={}, dominant_components={}, priority_score={}\n",
            edge.edge_index(),
            format_state_text(edge.from_state().lambdas()),
            format_state_text(edge.to_state().lambdas()),
            edge.overlap_min(),
            option_string(edge.relative_overlap()),
            edge.delta_f(),
            option_string(edge.uncertainty()),
            option_string(edge.relative_uncertainty()),
            option_string(edge.block_cv()),
            severity_name(edge.severity()),
            format_string_list(edge.dominant_components()),
            edge.priority_score(),
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
            if !suggestion.focus_components().is_empty() {
                output.push_str(&format!(
                    ", focus_components={}",
                    format_string_list(suggestion.focus_components())
                ));
            }
            if let Some(strategy) = suggestion.proposal_strategy() {
                output.push_str(&format!(
                    ", proposal_strategy={}",
                    proposal_strategy_name(strategy)
                ));
            }
            if let Some(state) = suggestion.proposed_state() {
                output.push_str(&format!(
                    ", proposed_lambda={}",
                    format_state_text(state.lambdas())
                ));
            }
            output.push_str(&format!(", priority_score={}", suggestion.priority_score()));
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
    let lambda_components = normalize_lambda_components(lambda_components);
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
                "relative_overlap": edge.relative_overlap(),
                "relative_uncertainty": edge.relative_uncertainty(),
                "dominant_components": edge.dominant_components(),
                "priority_score": edge.priority_score(),
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
                "focus_components": suggestion.focus_components(),
                "proposal_strategy": suggestion.proposal_strategy().map(proposal_strategy_name),
                "priority_score": suggestion.priority_score(),
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
            "relative_overlap",
            "delta_f",
            "uncertainty",
            "relative_uncertainty",
            "block_mean",
            "block_stddev",
            "block_cv",
            "dominant_components",
            "priority_score",
            "severity",
            "suggestion_kind",
            "focus_components",
            "proposal_strategy",
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
                option_string(edge.relative_overlap()),
                edge.delta_f().to_string(),
                option_string(edge.uncertainty()),
                option_string(edge.relative_uncertainty()),
                option_string(edge.block_mean()),
                option_string(edge.block_stddev()),
                option_string(edge.block_cv()),
                format_string_list(edge.dominant_components()),
                edge.priority_score().to_string(),
                severity_name(edge.severity()).to_string(),
                suggestion
                    .map(|suggestion| suggestion_name(suggestion.kind()).to_string())
                    .unwrap_or_else(|| suggestion_name(SuggestionKind::NoChange).to_string()),
                suggestion
                    .map(|suggestion| format_string_list(suggestion.focus_components()))
                    .unwrap_or_default(),
                suggestion
                    .and_then(|suggestion| suggestion.proposal_strategy())
                    .map(proposal_strategy_name)
                    .unwrap_or_default()
                    .to_string(),
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

fn render_html_report(
    advice: &ScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
    lambda_components: Option<Vec<String>>,
) -> Result<String, String> {
    let lambda_components = normalize_lambda_components(lambda_components);
    let max_priority = advice
        .edges()
        .iter()
        .map(|edge| edge.priority_score())
        .fold(0.0, f64::max)
        .max(
            advice
                .suggestions()
                .iter()
                .map(|suggestion| suggestion.priority_score())
                .fold(0.0, f64::max),
        )
        .max(1.0);

    let mut html = String::from(
        "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\
<title>alchemrs schedule advisor report</title><style>\
:root{--bg:#f7f4ec;--panel:#fffdf8;--ink:#1d1b19;--muted:#6a655e;--line:#d7cfc0;--good:#2f7d4a;--warn:#b9832f;--bad:#b5483d;--accent:#1f5e5b;}\
*{box-sizing:border-box}body{margin:0;font-family:Georgia,\"Times New Roman\",serif;background:radial-gradient(circle at top,#fffaf0,var(--bg));color:var(--ink)}\
.wrap{max-width:1180px;margin:0 auto;padding:32px 20px 56px}.hero{display:grid;gap:14px;margin-bottom:24px}.eyebrow{font:600 12px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.12em;text-transform:uppercase;color:var(--accent)}\
h1{margin:0;font-size:40px;line-height:1}.lede{max-width:72ch;color:var(--muted);font-size:16px;line-height:1.5}.grid{display:grid;gap:16px}.summary{grid-template-columns:repeat(auto-fit,minmax(180px,1fr));margin-bottom:20px}\
.card{background:var(--panel);border:1px solid var(--line);border-radius:18px;padding:16px 18px;box-shadow:0 10px 30px rgba(35,27,10,.06)}.label{font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}\
.value{margin-top:8px;font-size:28px;line-height:1.1}.sub{margin-top:6px;color:var(--muted);font-size:13px;line-height:1.4}.section{margin-top:28px}.section h2{margin:0 0 12px;font-size:22px}\
.stack{display:grid;gap:14px}.pill{display:inline-block;padding:4px 10px;border-radius:999px;font:600 12px/1.2 ui-monospace,Consolas,monospace;text-transform:uppercase;letter-spacing:.06em;background:#efe7d6;color:var(--ink)}\
.pill.healthy{background:rgba(47,125,74,.12);color:var(--good)}.pill.monitor{background:rgba(185,131,47,.14);color:var(--warn)}.pill.add_sampling,.pill.add_window{background:rgba(181,72,61,.12);color:var(--bad)}\
.badge-row{display:flex;flex-wrap:wrap;gap:8px;margin-top:10px}.badge{display:inline-block;padding:3px 8px;border-radius:999px;background:#efe7d6;color:var(--ink);font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.04em}.badge.rationale{background:rgba(31,94,91,.10);color:var(--accent)}\
.metric-row{display:grid;grid-template-columns:160px 1fr auto;gap:12px;align-items:center;margin-top:10px}.metric-name{font-size:13px;color:var(--muted)}\
.bar{height:10px;border-radius:999px;background:#eee6d8;overflow:hidden}.fill{height:100%;background:linear-gradient(90deg,var(--accent),#8db8a3)}.fill.bad{background:linear-gradient(90deg,#c46d52,var(--bad))}\
.mono{font:500 13px/1.45 ui-monospace,Consolas,monospace}.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px;margin-top:14px}.kv > div{padding-top:10px;border-top:1px solid var(--line)}\
.muted{color:var(--muted)}.empty{color:var(--muted);font-style:italic}.table{display:grid;gap:12px}.suggestion-title{display:none}.footer{margin-top:28px;color:var(--muted);font-size:13px}\
.viz{margin-top:14px;padding:12px;border:1px solid var(--line);border-radius:14px;background:linear-gradient(180deg,rgba(255,255,255,.7),rgba(239,231,214,.45))}.viz svg{width:100%;height:auto;display:block}\
.legend-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:12px}.legend-item{padding:12px;border:1px solid rgba(215,207,192,.9);border-radius:14px;background:rgba(255,253,248,.88)}.legend-title{font:600 12px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.06em;text-transform:uppercase;color:var(--muted)}\
.legend-row{display:flex;align-items:center;gap:10px;margin-top:8px}.swatch{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;flex:0 0 16px}.swatch.circle{border-radius:999px}.swatch.source{background:#1f5e5b}.swatch.target{background:#b5483d}.swatch.proposal{width:12px;height:12px;background:#b9832f;transform:rotate(45deg);border-radius:2px}.swatch.track{width:54px;height:10px;border-radius:999px;background:#eee6d8;overflow:hidden}.swatch.track::before{content:\"\";display:block;height:100%;width:75%;background:linear-gradient(90deg,#d8b56b,var(--warn))}.swatch.track.focus::before{background:linear-gradient(90deg,#5d9b98,var(--accent))}\
.queue-table{width:100%;border-collapse:collapse;margin-top:12px}.queue-table th,.queue-table td{padding:10px 12px;border-top:1px solid var(--line);text-align:left;vertical-align:top}.queue-table th{font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}.queue-table td{font-size:13px}.queue-table .mono{font-size:12px}.queue-link{color:var(--accent);text-decoration:none;border-bottom:1px solid rgba(31,94,91,.25)}.queue-link:hover{border-bottom-color:var(--accent)}\
.component-grid{display:grid;gap:8px;margin-top:12px}.component-row{display:grid;grid-template-columns:minmax(120px,1.3fr) repeat(3,minmax(72px,1fr)) minmax(140px,1.4fr) minmax(80px,1fr);gap:10px;align-items:center;padding:10px 12px;border:1px solid rgba(215,207,192,.9);border-radius:12px;background:rgba(255,253,248,.88)}\
.component-row.focus{border-color:rgba(31,94,91,.4);background:rgba(31,94,91,.06)}.component-head{font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}\
.component-name{font:600 13px/1.25 ui-monospace,Consolas,monospace}.component-value{font:500 12px/1.25 ui-monospace,Consolas,monospace;color:var(--ink)}.status{display:inline-block;padding:3px 8px;border-radius:999px;font:600 10px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.06em;text-transform:uppercase}\
.status.fixed{background:rgba(106,101,94,.12);color:var(--muted)}.status.changed{background:rgba(185,131,47,.14);color:var(--warn)}.status.dominant,.status.split{background:rgba(31,94,91,.12);color:var(--accent)}.status.held{background:rgba(181,72,61,.12);color:var(--bad)}\
.delta-track{height:10px;border-radius:999px;background:#eee6d8;overflow:hidden;position:relative}.delta-fill{height:100%;background:linear-gradient(90deg,#d8b56b,var(--warn))}.delta-fill.focus{background:linear-gradient(90deg,#5d9b98,var(--accent))}.delta-value{margin-top:4px;font:500 11px/1.2 ui-monospace,Consolas,monospace;color:var(--muted)}\
</style></head><body><div class=\"wrap\">",
    );

    html.push_str(
        "<header class=\"hero\"><div class=\"eyebrow\">alchemrs schedule advisor</div><h1>Lambda Schedule Report</h1>\
<div class=\"lede\">Threshold-driven local schedule diagnostics with neighbor-aware ranking and component-aware insertion proposals.</div></header>",
    );

    html.push_str("<section class=\"grid summary\">");
    html.push_str(&summary_card(
        "Estimator",
        estimator_name(run_options.estimator),
        "local fit used by this advisor",
    ));
    html.push_str(&summary_card(
        "Windows",
        &sample_counts.windows.to_string(),
        "input windows after loading",
    ));
    html.push_str(&summary_card(
        "Samples Kept",
        &sample_counts.samples_kept.to_string(),
        "after burn-in trimming and preprocessing",
    ));
    html.push_str(&summary_card(
        "Top Priority",
        &format!("{:.3}", max_priority),
        "maximum advisor priority score",
    ));
    html.push_str("</section>");

    html.push_str(
        "<section class=\"section\"><div class=\"card\"><h2>Configuration</h2><div class=\"kv\">",
    );
    html.push_str(&kv_html(
        "u_nk observable",
        input_options.u_nk_observable_name().unwrap_or(""),
    ));
    html.push_str(&kv_html(
        "decorrelate",
        &input_options.decorrelate.to_string(),
    ));
    html.push_str(&kv_html(
        "auto equilibrate",
        &input_options.auto_equilibrate.to_string(),
    ));
    html.push_str(&kv_html(
        "overlap min",
        &format_report_number(run_options.overlap_min),
    ));
    html.push_str(&kv_html(
        "block CV min",
        &format_report_number(run_options.block_cv_min),
    ));
    html.push_str(&kv_html("n blocks", &run_options.n_blocks.to_string()));
    if let Some(labels) = lambda_components
        .as_ref()
        .filter(|labels| !labels.is_empty())
    {
        html.push_str(&kv_html("lambda components", &labels.join(", ")));
    }
    html.push_str("</div></div></section>");

    let mut ranked_edges = advice.edges().iter().collect::<Vec<_>>();
    ranked_edges.sort_by(|left, right| {
        right
            .priority_score()
            .partial_cmp(&left.priority_score())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    html.push_str(
        "<section class=\"section\"><div class=\"card\"><h2>Priority Queue</h2><div class=\"sub\">Top-ranked edges sorted by advisor priority score for quick triage.</div><table class=\"queue-table\"><thead><tr><th>Edge</th><th>Priority</th><th>Severity</th><th>Weakest Component</th><th>Overlap</th><th>Endpoints</th></tr></thead><tbody>",
    );
    for edge in ranked_edges.into_iter().take(5) {
        html.push_str(&format!(
            "<tr><td class=\"mono\"><a class=\"queue-link\" href=\"#edge-{}\">edge {}</a></td><td class=\"mono\">{:.3}</td><td><span class=\"pill {}\">{}</span></td><td class=\"mono\">{}</td><td class=\"mono\">{}</td><td class=\"mono\">{} -&gt; {}</td></tr>",
            edge.edge_index(),
            edge.edge_index(),
            edge.priority_score(),
            severity_name(edge.severity()),
            escape_html(severity_name(edge.severity())),
            if edge.dominant_components().is_empty() {
                "&mdash;".to_string()
            } else {
                escape_html(&edge.dominant_components().join(", "))
            },
            format_report_number(edge.overlap_min()),
            escape_html(&format_report_state(edge.from_state().lambdas())),
            escape_html(&format_report_state(edge.to_state().lambdas()))
        ));
    }
    html.push_str("</tbody></table></div></section>");

    html.push_str("<section class=\"section\"><h2>Suggestions</h2><div class=\"stack\">");
    if advice.suggestions().is_empty() {
        html.push_str("<div class=\"card empty\">No schedule changes recommended for the current thresholds.</div>");
    } else {
        for suggestion in advice.suggestions() {
            let width = (suggestion.priority_score() / max_priority * 100.0).clamp(0.0, 100.0);
            html.push_str("<article class=\"card\">");
            html.push_str("<div>");
            html.push_str(&format!(
                "<span class=\"pill {}\">{}</span><a class=\"queue-link mono\" href=\"#edge-{}\"> edge {}: &lambda; {} -&gt; {}</a>",
                severity_class_from_suggestion(suggestion.kind()),
                escape_html(suggestion_name(suggestion.kind())),
                suggestion.edge_index(),
                suggestion.edge_index(),
                escape_html(&format_report_state(suggestion.from_state().lambdas())),
                escape_html(&format_report_state(suggestion.to_state().lambdas()))
            ));
            if let Some(strategy) = suggestion.proposal_strategy() {
                html.push_str(&format!(
                    "<span class=\"pill\">{}</span>",
                    escape_html(proposal_strategy_name(strategy))
                ));
            }
            html.push_str("</div>");
            html.push_str(&render_u_nk_suggestion_rationale_badges_html(
                advice,
                suggestion,
                run_options.overlap_min,
                run_options.block_cv_min,
            ));
            html.push_str(&format!(
                "<p>{}</p><div class=\"metric-row\"><div class=\"metric-name\">priority</div><div class=\"bar\"><div class=\"fill bad\" style=\"width:{:.2}%\"></div></div><div class=\"mono\">{:.3}</div></div>",
                escape_html(suggestion.reason()),
                width,
                suggestion.priority_score()
            ));
            html.push_str("<div class=\"kv\">");
            if !suggestion.focus_components().is_empty() {
                html.push_str(&kv_html(
                    "focus components",
                    &suggestion.focus_components().join(", "),
                ));
            }
            if let Some(state) = suggestion.proposed_state() {
                html.push_str(&kv_html(
                    "proposed lambda",
                    &format_report_state(state.lambdas()),
                ));
            }
            if let Some(edge) = advice.edges().get(suggestion.edge_index()) {
                html.push_str(&kv_html(
                    "overlap min",
                    &format_report_number(edge.overlap_min()),
                ));
                html.push_str(&kv_html("block CV", &report_option_string(edge.block_cv())));
                html.push_str(&kv_html(
                    "relative uncertainty",
                    &report_option_string(edge.relative_uncertainty()),
                ));
            }
            html.push_str("</div></article>");
        }
    }
    html.push_str("</div></section>");

    html.push_str("<section class=\"section\"><h2>Edges</h2><div class=\"table\">");
    for edge in advice.edges() {
        let width = (edge.priority_score() / max_priority * 100.0).clamp(0.0, 100.0);
        html.push_str(&format!(
            "<article class=\"card\" id=\"edge-{}\">",
            edge.edge_index()
        ));
        html.push_str(&format!(
            "<div><span class=\"pill {}\">{}</span><span class=\"mono\"> edge {}: &lambda; {} -&gt; {}</span></div>",
            severity_name(edge.severity()),
            escape_html(severity_name(edge.severity())),
            edge.edge_index(),
            escape_html(&format_report_state(edge.from_state().lambdas())),
            escape_html(&format_report_state(edge.to_state().lambdas()))
        ));
        html.push_str(&format!(
            "<div class=\"suggestion-title\"><span class=\"pill {}\">{}</span><span class=\"mono\">edge {}</span><span class=\"mono\">{} → {}</span></div>",
            severity_name(edge.severity()),
            escape_html(severity_name(edge.severity())),
            edge.edge_index(),
            escape_html(&format_report_state(edge.from_state().lambdas())),
            escape_html(&format_report_state(edge.to_state().lambdas()))
        ));
        html.push_str(&format!(
            "<div class=\"metric-row\"><div class=\"metric-name\">priority</div><div class=\"bar\"><div class=\"fill {}\" style=\"width:{:.2}%\"></div></div><div class=\"mono\">{:.3}</div></div>",
            if matches!(edge.severity(), EdgeSeverity::AddSampling | EdgeSeverity::AddWindow) { "bad" } else { "" },
            width,
            edge.priority_score()
        ));
        html.push_str("<div class=\"kv\">");
        html.push_str(&kv_html(
            "overlap min",
            &format_report_number(edge.overlap_min()),
        ));
        html.push_str(&kv_html(
            "relative overlap",
            &report_option_string(edge.relative_overlap()),
        ));
        html.push_str(&kv_html(
            "delta_f (kT)",
            &format_report_number(edge.delta_f()),
        ));
        html.push_str(&kv_html(
            "uncertainty",
            &report_option_string(edge.uncertainty()),
        ));
        html.push_str(&kv_html(
            "relative uncertainty",
            &report_option_string(edge.relative_uncertainty()),
        ));
        html.push_str(&kv_html("block CV", &report_option_string(edge.block_cv())));
        html.push_str(&kv_html(
            "dominant components",
            &edge.dominant_components().join(", "),
        ));
        html.push_str("</div></article>");
    }
    html.push_str("</div></section>");

    html.push_str(
        "<div class=\"footer\">Report generated by <span class=\"mono\">alchemrs advise-schedule</span>. Use the JSON output for exact machine-readable values.</div></div></body></html>",
    );
    Ok(html)
}

fn render_ti_html_report(
    advice: &TiScheduleAdvice,
    sample_counts: AnalysisSampleCounts,
    input_options: &AnalysisInputOptions,
    run_options: &AdviseRunOptions,
) -> Result<String, String> {
    let max_priority = advice
        .intervals()
        .iter()
        .map(|interval| interval.priority_score())
        .fold(0.0, f64::max)
        .max(
            advice
                .suggestions()
                .iter()
                .map(|suggestion| suggestion.priority_score())
                .fold(0.0, f64::max),
        )
        .max(1.0);
    let mut html = String::from(
        "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\
<title>alchemrs TI schedule advisor report</title><style>\
:root{--bg:#f7f4ec;--panel:#fffdf8;--ink:#1d1b19;--muted:#6a655e;--line:#d7cfc0;--good:#2f7d4a;--warn:#b9832f;--bad:#b5483d;--accent:#1f5e5b;}\
*{box-sizing:border-box}body{margin:0;font-family:Georgia,\"Times New Roman\",serif;background:radial-gradient(circle at top,#fffaf0,var(--bg));color:var(--ink)}\
.wrap{max-width:1180px;margin:0 auto;padding:32px 20px 56px}.hero{display:grid;gap:14px;margin-bottom:24px}.eyebrow{font:600 12px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.12em;text-transform:uppercase;color:var(--accent)}\
h1{margin:0;font-size:40px;line-height:1}.lede{max-width:72ch;color:var(--muted);font-size:16px;line-height:1.5}.grid{display:grid;gap:16px}.summary{grid-template-columns:repeat(auto-fit,minmax(180px,1fr));margin-bottom:20px}\
.plot-grid{grid-template-columns:repeat(auto-fit,minmax(360px,1fr));align-items:start}.plot-stack{display:grid;gap:12px}.plot-title{margin:0;font-size:18px}.plot-sub{margin:0;color:var(--muted);font-size:13px;line-height:1.45}.plot-frame{border:1px solid var(--line);border-radius:14px;background:#fcfaf4;padding:10px}.plot-empty{display:grid;place-items:center;min-height:220px;border:1px dashed var(--line);border-radius:12px;color:var(--muted);font-size:14px}\
.ti-series-plot{display:block;width:100%;height:auto}.axis{stroke:#a79d8e;stroke-width:1}.grid-line{stroke:#e8e1d3;stroke-width:1}.zero-line{stroke:#c9b8a8;stroke-width:1;stroke-dasharray:4 4}.series-fill{stroke:none}.series-fill.positive{fill:rgba(47,125,74,.18)}.series-fill.negative{fill:rgba(181,72,61,.18)}.series-line{fill:none;stroke:var(--accent);stroke-width:2.5}.series-line.curvature{stroke:#b5483d}.series-line.uncertainty{stroke:#b9832f}.series-line.method-trapezoidal{stroke:#1f5e5b}.series-line.method-simpson{stroke:#b5483d}.series-line.method-cubic-spline{stroke:#b9832f}.series-line.method-pchip{stroke:#2f7d4a}.series-line.method-akima{stroke:#6b7280}.series-point{fill:var(--accent);stroke:#fffaf0;stroke-width:1.5}.series-point.curvature{fill:#b5483d}.series-point.uncertainty{fill:#b9832f}.axis-label{fill:var(--muted);font:600 11px/1.2 ui-monospace,Consolas,monospace}.tick-label{fill:var(--muted);font:500 10px/1.2 ui-monospace,Consolas,monospace}\
.legend{display:flex;flex-wrap:wrap;gap:10px;margin-top:12px}.legend-item{display:flex;align-items:center;gap:8px;padding:6px 10px;border:1px solid var(--line);border-radius:999px;background:#fffaf0;color:var(--ink);font:500 12px/1.35 ui-monospace,Consolas,monospace}.legend-swatch{width:14px;height:3px;border-radius:999px;background:var(--accent);display:inline-block}.legend-swatch.method-trapezoidal{background:#1f5e5b}.legend-swatch.method-simpson{background:#b5483d}.legend-swatch.method-cubic-spline{background:#b9832f}.legend-swatch.method-pchip{background:#2f7d4a}.legend-swatch.method-akima{background:#6b7280}\
.card{background:var(--panel);border:1px solid var(--line);border-radius:18px;padding:16px 18px;box-shadow:0 10px 30px rgba(35,27,10,.06)}.label{font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}\
.value{margin-top:8px;font-size:28px;line-height:1.1}.sub{margin-top:6px;color:var(--muted);font-size:13px;line-height:1.4}.section{margin-top:28px}.section h2{margin:0 0 12px;font-size:22px}.stack{display:grid;gap:14px}\
.disclosure{padding:0;overflow:hidden}.disclosure-summary{display:flex;align-items:center;justify-content:space-between;gap:16px;list-style:none;cursor:pointer;padding:16px 18px;font-size:18px;font-weight:600}.disclosure-summary::-webkit-details-marker{display:none}.disclosure-summary::after{content:\"Show\";font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--accent)}.disclosure[open] .disclosure-summary::after{content:\"Hide\"}.disclosure-sub{margin:0;padding:0 18px 16px;color:var(--muted);font-size:13px;line-height:1.45}.disclosure-body{padding:0 18px 18px}\
.pill{display:inline-block;padding:4px 10px;border-radius:999px;font:600 12px/1.2 ui-monospace,Consolas,monospace;text-transform:uppercase;letter-spacing:.06em;background:#efe7d6;color:var(--ink)}\
.pill.healthy{background:rgba(47,125,74,.12);color:var(--good)}.pill.monitor{background:rgba(185,131,47,.14);color:var(--warn)}.pill.add_sampling,.pill.add_window,.pill.add_window_and_sampling{background:rgba(181,72,61,.12);color:var(--bad)}\
.badge-row{display:flex;flex-wrap:wrap;gap:8px;margin-top:10px}.badge{display:inline-block;padding:3px 8px;border-radius:999px;background:#efe7d6;color:var(--ink);font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.04em}.badge.rationale{background:rgba(31,94,91,.10);color:var(--accent)}\
.metric-row{display:grid;grid-template-columns:160px 1fr auto;gap:12px;align-items:center;margin-top:10px}.metric-name{font-size:13px;color:var(--muted)}.bar{height:10px;border-radius:999px;background:#eee6d8;overflow:hidden}.fill{height:100%;background:linear-gradient(90deg,var(--accent),#8db8a3)}.fill.bad{background:linear-gradient(90deg,#c46d52,var(--bad))}\
.mono{font:500 13px/1.45 ui-monospace,Consolas,monospace}.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px;margin-top:14px}.kv > div{padding-top:10px;border-top:1px solid var(--line)}\
.queue-table{width:100%;border-collapse:collapse;margin-top:12px}.queue-table th,.queue-table td{padding:10px 12px;border-top:1px solid var(--line);text-align:left;vertical-align:top}.queue-table th{font:600 11px/1.2 ui-monospace,Consolas,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}.queue-link{color:var(--accent);text-decoration:none;border-bottom:1px solid rgba(31,94,91,.25)}.queue-link:hover{border-bottom-color:var(--accent)}\
</style></head><body><div class=\"wrap\">",
    );
    html.push_str("<header class=\"hero\"><div class=\"eyebrow\">alchemrs schedule advisor</div><h1>TI Schedule Report</h1><div class=\"lede\">TI-native spacing diagnostics built from dH/dλ means, block stability, local slope, and curvature.</div></header>");
    html.push_str("<section class=\"grid summary\">");
    html.push_str(&summary_card(
        "Mode",
        "ti",
        "dH/dλ schedule diagnostics",
    ));
    html.push_str(&summary_card(
        "Windows",
        &sample_counts.windows.to_string(),
        "input windows after loading",
    ));
    html.push_str(&summary_card(
        "Samples Kept",
        &sample_counts.samples_kept.to_string(),
        "after burn-in trimming and preprocessing",
    ));
    html.push_str(&summary_card(
        "Top Priority",
        &format!("{max_priority:.3}"),
        "maximum advisor priority score",
    ));
    html.push_str("</section>");
    html.push_str(
        "<section class=\"section\"><div class=\"card\"><h2>Configuration</h2><div class=\"kv\">",
    );
    html.push_str(&kv_html(
        "decorrelate",
        &input_options.decorrelate.to_string(),
    ));
    html.push_str(&kv_html(
        "auto equilibrate",
        &input_options.auto_equilibrate.to_string(),
    ));
    html.push_str(&kv_html("n blocks", &run_options.n_blocks.to_string()));
    html.push_str(&kv_html(
        "block CV min",
        &format_report_number(run_options.block_cv_min),
    ));
    html.push_str("</div></div></section>");
    html.push_str("<section class=\"section\"><h2>Plots</h2><div class=\"grid plot-grid\">");
    html.push_str(&plot_card_html(
        "Mean dH/dλ",
        "Window means after preprocessing. Use this to spot steep regions and sign changes.",
        &render_ti_window_plot_svg(advice),
    ));
    html.push_str(&plot_card_html(
        "Curvature Magnitude",
        "Interval-local curvature diagnostic plotted at interval midpoints. Peaks indicate under-resolved TI regions.",
        &render_ti_curvature_plot_svg(advice),
    ));
    html.push_str(&plot_card_html(
        "Interval Uncertainty",
        "Propagated trapezoid uncertainty per interval midpoint. Peaks indicate where more TI sampling is likely needed.",
        &render_ti_interval_uncertainty_plot_svg(advice),
    ));
    html.push_str("</div></section>");
    html.push_str("<section class=\"section\"><h2>Method Differences</h2><div class=\"grid plot-grid\">");
    html.push_str(&render_ti_method_difference_plot_card_html(advice));
    html.push_str("</div></section>");
    html.push_str("<section class=\"section\"><details class=\"card disclosure\"><summary class=\"disclosure-summary\">Integration Method Curves</summary><p class=\"disclosure-sub\">Detailed per-method curve cards are available here when you want to inspect the raw interpolants, but they are hidden by default because the difference view is usually more informative.</p><div class=\"disclosure-body\"><div class=\"grid plot-grid\">");
    html.push_str(&render_ti_method_plot_cards_html(advice));
    html.push_str("</div></div></details></section>");
    let mut ranked = advice.intervals().iter().collect::<Vec<_>>();
    ranked.sort_by(|left, right| {
        right
            .priority_score()
            .partial_cmp(&left.priority_score())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    html.push_str("<section class=\"section\"><div class=\"card\"><h2>Priority Queue</h2><table class=\"queue-table\"><thead><tr><th>Interval</th><th>Endpoints</th><th>Slope</th><th>Curvature</th><th>Priority</th><th>Recommendation</th></tr></thead><tbody>");
    for interval in ranked.into_iter().take(5) {
        html.push_str(&format!(
            "<tr><td class=\"mono\"><a class=\"queue-link\" href=\"#interval-{}\">interval {}</a></td><td class=\"mono\">{} -&gt; {}</td><td class=\"mono\">{}</td><td class=\"mono\">{}</td><td class=\"mono\">{:.3}</td><td><span class=\"pill {}\">{}</span></td></tr>",
            interval.interval_index(),
            interval.interval_index(),
            escape_html(&format_report_state(interval.from_state().lambdas())),
            escape_html(&format_report_state(interval.to_state().lambdas())),
            format_report_number(interval.slope()),
            report_option_string(interval.curvature()),
            interval.priority_score(),
            ti_severity_name(interval.severity()),
            ti_severity_name(interval.severity())
        ));
    }
    html.push_str("</tbody></table></div></section>");
    html.push_str("<section class=\"section\"><h2>Suggestions</h2><div class=\"stack\">");
    if advice.suggestions().is_empty() {
        html.push_str("<div class=\"card\">No schedule changes recommended for the current TI thresholds.</div>");
    } else {
        for suggestion in advice.suggestions() {
            let width = (suggestion.priority_score() / max_priority * 100.0).clamp(0.0, 100.0);
            html.push_str("<article class=\"card\">");
            html.push_str(&format!(
                "<div><span class=\"pill {}\">{}</span><a class=\"queue-link mono\" href=\"#interval-{}\"> interval {}: λ {} -&gt; {}</a></div>",
                ti_suggestion_name(suggestion.kind()),
                ti_suggestion_name(suggestion.kind()),
                suggestion.interval_index(),
                suggestion.interval_index(),
                escape_html(&format_report_state(suggestion.from_state().lambdas())),
                escape_html(&format_report_state(suggestion.to_state().lambdas()))
            ));
            html.push_str(&render_ti_suggestion_rationale_badges_html(
                advice,
                suggestion,
                run_options.block_cv_min,
            ));
            html.push_str(&format!("<p>{}</p><div class=\"metric-row\"><div class=\"metric-name\">priority</div><div class=\"bar\"><div class=\"fill bad\" style=\"width:{width:.2}%\"></div></div><div class=\"mono\">{:.3}</div></div>", escape_html(suggestion.reason()), suggestion.priority_score()));
            html.push_str("<div class=\"kv\">");
            if let Some(state) = suggestion.proposed_state() {
                html.push_str(&kv_html(
                "proposed λ",
                    &format_report_state(state.lambdas()),
                ));
            }
            html.push_str("</div>");
            html.push_str(&render_ti_suggestion_evidence_html(
                advice,
                suggestion,
                run_options.block_cv_min,
            ));
            html.push_str("</article>");
        }
    }
    html.push_str("</div></section>");
    html.push_str("<section class=\"section\"><h2>Intervals</h2><div class=\"stack\">");
    for interval in advice.intervals() {
        let width = (interval.priority_score() / max_priority * 100.0).clamp(0.0, 100.0);
        html.push_str(&format!(
            "<article class=\"card\" id=\"interval-{}\">",
            interval.interval_index()
        ));
        html.push_str(&format!(
            "<div><span class=\"pill {}\">{}</span><span class=\"mono\"> interval {}: λ {} -&gt; {}</span></div>",
            ti_severity_name(interval.severity()),
            ti_severity_name(interval.severity()),
            interval.interval_index(),
            escape_html(&format_report_state(interval.from_state().lambdas())),
            escape_html(&format_report_state(interval.to_state().lambdas()))
        ));
        html.push_str(&format!("<div class=\"metric-row\"><div class=\"metric-name\">priority</div><div class=\"bar\"><div class=\"fill {}\" style=\"width:{width:.2}%\"></div></div><div class=\"mono\">{:.3}</div></div>", if matches!(interval.severity(), TiEdgeSeverity::AddSampling | TiEdgeSeverity::AddWindow | TiEdgeSeverity::AddWindowAndSampling) { "bad" } else { "" }, interval.priority_score()));
        html.push_str("<div class=\"kv\">");
        html.push_str(&kv_html(
            "delta λ",
            &format_report_number(interval.delta_lambda()),
        ));
        html.push_str(&kv_html(
            "left mean dH/dλ",
            &format_report_number(interval.left_mean_dhdl()),
        ));
        html.push_str(&kv_html(
            "right mean dH/dλ",
            &format_report_number(interval.right_mean_dhdl()),
        ));
        html.push_str(&kv_html("slope", &format_report_number(interval.slope())));
        html.push_str(&kv_html(
            "curvature",
            &report_option_string(interval.curvature()),
        ));
        html.push_str(&kv_html(
            "interval uncertainty",
            &report_option_string(interval.interval_uncertainty()),
        ));
        html.push_str("</div></article>");
    }
    html.push_str("</div></section><div class=\"footer\">Report generated by <span class=\"mono\">alchemrs advise-schedule</span>. Use the JSON output for exact machine-readable values.</div></div></body></html>");
    Ok(html)
}

fn plot_card_html(title: &str, subtitle: &str, body: &str) -> String {
    format!(
        "<div class=\"card plot-stack\"><h3 class=\"plot-title\">{}</h3><p class=\"plot-sub\">{}</p><div class=\"plot-frame\">{}</div></div>",
        escape_html(title),
        escape_html(subtitle),
        body
    )
}

fn render_u_nk_suggestion_rationale_badges_html(
    advice: &ScheduleAdvice,
    suggestion: &alchemrs::ScheduleSuggestion,
    overlap_min: f64,
    block_cv_min: f64,
) -> String {
    let Some(edge) = advice.edges().get(suggestion.edge_index()) else {
        return String::new();
    };
    let badges = u_nk_suggestion_rationale_badges(edge, suggestion, overlap_min, block_cv_min);
    if badges.is_empty() {
        return String::new();
    }

    let mut html = String::from("<div class=\"badge-row\">");
    for badge in badges {
        html.push_str(&format!(
            "<span class=\"badge rationale\">{}</span>",
            escape_html(&badge)
        ));
    }
    html.push_str("</div>");
    html
}

fn u_nk_suggestion_rationale_badges(
    edge: &alchemrs::AdjacentEdgeDiagnostic,
    suggestion: &alchemrs::ScheduleSuggestion,
    overlap_min: f64,
    block_cv_min: f64,
) -> Vec<String> {
    let mut badges = Vec::new();
    if edge.overlap_min() < overlap_min {
        badges.push("low_overlap".to_string());
    }
    if edge.relative_overlap().is_some_and(|value| value < 0.85) {
        badges.push("localized_overlap_gap".to_string());
    }
    if edge.block_cv().is_some_and(|value| value >= block_cv_min) {
        badges.push("high_block_cv".to_string());
    }
    if edge
        .relative_uncertainty()
        .is_some_and(|value| value >= 1.5)
    {
        badges.push("high_relative_uncertainty".to_string());
    }
    if suggestion.proposal_strategy() == Some(alchemrs::ProposalStrategy::FocusedSplit) {
        badges.push("focused_split".to_string());
    }

    match suggestion.kind() {
        SuggestionKind::InsertWindow => {
            if !badges.iter().any(|badge| badge == "low_overlap") {
                badges.push("low_overlap".to_string());
            }
        }
        SuggestionKind::ExtendSampling => {
            if !badges.iter().any(|badge| badge == "high_block_cv")
                && !badges
                    .iter()
                    .any(|badge| badge == "high_relative_uncertainty")
            {
                badges.push("high_relative_uncertainty".to_string());
            }
        }
        SuggestionKind::NoChange => {}
    }
    badges
}

fn render_ti_suggestion_evidence_html(
    advice: &TiScheduleAdvice,
    suggestion: &TiScheduleSuggestion,
    block_cv_min: f64,
) -> String {
    let Some(interval) = advice.intervals().get(suggestion.interval_index()) else {
        return String::new();
    };
    let left_window = advice.windows().get(suggestion.interval_index());
    let right_window = advice.windows().get(suggestion.interval_index() + 1);

    let mut html = String::from("<div class=\"kv\">");
    match suggestion.kind() {
        TiSuggestionKind::ExtendSampling => {
            html.push_str(&kv_html(
                "interval uncertainty",
                &report_option_string(interval.interval_uncertainty()),
            ));
            html.push_str(&kv_html(
                "left block CV",
                &format_ti_sampling_metric(
                    left_window.and_then(TiWindowDiagnostic::block_cv),
                    Some(block_cv_min),
                ),
            ));
            html.push_str(&kv_html(
                "right block CV",
                &format_ti_sampling_metric(
                    right_window.and_then(TiWindowDiagnostic::block_cv),
                    Some(block_cv_min),
                ),
            ));
            html.push_str(&kv_html(
                "left forward/reverse delta",
                &report_option_string(
                    left_window.and_then(TiWindowDiagnostic::forward_reverse_delta),
                ),
            ));
            html.push_str(&kv_html(
                "right forward/reverse delta",
                &report_option_string(
                    right_window.and_then(TiWindowDiagnostic::forward_reverse_delta),
                ),
            ));
        }
        TiSuggestionKind::InsertWindow => {
            html.push_str(&kv_html("slope", &format_report_number(interval.slope())));
            html.push_str(&kv_html(
                "curvature",
                &report_option_string(interval.curvature()),
            ));
            html.push_str(&kv_html(
                "delta λ",
                &format_report_number(interval.delta_lambda()),
            ));
        }
        TiSuggestionKind::InsertWindowAndExtendSampling => {
            html.push_str(&kv_html("slope", &format_report_number(interval.slope())));
            html.push_str(&kv_html(
                "curvature",
                &report_option_string(interval.curvature()),
            ));
            html.push_str(&kv_html(
                "interval uncertainty",
                &report_option_string(interval.interval_uncertainty()),
            ));
            html.push_str(&kv_html(
                "left block CV",
                &format_ti_sampling_metric(
                    left_window.and_then(TiWindowDiagnostic::block_cv),
                    Some(block_cv_min),
                ),
            ));
            html.push_str(&kv_html(
                "right block CV",
                &format_ti_sampling_metric(
                    right_window.and_then(TiWindowDiagnostic::block_cv),
                    Some(block_cv_min),
                ),
            ));
        }
        TiSuggestionKind::NoChange => {}
    }
    html.push_str("</div>");
    html
}

fn render_ti_suggestion_rationale_badges_html(
    advice: &TiScheduleAdvice,
    suggestion: &TiScheduleSuggestion,
    block_cv_min: f64,
) -> String {
    let badges = ti_suggestion_rationale_badges(advice, suggestion, block_cv_min);
    if badges.is_empty() {
        return String::new();
    }

    let mut html = String::from("<div class=\"badge-row\">");
    for badge in badges {
        html.push_str(&format!(
            "<span class=\"badge rationale\">{}</span>",
            escape_html(&badge)
        ));
    }
    html.push_str("</div>");
    html
}

fn ti_suggestion_rationale_badges(
    advice: &TiScheduleAdvice,
    suggestion: &TiScheduleSuggestion,
    block_cv_min: f64,
) -> Vec<String> {
    let Some(interval) = advice.intervals().get(suggestion.interval_index()) else {
        return Vec::new();
    };
    let left_window = advice.windows().get(suggestion.interval_index());
    let right_window = advice.windows().get(suggestion.interval_index() + 1);
    let defaults = TiScheduleAdvisorOptions::default();

    let slope_values = advice
        .intervals()
        .iter()
        .map(|item| item.abs_slope())
        .collect::<Vec<_>>();
    let slope_mean = mean_f64(&slope_values);
    let slope_std = stddev_f64(&slope_values, slope_mean);

    let curvature_values = advice
        .intervals()
        .iter()
        .filter_map(TiIntervalDiagnostic::curvature)
        .map(f64::abs)
        .collect::<Vec<_>>();
    let curvature_mean = mean_f64(&curvature_values);
    let curvature_std = stddev_f64(&curvature_values, curvature_mean);

    let uncertainty_values = advice
        .intervals()
        .iter()
        .filter_map(TiIntervalDiagnostic::interval_uncertainty)
        .collect::<Vec<_>>();
    let uncertainty_mean = mean_f64(&uncertainty_values);
    let uncertainty_std = stddev_f64(&uncertainty_values, uncertainty_mean);

    let slope_z = simple_z_score(interval.abs_slope(), slope_mean, slope_std);
    let curvature_z = interval
        .curvature()
        .map(f64::abs)
        .and_then(|value| simple_z_score(value, curvature_mean, curvature_std));
    let uncertainty_z = interval
        .interval_uncertainty()
        .and_then(|value| simple_z_score(value, uncertainty_mean, uncertainty_std));

    let mut badges = Vec::new();
    if slope_z.is_some_and(|value| value >= defaults.slope_z_min) {
        badges.push("high_slope".to_string());
    }
    if curvature_z.is_some_and(|value| value >= defaults.curvature_z_min) {
        badges.push("high_curvature".to_string());
    }
    if uncertainty_z.is_some_and(|value| value >= defaults.interval_uncertainty_z_min) {
        badges.push("high_interval_uncertainty".to_string());
    }
    if left_window
        .is_some_and(|window| window.block_cv().is_some_and(|value| value >= block_cv_min))
        || right_window
            .is_some_and(|window| window.block_cv().is_some_and(|value| value >= block_cv_min))
    {
        badges.push("high_block_cv".to_string());
    }
    if left_window.is_some_and(high_forward_reverse_delta)
        || right_window.is_some_and(high_forward_reverse_delta)
    {
        badges.push("high_forward_reverse_delta".to_string());
    }

    match suggestion.kind() {
        TiSuggestionKind::InsertWindow => {
            if badges.is_empty() {
                badges.push("high_curvature".to_string());
            }
        }
        TiSuggestionKind::ExtendSampling => {
            if badges.is_empty() {
                badges.push("high_block_cv".to_string());
            }
        }
        TiSuggestionKind::InsertWindowAndExtendSampling => {
            if !badges.iter().any(|badge| badge == "high_curvature") {
                badges.push("high_curvature".to_string());
            }
            if !badges.iter().any(|badge| badge == "high_block_cv")
                && !badges
                    .iter()
                    .any(|badge| badge == "high_interval_uncertainty")
            {
                badges.push("high_block_cv".to_string());
            }
        }
        TiSuggestionKind::NoChange => {}
    }
    badges
}

fn high_forward_reverse_delta(window: &TiWindowDiagnostic) -> bool {
    match (window.forward_reverse_delta(), window.sem_dhdl()) {
        (Some(delta), Some(sem)) => delta > sem * 2.0,
        _ => false,
    }
}

fn mean_f64(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        Some(values.iter().sum::<f64>() / values.len() as f64)
    }
}

fn stddev_f64(values: &[f64], mean: Option<f64>) -> Option<f64> {
    let mean = mean?;
    if values.len() < 2 {
        return None;
    }
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / values.len() as f64;
    Some(variance.sqrt())
}

fn simple_z_score(value: f64, mean: Option<f64>, stddev: Option<f64>) -> Option<f64> {
    let mean = mean?;
    let stddev = stddev?;
    if stddev <= 1.0e-12 {
        return None;
    }
    Some((value - mean) / stddev)
}

fn format_ti_sampling_metric(value: Option<f64>, threshold: Option<f64>) -> String {
    match (value, threshold) {
        (Some(value), Some(threshold)) => format!(
            "{} (threshold {})",
            format_report_number(value),
            format_report_number(threshold)
        ),
        (Some(value), None) => format_report_number(value),
        (None, _) => "n/a".to_string(),
    }
}

fn render_ti_window_plot_svg(advice: &TiScheduleAdvice) -> String {
    let points = advice
        .windows()
        .iter()
        .filter_map(|window| {
            window
                .state()
                .lambdas()
                .first()
                .copied()
                .map(|lambda| (lambda, window.mean_dhdl()))
        })
        .collect::<Vec<_>>();
    render_ti_zero_shaded_series_plot_svg(
        &points,
        &points,
        "λ",
        "mean dH/dλ",
        "series-line",
        "series-point",
        true,
    )
}

fn render_ti_curvature_plot_svg(advice: &TiScheduleAdvice) -> String {
    let points = advice
        .intervals()
        .iter()
        .filter_map(|interval| {
            let from = interval.from_state().lambdas().first().copied()?;
            let to = interval.to_state().lambdas().first().copied()?;
            interval
                .curvature()
                .map(|curvature| (0.5 * (from + to), curvature))
        })
        .collect::<Vec<_>>();
    render_ti_series_plot_svg(
        &points,
        &points,
        "λ midpoint",
        "|curvature|",
        true,
        "series-line curvature",
        "series-point curvature",
        true,
    )
}

fn render_ti_interval_uncertainty_plot_svg(advice: &TiScheduleAdvice) -> String {
    let points = advice
        .intervals()
        .iter()
        .filter_map(|interval| {
            let from = interval.from_state().lambdas().first().copied()?;
            let to = interval.to_state().lambdas().first().copied()?;
            interval
                .interval_uncertainty()
                .map(|uncertainty| (0.5 * (from + to), uncertainty))
        })
        .collect::<Vec<_>>();
    render_ti_series_plot_svg(
        &points,
        &points,
        "λ midpoint",
        "interval uncertainty",
        true,
        "series-line uncertainty",
        "series-point uncertainty",
        true,
    )
}

fn render_ti_method_plot_cards_html(advice: &TiScheduleAdvice) -> String {
    let Some((lambdas, values)) = ti_window_curve_inputs(advice) else {
        return "<div class=\"card\">Not enough TI windows to render integration-method plots.</div>"
            .to_string();
    };
    let window_points = lambdas
        .iter()
        .copied()
        .zip(values.iter().copied())
        .collect::<Vec<_>>();
    let method_curves = ti_method_curves(&lambdas, &values, 24);

    let mut html = String::new();
    for curve in method_curves {
        let draw_line = !matches!(curve.method, IntegrationMethod::GaussianQuadrature);
        let body = if matches!(curve.method, IntegrationMethod::GaussianQuadrature) {
            render_ti_series_plot_svg(
                &curve.points,
                &window_points,
                "λ",
                "method curve",
                true,
                "series-line",
                "series-point",
                draw_line,
            )
        } else {
            render_ti_zero_shaded_series_plot_svg(
                &curve.points,
                &window_points,
                "λ",
                "method curve",
                "series-line",
                "series-point",
                draw_line,
            )
        };
        html.push_str(&plot_card_html(
            ti_method_title(curve.method),
            ti_method_plot_subtitle(curve.method),
            &body,
        ));
    }

    if html.is_empty() {
        "<div class=\"card\">No applicable TI integration methods were available for the current schedule.</div>".to_string()
    } else {
        html
    }
}

fn render_ti_method_difference_plot_card_html(advice: &TiScheduleAdvice) -> String {
    let Some((lambdas, values)) = ti_window_curve_inputs(advice) else {
        return "<div class=\"card\">Not enough TI windows to compare integration-method curves.</div>"
            .to_string();
    };
    let method_curves = ti_method_curves(&lambdas, &values, 48);
    let Some(reference) = method_curves
        .iter()
        .find(|curve| curve.method == IntegrationMethod::Trapezoidal)
    else {
        return "<div class=\"card\">Trapezoidal reference curve was not available for method comparison.</div>"
            .to_string();
    };

    let mut difference_series = Vec::new();
    let mut legend_entries = Vec::new();
    for curve in method_curves.iter().filter(|curve| {
        !matches!(
            curve.method,
            IntegrationMethod::Trapezoidal | IntegrationMethod::GaussianQuadrature
        )
    }) {
        let Some(differences) = ti_curve_difference(reference, curve) else {
            continue;
        };
        let max_abs_delta = differences
            .iter()
            .map(|(_, value)| value.abs())
            .fold(0.0, f64::max);
        let signed_difference_integral = integrate_curve_trapezoidal(&differences);
        difference_series.push(TiPlotSeries {
            method: curve.method,
            points: differences,
        });
        legend_entries.push((curve.method, max_abs_delta, signed_difference_integral));
    }

    if difference_series.is_empty() {
        return "<div class=\"card\">Only one interpolating TI curve applies on this schedule, so there is no shape difference to plot.</div>".to_string();
    }

    let mut body = render_ti_multi_series_plot_svg(
        &difference_series,
        "λ",
        "delta from trapezoidal",
        true,
    );
    body.push_str(&render_ti_method_difference_legend_html(&legend_entries));

    plot_card_html(
        "Deviation From Trapezoidal",
        "All applicable interpolants pass through the same TI windows, so their raw curves often overlap. This zero-centered view magnifies the between-window shape differences.",
        &body,
    )
}

fn ti_window_curve_inputs(advice: &TiScheduleAdvice) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut pairs = advice
        .windows()
        .iter()
        .filter_map(|window| {
            window
                .state()
                .lambdas()
                .first()
                .copied()
                .map(|lambda| (lambda, window.mean_dhdl()))
        })
        .collect::<Vec<_>>();
    if pairs.len() < 2 {
        return None;
    }
    pairs.sort_by(|left, right| left.0.total_cmp(&right.0));
    let lambdas = pairs.iter().map(|(lambda, _)| *lambda).collect::<Vec<_>>();
    let values = pairs.iter().map(|(_, value)| *value).collect::<Vec<_>>();
    Some((lambdas, values))
}

fn ti_method_title(method: IntegrationMethod) -> &'static str {
    match method {
        IntegrationMethod::Trapezoidal => "Trapezoidal",
        IntegrationMethod::Simpson => "Simpson",
        IntegrationMethod::CubicSpline => "Cubic Spline",
        IntegrationMethod::Pchip => "PCHIP",
        IntegrationMethod::Akima => "Akima",
        IntegrationMethod::GaussianQuadrature => "Gaussian Quadrature",
    }
}

fn ti_method_plot_subtitle(method: IntegrationMethod) -> &'static str {
    match method {
        IntegrationMethod::Trapezoidal => {
            "Piecewise-linear interpolant between adjacent TI window means. Dots mark the original TI windows."
        }
        IntegrationMethod::Simpson => {
            "Composite quadratic interpolant across each pair of uniform λ intervals. Dots mark the original TI windows."
        }
        IntegrationMethod::CubicSpline => {
            "Natural cubic spline passing smoothly through the TI window means. Dots mark the original TI windows."
        }
        IntegrationMethod::Pchip => {
            "Shape-preserving cubic Hermite interpolant built from the TI window means. Dots mark the original TI windows."
        }
        IntegrationMethod::Akima => {
            "Local cubic interpolant with Akima slopes for nonuniform TI schedules. Dots mark the original TI windows."
        }
        IntegrationMethod::GaussianQuadrature => {
            "Supported Gauss-Legendre node set. This method does not define a unique interpolating curve, so the report shows the sampled nodes only."
        }
    }
}

#[derive(Clone)]
struct TiPlotSeries {
    method: IntegrationMethod,
    points: Vec<(f64, f64)>,
}

fn ti_method_curves(lambdas: &[f64], values: &[f64], samples_per_interval: usize) -> Vec<TiPlotSeries> {
    let mut curves = Vec::new();
    for method in [
        IntegrationMethod::Trapezoidal,
        IntegrationMethod::Simpson,
        IntegrationMethod::CubicSpline,
        IntegrationMethod::Pchip,
        IntegrationMethod::Akima,
        IntegrationMethod::GaussianQuadrature,
    ] {
        let Ok(points) = sample_ti_curve(lambdas, values, method, samples_per_interval) else {
            continue;
        };
        curves.push(TiPlotSeries { method, points });
    }
    curves
}

fn ti_curve_difference(reference: &TiPlotSeries, other: &TiPlotSeries) -> Option<Vec<(f64, f64)>> {
    let mut differences = Vec::with_capacity(reference.points.len());
    for (x, y_ref) in &reference.points {
        let y_other = interpolate_sampled_curve(&other.points, *x)?;
        differences.push((*x, y_other - y_ref));
    }
    Some(differences)
}

fn interpolate_sampled_curve(points: &[(f64, f64)], x: f64) -> Option<f64> {
    if points.is_empty() {
        return None;
    }
    if x < points.first()?.0 - 1.0e-9 || x > points.last()?.0 + 1.0e-9 {
        return None;
    }

    match points.binary_search_by(|(candidate_x, _)| candidate_x.total_cmp(&x)) {
        Ok(index) => Some(points[index].1),
        Err(0) => Some(points[0].1),
        Err(index) if index >= points.len() => Some(points.last()?.1),
        Err(index) => {
            let (x0, y0) = points[index - 1];
            let (x1, y1) = points[index];
            let dx = x1 - x0;
            if dx.abs() <= 1.0e-12 {
                Some(y0)
            } else {
                let t = (x - x0) / dx;
                Some(y0 + t * (y1 - y0))
            }
        }
    }
}

fn integrate_curve_trapezoidal(points: &[(f64, f64)]) -> f64 {
    if points.len() < 2 {
        return 0.0;
    }

    let mut integral = 0.0;
    for window in points.windows(2) {
        let (x0, y0) = window[0];
        let (x1, y1) = window[1];
        integral += (x1 - x0) * (y0 + y1) * 0.5;
    }
    integral
}

fn render_ti_method_difference_legend_html(entries: &[(IntegrationMethod, f64, f64)]) -> String {
    let mut html = String::from("<div class=\"legend\">");
    for (method, max_abs_delta, signed_difference_integral) in entries {
        html.push_str(&format!(
            "<div class=\"legend-item\"><span class=\"legend-swatch {}\"></span><span>{}: max |Δ| {}, ∫Δ dλ {}</span></div>",
            ti_method_css_class(*method),
            escape_html(ti_method_title(*method)),
            format_plot_number(*max_abs_delta),
            format_plot_number(*signed_difference_integral),
        ));
    }
    html.push_str("</div>");
    html
}

fn ti_method_css_class(method: IntegrationMethod) -> &'static str {
    match method {
        IntegrationMethod::Trapezoidal => "method-trapezoidal",
        IntegrationMethod::Simpson => "method-simpson",
        IntegrationMethod::CubicSpline => "method-cubic-spline",
        IntegrationMethod::Pchip => "method-pchip",
        IntegrationMethod::Akima => "method-akima",
        IntegrationMethod::GaussianQuadrature => "method-gaussian-quadrature",
    }
}

fn axis_tick_values(min: f64, max: f64, include_zero: bool) -> Vec<f64> {
    let mut ticks = (0..=4)
        .map(|tick| {
            let frac = tick as f64 / 4.0;
            max - frac * (max - min)
        })
        .collect::<Vec<_>>();

    if include_zero && min <= 0.0 && max >= 0.0 {
        if let Some((closest_index, _)) = ticks
            .iter()
            .enumerate()
            .map(|(index, value)| (index, value.abs()))
            .min_by(|left, right| left.1.total_cmp(&right.1))
        {
            ticks[closest_index] = 0.0;
            ticks.sort_by(|left, right| right.total_cmp(left));
        }
    }

    ticks
}

fn render_ti_series_plot_svg(
    line_points: &[(f64, f64)],
    marker_points: &[(f64, f64)],
    x_label: &str,
    y_label: &str,
    include_zero: bool,
    line_class: &str,
    point_class: &str,
    draw_line: bool,
) -> String {
    if (line_points.is_empty() && marker_points.is_empty()) || (draw_line && line_points.len() < 2)
    {
        return "<div class=\"plot-empty\">Not enough points to render plot.</div>".to_string();
    }

    let bounds_points = if marker_points.is_empty() {
        line_points
    } else {
        marker_points
    };

    let width = 520.0;
    let height = 240.0;
    let left = 58.0;
    let right = 18.0;
    let top = 16.0;
    let bottom = 40.0;
    let plot_width = width - left - right;
    let plot_height = height - top - bottom;

    let mut x_min = bounds_points
        .iter()
        .map(|(x, _)| *x)
        .fold(f64::INFINITY, f64::min);
    let mut x_max = bounds_points
        .iter()
        .map(|(x, _)| *x)
        .fold(f64::NEG_INFINITY, f64::max);
    if (x_max - x_min).abs() <= 1.0e-12 {
        x_min -= 0.5;
        x_max += 0.5;
    }

    let mut y_min = bounds_points
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::INFINITY, f64::min);
    let mut y_max = bounds_points
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::NEG_INFINITY, f64::max);
    if include_zero {
        y_min = y_min.min(0.0);
        y_max = y_max.max(0.0);
    }
    if (y_max - y_min).abs() <= 1.0e-12 {
        let pad = y_max.abs().max(1.0) * 0.25;
        y_min -= pad;
        y_max += pad;
    } else {
        let pad = (y_max - y_min) * 0.12;
        y_min -= pad;
        y_max += pad;
    }

    let map_x = |x: f64| left + (x - x_min) / (x_max - x_min) * plot_width;
    let map_y = |y: f64| top + (1.0 - (y - y_min) / (y_max - y_min)) * plot_height;

    let polyline = if draw_line {
        Some(
            line_points
                .iter()
                .map(|(x, y)| format!("{:.2},{:.2}", map_x(*x), map_y(*y)))
                .collect::<Vec<_>>()
                .join(" "),
        )
    } else {
        None
    };

    let mut svg = format!(
        "<svg class=\"ti-series-plot\" viewBox=\"0 0 {:.0} {:.0}\" xmlns=\"http://www.w3.org/2000/svg\" role=\"img\" aria-label=\"{} versus {}\">",
        width,
        height,
        escape_html(y_label),
        escape_html(x_label)
    );

    for tick in 0..=4 {
        let frac = tick as f64 / 4.0;
        let x = left + frac * plot_width;
        let x_value = x_min + frac * (x_max - x_min);
        svg.push_str(&format!(
            "<line class=\"grid-line\" x1=\"{x:.2}\" y1=\"{top:.2}\" x2=\"{x:.2}\" y2=\"{:.2}\" />",
            top + plot_height
        ));
        svg.push_str(&format!(
            "<text class=\"tick-label\" x=\"{x:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
            top + plot_height + 18.0,
            format_plot_number(x_value)
        ));
    }

    for y_value in axis_tick_values(y_min, y_max, include_zero) {
        let y = map_y(y_value);
        svg.push_str(&format!(
            "<line class=\"grid-line\" x1=\"{left:.2}\" y1=\"{y:.2}\" x2=\"{:.2}\" y2=\"{y:.2}\" />",
            left + plot_width
        ));
        svg.push_str(&format!(
            "<text class=\"tick-label\" x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\">{}</text>",
            left - 8.0,
            y + 3.0,
            format_plot_number(y_value)
        ));
    }

    if include_zero && y_min <= 0.0 && y_max >= 0.0 {
        let zero_y = map_y(0.0);
        svg.push_str(&format!(
            "<line class=\"zero-line\" x1=\"{left:.2}\" y1=\"{zero_y:.2}\" x2=\"{:.2}\" y2=\"{zero_y:.2}\" />",
            left + plot_width
        ));
    }

    svg.push_str(&format!(
        "<line class=\"axis\" x1=\"{left:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" />",
        top + plot_height,
        left + plot_width,
        top + plot_height
    ));
    svg.push_str(&format!(
        "<line class=\"axis\" x1=\"{left:.2}\" y1=\"{top:.2}\" x2=\"{left:.2}\" y2=\"{:.2}\" />",
        top + plot_height
    ));
    svg.push_str(&format!(
        "<text class=\"axis-label\" x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
        left + plot_width * 0.5,
        height - 8.0,
        escape_html(x_label)
    ));
    svg.push_str(&format!(
        "<text class=\"axis-label\" x=\"16\" y=\"{:.2}\" transform=\"rotate(-90 16 {:.2})\" text-anchor=\"middle\">{}</text>",
        top + plot_height * 0.5,
        top + plot_height * 0.5,
        escape_html(y_label)
    ));
    if let Some(polyline) = polyline.as_ref() {
        svg.push_str(&format!(
            "<polyline class=\"{}\" points=\"{}\" />",
            line_class, polyline
        ));
    }
    for (x, y) in marker_points {
        svg.push_str(&format!(
            "<circle class=\"{}\" cx=\"{:.2}\" cy=\"{:.2}\" r=\"3.5\" />",
            point_class,
            map_x(*x),
            map_y(*y)
        ));
    }
    svg.push_str("</svg>");
    svg
}

fn render_ti_zero_shaded_series_plot_svg(
    line_points: &[(f64, f64)],
    marker_points: &[(f64, f64)],
    x_label: &str,
    y_label: &str,
    line_class: &str,
    point_class: &str,
    draw_line: bool,
) -> String {
    if (line_points.is_empty() && marker_points.is_empty()) || (draw_line && line_points.len() < 2)
    {
        return "<div class=\"plot-empty\">Not enough points to render plot.</div>".to_string();
    }

    let bounds_points = if marker_points.is_empty() {
        line_points
    } else {
        marker_points
    };

    let width = 520.0;
    let height = 240.0;
    let left = 58.0;
    let right = 18.0;
    let top = 16.0;
    let bottom = 40.0;
    let plot_width = width - left - right;
    let plot_height = height - top - bottom;

    let mut x_min = bounds_points
        .iter()
        .map(|(x, _)| *x)
        .fold(f64::INFINITY, f64::min);
    let mut x_max = bounds_points
        .iter()
        .map(|(x, _)| *x)
        .fold(f64::NEG_INFINITY, f64::max);
    if (x_max - x_min).abs() <= 1.0e-12 {
        x_min -= 0.5;
        x_max += 0.5;
    }

    let mut y_min = bounds_points
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::INFINITY, f64::min);
    let mut y_max = bounds_points
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::NEG_INFINITY, f64::max);
    y_min = y_min.min(0.0);
    y_max = y_max.max(0.0);
    if (y_max - y_min).abs() <= 1.0e-12 {
        let pad = y_max.abs().max(1.0) * 0.25;
        y_min -= pad;
        y_max += pad;
    } else {
        let pad = (y_max - y_min) * 0.12;
        y_min -= pad;
        y_max += pad;
    }

    let map_x = |x: f64| left + (x - x_min) / (x_max - x_min) * plot_width;
    let map_y = |y: f64| top + (1.0 - (y - y_min) / (y_max - y_min)) * plot_height;
    let zero_y = map_y(0.0);

    let polyline = if draw_line {
        Some(
            line_points
                .iter()
                .map(|(x, y)| format!("{:.2},{:.2}", map_x(*x), map_y(*y)))
                .collect::<Vec<_>>()
                .join(" "),
        )
    } else {
        None
    };

    let mut svg = format!(
        "<svg class=\"ti-series-plot\" viewBox=\"0 0 {:.0} {:.0}\" xmlns=\"http://www.w3.org/2000/svg\" role=\"img\" aria-label=\"{} versus {}\">",
        width,
        height,
        escape_html(y_label),
        escape_html(x_label)
    );

    for tick in 0..=4 {
        let frac = tick as f64 / 4.0;
        let x = left + frac * plot_width;
        let x_value = x_min + frac * (x_max - x_min);
        svg.push_str(&format!(
            "<line class=\"grid-line\" x1=\"{x:.2}\" y1=\"{top:.2}\" x2=\"{x:.2}\" y2=\"{:.2}\" />",
            top + plot_height
        ));
        svg.push_str(&format!(
            "<text class=\"tick-label\" x=\"{x:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
            top + plot_height + 18.0,
            format_plot_number(x_value)
        ));
    }

    for y_value in axis_tick_values(y_min, y_max, true) {
        let y = map_y(y_value);
        svg.push_str(&format!(
            "<line class=\"grid-line\" x1=\"{left:.2}\" y1=\"{y:.2}\" x2=\"{:.2}\" y2=\"{y:.2}\" />",
            left + plot_width
        ));
        svg.push_str(&format!(
            "<text class=\"tick-label\" x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\">{}</text>",
            left - 8.0,
            y + 3.0,
            format_plot_number(y_value)
        ));
    }

    for polygon in zero_area_polygons(line_points, true) {
        svg.push_str(&format!(
            "<polygon class=\"series-fill positive\" points=\"{}\" />",
            polygon
                .iter()
                .map(|(x, y)| format!("{:.2},{:.2}", map_x(*x), map_y(*y)))
                .collect::<Vec<_>>()
                .join(" ")
        ));
    }
    for polygon in zero_area_polygons(line_points, false) {
        svg.push_str(&format!(
            "<polygon class=\"series-fill negative\" points=\"{}\" />",
            polygon
                .iter()
                .map(|(x, y)| format!("{:.2},{:.2}", map_x(*x), map_y(*y)))
                .collect::<Vec<_>>()
                .join(" ")
        ));
    }

    svg.push_str(&format!(
        "<line class=\"zero-line\" x1=\"{left:.2}\" y1=\"{zero_y:.2}\" x2=\"{:.2}\" y2=\"{zero_y:.2}\" />",
        left + plot_width
    ));
    svg.push_str(&format!(
        "<line class=\"axis\" x1=\"{left:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" />",
        top + plot_height,
        left + plot_width,
        top + plot_height
    ));
    svg.push_str(&format!(
        "<line class=\"axis\" x1=\"{left:.2}\" y1=\"{top:.2}\" x2=\"{left:.2}\" y2=\"{:.2}\" />",
        top + plot_height
    ));
    svg.push_str(&format!(
        "<text class=\"axis-label\" x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
        left + plot_width * 0.5,
        height - 8.0,
        escape_html(x_label)
    ));
    svg.push_str(&format!(
        "<text class=\"axis-label\" x=\"16\" y=\"{:.2}\" transform=\"rotate(-90 16 {:.2})\" text-anchor=\"middle\">{}</text>",
        top + plot_height * 0.5,
        top + plot_height * 0.5,
        escape_html(y_label)
    ));
    if let Some(polyline) = polyline.as_ref() {
        svg.push_str(&format!(
            "<polyline class=\"{}\" points=\"{}\" />",
            line_class, polyline
        ));
    }
    for (x, y) in marker_points {
        svg.push_str(&format!(
            "<circle class=\"{}\" cx=\"{:.2}\" cy=\"{:.2}\" r=\"3.5\" />",
            point_class,
            map_x(*x),
            map_y(*y)
        ));
    }
    svg.push_str("</svg>");
    svg
}

fn zero_area_polygons(points: &[(f64, f64)], positive: bool) -> Vec<Vec<(f64, f64)>> {
    let mut polygons = Vec::new();
    for window in points.windows(2) {
        let (x0, y0) = window[0];
        let (x1, y1) = window[1];
        if (x1 - x0).abs() <= 1.0e-12 {
            continue;
        }

        let crosses_zero = (y0 > 0.0 && y1 < 0.0) || (y0 < 0.0 && y1 > 0.0);
        let x_zero = if crosses_zero {
            Some(x0 + (x1 - x0) * (-y0) / (y1 - y0))
        } else {
            None
        };

        if positive {
            if y0 > 0.0 && y1 > 0.0 {
                polygons.push(vec![(x0, 0.0), (x0, y0), (x1, y1), (x1, 0.0)]);
            } else if y0 == 0.0 && y1 > 0.0 {
                polygons.push(vec![(x0, 0.0), (x1, y1), (x1, 0.0)]);
            } else if y0 > 0.0 && y1 == 0.0 {
                polygons.push(vec![(x0, 0.0), (x0, y0), (x1, 0.0)]);
            } else if y0 > 0.0 && y1 < 0.0 {
                polygons.push(vec![(x0, 0.0), (x0, y0), (x_zero.unwrap(), 0.0)]);
            } else if y0 < 0.0 && y1 > 0.0 {
                polygons.push(vec![(x_zero.unwrap(), 0.0), (x1, y1), (x1, 0.0)]);
            }
        } else if y0 < 0.0 && y1 < 0.0 {
            polygons.push(vec![(x0, 0.0), (x0, y0), (x1, y1), (x1, 0.0)]);
        } else if y0 == 0.0 && y1 < 0.0 {
            polygons.push(vec![(x0, 0.0), (x1, y1), (x1, 0.0)]);
        } else if y0 < 0.0 && y1 == 0.0 {
            polygons.push(vec![(x0, 0.0), (x0, y0), (x1, 0.0)]);
        } else if y0 < 0.0 && y1 > 0.0 {
            polygons.push(vec![(x0, 0.0), (x0, y0), (x_zero.unwrap(), 0.0)]);
        } else if y0 > 0.0 && y1 < 0.0 {
            polygons.push(vec![(x_zero.unwrap(), 0.0), (x1, y1), (x1, 0.0)]);
        }
    }
    polygons
}

fn render_ti_multi_series_plot_svg(
    series: &[TiPlotSeries],
    x_label: &str,
    y_label: &str,
    include_zero: bool,
) -> String {
    let points = series
        .iter()
        .flat_map(|item| item.points.iter().copied())
        .collect::<Vec<_>>();
    if points.is_empty() {
        return "<div class=\"plot-empty\">Not enough points to render plot.</div>".to_string();
    }

    let width = 520.0;
    let height = 240.0;
    let left = 58.0;
    let right = 18.0;
    let top = 16.0;
    let bottom = 40.0;
    let plot_width = width - left - right;
    let plot_height = height - top - bottom;

    let mut x_min = points.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
    let mut x_max = points
        .iter()
        .map(|(x, _)| *x)
        .fold(f64::NEG_INFINITY, f64::max);
    if (x_max - x_min).abs() <= 1.0e-12 {
        x_min -= 0.5;
        x_max += 0.5;
    }

    let mut y_min = points.iter().map(|(_, y)| *y).fold(f64::INFINITY, f64::min);
    let mut y_max = points
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::NEG_INFINITY, f64::max);
    if include_zero {
        y_min = y_min.min(0.0);
        y_max = y_max.max(0.0);
    }
    if (y_max - y_min).abs() <= 1.0e-12 {
        let pad = y_max.abs().max(1.0) * 0.25;
        y_min -= pad;
        y_max += pad;
    } else {
        let pad = (y_max - y_min) * 0.12;
        y_min -= pad;
        y_max += pad;
    }

    let map_x = |x: f64| left + (x - x_min) / (x_max - x_min) * plot_width;
    let map_y = |y: f64| top + (1.0 - (y - y_min) / (y_max - y_min)) * plot_height;

    let mut svg = format!(
        "<svg class=\"ti-series-plot\" viewBox=\"0 0 {:.0} {:.0}\" xmlns=\"http://www.w3.org/2000/svg\" role=\"img\" aria-label=\"{} versus {}\">",
        width,
        height,
        escape_html(y_label),
        escape_html(x_label)
    );

    for tick in 0..=4 {
        let frac = tick as f64 / 4.0;
        let x = left + frac * plot_width;
        let y = top + frac * plot_height;
        let x_value = x_min + frac * (x_max - x_min);
        let y_value = y_max - frac * (y_max - y_min);
        svg.push_str(&format!(
            "<line class=\"grid-line\" x1=\"{x:.2}\" y1=\"{top:.2}\" x2=\"{x:.2}\" y2=\"{:.2}\" />",
            top + plot_height
        ));
        svg.push_str(&format!(
            "<line class=\"grid-line\" x1=\"{left:.2}\" y1=\"{y:.2}\" x2=\"{:.2}\" y2=\"{y:.2}\" />",
            left + plot_width
        ));
        svg.push_str(&format!(
            "<text class=\"tick-label\" x=\"{x:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
            top + plot_height + 18.0,
            format_plot_number(x_value)
        ));
        svg.push_str(&format!(
            "<text class=\"tick-label\" x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\">{}</text>",
            left - 8.0,
            y + 3.0,
            format_plot_number(y_value)
        ));
    }

    if include_zero && y_min <= 0.0 && y_max >= 0.0 {
        let zero_y = map_y(0.0);
        svg.push_str(&format!(
            "<line class=\"zero-line\" x1=\"{left:.2}\" y1=\"{zero_y:.2}\" x2=\"{:.2}\" y2=\"{zero_y:.2}\" />",
            left + plot_width
        ));
    }

    svg.push_str(&format!(
        "<line class=\"axis\" x1=\"{left:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" />",
        top + plot_height,
        left + plot_width,
        top + plot_height
    ));
    svg.push_str(&format!(
        "<line class=\"axis\" x1=\"{left:.2}\" y1=\"{top:.2}\" x2=\"{left:.2}\" y2=\"{:.2}\" />",
        top + plot_height
    ));
    svg.push_str(&format!(
        "<text class=\"axis-label\" x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
        left + plot_width * 0.5,
        height - 8.0,
        escape_html(x_label)
    ));
    svg.push_str(&format!(
        "<text class=\"axis-label\" x=\"16\" y=\"{:.2}\" transform=\"rotate(-90 16 {:.2})\" text-anchor=\"middle\">{}</text>",
        top + plot_height * 0.5,
        top + plot_height * 0.5,
        escape_html(y_label)
    ));

    for item in series {
        let polyline = item
            .points
            .iter()
            .map(|(x, y)| format!("{:.2},{:.2}", map_x(*x), map_y(*y)))
            .collect::<Vec<_>>()
            .join(" ");
        svg.push_str(&format!(
            "<polyline class=\"series-line {}\" points=\"{}\" />",
            ti_method_css_class(item.method),
            polyline
        ));
    }

    svg.push_str("</svg>");
    svg
}

fn format_plot_number(value: f64) -> String {
    let formatted = if value.abs() >= 1000.0 {
        format!("{value:.0}")
    } else if value.abs() >= 100.0 {
        format!("{value:.1}")
    } else if value.abs() >= 10.0 {
        format!("{value:.2}")
    } else {
        format!("{value:.3}")
    };
    trim_trailing_zeros(&formatted)
}

fn estimator_name(estimator: AdvisorEstimatorArg) -> &'static str {
    match estimator {
        AdvisorEstimatorArg::Mbar => "mbar",
        AdvisorEstimatorArg::Bar => "bar",
    }
}

fn advise_input_kind_name(kind: AdviseInputKind) -> &'static str {
    match kind {
        AdviseInputKind::Auto => "auto",
        AdviseInputKind::UNk => "u_nk",
        AdviseInputKind::Dhdl => "dhdl",
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

fn ti_severity_name(severity: TiEdgeSeverity) -> &'static str {
    match severity {
        TiEdgeSeverity::Healthy => "healthy",
        TiEdgeSeverity::Monitor => "monitor",
        TiEdgeSeverity::AddSampling => "add_sampling",
        TiEdgeSeverity::AddWindow => "add_window",
        TiEdgeSeverity::AddWindowAndSampling => "add_window_and_sampling",
    }
}

fn suggestion_name(kind: SuggestionKind) -> &'static str {
    match kind {
        SuggestionKind::NoChange => "no_change",
        SuggestionKind::ExtendSampling => "extend_sampling",
        SuggestionKind::InsertWindow => "insert_window",
    }
}

fn ti_suggestion_name(kind: TiSuggestionKind) -> &'static str {
    match kind {
        TiSuggestionKind::NoChange => "no_change",
        TiSuggestionKind::ExtendSampling => "extend_sampling",
        TiSuggestionKind::InsertWindow => "insert_window",
        TiSuggestionKind::InsertWindowAndExtendSampling => "insert_window_and_extend_sampling",
    }
}

fn proposal_strategy_name(strategy: alchemrs::ProposalStrategy) -> &'static str {
    match strategy {
        alchemrs::ProposalStrategy::Midpoint => "midpoint",
        alchemrs::ProposalStrategy::FocusedSplit => "focused_split",
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

fn format_report_number(value: f64) -> String {
    trim_trailing_zeros(&format!("{value:.3}"))
}

fn report_option_string(value: Option<f64>) -> String {
    value.map(format_report_number).unwrap_or_default()
}

fn format_lambda_value(value: f64) -> String {
    let formatted = trim_trailing_zeros(&format!("{value:.3}"));
    if !formatted.contains('.') {
        format!("{formatted}.0")
    } else {
        formatted
    }
}

fn format_report_state(state: &[f64]) -> String {
    if state.len() == 1 {
        format_lambda_value(state[0])
    } else {
        format!(
            "[{}]",
            state
                .iter()
                .map(|value| format_lambda_value(*value))
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

fn trim_trailing_zeros(value: &str) -> String {
    let trimmed = value.trim_end_matches('0').trim_end_matches('.');
    match trimmed {
        "-0" | "" => "0".to_string(),
        other => other.to_string(),
    }
}

fn option_string(value: Option<f64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn format_string_list(values: &[String]) -> String {
    if values.is_empty() {
        String::new()
    } else {
        format!("[{}]", values.join(";"))
    }
}

fn normalize_lambda_components(lambda_components: Option<Vec<String>>) -> Option<Vec<String>> {
    let labels = lambda_components?
        .into_iter()
        .map(|label| label.trim().to_string())
        .filter(|label| !label.is_empty())
        .collect::<Vec<_>>();
    if labels.is_empty() {
        None
    } else {
        Some(labels)
    }
}

fn summary_card(label: &str, value: &str, sub: &str) -> String {
    format!(
        "<div class=\"card\"><div class=\"label\">{}</div><div class=\"value\">{}</div><div class=\"sub\">{}</div></div>",
        escape_html(label),
        escape_html(value),
        escape_html(sub)
    )
}

fn kv_html(label: &str, value: &str) -> String {
    format!(
        "<div><div class=\"label\">{}</div><div class=\"mono\">{}</div></div>",
        escape_html(label),
        if value.is_empty() {
            "&nbsp;".to_string()
        } else {
            escape_html(value)
        }
    )
}

fn escape_html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn severity_class_from_suggestion(kind: SuggestionKind) -> &'static str {
    match kind {
        SuggestionKind::NoChange => "healthy",
        SuggestionKind::ExtendSampling => "add_sampling",
        SuggestionKind::InsertWindow => "add_window",
    }
}

#[cfg(test)]
mod tests {
    use alchemrs::{
        advise_lambda_schedule, advise_ti_schedule, extract_u_nk, AdvisorEstimator, DhdlSeries,
        ScheduleAdvisorOptions, StatePoint, TiScheduleAdvisorOptions, UNkMatrix,
    };
    use serde_json::Value;

    use super::{
        axis_tick_values, render_advice, render_html_report, render_ti_html_report,
        render_ti_method_plot_cards_html, AdviseRunOptions,
        zero_area_polygons,
    };
    use crate::cli::input::{AnalysisInputOptions, AnalysisSampleCounts};
    use crate::cli::{AdviseInputKind, AdvisorEstimatorArg, OutputFormat, UNkObservable};

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

    fn sample_multidimensional_advice() -> alchemrs::ScheduleAdvice {
        let states = vec![
            StatePoint::new(vec![0.0, 0.0], 298.0).unwrap(),
            StatePoint::new(vec![0.8, 0.1], 298.0).unwrap(),
            StatePoint::new(vec![1.0, 1.0], 298.0).unwrap(),
        ];
        let labels = vec!["coul-lambda".to_string(), "vdw-lambda".to_string()];
        let windows = vec![
            build_window_with_labels(
                vec![0.0, 0.0],
                &[
                    [0.0, 2.0, 3.0],
                    [0.0, 2.1, 3.1],
                    [0.0, 1.9, 2.9],
                    [0.0, 2.0, 3.0],
                ],
                &states,
                labels.clone(),
            ),
            build_window_with_labels(
                vec![0.8, 0.1],
                &[
                    [2.0, 0.0, 0.1],
                    [2.1, 0.0, 0.1],
                    [1.9, 0.0, 0.2],
                    [2.0, 0.0, 0.1],
                ],
                &states,
                labels.clone(),
            ),
            build_window_with_labels(
                vec![1.0, 1.0],
                &[
                    [3.0, 0.1, 0.0],
                    [3.1, 0.1, 0.0],
                    [2.9, 0.2, 0.0],
                    [3.0, 0.1, 0.0],
                ],
                &states,
                labels.clone(),
            ),
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

    fn build_window_with_labels(
        sampled_lambdas: Vec<f64>,
        rows: &[[f64; 3]],
        evaluated_states: &[StatePoint],
        labels: Vec<String>,
    ) -> UNkMatrix {
        let n_samples = rows.len();
        let data = rows
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let time = (0..n_samples).map(|idx| idx as f64).collect::<Vec<_>>();
        UNkMatrix::new_with_labels(
            n_samples,
            3,
            data,
            time,
            Some(StatePoint::new(sampled_lambdas, 298.0).unwrap()),
            evaluated_states.to_vec(),
            Some(labels),
        )
        .unwrap()
    }

    fn sample_ti_advice() -> alchemrs::TiScheduleAdvice {
        let states = [
            StatePoint::new(vec![0.0], 298.0).unwrap(),
            StatePoint::new(vec![0.33], 298.0).unwrap(),
            StatePoint::new(vec![0.66], 298.0).unwrap(),
            StatePoint::new(vec![1.0], 298.0).unwrap(),
        ];
        let series = vec![
            DhdlSeries::new(
                states[0].clone(),
                vec![0.0, 1.0, 2.0, 3.0],
                vec![0.0, 0.0, 0.0, 0.0],
            )
            .unwrap(),
            DhdlSeries::new(
                states[1].clone(),
                vec![0.0, 1.0, 2.0, 3.0],
                vec![0.5, 0.5, 0.5, 0.5],
            )
            .unwrap(),
            DhdlSeries::new(
                states[2].clone(),
                vec![0.0, 1.0, 2.0, 3.0],
                vec![6.0, 6.0, 6.0, 6.0],
            )
            .unwrap(),
            DhdlSeries::new(
                states[3].clone(),
                vec![0.0, 1.0, 2.0, 3.0],
                vec![1.0, 1.0, 1.0, 1.0],
            )
            .unwrap(),
        ];
        advise_ti_schedule(
            &series,
            Some(TiScheduleAdvisorOptions {
                n_blocks: 2,
                block_cv_min: 1.0e30,
                curvature_z_min: 0.1,
                slope_z_min: 0.1,
                interval_uncertainty_z_min: 10.0,
                suggest_midpoints: true,
            }),
        )
        .unwrap()
    }

    fn sample_ti_sampling_advice() -> alchemrs::TiScheduleAdvice {
        let states = [
            StatePoint::new(vec![0.0], 298.0).unwrap(),
            StatePoint::new(vec![0.5], 298.0).unwrap(),
            StatePoint::new(vec![1.0], 298.0).unwrap(),
        ];
        let series = vec![
            DhdlSeries::new(
                states[0].clone(),
                vec![0.0, 1.0, 2.0, 3.0],
                vec![0.0, 0.0, 0.0, 0.0],
            )
            .unwrap(),
            DhdlSeries::new(
                states[1].clone(),
                vec![0.0, 1.0, 2.0, 3.0],
                vec![1.0, 1.0, 9.0, 9.0],
            )
            .unwrap(),
            DhdlSeries::new(
                states[2].clone(),
                vec![0.0, 1.0, 2.0, 3.0],
                vec![10.0, 10.0, 10.0, 10.0],
            )
            .unwrap(),
        ];
        advise_ti_schedule(
            &series,
            Some(TiScheduleAdvisorOptions {
                n_blocks: 2,
                block_cv_min: 0.15,
                curvature_z_min: 10.0,
                slope_z_min: 10.0,
                interval_uncertainty_z_min: 10.0,
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
                input_kind: AdviseInputKind::UNk,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 4,
                suggest_midpoints: true,
                output_format: OutputFormat::Json,
                output_path: None,
                report_path: None,
            },
            None,
        )
        .unwrap();

        let value: Value = serde_json::from_str(&output).unwrap();
        assert_eq!(value["edges"][0]["severity"], "add_window");
        assert_eq!(value["suggestions"][0]["kind"], "insert_window");
        assert_eq!(value["suggestions"][0]["proposal_strategy"], "midpoint");
        assert_eq!(
            value["suggestions"][0]["proposed_lambda"],
            serde_json::json!([0.0, 0.0, 0.025, 0.0, 0.0])
        );
    }

    #[test]
    fn renders_schedule_advice_as_html_report() {
        let output = render_html_report(
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
                input_kind: AdviseInputKind::UNk,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 4,
                suggest_midpoints: true,
                output_format: OutputFormat::Json,
                output_path: None,
                report_path: None,
            },
            None,
        )
        .unwrap();

        assert!(output.contains("Lambda Schedule Report"));
        assert!(output.contains("Priority Queue"));
        assert!(output.contains("Weakest Component"));
        assert!(output.contains("href=\"#edge-"));
        assert!(output.contains("id=\"edge-"));
        assert!(output.contains("<a class=\"queue-link mono\" href=\"#edge-"));
        assert!(output.contains("&lambda;"));
        assert!(output.contains("badge rationale"));
        assert!(output.contains("focused_split") || output.contains("midpoint"));
        assert!(output.contains("priority"));
        assert!(!output.contains("<svg"));
        assert!(!output.contains("Legend"));
    }

    #[test]
    fn omits_blank_lambda_components_from_html_report() {
        let output = render_html_report(
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
                input_kind: AdviseInputKind::UNk,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 4,
                suggest_midpoints: true,
                output_format: OutputFormat::Json,
                output_path: None,
                report_path: None,
            },
            Some(vec!["".to_string(), "   ".to_string()]),
        )
        .unwrap();

        assert!(!output.contains("lambda components"));
    }

    #[test]
    fn renders_multidimensional_component_breakdown_in_html_report() {
        let output = render_html_report(
            &sample_multidimensional_advice(),
            AnalysisSampleCounts {
                windows: 3,
                samples_in: 12,
                samples_after_burnin: 12,
                samples_kept: 12,
            },
            &AnalysisInputOptions {
                temperature: 300.0,
                decorrelate: false,
                remove_burnin: 0,
                auto_equilibrate: false,
                fast: false,
                conservative: true,
                nskip: 1,
                u_nk_observable: Some(UNkObservable::De),
            },
            &AdviseRunOptions {
                estimator: AdvisorEstimatorArg::Mbar,
                input_kind: AdviseInputKind::UNk,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 4,
                suggest_midpoints: true,
                output_format: OutputFormat::Json,
                output_path: None,
                report_path: None,
            },
            Some(vec!["coul-lambda".to_string(), "vdw-lambda".to_string()]),
        )
        .unwrap();

        assert!(output.contains("Priority Queue"));
        assert!(output.contains("href=\"#edge-0\"") || output.contains("href=\"#edge-1\""));
        assert!(output.contains("coul-lambda"));
        assert!(output.contains("vdw-lambda"));
        assert!(output.contains("lambda components"));
        assert!(output.contains("proposed lambda"));
        assert!(!output.contains("<div class=\"viz\">"));
        assert!(!output.contains("<svg"));
    }

    #[test]
    fn renders_ti_report_with_inline_plots() {
        let output = render_ti_html_report(
            &sample_ti_advice(),
            AnalysisSampleCounts {
                windows: 4,
                samples_in: 16,
                samples_after_burnin: 16,
                samples_kept: 16,
            },
            &AnalysisInputOptions {
                temperature: 300.0,
                decorrelate: false,
                remove_burnin: 0,
                auto_equilibrate: false,
                fast: false,
                conservative: true,
                nskip: 1,
                u_nk_observable: None,
            },
            &AdviseRunOptions {
                estimator: AdvisorEstimatorArg::Mbar,
                input_kind: AdviseInputKind::Dhdl,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 4,
                suggest_midpoints: true,
                output_format: OutputFormat::Json,
                output_path: None,
                report_path: None,
            },
        )
        .unwrap();

        assert!(output.contains("TI Schedule Report"));
        assert!(output.contains("Plots"));
        assert!(output.contains("Mean dH/dλ"));
        assert!(output.contains("Curvature Magnitude"));
        assert!(output.contains("Interval Uncertainty"));
        assert!(output.contains("Integration Method Curves"));
        assert!(output.contains("Method Differences"));
        assert!(output.contains("Deviation From Trapezoidal"));
        assert!(output.contains("∫Δ dλ"));
        assert!(output.contains("series-fill positive"));
        assert!(output.contains(">0</text>"));
        assert!(output.contains("Trapezoidal"));
        assert!(output.contains("Cubic Spline"));
        assert!(output.contains("PCHIP"));
        assert!(output.contains("Akima"));
        assert!(!output.contains("Gaussian Quadrature"));
        assert!(output.contains("ti-series-plot"));
        assert!(output.contains("λ midpoint"));
        assert!(output.contains("|curvature|"));
        assert!(output.contains("series-line curvature"));
        assert!(output.contains("series-line uncertainty"));
        assert!(output.contains("href=\"#interval-"));
        assert!(output.contains("id=\"interval-"));
        assert!(output.contains("<a class=\"queue-link mono\" href=\"#interval-"));
        assert!(output.contains("<th>Interval</th><th>Endpoints</th><th>Slope</th><th>Curvature</th><th>Priority</th><th>Recommendation</th>"));
        assert!(output.contains(": λ "));
        assert!(output.contains("<span class=\"mono\"> interval 0: λ "));
        assert!(output.contains("high_curvature"));
        assert!(output.contains("proposed λ</div><div class=\"mono\">0.495"));
        assert!(!output.contains("proposed λ</div><div class=\"mono\">n/a</div>"));
    }

    #[test]
    fn renders_ti_sampling_evidence_for_extend_sampling_suggestions() {
        let output = render_ti_html_report(
            &sample_ti_sampling_advice(),
            AnalysisSampleCounts {
                windows: 3,
                samples_in: 12,
                samples_after_burnin: 12,
                samples_kept: 12,
            },
            &AnalysisInputOptions {
                temperature: 300.0,
                decorrelate: false,
                remove_burnin: 0,
                auto_equilibrate: false,
                fast: false,
                conservative: true,
                nskip: 1,
                u_nk_observable: None,
            },
            &AdviseRunOptions {
                estimator: AdvisorEstimatorArg::Mbar,
                input_kind: AdviseInputKind::Dhdl,
                overlap_min: 0.03,
                block_cv_min: 0.15,
                n_blocks: 2,
                suggest_midpoints: true,
                output_format: OutputFormat::Json,
                output_path: None,
                report_path: None,
            },
        )
        .unwrap();

        assert!(output.contains("left block CV"));
        assert!(output.contains("right block CV"));
        assert!(output.contains("forward/reverse delta"));
        assert!(output.contains("Integration Method Curves"));
        assert!(output.contains("Method Differences"));
        assert!(output.contains("∫Δ dλ"));
        assert!(output.contains("Simpson"));
        assert!(output.contains("high_block_cv"));
        assert!(!output.contains("proposed λ</div><div class=\"mono\">n/a</div>"));
    }

    #[test]
    fn ti_method_plots_only_mark_original_window_means() {
        let output = render_ti_method_plot_cards_html(&sample_ti_advice());

        assert_eq!(output.matches("<polyline class=\"series-line\"").count(), 4);
        assert_eq!(output.matches("<circle class=\"series-point\"").count(), 16);
        assert!(output.contains("series-fill positive"));
    }

    #[test]
    fn axis_tick_values_include_zero_when_requested() {
        let ticks = axis_tick_values(-2.0, 5.0, true);

        assert!(ticks.contains(&0.0));
    }

    #[test]
    fn zero_area_polygons_split_sign_changes_at_crossings() {
        let points = vec![(0.0, -1.0), (1.0, 1.0)];

        let positive = zero_area_polygons(&points, true);
        let negative = zero_area_polygons(&points, false);

        assert_eq!(positive.len(), 1);
        assert_eq!(negative.len(), 1);
        assert!((positive[0][0].0 - 0.5).abs() < 1.0e-12);
        assert!((negative[0][2].0 - 0.5).abs() < 1.0e-12);
    }
}
