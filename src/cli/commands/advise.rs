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
    pub report_path: Option<PathBuf>,
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
        lambda_components.clone(),
    )
    .map_err(io::Error::other)?;

    if let Some(path) = run_options.report_path.as_deref() {
        let report = render_html_report(
            &advice,
            loaded.sample_counts,
            &input_options,
            &run_options,
            lambda_components,
        )
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
.metric-row{display:grid;grid-template-columns:160px 1fr auto;gap:12px;align-items:center;margin-top:10px}.metric-name{font-size:13px;color:var(--muted)}\
.bar{height:10px;border-radius:999px;background:#eee6d8;overflow:hidden}.fill{height:100%;background:linear-gradient(90deg,var(--accent),#8db8a3)}.fill.bad{background:linear-gradient(90deg,#c46d52,var(--bad))}\
.mono{font:500 13px/1.45 ui-monospace,Consolas,monospace}.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px;margin-top:14px}.kv div{padding-top:10px;border-top:1px solid var(--line)}\
.muted{color:var(--muted)}.empty{color:var(--muted);font-style:italic}.table{display:grid;gap:12px}.suggestion-title{display:flex;gap:10px;align-items:center;flex-wrap:wrap}.footer{margin-top:28px;color:var(--muted);font-size:13px}\
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

    html.push_str(&format!(
        "<header class=\"hero\"><div class=\"eyebrow\">alchemrs schedule advisor</div><h1>Lambda Schedule Report</h1>\
<div class=\"lede\">Threshold-driven adjacent-edge diagnostics with neighbor-aware ranking and component-aware insertion proposals.</div></header>"
    ));

    html.push_str("<section class=\"grid summary\">");
    html.push_str(&summary_card(
        "Estimator",
        estimator_name(run_options.estimator),
        "adjacent edge fit method",
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
        &run_options.overlap_min.to_string(),
    ));
    html.push_str(&kv_html(
        "block CV min",
        &run_options.block_cv_min.to_string(),
    ));
    html.push_str(&kv_html("n blocks", &run_options.n_blocks.to_string()));
    html.push_str(&kv_html(
        "lambda components",
        &lambda_components
            .as_ref()
            .map(|labels| labels.join(", "))
            .unwrap_or_default(),
    ));
    html.push_str("</div></div></section>");

    html.push_str(
        "<section class=\"section\"><div class=\"card\"><h2>Legend</h2><div class=\"legend-grid\">",
    );
    html.push_str(
        "<div class=\"legend-item\"><div class=\"legend-title\">Axis Markers</div>\
<div class=\"legend-row\"><span class=\"swatch circle source\"></span><span class=\"mono\">source sampled state</span></div>\
<div class=\"legend-row\"><span class=\"swatch circle target\"></span><span class=\"mono\">target sampled state</span></div>\
<div class=\"legend-row\"><span class=\"swatch proposal\"></span><span class=\"mono\">proposed insertion state</span></div></div>",
    );
    html.push_str(
        "<div class=\"legend-item\"><div class=\"legend-title\">Delta Bars</div>\
<div class=\"legend-row\"><span class=\"swatch track\"></span><span class=\"mono\">component delta</span></div>\
<div class=\"legend-row\"><span class=\"swatch track focus\"></span><span class=\"mono\">highlighted dominant/focus component</span></div></div>",
    );
    html.push_str(
        "<div class=\"legend-item\"><div class=\"legend-title\">Status Badges</div>\
<div class=\"legend-row\"><span class=\"status split\">bisect</span><span class=\"mono\">proposal moves this component</span></div>\
<div class=\"legend-row\"><span class=\"status held\">held</span><span class=\"mono\">proposal keeps source value</span></div>\
<div class=\"legend-row\"><span class=\"status dominant\">dominant</span><span class=\"mono\">largest edge-local change</span></div>\
<div class=\"legend-row\"><span class=\"status fixed\">fixed</span><span class=\"mono\">component does not change</span></div></div>",
    );
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
            "<tr><td class=\"mono\"><a class=\"queue-link\" href=\"#edge-{}\">edge {}</a></td><td class=\"mono\">{:.3}</td><td><span class=\"pill {}\">{}</span></td><td class=\"mono\">{}</td><td class=\"mono\">{:.6}</td><td class=\"mono\">{} -&gt; {}</td></tr>",
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
            edge.overlap_min(),
            escape_html(&format_state_text(edge.from_state().lambdas())),
            escape_html(&format_state_text(edge.to_state().lambdas()))
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
            html.push_str("<div class=\"suggestion-title\">");
            html.push_str(&format!(
                "<span class=\"pill {}\">{}</span><span class=\"mono\">edge {}</span>",
                severity_class_from_suggestion(suggestion.kind()),
                escape_html(suggestion_name(suggestion.kind())),
                suggestion.edge_index()
            ));
            if let Some(strategy) = suggestion.proposal_strategy() {
                html.push_str(&format!(
                    "<span class=\"pill\">{}</span>",
                    escape_html(proposal_strategy_name(strategy))
                ));
            }
            html.push_str("</div>");
            html.push_str(&format!(
                "<p>{}</p><div class=\"metric-row\"><div class=\"metric-name\">priority</div><div class=\"bar\"><div class=\"fill bad\" style=\"width:{:.2}%\"></div></div><div class=\"mono\">{:.3}</div></div>",
                escape_html(suggestion.reason()),
                width,
                suggestion.priority_score()
            ));
            html.push_str("<div class=\"viz\">");
            html.push_str(&render_lambda_axis_svg(
                suggestion.from_state().lambdas(),
                suggestion.to_state().lambdas(),
                suggestion.proposed_state().map(|state| state.lambdas()),
                lambda_components.as_deref(),
                suggestion.focus_components(),
            ));
            html.push_str(&render_component_breakdown(
                suggestion.from_state().lambdas(),
                suggestion.to_state().lambdas(),
                suggestion.proposed_state().map(|state| state.lambdas()),
                lambda_components.as_deref(),
                suggestion.focus_components(),
            ));
            html.push_str("</div>");
            html.push_str("<div class=\"kv\">");
            html.push_str(&kv_html(
                "from lambda",
                &format_state_text(suggestion.from_state().lambdas()),
            ));
            html.push_str(&kv_html(
                "to lambda",
                &format_state_text(suggestion.to_state().lambdas()),
            ));
            html.push_str(&kv_html(
                "focus components",
                &suggestion.focus_components().join(", "),
            ));
            html.push_str(&kv_html(
                "proposed lambda",
                &suggestion
                    .proposed_state()
                    .map(|state| format_state_text(state.lambdas()))
                    .unwrap_or_default(),
            ));
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
            "<div class=\"suggestion-title\"><span class=\"pill {}\">{}</span><span class=\"mono\">edge {}</span><span class=\"mono\">{} → {}</span></div>",
            severity_name(edge.severity()),
            escape_html(severity_name(edge.severity())),
            edge.edge_index(),
            escape_html(&format_state_text(edge.from_state().lambdas())),
            escape_html(&format_state_text(edge.to_state().lambdas()))
        ));
        html.push_str(&format!(
            "<div class=\"metric-row\"><div class=\"metric-name\">priority</div><div class=\"bar\"><div class=\"fill {}\" style=\"width:{:.2}%\"></div></div><div class=\"mono\">{:.3}</div></div>",
            if matches!(edge.severity(), EdgeSeverity::AddSampling | EdgeSeverity::AddWindow) { "bad" } else { "" },
            width,
            edge.priority_score()
        ));
        html.push_str("<div class=\"viz\">");
        html.push_str(&render_lambda_axis_svg(
            edge.from_state().lambdas(),
            edge.to_state().lambdas(),
            None,
            lambda_components.as_deref(),
            edge.dominant_components(),
        ));
        html.push_str(&render_component_breakdown(
            edge.from_state().lambdas(),
            edge.to_state().lambdas(),
            None,
            lambda_components.as_deref(),
            edge.dominant_components(),
        ));
        html.push_str("</div>");
        html.push_str("<div class=\"kv\">");
        html.push_str(&kv_html(
            "overlap min",
            &format!("{:.6}", edge.overlap_min()),
        ));
        html.push_str(&kv_html(
            "relative overlap",
            &option_string(edge.relative_overlap()),
        ));
        html.push_str(&kv_html("delta_f (kT)", &format!("{:.6}", edge.delta_f())));
        html.push_str(&kv_html("uncertainty", &option_string(edge.uncertainty())));
        html.push_str(&kv_html(
            "relative uncertainty",
            &option_string(edge.relative_uncertainty()),
        ));
        html.push_str(&kv_html("block CV", &option_string(edge.block_cv())));
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

fn format_string_list(values: &[String]) -> String {
    if values.is_empty() {
        String::new()
    } else {
        format!("[{}]", values.join(";"))
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

fn render_lambda_axis_svg(
    from: &[f64],
    to: &[f64],
    proposed: Option<&[f64]>,
    labels: Option<&[String]>,
    focus_components: &[String],
) -> String {
    let width = 420.0;
    let left = 120.0;
    let right = width - 18.0;
    let row_height = 34.0;
    let height = 24.0 + row_height * from.len() as f64;
    let mut svg = format!(
        "<svg viewBox=\"0 0 {width} {height}\" xmlns=\"http://www.w3.org/2000/svg\" role=\"img\" aria-label=\"lambda axis visualization\">"
    );

    for (index, (&from_value, &to_value)) in from.iter().zip(to.iter()).enumerate() {
        let y = 24.0 + index as f64 * row_height;
        let label = labels
            .and_then(|labels| labels.get(index))
            .cloned()
            .unwrap_or_else(|| format!("lambda[{index}]"));
        let is_focus = focus_components.iter().any(|value| value == &label);
        let proposed_value = proposed.and_then(|values| values.get(index)).copied();

        svg.push_str(&format!(
            "<text x=\"10\" y=\"{:.1}\" font-family=\"Consolas, monospace\" font-size=\"11\" fill=\"{}\">{}</text>",
            y + 4.0,
            if is_focus { "#1f5e5b" } else { "#6a655e" },
            escape_html(&label)
        ));
        svg.push_str(&format!(
            "<line x1=\"{left}\" y1=\"{y}\" x2=\"{right}\" y2=\"{y}\" stroke=\"{}\" stroke-width=\"3\" stroke-linecap=\"round\"/>",
            if is_focus { "#8db8a3" } else { "#d7cfc0" }
        ));
        svg.push_str(&format!(
            "<circle cx=\"{:.1}\" cy=\"{y}\" r=\"5\" fill=\"#1f5e5b\"/>",
            axis_position(from_value, left, right)
        ));
        svg.push_str(&format!(
            "<circle cx=\"{:.1}\" cy=\"{y}\" r=\"5\" fill=\"#b5483d\"/>",
            axis_position(to_value, left, right)
        ));
        if let Some(value) = proposed_value {
            let x = axis_position(value, left, right);
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"10\" height=\"10\" rx=\"2\" fill=\"#b9832f\" transform=\"rotate(45 {:.1} {:.1})\"/>",
                x - 5.0,
                y - 5.0,
                x,
                y
            ));
        }
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"Consolas, monospace\" font-size=\"10\" fill=\"#6a655e\">{:.3}</text>",
            left - 6.0,
            y - 9.0,
            from_value
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"Consolas, monospace\" font-size=\"10\" fill=\"#6a655e\" text-anchor=\"end\">{:.3}</text>",
            right + 6.0,
            y - 9.0,
            to_value
        ));
    }

    svg.push_str("</svg>");
    svg
}

fn render_component_breakdown(
    from: &[f64],
    to: &[f64],
    proposed: Option<&[f64]>,
    labels: Option<&[String]>,
    highlight_components: &[String],
) -> String {
    if from.len() <= 1 {
        return String::new();
    }
    let max_delta = from
        .iter()
        .zip(to.iter())
        .map(|(from_value, to_value)| (to_value - from_value).abs())
        .fold(0.0, f64::max)
        .max(1.0e-12);

    let mut html = String::from(
        "<div class=\"component-grid\"><div class=\"component-row component-head\"><div>component</div><div>from</div><div>proposal</div><div>to</div><div>delta</div><div>status</div></div>",
    );
    for (index, (&from_value, &to_value)) in from.iter().zip(to.iter()).enumerate() {
        let label = component_label(index, labels);
        let delta = (to_value - from_value).abs();
        let delta_width = (delta / max_delta * 100.0).clamp(0.0, 100.0);
        let proposed_value = proposed.and_then(|values| values.get(index)).copied();
        let is_highlight = highlight_components.iter().any(|value| value == &label);
        let (status_class, status_label) =
            component_status(from_value, to_value, proposed_value, is_highlight);
        html.push_str(&format!(
            "<div class=\"component-row{}\"><div class=\"component-name\">{}</div><div class=\"component-value\">{:.3}</div><div class=\"component-value\">{}</div><div class=\"component-value\">{:.3}</div><div><div class=\"delta-track\"><div class=\"delta-fill{}\" style=\"width:{:.2}%\"></div></div><div class=\"delta-value\">{:.3}</div></div><div><span class=\"status {}\">{}</span></div></div>",
            if is_highlight { " focus" } else { "" },
            escape_html(&label),
            from_value,
            proposed_value
                .map(|value| format!("{value:.3}"))
                .unwrap_or_else(|| "&mdash;".to_string()),
            to_value,
            if is_highlight { " focus" } else { "" },
            delta_width,
            delta,
            status_class,
            status_label,
        ));
    }
    html.push_str("</div>");
    html
}

fn axis_position(value: f64, left: f64, right: f64) -> f64 {
    let normalized = value.clamp(0.0, 1.0);
    left + normalized * (right - left)
}

fn component_label(index: usize, labels: Option<&[String]>) -> String {
    labels
        .and_then(|labels| labels.get(index))
        .cloned()
        .unwrap_or_else(|| format!("lambda[{index}]"))
}

fn component_status(
    from: f64,
    to: f64,
    proposed: Option<f64>,
    is_highlight: bool,
) -> (&'static str, &'static str) {
    let delta = (to - from).abs();
    if delta <= 1.0e-12 {
        return ("fixed", "fixed");
    }
    if let Some(proposed) = proposed {
        if (proposed - from).abs() <= 1.0e-12 {
            return ("held", "held");
        }
        return ("split", "bisect");
    }
    if is_highlight {
        ("dominant", "dominant")
    } else {
        ("changed", "changed")
    }
}

#[cfg(test)]
mod tests {
    use alchemrs::{
        advise_lambda_schedule, extract_u_nk, AdvisorEstimator, ScheduleAdvisorOptions, StatePoint,
        UNkMatrix,
    };
    use serde_json::Value;

    use super::{render_advice, render_html_report, AdviseRunOptions};
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
        assert!(output.contains("Legend"));
        assert!(output.contains("Priority Queue"));
        assert!(output.contains("Weakest Component"));
        assert!(output.contains("href=\"#edge-"));
        assert!(output.contains("id=\"edge-"));
        assert!(output.contains("source sampled state"));
        assert!(output.contains("focused_split") || output.contains("midpoint"));
        assert!(output.contains("priority"));
        assert!(output.contains("<svg"));
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

        assert!(output.contains("component-grid"));
        assert!(output.contains("highlighted dominant/focus component"));
        assert!(output.contains("Priority Queue"));
        assert!(output.contains("href=\"#edge-0\"") || output.contains("href=\"#edge-1\""));
        assert!(output.contains("coul-lambda"));
        assert!(output.contains("vdw-lambda"));
        assert!(output.contains("delta-track"));
        assert!(output.contains("0.800"));
        assert!(output.contains("bisect"));
        assert!(output.contains("held"));
    }
}
