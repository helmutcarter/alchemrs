use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use alchemrs::parse::amber::{extract_dhdl, extract_u_nk, extract_u_nk_with_potential};
use alchemrs::{
    advise_lambda_schedule, advise_ti_schedule, decorrelate_dhdl, decorrelate_u_nk,
    decorrelate_u_nk_with_observable, detect_equilibration_dhdl, detect_equilibration_u_nk,
    recommend_ti_method, AdvisorEstimator, DecorrelationOptions, ScheduleAdvisorOptions,
    TiScheduleAdvisorOptions, UNkSeriesMethod,
};
use alchemrs::{
    overlap_eigenvalues, overlap_matrix, BarEstimator, BarMethod, BarOptions, DhdlSeries,
    ExpEstimator, ExpOptions, MbarEstimator, MbarOptions, TiEstimator, TiOptions, UNkMatrix,
};
use serde_json::Value;

const TEMPERATURE_K: f64 = 300.0;

#[test]
fn advise_schedule_cli_outputs_expected_json_when_decorrelating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "advise-schedule",
            "--decorrelate",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows(&inputs);
    let advice = advise_lambda_schedule(
        &windows,
        Some(ScheduleAdvisorOptions {
            estimator: AdvisorEstimator::Mbar,
            overlap_min: 0.03,
            block_cv_min: 0.15,
            n_blocks: 4,
            suggest_midpoints: true,
        }),
    )
    .expect("compute schedule advice");
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    let sample_counts = payload["sample_counts"].as_object().expect("sample_counts");
    assert_eq!(
        sample_counts["windows"].as_u64(),
        Some(expected_counts.windows as u64)
    );
    assert_eq!(
        sample_counts["samples_in"].as_u64(),
        Some(expected_counts.samples_in as u64)
    );
    assert_eq!(
        sample_counts["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        sample_counts["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );

    assert_eq!(
        payload["provenance"]["advisor_estimator"].as_str(),
        Some("mbar")
    );
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
    assert_close(
        payload["provenance"]["overlap_min"]
            .as_f64()
            .expect("overlap_min"),
        0.03,
    );
    assert_close(
        payload["provenance"]["block_cv_min"]
            .as_f64()
            .expect("block_cv_min"),
        0.15,
    );
    assert_eq!(payload["provenance"]["n_blocks"].as_u64(), Some(4));

    let edges = payload["edges"].as_array().expect("edges");
    assert_eq!(edges.len(), advice.edges().len());
    for (actual, expected) in edges.iter().zip(advice.edges().iter()) {
        assert_eq!(
            actual["edge_index"].as_u64(),
            Some(expected.edge_index() as u64)
        );
        assert_eq!(
            actual["from_lambda"].as_f64(),
            Some(expected.from_state().lambdas()[0])
        );
        assert_eq!(
            actual["to_lambda"].as_f64(),
            Some(expected.to_state().lambdas()[0])
        );
        assert_close(
            actual["overlap_min"].as_f64().expect("overlap_min"),
            expected.overlap_min(),
        );
        if let Some(value) = expected.relative_overlap() {
            assert_close(
                actual["relative_overlap"]
                    .as_f64()
                    .expect("relative_overlap"),
                value,
            );
        } else {
            assert!(actual["relative_overlap"].is_null());
        }
        assert_close(
            actual["delta_f"].as_f64().expect("delta_f"),
            expected.delta_f(),
        );
        if let Some(value) = expected.relative_uncertainty() {
            assert_close(
                actual["relative_uncertainty"]
                    .as_f64()
                    .expect("relative_uncertainty"),
                value,
            );
        } else {
            assert!(actual["relative_uncertainty"].is_null());
        }
        assert_close(
            actual["priority_score"].as_f64().expect("priority_score"),
            expected.priority_score(),
        );
        assert_eq!(
            actual["severity"].as_str(),
            Some(schedule_severity_name(expected.severity()))
        );
        assert_eq!(
            actual["dominant_components"],
            Value::Array(
                expected
                    .dominant_components()
                    .iter()
                    .cloned()
                    .map(Value::from)
                    .collect()
            )
        );
    }

    let suggestions = payload["suggestions"].as_array().expect("suggestions");
    assert_eq!(suggestions.len(), advice.suggestions().len());
    for (actual, expected) in suggestions.iter().zip(advice.suggestions().iter()) {
        assert_eq!(
            actual["edge_index"].as_u64(),
            Some(expected.edge_index() as u64)
        );
        assert_eq!(
            actual["kind"].as_str(),
            Some(schedule_suggestion_name(expected.kind()))
        );
        assert_eq!(actual["reason"].as_str(), Some(expected.reason()));
        assert_close(
            actual["priority_score"].as_f64().expect("priority_score"),
            expected.priority_score(),
        );
        assert_eq!(
            actual["focus_components"],
            Value::Array(
                expected
                    .focus_components()
                    .iter()
                    .cloned()
                    .map(Value::from)
                    .collect()
            )
        );
        assert_eq!(
            actual["proposal_strategy"],
            expected
                .proposal_strategy()
                .map(schedule_proposal_strategy_name)
                .map(Value::from)
                .unwrap_or(Value::Null)
        );
        assert_eq!(
            actual["proposed_lambda"],
            expected
                .proposed_state()
                .map(|state| Value::from(state.lambdas()[0]))
                .unwrap_or(Value::Null)
        );
    }
}

#[test]
fn advise_schedule_cli_outputs_expected_ti_json_when_forced_to_dhdl() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "advise-schedule",
            "--input-kind",
            "dhdl",
            "--decorrelate",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let series = load_decorrelated_dhdl_series(&inputs);
    let advice = advise_ti_schedule(
        &series,
        Some(TiScheduleAdvisorOptions {
            n_blocks: 4,
            block_cv_min: 0.15,
            ..TiScheduleAdvisorOptions::default()
        }),
    )
    .expect("compute TI schedule advice");
    let expected_counts = expected_dhdl_counts(&inputs, &series);

    assert_eq!(payload["advisor_mode"].as_str(), Some("ti"));
    assert_eq!(payload["provenance"]["input_kind"].as_str(), Some("dhdl"));
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert_close(
        payload["provenance"]["block_cv_min"]
            .as_f64()
            .expect("block_cv_min"),
        0.15,
    );
    assert_eq!(payload["provenance"]["n_blocks"].as_u64(), Some(4));

    let sample_counts = payload["sample_counts"].as_object().expect("sample_counts");
    assert_eq!(
        sample_counts["windows"].as_u64(),
        Some(expected_counts.windows as u64)
    );
    assert_eq!(
        sample_counts["samples_in"].as_u64(),
        Some(expected_counts.samples_in as u64)
    );
    assert_eq!(
        sample_counts["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        sample_counts["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );

    let windows = payload["windows"].as_array().expect("windows");
    assert_eq!(windows.len(), advice.windows().len());
    for (actual, expected) in windows.iter().zip(advice.windows().iter()) {
        assert_eq!(
            actual["window_index"].as_u64(),
            Some(expected.window_index() as u64)
        );
        assert_eq!(actual["lambda"], json_state(expected.state().lambdas()));
        assert_close(
            actual["mean_dhdl"].as_f64().expect("mean_dhdl"),
            expected.mean_dhdl(),
        );
        assert_json_optional_f64_close(&actual["sem_dhdl"], expected.sem_dhdl());
        assert_json_optional_f64_close(&actual["block_cv"], expected.block_cv());
        assert_json_optional_f64_close(
            &actual["split_half_dhdl_delta"],
            expected.split_half_dhdl_delta(),
        );
    }

    let intervals = payload["intervals"].as_array().expect("intervals");
    assert_eq!(intervals.len(), advice.intervals().len());
    for (actual, expected) in intervals.iter().zip(advice.intervals().iter()) {
        assert_eq!(
            actual["interval_index"].as_u64(),
            Some(expected.interval_index() as u64)
        );
        assert_eq!(
            actual["from_lambda"],
            json_state(expected.from_state().lambdas())
        );
        assert_eq!(
            actual["to_lambda"],
            json_state(expected.to_state().lambdas())
        );
        assert_close(
            actual["delta_lambda"].as_f64().expect("delta_lambda"),
            expected.delta_lambda(),
        );
        assert_close(
            actual["left_mean_dhdl"].as_f64().expect("left_mean_dhdl"),
            expected.left_mean_dhdl(),
        );
        assert_close(
            actual["right_mean_dhdl"].as_f64().expect("right_mean_dhdl"),
            expected.right_mean_dhdl(),
        );
        assert_close(
            actual["trapezoid_contribution"]
                .as_f64()
                .expect("trapezoid_contribution"),
            expected.trapezoid_contribution(),
        );
        assert_close(actual["slope"].as_f64().expect("slope"), expected.slope());
        assert_close(
            actual["abs_slope"].as_f64().expect("abs_slope"),
            expected.abs_slope(),
        );
        assert_json_optional_f64_close(&actual["curvature"], expected.curvature());
        assert_json_optional_f64_close(
            &actual["interval_uncertainty"],
            expected.interval_uncertainty(),
        );
        assert_eq!(
            actual["severity"].as_str(),
            Some(ti_schedule_severity_name(expected.severity()))
        );
        assert_close(
            actual["priority_score"].as_f64().expect("priority_score"),
            expected.priority_score(),
        );
    }

    let suggestions = payload["suggestions"].as_array().expect("suggestions");
    assert_eq!(suggestions.len(), advice.suggestions().len());
    for (actual, expected) in suggestions.iter().zip(advice.suggestions().iter()) {
        assert_eq!(
            actual["interval_index"].as_u64(),
            Some(expected.interval_index() as u64)
        );
        assert_eq!(
            actual["kind"].as_str(),
            Some(ti_schedule_suggestion_name(expected.kind()))
        );
        assert_eq!(
            actual["from_lambda"],
            json_state(expected.from_state().lambdas())
        );
        assert_eq!(
            actual["to_lambda"],
            json_state(expected.to_state().lambdas())
        );
        assert_eq!(
            actual["proposed_lambda"],
            expected
                .proposed_state()
                .map(|state| json_state(state.lambdas()))
                .unwrap_or(Value::Null)
        );
        assert_close(
            actual["priority_score"].as_f64().expect("priority_score"),
            expected.priority_score(),
        );
        assert_eq!(actual["reason"].as_str(), Some(expected.reason()));
    }
}

#[test]
fn advise_schedule_cli_writes_html_report() {
    let inputs = acetamide_inputs();
    let report_path = unique_output_path("advise-schedule-report.html");
    let report_path_string = report_path.to_string_lossy().to_string();

    let output = run_cli(
        &[
            "advise-schedule",
            "--decorrelate",
            "--report",
            &report_path_string,
            "--output-format",
            "json",
        ],
        &inputs,
    );

    let payload = parse_json_output(&output);
    assert!(payload["edges"].is_array());

    let written = fs::read_to_string(&report_path).expect("read HTML report");
    assert!(written.contains("Lambda Schedule Report"));
    assert!(written.contains("Suggestions"));
    assert!(written.contains("Edges"));
    assert!(written.contains("priority"));
    assert!(!written.contains("<svg"));

    fs::remove_file(&report_path).expect("remove HTML report");
}

fn schedule_severity_name(severity: alchemrs::EdgeSeverity) -> &'static str {
    match severity {
        alchemrs::EdgeSeverity::Healthy => "healthy",
        alchemrs::EdgeSeverity::Monitor => "monitor",
        alchemrs::EdgeSeverity::AddSampling => "add_sampling",
        alchemrs::EdgeSeverity::AddWindow => "add_window",
    }
}

fn schedule_suggestion_name(kind: alchemrs::SuggestionKind) -> &'static str {
    match kind {
        alchemrs::SuggestionKind::NoChange => "no_change",
        alchemrs::SuggestionKind::ExtendSampling => "extend_sampling",
        alchemrs::SuggestionKind::InsertWindow => "insert_window",
    }
}

fn schedule_proposal_strategy_name(strategy: alchemrs::ProposalStrategy) -> &'static str {
    match strategy {
        alchemrs::ProposalStrategy::Midpoint => "midpoint",
        alchemrs::ProposalStrategy::FocusedSplit => "focused_split",
    }
}

fn ti_schedule_severity_name(severity: alchemrs::TiEdgeSeverity) -> &'static str {
    match severity {
        alchemrs::TiEdgeSeverity::Healthy => "healthy",
        alchemrs::TiEdgeSeverity::Monitor => "monitor",
        alchemrs::TiEdgeSeverity::AddSampling => "add_sampling",
        alchemrs::TiEdgeSeverity::AddWindow => "add_window",
        alchemrs::TiEdgeSeverity::AddWindowAndSampling => "add_window_and_sampling",
    }
}

fn ti_schedule_suggestion_name(kind: alchemrs::TiSuggestionKind) -> &'static str {
    match kind {
        alchemrs::TiSuggestionKind::NoChange => "no_change",
        alchemrs::TiSuggestionKind::ExtendSampling => "extend_sampling",
        alchemrs::TiSuggestionKind::InsertWindow => "insert_window",
        alchemrs::TiSuggestionKind::InsertWindowAndExtendSampling => {
            "insert_window_and_extend_sampling"
        }
    }
}

#[test]
fn mbar_cli_outputs_expected_json_with_overlap_summary_when_decorrelating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "mbar",
            "--decorrelate",
            "--output-format",
            "json",
            "--overlap-summary",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows(&inputs);
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: 10_000,
        tolerance: 1.0e-7,
        parallel: false,
        ..MbarOptions::default()
    });
    let result = estimator
        .estimate_with_uncertainty(&windows)
        .expect("fit MBAR");
    let delta_index = result.n_states() - 1;
    let expected_overlap = compute_overlap_summary(&windows);
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_close(payload["from_lambda"].as_f64().expect("from_lambda"), 0.0);
    assert_close(payload["to_lambda"].as_f64().expect("to_lambda"), 1.0);
    assert_eq!(payload["units"].as_str(), Some("kT"));
    assert_eq!(payload["provenance"]["estimator"].as_str(), Some("mbar"));
    assert_close(
        payload["provenance"]["temperature_k"]
            .as_f64()
            .expect("temperature_k"),
        TEMPERATURE_K,
    );
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["remove_burnin"].as_u64(), Some(0));
    assert_eq!(
        payload["provenance"]["auto_equilibrate"].as_bool(),
        Some(false)
    );
    assert_eq!(payload["provenance"]["fast"].as_bool(), Some(false));
    assert_eq!(payload["provenance"]["conservative"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["nskip"].as_u64(), Some(1));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
    assert_eq!(
        payload["provenance"]["windows"].as_u64(),
        Some(expected_counts.windows as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_in"].as_u64(),
        Some(expected_counts.samples_in as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );

    let overlap = payload["overlap"].as_object().expect("overlap object");
    assert_close(
        overlap["scalar"].as_f64().expect("overlap scalar"),
        expected_overlap.0,
    );
    let eigenvalues = overlap["eigenvalues"]
        .as_array()
        .expect("overlap eigenvalues");
    assert_eq!(eigenvalues.len(), expected_overlap.1.len());
    for (actual, expected) in eigenvalues.iter().zip(expected_overlap.1.iter()) {
        assert_close(actual.as_f64().expect("eigenvalue"), *expected);
    }
}

#[test]
fn mbar_cli_outputs_expected_json_when_auto_equilibrating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &["mbar", "--auto-equilibrate", "--output-format", "json"],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_auto_equilibrated_windows(&inputs);
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: 10_000,
        tolerance: 1.0e-7,
        parallel: false,
        ..MbarOptions::default()
    });
    let result = estimator
        .estimate_with_uncertainty(&windows)
        .expect("fit MBAR");
    let delta_index = result.n_states() - 1;
    let expected_counts = expected_u_nk_counts_after_auto_equilibration(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_eq!(
        payload["provenance"]["auto_equilibrate"].as_bool(),
        Some(true)
    );
    assert_eq!(payload["provenance"]["fast"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["conservative"].as_bool(), Some(false));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
    assert_eq!(
        payload["provenance"]["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
}

#[test]
fn mbar_cli_accepts_nonconservative_decorrelation() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "mbar",
            "--decorrelate",
            "--conservative=false",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows_nonconservative(&inputs);
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: 10_000,
        tolerance: 1.0e-7,
        parallel: false,
        ..MbarOptions::default()
    });
    let result = estimator
        .estimate_with_uncertainty(&windows)
        .expect("fit MBAR");
    let delta_index = result.n_states() - 1;
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["conservative"].as_bool(), Some(false));
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
}

#[test]
fn bar_cli_outputs_expected_json_with_overlap_summary_when_decorrelating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "bar",
            "--decorrelate",
            "--output-format",
            "json",
            "--overlap-summary",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows(&inputs);
    let estimator = BarEstimator::new(BarOptions {
        method: BarMethod::FalsePosition,
        parallel: false,
        ..BarOptions::default()
    });
    let fit = estimator.fit(&windows).expect("fit BAR");
    let result = fit.result().expect("BAR result");
    let delta_index = result.n_states() - 1;
    let expected_overlap = compute_overlap_summary(&windows);
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    if let Some(uncertainties) = result.uncertainties() {
        let sigma = uncertainties[delta_index];
        if sigma.is_finite() {
            assert_close(payload["uncertainty"].as_f64().expect("uncertainty"), sigma);
        } else {
            assert!(
                payload["uncertainty"].is_null(),
                "expected null uncertainty"
            );
        }
    } else {
        assert!(
            payload["uncertainty"].is_null(),
            "expected null uncertainty"
        );
    }
    assert_close(payload["from_lambda"].as_f64().expect("from_lambda"), 0.0);
    assert_close(payload["to_lambda"].as_f64().expect("to_lambda"), 1.0);
    assert_eq!(payload["units"].as_str(), Some("kT"));
    assert_eq!(payload["provenance"]["estimator"].as_str(), Some("bar"));
    assert_close(
        payload["provenance"]["temperature_k"]
            .as_f64()
            .expect("temperature_k"),
        TEMPERATURE_K,
    );
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["remove_burnin"].as_u64(), Some(0));
    assert_eq!(
        payload["provenance"]["auto_equilibrate"].as_bool(),
        Some(false)
    );
    assert_eq!(payload["provenance"]["fast"].as_bool(), Some(false));
    assert_eq!(payload["provenance"]["conservative"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["nskip"].as_u64(), Some(1));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
    assert_eq!(
        payload["provenance"]["windows"].as_u64(),
        Some(expected_counts.windows as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_in"].as_u64(),
        Some(expected_counts.samples_in as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );

    let overlap = payload["overlap"].as_object().expect("overlap object");
    assert_close(
        overlap["scalar"].as_f64().expect("overlap scalar"),
        expected_overlap.0,
    );
    let eigenvalues = overlap["eigenvalues"]
        .as_array()
        .expect("overlap eigenvalues");
    assert_eq!(eigenvalues.len(), expected_overlap.1.len());
    for (actual, expected) in eigenvalues.iter().zip(expected_overlap.1.iter()) {
        assert_close(actual.as_f64().expect("eigenvalue"), *expected);
    }
}

#[test]
fn ti_cli_outputs_expected_json_when_decorrelating() {
    let inputs = acetamide_inputs();
    let output = run_cli(&["ti", "--decorrelate", "--output-format", "json"], &inputs);
    let payload = parse_json_output(&output);

    let series = load_decorrelated_dhdl_series(&inputs);
    let estimator = TiEstimator::new(TiOptions {
        parallel: false,
        ..TiOptions::default()
    });
    let fit = estimator.fit(&series).expect("fit TI");
    let result = fit.result().expect("TI result");
    let expected_counts = expected_dhdl_counts(&inputs, &series);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.delta_f(),
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainty().expect("uncertainty"),
    );
    assert_close(
        payload["from_lambda"].as_f64().expect("from_lambda"),
        result.from_state().lambdas()[0],
    );
    assert_close(
        payload["to_lambda"].as_f64().expect("to_lambda"),
        result.to_state().lambdas()[0],
    );
    assert!(payload["overlap"].is_null(), "expected null overlap");
    assert_eq!(payload["provenance"]["estimator"].as_str(), Some("ti"));
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert!(payload["provenance"]["u_nk_observable"].is_null());
    assert_eq!(
        payload["provenance"]["ti_method"].as_str(),
        Some("trapezoidal")
    );
    assert!(payload["provenance"]["ti_method_reason"].is_null());
    assert_eq!(
        payload["provenance"]["windows"].as_u64(),
        Some(expected_counts.windows as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_in"].as_u64(),
        Some(expected_counts.samples_in as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
}

#[test]
fn ti_cli_outputs_expected_json_when_auto_equilibrating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &["ti", "--auto-equilibrate", "--output-format", "json"],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let series = load_auto_equilibrated_dhdl_series(&inputs);
    let estimator = TiEstimator::new(TiOptions {
        parallel: false,
        ..TiOptions::default()
    });
    let fit = estimator.fit(&series).expect("fit TI");
    let result = fit.result().expect("TI result");
    let expected_counts = expected_dhdl_counts_after_auto_equilibration(&inputs, &series);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.delta_f(),
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainty().expect("uncertainty"),
    );
    assert_eq!(
        payload["provenance"]["auto_equilibrate"].as_bool(),
        Some(true)
    );
    assert_eq!(payload["provenance"]["fast"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["conservative"].as_bool(), Some(false));
    assert!(payload["provenance"]["u_nk_observable"].is_null());
    assert_eq!(
        payload["provenance"]["ti_method"].as_str(),
        Some("trapezoidal")
    );
    assert!(payload["provenance"]["ti_method_reason"].is_null());
    assert_eq!(
        payload["provenance"]["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
}

#[test]
fn ti_cli_rejects_u_nk_observable_with_explanatory_error() {
    let inputs = acetamide_inputs();
    let output = run_cli_failure(
        &["ti", "--u-nk-observable", "de", "--output-format", "json"],
        &inputs,
    );

    assert!(!output.status.success(), "expected command to fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(
            "--u-nk-observable is only valid for bar, mbar, exp, and dexp; ti uses dH/dlambda for preprocessing."
        ),
        "unexpected stderr:\n{stderr}"
    );
}

#[test]
fn ti_cli_rejects_nonquadrature_schedule_for_gaussian_quadrature() {
    let inputs = acetamide_inputs();
    let output = run_cli_failure(
        &[
            "ti",
            "--method",
            "gaussian-quadrature",
            "--output-format",
            "json",
        ],
        &inputs,
    );

    assert!(!output.status.success(), "expected command to fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("use trapezoidal TI for arbitrary schedules"),
        "unexpected stderr:\n{stderr}"
    );
}

#[test]
fn ti_cli_auto_reports_recommended_method_in_json() {
    let inputs = acetamide_inputs();
    let output = run_cli(&["ti", "--method", "auto", "--output-format", "json"], &inputs);
    let payload = parse_json_output(&output);

    let series = inputs
        .iter()
        .map(|path| extract_dhdl(path, TEMPERATURE_K).expect("extract dhdl"))
        .collect::<Vec<_>>();
    let recommendation = recommend_ti_method(&series, None).expect("TI recommendation");
    let estimator = TiEstimator::new(TiOptions {
        method: recommendation.recommended_method(),
        parallel: false,
    });
    let fit = estimator.fit(&series).expect("fit TI");
    let result = fit.result().expect("TI result");

    assert_eq!(
        payload["provenance"]["ti_method"].as_str(),
        Some(recommendation.recommended_method().as_str())
    );
    assert_eq!(
        payload["provenance"]["ti_method_reason"].as_str(),
        Some(recommendation.reason())
    );
    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.delta_f(),
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainty().expect("uncertainty"),
    );
}

#[test]
fn exp_cli_outputs_expected_json_with_overlap_summary_when_decorrelating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "exp",
            "--decorrelate",
            "--output-format",
            "json",
            "--overlap-summary",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows(&inputs);
    let estimator = ExpEstimator::new(ExpOptions { parallel: false });
    let fit = estimator.fit(&windows).expect("fit EXP");
    let result = fit
        .result_with_uncertainty()
        .expect("EXP result with uncertainty");
    let delta_index = result.n_states() - 1;
    let expected_overlap = compute_overlap_summary(&windows);
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_close(payload["from_lambda"].as_f64().expect("from_lambda"), 0.0);
    assert_close(payload["to_lambda"].as_f64().expect("to_lambda"), 1.0);
    assert_eq!(payload["provenance"]["estimator"].as_str(), Some("exp"));
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
    assert_close(
        payload["overlap"]["scalar"]
            .as_f64()
            .expect("overlap scalar"),
        expected_overlap.0,
    );
}

#[test]
fn dexp_cli_outputs_expected_json_with_overlap_summary_when_decorrelating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "dexp",
            "--decorrelate",
            "--output-format",
            "json",
            "--overlap-summary",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows(&inputs);
    let estimator = ExpEstimator::new(ExpOptions { parallel: false });
    let fit = estimator.fit(&windows).expect("fit DEXP");
    let result = fit
        .result_with_uncertainty()
        .expect("DEXP result with uncertainty");
    let n_states = result.n_states();
    let delta_index = (n_states - 1) * n_states;
    let expected_overlap = compute_overlap_summary(&windows);
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_close(payload["from_lambda"].as_f64().expect("from_lambda"), 1.0);
    assert_close(payload["to_lambda"].as_f64().expect("to_lambda"), 0.0);
    assert_eq!(payload["provenance"]["estimator"].as_str(), Some("dexp"));
    assert_eq!(payload["provenance"]["decorrelate"].as_bool(), Some(true));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
    assert_close(
        payload["overlap"]["scalar"]
            .as_f64()
            .expect("overlap scalar"),
        expected_overlap.0,
    );
}

#[test]
fn mbar_cli_writes_json_to_output_file() {
    let inputs = acetamide_inputs();
    let output_path = unique_output_path("mbar-output.json");
    let output_path_string = output_path.to_string_lossy().to_string();
    let output = run_cli(
        &[
            "mbar",
            "--decorrelate",
            "--output-format",
            "json",
            "--output",
            &output_path_string,
        ],
        &inputs,
    );

    assert!(
        output.stdout.is_empty(),
        "expected stdout to be empty when writing to a file"
    );

    let written = fs::read(&output_path).expect("read CLI output file");
    let payload: Value = serde_json::from_slice(&written).expect("parse CLI JSON file");
    assert!(payload["delta_f"].as_f64().expect("delta_f").is_finite());
    assert_eq!(payload["units"].as_str(), Some("kT"));
    assert_eq!(payload["provenance"]["estimator"].as_str(), Some("mbar"));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
    assert!(
        payload["provenance"]["samples_kept"]
            .as_u64()
            .expect("samples_kept")
            > 0
    );

    fs::remove_file(&output_path).expect("remove CLI output file");
}

#[test]
fn mbar_cli_reports_effective_auto_equilibrate_settings() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "mbar",
            "--decorrelate",
            "--auto-equilibrate",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    assert_eq!(
        payload["provenance"]["auto_equilibrate"].as_bool(),
        Some(true)
    );
    assert_eq!(payload["provenance"]["fast"].as_bool(), Some(true));
    assert_eq!(payload["provenance"]["conservative"].as_bool(), Some(false));
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("de")
    );
}

#[test]
fn mbar_cli_auto_equilibrates_before_decorrelating() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "mbar",
            "--auto-equilibrate",
            "--decorrelate",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let auto_windows = load_auto_equilibrated_windows(&inputs);
    let windows = load_auto_equilibrated_decorrelated_windows(&inputs);
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: 10_000,
        tolerance: 1.0e-7,
        parallel: false,
        ..MbarOptions::default()
    });
    let result = estimator
        .estimate_with_uncertainty(&windows)
        .expect("fit MBAR");
    let delta_index = result.n_states() - 1;
    let expected_counts = expected_u_nk_counts_after_auto_equilibration_and_decorrelation(
        &inputs,
        &auto_windows,
        &windows,
    );

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_eq!(
        payload["provenance"]["samples_after_burnin"].as_u64(),
        Some(expected_counts.samples_after_burnin as u64)
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
}

#[test]
fn mbar_cli_supports_epot_observable_for_u_nk_preprocessing() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "mbar",
            "--decorrelate",
            "--u-nk-observable",
            "epot",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows_epot(&inputs);
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: 10_000,
        tolerance: 1.0e-7,
        parallel: false,
        ..MbarOptions::default()
    });
    let result = estimator
        .estimate_with_uncertainty(&windows)
        .expect("fit MBAR");
    let delta_index = result.n_states() - 1;
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("epot")
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
}

#[test]
fn mbar_cli_supports_all_observable_for_u_nk_preprocessing() {
    let inputs = acetamide_inputs();
    let output = run_cli(
        &[
            "mbar",
            "--decorrelate",
            "--u-nk-observable",
            "all",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let windows = load_decorrelated_windows_all(&inputs);
    let estimator = MbarEstimator::new(MbarOptions {
        max_iterations: 10_000,
        tolerance: 1.0e-7,
        parallel: false,
        ..MbarOptions::default()
    });
    let result = estimator
        .estimate_with_uncertainty(&windows)
        .expect("fit MBAR");
    let delta_index = result.n_states() - 1;
    let expected_counts = expected_u_nk_counts(&inputs, &windows);

    assert_close(
        payload["delta_f"].as_f64().expect("delta_f"),
        result.values()[delta_index],
    );
    assert_close(
        payload["uncertainty"].as_f64().expect("uncertainty"),
        result.uncertainties().expect("uncertainties")[delta_index],
    );
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("all")
    );
    assert_eq!(
        payload["provenance"]["samples_kept"].as_u64(),
        Some(expected_counts.samples_kept as u64)
    );
}

#[test]
fn mbar_cli_reports_nonfinite_de_observable_error() {
    let input_path = write_nonfinite_u_nk_input();
    let output = run_cli_failure(
        &[
            "mbar",
            "--decorrelate",
            "--u-nk-observable",
            "de",
            "--output-format",
            "json",
        ],
        std::slice::from_ref(&input_path),
    );

    assert!(!output.status.success(), "expected command to fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("decorrelation is unsupported for +inf reduced energies"),
        "unexpected stderr:\n{stderr}"
    );

    fs::remove_file(&input_path).expect("remove temporary AMBER input");
}

#[test]
fn mbar_cli_outputs_multidimensional_gromacs_json_with_lambda_components() {
    let inputs = gromacs_lambda15_input();
    let output = run_cli(
        &[
            "mbar",
            "--temperature",
            "298",
            "--u-nk-observable",
            "all",
            "--output-format",
            "json",
        ],
        &inputs,
    );
    let payload = parse_json_output(&output);

    let from_lambda = payload["from_lambda"]
        .as_array()
        .expect("from_lambda array");
    let to_lambda = payload["to_lambda"].as_array().expect("to_lambda array");
    assert_eq!(from_lambda.len(), 5);
    assert_eq!(to_lambda.len(), 5);
    assert_close(from_lambda[0].as_f64().expect("from[0]"), 0.0);
    assert_close(from_lambda[1].as_f64().expect("from[1]"), 0.0);
    assert_close(from_lambda[2].as_f64().expect("from[2]"), 0.7);
    assert_close(from_lambda[3].as_f64().expect("from[3]"), 0.0);
    assert_close(from_lambda[4].as_f64().expect("from[4]"), 0.0);
    assert_close(to_lambda[0].as_f64().expect("to[0]"), 0.0);
    assert_close(to_lambda[1].as_f64().expect("to[1]"), 0.0);
    assert_close(to_lambda[2].as_f64().expect("to[2]"), 0.9);
    assert_close(to_lambda[3].as_f64().expect("to[3]"), 0.0);
    assert_close(to_lambda[4].as_f64().expect("to[4]"), 0.0);
    assert_eq!(
        payload["provenance"]["lambda_components"]
            .as_array()
            .expect("lambda_components"),
        &vec![
            Value::from("mass-lambda"),
            Value::from("coul-lambda"),
            Value::from("vdw-lambda"),
            Value::from("bonded-lambda"),
            Value::from("restraint-lambda"),
        ]
    );
    assert_eq!(
        payload["provenance"]["u_nk_observable"].as_str(),
        Some("all")
    );
    assert_eq!(payload["provenance"]["windows"].as_u64(), Some(1));
    assert_eq!(payload["provenance"]["samples_kept"].as_u64(), Some(200));
}

#[test]
fn mbar_cli_outputs_multidimensional_gromacs_text_with_lambda_components() {
    let inputs = gromacs_lambda15_input();
    let output = run_cli(
        &[
            "mbar",
            "--temperature",
            "298",
            "--u-nk-observable",
            "all",
            "--output-format",
            "text",
        ],
        &inputs,
    );
    let stdout = String::from_utf8(output.stdout).expect("utf8 stdout");

    assert!(stdout.contains("from_lambda: [0, 0, 0.7, 0, 0]"));
    assert!(stdout.contains("to_lambda: [0, 0, 0.9, 0, 0]"));
    assert!(stdout.contains(
        "lambda_components: mass-lambda, coul-lambda, vdw-lambda, bonded-lambda, restraint-lambda"
    ));
    assert!(stdout.contains("u_nk_observable: all"));
}

#[test]
fn mbar_cli_outputs_multidimensional_gromacs_csv_with_lambda_components() {
    let inputs = gromacs_lambda15_input();
    let output = run_cli(
        &[
            "mbar",
            "--temperature",
            "298",
            "--u-nk-observable",
            "epot",
            "--output-format",
            "csv",
        ],
        &inputs,
    );
    let stdout = String::from_utf8(output.stdout).expect("utf8 stdout");

    assert!(stdout.contains("lambda_components"));
    assert!(stdout.contains("[0;0;0.7;0;0]"));
    assert!(stdout.contains("[0;0;0.9;0;0]"));
    assert!(stdout.contains("[mass-lambda;coul-lambda;vdw-lambda;bonded-lambda;restraint-lambda]"));
    assert!(stdout.contains(",epot,"));
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn acetamide_inputs() -> Vec<PathBuf> {
    let base = repo_root()
        .join("fixtures")
        .join("amber")
        .join("acetamide_tiny");
    [
        "0.0", "0.05", "0.1", "0.15", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.85",
        "0.90", "0.95", "1.0",
    ]
    .into_iter()
    .map(|lambda| base.join(lambda).join("acetamide.prod.out"))
    .collect()
}

fn gromacs_lambda15_input() -> Vec<PathBuf> {
    vec![repo_root()
        .join("fixtures")
        .join("gromacs")
        .join("lambda_15.xvg")]
}

fn run_cli(args: &[&str], inputs: &[PathBuf]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_alchemrs"));
    command.args(args);
    for input in inputs {
        command.arg(input);
    }
    let output = command.output().expect("run alchemrs");
    assert!(
        output.status.success(),
        "command failed with stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_cli_failure(args: &[&str], inputs: &[PathBuf]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_alchemrs"));
    command.args(args);
    for input in inputs {
        command.arg(input);
    }
    command.output().expect("run alchemrs")
}

fn parse_json_output(output: &Output) -> Value {
    serde_json::from_slice(&output.stdout).expect("parse CLI JSON output")
}

fn json_state(state: &[f64]) -> Value {
    if state.len() == 1 {
        Value::from(state[0])
    } else {
        Value::Array(state.iter().copied().map(Value::from).collect())
    }
}

fn load_decorrelated_windows(inputs: &[PathBuf]) -> Vec<UNkMatrix> {
    inputs
        .iter()
        .map(|path| load_decorrelated_window(path))
        .collect()
}

fn load_decorrelated_windows_nonconservative(inputs: &[PathBuf]) -> Vec<UNkMatrix> {
    inputs
        .iter()
        .map(|path| load_decorrelated_window_nonconservative(path))
        .collect()
}

fn load_decorrelated_windows_epot(inputs: &[PathBuf]) -> Vec<UNkMatrix> {
    inputs
        .iter()
        .map(|path| load_decorrelated_window_epot(path))
        .collect()
}

fn load_decorrelated_windows_all(inputs: &[PathBuf]) -> Vec<UNkMatrix> {
    inputs
        .iter()
        .map(|path| load_decorrelated_window_all(path))
        .collect()
}

fn load_auto_equilibrated_windows(inputs: &[PathBuf]) -> Vec<UNkMatrix> {
    inputs
        .iter()
        .map(|path| load_auto_equilibrated_window(path))
        .collect()
}

fn load_auto_equilibrated_decorrelated_windows(inputs: &[PathBuf]) -> Vec<UNkMatrix> {
    inputs
        .iter()
        .map(|path| load_auto_equilibrated_decorrelated_window(path))
        .collect()
}

fn load_decorrelated_dhdl_series(inputs: &[PathBuf]) -> Vec<DhdlSeries> {
    inputs
        .iter()
        .map(|path| {
            let series = extract_dhdl(path, TEMPERATURE_K).expect("parse dhdl");
            decorrelate_dhdl(&series, &DecorrelationOptions::default()).expect("decorrelate dhdl")
        })
        .collect()
}

fn load_auto_equilibrated_dhdl_series(inputs: &[PathBuf]) -> Vec<DhdlSeries> {
    inputs
        .iter()
        .map(|path| {
            let series = extract_dhdl(path, TEMPERATURE_K).expect("parse dhdl");
            let equilibration = detect_equilibration_dhdl(&series, &auto_equilibrate_options())
                .expect("detect equilibration");
            trim_dhdl_series(series, equilibration.t0)
        })
        .collect()
}

fn load_decorrelated_window(path: &Path) -> UNkMatrix {
    let u_nk = extract_u_nk(path, TEMPERATURE_K).expect("parse u_nk");
    decorrelate_u_nk(&u_nk, UNkSeriesMethod::DE, &DecorrelationOptions::default())
        .expect("decorrelate u_nk")
}

fn load_decorrelated_window_nonconservative(path: &Path) -> UNkMatrix {
    let u_nk = extract_u_nk(path, TEMPERATURE_K).expect("parse u_nk");
    decorrelate_u_nk(
        &u_nk,
        UNkSeriesMethod::DE,
        &nonconservative_decorrelation_options(),
    )
    .expect("decorrelate u_nk")
}

fn load_decorrelated_window_epot(path: &Path) -> UNkMatrix {
    let (u_nk, potential) =
        extract_u_nk_with_potential(path, TEMPERATURE_K).expect("parse u_nk with potential");
    decorrelate_u_nk_with_observable(&u_nk, &potential, &DecorrelationOptions::default())
        .expect("decorrelate u_nk with EPtot")
}

fn load_decorrelated_window_all(path: &Path) -> UNkMatrix {
    let u_nk = extract_u_nk(path, TEMPERATURE_K).expect("parse u_nk");
    decorrelate_u_nk(
        &u_nk,
        UNkSeriesMethod::All,
        &DecorrelationOptions::default(),
    )
    .expect("decorrelate u_nk with all")
}

fn load_auto_equilibrated_window(path: &Path) -> UNkMatrix {
    let u_nk = extract_u_nk(path, TEMPERATURE_K).expect("parse u_nk");
    let equilibration =
        detect_equilibration_u_nk(&u_nk, UNkSeriesMethod::DE, &auto_equilibrate_options())
            .expect("detect equilibration on dE");
    trim_u_nk_matrix(u_nk, equilibration.t0)
}

fn load_auto_equilibrated_decorrelated_window(path: &Path) -> UNkMatrix {
    let u_nk = extract_u_nk(path, TEMPERATURE_K).expect("parse u_nk");
    let equilibration =
        detect_equilibration_u_nk(&u_nk, UNkSeriesMethod::DE, &auto_equilibrate_options())
            .expect("detect equilibration on dE");
    let u_nk = trim_u_nk_matrix(u_nk, equilibration.t0);
    decorrelate_u_nk(&u_nk, UNkSeriesMethod::DE, &auto_equilibrate_options())
        .expect("decorrelate u_nk after auto equilibration")
}

fn expected_u_nk_counts(
    inputs: &[PathBuf],
    windows: &[UNkMatrix],
) -> alchemrs_cli_input_like::AnalysisSampleCountsLike {
    let samples_in = inputs
        .iter()
        .map(|path| {
            extract_u_nk(path, TEMPERATURE_K)
                .expect("parse raw u_nk")
                .n_samples()
        })
        .sum();
    let samples_kept = windows.iter().map(|window| window.n_samples()).sum();
    alchemrs_cli_input_like::AnalysisSampleCountsLike {
        windows: windows.len(),
        samples_in,
        samples_after_burnin: samples_in,
        samples_kept,
    }
}

fn expected_u_nk_counts_after_auto_equilibration(
    inputs: &[PathBuf],
    windows: &[UNkMatrix],
) -> alchemrs_cli_input_like::AnalysisSampleCountsLike {
    let samples_in = inputs
        .iter()
        .map(|path| {
            extract_u_nk(path, TEMPERATURE_K)
                .expect("parse raw u_nk")
                .n_samples()
        })
        .sum();
    let samples_after_burnin = windows.iter().map(|window| window.n_samples()).sum();
    alchemrs_cli_input_like::AnalysisSampleCountsLike {
        windows: windows.len(),
        samples_in,
        samples_after_burnin,
        samples_kept: samples_after_burnin,
    }
}

fn expected_u_nk_counts_after_auto_equilibration_and_decorrelation(
    inputs: &[PathBuf],
    auto_windows: &[UNkMatrix],
    windows: &[UNkMatrix],
) -> alchemrs_cli_input_like::AnalysisSampleCountsLike {
    let samples_in = inputs
        .iter()
        .map(|path| {
            extract_u_nk(path, TEMPERATURE_K)
                .expect("parse raw u_nk")
                .n_samples()
        })
        .sum();
    let samples_after_burnin = auto_windows.iter().map(|window| window.n_samples()).sum();
    let samples_kept = windows.iter().map(|window| window.n_samples()).sum();
    alchemrs_cli_input_like::AnalysisSampleCountsLike {
        windows: windows.len(),
        samples_in,
        samples_after_burnin,
        samples_kept,
    }
}

fn expected_dhdl_counts(
    inputs: &[PathBuf],
    series: &[DhdlSeries],
) -> alchemrs_cli_input_like::AnalysisSampleCountsLike {
    let samples_in = inputs
        .iter()
        .map(|path| {
            extract_dhdl(path, TEMPERATURE_K)
                .expect("parse raw dhdl")
                .values()
                .len()
        })
        .sum();
    let samples_kept = series.iter().map(|item| item.values().len()).sum();
    alchemrs_cli_input_like::AnalysisSampleCountsLike {
        windows: series.len(),
        samples_in,
        samples_after_burnin: samples_in,
        samples_kept,
    }
}

fn expected_dhdl_counts_after_auto_equilibration(
    inputs: &[PathBuf],
    series: &[DhdlSeries],
) -> alchemrs_cli_input_like::AnalysisSampleCountsLike {
    let samples_in = inputs
        .iter()
        .map(|path| {
            extract_dhdl(path, TEMPERATURE_K)
                .expect("parse raw dhdl")
                .values()
                .len()
        })
        .sum();
    let samples_after_burnin = series.iter().map(|item| item.values().len()).sum();
    alchemrs_cli_input_like::AnalysisSampleCountsLike {
        windows: series.len(),
        samples_in,
        samples_after_burnin,
        samples_kept: samples_after_burnin,
    }
}

fn auto_equilibrate_options() -> DecorrelationOptions {
    DecorrelationOptions {
        fast: true,
        conservative: false,
        ..DecorrelationOptions::default()
    }
}

fn nonconservative_decorrelation_options() -> DecorrelationOptions {
    DecorrelationOptions {
        conservative: false,
        ..DecorrelationOptions::default()
    }
}

fn trim_u_nk_matrix(u_nk: UNkMatrix, remove_burnin: usize) -> UNkMatrix {
    if remove_burnin == 0 {
        return u_nk;
    }

    let n_samples = u_nk.n_samples();
    let n_states = u_nk.n_states();
    let data = u_nk.data()[(remove_burnin * n_states)..].to_vec();
    let time = u_nk.time_ps()[remove_burnin..].to_vec();
    UNkMatrix::new(
        n_samples - remove_burnin,
        n_states,
        data,
        time,
        u_nk.sampled_state().cloned(),
        u_nk.evaluated_states().to_vec(),
    )
    .expect("trim u_nk")
}

fn trim_dhdl_series(dhdl: DhdlSeries, remove_burnin: usize) -> DhdlSeries {
    if remove_burnin == 0 {
        return dhdl;
    }

    DhdlSeries::new(
        dhdl.state().clone(),
        dhdl.time_ps()[remove_burnin..].to_vec(),
        dhdl.values()[remove_burnin..].to_vec(),
    )
    .expect("trim dhdl")
}

fn compute_overlap_summary(windows: &[UNkMatrix]) -> (f64, Vec<f64>) {
    let overlap = overlap_matrix(
        windows,
        Some(MbarOptions {
            parallel: false,
            ..MbarOptions::default()
        }),
    )
    .expect("compute overlap matrix");
    let eigenvalues = overlap_eigenvalues(&overlap).expect("compute overlap eigenvalues");
    (1.0 - eigenvalues[1], eigenvalues)
}

fn unique_output_path(file_name: &str) -> PathBuf {
    let suffix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("alchemrs-{suffix}-{file_name}"))
}

fn write_nonfinite_u_nk_input() -> PathBuf {
    let path = unique_output_path("nonfinite-u-nk.out");
    let content = r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.0000
FEP MBAR options:
    bar_intervall = 10
    MBAR - lambda values considered:
      total   0.0000 0.1000
    Extra
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
MBAR Energy analysis:
Energy at 0.0000 =    -10.0
Energy at 0.1000 = **********
 ---
MBAR Energy analysis:
Energy at 0.0000 =    -11.0
Energy at 0.1000 =     -8.0
 ---
"#;
    fs::write(&path, content).expect("write temporary AMBER input");
    path
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-9;
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected}, got {actual}"
    );
}

fn assert_json_optional_f64_close(actual: &Value, expected: Option<f64>) {
    match expected {
        Some(expected) => assert_close(actual.as_f64().expect("optional float"), expected),
        None => assert!(actual.is_null(), "expected null, got {actual:?}"),
    }
}

mod alchemrs_cli_input_like {
    pub struct AnalysisSampleCountsLike {
        pub windows: usize,
        pub samples_in: usize,
        pub samples_after_burnin: usize,
        pub samples_kept: usize,
    }
}
