use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use alchemrs_analysis::{overlap_eigenvalues, overlap_matrix};
use alchemrs_estimators::{BarEstimator, BarMethod, BarOptions, MbarEstimator, MbarOptions};
use alchemrs_parse::amber::extract_u_nk_with_potential;
use alchemrs_prep::{decorrelate_u_nk_with_observable, DecorrelationOptions};
use serde_json::Value;

const TEMPERATURE_K: f64 = 300.0;

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
        compute_uncertainty: true,
        parallel: false,
        ..MbarOptions::default()
    });
    let result = estimator.fit(&windows).expect("fit MBAR");
    let delta_index = result.n_states() - 1;
    let expected_overlap = compute_overlap_summary(&windows);

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
    let result = estimator.fit(&windows).expect("fit BAR");
    let delta_index = result.n_states() - 1;
    let expected_overlap = compute_overlap_summary(&windows);

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

    fs::remove_file(&output_path).expect("remove CLI output file");
}

fn workspace_root() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .expect("crates dir")
        .parent()
        .expect("workspace root")
        .to_path_buf()
}

fn acetamide_inputs() -> Vec<PathBuf> {
    let base = workspace_root()
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

fn run_cli(args: &[&str], inputs: &[PathBuf]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_alchemrs-cli"));
    command.args(args);
    for input in inputs {
        command.arg(input);
    }
    let output = command.output().expect("run alchemrs-cli");
    assert!(
        output.status.success(),
        "command failed with stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn parse_json_output(output: &Output) -> Value {
    serde_json::from_slice(&output.stdout).expect("parse CLI JSON output")
}

fn load_decorrelated_windows(inputs: &[PathBuf]) -> Vec<alchemrs_core::UNkMatrix> {
    inputs
        .iter()
        .map(|path| load_decorrelated_window(path))
        .collect()
}

fn load_decorrelated_window(path: &Path) -> alchemrs_core::UNkMatrix {
    let (u_nk, potential) =
        extract_u_nk_with_potential(path, TEMPERATURE_K).expect("parse u_nk with potential");
    decorrelate_u_nk_with_observable(&u_nk, &potential, &DecorrelationOptions::default())
        .expect("decorrelate u_nk")
}

fn compute_overlap_summary(windows: &[alchemrs_core::UNkMatrix]) -> (f64, Vec<f64>) {
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
    std::env::temp_dir().join(format!("alchemrs-cli-{suffix}-{file_name}"))
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-9;
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected}, got {actual}"
    );
}
