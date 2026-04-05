use std::fs;
use std::path::PathBuf;

use alchemrs::{extract_u_nk, ExpEstimator};

fn read_matrix(path: &str) -> Vec<f64> {
    let content = fs::read_to_string(path).expect("read expected matrix");
    let mut values = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        for entry in line.split(',') {
            let val = entry.trim().parse::<f64>().expect("parse value");
            values.push(val);
        }
    }
    values
}

#[test]
fn exp_matches_pymbar_matrix() {
    let base = env!("CARGO_MANIFEST_DIR");
    let mut paths: Vec<PathBuf> = fs::read_dir(format!("{base}/fixtures/amber/acetamide_tiny"))
        .expect("read fixture directory")
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let path = entry.path();
            if path.is_dir() {
                Some(path.join("acetamide.prod.out"))
            } else {
                None
            }
        })
        .collect();
    paths.sort_by(|a, b| {
        let la = a
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|s| s.to_str())
            .unwrap();
        let lb = b
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|s| s.to_str())
            .unwrap();
        la.parse::<f64>()
            .unwrap()
            .partial_cmp(&lb.parse::<f64>().unwrap())
            .unwrap()
    });

    let mut windows = Vec::new();
    for path in paths {
        windows.push(extract_u_nk(path, 300.0).expect("parse AMBER output"));
    }

    let estimator = ExpEstimator::default();
    let fit = estimator.fit(&windows).expect("EXP fit");
    let result = fit
        .result_with_uncertainty()
        .expect("EXP result with uncertainty");

    let expected_delta_path = format!("{base}/fixtures/amber/acetamide_tiny/exp_matrix.txt");
    let expected_sigma_path = format!("{base}/fixtures/amber/acetamide_tiny/exp_sigma.txt");
    let expected_delta = read_matrix(&expected_delta_path);
    let expected_sigma = read_matrix(&expected_sigma_path);

    let values = result.values();
    let sigma = result.uncertainties().expect("uncertainties missing");
    assert_eq!(values.len(), expected_delta.len());
    assert_eq!(sigma.len(), expected_sigma.len());

    for (idx, (value, exp)) in values.iter().zip(expected_delta.iter()).enumerate() {
        assert!(
            (value - exp).abs() < 1e-6,
            "delta_f mismatch at {idx}: {value} vs {exp}"
        );
    }

    for (idx, (value, exp)) in sigma.iter().zip(expected_sigma.iter()).enumerate() {
        assert!(
            (value - exp).abs() < 1e-6,
            "sigma mismatch at {idx}: {value} vs {exp}"
        );
    }
}
