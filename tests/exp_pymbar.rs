use std::fs;
use std::path::PathBuf;

use alchemrs::{extract_u_nk, IexpEstimator};

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
fn exp_matches_pymbar_adjacent_fep_chain() {
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

    let estimator = IexpEstimator::default();
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

    let n_states = result.n_states();
    let mut expected_forward_endpoint = 0.0;
    let mut expected_reverse_endpoint = 0.0;
    let mut expected_forward_var = 0.0;
    let mut expected_reverse_var = 0.0;

    for state_idx in 0..(n_states - 1) {
        let upper = state_idx * n_states + state_idx + 1;
        let lower = (state_idx + 1) * n_states + state_idx;

        assert!(
            (values[upper] - expected_delta[upper]).abs() < 1e-6,
            "forward adjacent delta_f mismatch at {upper}: {} vs {}",
            values[upper],
            expected_delta[upper]
        );
        assert!(
            (values[lower] - expected_delta[lower]).abs() < 1e-6,
            "reverse adjacent delta_f mismatch at {lower}: {} vs {}",
            values[lower],
            expected_delta[lower]
        );
        assert!(
            (sigma[upper] - expected_sigma[upper]).abs() < 1e-6,
            "forward adjacent sigma mismatch at {upper}: {} vs {}",
            sigma[upper],
            expected_sigma[upper]
        );
        assert!(
            (sigma[lower] - expected_sigma[lower]).abs() < 1e-6,
            "reverse adjacent sigma mismatch at {lower}: {} vs {}",
            sigma[lower],
            expected_sigma[lower]
        );

        expected_forward_endpoint += expected_delta[upper];
        expected_reverse_endpoint += expected_delta[lower];
        expected_forward_var += expected_sigma[upper] * expected_sigma[upper];
        expected_reverse_var += expected_sigma[lower] * expected_sigma[lower];
    }

    assert!(
        (values[n_states - 1] - expected_forward_endpoint).abs() < 1e-6,
        "forward endpoint delta_f mismatch: {} vs {}",
        values[n_states - 1],
        expected_forward_endpoint
    );
    assert!(
        (values[(n_states - 1) * n_states] - expected_reverse_endpoint).abs() < 1e-6,
        "reverse endpoint delta_f mismatch: {} vs {}",
        values[(n_states - 1) * n_states],
        expected_reverse_endpoint
    );
    assert!(
        (sigma[n_states - 1] - expected_forward_var.sqrt()).abs() < 1e-6,
        "forward endpoint sigma mismatch: {} vs {}",
        sigma[n_states - 1],
        expected_forward_var.sqrt()
    );
    assert!(
        (sigma[(n_states - 1) * n_states] - expected_reverse_var.sqrt()).abs() < 1e-6,
        "reverse endpoint sigma mismatch: {} vs {}",
        sigma[(n_states - 1) * n_states],
        expected_reverse_var.sqrt()
    );
}
