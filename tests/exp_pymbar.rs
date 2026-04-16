use std::fs;
use std::path::PathBuf;

use alchemrs::{extract_u_nk, IexpEstimator, StatePoint, UNkMatrix};

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

fn state_index(window: &UNkMatrix, state: &StatePoint) -> usize {
    window
        .evaluated_states()
        .iter()
        .position(|candidate| candidate == state)
        .expect("state present in evaluated_states")
}

fn work_values(window: &UNkMatrix, from_idx: usize, to_idx: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(window.n_samples());
    for sample_idx in 0..window.n_samples() {
        let offset = sample_idx * window.n_states();
        out.push(window.data()[offset + to_idx] - window.data()[offset + from_idx]);
    }
    out
}

fn exp_sigma_from_work(work: &[f64]) -> f64 {
    let n = work.len() as f64;
    let max_arg = work
        .iter()
        .map(|value| -value)
        .fold(f64::NEG_INFINITY, f64::max);
    let x = work
        .iter()
        .map(|value| (-value - max_arg).exp())
        .collect::<Vec<_>>();
    let mean = x.iter().sum::<f64>() / n;
    let variance = x
        .iter()
        .map(|value| {
            let diff = value - mean;
            diff * diff
        })
        .sum::<f64>()
        / (x.len() - 1) as f64;
    variance.sqrt() / (n.sqrt() * mean)
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
    let expected_delta = read_matrix(&expected_delta_path);

    let values = result.values();
    let sigma = result.uncertainties().expect("uncertainties missing");
    assert_eq!(values.len(), expected_delta.len());

    let n_states = result.n_states();
    let mut expected_forward_endpoint = 0.0;
    let mut expected_reverse_endpoint = 0.0;
    let mut expected_forward_var = 0.0;
    let mut expected_reverse_var = 0.0;

    for state_idx in 0..(n_states - 1) {
        let upper = state_idx * n_states + state_idx + 1;
        let lower = (state_idx + 1) * n_states + state_idx;
        let from_state = windows[state_idx].sampled_state().expect("sampled_state");
        let to_state = windows[state_idx + 1]
            .sampled_state()
            .expect("sampled_state");
        let forward_sigma = exp_sigma_from_work(&work_values(
            &windows[state_idx],
            state_index(&windows[state_idx], from_state),
            state_index(&windows[state_idx], to_state),
        ));
        let reverse_sigma = exp_sigma_from_work(&work_values(
            &windows[state_idx + 1],
            state_index(&windows[state_idx + 1], to_state),
            state_index(&windows[state_idx + 1], from_state),
        ));

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
            (sigma[upper] - forward_sigma).abs() < 1e-6,
            "forward adjacent sigma mismatch at {upper}: {} vs {}",
            sigma[upper],
            forward_sigma
        );
        assert!(
            (sigma[lower] - reverse_sigma).abs() < 1e-6,
            "reverse adjacent sigma mismatch at {lower}: {} vs {}",
            sigma[lower],
            reverse_sigma
        );

        expected_forward_endpoint += expected_delta[upper];
        expected_reverse_endpoint += expected_delta[lower];
        expected_forward_var += forward_sigma * forward_sigma;
        expected_reverse_var += reverse_sigma * reverse_sigma;
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
