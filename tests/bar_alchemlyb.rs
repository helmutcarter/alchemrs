use std::fs;
use std::path::PathBuf;

use alchemrs::{extract_u_nk, BarEstimator};

#[test]
fn bar_matches_alchemlyb_all_windows() {
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

    let estimator = BarEstimator::default();
    let fit = estimator.fit(&windows).expect("BAR fit");
    let result = fit.result().expect("BAR result");

    let expected_delta = 37.39789124;
    let expected_sigma = 0.44660501995747376;

    let n = result.n_states();
    let delta_index = n - 1;
    let delta = result.values()[delta_index];
    assert!(
        (delta - expected_delta).abs() < 1e-6,
        "delta_f mismatch: {} vs {}",
        delta,
        expected_delta
    );

    let sigma = result.uncertainties().expect("uncertainties missing")[delta_index];
    assert!(
        (sigma - expected_sigma).abs() < 1e-6,
        "uncertainty mismatch: {} vs {}",
        sigma,
        expected_sigma
    );
}
