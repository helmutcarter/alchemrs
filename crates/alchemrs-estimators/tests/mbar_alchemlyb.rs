use std::fs;
use std::path::PathBuf;

use alchemrs_estimators::MbarEstimator;
use alchemrs_parse::amber::extract_u_nk;

fn read_expected(path: &str) -> (f64, f64) {
    let content = fs::read_to_string(path).expect("read expected mbar output");
    let mut parts = content.trim().split(',');
    let delta = parts
        .next()
        .expect("delta")
        .trim()
        .parse::<f64>()
        .expect("parse expected delta_f");
    let sigma = parts
        .next()
        .expect("sigma")
        .trim()
        .parse::<f64>()
        .expect("parse expected sigma");
    (delta, sigma)
}

#[test]
fn mbar_matches_alchemlyb_all_windows() {
    let base = env!("CARGO_MANIFEST_DIR");
    let mut paths: Vec<PathBuf> = fs::read_dir(format!("{base}/../../fixtures/amber/acetamide_tiny"))
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
        let la = a.parent().and_then(|p| p.file_name()).and_then(|s| s.to_str()).unwrap();
        let lb = b.parent().and_then(|p| p.file_name()).and_then(|s| s.to_str()).unwrap();
        la.parse::<f64>().unwrap().partial_cmp(&lb.parse::<f64>().unwrap()).unwrap()
    });

    let mut windows = Vec::new();
    for path in paths {
        windows.push(extract_u_nk(path, 300.0).expect("parse AMBER output"));
    }

    let estimator = MbarEstimator::default();
    let result = estimator.fit(&windows).expect("MBAR fit");

    let expected_path = format!(
        "{base}/../../fixtures/amber/acetamide_tiny/mbar_0.0_1.0.delta_f_sigma.txt"
    );
    let (expected_delta, expected_sigma) = read_expected(&expected_path);

    let n = result.n_states();
    let delta = result.values()[0 * n + (n - 1)];
    assert!(
        (delta - expected_delta).abs() < 1e-5,
        "delta_f mismatch: {} vs {}",
        delta,
        expected_delta
    );

    let sigma = result
        .uncertainties()
        .expect("uncertainties missing")[0 * n + (n - 1)];
    assert!(
        (sigma - expected_sigma).abs() < 1e-5,
        "sigma mismatch: {} vs {}",
        sigma,
        expected_sigma
    );
}
