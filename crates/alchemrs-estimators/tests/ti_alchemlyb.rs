use alchemrs_estimators::TiEstimator;
use alchemrs_parse::amber::extract_dhdl;
use std::fs;
use std::path::PathBuf;

fn read_expected(path: &str) -> (f64, f64) {
    let content = fs::read_to_string(path).expect("read expected delta_f");
    let mut parts = content.trim().split(',');
    let delta = parts
        .next()
        .expect("delta")
        .trim()
        .parse::<f64>()
        .expect("parse delta");
    let sigma = parts
        .next()
        .expect("sigma")
        .trim()
        .parse::<f64>()
        .expect("parse sigma");
    (delta, sigma)
}

#[test]
fn ti_matches_alchemlyb_all_windows() {
    let base = env!("CARGO_MANIFEST_DIR");
    let mut paths: Vec<PathBuf> =
        fs::read_dir(format!("{base}/../../fixtures/amber/acetamide_tiny"))
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

    let mut series = Vec::new();
    for path in paths {
        series.push(extract_dhdl(path, 300.0).expect("parse AMBER output"));
    }

    let estimator = TiEstimator::default();
    let result = estimator.fit(&series).expect("TI fit");

    let expected_path =
        format!("{base}/../../fixtures/amber/acetamide_tiny/ti_0.0_1.0.delta_f_sigma.txt");
    let (expected_delta, expected_sigma) = read_expected(&expected_path);
    assert!(
        (result.delta_f() - expected_delta).abs() < 1e-6,
        "delta_f mismatch: {} vs {}",
        result.delta_f(),
        expected_delta
    );
    let sigma = result.uncertainty().expect("uncertainty missing");
    assert!(
        (sigma - expected_sigma).abs() < 1e-6,
        "uncertainty mismatch: {} vs {}",
        sigma,
        expected_sigma
    );
}
