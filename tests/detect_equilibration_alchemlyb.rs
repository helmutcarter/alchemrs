use std::fs;

use alchemrs::{detect_equilibration_dhdl, extract_dhdl, DecorrelationOptions};

fn read_expected(path: &str) -> (usize, f64, f64) {
    let content = fs::read_to_string(path).expect("read expected equilibration");
    let parts: Vec<&str> = content.trim().split(',').collect();
    let t0 = parts[0].trim().parse::<usize>().expect("parse t0");
    let g = parts[1].trim().parse::<f64>().expect("parse g");
    let neff = parts[2].trim().parse::<f64>().expect("parse neff");
    (t0, g, neff)
}

#[test]
fn detect_equilibration_matches_repository_reference() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}/fixtures/amber/acetamide_tiny/0.1/acetamide.prod.out");
    let series = extract_dhdl(path, 300.0).expect("parse AMBER output");
    let expected_path =
        format!("{base}/fixtures/amber/acetamide_tiny/detect_equilibration_0.1.txt");
    let (t0, g, neff) = read_expected(&expected_path);

    let options = DecorrelationOptions::default();
    let result = detect_equilibration_dhdl(&series, &options).expect("detect equilibration");

    assert_eq!(result.t0, t0);
    assert!(
        (result.g - g).abs() < 1e-6,
        "g mismatch: {} vs {}",
        result.g,
        g
    );
    assert!(
        (result.neff_max - neff).abs() < 1e-6,
        "neff mismatch: {} vs {}",
        result.neff_max,
        neff
    );
}
