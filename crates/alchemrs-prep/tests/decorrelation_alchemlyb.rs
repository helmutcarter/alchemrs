use std::fs;

use alchemrs_parse::amber::extract_dhdl;
use alchemrs_prep::{decorrelate_dhdl, DecorrelationOptions};

fn load_expected(path: &str) -> Vec<(f64, f64)> {
    let content = fs::read_to_string(path).expect("read expected decorrelation file");
    let mut out = Vec::new();
    let mut header_cols = None;
    for (idx, line) in content.lines().enumerate() {
        if idx == 0 {
            header_cols = Some(line.split(',').count());
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        let cols = header_cols.unwrap_or(parts.len());
        let time = parts.first()
            .expect("time")
            .parse::<f64>()
            .expect("time parse");
        let dhdl_index = if cols >= 3 { 2 } else { 1 };
        let dhdl = parts
            .get(dhdl_index)
            .expect("dhdl")
            .parse::<f64>()
            .expect("dhdl parse");
        out.push((time, dhdl));
    }
    out
}

#[test]
fn decorrelate_dhdl_matches_alchemlyb() {
    let base = env!("CARGO_MANIFEST_DIR");
    let input = format!("{base}/../../fixtures/amber/acetamide_tiny/0.1/acetamide.prod.out");
    let expected_path =
        format!("{base}/../../fixtures/amber/acetamide_tiny/0.1/dhdl.decorrelated.csv");

    let series = extract_dhdl(input, 300.0).expect("parse AMBER output");
    let result = decorrelate_dhdl(&series, &DecorrelationOptions::default())
        .expect("decorrelate dhdl");

    let expected = load_expected(&expected_path);
    assert_eq!(result.values().len(), expected.len());
    for (idx, (time, dhdl)) in expected.iter().enumerate() {
        let rt = result.time_ps()[idx];
        let rv = result.values()[idx];
        assert!(
            (rt - time).abs() < 1e-6,
            "time mismatch at {idx}: {rt} vs {time}"
        );
        assert!(
            (rv - dhdl).abs() < 1e-6,
            "dhdl mismatch at {idx}: {rv} vs {dhdl}"
        );
    }
}
