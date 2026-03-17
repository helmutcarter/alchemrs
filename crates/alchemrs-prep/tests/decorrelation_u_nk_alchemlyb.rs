use std::fs;

use alchemrs_parse::amber::extract_u_nk;
use alchemrs_prep::{decorrelate_u_nk, DecorrelationOptions, UNkSeriesMethod};

fn load_expected(path: &str) -> (Vec<f64>, Vec<f64>) {
    let content = fs::read_to_string(path).expect("read expected u_nk file");
    let mut lines = content.lines();
    let header = lines.next().expect("header");
    let cols: Vec<&str> = header.split(',').collect();
    if cols.len() < 3 {
        panic!("expected at least time, lambdas, and one state column");
    }
    let n_states = cols.len() - 2;

    let mut times = Vec::new();
    let mut data = Vec::new();
    for line in lines {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() != cols.len() {
            panic!("unexpected column count in row");
        }
        let time = parts[0].parse::<f64>().expect("time parse");
        times.push(time);
        for value in &parts[2..] {
            data.push(value.parse::<f64>().expect("value parse"));
        }
    }

    let expected_len = times.len() * n_states;
    assert_eq!(data.len(), expected_len);
    (times, data)
}

#[test]
fn decorrelate_u_nk_matches_alchemlyb() {
    let base = env!("CARGO_MANIFEST_DIR");
    let input = format!("{base}/../../fixtures/amber/acetamide_tiny/0.1/acetamide.prod.out");
    let expected_path =
        format!("{base}/../../fixtures/amber/acetamide_tiny/0.1/u_nk.decorrelated.csv");

    let u_nk = extract_u_nk(input, 300.0).expect("parse AMBER output");
    let result = decorrelate_u_nk(&u_nk, UNkSeriesMethod::DE, &DecorrelationOptions::default())
        .expect("decorrelate u_nk");

    let (times, data) = load_expected(&expected_path);
    assert_eq!(result.time_ps().len(), times.len());
    assert_eq!(result.data().len(), data.len());

    for (idx, time) in times.iter().enumerate() {
        let rt = result.time_ps()[idx];
        assert!(
            (rt - time).abs() < 1e-6,
            "time mismatch at {idx}: {rt} vs {time}"
        );
    }

    for (idx, value) in data.iter().enumerate() {
        let rv = result.data()[idx];
        assert!(
            (rv - value).abs() < 1e-6,
            "u_nk mismatch at {idx}: {rv} vs {value}"
        );
    }
}
