use std::fs;
use std::path::PathBuf;

use alchemrs_analysis::{overlap_eigenvalues, overlap_matrix, overlap_scalar};
use alchemrs_parse::amber::extract_u_nk;

fn read_expected(path: &str) -> Vec<f64> {
    let content = fs::read_to_string(path).expect("read expected overlap matrix");
    let mut values = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        for entry in line.split(',') {
            let val = entry.trim().parse::<f64>().expect("parse overlap value");
            values.push(val);
        }
    }
    values
}

fn read_scalar(path: &str) -> f64 {
    let content = fs::read_to_string(path).expect("read expected overlap scalar");
    content.trim().parse::<f64>().expect("parse scalar")
}

#[test]
fn overlap_matches_alchemlyb() {
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

    let mut windows = Vec::new();
    for path in paths {
        windows.push(extract_u_nk(path, 300.0).expect("parse AMBER output"));
    }

    let overlap = overlap_matrix(&windows, None).expect("overlap matrix");
    let expected_path = format!("{base}/../../fixtures/amber/acetamide_tiny/overlap_matrix.txt");
    let expected = read_expected(&expected_path);

    let values = overlap.values();
    assert_eq!(values.len(), expected.len());
    for (idx, (value, exp)) in values.iter().zip(expected.iter()).enumerate() {
        assert!(
            (value - exp).abs() < 1e-6,
            "overlap mismatch at {idx}: {value} vs {exp}"
        );
    }

    let eigen_path = format!("{base}/../../fixtures/amber/acetamide_tiny/overlap_eigenvalues.txt");
    let scalar_path = format!("{base}/../../fixtures/amber/acetamide_tiny/overlap_scalar.txt");
    let expected_eigen = read_expected(&eigen_path);
    let expected_scalar = read_scalar(&scalar_path);

    let eigen = overlap_eigenvalues(&overlap).expect("eigenvalues");
    assert_eq!(eigen.len(), expected_eigen.len());
    for (idx, (value, exp)) in eigen.iter().zip(expected_eigen.iter()).enumerate() {
        assert!(
            (value - exp).abs() < 1e-6,
            "eigenvalue mismatch at {idx}: {value} vs {exp}"
        );
    }

    let scalar = overlap_scalar(&overlap).expect("scalar");
    assert!(
        (scalar - expected_scalar).abs() < 1e-6,
        "scalar mismatch: {scalar} vs {expected_scalar}"
    );
}
