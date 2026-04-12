use std::fs;
use std::path::PathBuf;

use alchemrs::{extract_u_nk, UwhamEstimator};

fn amber_fixture_paths() -> Vec<PathBuf> {
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
    paths
}

fn main() {
    let output_dir = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("target")
                .join("uwham-reference")
                .join("acetamide_tiny")
        });

    let mut windows = Vec::new();
    for path in amber_fixture_paths() {
        windows.push(extract_u_nk(path, 300.0).expect("parse AMBER output"));
    }

    let fit = UwhamEstimator::default().fit(&windows).expect("UWHAM fit");
    fit.write_reference_inputs(&output_dir)
        .expect("write UWHAM reference inputs");

    println!("wrote UWHAM reference inputs to {}", output_dir.display());
}
