use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

use alchemrs::{extract_u_nk, UNkMatrix, UwhamEstimator, UwhamFit};
use tempfile::tempdir;

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

fn reduced_fixture_paths() -> Vec<PathBuf> {
    amber_fixture_paths()
        .into_iter()
        .filter(|path| {
            path.parent()
                .and_then(|p| p.file_name())
                .and_then(|s| s.to_str())
                .is_some_and(|name| matches!(name, "0.2" | "0.3" | "0.4" | "0.5" | "0.6"))
        })
        .collect()
}

fn amber_fixture_windows() -> &'static [UNkMatrix] {
    static WINDOWS: OnceLock<Vec<UNkMatrix>> = OnceLock::new();
    WINDOWS.get_or_init(|| {
        reduced_fixture_paths()
            .into_iter()
            .map(|path| extract_u_nk(path, 300.0).expect("parse AMBER output"))
            .collect()
    })
}

fn amber_fixture_uwham_fit() -> &'static UwhamFit {
    static FIT: OnceLock<UwhamFit> = OnceLock::new();
    FIT.get_or_init(|| {
        UwhamEstimator::default()
            .fit(amber_fixture_windows())
            .expect("UWHAM fit")
    })
}

fn read_csv_vector(path: &str) -> Vec<f64> {
    fs::read_to_string(path)
        .expect("read csv vector")
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.parse::<f64>().expect("parse vector value"))
        .collect()
}

fn read_csv_matrix(path: &str) -> Vec<f64> {
    fs::read_to_string(path)
        .expect("read csv matrix")
        .lines()
        .filter(|line| !line.trim().is_empty())
        .flat_map(|line| {
            line.split(',')
                .map(|value| value.parse::<f64>().expect("parse matrix value"))
                .collect::<Vec<_>>()
        })
        .collect()
}

#[test]
fn uwham_produces_finite_result_for_amber_fixture() {
    let uwham = amber_fixture_uwham_fit().result().expect("UWHAM delta_f");

    assert!(uwham.values().iter().all(|value| value.is_finite()));
    for i in 0..uwham.n_states() {
        assert_eq!(uwham.values()[i * uwham.n_states() + i], 0.0);
    }
}

#[test]
fn uwham_weights_satisfy_sampled_state_check_on_fixture() {
    let uwham = amber_fixture_uwham_fit();
    assert!(uwham.free_energies().iter().all(|value| value.is_finite()));

    let n_observations = uwham.n_observations() as f64;
    for state_idx in 0..uwham.n_states() {
        let mean_weight = uwham
            .weights()
            .chunks(uwham.n_states())
            .map(|row| row[state_idx])
            .sum::<f64>()
            / n_observations;
        assert!(
            (mean_weight - 1.0).abs() < 1e-6,
            "sampled-state check failed at {state_idx}: {mean_weight}"
        );
    }
}

#[test]
fn uwham_matches_committed_r_reference_on_fixture() {
    let uwham = amber_fixture_uwham_fit();
    let ze = read_csv_vector("fixtures/uwham/ze.csv");
    let delta_f = read_csv_matrix("fixtures/uwham/delta_f.csv");

    assert_eq!(ze.len(), uwham.free_energies().len());
    let gauge_shift = uwham.free_energies()[uwham.n_states() - 1] - ze[ze.len() - 1];
    for (idx, (lhs, rhs)) in ze.iter().zip(uwham.free_energies().iter()).enumerate() {
        let shifted_rhs = rhs - gauge_shift;
        assert!(
            (lhs - shifted_rhs).abs() < 1e-5,
            "ze mismatch at {idx}: {lhs} vs {shifted_rhs}"
        );
    }

    let result = uwham.result().expect("UWHAM delta_f");
    assert_eq!(delta_f.len(), result.values().len());
    for (idx, (lhs, rhs)) in delta_f.iter().zip(result.values().iter()).enumerate() {
        assert!(
            (lhs - rhs).abs() < 1e-5,
            "delta_f mismatch at {idx}: {lhs} vs {rhs}"
        );
    }
}

#[test]
fn uwham_exports_reference_inputs() {
    let fit = amber_fixture_uwham_fit();
    let dir = tempdir().expect("tempdir");
    fit.write_reference_inputs(dir.path())
        .expect("write UWHAM reference inputs");

    for name in [
        "logQ.csv",
        "size.csv",
        "delta_f.csv",
        "ze.csv",
        "README.txt",
    ] {
        let path = dir.path().join(name);
        assert!(path.exists(), "missing exported file: {}", path.display());
        let content = fs::read_to_string(&path).expect("read exported file");
        assert!(!content.trim().is_empty(), "empty exported file: {name}");
    }
}
