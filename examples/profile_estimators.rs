use std::env;
use std::fs;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::Instant;

use alchemrs::{
    decorrelate_u_nk_with_observable, extract_dhdl, extract_u_nk_with_potential, BarEstimator,
    BarOptions, DecorrelationOptions, DhdlSeries, IexpEstimator, IexpOptions, MbarEstimator,
    MbarOptions, TiEstimator, TiOptions, UNkMatrix, UwhamEstimator, UwhamOptions,
};

type DynError = Box<dyn std::error::Error>;

fn main() -> Result<(), DynError> {
    let mut args = env::args().skip(1);
    let workload = args.next().ok_or_else(usage)?;
    let iterations = args
        .next()
        .map(|value| value.parse::<usize>())
        .transpose()?
        .unwrap_or(1);
    if iterations == 0 {
        return Err("iterations must be positive".into());
    }

    match workload.as_str() {
        "amber-parse-u-nk" => measure(&workload, iterations, bench_amber_parse_u_nk)?,
        "amber-parse-dhdl" => measure(&workload, iterations, bench_amber_parse_dhdl)?,
        "amber-decorrelate-u-nk" => {
            let parsed = load_amber_u_nk_with_potential()?;
            measure(&workload, iterations, || bench_decorrelate_u_nk(&parsed))?;
        }
        "amber-mbar-serial" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = MbarEstimator::new(MbarOptions {
                parallel: false,
                ..MbarOptions::default()
            });
            measure(&workload, iterations, || bench_mbar(&estimator, &windows))?;
        }
        "amber-mbar-parallel" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = MbarEstimator::new(MbarOptions {
                parallel: true,
                ..MbarOptions::default()
            });
            measure(&workload, iterations, || bench_mbar(&estimator, &windows))?;
        }
        "amber-uwham-serial" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = UwhamEstimator::new(UwhamOptions {
                parallel: false,
                ..UwhamOptions::default()
            });
            measure(&workload, iterations, || bench_uwham(&estimator, &windows))?;
        }
        "amber-uwham-parallel" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = UwhamEstimator::new(UwhamOptions {
                parallel: true,
                ..UwhamOptions::default()
            });
            measure(&workload, iterations, || bench_uwham(&estimator, &windows))?;
        }
        "amber-bar-serial" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = BarEstimator::new(BarOptions {
                parallel: false,
                ..BarOptions::default()
            });
            measure(&workload, iterations, || bench_bar(&estimator, &windows))?;
        }
        "amber-bar-parallel" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = BarEstimator::new(BarOptions {
                parallel: true,
                ..BarOptions::default()
            });
            measure(&workload, iterations, || bench_bar(&estimator, &windows))?;
        }
        "amber-exp-serial" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = IexpEstimator::new(IexpOptions { parallel: false });
            measure(&workload, iterations, || bench_exp(&estimator, &windows))?;
        }
        "amber-exp-parallel" => {
            let windows = amber_windows_only(load_amber_u_nk_with_potential()?);
            let estimator = IexpEstimator::new(IexpOptions { parallel: true });
            measure(&workload, iterations, || bench_exp(&estimator, &windows))?;
        }
        "amber-ti-serial" => {
            let series = load_amber_dhdl_series()?;
            let estimator = TiEstimator::new(TiOptions {
                parallel: false,
                ..TiOptions::default()
            });
            measure(&workload, iterations, || bench_ti(&estimator, &series))?;
        }
        "amber-ti-parallel" => {
            let series = load_amber_dhdl_series()?;
            let estimator = TiEstimator::new(TiOptions {
                parallel: true,
                ..TiOptions::default()
            });
            measure(&workload, iterations, || bench_ti(&estimator, &series))?;
        }
        "gromacs-parse-u-nk" => measure(&workload, iterations, bench_gromacs_parse_u_nk)?,
        "gromacs-decorrelate-u-nk" => {
            let parsed = load_gromacs_u_nk_with_potential()?;
            measure(&workload, iterations, || bench_decorrelate_u_nk(&parsed))?;
        }
        "gromacs-bar-serial" => {
            let windows = prepare_gromacs_decorrelated_windows()?;
            let estimator = BarEstimator::new(BarOptions {
                parallel: false,
                ..BarOptions::default()
            });
            measure(&workload, iterations, || bench_bar(&estimator, &windows))?;
        }
        "gromacs-bar-parallel" => {
            let windows = prepare_gromacs_decorrelated_windows()?;
            let estimator = BarEstimator::new(BarOptions {
                parallel: true,
                ..BarOptions::default()
            });
            measure(&workload, iterations, || bench_bar(&estimator, &windows))?;
        }
        "gromacs-exp-serial" => {
            let windows = prepare_gromacs_decorrelated_windows()?;
            let estimator = IexpEstimator::new(IexpOptions { parallel: false });
            measure(&workload, iterations, || bench_exp(&estimator, &windows))?;
        }
        "gromacs-exp-parallel" => {
            let windows = prepare_gromacs_decorrelated_windows()?;
            let estimator = IexpEstimator::new(IexpOptions { parallel: true });
            measure(&workload, iterations, || bench_exp(&estimator, &windows))?;
        }
        _ => return Err(usage().into()),
    }

    Ok(())
}

fn usage() -> String {
    "usage: profile_estimators <workload> [iterations]".to_string()
}

fn measure(
    workload: &str,
    iterations: usize,
    mut f: impl FnMut() -> Result<usize, DynError>,
) -> Result<(), DynError> {
    let mut checksum = 0usize;
    let start = Instant::now();
    for _ in 0..iterations {
        checksum = checksum.wrapping_add(f()?);
    }
    let elapsed = start.elapsed();
    black_box(checksum);
    println!(
        "workload={workload} iterations={iterations} total_s={:.6} per_iter_ms={:.3} checksum={checksum}",
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1000.0 / iterations as f64
    );
    Ok(())
}

fn bench_amber_parse_u_nk() -> Result<usize, DynError> {
    parse_u_nk_with_potential(&amber_paths()?, 300.0)
}

fn bench_gromacs_parse_u_nk() -> Result<usize, DynError> {
    parse_u_nk_with_potential(&gromacs_paths()?, 298.0)
}

fn bench_amber_parse_dhdl() -> Result<usize, DynError> {
    let mut checksum = 0usize;
    for path in amber_paths()? {
        let series = extract_dhdl(path, 300.0)?;
        checksum = checksum.wrapping_add(series.values().len());
    }
    Ok(checksum)
}

fn bench_decorrelate_u_nk(parsed: &[(UNkMatrix, Vec<f64>)]) -> Result<usize, DynError> {
    let mut checksum = 0usize;
    for (window, epot) in parsed {
        let decorrelated =
            decorrelate_u_nk_with_observable(window, epot, &DecorrelationOptions::default())?;
        checksum = checksum.wrapping_add(decorrelated.n_samples());
    }
    Ok(checksum)
}

fn bench_mbar(estimator: &MbarEstimator, windows: &[UNkMatrix]) -> Result<usize, DynError> {
    let result = estimator.estimate_with_uncertainty(windows)?;
    let delta = result.values()[result.n_states() - 1];
    Ok(result.values().len().wrapping_add(delta.to_bits() as usize))
}

fn bench_uwham(estimator: &UwhamEstimator, windows: &[UNkMatrix]) -> Result<usize, DynError> {
    let fit = estimator.fit(windows)?;
    let result = fit.result_with_uncertainty()?;
    let delta = result.values()[result.n_states() - 1];
    Ok(result.values().len().wrapping_add(delta.to_bits() as usize))
}

fn bench_bar(estimator: &BarEstimator, windows: &[UNkMatrix]) -> Result<usize, DynError> {
    let result = estimator.estimate(windows)?;
    let delta = result.values()[result.n_states() - 1];
    Ok(result.values().len().wrapping_add(delta.to_bits() as usize))
}

fn bench_exp(estimator: &IexpEstimator, windows: &[UNkMatrix]) -> Result<usize, DynError> {
    let result = estimator.estimate_with_uncertainty(windows)?;
    let delta = result.values()[result.n_states() - 1];
    Ok(result.values().len().wrapping_add(delta.to_bits() as usize))
}

fn bench_ti(estimator: &TiEstimator, series: &[DhdlSeries]) -> Result<usize, DynError> {
    let result = estimator.fit(series)?.result()?;
    Ok(result.delta_f().to_bits() as usize)
}

fn load_amber_u_nk_with_potential() -> Result<Vec<(UNkMatrix, Vec<f64>)>, DynError> {
    let paths = amber_paths()?;
    let mut parsed = Vec::with_capacity(paths.len());
    for path in &paths {
        parsed.push(extract_u_nk_with_potential(path, 300.0)?);
    }
    Ok(parsed)
}

fn load_gromacs_u_nk_with_potential() -> Result<Vec<(UNkMatrix, Vec<f64>)>, DynError> {
    let paths = gromacs_paths()?;
    let mut parsed = Vec::with_capacity(paths.len());
    for path in &paths {
        parsed.push(extract_u_nk_with_potential(path, 298.0)?);
    }
    Ok(parsed)
}

fn prepare_gromacs_decorrelated_windows() -> Result<Vec<UNkMatrix>, DynError> {
    let parsed = load_gromacs_u_nk_with_potential()?;
    let mut windows = Vec::with_capacity(parsed.len());
    for (window, epot) in &parsed {
        windows.push(decorrelate_u_nk_with_observable(
            window,
            epot,
            &DecorrelationOptions::default(),
        )?);
    }
    Ok(windows)
}

fn load_amber_dhdl_series() -> Result<Vec<DhdlSeries>, DynError> {
    let paths = amber_paths()?;
    let mut series = Vec::with_capacity(paths.len());
    for path in &paths {
        series.push(extract_dhdl(path, 300.0)?);
    }
    Ok(series)
}

fn amber_windows_only(parsed: Vec<(UNkMatrix, Vec<f64>)>) -> Vec<UNkMatrix> {
    parsed.into_iter().map(|(window, _)| window).collect()
}

fn parse_u_nk_with_potential(paths: &[PathBuf], temperature_k: f64) -> Result<usize, DynError> {
    let mut checksum = 0usize;
    for path in paths {
        let (window, epot) = extract_u_nk_with_potential(path, temperature_k)?;
        checksum = checksum
            .wrapping_add(window.n_samples())
            .wrapping_add(window.n_states())
            .wrapping_add(epot.len());
    }
    Ok(checksum)
}

fn amber_paths() -> Result<Vec<PathBuf>, DynError> {
    let base = manifest_dir().join("fixtures/amber/acetamide_tiny");
    let mut paths = fs::read_dir(base)?
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let path = entry.path();
            if path.is_dir() {
                Some(path.join("acetamide.prod.out"))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    paths.sort_by(|left, right| {
        amber_lambda(left)
            .partial_cmp(&amber_lambda(right))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(paths)
}

fn amber_lambda(path: &Path) -> f64 {
    path.parent()
        .and_then(|parent| parent.file_name())
        .and_then(|name| name.to_str())
        .and_then(|name| name.parse::<f64>().ok())
        .expect("fixture path should encode lambda as folder name")
}

fn gromacs_paths() -> Result<Vec<PathBuf>, DynError> {
    let base = manifest_dir().join("fixtures/gromacs/1k_bar_samples");
    let mut paths = fs::read_dir(base)?
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let path = entry.path();
            if path.is_dir() {
                Some(path.join("dhdl.xvg"))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    paths.sort_by_key(|path| gromacs_lambda_index(path));
    Ok(paths)
}

fn gromacs_lambda_index(path: &Path) -> usize {
    path.parent()
        .and_then(|parent| parent.file_name())
        .and_then(|name| name.to_str())
        .and_then(|name| name.strip_prefix("lambda-"))
        .and_then(|name| name.parse::<usize>().ok())
        .expect("fixture path should encode lambda index as lambda-N")
}

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
