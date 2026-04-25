use std::env;
use std::fs;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::Instant;

use alchemrs::{
    decorrelate_u_nk_with_observable, extract_dhdl, extract_u_nk_with_potential, BarEstimator,
    BarOptions, DecorrelationOptions, DhdlSeries, IexpEstimator, IexpOptions, MbarEstimator,
    MbarOptions, MbarSolver, StatePoint, TiEstimator, TiOptions, UNkMatrix, UwhamEstimator,
    UwhamOptions, UwhamSolver,
};

type DynError = Box<dyn std::error::Error>;

fn main() -> Result<(), DynError> {
    let mut args = env::args().skip(1);
    let workload = args.next().ok_or_else(usage)?;
    let (iterations, scale) = parse_iterations_and_scale(args)?;
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
        "synthetic-mbar-fixed" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = MbarEstimator::new(MbarOptions {
                solver: MbarSolver::FixedPoint,
                ..MbarOptions::default()
            });
            measure_with_scale(&workload, scale, iterations, || {
                bench_mbar_fit_only(&estimator, &windows)
            })?;
        }
        "synthetic-mbar-lbfgs" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = MbarEstimator::new(MbarOptions {
                solver: MbarSolver::Lbfgs,
                ..MbarOptions::default()
            });
            measure_with_scale(&workload, scale, iterations, || {
                bench_mbar_fit_only(&estimator, &windows)
            })?;
        }
        "synthetic-mbar-default" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = MbarEstimator::new(MbarOptions::default());
            measure_with_scale(&workload, scale, iterations, || {
                bench_mbar_fit_only(&estimator, &windows)
            })?;
        }
        "synthetic-mbar-uncertainty" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = MbarEstimator::new(MbarOptions {
                solver: MbarSolver::FixedPoint,
                ..MbarOptions::default()
            });
            measure_with_scale(&workload, scale, iterations, || {
                bench_mbar(&estimator, &windows)
            })?;
        }
        "synthetic-mbar-lbfgs-uncertainty" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = MbarEstimator::new(MbarOptions {
                solver: MbarSolver::Lbfgs,
                ..MbarOptions::default()
            });
            measure_with_scale(&workload, scale, iterations, || {
                bench_mbar(&estimator, &windows)
            })?;
        }
        "synthetic-uwham" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = UwhamEstimator::new(UwhamOptions::default());
            measure_with_scale(&workload, scale, iterations, || {
                bench_uwham(&estimator, &windows)
            })?;
        }
        "synthetic-uwham-newton" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = UwhamEstimator::new(UwhamOptions {
                solver: UwhamSolver::Newton,
                ..UwhamOptions::default()
            });
            measure_with_scale(&workload, scale, iterations, || {
                bench_uwham(&estimator, &windows)
            })?;
        }
        "synthetic-uwham-lbfgs" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = UwhamEstimator::new(UwhamOptions {
                solver: UwhamSolver::Lbfgs,
                ..UwhamOptions::default()
            });
            measure_with_scale(&workload, scale, iterations, || {
                bench_uwham(&estimator, &windows)
            })?;
        }
        "synthetic-uwham-newton-stats" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = UwhamEstimator::new(UwhamOptions {
                solver: UwhamSolver::Newton,
                ..UwhamOptions::default()
            });
            print_uwham_stats(&workload, scale, &estimator, &windows)?;
        }
        "synthetic-uwham-lbfgs-stats" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = UwhamEstimator::new(UwhamOptions {
                solver: UwhamSolver::Lbfgs,
                ..UwhamOptions::default()
            });
            print_uwham_stats(&workload, scale, &estimator, &windows)?;
        }
        "synthetic-uwham-wide-scan" => {
            bench_uwham_wide_scan()?;
        }
        "synthetic-bar" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = BarEstimator::new(BarOptions::default());
            measure_with_scale(&workload, scale, iterations, || {
                bench_bar(&estimator, &windows)
            })?;
        }
        "synthetic-exp" => {
            let windows = synthetic_u_nk_windows(scale)?;
            let estimator = IexpEstimator::new(IexpOptions::default());
            measure_with_scale(&workload, scale, iterations, || {
                bench_exp(&estimator, &windows)
            })?;
        }
        "synthetic-prep-u-nk" => {
            let parsed = synthetic_u_nk_with_observable(scale)?;
            measure_with_scale(&workload, scale, iterations, || {
                bench_decorrelate_u_nk(&parsed)
            })?;
        }
        "synthetic-ti" => {
            let series = synthetic_dhdl_series(scale)?;
            let estimator = TiEstimator::new(TiOptions::default());
            measure_with_scale(&workload, scale, iterations, || {
                bench_ti(&estimator, &series)
            })?;
        }
        _ => return Err(usage().into()),
    }

    Ok(())
}

fn usage() -> String {
    "usage: profile_estimators <workload> [iterations] [scale]\n\
     scales: smoke, medium, large, stress\n\
     synthetic workloads: synthetic-mbar-default, synthetic-mbar-fixed, synthetic-mbar-lbfgs, synthetic-mbar-uncertainty, synthetic-mbar-lbfgs-uncertainty, synthetic-uwham, synthetic-uwham-newton, synthetic-uwham-lbfgs, synthetic-uwham-newton-stats, synthetic-uwham-lbfgs-stats, synthetic-uwham-wide-scan, synthetic-bar, synthetic-exp, synthetic-prep-u-nk, synthetic-ti"
        .to_string()
}

fn parse_iterations_and_scale(
    mut args: impl Iterator<Item = String>,
) -> Result<(usize, SyntheticScale), DynError> {
    let Some(first) = args.next() else {
        return Ok((1, SyntheticScale::Smoke));
    };

    match first.parse::<usize>() {
        Ok(iterations) => {
            let scale = args
                .next()
                .map(|value| SyntheticScale::parse(&value))
                .transpose()?
                .unwrap_or(SyntheticScale::Smoke);
            Ok((iterations, scale))
        }
        Err(_) => Ok((1, SyntheticScale::parse(&first)?)),
    }
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

fn measure_with_scale(
    workload: &str,
    scale: SyntheticScale,
    iterations: usize,
    f: impl FnMut() -> Result<usize, DynError>,
) -> Result<(), DynError> {
    println!(
        "synthetic_scale={} states={} samples_per_window={}",
        scale.name(),
        scale.n_states(),
        scale.samples_per_window()
    );
    measure(workload, iterations, f)
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

fn bench_mbar_fit_only(
    estimator: &MbarEstimator,
    windows: &[UNkMatrix],
) -> Result<usize, DynError> {
    let fit = estimator.fit(windows)?;
    let delta = fit.free_energies()[fit.n_states() - 1];
    Ok(fit.n_states().wrapping_add(delta.to_bits() as usize))
}

fn bench_uwham(estimator: &UwhamEstimator, windows: &[UNkMatrix]) -> Result<usize, DynError> {
    let fit = estimator.fit(windows)?;
    let result = fit.result_with_uncertainty()?;
    let delta = result.values()[result.n_states() - 1];
    Ok(result.values().len().wrapping_add(delta.to_bits() as usize))
}

fn print_uwham_stats(
    workload: &str,
    scale: SyntheticScale,
    estimator: &UwhamEstimator,
    windows: &[UNkMatrix],
) -> Result<(), DynError> {
    println!(
        "synthetic_scale={} states={} samples_per_window={}",
        scale.name(),
        scale.n_states(),
        scale.samples_per_window()
    );
    let start = Instant::now();
    let fit = estimator.fit(windows)?;
    let result = fit.result_with_uncertainty()?;
    let elapsed = start.elapsed();
    let stats = fit.solver_stats();
    let delta = result.values()[result.n_states() - 1];
    let checksum = result.values().len().wrapping_add(delta.to_bits() as usize);
    black_box(checksum);
    println!(
        "workload={workload} solver={:?} total_s={:.6} checksum={checksum} iterations={} objective_evaluations={} gradient_evaluations={} hessian_evaluations={} line_search_evaluations={} accepted_newton_steps={} accepted_gradient_steps={} accepted_lbfgs_steps={} fallback_direction_resets={}",
        stats.solver,
        elapsed.as_secs_f64(),
        stats.iterations,
        stats.objective_evaluations,
        stats.gradient_evaluations,
        stats.hessian_evaluations,
        stats.line_search_evaluations,
        stats.accepted_newton_steps,
        stats.accepted_gradient_steps,
        stats.accepted_lbfgs_steps,
        stats.fallback_direction_resets
    );
    Ok(())
}

fn bench_uwham_wide_scan() -> Result<(), DynError> {
    for &(n_states, n_samples) in &[(16, 20), (32, 20), (64, 20), (96, 10), (128, 10), (160, 8)] {
        let windows = synthetic_u_nk_windows_custom(n_states, n_samples)?;
        for solver in [UwhamSolver::Newton, UwhamSolver::Lbfgs] {
            let estimator = UwhamEstimator::new(UwhamOptions {
                solver,
                ..UwhamOptions::default()
            });
            let start = Instant::now();
            let fit = estimator.fit(&windows)?;
            let result = fit.result()?;
            let elapsed = start.elapsed();
            let stats = fit.solver_stats();
            let delta = result.values()[result.n_states() - 1];
            let checksum = result.values().len().wrapping_add(delta.to_bits() as usize);
            black_box(checksum);
            println!(
                "workload=synthetic-uwham-wide-scan solver={solver:?} states={n_states} samples_per_window={n_samples} observations={} total_ms={:.3} checksum={checksum} iterations={} objective_evaluations={} gradient_evaluations={} hessian_evaluations={} line_search_evaluations={}",
                n_states * n_samples,
                elapsed.as_secs_f64() * 1000.0,
                stats.iterations,
                stats.objective_evaluations,
                stats.gradient_evaluations,
                stats.hessian_evaluations,
                stats.line_search_evaluations
            );
        }
    }
    Ok(())
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

#[derive(Debug, Clone, Copy)]
enum SyntheticScale {
    Smoke,
    Medium,
    Large,
    Stress,
}

impl SyntheticScale {
    fn parse(value: &str) -> Result<Self, DynError> {
        match value {
            "smoke" => Ok(Self::Smoke),
            "medium" => Ok(Self::Medium),
            "large" => Ok(Self::Large),
            "stress" => Ok(Self::Stress),
            _ => Err(format!("unknown synthetic scale '{value}'").into()),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Smoke => "smoke",
            Self::Medium => "medium",
            Self::Large => "large",
            Self::Stress => "stress",
        }
    }

    fn n_states(self) -> usize {
        match self {
            Self::Smoke => 8,
            Self::Medium => 16,
            Self::Large => 32,
            Self::Stress => 64,
        }
    }

    fn samples_per_window(self) -> usize {
        match self {
            Self::Smoke => 500,
            Self::Medium => 2_000,
            Self::Large => 5_000,
            Self::Stress => 10_000,
        }
    }
}

fn synthetic_u_nk_with_observable(
    scale: SyntheticScale,
) -> Result<Vec<(UNkMatrix, Vec<f64>)>, DynError> {
    let states = synthetic_states(scale.n_states())?;
    let mut windows = Vec::with_capacity(states.len());
    for sampled_index in 0..states.len() {
        let window = synthetic_u_nk_window(&states, sampled_index, scale.samples_per_window())?;
        let observable = synthetic_observable(sampled_index, scale.samples_per_window());
        windows.push((window, observable));
    }
    Ok(windows)
}

fn synthetic_u_nk_windows(scale: SyntheticScale) -> Result<Vec<UNkMatrix>, DynError> {
    synthetic_u_nk_windows_custom(scale.n_states(), scale.samples_per_window())
}

fn synthetic_u_nk_windows_custom(
    n_states: usize,
    samples_per_window: usize,
) -> Result<Vec<UNkMatrix>, DynError> {
    let states = synthetic_states(n_states)?;
    (0..states.len())
        .map(|sampled_index| synthetic_u_nk_window(&states, sampled_index, samples_per_window))
        .collect()
}

fn synthetic_u_nk_window(
    states: &[StatePoint],
    sampled_index: usize,
    n_samples: usize,
) -> Result<UNkMatrix, DynError> {
    let n_states = states.len();
    let sampled_lambda = states[sampled_index].lambdas()[0];
    let mut data = Vec::with_capacity(n_samples * n_states);
    let mut time_ps = Vec::with_capacity(n_samples);

    for sample_index in 0..n_samples {
        time_ps.push(sample_index as f64);
        let phase = sample_index as f64 * 0.017 + sampled_index as f64 * 0.31;
        let coordinate = sampled_lambda + 0.08 * phase.sin() + 0.03 * (phase * 0.37).cos();
        let common = 0.02 * (sample_index as f64 * 0.011).sin();

        for state in states {
            let lambda = state.lambdas()[0];
            let displacement = lambda - coordinate;
            let state_bias = 0.25 * lambda + 0.05 * lambda * lambda;
            data.push(18.0 * displacement * displacement + state_bias + common);
        }
    }

    UNkMatrix::new_with_labels(
        n_samples,
        n_states,
        data,
        time_ps,
        Some(states[sampled_index].clone()),
        states.to_vec(),
        Some(vec!["lambda".to_string()]),
    )
    .map_err(Into::into)
}

fn synthetic_dhdl_series(scale: SyntheticScale) -> Result<Vec<DhdlSeries>, DynError> {
    let states = synthetic_states(scale.n_states())?;
    let n_samples = scale.samples_per_window();
    let mut series = Vec::with_capacity(states.len());

    for (state_index, state) in states.into_iter().enumerate() {
        let lambda = state.lambdas()[0];
        let mut time_ps = Vec::with_capacity(n_samples);
        let mut values = Vec::with_capacity(n_samples);
        for sample_index in 0..n_samples {
            time_ps.push(sample_index as f64);
            let phase = sample_index as f64 * 0.013 + state_index as f64 * 0.19;
            values.push(0.4 + 1.2 * lambda + 0.08 * phase.sin() + 0.03 * (phase * 0.41).cos());
        }
        series.push(DhdlSeries::new(state, time_ps, values)?);
    }

    Ok(series)
}

fn synthetic_states(n_states: usize) -> Result<Vec<StatePoint>, DynError> {
    let denom = n_states.saturating_sub(1) as f64;
    (0..n_states)
        .map(|index| {
            let lambda = if denom == 0.0 {
                0.0
            } else {
                index as f64 / denom
            };
            StatePoint::new(vec![lambda], 300.0).map_err(Into::into)
        })
        .collect()
}

fn synthetic_observable(sampled_index: usize, n_samples: usize) -> Vec<f64> {
    (0..n_samples)
        .map(|sample_index| {
            let phase = sample_index as f64 * 0.021 + sampled_index as f64 * 0.17;
            0.5 + 0.1 * sampled_index as f64 + phase.sin() + 0.2 * (phase * 0.23).cos()
        })
        .collect()
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
