use std::collections::BTreeMap;
use std::env;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use alchemrs::estimators::{IntegrationMethod, TiEstimator, TiOptions};
use alchemrs::parse::amber::extract_dhdl;
use alchemrs::{advise_ti_schedule, DhdlSeries, TiScheduleAdvisorOptions, TiSuggestionKind};
use rayon::prelude::*;

const DEFAULT_TEMPERATURE_K: f64 = 300.0;
const DEFAULT_ERROR_TOLERANCE_KCAL_MOL: f64 = 1.0;
const DEFAULT_OUTPUT_DIR: &str = "target/ti-advisor-calibration";
const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;

#[derive(Debug, Clone)]
struct Config {
    dataset_root: PathBuf,
    output_dir: PathBuf,
    temperature_k: f64,
    error_tolerance_reduced: f64,
    error_tolerance_kcal_mol: f64,
    exhaustive_drops: bool,
}

#[derive(Debug, Clone)]
struct Replicate {
    system: String,
    replicate: String,
    path: PathBuf,
    series: Vec<DhdlSeries>,
    delta_f: f64,
}

#[derive(Debug, Clone)]
struct CalibrationCase {
    system: String,
    replicate: String,
    case_kind: String,
    n_windows: usize,
    samples_per_window_min: usize,
    reference_delta_f: f64,
    case_delta_f: f64,
    abs_error: f64,
    true_window_issue: bool,
    true_sampling_issue: bool,
    series: Vec<DhdlSeries>,
}

#[derive(Debug, Clone, Copy)]
struct Prediction {
    window_issue: bool,
    sampling_issue: bool,
}

#[derive(Debug, Clone)]
struct SweepResult {
    options: TiScheduleAdvisorOptions,
    window_precision: f64,
    window_recall: f64,
    window_f1: f64,
    sampling_precision: f64,
    sampling_recall: f64,
    sampling_f1: f64,
    macro_f1: f64,
    false_negatives: usize,
    false_positives: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    let config = parse_args()?;
    fs::create_dir_all(&config.output_dir)?;

    let replicate_dirs = find_replicate_dirs(&config.dataset_root)?;
    if replicate_dirs.is_empty() {
        return Err(format!(
            "no replicate directories found under {}",
            config.dataset_root.display()
        )
        .into());
    }

    let replicates = load_replicates(&replicate_dirs, config.temperature_k)?;
    if replicates.is_empty() {
        return Err("no complete TI replicate data could be parsed".into());
    }

    let references = reference_delta_f_by_system(&replicates);
    if config.exhaustive_drops {
        run_exhaustive_drop_scan(&config, replicate_dirs.len(), &replicates, &references)?;
        return Ok(());
    }

    let cases = build_cases(&replicates, &references, config.error_tolerance_reduced)?;
    if cases.is_empty() {
        return Err("no calibration cases were generated".into());
    }

    write_case_csv(&config.output_dir.join("cases.csv"), &cases)?;

    let sweep = sweep_options(&cases)?;
    write_sweep_csv(&config.output_dir.join("threshold-sweep.csv"), &sweep)?;
    write_report(
        &config.output_dir.join("report.md"),
        &config,
        replicate_dirs.len(),
        &replicates,
        &cases,
        &sweep,
    )?;

    if let Some(best) = sweep.first() {
        println!("replicate_dirs={}", replicate_dirs.len());
        println!("parsed_replicates={}", replicates.len());
        println!("calibration_cases={}", cases.len());
        println!(
            "best block_cv_min={:.3} slope_z_min={:.2} curvature_z_min={:.2} interval_uncertainty_z_min={:.2} macro_f1={:.3}",
            best.options.block_cv_min,
            best.options.slope_z_min,
            best.options.curvature_z_min,
            best.options.interval_uncertainty_z_min,
            best.macro_f1
        );
        println!("wrote {}", config.output_dir.display());
    }

    Ok(())
}

fn parse_args() -> Result<Config, Box<dyn Error>> {
    let mut dataset_root = None;
    let mut output_dir = PathBuf::from(DEFAULT_OUTPUT_DIR);
    let mut temperature_k = DEFAULT_TEMPERATURE_K;
    let mut error_tolerance_reduced = None;
    let mut error_tolerance_kcal_mol = DEFAULT_ERROR_TOLERANCE_KCAL_MOL;
    let mut exhaustive_drops = false;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output-dir" => {
                output_dir = PathBuf::from(next_arg(&mut args, "--output-dir")?);
            }
            "--temperature-k" => {
                temperature_k = next_arg(&mut args, "--temperature-k")?.parse()?;
            }
            "--error-tolerance" => {
                error_tolerance_reduced = Some(next_arg(&mut args, "--error-tolerance")?.parse()?);
            }
            "--error-tolerance-kcal" => {
                error_tolerance_kcal_mol =
                    next_arg(&mut args, "--error-tolerance-kcal")?.parse()?;
            }
            "--exhaustive-drops" => {
                exhaustive_drops = true;
            }
            "--help" | "-h" => {
                print_usage();
                std::process::exit(0);
            }
            _ if arg.starts_with('-') => {
                return Err(format!("unknown option: {arg}").into());
            }
            _ => {
                if dataset_root.replace(PathBuf::from(arg)).is_some() {
                    return Err("expected a single dataset root".into());
                }
            }
        }
    }

    let dataset_root = dataset_root.ok_or_else(|| {
        "usage: cargo run --example calibrate_ti_advisor -- <dataset-root> [--output-dir DIR] [--temperature-k K] [--error-tolerance-kcal KCAL] [--error-tolerance REDUCED_DG]"
            .to_string()
    })?;

    let converted_tolerance = error_tolerance_kcal_mol / (K_B_KCAL_PER_MOL_K * temperature_k);
    let error_tolerance_reduced = error_tolerance_reduced.unwrap_or(converted_tolerance);
    let error_tolerance_kcal_mol = error_tolerance_reduced * K_B_KCAL_PER_MOL_K * temperature_k;

    Ok(Config {
        dataset_root,
        output_dir,
        temperature_k,
        error_tolerance_reduced,
        error_tolerance_kcal_mol,
        exhaustive_drops,
    })
}

fn next_arg(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, Box<dyn Error>> {
    args.next()
        .ok_or_else(|| format!("{flag} requires a value").into())
}

fn print_usage() {
    eprintln!(
        "usage: cargo run --example calibrate_ti_advisor -- <dataset-root> [--output-dir DIR] [--temperature-k K] [--error-tolerance-kcal KCAL] [--error-tolerance REDUCED_DG] [--exhaustive-drops]"
    );
}

fn find_replicate_dirs(root: &Path) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let mut out = Vec::new();
    for system_entry in fs::read_dir(root)? {
        let flexible_ti = system_entry?.path().join("flexible_TI");
        if !flexible_ti.is_dir() {
            continue;
        }
        for replicate_entry in fs::read_dir(flexible_ti)? {
            let replicate_path = replicate_entry?.path();
            if replicate_path
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with("replicate_"))
                && replicate_path.is_dir()
            {
                out.push(replicate_path);
            }
        }
    }
    out.sort();
    Ok(out)
}

fn load_replicates(
    paths: &[PathBuf],
    temperature_k: f64,
) -> Result<Vec<Replicate>, Box<dyn Error>> {
    let loaded = paths
        .par_iter()
        .map(|path| parse_replicate(path, temperature_k))
        .collect::<Vec<_>>();

    let mut replicates = Vec::new();
    for item in loaded {
        if let Some(replicate) = item.map_err(|err| -> Box<dyn Error> { err.into() })? {
            replicates.push(replicate);
        }
    }
    replicates.sort_by(|left, right| {
        left.system
            .cmp(&right.system)
            .then_with(|| left.replicate.cmp(&right.replicate))
    });

    Ok(replicates)
}

fn parse_replicate(path: &Path, temperature_k: f64) -> Result<Option<Replicate>, String> {
    let estimator = TiEstimator::new(TiOptions {
        method: IntegrationMethod::Trapezoidal,
        parallel: false,
    });

    let inputs = amber_outputs_for_replicate(path).map_err(|err| err.to_string())?;
    if inputs.len() < 2 {
        return Ok(None);
    }
    let mut series = Vec::with_capacity(inputs.len());
    for input in inputs {
        series.push(extract_dhdl(&input, temperature_k).map_err(|err| err.to_string())?);
    }
    series.sort_by(|a, b| lambda_of(a).total_cmp(&lambda_of(b)));
    let delta_f = estimator
        .fit(&series)
        .map_err(|err| err.to_string())?
        .delta_f();
    Ok(Some(Replicate {
        system: system_name(path).map_err(|err| err.to_string())?,
        replicate: path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("replicate")
            .to_string(),
        path: path.to_path_buf(),
        series,
        delta_f,
    }))
}

fn amber_outputs_for_replicate(path: &Path) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let mut out = Vec::new();
    for lambda_entry in fs::read_dir(path)? {
        let lambda_path = lambda_entry?.path();
        if !lambda_path.is_dir() {
            continue;
        }
        for entry in fs::read_dir(lambda_path)? {
            let entry_path = entry?.path();
            if entry_path
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.ends_with(".prod.out"))
            {
                out.push(entry_path);
            }
        }
    }
    out.sort_by(|a, b| path_lambda(a).total_cmp(&path_lambda(b)));
    Ok(out)
}

fn path_lambda(path: &Path) -> f64 {
    path.parent()
        .and_then(Path::file_name)
        .and_then(|name| name.to_str())
        .and_then(|name| name.parse::<f64>().ok())
        .unwrap_or(f64::INFINITY)
}

fn system_name(replicate_path: &Path) -> Result<String, Box<dyn Error>> {
    Ok(replicate_path
        .parent()
        .and_then(Path::parent)
        .and_then(Path::file_name)
        .and_then(|name| name.to_str())
        .ok_or_else(|| {
            format!(
                "could not determine system for {}",
                replicate_path.display()
            )
        })?
        .to_string())
}

fn lambda_of(series: &DhdlSeries) -> f64 {
    series
        .state()
        .lambdas()
        .first()
        .copied()
        .unwrap_or(f64::NAN)
}

fn reference_delta_f_by_system(replicates: &[Replicate]) -> BTreeMap<String, f64> {
    let mut grouped: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for replicate in replicates {
        grouped
            .entry(replicate.system.clone())
            .or_default()
            .push(replicate.delta_f);
    }
    grouped
        .into_iter()
        .map(|(system, values)| {
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            (system, mean)
        })
        .collect()
}

fn build_cases(
    replicates: &[Replicate],
    references: &BTreeMap<String, f64>,
    error_tolerance: f64,
) -> Result<Vec<CalibrationCase>, Box<dyn Error>> {
    let mut cases = Vec::new();
    let estimator = TiEstimator::new(TiOptions::default());

    for replicate in replicates {
        let reference_delta_f = references[&replicate.system];
        cases.push(make_case(
            replicate,
            "full",
            replicate.series.clone(),
            reference_delta_f,
            error_tolerance,
            false,
            false,
            &estimator,
        )?);

        for fraction in [0.25, 0.50, 0.75] {
            let truncated = replicate
                .series
                .iter()
                .map(|series| truncate_series(series, fraction))
                .collect::<Result<Vec<_>, _>>()?;
            cases.push(make_case(
                replicate,
                &format!("truncate_{:.0}pct", fraction * 100.0),
                truncated,
                reference_delta_f,
                error_tolerance,
                false,
                true,
                &estimator,
            )?);
        }

        if replicate.series.len() >= 5 {
            let thinned = replicate
                .series
                .iter()
                .enumerate()
                .filter_map(|(idx, series)| {
                    let keep_endpoint = idx == 0 || idx + 1 == replicate.series.len();
                    if keep_endpoint || idx % 2 == 0 {
                        Some(series.clone())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            cases.push(make_case(
                replicate,
                "drop_alternating",
                thinned,
                reference_delta_f,
                error_tolerance,
                true,
                false,
                &estimator,
            )?);
        }

        for dropped in representative_drop_indices(replicate.series.len()) {
            let missing_one = replicate
                .series
                .iter()
                .enumerate()
                .filter_map(|(idx, series)| {
                    if idx == dropped {
                        None
                    } else {
                        Some(series.clone())
                    }
                })
                .collect::<Vec<_>>();
            cases.push(make_case(
                replicate,
                &format!("drop_lambda_{:.4}", lambda_of(&replicate.series[dropped])),
                missing_one,
                reference_delta_f,
                error_tolerance,
                true,
                false,
                &estimator,
            )?);
        }
    }

    Ok(cases)
}

fn representative_drop_indices(n_windows: usize) -> Vec<usize> {
    if n_windows <= 2 {
        return Vec::new();
    }
    let mut indices = vec![n_windows / 4, n_windows / 2, (3 * n_windows) / 4];
    indices.retain(|idx| *idx > 0 && *idx + 1 < n_windows);
    indices.sort_unstable();
    indices.dedup();
    indices
}

fn truncate_series(series: &DhdlSeries, fraction: f64) -> Result<DhdlSeries, Box<dyn Error>> {
    let keep = ((series.values().len() as f64) * fraction).round() as usize;
    let keep = keep.clamp(2, series.values().len());
    Ok(DhdlSeries::new(
        series.state().clone(),
        series.time_ps()[..keep].to_vec(),
        series.values()[..keep].to_vec(),
    )?)
}

#[derive(Debug, Clone)]
struct ExhaustiveWorstCase {
    system: String,
    replicate: String,
    truncation_fraction: f64,
    n_windows: usize,
    dropped_windows: usize,
    case_delta_f: f64,
    abs_error_reduced: f64,
}

#[derive(Debug, Clone)]
struct ExhaustiveSummary {
    replicate_dir_count: usize,
    parsed_replicates: usize,
    systems: usize,
    total_cases: usize,
    failures: usize,
    worst: ExhaustiveWorstCase,
    window_importance: Vec<WindowImportance>,
}

#[derive(Debug, Clone)]
struct WindowImportance {
    lambda: f64,
    paired_comparisons: usize,
    failures_introduced: usize,
    failures_removed: usize,
    mean_abs_error_delta_reduced: f64,
    mean_abs_error_delta_kcal_mol: f64,
    max_abs_error_when_dropped_kcal_mol: f64,
}

fn run_exhaustive_drop_scan(
    config: &Config,
    replicate_dir_count: usize,
    replicates: &[Replicate],
    references: &BTreeMap<String, f64>,
) -> Result<(), Box<dyn Error>> {
    let csv_path = config.output_dir.join("exhaustive-drop-cases.csv");
    let mut writer = csv::Writer::from_path(&csv_path)?;
    writer.write_record([
        "system",
        "replicate",
        "truncation_fraction",
        "n_windows",
        "dropped_windows",
        "reference_delta_f",
        "case_delta_f",
        "abs_error_reduced",
        "abs_error_kcal_mol",
        "exceeds_tolerance",
    ])?;

    let mut total_cases = 0usize;
    let mut failures = 0usize;
    let mut worst: Option<ExhaustiveWorstCase> = None;
    let mut window_importance = Vec::new();

    for replicate in replicates {
        let reference_delta_f = references[&replicate.system];
        for truncation_fraction in [1.0, 0.75, 0.50, 0.25] {
            let reduced = precompute_truncated_windows(&replicate.series, truncation_fraction)?;
            let n_windows = reduced.len();
            if n_windows < 2 {
                continue;
            }
            let interior_count = n_windows - 2;
            let subset_count = 1usize << interior_count;
            let mut case_errors = vec![0.0; subset_count];
            let mut case_failures = vec![false; subset_count];
            for mask in 0..subset_count {
                let case_delta_f = integrate_subset_trapezoidal(&reduced, mask);
                let abs_error_reduced = (case_delta_f - reference_delta_f).abs();
                let abs_error_kcal_mol =
                    abs_error_reduced * K_B_KCAL_PER_MOL_K * config.temperature_k;
                let exceeds_tolerance = abs_error_reduced >= config.error_tolerance_reduced;
                let kept_interior = mask.count_ones() as usize;
                let kept_windows = kept_interior + 2;
                let dropped_windows = n_windows - kept_windows;

                total_cases += 1;
                if exceeds_tolerance {
                    failures += 1;
                }
                case_errors[mask] = abs_error_reduced;
                case_failures[mask] = exceeds_tolerance;
                if worst
                    .as_ref()
                    .is_none_or(|item| abs_error_reduced > item.abs_error_reduced)
                {
                    worst = Some(ExhaustiveWorstCase {
                        system: replicate.system.clone(),
                        replicate: replicate.replicate.clone(),
                        truncation_fraction,
                        n_windows: kept_windows,
                        dropped_windows,
                        case_delta_f,
                        abs_error_reduced,
                    });
                }

                writer.write_record([
                    replicate.system.as_str(),
                    replicate.replicate.as_str(),
                    &format!("{truncation_fraction:.2}"),
                    &kept_windows.to_string(),
                    &dropped_windows.to_string(),
                    &format!("{reference_delta_f:.8}"),
                    &format!("{case_delta_f:.8}"),
                    &format!("{abs_error_reduced:.8}"),
                    &format!("{abs_error_kcal_mol:.8}"),
                    &exceeds_tolerance.to_string(),
                ])?;
            }
            accumulate_window_importance(
                &reduced,
                &case_errors,
                &case_failures,
                config.temperature_k,
                &mut window_importance,
            );
        }
    }
    writer.flush()?;

    let summary = ExhaustiveSummary {
        replicate_dir_count,
        parsed_replicates: replicates.len(),
        systems: references.len(),
        total_cases,
        failures,
        worst: worst.ok_or("no exhaustive cases were generated")?,
        window_importance: summarize_window_importance(window_importance),
    };
    write_window_importance_csv(
        &config.output_dir.join("window-importance.csv"),
        &summary.window_importance,
    )?;
    write_exhaustive_report(
        &config.output_dir.join("exhaustive-drop-report.md"),
        config,
        &summary,
    )?;

    println!("replicate_dirs={}", summary.replicate_dir_count);
    println!("parsed_replicates={}", summary.parsed_replicates);
    println!("exhaustive_cases={}", summary.total_cases);
    println!("failures_over_tolerance={}", summary.failures);
    println!(
        "max_abs_error_reduced={:.6} max_abs_error_kcal={:.6}",
        summary.worst.abs_error_reduced,
        summary.worst.abs_error_reduced * K_B_KCAL_PER_MOL_K * config.temperature_k
    );
    println!("wrote {}", config.output_dir.display());

    Ok(())
}

fn accumulate_window_importance(
    windows: &[ReducedWindow],
    case_errors: &[f64],
    case_failures: &[bool],
    temperature_k: f64,
    out: &mut Vec<WindowImportance>,
) {
    let interior_count = windows.len().saturating_sub(2);
    for interior_idx in 0..interior_count {
        let bit = 1usize << interior_idx;
        let lambda = windows[interior_idx + 1].lambda;
        let mut paired_comparisons = 0usize;
        let mut failures_introduced = 0usize;
        let mut failures_removed = 0usize;
        let mut delta_sum = 0.0;
        let mut max_abs_error_when_dropped: f64 = 0.0;

        for mask_with_window in 0..case_errors.len() {
            if (mask_with_window & bit) == 0 {
                continue;
            }
            let mask_without_window = mask_with_window & !bit;
            let kept_error = case_errors[mask_with_window];
            let dropped_error = case_errors[mask_without_window];
            let kept_failure = case_failures[mask_with_window];
            let dropped_failure = case_failures[mask_without_window];

            paired_comparisons += 1;
            delta_sum += dropped_error - kept_error;
            max_abs_error_when_dropped = max_abs_error_when_dropped.max(dropped_error);
            match (kept_failure, dropped_failure) {
                (false, true) => failures_introduced += 1,
                (true, false) => failures_removed += 1,
                _ => {}
            }
        }

        out.push(WindowImportance {
            lambda,
            paired_comparisons,
            failures_introduced,
            failures_removed,
            mean_abs_error_delta_reduced: delta_sum / paired_comparisons as f64,
            mean_abs_error_delta_kcal_mol: delta_sum / paired_comparisons as f64
                * K_B_KCAL_PER_MOL_K
                * temperature_k,
            max_abs_error_when_dropped_kcal_mol: max_abs_error_when_dropped
                * K_B_KCAL_PER_MOL_K
                * temperature_k,
        });
    }
}

fn summarize_window_importance(items: Vec<WindowImportance>) -> Vec<WindowImportance> {
    let mut grouped: BTreeMap<String, Vec<WindowImportance>> = BTreeMap::new();
    for item in items {
        grouped
            .entry(format!("{:.8}", item.lambda))
            .or_default()
            .push(item);
    }
    let mut out = grouped
        .into_iter()
        .map(|(lambda, values)| {
            let paired_comparisons = values
                .iter()
                .map(|item| item.paired_comparisons)
                .sum::<usize>();
            let weighted_delta_reduced = values
                .iter()
                .map(|item| item.mean_abs_error_delta_reduced * item.paired_comparisons as f64)
                .sum::<f64>()
                / paired_comparisons as f64;
            let weighted_delta_kcal = values
                .iter()
                .map(|item| item.mean_abs_error_delta_kcal_mol * item.paired_comparisons as f64)
                .sum::<f64>()
                / paired_comparisons as f64;
            WindowImportance {
                lambda: lambda.parse().expect("lambda key is numeric"),
                paired_comparisons,
                failures_introduced: values.iter().map(|item| item.failures_introduced).sum(),
                failures_removed: values.iter().map(|item| item.failures_removed).sum(),
                mean_abs_error_delta_reduced: weighted_delta_reduced,
                mean_abs_error_delta_kcal_mol: weighted_delta_kcal,
                max_abs_error_when_dropped_kcal_mol: values
                    .iter()
                    .map(|item| item.max_abs_error_when_dropped_kcal_mol)
                    .fold(0.0, f64::max),
            }
        })
        .collect::<Vec<_>>();
    out.sort_by(|left, right| {
        right
            .failures_introduced
            .cmp(&left.failures_introduced)
            .then_with(|| {
                right
                    .mean_abs_error_delta_kcal_mol
                    .total_cmp(&left.mean_abs_error_delta_kcal_mol)
            })
    });
    out
}

fn write_window_importance_csv(
    path: &Path,
    items: &[WindowImportance],
) -> Result<(), Box<dyn Error>> {
    let mut writer = csv::Writer::from_path(path)?;
    writer.write_record([
        "lambda",
        "paired_comparisons",
        "failures_introduced",
        "failures_removed",
        "net_failures_introduced",
        "mean_abs_error_delta_reduced",
        "mean_abs_error_delta_kcal_mol",
        "max_abs_error_when_dropped_kcal_mol",
    ])?;
    for item in items {
        writer.write_record([
            &format!("{:.4}", item.lambda),
            &item.paired_comparisons.to_string(),
            &item.failures_introduced.to_string(),
            &item.failures_removed.to_string(),
            &(item.failures_introduced as isize - item.failures_removed as isize).to_string(),
            &format!("{:.8}", item.mean_abs_error_delta_reduced),
            &format!("{:.8}", item.mean_abs_error_delta_kcal_mol),
            &format!("{:.8}", item.max_abs_error_when_dropped_kcal_mol),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct ReducedWindow {
    lambda: f64,
    mean_dhdl: f64,
}

fn precompute_truncated_windows(
    series: &[DhdlSeries],
    fraction: f64,
) -> Result<Vec<ReducedWindow>, Box<dyn Error>> {
    series
        .iter()
        .map(|item| {
            let keep = ((item.values().len() as f64) * fraction).round() as usize;
            let keep = keep.clamp(2, item.values().len());
            let mean_dhdl = item.values()[..keep].iter().sum::<f64>() / keep as f64;
            Ok(ReducedWindow {
                lambda: lambda_of(item),
                mean_dhdl,
            })
        })
        .collect()
}

fn integrate_subset_trapezoidal(windows: &[ReducedWindow], interior_mask: usize) -> f64 {
    let mut previous = windows[0];
    let mut delta_f = 0.0;
    for interior_idx in 0..windows.len().saturating_sub(2) {
        if (interior_mask & (1usize << interior_idx)) != 0 {
            let current = windows[interior_idx + 1];
            delta_f += trapezoid(previous, current);
            previous = current;
        }
    }
    delta_f + trapezoid(previous, *windows.last().expect("at least two windows"))
}

fn trapezoid(left: ReducedWindow, right: ReducedWindow) -> f64 {
    0.5 * (left.mean_dhdl + right.mean_dhdl) * (right.lambda - left.lambda)
}

fn write_exhaustive_report(
    path: &Path,
    config: &Config,
    summary: &ExhaustiveSummary,
) -> Result<(), Box<dyn Error>> {
    let failure_rate = if summary.total_cases == 0 {
        0.0
    } else {
        summary.failures as f64 / summary.total_cases as f64
    };
    let worst_error_kcal =
        summary.worst.abs_error_reduced * K_B_KCAL_PER_MOL_K * config.temperature_k;
    let report = format!(
        "# Exhaustive TI Window-Drop Stress Test\n\n\
         ## Dataset\n\n\
         - Root: `{}`\n\
         - Strict replicate directories found: {}\n\
         - Parsed replicates: {}\n\
         - Systems represented: {}\n\
         - Temperature: {:.3} K\n\
         - Error tolerance: {:.6} reduced units\n\
         - Error tolerance: {:.3} kcal/mol\n\n\
         ## Exhaustive Search\n\n\
         - Truncation fractions tested: 1.00, 0.75, 0.50, 0.25\n\
         - Window subsets: every interior lambda-window subset, endpoints retained\n\
         - Total cases: {}\n\
         - Cases exceeding tolerance: {}\n\
         - Failure rate: {:.6}\n\n\
         ## Worst Case\n\n\
         - System: `{}`\n\
         - Replicate: `{}`\n\
         - Truncation fraction: {:.2}\n\
         - Kept windows: {}\n\
         - Dropped windows: {}\n\
         - Case delta_f: {:.8} reduced units\n\
         - Absolute error: {:.8} reduced units\n\
         - Absolute error: {:.8} kcal/mol\n\n\
         ## Most Important Windows\n\n\
         Ranked by paired failures introduced when dropping exactly that interior window from otherwise identical subsets:\n\n{}\n\n\
         See `exhaustive-drop-cases.csv` for all cases and `window-importance.csv` for the full marginal ranking.\n",
        config.dataset_root.display(),
        summary.replicate_dir_count,
        summary.parsed_replicates,
        summary.systems,
        config.temperature_k,
        config.error_tolerance_reduced,
        config.error_tolerance_kcal_mol,
        summary.total_cases,
        summary.failures,
        failure_rate,
        summary.worst.system,
        summary.worst.replicate,
        summary.worst.truncation_fraction,
        summary.worst.n_windows,
        summary.worst.dropped_windows,
        summary.worst.case_delta_f,
        summary.worst.abs_error_reduced,
        worst_error_kcal,
        format_window_importance_table(&summary.window_importance),
    );
    fs::write(path, report)?;
    Ok(())
}

fn format_window_importance_table(items: &[WindowImportance]) -> String {
    let mut lines = vec![
        "| lambda | failures introduced | failures removed | mean error delta kcal/mol | max error when dropped kcal/mol |".to_string(),
        "| ---: | ---: | ---: | ---: | ---: |".to_string(),
    ];
    for item in items.iter().take(8) {
        lines.push(format!(
            "| {:.4} | {} | {} | {:.6} | {:.6} |",
            item.lambda,
            item.failures_introduced,
            item.failures_removed,
            item.mean_abs_error_delta_kcal_mol,
            item.max_abs_error_when_dropped_kcal_mol
        ));
    }
    lines.join("\n")
}

#[allow(clippy::too_many_arguments)]
fn make_case(
    replicate: &Replicate,
    case_kind: &str,
    series: Vec<DhdlSeries>,
    reference_delta_f: f64,
    error_tolerance: f64,
    candidate_window_issue: bool,
    candidate_sampling_issue: bool,
    estimator: &TiEstimator,
) -> Result<CalibrationCase, Box<dyn Error>> {
    let case_delta_f = estimator.fit(&series)?.delta_f();
    let abs_error = (case_delta_f - reference_delta_f).abs();
    let true_issue = abs_error >= error_tolerance;
    Ok(CalibrationCase {
        system: replicate.system.clone(),
        replicate: replicate.replicate.clone(),
        case_kind: case_kind.to_string(),
        n_windows: series.len(),
        samples_per_window_min: series
            .iter()
            .map(|item| item.values().len())
            .min()
            .unwrap_or_default(),
        reference_delta_f,
        case_delta_f,
        abs_error,
        true_window_issue: candidate_window_issue && true_issue,
        true_sampling_issue: candidate_sampling_issue && true_issue,
        series,
    })
}

fn sweep_options(cases: &[CalibrationCase]) -> Result<Vec<SweepResult>, Box<dyn Error>> {
    let block_cv_values = [0.05, 0.10, 0.15, 0.20, 0.25];
    let z_values = [1.00, 1.50, 2.00];
    let mut results = Vec::new();

    for block_cv_min in block_cv_values {
        for slope_z_min in z_values {
            for curvature_z_min in z_values {
                for interval_uncertainty_z_min in z_values {
                    let options = TiScheduleAdvisorOptions {
                        n_blocks: 4,
                        block_cv_min,
                        curvature_z_min,
                        slope_z_min,
                        interval_uncertainty_z_min,
                        suggest_midpoints: true,
                    };
                    results.push(score_options(cases, options)?);
                }
            }
        }
    }

    results.sort_by(|left, right| {
        right
            .macro_f1
            .total_cmp(&left.macro_f1)
            .then_with(|| left.false_negatives.cmp(&right.false_negatives))
            .then_with(|| left.false_positives.cmp(&right.false_positives))
    });
    Ok(results)
}

fn score_options(
    cases: &[CalibrationCase],
    options: TiScheduleAdvisorOptions,
) -> Result<SweepResult, Box<dyn Error>> {
    let mut window_counts = Counts::default();
    let mut sampling_counts = Counts::default();

    for case in cases {
        let prediction = predict(&case.series, options)?;
        window_counts.add(case.true_window_issue, prediction.window_issue);
        sampling_counts.add(case.true_sampling_issue, prediction.sampling_issue);
    }

    let window_precision = window_counts.precision();
    let window_recall = window_counts.recall();
    let window_f1 = window_counts.f1();
    let sampling_precision = sampling_counts.precision();
    let sampling_recall = sampling_counts.recall();
    let sampling_f1 = sampling_counts.f1();

    Ok(SweepResult {
        options,
        window_precision,
        window_recall,
        window_f1,
        sampling_precision,
        sampling_recall,
        sampling_f1,
        macro_f1: 0.5 * (window_f1 + sampling_f1),
        false_negatives: window_counts.fn_ + sampling_counts.fn_,
        false_positives: window_counts.fp + sampling_counts.fp,
    })
}

fn predict(
    series: &[DhdlSeries],
    options: TiScheduleAdvisorOptions,
) -> Result<Prediction, Box<dyn Error>> {
    let advice = advise_ti_schedule(series, Some(options))?;
    let mut prediction = Prediction {
        window_issue: false,
        sampling_issue: false,
    };
    for suggestion in advice.suggestions() {
        match suggestion.kind() {
            TiSuggestionKind::InsertWindow => prediction.window_issue = true,
            TiSuggestionKind::ExtendSampling => prediction.sampling_issue = true,
            TiSuggestionKind::InsertWindowAndExtendSampling => {
                prediction.window_issue = true;
                prediction.sampling_issue = true;
            }
            TiSuggestionKind::NoChange => {}
        }
    }
    Ok(prediction)
}

#[derive(Debug, Default, Clone, Copy)]
struct Counts {
    tp: usize,
    fp: usize,
    fn_: usize,
}

impl Counts {
    fn add(&mut self, actual: bool, predicted: bool) {
        match (actual, predicted) {
            (true, true) => self.tp += 1,
            (false, true) => self.fp += 1,
            (true, false) => self.fn_ += 1,
            (false, false) => {}
        }
    }

    fn precision(self) -> f64 {
        let denominator = self.tp + self.fp;
        if denominator == 0 {
            0.0
        } else {
            self.tp as f64 / denominator as f64
        }
    }

    fn recall(self) -> f64 {
        let denominator = self.tp + self.fn_;
        if denominator == 0 {
            0.0
        } else {
            self.tp as f64 / denominator as f64
        }
    }

    fn f1(self) -> f64 {
        let precision = self.precision();
        let recall = self.recall();
        if precision + recall == 0.0 {
            0.0
        } else {
            2.0 * precision * recall / (precision + recall)
        }
    }
}

fn write_case_csv(path: &Path, cases: &[CalibrationCase]) -> Result<(), Box<dyn Error>> {
    let mut writer = csv::Writer::from_path(path)?;
    writer.write_record([
        "system",
        "replicate",
        "case_kind",
        "n_windows",
        "samples_per_window_min",
        "reference_delta_f",
        "case_delta_f",
        "abs_error",
        "true_window_issue",
        "true_sampling_issue",
    ])?;
    for case in cases {
        writer.write_record([
            case.system.as_str(),
            case.replicate.as_str(),
            case.case_kind.as_str(),
            &case.n_windows.to_string(),
            &case.samples_per_window_min.to_string(),
            &format!("{:.8}", case.reference_delta_f),
            &format!("{:.8}", case.case_delta_f),
            &format!("{:.8}", case.abs_error),
            &case.true_window_issue.to_string(),
            &case.true_sampling_issue.to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_sweep_csv(path: &Path, sweep: &[SweepResult]) -> Result<(), Box<dyn Error>> {
    let mut writer = csv::Writer::from_path(path)?;
    writer.write_record([
        "rank",
        "block_cv_min",
        "slope_z_min",
        "curvature_z_min",
        "interval_uncertainty_z_min",
        "window_precision",
        "window_recall",
        "window_f1",
        "sampling_precision",
        "sampling_recall",
        "sampling_f1",
        "macro_f1",
        "false_negatives",
        "false_positives",
    ])?;
    for (idx, result) in sweep.iter().enumerate() {
        writer.write_record([
            &(idx + 1).to_string(),
            &format!("{:.3}", result.options.block_cv_min),
            &format!("{:.2}", result.options.slope_z_min),
            &format!("{:.2}", result.options.curvature_z_min),
            &format!("{:.2}", result.options.interval_uncertainty_z_min),
            &format!("{:.4}", result.window_precision),
            &format!("{:.4}", result.window_recall),
            &format!("{:.4}", result.window_f1),
            &format!("{:.4}", result.sampling_precision),
            &format!("{:.4}", result.sampling_recall),
            &format!("{:.4}", result.sampling_f1),
            &format!("{:.4}", result.macro_f1),
            &result.false_negatives.to_string(),
            &result.false_positives.to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_report(
    path: &Path,
    config: &Config,
    replicate_dir_count: usize,
    replicates: &[Replicate],
    cases: &[CalibrationCase],
    sweep: &[SweepResult],
) -> Result<(), Box<dyn Error>> {
    let best = sweep.first().ok_or("empty sweep")?;
    let systems = reference_delta_f_by_system(replicates).len();
    let mean_abs_full_error = cases
        .iter()
        .filter(|case| case.case_kind == "full")
        .map(|case| case.abs_error)
        .sum::<f64>()
        / replicates.len() as f64;
    let parsed_paths = replicates
        .iter()
        .take(5)
        .map(|replicate| format!("- `{}`", replicate.path.display()))
        .collect::<Vec<_>>()
        .join("\n");

    let report = format!(
        "# TI Advisor Calibration Report\n\n\
         ## Dataset\n\n\
         - Root: `{}`\n\
         - Strict replicate directories found: {}\n\
         - Parsed replicates: {}\n\
         - Systems represented: {}\n\
         - Generated calibration cases: {}\n\
         - Temperature: {:.3} K\n\
         - Error tolerance: {:.6} reduced units\n\
         - Error tolerance: {:.3} kcal/mol\n\
         - Mean full-replicate absolute error against system reference: {:.6}\n\n\
         Example parsed replicate paths:\n\n{}\n\n\
         ## Best Threshold Set\n\n\
         - `block_cv_min`: {:.3}\n\
         - `slope_z_min`: {:.2}\n\
         - `curvature_z_min`: {:.2}\n\
         - `interval_uncertainty_z_min`: {:.2}\n\
         - Window precision/recall/F1: {:.3} / {:.3} / {:.3}\n\
         - Sampling precision/recall/F1: {:.3} / {:.3} / {:.3}\n\
         - Macro F1: {:.3}\n\
         - False negatives: {}\n\
         - False positives: {}\n\n\
         ## Interpretation\n\n\
         The labels are empirical, not absolute truth: full 15-window replicate means are averaged by system and used as the reference. \
         Dropped-window cases test whether schedule coarsening creates integration error. Truncated-trajectory cases test whether reduced sampling creates error. \
         This is appropriate for calibrating defaults to this AMBER flexible-TI protocol, but the recommended values should be validated on a held-out protocol before treating them as universal defaults.\n\n\
         See `cases.csv` for generated labels and `threshold-sweep.csv` for all threshold combinations.\n",
        config.dataset_root.display(),
        replicate_dir_count,
        replicates.len(),
        systems,
        cases.len(),
        config.temperature_k,
        config.error_tolerance_reduced,
        config.error_tolerance_kcal_mol,
        mean_abs_full_error,
        parsed_paths,
        best.options.block_cv_min,
        best.options.slope_z_min,
        best.options.curvature_z_min,
        best.options.interval_uncertainty_z_min,
        best.window_precision,
        best.window_recall,
        best.window_f1,
        best.sampling_precision,
        best.sampling_recall,
        best.sampling_f1,
        best.macro_f1,
        best.false_negatives,
        best.false_positives
    );
    fs::write(path, report)?;
    Ok(())
}
