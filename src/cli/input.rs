use std::path::PathBuf;

use alchemrs::{
    extract_dhdl, extract_u_nk, extract_u_nk_with_potential,
    decorrelate_dhdl, decorrelate_u_nk, decorrelate_u_nk_with_observable,
    detect_equilibration_dhdl, detect_equilibration_observable, detect_equilibration_u_nk,
    DecorrelationOptions, DhdlSeries, UNkMatrix, UNkSeriesMethod,
};

use crate::cli::UNkObservable;
use crate::CliResult;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnalysisSampleCounts {
    pub windows: usize,
    pub samples_in: usize,
    pub samples_after_burnin: usize,
    pub samples_kept: usize,
}

pub struct LoadedWindows {
    pub windows: Vec<UNkMatrix>,
    pub sample_counts: AnalysisSampleCounts,
}

pub struct LoadedDhdlSeries {
    pub series: Vec<DhdlSeries>,
    pub sample_counts: AnalysisSampleCounts,
}

#[derive(Debug, Clone, Copy)]
pub struct AnalysisInputOptions {
    pub temperature: f64,
    pub decorrelate: bool,
    pub remove_burnin: usize,
    pub auto_equilibrate: bool,
    pub fast: bool,
    pub conservative: bool,
    pub nskip: usize,
    pub u_nk_observable: Option<UNkObservable>,
}

impl AnalysisInputOptions {
    fn normalized_equilibration_flags(self) -> (bool, bool) {
        if self.auto_equilibrate {
            (true, false)
        } else {
            (self.fast, self.conservative)
        }
    }

    fn decorrelation_options(self) -> DecorrelationOptions {
        let (fast, conservative) = self.normalized_equilibration_flags();
        DecorrelationOptions {
            remove_burnin: false,
            fast,
            conservative,
            nskip: self.nskip,
            ..DecorrelationOptions::default()
        }
    }

    fn equilibration_options(self) -> DecorrelationOptions {
        let (fast, conservative) = self.normalized_equilibration_flags();
        DecorrelationOptions {
            remove_burnin: false,
            fast,
            conservative,
            nskip: self.nskip,
            ..DecorrelationOptions::default()
        }
    }

    pub fn effective_fast(self) -> bool {
        self.normalized_equilibration_flags().0
    }

    pub fn effective_conservative(self) -> bool {
        self.normalized_equilibration_flags().1
    }

    pub fn u_nk_observable_name(self) -> Option<&'static str> {
        self.u_nk_observable.map(UNkObservable::as_str)
    }
}

fn map_cli_dhdl_error(error: Box<dyn std::error::Error>) -> Box<dyn std::error::Error> {
    let message = error.to_string();
    if message.contains("scalar dH/dlambda parsing is unsupported") {
        return "multidimensional GROMACS dH/dlambda inputs are not supported by the CLI yet; use u_nk-based workflows or provide a one-dimensional schedule".into();
    }
    error
}

fn map_cli_u_nk_error(error: Box<dyn std::error::Error>) -> Box<dyn std::error::Error> {
    let message = error.to_string();
    if message == "DE u_nk preprocessing requires one-dimensional lambda states" {
        return "the CLI `de` observable is only supported for one-dimensional lambda schedules; use `--u-nk-observable all` or `epot`, or analyze a one-dimensional schedule".into();
    }
    error
}

pub fn load_windows(
    inputs: Vec<PathBuf>,
    options: AnalysisInputOptions,
) -> CliResult<LoadedWindows> {
    let mut windows = Vec::with_capacity(inputs.len());
    let mut samples_in = 0;
    let mut samples_after_burnin = 0;
    let mut samples_kept = 0;
    for path in inputs {
        if options.decorrelate || options.auto_equilibrate {
            let observable = options
                .u_nk_observable
                .ok_or("u_nk observable is required for u_nk estimators")?;
            let u_nk = match observable {
                UNkObservable::Epot => {
                    let (mut u_nk, potential) = extract_u_nk_with_potential(&path, options.temperature)
                        .map_err(|err| map_cli_u_nk_error(err.into()))?;
                    let mut potential = potential;
                    samples_in += u_nk.n_samples();
                    u_nk = trim_u_nk(u_nk, options.remove_burnin)?;
                    potential = trim_values(potential, options.remove_burnin)?;
                    if options.auto_equilibrate {
                        let equilibration = detect_equilibration_observable(
                            &potential,
                            options.effective_fast(),
                            options.nskip,
                        )
                        .map_err(|err| map_cli_u_nk_error(err.into()))?;
                        u_nk = trim_u_nk(u_nk, equilibration.t0)?;
                        potential = trim_values(potential, equilibration.t0)?;
                    }
                    samples_after_burnin += u_nk.n_samples();
                    if options.decorrelate {
                        u_nk = decorrelate_u_nk_with_observable(
                            &u_nk,
                            &potential,
                            &options.decorrelation_options(),
                        )
                        .map_err(|err| map_cli_u_nk_error(err.into()))?;
                    }
                    u_nk
                }
                UNkObservable::De | UNkObservable::All => {
                    let method = match observable {
                        UNkObservable::De => UNkSeriesMethod::DE,
                        UNkObservable::All => UNkSeriesMethod::All,
                        UNkObservable::Epot => unreachable!(),
                    };
                    let mut u_nk = extract_u_nk(&path, options.temperature)
                        .map_err(|err| map_cli_u_nk_error(err.into()))?;
                    samples_in += u_nk.n_samples();
                    u_nk = trim_u_nk(u_nk, options.remove_burnin)?;
                    if options.auto_equilibrate {
                        let equilibration = detect_equilibration_u_nk(
                            &u_nk,
                            method,
                            &options.equilibration_options(),
                        )
                        .map_err(|err| map_cli_u_nk_error(err.into()))?;
                        u_nk = trim_u_nk(u_nk, equilibration.t0)?;
                    }
                    samples_after_burnin += u_nk.n_samples();
                    if options.decorrelate {
                        u_nk = decorrelate_u_nk(&u_nk, method, &options.decorrelation_options())
                            .map_err(|err| map_cli_u_nk_error(err.into()))?;
                    }
                    u_nk
                }
            };
            samples_kept += u_nk.n_samples();
            windows.push(u_nk);
            continue;
        }
        let mut u_nk = extract_u_nk(path, options.temperature)
            .map_err(|err| map_cli_u_nk_error(err.into()))?;
        samples_in += u_nk.n_samples();
        u_nk = trim_u_nk(u_nk, options.remove_burnin)?;
        samples_after_burnin += u_nk.n_samples();
        samples_kept += u_nk.n_samples();
        windows.push(u_nk);
    }
    Ok(LoadedWindows {
        sample_counts: AnalysisSampleCounts {
            windows: windows.len(),
            samples_in,
            samples_after_burnin,
            samples_kept,
        },
        windows,
    })
}

pub fn load_dhdl_series(
    inputs: Vec<PathBuf>,
    options: AnalysisInputOptions,
) -> CliResult<LoadedDhdlSeries> {
    let mut series = Vec::with_capacity(inputs.len());
    let mut samples_in = 0;
    let mut samples_after_burnin = 0;
    let mut samples_kept = 0;
    for path in inputs {
        let mut dhdl = extract_dhdl(path, options.temperature)
            .map_err(|err| map_cli_dhdl_error(err.into()))?;
        samples_in += dhdl.values().len();
        dhdl = trim_dhdl(dhdl, options.remove_burnin)?;
        if options.auto_equilibrate {
            let equilibration = detect_equilibration_dhdl(&dhdl, &options.equilibration_options())?;
            dhdl = trim_dhdl(dhdl, equilibration.t0)?;
        }
        samples_after_burnin += dhdl.values().len();
        if options.decorrelate {
            dhdl = decorrelate_dhdl(&dhdl, &options.decorrelation_options())?;
        }
        samples_kept += dhdl.values().len();
        series.push(dhdl);
    }
    Ok(LoadedDhdlSeries {
        sample_counts: AnalysisSampleCounts {
            windows: series.len(),
            samples_in,
            samples_after_burnin,
            samples_kept,
        },
        series,
    })
}

fn trim_u_nk(u_nk: UNkMatrix, remove_burnin: usize) -> CliResult<UNkMatrix> {
    if remove_burnin == 0 {
        return Ok(u_nk);
    }

    let n_samples = u_nk.n_samples();
    if remove_burnin >= n_samples {
        return Err("remove-burnin exceeds series length".into());
    }

    let n_states = u_nk.n_states();
    let data = u_nk.data()[(remove_burnin * n_states)..].to_vec();
    let time = u_nk.time_ps()[remove_burnin..].to_vec();
    Ok(UNkMatrix::new_with_labels(
        n_samples - remove_burnin,
        n_states,
        data,
        time,
        u_nk.sampled_state().cloned(),
        u_nk.evaluated_states().to_vec(),
        u_nk.lambda_labels().map(|labels| labels.to_vec()),
    )?)
}

fn trim_dhdl(dhdl: DhdlSeries, remove_burnin: usize) -> CliResult<DhdlSeries> {
    if remove_burnin == 0 {
        return Ok(dhdl);
    }

    let values = dhdl.values();
    if remove_burnin >= values.len() {
        return Err("remove-burnin exceeds series length".into());
    }

    let new_time = dhdl.time_ps()[remove_burnin..].to_vec();
    let new_values = values[remove_burnin..].to_vec();
    Ok(DhdlSeries::new(dhdl.state().clone(), new_time, new_values)?)
}

fn trim_values(values: Vec<f64>, remove_burnin: usize) -> CliResult<Vec<f64>> {
    if remove_burnin == 0 {
        return Ok(values);
    }

    if remove_burnin >= values.len() {
        return Err("remove-burnin exceeds series length".into());
    }

    Ok(values[remove_burnin..].to_vec())
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use crate::cli::UNkObservable;

    use super::{load_dhdl_series, load_windows, AnalysisInputOptions};

    #[test]
    fn load_windows_accepts_multidimensional_gromacs_states_for_cli_estimators() {
        let base = env!("CARGO_MANIFEST_DIR");
        let path = PathBuf::from(format!("{base}/fixtures/gromacs/lambda_15.xvg"));
        let loaded = load_windows(
            vec![path],
            AnalysisInputOptions {
                temperature: 298.0,
                decorrelate: false,
                remove_burnin: 0,
                auto_equilibrate: false,
                fast: false,
                conservative: true,
                nskip: 1,
                u_nk_observable: Some(UNkObservable::All),
            },
        )
        .expect("expected multidimensional CLI u_nk load to succeed");

        assert_eq!(loaded.windows.len(), 1);
        assert_eq!(
            loaded.windows[0].sampled_state().unwrap().lambdas(),
            &[0.0, 0.0, 0.8, 0.0, 0.0]
        );
        assert_eq!(
            loaded.windows[0].lambda_labels().unwrap(),
            &[
                "mass-lambda",
                "coul-lambda",
                "vdw-lambda",
                "bonded-lambda",
                "restraint-lambda",
            ]
        );
    }

    #[test]
    fn load_dhdl_series_reports_multidimensional_gromacs_ti_limit() {
        let base = env!("CARGO_MANIFEST_DIR");
        let path = PathBuf::from(format!("{base}/fixtures/gromacs/lambda_15.xvg"));
        let err = match load_dhdl_series(
            vec![path],
            AnalysisInputOptions {
                temperature: 298.0,
                decorrelate: false,
                remove_burnin: 0,
                auto_equilibrate: false,
                fast: false,
                conservative: true,
                nskip: 1,
                u_nk_observable: None,
            },
        ) {
            Ok(_) => panic!("expected multidimensional CLI dH/dlambda load to fail"),
            Err(err) => err,
        };

        assert_eq!(
            err.to_string(),
            "multidimensional GROMACS dH/dlambda inputs are not supported by the CLI yet; use u_nk-based workflows or provide a one-dimensional schedule"
        );
    }

    #[test]
    fn auto_equilibrate_overrides_decorrelation_flags() {
        let decorrelation = AnalysisInputOptions {
            temperature: 300.0,
            decorrelate: true,
            remove_burnin: 0,
            auto_equilibrate: true,
            fast: false,
            conservative: true,
            nskip: 7,
            u_nk_observable: Some(UNkObservable::De),
        }
        .decorrelation_options();

        assert!(!decorrelation.remove_burnin);
        assert!(decorrelation.fast);
        assert!(!decorrelation.conservative);
        assert_eq!(decorrelation.nskip, 7);
    }

    #[test]
    fn auto_equilibrate_uses_effective_flags_for_equilibration() {
        let equilibration = AnalysisInputOptions {
            temperature: 300.0,
            decorrelate: false,
            remove_burnin: 0,
            auto_equilibrate: true,
            fast: false,
            conservative: true,
            nskip: 7,
            u_nk_observable: Some(UNkObservable::De),
        }
        .equilibration_options();

        assert!(!equilibration.remove_burnin);
        assert!(equilibration.fast);
        assert!(!equilibration.conservative);
        assert_eq!(equilibration.nskip, 7);
    }
}
