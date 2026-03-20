use std::path::PathBuf;

use alchemrs_core::{DhdlSeries, UNkMatrix};
use alchemrs_parse::amber::{extract_dhdl, extract_u_nk, extract_u_nk_with_potential};
use alchemrs_prep::{decorrelate_dhdl, decorrelate_u_nk_with_observable, DecorrelationOptions};

use crate::CliResult;

#[derive(Debug, Clone, Copy)]
pub struct AnalysisInputOptions {
    pub temperature: f64,
    pub decorrelate: bool,
    pub remove_burnin: usize,
    pub auto_equilibrate: bool,
    pub fast: bool,
    pub conservative: bool,
    pub nskip: usize,
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
            remove_burnin: self.auto_equilibrate,
            fast,
            conservative,
            nskip: self.nskip,
            ..DecorrelationOptions::default()
        }
    }
}

pub fn load_windows(
    inputs: Vec<PathBuf>,
    options: AnalysisInputOptions,
) -> CliResult<Vec<UNkMatrix>> {
    let mut windows = Vec::with_capacity(inputs.len());
    for path in inputs {
        if options.decorrelate {
            let (mut u_nk, potential) = extract_u_nk_with_potential(&path, options.temperature)?;
            u_nk = trim_u_nk(u_nk, options.remove_burnin)?;
            let potential = trim_values(potential, options.remove_burnin)?;
            u_nk = decorrelate_u_nk_with_observable(
                &u_nk,
                &potential,
                &options.decorrelation_options(),
            )?;
            windows.push(u_nk);
            continue;
        }
        let mut u_nk = extract_u_nk(path, options.temperature)?;
        u_nk = trim_u_nk(u_nk, options.remove_burnin)?;
        windows.push(u_nk);
    }
    Ok(windows)
}

pub fn load_dhdl_series(
    inputs: Vec<PathBuf>,
    options: AnalysisInputOptions,
) -> CliResult<Vec<DhdlSeries>> {
    let mut series = Vec::with_capacity(inputs.len());
    for path in inputs {
        let mut dhdl = extract_dhdl(path, options.temperature)?;
        dhdl = trim_dhdl(dhdl, options.remove_burnin)?;
        if options.decorrelate {
            dhdl = decorrelate_dhdl(&dhdl, &options.decorrelation_options())?;
        }
        series.push(dhdl);
    }
    Ok(series)
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
    Ok(UNkMatrix::new(
        n_samples - remove_burnin,
        n_states,
        data,
        time,
        u_nk.sampled_state().cloned(),
        u_nk.evaluated_states().to_vec(),
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
