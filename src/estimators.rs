use crate::core::{CoreError, DeltaFMatrix, DhdlSeries, FreeEnergyEstimate, Result};
use crate::core::{StatePoint, UNkMatrix};

type CombinedWindows = (Vec<Vec<f64>>, Vec<f64>, Vec<StatePoint>);
type PairEstimate = (f64, f64);
type ExpRow = (usize, Vec<f64>, Vec<f64>);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrationMethod {
    Trapezoidal,
    Simpson,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TiOptions {
    pub method: IntegrationMethod,
    pub parallel: bool,
}

impl Default for TiOptions {
    fn default() -> Self {
        Self {
            method: IntegrationMethod::Trapezoidal,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct TiEstimator {
    pub options: TiOptions,
}

impl TiEstimator {
    pub fn new(options: TiOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, series: &[DhdlSeries]) -> Result<FreeEnergyEstimate> {
        if series.len() < 2 {
            return Err(CoreError::InvalidShape {
                expected: 2,
                found: series.len(),
            });
        }

        let mut points: Vec<(f64, f64, f64, StatePoint)> = if self.options.parallel {
            use rayon::prelude::*;
            series
                .par_iter()
                .map(|item| {
                    let lambda = extract_lambda(item.state())?;
                    let mean = mean_values(item.values())?;
                    let sem2 = sem2_values(item.values())?;
                    Ok((lambda, mean, sem2, item.state().clone()))
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            let mut out = Vec::with_capacity(series.len());
            for item in series {
                let lambda = extract_lambda(item.state())?;
                let mean = mean_values(item.values())?;
                let sem2 = sem2_values(item.values())?;
                out.push((lambda, mean, sem2, item.state().clone()));
            }
            out
        };
        points.sort_by(|a, b| a.0.total_cmp(&b.0));

        let lambdas: Vec<f64> = points.iter().map(|(l, _, _, _)| *l).collect();
        let values: Vec<f64> = points.iter().map(|(_, v, _, _)| *v).collect();
        let sem2_values: Vec<f64> = points.iter().map(|(_, _, s, _)| *s).collect();

        let delta_f = match self.options.method {
            IntegrationMethod::Trapezoidal => integrate_trapezoidal(&lambdas, &values)?,
            IntegrationMethod::Simpson => integrate_simpson(&lambdas, &values)?,
        };

        let uncertainty = match self.options.method {
            IntegrationMethod::Trapezoidal => {
                Some(trapezoidal_uncertainty(&lambdas, &sem2_values)?)
            }
            IntegrationMethod::Simpson => None,
        };

        let from_state = points.first().unwrap().3.clone();
        let to_state = points.last().unwrap().3.clone();
        FreeEnergyEstimate::new(delta_f, uncertainty, from_state, to_state)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarMethod {
    FalsePosition,
    SelfConsistentIteration,
    Bisection,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarUncertainty {
    Bar,
}

#[derive(Debug, Clone)]
pub struct BarOptions {
    pub maximum_iterations: usize,
    pub relative_tolerance: f64,
    pub method: BarMethod,
    pub uncertainty: BarUncertainty,
    pub parallel: bool,
}

impl Default for BarOptions {
    fn default() -> Self {
        Self {
            maximum_iterations: 10_000,
            relative_tolerance: 1.0e-7,
            method: BarMethod::FalsePosition,
            uncertainty: BarUncertainty::Bar,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct BarEstimator {
    pub options: BarOptions,
}

impl BarEstimator {
    pub fn new(options: BarOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        if windows.len() < 2 {
            return Err(CoreError::InvalidShape {
                expected: 2,
                found: windows.len(),
            });
        }
        let eval_states = ensure_consistent_states(windows)?;
        let lambdas: Vec<f64> = eval_states.iter().map(|s| s.lambdas()[0]).collect();
        let mut window_by_index: Vec<Option<&UNkMatrix>> = vec![None; lambdas.len()];
        for window in windows {
            let sampled = window.sampled_state().ok_or_else(|| {
                CoreError::InvalidState("sampled_state required for BAR".to_string())
            })?;
            let lambda = sampled.lambdas()[0];
            let idx = find_lambda_index(&lambdas, lambda)?;
            if window_by_index[idx].is_some() {
                return Err(CoreError::InvalidState(format!(
                    "multiple windows for sampled_state lambda {lambda}"
                )));
            }
            window_by_index[idx] = Some(window);
        }

        let pair_results: Vec<Result<PairEstimate>> = if self.options.parallel {
            use rayon::prelude::*;
            (0..(lambdas.len() - 1))
                .into_par_iter()
                .map(|idx| {
                    let lambda = lambdas[idx];
                    let lambda_next = lambdas[idx + 1];
                    let win_f = window_by_index[idx].ok_or_else(|| {
                        CoreError::InvalidState(format!("missing window for lambda {lambda}"))
                    })?;
                    let win_r = window_by_index[idx + 1].ok_or_else(|| {
                        CoreError::InvalidState(format!("missing window for lambda {lambda_next}"))
                    })?;
                    let w_f = work_values(win_f, idx, idx + 1)?;
                    let w_r = work_values(win_r, idx, idx + 1)?
                        .into_iter()
                        .map(|w| -w)
                        .collect::<Vec<_>>();
                    let (df, ddf) = bar_estimate(
                        &w_f,
                        &w_r,
                        self.options.method,
                        self.options.maximum_iterations,
                        self.options.relative_tolerance,
                    )?;
                    Ok((df, ddf * ddf))
                })
                .collect::<Vec<_>>()
        } else {
            let mut out = Vec::new();
            for idx in 0..(lambdas.len() - 1) {
                let lambda = lambdas[idx];
                let lambda_next = lambdas[idx + 1];
                let win_f = window_by_index[idx].ok_or_else(|| {
                    CoreError::InvalidState(format!("missing window for lambda {lambda}"))
                })?;
                let win_r = window_by_index[idx + 1].ok_or_else(|| {
                    CoreError::InvalidState(format!("missing window for lambda {lambda_next}"))
                })?;
                let w_f = work_values(win_f, idx, idx + 1)?;
                let w_r = work_values(win_r, idx, idx + 1)?
                    .into_iter()
                    .map(|w| -w)
                    .collect::<Vec<_>>();
                let (df, ddf) = bar_estimate(
                    &w_f,
                    &w_r,
                    self.options.method,
                    self.options.maximum_iterations,
                    self.options.relative_tolerance,
                )?;
                out.push(Ok((df, ddf * ddf)));
            }
            out
        };

        let mut deltas = Vec::with_capacity(pair_results.len());
        let mut d_deltas = Vec::with_capacity(pair_results.len());
        for result in pair_results {
            let (delta, d_delta) = result?;
            deltas.push(delta);
            d_deltas.push(d_delta);
        }

        let n_states = lambdas.len();
        let mut adelta = vec![0.0; n_states * n_states];
        let mut ad_delta = vec![f64::NAN; n_states * n_states];

        for j in 0..deltas.len() {
            let mut out = Vec::new();
            let mut dout = Vec::new();
            for i in 0..(deltas.len() - j) {
                out.push(deltas[i..=i + j].iter().sum::<f64>());
                if j == 0 {
                    dout.push(d_deltas[i..=i + j].iter().sum::<f64>());
                } else {
                    dout.push(f64::NAN);
                }
            }
            for (i, value) in out.into_iter().enumerate() {
                let row = i;
                let col = i + j + 1;
                adelta[row * n_states + col] = value;
            }
            for (i, value) in dout.into_iter().enumerate() {
                let row = i;
                let col = i + j + 1;
                ad_delta[row * n_states + col] = value;
            }
        }

        for i in 0..n_states {
            for j in (i + 1)..n_states {
                let val = adelta[i * n_states + j];
                adelta[j * n_states + i] = -val;
                let unc = ad_delta[i * n_states + j];
                ad_delta[j * n_states + i] = unc;
            }
        }

        DeltaFMatrix::new(adelta, Some(ad_delta), n_states, eval_states)
    }
}

#[derive(Debug, Clone)]
pub struct MbarOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub initial_f_k: Option<Vec<f64>>,
    pub compute_uncertainty: bool,
    pub parallel: bool,
}

impl Default for MbarOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
            tolerance: 1.0e-7,
            initial_f_k: None,
            compute_uncertainty: true,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MbarEstimator {
    pub options: MbarOptions,
}

impl MbarEstimator {
    pub fn new(options: MbarOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        if windows.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        let (u_kn, n_k, states) = combine_windows(windows)?;
        let mut f_k = initial_f_k(&self.options, states.len())?;

        mbar_solve(
            &u_kn,
            &n_k,
            &mut f_k,
            self.options.tolerance,
            self.options.max_iterations,
        )?;

        let n_states = states.len();
        let mut values = vec![0.0; n_states * n_states];
        for i in 0..n_states {
            for j in 0..n_states {
                values[i * n_states + j] = f_k[j] - f_k[i];
            }
        }
        let uncertainties = if self.options.compute_uncertainty {
            Some(mbar_uncertainty(&u_kn, &n_k, &f_k)?)
        } else {
            None
        };
        DeltaFMatrix::new(values, uncertainties, n_states, states)
    }
}

pub fn mbar_log_weights_from_windows(
    windows: &[UNkMatrix],
    options: &MbarOptions,
) -> Result<(Vec<f64>, Vec<f64>, Vec<StatePoint>)> {
    if windows.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let (u_kn, n_k, states) = combine_windows(windows)?;
    let mut f_k = initial_f_k(options, states.len())?;
    mbar_solve(
        &u_kn,
        &n_k,
        &mut f_k,
        options.tolerance,
        options.max_iterations,
    )?;
    let log_weights = mbar_log_weights(&u_kn, &n_k, &f_k)?;
    Ok((log_weights, n_k, states))
}

fn initial_f_k(options: &MbarOptions, n_states: usize) -> Result<Vec<f64>> {
    if let Some(ref initial) = options.initial_f_k {
        if initial.len() != n_states {
            return Err(CoreError::InvalidShape {
                expected: n_states,
                found: initial.len(),
            });
        }
        Ok(initial.clone())
    } else {
        Ok(vec![0.0; n_states])
    }
}

fn ensure_consistent_states(windows: &[UNkMatrix]) -> Result<Vec<StatePoint>> {
    let first = windows[0].evaluated_states();
    let first_temperature = first[0].temperature_k();
    for state in first {
        if state.lambdas().len() != 1 {
            return Err(CoreError::Unsupported(
                "estimators require one-dimensional lambda states".to_string(),
            ));
        }
        if (state.temperature_k() - first_temperature).abs() > 1e-6 {
            return Err(CoreError::InvalidState(
                "evaluated_states temperatures differ within a window".to_string(),
            ));
        }
    }
    for window in windows.iter().skip(1) {
        let states = window.evaluated_states();
        if states.len() != first.len() {
            return Err(CoreError::InvalidShape {
                expected: first.len(),
                found: states.len(),
            });
        }
        for (a, b) in first.iter().zip(states.iter()) {
            if a.lambdas().len() != b.lambdas().len() {
                return Err(CoreError::InvalidShape {
                    expected: a.lambdas().len(),
                    found: b.lambdas().len(),
                });
            }
            if a.lambdas().len() != 1 {
                return Err(CoreError::Unsupported(
                    "estimators require one-dimensional lambda states".to_string(),
                ));
            }
            if (a.temperature_k() - b.temperature_k()).abs() > 1e-6 {
                return Err(CoreError::InvalidState(
                    "evaluated_states temperatures differ between windows".to_string(),
                ));
            }
            if (a.lambdas()[0] - b.lambdas()[0]).abs() > 1e-6 {
                return Err(CoreError::InvalidState(
                    "evaluated_states differ between windows".to_string(),
                ));
            }
        }
    }
    for window in windows {
        let sampled = window.sampled_state().ok_or_else(|| {
            CoreError::InvalidState("sampled_state required for estimator".to_string())
        })?;
        if sampled.lambdas().len() != 1 {
            return Err(CoreError::Unsupported(
                "estimators require one-dimensional lambda states".to_string(),
            ));
        }
        if (sampled.temperature_k() - first_temperature).abs() > 1e-6 {
            return Err(CoreError::InvalidState(
                "sampled_state temperature differs from evaluated_states".to_string(),
            ));
        }
        if !first
            .iter()
            .any(|state| (state.lambdas()[0] - sampled.lambdas()[0]).abs() < 1e-6)
        {
            return Err(CoreError::InvalidState(
                "sampled_state not found in evaluated_states".to_string(),
            ));
        }
    }
    Ok(first.to_vec())
}

fn find_lambda_index(lambdas: &[f64], target: f64) -> Result<usize> {
    for (idx, value) in lambdas.iter().enumerate() {
        if (*value - target).abs() < 1e-6 {
            return Ok(idx);
        }
    }
    Err(CoreError::InvalidState(
        "sampled_state not found in evaluated_states".to_string(),
    ))
}

fn work_values(window: &UNkMatrix, idx_a: usize, idx_b: usize) -> Result<Vec<f64>> {
    let n_states = window.n_states();
    let n_samples = window.n_samples();
    if idx_b >= n_states {
        return Err(CoreError::InvalidShape {
            expected: n_states,
            found: idx_b + 1,
        });
    }
    let data = window.data();
    let mut out = Vec::with_capacity(n_samples);
    for sample_idx in 0..n_samples {
        let offset = sample_idx * n_states;
        let a = data[offset + idx_a];
        let b = data[offset + idx_b];
        if a.is_nan() || b.is_nan() {
            return Err(CoreError::NonFiniteValue(format!(
                "u_nk work input for states ({idx_a}, {idx_b}) at sample {sample_idx} must not be NaN"
            )));
        }
        let delta = b - a;
        if delta.is_nan() {
            return Err(CoreError::NonFiniteValue(format!(
                "work value for states ({idx_a}, {idx_b}) at sample {sample_idx} must not be NaN"
            )));
        }
        out.push(delta);
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct ExpOptions {
    pub compute_uncertainty: bool,
    pub parallel: bool,
}

impl Default for ExpOptions {
    fn default() -> Self {
        Self {
            compute_uncertainty: true,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ExpEstimator {
    pub options: ExpOptions,
}

impl ExpEstimator {
    pub fn new(options: ExpOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, windows: &[UNkMatrix]) -> Result<DeltaFMatrix> {
        if windows.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        let states = ensure_consistent_states(windows)?;
        let n_states = states.len();
        let lambdas: Vec<f64> = states.iter().map(|s| s.lambdas()[0]).collect();
        let mut window_map: Vec<Option<&UNkMatrix>> = vec![None; n_states];

        for window in windows {
            let sampled = window.sampled_state().ok_or_else(|| {
                CoreError::InvalidState("sampled_state required for EXP".to_string())
            })?;
            let idx = find_lambda_index(&lambdas, sampled.lambdas()[0])?;
            if window_map[idx].is_some() {
                return Err(CoreError::InvalidState(
                    "multiple windows for same sampled_state".to_string(),
                ));
            }
            window_map[idx] = Some(window);
        }

        for (idx, entry) in window_map.iter().enumerate() {
            if entry.is_none() {
                return Err(CoreError::InvalidState(format!(
                    "missing window for state {idx}"
                )));
            }
        }

        let mut values = vec![0.0; n_states * n_states];
        let mut uncertainties = if self.options.compute_uncertainty {
            Some(vec![0.0; n_states * n_states])
        } else {
            None
        };

        if self.options.parallel {
            use rayon::prelude::*;
            let rows: Vec<Result<ExpRow>> = (0..n_states)
                .into_par_iter()
                .map(|i| {
                    let window = window_map[i].expect("window present");
                    let mut row = Vec::with_capacity(n_states);
                    let mut row_unc = Vec::with_capacity(n_states);
                    for j in 0..n_states {
                        let work = work_values(window, i, j)?;
                        let delta = exp_delta_f(&work)?;
                        row.push(delta);
                        if self.options.compute_uncertainty {
                            row_unc.push(exp_uncertainty(&work)?);
                        }
                    }
                    Ok((i, row, row_unc))
                })
                .collect::<Vec<_>>();
            for entry in rows {
                let (i, row, row_unc) = entry?;
                for (j, value) in row.into_iter().enumerate() {
                    values[i * n_states + j] = value;
                }
                if let Some(ref mut unc) = uncertainties {
                    for (j, value) in row_unc.into_iter().enumerate() {
                        unc[i * n_states + j] = value;
                    }
                }
            }
        } else {
            for i in 0..n_states {
                let window = window_map[i].expect("window present");
                for j in 0..n_states {
                    let work = work_values(window, i, j)?;
                    let delta = exp_delta_f(&work)?;
                    values[i * n_states + j] = delta;
                    if let Some(ref mut unc) = uncertainties {
                        unc[i * n_states + j] = exp_uncertainty(&work)?;
                    }
                }
            }
        }

        DeltaFMatrix::new(values, uncertainties, n_states, states)
    }
}

fn combine_windows(windows: &[UNkMatrix]) -> Result<CombinedWindows> {
    let states = ensure_consistent_states(windows)?;
    let n_states = states.len();
    let mut n_k = vec![0.0; n_states];
    let mut offsets = Vec::with_capacity(windows.len());
    let mut total_samples = 0usize;

    let lambdas: Vec<f64> = states.iter().map(|s| s.lambdas()[0]).collect();

    for window in windows {
        let sampled = window.sampled_state().ok_or_else(|| {
            CoreError::InvalidState("sampled_state required for MBAR".to_string())
        })?;
        let lambda = sampled.lambdas()[0];
        let idx = find_lambda_index(&lambdas, lambda)?;
        n_k[idx] += window.n_samples() as f64;
        offsets.push(total_samples);
        total_samples += window.n_samples();
    }

    let mut u_kn = vec![vec![0.0; total_samples]; n_states];
    for (win_idx, window) in windows.iter().enumerate() {
        let offset = offsets[win_idx];
        let data = window.data();
        for sample_idx in 0..window.n_samples() {
            let row_offset = sample_idx * n_states;
            let out_idx = offset + sample_idx;
            for k in 0..n_states {
                u_kn[k][out_idx] = data[row_offset + k];
            }
        }
    }

    validate_mbar_input(&u_kn)?;

    Ok((u_kn, n_k, states))
}

fn mbar_solve(
    u_kn: &[Vec<f64>],
    n_k: &[f64],
    f_k: &mut [f64],
    tolerance: f64,
    max_iterations: usize,
) -> Result<()> {
    let n_states = u_kn.len();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let n_samples = u_kn[0].len();
    let ln_n_k: Vec<f64> = n_k
        .iter()
        .map(|n| if *n > 0.0 { n.ln() } else { f64::NEG_INFINITY })
        .collect();

    let mut log_denominator = vec![0.0; n_samples];
    for _ in 0..max_iterations {
        for n in 0..n_samples {
            let mut max_arg = f64::NEG_INFINITY;
            for k in 0..n_states {
                if n_k[k] == 0.0 {
                    continue;
                }
                let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
                if val > max_arg {
                    max_arg = val;
                }
            }
            let mut sum = 0.0;
            for k in 0..n_states {
                if n_k[k] == 0.0 {
                    continue;
                }
                let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
                sum += (val - max_arg).exp();
            }
            log_denominator[n] = max_arg + sum.ln();
            if !log_denominator[n].is_finite() {
                return Err(CoreError::NonFiniteValue(format!(
                    "MBAR denominator became non-finite at sample {n}"
                )));
            }
        }

        let mut f_new = vec![0.0; n_states];
        for k in 0..n_states {
            let mut max_arg = f64::NEG_INFINITY;
            for n in 0..n_samples {
                let val = -log_denominator[n] - u_kn[k][n];
                if val > max_arg {
                    max_arg = val;
                }
            }
            let mut sum = 0.0;
            for n in 0..n_samples {
                let val = -log_denominator[n] - u_kn[k][n];
                sum += (val - max_arg).exp();
            }
            f_new[k] = -(max_arg + sum.ln());
            if !f_new[k].is_finite() {
                return Err(CoreError::NonFiniteValue(format!(
                    "MBAR free energy became non-finite for state {k}"
                )));
            }
        }
        let shift = f_new[0];
        for value in f_new.iter_mut().take(n_states) {
            *value -= shift;
        }

        let mut max_delta = 0.0;
        for k in 0..n_states {
            let delta = (f_new[k] - f_k[k]).abs();
            if delta > max_delta {
                max_delta = delta;
            }
        }
        f_k.copy_from_slice(&f_new);
        if max_delta < tolerance {
            return Ok(());
        }
    }

    Err(CoreError::ConvergenceFailure)
}

fn mbar_uncertainty(u_kn: &[Vec<f64>], n_k: &[f64], f_k: &[f64]) -> Result<Vec<f64>> {
    use nalgebra::{DMatrix, SymmetricEigen};

    let n_states = u_kn.len();
    let n_samples = u_kn[0].len();

    let log_w = mbar_log_weights(u_kn, n_k, f_k)?;
    let mut w_data = Vec::with_capacity(n_samples * n_states);
    for n in 0..n_samples {
        for k in 0..n_states {
            w_data.push(log_w[n * n_states + k].exp());
        }
    }
    let w = DMatrix::from_row_slice(n_samples, n_states, &w_data);

    let wtw = &w.transpose() * &w;
    let eigen = SymmetricEigen::new(wtw);
    let mut s2 = eigen.eigenvalues;
    for value in s2.iter_mut() {
        if *value < 0.0 {
            *value = 0.0;
        }
    }
    let sigma = DMatrix::from_diagonal(&s2.map(|v| v.sqrt()));
    let v = eigen.eigenvectors;

    let mut ndiag = DMatrix::zeros(n_states, n_states);
    for i in 0..n_states {
        ndiag[(i, i)] = n_k[i];
    }

    let identity = DMatrix::identity(n_states, n_states);
    let a = &identity - &sigma * &v.transpose() * &ndiag * &v * &sigma;
    let a_inv = pseudoinverse(&a, 1.0e-10)?;
    let theta = &v * &sigma * a_inv * &sigma * v.transpose();

    let mut uncertainties = vec![0.0; n_states * n_states];
    for i in 0..n_states {
        for j in 0..n_states {
            let val = theta[(i, i)] + theta[(j, j)] - 2.0 * theta[(i, j)];
            uncertainties[i * n_states + j] = if val < 0.0 && val > -1e-10 {
                0.0
            } else if val < 0.0 {
                val.abs().sqrt()
            } else {
                val.sqrt()
            };
        }
    }
    Ok(uncertainties)
}

fn mbar_log_weights(u_kn: &[Vec<f64>], n_k: &[f64], f_k: &[f64]) -> Result<Vec<f64>> {
    let n_states = u_kn.len();
    let n_samples = u_kn[0].len();
    let ln_n_k: Vec<f64> = n_k
        .iter()
        .map(|n| if *n > 0.0 { n.ln() } else { f64::NEG_INFINITY })
        .collect();

    let mut log_w = vec![0.0; n_samples * n_states];
    for n in 0..n_samples {
        let mut max_arg = f64::NEG_INFINITY;
        for k in 0..n_states {
            if n_k[k] == 0.0 {
                continue;
            }
            let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
            if val > max_arg {
                max_arg = val;
            }
        }
        let mut sum = 0.0;
        for k in 0..n_states {
            if n_k[k] == 0.0 {
                continue;
            }
            let val = f_k[k] - u_kn[k][n] + ln_n_k[k];
            sum += (val - max_arg).exp();
        }
        let log_denom = max_arg + sum.ln();
        for k in 0..n_states {
            log_w[n * n_states + k] = f_k[k] - u_kn[k][n] - log_denom;
        }
    }
    Ok(log_w)
}

fn pseudoinverse(matrix: &nalgebra::DMatrix<f64>, tol: f64) -> Result<nalgebra::DMatrix<f64>> {
    let svd = matrix.clone().svd(true, true);
    let u = svd.u.ok_or(CoreError::ConvergenceFailure)?;
    let v_t = svd.v_t.ok_or(CoreError::ConvergenceFailure)?;
    let mut sigma_inv = nalgebra::DMatrix::zeros(matrix.nrows(), matrix.ncols());
    for i in 0..svd.singular_values.len() {
        let value = svd.singular_values[i];
        if value > tol {
            sigma_inv[(i, i)] = 1.0 / value;
        }
    }
    Ok(v_t.transpose() * sigma_inv * u.transpose())
}

fn logsumexp(values: &[f64]) -> f64 {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    let sum = values.iter().map(|v| (v - max).exp()).sum::<f64>();
    max + sum.ln()
}

fn exp_delta_f(work: &[f64]) -> Result<f64> {
    if work
        .iter()
        .any(|value| *value == f64::NEG_INFINITY || value.is_nan())
    {
        return Err(CoreError::NonFiniteValue(
            "EXP work values must be finite or +inf".to_string(),
        ));
    }
    let t = work.len() as f64;
    let neg: Vec<f64> = work.iter().map(|w| -w).collect();
    let delta_f = -(logsumexp(&neg) - t.ln());
    if !delta_f.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "EXP delta_f became non-finite; no finite Boltzmann weight remained".to_string(),
        ));
    }
    Ok(delta_f)
}

fn exp_uncertainty(work: &[f64]) -> Result<f64> {
    if work.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    if work
        .iter()
        .any(|value| *value == f64::NEG_INFINITY || value.is_nan())
    {
        return Err(CoreError::NonFiniteValue(
            "EXP work values must be finite or +inf".to_string(),
        ));
    }
    let t = work.len() as f64;
    let mut max_arg = f64::NEG_INFINITY;
    for value in work {
        let arg = -value;
        if arg > max_arg {
            max_arg = arg;
        }
    }
    let mut x = Vec::with_capacity(work.len());
    for value in work {
        x.push((-value - max_arg).exp());
    }
    let mean = x.iter().sum::<f64>() / t;
    if mean == 0.0 || !mean.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "EXP uncertainty became non-finite; no finite Boltzmann weight remained".to_string(),
        ));
    }
    let mut var = 0.0;
    for value in &x {
        let diff = value - mean;
        var += diff * diff;
    }
    var /= t;
    let std = var.sqrt();
    let dx = std / (t.sqrt());
    let uncertainty = dx / mean;
    if !uncertainty.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "EXP uncertainty became non-finite".to_string(),
        ));
    }
    Ok(uncertainty)
}

fn exp_delta(w_f: &[f64]) -> Result<f64> {
    let log_sum = logsumexp(&w_f.iter().map(|w| -w).collect::<Vec<_>>());
    let t = w_f.len() as f64;
    let delta = -(log_sum - t.ln());
    if !delta.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "BAR initial estimate became non-finite".to_string(),
        ));
    }
    Ok(delta)
}

fn bar_zero(w_f: &[f64], w_r: &[f64], delta_f: f64) -> f64 {
    let t_f = w_f.len() as f64;
    let t_r = w_r.len() as f64;
    let m = (t_f / t_r).ln();

    let log_f_f: Vec<f64> = w_f.iter().map(|w| neg_log1pexp(m + w - delta_f)).collect();
    let log_numer = logsumexp(&log_f_f);

    let log_f_r: Vec<f64> = w_r
        .iter()
        .map(|w| neg_log1pexp(-(m - w - delta_f)))
        .collect();
    let log_denom = logsumexp(&log_f_r);

    log_numer - log_denom
}

fn neg_log1pexp(value: f64) -> f64 {
    if value == f64::INFINITY {
        f64::NEG_INFINITY
    } else if value == f64::NEG_INFINITY {
        0.0
    } else if value > 0.0 {
        -value - (-value).exp().ln_1p()
    } else {
        -value.exp().ln_1p()
    }
}

fn bar_estimate(
    w_f: &[f64],
    w_r: &[f64],
    method: BarMethod,
    maximum_iterations: usize,
    relative_tolerance: f64,
) -> Result<(f64, f64)> {
    if w_f.is_empty() || w_r.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    if w_f.iter().any(|v| v.is_nan()) || w_r.iter().any(|v| v.is_nan()) {
        return Err(CoreError::NonFiniteValue(
            "BAR work values must not contain NaN".to_string(),
        ));
    }
    if !w_f.iter().any(|v| v.is_finite()) || !w_r.iter().any(|v| v.is_finite()) {
        return Err(CoreError::NonFiniteValue(
            "BAR requires at least one finite work value in each direction".to_string(),
        ));
    }
    let mut delta_f = 0.0;
    let initial_upper = exp_delta(w_f)?;
    let initial_lower = -exp_delta(w_r)?;
    let (mut lower, mut upper, mut f_lower, mut f_upper) =
        bracket_bar_root(w_f, w_r, initial_lower, initial_upper, maximum_iterations)?;

    if !f_upper.is_finite() || !f_lower.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "BAR overlap function became non-finite while bracketing the solution".to_string(),
        ));
    }

    for _ in 0..=maximum_iterations {
        let old = delta_f;
        match method {
            BarMethod::FalsePosition => {
                delta_f = if upper == 0.0 && lower == 0.0 {
                    0.0
                } else {
                    upper - f_upper * (upper - lower) / (f_upper - f_lower)
                };
            }
            BarMethod::Bisection => {
                delta_f = (upper + lower) / 2.0;
            }
            BarMethod::SelfConsistentIteration => {
                delta_f = -bar_zero(w_f, w_r, delta_f) + delta_f;
            }
        }
        let f_new = bar_zero(w_f, w_r, delta_f);
        if delta_f != 0.0 {
            let relative_change = ((delta_f - old) / delta_f).abs();
            if relative_change < relative_tolerance {
                break;
            }
        }

        if matches!(method, BarMethod::FalsePosition | BarMethod::Bisection) {
            if f_upper * f_new < 0.0 {
                lower = delta_f;
                f_lower = f_new;
            } else if f_lower * f_new <= 0.0 {
                upper = delta_f;
                f_upper = f_new;
            } else {
                return Err(CoreError::ConvergenceFailure);
            }
        }
    }

    let d_delta_f = bar_uncertainty(w_f, w_r, delta_f)?;
    Ok((delta_f, d_delta_f))
}

fn bracket_bar_root(
    w_f: &[f64],
    w_r: &[f64],
    initial_lower: f64,
    initial_upper: f64,
    maximum_iterations: usize,
) -> Result<(f64, f64, f64, f64)> {
    let mut lower = initial_lower.min(initial_upper);
    let mut upper = initial_lower.max(initial_upper);
    let mut f_lower = bar_zero(w_f, w_r, lower);
    let mut f_upper = bar_zero(w_f, w_r, upper);

    if !f_lower.is_finite() || !f_upper.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "BAR overlap function became non-finite while bracketing the solution".to_string(),
        ));
    }
    if f_lower == 0.0 || f_upper == 0.0 || f_lower * f_upper < 0.0 {
        return Ok((lower, upper, f_lower, f_upper));
    }

    let mut expansion = (upper - lower).abs().max(0.1);
    for _ in 0..maximum_iterations {
        lower -= expansion;
        upper += expansion;
        f_lower = bar_zero(w_f, w_r, lower);
        f_upper = bar_zero(w_f, w_r, upper);
        if !f_lower.is_finite() || !f_upper.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "BAR overlap function became non-finite while bracketing the solution".to_string(),
            ));
        }
        if f_lower == 0.0 || f_upper == 0.0 || f_lower * f_upper < 0.0 {
            return Ok((lower, upper, f_lower, f_upper));
        }
        expansion *= 2.0;
    }

    Err(CoreError::ConvergenceFailure)
}

fn bar_uncertainty(w_f: &[f64], w_r: &[f64], delta_f: f64) -> Result<f64> {
    if w_f.iter().any(|v| !v.is_finite()) || w_r.iter().any(|v| !v.is_finite()) {
        return Ok(f64::NAN);
    }
    let t_f = w_f.len() as f64;
    let t_r = w_r.len() as f64;
    let m = (t_f / t_r).ln();
    let c = m - delta_f;

    let exp_arg_f: Vec<f64> = w_f.iter().map(|w| w + c).collect();
    let max_arg_f = exp_arg_f.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let log_f_f: Vec<f64> = exp_arg_f
        .iter()
        .map(|arg| -(((-max_arg_f).exp() + (arg - max_arg_f).exp()).ln()))
        .collect();
    let af_f = (logsumexp(&log_f_f) - max_arg_f).exp() / t_f;
    let af_f2 = (logsumexp(&log_f_f.iter().map(|v| 2.0 * v).collect::<Vec<_>>()) - 2.0 * max_arg_f)
        .exp()
        / t_f;

    let exp_arg_r: Vec<f64> = w_r.iter().map(|w| w - c).collect();
    let max_arg_r = exp_arg_r.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let log_f_r: Vec<f64> = exp_arg_r
        .iter()
        .map(|arg| -(((-max_arg_r).exp() + (arg - max_arg_r).exp()).ln()))
        .collect();
    let af_r = (logsumexp(&log_f_r) - max_arg_r).exp() / t_r;
    let af_r2 = (logsumexp(&log_f_r.iter().map(|v| 2.0 * v).collect::<Vec<_>>()) - 2.0 * max_arg_r)
        .exp()
        / t_r;

    let nrat = (t_f + t_r) / (t_f * t_r);
    let variance = (af_f2 / (af_f * af_f)) / t_f + (af_r2 / (af_r * af_r)) / t_r - nrat;
    Ok(variance.max(0.0).sqrt())
}

fn validate_mbar_input(u_kn: &[Vec<f64>]) -> Result<()> {
    if u_kn.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let n_samples = u_kn[0].len();
    for (state_idx, row) in u_kn.iter().enumerate() {
        if row.len() != n_samples {
            return Err(CoreError::InvalidShape {
                expected: n_samples,
                found: row.len(),
            });
        }
        let has_finite = row.iter().any(|value| value.is_finite());
        if !has_finite {
            return Err(CoreError::NonFiniteValue(format!(
                "MBAR state {state_idx} has no finite reduced potentials"
            )));
        }
        for (sample_idx, value) in row.iter().enumerate() {
            if value.is_nan() || *value == f64::NEG_INFINITY {
                return Err(CoreError::NonFiniteValue(format!(
                    "MBAR reduced potential for state {state_idx} sample {sample_idx} must be finite or +inf"
                )));
            }
        }
    }
    for sample_idx in 0..n_samples {
        if !u_kn.iter().any(|row| row[sample_idx].is_finite()) {
            return Err(CoreError::NonFiniteValue(format!(
                "MBAR sample {sample_idx} has no finite reduced potentials"
            )));
        }
    }
    Ok(())
}

fn extract_lambda(state: &StatePoint) -> Result<f64> {
    let lambdas = state.lambdas();
    if lambdas.len() != 1 {
        return Err(CoreError::InvalidState(
            "TI requires exactly one lambda dimension".to_string(),
        ));
    }
    Ok(lambdas[0])
}

fn mean_values(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let sum: f64 = values.iter().sum();
    Ok(sum / (values.len() as f64))
}

fn sem2_values(values: &[f64]) -> Result<f64> {
    if values.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: values.len(),
        });
    }
    let mean = mean_values(values)?;
    let mut sum = 0.0;
    for value in values {
        let diff = value - mean;
        sum += diff * diff;
    }
    let variance = sum / ((values.len() - 1) as f64);
    Ok(variance / (values.len() as f64))
}

fn integrate_trapezoidal(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    let mut total = 0.0;
    for i in 0..(lambdas.len() - 1) {
        let dx = lambdas[i + 1] - lambdas[i];
        total += dx * (values[i] + values[i + 1]) * 0.5;
    }
    Ok(total)
}

fn integrate_simpson(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    if (lambdas.len() & 1) == 0 {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len() + 1,
            found: lambdas.len(),
        });
    }
    let n = lambdas.len();
    let h = (lambdas[n - 1] - lambdas[0]) / ((n - 1) as f64);
    if h == 0.0 {
        return Err(CoreError::InvalidState(
            "lambda spacing must be non-zero".to_string(),
        ));
    }
    let tol = h.abs() * 1e-8;
    for i in 1..n {
        let expected = lambdas[0] + (i as f64) * h;
        if (lambdas[i] - expected).abs() > tol {
            return Err(CoreError::Unsupported(
                "Simpson integration requires uniform lambda spacing".to_string(),
            ));
        }
    }

    let mut total = values[0] + values[n - 1];
    for (i, value) in values.iter().enumerate().take(n - 1).skip(1) {
        if i % 2 == 0 {
            total += 2.0 * value;
        } else {
            total += 4.0 * value;
        }
    }
    Ok(total * h / 3.0)
}

fn trapezoidal_uncertainty(lambdas: &[f64], sem2: &[f64]) -> Result<f64> {
    if lambdas.len() != sem2.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: sem2.len(),
        });
    }
    if lambdas.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: lambdas.len(),
        });
    }
    let mut variance = 0.0;
    for i in 0..lambdas.len() {
        let dl_prev = if i == 0 {
            0.0
        } else {
            lambdas[i] - lambdas[i - 1]
        };
        let dl_next = if i + 1 < lambdas.len() {
            lambdas[i + 1] - lambdas[i]
        } else {
            0.0
        };
        let coeff = dl_prev + dl_next;
        variance += (coeff * coeff) * sem2[i] * 0.25;
    }
    Ok(variance.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::StatePoint;

    fn make_two_state_windows() -> Vec<UNkMatrix> {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let states = vec![s0.clone(), s1.clone()];

        let data0 = vec![0.0, 0.2, 0.1, 0.3, 0.2, 0.4];
        let data1 = vec![0.1, 0.0, 0.2, 0.1, 0.3, 0.2];
        let time = vec![0.0, 1.0, 2.0];

        vec![
            UNkMatrix::new(3, 2, data0, time.clone(), Some(s0), states.clone()).unwrap(),
            UNkMatrix::new(3, 2, data1, time, Some(s1), states).unwrap(),
        ]
    }

    fn make_window(
        sampled_lambda: f64,
        sampled_temperature: f64,
        evaluated_lambdas: [f64; 2],
        evaluated_temperature: f64,
    ) -> UNkMatrix {
        let sampled_state = StatePoint::new(vec![sampled_lambda], sampled_temperature).unwrap();
        let evaluated_states = evaluated_lambdas
            .into_iter()
            .map(|lambda| StatePoint::new(vec![lambda], evaluated_temperature).unwrap())
            .collect::<Vec<_>>();
        let data = vec![0.0, 0.2, 0.1, 0.3, 0.2, 0.4];
        let time = vec![0.0, 1.0, 2.0];
        UNkMatrix::new(3, 2, data, time, Some(sampled_state), evaluated_states).unwrap()
    }

    fn make_two_state_windows_with_positive_infinity() -> Vec<UNkMatrix> {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let states = vec![s0.clone(), s1.clone()];

        let data0 = vec![0.0, 1.0, 0.0, 2.0, 0.0, f64::INFINITY];
        let data1 = vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0];
        let time = vec![0.0, 1.0, 2.0];

        vec![
            UNkMatrix::new(3, 2, data0, time.clone(), Some(s0), states.clone()).unwrap(),
            UNkMatrix::new(3, 2, data1, time, Some(s1), states).unwrap(),
        ]
    }

    fn assert_vec_eq_with_nan(left: Option<&[f64]>, right: Option<&[f64]>) {
        match (left, right) {
            (None, None) => {}
            (Some(a), Some(b)) => {
                assert_eq!(a.len(), b.len());
                for (idx, (x, y)) in a.iter().zip(b.iter()).enumerate() {
                    if x.is_nan() && y.is_nan() {
                        continue;
                    }
                    assert!((*x - *y).abs() < 1e-12, "mismatch at {idx}: {x} vs {y}");
                }
            }
            _ => panic!("option mismatch"),
        }
    }

    #[test]
    fn ti_fit_requires_two_windows() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let series = DhdlSeries::new(state, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let estimator = TiEstimator::default();
        let err = estimator.fit(&[series]).unwrap_err();
        assert!(matches!(err, CoreError::InvalidShape { .. }));
    }

    #[test]
    fn ti_trapezoidal_two_points() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap();
        let estimator = TiEstimator::default();
        let result = estimator.fit(&[d0, d1]).unwrap();
        assert!((result.delta_f() - 1.0).abs() < 1e-12);
        assert_eq!(result.uncertainty(), Some(0.0));
    }

    #[test]
    fn ti_simpson_three_points() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d2 = DhdlSeries::new(s2, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Simpson,
            parallel: false,
        });
        let result = estimator.fit(&[d0, d1, d2]).unwrap();
        assert!((result.delta_f() - 1.0).abs() < 1e-12);
        assert_eq!(result.uncertainty(), None);
    }

    #[test]
    fn ti_simpson_rejects_nonuniform_spacing() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.3], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d2 = DhdlSeries::new(s2, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Simpson,
            parallel: false,
        });
        let err = estimator.fit(&[d0, d1, d2]).unwrap_err();
        assert!(matches!(err, CoreError::Unsupported(_)));
    }

    #[test]
    fn mbar_requires_windows() {
        let estimator = MbarEstimator::default();
        let err = estimator.fit(&[]).unwrap_err();
        assert!(matches!(err, CoreError::InvalidShape { .. }));
    }

    #[test]
    fn ti_parallel_matches_serial() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0, 2.0], vec![0.0, 0.1, 0.2]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0, 2.0], vec![1.0, 1.1, 1.2]).unwrap();

        let serial = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Trapezoidal,
            parallel: false,
        })
        .fit(&[d0.clone(), d1.clone()])
        .unwrap();
        let parallel = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Trapezoidal,
            parallel: true,
        })
        .fit(&[d0, d1])
        .unwrap();

        assert!((serial.delta_f() - parallel.delta_f()).abs() < 1e-12);
        assert_eq!(serial.uncertainty(), parallel.uncertainty());
    }

    #[test]
    fn bar_parallel_matches_serial() {
        let windows = make_two_state_windows();
        let serial = BarEstimator::new(BarOptions {
            parallel: false,
            ..BarOptions::default()
        })
        .fit(&windows)
        .unwrap();
        let parallel = BarEstimator::new(BarOptions {
            parallel: true,
            ..BarOptions::default()
        })
        .fit(&windows)
        .unwrap();

        assert_eq!(serial.values(), parallel.values());
        assert_vec_eq_with_nan(serial.uncertainties(), parallel.uncertainties());
    }

    #[test]
    fn exp_parallel_matches_serial() {
        let windows = make_two_state_windows();
        let serial = ExpEstimator::new(ExpOptions {
            parallel: false,
            ..ExpOptions::default()
        })
        .fit(&windows)
        .unwrap();
        let parallel = ExpEstimator::new(ExpOptions {
            parallel: true,
            ..ExpOptions::default()
        })
        .fit(&windows)
        .unwrap();

        assert_eq!(serial.values(), parallel.values());
        assert_vec_eq_with_nan(serial.uncertainties(), parallel.uncertainties());
    }

    #[test]
    fn mbar_parallel_matches_serial() {
        let windows = make_two_state_windows();
        let serial = MbarEstimator::new(MbarOptions {
            parallel: false,
            ..MbarOptions::default()
        })
        .fit(&windows)
        .unwrap();
        let parallel = MbarEstimator::new(MbarOptions {
            parallel: true,
            ..MbarOptions::default()
        })
        .fit(&windows)
        .unwrap();

        assert_eq!(serial.values(), parallel.values());
        assert_vec_eq_with_nan(serial.uncertainties(), parallel.uncertainties());
    }

    #[test]
    fn bar_rejects_mismatched_evaluated_state_grid() {
        let windows = vec![
            make_window(0.0, 300.0, [0.0, 1.0], 300.0),
            make_window(1.0, 300.0, [0.0, 2.0], 300.0),
        ];

        let err = BarEstimator::default().fit(&windows).unwrap_err();
        assert!(
            matches!(err, CoreError::InvalidState(message) if message == "evaluated_states differ between windows")
        );
    }

    #[test]
    fn mbar_rejects_sampled_state_temperature_mismatch() {
        let windows = vec![
            make_window(0.0, 300.0, [0.0, 1.0], 300.0),
            make_window(1.0, 310.0, [0.0, 1.0], 300.0),
        ];

        let err = MbarEstimator::default().fit(&windows).unwrap_err();
        assert!(
            matches!(err, CoreError::InvalidState(message) if message == "sampled_state temperature differs from evaluated_states")
        );
    }

    #[test]
    fn exp_rejects_sampled_state_missing_from_grid() {
        let windows = vec![
            make_window(0.0, 300.0, [0.0, 1.0], 300.0),
            make_window(0.5, 300.0, [0.0, 1.0], 300.0),
        ];

        let err = ExpEstimator::default().fit(&windows).unwrap_err();
        assert!(
            matches!(err, CoreError::InvalidState(message) if message == "sampled_state not found in evaluated_states")
        );
    }

    #[test]
    fn bar_rejects_duplicate_windows_for_same_sampled_state() {
        let windows = vec![
            make_window(0.0, 300.0, [0.0, 1.0], 300.0),
            make_window(0.0, 300.0, [0.0, 1.0], 300.0),
        ];

        let err = BarEstimator::default().fit(&windows).unwrap_err();
        assert!(
            matches!(err, CoreError::InvalidState(message) if message == "multiple windows for sampled_state lambda 0")
        );
    }

    #[test]
    fn exp_supports_positive_infinity_work_values() {
        let windows = make_two_state_windows_with_positive_infinity();
        let result = ExpEstimator::default().fit(&windows).unwrap();
        assert!(result.values().iter().all(|value| value.is_finite()));
    }

    #[test]
    fn bar_supports_positive_infinity_work_values() {
        let windows = make_two_state_windows_with_positive_infinity();
        let result = BarEstimator::default().fit(&windows).unwrap();
        assert!(result.values().iter().all(|value| value.is_finite()));
        assert!(result
            .uncertainties()
            .expect("uncertainties")
            .iter()
            .any(|value| value.is_nan()));
    }

    #[test]
    fn bar_bracketing_respects_iteration_limit() {
        let err = bar_estimate(
            &[5.0, 5.0],
            &[-10.0, 5.0],
            BarMethod::FalsePosition,
            0,
            1.0e-7,
        )
        .unwrap_err();

        assert!(matches!(err, CoreError::ConvergenceFailure));
    }

    #[test]
    fn mbar_supports_positive_infinity_reduced_potentials() {
        let windows = make_two_state_windows_with_positive_infinity();
        let result = MbarEstimator::new(MbarOptions {
            compute_uncertainty: false,
            ..MbarOptions::default()
        })
        .fit(&windows)
        .unwrap();
        assert!(result.values().iter().all(|value| value.is_finite()));
    }
}
