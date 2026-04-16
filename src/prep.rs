use crate::data::{find_scalar_lambda_state_index_exact, DhdlSeries, UNkMatrix};
use crate::error::{CoreError, Result};
use rustfft::{num_complex::Complex, FftPlanner};
use std::borrow::Cow;
use std::cell::RefCell;

type SeriesCowPair<'a> = (Cow<'a, [f64]>, Cow<'a, [f64]>);
const DIRECT_AUTOCORRELATION_TRIAL_STEPS: usize = 64;

thread_local! {
    static FFT_PLANNER: RefCell<FftPlanner<f64>> = RefCell::new(FftPlanner::<f64>::new());
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UNkSeriesMethod {
    All,
    DE,
}

#[derive(Debug, Clone)]
pub struct DecorrelationOptions {
    pub drop_duplicates: bool,
    pub sort: bool,
    pub conservative: bool,
    pub remove_burnin: bool,
    pub fast: bool,
    pub nskip: usize,
    pub lower: Option<f64>,
    pub upper: Option<f64>,
    pub step: Option<usize>,
}

impl Default for DecorrelationOptions {
    fn default() -> Self {
        Self {
            drop_duplicates: true,
            sort: true,
            conservative: true,
            remove_burnin: false,
            fast: false,
            nskip: 1,
            lower: None,
            upper: None,
            step: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibrationResult {
    pub t0: usize,
    pub g: f64,
    pub neff_max: f64,
}

pub fn decorrelate_dhdl(series: &DhdlSeries, options: &DecorrelationOptions) -> Result<DhdlSeries> {
    let (time, values) = prepare_series_view(
        series.time_ps(),
        series.values(),
        options.drop_duplicates,
        options.sort,
    )?;
    let (time, values) = apply_time_slice_view(time.as_ref(), values.as_ref(), options)?;
    let indices = subsample_indices(values.as_ref(), options)?;
    let (time, values) = apply_indices(time.as_ref(), values.as_ref(), &indices);
    DhdlSeries::new(series.state().clone(), time, values)
}

pub fn detect_equilibration_dhdl(
    series: &DhdlSeries,
    options: &DecorrelationOptions,
) -> Result<EquilibrationResult> {
    let (time, values) = prepare_series_view(
        series.time_ps(),
        series.values(),
        options.drop_duplicates,
        options.sort,
    )?;
    let (_time, values) = apply_time_slice_view(time.as_ref(), values.as_ref(), options)?;
    detect_equilibration(values.as_ref(), options.fast, options.nskip)
}

pub fn decorrelate_u_nk(
    u_nk: &UNkMatrix,
    method: UNkSeriesMethod,
    options: &DecorrelationOptions,
) -> Result<UNkMatrix> {
    let (time, data) = prepare_u_nk(
        u_nk.time_ps(),
        u_nk.data(),
        u_nk.n_states(),
        options.drop_duplicates,
        options.sort,
    )?;

    let series = u_nk_series(u_nk, &data, method)?;
    let (time, data, series) =
        apply_time_slice_u_nk(&time, &data, u_nk.n_states(), &series, options)?;
    let indices = subsample_indices(&series, options)?;

    let (time, data) = apply_indices_u_nk(&time, &data, u_nk.n_states(), &indices);

    UNkMatrix::new_with_labels(
        indices.len(),
        u_nk.n_states(),
        data,
        time,
        u_nk.sampled_state().cloned(),
        u_nk.evaluated_states().to_vec(),
        u_nk.lambda_labels().map(|labels| labels.to_vec()),
    )
}

pub fn decorrelate_u_nk_with_observable(
    u_nk: &UNkMatrix,
    observable: &[f64],
    options: &DecorrelationOptions,
) -> Result<UNkMatrix> {
    let (time, data, observable) = prepare_u_nk_with_observable(
        u_nk.time_ps(),
        u_nk.data(),
        u_nk.n_states(),
        observable,
        options.drop_duplicates,
        options.sort,
    )?;
    let (time, data, observable) =
        apply_time_slice_u_nk(&time, &data, u_nk.n_states(), &observable, options)?;
    let indices = subsample_indices(&observable, options)?;
    let (time, data) = apply_indices_u_nk(&time, &data, u_nk.n_states(), &indices);

    UNkMatrix::new_with_labels(
        indices.len(),
        u_nk.n_states(),
        data,
        time,
        u_nk.sampled_state().cloned(),
        u_nk.evaluated_states().to_vec(),
        u_nk.lambda_labels().map(|labels| labels.to_vec()),
    )
}

pub fn detect_equilibration_u_nk(
    u_nk: &UNkMatrix,
    method: UNkSeriesMethod,
    options: &DecorrelationOptions,
) -> Result<EquilibrationResult> {
    let (time, data) = prepare_u_nk(
        u_nk.time_ps(),
        u_nk.data(),
        u_nk.n_states(),
        options.drop_duplicates,
        options.sort,
    )?;
    let series = u_nk_series(u_nk, &data, method)?;
    let (time, _data, series) =
        apply_time_slice_u_nk(&time, &data, u_nk.n_states(), &series, options)?;
    let _ = time;
    detect_equilibration(&series, options.fast, options.nskip)
}

pub fn detect_equilibration_observable(
    observable: &[f64],
    fast: bool,
    nskip: usize,
) -> Result<EquilibrationResult> {
    if observable.iter().any(|value| !value.is_finite()) {
        return Err(CoreError::NonFiniteValue(
            "equilibration observable must be finite".to_string(),
        ));
    }
    detect_equilibration(observable, fast, nskip)
}

fn prepare_series_view<'a>(
    time: &'a [f64],
    values: &'a [f64],
    drop_duplicates: bool,
    sort: bool,
) -> Result<SeriesCowPair<'a>> {
    if time.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: time.len(),
            found: values.len(),
        });
    }

    let has_duplicates = has_duplicates(time);
    if has_duplicates && !drop_duplicates {
        return Err(CoreError::InvalidState(
            "duplicate time values found".to_string(),
        ));
    }

    let sorted = is_sorted(time);
    if !sorted && !sort {
        return Err(CoreError::InvalidTimeOrder(
            "time values are not sorted".to_string(),
        ));
    }

    if !has_duplicates && sorted {
        return Ok((Cow::Borrowed(time), Cow::Borrowed(values)));
    }

    let mut time = time.to_vec();
    let mut values = values.to_vec();
    if has_duplicates {
        (time, values) = drop_duplicate_times(&time, &values);
    }
    if !sorted {
        (time, values) = sort_by_time(&time, &values);
    }
    Ok((Cow::Owned(time), Cow::Owned(values)))
}

fn prepare_u_nk(
    time: &[f64],
    data: &[f64],
    n_states: usize,
    drop_duplicates: bool,
    sort: bool,
) -> Result<(Vec<f64>, Vec<f64>)> {
    if time.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let mut time = time.to_vec();
    let mut data = data.to_vec();
    if has_duplicates(&time) {
        if drop_duplicates {
            let (new_time, indices) = unique_time_indices(&time);
            time = new_time;
            data = select_u_nk_rows(&data, n_states, &indices);
        } else {
            return Err(CoreError::InvalidState(
                "duplicate time values found".to_string(),
            ));
        }
    }
    if !is_sorted(&time) {
        if sort {
            let (sorted_time, indices) = sorted_time_indices(&time);
            time = sorted_time;
            data = select_u_nk_rows(&data, n_states, &indices);
        } else {
            return Err(CoreError::InvalidTimeOrder(
                "time values are not sorted".to_string(),
            ));
        }
    }
    Ok((time, data))
}

fn prepare_u_nk_with_observable(
    time: &[f64],
    data: &[f64],
    n_states: usize,
    observable: &[f64],
    drop_duplicates: bool,
    sort: bool,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    if time.len() != observable.len() {
        return Err(CoreError::InvalidShape {
            expected: time.len(),
            found: observable.len(),
        });
    }
    if observable.iter().any(|value| !value.is_finite()) {
        return Err(CoreError::NonFiniteValue(
            "decorrelation observable must be finite".to_string(),
        ));
    }
    let mut time = time.to_vec();
    let mut data = data.to_vec();
    let mut observable = observable.to_vec();
    if has_duplicates(&time) {
        if drop_duplicates {
            let (new_time, indices) = unique_time_indices(&time);
            time = new_time;
            data = select_u_nk_rows(&data, n_states, &indices);
            observable = select_values(&observable, &indices);
        } else {
            return Err(CoreError::InvalidState(
                "duplicate time values found".to_string(),
            ));
        }
    }
    if !is_sorted(&time) {
        if sort {
            let (sorted_time, indices) = sorted_time_indices(&time);
            time = sorted_time;
            data = select_u_nk_rows(&data, n_states, &indices);
            observable = select_values(&observable, &indices);
        } else {
            return Err(CoreError::InvalidTimeOrder(
                "time values are not sorted".to_string(),
            ));
        }
    }
    Ok((time, data, observable))
}

fn apply_time_slice_view<'a>(
    time: &'a [f64],
    values: &'a [f64],
    options: &DecorrelationOptions,
) -> Result<SeriesCowPair<'a>> {
    if options.lower.is_none() && options.upper.is_none() && options.step.is_none() {
        if time.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        return Ok((Cow::Borrowed(time), Cow::Borrowed(values)));
    }

    let indices = slice_time_indices(time, options)?;
    let (time, values) = apply_indices(time, values, &indices);
    Ok((Cow::Owned(time), Cow::Owned(values)))
}

fn apply_time_slice_u_nk(
    time: &[f64],
    data: &[f64],
    n_states: usize,
    series: &[f64],
    options: &DecorrelationOptions,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let indices = slice_time_indices(time, options)?;
    let (time, series) = apply_indices(time, series, &indices);
    let data = select_u_nk_rows(data, n_states, &indices);
    Ok((time, data, series))
}

fn subsample_indices(values: &[f64], options: &DecorrelationOptions) -> Result<Vec<usize>> {
    if options.remove_burnin {
        let result = detect_equilibration(values, options.fast, options.nskip)?;
        let t = result.t0;
        let g = result.g;
        let values = &values[t..];
        let indices =
            subsample_correlated_data(values, Some(g), options.fast, options.conservative)?;
        Ok(indices.into_iter().map(|idx| idx + t).collect())
    } else {
        subsample_correlated_data(values, None, options.fast, options.conservative)
    }
}

fn apply_indices(time: &[f64], values: &[f64], indices: &[usize]) -> (Vec<f64>, Vec<f64>) {
    let mut out_time = Vec::with_capacity(indices.len());
    let mut out_values = Vec::with_capacity(indices.len());
    for &idx in indices {
        out_time.push(time[idx]);
        out_values.push(values[idx]);
    }
    (out_time, out_values)
}

fn apply_indices_u_nk(
    time: &[f64],
    data: &[f64],
    n_states: usize,
    indices: &[usize],
) -> (Vec<f64>, Vec<f64>) {
    let mut out_time = Vec::with_capacity(indices.len());
    let mut out_data = Vec::with_capacity(indices.len() * n_states);
    for &idx in indices {
        out_time.push(time[idx]);
        let offset = idx * n_states;
        out_data.extend_from_slice(&data[offset..offset + n_states]);
    }
    (out_time, out_data)
}

fn slice_time_indices(time: &[f64], options: &DecorrelationOptions) -> Result<Vec<usize>> {
    let mut indices = Vec::new();
    let lower = options.lower.unwrap_or(f64::NEG_INFINITY);
    let upper = options.upper.unwrap_or(f64::INFINITY);
    if lower > upper {
        return Err(CoreError::InvalidState(
            "lower bound is greater than upper bound".to_string(),
        ));
    }
    for (idx, &t) in time.iter().enumerate() {
        if t < lower || t > upper {
            continue;
        }
        indices.push(idx);
    }
    if let Some(step) = options.step {
        if step == 0 {
            return Err(CoreError::InvalidState(
                "step must be greater than zero".to_string(),
            ));
        }
        indices = indices.into_iter().step_by(step).collect();
    }
    if indices.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    Ok(indices)
}

fn u_nk_series(u_nk: &UNkMatrix, data: &[f64], method: UNkSeriesMethod) -> Result<Vec<f64>> {
    let n_states = u_nk.n_states();
    let n_samples = data.len() / n_states;
    let mut series = Vec::with_capacity(n_samples);
    match method {
        UNkSeriesMethod::All => {
            for sample_idx in 0..n_samples {
                let offset = sample_idx * n_states;
                let sum: f64 = data[offset..offset + n_states].iter().sum();
                series.push(sum);
            }
        }
        UNkSeriesMethod::DE => {
            let sampled = u_nk.sampled_state().ok_or_else(|| {
                CoreError::InvalidState("sampled_state must be set for dE".to_string())
            })?;
            if sampled.lambdas().len() != 1
                || u_nk
                    .evaluated_states()
                    .iter()
                    .any(|state| state.lambdas().len() != 1)
            {
                return Err(CoreError::RequiresOneDimensionalLambda {
                    operation: "DE u_nk preprocessing",
                });
            }
            if n_states < 2 {
                return Err(CoreError::InvalidShape {
                    expected: 2,
                    found: n_states,
                });
            }
            let sampled_lambda = sampled.lambdas()[0];
            let index =
                find_scalar_lambda_state_index_exact(u_nk.evaluated_states(), sampled_lambda)?;
            let other_index = if index + 1 < n_states {
                index + 1
            } else {
                index - 1
            };
            for sample_idx in 0..n_samples {
                let offset = sample_idx * n_states;
                let current = data[offset + index];
                let other = data[offset + other_index];
                series.push(other - current);
            }
        }
    }
    if let Some((idx, _)) = series
        .iter()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(CoreError::NonFiniteValue(format!(
            "u_nk decorrelation series contains a non-finite value at sample {idx}; decorrelation is unsupported for +inf reduced energies"
        )));
    }
    Ok(series)
}

fn has_duplicates(time: &[f64]) -> bool {
    let mut seen = std::collections::HashSet::new();
    for &value in time {
        if !seen.insert(value.to_bits()) {
            return true;
        }
    }
    false
}

fn drop_duplicate_times(time: &[f64], values: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut seen = std::collections::HashSet::new();
    let mut out_time = Vec::new();
    let mut out_values = Vec::new();
    for (&t, &v) in time.iter().zip(values.iter()) {
        if seen.insert(t.to_bits()) {
            out_time.push(t);
            out_values.push(v);
        }
    }
    (out_time, out_values)
}

fn unique_time_indices(time: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut seen = std::collections::HashSet::new();
    let mut out_time = Vec::new();
    let mut indices = Vec::new();
    for (idx, &t) in time.iter().enumerate() {
        if seen.insert(t.to_bits()) {
            out_time.push(t);
            indices.push(idx);
        }
    }
    (out_time, indices)
}

fn sorted_time_indices(time: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut pairs: Vec<(usize, f64)> = time.iter().cloned().enumerate().collect();
    pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let sorted_time: Vec<f64> = pairs.iter().map(|(_, t)| *t).collect();
    let indices: Vec<usize> = pairs.iter().map(|(idx, _)| *idx).collect();
    (sorted_time, indices)
}

fn sort_by_time(time: &[f64], values: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let (sorted_time, indices) = sorted_time_indices(time);
    let mut sorted_values = Vec::with_capacity(values.len());
    for idx in indices {
        sorted_values.push(values[idx]);
    }
    (sorted_time, sorted_values)
}

fn is_sorted(time: &[f64]) -> bool {
    time.windows(2).all(|pair| pair[1] >= pair[0])
}

fn select_u_nk_rows(data: &[f64], n_states: usize, indices: &[usize]) -> Vec<f64> {
    let mut out = Vec::with_capacity(indices.len() * n_states);
    for &idx in indices {
        let offset = idx * n_states;
        out.extend_from_slice(&data[offset..offset + n_states]);
    }
    out
}

fn select_values(values: &[f64], indices: &[usize]) -> Vec<f64> {
    let mut out = Vec::with_capacity(indices.len());
    for &idx in indices {
        out.push(values[idx]);
    }
    out
}

fn statistical_inefficiency(values: &[f64], fast: bool) -> Result<f64> {
    if values.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let mut d_values = Vec::with_capacity(values.len());
    for &value in values {
        d_values.push(value - mean);
    }
    let sigma2 = d_values.iter().map(|v| v * v).sum::<f64>() / values.len() as f64;
    if sigma2 == 0.0 {
        return Err(CoreError::InvalidState(
            "sample covariance sigma^2 = 0".to_string(),
        ));
    }

    if let Some(g) = statistical_inefficiency_direct_limited(
        &d_values,
        sigma2,
        fast,
        DIRECT_AUTOCORRELATION_TRIAL_STEPS,
    )? {
        return Ok(g);
    }
    statistical_inefficiency_fft(&d_values, sigma2, fast)
}

#[cfg(test)]
fn statistical_inefficiency_direct(d_values: &[f64], sigma2: f64, fast: bool) -> Result<f64> {
    statistical_inefficiency_direct_inner(d_values, sigma2, fast, None).map(|result| {
        result.expect("unbounded direct statistical inefficiency should always complete")
    })
}

fn statistical_inefficiency_direct_limited(
    d_values: &[f64],
    sigma2: f64,
    fast: bool,
    max_steps: usize,
) -> Result<Option<f64>> {
    statistical_inefficiency_direct_inner(d_values, sigma2, fast, Some(max_steps))
}

fn statistical_inefficiency_direct_inner(
    d_values: &[f64],
    sigma2: f64,
    fast: bool,
    max_steps: Option<usize>,
) -> Result<Option<f64>> {
    let n = d_values.len();
    let mut g = 1.0;
    let mut t = 1usize;
    let mut increment = 1usize;
    let mintime = 3usize;
    let mut steps = 0usize;
    while t < n - 1 {
        if max_steps.is_some_and(|limit| steps >= limit) {
            return Ok(None);
        }
        let mut sum = 0.0;
        let denom = (n - t) as f64;
        for i in 0..(n - t) {
            sum += d_values[i] * d_values[i + t] + d_values[i + t] * d_values[i];
        }
        let c = sum / (2.0 * denom * sigma2);
        if c <= 0.0 && t > mintime {
            break;
        }
        g += 2.0 * c * (1.0 - (t as f64) / (n as f64)) * (increment as f64);
        t += increment;
        if fast {
            increment += 1;
        }
        steps += 1;
    }
    if g < 1.0 {
        g = 1.0;
    }
    Ok(Some(g))
}

fn statistical_inefficiency_fft(d_values: &[f64], sigma2: f64, fast: bool) -> Result<f64> {
    let n = d_values.len();
    let autocov = autocovariance_fft(d_values);
    let mut g = 1.0;
    let mut t = 1usize;
    let mut increment = 1usize;
    let mintime = 3usize;
    while t < n - 1 {
        let c = autocov[t] / ((n - t) as f64 * sigma2);
        if c <= 0.0 && t > mintime {
            break;
        }
        g += 2.0 * c * (1.0 - (t as f64) / (n as f64)) * (increment as f64);
        t += increment;
        if fast {
            increment += 1;
        }
    }
    if g < 1.0 {
        g = 1.0;
    }
    Ok(g)
}

fn autocovariance_fft(d_values: &[f64]) -> Vec<f64> {
    let n = d_values.len();
    let fft_len = (n * 2).next_power_of_two();
    let (forward, inverse) = FFT_PLANNER.with(|planner| {
        let mut planner = planner.borrow_mut();
        (
            planner.plan_fft_forward(fft_len),
            planner.plan_fft_inverse(fft_len),
        )
    });

    let mut buffer = vec![Complex::new(0.0, 0.0); fft_len];
    for (slot, &value) in buffer.iter_mut().zip(d_values.iter()) {
        slot.re = value;
    }

    forward.process(&mut buffer);
    for value in &mut buffer {
        *value *= value.conj();
    }
    inverse.process(&mut buffer);

    let scale = 1.0 / fft_len as f64;
    buffer
        .into_iter()
        .take(n)
        .map(|value| value.re * scale)
        .collect()
}

fn subsample_correlated_data(
    values: &[f64],
    g: Option<f64>,
    fast: bool,
    conservative: bool,
) -> Result<Vec<usize>> {
    let g = match g {
        Some(value) => value,
        None => statistical_inefficiency(values, fast)?,
    };
    if conservative {
        let stride = g.ceil() as usize;
        Ok((0..values.len()).step_by(stride.max(1)).collect())
    } else {
        let mut indices = Vec::new();
        let mut n = 0usize;
        loop {
            let t = (n as f64 * g).round() as usize;
            if t >= values.len() {
                break;
            }
            if indices.last().copied() != Some(t) {
                indices.push(t);
            }
            n += 1;
        }
        Ok(indices)
    }
}

fn detect_equilibration(values: &[f64], fast: bool, nskip: usize) -> Result<EquilibrationResult> {
    if values.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|v| {
            let d = v - mean;
            d * d
        })
        .sum::<f64>()
        / values.len() as f64;
    if variance == 0.0 {
        return Ok(EquilibrationResult {
            t0: 0,
            g: 1.0,
            neff_max: 1.0,
        });
    }

    let step = nskip.max(1);
    let mut best_result = None;
    for t in (0..values.len().saturating_sub(1)).step_by(step) {
        let slice = &values[t..];
        let g = statistical_inefficiency(slice, fast).unwrap_or(values.len() as f64);
        // Use the actual retained sample count for the suffix `values[t..]`.
        // This intentionally differs from pymbar/alchemlyb, which use `T - t + 1`;
        // on the bundled equilibration fixture that changes `neff_max` from 21 to 20.
        let neff = (values.len() - t) as f64 / g;
        let candidate = EquilibrationResult {
            t0: t,
            g,
            neff_max: neff,
        };
        // When `nskip > 1`, only the sampled time origins are valid candidates.
        // This intentionally deviates from pymbar's current implementation, which
        // pre-fills arrays for all indices and then takes `argmax()` over those
        // placeholders as well. The documented meaning of `nskip` is to try every
        // `nskip`-th sample as a potential time origin, so skipped indices are
        // treated as ineligible here.
        best_result = Some(best_equilibration_result(best_result, candidate));
    }
    best_result.ok_or(CoreError::InvalidShape {
        expected: 2,
        found: values.len(),
    })
}

fn best_equilibration_result(
    best: Option<EquilibrationResult>,
    candidate: EquilibrationResult,
) -> EquilibrationResult {
    match best {
        Some(current) if current.neff_max >= candidate.neff_max => current,
        _ => candidate,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::StatePoint;
    use std::time::Instant;

    #[test]
    fn decorrelate_dhdl_drops_duplicates() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let series = DhdlSeries::new(state, vec![0.0, 0.0, 1.0], vec![1.0, 2.0, 3.0]).unwrap();
        let result = decorrelate_dhdl(&series, &DecorrelationOptions::default()).unwrap();
        assert_eq!(result.time_ps().len(), 2);
    }

    #[test]
    fn decorrelate_u_nk_all_keeps_shape() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let u_nk = UNkMatrix::new(
            3,
            2,
            vec![1.0, 2.0, 1.1, 2.1, 1.2, 2.2],
            vec![0.0, 1.0, 2.0],
            Some(state.clone()),
            vec![state.clone(), StatePoint::new(vec![1.0], 300.0).unwrap()],
        )
        .unwrap();
        let result = decorrelate_u_nk(
            &u_nk,
            UNkSeriesMethod::All,
            &DecorrelationOptions::default(),
        )
        .unwrap();
        assert!(result.n_samples() > 0);
        assert_eq!(result.n_states(), 2);
    }

    #[test]
    fn detect_equilibration_constant_series() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let series = DhdlSeries::new(state, vec![0.0, 1.0, 2.0], vec![1.0, 1.0, 1.0]).unwrap();
        let result = detect_equilibration_dhdl(&series, &DecorrelationOptions::default()).unwrap();
        assert_eq!(result.t0, 0);
        assert_eq!(result.g, 1.0);
        assert_eq!(result.neff_max, 1.0);
    }

    #[test]
    fn detect_equilibration_observable_rejects_nonfinite_values() {
        let err = detect_equilibration_observable(&[1.0, f64::INFINITY], false, 1).unwrap_err();
        assert!(matches!(err, CoreError::NonFiniteValue(_)));
    }

    #[test]
    fn best_equilibration_result_uses_only_explicit_candidates() {
        let best = best_equilibration_result(
            Some(EquilibrationResult {
                t0: 0,
                g: 3.0,
                neff_max: 2.0,
            }),
            EquilibrationResult {
                t0: 4,
                g: 2.0,
                neff_max: 3.0,
            },
        );

        assert_eq!(best.t0, 4);
        assert_eq!(best.g, 2.0);
        assert_eq!(best.neff_max, 3.0);
    }

    #[test]
    fn decorrelate_u_nk_de_rejects_multidimensional_states() {
        let sampled = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
        let neighbor = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
        let u_nk = UNkMatrix::new(
            2,
            2,
            vec![0.0, 1.0, 0.0, 2.0],
            vec![0.0, 1.0],
            Some(sampled.clone()),
            vec![sampled, neighbor],
        )
        .unwrap();

        let err = decorrelate_u_nk(&u_nk, UNkSeriesMethod::DE, &DecorrelationOptions::default())
            .unwrap_err();
        assert!(matches!(
            err,
            CoreError::RequiresOneDimensionalLambda {
                operation: "DE u_nk preprocessing"
            }
        ));
    }

    #[test]
    fn decorrelate_u_nk_de_requires_two_states() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let u_nk = UNkMatrix::new(
            2,
            1,
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            Some(state.clone()),
            vec![state],
        )
        .unwrap();

        let err = decorrelate_u_nk(&u_nk, UNkSeriesMethod::DE, &DecorrelationOptions::default())
            .unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 2,
                found: 1
            }
        ));
    }

    #[test]
    fn detect_equilibration_u_nk_de_requires_two_states() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let u_nk = UNkMatrix::new(
            2,
            1,
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            Some(state.clone()),
            vec![state],
        )
        .unwrap();

        let err =
            detect_equilibration_u_nk(&u_nk, UNkSeriesMethod::DE, &DecorrelationOptions::default())
                .unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 2,
                found: 1
            }
        ));
    }

    #[test]
    fn decorrelate_u_nk_rejects_positive_infinity_series() {
        let state0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let u_nk = UNkMatrix::new(
            2,
            2,
            vec![0.0, f64::INFINITY, 0.0, 1.0],
            vec![0.0, 1.0],
            Some(state0.clone()),
            vec![state0, state1],
        )
        .unwrap();

        let err = decorrelate_u_nk(&u_nk, UNkSeriesMethod::DE, &DecorrelationOptions::default())
            .unwrap_err();
        assert!(matches!(err, CoreError::NonFiniteValue(_)));
    }

    #[test]
    fn decorrelate_u_nk_with_observable_accepts_positive_infinity_matrix() {
        let state0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let u_nk = UNkMatrix::new(
            4,
            2,
            vec![0.0, 1.0, 0.0, f64::INFINITY, 0.0, 2.0, 0.0, 3.0],
            vec![0.0, 1.0, 2.0, 3.0],
            Some(state0.clone()),
            vec![state0, state1],
        )
        .unwrap();

        let result = decorrelate_u_nk_with_observable(
            &u_nk,
            &[10.0, 11.0, 12.0, 13.0],
            &DecorrelationOptions::default(),
        )
        .unwrap();

        assert_eq!(result.n_samples(), 4);
        assert_eq!(result.data()[3], f64::INFINITY);
    }

    #[test]
    fn statistical_inefficiency_fft_matches_direct_for_long_series() {
        let values: Vec<f64> = (0..4096)
            .map(|idx| {
                let x = idx as f64;
                (x / 17.0).sin() + 0.3 * (x / 31.0).cos() + 0.05 * (x / 7.0).sin()
            })
            .collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let centered: Vec<f64> = values.iter().map(|value| value - mean).collect();
        let sigma2 =
            centered.iter().map(|value| value * value).sum::<f64>() / centered.len() as f64;

        let direct = statistical_inefficiency_direct(&centered, sigma2, false).unwrap();
        let fft = statistical_inefficiency_fft(&centered, sigma2, false).unwrap();

        assert!((fft - direct).abs() < 1.0e-9, "{fft} vs {direct}");
    }

    #[test]
    #[ignore = "benchmark"]
    fn benchmark_statistical_inefficiency_backends() {
        let values: Vec<f64> = (0..262_144)
            .map(|idx| {
                let x = idx as f64;
                (x / 17.0).sin() + 0.3 * (x / 31.0).cos() + 0.05 * (x / 7.0).sin()
            })
            .collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let centered: Vec<f64> = values.iter().map(|value| value - mean).collect();
        let sigma2 =
            centered.iter().map(|value| value * value).sum::<f64>() / centered.len() as f64;

        let start = Instant::now();
        let direct = statistical_inefficiency_direct(&centered, sigma2, false).unwrap();
        let direct_elapsed = start.elapsed();

        let start = Instant::now();
        let fft = statistical_inefficiency_fft(&centered, sigma2, false).unwrap();
        let fft_elapsed = start.elapsed();

        eprintln!(
            "direct={direct_elapsed:?} fft={fft_elapsed:?} speedup={:.2}x",
            direct_elapsed.as_secs_f64() / fft_elapsed.as_secs_f64()
        );
        assert!((fft - direct).abs() < 1.0e-9, "{fft} vs {direct}");
    }
}
