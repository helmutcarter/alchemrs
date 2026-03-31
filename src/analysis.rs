use crate::data::{DhdlSeries, DeltaFMatrix, FreeEnergyEstimate, OverlapMatrix, StatePoint, UNkMatrix};
use crate::error::{CoreError, Result};
use crate::estimators::{
    mbar_log_weights_from_windows, BarEstimator, BarOptions, ExpEstimator, ExpOptions,
    MbarEstimator, MbarOptions, TiEstimator, TiOptions,
};
use nalgebra::{DMatrix, Schur};

#[derive(Debug, Clone, PartialEq)]
pub struct ConvergencePoint {
    n_windows: usize,
    delta_f: f64,
    uncertainty: Option<f64>,
    from_state: StatePoint,
    to_state: StatePoint,
    lambda_labels: Option<Vec<String>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlockEstimate {
    block_index: usize,
    n_blocks: usize,
    delta_f: f64,
    uncertainty: Option<f64>,
    from_state: StatePoint,
    to_state: StatePoint,
    lambda_labels: Option<Vec<String>>,
}

impl BlockEstimate {
    pub fn new(
        block_index: usize,
        n_blocks: usize,
        delta_f: f64,
        uncertainty: Option<f64>,
        from_state: StatePoint,
        to_state: StatePoint,
        lambda_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        if n_blocks == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        if block_index >= n_blocks {
            return Err(CoreError::InvalidShape {
                expected: n_blocks,
                found: block_index + 1,
            });
        }
        if !delta_f.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "delta_f must be finite".to_string(),
            ));
        }
        if let Some(value) = uncertainty {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "uncertainty must be finite".to_string(),
                ));
            }
        }
        Ok(Self {
            block_index,
            n_blocks,
            delta_f,
            uncertainty,
            from_state,
            to_state,
            lambda_labels,
        })
    }

    pub fn block_index(&self) -> usize {
        self.block_index
    }

    pub fn n_blocks(&self) -> usize {
        self.n_blocks
    }

    pub fn delta_f(&self) -> f64 {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

impl ConvergencePoint {
    pub fn new(
        n_windows: usize,
        delta_f: f64,
        uncertainty: Option<f64>,
        from_state: StatePoint,
        to_state: StatePoint,
        lambda_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        if n_windows == 0 {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }
        if !delta_f.is_finite() {
            return Err(CoreError::NonFiniteValue(
                "delta_f must be finite".to_string(),
            ));
        }
        if let Some(value) = uncertainty {
            if !value.is_finite() {
                return Err(CoreError::NonFiniteValue(
                    "uncertainty must be finite".to_string(),
                ));
            }
        }
        Ok(Self {
            n_windows,
            delta_f,
            uncertainty,
            from_state,
            to_state,
            lambda_labels,
        })
    }

    pub fn n_windows(&self) -> usize {
        self.n_windows
    }

    pub fn delta_f(&self) -> f64 {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn from_state(&self) -> &StatePoint {
        &self.from_state
    }

    pub fn to_state(&self) -> &StatePoint {
        &self.to_state
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.lambda_labels.as_deref()
    }
}

pub fn overlap_matrix(
    windows: &[UNkMatrix],
    options: Option<MbarOptions>,
) -> Result<OverlapMatrix> {
    let options = options.unwrap_or_default();
    let lambda_labels = windows
        .first()
        .and_then(|window| window.lambda_labels().map(|labels| labels.to_vec()));
    let (log_w, n_k, states) = mbar_log_weights_from_windows(windows, &options)?;
    let n_states = states.len();
    if log_w.len() % n_states != 0 {
        return Err(CoreError::InvalidShape {
            expected: n_states,
            found: log_w.len(),
        });
    }
    let n_samples = log_w.len() / n_states;

    let mut wtw = vec![0.0; n_states * n_states];
    let mut weights = vec![0.0; n_states];
    for n in 0..n_samples {
        let base = n * n_states;
        for i in 0..n_states {
            weights[i] = log_w[base + i].exp();
        }
        for i in 0..n_states {
            let wi = weights[i];
            for j in 0..n_states {
                wtw[i * n_states + j] += wi * weights[j];
            }
        }
    }

    let mut values = vec![0.0; n_states * n_states];
    for i in 0..n_states {
        for j in 0..n_states {
            values[i * n_states + j] = wtw[i * n_states + j] * n_k[j];
        }
    }

    OverlapMatrix::new_with_labels(values, n_states, states, lambda_labels)
}

pub fn overlap_eigenvalues(overlap: &OverlapMatrix) -> Result<Vec<f64>> {
    let n_states = overlap.n_states();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let matrix = DMatrix::from_row_slice(n_states, n_states, overlap.values());
    let eigenvalues = Schur::new(matrix).complex_eigenvalues();
    let mut values = Vec::with_capacity(n_states);
    for eig in eigenvalues.iter() {
        if eig.im.abs() > 1e-8 {
            return Err(CoreError::InvalidState(
                "overlap eigenvalue has significant imaginary component".to_string(),
            ));
        }
        values.push(eig.re);
    }
    values.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    Ok(values)
}

pub fn overlap_scalar(overlap: &OverlapMatrix) -> Result<f64> {
    if overlap.n_states() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: overlap.n_states(),
        });
    }
    let eigenvalues = overlap_eigenvalues(overlap)?;
    Ok(1.0 - eigenvalues[1])
}

pub fn ti_convergence(series: &[DhdlSeries], options: Option<TiOptions>) -> Result<Vec<ConvergencePoint>> {
    if series.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: series.len(),
        });
    }
    let estimator = TiEstimator::new(options.unwrap_or_default());
    let mut points = Vec::with_capacity(series.len() - 1);
    for end in 2..=series.len() {
        let result = estimator.fit(&series[..end])?;
        points.push(convergence_point_from_scalar(end, &result)?);
    }
    Ok(points)
}

pub fn bar_convergence(windows: &[UNkMatrix], options: Option<BarOptions>) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| BarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub fn mbar_convergence(
    windows: &[UNkMatrix],
    options: Option<MbarOptions>,
) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| MbarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        1,
    )
}

pub fn exp_convergence(windows: &[UNkMatrix], options: Option<ExpOptions>) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub fn dexp_convergence(
    windows: &[UNkMatrix],
    options: Option<ExpOptions>,
) -> Result<Vec<ConvergencePoint>> {
    convergence_from_matrix_windows(
        windows,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        true,
        2,
    )
}

pub(crate) fn ti_block_average(
    series: &[DhdlSeries],
    n_blocks: usize,
    options: Option<TiOptions>,
) -> Result<Vec<BlockEstimate>> {
    if series.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: series.len(),
        });
    }
    if n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let estimator = TiEstimator::new(options.unwrap_or_default());
    let blocked = series
        .iter()
        .map(|item| split_dhdl_series(item, n_blocks))
        .collect::<Result<Vec<_>>>()?;

    let mut points = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let block = blocked
            .iter()
            .map(|chunks| chunks[block_index].clone())
            .collect::<Vec<_>>();
        let result = estimator.fit(&block)?;
        points.push(block_estimate_from_scalar(block_index, n_blocks, &result)?);
    }
    Ok(points)
}

pub(crate) fn mbar_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<MbarOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| MbarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        1,
    )
}

pub(crate) fn bar_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<BarOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| BarEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub(crate) fn exp_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<ExpOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        false,
        2,
    )
}

pub(crate) fn dexp_block_average(
    windows: &[UNkMatrix],
    n_blocks: usize,
    options: Option<ExpOptions>,
) -> Result<Vec<BlockEstimate>> {
    block_average_from_windows(
        windows,
        n_blocks,
        |subset| ExpEstimator::new(options.clone().unwrap_or_default()).fit(subset),
        true,
        2,
    )
}

fn convergence_from_matrix_windows<F>(
    windows: &[UNkMatrix],
    fit: F,
    reverse: bool,
    minimum_windows: usize,
) -> Result<Vec<ConvergencePoint>>
where
    F: Fn(&[UNkMatrix]) -> Result<DeltaFMatrix>,
{
    if windows.len() < minimum_windows {
        return Err(CoreError::InvalidShape {
            expected: minimum_windows,
            found: windows.len(),
        });
    }
    let mut points = Vec::with_capacity(windows.len() - minimum_windows + 1);
    for count in minimum_windows..=windows.len() {
        let subset = if reverse {
            &windows[windows.len() - count..]
        } else {
            &windows[..count]
        };
        let trimmed = trim_windows_to_sampled_states(subset)?;
        let result = fit(&trimmed)?;
        points.push(convergence_point_from_matrix(count, &result, reverse)?);
    }
    Ok(points)
}

fn block_average_from_windows<F>(
    windows: &[UNkMatrix],
    n_blocks: usize,
    fit: F,
    reverse: bool,
    minimum_windows: usize,
) -> Result<Vec<BlockEstimate>>
where
    F: Fn(&[UNkMatrix]) -> Result<DeltaFMatrix>,
{
    if windows.len() < minimum_windows {
        return Err(CoreError::InvalidShape {
            expected: minimum_windows,
            found: windows.len(),
        });
    }
    if n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let blocked = windows
        .iter()
        .map(|window| split_u_nk_window(window, n_blocks))
        .collect::<Result<Vec<_>>>()?;

    let mut points = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let subset = blocked
            .iter()
            .map(|chunks| chunks[block_index].clone())
            .collect::<Vec<_>>();
        let trimmed = trim_windows_to_sampled_states(&subset)?;
        let result = fit(&trimmed)?;
        points.push(block_estimate_from_matrix(block_index, n_blocks, &result, reverse)?);
    }
    Ok(points)
}

fn trim_windows_to_sampled_states(windows: &[UNkMatrix]) -> Result<Vec<UNkMatrix>> {
    let states = windows
        .iter()
        .map(|window| {
            window.sampled_state().cloned().ok_or_else(|| {
                CoreError::InvalidState(
                    "sampled_state required for convergence analysis".to_string(),
                )
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut trimmed = Vec::with_capacity(windows.len());
    for window in windows {
        let indices = states
            .iter()
            .map(|state| {
                window
                    .evaluated_states()
                    .iter()
                    .position(|candidate| candidate == state)
                    .ok_or_else(|| {
                        CoreError::InvalidState(
                            "sampled_state not found in evaluated_states".to_string(),
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;

        let mut data = Vec::with_capacity(window.n_samples() * states.len());
        for sample_idx in 0..window.n_samples() {
            let row_offset = sample_idx * window.n_states();
            for &idx in &indices {
                data.push(window.data()[row_offset + idx]);
            }
        }

        trimmed.push(UNkMatrix::new_with_labels(
            window.n_samples(),
            states.len(),
            data,
            window.time_ps().to_vec(),
            window.sampled_state().cloned(),
            states.clone(),
            window.lambda_labels().map(|labels| labels.to_vec()),
        )?);
    }

    Ok(trimmed)
}

fn split_dhdl_series(series: &DhdlSeries, n_blocks: usize) -> Result<Vec<DhdlSeries>> {
    let block_len = block_length(series.values().len(), n_blocks)?;
    let mut out = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let start = block_index * block_len;
        let end = start + block_len;
        out.push(DhdlSeries::new(
            series.state().clone(),
            series.time_ps()[start..end].to_vec(),
            series.values()[start..end].to_vec(),
        )?);
    }
    Ok(out)
}

fn split_u_nk_window(window: &UNkMatrix, n_blocks: usize) -> Result<Vec<UNkMatrix>> {
    let block_len = block_length(window.n_samples(), n_blocks)?;
    let row_width = window.n_states();
    let mut out = Vec::with_capacity(n_blocks);
    for block_index in 0..n_blocks {
        let start = block_index * block_len;
        let end = start + block_len;
        let mut data = Vec::with_capacity(block_len * row_width);
        for sample_idx in start..end {
            let row_start = sample_idx * row_width;
            data.extend_from_slice(&window.data()[row_start..row_start + row_width]);
        }
        out.push(UNkMatrix::new_with_labels(
            block_len,
            row_width,
            data,
            window.time_ps()[start..end].to_vec(),
            window.sampled_state().cloned(),
            window.evaluated_states().to_vec(),
            window.lambda_labels().map(|labels| labels.to_vec()),
        )?);
    }
    Ok(out)
}

fn block_length(n_samples: usize, n_blocks: usize) -> Result<usize> {
    if n_blocks == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let block_len = n_samples / n_blocks;
    if block_len == 0 {
        return Err(CoreError::InvalidShape {
            expected: n_blocks,
            found: n_samples,
        });
    }
    Ok(block_len)
}

fn convergence_point_from_scalar(n_windows: usize, result: &FreeEnergyEstimate) -> Result<ConvergencePoint> {
    ConvergencePoint::new(
        n_windows,
        result.delta_f(),
        result.uncertainty(),
        result.from_state().clone(),
        result.to_state().clone(),
        None,
    )
}

fn block_estimate_from_scalar(
    block_index: usize,
    n_blocks: usize,
    result: &FreeEnergyEstimate,
) -> Result<BlockEstimate> {
    BlockEstimate::new(
        block_index,
        n_blocks,
        result.delta_f(),
        result.uncertainty(),
        result.from_state().clone(),
        result.to_state().clone(),
        None,
    )
}

fn convergence_point_from_matrix(
    n_windows: usize,
    result: &DeltaFMatrix,
    reverse: bool,
) -> Result<ConvergencePoint> {
    let n_states = result.n_states();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let index = if reverse {
        (n_states - 1) * n_states
    } else {
        n_states - 1
    };
    let (from_state, to_state) = if reverse {
        (
            result.states().last().unwrap().clone(),
            result.states().first().unwrap().clone(),
        )
    } else {
        (
            result.states().first().unwrap().clone(),
            result.states().last().unwrap().clone(),
        )
    };
    ConvergencePoint::new(
        n_windows,
        result.values()[index],
        result.uncertainties().map(|values| values[index]),
        from_state,
        to_state,
        result.lambda_labels().map(|labels| labels.to_vec()),
    )
}

fn block_estimate_from_matrix(
    block_index: usize,
    n_blocks: usize,
    result: &DeltaFMatrix,
    reverse: bool,
) -> Result<BlockEstimate> {
    let n_states = result.n_states();
    if n_states == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    let index = if reverse {
        (n_states - 1) * n_states
    } else {
        n_states - 1
    };
    BlockEstimate::new(
        block_index,
        n_blocks,
        result.values()[index],
        result.uncertainties().map(|values| values[index]),
        if reverse {
            result.states().last().unwrap().clone()
        } else {
            result.states().first().unwrap().clone()
        },
        if reverse {
            result.states().first().unwrap().clone()
        } else {
            result.states().last().unwrap().clone()
        },
        result.lambda_labels().map(|labels| labels.to_vec()),
    )
}
