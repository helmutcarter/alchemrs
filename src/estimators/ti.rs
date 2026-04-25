use crate::analysis::{self, BlockEstimate};
use crate::data::{DhdlSeries, FreeEnergyEstimate, StatePoint};
use crate::error::{CoreError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrationMethod {
    Trapezoidal,
    Simpson,
    CubicSpline,
    Pchip,
    Akima,
    GaussianQuadrature,
}

impl IntegrationMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Trapezoidal => "trapezoidal",
            Self::Simpson => "simpson",
            Self::CubicSpline => "cubic-spline",
            Self::Pchip => "pchip",
            Self::Akima => "akima",
            Self::GaussianQuadrature => "gaussian-quadrature",
        }
    }
}

const GAUSSIAN_QUADRATURE_LAMBDA_ABS_TOL: f64 = 1.0e-5;
const NONLINEAR_TI_UNCERTAINTY_MIN_STEP: f64 = 1.0e-6;

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

#[derive(Debug, Clone)]
pub struct TiFit {
    delta_f: f64,
    uncertainty: Option<f64>,
    from_state: StatePoint,
    to_state: StatePoint,
    method: IntegrationMethod,
}

impl TiEstimator {
    pub fn new(options: TiOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, series: &[DhdlSeries]) -> Result<TiFit> {
        let min_windows = minimum_ti_windows(self.options.method);
        if series.len() < min_windows {
            return Err(CoreError::InvalidShape {
                expected: min_windows,
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
        validate_strictly_increasing_lambdas(&lambdas)?;

        let delta_f = match self.options.method {
            IntegrationMethod::Trapezoidal => integrate_trapezoidal(&lambdas, &values)?,
            IntegrationMethod::Simpson => integrate_simpson(&lambdas, &values)?,
            IntegrationMethod::CubicSpline => integrate_cubic_spline(&lambdas, &values)?,
            IntegrationMethod::Pchip => integrate_pchip(&lambdas, &values)?,
            IntegrationMethod::Akima => integrate_akima(&lambdas, &values)?,
            IntegrationMethod::GaussianQuadrature => {
                integrate_gaussian_quadrature(&lambdas, &values)?
            }
        };

        let uncertainty = match self.options.method {
            IntegrationMethod::Trapezoidal => {
                Some(trapezoidal_uncertainty(&lambdas, &sem2_values)?)
            }
            IntegrationMethod::Simpson => Some(simpson_uncertainty(&lambdas, &sem2_values)?),
            IntegrationMethod::CubicSpline => {
                Some(cubic_spline_uncertainty(&lambdas, &sem2_values)?)
            }
            IntegrationMethod::Pchip => Some(pchip_uncertainty(&lambdas, &values, &sem2_values)?),
            IntegrationMethod::Akima => Some(akima_uncertainty(&lambdas, &values, &sem2_values)?),
            IntegrationMethod::GaussianQuadrature => {
                Some(gaussian_quadrature_uncertainty(&lambdas, &sem2_values)?)
            }
        };

        let (from_state, to_state) = match self.options.method {
            IntegrationMethod::GaussianQuadrature => {
                // Gauss-Legendre TI integrates the full [0, 1] interval from interior nodes,
                // so the reported endpoints are synthetic lambda=0 and lambda=1 states.
                gaussian_quadrature_endpoint_states(points.first().unwrap().3.temperature_k())?
            }
            _ => (
                points.first().unwrap().3.clone(),
                points.last().unwrap().3.clone(),
            ),
        };
        Ok(TiFit {
            delta_f,
            uncertainty,
            from_state,
            to_state,
            method: self.options.method,
        })
    }

    pub fn estimate(&self, series: &[DhdlSeries]) -> Result<FreeEnergyEstimate> {
        self.fit(series)?.result()
    }

    pub fn block_average(
        &self,
        series: &[DhdlSeries],
        n_blocks: usize,
    ) -> Result<Vec<BlockEstimate>> {
        analysis::ti_block_average(series, n_blocks, Some(self.options.clone()))
    }
}

pub fn sample_ti_curve(
    lambdas: &[f64],
    values: &[f64],
    method: IntegrationMethod,
    samples_per_interval: usize,
) -> Result<Vec<(f64, f64)>> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if samples_per_interval == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    match method {
        IntegrationMethod::Trapezoidal => {
            sample_trapezoidal_curve(lambdas, values, samples_per_interval)
        }
        IntegrationMethod::Simpson => sample_simpson_curve(lambdas, values, samples_per_interval),
        IntegrationMethod::CubicSpline => {
            let second = natural_cubic_second_derivatives(lambdas, values)?;
            sample_cubic_spline_curve(lambdas, values, &second, samples_per_interval)
        }
        IntegrationMethod::Pchip => {
            let derivatives = pchip_derivatives(lambdas, values)?;
            sample_cubic_hermite_curve(lambdas, values, &derivatives, samples_per_interval, "PCHIP")
        }
        IntegrationMethod::Akima => {
            let derivatives = akima_derivatives(lambdas, values)?;
            sample_cubic_hermite_curve(lambdas, values, &derivatives, samples_per_interval, "Akima")
        }
        IntegrationMethod::GaussianQuadrature => {
            gaussian_quadrature_rule(lambdas)?;
            Ok(lambdas
                .iter()
                .copied()
                .zip(values.iter().copied())
                .collect())
        }
    }
}

impl TiFit {
    pub fn method(&self) -> IntegrationMethod {
        self.method
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

    pub fn result(&self) -> Result<FreeEnergyEstimate> {
        FreeEnergyEstimate::new(
            self.delta_f,
            self.uncertainty,
            self.from_state.clone(),
            self.to_state.clone(),
        )
    }
}

fn minimum_ti_windows(method: IntegrationMethod) -> usize {
    match method {
        IntegrationMethod::Trapezoidal => 2,
        IntegrationMethod::Simpson => 3,
        IntegrationMethod::CubicSpline => 2,
        IntegrationMethod::Pchip => 2,
        IntegrationMethod::Akima => 2,
        IntegrationMethod::GaussianQuadrature => 1,
    }
}

fn extract_lambda(state: &StatePoint) -> Result<f64> {
    let lambdas = state.lambdas();
    if lambdas.len() != 1 {
        return Err(CoreError::RequiresOneDimensionalLambda {
            operation: "estimators",
        });
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

fn validate_strictly_increasing_lambdas(lambdas: &[f64]) -> Result<()> {
    for pair in lambdas.windows(2) {
        if pair[1] <= pair[0] {
            return Err(CoreError::Unsupported(
                "TI integration requires strictly increasing lambda spacing".to_string(),
            ));
        }
    }
    Ok(())
}

fn integrate_trapezoidal(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    if lambdas.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: lambdas.len(),
        });
    }
    let mut total = 0.0;
    for i in 0..(lambdas.len() - 1) {
        let dx = lambdas[i + 1] - lambdas[i];
        total += dx * (values[i] + values[i + 1]) * 0.5;
    }
    Ok(total)
}

fn integrate_simpson(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    if lambdas.len() < 3 {
        return Err(CoreError::InvalidShape {
            expected: 3,
            found: lambdas.len(),
        });
    }
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

fn integrate_gaussian_quadrature(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    let rule = gaussian_quadrature_rule(lambdas)?;
    Ok(rule
        .weights
        .iter()
        .zip(values.iter())
        .map(|(weight, value)| weight * value)
        .sum())
}

fn sample_trapezoidal_curve(
    lambdas: &[f64],
    values: &[f64],
    samples_per_interval: usize,
) -> Result<Vec<(f64, f64)>> {
    interval_widths_and_slopes(lambdas, values, "trapezoidal")?;
    let mut points = Vec::with_capacity((lambdas.len() - 1) * samples_per_interval + 1);
    for interval_index in 0..(lambdas.len() - 1) {
        let x0 = lambdas[interval_index];
        let x1 = lambdas[interval_index + 1];
        let y0 = values[interval_index];
        let y1 = values[interval_index + 1];
        let steps = samples_per_interval;
        for step in 0..=steps {
            if interval_index > 0 && step == 0 {
                continue;
            }
            let t = step as f64 / (steps as f64);
            let x = x0 + (x1 - x0) * t;
            let y = y0 + (y1 - y0) * t;
            points.push((x, y));
        }
    }
    Ok(points)
}

fn sample_simpson_curve(
    lambdas: &[f64],
    values: &[f64],
    samples_per_interval: usize,
) -> Result<Vec<(f64, f64)>> {
    simpson_weights(lambdas)?;
    let mut points =
        Vec::with_capacity(((lambdas.len() - 1) * samples_per_interval).saturating_add(1));
    for start in (0..(lambdas.len() - 2)).step_by(2) {
        let x0 = lambdas[start];
        let x1 = lambdas[start + 1];
        let x2 = lambdas[start + 2];
        let y0 = values[start];
        let y1 = values[start + 1];
        let y2 = values[start + 2];
        let steps = samples_per_interval * 2;
        for step in 0..=steps {
            if start > 0 && step == 0 {
                continue;
            }
            let t = step as f64 / (steps as f64);
            let x = x0 + (x2 - x0) * t;
            let y = lagrange_quadratic(x0, y0, x1, y1, x2, y2, x);
            points.push((x, y));
        }
    }
    Ok(points)
}

fn simpson_uncertainty(lambdas: &[f64], sem2: &[f64]) -> Result<f64> {
    if lambdas.len() != sem2.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: sem2.len(),
        });
    }

    let weights = simpson_weights(lambdas)?;
    let variance: f64 = weights
        .iter()
        .zip(sem2.iter())
        .map(|(weight, sem2)| weight * weight * sem2)
        .sum();
    Ok(variance.sqrt())
}

fn integrate_cubic_spline(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if lambdas.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: lambdas.len(),
        });
    }

    // Natural cubic splines set the endpoint second derivatives to zero and
    // integrate the resulting piecewise cubic exactly on each interval.
    let second_derivatives = natural_cubic_second_derivatives(lambdas, values)?;
    let mut total = 0.0;
    for i in 0..(lambdas.len() - 1) {
        let h = lambdas[i + 1] - lambdas[i];
        total += h * (values[i] + values[i + 1]) * 0.5
            - (h * h * h) * (second_derivatives[i] + second_derivatives[i + 1]) / 24.0;
    }
    Ok(total)
}

fn sample_cubic_spline_curve(
    lambdas: &[f64],
    values: &[f64],
    second_derivatives: &[f64],
    samples_per_interval: usize,
) -> Result<Vec<(f64, f64)>> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if lambdas.len() != second_derivatives.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: second_derivatives.len(),
        });
    }
    let mut points = Vec::with_capacity((lambdas.len() - 1) * samples_per_interval + 1);
    for interval_index in 0..(lambdas.len() - 1) {
        let x0 = lambdas[interval_index];
        let x1 = lambdas[interval_index + 1];
        let y0 = values[interval_index];
        let y1 = values[interval_index + 1];
        let m0 = second_derivatives[interval_index];
        let m1 = second_derivatives[interval_index + 1];
        let h = x1 - x0;
        if h <= 0.0 {
            return Err(CoreError::Unsupported(
                "Cubic spline integration requires strictly increasing lambda spacing".to_string(),
            ));
        }
        for step in 0..=samples_per_interval {
            if interval_index > 0 && step == 0 {
                continue;
            }
            let t = step as f64 / (samples_per_interval as f64);
            let x = x0 + h * t;
            let a = (x1 - x) / h;
            let b = (x - x0) / h;
            let y =
                a * y0 + b * y1 + (((a * a * a) - a) * m0 + ((b * b * b) - b) * m1) * (h * h) / 6.0;
            points.push((x, y));
        }
    }
    Ok(points)
}

fn integrate_pchip(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    let derivatives = pchip_derivatives(lambdas, values)?;
    integrate_cubic_hermite(lambdas, values, &derivatives, "PCHIP")
}

fn integrate_akima(lambdas: &[f64], values: &[f64]) -> Result<f64> {
    let derivatives = akima_derivatives(lambdas, values)?;
    integrate_cubic_hermite(lambdas, values, &derivatives, "Akima")
}

fn cubic_spline_uncertainty(lambdas: &[f64], sem2: &[f64]) -> Result<f64> {
    if lambdas.len() != sem2.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: sem2.len(),
        });
    }

    // For a fixed lambda grid, the natural cubic-spline integral is a linear
    // functional of the window means. Reusing those implied weights lets us
    // propagate the per-window SEMs with the same independence assumption used
    // by the other TI uncertainty estimates.
    let weights = cubic_spline_weights(lambdas)?;
    let variance: f64 = weights
        .iter()
        .zip(sem2.iter())
        .map(|(weight, sem2)| weight * weight * sem2)
        .sum();
    Ok(variance.sqrt())
}

fn pchip_uncertainty(lambdas: &[f64], values: &[f64], sem2: &[f64]) -> Result<f64> {
    nonlinear_ti_uncertainty(lambdas, values, sem2, "PCHIP", integrate_pchip)
}

fn akima_uncertainty(lambdas: &[f64], values: &[f64], sem2: &[f64]) -> Result<f64> {
    nonlinear_ti_uncertainty(lambdas, values, sem2, "Akima", integrate_akima)
}

fn integrate_cubic_hermite(
    lambdas: &[f64],
    values: &[f64],
    derivatives: &[f64],
    method_name: &str,
) -> Result<f64> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if lambdas.len() != derivatives.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: derivatives.len(),
        });
    }

    let mut total = 0.0;
    for i in 0..(lambdas.len() - 1) {
        let h = lambdas[i + 1] - lambdas[i];
        if h <= 0.0 {
            return Err(CoreError::Unsupported(format!(
                "{method_name} integration requires strictly increasing lambda spacing"
            )));
        }
        total += h * (values[i] + values[i + 1]) * 0.5
            + (h * h) * (derivatives[i] - derivatives[i + 1]) / 12.0;
    }
    Ok(total)
}

fn sample_cubic_hermite_curve(
    lambdas: &[f64],
    values: &[f64],
    derivatives: &[f64],
    samples_per_interval: usize,
    method_name: &str,
) -> Result<Vec<(f64, f64)>> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if lambdas.len() != derivatives.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: derivatives.len(),
        });
    }
    let mut points = Vec::with_capacity((lambdas.len() - 1) * samples_per_interval + 1);
    for interval_index in 0..(lambdas.len() - 1) {
        let x0 = lambdas[interval_index];
        let x1 = lambdas[interval_index + 1];
        let y0 = values[interval_index];
        let y1 = values[interval_index + 1];
        let d0 = derivatives[interval_index];
        let d1 = derivatives[interval_index + 1];
        let h = x1 - x0;
        if h <= 0.0 {
            return Err(CoreError::Unsupported(format!(
                "{method_name} integration requires strictly increasing lambda spacing"
            )));
        }
        for step in 0..=samples_per_interval {
            if interval_index > 0 && step == 0 {
                continue;
            }
            let t = step as f64 / (samples_per_interval as f64);
            let x = x0 + h * t;
            let h00 = 2.0 * t * t * t - 3.0 * t * t + 1.0;
            let h10 = t * t * t - 2.0 * t * t + t;
            let h01 = -2.0 * t * t * t + 3.0 * t * t;
            let h11 = t * t * t - t * t;
            let y = h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1;
            points.push((x, y));
        }
    }
    Ok(points)
}

fn lagrange_quadratic(x0: f64, y0: f64, x1: f64, y1: f64, x2: f64, y2: f64, x: f64) -> f64 {
    let l0 = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2));
    let l1 = ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2));
    let l2 = ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1));
    y0 * l0 + y1 * l1 + y2 * l2
}

fn gaussian_quadrature_uncertainty(lambdas: &[f64], sem2: &[f64]) -> Result<f64> {
    if lambdas.len() != sem2.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: sem2.len(),
        });
    }
    let rule = gaussian_quadrature_rule(lambdas)?;
    let variance: f64 = rule
        .weights
        .iter()
        .zip(sem2.iter())
        .map(|(weight, sem2)| weight * weight * sem2)
        .sum();
    Ok(variance.sqrt())
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

fn natural_cubic_second_derivatives(lambdas: &[f64], values: &[f64]) -> Result<Vec<f64>> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if lambdas.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: lambdas.len(),
        });
    }

    let n = lambdas.len();
    let mut intervals = Vec::with_capacity(n - 1);
    for i in 0..(n - 1) {
        let width = lambdas[i + 1] - lambdas[i];
        if width <= 0.0 {
            return Err(CoreError::Unsupported(
                "Cubic spline integration requires strictly increasing lambda spacing".to_string(),
            ));
        }
        intervals.push(width);
    }

    if n == 2 {
        return Ok(vec![0.0, 0.0]);
    }

    let interior = n - 2;
    let mut lower = vec![0.0; interior];
    let mut diag = vec![0.0; interior];
    let mut upper = vec![0.0; interior];
    let mut rhs = vec![0.0; interior];

    for row in 0..interior {
        let i = row + 1;
        let h_prev = intervals[i - 1];
        let h_next = intervals[i];
        lower[row] = if row == 0 { 0.0 } else { h_prev };
        diag[row] = 2.0 * (h_prev + h_next);
        upper[row] = if row + 1 == interior { 0.0 } else { h_next };
        rhs[row] =
            6.0 * (((values[i + 1] - values[i]) / h_next) - ((values[i] - values[i - 1]) / h_prev));
    }

    for row in 1..interior {
        let factor = lower[row] / diag[row - 1];
        diag[row] -= factor * upper[row - 1];
        rhs[row] -= factor * rhs[row - 1];
    }

    let mut second = vec![0.0; n];
    second[n - 2] = rhs[interior - 1] / diag[interior - 1];
    for row in (0..(interior - 1)).rev() {
        second[row + 1] = (rhs[row] - upper[row] * second[row + 2]) / diag[row];
    }
    Ok(second)
}

fn pchip_derivatives(lambdas: &[f64], values: &[f64]) -> Result<Vec<f64>> {
    let (intervals, slopes) = interval_widths_and_slopes(lambdas, values, "PCHIP")?;
    if values.len() == 2 {
        return Ok(vec![slopes[0], slopes[0]]);
    }

    let mut derivatives = vec![0.0; values.len()];
    for k in 1..(values.len() - 1) {
        let left = slopes[k - 1];
        let right = slopes[k];
        if left == 0.0 || right == 0.0 || left.signum() != right.signum() {
            derivatives[k] = 0.0;
        } else {
            let w1 = 2.0 * intervals[k] + intervals[k - 1];
            let w2 = intervals[k] + 2.0 * intervals[k - 1];
            derivatives[k] = (w1 + w2) / (w1 / left + w2 / right);
        }
    }

    derivatives[0] = pchip_edge_derivative(intervals[0], intervals[1], slopes[0], slopes[1]);
    let last = values.len() - 1;
    derivatives[last] = pchip_edge_derivative(
        intervals[last - 1],
        intervals[last - 2],
        slopes[last - 1],
        slopes[last - 2],
    );
    Ok(derivatives)
}

fn pchip_edge_derivative(h0: f64, h1: f64, m0: f64, m1: f64) -> f64 {
    let derivative = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
    if derivative.signum() != m0.signum() {
        0.0
    } else if m0.signum() != m1.signum() && derivative.abs() > 3.0 * m0.abs() {
        3.0 * m0
    } else {
        derivative
    }
}

fn akima_derivatives(lambdas: &[f64], values: &[f64]) -> Result<Vec<f64>> {
    let (_intervals, slopes) = interval_widths_and_slopes(lambdas, values, "Akima")?;
    if values.len() == 2 {
        return Ok(vec![slopes[0], slopes[0]]);
    }

    let n = values.len();
    let mut extended = vec![0.0; n + 3];
    for (idx, slope) in slopes.iter().enumerate() {
        extended[idx + 2] = *slope;
    }
    extended[1] = 2.0 * extended[2] - extended[3];
    extended[0] = 2.0 * extended[1] - extended[2];
    extended[n + 1] = 2.0 * extended[n] - extended[n - 1];
    extended[n + 2] = 2.0 * extended[n + 1] - extended[n];

    let mut derivatives = (0..n)
        .map(|idx| 0.5 * (extended[idx + 3] + extended[idx]))
        .collect::<Vec<_>>();
    let dm = extended
        .windows(2)
        .map(|pair| (pair[1] - pair[0]).abs())
        .collect::<Vec<_>>();
    let f1 = &dm[2..];
    let f2 = &dm[..n];
    let break_mult = 1.0e-9;
    let f12_max = f1
        .iter()
        .zip(f2.iter())
        .map(|(left, right)| left + right)
        .fold(f64::NEG_INFINITY, f64::max);

    for idx in 0..n {
        let f12 = f1[idx] + f2[idx];
        if f12 > break_mult * f12_max {
            derivatives[idx] =
                extended[idx + 1] + (f2[idx] / f12) * (extended[idx + 2] - extended[idx + 1]);
        }
    }

    Ok(derivatives)
}

fn cubic_spline_weights(lambdas: &[f64]) -> Result<Vec<f64>> {
    if lambdas.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: lambdas.len(),
        });
    }

    let mut weights = Vec::with_capacity(lambdas.len());
    let mut basis = vec![0.0; lambdas.len()];
    for basis_index in 0..lambdas.len() {
        basis.fill(0.0);
        basis[basis_index] = 1.0;
        weights.push(integrate_cubic_spline(lambdas, &basis)?);
    }
    Ok(weights)
}

fn simpson_weights(lambdas: &[f64]) -> Result<Vec<f64>> {
    if lambdas.len() < 3 {
        return Err(CoreError::InvalidShape {
            expected: 3,
            found: lambdas.len(),
        });
    }
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
    let mut weights = Vec::with_capacity(n);
    for i in 0..n {
        let expected = lambdas[0] + (i as f64) * h;
        if (lambdas[i] - expected).abs() > tol {
            return Err(CoreError::Unsupported(
                "Simpson integration requires uniform lambda spacing".to_string(),
            ));
        }
        let coeff = if i == 0 || i + 1 == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        weights.push(coeff * h / 3.0);
    }
    Ok(weights)
}

fn nonlinear_ti_uncertainty(
    lambdas: &[f64],
    values: &[f64],
    sem2: &[f64],
    method_name: &str,
    integrator: fn(&[f64], &[f64]) -> Result<f64>,
) -> Result<f64> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if values.len() != sem2.len() {
        return Err(CoreError::InvalidShape {
            expected: values.len(),
            found: sem2.len(),
        });
    }

    // PCHIP and Akima depend nonlinearly on the window means, so propagate the
    // SEMs through a numerical Jacobian of the integrated estimate instead of
    // reusing the linear-weight shortcut used by trapezoidal/GQ/natural spline TI.
    let mut plus = values.to_vec();
    let mut minus = values.to_vec();
    let mut variance = 0.0;
    for index in 0..values.len() {
        let step = (f64::EPSILON.sqrt() * (1.0 + values[index].abs()))
            .max(NONLINEAR_TI_UNCERTAINTY_MIN_STEP);
        plus[index] = values[index] + step;
        minus[index] = values[index] - step;

        let upper = integrator(lambdas, &plus)?;
        let lower = integrator(lambdas, &minus)?;
        let derivative = (upper - lower) / (2.0 * step);

        if !derivative.is_finite() {
            return Err(CoreError::NonFiniteValue(format!(
                "{method_name} uncertainty Jacobian must be finite"
            )));
        }
        variance += derivative * derivative * sem2[index];
        plus[index] = values[index];
        minus[index] = values[index];
    }

    Ok(variance.sqrt())
}

fn interval_widths_and_slopes(
    lambdas: &[f64],
    values: &[f64],
    method_name: &str,
) -> Result<(Vec<f64>, Vec<f64>)> {
    if lambdas.len() != values.len() {
        return Err(CoreError::InvalidShape {
            expected: lambdas.len(),
            found: values.len(),
        });
    }
    if lambdas.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: lambdas.len(),
        });
    }

    let mut intervals = Vec::with_capacity(lambdas.len() - 1);
    let mut slopes = Vec::with_capacity(lambdas.len() - 1);
    for i in 0..(lambdas.len() - 1) {
        let width = lambdas[i + 1] - lambdas[i];
        if width <= 0.0 {
            return Err(CoreError::Unsupported(format!(
                "{method_name} integration requires strictly increasing lambda spacing"
            )));
        }
        intervals.push(width);
        slopes.push((values[i + 1] - values[i]) / width);
    }
    Ok((intervals, slopes))
}

fn gaussian_quadrature_rule(lambdas: &[f64]) -> Result<&'static QuadratureRule> {
    let Some(rule) = gaussian_quadrature_rule_by_count(lambdas.len()) else {
        return Err(CoreError::Unsupported(format!(
            "Gaussian quadrature currently supports 1 through 16 lambda windows; found {}",
            lambdas.len()
        )));
    };

    for (actual, expected) in lambdas.iter().zip(rule.lambdas.iter()) {
        if (actual - expected).abs() > GAUSSIAN_QUADRATURE_LAMBDA_ABS_TOL {
            return Err(CoreError::Unsupported(
                "Gaussian quadrature requires the supported Gauss-Legendre lambda nodes on [0, 1]; use trapezoidal TI for arbitrary schedules".to_string(),
            ));
        }
    }

    Ok(rule)
}

fn gaussian_quadrature_endpoint_states(temperature_k: f64) -> Result<(StatePoint, StatePoint)> {
    Ok((
        StatePoint::new(vec![0.0], temperature_k)?,
        StatePoint::new(vec![1.0], temperature_k)?,
    ))
}

#[derive(Debug)]
struct QuadratureRule {
    lambdas: &'static [f64],
    weights: &'static [f64],
}

fn gaussian_quadrature_rule_by_count(count: usize) -> Option<&'static QuadratureRule> {
    match count {
        1 => Some(&GQ_RULE_1),
        2 => Some(&GQ_RULE_2),
        3 => Some(&GQ_RULE_3),
        4 => Some(&GQ_RULE_4),
        5 => Some(&GQ_RULE_5),
        6 => Some(&GQ_RULE_6),
        7 => Some(&GQ_RULE_7),
        8 => Some(&GQ_RULE_8),
        9 => Some(&GQ_RULE_9),
        10 => Some(&GQ_RULE_10),
        11 => Some(&GQ_RULE_11),
        12 => Some(&GQ_RULE_12),
        13 => Some(&GQ_RULE_13),
        14 => Some(&GQ_RULE_14),
        15 => Some(&GQ_RULE_15),
        16 => Some(&GQ_RULE_16),
        _ => None,
    }
}

// These are Gauss-Legendre nodes and weights transformed from [-1, 1] to [0, 1].
// The lambda check tolerates normal decimal rounding when users specify irrational nodes.
const GQ_RULE_1: QuadratureRule = QuadratureRule {
    lambdas: &[0.5],
    weights: &[1.0],
};

const GQ_RULE_2: QuadratureRule = QuadratureRule {
    lambdas: &[0.21132486540518713, 0.7886751345948129],
    weights: &[0.5, 0.5],
};

const GQ_RULE_3: QuadratureRule = QuadratureRule {
    lambdas: &[0.1127016653792583, 0.5, 0.8872983346207417],
    weights: &[0.27777777777777785, 0.4444444444444444, 0.27777777777777785],
};

const GQ_RULE_4: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.06943184420297371,
        0.33000947820757187,
        0.6699905217924281,
        0.9305681557970262,
    ],
    weights: &[
        0.17392742256872684,
        0.3260725774312731,
        0.3260725774312731,
        0.17392742256872684,
    ],
};

const GQ_RULE_5: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.04691007703066802,
        0.23076534494715845,
        0.5,
        0.7692346550528415,
        0.9530899229693319,
    ],
    weights: &[
        0.11846344252809471,
        0.2393143352496831,
        0.2844444444444445,
        0.2393143352496831,
        0.11846344252809471,
    ],
};

const GQ_RULE_6: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.033765242898423975,
        0.16939530676686776,
        0.3806904069584015,
        0.6193095930415985,
        0.8306046932331322,
        0.966234757101576,
    ],
    weights: &[
        0.08566224618958487,
        0.18038078652406947,
        0.2339569672863457,
        0.2339569672863457,
        0.18038078652406947,
        0.08566224618958487,
    ],
};

const GQ_RULE_7: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.025446043828620757,
        0.12923440720030277,
        0.2970774243113014,
        0.5,
        0.7029225756886985,
        0.8707655927996972,
        0.9745539561713792,
    ],
    weights: &[
        0.06474248308443532,
        0.1398526957446383,
        0.19091502525255916,
        0.20897959183673448,
        0.19091502525255916,
        0.1398526957446383,
        0.06474248308443532,
    ],
};

const GQ_RULE_8: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.019855071751231912,
        0.10166676129318664,
        0.2372337950418355,
        0.4082826787521751,
        0.5917173212478248,
        0.7627662049581645,
        0.8983332387068134,
        0.9801449282487681,
    ],
    weights: &[
        0.050614268145188344,
        0.11119051722668717,
        0.15685332293894352,
        0.18134189168918088,
        0.18134189168918088,
        0.15685332293894352,
        0.11119051722668717,
        0.050614268145188344,
    ],
};

const GQ_RULE_9: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.015919880246186957,
        0.08198444633668212,
        0.19331428364970482,
        0.33787328829809554,
        0.5,
        0.6621267117019045,
        0.8066857163502952,
        0.9180155536633179,
        0.984080119753813,
    ],
    weights: &[
        0.04063719418078736,
        0.09032408034742856,
        0.13030534820146783,
        0.1561735385200014,
        0.16511967750062984,
        0.1561735385200014,
        0.13030534820146783,
        0.09032408034742856,
        0.04063719418078736,
    ],
};

const GQ_RULE_10: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.013046735741414128,
        0.06746831665550773,
        0.16029521585048778,
        0.2833023029353764,
        0.4255628305091844,
        0.5744371694908156,
        0.7166976970646236,
        0.8397047841495122,
        0.9325316833444923,
        0.9869532642585859,
    ],
    weights: &[
        0.033335672154344034,
        0.07472567457529018,
        0.109_543_181_257_991,
        0.13463335965499826,
        0.147_762_112_357_376_5,
        0.147_762_112_357_376_5,
        0.13463335965499826,
        0.109_543_181_257_991,
        0.07472567457529018,
        0.033335672154344034,
    ],
};

const GQ_RULE_11: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.010885670926971514,
        0.05646870011595234,
        0.13492399721297532,
        0.2404519353965941,
        0.36522842202382755,
        0.5,
        0.6347715779761725,
        0.7595480646034058,
        0.8650760027870247,
        0.9435312998840477,
        0.9891143290730284,
    ],
    weights: &[
        0.027834283558086582,
        0.06279018473245235,
        0.09314510546386721,
        0.11659688229599534,
        0.13140227225512338,
        0.13646254338895045,
        0.13140227225512338,
        0.11659688229599534,
        0.09314510546386721,
        0.06279018473245235,
        0.027834283558086582,
    ],
};

const GQ_RULE_12: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.009219682876640378,
        0.0479413718147626,
        0.11504866290284765,
        0.20634102285669126,
        0.31608425050090994,
        0.43738329574426554,
        0.5626167042557344,
        0.68391574949909,
        0.7936589771433087,
        0.8849513370971524,
        0.9520586281852375,
        0.9907803171233596,
    ],
    weights: &[
        0.02358766819325601,
        0.05346966299765944,
        0.08003916427167306,
        0.10158371336153282,
        0.11674626826917732,
        0.12457352290670134,
        0.12457352290670134,
        0.11674626826917732,
        0.10158371336153282,
        0.08003916427167306,
        0.05346966299765944,
        0.02358766819325601,
    ],
};

const GQ_RULE_13: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.007908472640705933,
        0.04120080038851104,
        0.09921095463334506,
        0.1788253302798299,
        0.2757536244817766,
        0.3847708420224326,
        0.5,
        0.6152291579775674,
        0.7242463755182233,
        0.8211746697201701,
        0.9007890453666549,
        0.958799199611489,
        0.9920915273592941,
    ],
    weights: &[
        0.02024200238265794,
        0.0460607499188643,
        0.06943675510989368,
        0.08907299038097276,
        0.10390802376844428,
        0.11314159013144857,
        0.11627577661543695,
        0.11314159013144857,
        0.10390802376844428,
        0.08907299038097276,
        0.06943675510989368,
        0.0460607499188643,
        0.02024200238265794,
    ],
};

const GQ_RULE_14: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.006858095651593843,
        0.03578255816821324,
        0.08639934246511749,
        0.15635354759415726,
        0.24237568182092295,
        0.340_443_815_536_055_1,
        0.44597252564632817,
        0.5540274743536718,
        0.6595561844639448,
        0.757624318179077,
        0.8436464524058427,
        0.9136006575348825,
        0.9642174418317868,
        0.9931419043484062,
    ],
    weights: &[
        0.017559730165876187,
        0.04007904357988015,
        0.06075928534395148,
        0.0786015835790967,
        0.09276919873896881,
        0.10259923186064777,
        0.10763192673157883,
        0.10763192673157883,
        0.10259923186064777,
        0.09276919873896881,
        0.0786015835790967,
        0.06075928534395148,
        0.04007904357988015,
        0.017559730165876187,
    ],
};

const GQ_RULE_15: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.006003740989757311,
        0.031363303799647024,
        0.0758967082947864,
        0.13779113431991497,
        0.21451391369573058,
        0.3029243264612183,
        0.39940295300128276,
        0.5,
        0.6005970469987173,
        0.6970756735387818,
        0.7854860863042694,
        0.862208865680085,
        0.9241032917052137,
        0.968636696200353,
        0.9939962590102427,
    ],
    weights: &[
        0.015376620998059323,
        0.035183023744054034,
        0.053579610233585886,
        0.06978533896307695,
        0.08313460290849689,
        0.09308050000778094,
        0.09921574266355562,
        0.10128912096278045,
        0.09921574266355562,
        0.09308050000778094,
        0.08313460290849689,
        0.06978533896307695,
        0.053579610233585886,
        0.035183023744054034,
        0.015376620998059323,
    ],
};

const GQ_RULE_16: QuadratureRule = QuadratureRule {
    lambdas: &[
        0.005299532504175031,
        0.0277124884633837,
        0.06718439880608412,
        0.1222977958224985,
        0.1910618777986781,
        0.2709916111713863,
        0.35919822461037054,
        0.4524937450811813,
        0.5475062549188188,
        0.6408017753896295,
        0.7290083888286136,
        0.8089381222013219,
        0.8777022041775016,
        0.9328156011939159,
        0.9722875115366163,
        0.994700467495825,
    ],
    weights: &[
        0.013576229705877019,
        0.031126761969323853,
        0.047579255841246296,
        0.062314485627767015,
        0.07479799440828838,
        0.08457825969750131,
        0.0913017075224618,
        0.0947253052275343,
        0.0947253052275343,
        0.0913017075224618,
        0.08457825969750131,
        0.07479799440828838,
        0.062314485627767015,
        0.047579255841246296,
        0.031126761969323853,
        0.013576229705877019,
    ],
};
