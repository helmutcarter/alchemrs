use crate::analysis::{self, BlockEstimate};
use crate::data::{DhdlSeries, FreeEnergyEstimate, StatePoint};
use crate::error::{CoreError, Result};

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
