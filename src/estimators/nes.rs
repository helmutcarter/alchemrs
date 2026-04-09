use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use crate::data::{FreeEnergyEstimate, SwitchingTrajectory};
use crate::error::{CoreError, Result};

#[derive(Debug, Clone, Default)]
pub struct NesOptions {
    pub n_bootstrap: usize,
    pub seed: u64,
}

#[derive(Debug, Clone, Default)]
pub struct NesEstimator {
    pub options: NesOptions,
}

#[derive(Debug, Clone)]
pub struct NesFit {
    delta_f: f64,
    uncertainty: Option<f64>,
    from_state: crate::data::StatePoint,
    to_state: crate::data::StatePoint,
}

impl NesEstimator {
    pub fn new(options: NesOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, trajectories: &[SwitchingTrajectory]) -> Result<NesFit> {
        if trajectories.is_empty() {
            return Err(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            });
        }

        let first = &trajectories[0];
        for trajectory in trajectories.iter().skip(1) {
            if trajectory.initial_state() != first.initial_state()
                || trajectory.final_state() != first.final_state()
            {
                return Err(CoreError::InvalidState(
                    "all switching trajectories must share the same initial and final states"
                        .to_string(),
                ));
            }
        }

        let work: Vec<f64> = trajectories.iter().map(|t| t.reduced_work()).collect();
        let delta_f = jarzynski_delta_f(&work)?;
        let uncertainty = if self.options.n_bootstrap == 0 {
            Some(analytic_uncertainty(&work)?)
        } else {
            Some(bootstrap_uncertainty(
                &work,
                self.options.n_bootstrap,
                self.options.seed,
            )?)
        };

        Ok(NesFit {
            delta_f,
            uncertainty,
            from_state: first.initial_state().clone(),
            to_state: first.final_state().clone(),
        })
    }

    pub fn estimate(&self, trajectories: &[SwitchingTrajectory]) -> Result<FreeEnergyEstimate> {
        self.fit(trajectories)?.result()
    }
}

impl NesFit {
    pub fn result(&self) -> Result<FreeEnergyEstimate> {
        FreeEnergyEstimate::new(
            self.delta_f,
            self.uncertainty,
            self.from_state.clone(),
            self.to_state.clone(),
        )
    }
}

pub(crate) fn jarzynski_delta_f(work: &[f64]) -> Result<f64> {
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
            "NES work values must be finite or +inf".to_string(),
        ));
    }
    let n = work.len() as f64;
    let max_arg = work
        .iter()
        .map(|value| -value)
        .fold(f64::NEG_INFINITY, f64::max);
    let sum = work
        .iter()
        .map(|value| (-value - max_arg).exp())
        .sum::<f64>();
    let delta_f = -(max_arg + sum.ln() - n.ln());
    if !delta_f.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "NES delta_f became non-finite; no finite Boltzmann weight remained".to_string(),
        ));
    }
    Ok(delta_f)
}

pub(crate) fn analytic_uncertainty(work: &[f64]) -> Result<f64> {
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
            "NES work values must be finite or +inf".to_string(),
        ));
    }
    let n = work.len() as f64;
    let max_arg = work
        .iter()
        .map(|value| -value)
        .fold(f64::NEG_INFINITY, f64::max);
    let x: Vec<f64> = work.iter().map(|value| (-value - max_arg).exp()).collect();
    let mean = x.iter().sum::<f64>() / n;
    if mean == 0.0 || !mean.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "NES analytic uncertainty became non-finite; no finite Boltzmann weight remained"
                .to_string(),
        ));
    }
    let variance = x
        .iter()
        .map(|value| {
            let diff = value - mean;
            diff * diff
        })
        .sum::<f64>()
        / n;
    let uncertainty = variance.sqrt() / (n.sqrt() * mean);
    if !uncertainty.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "NES analytic uncertainty became non-finite".to_string(),
        ));
    }
    Ok(uncertainty)
}

fn bootstrap_uncertainty(work: &[f64], n_bootstrap: usize, seed: u64) -> Result<f64> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut samples = Vec::with_capacity(n_bootstrap);
    let mut resample = Vec::with_capacity(work.len());
    for _ in 0..n_bootstrap {
        resample.clear();
        for _ in 0..work.len() {
            let idx = rng.gen_range(0..work.len());
            resample.push(work[idx]);
        }
        samples.push(jarzynski_delta_f(&resample)?);
    }

    let mean = samples.iter().sum::<f64>() / samples.len() as f64;
    let variance = samples
        .iter()
        .map(|value| {
            let diff = value - mean;
            diff * diff
        })
        .sum::<f64>()
        / samples.len() as f64;
    Ok(variance.sqrt())
}
