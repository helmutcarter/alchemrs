use crate::data::{AtmLogQMatrix, AtmSampleSet, DeltaFMatrix, FreeEnergyEstimate, StatePoint};
use crate::error::{CoreError, Result};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use super::uwham::{fit_log_q_input, UwhamLogQInput, UwhamOptions};

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum AtmUncertaintyMethod {
    #[default]
    Analytical,
    Bootstrap {
        n_bootstrap: usize,
        seed: u64,
    },
    None,
}

#[derive(Debug, Clone)]
pub struct AtmOptions {
    pub uwham: UwhamOptions,
    pub uncertainty: AtmUncertaintyMethod,
}

impl Default for AtmOptions {
    fn default() -> Self {
        Self {
            uwham: UwhamOptions::default(),
            uncertainty: AtmUncertaintyMethod::Analytical,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct AtmEstimator {
    pub options: AtmOptions,
}

#[derive(Debug)]
pub struct AtmFit {
    inner: super::uwham::UwhamFit,
    leg_uncertainty: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtmBindingEstimate {
    delta_f: f64,
    uncertainty: Option<f64>,
    leg1: FreeEnergyEstimate,
    leg2: FreeEnergyEstimate,
}

impl AtmEstimator {
    pub fn new(options: AtmOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, data: &AtmLogQMatrix) -> Result<AtmFit> {
        if matches!(
            self.options.uncertainty,
            AtmUncertaintyMethod::Bootstrap { .. }
        ) {
            return Err(CoreError::Unsupported(
                "ATM bootstrap uncertainty requires AtmSampleSet input".to_string(),
            ));
        }
        let input = UwhamLogQInput::from_atm(data);
        let inner = fit_log_q_input(input, &self.options.uwham)?;
        let leg_uncertainty = match self.options.uncertainty {
            AtmUncertaintyMethod::Analytical => Some(endpoint_uncertainty(&inner)?),
            AtmUncertaintyMethod::Bootstrap { .. } | AtmUncertaintyMethod::None => None,
        };
        Ok(AtmFit {
            inner,
            leg_uncertainty,
        })
    }

    pub fn fit_leg(&self, samples: &AtmSampleSet) -> Result<AtmFit> {
        let matrix = samples.to_log_q_matrix()?;
        let input = UwhamLogQInput::from_atm(&matrix);
        let inner = fit_log_q_input(input, &self.options.uwham)?;
        let leg_uncertainty = match self.options.uncertainty {
            AtmUncertaintyMethod::Analytical => Some(endpoint_uncertainty(&inner)?),
            AtmUncertaintyMethod::Bootstrap { n_bootstrap, seed } => {
                Some(bootstrap_leg_uncertainty(self, samples, n_bootstrap, seed)?)
            }
            AtmUncertaintyMethod::None => None,
        };
        Ok(AtmFit {
            inner,
            leg_uncertainty,
        })
    }

    pub fn estimate(&self, data: &AtmLogQMatrix) -> Result<DeltaFMatrix> {
        self.fit(data)?.result()
    }

    pub fn estimate_leg(&self, samples: &AtmSampleSet) -> Result<FreeEnergyEstimate> {
        self.fit_leg(samples)?.leg_result()
    }

    pub fn estimate_binding(
        &self,
        leg1: &AtmSampleSet,
        leg2: &AtmSampleSet,
    ) -> Result<AtmBindingEstimate> {
        if leg1.schedule().direction() == leg2.schedule().direction() {
            return Err(CoreError::InvalidState(
                "ATM binding analysis requires legs with opposite directions".to_string(),
            ));
        }
        if (leg1.schedule().temperature_k() - leg2.schedule().temperature_k()).abs() > 1.0e-9 {
            return Err(CoreError::InvalidState(
                "ATM binding analysis requires both legs to use the same temperature".to_string(),
            ));
        }

        let leg1_estimate = self.estimate_leg(leg1)?;
        let leg2_estimate = self.estimate_leg(leg2)?;
        let delta_f = leg1_estimate.delta_f() - leg2_estimate.delta_f();
        let uncertainty = match (leg1_estimate.uncertainty(), leg2_estimate.uncertainty()) {
            (Some(left), Some(right)) => Some((left * left + right * right).sqrt()),
            _ => None,
        };
        Ok(AtmBindingEstimate {
            delta_f,
            uncertainty,
            leg1: leg1_estimate,
            leg2: leg2_estimate,
        })
    }

    pub fn estimate_rbfe(
        &self,
        leg1: &AtmSampleSet,
        leg2: &AtmSampleSet,
    ) -> Result<AtmBindingEstimate> {
        self.estimate_binding(leg1, leg2)
    }

    pub fn estimate_abfe(
        &self,
        leg1: &AtmSampleSet,
        leg2: &AtmSampleSet,
    ) -> Result<AtmBindingEstimate> {
        self.estimate_binding(leg1, leg2)
    }
}

impl AtmFit {
    pub fn n_observations(&self) -> usize {
        self.inner.n_observations()
    }

    pub fn n_states(&self) -> usize {
        self.inner.n_states()
    }

    pub fn states(&self) -> &[StatePoint] {
        self.inner.states()
    }

    pub fn lambda_labels(&self) -> Option<&[String]> {
        self.inner.lambda_labels()
    }

    pub fn free_energies(&self) -> &[f64] {
        self.inner.free_energies()
    }

    pub fn weights(&self) -> &[f64] {
        self.inner.weights()
    }

    pub fn result(&self) -> Result<DeltaFMatrix> {
        self.inner.result()
    }

    pub fn leg_uncertainty(&self) -> Option<f64> {
        self.leg_uncertainty
    }

    pub fn leg_result(&self) -> Result<FreeEnergyEstimate> {
        let from_state = self
            .states()
            .first()
            .ok_or(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            })?
            .clone();
        let to_state = self
            .states()
            .last()
            .ok_or(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            })?
            .clone();
        let delta_f = self
            .free_energies()
            .first()
            .zip(self.free_energies().last())
            .map(|(first, last)| first - last)
            .ok_or(CoreError::InvalidShape {
                expected: 1,
                found: 0,
            })?;
        FreeEnergyEstimate::new(delta_f, self.leg_uncertainty, from_state, to_state)
    }
}

impl AtmBindingEstimate {
    pub fn delta_f(&self) -> f64 {
        self.delta_f
    }

    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty
    }

    pub fn leg1(&self) -> &FreeEnergyEstimate {
        &self.leg1
    }

    pub fn leg2(&self) -> &FreeEnergyEstimate {
        &self.leg2
    }
}

fn endpoint_uncertainty(fit: &super::uwham::UwhamFit) -> Result<f64> {
    let n_states = fit.n_states();
    let covariance = fit.covariance()?;
    let variance = covariance[0] + covariance[(n_states - 1) * n_states + (n_states - 1)]
        - 2.0 * covariance[n_states - 1];
    if !variance.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "ATM analytical uncertainty became non-finite".to_string(),
        ));
    }
    Ok(variance.max(0.0).sqrt())
}

fn bootstrap_leg_uncertainty(
    estimator: &AtmEstimator,
    samples: &AtmSampleSet,
    n_bootstrap: usize,
    seed: u64,
) -> Result<f64> {
    if n_bootstrap == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let mut grouped = std::collections::BTreeMap::<usize, Vec<_>>::new();
    for sample in samples.samples() {
        grouped
            .entry(sample.state_id())
            .or_default()
            .push(sample.clone());
    }

    let mut rng = SmallRng::seed_from_u64(seed);
    let mut estimates = Vec::with_capacity(n_bootstrap);
    for _ in 0..n_bootstrap {
        let mut resampled = Vec::with_capacity(samples.samples().len());
        for state in samples.schedule().states() {
            let Some(group) = grouped.get(&state.state_id()) else {
                continue;
            };
            for _ in 0..group.len() {
                let idx = rng.gen_range(0..group.len());
                resampled.push(group[idx].clone());
            }
        }
        let bootstrap_set = AtmSampleSet::new(samples.schedule().clone(), resampled)?;
        let input = UwhamLogQInput::from_atm(&bootstrap_set.to_log_q_matrix()?);
        let inner = fit_log_q_input(input, &estimator.options.uwham)?;
        estimates.push(leg_delta_f(&inner)?);
    }

    let mean = estimates.iter().sum::<f64>() / estimates.len() as f64;
    let variance = estimates
        .iter()
        .map(|value| {
            let diff = value - mean;
            diff * diff
        })
        .sum::<f64>()
        / estimates.len() as f64;
    if !variance.is_finite() {
        return Err(CoreError::NonFiniteValue(
            "ATM bootstrap uncertainty became non-finite".to_string(),
        ));
    }
    Ok(variance.sqrt())
}

fn leg_delta_f(fit: &super::uwham::UwhamFit) -> Result<f64> {
    fit.free_energies()
        .first()
        .zip(fit.free_energies().last())
        .map(|(first, last)| first - last)
        .ok_or(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        })
}
