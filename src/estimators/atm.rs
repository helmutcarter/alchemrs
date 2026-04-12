use crate::data::{AtmLogQMatrix, DeltaFMatrix, StatePoint};
use crate::error::Result;

use super::uwham::{fit_log_q_input, UwhamLogQInput, UwhamOptions};

#[derive(Debug, Clone, Default)]
pub struct AtmOptions {
    pub uwham: UwhamOptions,
}

#[derive(Debug, Clone, Default)]
pub struct AtmEstimator {
    pub options: AtmOptions,
}

#[derive(Debug)]
pub struct AtmFit {
    inner: super::uwham::UwhamFit,
}

impl AtmEstimator {
    pub fn new(options: AtmOptions) -> Self {
        Self { options }
    }

    pub fn fit(&self, data: &AtmLogQMatrix) -> Result<AtmFit> {
        let input = UwhamLogQInput::from_atm(data);
        let inner = fit_log_q_input(input, &self.options.uwham)?;
        Ok(AtmFit { inner })
    }

    pub fn estimate(&self, data: &AtmLogQMatrix) -> Result<DeltaFMatrix> {
        self.fit(data)?.result()
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
}
