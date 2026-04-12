mod atm;
mod bar;
mod common;
mod exp;
mod mbar;
mod nes;
mod ti;
mod uwham;

pub use atm::{AtmBindingEstimate, AtmEstimator, AtmFit, AtmOptions};
pub use bar::{BarEstimator, BarFit, BarMethod, BarOptions};
pub use exp::{IexpEstimator, IexpFit, IexpOptions};
pub use mbar::{MbarEstimator, MbarFit, MbarOptions, MbarSolver};
pub use nes::{NesEstimator, NesFit, NesOptions};
pub use ti::{sample_ti_curve, IntegrationMethod, TiEstimator, TiFit, TiOptions};
pub use uwham::{UwhamEstimator, UwhamFit, UwhamOptions};

#[cfg(test)]
mod tests {
    use crate::data::{DhdlSeries, StatePoint, SwitchingTrajectory, UNkMatrix};
    use crate::error::CoreError;

    use super::bar::{bar_estimate, BarEstimator, BarMethod, BarOptions};
    use super::common::work_values;
    use super::exp::{IexpEstimator, IexpOptions};
    use super::mbar::{MbarEstimator, MbarOptions, MbarSolver};
    use super::nes::{analytic_uncertainty, jarzynski_delta_f, NesEstimator, NesOptions};
    use super::ti::{IntegrationMethod, TiEstimator, TiOptions};
    use super::uwham::{UwhamEstimator, UwhamOptions};

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
        try_make_window(
            sampled_lambda,
            sampled_temperature,
            evaluated_lambdas,
            evaluated_temperature,
        )
        .unwrap()
    }

    fn try_make_window(
        sampled_lambda: f64,
        sampled_temperature: f64,
        evaluated_lambdas: [f64; 2],
        evaluated_temperature: f64,
    ) -> crate::error::Result<UNkMatrix> {
        let sampled_state = StatePoint::new(vec![sampled_lambda], sampled_temperature).unwrap();
        let evaluated_states = evaluated_lambdas
            .into_iter()
            .map(|lambda| StatePoint::new(vec![lambda], evaluated_temperature).unwrap())
            .collect::<Vec<_>>();
        let data = vec![0.0, 0.2, 0.1, 0.3, 0.2, 0.4];
        let time = vec![0.0, 1.0, 2.0];
        UNkMatrix::new(3, 2, data, time, Some(sampled_state), evaluated_states)
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

    fn make_multidimensional_windows() -> Vec<UNkMatrix> {
        let s0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
        let states = vec![s0.clone(), s1.clone()];
        let labels = Some(vec!["coul-lambda".to_string(), "vdw-lambda".to_string()]);

        let data0 = vec![0.0, 0.2, 0.1, 0.3, 0.2, 0.4];
        let data1 = vec![0.1, 0.0, 0.2, 0.1, 0.3, 0.2];
        let time = vec![0.0, 1.0, 2.0];

        vec![
            UNkMatrix::new_with_labels(
                3,
                2,
                data0,
                time.clone(),
                Some(s0),
                states.clone(),
                labels.clone(),
            )
            .unwrap(),
            UNkMatrix::new_with_labels(3, 2, data1, time, Some(s1), states, labels).unwrap(),
        ]
    }

    fn make_multidimensional_dhdl_series() -> Vec<DhdlSeries> {
        vec![
            DhdlSeries::new(
                StatePoint::new(vec![0.0, 0.0], 300.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![1.0, 1.1, 1.2],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![1.0, 0.0], 300.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![2.0, 2.1, 2.2],
            )
            .unwrap(),
        ]
    }

    fn make_three_state_local_bar_windows() -> Vec<UNkMatrix> {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let time = vec![0.0, 1.0, 2.0];

        let w0 = UNkMatrix::new(
            3,
            2,
            vec![0.0, 0.2, 0.1, 0.3, 0.2, 0.4],
            time.clone(),
            Some(s0.clone()),
            vec![s0.clone(), s1.clone()],
        )
        .unwrap();
        let w1 = UNkMatrix::new(
            3,
            3,
            vec![0.2, 0.0, 0.3, 0.3, 0.1, 0.2, 0.4, 0.2, 0.1],
            time.clone(),
            Some(s1.clone()),
            vec![s0, s1.clone(), s2.clone()],
        )
        .unwrap();
        let w2 = UNkMatrix::new(
            3,
            2,
            vec![0.3, 0.0, 0.2, 0.1, 0.4, 0.2],
            time,
            Some(s2),
            vec![s1, StatePoint::new(vec![1.0], 300.0).unwrap()],
        )
        .unwrap();

        vec![w0, w1, w2]
    }

    fn make_switching_trajectories() -> Vec<SwitchingTrajectory> {
        let from = StatePoint::new(vec![0.0], 300.0).unwrap();
        let to = StatePoint::new(vec![1.0], 300.0).unwrap();
        vec![
            SwitchingTrajectory::new(from.clone(), to.clone(), 0.0).unwrap(),
            SwitchingTrajectory::new(from.clone(), to.clone(), 1.0).unwrap(),
            SwitchingTrajectory::new(from, to, 2.0).unwrap(),
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

    fn assert_slice_close(left: &[f64], right: &[f64], tolerance: f64) {
        assert_eq!(left.len(), right.len());
        for (idx, (lhs, rhs)) in left.iter().zip(right.iter()).enumerate() {
            assert!(
                (*lhs - *rhs).abs() < tolerance,
                "mismatch at {idx}: {lhs} vs {rhs}"
            );
        }
    }

    fn assert_option_slice_close(left: Option<&[f64]>, right: Option<&[f64]>, tolerance: f64) {
        match (left, right) {
            (None, None) => {}
            (Some(a), Some(b)) => assert_slice_close(a, b, tolerance),
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
        let result = estimator.fit(&[d0, d1]).unwrap().result().unwrap();
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
        let result = estimator.fit(&[d0, d1, d2]).unwrap().result().unwrap();
        assert!((result.delta_f() - 1.0).abs() < 1e-12);
        assert_eq!(result.uncertainty(), Some(0.0));
    }

    #[test]
    fn ti_simpson_propagates_uncertainty_from_rule_weights() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0, 2.0], vec![0.0, 0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0, 2.0], vec![1.0, 2.0, 0.0]).unwrap();
        let d2 = DhdlSeries::new(s2, vec![0.0, 1.0, 2.0], vec![2.0, 2.0, 2.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Simpson,
            parallel: false,
        });
        let result = estimator.fit(&[d0, d1, d2]).unwrap().result().unwrap();
        let expected_sigma = (4.0_f64 / 27.0_f64).sqrt();
        assert!((result.uncertainty().unwrap() - expected_sigma).abs() < 1e-12);
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
    fn ti_cubic_spline_integrates_nonuniform_schedule() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d2 = DhdlSeries::new(s2, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::CubicSpline,
            parallel: false,
        });
        let result = estimator.fit(&[d0, d1, d2]).unwrap().result().unwrap();
        assert!((result.delta_f() - 0.625).abs() < 1e-12);
        assert_eq!(result.uncertainty(), Some(0.0));
    }

    #[test]
    fn ti_cubic_spline_propagates_uncertainty_from_spline_weights() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0, 2.0], vec![0.0, 0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0, 2.0], vec![1.0, 2.0, 0.0]).unwrap();
        let d2 = DhdlSeries::new(s2, vec![0.0, 1.0, 2.0], vec![0.0, 0.0, 0.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::CubicSpline,
            parallel: false,
        });
        let result = estimator.fit(&[d0, d1, d2]).unwrap().result().unwrap();
        let expected_sigma = 0.625_f64 / 3.0_f64.sqrt();
        assert!((result.uncertainty().unwrap() - expected_sigma).abs() < 1e-12);
    }

    #[test]
    fn ti_cubic_spline_rejects_duplicate_lambda_points() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::CubicSpline,
            parallel: false,
        });
        let err = estimator.fit(&[d0, d1]).unwrap_err();
        assert!(matches!(err, CoreError::Unsupported(_)));
    }

    #[test]
    fn ti_pchip_integrates_piecewise_cubic_hermite_schedule() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d2 = DhdlSeries::new(s2, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Pchip,
            parallel: false,
        });
        let result = estimator.fit(&[d0, d1, d2]).unwrap().result().unwrap();
        assert!((result.delta_f() - (2.0 / 3.0)).abs() < 1e-12);
        assert_eq!(result.uncertainty(), Some(0.0));
    }

    #[test]
    fn ti_akima_integrates_piecewise_cubic_hermite_schedule() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d2 = DhdlSeries::new(s2, vec![0.0, 1.0], vec![0.0, 0.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Akima,
            parallel: false,
        });
        let result = estimator.fit(&[d0, d1, d2]).unwrap().result().unwrap();
        assert!((result.delta_f() - (2.0 / 3.0)).abs() < 1e-12);
        assert_eq!(result.uncertainty(), Some(0.0));
    }

    #[test]
    fn ti_pchip_rejects_duplicate_lambda_points() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Pchip,
            parallel: false,
        });
        let err = estimator.fit(&[d0, d1]).unwrap_err();
        assert!(matches!(err, CoreError::Unsupported(_)));
    }

    #[test]
    fn ti_akima_rejects_duplicate_lambda_points() {
        let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Akima,
            parallel: false,
        });
        let err = estimator.fit(&[d0, d1]).unwrap_err();
        assert!(matches!(err, CoreError::Unsupported(_)));
    }

    #[test]
    fn ti_gaussian_quadrature_integrates_over_full_unit_interval() {
        let l0 = 0.21132486540518713;
        let l1 = 0.7886751345948129;
        let s0 = StatePoint::new(vec![l0], 300.0).unwrap();
        let s1 = StatePoint::new(vec![l1], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![l0.powi(3), l0.powi(3)]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![l1.powi(3), l1.powi(3)]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::GaussianQuadrature,
            parallel: false,
        });
        let result = estimator.fit(&[d0, d1]).unwrap().result().unwrap();
        assert!((result.delta_f() - 0.25).abs() < 1e-12);
        assert_eq!(result.uncertainty(), Some(0.0));
        assert_eq!(result.from_state().lambdas(), &[0.0]);
        assert_eq!(result.to_state().lambdas(), &[1.0]);
    }

    #[test]
    fn ti_gaussian_quadrature_rejects_nonquadrature_schedule() {
        let s0 = StatePoint::new(vec![0.25], 300.0).unwrap();
        let s1 = StatePoint::new(vec![0.75], 300.0).unwrap();
        let d0 = DhdlSeries::new(s0, vec![0.0, 1.0], vec![1.0, 1.0]).unwrap();
        let d1 = DhdlSeries::new(s1, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::GaussianQuadrature,
            parallel: false,
        });
        let err = estimator.fit(&[d0, d1]).unwrap_err();
        assert!(matches!(err, CoreError::Unsupported(_)));
    }

    #[test]
    fn ti_gaussian_quadrature_rejects_unsupported_window_count() {
        let series = (0..17)
            .map(|idx| {
                DhdlSeries::new(
                    StatePoint::new(vec![idx as f64 / 16.0], 300.0).unwrap(),
                    vec![0.0, 1.0],
                    vec![idx as f64, idx as f64],
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        let estimator = TiEstimator::new(TiOptions {
            method: IntegrationMethod::GaussianQuadrature,
            parallel: false,
        });
        let err = estimator.fit(&series).unwrap_err();
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
        .estimate(&[d0.clone(), d1.clone()])
        .unwrap();
        let parallel = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Trapezoidal,
            parallel: true,
        })
        .estimate(&[d0, d1])
        .unwrap();

        assert!((serial.delta_f() - parallel.delta_f()).abs() < 1e-12);
        assert_eq!(serial.uncertainty(), parallel.uncertainty());
    }

    #[test]
    fn ti_cubic_spline_parallel_matches_serial() {
        let d0 = DhdlSeries::new(
            StatePoint::new(vec![0.0], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![0.0, 0.1, 0.2],
        )
        .unwrap();
        let d1 = DhdlSeries::new(
            StatePoint::new(vec![0.5], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![1.0, 1.1, 1.2],
        )
        .unwrap();
        let d2 = DhdlSeries::new(
            StatePoint::new(vec![1.0], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![0.0, 0.1, 0.2],
        )
        .unwrap();

        let serial = TiEstimator::new(TiOptions {
            method: IntegrationMethod::CubicSpline,
            parallel: false,
        })
        .estimate(&[d0.clone(), d1.clone(), d2.clone()])
        .unwrap();
        let parallel = TiEstimator::new(TiOptions {
            method: IntegrationMethod::CubicSpline,
            parallel: true,
        })
        .estimate(&[d0, d1, d2])
        .unwrap();

        assert!((serial.delta_f() - parallel.delta_f()).abs() < 1e-12);
        assert_eq!(serial.uncertainty(), parallel.uncertainty());
    }

    #[test]
    fn ti_pchip_parallel_matches_serial() {
        let d0 = DhdlSeries::new(
            StatePoint::new(vec![0.0], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![0.0, 0.1, 0.2],
        )
        .unwrap();
        let d1 = DhdlSeries::new(
            StatePoint::new(vec![0.5], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![1.0, 1.1, 1.2],
        )
        .unwrap();
        let d2 = DhdlSeries::new(
            StatePoint::new(vec![1.0], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![0.0, 0.1, 0.2],
        )
        .unwrap();

        let serial = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Pchip,
            parallel: false,
        })
        .estimate(&[d0.clone(), d1.clone(), d2.clone()])
        .unwrap();
        let parallel = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Pchip,
            parallel: true,
        })
        .estimate(&[d0, d1, d2])
        .unwrap();

        assert!((serial.delta_f() - parallel.delta_f()).abs() < 1e-12);
        assert_eq!(serial.uncertainty(), parallel.uncertainty());
    }

    #[test]
    fn ti_akima_parallel_matches_serial() {
        let d0 = DhdlSeries::new(
            StatePoint::new(vec![0.0], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![0.0, 0.1, 0.2],
        )
        .unwrap();
        let d1 = DhdlSeries::new(
            StatePoint::new(vec![0.5], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![1.0, 1.1, 1.2],
        )
        .unwrap();
        let d2 = DhdlSeries::new(
            StatePoint::new(vec![1.0], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![0.0, 0.1, 0.2],
        )
        .unwrap();

        let serial = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Akima,
            parallel: false,
        })
        .estimate(&[d0.clone(), d1.clone(), d2.clone()])
        .unwrap();
        let parallel = TiEstimator::new(TiOptions {
            method: IntegrationMethod::Akima,
            parallel: true,
        })
        .estimate(&[d0, d1, d2])
        .unwrap();

        assert!((serial.delta_f() - parallel.delta_f()).abs() < 1e-12);
        assert_eq!(serial.uncertainty(), parallel.uncertainty());
    }

    #[test]
    fn ti_gaussian_quadrature_parallel_matches_serial() {
        let l0 = 0.21132486540518713;
        let l1 = 0.7886751345948129;
        let d0 = DhdlSeries::new(
            StatePoint::new(vec![l0], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![l0, l0 + 0.1, l0 + 0.2],
        )
        .unwrap();
        let d1 = DhdlSeries::new(
            StatePoint::new(vec![l1], 300.0).unwrap(),
            vec![0.0, 1.0, 2.0],
            vec![l1, l1 + 0.1, l1 + 0.2],
        )
        .unwrap();

        let serial = TiEstimator::new(TiOptions {
            method: IntegrationMethod::GaussianQuadrature,
            parallel: false,
        })
        .estimate(&[d0.clone(), d1.clone()])
        .unwrap();
        let parallel = TiEstimator::new(TiOptions {
            method: IntegrationMethod::GaussianQuadrature,
            parallel: true,
        })
        .estimate(&[d0, d1])
        .unwrap();

        assert!((serial.delta_f() - parallel.delta_f()).abs() < 1e-12);
        assert_eq!(serial.uncertainty(), parallel.uncertainty());
        assert_eq!(serial.from_state().lambdas(), &[0.0]);
        assert_eq!(serial.to_state().lambdas(), &[1.0]);
    }

    #[test]
    fn bar_parallel_matches_serial() {
        let windows = make_two_state_windows();
        let serial = BarEstimator::new(BarOptions {
            parallel: false,
            ..BarOptions::default()
        })
        .estimate(&windows)
        .unwrap();
        let parallel = BarEstimator::new(BarOptions {
            parallel: true,
            ..BarOptions::default()
        })
        .estimate(&windows)
        .unwrap();

        assert_eq!(serial.values(), parallel.values());
        assert_vec_eq_with_nan(serial.uncertainties(), parallel.uncertainties());
    }

    #[test]
    fn exp_parallel_matches_serial() {
        let windows = make_two_state_windows();
        let serial = IexpEstimator::new(IexpOptions { parallel: false })
            .estimate_with_uncertainty(&windows)
            .unwrap();
        let parallel = IexpEstimator::new(IexpOptions { parallel: true })
            .estimate_with_uncertainty(&windows)
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
        .estimate_with_uncertainty(&windows)
        .unwrap();
        let parallel = MbarEstimator::new(MbarOptions {
            parallel: true,
            ..MbarOptions::default()
        })
        .estimate_with_uncertainty(&windows)
        .unwrap();

        assert_eq!(serial.values(), parallel.values());
        assert_vec_eq_with_nan(serial.uncertainties(), parallel.uncertainties());
    }

    #[test]
    fn mbar_parallel_fit_derived_outputs_match_serial() {
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

        assert_eq!(serial.free_energies(), parallel.free_energies());
        assert_eq!(serial.log_weights(), parallel.log_weights());

        let serial_overlap = serial.overlap_matrix().unwrap();
        let parallel_overlap = parallel.overlap_matrix().unwrap();
        assert_eq!(serial_overlap.values(), parallel_overlap.values());
    }

    #[test]
    fn mbar_lbfgs_matches_fixed_point() {
        let windows = make_two_state_windows();
        let fixed_point = MbarEstimator::new(MbarOptions {
            solver: MbarSolver::FixedPoint,
            ..MbarOptions::default()
        })
        .fit(&windows)
        .unwrap();
        let lbfgs = MbarEstimator::new(MbarOptions {
            solver: MbarSolver::Lbfgs,
            ..MbarOptions::default()
        })
        .fit(&windows)
        .unwrap();

        assert_slice_close(fixed_point.free_energies(), lbfgs.free_energies(), 1e-8);
        assert_slice_close(fixed_point.log_weights(), lbfgs.log_weights(), 1e-8);

        let fixed_point_result = fixed_point.result_with_uncertainty().unwrap();
        let lbfgs_result = lbfgs.result_with_uncertainty().unwrap();
        assert_slice_close(fixed_point_result.values(), lbfgs_result.values(), 1e-8);
        assert_option_slice_close(
            fixed_point_result.uncertainties(),
            lbfgs_result.uncertainties(),
            1e-8,
        );
    }

    #[test]
    fn mbar_lbfgs_parallel_matches_serial() {
        let windows = make_two_state_windows();
        let serial = MbarEstimator::new(MbarOptions {
            solver: MbarSolver::Lbfgs,
            parallel: false,
            ..MbarOptions::default()
        })
        .fit(&windows)
        .unwrap();
        let parallel = MbarEstimator::new(MbarOptions {
            solver: MbarSolver::Lbfgs,
            parallel: true,
            ..MbarOptions::default()
        })
        .fit(&windows)
        .unwrap();

        assert_slice_close(serial.free_energies(), parallel.free_energies(), 1e-12);
        assert_slice_close(serial.log_weights(), parallel.log_weights(), 1e-12);

        let serial_result = serial.result_with_uncertainty().unwrap();
        let parallel_result = parallel.result_with_uncertainty().unwrap();
        assert_slice_close(serial_result.values(), parallel_result.values(), 1e-12);
        assert_option_slice_close(
            serial_result.uncertainties(),
            parallel_result.uncertainties(),
            1e-12,
        );
    }

    #[test]
    fn bar_rejects_mismatched_evaluated_state_grid() {
        let windows = vec![
            make_window(0.0, 300.0, [0.0, 1.0], 300.0),
            make_window(2.0, 300.0, [0.0, 2.0], 300.0),
        ];

        let err = BarEstimator::default().fit(&windows).unwrap_err();
        assert!(
            matches!(err, CoreError::InvalidState(message) if message == "window for state 0 is missing adjacent state 1 in evaluated_states")
        );
    }

    #[test]
    fn bar_supports_local_neighbor_evaluated_grids() {
        let fit = BarEstimator::default()
            .fit(&make_three_state_local_bar_windows())
            .unwrap();
        let result = fit.result().unwrap();
        assert_eq!(fit.n_states(), 3);
        assert_eq!(fit.states()[0].lambdas(), &[0.0]);
        assert_eq!(fit.states()[1].lambdas(), &[0.5]);
        assert_eq!(fit.states()[2].lambdas(), &[1.0]);
        assert_eq!(result.n_states(), 3);
        assert!(result.values()[1].is_finite());
        assert!(result.values()[5].is_finite());
    }

    #[test]
    fn exp_supports_local_neighbor_evaluated_grids() {
        let fit = IexpEstimator::default()
            .fit(&make_three_state_local_bar_windows())
            .unwrap();
        let result = fit.result_with_uncertainty().unwrap();
        assert_eq!(fit.n_states(), 3);
        assert_eq!(fit.states()[0].lambdas(), &[0.0]);
        assert_eq!(fit.states()[1].lambdas(), &[0.5]);
        assert_eq!(fit.states()[2].lambdas(), &[1.0]);
        assert_eq!(result.n_states(), 3);
        assert!(result.values()[1].is_finite());
        assert!(result.values()[2].is_finite());
        assert!(result.values()[3].is_finite());
        assert!(result.values()[5].is_finite());
        assert!(result
            .uncertainties()
            .expect("uncertainties")
            .iter()
            .all(|value| value.is_finite()));
    }

    #[test]
    fn mbar_reports_local_neighbor_grid_error_informatively() {
        let err = MbarEstimator::default()
            .fit(&make_three_state_local_bar_windows())
            .unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidState(message)
                if message.contains("local-neighbor Delta H evaluations")
                    && message.contains("Use BAR")
                    && message.contains("MBAR requires a common full state grid")
        ));
    }

    #[test]
    fn window_construction_rejects_sampled_state_temperature_mismatch() {
        let err = try_make_window(1.0, 310.0, [0.0, 1.0], 300.0).unwrap_err();
        assert!(
            matches!(err, CoreError::InvalidState(message) if message == "sampled_state not found in evaluated_states")
        );
    }

    #[test]
    fn window_construction_rejects_sampled_state_missing_from_grid() {
        let err = try_make_window(0.5, 300.0, [0.0, 1.0], 300.0).unwrap_err();
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
            matches!(err, CoreError::InvalidState(message) if message == "multiple windows for same sampled_state")
        );
    }

    #[test]
    fn ti_rejects_multidimensional_lambda_states() {
        let err = TiEstimator::default()
            .fit(&make_multidimensional_dhdl_series())
            .unwrap_err();
        assert!(matches!(
            err,
            CoreError::RequiresOneDimensionalLambda {
                operation: "estimators"
            }
        ));
    }

    #[test]
    fn mbar_supports_multidimensional_lambda_states() {
        let fit = MbarEstimator::default()
            .fit(&make_multidimensional_windows())
            .unwrap();
        let result = fit.result_with_uncertainty().unwrap();
        assert_eq!(fit.n_states(), 2);
        assert_eq!(fit.states()[0].lambdas(), &[0.0, 0.0]);
        assert_eq!(fit.states()[1].lambdas(), &[1.0, 0.0]);
        assert_eq!(fit.lambda_labels().unwrap(), &["coul-lambda", "vdw-lambda"]);
        assert_eq!(result.n_states(), 2);
    }

    #[test]
    fn uwham_supports_multidimensional_lambda_states() {
        let fit = UwhamEstimator::default()
            .fit(&make_multidimensional_windows())
            .unwrap();
        let result = fit.result().unwrap();
        assert_eq!(fit.n_states(), 2);
        assert_eq!(fit.states()[0].lambdas(), &[0.0, 0.0]);
        assert_eq!(fit.states()[1].lambdas(), &[1.0, 0.0]);
        assert_eq!(fit.lambda_labels().unwrap(), &["coul-lambda", "vdw-lambda"]);
        assert_eq!(result.n_states(), 2);
    }

    #[test]
    fn exp_supports_multidimensional_lambda_states() {
        let fit = IexpEstimator::default()
            .fit(&make_multidimensional_windows())
            .unwrap();
        let result = fit.result_with_uncertainty().unwrap();
        assert_eq!(fit.n_states(), 2);
        assert_eq!(fit.states()[0].lambdas(), &[0.0, 0.0]);
        assert_eq!(fit.states()[1].lambdas(), &[1.0, 0.0]);
        assert_eq!(result.n_states(), 2);
        assert_eq!(result.states()[0].lambdas(), &[0.0, 0.0]);
        assert_eq!(result.states()[1].lambdas(), &[1.0, 0.0]);
        assert_eq!(
            result.lambda_labels().unwrap(),
            &["coul-lambda", "vdw-lambda"]
        );
    }

    #[test]
    fn bar_supports_multidimensional_lambda_states() {
        let fit = BarEstimator::default()
            .fit(&make_multidimensional_windows())
            .unwrap();
        let result = fit.result().unwrap();
        assert_eq!(fit.n_states(), 2);
        assert_eq!(fit.states()[0].lambdas(), &[0.0, 0.0]);
        assert_eq!(fit.states()[1].lambdas(), &[1.0, 0.0]);
        assert_eq!(result.n_states(), 2);
        assert_eq!(result.states()[0].lambdas(), &[0.0, 0.0]);
        assert_eq!(result.states()[1].lambdas(), &[1.0, 0.0]);
        assert_eq!(
            result.lambda_labels().unwrap(),
            &["coul-lambda", "vdw-lambda"]
        );
    }

    #[test]
    fn exp_supports_positive_infinity_work_values() {
        let windows = make_two_state_windows_with_positive_infinity();
        let result = IexpEstimator::default()
            .estimate_with_uncertainty(&windows)
            .unwrap();
        assert!(result.values().iter().all(|value| value.is_finite()));
    }

    #[test]
    fn uwham_supports_positive_infinity_u_nk_values() {
        let windows = make_two_state_windows_with_positive_infinity();
        let fit = UwhamEstimator::default().fit(&windows).unwrap();
        let result = fit.result().unwrap();
        assert!(result.values().iter().all(|value| value.is_finite()));
        assert!(fit.weights().iter().all(|value| value.is_finite()));
    }

    #[test]
    fn uwham_returns_antisymmetric_delta_f_matrix() {
        let windows = make_two_state_windows();
        let uwham = UwhamEstimator::default().estimate(&windows).unwrap();
        assert_eq!(uwham.values()[0], 0.0);
        assert_eq!(uwham.values()[3], 0.0);
        assert!((uwham.values()[1] + uwham.values()[2]).abs() < 1e-12);
    }

    #[test]
    fn uwham_parallel_matches_serial() {
        let windows = make_two_state_windows();
        let serial = UwhamEstimator::new(UwhamOptions {
            parallel: false,
            ..UwhamOptions::default()
        })
        .fit(&windows)
        .unwrap();
        let parallel = UwhamEstimator::new(UwhamOptions {
            parallel: true,
            ..UwhamOptions::default()
        })
        .fit(&windows)
        .unwrap();
        assert_slice_close(serial.free_energies(), parallel.free_energies(), 1e-12);
        assert_slice_close(serial.weights(), parallel.weights(), 1e-12);
        assert_slice_close(
            serial.result().unwrap().values(),
            parallel.result().unwrap().values(),
            1e-12,
        );
    }

    #[test]
    fn uwham_respects_iteration_limit() {
        let windows = make_two_state_windows();
        let err = UwhamEstimator::new(UwhamOptions {
            max_iterations: 0,
            tolerance: 1.0e-7,
            parallel: false,
        })
        .fit(&windows)
        .unwrap_err();
        assert!(matches!(err, CoreError::ConvergenceFailure));
    }

    #[test]
    fn bar_supports_positive_infinity_work_values() {
        let windows = make_two_state_windows_with_positive_infinity();
        let result = BarEstimator::default().estimate(&windows).unwrap();
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
    fn bar_returns_error_when_iteration_budget_is_exhausted() {
        let err = bar_estimate(
            &[-5.0, -5.0],
            &[-5.0, -2.0],
            BarMethod::FalsePosition,
            0,
            1.0e-7,
        )
        .unwrap_err();

        assert!(matches!(err, CoreError::ConvergenceFailure));
    }

    #[test]
    fn bar_adjacent_uncertainty_matches_solver_sigma() {
        let windows = make_two_state_windows();
        let result = BarEstimator::default().estimate(&windows).unwrap();

        let w_f = work_values(&windows[0], 0, 1).unwrap();
        let w_r = work_values(&windows[1], 0, 1)
            .unwrap()
            .into_iter()
            .map(|w| -w)
            .collect::<Vec<_>>();
        let (_, expected_sigma) = bar_estimate(
            &w_f,
            &w_r,
            BarMethod::FalsePosition,
            BarOptions::default().maximum_iterations,
            BarOptions::default().relative_tolerance,
        )
        .unwrap();

        let sigma = result.uncertainties().expect("uncertainties")[1];
        assert!((sigma - expected_sigma).abs() < 1e-12);
        assert!((result.uncertainties().expect("uncertainties")[2] - expected_sigma).abs() < 1e-12);
    }

    #[test]
    fn bar_cumulative_uncertainty_propagates_adjacent_sigmas() {
        let windows = make_three_state_local_bar_windows();
        let result = BarEstimator::default().estimate(&windows).unwrap();
        let uncertainties = result.uncertainties().expect("uncertainties");
        let expected =
            (uncertainties[1] * uncertainties[1] + uncertainties[5] * uncertainties[5]).sqrt();
        assert!((uncertainties[2] - expected).abs() < 1e-12);
        assert!((uncertainties[6] - expected).abs() < 1e-12);
    }

    #[test]
    fn mbar_supports_positive_infinity_reduced_potentials() {
        let windows = make_two_state_windows_with_positive_infinity();
        let result = MbarEstimator::new(MbarOptions::default())
            .estimate(&windows)
            .unwrap();
        assert!(result.values().iter().all(|value| value.is_finite()));
    }

    #[test]
    fn nes_matches_jarzynski_average() {
        let trajectories = make_switching_trajectories();
        let expected = jarzynski_delta_f(&[0.0, 1.0, 2.0]).unwrap();
        let result = NesEstimator::default().estimate(&trajectories).unwrap();
        assert!((result.delta_f() - expected).abs() < 1e-12);
    }

    #[test]
    fn nes_default_uncertainty_matches_analytic_delta_method() {
        let trajectories = make_switching_trajectories();
        let expected = analytic_uncertainty(&[0.0, 1.0, 2.0]).unwrap();
        let result = NesEstimator::default().estimate(&trajectories).unwrap();
        assert!((result.uncertainty().expect("uncertainty") - expected).abs() < 1e-12);
    }

    #[test]
    fn nes_bootstrap_is_deterministic_for_fixed_seed() {
        let trajectories = make_switching_trajectories();
        let estimator = NesEstimator::new(NesOptions {
            n_bootstrap: 32,
            seed: 17,
        });
        let left = estimator.estimate(&trajectories).unwrap();
        let right = estimator.estimate(&trajectories).unwrap();
        assert_eq!(left.uncertainty(), right.uncertainty());
        assert!(left.uncertainty().expect("uncertainty") >= 0.0);
    }
}
