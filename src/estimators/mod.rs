mod bar;
mod common;
mod exp;
mod mbar;
mod ti;

pub use bar::{BarEstimator, BarMethod, BarOptions, BarUncertainty};
pub use exp::{ExpEstimator, ExpOptions};
pub use mbar::{mbar_log_weights_from_windows, MbarEstimator, MbarOptions};
pub use ti::{IntegrationMethod, TiEstimator, TiOptions};

#[cfg(test)]
mod tests {
    use crate::data::{DhdlSeries, StatePoint, UNkMatrix};
    use crate::error::CoreError;

    use super::bar::{bar_estimate, BarEstimator, BarMethod, BarOptions};
    use super::common::work_values;
    use super::exp::{ExpEstimator, ExpOptions};
    use super::mbar::{MbarEstimator, MbarOptions};
    use super::ti::{IntegrationMethod, TiEstimator, TiOptions};

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
            make_window(2.0, 300.0, [0.0, 2.0], 300.0),
        ];

        let err = BarEstimator::default().fit(&windows).unwrap_err();
        assert!(
            matches!(err, CoreError::InvalidState(message) if message == "evaluated_states differ between windows")
        );
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
        let result = MbarEstimator::default()
            .fit(&make_multidimensional_windows())
            .unwrap();
        assert_eq!(result.n_states(), 2);
        assert_eq!(result.states()[0].lambdas(), &[0.0, 0.0]);
        assert_eq!(result.states()[1].lambdas(), &[1.0, 0.0]);
        assert_eq!(
            result.lambda_labels().unwrap(),
            &["coul-lambda", "vdw-lambda"]
        );
    }

    #[test]
    fn exp_supports_multidimensional_lambda_states() {
        let result = ExpEstimator::default()
            .fit(&make_multidimensional_windows())
            .unwrap();
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
        let result = BarEstimator::default()
            .fit(&make_multidimensional_windows())
            .unwrap();
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
        let result = BarEstimator::default().fit(&windows).unwrap();

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
