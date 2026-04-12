use alchemrs::{
    AtmDirection, AtmEstimator, AtmLogQMatrix, AtmSample, AtmSampleSet, AtmSchedule, AtmState,
    CoreError, StatePoint, UwhamOptions,
};

fn two_states() -> Vec<StatePoint> {
    vec![
        StatePoint::new(vec![0.0], 300.0).unwrap(),
        StatePoint::new(vec![1.0], 300.0).unwrap(),
    ]
}

#[test]
fn atm_matrix_rejects_zero_observations() {
    let err = AtmLogQMatrix::new(0, 2, Vec::new(), two_states(), vec![0, 0]).unwrap_err();
    assert!(matches!(
        err,
        CoreError::InvalidShape {
            expected: 1,
            found: 0
        }
    ));
}

#[test]
fn atm_matrix_rejects_count_length_mismatch() {
    let err =
        AtmLogQMatrix::new(2, 2, vec![0.0, -1.0, -2.0, -3.0], two_states(), vec![2]).unwrap_err();
    assert!(matches!(err, CoreError::InvalidShape { .. }));
}

#[test]
fn atm_matrix_accepts_valid_log_q_input() {
    let matrix = AtmLogQMatrix::new_with_labels(
        2,
        2,
        vec![0.0, -1.0, -2.0, f64::NEG_INFINITY],
        two_states(),
        vec![1, 1],
        Some(vec!["lambda".to_string()]),
    )
    .unwrap();

    assert_eq!(matrix.n_observations(), 2);
    assert_eq!(matrix.n_states(), 2);
    assert_eq!(matrix.sampled_counts(), &[1, 1]);
    assert_eq!(matrix.lambda_labels().unwrap(), &["lambda"]);
}

#[test]
fn atm_estimator_produces_antisymmetric_delta_f_matrix() {
    let matrix = AtmLogQMatrix::new_with_labels(
        4,
        2,
        vec![0.0, -0.2, -0.1, -0.3, -0.2, 0.0, -0.3, -0.1],
        two_states(),
        vec![2, 2],
        Some(vec!["lambda".to_string()]),
    )
    .unwrap();

    let fit = AtmEstimator::default().fit(&matrix).unwrap();
    let result = fit.result().unwrap();

    assert_eq!(fit.n_observations(), 4);
    assert_eq!(fit.n_states(), 2);
    assert_eq!(fit.lambda_labels().unwrap(), &["lambda"]);
    assert_eq!(result.values()[0], 0.0);
    assert_eq!(result.values()[3], 0.0);
    assert!((result.values()[1] + result.values()[2]).abs() < 1e-12);
}

#[test]
fn atm_estimator_parallel_matches_serial() {
    let matrix = AtmLogQMatrix::new(
        4,
        2,
        vec![0.0, -0.2, -0.1, -0.3, -0.2, 0.0, -0.3, -0.1],
        two_states(),
        vec![2, 2],
    )
    .unwrap();

    let serial = AtmEstimator::default().fit(&matrix).unwrap();
    let parallel = AtmEstimator::new(alchemrs::AtmOptions {
        uwham: UwhamOptions {
            parallel: true,
            ..UwhamOptions::default()
        },
    })
    .fit(&matrix)
    .unwrap();

    assert_eq!(serial.n_states(), parallel.n_states());
    for (left, right) in serial
        .free_energies()
        .iter()
        .zip(parallel.free_energies().iter())
    {
        assert!((left - right).abs() < 1e-12);
    }
    for (left, right) in serial.weights().iter().zip(parallel.weights().iter()) {
        assert!((left - right).abs() < 1e-12);
    }
}

#[test]
fn atm_sample_set_to_log_q_matrix_matches_expected_standard_softplus_limit() {
    let schedule = AtmSchedule::new(vec![
        AtmState::new(
            0,
            AtmDirection::Forward,
            0.0,
            0.0,
            None,
            0.0,
            0.0,
            None,
            0.0,
            300.0,
        )
        .unwrap(),
        AtmState::new(
            1,
            AtmDirection::Forward,
            1.0,
            1.0,
            None,
            0.0,
            0.0,
            None,
            0.0,
            300.0,
        )
        .unwrap(),
    ])
    .unwrap();
    let samples = vec![
        AtmSample::new(0, 10.0, 2.0).unwrap(),
        AtmSample::new(1, 12.0, 2.0).unwrap(),
    ];
    let matrix = AtmSampleSet::new(schedule, samples)
        .unwrap()
        .to_log_q_matrix()
        .unwrap();

    let beta = 1.0 / (0.001_986_209 * 300.0);
    let expected = [-beta * 10.0, -beta * 12.0, -beta * 10.0, -beta * 12.0];
    for (lhs, rhs) in matrix.log_q().iter().zip(expected.iter()) {
        assert!((lhs - rhs).abs() < 1.0e-12);
    }
    assert_eq!(matrix.sampled_counts(), &[1, 1]);
}

#[test]
fn atm_sample_set_reverses_state_order_for_reverse_leg() {
    let schedule = AtmSchedule::new(vec![
        AtmState::new(
            0,
            AtmDirection::Reverse,
            0.0,
            0.0,
            None,
            0.0,
            0.0,
            None,
            0.0,
            300.0,
        )
        .unwrap(),
        AtmState::new(
            1,
            AtmDirection::Reverse,
            1.0,
            1.0,
            None,
            0.0,
            0.0,
            None,
            0.0,
            300.0,
        )
        .unwrap(),
    ])
    .unwrap();
    let matrix = AtmSampleSet::new(
        schedule,
        vec![
            AtmSample::new(0, 10.0, 2.0).unwrap(),
            AtmSample::new(1, 12.0, 2.0).unwrap(),
        ],
    )
    .unwrap()
    .to_log_q_matrix()
    .unwrap();

    assert_eq!(matrix.states()[0].lambdas()[0], 1.0);
    assert_eq!(matrix.states()[1].lambdas()[0], 0.0);
}
