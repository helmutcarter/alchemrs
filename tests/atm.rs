use alchemrs::{AtmLogQMatrix, CoreError, StatePoint};

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
