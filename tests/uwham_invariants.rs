use alchemrs::{CoreError, StatePoint, UNkMatrix, UwhamEstimator};

fn make_two_state_windows() -> Vec<UNkMatrix> {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let states = vec![s0.clone(), s1.clone()];
    let time = vec![0.0, 1.0, 2.0];

    let data0 = vec![0.0, 0.2, 0.1, 0.3, 0.2, 0.4];
    let data1 = vec![0.1, 0.0, 0.2, 0.1, 0.3, 0.2];

    vec![
        UNkMatrix::new(3, 2, data0, time.clone(), Some(s0), states.clone()).unwrap(),
        UNkMatrix::new(3, 2, data1, time, Some(s1), states).unwrap(),
    ]
}

fn make_two_state_windows_with_positive_infinity() -> Vec<UNkMatrix> {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let states = vec![s0.clone(), s1.clone()];
    let time = vec![0.0, 1.0, 2.0];

    let data0 = vec![0.0, 1.0, 0.0, 2.0, 0.0, f64::INFINITY];
    let data1 = vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0];

    vec![
        UNkMatrix::new(3, 2, data0, time.clone(), Some(s0), states.clone()).unwrap(),
        UNkMatrix::new(3, 2, data1, time, Some(s1), states).unwrap(),
    ]
}

fn reverse_samples(window: &UNkMatrix) -> UNkMatrix {
    let mut data = Vec::with_capacity(window.data().len());
    for sample_idx in (0..window.n_samples()).rev() {
        let start = sample_idx * window.n_states();
        data.extend_from_slice(&window.data()[start..start + window.n_states()]);
    }

    UNkMatrix::new_with_labels(
        window.n_samples(),
        window.n_states(),
        data,
        window.time_ps().to_vec(),
        window.sampled_state().cloned(),
        window.evaluated_states().to_vec(),
        window.lambda_labels().map(|labels| labels.to_vec()),
    )
    .unwrap()
}

#[test]
fn uwham_diagonal_is_zero_and_matrix_is_antisymmetric() {
    let result = UwhamEstimator::default()
        .estimate(&make_two_state_windows())
        .unwrap();

    let n_states = result.n_states();
    for i in 0..n_states {
        assert_eq!(result.values()[i * n_states + i], 0.0);
        for j in 0..n_states {
            let lhs = result.values()[i * n_states + j];
            let rhs = result.values()[j * n_states + i];
            assert!(
                (lhs + rhs).abs() < 1e-12,
                "antisymmetry failed at ({i}, {j})"
            );
        }
    }
}

#[test]
fn uwham_is_invariant_to_window_order() {
    let windows = make_two_state_windows();
    let original = UwhamEstimator::default().estimate(&windows).unwrap();

    let mut swapped = windows;
    swapped.reverse();
    let swapped_result = UwhamEstimator::default().estimate(&swapped).unwrap();

    assert_eq!(original.states(), swapped_result.states());
    for (lhs, rhs) in original.values().iter().zip(swapped_result.values().iter()) {
        assert!((lhs - rhs).abs() < 1e-12);
    }
}

#[test]
fn uwham_is_invariant_to_sample_order_within_windows() {
    let windows = make_two_state_windows();
    let original = UwhamEstimator::default().estimate(&windows).unwrap();

    let reversed = windows.iter().map(reverse_samples).collect::<Vec<_>>();
    let reversed_result = UwhamEstimator::default().estimate(&reversed).unwrap();

    assert_eq!(original.states(), reversed_result.states());
    for (lhs, rhs) in original
        .values()
        .iter()
        .zip(reversed_result.values().iter())
    {
        assert!((lhs - rhs).abs() < 1e-12);
    }
}

#[test]
fn uwham_supports_positive_infinity_inputs() {
    let fit = UwhamEstimator::default()
        .fit(&make_two_state_windows_with_positive_infinity())
        .unwrap();

    assert!(fit.free_energies().iter().all(|value| value.is_finite()));
    assert!(fit.weights().iter().all(|value| value.is_finite()));
}

#[test]
fn uwham_rejects_empty_input() {
    let err = UwhamEstimator::default().fit(&[]).unwrap_err();
    assert!(matches!(
        err,
        CoreError::InvalidShape {
            expected: 1,
            found: 0
        }
    ));
}
