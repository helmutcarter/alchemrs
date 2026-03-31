use alchemrs::{overlap_matrix, StatePoint, UNkMatrix};

fn build_window(sampled: StatePoint, evaluated: &[StatePoint]) -> UNkMatrix {
    let data = vec![0.0, 0.0];
    let time_ps = vec![0.0];
    UNkMatrix::new(
        1,
        evaluated.len(),
        data,
        time_ps,
        Some(sampled),
        evaluated.to_vec(),
    )
    .expect("window")
}

fn build_window_with_labels(sampled: StatePoint, evaluated: &[StatePoint]) -> UNkMatrix {
    let data = vec![0.0, 0.0];
    let time_ps = vec![0.0];
    UNkMatrix::new_with_labels(
        1,
        evaluated.len(),
        data,
        time_ps,
        Some(sampled),
        evaluated.to_vec(),
        Some(vec!["coul-lambda".to_string(), "vdw-lambda".to_string()]),
    )
    .expect("window")
}

#[test]
fn overlap_matrix_rows_sum_to_one() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];

    let w0 = build_window(s0, &evaluated);
    let w1 = build_window(s1, &evaluated);

    let overlap = overlap_matrix(&[w0, w1], None).expect("overlap");
    let n = overlap.n_states();
    let values = overlap.values();

    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += values[i * n + j];
        }
        assert!((sum - 1.0).abs() < 1e-12, "row {i} sums to {sum}");
    }
}

#[test]
fn overlap_matrix_preserves_lambda_labels() {
    let s0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];

    let w0 = build_window_with_labels(s0, &evaluated);
    let w1 = build_window_with_labels(s1, &evaluated);

    let overlap = overlap_matrix(&[w0, w1], None).expect("overlap");
    assert_eq!(overlap.lambda_labels().unwrap(), &["coul-lambda", "vdw-lambda"]);
}
