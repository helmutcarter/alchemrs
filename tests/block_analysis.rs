use alchemrs::{
    BarEstimator, CoreError, DhdlSeries, IexpEstimator, MbarEstimator, StatePoint, TiEstimator,
    UNkMatrix,
};
fn build_window(sampled: StatePoint, evaluated: &[StatePoint], n_samples: usize) -> UNkMatrix {
    let mut data = Vec::with_capacity(n_samples * evaluated.len());
    let mut time_ps = Vec::with_capacity(n_samples);
    for sample_idx in 0..n_samples {
        time_ps.push(sample_idx as f64);
        for state_idx in 0..evaluated.len() {
            data.push((sample_idx + state_idx) as f64 * 0.1);
        }
    }
    UNkMatrix::new(
        n_samples,
        evaluated.len(),
        data,
        time_ps,
        Some(sampled),
        evaluated.to_vec(),
    )
    .expect("window")
}

fn build_known_bar_windows() -> [UNkMatrix; 2] {
    let s0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];
    let time = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];

    let w0 = UNkMatrix::new(
        6,
        2,
        vec![0.0, 0.2, 0.1, 0.3, 0.2, 0.4, 0.0, 0.2, 0.1, 0.3, 0.2, 0.4],
        time.clone(),
        Some(s0),
        evaluated.clone(),
    )
    .unwrap();
    let w1 = UNkMatrix::new(
        6,
        2,
        vec![0.1, 0.0, 0.2, 0.1, 0.3, 0.2, 0.1, 0.0, 0.2, 0.1, 0.3, 0.2],
        time,
        Some(s1),
        evaluated,
    )
    .unwrap();

    [w0, w1]
}

fn trim_window_to_states(window: &UNkMatrix, state_indices: &[usize]) -> UNkMatrix {
    let mut data = Vec::with_capacity(window.n_samples() * state_indices.len());
    for sample_idx in 0..window.n_samples() {
        let row_offset = sample_idx * window.n_states();
        for &idx in state_indices {
            data.push(window.data()[row_offset + idx]);
        }
    }
    UNkMatrix::new_with_labels(
        window.n_samples(),
        state_indices.len(),
        data,
        window.time_ps().to_vec(),
        window.sampled_state().cloned(),
        state_indices
            .iter()
            .map(|&idx| window.evaluated_states()[idx].clone())
            .collect(),
        window.lambda_labels().map(|labels| labels.to_vec()),
    )
    .expect("trimmed window")
}

#[test]
fn ti_block_average_returns_requested_blocks() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let series = vec![
        DhdlSeries::new(
            s0.clone(),
            vec![0.0, 1.0, 2.0, 3.0],
            vec![1.0, 1.2, 0.8, 1.1],
        )
        .unwrap(),
        DhdlSeries::new(
            s1.clone(),
            vec![0.0, 1.0, 2.0, 3.0],
            vec![2.0, 2.2, 1.8, 2.1],
        )
        .unwrap(),
    ];

    let blocks = TiEstimator::default()
        .block_average(&series, 2)
        .expect("ti blocks");
    assert_eq!(blocks.len(), 2);
    assert_eq!(blocks[0].block_index(), 0);
    assert_eq!(blocks[1].block_index(), 1);
    assert_eq!(blocks[0].n_blocks(), 2);
    assert_eq!(blocks[0].from_state().lambdas(), s0.lambdas());
    assert_eq!(blocks[0].to_state().lambdas(), s1.lambdas());
}

#[test]
fn ti_block_average_requires_enough_samples_per_block() {
    let state = StatePoint::new(vec![0.0], 300.0).unwrap();
    let series = vec![
        DhdlSeries::new(state.clone(), vec![0.0, 1.0], vec![1.0, 1.1]).unwrap(),
        DhdlSeries::new(
            StatePoint::new(vec![1.0], 300.0).unwrap(),
            vec![0.0, 1.0],
            vec![2.0, 2.1],
        )
        .unwrap(),
    ];

    let err = TiEstimator::default()
        .block_average(&series, 3)
        .unwrap_err();
    assert!(matches!(
        err,
        CoreError::InvalidShape {
            expected: 3,
            found: 2
        }
    ));
}

#[test]
fn mbar_block_average_preserves_labels() {
    let s0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];
    let w0 = UNkMatrix::new_with_labels(
        4,
        2,
        vec![0.0, 0.2, 0.1, 0.3, 0.0, 0.2, 0.1, 0.3],
        vec![0.0, 1.0, 2.0, 3.0],
        Some(s0.clone()),
        evaluated.clone(),
        Some(vec!["coul-lambda".to_string(), "vdw-lambda".to_string()]),
    )
    .unwrap();
    let w1 = UNkMatrix::new_with_labels(
        4,
        2,
        vec![0.2, 0.0, 0.3, 0.1, 0.2, 0.0, 0.3, 0.1],
        vec![0.0, 1.0, 2.0, 3.0],
        Some(s1.clone()),
        evaluated,
        Some(vec!["coul-lambda".to_string(), "vdw-lambda".to_string()]),
    )
    .unwrap();

    let blocks = MbarEstimator::default()
        .block_average(&[w0, w1], 2)
        .expect("mbar blocks");
    assert_eq!(blocks.len(), 2);
    assert_eq!(
        blocks[0].lambda_labels().unwrap(),
        &["coul-lambda", "vdw-lambda"]
    );
    assert_eq!(blocks[0].from_state().lambdas(), s0.lambdas());
    assert_eq!(blocks[0].to_state().lambdas(), s1.lambdas());
}

#[test]
fn bar_block_average_returns_requested_blocks() {
    let [w0, w1] = build_known_bar_windows();
    let s0 = w0.sampled_state().unwrap().clone();
    let s1 = w1.sampled_state().unwrap().clone();

    let blocks = BarEstimator::default()
        .block_average(&[w0, w1], 2)
        .expect("bar blocks");
    assert_eq!(blocks.len(), 2);
    assert_eq!(blocks[0].block_index(), 0);
    assert_eq!(blocks[1].block_index(), 1);
    assert_eq!(blocks[0].from_state().lambdas(), s0.lambdas());
    assert_eq!(blocks[0].to_state().lambdas(), s1.lambdas());
}

#[test]
fn exp_block_average_returns_requested_blocks() {
    let s0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];
    let w0 = build_window(s0.clone(), &evaluated, 4);
    let w1 = build_window(s1.clone(), &evaluated, 4);

    let blocks = IexpEstimator::default()
        .block_average(&[w0, w1], 2)
        .expect("exp blocks");
    assert_eq!(blocks.len(), 2);
    assert_eq!(blocks[0].block_index(), 0);
    assert_eq!(blocks[1].block_index(), 1);
    assert_eq!(blocks[0].from_state().lambdas(), s0.lambdas());
    assert_eq!(blocks[0].to_state().lambdas(), s1.lambdas());
}

#[test]
fn dexp_block_average_uses_reverse_endpoints() {
    let s0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];
    let w0 = build_window(s0.clone(), &evaluated, 4);
    let w1 = build_window(s1.clone(), &evaluated, 4);

    let blocks = IexpEstimator::default()
        .reverse_block_average(&[w0, w1], 2)
        .expect("dexp blocks");
    assert_eq!(blocks.len(), 2);
    assert_eq!(blocks[0].from_state().lambdas(), s1.lambdas());
    assert_eq!(blocks[0].to_state().lambdas(), s0.lambdas());
}

#[test]
fn mbar_block_average_requires_nonzero_blocks() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone()];
    let window = build_window(s0, &evaluated, 2);

    let err = MbarEstimator::default()
        .block_average(&[window], 0)
        .unwrap_err();
    assert!(matches!(
        err,
        CoreError::InvalidShape {
            expected: 1,
            found: 0
        }
    ));
}

#[test]
fn real_gromacs_windows_support_mbar_block_average_after_grid_intersection() {
    let base = env!("CARGO_MANIFEST_DIR");
    let w2 = alchemrs::extract_u_nk(
        format!("{base}/fixtures/gromacs/1k_bar_samples/lambda-2/dhdl.xvg"),
        298.0,
    )
    .expect("parse lambda-2");
    let w3 = alchemrs::extract_u_nk(
        format!("{base}/fixtures/gromacs/1k_bar_samples/lambda-3/dhdl.xvg"),
        298.0,
    )
    .expect("parse lambda-3");

    let w2 = trim_window_to_states(&w2, &[1, 2]);
    let w3 = trim_window_to_states(&w3, &[0, 1]);

    assert_eq!(
        w2.evaluated_states()[0].lambdas(),
        &[0.0, 0.0, 0.1, 0.0, 0.0]
    );
    assert_eq!(
        w2.evaluated_states()[1].lambdas(),
        &[0.0, 0.0, 0.15, 0.0, 0.0]
    );
    assert_eq!(w3.evaluated_states(), w2.evaluated_states());

    let blocks = MbarEstimator::default()
        .block_average(&[w2, w3], 4)
        .expect("real gromacs mbar blocks");
    assert_eq!(blocks.len(), 4);
    assert_eq!(
        blocks[0].lambda_labels().unwrap(),
        &[
            "mass-lambda",
            "coul-lambda",
            "vdw-lambda",
            "bonded-lambda",
            "restraint-lambda",
        ]
    );
    assert_eq!(blocks[0].from_state().lambdas(), &[0.0, 0.0, 0.1, 0.0, 0.0]);
    assert_eq!(blocks[0].to_state().lambdas(), &[0.0, 0.0, 0.15, 0.0, 0.0]);
}
