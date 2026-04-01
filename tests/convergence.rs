use alchemrs::{
    bar_convergence, dexp_convergence, mbar_convergence, ti_convergence, CoreError, DhdlSeries,
    StatePoint, UNkMatrix,
};

fn build_window(sampled: StatePoint, evaluated: &[StatePoint]) -> UNkMatrix {
    let data = vec![0.0; evaluated.len()];
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
    let data = vec![0.0; evaluated.len()];
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
fn mbar_convergence_tracks_prefix_count_and_labels() {
    let s0 = StatePoint::new(vec![0.0, 0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![0.5, 0.0], 300.0).unwrap();
    let s2 = StatePoint::new(vec![1.0, 0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone(), s2.clone()];
    let windows = vec![
        build_window_with_labels(s0.clone(), &evaluated),
        build_window_with_labels(s1, &evaluated),
        build_window_with_labels(s2.clone(), &evaluated),
    ];

    let points = mbar_convergence(&windows, None).expect("mbar convergence");
    assert_eq!(points.len(), 3);
    assert_eq!(points[0].n_windows(), 1);
    assert_eq!(points[2].n_windows(), 3);
    assert_eq!(points[2].from_state().lambdas(), s0.lambdas());
    assert_eq!(points[2].to_state().lambdas(), s2.lambdas());
    assert_eq!(
        points[2].lambda_labels().unwrap(),
        &["coul-lambda", "vdw-lambda"]
    );
}

#[test]
fn dexp_convergence_uses_reverse_endpoints() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];
    let windows = vec![
        build_window(s0.clone(), &evaluated),
        build_window(s1.clone(), &evaluated),
    ];

    let points = dexp_convergence(&windows, None).expect("dexp convergence");
    assert_eq!(points.len(), 1);
    assert_eq!(points[0].from_state().lambdas(), s1.lambdas());
    assert_eq!(points[0].to_state().lambdas(), s0.lambdas());
}

#[test]
fn bar_convergence_requires_two_windows() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let evaluated = vec![s0.clone()];
    let windows = vec![build_window(s0, &evaluated)];

    let err = bar_convergence(&windows, None).unwrap_err();
    assert!(matches!(
        err,
        CoreError::InvalidShape {
            expected: 2,
            found: 1
        }
    ));
}

#[test]
fn ti_convergence_returns_prefix_series() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![0.5], 300.0).unwrap();
    let s2 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let series = vec![
        DhdlSeries::new(s0.clone(), vec![0.0, 1.0], vec![1.0, 1.0]).unwrap(),
        DhdlSeries::new(s1, vec![0.0, 1.0], vec![2.0, 2.0]).unwrap(),
        DhdlSeries::new(s2.clone(), vec![0.0, 1.0], vec![3.0, 3.0]).unwrap(),
    ];

    let points = ti_convergence(&series, None).expect("ti convergence");
    assert_eq!(points.len(), 2);
    assert_eq!(points[0].n_windows(), 2);
    assert_eq!(points[1].n_windows(), 3);
    assert_eq!(points[1].from_state().lambdas(), s0.lambdas());
    assert_eq!(points[1].to_state().lambdas(), s2.lambdas());
    assert!(points[1].lambda_labels().is_none());
}
