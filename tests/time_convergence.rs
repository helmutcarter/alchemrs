use alchemrs::{
    mbar_time_convergence, ti_time_convergence, DhdlSeries, StatePoint, TimeConvergenceOptions,
    UNkMatrix,
};

#[test]
fn ti_time_convergence_uses_shared_elapsed_time() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let series = vec![
        DhdlSeries::new(s0, vec![10.0, 12.0, 14.0, 16.0], vec![0.0, 0.0, 0.0, 0.0]).unwrap(),
        DhdlSeries::new(s1, vec![20.0, 21.0, 22.0, 23.0], vec![2.0, 2.0, 4.0, 4.0]).unwrap(),
    ];

    let points = ti_time_convergence(
        &series,
        None,
        Some(TimeConvergenceOptions {
            n_points: 3,
            min_samples_per_window: 2,
        }),
    )
    .unwrap();

    assert_eq!(points.len(), 3);
    assert_eq!(points[0].elapsed_time_ps(), 2.0);
    assert_eq!(points[2].elapsed_time_ps(), 3.0);
    assert!((points[0].delta_f() - 4.0 / 3.0).abs() < 1.0e-12);
    assert!((points[2].delta_f() - 1.5).abs() < 1.0e-12);
}

#[test]
fn mbar_time_convergence_returns_shared_elapsed_points() {
    let s0 = StatePoint::new(vec![0.0], 300.0).unwrap();
    let s1 = StatePoint::new(vec![1.0], 300.0).unwrap();
    let evaluated = vec![s0.clone(), s1.clone()];
    let w0 = UNkMatrix::new(
        4,
        2,
        vec![0.0, 0.2, 0.0, 0.3, 0.0, 0.4, 0.0, 0.5],
        vec![10.0, 11.0, 12.0, 13.0],
        Some(s0),
        evaluated.clone(),
    )
    .unwrap();
    let w1 = UNkMatrix::new(
        4,
        2,
        vec![0.2, 0.0, 0.3, 0.0, 0.4, 0.0, 0.5, 0.0],
        vec![20.0, 21.0, 22.0, 23.0],
        Some(s1),
        evaluated,
    )
    .unwrap();

    let points = mbar_time_convergence(
        &[w0, w1],
        None,
        Some(TimeConvergenceOptions {
            n_points: 2,
            min_samples_per_window: 2,
        }),
    )
    .unwrap();

    assert_eq!(points.len(), 2);
    assert_eq!(points[0].elapsed_time_ps(), 1.0);
    assert_eq!(points[1].elapsed_time_ps(), 3.0);
    assert!(points.iter().all(|point| point.delta_f().is_finite()));
}
