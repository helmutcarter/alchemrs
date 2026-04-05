use alchemrs::{DhdlSeries, IntegrationMethod, StatePoint, TiEstimator, TiOptions};

fn build_series(lambdas: &[f64], means: &[f64], sem2: &[f64]) -> Vec<DhdlSeries> {
    lambdas
        .iter()
        .zip(means.iter())
        .zip(sem2.iter())
        .map(|((&lambda, &mean), &target_sem2)| {
            let spread = (3.0 * target_sem2).sqrt();
            DhdlSeries::new(
                StatePoint::new(vec![lambda], 300.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![mean - spread, mean, mean + spread],
            )
            .unwrap()
        })
        .collect()
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected}, got {actual}"
    );
}

// Reference values below were generated with SciPy 1.8.1 using:
// - scipy.interpolate.PchipInterpolator(...).integrate(x[0], x[-1])
// - scipy.interpolate.Akima1DInterpolator(...).integrate(x[0], x[-1])
// and sigma from the same central-difference Jacobian propagation used here.

#[test]
fn pchip_matches_scipy_reference_case1() {
    let series = build_series(
        &[0.0, 0.2, 0.55, 1.0],
        &[0.3, 1.1, -0.4, 0.8],
        &[0.01, 0.04, 0.09, 0.16],
    );
    let result = TiEstimator::new(TiOptions {
        method: IntegrationMethod::Pchip,
        parallel: false,
    })
    .fit(&series)
    .unwrap()
    .result()
    .unwrap();

    assert_close(result.delta_f(), 0.264_883_319_805_195, 1.0e-12);
    assert_close(
        result.uncertainty().expect("pchip uncertainty"),
        0.168_599_525_547_298_4,
        1.0e-10,
    );
}

#[test]
fn akima_matches_scipy_reference_case1() {
    let series = build_series(
        &[0.0, 0.2, 0.55, 1.0],
        &[0.3, 1.1, -0.4, 0.8],
        &[0.01, 0.04, 0.09, 0.16],
    );
    let result = TiEstimator::new(TiOptions {
        method: IntegrationMethod::Akima,
        parallel: false,
    })
    .fit(&series)
    .unwrap()
    .result()
    .unwrap();

    assert_close(result.delta_f(), 0.269_138_764_880_952_6, 1.0e-12);
    assert_close(
        result.uncertainty().expect("akima uncertainty"),
        0.171_382_858_364_074_74,
        1.0e-10,
    );
}

#[test]
fn pchip_matches_scipy_reference_case2() {
    let series = build_series(
        &[0.0, 0.15, 0.5, 0.85, 1.0],
        &[0.0, 0.7, 1.4, 1.0, 0.2],
        &[0.0025, 0.01, 0.0225, 0.04, 0.0625],
    );
    let result = TiEstimator::new(TiOptions {
        method: IntegrationMethod::Pchip,
        parallel: false,
    })
    .fit(&series)
    .unwrap()
    .result()
    .unwrap();

    assert_close(result.delta_f(), 0.994_422_412_155_335_4, 1.0e-12);
    assert_close(
        result.uncertainty().expect("pchip uncertainty"),
        0.082_405_046_598_789_75,
        1.0e-10,
    );
}

#[test]
fn akima_matches_scipy_reference_case2() {
    let series = build_series(
        &[0.0, 0.15, 0.5, 0.85, 1.0],
        &[0.0, 0.7, 1.4, 1.0, 0.2],
        &[0.0025, 0.01, 0.0225, 0.04, 0.0625],
    );
    let result = TiEstimator::new(TiOptions {
        method: IntegrationMethod::Akima,
        parallel: false,
    })
    .fit(&series)
    .unwrap()
    .result()
    .unwrap();

    assert_close(result.delta_f(), 1.008_356_891_937_102_7, 1.0e-12);
    assert_close(
        result.uncertainty().expect("akima uncertainty"),
        0.083_717_122_423_865_17,
        1.0e-10,
    );
}
