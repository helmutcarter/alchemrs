use alchemrs::parse::gromacs::{extract_dhdl as extract_gromacs_dhdl, extract_u_nk as extract_gromacs_u_nk};
use alchemrs::{extract_dhdl, extract_u_nk, extract_u_nk_with_potential, CoreError};

const REAL_FIXTURE: &str = "/fixtures/gromacs/lambda15.dhdl.xvg";
const WINDOW_FIXTURES: &[(&str, f64, &[f64])] = &[
    ("/fixtures/gromacs/lambda_0_dhdl.xvg", 0.0, &[0.0, 0.0, 0.0, 0.0, 0.0]),
    ("/fixtures/gromacs/lambda_1_dhdl.xvg", 0.05, &[0.0, 0.0, 0.05, 0.0, 0.0]),
    ("/fixtures/gromacs/lambda_2_dhdl.xvg", 0.1, &[0.0, 0.0, 0.1, 0.0, 0.0]),
    ("/fixtures/gromacs/lambda_3_dhdl.xvg", 0.15, &[0.0, 0.0, 0.15, 0.0, 0.0]),
    ("/fixtures/gromacs/lambda_15_dhdl.xvg", 0.8, &[0.0, 0.0, 0.8, 0.0, 0.0]),
];

#[test]
fn gromacs_extract_u_nk_from_real_multidimensional_file() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}{REAL_FIXTURE}");
    let u_nk = extract_gromacs_u_nk(path, 298.0).expect("parse real GROMACS dhdl.xvg");

    assert_eq!(u_nk.sampled_state().unwrap().lambdas(), &[0.0, 0.0, 0.8, 0.0, 0.0]);
    assert_eq!(
        u_nk.lambda_labels().unwrap(),
        &[
            "mass-lambda",
            "coul-lambda",
            "vdw-lambda",
            "bonded-lambda",
            "restraint-lambda",
        ]
    );
    assert_eq!(u_nk.n_states(), 3);
    assert!(u_nk.n_samples() > 10_000);
    assert_eq!(u_nk.evaluated_states().first().unwrap().lambdas(), &[0.0, 0.0, 0.7, 0.0, 0.0]);
    assert_eq!(u_nk.evaluated_states()[1].lambdas(), &[0.0, 0.0, 0.8, 0.0, 0.0]);
    assert_eq!(u_nk.evaluated_states().last().unwrap().lambdas(), &[0.0, 0.0, 0.9, 0.0, 0.0]);
}

#[test]
fn top_level_dispatch_supports_real_gromacs_u_nk_fixture() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}{REAL_FIXTURE}");

    let u_nk = extract_u_nk(&path, 298.0).expect("dispatch u_nk parser");
    let (_u_nk_with_potential, potential) =
        extract_u_nk_with_potential(&path, 298.0).expect("dispatch u_nk+potential parser");

    assert_eq!(u_nk.n_states(), 3);
    assert_eq!(potential.len(), u_nk.n_samples());
    assert!((potential[0] - 25663.127).abs() < 1e-6);
}

#[test]
fn real_multidimensional_gromacs_dhdl_rejects_scalar_parser() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}{REAL_FIXTURE}");

    let err = extract_dhdl(&path, 298.0).unwrap_err();
    assert!(matches!(err, CoreError::Parse(message) if message.contains("scalar dH/dlambda parsing is unsupported")));

    let err = extract_gromacs_dhdl(&path, 298.0).unwrap_err();
    let _ = extract_gromacs_u_nk(path, 298.0).expect("direct GROMACS u_nk parser still works");
    assert!(err.to_string().contains("scalar dH/dlambda parsing is unsupported"));
}

#[test]
fn trimmed_real_gromacs_windows_parse_across_multiple_lambda_states() {
    let base = env!("CARGO_MANIFEST_DIR");

    for (fixture, vdw_lambda, sampled) in WINDOW_FIXTURES {
        let path = format!("{base}{fixture}");
        let u_nk = extract_gromacs_u_nk(path, 298.0).expect("parse real trimmed GROMACS window");
        assert_eq!(u_nk.n_samples(), 200);
        assert_eq!(u_nk.sampled_state().unwrap().lambdas(), *sampled);
        assert_eq!(u_nk.sampled_state().unwrap().lambdas()[2], *vdw_lambda);
        assert_eq!(
            u_nk.lambda_labels().unwrap(),
            &[
                "mass-lambda",
                "coul-lambda",
                "vdw-lambda",
                "bonded-lambda",
                "restraint-lambda",
            ]
        );
    }
}
