use alchemrs::parse::gromacs::{extract_dhdl as extract_gromacs_dhdl, extract_u_nk as extract_gromacs_u_nk};
use alchemrs::{extract_dhdl, extract_u_nk, extract_u_nk_with_potential};

#[test]
fn gromacs_extract_dhdl_from_real_file() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}/fixtures/gromacs/simple/dhdl.xvg");
    let series = extract_gromacs_dhdl(path, 300.0).expect("parse GROMACS dhdl.xvg");
    assert_eq!(series.time_ps(), &[0.0, 2.0, 4.0]);
    assert_eq!(series.values(), &[1.0, 2.0, 3.0]);
    assert_eq!(series.state().lambdas(), &[0.1]);
}

#[test]
fn gromacs_extract_u_nk_from_real_file() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}/fixtures/gromacs/simple/dhdl.xvg");
    let u_nk = extract_gromacs_u_nk(path, 300.0).expect("parse GROMACS dhdl.xvg");
    assert_eq!(u_nk.n_states(), 3);
    assert_eq!(u_nk.n_samples(), 3);
    assert_eq!(u_nk.time_ps(), &[0.0, 2.0, 4.0]);
    assert_eq!(u_nk.evaluated_states()[0].lambdas(), &[0.0]);
    assert_eq!(u_nk.evaluated_states()[1].lambdas(), &[0.1]);
    assert_eq!(u_nk.evaluated_states()[2].lambdas(), &[0.2]);
    assert_eq!(u_nk.data()[4], 0.0);
    assert_eq!(u_nk.data()[5], f64::INFINITY);
}

#[test]
fn top_level_dispatch_supports_gromacs_fixture() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}/fixtures/gromacs/simple/dhdl.xvg");

    let dhdl = extract_dhdl(&path, 300.0).expect("dispatch dhdl parser");
    let u_nk = extract_u_nk(&path, 300.0).expect("dispatch u_nk parser");
    let (_u_nk_with_potential, potential) =
        extract_u_nk_with_potential(&path, 300.0).expect("dispatch u_nk+potential parser");

    assert_eq!(dhdl.values().len(), 3);
    assert_eq!(u_nk.n_samples(), 3);
    assert_eq!(potential, vec![-10.0, -11.0, -12.0]);
}
