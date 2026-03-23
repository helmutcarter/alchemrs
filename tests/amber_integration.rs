use alchemrs::parse::amber::{extract_dhdl, extract_u_nk};

#[test]
fn amber_extract_dhdl_from_real_file() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}/fixtures/amber/acetamide_tiny/0.1/acetamide.prod.out");
    let series = extract_dhdl(path, 300.0).expect("parse AMBER output");
    assert!(!series.values().is_empty());
    assert_eq!(series.values().len(), series.time_ps().len());
}

#[test]
fn amber_extract_u_nk_from_real_file() {
    let base = env!("CARGO_MANIFEST_DIR");
    let path = format!("{base}/fixtures/amber/acetamide_tiny/0.1/acetamide.prod.out");
    let u_nk = extract_u_nk(path, 300.0).expect("parse AMBER output");
    assert!(u_nk.n_states() > 0);
    assert!(u_nk.n_samples() > 0);
}
