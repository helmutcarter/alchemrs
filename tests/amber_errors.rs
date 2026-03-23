use std::io::Write;

use alchemrs::parse::amber::{extract_dhdl, extract_u_nk, AmberParseError};

fn write_temp_amber(content: &str) -> tempfile::NamedTempFile {
    let mut file = tempfile::NamedTempFile::new().expect("create temp file");
    file.write_all(content.as_bytes())
        .expect("write temp AMBER output");
    file
}

#[test]
fn amber_dhdl_reports_missing_temp0() {
    let file = write_temp_amber(
        r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Nature and format of output:
 ntpr =       10
Molecular dynamics:
 dt = 0.002
Free energy options:
 clambda = 0.5000
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
 NSTEP =       10  TIME(PS) =       0.02  DV/DL =     1.0000
 ---
"#,
    );

    let error = extract_dhdl(file.path(), 300.0).expect_err("temp0 should be required");
    assert_eq!(error, AmberParseError::MissingField { field: "temp0" });
}

#[test]
fn amber_u_nk_reports_invalid_header_value() {
    let file = write_temp_amber(
        r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Molecular dynamics:
 dt = nope
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.0000
FEP MBAR options:
    bar_intervall = 10
    MBAR - lambda values considered:
      total   0.0000 0.1000
    Extra
   3.  ATOMIC
 begin time coords = 0.0
"#,
    );

    let error = extract_u_nk(file.path(), 300.0).expect_err("invalid dt should fail");
    assert_eq!(
        error,
        AmberParseError::InvalidField {
            field: "dt",
            value: "nope".to_string(),
        }
    );
}

#[test]
fn amber_u_nk_reports_incomplete_mbar_block() {
    let file = write_temp_amber(
        r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.0000
FEP MBAR options:
    bar_intervall = 10
    MBAR - lambda values considered:
      total   0.0000 0.1000
    Extra
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
MBAR Energy analysis:
Energy at 0.0000 =    -10.0
 ---
"#,
    );

    let error = extract_u_nk(file.path(), 300.0).expect_err("missing state should fail");
    assert_eq!(
        error,
        AmberParseError::IncompleteMbarBlock {
            sample_idx: 0,
            missing_lambdas: vec![0.1],
        }
    );
}

#[test]
fn amber_u_nk_reports_duplicate_mbar_energy_entries() {
    let file = write_temp_amber(
        r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.0000
FEP MBAR options:
    bar_intervall = 10
    MBAR - lambda values considered:
      total   0.0000 0.1000
    Extra
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
MBAR Energy analysis:
Energy at 0.0000 =    -10.0
Energy at 0.1000 =     -9.0
Energy at 0.1000 =     -8.0
 ---
"#,
    );

    let error = extract_u_nk(file.path(), 300.0).expect_err("duplicate lambda should fail");
    assert_eq!(
        error,
        AmberParseError::DuplicateMbarEnergy {
            sample_idx: 0,
            lambda: 0.1,
        }
    );
}
