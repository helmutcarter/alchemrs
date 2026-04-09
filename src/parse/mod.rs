use crate::data::{DhdlSeries, SwitchingTrajectory, UNkMatrix};
use crate::error::{CoreError, Result as CoreResult};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub mod amber;
pub mod gromacs;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ParseFormat {
    Amber,
    Gromacs,
}

pub fn extract_dhdl(path: impl AsRef<Path>, temperature_k: f64) -> CoreResult<DhdlSeries> {
    match detect_format(path.as_ref())? {
        ParseFormat::Amber => amber::extract_dhdl(path, temperature_k).map_err(Into::into),
        ParseFormat::Gromacs => gromacs::extract_dhdl(path, temperature_k).map_err(Into::into),
    }
}

pub fn infer_temperature(path: impl AsRef<Path>) -> CoreResult<f64> {
    let path = path.as_ref();
    match detect_format(path)? {
        ParseFormat::Amber => amber::extract_temperature(path).map_err(Into::into),
        ParseFormat::Gromacs => gromacs::extract_temperature(path).map_err(Into::into),
    }
}

pub fn extract_u_nk(path: impl AsRef<Path>, temperature_k: f64) -> CoreResult<UNkMatrix> {
    match detect_format(path.as_ref())? {
        ParseFormat::Amber => amber::extract_u_nk(path, temperature_k).map_err(Into::into),
        ParseFormat::Gromacs => gromacs::extract_u_nk(path, temperature_k).map_err(Into::into),
    }
}

pub fn extract_u_nk_with_potential(
    path: impl AsRef<Path>,
    temperature_k: f64,
) -> CoreResult<(UNkMatrix, Vec<f64>)> {
    match detect_format(path.as_ref())? {
        ParseFormat::Amber => {
            amber::extract_u_nk_with_potential(path, temperature_k).map_err(Into::into)
        }
        ParseFormat::Gromacs => {
            gromacs::extract_u_nk_with_potential(path, temperature_k).map_err(Into::into)
        }
    }
}

pub fn extract_nes_trajectory(
    path: impl AsRef<Path>,
    temperature_k: f64,
) -> CoreResult<SwitchingTrajectory> {
    match detect_format(path.as_ref())? {
        ParseFormat::Amber => {
            amber::extract_nes_trajectory(path, temperature_k).map_err(Into::into)
        }
        ParseFormat::Gromacs => Err(CoreError::Unsupported(
            "NES trajectory parsing is currently only supported for AMBER outputs".to_string(),
        )),
    }
}

fn detect_format(path: &Path) -> CoreResult<ParseFormat> {
    let file = File::open(path).map_err(|err| {
        CoreError::Parse(format!("failed to open input for format detection: {err}"))
    })?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let is_xvg = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("xvg"));

    for _ in 0..128 {
        line.clear();
        let bytes = reader.read_line(&mut line).map_err(|err| {
            CoreError::Parse(format!("failed to read input for format detection: {err}"))
        })?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.starts_with('@') || trimmed.starts_with('#') {
            return Ok(ParseFormat::Gromacs);
        }
        if trimmed.contains("CONTROL  DATA  FOR  THE  RUN")
            || trimmed.contains("MBAR Energy analysis:")
            || trimmed.contains("DV/DL")
            || trimmed.contains("clambda")
        {
            return Ok(ParseFormat::Amber);
        }
    }

    if is_xvg {
        Ok(ParseFormat::Gromacs)
    } else {
        Ok(ParseFormat::Amber)
    }
}

#[cfg(test)]
mod tests {
    use super::amber::extract_nes_trajectory;
    use super::amber::extract_u_nk;
    use super::amber::extract_u_nk_with_potential;
    use super::amber::{extract_dhdl, extract_dhdl_with_options, AmberDhdlOptions};
    use std::io::Write;

    #[test]
    fn parse_simple_amber_dhdl() {
        let content = r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Nature and format of output:
 ntpr =       10
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.5000
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
Summary of dvdl values over  2 steps:
    1.0000
    2.0000
End of dvdl summary
   5.  TIMINGS
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        let series = extract_dhdl(file.path(), 300.0).unwrap();
        let beta = 1.0 / (0.00198720425864083 * 300.0);
        let expected = vec![1.0 * beta, 2.0 * beta];
        assert_eq!(series.values(), expected.as_slice());
        assert_eq!(series.time_ps(), &[0.02, 0.04]);
        assert_eq!(series.state().lambdas(), &[0.5]);
    }

    #[test]
    fn parse_dense_amber_dhdl_summary() {
        let content = r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Nature and format of output:
 ntpr =       10
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.5000
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
Summary of dvdl values over  21 steps:
    0.0000
    1.0000
    2.0000
    3.0000
    4.0000
    5.0000
    6.0000
    7.0000
    8.0000
    9.0000
    10.0000
    11.0000
    12.0000
    13.0000
    14.0000
    15.0000
    16.0000
    17.0000
    18.0000
    19.0000
    20.0000
End of dvdl summary
   5.  TIMINGS
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        let series = extract_dhdl(file.path(), 300.0).unwrap();
        let beta = 1.0 / (0.00198720425864083 * 300.0);
        let expected = vec![0.0 * beta, 10.0 * beta, 20.0 * beta];
        assert_eq!(series.values(), expected.as_slice());
        assert_eq!(series.time_ps(), &[0.0, 0.02, 0.04]);
        assert_eq!(series.state().lambdas(), &[0.5]);
    }

    #[test]
    fn parse_dense_amber_dhdl_summary_in_full_mode() {
        let content = r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Nature and format of output:
 ntpr =       10
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.5000
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
Summary of dvdl values over  21 steps:
    0.0000
    1.0000
    2.0000
    3.0000
    4.0000
    5.0000
    6.0000
    7.0000
    8.0000
    9.0000
    10.0000
    11.0000
    12.0000
    13.0000
    14.0000
    15.0000
    16.0000
    17.0000
    18.0000
    19.0000
    20.0000
End of dvdl summary
   5.  TIMINGS
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        let series = extract_dhdl_with_options(
            file.path(),
            300.0,
            AmberDhdlOptions {
                input_stride: Some(1),
            },
        )
        .unwrap();
        let beta = 1.0 / (0.00198720425864083 * 300.0);
        let expected: Vec<f64> = (0..=20).map(|value| value as f64 * beta).collect();
        let expected_time: Vec<f64> = (0..=20).map(|idx| idx as f64 * 0.002).collect();
        assert_eq!(series.values(), expected.as_slice());
        assert_eq!(series.time_ps(), expected_time.as_slice());
        assert_eq!(series.state().lambdas(), &[0.5]);
    }

    #[test]
    fn parse_simple_amber_u_nk() {
        let content = r#"
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
 ---
MBAR Energy analysis:
Energy at 0.0000 =    -11.0
Energy at 0.1000 =     -8.0
 ---
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let u_nk = extract_u_nk(file.path(), 300.0).unwrap();
        assert_eq!(u_nk.n_states(), 2);
        assert_eq!(u_nk.n_samples(), 2);
        assert_eq!(u_nk.data().len(), 4);
    }

    #[test]
    fn parse_simple_amber_nes_trajectory() {
        let content = r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.5000
 dynlmb = 0.2000
 ntave = 10
   4.  RESULTS
Summary of dvdl values over  4 steps:
    1.0000
    2.0000
    3.0000
    4.0000
End of dvdl summary
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        let trajectory = extract_nes_trajectory(file.path(), 300.0).unwrap();
        let beta = 1.0 / (0.00198720425864083 * 300.0);
        let expected_work = (1.0 + 2.0 + 3.0 + 4.0) * (0.2 / 10.0) * beta;
        assert!((trajectory.reduced_work() - expected_work).abs() < 1e-12);
        assert_eq!(trajectory.initial_state().lambdas(), &[0.5]);
        assert_eq!(trajectory.final_state().lambdas(), &[0.58]);
        assert_eq!(trajectory.lambda_path(), &[0.5, 0.52, 0.54, 0.56]);
        assert_eq!(
            trajectory.dvdl_path(),
            &[1.0 * beta, 2.0 * beta, 3.0 * beta, 4.0 * beta]
        );
    }

    #[test]
    fn parse_amber_u_nk_retains_positive_infinity_samples() {
        let content = r#"
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
Energy at 0.1000 = **********
 ---
MBAR Energy analysis:
Energy at 0.0000 =    -11.0
Energy at 0.1000 =     -8.0
 ---
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let u_nk = extract_u_nk(file.path(), 300.0).unwrap();
        assert_eq!(u_nk.n_samples(), 2);
        assert_eq!(u_nk.time_ps(), &[0.02, 0.04]);
        assert_eq!(u_nk.data()[1], f64::INFINITY);
    }

    #[test]
    fn parse_amber_u_nk_handles_final_block_at_eof() {
        let content = r#"
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
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let u_nk = extract_u_nk(file.path(), 300.0).unwrap();
        assert_eq!(u_nk.n_samples(), 1);
        assert_eq!(u_nk.data().len(), 2);
    }

    #[test]
    fn parse_amber_u_nk_with_potential_samples() {
        let content = r#"
   2.  CONTROL  DATA  FOR  THE  RUN
Nature and format of output:
 ntpr =       10
Molecular dynamics:
 dt = 0.002
temperature regulation:
 temp0 = 300.0
Free energy options:
 clambda = 0.1000
FEP MBAR options:
    ifmbar = 1, bar_intervall = 10
    mbar_states = 2
    MBAR - lambda values considered:
      2 total:  0.0000 0.1000
    Extra energies will be computed 2 times.
   3.  ATOMIC
 begin time coords = 0.0
   4.  RESULTS
MBAR Energy analysis:
Energy at 0.0000 =    -10.0
Energy at 0.1000 =     -9.0
 ---
MBAR Energy analysis:
Energy at 0.0000 =    -11.0
Energy at 0.1000 =     -8.0
 ---
 NSTEP =       10  TIME(PS) =       0.02  EPtot =   -9.0000  DV/DL =     1.0000
 ---
 NSTEP =       20  TIME(PS) =       0.04  EPtot =   -8.0000  DV/DL =     2.0000
 ---
   5.  TIMINGS
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        let (u_nk, potential) = extract_u_nk_with_potential(file.path(), 300.0).unwrap();
        assert_eq!(u_nk.n_samples(), 2);
        assert_eq!(potential, vec![-9.0, -8.0]);
    }
}
