use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::OnceLock;

use alchemrs_core::{CoreError, DhdlSeries, StatePoint};
use thiserror::Error;

pub mod amber {
    use super::*;

    const FP_RE: &str = r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?";
    const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;

    pub type Result<T> = std::result::Result<T, AmberParseError>;

    #[derive(Debug, Error, Clone, PartialEq)]
    pub enum AmberParseError {
        #[error("failed to {operation} AMBER output: {message}")]
        Io {
            operation: &'static str,
            message: String,
        },
        #[error("failed to locate {field} in AMBER output")]
        MissingField { field: &'static str },
        #[error("invalid {field} value `{value}` in AMBER output")]
        InvalidField { field: &'static str, value: String },
        #[error("missing value for {field} in AMBER output")]
        MissingFieldValue { field: &'static str },
        #[error(
            "temperature mismatch: file {file_temperature_k:.2} K vs input {input_temperature_k:.2} K"
        )]
        TemperatureMismatch {
            file_temperature_k: f64,
            input_temperature_k: f64,
        },
        #[error("missing {field} in NSTEP block")]
        MissingGradientField { field: &'static str },
        #[error("missing {field} in NSTEP block")]
        MissingPotentialField { field: &'static str },
        #[error("no DV/DL gradients found in AMBER output")]
        NoGradients,
        #[error("no EPtot values found in AMBER output")]
        NoPotentialSamples,
        #[error("no MBAR lambda values found in AMBER output")]
        NoMbarLambdas,
        #[error("clambda {clambda:.4} not found in MBAR lambda list; cannot build u_nk")]
        LambdaNotInGrid { clambda: f64 },
        #[error("no MBAR energy samples found in AMBER output")]
        NoMbarSamples,
        #[error("all MBAR samples contained non-finite values")]
        NoFiniteMbarSamples,
        #[error(
            "EPtot sample count {potential_samples} does not match MBAR sample count {u_nk_samples}"
        )]
        PotentialSampleCountMismatch {
            u_nk_samples: usize,
            potential_samples: usize,
        },
        #[error(
            "missing MBAR energy entries in block {sample_idx} for lambdas {missing_lambdas:?}"
        )]
        IncompleteMbarBlock {
            sample_idx: usize,
            missing_lambdas: Vec<f64>,
        },
        #[error("duplicate MBAR energy entry in block {sample_idx} for lambda {lambda:.4}")]
        DuplicateMbarEnergy { sample_idx: usize, lambda: f64 },
        #[error("invalid MBAR energy line `{line}`")]
        InvalidMbarEnergyLine { line: String },
        #[error("invalid MBAR energy value `{value}` in line `{line}`")]
        InvalidMbarEnergyValue { line: String, value: String },
        #[error("invalid parsed data: {message}")]
        InvalidParsedData { message: String },
    }

    impl From<AmberParseError> for CoreError {
        fn from(error: AmberParseError) -> Self {
            CoreError::Parse(error.to_string())
        }
    }

    impl From<CoreError> for AmberParseError {
        fn from(error: CoreError) -> Self {
            AmberParseError::InvalidParsedData {
                message: error.to_string(),
            }
        }
    }

    pub fn extract_dhdl(path: impl AsRef<Path>, temperature_k: f64) -> Result<DhdlSeries> {
        let file = File::open(path.as_ref()).map_err(|err| AmberParseError::Io {
            operation: "open",
            message: err.to_string(),
        })?;
        let mut reader = BufReader::new(file);

        let mut temp0: Option<f64> = None;
        let mut clambda: Option<f64> = None;
        let mut dt: Option<f64> = None;
        let mut ntpr: Option<f64> = None;
        let mut t0: Option<f64> = None;

        let mut gradients: Vec<f64> = Vec::new();
        let mut last_nstep: Option<i64> = None;
        let mut pending_block: Option<String> = None;
        let mut line = String::new();

        loop {
            line.clear();
            let bytes = reader
                .read_line(&mut line)
                .map_err(|err| AmberParseError::Io {
                    operation: "read",
                    message: err.to_string(),
                })?;
            if bytes == 0 {
                break;
            }
            let line = line.trim_end_matches(&['\r', '\n'][..]);

            if temp0.is_none() {
                if let Some(value) = capture_field(line, "temp0")? {
                    temp0 = Some(value);
                }
            }
            if clambda.is_none() {
                if let Some(value) = capture_field(line, "clambda")? {
                    clambda = Some(value);
                }
            }
            if dt.is_none() {
                if let Some(value) = capture_field(line, "dt")? {
                    dt = Some(value);
                }
            }
            if ntpr.is_none() {
                if let Some(value) = capture_field(line, "ntpr")? {
                    ntpr = Some(value);
                }
            }
            if t0.is_none() && line.contains("begin time") {
                if let Some(value) = capture_field(line, "coords")? {
                    t0 = Some(value);
                }
            }

            if let Some(mut block) = pending_block.take() {
                block.push(' ');
                block.push_str(line);
                if line.trim_start().starts_with("---") {
                    handle_gradient_block(&block, &mut last_nstep, &mut gradients)?;
                } else {
                    pending_block = Some(block);
                }
                continue;
            }

            if line.starts_with(" NSTEP") || line.starts_with("NSTEP") {
                let block = line.to_string();
                if line.trim_start().starts_with("---") {
                    handle_gradient_block(&block, &mut last_nstep, &mut gradients)?;
                } else {
                    pending_block = Some(block);
                }
            }
        }

        if let Some(block) = pending_block {
            handle_gradient_block(&block, &mut last_nstep, &mut gradients)?;
        }

        let temp0 = temp0.ok_or(AmberParseError::MissingField { field: "temp0" })?;
        if (temperature_k - temp0).abs() > 1e-2 {
            return Err(AmberParseError::TemperatureMismatch {
                file_temperature_k: temp0,
                input_temperature_k: temperature_k,
            });
        }
        let clambda = clambda.ok_or(AmberParseError::MissingField { field: "clambda" })?;
        let dt = dt.ok_or(AmberParseError::MissingField { field: "dt" })?;
        let ntpr = ntpr.ok_or(AmberParseError::MissingField { field: "ntpr" })?;
        let t0 = t0.ok_or(AmberParseError::MissingField { field: "t0" })?;

        if gradients.is_empty() {
            return Err(AmberParseError::NoGradients);
        }

        let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);
        for value in gradients.iter_mut() {
            *value *= beta;
        }

        let time_ps: Vec<f64> = (0..gradients.len())
            .map(|idx| t0 + ((idx + 1) as f64) * dt * ntpr)
            .collect();

        let state = StatePoint::new(vec![clambda], temperature_k)?;
        DhdlSeries::new(state, time_ps, gradients).map_err(Into::into)
    }

    pub fn extract_u_nk(
        path: impl AsRef<Path>,
        temperature_k: f64,
    ) -> Result<alchemrs_core::UNkMatrix> {
        let (u_nk, _potential) = extract_u_nk_internal(path.as_ref(), temperature_k, false)?;
        Ok(u_nk)
    }

    pub fn extract_u_nk_with_potential(
        path: impl AsRef<Path>,
        temperature_k: f64,
    ) -> Result<(alchemrs_core::UNkMatrix, Vec<f64>)> {
        let (u_nk, potential) = extract_u_nk_internal(path.as_ref(), temperature_k, true)?;
        Ok((u_nk, potential.expect("potential samples requested")))
    }

    fn extract_u_nk_internal(
        path: &Path,
        temperature_k: f64,
        include_potential: bool,
    ) -> Result<(alchemrs_core::UNkMatrix, Option<Vec<f64>>)> {
        let header = read_u_nk_header(path)?;
        if (temperature_k - header.temp0).abs() > 1e-2 {
            return Err(AmberParseError::TemperatureMismatch {
                file_temperature_k: header.temp0,
                input_temperature_k: temperature_k,
            });
        }

        if header.mbar_lambdas.is_empty() {
            return Err(AmberParseError::NoMbarLambdas);
        }

        let lambda_index = build_lambda_index(&header.mbar_lambdas);
        let reference_idx = lambda_index
            .get(&lambda_key(header.clambda))
            .copied()
            .ok_or(AmberParseError::LambdaNotInGrid {
                clambda: header.clambda,
            })?;

        let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);
        let file = File::open(path).map_err(|err| AmberParseError::Io {
            operation: "open",
            message: err.to_string(),
        })?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let mut in_mbar_block = false;
        let mut encountered_blocks = 0usize;
        let n_states = header.mbar_lambdas.len();
        let mut block_energies = vec![None; n_states];
        let mut reduced_row = Vec::with_capacity(n_states);
        let mut time_ps = Vec::new();
        let mut data = Vec::new();
        let mut potential_samples = Vec::new();
        let mut pending_nstep_block: Option<String> = None;
        let mut last_potential_nstep: Option<i64> = None;

        loop {
            line.clear();
            let bytes = reader
                .read_line(&mut line)
                .map_err(|err| AmberParseError::Io {
                    operation: "read",
                    message: err.to_string(),
                })?;
            if bytes == 0 {
                break;
            }
            let line = line.trim_end_matches(&['\r', '\n'][..]);

            if include_potential {
                if let Some(mut block) = pending_nstep_block.take() {
                    block.push(' ');
                    block.push_str(line);
                    if line.trim_start().starts_with("---") {
                        handle_potential_block(
                            &block,
                            &mut last_potential_nstep,
                            &mut potential_samples,
                        )?;
                    } else {
                        pending_nstep_block = Some(block);
                    }
                    continue;
                }
            }

            if include_potential && (line.starts_with(" NSTEP") || line.starts_with("NSTEP")) {
                pending_nstep_block = Some(line.to_string());
                continue;
            }

            if is_mbar_energy_start(line) {
                in_mbar_block = true;
                block_energies.fill(None);
                continue;
            }
            if !in_mbar_block {
                continue;
            }
            if line.trim_start().starts_with("---") {
                finalize_mbar_block(
                    &block_energies,
                    &header.mbar_lambdas,
                    reference_idx,
                    beta,
                    encountered_blocks,
                    header.t0,
                    header.dt,
                    header.bar_intervall,
                    &mut reduced_row,
                    &mut time_ps,
                    &mut data,
                )?;
                encountered_blocks += 1;
                in_mbar_block = false;
                continue;
            }
            parse_mbar_energy_line(line, encountered_blocks, &lambda_index, &mut block_energies)?;
        }

        if let Some(block) = pending_nstep_block {
            if include_potential {
                handle_potential_block(&block, &mut last_potential_nstep, &mut potential_samples)?;
            }
        }

        if in_mbar_block {
            finalize_mbar_block(
                &block_energies,
                &header.mbar_lambdas,
                reference_idx,
                beta,
                encountered_blocks,
                header.t0,
                header.dt,
                header.bar_intervall,
                &mut reduced_row,
                &mut time_ps,
                &mut data,
            )?;
            encountered_blocks += 1;
        }

        if encountered_blocks == 0 {
            return Err(AmberParseError::NoMbarSamples);
        }
        if time_ps.is_empty() {
            return Err(AmberParseError::NoFiniteMbarSamples);
        }

        let sampled_state = Some(StatePoint::new(vec![header.clambda], temperature_k)?);
        let evaluated_states = header
            .mbar_lambdas
            .iter()
            .map(|value| StatePoint::new(vec![*value], temperature_k))
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(AmberParseError::from)?;

        let u_nk = alchemrs_core::UNkMatrix::new(
            time_ps.len(),
            n_states,
            data,
            time_ps,
            sampled_state,
            evaluated_states,
        )
        .map_err(AmberParseError::from)?;

        if !include_potential {
            return Ok((u_nk, None));
        }

        if potential_samples.is_empty() {
            return Err(AmberParseError::NoPotentialSamples);
        }
        if potential_samples.len() != u_nk.n_samples() {
            return Err(AmberParseError::PotentialSampleCountMismatch {
                u_nk_samples: u_nk.n_samples(),
                potential_samples: potential_samples.len(),
            });
        }

        Ok((u_nk, Some(potential_samples)))
    }

    struct UNkHeader {
        temp0: f64,
        clambda: f64,
        dt: f64,
        bar_intervall: f64,
        t0: f64,
        mbar_lambdas: Vec<f64>,
    }

    fn capture_field(line: &str, field: &'static str) -> Result<Option<f64>> {
        let mut iter = line.split_whitespace().peekable();
        while let Some(token) = iter.next() {
            if token == field {
                if let Some(next) = iter.next() {
                    let value = if next == "=" {
                        iter.next()
                            .ok_or(AmberParseError::MissingFieldValue { field })?
                    } else {
                        next
                    };
                    return parse_numeric_field(field, value).map(Some);
                }
            } else if token.starts_with(field) {
                if let Some(eq_idx) = token.find('=') {
                    let value = &token[(eq_idx + 1)..];
                    if !value.is_empty() {
                        return parse_numeric_field(field, value).map(Some);
                    }
                }
            }
        }
        Ok(None)
    }

    fn parse_numeric_field(field: &'static str, raw_value: &str) -> Result<f64> {
        let trimmed = raw_value.trim_matches(|c: char| {
            !c.is_ascii_digit() && c != '.' && c != '-' && c != '+' && c != 'e' && c != 'E'
        });
        if trimmed.is_empty() {
            return Err(AmberParseError::InvalidField {
                field,
                value: raw_value.to_string(),
            });
        }
        trimmed
            .parse::<f64>()
            .map_err(|_| AmberParseError::InvalidField {
                field,
                value: raw_value.to_string(),
            })
    }

    fn read_u_nk_header(path: &Path) -> Result<UNkHeader> {
        let file = File::open(path).map_err(|err| AmberParseError::Io {
            operation: "open",
            message: err.to_string(),
        })?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let mut temp0: Option<f64> = None;
        let mut clambda: Option<f64> = None;
        let mut dt: Option<f64> = None;
        let mut bar_intervall: Option<f64> = None;
        let mut t0: Option<f64> = None;
        let mut in_mbar = false;
        let mut mbar_lambdas = Vec::new();

        loop {
            line.clear();
            let bytes = reader
                .read_line(&mut line)
                .map_err(|err| AmberParseError::Io {
                    operation: "read",
                    message: err.to_string(),
                })?;
            if bytes == 0 {
                break;
            }
            let line = line.trim_end_matches(&['\r', '\n'][..]);

            if temp0.is_none() {
                temp0 = capture_field(line, "temp0")?;
            }
            if clambda.is_none() {
                clambda = capture_field(line, "clambda")?;
            }
            if dt.is_none() {
                dt = capture_field(line, "dt")?;
            }
            if bar_intervall.is_none() {
                bar_intervall = capture_field(line, "bar_intervall")?;
            }
            if t0.is_none() && line.contains("begin time") {
                t0 = capture_field(line, "coords")?;
            }

            if line.starts_with("    MBAR - lambda values considered:") {
                in_mbar = true;
                continue;
            }
            if in_mbar {
                if line.starts_with("    Extra") {
                    in_mbar = false;
                    continue;
                }
                parse_mbar_lambda_line(line, &mut mbar_lambdas)?;
            }
        }

        Ok(UNkHeader {
            temp0: temp0.ok_or(AmberParseError::MissingField { field: "temp0" })?,
            clambda: clambda.ok_or(AmberParseError::MissingField { field: "clambda" })?,
            dt: dt.ok_or(AmberParseError::MissingField { field: "dt" })?,
            bar_intervall: bar_intervall.ok_or(AmberParseError::MissingField {
                field: "bar_intervall",
            })?,
            t0: t0.ok_or(AmberParseError::MissingField { field: "t0" })?,
            mbar_lambdas,
        })
    }

    fn is_mbar_energy_start(line: &str) -> bool {
        line.trim_start().starts_with("MBAR Energy analysis:")
    }

    fn parse_mbar_lambda_line(line: &str, mbar_lambdas: &mut Vec<f64>) -> Result<()> {
        for token in line.split_whitespace() {
            if token.contains('.') {
                if let Ok(value) = token.parse::<f64>() {
                    mbar_lambdas.push(value);
                }
            }
        }
        Ok(())
    }

    fn float_regex() -> &'static regex::Regex {
        static FLOAT_REGEX: OnceLock<regex::Regex> = OnceLock::new();
        FLOAT_REGEX.get_or_init(|| regex::Regex::new(FP_RE).expect("valid float regex"))
    }

    fn lambda_key(value: f64) -> i64 {
        (value * 10_000.0).round() as i64
    }

    fn build_lambda_index(mbar_lambdas: &[f64]) -> std::collections::HashMap<i64, usize> {
        let mut index = std::collections::HashMap::with_capacity(mbar_lambdas.len());
        for (idx, lambda) in mbar_lambdas.iter().enumerate() {
            index.insert(lambda_key(*lambda), idx);
        }
        index
    }

    fn parse_mbar_energy_line(
        line: &str,
        sample_idx: usize,
        lambda_index: &std::collections::HashMap<i64, usize>,
        energies: &mut [Option<f64>],
    ) -> Result<()> {
        if !line.contains("Energy at") {
            return Ok(());
        }

        let mut numbers = float_regex().find_iter(line);
        let Some(lambda_match) = numbers.next() else {
            return Err(AmberParseError::InvalidMbarEnergyLine {
                line: line.to_string(),
            });
        };
        let lambda_token = lambda_match.as_str();
        let lambda =
            lambda_token
                .parse::<f64>()
                .map_err(|_| AmberParseError::InvalidMbarEnergyValue {
                    line: line.to_string(),
                    value: lambda_token.to_string(),
                })?;
        let energy = if let Some(energy_match) = numbers.next() {
            let energy_token = energy_match.as_str();
            energy_token
                .parse::<f64>()
                .map_err(|_| AmberParseError::InvalidMbarEnergyValue {
                    line: line.to_string(),
                    value: energy_token.to_string(),
                })?
        } else if line.contains('*') {
            f64::INFINITY
        } else {
            return Err(AmberParseError::InvalidMbarEnergyLine {
                line: line.to_string(),
            });
        };

        if let Some(&idx) = lambda_index.get(&lambda_key(lambda)) {
            if energies[idx].is_some() {
                return Err(AmberParseError::DuplicateMbarEnergy { sample_idx, lambda });
            }
            energies[idx] = Some(energy);
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn finalize_mbar_block(
        energies: &[Option<f64>],
        lambdas: &[f64],
        reference_idx: usize,
        beta: f64,
        sample_idx: usize,
        t0: f64,
        dt: f64,
        bar_intervall: f64,
        reduced_row: &mut Vec<f64>,
        time_ps: &mut Vec<f64>,
        data: &mut Vec<f64>,
    ) -> Result<()> {
        if energies.iter().any(|value| value.is_none()) {
            let missing_lambdas = energies
                .iter()
                .zip(lambdas.iter().copied())
                .filter_map(|(value, lambda)| value.is_none().then_some(lambda))
                .collect();
            return Err(AmberParseError::IncompleteMbarBlock {
                sample_idx,
                missing_lambdas,
            });
        }

        let reference = energies[reference_idx].unwrap();
        if !reference.is_finite() {
            return Ok(());
        }
        reduced_row.clear();
        reduced_row.reserve(energies.len());
        for energy in energies.iter().map(|value| value.unwrap()) {
            let reduced = beta * (energy - reference);
            if reduced.is_nan() || reduced == f64::NEG_INFINITY {
                return Ok(());
            }
            reduced_row.push(reduced);
        }

        time_ps.push(t0 + ((sample_idx + 1) as f64) * dt * bar_intervall);
        data.extend(reduced_row.iter().copied());
        Ok(())
    }

    fn handle_gradient_block(
        block: &str,
        last_nstep: &mut Option<i64>,
        gradients: &mut Vec<f64>,
    ) -> Result<()> {
        let nstep = capture_field(block, "NSTEP")?
            .map(|value| value as i64)
            .ok_or(AmberParseError::MissingGradientField { field: "NSTEP" })?;
        if let Some(prev) = last_nstep {
            if *prev == nstep {
                return Ok(());
            }
        }
        let dvdl = capture_field(block, "DV/DL")?
            .ok_or(AmberParseError::MissingGradientField { field: "DV/DL" })?;
        gradients.push(dvdl);
        *last_nstep = Some(nstep);
        Ok(())
    }

    fn handle_potential_block(
        block: &str,
        last_nstep: &mut Option<i64>,
        potential_samples: &mut Vec<f64>,
    ) -> Result<()> {
        let nstep = capture_field(block, "NSTEP")?
            .map(|value| value as i64)
            .ok_or(AmberParseError::MissingPotentialField { field: "NSTEP" })?;
        if let Some(prev) = last_nstep {
            if *prev == nstep {
                return Ok(());
            }
        }
        let epot = capture_field(block, "EPtot")?
            .ok_or(AmberParseError::MissingPotentialField { field: "EPtot" })?;
        potential_samples.push(epot);
        *last_nstep = Some(nstep);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::amber::extract_dhdl;
    use super::amber::extract_u_nk;
    use super::amber::extract_u_nk_with_potential;
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
 NSTEP =       10  TIME(PS) =       0.02  DV/DL =     1.0000
 ---
 NSTEP =       20  TIME(PS) =       0.04  DV/DL =     2.0000
 ---
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
