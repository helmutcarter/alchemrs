use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use alchemrs_core::{CoreError, DhdlSeries, Result, StatePoint};

pub mod amber {
    use super::*;

    const FP_RE: &str = r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?";
    const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;

    pub fn extract_dhdl(path: impl AsRef<Path>, temperature_k: f64) -> Result<DhdlSeries> {
        let file = File::open(path.as_ref())
            .map_err(|err| CoreError::Parse(format!("failed to open file: {err}")))?;
        let reader = BufReader::new(file);

        let mut temp0: Option<f64> = None;
        let mut clambda: Option<f64> = None;
        let mut dt: Option<f64> = None;
        let mut ntpr: Option<f64> = None;
        let mut t0: Option<f64> = None;

        let mut gradients: Vec<f64> = Vec::new();
        let mut last_nstep: Option<i64> = None;
        let mut pending_block: Option<String> = None;

        for line in reader.lines() {
            let line = line.map_err(|err| CoreError::Parse(format!("read error: {err}")))?;

            if temp0.is_none() {
                if let Some(value) = capture_field(&line, "temp0")? {
                    temp0 = Some(value);
                }
            }
            if clambda.is_none() {
                if let Some(value) = capture_field(&line, "clambda")? {
                    clambda = Some(value);
                }
            }
            if dt.is_none() {
                if let Some(value) = capture_field(&line, "dt")? {
                    dt = Some(value);
                }
            }
            if ntpr.is_none() {
                if let Some(value) = capture_field(&line, "ntpr")? {
                    ntpr = Some(value);
                }
            }
            if t0.is_none() && line.contains("begin time") {
                if let Some(value) = capture_field(&line, "coords")? {
                    t0 = Some(value);
                }
            }

            if let Some(mut block) = pending_block.take() {
                block.push(' ');
                block.push_str(&line);
                if line.trim_start().starts_with("---") {
                    handle_gradient_block(&block, &mut last_nstep, &mut gradients)?;
                } else {
                    pending_block = Some(block);
                }
                continue;
            }

            if line.starts_with(" NSTEP") || line.starts_with("NSTEP") {
                let block = line.clone();
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

        let temp0 = temp0.ok_or_else(|| {
            CoreError::Parse("failed to locate temp0 in AMBER output".to_string())
        })?;
        if (temperature_k - temp0).abs() > 1e-2 {
            return Err(CoreError::Parse(format!(
                "temperature mismatch: file {temp0:.2} K vs input {temperature_k:.2} K"
            )));
        }
        let clambda = clambda.ok_or_else(|| {
            CoreError::Parse("failed to locate clambda in AMBER output".to_string())
        })?;
        let dt = dt.ok_or_else(|| CoreError::Parse("failed to locate dt".to_string()))?;
        let ntpr =
            ntpr.ok_or_else(|| CoreError::Parse("failed to locate ntpr".to_string()))?;
        let t0 = t0.ok_or_else(|| CoreError::Parse("failed to locate t0".to_string()))?;

        if gradients.is_empty() {
            return Err(CoreError::Parse(
                "no DV/DL gradients found in AMBER output".to_string(),
            ));
        }

        let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);
        for value in gradients.iter_mut() {
            *value *= beta;
        }

        let time_ps: Vec<f64> = (0..gradients.len())
            .map(|idx| t0 + ((idx + 1) as f64) * dt * ntpr)
            .collect();

        let state = StatePoint::new(vec![clambda], temperature_k)?;
        DhdlSeries::new(state, time_ps, gradients)
    }

    pub fn extract_u_nk(path: impl AsRef<Path>, temperature_k: f64) -> Result<alchemrs_core::UNkMatrix> {
        let file = File::open(path.as_ref())
            .map_err(|err| CoreError::Parse(format!("failed to open file: {err}")))?;
        let reader = BufReader::new(file);
        let lines: Vec<String> = reader
            .lines()
            .collect::<std::result::Result<_, _>>()
            .map_err(|err| CoreError::Parse(format!("read error: {err}")))?;

        let temp0 = find_first_field(&lines, "temp0")?
            .ok_or_else(|| CoreError::Parse("failed to locate temp0 in AMBER output".to_string()))?;
        if (temperature_k - temp0).abs() > 1e-2 {
            return Err(CoreError::Parse(format!(
                "temperature mismatch: file {temp0:.2} K vs input {temperature_k:.2} K"
            )));
        }

        let clambda = find_first_field(&lines, "clambda")?
            .ok_or_else(|| CoreError::Parse("failed to locate clambda in AMBER output".to_string()))?;
        let dt = find_first_field(&lines, "dt")?
            .ok_or_else(|| CoreError::Parse("failed to locate dt".to_string()))?;
        let bar_intervall = find_first_field(&lines, "bar_intervall")?
            .ok_or_else(|| CoreError::Parse("failed to locate bar_intervall".to_string()))?;
        let t0 = find_begin_time(&lines)?
            .ok_or_else(|| CoreError::Parse("failed to locate t0".to_string()))?;

        let mbar_lambdas = extract_mbar_lambdas(&lines)?;
        if mbar_lambdas.is_empty() {
            return Err(CoreError::Parse(
                "no MBAR lambda values found in AMBER output".to_string(),
            ));
        }

        let clambda_key = format!("{clambda:.4}");
        let reference_idx = mbar_lambdas
            .iter()
            .position(|value| value == &clambda_key)
            .ok_or_else(|| {
                CoreError::Parse(
                    "clambda not found in MBAR lambda list; cannot build u_nk".to_string(),
                )
            })?;

        let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);
        let mut mbar_energies: Vec<Vec<f64>> = vec![Vec::new(); mbar_lambdas.len()];

        let mut idx = 0;
        while idx < lines.len() {
            let line = &lines[idx];
            if is_mbar_energy_start(line) {
                let mut block_lines = Vec::new();
                idx += 1;
                while idx < lines.len() {
                    let next = &lines[idx];
                    if next.trim_start().starts_with("---") {
                        break;
                    }
                    block_lines.push(next.clone());
                    idx += 1;
                }
                let energies = parse_mbar_block(&block_lines, &mbar_lambdas)?;
                let reference = energies[reference_idx];
                for (state_idx, energy) in energies.into_iter().enumerate() {
                    mbar_energies[state_idx].push(beta * (energy - reference));
                }
            }
            idx += 1;
        }

        if mbar_energies.iter().all(|values| values.is_empty()) {
            return Err(CoreError::Parse(
                "no MBAR energy samples found in AMBER output".to_string(),
            ));
        }

        let n_states = mbar_lambdas.len();
        let n_samples = mbar_energies[0].len();
        for (idx, values) in mbar_energies.iter().enumerate() {
            if values.len() != n_samples {
                return Err(CoreError::InvalidShape {
                    expected: n_samples,
                    found: values.len(),
                });
            }
            if values.is_empty() {
                return Err(CoreError::Parse(format!(
                    "no MBAR samples collected for state {idx}"
                )));
            }
        }

        let mut valid_indices = Vec::new();
        (0..n_samples).for_each(|sample_idx| {
            if (0..n_states).all(|state_idx| {
                mbar_energies[state_idx][sample_idx].is_finite()
            }) {
                valid_indices.push(sample_idx);
            }
        });

        if valid_indices.is_empty() {
            return Err(CoreError::Parse(
                "all MBAR samples contained non-finite values".to_string(),
            ));
        }

        let n_samples = valid_indices.len();
        let time_ps: Vec<f64> = valid_indices
            .iter()
            .map(|idx| t0 + ((*idx + 1) as f64) * dt * bar_intervall)
            .collect();

        let mut data = Vec::with_capacity(valid_indices.len() * n_states);
        for &sample_idx in &valid_indices {
            (0..n_states).for_each(|state_idx| {
                data.push(mbar_energies[state_idx][sample_idx]);
            });
        }

        let sampled_state = Some(StatePoint::new(vec![clambda], temperature_k)?);
        let evaluated_states = mbar_lambdas
            .iter()
            .map(|value| {
                let lambda = value
                    .trim()
                    .parse::<f64>()
                    .map_err(|err| CoreError::Parse(format!("invalid MBAR lambda: {err}")))?;
                StatePoint::new(vec![lambda], temperature_k)
            })
            .collect::<Result<Vec<_>>>()?;

        alchemrs_core::UNkMatrix::new(
            n_samples,
            n_states,
            data,
            time_ps,
            sampled_state,
            evaluated_states,
        )
    }

    fn capture_field(line: &str, field: &str) -> Result<Option<f64>> {
        let mut iter = line.split_whitespace().peekable();
        while let Some(token) = iter.next() {
            if token == field {
                if let Some(next) = iter.next() {
                    let value = if next == "=" {
                        iter.next()
                            .ok_or_else(|| CoreError::Parse(format!("missing {field} value")))? 
                    } else {
                        next
                    };
                    let value = value.trim_matches(|c: char| !c.is_ascii_digit() && c != '.' && c != '-' && c != '+' && c != 'e' && c != 'E');
                    let parsed = value
                        .parse::<f64>()
                        .map_err(|err| CoreError::Parse(format!("invalid {field}: {err}")))?;
                    return Ok(Some(parsed));
                }
            } else if token.starts_with(field) {
                if let Some(eq_idx) = token.find('=') {
                    let value = &token[(eq_idx + 1)..];
                    if !value.is_empty() {
                        let value = value.trim_matches(|c: char| !c.is_ascii_digit() && c != '.' && c != '-' && c != '+' && c != 'e' && c != 'E');
                        let parsed = value
                            .parse::<f64>()
                            .map_err(|err| CoreError::Parse(format!("invalid {field}: {err}")))?;
                        return Ok(Some(parsed));
                    }
                }
            }
        }
        Ok(None)
    }

    fn find_first_field(lines: &[String], field: &str) -> Result<Option<f64>> {
        for line in lines {
            if let Some(value) = capture_field(line, field)? {
                return Ok(Some(value));
            }
        }
        Ok(None)
    }

    fn find_begin_time(lines: &[String]) -> Result<Option<f64>> {
        for line in lines {
            if line.contains("begin time") {
                if let Some(value) = capture_field(line, "coords")? {
                    return Ok(Some(value));
                }
            }
        }
        Ok(None)
    }

    fn extract_mbar_lambdas(lines: &[String]) -> Result<Vec<String>> {
        let mut in_mbar = false;
        let mut lambdas = Vec::new();
        for line in lines {
            if line.starts_with("    MBAR - lambda values considered:") {
                in_mbar = true;
                continue;
            }
            if in_mbar {
                if line.starts_with("    Extra") {
                    break;
                }
                for token in line.split_whitespace() {
                    if token.contains('.') && token.parse::<f64>().is_ok() {
                        lambdas.push(token.to_string());
                    }
                }
            }
        }
        Ok(lambdas)
    }

    fn is_mbar_energy_start(line: &str) -> bool {
        line.trim_start().starts_with("MBAR Energy analysis:")
    }

    fn parse_mbar_block(lines: &[String], mbar_lambdas: &[String]) -> Result<Vec<f64>> {
        let mut energies: Vec<Option<f64>> = vec![None; mbar_lambdas.len()];
        let mbar_regex = regex::Regex::new(FP_RE)
            .map_err(|err| CoreError::Parse(format!("regex error: {err}")))?;

        for line in lines {
            if !line.contains("Energy at") {
                continue;
            }
            let mut numbers = Vec::new();
            for caps in mbar_regex.find_iter(line) {
                numbers.push(
                    caps.as_str()
                        .parse::<f64>()
                        .map_err(|err| CoreError::Parse(format!("invalid MBAR energy: {err}")))?,
                );
            }
            if numbers.is_empty() {
                continue;
            }
            let lambda = numbers[0];
            let energy = if numbers.len() >= 2 {
                numbers[1]
            } else if line.contains('*') {
                f64::INFINITY
            } else {
                continue;
            };
            let lambda_key = format!("{lambda:.4}");
            if let Some(idx) = mbar_lambdas.iter().position(|val| val == &lambda_key) {
                energies[idx] = Some(energy);
            }
        }

        if energies.iter().any(|value| value.is_none()) {
            return Err(CoreError::Parse(
                "missing MBAR energy entries in block".to_string(),
            ));
        }
        Ok(energies.into_iter().map(|value| value.unwrap()).collect())
    }

    fn handle_gradient_block(
        block: &str,
        last_nstep: &mut Option<i64>,
        gradients: &mut Vec<f64>,
    ) -> Result<()> {
        let nstep = capture_field(block, "NSTEP")?
            .map(|value| value as i64)
            .ok_or_else(|| CoreError::Parse("missing NSTEP in block".to_string()))?;
        if let Some(prev) = last_nstep {
            if *prev == nstep {
                return Ok(());
            }
        }
        let dvdl = capture_field(block, "DV/DL")?.ok_or_else(|| {
            CoreError::Parse("missing DV/DL in NSTEP block".to_string())
        })?;
        gradients.push(dvdl);
        *last_nstep = Some(nstep);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::amber::extract_dhdl;
    use super::amber::extract_u_nk;
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
}
