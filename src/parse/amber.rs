use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::OnceLock;

use crate::data::{DhdlSeries, StatePoint, SwitchingTrajectory, UNkMatrix};
use crate::error::CoreError;
use thiserror::Error;

const FP_RE: &str = r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?";
const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;

pub type Result<T> = std::result::Result<T, AmberParseError>;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AmberDhdlOptions {
    pub input_stride: Option<usize>,
}

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
    MissingPotentialField { field: &'static str },
    #[error("no DV/DL gradients found in AMBER output")]
    NoGradients,
    #[error("no DV/DL summary block found in AMBER output")]
    NoDvdlSummary,
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
    #[error("missing MBAR energy entries in block {sample_idx} for lambdas {missing_lambdas:?}")]
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
    extract_dhdl_with_options(path, temperature_k, AmberDhdlOptions::default())
}

pub fn extract_dhdl_with_options(
    path: impl AsRef<Path>,
    temperature_k: f64,
    options: AmberDhdlOptions,
) -> Result<DhdlSeries> {
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
    let mut time_ps: Vec<f64> = Vec::new();
    let mut summary_state: Option<DvdlSummaryState> = None;
    let mut summary_dt: Option<f64> = None;
    let mut summary_t0: Option<f64> = None;
    let mut line = Vec::new();

    loop {
        line.clear();
        let bytes = reader
            .read_until(b'\n', &mut line)
            .map_err(|err| AmberParseError::Io {
                operation: "read",
                message: err.to_string(),
            })?;
        if bytes == 0 {
            break;
        }
        let line = trim_line_end_bytes(&line);

        if let Some(state) = summary_state.as_mut() {
            if is_dvdl_summary_end(line) {
                break;
            }
            if line.is_empty() {
                continue;
            }
            if let Some(step) = state.next_retained_step() {
                let value = parse_dvdl_summary_value_bytes(line)?;
                gradients.push(value);
                time_ps.push(
                    summary_t0.expect("summary t0 set")
                        + (step as f64) * summary_dt.expect("summary dt set"),
                );
            }
            continue;
        }

        if is_dvdl_summary_start(line) {
            let ntpr = ntpr.ok_or(AmberParseError::MissingField { field: "ntpr" })?;
            let dt = dt.ok_or(AmberParseError::MissingField { field: "dt" })?;
            let t0 = t0.ok_or(AmberParseError::MissingField { field: "t0" })?;
            let summary = DvdlSummaryState::new(
                parse_dvdl_summary_count_bytes(line)?,
                summary_stride(ntpr),
                options.input_stride,
            );
            gradients.reserve(summary.retained_len());
            time_ps.reserve(summary.retained_len());
            summary_dt = Some(dt);
            summary_t0 = Some(t0);
            summary_state = Some(summary);
            continue;
        }

        let needs_header_parse =
            temp0.is_none() || clambda.is_none() || dt.is_none() || ntpr.is_none() || t0.is_none();
        if !needs_header_parse {
            continue;
        }

        let line = parse_line_bytes(line)?;
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
        if t0.is_none() && line_contains_begin_time(line.as_bytes()) {
            if let Some(value) = capture_field(line, "coords")? {
                t0 = Some(value);
            }
        }
    }

    let temp0 = temp0.ok_or(AmberParseError::MissingField { field: "temp0" })?;
    if (temperature_k - temp0).abs() > 1e-2 {
        return Err(AmberParseError::TemperatureMismatch {
            file_temperature_k: temp0,
            input_temperature_k: temperature_k,
        });
    }
    let clambda = clambda.ok_or(AmberParseError::MissingField { field: "clambda" })?;
    let _dt = dt.ok_or(AmberParseError::MissingField { field: "dt" })?;
    ntpr.ok_or(AmberParseError::MissingField { field: "ntpr" })?;
    let _t0 = t0.ok_or(AmberParseError::MissingField { field: "t0" })?;

    if gradients.is_empty() {
        return Err(AmberParseError::NoGradients);
    }

    let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);
    for value in gradients.iter_mut() {
        *value *= beta;
    }

    let state = StatePoint::new(vec![clambda], temperature_k)?;
    DhdlSeries::new(state, time_ps, gradients).map_err(Into::into)
}

pub fn extract_temperature(path: impl AsRef<Path>) -> Result<f64> {
    read_temperature(path.as_ref())
}

pub fn extract_nes_trajectory(
    path: impl AsRef<Path>,
    temperature_k: f64,
) -> Result<SwitchingTrajectory> {
    let file = File::open(path.as_ref()).map_err(|err| AmberParseError::Io {
        operation: "open",
        message: err.to_string(),
    })?;
    let mut reader = BufReader::new(file);

    let mut temp0: Option<f64> = None;
    let mut clambda: Option<f64> = None;
    let mut dynlmb: Option<f64> = None;
    let mut ntave: Option<f64> = None;
    let mut summary_expected_steps: Option<usize> = None;
    let mut summary_steps = 0usize;
    let mut summary_dvdl_sum = 0.0f64;
    let mut in_summary = false;
    let mut line = Vec::new();

    loop {
        line.clear();
        let bytes = reader
            .read_until(b'\n', &mut line)
            .map_err(|err| AmberParseError::Io {
                operation: "read",
                message: err.to_string(),
            })?;
        if bytes == 0 {
            break;
        }
        let line = trim_line_end_bytes(&line);

        if in_summary {
            if is_dvdl_summary_end(line) {
                break;
            }
            if line.is_empty() {
                continue;
            }
            summary_dvdl_sum += parse_dvdl_summary_value_bytes(line)?;
            summary_steps += 1;
            continue;
        }

        if is_dvdl_summary_start(line) {
            summary_expected_steps = Some(parse_dvdl_summary_count_bytes(line)?);
            in_summary = true;
            continue;
        }

        if temp0.is_some() && clambda.is_some() && dynlmb.is_some() && ntave.is_some() {
            continue;
        }

        let line = parse_line_bytes(line)?;
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
        if dynlmb.is_none() {
            if let Some(value) = capture_field(line, "dynlmb")? {
                dynlmb = Some(value);
            }
        }
        if ntave.is_none() {
            if let Some(value) = capture_field(line, "ntave")? {
                ntave = Some(value);
            }
        }
    }

    let temp0 = temp0.ok_or(AmberParseError::MissingField { field: "temp0" })?;
    if (temperature_k - temp0).abs() > 1e-2 {
        return Err(AmberParseError::TemperatureMismatch {
            file_temperature_k: temp0,
            input_temperature_k: temperature_k,
        });
    }
    let clambda = clambda.ok_or(AmberParseError::MissingField { field: "clambda" })?;
    let dynlmb = dynlmb.ok_or(AmberParseError::MissingField { field: "dynlmb" })?;
    let ntave = ntave.ok_or(AmberParseError::MissingField { field: "ntave" })?;
    if summary_expected_steps.is_none() {
        return Err(AmberParseError::NoDvdlSummary);
    }
    if summary_steps == 0 {
        return Err(AmberParseError::NoGradients);
    }
    if summary_expected_steps != Some(summary_steps) {
        return Err(AmberParseError::InvalidParsedData {
            message: format!(
                "DV/DL summary count mismatch: expected {} values but found {}",
                summary_expected_steps.expect("checked above"),
                summary_steps
            ),
        });
    }

    let step_lambda = dynlmb / ntave;
    let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);
    let reduced_work = summary_dvdl_sum * step_lambda * beta;
    let initial_state = StatePoint::new(vec![clambda], temperature_k)?;
    let final_state = StatePoint::new(
        vec![clambda + step_lambda * summary_steps as f64],
        temperature_k,
    )?;
    SwitchingTrajectory::new(initial_state, final_state, reduced_work).map_err(Into::into)
}

pub fn extract_u_nk(path: impl AsRef<Path>, temperature_k: f64) -> Result<UNkMatrix> {
    let (u_nk, _potential) = extract_u_nk_internal(path.as_ref(), temperature_k, false)?;
    Ok(u_nk)
}

pub fn extract_u_nk_with_potential(
    path: impl AsRef<Path>,
    temperature_k: f64,
) -> Result<(UNkMatrix, Vec<f64>)> {
    let (u_nk, potential) = extract_u_nk_internal(path.as_ref(), temperature_k, true)?;
    Ok((u_nk, potential.expect("potential samples requested")))
}

fn extract_u_nk_internal(
    path: &Path,
    temperature_k: f64,
    include_potential: bool,
) -> Result<(UNkMatrix, Option<Vec<f64>>)> {
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

    let u_nk = UNkMatrix::new(
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

fn read_temperature(path: &Path) -> Result<f64> {
    let file = File::open(path).map_err(|err| AmberParseError::Io {
        operation: "open",
        message: err.to_string(),
    })?;
    let mut reader = BufReader::new(file);
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
        if let Some(temp0) = capture_field(line, "temp0")? {
            return Ok(temp0);
        }
    }

    Err(AmberParseError::MissingField { field: "temp0" })
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

fn trim_line_end_bytes(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && matches!(line[end - 1], b'\n' | b'\r') {
        end -= 1;
    }
    &line[..end]
}

fn parse_line_bytes(line: &[u8]) -> Result<&str> {
    std::str::from_utf8(line).map_err(|err| AmberParseError::Io {
        operation: "read",
        message: format!("invalid UTF-8 in AMBER output: {err}"),
    })
}

fn trim_ascii_start_bytes(line: &[u8]) -> &[u8] {
    let start = line
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .unwrap_or(line.len());
    &line[start..]
}

fn line_contains_begin_time(line: &[u8]) -> bool {
    line.windows(b"begin time".len())
        .any(|window| window == b"begin time")
}

fn is_dvdl_summary_start(line: &[u8]) -> bool {
    trim_ascii_start_bytes(line).starts_with(b"Summary of dvdl values over")
}

fn is_dvdl_summary_end(line: &[u8]) -> bool {
    trim_ascii_start_bytes(line).starts_with(b"End of dvdl summary")
}

fn trim_ascii_bytes(line: &[u8]) -> &[u8] {
    let mut start = 0;
    let mut end = line.len();
    while start < end && line[start].is_ascii_whitespace() {
        start += 1;
    }
    while end > start && line[end - 1].is_ascii_whitespace() {
        end -= 1;
    }
    &line[start..end]
}

fn parse_dvdl_summary_count_bytes(line: &[u8]) -> Result<usize> {
    let trimmed = trim_ascii_bytes(line);
    let mut idx = 0;
    while idx < trimmed.len() {
        if trimmed[idx].is_ascii_digit() {
            let start = idx;
            idx += 1;
            while idx < trimmed.len() && trimmed[idx].is_ascii_digit() {
                idx += 1;
            }
            return parse_ascii_usize(&trimmed[start..idx]).ok_or_else(|| {
                AmberParseError::InvalidField {
                    field: "dvdl_summary_steps",
                    value: String::from_utf8_lossy(trimmed).into_owned(),
                }
            });
        }
        idx += 1;
    }
    Err(AmberParseError::InvalidField {
        field: "dvdl_summary_steps",
        value: String::from_utf8_lossy(trimmed).into_owned(),
    })
}

fn parse_dvdl_summary_value_bytes(line: &[u8]) -> Result<f64> {
    let trimmed = trim_ascii_bytes(line);
    parse_ascii_f64(trimmed).ok_or_else(|| AmberParseError::InvalidField {
        field: "DV/DL",
        value: String::from_utf8_lossy(trimmed).into_owned(),
    })
}

fn parse_ascii_usize(bytes: &[u8]) -> Option<usize> {
    if bytes.is_empty() {
        return None;
    }
    let mut value = 0usize;
    for byte in bytes {
        if !byte.is_ascii_digit() {
            return None;
        }
        value = value.checked_mul(10)?.checked_add((byte - b'0') as usize)?;
    }
    Some(value)
}

fn parse_ascii_f64(bytes: &[u8]) -> Option<f64> {
    if bytes.is_empty() {
        return None;
    }

    let mut idx = 0usize;
    let mut sign = 1.0f64;
    match bytes[idx] {
        b'+' => idx += 1,
        b'-' => {
            sign = -1.0;
            idx += 1;
        }
        _ => {}
    }
    if idx >= bytes.len() {
        return None;
    }

    let mut value = 0.0f64;
    let mut saw_digit = false;
    while idx < bytes.len() && bytes[idx].is_ascii_digit() {
        value = value * 10.0 + (bytes[idx] - b'0') as f64;
        idx += 1;
        saw_digit = true;
    }

    if idx < bytes.len() && bytes[idx] == b'.' {
        idx += 1;
        let mut scale = 0.1f64;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            value += (bytes[idx] - b'0') as f64 * scale;
            scale *= 0.1;
            idx += 1;
            saw_digit = true;
        }
    }

    if !saw_digit {
        return None;
    }

    if idx < bytes.len() && matches!(bytes[idx], b'e' | b'E') {
        idx += 1;
        if idx >= bytes.len() {
            return None;
        }
        let mut exponent_sign = 1i32;
        match bytes[idx] {
            b'+' => idx += 1,
            b'-' => {
                exponent_sign = -1;
                idx += 1;
            }
            _ => {}
        }
        if idx >= bytes.len() || !bytes[idx].is_ascii_digit() {
            return None;
        }
        let mut exponent = 0i32;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            exponent = exponent
                .checked_mul(10)?
                .checked_add((bytes[idx] - b'0') as i32)?;
            idx += 1;
        }
        value *= 10f64.powi(exponent_sign * exponent);
    }

    if idx != bytes.len() {
        return None;
    }

    Some(sign * value)
}

fn summary_stride(ntpr: f64) -> usize {
    ntpr.round().max(1.0) as usize
}

#[derive(Debug, Clone, Copy)]
enum DvdlSummaryMode {
    Dense { stride: usize },
    Compact { stride: usize },
}

#[derive(Debug, Clone, Copy)]
struct DvdlSummaryState {
    mode: DvdlSummaryMode,
    total_steps: usize,
    index: usize,
}

impl DvdlSummaryState {
    fn new(total_steps: usize, default_stride: usize, input_stride: Option<usize>) -> Self {
        let mode = if input_stride.is_none() && total_steps <= default_stride {
            DvdlSummaryMode::Compact {
                stride: default_stride,
            }
        } else {
            DvdlSummaryMode::Dense {
                stride: input_stride.unwrap_or(default_stride),
            }
        };
        Self {
            mode,
            total_steps,
            index: 0,
        }
    }

    fn next_retained_step(&mut self) -> Option<usize> {
        let retained = match self.mode {
            DvdlSummaryMode::Dense { stride } => (self.index % stride == 0).then_some(self.index),
            DvdlSummaryMode::Compact { stride } => Some((self.index + 1) * stride),
        };
        self.index += 1;
        retained
    }

    fn retained_len(&self) -> usize {
        match self.mode {
            DvdlSummaryMode::Dense { stride } => self.total_steps.div_ceil(stride),
            DvdlSummaryMode::Compact { .. } => self.total_steps,
        }
    }
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
