use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::OnceLock;

use thiserror::Error;

use crate::data::{DhdlSeries, StatePoint, UNkMatrix};
use crate::error::CoreError;

const FP_RE: &str = r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?";
const K_B_KCAL_PER_MOL_K: f64 = 0.00198720425864083;

pub type Result<T> = std::result::Result<T, GromacsParseError>;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum GromacsParseError {
    #[error("failed to {operation} GROMACS output: {message}")]
    Io {
        operation: &'static str,
        message: String,
    },
    #[error("failed to locate {field} in GROMACS dhdl.xvg output")]
    MissingField { field: &'static str },
    #[error("invalid {field} value `{value}` in GROMACS dhdl.xvg output")]
    InvalidField { field: &'static str, value: String },
    #[error(
        "temperature mismatch: file {file_temperature_k:.2} K vs input {input_temperature_k:.2} K"
    )]
    TemperatureMismatch {
        file_temperature_k: f64,
        input_temperature_k: f64,
    },
    #[error("missing dH/dlambda legend in GROMACS dhdl.xvg output")]
    MissingDhdlLegend,
    #[error("no Delta H legends found in GROMACS dhdl.xvg output")]
    NoDeltaHColumns,
    #[error("no potential-energy legend found in GROMACS dhdl.xvg output")]
    NoPotentialSamples,
    #[error("duplicate lambda legend for state {lambda:.4}")]
    DuplicateLambda { lambda: f64 },
    #[error("invalid legend line `{line}`")]
    InvalidLegendLine { line: String },
    #[error("invalid data line `{line}`")]
    InvalidDataLine { line: String },
    #[error("invalid data value `{value}` in line `{line}`")]
    InvalidDataValue { line: String, value: String },
    #[error("no numerical samples found in GROMACS dhdl.xvg output")]
    NoSamples,
    #[error("invalid parsed data: {message}")]
    InvalidParsedData { message: String },
}

impl From<GromacsParseError> for CoreError {
    fn from(error: GromacsParseError) -> Self {
        CoreError::Parse(error.to_string())
    }
}

impl From<CoreError> for GromacsParseError {
    fn from(error: CoreError) -> Self {
        GromacsParseError::InvalidParsedData {
            message: error.to_string(),
        }
    }
}

#[derive(Debug, Clone)]
struct Header {
    temperature_k: f64,
    sampled_lambda: f64,
    dhdl_column: usize,
    delta_h_columns: Vec<(usize, f64)>,
    potential_column: Option<usize>,
}

pub fn extract_dhdl(path: impl AsRef<Path>, temperature_k: f64) -> Result<DhdlSeries> {
    let header = read_header(path.as_ref())?;
    validate_temperature(header.temperature_k, temperature_k)?;

    let mut time_ps = Vec::new();
    let mut values = Vec::new();
    for row in read_data_rows(path.as_ref())? {
        time_ps.push(row.time_ps);
        values.push(
            *row.values
                .get(header.dhdl_column)
                .ok_or_else(|| GromacsParseError::InvalidDataLine {
                    line: row.raw_line.clone(),
                })?,
        );
    }

    if values.is_empty() {
        return Err(GromacsParseError::NoSamples);
    }

    let state = StatePoint::new(vec![header.sampled_lambda], temperature_k)?;
    DhdlSeries::new(state, time_ps, values).map_err(Into::into)
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
    let header = read_header(path)?;
    validate_temperature(header.temperature_k, temperature_k)?;

    if header.delta_h_columns.is_empty() {
        return Err(GromacsParseError::NoDeltaHColumns);
    }
    if include_potential && header.potential_column.is_none() {
        return Err(GromacsParseError::NoPotentialSamples);
    }

    let mut ordered_lambdas = vec![header.sampled_lambda];
    ordered_lambdas.extend(header.delta_h_columns.iter().map(|(_, lambda)| *lambda));
    ordered_lambdas.sort_by(|left, right| left.partial_cmp(right).expect("finite lambda"));

    let mut lambda_to_column = HashMap::with_capacity(ordered_lambdas.len());
    for (idx, lambda) in ordered_lambdas.iter().copied().enumerate() {
        if lambda_to_column.insert(lambda_key(lambda), idx).is_some() {
            return Err(GromacsParseError::DuplicateLambda { lambda });
        }
    }

    let sampled_idx = *lambda_to_column
        .get(&lambda_key(header.sampled_lambda))
        .expect("sampled lambda indexed");
    let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);

    let mut time_ps = Vec::new();
    let mut data = Vec::new();
    let mut potential = include_potential.then(Vec::new);

    for row in read_data_rows(path)? {
        let mut reduced = vec![0.0; ordered_lambdas.len()];
        reduced[sampled_idx] = 0.0;
        for (series_idx, foreign_lambda) in &header.delta_h_columns {
            let raw = *row.values.get(*series_idx).ok_or_else(|| {
                GromacsParseError::InvalidDataLine {
                    line: row.raw_line.clone(),
                }
            })?;
            let target_idx = *lambda_to_column
                .get(&lambda_key(*foreign_lambda))
                .expect("foreign lambda indexed");
            reduced[target_idx] = beta * raw;
        }
        if let Some(potential_column) = header.potential_column {
            if let Some(samples) = potential.as_mut() {
                samples.push(*row.values.get(potential_column).ok_or_else(|| {
                    GromacsParseError::InvalidDataLine {
                        line: row.raw_line.clone(),
                    }
                })?);
            }
        }
        time_ps.push(row.time_ps);
        data.extend(reduced);
    }

    if time_ps.is_empty() {
        return Err(GromacsParseError::NoSamples);
    }

    let sampled_state = Some(StatePoint::new(vec![header.sampled_lambda], temperature_k)?);
    let evaluated_states = ordered_lambdas
        .iter()
        .copied()
        .map(|lambda| StatePoint::new(vec![lambda], temperature_k))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(GromacsParseError::from)?;

    let u_nk = UNkMatrix::new(
        time_ps.len(),
        evaluated_states.len(),
        data,
        time_ps,
        sampled_state,
        evaluated_states,
    )
    .map_err(GromacsParseError::from)?;

    Ok((u_nk, potential))
}

fn validate_temperature(file_temperature_k: f64, temperature_k: f64) -> Result<()> {
    if (temperature_k - file_temperature_k).abs() > 1e-2 {
        return Err(GromacsParseError::TemperatureMismatch {
            file_temperature_k,
            input_temperature_k: temperature_k,
        });
    }
    Ok(())
}

fn read_header(path: &Path) -> Result<Header> {
    let file = File::open(path).map_err(|err| GromacsParseError::Io {
        operation: "open",
        message: err.to_string(),
    })?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();

    let mut temperature_k = None;
    let mut sampled_lambda = None;
    let mut dhdl_column = None;
    let mut delta_h_columns = Vec::new();
    let mut potential_column = None;

    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|err| GromacsParseError::Io {
                operation: "read",
                message: err.to_string(),
            })?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
        if !trimmed.starts_with('@') && !trimmed.starts_with('#') {
            continue;
        }

        if temperature_k.is_none() {
            temperature_k = parse_temperature(trimmed)?;
        }
        if sampled_lambda.is_none() {
            sampled_lambda = parse_header_lambda(trimmed)?;
        }

        let Some((series_idx, legend)) = parse_legend(trimmed)? else {
            continue;
        };

        if is_dhdl_legend(&legend) {
            dhdl_column = Some(series_idx);
            if sampled_lambda.is_none() {
                sampled_lambda = parse_lambda_from_text(&legend);
            }
            continue;
        }

        if is_delta_h_legend(&legend) {
            let lambda = parse_lambda_from_text(&legend).ok_or(GromacsParseError::InvalidField {
                field: "foreign lambda",
                value: legend.clone(),
            })?;
            delta_h_columns.push((series_idx, lambda));
            continue;
        }

        if potential_column.is_none() && is_potential_legend(&legend) {
            potential_column = Some(series_idx);
        }
    }

    Ok(Header {
        temperature_k: temperature_k.ok_or(GromacsParseError::MissingField {
            field: "temperature",
        })?,
        sampled_lambda: sampled_lambda.ok_or(GromacsParseError::MissingField {
            field: "lambda",
        })?,
        dhdl_column: dhdl_column.ok_or(GromacsParseError::MissingDhdlLegend)?,
        delta_h_columns,
        potential_column,
    })
}

#[derive(Debug, Clone)]
struct DataRow {
    raw_line: String,
    time_ps: f64,
    values: Vec<f64>,
}

fn read_data_rows(path: &Path) -> Result<Vec<DataRow>> {
    let file = File::open(path).map_err(|err| GromacsParseError::Io {
        operation: "open",
        message: err.to_string(),
    })?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut rows = Vec::new();

    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|err| GromacsParseError::Io {
                operation: "read",
                message: err.to_string(),
            })?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]).trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with('@')
            || trimmed.starts_with('&')
        {
            continue;
        }

        let mut fields = trimmed.split_whitespace();
        let time_token = fields.next().ok_or_else(|| GromacsParseError::InvalidDataLine {
            line: trimmed.to_string(),
        })?;
        let time_ps = parse_data_value(trimmed, time_token)?;
        let values = fields
            .map(|token| parse_data_value(trimmed, token))
            .collect::<Result<Vec<_>>>()?;

        rows.push(DataRow {
            raw_line: trimmed.to_string(),
            time_ps,
            values,
        });
    }

    Ok(rows)
}

fn parse_data_value(line: &str, token: &str) -> Result<f64> {
    match token {
        "inf" | "+inf" | "INF" | "+INF" => Ok(f64::INFINITY),
        "-inf" | "-INF" => Ok(f64::NEG_INFINITY),
        "nan" | "NaN" | "NAN" => Ok(f64::NAN),
        _ => token
            .parse::<f64>()
            .map_err(|_| GromacsParseError::InvalidDataValue {
                line: line.to_string(),
                value: token.to_string(),
            }),
    }
}

fn parse_temperature(line: &str) -> Result<Option<f64>> {
    let Some(idx) = line.find("T =") else {
        return Ok(None);
    };
    let value = first_float(&line[(idx + 3)..]).ok_or(GromacsParseError::InvalidField {
        field: "temperature",
        value: line.to_string(),
    })?;
    Ok(Some(value))
}

fn parse_header_lambda(line: &str) -> Result<Option<f64>> {
    let lower = line.to_ascii_lowercase();
    if !lower.contains("lambda") {
        return Ok(None);
    }
    Ok(parse_lambda_from_text(line))
}

fn parse_legend(line: &str) -> Result<Option<(usize, String)>> {
    if !line.starts_with('@') || !line.contains(" legend ") {
        return Ok(None);
    }

    let mut parts = line.splitn(4, ' ');
    let _at = parts.next();
    let series = parts
        .next()
        .ok_or_else(|| GromacsParseError::InvalidLegendLine {
            line: line.to_string(),
        })?;
    let keyword = parts
        .next()
        .ok_or_else(|| GromacsParseError::InvalidLegendLine {
            line: line.to_string(),
        })?;
    let remainder = parts
        .next()
        .ok_or_else(|| GromacsParseError::InvalidLegendLine {
            line: line.to_string(),
        })?;

    if keyword != "legend" || !series.starts_with('s') {
        return Ok(None);
    }

    let series_idx = series[1..]
        .parse::<usize>()
        .map_err(|_| GromacsParseError::InvalidLegendLine {
            line: line.to_string(),
        })?;
    let legend = extract_quoted_text(remainder).ok_or_else(|| GromacsParseError::InvalidLegendLine {
        line: line.to_string(),
    })?;
    Ok(Some((series_idx, legend)))
}

fn extract_quoted_text(text: &str) -> Option<String> {
    let start = text.find('"')?;
    let end = text.rfind('"')?;
    (end > start).then(|| text[(start + 1)..end].to_string())
}

fn is_dhdl_legend(legend: &str) -> bool {
    legend.contains("dH")
}

fn is_delta_h_legend(legend: &str) -> bool {
    legend.contains("Delta H")
        || legend.contains('Δ')
        || (legend.contains('D') && legend.contains('H') && !legend.contains("dH"))
}

fn is_potential_legend(legend: &str) -> bool {
    let lower = legend.to_ascii_lowercase();
    lower.contains("potential energy") || lower.contains("total energy")
}

fn parse_lambda_from_text(text: &str) -> Option<f64> {
    last_float(text)
}

fn float_regex() -> &'static regex::Regex {
    static FLOAT_REGEX: OnceLock<regex::Regex> = OnceLock::new();
    FLOAT_REGEX.get_or_init(|| regex::Regex::new(FP_RE).expect("valid float regex"))
}

fn first_float(text: &str) -> Option<f64> {
    float_regex().find(text)?.as_str().parse::<f64>().ok()
}

fn last_float(text: &str) -> Option<f64> {
    float_regex()
        .find_iter(text)
        .last()
        .and_then(|m| m.as_str().parse::<f64>().ok())
}

fn lambda_key(value: f64) -> i64 {
    (value * 10_000.0).round() as i64
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::{extract_dhdl, extract_u_nk, extract_u_nk_with_potential, GromacsParseError};

    #[test]
    fn parse_simple_gromacs_dhdl() {
        let content = r#"
# This file was created by GROMACS
@ subtitle "T = 300.0 lambda = 0.1000"
@ s0 legend "dH/dlambda = 0.1000"
@ s1 legend "Delta H to 0.0000"
@ s2 legend "Delta H to 0.2000"
0.0 1.0 -0.5 0.25
2.0 2.0 -0.4 0.30
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let series = extract_dhdl(file.path(), 300.0).unwrap();
        assert_eq!(series.time_ps(), &[0.0, 2.0]);
        assert_eq!(series.values(), &[1.0, 2.0]);
        assert_eq!(series.state().lambdas(), &[0.1]);
    }

    #[test]
    fn parse_simple_gromacs_u_nk() {
        let content = r#"
@ subtitle "T = 300.0 lambda = 0.1000"
@ s0 legend "dH/dlambda = 0.1000"
@ s1 legend "Delta H to 0.0000"
@ s2 legend "Delta H to 0.2000"
0.0 1.0 -0.5 0.25
2.0 2.0 -0.4 inf
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let u_nk = extract_u_nk(file.path(), 300.0).unwrap();
        let beta = 1.0 / (0.00198720425864083 * 300.0);
        assert_eq!(u_nk.n_samples(), 2);
        assert_eq!(u_nk.n_states(), 3);
        assert_eq!(u_nk.time_ps(), &[0.0, 2.0]);
        assert_eq!(u_nk.evaluated_states()[0].lambdas(), &[0.0]);
        assert_eq!(u_nk.evaluated_states()[1].lambdas(), &[0.1]);
        assert_eq!(u_nk.evaluated_states()[2].lambdas(), &[0.2]);
        assert!((u_nk.data()[0] - (-0.5 * beta)).abs() < 1e-12);
        assert_eq!(u_nk.data()[1], 0.0);
        assert!((u_nk.data()[2] - (0.25 * beta)).abs() < 1e-12);
        assert_eq!(u_nk.data()[5], f64::INFINITY);
    }

    #[test]
    fn parse_gromacs_u_nk_with_potential() {
        let content = r#"
@ subtitle "T = 300.0 lambda = 0.1000"
@ s0 legend "dH/dlambda = 0.1000"
@ s1 legend "Delta H to 0.0000"
@ s2 legend "Delta H to 0.2000"
@ s3 legend "Potential Energy"
0.0 1.0 -0.5 0.25 -10.0
2.0 2.0 -0.4 0.30 -11.0
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let (u_nk, potential) = extract_u_nk_with_potential(file.path(), 300.0).unwrap();
        assert_eq!(u_nk.n_samples(), 2);
        assert_eq!(potential, vec![-10.0, -11.0]);
    }

    #[test]
    fn gromacs_u_nk_requires_delta_h_columns() {
        let content = r#"
@ subtitle "T = 300.0 lambda = 0.1000"
@ s0 legend "dH/dlambda = 0.1000"
0.0 1.0
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let error = extract_u_nk(file.path(), 300.0).unwrap_err();
        assert_eq!(error, GromacsParseError::NoDeltaHColumns);
    }
}
