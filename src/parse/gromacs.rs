use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
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
    #[error("GROMACS dhdl.xvg contains {count} dH/dlambda components; scalar dH/dlambda parsing is unsupported for multidimensional lambda schedules")]
    MultipleDhdlComponents { count: usize },
    #[error("no Delta H legends found in GROMACS dhdl.xvg output")]
    NoDeltaHColumns,
    #[error("no potential-energy legend found in GROMACS dhdl.xvg output")]
    NoPotentialSamples,
    #[error("duplicate evaluated state {state} in GROMACS dhdl.xvg output")]
    DuplicateState { state: String },
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
    sampled_state: Vec<f64>,
    lambda_labels: Option<Vec<String>>,
    dhdl_columns: Vec<usize>,
    delta_h_columns: Vec<(usize, Vec<f64>)>,
    potential_column: Option<usize>,
}

pub fn extract_dhdl(path: impl AsRef<Path>, temperature_k: f64) -> Result<DhdlSeries> {
    let header = read_header(path.as_ref())?;
    validate_temperature(header.temperature_k, temperature_k)?;

    if header.dhdl_columns.is_empty() {
        return Err(GromacsParseError::MissingDhdlLegend);
    }
    if header.dhdl_columns.len() != 1 {
        return Err(GromacsParseError::MultipleDhdlComponents {
            count: header.dhdl_columns.len(),
        });
    }

    let column = header.dhdl_columns[0];
    let mut time_ps = Vec::new();
    let mut values = Vec::new();
    for row in read_data_rows(path.as_ref())? {
        time_ps.push(row.time_ps);
        values.push(
            *row.values
                .get(column)
                .ok_or_else(|| GromacsParseError::InvalidDataLine {
                    line: row.raw_line.clone(),
                })?,
        );
    }

    if values.is_empty() {
        return Err(GromacsParseError::NoSamples);
    }

    let state = StatePoint::new(header.sampled_state, temperature_k)?;
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

    let sampled_key = state_key(&header.sampled_state);
    let mut seen_targets = HashSet::new();
    for (_, state) in &header.delta_h_columns {
        let key = state_key(state);
        if key != sampled_key && !seen_targets.insert(key) {
            return Err(GromacsParseError::DuplicateState {
                state: format_state(state),
            });
        }
    }

    let mut evaluated_states = vec![header.sampled_state.clone()];
    let mut evaluated_keys = HashSet::from([sampled_key.clone()]);
    for (_, state) in &header.delta_h_columns {
        let key = state_key(state);
        if evaluated_keys.insert(key) {
            evaluated_states.push(state.clone());
        }
    }
    evaluated_states.sort_by(|left, right| compare_state_vectors(left, right));

    let mut state_to_column = HashMap::with_capacity(evaluated_states.len());
    for (idx, state) in evaluated_states.iter().enumerate() {
        state_to_column.insert(state_key(state), idx);
    }

    let sampled_idx = *state_to_column
        .get(&sampled_key)
        .expect("sampled state indexed");
    let beta = 1.0 / (K_B_KCAL_PER_MOL_K * temperature_k);

    let mut time_ps = Vec::new();
    let mut data = Vec::new();
    let mut potential = include_potential.then(Vec::new);

    for row in read_data_rows(path)? {
        let mut reduced = vec![0.0; evaluated_states.len()];
        for (series_idx, target_state) in &header.delta_h_columns {
            let raw = *row.values.get(*series_idx).ok_or_else(|| {
                GromacsParseError::InvalidDataLine {
                    line: row.raw_line.clone(),
                }
            })?;
            let target_idx = *state_to_column
                .get(&state_key(target_state))
                .expect("target state indexed");
            if target_idx != sampled_idx {
                reduced[target_idx] = beta * raw;
            }
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

    let sampled_state = Some(StatePoint::new(header.sampled_state, temperature_k)?);
    let lambda_labels = header.lambda_labels.clone();
    let evaluated_states = evaluated_states
        .into_iter()
        .map(|state| StatePoint::new(state, temperature_k))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(GromacsParseError::from)?;

    let u_nk = UNkMatrix::new_with_labels(
        time_ps.len(),
        evaluated_states.len(),
        data,
        time_ps,
        sampled_state,
        evaluated_states,
        lambda_labels,
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
    let mut sampled_state = None;
    let mut lambda_labels = None;
    let mut dhdl_columns = Vec::new();
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
        if sampled_state.is_none() && trimmed.starts_with("@") && trimmed.contains("subtitle") {
            sampled_state = parse_state_vector(trimmed);
            lambda_labels = parse_state_labels(trimmed);
        }

        let Some((series_idx, legend)) = parse_legend(trimmed)? else {
            continue;
        };

        if is_dhdl_legend(&legend) {
            dhdl_columns.push(series_idx);
            continue;
        }

        if is_delta_h_legend(&legend) {
            let state = parse_state_vector(&legend).ok_or(GromacsParseError::InvalidField {
                field: "foreign lambda",
                value: legend.clone(),
            })?;
            delta_h_columns.push((series_idx, state));
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
        sampled_state: sampled_state.ok_or(GromacsParseError::MissingField { field: "lambda" })?,
        lambda_labels,
        dhdl_columns,
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

fn parse_legend(line: &str) -> Result<Option<(usize, String)>> {
    if !line.starts_with('@') || !line.contains(" legend ") {
        return Ok(None);
    }

    let mut parts = line.split_whitespace();
    let marker = parts.next();
    let series = parts.next();
    let keyword = parts.next();
    if marker != Some("@") || keyword != Some("legend") {
        return Ok(None);
    }
    let Some(series) = series else {
        return Err(GromacsParseError::InvalidLegendLine {
            line: line.to_string(),
        });
    };
    if !series.starts_with('s') {
        return Ok(None);
    }

    let series_idx = series[1..]
        .parse::<usize>()
        .map_err(|_| GromacsParseError::InvalidLegendLine {
            line: line.to_string(),
        })?;
    let legend = extract_quoted_text(line).ok_or_else(|| GromacsParseError::InvalidLegendLine {
        line: line.to_string(),
    })?;
    Ok(Some((series_idx, legend)))
}

fn extract_quoted_text(text: &str) -> Option<String> {
    let start = text.find('"')?;
    let end = text.rfind('"')?;
    (end > start).then(|| text[(start + 1)..end].to_string())
}

fn normalize_xvg_markup(text: &str) -> String {
    text.replace(r"\xD\f{}H", "Delta H")
        .replace(r"\xD", "Delta")
        .replace(r"\xl\f{}", "lambda")
        .replace(r"\xl", "lambda")
        .replace(r"\f{}", "")
        .replace(r"\S-1\N", "")
        .replace(r"\N", "")
}

fn is_dhdl_legend(legend: &str) -> bool {
    normalize_xvg_markup(legend)
        .to_ascii_lowercase()
        .contains("dh/dlambda")
}

fn is_delta_h_legend(legend: &str) -> bool {
    normalize_xvg_markup(legend)
        .to_ascii_lowercase()
        .contains("delta h")
}

fn is_potential_legend(legend: &str) -> bool {
    let lower = normalize_xvg_markup(legend).to_ascii_lowercase();
    lower.contains("potential energy") || lower.contains("total energy")
}

fn parse_state_vector(text: &str) -> Option<Vec<f64>> {
    if let Some(values) = parse_parenthesized_vector(text) {
        return Some(values);
    }
    last_float(text).map(|value| vec![value])
}

fn parse_state_labels(text: &str) -> Option<Vec<String>> {
    let end = text.rfind(')')?;
    let before_values = &text[..end];
    let labels_end = before_values.rfind(") =")?;
    let labels_start = before_values[..labels_end].rfind('(')?;
    let labels = before_values[(labels_start + 1)..labels_end]
        .split(',')
        .map(|label| label.trim().to_string())
        .filter(|label| !label.is_empty())
        .collect::<Vec<_>>();
    (!labels.is_empty()).then_some(labels)
}

fn parse_parenthesized_vector(text: &str) -> Option<Vec<f64>> {
    let end = text.rfind(')')?;
    let start = text[..end].rfind('(')?;
    let values = float_regex()
        .find_iter(&text[(start + 1)..end])
        .filter_map(|m| m.as_str().parse::<f64>().ok())
        .collect::<Vec<_>>();
    (!values.is_empty()).then_some(values)
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

fn compare_state_vectors(left: &[f64], right: &[f64]) -> Ordering {
    for (lhs, rhs) in left.iter().zip(right.iter()) {
        let ordering = lhs.total_cmp(rhs);
        if ordering != Ordering::Equal {
            return ordering;
        }
    }
    left.len().cmp(&right.len())
}

fn state_key(values: &[f64]) -> Vec<i64> {
    values.iter().map(|value| (value * 10_000.0).round() as i64).collect()
}

fn format_state(values: &[f64]) -> String {
    let formatted = values
        .iter()
        .map(|value| format!("{value:.4}"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("({formatted})")
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
    fn parse_multidimensional_gromacs_u_nk() {
        let content = r#"
@ subtitle "T = 300.0 lambda state 9: (coul-lambda, vdw-lambda) = (0.4000, 1.0000)"
@ s0 legend "Total Energy (kJ/mol)"
@ s1 legend "dH/d\xl\f{} coul-lambda = 0.4000"
@ s2 legend "dH/d\xl\f{} vdw-lambda = 1.0000"
@ s3 legend "\xD\f{}H \xl\f{} to (0.3000, 1.0000)"
@ s4 legend "\xD\f{}H \xl\f{} to (0.4000, 1.0000)"
@ s5 legend "\xD\f{}H \xl\f{} to (0.5000, 1.0000)"
0.0 -10.0 1.0 2.0 -5.0 0.0 5.0
2.0 -11.0 1.5 2.5 -4.0 0.0 inf
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let (u_nk, potential) = extract_u_nk_with_potential(file.path(), 300.0).unwrap();
        assert_eq!(u_nk.sampled_state().unwrap().lambdas(), &[0.4, 1.0]);
        assert_eq!(
            u_nk.lambda_labels().unwrap(),
            &["coul-lambda", "vdw-lambda"]
        );
        assert_eq!(u_nk.n_states(), 3);
        assert_eq!(u_nk.evaluated_states()[0].lambdas(), &[0.3, 1.0]);
        assert_eq!(u_nk.evaluated_states()[1].lambdas(), &[0.4, 1.0]);
        assert_eq!(u_nk.evaluated_states()[2].lambdas(), &[0.5, 1.0]);
        assert_eq!(potential, vec![-10.0, -11.0]);
    }

    #[test]
    fn multidimensional_gromacs_dhdl_is_unsupported() {
        let content = r#"
@ subtitle "T = 300.0 lambda state 9: (coul-lambda, vdw-lambda) = (0.4000, 1.0000)"
@ s0 legend "dH/d\xl\f{} coul-lambda = 0.4000"
@ s1 legend "dH/d\xl\f{} vdw-lambda = 1.0000"
0.0 1.0 2.0
"#;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let error = extract_dhdl(file.path(), 300.0).unwrap_err();
        assert_eq!(error, GromacsParseError::MultipleDhdlComponents { count: 2 });
    }
}
