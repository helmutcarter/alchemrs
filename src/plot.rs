use crate::analysis::{BlockEstimate, ConvergencePoint};
use crate::data::{DeltaFMatrix, DhdlSeries, OverlapMatrix, StatePoint};
use crate::error::{CoreError, Result};
use plotters::prelude::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConvergencePlotOptions {
    pub width: u32,
    pub height: u32,
    pub title: String,
    pub x_label: String,
    pub y_label: String,
}

impl Default for ConvergencePlotOptions {
    fn default() -> Self {
        Self {
            width: 800,
            height: 500,
            title: "Free Energy Convergence".to_string(),
            x_label: "Windows Included".to_string(),
            y_label: "Delta F (kT)".to_string(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OverlapPlotOptions {
    pub width: u32,
    pub height: u32,
    pub title: String,
}

impl Default for OverlapPlotOptions {
    fn default() -> Self {
        Self {
            width: 700,
            height: 700,
            title: "MBAR Overlap Matrix".to_string(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeltaFStatePlotOptions {
    pub width: u32,
    pub height: u32,
    pub title: String,
    pub y_label: String,
}

impl Default for DeltaFStatePlotOptions {
    fn default() -> Self {
        Self {
            width: 900,
            height: 500,
            title: "Adjacent-State Free Energies".to_string(),
            y_label: "Delta F (kT)".to_string(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TiDhdlPlotOptions {
    pub width: u32,
    pub height: u32,
    pub title: String,
    pub x_label: String,
    pub y_label: String,
}

impl Default for TiDhdlPlotOptions {
    fn default() -> Self {
        Self {
            width: 800,
            height: 500,
            title: "TI dH/dlambda".to_string(),
            x_label: "Lambda".to_string(),
            y_label: "<dH/dlambda> (kT)".to_string(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockAveragePlotOptions {
    pub width: u32,
    pub height: u32,
    pub title: String,
    pub x_label: String,
    pub y_label: String,
}

impl Default for BlockAveragePlotOptions {
    fn default() -> Self {
        Self {
            width: 800,
            height: 500,
            title: "Block Average".to_string(),
            x_label: "Block Index".to_string(),
            y_label: "Delta F (kT)".to_string(),
        }
    }
}

pub fn render_convergence_svg(
    points: &[ConvergencePoint],
    options: Option<ConvergencePlotOptions>,
) -> Result<String> {
    if points.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let options = options.unwrap_or_default();
    let mut svg = String::new();

    {
        let backend = SVGBackend::with_string(&mut svg, (options.width, options.height));
        let root = backend.into_drawing_area();
        root.fill(&WHITE)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        let (x_min, x_max) = x_bounds(points);
        let (y_min, y_max) = y_bounds(points);

        let mut chart = ChartBuilder::on(&root)
            .caption(options.title, ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .configure_mesh()
            .x_desc(options.x_label)
            .y_desc(options.y_label)
            .draw()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(LineSeries::new(
                points
                    .iter()
                    .map(|point| (point.n_windows() as i32, point.delta_f())),
                &BLUE,
            ))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(points.iter().map(|point| {
                Circle::new(
                    (point.n_windows() as i32, point.delta_f()),
                    4,
                    BLUE.filled(),
                )
            }))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        for point in points {
            if let Some(uncertainty) = point.uncertainty() {
                let x = point.n_windows() as i32;
                let y0 = point.delta_f() - uncertainty;
                let y1 = point.delta_f() + uncertainty;
                chart
                    .draw_series(std::iter::once(PathElement::new(
                        vec![(x, y0), (x, y1)],
                        BLACK,
                    )))
                    .map_err(|err| CoreError::InvalidState(err.to_string()))?;
            }
        }

        root.present()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;
    }

    Ok(svg)
}

pub fn render_overlap_matrix_svg(
    overlap: &OverlapMatrix,
    options: Option<OverlapPlotOptions>,
) -> Result<String> {
    if overlap.n_states() == 0 {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let options = options.unwrap_or_default();
    let mut svg = String::new();

    {
        let backend = SVGBackend::with_string(&mut svg, (options.width, options.height));
        let root = backend.into_drawing_area();
        root.fill(&WHITE)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        let n_states = overlap.n_states();
        let labels = overlap
            .states()
            .iter()
            .map(format_state_label)
            .collect::<Vec<_>>();

        let mut chart = ChartBuilder::on(&root)
            .caption(options.title, ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(120)
            .y_label_area_size(120)
            .build_cartesian_2d(0..n_states as i32, 0..n_states as i32)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .configure_mesh()
            .disable_mesh()
            .x_desc("To State")
            .y_desc("From State")
            .x_labels(n_states)
            .y_labels(n_states)
            .x_label_formatter(&|value| {
                labels
                    .get(*value as usize)
                    .cloned()
                    .unwrap_or_else(String::new)
            })
            .y_label_formatter(&|value| {
                labels
                    .get(*value as usize)
                    .cloned()
                    .unwrap_or_else(String::new)
            })
            .draw()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series((0..n_states).flat_map(|row| {
                (0..n_states).map(move |col| {
                    let value = overlap.values()[row * n_states + col];
                    let color = overlap_color(value);
                    Rectangle::new(
                        [(col as i32, row as i32), (col as i32 + 1, row as i32 + 1)],
                        color.filled(),
                    )
                })
            }))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series((0..n_states).flat_map(|row| {
                (0..n_states).map(move |col| {
                    let value = overlap.values()[row * n_states + col];
                    Text::new(
                        format!("{value:.2}"),
                        (col as i32, row as i32),
                        ("sans-serif", 16).into_font().color(&BLACK),
                    )
                })
            }))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        root.present()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;
    }

    Ok(svg)
}

pub fn render_delta_f_state_svg(
    delta_f: &DeltaFMatrix,
    options: Option<DeltaFStatePlotOptions>,
) -> Result<String> {
    if delta_f.n_states() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: delta_f.n_states(),
        });
    }

    let options = options.unwrap_or_default();
    let mut svg = String::new();
    let adjacent = adjacent_state_values(delta_f)?;
    let labels = adjacent
        .iter()
        .map(|(from, to, _, _)| format!("{}→{}", format_state_label(from), format_state_label(to)))
        .collect::<Vec<_>>();
    let (y_min, y_max) = delta_f_bounds(&adjacent);

    {
        let backend = SVGBackend::with_string(&mut svg, (options.width, options.height));
        let root = backend.into_drawing_area();
        root.fill(&WHITE)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        let mut chart = ChartBuilder::on(&root)
            .caption(options.title, ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(80)
            .y_label_area_size(60)
            .build_cartesian_2d(0..adjacent.len() as i32, y_min..y_max)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .configure_mesh()
            .x_desc("Adjacent State Pair")
            .y_desc(options.y_label)
            .x_labels(adjacent.len())
            .x_label_formatter(&|value| {
                labels
                    .get(*value as usize)
                    .cloned()
                    .unwrap_or_else(String::new)
            })
            .draw()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(adjacent.iter().enumerate().map(|(idx, (_, _, value, _))| {
                let x0 = idx as i32;
                let x1 = x0 + 1;
                let style = if *value >= 0.0 {
                    RGBColor(44, 120, 115).filled()
                } else {
                    RGBColor(188, 80, 48).filled()
                };
                Rectangle::new([(x0, 0.0), (x1, *value)], style)
            }))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        for (idx, (_, _, value, uncertainty)) in adjacent.iter().enumerate() {
            if let Some(uncertainty) = uncertainty {
                let x = idx as i32;
                chart
                    .draw_series(std::iter::once(PathElement::new(
                        vec![(x, value - uncertainty), (x, value + uncertainty)],
                        BLACK,
                    )))
                    .map_err(|err| CoreError::InvalidState(err.to_string()))?;
            }
        }

        root.present()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;
    }

    Ok(svg)
}

pub fn render_ti_dhdl_svg(
    series: &[DhdlSeries],
    options: Option<TiDhdlPlotOptions>,
) -> Result<String> {
    if series.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let options = options.unwrap_or_default();
    let mut points = Vec::with_capacity(series.len());
    for item in series {
        let lambda = extract_scalar_lambda(item.state())?;
        let mean = mean(item.values())?;
        let sem = sem(item.values())?;
        points.push((lambda, mean, sem));
    }
    points.sort_by(|a, b| a.0.total_cmp(&b.0));

    let (x_min, x_max) = ti_x_bounds(&points);
    let (y_min, y_max) = ti_y_bounds(&points);
    let mut svg = String::new();

    {
        let backend = SVGBackend::with_string(&mut svg, (options.width, options.height));
        let root = backend.into_drawing_area();
        root.fill(&WHITE)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        let mut chart = ChartBuilder::on(&root)
            .caption(options.title, ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .configure_mesh()
            .x_desc(options.x_label)
            .y_desc(options.y_label)
            .draw()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(AreaSeries::new(
                points.iter().map(|(lambda, value, _)| (*lambda, *value)),
                0.0,
                RGBColor(111, 168, 220).mix(0.25),
            ))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(LineSeries::new(
                points.iter().map(|(lambda, value, _)| (*lambda, *value)),
                &RGBColor(44, 120, 115),
            ))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(points.iter().map(|(lambda, value, _)| {
                Circle::new((*lambda, *value), 4, RGBColor(44, 120, 115).filled())
            }))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        for (lambda, value, sem) in &points {
            chart
                .draw_series(std::iter::once(PathElement::new(
                    vec![(*lambda, value - sem), (*lambda, value + sem)],
                    BLACK,
                )))
                .map_err(|err| CoreError::InvalidState(err.to_string()))?;
        }

        root.present()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;
    }

    Ok(svg)
}

pub fn render_block_average_svg(
    blocks: &[BlockEstimate],
    options: Option<BlockAveragePlotOptions>,
) -> Result<String> {
    if blocks.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }

    let options = options.unwrap_or_default();
    let x_min = 0;
    let x_max = blocks.len() as i32;
    let (y_min, y_max) = block_average_bounds(blocks);
    let mut svg = String::new();

    {
        let backend = SVGBackend::with_string(&mut svg, (options.width, options.height));
        let root = backend.into_drawing_area();
        root.fill(&WHITE)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        let mut chart = ChartBuilder::on(&root)
            .caption(options.title, ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .configure_mesh()
            .x_desc(options.x_label)
            .y_desc(options.y_label)
            .x_labels(blocks.len())
            .x_label_formatter(&|value| format!("{}", value + 1))
            .draw()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(LineSeries::new(
                blocks
                    .iter()
                    .map(|block| (block.block_index() as i32 + 1, block.delta_f())),
                &RGBColor(117, 112, 179),
            ))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        chart
            .draw_series(blocks.iter().map(|block| {
                Circle::new(
                    (block.block_index() as i32 + 1, block.delta_f()),
                    4,
                    RGBColor(117, 112, 179).filled(),
                )
            }))
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;

        for block in blocks {
            if let Some(uncertainty) = block.uncertainty() {
                let x = block.block_index() as i32 + 1;
                chart
                    .draw_series(std::iter::once(PathElement::new(
                        vec![
                            (x, block.delta_f() - uncertainty),
                            (x, block.delta_f() + uncertainty),
                        ],
                        BLACK,
                    )))
                    .map_err(|err| CoreError::InvalidState(err.to_string()))?;
            }
        }

        root.present()
            .map_err(|err| CoreError::InvalidState(err.to_string()))?;
    }

    Ok(svg)
}

fn x_bounds(points: &[ConvergencePoint]) -> (i32, i32) {
    let min = points.iter().map(|point| point.n_windows()).min().unwrap() as i32;
    let max = points.iter().map(|point| point.n_windows()).max().unwrap() as i32;
    if min == max {
        (min - 1, max + 1)
    } else {
        (min, max)
    }
}

fn y_bounds(points: &[ConvergencePoint]) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for point in points {
        let delta_f = point.delta_f();
        min = min.min(delta_f);
        max = max.max(delta_f);
        if let Some(uncertainty) = point.uncertainty() {
            min = min.min(delta_f - uncertainty);
            max = max.max(delta_f + uncertainty);
        }
    }
    if min == max {
        (min - 1.0, max + 1.0)
    } else {
        let padding = (max - min) * 0.1;
        (min - padding, max + padding)
    }
}

fn overlap_color(value: f64) -> RGBColor {
    let clamped = value.clamp(0.0, 1.0);
    let start = (255.0, 255.0, 255.0);
    let end = (38.0, 93.0, 171.0);
    RGBColor(
        lerp_channel(start.0, end.0, clamped),
        lerp_channel(start.1, end.1, clamped),
        lerp_channel(start.2, end.2, clamped),
    )
}

type AdjacentStateValue<'a> = (&'a StatePoint, &'a StatePoint, f64, Option<f64>);

fn adjacent_state_values(delta_f: &DeltaFMatrix) -> Result<Vec<AdjacentStateValue<'_>>> {
    let n_states = delta_f.n_states();
    let uncertainties = delta_f.uncertainties();
    let mut values = Vec::with_capacity(n_states - 1);
    for idx in 0..(n_states - 1) {
        let matrix_idx = idx * n_states + idx + 1;
        values.push((
            &delta_f.states()[idx],
            &delta_f.states()[idx + 1],
            delta_f.values()[matrix_idx],
            uncertainties
                .map(|values| values[matrix_idx])
                .filter(|value| !value.is_nan()),
        ));
    }
    Ok(values)
}

fn delta_f_bounds(adjacent: &[AdjacentStateValue<'_>]) -> (f64, f64) {
    let mut min: f64 = 0.0;
    let mut max: f64 = 0.0;
    for (_, _, value, uncertainty) in adjacent {
        min = min.min(*value);
        max = max.max(*value);
        if let Some(uncertainty) = uncertainty {
            min = min.min(*value - uncertainty);
            max = max.max(*value + uncertainty);
        }
    }
    if min == max {
        (min - 1.0, max + 1.0)
    } else {
        let padding = (max - min) * 0.1;
        (min - padding, max + padding)
    }
}

fn extract_scalar_lambda(state: &StatePoint) -> Result<f64> {
    if state.lambdas().len() != 1 {
        return Err(CoreError::Unsupported(
            "TI plotting requires one-dimensional lambda states".to_string(),
        ));
    }
    Ok(state.lambdas()[0])
}

fn mean(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        return Err(CoreError::InvalidShape {
            expected: 1,
            found: 0,
        });
    }
    Ok(values.iter().sum::<f64>() / values.len() as f64)
}

fn sem(values: &[f64]) -> Result<f64> {
    if values.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: values.len(),
        });
    }
    let average = mean(values)?;
    let variance = values
        .iter()
        .map(|value| {
            let diff = value - average;
            diff * diff
        })
        .sum::<f64>()
        / ((values.len() - 1) as f64);
    Ok((variance / values.len() as f64).sqrt())
}

fn ti_x_bounds(points: &[(f64, f64, f64)]) -> (f64, f64) {
    let min = points
        .iter()
        .map(|(lambda, _, _)| *lambda)
        .fold(f64::INFINITY, f64::min);
    let max = points
        .iter()
        .map(|(lambda, _, _)| *lambda)
        .fold(f64::NEG_INFINITY, f64::max);
    if min == max {
        (min - 0.1, max + 0.1)
    } else {
        (min, max)
    }
}

fn ti_y_bounds(points: &[(f64, f64, f64)]) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for (_, value, sem) in points {
        min = min.min(*value - *sem);
        max = max.max(*value + *sem);
    }
    if min == max {
        (min - 1.0, max + 1.0)
    } else {
        let padding = (max - min) * 0.1;
        (min - padding, max + padding)
    }
}

fn block_average_bounds(blocks: &[BlockEstimate]) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for block in blocks {
        min = min.min(block.delta_f());
        max = max.max(block.delta_f());
        if let Some(uncertainty) = block.uncertainty() {
            min = min.min(block.delta_f() - uncertainty);
            max = max.max(block.delta_f() + uncertainty);
        }
    }
    if min == max {
        (min - 1.0, max + 1.0)
    } else {
        let padding = (max - min) * 0.1;
        (min - padding, max + padding)
    }
}

fn lerp_channel(start: f64, end: f64, t: f64) -> u8 {
    (start + (end - start) * t).round() as u8
}

fn format_state_label(state: &StatePoint) -> String {
    let lambdas = state
        .lambdas()
        .iter()
        .map(|value| format!("{value:.3}"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{lambdas}]")
}

#[cfg(test)]
mod tests {
    use super::{
        render_block_average_svg, render_convergence_svg, render_delta_f_state_svg,
        render_overlap_matrix_svg, render_ti_dhdl_svg, BlockAveragePlotOptions,
        ConvergencePlotOptions, DeltaFStatePlotOptions, OverlapPlotOptions, TiDhdlPlotOptions,
    };
    use crate::analysis::{BlockEstimate, ConvergencePoint};
    use crate::data::{DeltaFMatrix, DhdlSeries, OverlapMatrix, StatePoint};
    use crate::error::CoreError;

    fn point(n_windows: usize, delta_f: f64, uncertainty: Option<f64>) -> ConvergencePoint {
        let from = StatePoint::new(vec![0.0], 300.0).unwrap();
        let to = StatePoint::new(vec![1.0], 300.0).unwrap();
        ConvergencePoint::new(n_windows, delta_f, uncertainty, from, to, None).unwrap()
    }

    #[test]
    fn render_convergence_svg_rejects_empty_input() {
        let err = render_convergence_svg(&[], None).unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 1,
                found: 0
            }
        ));
    }

    #[test]
    fn render_convergence_svg_returns_svg_document() {
        let points = vec![point(1, 0.0, Some(0.1)), point(2, 1.0, Some(0.2))];
        let svg = render_convergence_svg(
            &points,
            Some(ConvergencePlotOptions {
                title: "MBAR Convergence".to_string(),
                ..ConvergencePlotOptions::default()
            }),
        )
        .unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("MBAR Convergence"));
    }

    #[test]
    fn render_overlap_matrix_svg_rejects_empty_matrix() {
        let overlap = OverlapMatrix::new(vec![], 0, vec![]).unwrap();
        let err = render_overlap_matrix_svg(&overlap, None).unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 1,
                found: 0
            }
        ));
    }

    #[test]
    fn render_overlap_matrix_svg_returns_svg_document() {
        let state0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let overlap =
            OverlapMatrix::new(vec![1.0, 0.25, 0.25, 1.0], 2, vec![state0, state1]).unwrap();
        let svg = render_overlap_matrix_svg(
            &overlap,
            Some(OverlapPlotOptions {
                title: "Overlap Heatmap".to_string(),
                ..OverlapPlotOptions::default()
            }),
        )
        .unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Overlap Heatmap"));
        assert!(svg.contains("[0.000]"));
        assert!(svg.contains("[1.000]"));
        assert!(svg.contains("1.00"));
        assert!(svg.contains("0.25"));
    }

    #[test]
    fn render_delta_f_state_svg_rejects_single_state() {
        let state = StatePoint::new(vec![0.0], 300.0).unwrap();
        let matrix = DeltaFMatrix::new(vec![0.0], None, 1, vec![state]).unwrap();
        let err = render_delta_f_state_svg(&matrix, None).unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 2,
                found: 1
            }
        ));
    }

    #[test]
    fn render_delta_f_state_svg_returns_svg_document() {
        let state0 = StatePoint::new(vec![0.0], 300.0).unwrap();
        let state1 = StatePoint::new(vec![0.5], 300.0).unwrap();
        let state2 = StatePoint::new(vec![1.0], 300.0).unwrap();
        let matrix = DeltaFMatrix::new(
            vec![0.0, 0.8, 1.5, -0.8, 0.0, 0.7, -1.5, -0.7, 0.0],
            Some(vec![0.0, 0.1, f64::NAN, 0.1, 0.0, 0.2, f64::NAN, 0.2, 0.0]),
            3,
            vec![state0, state1, state2],
        )
        .unwrap();
        let svg = render_delta_f_state_svg(
            &matrix,
            Some(DeltaFStatePlotOptions {
                title: "Adjacent Delta F".to_string(),
                ..DeltaFStatePlotOptions::default()
            }),
        )
        .unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Adjacent Delta F"));
        assert!(svg.contains("[0.000]→[0.500]"));
        assert!(svg.contains("[0.500]→[1.000]"));
    }

    #[test]
    fn render_ti_dhdl_svg_rejects_empty_input() {
        let err = render_ti_dhdl_svg(&[], None).unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 1,
                found: 0
            }
        ));
    }

    #[test]
    fn render_ti_dhdl_svg_rejects_multidimensional_states() {
        let state = StatePoint::new(vec![0.0, 0.5], 300.0).unwrap();
        let series = vec![DhdlSeries::new(state, vec![0.0, 1.0], vec![1.0, 1.2]).unwrap()];
        let err = render_ti_dhdl_svg(&series, None).unwrap_err();
        assert!(matches!(
            err,
            CoreError::Unsupported(message) if message == "TI plotting requires one-dimensional lambda states"
        ));
    }

    #[test]
    fn render_ti_dhdl_svg_returns_svg_document() {
        let series = vec![
            DhdlSeries::new(
                StatePoint::new(vec![0.0], 300.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![3.0, 3.2, 2.8],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![0.5], 300.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![1.8, 2.0, 1.9],
            )
            .unwrap(),
            DhdlSeries::new(
                StatePoint::new(vec![1.0], 300.0).unwrap(),
                vec![0.0, 1.0, 2.0],
                vec![0.2, 0.3, 0.1],
            )
            .unwrap(),
        ];
        let svg = render_ti_dhdl_svg(
            &series,
            Some(TiDhdlPlotOptions {
                title: "TI dH/dlambda".to_string(),
                ..TiDhdlPlotOptions::default()
            }),
        )
        .unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("TI dH/dlambda"));
    }

    #[test]
    fn render_block_average_svg_rejects_empty_input() {
        let err = render_block_average_svg(&[], None).unwrap_err();
        assert!(matches!(
            err,
            CoreError::InvalidShape {
                expected: 1,
                found: 0
            }
        ));
    }

    #[test]
    fn render_block_average_svg_returns_svg_document() {
        let from = StatePoint::new(vec![0.0], 300.0).unwrap();
        let to = StatePoint::new(vec![1.0], 300.0).unwrap();
        let blocks = vec![
            BlockEstimate::new(0, 3, -2.0, Some(0.2), from.clone(), to.clone(), None).unwrap(),
            BlockEstimate::new(1, 3, -2.2, Some(0.15), from.clone(), to.clone(), None).unwrap(),
            BlockEstimate::new(2, 3, -2.1, Some(0.1), from, to, None).unwrap(),
        ];
        let svg = render_block_average_svg(
            &blocks,
            Some(BlockAveragePlotOptions {
                title: "MBAR Block Average".to_string(),
                ..BlockAveragePlotOptions::default()
            }),
        )
        .unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("MBAR Block Average"));
    }
}
