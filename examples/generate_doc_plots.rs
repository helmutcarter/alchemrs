use alchemrs::{
    BlockAveragePlotOptions, BlockEstimate, ConvergencePoint, ConvergencePlotOptions,
    DeltaFMatrix, DeltaFStatePlotOptions, DhdlSeries, OverlapMatrix, OverlapPlotOptions,
    StatePoint, TiDhdlPlotOptions, render_block_average_svg, render_convergence_svg,
    render_delta_f_state_svg, render_overlap_matrix_svg, render_ti_dhdl_svg,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("docs/src/assets/plots")?;

    let convergence_svg = render_convergence_svg(
        &convergence_points()?,
        Some(ConvergencePlotOptions {
            title: "MBAR Convergence Example".to_string(),
            ..ConvergencePlotOptions::default()
        }),
    )?;
    std::fs::write("docs/src/assets/plots/mbar-convergence.svg", convergence_svg)?;

    let overlap_svg = render_overlap_matrix_svg(
        &overlap_matrix()?,
        Some(OverlapPlotOptions {
            title: "MBAR Overlap Example".to_string(),
            ..OverlapPlotOptions::default()
        }),
    )?;
    std::fs::write("docs/src/assets/plots/mbar-overlap.svg", overlap_svg)?;

    let delta_f_svg = render_delta_f_state_svg(
        &delta_f_matrix()?,
        Some(DeltaFStatePlotOptions {
            title: "Adjacent-State Delta F Example".to_string(),
            ..DeltaFStatePlotOptions::default()
        }),
    )?;
    std::fs::write("docs/src/assets/plots/deltaf-state.svg", delta_f_svg)?;

    let ti_svg = render_ti_dhdl_svg(
        &ti_dhdl_series()?,
        Some(TiDhdlPlotOptions {
            title: "TI dH/dlambda Example".to_string(),
            ..TiDhdlPlotOptions::default()
        }),
    )?;
    std::fs::write("docs/src/assets/plots/ti-dhdl.svg", ti_svg)?;

    let block_svg = render_block_average_svg(
        &block_estimates()?,
        Some(BlockAveragePlotOptions {
            title: "MBAR Block Average Example".to_string(),
            ..BlockAveragePlotOptions::default()
        }),
    )?;
    std::fs::write("docs/src/assets/plots/block-average.svg", block_svg)?;

    Ok(())
}

fn convergence_points() -> Result<Vec<ConvergencePoint>, Box<dyn std::error::Error>> {
    let from = StatePoint::new(vec![0.0], 300.0)?;
    let to = StatePoint::new(vec![1.0], 300.0)?;
    Ok(vec![
        ConvergencePoint::new(2, -1.80, Some(0.35), from.clone(), to.clone(), None)?,
        ConvergencePoint::new(3, -2.05, Some(0.24), from.clone(), to.clone(), None)?,
        ConvergencePoint::new(4, -2.18, Some(0.17), from.clone(), to.clone(), None)?,
        ConvergencePoint::new(5, -2.25, Some(0.12), from.clone(), to.clone(), None)?,
        ConvergencePoint::new(6, -2.28, Some(0.09), from, to, None)?,
    ])
}

fn overlap_matrix() -> Result<OverlapMatrix, Box<dyn std::error::Error>> {
    let s0 = StatePoint::new(vec![0.0], 300.0)?;
    let s1 = StatePoint::new(vec![0.5], 300.0)?;
    let s2 = StatePoint::new(vec![1.0], 300.0)?;
    Ok(OverlapMatrix::new(
        vec![
            1.00, 0.32, 0.08,
            0.32, 1.00, 0.29,
            0.08, 0.29, 1.00,
        ],
        3,
        vec![s0, s1, s2],
    )?)
}

fn delta_f_matrix() -> Result<DeltaFMatrix, Box<dyn std::error::Error>> {
    let s0 = StatePoint::new(vec![0.0], 300.0)?;
    let s1 = StatePoint::new(vec![0.5], 300.0)?;
    let s2 = StatePoint::new(vec![1.0], 300.0)?;
    Ok(DeltaFMatrix::new(
        vec![
            0.0, 0.82, 1.46,
            -0.82, 0.0, 0.64,
            -1.46, -0.64, 0.0,
        ],
        Some(vec![
            0.0, 0.12, f64::NAN,
            0.12, 0.0, 0.09,
            f64::NAN, 0.09, 0.0,
        ]),
        3,
        vec![s0, s1, s2],
    )?)
}

fn ti_dhdl_series() -> Result<Vec<DhdlSeries>, Box<dyn std::error::Error>> {
    Ok(vec![
        DhdlSeries::new(
            StatePoint::new(vec![0.0], 300.0)?,
            vec![0.0, 1.0, 2.0, 3.0],
            vec![3.0, 3.1, 2.9, 3.0],
        )?,
        DhdlSeries::new(
            StatePoint::new(vec![0.5], 300.0)?,
            vec![0.0, 1.0, 2.0, 3.0],
            vec![1.8, 2.0, 1.9, 1.7],
        )?,
        DhdlSeries::new(
            StatePoint::new(vec![1.0], 300.0)?,
            vec![0.0, 1.0, 2.0, 3.0],
            vec![0.3, 0.2, 0.1, 0.2],
        )?,
    ])
}

fn block_estimates() -> Result<Vec<BlockEstimate>, Box<dyn std::error::Error>> {
    let from = StatePoint::new(vec![0.0], 300.0)?;
    let to = StatePoint::new(vec![1.0], 300.0)?;
    Ok(vec![
        BlockEstimate::new(0, 4, -2.05, Some(0.18), from.clone(), to.clone(), None)?,
        BlockEstimate::new(1, 4, -2.12, Some(0.12), from.clone(), to.clone(), None)?,
        BlockEstimate::new(2, 4, -2.18, Some(0.09), from.clone(), to.clone(), None)?,
        BlockEstimate::new(3, 4, -2.09, Some(0.11), from, to, None)?,
    ])
}
