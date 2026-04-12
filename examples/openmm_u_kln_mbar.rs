use alchemrs::{MbarEstimator, MbarOptions, StatePoint, UNkMatrix};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // OpenMM workflows commonly accumulate reduced potentials as u_kln[k][l][n]:
    // sampled state k, evaluated state l, sample n.
    //
    // This example shows how to convert that tensor into the per-window UNkMatrix
    // objects expected by alchemrs so the existing MBAR APIs can be used directly.
    let lambdas = vec![1.0, 0.5, 0.0];
    let temperature_k = 120.0;
    let u_kln = vec![
        vec![
            vec![0.0, 0.1, 0.2, 0.1],
            vec![0.6, 0.7, 0.8, 0.7],
            vec![1.5, 1.7, 1.8, 1.6],
        ],
        vec![
            vec![0.5, 0.6, 0.7, 0.6],
            vec![0.0, 0.1, 0.1, 0.0],
            vec![0.6, 0.7, 0.8, 0.7],
        ],
        vec![
            vec![1.3, 1.4, 1.5, 1.4],
            vec![0.5, 0.6, 0.7, 0.6],
            vec![0.0, 0.1, 0.2, 0.1],
        ],
    ];

    let windows = windows_from_openmm_u_kln(&u_kln, &lambdas, temperature_k)?;
    let fit = MbarEstimator::new(MbarOptions::default()).fit(&windows)?;
    let result = fit.result_with_uncertainty()?;
    let delta_index = result.n_states() - 1;

    println!(
        "MBAR dG = {:.6} +/- {:.6} kT",
        result.values()[delta_index],
        result.uncertainties().unwrap()[delta_index]
    );
    println!("overlap scalar = {:.6}", fit.overlap_scalar()?);

    Ok(())
}

fn windows_from_openmm_u_kln(
    u_kln: &[Vec<Vec<f64>>],
    lambdas: &[f64],
    temperature_k: f64,
) -> Result<Vec<UNkMatrix>, Box<dyn std::error::Error>> {
    if u_kln.is_empty() {
        return Err("u_kln must contain at least one sampled state".into());
    }
    if u_kln.len() != lambdas.len() {
        return Err("u_kln sampled-state dimension must match lambdas".into());
    }

    let evaluated_states = lambdas
        .iter()
        .map(|&lambda| StatePoint::new(vec![lambda], temperature_k))
        .collect::<alchemrs::Result<Vec<_>>>()?;
    let n_states = evaluated_states.len();

    let mut windows = Vec::with_capacity(u_kln.len());
    for (sampled_idx, evaluated_by_state) in u_kln.iter().enumerate() {
        if evaluated_by_state.len() != n_states {
            return Err(format!(
                "u_kln[{sampled_idx}] has {} evaluated states, expected {n_states}",
                evaluated_by_state.len()
            )
            .into());
        }

        let n_samples = evaluated_by_state[0].len();
        if n_samples == 0 {
            return Err(format!("u_kln[{sampled_idx}] has no samples").into());
        }
        if evaluated_by_state
            .iter()
            .any(|series| series.len() != n_samples)
        {
            return Err(format!("u_kln[{sampled_idx}] is ragged across evaluated states").into());
        }

        let mut sample_major = Vec::with_capacity(n_samples * n_states);
        for sample_idx in 0..n_samples {
            for evaluated_idx in 0..n_states {
                sample_major.push(evaluated_by_state[evaluated_idx][sample_idx]);
            }
        }

        let time_ps = (0..n_samples).map(|idx| idx as f64).collect::<Vec<_>>();
        let sampled_state = evaluated_states[sampled_idx].clone();
        windows.push(UNkMatrix::new(
            n_samples,
            n_states,
            sample_major,
            time_ps,
            Some(sampled_state),
            evaluated_states.clone(),
        )?);
    }

    Ok(windows)
}
