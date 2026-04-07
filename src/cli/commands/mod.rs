mod advise;
mod bar;
mod exp;
mod mbar;
mod ti;

use crate::cli::input::{resolve_input_temperature, AnalysisInputOptions};
use crate::cli::Command;
use crate::CliResult;

pub fn run(command: Command) -> CliResult<()> {
    match command {
        Command::AdviseSchedule {
            inputs,
            temperature,
            output_units,
            estimator,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            u_nk_observable,
            input_kind,
            overlap_min,
            block_cv_min,
            n_blocks,
            no_midpoints,
            output_format,
            output,
            report,
        } => {
            let temperature = resolve_input_temperature(&inputs, temperature)?;
            advise::run(
                inputs,
                AnalysisInputOptions {
                    temperature,
                    decorrelate,
                    remove_burnin,
                    auto_equilibrate,
                    fast,
                    conservative,
                    nskip,
                    u_nk_observable: Some(u_nk_observable),
                },
                advise::AdviseRunOptions {
                    output_units,
                    estimator,
                    input_kind,
                    overlap_min,
                    block_cv_min,
                    n_blocks,
                    suggest_midpoints: !no_midpoints,
                    output_format,
                    output_path: output,
                    report_path: report,
                },
            )
        }
        Command::Ti {
            inputs,
            temperature,
            method,
            output_units,
            output_format,
            output,
            parallel,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            u_nk_observable,
        } => {
            let temperature = resolve_input_temperature(&inputs, temperature)?;
            if u_nk_observable.is_some() {
                return Err(
                    "--u-nk-observable is only valid for bar, mbar, iexp, and dexp; ti uses dH/dlambda for preprocessing."
                        .into(),
                );
            }
            ti::run(
                inputs,
                AnalysisInputOptions {
                    temperature,
                    decorrelate,
                    remove_burnin,
                    auto_equilibrate,
                    fast,
                    conservative,
                    nskip,
                    u_nk_observable: None,
                },
                method,
                output_units,
                output_format,
                output,
                parallel,
            )
        }
        Command::Bar {
            inputs,
            temperature,
            method,
            output_units,
            output_format,
            output,
            overlap_summary,
            parallel,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            u_nk_observable,
        } => {
            let temperature = resolve_input_temperature(&inputs, temperature)?;
            bar::run(
                inputs,
                AnalysisInputOptions {
                    temperature,
                    decorrelate,
                    remove_burnin,
                    auto_equilibrate,
                    fast,
                    conservative,
                    nskip,
                    u_nk_observable: Some(u_nk_observable),
                },
                bar::BarRunOptions {
                    method,
                    output_units,
                    output_format,
                    output_path: output,
                    overlap_summary,
                    parallel,
                },
            )
        }
        Command::Iexp {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            u_nk_observable,
            no_uncertainty,
            output_units,
            output_format,
            output,
            overlap_summary,
            parallel,
        } => {
            let temperature = resolve_input_temperature(&inputs, temperature)?;
            exp::run_forward(
                inputs,
                AnalysisInputOptions {
                    temperature,
                    decorrelate,
                    remove_burnin,
                    auto_equilibrate,
                    fast,
                    conservative,
                    nskip,
                    u_nk_observable: Some(u_nk_observable),
                },
                exp::ExpRunOptions {
                    no_uncertainty,
                    output_units,
                    output_format,
                    output_path: output,
                    overlap_summary,
                    parallel,
                },
            )
        }
        Command::Dexp {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            u_nk_observable,
            no_uncertainty,
            output_units,
            output_format,
            output,
            overlap_summary,
            parallel,
        } => {
            let temperature = resolve_input_temperature(&inputs, temperature)?;
            exp::run_reverse(
                inputs,
                AnalysisInputOptions {
                    temperature,
                    decorrelate,
                    remove_burnin,
                    auto_equilibrate,
                    fast,
                    conservative,
                    nskip,
                    u_nk_observable: Some(u_nk_observable),
                },
                exp::ExpRunOptions {
                    no_uncertainty,
                    output_units,
                    output_format,
                    output_path: output,
                    overlap_summary,
                    parallel,
                },
            )
        }
        Command::Mbar {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            u_nk_observable,
            max_iterations,
            tolerance,
            fast_mbar,
            no_uncertainty,
            output_units,
            output_format,
            output,
            overlap_summary,
            parallel,
        } => {
            let temperature = resolve_input_temperature(&inputs, temperature)?;
            mbar::run(
                inputs,
                AnalysisInputOptions {
                    temperature,
                    decorrelate,
                    remove_burnin,
                    auto_equilibrate,
                    fast,
                    conservative,
                    nskip,
                    u_nk_observable: Some(u_nk_observable),
                },
                mbar::MbarRunOptions {
                    max_iterations,
                    tolerance,
                    fast_mbar,
                    no_uncertainty,
                    output_units,
                    output_format,
                    output_path: output,
                    overlap_summary,
                    parallel,
                },
            )
        }
    }
}
