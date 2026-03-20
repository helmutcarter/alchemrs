mod bar;
mod exp;
mod mbar;
mod ti;

use crate::cli::Command;
use crate::input::AnalysisInputOptions;
use crate::CliResult;

pub fn run(command: Command) -> CliResult<()> {
    match command {
        Command::Ti {
            inputs,
            temperature,
            method,
            output_units,
            parallel,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
        } => ti::run(
            inputs,
            input_options(
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            ),
            method,
            output_units,
            parallel,
        ),
        Command::Bar {
            inputs,
            temperature,
            method,
            output_units,
            parallel,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
        } => bar::run(
            inputs,
            input_options(
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            ),
            method,
            output_units,
            parallel,
        ),
        Command::Exp {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            no_uncertainty,
            output_units,
            parallel,
        } => exp::run_forward(
            inputs,
            input_options(
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            ),
            no_uncertainty,
            output_units,
            parallel,
        ),
        Command::Dexp {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            no_uncertainty,
            output_units,
            parallel,
        } => exp::run_reverse(
            inputs,
            input_options(
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            ),
            no_uncertainty,
            output_units,
            parallel,
        ),
        Command::Mbar {
            inputs,
            temperature,
            decorrelate,
            remove_burnin,
            auto_equilibrate,
            fast,
            conservative,
            nskip,
            max_iterations,
            tolerance,
            no_uncertainty,
            output_units,
            parallel,
        } => mbar::run(
            inputs,
            input_options(
                temperature,
                decorrelate,
                remove_burnin,
                auto_equilibrate,
                fast,
                conservative,
                nskip,
            ),
            max_iterations,
            tolerance,
            no_uncertainty,
            output_units,
            parallel,
        ),
    }
}

fn input_options(
    temperature: f64,
    decorrelate: bool,
    remove_burnin: usize,
    auto_equilibrate: bool,
    fast: bool,
    conservative: bool,
    nskip: usize,
) -> AnalysisInputOptions {
    AnalysisInputOptions {
        temperature,
        decorrelate,
        remove_burnin,
        auto_equilibrate,
        fast,
        conservative,
        nskip,
    }
}
