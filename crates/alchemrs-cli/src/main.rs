pub type CliResult<T> = Result<T, Box<dyn std::error::Error>>;

mod cli;
mod commands;
mod input;
mod output;
mod overlap;

use clap::Parser;

use crate::cli::Cli;

fn main() -> CliResult<()> {
    let cli = Cli::parse();
    commands::run(cli.command)
}
