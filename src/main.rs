pub type CliResult<T> = Result<T, Box<dyn std::error::Error>>;

mod cli;

use clap::Parser;

use crate::cli::Cli;

fn main() -> CliResult<()> {
    let cli = Cli::parse();
    crate::cli::commands::run(cli.command)
}
