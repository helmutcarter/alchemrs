pub type CliResult<T> = Result<T, Box<dyn std::error::Error>>;

mod cli;

use clap::Parser;
use rayon::ThreadPoolBuilder;

use crate::cli::Cli;

fn main() -> CliResult<()> {
    let cli = Cli::parse();
    if let Some(threads) = cli.threads {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|error| format!("failed to initialize Rayon thread pool: {error}"))?;
    }
    crate::cli::commands::run(cli.command)
}
