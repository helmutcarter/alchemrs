use alchemrs::{overlap_eigenvalues, overlap_matrix, CoreError, MbarOptions, UNkMatrix};

use crate::cli::output::OverlapSummary;
use crate::CliResult;

pub fn summarize_overlap(
    windows: &[UNkMatrix],
    options: Option<MbarOptions>,
) -> CliResult<OverlapSummary> {
    let overlap = overlap_matrix(windows, options)?;
    let eigenvalues = overlap_eigenvalues(&overlap)?;
    if eigenvalues.len() < 2 {
        return Err(CoreError::InvalidShape {
            expected: 2,
            found: eigenvalues.len(),
        }
        .into());
    }
    Ok(OverlapSummary {
        scalar: 1.0 - eigenvalues[1],
        eigenvalues,
    })
}
