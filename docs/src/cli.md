# CLI Guide

The CLI is the primary entry point for `alchemrs`. It ships as the `alchemrs` binary in the same package as the Rust library crate.

Top-level commands:

- `advise-schedule`
- `ti`
- `bar`
- `mbar`
- `iexp`
- `dexp`

## Build

```bash
cargo build --release
```

## Top-Level Usage

```text
alchemrs <COMMAND>
```

Top-level options:

- `-h`, `--help`
  - Print top-level help.

## Common CLI Behavior

Most commands follow the same overall workflow:

1. load one simulation output per lambda window
2. infer temperature from input files unless `--temperature` is provided
3. optionally trim initial samples with `--remove-burnin`
4. optionally run equilibration detection with `--auto-equilibrate`
5. optionally decorrelate retained samples with `--decorrelate`
6. fit the requested estimator or advisor
7. write text, JSON, CSV, or an HTML report depending on command and flags

Input files are positional arguments:

- `INPUTS...`
  - Required for every subcommand.
  - Accepted file types are AMBER `.out` and GROMACS `dhdl.xvg`.
  - Pass one file per lambda window.

## Shared Option Semantics

These options recur across commands.

- `--temperature <TEMPERATURE>`
  - Temperature in kelvin.
  - If omitted, the CLI attempts to infer it from the input files.

- `--decorrelate`
  - Apply decorrelation after burn-in removal / equilibration detection.
  - `ti` decorrelates the `dH/dλ` series.
  - `bar`, `mbar`, `iexp`, `dexp`, and `advise-schedule` decorrelate using the selected `u_nk` observable when operating on `u_nk` data.

- `--remove-burnin <REMOVE_BURNIN>`
  - Skip this many initial samples before any other analysis.
  - Default: `0`

- `--auto-equilibrate`
  - Automatically detect equilibration and remove burn-in.

- `--fast`
  - Use the fast statistical inefficiency estimate.

- `--conservative[=<CONSERVATIVE>]`
  - Control conservative subsampling.
  - Default: `true`
  - Allowed values: `true`, `false`
  - Examples:
    - `--conservative`
    - `--conservative=true`
    - `--conservative=false`

- `--nskip <NSKIP>`
  - Stride used during equilibration detection.
  - Default: `1`

- `--output-units <OUTPUT_UNITS>`
  - Available on estimator commands and `advise-schedule`.
  - Allowed values:
    - `kt`
    - `kcal`
    - `kj`
  - Default: `kt`

- `--output-format <OUTPUT_FORMAT>`
  - Allowed values:
    - `text`
    - `json`
    - `csv`
  - Default: `text`

- `--output <OUTPUT>`
  - Write command output to a file instead of stdout.

- `--parallel`
  - Enable parallel processing.
  - Available on `ti`, `bar`, `mbar`, `iexp`, and `dexp`.

- `--u-nk-observable <U_NK_OBSERVABLE>`
  - Available on `bar`, `mbar`, `iexp`, `dexp`, and `advise-schedule`.
  - Hidden-but-recognized on `ti` only to provide a clearer error message.
  - Allowed values:
    - `epot`
    - `de`
    - `all`
  - Default: `de`

- `--overlap-summary`
  - Include overlap scalar and overlap eigenvalues in estimator output.
  - Available on `bar`, `mbar`, `iexp`, and `dexp`.

## Observable Selection

For `u_nk`-based estimators:

- `de`
  - Default.
  - Matches the adjacent-state `ΔE` style observable used by `alchemlyb`.

- `all`
  - Uses the full `u_nk` row sum.

- `epot`
  - Uses an engine-provided potential-energy observable.
  - AMBER: `EPtot`
  - GROMACS `dhdl.xvg`: `Potential Energy`

Use `epot` when:

- you want the CLI's external-observable path
- the `u_nk` matrix contains positive infinity values that make `de` unusable
- you need a preprocessing observable that does not depend on adjacent-state `ΔE` semantics

For multidimensional `u_nk` schedules:

- `bar`, `mbar`, `iexp`, and `dexp` accept them
- CLI output renders lambda states as JSON arrays or bracketed values in text / CSV output
- parser-derived lambda component labels are included in provenance when available
- `de` remains one-dimensional, so multidimensional schedules generally require `all` or `epot`

## TI-Specific Behavior

Current TI-specific behavior:

- `ti --method auto` automatically selects an integration method and records both `ti_method` and `ti_method_reason` in provenance.
- For GROMACS files with multidimensional lambda schedules, the CLI currently rejects TI input because multiple `dH/dλ` components cannot yet be collapsed into one scalar TI series safely.
- `--u-nk-observable` is intentionally not part of the public `ti --help` surface. If supplied, it is only accepted so the command can emit a domain-specific error instead of a generic unknown-flag parse error.

## `advise-schedule`

Usage:

```text
alchemrs advise-schedule [OPTIONS] <INPUTS>...
```

Purpose:

- Run schedule diagnostics instead of computing a final scalar free-energy estimate.
- Supports both `u_nk`-based schedule analysis and TI-style `dH/dλ` schedule analysis.
- In `u_nk` report mode, includes an MBAR-derived overlap-matrix view inspired by the overlap-matrix diagnostic discussed in Klimovich, Shirts, and Mobley, "Guidelines for the analysis of free energy calculations", J Comput Aided Mol Des 29, 397-411 (2015), doi:10.1007/s10822-015-9840-9.

Arguments:

- `INPUTS...`
  - Required.
  - Simulation output files, one per lambda window.

Options:

- `--temperature <TEMPERATURE>`
  - Temperature in kelvin.
  - Inferred from inputs when omitted.

- `--estimator <ESTIMATOR>`
  - Estimator used for adjacent-edge estimates in `u_nk` advisor mode.
  - Allowed values:
    - `mbar`
    - `bar`
  - Default: `mbar`

- `--output-units <OUTPUT_UNITS>`
  - Output units for energy-valued advisor diagnostics.
  - Allowed values:
    - `kt`
    - `kcal`
    - `kj`
  - Default: `kt`

- `--decorrelate`
  - Apply decorrelation to each window using the selected observable.

- `--remove-burnin <REMOVE_BURNIN>`
  - Skip this many initial samples before analysis.
  - Default: `0`

- `--auto-equilibrate`
  - Detect equilibration and remove burn-in automatically.

- `--fast`
  - Use the fast statistical inefficiency estimate.

- `--conservative[=<CONSERVATIVE>]`
  - Use conservative subsampling.
  - Default: `true`
  - Allowed values: `true`, `false`

- `--nskip <NSKIP>`
  - Equilibration-detection stride.
  - Default: `1`

- `--u-nk-observable <U_NK_OBSERVABLE>`
  - Observable used for `u_nk` auto-equilibration and decorrelation.
  - Allowed values:
    - `epot`
    - `de`
    - `all`
  - Default: `de`

- `--input-kind <INPUT_KIND>`
  - Force the advisor to treat inputs as `u_nk` or `dH/dλ` data.
  - Allowed values:
    - `auto`
    - `u-nk`
    - `dhdl`
  - Default: `auto`

- `--overlap-min <OVERLAP_MIN>`
  - Minimum adjacent overlap before suggesting a new window.
  - Default: `0.03`

- `--block-cv-min <BLOCK_CV_MIN>`
  - Minimum block coefficient of variation before suggesting more sampling.
  - Default: `0.15`

- `--n-blocks <N_BLOCKS>`
  - Number of blocks used for block averaging.
  - Default: `4`

- `--no-midpoints`
  - Disable midpoint proposals for insert-window suggestions.

- `--output-format <OUTPUT_FORMAT>`
  - Output format.
  - Allowed values:
    - `text`
    - `json`
    - `csv`
  - Default: `text`

- `--output <OUTPUT>`
  - Write advisor output to a file.

- `--report <REPORT>`
  - Write a standalone HTML advisor report to a file.

Behavior notes:

- `u_nk` mode reports overlap-driven diagnostics and window insertion / sampling suggestions.
- `dhdl` mode reports TI spacing diagnostics such as means, SEMs, block CV, split-half drift, slope, curvature, trapezoid contributions, and interval uncertainty.
- When `--report` is provided, the command writes a standalone HTML report in addition to normal structured output.

Examples:

```bash
cargo run --release -- advise-schedule \
  --temperature 300 \
  --decorrelate \
  --u-nk-observable de \
  --report schedule-report.html \
  --output-format json \
  path/to/*/prod.out
```

Force TI-style schedule diagnostics:

```bash
cargo run --release -- advise-schedule \
  --temperature 300 \
  --input-kind dhdl \
  --decorrelate \
  --report ti-schedule-report.html \
  --output-format json \
  path/to/*/prod.out
```

## `ti`

Usage:

```text
alchemrs ti [OPTIONS] <INPUTS>...
```

Purpose:

- Perform thermodynamic integration on `dH/dλ` windows.

Arguments:

- `INPUTS...`
  - Required.
  - Simulation output files, one per lambda window.

Options:

- `--temperature <TEMPERATURE>`
  - Temperature in kelvin.
  - Inferred from inputs when omitted.

- `--method <METHOD>`
  - Integration method.
  - Allowed values:
    - `auto`
    - `trapezoidal`
    - `simpson`
    - `cubic-spline`
    - `pchip`
    - `akima`
    - `gaussian-quadrature`
  - Default: `trapezoidal`

- `--output-units <OUTPUT_UNITS>`
  - Output units.
  - Allowed values:
    - `kt`
    - `kcal`
    - `kj`
  - Default: `kt`

- `--output-format <OUTPUT_FORMAT>`
  - Output format.
  - Allowed values:
    - `text`
    - `json`
    - `csv`
  - Default: `text`

- `--output <OUTPUT>`
  - Write output to a file.

- `--parallel`
  - Enable parallel processing.

- `--decorrelate`
  - Apply decorrelation to each window using the `dH/dλ` series.

- `--remove-burnin <REMOVE_BURNIN>`
  - Skip this many initial samples before analysis.
  - Default: `0`

- `--auto-equilibrate`
  - Detect equilibration and remove burn-in automatically.

- `--fast`
  - Use the fast statistical inefficiency estimate.

- `--conservative[=<CONSERVATIVE>]`
  - Use conservative subsampling.
  - Default: `true`
  - Allowed values: `true`, `false`

- `--nskip <NSKIP>`
  - Equilibration-detection stride.
  - Default: `1`

Notes:

- `ti` does not publicly expose `--u-nk-observable`.
- `--method auto` chooses a supported TI integration method after preprocessing.

Example:

```bash
cargo run --release -- ti \
  --temperature 300 \
  --method auto \
  --decorrelate \
  path/to/*/prod.out
```

## `bar`

Usage:

```text
alchemrs bar [OPTIONS] <INPUTS>...
```

Purpose:

- Compute free-energy differences with BAR.

Arguments:

- `INPUTS...`
  - Required.
  - Simulation output files, one per lambda window.

Options:

- `--temperature <TEMPERATURE>`
  - Temperature in kelvin.
  - Inferred from inputs when omitted.

- `--method <METHOD>`
  - BAR root-finding method.
  - Allowed values:
    - `false-position`
    - `self-consistent-iteration`
    - `bisection`
  - Default: `false-position`

- `--output-units <OUTPUT_UNITS>`
  - Output units.
  - Allowed values:
    - `kt`
    - `kcal`
    - `kj`
  - Default: `kt`

- `--output-format <OUTPUT_FORMAT>`
  - Output format.
  - Allowed values:
    - `text`
    - `json`
    - `csv`
  - Default: `text`

- `--output <OUTPUT>`
  - Write output to a file.

- `--overlap-summary`
  - Include overlap scalar and eigenvalues in the output.

- `--parallel`
  - Enable parallel processing.

- `--decorrelate`
  - Apply decorrelation using the selected `u_nk` observable.

- `--remove-burnin <REMOVE_BURNIN>`
  - Skip this many initial samples before analysis.
  - Default: `0`

- `--auto-equilibrate`
  - Detect equilibration and remove burn-in automatically.

- `--fast`
  - Use the fast statistical inefficiency estimate.

- `--conservative[=<CONSERVATIVE>]`
  - Use conservative subsampling.
  - Default: `true`
  - Allowed values: `true`, `false`

- `--nskip <NSKIP>`
  - Equilibration-detection stride.
  - Default: `1`

- `--u-nk-observable <U_NK_OBSERVABLE>`
  - Observable used for `u_nk` auto-equilibration and decorrelation.
  - Allowed values:
    - `epot`
    - `de`
    - `all`
  - Default: `de`

Example:

```bash
cargo run --release -- bar \
  --temperature 300 \
  --decorrelate \
  --u-nk-observable de \
  --overlap-summary \
  path/to/*/prod.out
```

## `mbar`

Usage:

```text
alchemrs mbar [OPTIONS] <INPUTS>...
```

Purpose:

- Compute free-energy differences with MBAR.

Arguments:

- `INPUTS...`
  - Required.
  - Simulation output files, one per lambda window.

Options:

- `--temperature <TEMPERATURE>`
  - Temperature in kelvin.
  - Inferred from inputs when omitted.

- `--decorrelate`
  - Apply decorrelation using the selected `u_nk` observable.

- `--remove-burnin <REMOVE_BURNIN>`
  - Skip this many initial samples before analysis.
  - Default: `0`

- `--auto-equilibrate`
  - Detect equilibration and remove burn-in automatically.

- `--fast`
  - Use the fast statistical inefficiency estimate.

- `--conservative[=<CONSERVATIVE>]`
  - Use conservative subsampling.
  - Default: `true`
  - Allowed values: `true`, `false`

- `--nskip <NSKIP>`
  - Equilibration-detection stride.
  - Default: `1`

- `--u-nk-observable <U_NK_OBSERVABLE>`
  - Observable used for `u_nk` auto-equilibration and decorrelation.
  - Allowed values:
    - `epot`
    - `de`
    - `all`
  - Default: `de`

- `--max-iterations <MAX_ITERATIONS>`
  - Maximum number of MBAR iterations.
  - Default: `10000`

- `--tolerance <TOLERANCE>`
  - Relative convergence tolerance for MBAR.
  - Default: `1.0e-7`

- `--fast-mbar`
  - Use the fast L-BFGS MBAR backend.

- `--no-uncertainty`
  - Disable uncertainty estimation.

- `--output-units <OUTPUT_UNITS>`
  - Output units.
  - Allowed values:
    - `kt`
    - `kcal`
    - `kj`
  - Default: `kt`

- `--output-format <OUTPUT_FORMAT>`
  - Output format.
  - Allowed values:
    - `text`
    - `json`
    - `csv`
  - Default: `text`

- `--output <OUTPUT>`
  - Write output to a file.

- `--overlap-summary`
  - Include overlap scalar and eigenvalues in the output.

- `--parallel`
  - Enable parallel processing.

Example:

```bash
cargo run --release -- mbar \
  --temperature 300 \
  --auto-equilibrate \
  --decorrelate \
  --u-nk-observable epot \
  --output-format json \
  path/to/*/prod.out
```

## FEP

See [Estimators](estimators.md#iexp-and-dexp) for background on naming conventions

### `iexp`

Usage:

```text
alchemrs iexp [OPTIONS] <INPUTS>...
```

Purpose:

- Compute exponential averaging estimates, typically in the direction of particle insertion.

Arguments:

- `INPUTS...`
  - Required.
  - Simulation output files, one per lambda window.

Options:

- `--temperature <TEMPERATURE>`
  - Temperature in kelvin.
  - Inferred from inputs when omitted.

- `--decorrelate`
  - Apply decorrelation using the selected `u_nk` observable.

- `--remove-burnin <REMOVE_BURNIN>`
  - Skip this many initial samples before analysis.
  - Default: `0`

- `--auto-equilibrate`
  - Detect equilibration and remove burn-in automatically.

- `--fast`
  - Use the fast statistical inefficiency estimate.

- `--conservative[=<CONSERVATIVE>]`
  - Use conservative subsampling.
  - Default: `true`
  - Allowed values: `true`, `false`

- `--nskip <NSKIP>`
  - Equilibration-detection stride.
  - Default: `1`

- `--u-nk-observable <U_NK_OBSERVABLE>`
  - Observable used for `u_nk` auto-equilibration and decorrelation.
  - Allowed values:
    - `epot`
    - `de`
    - `all`
  - Default: `de`

- `--no-uncertainty`
  - Disable uncertainty estimation.

- `--output-units <OUTPUT_UNITS>`
  - Output units.
  - Allowed values:
    - `kt`
    - `kcal`
    - `kj`
  - Default: `kt`

- `--output-format <OUTPUT_FORMAT>`
  - Output format.
  - Allowed values:
    - `text`
    - `json`
    - `csv`
  - Default: `text`

- `--output <OUTPUT>`
  - Write output to a file.

- `--overlap-summary`
  - Include overlap scalar and eigenvalues in the output.

- `--parallel`
  - Enable parallel processing.

Example:

```bash
cargo run --release -- iexp --temperature 300 path/to/*/prod.out
```

### `dexp`

Usage:

```text
alchemrs dexp [OPTIONS] <INPUTS>...
```

Purpose:

- Compute exponential averaging estimates, typically in the direction of particle deletion.

Arguments:

- `INPUTS...`
  - Required.
  - Simulation output files, one per lambda window.

Options:

- `--temperature <TEMPERATURE>`
  - Temperature in kelvin.
  - Inferred from inputs when omitted.

- `--decorrelate`
  - Apply decorrelation using the selected `u_nk` observable.

- `--remove-burnin <REMOVE_BURNIN>`
  - Skip this many initial samples before analysis.
  - Default: `0`

- `--auto-equilibrate`
  - Detect equilibration and remove burn-in automatically.

- `--fast`
  - Use the fast statistical inefficiency estimate.

- `--conservative[=<CONSERVATIVE>]`
  - Use conservative subsampling.
  - Default: `true`
  - Allowed values: `true`, `false`

- `--nskip <NSKIP>`
  - Equilibration-detection stride.
  - Default: `1`

- `--u-nk-observable <U_NK_OBSERVABLE>`
  - Observable used for `u_nk` auto-equilibration and decorrelation.
  - Allowed values:
    - `epot`
    - `de`
    - `all`
  - Default: `de`

- `--no-uncertainty`
  - Disable uncertainty estimation.

- `--output-units <OUTPUT_UNITS>`
  - Output units.
  - Allowed values:
    - `kt`
    - `kcal`
    - `kj`
  - Default: `kt`

- `--output-format <OUTPUT_FORMAT>`
  - Output format.
  - Allowed values:
    - `text`
    - `json`
    - `csv`
  - Default: `text`

- `--output <OUTPUT>`
  - Write output to a file.

- `--overlap-summary`
  - Include overlap scalar and eigenvalues in the output.

- `--parallel`
  - Enable parallel processing.

Example:

```bash
cargo run --release -- dexp --temperature 300 path/to/*/prod.out
```

## Effective Settings with Auto-Equilibration

When `--auto-equilibrate` is enabled, the CLI reports the effective preprocessing policy in provenance:

- `fast = true`
- `conservative = false`

That override is deliberate and mirrors the `alchemlyb` / `pymbar` equilibration-detection workflow.
