# Overview
This project aims to create a Rust-native library for alchemical free energy analysis that can be used directly by Rust crates, while later supporting Python through thin bindings.

The goal is not to optimize isolated Python hot paths. The goal is to make Rust the canonical implementation and API surface.

The library should eventually support:

- ingestion of alchemical analysis data from MD engines
- typed internal representations of u_nk, dH/dλ, state metadata, and analysis results
- estimators such as TI, BAR, and MBAR
- preprocessing operations such as equilibration trimming and decorrelation
- convergence and overlap analysis
- a Python compatibility layer built on top of the Rust core

This document defines the architecture, scope, data model, public API direction, milestones, and implementation constraints.

# Goals
Primary goals:
- First-class support for alchemical free energy methods in Rust
- All core algorithms, parsers, and data structures live in Rust.
    - Python is a client binding, not the implementation center.
- Provide a clean Rust API
    - The library should feel like an idiomatic Rust crate.
    - Internal design should not be shaped around pandas or Python objects.
- Preserve scientific correctness
    - Results should reproduce trusted reference outputs within reasonable numerical tolerance, enforced with testing.
- Support incremental adoption
    - The project should deliver useful standalone milestones
- Enable Python compatibility later
    - The Rust core should be easy to wrap using PyO3 or similar bindings.

# Non-goals

These are explicitly out of scope for the initial design:
- one-to-one reproduction of `alchemlyb`
- first-class plotting in the Rust core
- support for every MD engine immediately
- preserving pandas as the central abstraction
- implementing every convenience helper before core estimators are stable
- designing primarily for notebook ergonomics

# Product vision

The final system should support a workflow like:
1) Parse engine output into typed Rust structures
2) Validate and normalize the data
3) Trim and decorrelate time series if needed
4) Run TI, BAR, or MBAR
5) Compute uncertainty, overlap, and convergence diagnostics
6) Expose the same core functionality through python bindings and a CLI

Example future Rust usage:
~~~rust
use alchemrs::parse::amber;
use alchemrs::prep;
use alchemrs::estimators::MBAR;

fn main() -> Result<(), alchemrs::Error> {
    let u_nk = amber::extract_u_nk("prod.out", 298.15)?;
    let u_nk = prep::decorrelate_u_nk(&u_nk)?;
    let result = MBAR::new().fit(&u_nk)?;
    println!("ΔG = {} ± {}", result.delta_f, result.uncertainty);
    Ok(())
}
~~~
# Design principles
## Rust-first, not Python-first
The internal architecture should remain valid even if Python bindings never exist.

## Strong typing over implicit conventions
State metadata, dimensions, and measurement types should be encoded in types and validated structures, not inferred from DataFrame column names.

## Small stable core
Keep the initial public API narrow. Expand only after real usage reveals what should be stabilized.

## Clear separation of concerns
Parsing, data representation, preprocessing, estimation, and bindings must remain separate layers.

## Scientific transparency
Transformations to data should be explicit and inspectable. Silent coercions should be avoided.

## Reproducibility
Outputs must be deterministic when possible and accompanied by sufficient metadata.

# Proposed architecture
## Initial repository layout
~~~text
alchemrs/
  Cargo.toml
  crates/
    alchemrs-core/
    alchemrs-estimators/
    alchemrs-prep/
    alchemrs-parse/
    alchemrs-analysis/
    alchemrs-python/        # later
  examples/
  tests/
  fixtures/
  docs/
~~~
A simpler starting point is acceptable:
~~~text
alchemrs/
  crates/
    core/
    python/
~~~
but the architecture should conceptually preserve the following separation.

## Crate responsibilities
### `alchemrs-core`
Defines canonical domain types and shared errors.

- lambda/state metadata
- matrix and timeseries wrappers
- free energy result objects
- overlap and convergence result types
- validation logic
- common traits
### `alchemrs-estimators`
Implements TI, BAR, and MBAR.

- estimator algorithms
- uncertainty estimation
- estimator-specific configuration
- fit result types
### `alchemrs-prep`
Implements preprocessing operations.

- equilibration trimming
- subsampling/decorrelation
- timeseries slicing
- state-wise normalization helpers
### `alchemrs-parse`
Parses engine-specific output into core types.

- AMBER parser
- GROMACS parser
- NAMD parser
- shared parser utilities
- validation of parsed metadata
### `alchemrs-analysis`
Implements diagnostics and analysis.

- overlap matrices
- convergence analysis
- forward/reverse comparisons
- statistical summaries
### `alchemrs-python`
Thin Python wrapper over the Rust core.

- conversion between Python objects and Rust types
- exception mapping
- compatibility-facing convenience functions
- This crate should contain as little scientific logic as possible.

# Core domain model
The most important design decision is the canonical representation of alchemical data.

## Requirements for the data model
The internal data model must:
- represent values and metadata explicitly
- encode dimensions and invariants clearly
- support efficient numerical computation
- avoid dependence on Python objects or DataFrames
- be usable across parsers and estimators

## Fundamental concepts
### State point
- Represents the thermodynamic or alchemical state associated with samples or columns.
Possible fields:
- lambda vector
- temperature
- pressure if needed
- engine-specific identifiers only if required

Example:
~~~rust
pub struct StatePoint {
    pub lambdas: Vec<f64>,
    pub temperature_k: f64,
}
~~~
Longer term, this may become more structured, for example separate named lambda dimensions.
### Time series
- Represents time-indexed scalar or vector observations.
Example:
~~~rust
pub struct TimeSeries<T> {
    pub time: Vec<f64>,
    pub values: Vec<T>,
}
~~~
### DhdlSeries
- Represents dH/dλ samples associated with a specific state.
Example:
~~~rust
pub struct DhdlSeries {
    pub state: StatePoint,
    pub time_ps: Vec<f64>,
    pub values: Vec<f64>,
}
~~~
### UNkMatrix
- Represents reduced potentials for samples evaluated across states.
Conceptually:
- rows = samples
- columns = evaluated states
- plus metadata identifying sampled state and evaluated states
Example:
~~~rust
pub struct UNkMatrix {
    pub n_samples: usize,
    pub n_states: usize,
    pub data: Vec<f64>, // row-major
    pub sampled_state: Option<StatePoint>,
    pub evaluated_states: Vec<StatePoint>,
}
~~~
For real implementation, backing this with ndarray may be better than raw Vec<f64>.

### Result objects
~~~rust
pub struct FreeEnergyEstimate {
    pub delta_f: f64,
    pub uncertainty: Option<f64>,
    pub from_state: StatePoint,
    pub to_state: StatePoint,
}
~~~
For multi-state estimators, results may include full pairwise matrices:
~~~rust
pub struct DeltaFMatrix {
    pub values: Vec<f64>,
    pub uncertainties: Option<Vec<f64>>,
    pub n_states: usize,
    pub states: Vec<StatePoint>,
}
~~~
### Diagnostics
~~~rust
pub struct OverlapMatrix {
    pub values: Vec<f64>,
    pub n_states: usize,
    pub states: Vec<StatePoint>,
}
~~~
# Invariants and validation
Every core type should have explicit invariants.

## StatePoint
- lambda values must be finite
- temperature must be finite and positive
## DhdlSeries
`time_ps.len() == values.len()`
- time must be monotonic nondecreasing
- values must be finite
- associated StatePoint must be valid
## UNkMatrix
`data.len() == n_samples * n_states`
`evaluated_states.len() == n_states`
- all entries finite unless missing values are deliberately supported
- sampled-state metadata consistent with parser output
Validation should occur at construction time through smart constructors rather than relying on public mutable fields.
Example:
~~~rust
impl DhdlSeries {
    pub fn new(state: StatePoint, time_ps: Vec<f64>, values: Vec<f64>) -> Result<Self, Error> {
        // validate
        Ok(Self { state, time_ps, values })
    }
}
~~~
# Public API philosophy
## Avoid mirroring Python dynamically

The Rust API should not expose string-based mode switches everywhere.
Prefer this:
~~~rust
pub enum DecorrelationMethod {
    StatisticalInefficiency,
    EveryNth(usize),
}
~~~
instead of:
~~~rust
fn decorrelate(data: &Data, method: &str) -> ...
~~~

## Prefer typed configuration structs
If an estimator has multiple options, use config structs or builders.
Example:
~~~rust
pub struct MbarOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub initialize_from_bar: bool,
}
~~~

## Keep traits minimal
Traits should exist only where there is a real abstraction benefit.
Possible trait:
~~~rust
pub trait Estimator<Input> {
    type Output;
    fn fit(&self, input: &Input) -> Result<Self::Output, Error>;
}
~~~
Do not introduce excessive generic machinery early.

# Error model

The Rust core must use typed errors.
Example:
~~~rust
#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("invalid shape: expected {expected}, found {found}")]
    InvalidShape { expected: usize, found: usize },

    #[error("invalid state metadata: {0}")]
    InvalidState(String),

    #[error("non-finite numerical value encountered")]
    NonFiniteValue,

    #[error("parser error: {0}")]
    Parse(String),

    #[error("estimator failed to converge")]
    ConvergenceFailure,

    #[error("unsupported input: {0}")]
    Unsupported(String),
}
~~~
Principles:
- no panics for user-facing invalid input
- avoid opaque catch-all errors unless wrapping external crates
- include enough context to debug scientific failures

Python bindings will later map these into Python exceptions.

# Numerical backend choices
## Array representation

Recommended approach:
- use ndarray internally for multi-dimensional numerical work
- expose stable wrapper structs in the public API
- avoid exposing raw ndarray types everywhere unless necessary

This keeps implementation flexible while preserving ergonomics.

## Floating-point precision

Use f64 by default throughout the core. Use fused-multiply add where possible.
Reasoning:
- scientific analysis requires stable precision
- alignment with existing scientific Python expectations
- avoids complexity of making everything generic over float types

## Linear algebra

If MBAR requires matrix operations or solvers, evaluate:
- ndarray
- ndarray-linalg
- custom iterative implementations where practical

The abstraction should not use a heavyweight backend unless required.

# Estimator design
## TI
Inputs
- a collection of DhdlSeries values ordered by lambda state
- integration configuration

Outputs
- overall free energy estimate
- uncertainty 
- per-window summaries if useful

Initial API sketch
~~~rust
pub struct TiOptions {
    pub method: IntegrationMethod,
}

pub enum IntegrationMethod {
    Trapezoidal,
    Simpson,
}

pub struct TiEstimator {
    pub options: TiOptions,
}

impl TiEstimator {
    pub fn fit(&self, series: &[DhdlSeries]) -> Result<FreeEnergyEstimate, Error> {
        todo!()
    }
}
~~~
## BAR
Inputs
- forward and reverse work values, or adjacent-state reduced potential differences depending on representation chosen

Outputs
- pairwise free energy estimate
- uncertainty
- convergence metadata if available

Need to decide whether BAR should operate on:
- raw work series
- adjacent u_nk slices
- both, via conversion helpers

Recommendation:
- begin with one clean input type
- add conversion helpers later

## MBAR
- MBAR is the most complex estimator and should come later.
Requirements
- support multi-state reduced potential data
- stable convergence behavior
- explicit options for iteration and tolerance
- pairwise free energy outputs
- overlap diagnostics where possible

Initial API sketch
~~~rust
pub struct MbarOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub initial_f_k: Option<Vec<f64>>,
}

pub struct MbarEstimator {
    pub options: MbarOptions,
}

impl MbarEstimator {
    pub fn fit(&self, u_nk: &UNkMatrix) -> Result<DeltaFMatrix, Error> {
        todo!()
    }
}
~~~
# Preprocessing design

Preprocessing must be explicit because it changes the data used by estimators.

## Equilibration trimming

Functions should return a new structure or a view-like wrapper.

Example:
~~~rust
pub fn remove_burnin(series: &DhdlSeries, start_index: usize) -> Result<DhdlSeries, Error>;
~~~
## Decorrelation

Possible methods:
- statistical inefficiency-based subsampling
- user-specified stride
- no decorrelation
- should match `alchemlyb`/`pymbar` to start with

Example:
~~~rust
pub fn decorrelate_dhdl(
    series: &DhdlSeries,
    method: DecorrelationMethod,
) -> Result<DhdlSeries, Error>;
~~~
For UNkMatrix, decorrelation may require preserving row/state relationships carefully.

## Auditability

Every preprocessing step should be easy to inspect. The library should not silently trim or decorrelate without clear user intent.

# Parsing design

Parsers should convert engine outputs into canonical core types. They should not return engine-specific ad hoc structures unless there is a strong reason.

## Parser goals

- robust extraction of alchemical observables
- clear error messages
- precise mapping from engine-specific fields to canonical types
- minimal leakage of engine-specific quirks into the core API

## Parser module layout
~~~text
alchemrs-parse/
  amber/
  gromacs/
  namd/
  common/
~~~
## Parser API sketch
~~~rust
pub mod amber {
    pub fn extract_dhdl(path: impl AsRef<std::path::Path>, temperature_k: f64)
        -> Result<DhdlSeries, Error>;

    pub fn extract_u_nk(path: impl AsRef<std::path::Path>, temperature_k: f64)
        -> Result<UNkMatrix, Error>;
}
~~~
## Parser strategy

Implement parser families one at a time. Do not attempt all engines in the first milestone.
Recommendation:
- start with AMBER (`.rst7` + `.parm7`) because I know it best, and it'll be most useful to me
- prioritize fixture coverage over theoretical completeness

# Analysis and diagnostics
These components are important, but they should follow the estimators rather than block them. Performance concerns should be secondary.

## Potential analysis functions:
- overlap matrix computation
- per-window contribution summary
- cumulative convergence analysis
- forward vs reverse estimate comparison
- uncertainty breakdowns
- API should remain result-centric and typed.
Example:
~~~rust
pub fn overlap_matrix(u_nk: &UNkMatrix) -> Result<OverlapMatrix, Error>;
~~~
# Python binding strategy
Python support should be added after the Rust core is coherent.
## Principles
- binding layer is thin
- no core scientific logic in binding code
- Python compatibility is an adapter, not the design center
## Responsibilities
- convert Python arrays/DataFrames into Rust core types
- call Rust estimators/parsers
- convert results back into Python-friendly objects
- map Rust errors to Python exceptions

## Likely stack
- PyO3
- maturin
The Python package should be treated as a separate product layer built on the Rust core.

# Testing strategy
Testing is central to this project. It will ensure first-class support for MD engines and trajectory formats that I have less familiarity with.

## Unit tests

Every core type and estimator should have unit tests for:
- shape validation
- invalid input rejection
- simple known-result cases
- edge cases such as empty or singular inputs

## Fixture tests
Use minimal engine outputs and reference results.

Tests should verify:
- parser correctness
- TI/BAR/MBAR numerical output
- metadata preservation
- preprocessing behavior

## Cross-validation against trusted outputs
A critical part of development will be comparing Rust outputs against established reference calculations.

Comparison should include:
- free energy estimates
- uncertainties
- pairwise matrices
- overlap diagnostics where applicable

## Property testing
Where useful, add property tests for:
- shape-preserving transformations
- state ordering behavior
- invariants under slicing/trimming

# Documentation strategy
The Rust docs must stand on their own.

## Required documentation
- crate-level overview
- type-level documentation
- invariants on constructors
- examples for core workflows
- estimator assumptions
- parser limitations

## Example programs
Include examples for:
- TI from canonical dH/dλ
- BAR from pairwise data
- MBAR from u_nk
- parsing one engine format
- preprocessing plus estimation workflow

# Versioning and API stability

This project should be conservative about public API exposure.
## Guidelines
- keep most fields private
- use constructors and accessors
- avoid exposing internal implementation types
- prefer a narrow pub surface until design stabilizes
- do not promise semver stability too early

# Milestones
## Milestone 1: Core data model
### Deliverables:
- StatePoint
- DhdlSeries
- UNkMatrix
- `result` types
- typed error enum
- validation constructors
- unit tests
### Success criteria:
- canonical core types compile cleanly
- invariants enforced
- examples demonstrate construction and inspection

## Milestone 2: TI
### Deliverables:
- TI estimator
- integration methods
- fixture tests
- documentation examples
### Success criteria:
- reproduces trusted TI values on reference data
- handles invalid or incomplete inputs robustly

## Milestone 3: BAR
### Deliverables:
- BAR estimator
- uncertainty handling
- simple diagnostics
- reference tests
### Success criteria:
- stable convergence on test cases
- clear API for pairwise free energy estimation

## Milestone 4: MBAR
### Deliverables:
- MBAR estimator
- convergence controls
- pairwise matrix outputs
- overlap-related hooks
### Success criteria:
- numerically reliable on representative test cases
- practical performance on moderate datasets
## Milestone 5: First parser family
### Deliverables:
- one complete engine parser module
- fixtures
- parse-to-estimate example
### Success criteria:
- end-to-end workflow from engine output to ΔG

## Milestone 6: Preprocessing
### Deliverables:
- equilibration trimming
- decorrelation/subsampling
- transformation auditability
### Success criteria:
- preprocessing can be composed safely before estimators

## Milestone 7: Python bindings
### Deliverables:
- thin PyO3 wrapper
- packaging setup
- basic compatibility examples
### Success criteria:
- Python can call the Rust core successfully
- Rust remains the implementation center

# Recommended implementation order
Follow this order:
1) core types
2) TI
3) BAR
4) MBAR
5) first parser
6) preprocessing
7) diagnostics
8) Python bindings

This order minimizes architectural thrash and gets a useful Rust library working early.

# Key risks
## Data model mismatch
The chosen canonical representation may not cover real parser and estimator needs cleanly.
Mitigation:
- keep initial API small
- test against real fixtures early
- revise core types before too much parser code exists

## MBAR complexity
MBAR implementation may dominate time and complexity.
Mitigation:
- do not begin with MBAR
- validate core types through TI and BAR first
- isolate numerical solver design

## Parser sprawl
Parser support expands too quickly and fragments the project.
Mitigation:
- support one engine family first
- require fixtures for every parser addition

## Python compatibility pressure
Trying to preserve Python ergonomics too early may distort the Rust core.
Mitigation:
- postpone bindings
- treat compatibility as a wrapper concern

# Open design questions
These should be resolved during implementation:
- Should lambda dimensions be named explicitly instead of stored as plain vectors?
- Should UNkMatrix expose ndarray publicly or keep it internal?
- What is the cleanest canonical input representation for BAR?
- How should uncertainty be represented for estimators with matrix outputs?
- Should time be stored as `f64` in picoseconds, or as a stronger unit-aware type later? What about `usize` in femtoseconds?
- Should missing values be supported in parsed data, or rejected at construction?
- How are infintite potential energy values handled by different engines? How should I treat them?

