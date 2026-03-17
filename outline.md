## Implementation Outline

### 1) Repository layout and workspace setup
- Create a Cargo workspace with crates aligned to the design doc:
  - `crates/alchemrs-core` (core types + errors)
  - `crates/alchemrs-estimators` (TI first)
  - `crates/alchemrs-parse` (AMBER parser first)
  - Add a top-level `examples/` folder for end-to-end workflows.
- Decide whether to start minimal (core + estimators + parse) and expand later.

### 2) Milestone 1: Core data model (alchemrs-core)
- Implement typed error enum with context (`InvalidShape`, `InvalidState`, `NonFiniteValue`, `Parse`, `Unsupported`, etc.).
- Implement `StatePoint` with validation (finite lambdas, positive temperature).
- Implement `DhdlSeries` with validation (lengths match, monotonic time, finite values).
- Implement `UNkMatrix` with validation (shape, finite data, state metadata length).
- Implement `FreeEnergyEstimate` and `DeltaFMatrix` result types.
- Add unit tests for all invariants and error cases.

### 3) Milestone 2: TI estimator (alchemrs-estimators)
- Define `IntegrationMethod` enum (Trapezoidal, Simpson).
- Implement `TiEstimator::fit(&[DhdlSeries]) -> FreeEnergyEstimate`.
- Ensure ordering by lambda state is explicit or validated.
- Add unit tests on small synthetic datasets with known results.

### 4) Milestone 3: AMBER parser (alchemrs-parse)
- Implement `amber::extract_dhdl(path, temperature_k) -> DhdlSeries`.
- Define parser utilities for reading AMBER output blocks robustly.
- Map AMBER fields to canonical `StatePoint` (lambda vector + temperature).
- Add fixture tests using local AMBER outputs (read-only fixtures).

### 5) End-to-end example
- Add an example showing: parse AMBER -> TI estimate -> print ΔG.
- Use a small fixture to ensure deterministic output.

### 6) Validation and cross-checking
- Compare TI output with alchemlyb on the same fixture.
- Document tolerances and any expected deviations.

### 7) Follow-on milestones (post-initial delivery)
- BAR estimator, then MBAR (with careful solver design).
- Preprocessing (equilibration trimming, decorrelation).
- Diagnostics (overlap, convergence).
- Python bindings after the Rust core is stable.
