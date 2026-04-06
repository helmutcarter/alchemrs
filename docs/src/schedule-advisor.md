# Schedule Advisor

The lambda-schedule advisor is a diagnostic CLI workflow exposed through `advise-schedule`. It supports both the existing `u_nk` overlap-based path and a TI-native `dH/dlambda` spacing path.

Command:

- `advise-schedule`

It does not fit a final production free energy estimate. Instead, it inspects the sampled schedule and reports where the current setup looks weak.

## What it uses

In `u_nk` mode, the current MVP combines:

- adjacent overlap from the MBAR-derived overlap matrix
- adjacent BAR or MBAR edge estimates
- edge-local block-average variability

In TI mode (`--input-kind dhdl`), it combines:

- per-window `mean_dhdl` and SEM
- window-local block-average variability
- forward/reverse `dH/dlambda` stability within each window
- interval-local slope, curvature, and trapezoid contribution
- interval uncertainty propagated from neighboring window SEMs

The output is deterministic and threshold-driven in both modes.

## Typical usage

```bash
cargo run --release -- advise-schedule \
  --temperature 300 \
  --decorrelate \
  --u-nk-observable de \
  --output-format json \
  path/to/*/prod.out
```

Useful flags:

- `--input-kind <auto|u-nk|dhdl>`
- `--estimator <mbar|bar>`
- `--overlap-min <VALUE>`
- `--block-cv-min <VALUE>`
- `--n-blocks <N>`
- `--no-midpoints`
- `--report <PATH>`

TI usage:

```bash
cargo run --release -- advise-schedule \
  --temperature 300 \
  --input-kind dhdl \
  --decorrelate \
  --output-format json \
  path/to/*/prod.out
```

## Output

JSON output includes:

- `sample_counts`
- `provenance`
- `edges` in `u_nk` mode
- `windows` and `intervals` in TI mode
- `suggestions`

If `--report` is provided, the CLI also writes a standalone HTML report with:

- a configuration summary
- a top priority queue of the highest-risk edges, including their weakest components
- ranked schedule suggestions
- edge-level diagnostics with severity badges
- for TI mode, an integration-method shape gallery for each applicable TI method on the current lambda grid
- inline priority bars for quick scanning
- inline SVG lambda-axis visuals for each edge and proposal
- an in-report legend for source, target, proposal, delta-bar, and status semantics
- per-component rows for multidimensional schedules showing which lambda components are bisected, held at the source state, dominant, or fixed
- normalized delta bars so the largest component jump on each edge is immediately visible

Each edge records:

- adjacent lambda endpoints
- lambda delta
- forward and reverse overlap values
- minimum directional overlap
- optional relative overlap versus neighboring edges
- edge `delta_f`
- optional uncertainty
- optional relative uncertainty versus neighboring edges
- optional block mean / standard deviation / coefficient of variation
- dominant changing lambda components
- priority score
- severity classification

Each suggestion records:

- suggestion kind
- target edge index
- lambda endpoints
- focus components
- proposal strategy
- optional midpoint proposal
- priority score
- reason string

In TI mode, each window records:

- lambda state
- `mean_dhdl`
- optional `sem_dhdl`
- optional block mean / standard deviation / coefficient of variation
- optional forward/reverse `dH/dlambda` delta

Each TI interval records:

- lambda endpoints
- interval width
- left and right `mean_dhdl`
- trapezoid contribution
- slope and absolute slope
- optional curvature
- optional interval uncertainty
- priority score
- severity classification

## Current heuristics

The current `u_nk` advisor classifies edges as:

- `healthy`
- `monitor`
- `add_sampling`
- `add_window`

The current pass uses:

- `overlap_min < --overlap-min` to suggest a new window
- `block_cv >= --block-cv-min` to suggest more sampling
- neighbor-aware comparisons so a uniformly weak schedule is not treated the same as one isolated bad edge
- component-level labeling so multidimensional schedules can report which lambda components dominate the jump

If midpoint proposals are enabled, `insert_window` suggestions return one of:

- `midpoint`: the component-wise midpoint between the two edge endpoints
- `focused_split`: midpoint only on the dominant changing component while holding the others at the source state

The current TI advisor classifies intervals as:

- `healthy`
- `monitor`
- `add_sampling`
- `add_window`
- `add_window_and_sampling`

The TI pass uses:

- high slope / curvature z-scores to suggest more lambda resolution
- high block CV or interval uncertainty to suggest more sampling
- midpoint proposals for inserted TI windows
- one-dimensional lambda schedules only

## Limits

- This is currently a local adjacent-edge advisor, not a full workflow planner.
- The recommendation engine is intentionally simple and should be treated as a diagnostic starting point.
- Multidimensional schedules are supported through vector endpoints and preserved lambda-component labels, but the midpoint proposal is still purely geometric.
- TI mode is currently limited to one-dimensional `dH/dlambda` schedules.
- `--input-kind auto` keeps the existing `u_nk` behavior; use `--input-kind dhdl` to force TI diagnostics through `advise-schedule`.
