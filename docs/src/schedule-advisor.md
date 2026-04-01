# Schedule Advisor

The lambda-schedule advisor is a diagnostic CLI workflow for `u_nk` inputs.

Command:

- `advise-schedule`

It does not fit a final production free energy estimate. Instead, it inspects adjacent sampled-state edges and reports where the current schedule looks weak.

## What it uses

The current MVP combines:

- adjacent overlap from the MBAR-derived overlap matrix
- adjacent BAR or MBAR edge estimates
- edge-local block-average variability

The output is deterministic and threshold-driven.

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

- `--estimator <mbar|bar>`
- `--overlap-min <VALUE>`
- `--block-cv-min <VALUE>`
- `--n-blocks <N>`
- `--no-midpoints`
- `--report <PATH>`

## Output

JSON output includes:

- `sample_counts`
- `provenance`
- `edges`
- `suggestions`

If `--report` is provided, the CLI also writes a standalone HTML report with:

- a configuration summary
- a top priority queue of the highest-risk edges, including their weakest components
- ranked schedule suggestions
- edge-level diagnostics with severity badges
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

## Current heuristics

The current advisor classifies edges as:

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

## Limits

- This is currently a local adjacent-edge advisor, not a full workflow planner.
- The recommendation engine is intentionally simple and should be treated as a diagnostic starting point.
- Multidimensional schedules are supported through vector endpoints and preserved lambda-component labels, but the midpoint proposal is still purely geometric.
