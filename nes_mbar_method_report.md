# NES+MBAR Method Report

## Scope

This document describes the nonequilibrium-switching-plus-state-reweighting
(`nes-mbar`) method currently implemented in `alchemrs`, the statistical
mechanical theory that motivates it, the exact data representation used by the
code, the estimator currently implemented, and the main limitations that remain
before this should be treated as a production-quality method.

The current implementation is **forward-only** and is intended as a
**research-stage estimator**. It is not ordinary equilibrium MBAR applied to
switching trajectories. Instead, it is a path-reweighted estimator that uses:

- nonequilibrium work accumulated along a driven protocol, and
- reduced energies `u_k(x)` evaluated at a discrete set of target lambda states

to estimate equilibrium free energies on that discrete target grid.

## Scientific Motivation

Ordinary MBAR assumes that samples are drawn from equilibrium distributions at
fixed thermodynamic states. Nonequilibrium switching trajectories violate that
assumption: each saved configuration is drawn from a driven path ensemble, not
from the equilibrium distribution associated with the instantaneous protocol
lambda.

The NES+MBAR design here is motivated by a standard chain of statistical
mechanical results:

1. Jarzynski's equality shows that equilibrium free energy differences can be
   obtained from an ensemble of nonequilibrium work measurements:

   \[
   e^{-\beta \Delta F} = \left\langle e^{-\beta W} \right\rangle
   \]

2. Hummer-Szabo-type relations show that snapshots drawn at intermediate times
   along a driven protocol can be reweighted back to equilibrium at that
   intermediate protocol value using the nonequilibrium work accumulated up to
   that time.

3. Multiple-time-slice estimators generalize this idea and pool information from
   many protocol times, rather than using only the terminal work.

This is the conceptual basis for combining nonequilibrium switching with
discrete-state energy reweighting: each driven snapshot is first corrected for
nonequilibrium bias by a work factor, and then reweighted from the instantaneous
protocol state to one of the discrete target states using the evaluated
`u_k(x)`.

## Primary Theoretical References

- Jarzynski, C. "Nonequilibrium Equality for Free Energy Differences,"
  *Phys. Rev. Lett.* **78**, 2690-2693 (1997).
  https://www.physics.uci.edu/~tritz/BP/Jarz97.pdf

- Hummer, G.; Szabo, A. "Free energy reconstruction from nonequilibrium
  single-molecule pulling experiments," *PNAS* **98**, 3658-3661 (2001).
  https://projects.h-its.org/dbase/js/2004-06-15/FreeEnNoneq1MolPullExp_PNAS2001.pdf

- Minh, D. D. L.; Chodera, J. D. "Optimal estimators and asymptotic variances
  for nonequilibrium path-ensemble averages," *J. Chem. Phys.* **131**, 134110
  (2009). arXiv summary:
  https://arxiv.org/abs/0907.4776

- Minh, D. D. L.; Chodera, J. D. "Estimating equilibrium ensemble averages
  using multiple time slices from driven nonequilibrium processes," *J. Chem.
  Phys.* **134**, 024111 (2011). arXiv summary:
  https://arxiv.org/abs/1010.5467

## Statistical Mechanical Basis

### 1. Jarzynski equality

Jarzynski's result provides the starting point. For trajectories initialized
from equilibrium at the initial state,

\[
\Delta F = -\beta^{-1} \ln \left\langle e^{-\beta W} \right\rangle
\]

where the average is over repeated nonequilibrium switches performed under the
same protocol.

This is the scalar endpoint result currently used in the existing `nes`
estimator in `alchemrs`.

### 2. Reweighting intermediate-time configurations

Hummer and Szabo show that if a trajectory is started from equilibrium and then
driven by a time-dependent Hamiltonian, snapshots collected at time `t` can be
reweighted with a factor involving the accumulated work up to time `t`.

In the notation relevant here, let:

- `x_{m,j}` be the configuration from trajectory `m` at saved protocol slice `j`
- `\lambda_j` be the instantaneous protocol lambda at that slice
- `w_{m,j}` be the **reduced** accumulated work up to that slice:

\[
w_{m,j} = \beta W_{m,j}
\]

- `u_{\lambda_j}(x_{m,j})` be the reduced potential at the instantaneous
  protocol state
- `u_k(x_{m,j})` be the reduced potential at discrete target state `k`

Then, for a target equilibrium state `k`, the natural nonequilibrium
reweighting factor for that sample is

\[
g_{m,j \to k}
=
\exp\left[-w_{m,j} - \bigl(u_k(x_{m,j}) - u_{\lambda_j}(x_{m,j})\bigr)\right]
\]

This combines:

- a **path reweighting factor** `exp(-w)` that corrects for the driven protocol,
  and
- a **state reweighting factor** `exp(-(u_k - u_{\lambda_j}))` that transfers
  the sample from the protocol state to target state `k`.

### 3. Multiple-time-slice pooling

The multi-time-slice theory of Minh and Chodera shows that information from
many protocol times can be pooled. The most general formulation uses
time-slice-specific weights or extended bridge sampling / EBS-style optimal
weights.

The current `alchemrs` implementation uses the simplest pooled forward-only
version: it averages all retained samples with equal weight after optional
postprocessing stride.

## What Data the Current Method Uses

The implementation is designed for AMBER nonequilibrium outputs that contain:

- protocol metadata (`temp0`, `clambda`, `dynlmb`, `ntave`)
- repeated switching blocks with `DV/DL, AVERAGES OVER ... STEPS`
- explicit `Dynamically changing lambda` events
- `MBAR Energy analysis:` blocks associated with those lambda jumps

The real parser entry point is:

- [src/parse/amber.rs](/gibbs/helmut/code/rust/alchemrs/src/parse/amber.rs)
  `extract_nes_mbar_block_trajectory(...)`

The public re-export is:

- [src/lib.rs](/gibbs/helmut/code/rust/alchemrs/src/lib.rs)

### Parsed block representation

Each retained lambda block is represented by
[NesMbarBlockSample](/gibbs/helmut/code/rust/alchemrs/src/data/timeseries.rs):

- `block_index`
- `time_ps`
- `lambda_before`
- `lambda_after`
- `reduced_work_after`
- `reduced_dvdl_avg`
- optional `reduced_rms_dvdl`
- optional `reduced_energy_protocol`
- `reduced_energies_states`

The full switching run is represented by
[NesMbarBlockTrajectory](/gibbs/helmut/code/rust/alchemrs/src/data/timeseries.rs),
which stores:

- `initial_state`
- `final_state`
- `target_states`
- `blocks`

### Why the block representation was introduced

Earlier work in this repository treated every saved `MBAR Energy analysis`
section as if it defined an independent intermediate-time sample with a
continuously updated protocol lambda and accumulated work. Careful inspection of
the AMBER output ordering showed that this assumption was too optimistic:

- `DV/DL, AVERAGES OVER 10 STEPS` defines a natural work increment,
- `Dynamically changing lambda` marks the actual protocol jump,
- the associated `MBAR Energy analysis` block comes after that jump,
- and the protocol-state energy is not always unambiguously available for
  off-grid instantaneous lambda values.

Because the multi-time-slice estimator only makes sense when
`(x_j, \lambda_j, w_j, u_{\lambda_j}(x_j))` refer to the same physical point on
the driven path, the parser was changed to use one retained sample per
lambda-jump block rather than one retained sample per printed MBAR block.

### Protocol energy

The block parser only assigns `reduced_energy_protocol` when the post-jump
protocol lambda is exactly present in the discrete MBAR grid. If the protocol
lambda is off-grid, that field is left unset instead of inventing a surrogate
energy from an ambiguous `EPtot`.

This is deliberately conservative: it reduces the number of usable retained
samples, but it avoids assigning a nonequilibrium reweighting factor whose
protocol-state reference energy is not clearly defined by the file.

### Discrete target energies

The `MBAR Energy analysis:` block is parsed into the per-block vector
`u_k(x_{m,b})` over the chosen target lambda grid.

### Accumulated work

The block parser updates reduced work only at lambda-jump events:

\[
w \leftarrow w + (\beta \, \langle dV/d\lambda \rangle_{\mathrm{block}})\,\Delta \lambda_{\mathrm{block}}
\]

where `\Delta \lambda_{\mathrm{block}}` is the actual printed lambda jump from
the `Dynamically changing lambda` line.

This is the first statistically defensible path bookkeeping used by the
implementation. It still requires validation, but it is substantially more
faithful to the AMBER switching output than the earlier per-step reconstruction.

## Current Estimator in `alchemrs`

The estimator is implemented in
[src/estimators/nes_mbar.rs](/gibbs/helmut/code/rust/alchemrs/src/estimators/nes_mbar.rs)
as:

- `NesMbarOptions`
- `NesMbarEstimator`
- `NesMbarFit`

### State set

The estimator state list is:

1. the initial switching state
2. all discrete target states from the MBAR energy grid

with duplicate removal if the initial state is already present in the discrete
grid.

### Point estimator

For each target state `k`, the current implementation computes

\[
\hat c_k
=
\frac{1}{N_k}
\sum_{m,b \in \text{retained finite blocks}}
\exp\left[-w_{m,b} - \bigl(u_k(x_{m,b}) - u_{\lambda_b}(x_{m,b})\bigr)\right]
\]

and returns

\[
\hat f_k - f_0 = -\ln \hat c_k
\]

where `N_k` is the number of retained finite block contributions.

In implementation terms:

- all trajectories are pooled
- all retained protocol blocks are pooled
- the log-sum-exp trick is used for numerical stability
- the initial state is fixed at relative free energy `0`

This is still a direct forward-only pooled path-reweighting estimator. The
block-based parser correction changed the sampling unit, but it did **not**
change the basic pooled estimator structure.

## Current Validation Status

The block-based correction matters scientifically, but it does not by itself
make the method production-ready.

On the real AMBER `10/40/80 ps` forward switching ensembles studied during
development:

- the old per-step interpretation was too optimistic about the independence and
  meaning of intermediate samples,
- the corrected block-based implementation is more expensive and more
  conservative,
- but the corrected endpoint estimates are still not consistently better than
  plain NES.

Representative ensemble results after the block-based correction:

- `10 ps`: plain NES `45.4739 kT`, block-based `nes-mbar` `50.1889 kT`
- `40 ps`: plain NES `35.6568 kT`, block-based `nes-mbar` `37.3642 kT`
- `80 ps`: plain NES `27.5817 kT`, block-based `nes-mbar` `29.6091 kT`

So the bookkeeping fix did not rescue the estimator. It mainly clarified that
the remaining problem lies in the **estimator design and weighting**, not just
in parser mechanics.

## Revised Main Conclusion

The current block-based `nes-mbar` implementation should be regarded as a
better-defined research baseline, not as a validated production method.

The evidence now supports these conclusions:

1. the earlier per-step interpretation was not statistically trustworthy,
2. block-level path bookkeeping is the right basis for future work,
3. the current pooled forward-only estimator remains too crude even with the
   corrected parser,
4. further heuristic weighting alone is unlikely to be enough.

The next serious method-development step should be a redesign of the estimator
around the block objects themselves, with a more principled multiple-time-slice
or bridge-sampling weighting model rather than more ad hoc filtering.

### `sample_stride`

The Rust API and CLI both expose `sample_stride`.

This is important for the intended research program because the raw data may be
written with `ntpr=1`, which produces extremely dense and strongly correlated
snapshots. The current estimator therefore supports a postprocessing stride that
retains every `n`th sample:

- `sample_stride = 1`: use all saved snapshots
- `sample_stride = 10`: use every 10th saved snapshot
- etc.

The current implementation does **not** yet support block averaging within
switching windows, but it is designed so that such a comparison can be added.

### Matrix assembly

The estimator computes a vector of relative free energies with respect to the
initial state and then converts them into a full pairwise
[DeltaFMatrix](/gibbs/helmut/code/rust/alchemrs/src/data/results.rs)
by subtraction:

\[
\Delta f_{ij} = f_j - f_i
\]

## Current Uncertainty Treatment

### What is implemented

The current estimator supports **trajectory bootstrap** uncertainty:

- whole trajectories are resampled with replacement
- the full NES+MBAR estimate is recomputed for each bootstrap replicate
- the standard deviation across replicates is reported as the statewise
  uncertainty

This is the correct bootstrap unit for the current data model, because protocol
slices within a single switching trajectory are strongly correlated.

### What is not implemented

The following are not yet implemented:

- analytic asymptotic variance based on EBS / optimal path-average theory
- covariance-aware uncertainty propagation for the full `DeltaFMatrix`
- time-slice-aware uncertainty decomposition

The current `DeltaFMatrix` uncertainty assembly uses simple quadrature from the
statewise bootstrap standard deviations:

\[
\sigma_{ij} \approx \sqrt{\sigma_i^2 + \sigma_j^2}
\]

This is a pragmatic first step, not a fully covariance-aware matrix uncertainty
estimate.

### CLI default

The current CLI defaults to no uncertainty unless `--n-bootstrap` is provided.
That is appropriate for now because the research-stage method should not imply a
false sense of finished uncertainty theory.

## CLI Surface

The command is:

```bash
alchemrs nes-mbar --temperature 300 path/to/run_*/fwd.out
```

The implementation is in
[src/cli/commands/nes_mbar.rs](/gibbs/helmut/code/rust/alchemrs/src/cli/commands/nes_mbar.rs),
with command wiring in:

- [src/cli/mod.rs](/gibbs/helmut/code/rust/alchemrs/src/cli/mod.rs)
- [src/cli/commands/mod.rs](/gibbs/helmut/code/rust/alchemrs/src/cli/commands/mod.rs)

Current CLI options include:

- `--temperature`
- `--n-bootstrap`
- `--seed`
- `--sample-stride`
- standard output units / output format / output file options

The current CLI reports the scalar endpoint result from the first state to the
last state in the assembled free energy matrix, along with sample counts.

## Why This Is Not Ordinary MBAR

It is important to state this explicitly.

The current method is **not**:

- equilibrium MBAR on nonequilibrium samples
- or a direct replacement for the existing `mbar` estimator

Ordinary MBAR would assume that the saved configurations at protocol slice `j`
were equilibrium samples from state `\lambda_j`. They are not. The factor
`exp(-w_{m,j})` is what corrects for that nonequilibrium bias.

So the present method is best thought of as:

1. nonequilibrium path reweighting to the instantaneous protocol state
2. state reweighting from the protocol state to a discrete target grid
3. pooling over many protocol slices and trajectories

This is why the method belongs in its own estimator family instead of being
bolted onto equilibrium MBAR semantics.

## What the Current Implementation Gets Right

- It uses trajectory-level path information, not just endpoint work.
- It uses the explicit target-state reduced energies `u_k(x)` written by AMBER.
- It distinguishes the protocol energy from the discrete target-state energies.
- It avoids the invalid assumption that raw switching snapshots are equilibrium
  MBAR samples.
- It treats the switching trajectory, not the individual time slice, as the
  natural bootstrap unit.
- It exposes `sample_stride`, which is exactly the right research knob for the
  current `ntpr=1` design.

## Current Limitations and Open Methodological Questions

### 1. Uniform pooling across time slices

At present, all retained protocol slices enter with equal weight. This is a
useful baseline, but it is not obviously optimal.

The more complete theory suggests replacing uniform pooling with:

- explicit time-slice weights `\alpha_j`
- or an extended-bridge-sampling / EBS-derived optimal combination rule

This is likely the single biggest theoretical improvement still missing.

### 2. Correlation between adjacent time slices

With `ntpr=1`, adjacent samples are extremely correlated. The estimator
currently addresses this only through `sample_stride`. That is sufficient for
research comparisons, but not yet a principled decorrelation scheme.

The obvious follow-up experiments are:

- pure stride
- block-end snapshots only
- block-average representative samples
- different work-integration granularities

### 3. Work reconstruction details

The current work accumulation rule should be validated carefully against the
AMBER switching protocol semantics and, ideally, against an independently
computed reference from the raw simulation outputs.

### 4. Uncertainty theory

Bootstrap is scientifically defensible and useful, but it is not the end state.
For a mature version of this method, one would want:

- asymptotic variance expressions from the optimal path-average framework
- covariance-aware uncertainty propagation for pairwise free energies
- possibly diagnostics for effective sample size of the nonequilibrium weights

### 5. Bias and rare-event dominance

As with Jarzynski-type estimators, large dissipation and rare-event dominance
can still create high bias or slow convergence. Pooling intermediate slices may
help substantially, but does not eliminate the need for diagnostic work.

## Recommended Research Agenda

The present implementation is strong enough to support a systematic research
study. The most useful next experiments are:

1. **Stride sensitivity**
   Compare estimates and bootstrap uncertainty for `sample_stride = 1, 2, 5,
   10, 20, ...`.

2. **Block averaging**
   Compare per-step inclusion with block-averaged representatives over the
   natural AMBER averaging window.

3. **Work reconstruction validation**
   Confirm the accumulated-work reconstruction against an independent reference
   implementation.

4. **Uniform vs. time-slice-weighted pooling**
   Test whether simple nonuniform weights improve variance and stability.

5. **Theory upgrade**
   Replace the current equal-weight pooled estimator with a multiple-time-slice
   EBS-style estimator.

6. **Uncertainty upgrade**
   Add asymptotic variance and covariance estimates once the estimator form is
   stabilized.

The likely end product of this research phase is:

- a recommended production `ntpr`
- a recommended postprocessing stride or block size
- a recommended uncertainty method
- and a better-supported production `nes-mbar` default estimator

## Bottom Line

The new `alchemrs` NES+MBAR method is a scientifically motivated,
nonequilibrium path-reweighting estimator that combines:

- forward switching work,
- intermediate nonequilibrium snapshots,
- and per-snapshot reduced energies on a discrete target lambda grid.

Its current implementation is a sound **first research implementation**:

- forward-only
- pooled over protocol slices
- optionally trajectory-bootstrapped
- intentionally equipped with `sample_stride` for methodological exploration

It should not yet be described as final or optimal. The correct framing is that
`alchemrs` now has a working experimental platform for studying a genuine
nonequilibrium state-reweighting estimator, with clear upgrade paths toward a
more rigorous multiple-time-slice and variance-aware formulation.
