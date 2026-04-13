# NES+MBAR Methods Appendix

## Overview

The NES+MBAR method implemented in `alchemrs` is a forward-only nonequilibrium
path-reweighting estimator for discrete target lambda states. The current
research implementation now uses **lambda-jump blocks** rather than raw
per-step saved samples. It combines:

- nonequilibrium switching trajectories initiated from equilibrium at the
  initial state,
- the accumulated switching work up to retained lambda blocks, and
- reduced energies evaluated on a discrete set of target lambda states for one
  representative configuration associated with each retained block.

The method is designed for AMBER nonequilibrium switching outputs containing
block-average `DV/DL` data, explicit lambda-jump messages, and associated
`MBAR Energy analysis` blocks.

## Notation

Let:

- `m = 1, ..., M` index nonequilibrium switching trajectories
- `b = 1, ..., B_m` index retained lambda-jump blocks within trajectory `m`
- `x_{m,b}` be the representative configuration attached to block `b`
- `\lambda_{m,b}^{-}` and `\lambda_{m,b}^{+}` be the protocol lambda before and
  after that jump
- `w_{m,b} = \beta W_{m,b}` be the cumulative **reduced** work after block `b`
- `u_{\lambda_{m,b}^{+}}(x_{m,b})` be the reduced potential at the post-jump
  protocol state, when that state exists explicitly in the discrete MBAR grid
- `u_k(x_{m,b})` be the reduced potential at target discrete state `k`

All reduced potentials are defined in units of `kT`:

\[
u(x;\lambda) = \beta U(x;\lambda), \qquad \beta = (k_B T)^{-1}
\]

## Estimator

For each target state `k`, the current implementation constructs the
nonequilibrium reweighting factor

\[
g_{m,b \to k}
=
\exp\left[
-w_{m,b}
-\bigl(u_k(x_{m,b}) - u_{\lambda_{m,b}^{+}}(x_{m,b})\bigr)
\right]
\]

and estimates the free energy of state `k` relative to the initial state by
pooling all retained finite block contributions:

\[
\hat c_k
=
\frac{1}{N_k}
\sum_{m,b \in \mathrm{retained}} g_{m,b \to k}
\]

\[
\hat f_k - f_0 = -\ln \hat c_k
\]

where `N_k` is the number of retained finite-weight contributions.

Pairwise free energies are then obtained by subtraction:

\[
\Delta \hat f_{ij} = \hat f_j - \hat f_i
\]

## Work accumulation

In the present implementation, accumulated reduced work is reconstructed
incrementally from the AMBER-reported block-average `DV/DL` values as

\[
w \leftarrow w + (\beta \, \langle dV/d\lambda \rangle_{\mathrm{block}})\,\Delta\lambda_{\mathrm{block}}
\]

once per printed lambda jump.

This was introduced specifically to avoid the earlier, less defensible
per-step reconstruction.

## Protocol-state reduced energy

The protocol-state reduced energy is only retained when the post-jump protocol
lambda is exactly one of the discrete MBAR states printed by AMBER. If the
protocol state is off-grid, the current implementation leaves that quantity
undefined and the corresponding block does not contribute a finite weight.

## Sample retention and stride

Because nonequilibrium outputs may be written with `ntpr = 1`, the implementation
still supports a postprocessing stride over retained blocks:

- `sample_stride = 1`: retain every parsed block
- `sample_stride = n`: retain every `n`th parsed block

This stride is applied after parsing and before estimation.

## Uncertainty estimation

The current implementation supports trajectory bootstrap uncertainty. Whole
switching trajectories are resampled with replacement, the full estimator is
recomputed for each bootstrap replicate, and the standard deviation across
replicates is used as the statewise uncertainty estimate.

This bootstrap is performed at the trajectory level rather than the per-slice
level because slices within a trajectory are correlated.

## Assumptions

The present formulation assumes:

1. each switching trajectory is initialized from equilibrium at the initial
   state,
2. the AMBER output provides a correctly matched MBAR energy block for each
   retained lambda jump,
3. the retained post-jump protocol state is on the discrete MBAR grid whenever
   a finite block weight is used,
4. the accumulated block work reconstruction is consistent with the actual switching
   protocol,
5. trajectory-level bootstrap is an adequate first uncertainty estimate for the
   current research-stage implementation.

## Current Status

This should be regarded as a research implementation of a
better-defined block-based nonequilibrium state-reweighting estimator, not as a
final production method. The main remaining methodological upgrades are:

- a redesigned block-based estimator rather than the current simple pooled form
- optimized multiple-time-slice weighting rather than uniform pooling
- more complete uncertainty theory
- validation of protocol-energy matching and postprocessing choices such as
  stride or block averaging

## References

- Jarzynski, C. *Phys. Rev. Lett.* **78**, 2690-2693 (1997).
- Hummer, G.; Szabo, A. *PNAS* **98**, 3658-3661 (2001).
- Minh, D. D. L.; Chodera, J. D. *J. Chem. Phys.* **131**, 134110 (2009).
- Minh, D. D. L.; Chodera, J. D. *J. Chem. Phys.* **134**, 024111 (2011).
