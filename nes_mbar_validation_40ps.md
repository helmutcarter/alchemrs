# NES+MBAR Validation on 40 ps AMBER Ensemble

## Dataset

Validation was performed on the forward nonequilibrium AMBER ensemble:

`/gibbs/helmut/projects/flexible_ligands_GIST/steve/solvation_free_energy/flexible/nes_experiment/40ps/nes_mbar/run_*/fwd.out`

This ensemble contains 100 switching trajectories with:

- `ntpr = 1`
- dense `MBAR Energy analysis` output
- dense nonequilibrium switching samples appropriate for NES+MBAR testing

## Validation Goals

The first validation questions were:

1. Does the current NES+MBAR work reconstruction agree with the existing NES
   parser?
2. Does simple postprocessing of the dense `ntpr = 1` data improve the
   estimator?
3. Is temporal redundancy likely to be the dominant source of the current
   estimator offset?

## Reference Values

### Plain NES anchor

Using the existing scalar forward Jarzynski estimator:

```bash
cargo run --release -- nes --temperature 300 \
  /gibbs/helmut/projects/flexible_ligands_GIST/steve/solvation_free_energy/flexible/nes_experiment/40ps/nes_mbar/run_*/fwd.out
```

the endpoint estimate is:

- `delta_f = 35.656789577222 kT`
- `uncertainty = 0.5894804935449011 kT`

### Current NES+MBAR baseline

Using the current pooled forward-only NES+MBAR estimator with no sample stride:

```bash
cargo run --release -- nes-mbar --temperature 300 \
  /gibbs/helmut/projects/flexible_ligands_GIST/steve/solvation_free_energy/flexible/nes_experiment/40ps/nes_mbar/run_*/fwd.out
```

the endpoint estimate is:

- `delta_f = 37.382581880584 kT`

So the current NES+MBAR estimate lies about:

- `+1.725792303362 kT`

above the plain NES anchor.

## Work Reconstruction Check

To test whether the offset is dominated by a work-reconstruction mismatch, the
terminal accumulated work from the new NES+MBAR parser was compared against the
terminal reduced work from the existing NES parser trajectory by trajectory.

Measured differences:

- mean absolute terminal work difference: `0.009606671321 kT`
- maximum absolute terminal work difference: `0.013404917680 kT`

This is small compared with the `~1.7 kT` estimator offset. Therefore, the
current discrepancy is **not** primarily explained by a large inconsistency in
the accumulated work reconstruction between the two parsers.

## Postprocessing Sweep

The following variants were compared using the validation example:

```bash
cargo run --release --example nes_mbar_validation -- --temperature 300 \
  /gibbs/helmut/projects/flexible_ligands_GIST/steve/solvation_free_energy/flexible/nes_experiment/40ps/nes_mbar/run_*/fwd.out
```

### Stride variants

| Variant | Parameter | `delta_f` (kT) | `delta_f - NES` (kT) | Samples kept |
|---|---:|---:|---:|---:|
| stride | 1 | 37.382581880584 | 1.725792303362 | 2,000,000 |
| stride | 2 | 37.383263472018 | 1.726473894796 | 1,000,000 |
| stride | 5 | 37.384846102356 | 1.728056525134 | 400,000 |
| stride | 10 | 37.371260841767 | 1.714471264545 | 200,000 |
| stride | 20 | 37.279794563961 | 1.623004986739 | 100,000 |
| stride | 50 | 37.331361713391 | 1.674572136169 | 40,000 |
| stride | 100 | 37.459650653085 | 1.802861075863 | 20,000 |

### Block-end representatives

Each block-end variant keeps the last sample in each consecutive block.

| Variant | Parameter | `delta_f` (kT) | `delta_f - NES` (kT) | Samples kept |
|---|---:|---:|---:|---:|
| block_end | 5 | 37.377812532309 | 1.721022955087 | 400,000 |
| block_end | 10 | 37.358695909710 | 1.701906332488 | 200,000 |
| block_end | 20 | 37.260496375768 | 1.603706798546 | 100,000 |
| block_end | 50 | 37.316502497164 | 1.659712919942 | 40,000 |

### Block-averaged representatives

Each block-average variant replaces consecutive samples by a representative
sample with:

- block-averaged `lambda_protocol`
- block-averaged `u_protocol`
- block-averaged `u_k`
- terminal work from the end of the block

| Variant | Parameter | `delta_f` (kT) | `delta_f - NES` (kT) | Samples kept |
|---|---:|---:|---:|---:|
| block_avg | 5 | 37.391931159907 | 1.735141582685 | 400,000 |
| block_avg | 10 | 37.429290756927 | 1.772501179705 | 200,000 |
| block_avg | 20 | 37.484568761779 | 1.827779184557 | 100,000 |
| block_avg | 50 | 37.557221701083 | 1.900432123861 | 40,000 |

## Interpretation

### 1. Simple stride does not solve the offset

The stride sweep changes the estimate only modestly:

- best tested stride result: `37.279794563961 kT` at `stride = 20`
- improvement relative to stride 1: only about `0.10 kT`

This is much smaller than the full `~1.7 kT` offset from the plain NES anchor.

### 2. Block-end selection helps slightly

Among the tested postprocessing variants, `block_end` performs better than
`block_avg`:

- best tested block-end result: `37.260496375768 kT` at block size `20`
- this is the closest tested NES+MBAR variant to the NES anchor

This suggests that reducing temporal redundancy can help slightly, but the
effect is modest.

### 3. Block averaging does not improve accuracy here

All tested block-average variants moved the estimate *farther* from the NES
anchor than the corresponding block-end variants. On this dataset, naive
block-averaging of the energies within the switching window does not appear to
be the right postprocessing choice.

### 4. Redundancy is not the dominant problem

Because:

- work reconstruction already matches the plain NES parser closely, and
- stride/block processing only changes the estimate by at most `~0.1-0.12 kT`

the dominant source of the remaining offset is unlikely to be simple temporal
redundancy in the `ntpr = 1` data.

## Main Conclusion

The first 40 ps validation indicates that:

1. the parser and accumulated-work handling are internally consistent to within
   `~0.01 kT` at the per-trajectory terminal-work level,
2. simple stride or naive block averaging is **not** sufficient to remove the
   remaining estimator offset,
3. the next methodological improvement should target the **estimator weights**,
   not just denser or coarser postprocessing.

## Rejected Weighting Experiment

An experimental inverse-variance slice combiner was also implemented and tested
on the same 40 ps ensemble. In practice it strongly over-weighted apparently
low-variance slices and drove the endpoint estimate to an unphysical value of
about `135.8 kT`.

That experiment was therefore **rejected as a default estimator change**. The
result is useful because it shows that simply weighting slices by an
independent-slice inverse variance heuristic is not adequate here; the slice
statistics are too heterogeneous and too strongly coupled to support that
approximation cleanly.

## Recommended Next Step

The most important next accuracy improvement is to replace uniform pooling over
all retained time slices with a more principled multiple-time-slice weighting
scheme.

In practical terms, the next steps should be:

1. keep `sample_stride` support for research sweeps,
2. prefer block-end style representatives over naive block averaging if a
   coarse variant is needed,
3. implement time-slice weighting beyond uniform pooling,
4. add diagnostics showing which slices and trajectories dominate the final
   NES+MBAR weights.

The current results suggest that **estimator structure**, not merely
sample-thinning strategy, is now the limiting factor.

## Later Correction: Block-Based Parser/Estimator

After the original 40 ps validation above, the `nes-mbar` implementation was
corrected to use **lambda-jump blocks** rather than raw per-step saved samples.
This change was motivated by inspection of the AMBER output ordering:

- `DV/DL, AVERAGES OVER ... STEPS` defines a natural work increment,
- `Dynamically changing lambda` marks the actual protocol jump,
- the associated `MBAR Energy analysis` block follows that jump,
- and per-step treatment was too optimistic about the meaning of intermediate
  samples.

With the corrected block-based parser and the same default pooled estimator, the
40 ps ensemble result became:

- plain NES: `35.656789577222 kT`
- block-based `nes-mbar`: `37.364218258206 kT`

So the parser/bookkeeping correction did **not** qualitatively change the 40 ps
picture. It made the method more statistically defensible, but the estimator
still remains about:

- `+1.707428680984 kT`

above the plain NES anchor on this ensemble.

The practical consequence is important: for 40 ps, the earlier conclusion still
holds after the bookkeeping fix. The limiting issue is the estimator itself,
not just temporal redundancy or parser mechanics.
