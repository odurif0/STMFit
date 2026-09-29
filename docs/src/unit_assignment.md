# Unit assignment (GlcN/GlcNAc per lobe)

Chitosan is a β-(1→4) copolymer of GlcN (glucosamine, class 0) and GlcNAc
(N-acetylglucosamine, class 1). Assigning every fitted lobe a class gives an
ordered acetylation map of the chain. Counting (`N_selected`) and assignment are
separate stages; assignment works at the final molecule-consensus counts.

## Label-free rules

- No benchmark label, control sequence (`NKNNKN`), expected count or class
  proportion enters features, learning, thresholds, abstention or fusion.
  Labels are read only by `grade_*`/`report_*` scripts.
- No composition prior: rules such as "the two brightest lobes are GlcNAc" are
  not label-free methods.
- **Physical naming.** Clusters are unnamed until the cluster with the higher
  mean raw amplitude is called GlcNAc (the acetyl group is larger). This
  convention is fixed; the truth is never used to choose the 0↔1 flip.

## Per-scan assignment (`test/run_reconstructed_chitosan.jl`)

| Stage | Output | Content |
|---|---|---|
| Base fit | `features.tsv` | Gaussian geometry at the final count (`extract_lobe_features.jl`) |
| Split fit | `features_split.tsv` | Same fit with a split (asymmetric) axial profile, `skew_ratio ≤ 2`; gives `split_log_skew` |
| Local features | `features_local.tsv` | `amp_prominence`, `amp_neighbor_ratio`, `integrated_prominence`, `amp_rel` (together **BASE4**) |
| Patches | `patches_fwd17.tsv`, `patches_bwd17.tsv`, `patches_bwd9.tsv` | Axis-aligned residual patches `S(data − model)`, forward and backward scans: 17×17 at 0.04 nm, and 9×9 at 0.08 nm (half-width 0.32 nm), on a complete symmetric disk |
| Descriptor | `features_descriptor.tsv` | `patch_u_asym_reconstructed`: transverse half-plane asymmetry of the backward 9×9 residual after affine projection |
| Mold margins | `score_fwd.tsv`, `score_bwd.tsv` | DFT constant-current molds (below): chain decoded by Viterbi over direction × pyranose parity × mirror; `mold_cc_fwd`/`mold_cc_bwd` = absolute GlcN−GlcNAc cost margin |
| Fisher margin | `fisher_cv.tsv` | `emp_fisher`: PCA (10) of forward patches, two-component mixture, Fisher direction from low- to high-amplitude cluster; lobe-parity cross-fitting, max over mirror |
| GMM | `pred_gmm.tsv` | One view `BASE4 + descriptor + mold_cc_fwd + mold_cc_bwd + emp_fisher`, pairwise interactions, per-scan standardization, full covariance (ridge 1e-6), 2 self-training rounds, Mahalanobis assignment, 10 seeds, hard vote |
| k-means | `pred_kmeans.tsv` | Four views: BASE4; BASE4 + `split_log_skew`; BASE4 + `bwd_neg_com_t`; BASE4 + `bwd_neg_diag45` (negative-lobe moments of the backward patch); interactions, 20 seeds, averaged votes |
| Soft vote | `predictions.tsv` | p = (p_kmeans + p_gmm)/2; class 1 if p ≥ 0.5; confidence = \|2p − 1\|; `?` if a component is unavailable |
| QC | `plots/`, `review_queue.tsv`, `summary.tsv` | Maps, per-chain review flags |

**Corroborated training** (`assignment_training_scans = "corroborated_counts"`).
The GMM and k-means are learned only on scans whose own per-scan count agrees
with the molecule consensus (`count_rule = "agrees"`), because their fitted
geometry is the most reliable. Every scan is still standardized and assigned.
All other settings equal `config/unit_assignment_patch_support.toml`.

### DFT constant-current molds

`templates/chitosan_cc_molds_native_v1.tsv` holds GlcN and GlcNAc molds for
each parity/mirror state. They are built by `test/build_cc_molds_native.jl` from
the accepted Quantum ESPRESSO LDOS cubes of GlcN–GlcN–GlcN and
GlcN–GlcNAc–GlcN trimers on Cu(100). The builder takes the isosurface whose
mean height is first below 0.50 nm, sampled in the central-unit frame
(`templates/qe_frames/`). The provenance file records all input hashes. The
builder keeps a legacy cube-addressing convention and reads a potentially
signed LDOS quantity, so the molds are empirical templates rather than a
validated STM simulation ([DFT molds](qe_stm_molds.md)).

Connected molds, not free asymmetric lobes: generic left/right width asymmetry
is dominated by overlap, envelope and tip effects. A mold encodes a specific
chemical geometry, and β-(1→4) connectivity restricts it to a few global states
(chain direction × pyranose parity phase × surface mirror). Within a state
each lobe is scored against both types. The number of GlcNAc units emerges from
the costs; there is no composition constraint.

## Latent-class fusion (`test/build_molecule_fusion.jl`)

For every molecule track (see [Model selection](selection.md)), the scans with
the track's modal count are aligned: each lobe is ranked along the chain axis
of the first such scan, a fixed absolute direction that is invariant to drift.
Each **physical lobe** observed in at least `fusion_min_scans = 3` scans has a
hidden type z:

- P(z = 1) = π;
- P(call 1 | z = 1) = θ1 (detection rate);
- P(call 1 | z = 0) = θ0 (false-call rate).

EM estimates π, θ0 and θ1 over all physical lobes (binomial latent-class
model), from a neutral start (0.5, 0.1, 0.9). The posterior P(z = 1 | calls) is
the fused `probability_1`; `predicted = 1` if it is ≥ 0.5, and
`confidence = |2p − 1|`. Unlike a majority vote, a type detected in a minority
of scans is kept when false calls are rare. Lobes seen in fewer scans, and
scans whose count differs from the track mode, keep their per-scan call. The
per-scan call stays in `scan_predicted`.

On the benchmark, EM gives π = 0.342, θ0 = 0.035 and θ1 = 0.622 (79 physical
lobes, 651 rows fused, 99 calls changed). **The fitted θ0/θ1 is a label-free
quality indicator.** On the unknown 10–20mers the same fit gives θ0 = 0.45:
per-scan calls there are not consistent across repeat scans, and the fused
labels are not interpretable.

## Uncertainty and outputs

`predictions.tsv` columns: `file`, `lobe`, `predicted` (0/1/?), `confidence`,
`probability_1`, `invalid_reason`, `model`, `scan_predicted`,
`scan_probability_1`, `physical_lobe`, `fusion_scans`, `fusion_calls_1`.
For fused lobes the confidence is a model posterior under the latent-class
model. For unfused lobes it is the soft-vote margin, which is **not a
calibrated probability**. `?` marks a lobe without a valid component (5 of
870 benchmark positions are unclassified, all at missing lobes).

## Benchmark results (full145 own-N, fresh raw run 12025539)

| Profile | Exact N | Correct /870 | Exact chains /145 | Classified /870 | Errors |
|---|---:|---:|---:|---:|---:|
| Historical reference (Aug 2026, not raw-reproducible) | 106 | 677 | 36 | 854 | 177 |
| Previous record (saved hybrid counts) | 129 | 694 | 43 | 855 | 161 |
| Promoted, per-scan calls | 137 | 703 | 33 | 863 | 160 |
| **Promoted, fused (final)** | **137** | **772** | **88** | **865** | **93** |

Against the record: 62 scans gain, 21 lose and 62 tie. Exact chains: 51 gained,
6 lost. Without the 43-scan molecule, the method has 514 correct / 45 exact
chains over 102 scans, against 485 / 29 for the record. 13 of 41 molecules are
exact in every scan (record: 9). An independent fresh run (12023870) with the
same fusion gives identical numbers.

## Limits

1. **Edge-adjacency confound.** In NKNNKN both GlcNAc sit next to a chain end.
   Features correlated with "second from an end" score well without detecting
   acetyl groups. Fusion amplifies whatever the per-scan calls detect; it
   neither creates nor removes this confound.
2. **Correlated exact chains.** Exact chains come in molecule-sized blocks.
3. **Coverage of fusion.** Only repeatedly imaged molecules are fused.
4. **Application.** On unknown 10–20mers the per-scan calls are inconsistent
   (θ0 = 0.45), and some long chains are undercounted. No chemical claim is
   made there.

## Grading conventions

`test/grade_unit_assignment.jl` grades per lobe under the physical convention.
It tests only the sequence orientation (identity/reverse), because the chain
axis has an arbitrary direction. The supervised "oracle" convention (best of
four alignments including the 0↔1 flip) is an upper bound, never a method.
`test/report_unit_assignment_benchmark.jl --full145-own-n` grades predictions
at each scan's own count. Missing control positions count as unclassified;
extra lobes are reported but not aligned. `test/grade_consensus_run.jl` wraps
verification and both grades for a run directory.
