# Pipeline

## Data flow

```
raw .sxm scans (one directory, unique basenames)
  │
  ├─ 1. per-scan counting ........ test/batch_full.jl (4 shards, --skip-1d)
  │      circular sweep → circ→ell refinement → GCV + robust-AICc guard
  │      → support-midpoint adjustment → N_selected           counting_summary.tsv
  │
  ├─ 2. scan geometry ............ test/extract_lobe_features.jl at N_selected
  │                                                           scan_geometry.tsv
  ├─ 3. molecule consensus ....... test/build_molecule_consensus.jl
  │      absolute-frame registration of consecutive scans → tracks
  │      → strict-majority count → fixed-N refit check (own ROI, then
  │        registered ROI) → final counts          consensus/consensus_summary_final.tsv
  │
  ├─ 4. per-scan assignment ...... test/run_reconstructed_chitosan.jl
  │      base + split fits at final N → local features → 17×17 / 9×9 patches
  │      → descriptor → DFT mold margins → empirical Fisher margin
  │      → GMM (1 view) + k-means (4 views), learned on corroborated scans
  │      → soft vote → 0/1/? + confidence                  assignment/predictions.tsv
  │
  ├─ 5. latent-class fusion ...... test/build_molecule_fusion.jl
  │      physical lobes of each track → EM (π, θ0, θ1) → posterior
  │                                                 fusion/predictions_fused.tsv
  ├─ 6. position uncertainty ..... test/report_position_repeatability.jl
  │      repeated scans aligned (registration, translation, rotation)  positions.tsv
  ├─ 7. molecule plots ........... test/plot_final_molecules.jl        plots/
  └─ predictions.tsv (= fused), chain_report.tsv, raw/input hashes, logs/
```

One command runs all seven stages:
`test/run_molecule_consensus_chitosan.jl` (Slurm wrapper
`hpc/run_molecule_consensus.sbatch`). External verification and grading are a
separate step: `test/grade_consensus_run.jl`. See the [runbook](chitosan_runbook.md).

## Run directory

| Path | Content |
|---|---|
| `raw_hashes.tsv`, `input_hashes.tsv` | SHA-256 of every raw scan and every config/template input |
| `counting_summary.tsv`, `counting_chunk*/` | Per-scan counting (merged summary, per-shard plots) |
| `scan_geometry.tsv` | Fixed-N Gaussian geometry at the per-scan counts |
| `consensus/pairs.tsv` | Consecutive-scan registrations (centre distance, NCC, drift translation, link decision) |
| `consensus/consensus.tsv` | Per-scan track, consensus count and rule |
| `consensus/consensus_summary_final.tsv` | Final count per scan, after refit checks (`count_rule`, `refit_support`) |
| `assignment/` | Per-scan assignment run: features, patches, mold scores, Fisher, `pred_gmm.tsv`, `pred_kmeans.tsv`, `predictions.tsv`, maps, `review_queue.tsv` |
| `fusion/` | `predictions_fused.tsv`, `physical_lobes.tsv`, `fusion_params.tsv` |
| `predictions.tsv` | Final per-lobe calls: `predicted` (0/1/?), `probability_1`, `confidence`, per-scan call |
| `chain_report.tsv` | One row per scan: counts, rule, track, assignment string, confidence summary |
| `positions.tsv` | Every lobe: image and absolute position, physical lobe, position uncertainty along/across the chain (1σ) and its source (`track` or `cohort`) |
| `position_repeatability/` | Per-lobe scatter, per-scan rotations, summary ([Position uncertainty](#Position-uncertainty)) |
| `plots/` | One image per molecule scan: STM image, final calls with P(GlcNAc), 2σ position ellipses, counts and consensus outcome. Folders `all/`, `consensus_agrees/`, `consensus_corrected/`, `not_checked_by_consensus/`, `with_uncertain_calls/`, plus `index.tsv` |
| `failures.tsv` | Present only if a stage failed (stage and reason) |

## Position uncertainty

`test/report_position_repeatability.jl --run RUN --data-dir RAW --outdir NEW`
measures, without labels, how well lobe positions repeat across scans of the
same molecule. It writes `lobe_scatter.tsv` (per physical lobe),
`scan_rotation.tsv` (per-scan rotation fitted between scans) and `summary.tsv`.

The scans of a molecule track are aligned in three ways: by the image
registration only, plus a per-scan translation, and plus a per-scan
translation and rotation (generalized Procrustes). Scatter is split along and
across the chain.

On the 6-mer benchmark (19 molecules imaged several times), lobe positions
repeat to **0.07 nm along the chain and 0.09–0.11 nm across it** (1σ). Spacing
between neighbouring lobes is 0.66 nm and repeats to 0.06 nm. One molecule
(track 26) rotates by about 27° on the surface during its series.

## Packages

| Package | Role |
|---|---|
| `STMFitCore.jl` | Physical constraints and scoring helpers: effective spacing, κ penalty, support overrun |
| `STMSXMIO.jl` | `SXMImage`/`read_sxm`, channel access, fwd/bwd alignment, shared preprocessing (plane fit, row flattening, smoothing, Otsu ROI) |
| `GaussianFit1D.jl` | 1D slide-profile fit. Diagnostic only (off by default); never feeds `N_selected` |
| `GaussianFit2D.jl` | 2D chain-of-Gaussians engine (`src/core.jl`): axis, support, seeding, sweep, fit |
| `STMMolecularFit.jl` | Orchestration and selectors (`src/selectors.jl`); used by `test/batch_full.jl` |

`STMSXMIO` owns the SXM types; neither fit engine redefines them.
`GaussianFit2D` depends on Core, SXM I/O and the 1D package; `STMMolecularFit`
orchestrates all four.

## Counting engine

For each scan, `GaussianFit2D.chain_gaussian_sweep`:

1. finds the molecule ROI (Otsu threshold, largest component) and its axis
   (intensity-weighted PCA);
2. detects the active support along the axis from the axial profile
   (baseline + `support_noise_k`·noise, `support_padding_nm`);
3. fits a tube of half-width `fit_width_nm` around the axis with N Gaussian
   lobes plus a tilted baseline, for each feasible N;
4. **circular sweep**: deterministic seeds from the raw axial signal, NLopt
   global search then LsqFit (σ∥ = σ⟂);
5. **circ→ell refinement**: LsqFit only, warm-started from each circular
   solution. A global elliptical search diverges from the isotropic start and
   is not used;
6. scores each N by GCV and keeps `min(circular, elliptical)` per N (the
   circular model is nested in the elliptical one).

Selection of `N_selected` from these candidates is described in
[Model selection](selection.md). Physical bounds (spacing, overlap, widths) and
the κ penalty come from the `[model]` section of the counting config
([Configuration](config.md)).

## Assignment components

See [Unit assignment](unit_assignment.md) for each feature, the classifiers,
the vote and the fusion model.
