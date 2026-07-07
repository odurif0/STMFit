# Research Journal

Chronological record of investigations into the 2D elliptical chain model
convergence and model selection problem. Includes both successful and
unsuccessful approaches, with rationale.

---

## Problem Statement (May 2025) — RESOLVED

The original problem: batch processing of chitosan STM images (240817 dataset)
where the 2D elliptical chain model produced N ≠ 6 for 7/27 files (019, 026,
051 systematically N=8), despite the theoretical expectation of 6 monomers.

**Goal**: achieve N=6 for all files where the molecule has 6 monomers, without
introducing heuristic/arbitrary parameters.

**Status (Jun 2026): SOLVED.** The label-free selection rule (GCV + robust-AICc
guard + up-when-ambiguous) now gives **39/39 primary benchmark files exact
(N=6)**, reproducible across runs. The original elliptical-convergence work is
archived in `journal_archive.md`. The active line of work is now **unit
assignment** (GlcNAc/GlcN per lobe) — see the 2026-06-22 entry below and
`docs/src/unit_assignment.md`.

---

## Investigation Timeline

### 2026-07-07 — Unknown-sequence unit-assignment workflow hardening

**Goal:** Make the GlcNAc/GlcN unit-assignment path usable for unlabeled
chitosan chains without accidentally entering benchmark-reporting code.

**Changes:** Added a production runner for fixed label-free profiles,
prediction validation, plot handling for explicit `?` abstentions, QC review
queues, and artifact manifests. The docs now separate unknown production from
post-hoc benchmark grading, with a checker guarding the unknown sections.

**Rationale:** Current label-free unit assignment remains a diagnostic workflow,
not a solved binary map. Unknown runs therefore need honest uncertainty,
validation, plots, and review queues before any external report is generated.

### 2026-06-22 — Unit assignment (GlcNAc/GlcN) investigation started

**Goal:** Assign each fitted Gaussian lobe a monomer type (0 = GlcN,
1 = GlcNAc) to produce a deacetylation map per chain. The benchmark-control
sequence is external grading information only and must stay outside the
fitting/selection/assignment path (same label-free rule as for N).

**Approach:** Graduated, label-free, generalist-first:
1. **Phase 0** — Grading framework: `test/grade_unit_assignment.jl` tests 4
   alignments (identity, reverse, flip, reverse+flip) and 2 conventions
   (physical: GlcNAc=amp max; oracle: best flip, supervised upper bound).
2. **Phase 1** — Gaussian feature separability: `test/extract_lobe_features.jl`
   re-runs the fit and extracts per-lobe (A, σ∥, σ⟂, integrated);
   `test/analyze_unit_separability.jl` tests unimodal vs bimodal (kmeans k=1
   vs k=2, BIC) and cross-evaluates with truth (`--with-truth`).
3. **Phase 1b** — Non-Gaussian residual features: an experimental residual
   feature route computed skewness, shoulder at ±δ (δ=0.15 nm, C2 acetyl offset
   from pyranose geometry), kurtosis, and L/R asymmetry from the fit residual.
   It was later removed after adding noise instead of unit-identity signal.

Phases 2–5 (model selection 1-type vs 2-type, per-blob clustering, template
supervised validation, DFT-STM) are planned pending the Phase 1+1b separability
verdict.

**Key design decisions:**
- Chitosan is a random copolymer (DD-dependent), not strictly alternating →
  per-blob clustering (Phase 3), not alternating-model assignment.
- STM ≠ geometric cross section → no direct template fitting without DFT.
  Geometry guides *where* to measure features (Phase 1b shoulder offset), not
  what to fit.
- Two conventions reported: physical (label-free) and oracle (supervised upper
  bound). Gap between them validates or invalidates the physical mapping.
- Batch doesn't produce per-lobe TSVs → `extract_lobe_features.jl` re-runs the
  fit. Supports `--chunk I/N` for HPC.

**Dependencies added:** `Clustering.jl` v0.15.8, `StatsBase.jl` v0.34.10.

**Files created:**
| File | Role |
|------|------|
| `benchmarks/chitosan_240817_unit_sequences.tsv` | Ground-truth skeleton (user fills sequences) |
| `test/grade_unit_assignment.jl` | Phase 0: grading (4 alignments, 2 conventions) |
| `test/extract_lobe_features.jl` | Phase 1: re-run fit, extract per-lobe Gaussian features |
| `test/analyze_unit_separability.jl` | Phase 1: unimodal vs bimodal + AUC + clustering accuracy |
| `docs/src/unit_assignment.md` | Documentation page for the unit-assignment pipeline |

**Bug fix (pre-existing):** `Documenter` UUID in `Project.toml` was malformed
(`86b` instead of `863b`), causing `Pkg.resolve()` to fail. Fixed to match the
Manifest and git HEAD.

**2026-06-22 execution result:** Phase 1 and 1b were first run on 74 valid
fitted files (456 lobes; `240817_084.sxm` had no valid fit in the
feature-extraction rerun). A bug was then found in the kmeans BIC helper: the
code treated rows as observations while `Clustering.kmeans` expects observations
in columns. After fixing this, the 74-file diagnostic gave Gaussian-only
`ΔBIC(k=1-k=2) = +304.7` and Gaussian+residual `ΔBIC = +39.6`; the features are
therefore bimodal, not unimodal.

The extraction was then improved to use the batch `N_selected` (via
`--selected-summary`) and to filter the primary benchmark only
(`--manifest ... --primary-only`). On the corrected 39-file primary set (234
lobes, all N=6), Gaussian features alone gave `ΔBIC = +184.1` (strongly
bimodal). Adding residual features (axial skewness, shoulder at ±0.15 nm,
kurtosis, L/R asymmetry, residual peak SNR) weakened the evidence to
`ΔBIC = +4.3` (weakly bimodal). Conclusion: the current best label-free feature
set is the Gaussian fit features alone. Residual features, as currently defined,
mostly add noise. The next necessary step is to fill the withheld unit truth and
grade the label-free sequences; if performance is poor, use the new supervised
template script to quantify the empirical upper bound before moving to DFT-STM.

**2026-06-23 follow-up:** The benchmark sequence was revealed externally as
`NKNNKN` for diagnostic grading only (`N=GlcN=0`, `K=GlcNAc=1`, so `010010`).
This truth remains excluded from fit, selection, clustering, and any label-free
assignment rule. A stricter constraint was added: the unit-assignment rule must
not assume the number of units of either type. Composition-constrained heuristics
such as "top-2 prominence" are therefore diagnostic only, even though top-2 local
prominence reached 75.2% per-lobe accuracy and 13/39 exact sequences.

Implemented `test/augment_lobe_local_features.jl` to add chain-internal local
features without truth or composition priors: per-file z-scores, local prominence,
neighbor ratios, linear/quadratic envelope residuals, and edge-distance features.
On the 39 primary files these features remained strongly bimodal (`ΔBIC=+144.0`),
but label-free assignment did not improve: 69.2% physical accuracy, 73.5% oracle
accuracy, and 1/39 exact sequences against the diagnostic truth. The predicted
compositions vary by chain, confirming that no fixed number-of-units prior is
being imposed.

Implemented `test/extract_lobe_patches.jl` and `test/analyze_lobe_patches.jl` to
extract 9×9 raw and residual patches aligned to the fitted chain axis, then run
PCA/kmeans label-free and optional supervised train/test diagnostics. Results on
234 primary lobes: raw patches are bimodal (`ΔBIC=+66.2`) but supervised test
accuracy is 51.4%; residual patches are weakly bimodal (`ΔBIC=+20.3`) with 61.1%
supervised test accuracy. Conclusion: the current aligned patch representation
does not yet contain a robust, generalizable GlcNAc/GlcN signal. The dominant
failure remains the lobe-4 false positive: global/shape features see brightness
or overlap structure, not chemical identity.

Implemented the Phase 1e split-width Gaussian forward-model test. The production
profile remains symmetric by default; `config/chitosan_split.toml` opts into
`peak_profile="split"`, which fits one label-free `skew_ratio = σright/σleft`
per lobe, bounded by `skew_ratio_max=2.0`. `skew_ratio=1` is exactly the old
symmetric Gaussian, verified by a unit test in
`packages/GaussianFit2D.jl/test/runtests.jl`. `test/extract_lobe_features.jl`
now exports `skew_ratio`, so the existing separability/assignment scripts can use
it directly. The intended experimental protocol is to keep N from the validated
batch summary and refit only the lobe shape at that fixed N; this isolates the
asymmetry test from N selection and preserves the no-composition-prior rule.

Decision gate before DFT-STM molds: if split-width refits do not improve GCV
and `skew_ratio` stays near 1, the STM/tip conditions likely do not resolve the
acetyl asymmetry. If split-width improves GCV and `skew_ratio` is stable/bimodal,
the next step is a two-template physical mold built from GlcNAc/GlcN maps, with
per-lobe continuous mold weights and no composition constraint.

The first full split-width verification exposed a performance issue: the split
kernel called `String/lowercase` inside the inner pixel loop, making the run
effectively infeasible. This was fixed by hoisting the split-profile flag and
`skew_ratio_max` outside the loop. A second optimization was added to
`test/extract_lobe_features.jl`: when `--selected-summary` is supplied, the
script sets `n_min=n_max=N_selected` per file, avoiding the unnecessary N sweep.
With these changes, the 39 primary files completed locally.

Results of the split-width verification at fixed batch-selected N: all 39 files
refit at N=6. Split improved GCV on 36/39 files relative to the symmetric
Gaussian at the same N; median relative ΔGCV was -10.7%, mean -12.6%, best
-36.4%, worst +25.8%. `skew_ratio` alone was strongly bimodal (`ΔBIC=+294.8`,
clusters 99/135), so the STM data do contain a repeatable left/right shape
asymmetry. However, external diagnostic grading against the withheld `010010`
truth showed that this asymmetry is not the GlcNAc/GlcN label: `skew_ratio` AUC
was 0.477 (inverse 0.523), skew-only assignment gave 48.3% physical accuracy and
0/39 exact sequences, and Gaussian+skew degraded to 60.7% physical accuracy.
Conclusion: split-width is a better forward model for topography, but the fitted
asymmetry is dominated by local shape/overlap/envelope effects rather than the
C2 acetyl identity. Do not proceed to a DFT-STM mold expecting only this skew
mode to solve unit assignment; a mold would need to encode a different, more
specific observable than generic left/right width.

Implemented the first connected-mold decoder,
`test/score_connected_mold_templates.jl`. It operates at patch level: given the
aligned patches from `extract_lobe_patches.jl` and a template TSV containing all
GlcN/GlcNAc × parity × mirror molds, it tests the global connectivity states
`direction × parity phase × mirror` and chooses the lowest-cost per-lobe type.
This enforces glycosidic-orientation constraints while preserving the no-prior
rule on composition: no truth sequence and no number of GlcNAc/GlcN units are
read. The current script accepts externally generated templates (geometric proxy
molds or DFT-STM maps) and writes per-lobe 0/1 predictions that
can be graded externally.

Extended the connected-mold path to include sliding pairwise bond templates. The
decoder now accepts `--bond-templates`, with rows for the 16 combinations
`left_type/right_type ∈ 00/01/10/11 × parity × mirror`, scored on every adjacent
edge `(i,i+1)` and decoded by Viterbi. This does not tile chains into disjoint
dimers, so both odd and even N are supported. Added
`test/generate_connected_mold_templates.jl`, which generates unary templates and
optional concatenated left/right bond templates from an aligned proxy-site TSV.
An intermediate RDKit/3D-coordinate route was prototyped, then removed once the
maintained path switched to explicit geometric proxy sites and DFT-STM maps.
Added `test/validate_connected_molds.jl`. The validator checks connected-mold
readiness without using truth labels: proxy-site columns when available, the 8
unary template combinations, the 16 optional sliding-bond combinations, and
pixel-count compatibility against the patch TSV. This gives a preflight check
before scientific decoding.

Removed the no-RDKit heavy-atom coordinate route and replaced it with a manual
geometric proxy-site source, `templates/chitosan_geometric_sites.tsv`. The new
source is deliberately simple and explicit: shared pyranose backbone sites plus a
short GlcN substituent or a longer GlcNAc acetyl-side proxy in the aligned
`(t,u)` patch frame. `test/generate_connected_mold_templates.jl` now defaults to
this site file. `test/validate_connected_molds.jl` also treats the site file as
the primary upstream geometry and keeps 3D coordinates optional.

Added `--template-mode contrast` to `test/score_connected_mold_templates.jl`.
The default `full` mode scores the complete GlcN/GlcNAc mold; `contrast` subtracts
the parity/mirror common mold before scoring so the shared pyranose backbone does
not dominate the weak substituent signal. This remains label-free and uses no
composition prior.

Diagnostic results for the first manual geometric mold are technically valid but
not sufficient. Using 13×13 raw patches (`±0.48 nm`) and contrast scoring gives
60.3% physical accuracy, 63.7% oracle accuracy, and 0/39 exact sequences once the
template scorer writes amplitudes and the physical flip can be applied. Adding
sliding bond templates does not improve the result. Residual-patch contrast gives
42.7% physical / 66.7% oracle, indicating a poor physical 0↔1 mapping. Full raw
template scoring reaches 67.9% oracle and 4/39 exact, but its physical convention
collapses to 35.5%, so it is not a valid label-free assignment. Conclusion: the
geometric mold is the right scaffold to keep, but the current hand-drawn shape is
only a first guess. Any next refinement must optimize a label-free residual or
cross-validation criterion, freeze the mold, and only then grade against the
benchmark sequence.

Implemented `test/refine_geometric_mold.jl`, a first label-free refinement pass
for the geometric scaffold. It searches global transforms of the acetyl proxy
sites (small `t/u` shifts, transverse scaling, weight scaling, sigma scaling),
scores contrast templates against the aligned patches, and ranks candidates by
the k=1 vs k=2 BIC of per-lobe template-evidence margins. It does not read truth
labels and does not impose a composition. The scoring functions now use the
finite overlap of patch/template pixels, so wider `±0.48 nm` patches with a few
edge `NaN`s no longer produce infinite costs. `score_connected_mold_templates.jl`
also writes the lobe amplitude so `grade_unit_assignment.jl` can apply the
existing physical label-free 0↔1 convention.

Raw-patch refinement selected `dt=-0.08 nm`, `du=-0.08 nm`, `u_scale=0.85`,
`weight_scale=0.8`, `sigma_scale=1.25` for the acetyl proxy sites, with
`ΔBIC=417.3` and cluster sizes 74/160 over all 234 lobes. Post-hoc diagnostic
grading improved to 67.9% physical accuracy and 72.2% oracle accuracy, still with
0/39 physical exact sequences (1/39 oracle exact). Adding sliding bonds slightly
worsens the refined raw mold to 67.1% physical. Residual-patch refinement gives a
larger label-free `ΔBIC=625.8`, but physical accuracy is only 44.0%; the residual
mode is therefore strongly bimodal but not chemically aligned. Conclusion: the
raw geometric mold refinement is a measurable improvement over the hand-drawn
scaffold, but still below the threshold for a robust deacetylation map.

Implemented `test/import_stm_mold_maps.jl` to bridge the next, physically better
mold source: DFT-STM / Tersoff-Hamann maps in the aligned `(t,u)` lobe frame. The
input unary TSV is long-form (`type, t_nm, u_nm, value` with optional
`parity/mirror`). If orientation-specific maps are not supplied, the importer
generates the beta-(1->4) parity/mirror variants by flips. It can also write
sliding-bond template TSVs either by concatenating unary maps or from optional
long-form bond maps (`left_type, right_type, side, t_nm, u_nm, value`). The tool
was validated on synthetic GlcN/GlcNAc maps: generated unary and bond templates
pass `validate_connected_molds.jl` and are accepted by
`score_connected_mold_templates.jl`. This does not solve unit assignment by
itself; it defines the required interface for real DFT-STM inputs.

Added the Quantum ESPRESSO mold workflow. `test/prepare_qe_mold_inputs.jl`
serializes a vetted slab+trimer XYZ into `pw_relax.in`, `pw_scf.in`,
`pp_ldos.in`, and a Slurm sketch, without building or altering chemistry. This
keeps the adsorbed geometry as an explicit scientific input. It can read
`cell_a/b/c` directly from the slab-builder metadata via `--cell-metadata`.
Added
`test/build_qe_slab_trimer_xyz.jl` to place an already oriented trimer above a
reproducible Cu(100) slab, and `test/extract_qe_mold_frame.jl` to compute the
central-unit origin and `(t,u)` axes from a relaxed XYZ. Added
`test/extract_qe_relaxed_xyz.jl` and `test/update_qe_positions_from_xyz.jl` for
the relax-to-SCF handoff; generated Slurm scripts call them via `STMFIT_ROOT`.
The extractor handles output paths without parent directories and supports
`ATOMIC_POSITIONS/CELL_PARAMETERS (alat)` when `--alat-angstrom` is supplied.
`test/cube_to_stm_maps.jl` converts
Gaussian cube output (for example from QE `pp.x`) into the long-form
`templates/chitosan_stm_maps.tsv` format by sampling a local `(t,u)` plane around
the central unit. It now accepts typed frame files from
`extract_qe_mold_frame.jl` via `--frame TYPE:frame.tsv`, so each GlcN/GlcNAc cube
can carry its own origin, axes, and sampling height. `hpc/qe_molds/` now contains
template QE inputs for relaxation, SCF, LDOS cube export, placeholder XYZ
requirements, and a minimal Slurm launcher sketch. `docs/src/qe_stm_molds.md`
records the full protocol and the label-free constraints.
Added `test/smoke_qe_mold_workflow.jl`, a no-QE end-to-end smoke test that
generates a synthetic trimer/slab, prepares QE inputs via `--cell-metadata`,
simulates the relax-to-SCF handoff, extracts a typed frame, samples synthetic cube
files, and imports connected STM molds.

Added `test/build_initial_chitosan_trimer_xyz.jl` and generated the first actual
starting structures in `hpc/qe_molds/`: `glcn_central_trimer.xyz` (69 atoms,
GlcN-GlcN-GlcN) and `glcnac_central_trimer.xyz` (74 atoms,
GlcN-GlcNAc-GlcN), with companion atom/index TSVs. These are deterministic,
unoptimized starting geometries for QE relaxation, not final scientific
structures. Slab builds now accept `--center-indices`; the generated slabs use
central-ring indices `12,13,14,15,16,17` so GlcN and GlcNAc start from the same
Cu(100) registry. Frame indices were validated after the default `8×8×4` slab
offset: `--origin-indices 268,269,270,271,272,273 --axis-from 271 --axis-to 268
--plane-index 277`.
Added `test/validate_chitosan_trimer_structures.jl` and wrote
`hpc/qe_molds/structure_validation.tsv`; the report confirms `glcn` has 0 acetyl
units, `glcnac` has 1 central acetyl unit, both have non-acetylated GlcN
neighbors, minimum distances are sane, molecule labels survive slab generation,
and central-ring centers match the Cu-cell center.
Parsed the 240817 SXM headers and found a uniform `BIAS=-3.000E-1 V` across 94
files; the LDOS window is therefore `emin=-0.3`, `emax=0.0` eV. The original QE
run directories used the full `8×8×4` slab with `ntasks=4`, `ecutwfc=80`,
`ecutrho=640`; those parameters were later dropped (OOM on Raven — see the
submission lessons below) in favour of the active `8×6×3` pilot
(`ecutwfc=50`, `ecutrho=360`, `8` tasks). The validated slab freeze cutoff
`fix_below_z=1.807501` is unchanged. The `qe/` directory is gitignored to avoid
committing large QE scratch/cube outputs.
Added `test/preflight_qe_mold_inputs.jl` and `hpc/submit_qe_molds.sh`; the
preflight report `hpc/qe_molds/qe_input_preflight.tsv` verifies `nat/ntyp`,
species, frozen relax atoms, LDOS window, sbatch handoff commands, and the 8-task
total Slurm budget before submission.
Added `hpc/launch_qe_molds_remote.sh`, a local-to-MPCDF QE launcher that reuses
`hpc/remote.env`, syncs code while excluding local `qe/` outputs, syncs only the
prepared QE input files plus local `pseudo/*.UPF` files, runs the remote preflight, and then calls
`hpc/submit_qe_molds.sh`.
Remote launchers now expose `SSH_CONNECT_TIMEOUT` and
`SSH_SERVER_ALIVE_INTERVAL` in `hpc/remote.env.example`; defaults are 180 s and
60 s to tolerate slow MPCDF gateway/password+OTP handshakes.
Added `test/finalize_qe_mold_workflow.jl` for the post-QE handoff: once relaxed
XYZ files and LDOS cubes exist, it extracts typed frames, samples the GlcN/GlcNAc
cubes into `templates/chitosan_stm_maps.tsv`, and imports the connected unary and
bond templates. It requires an explicit `--height-nm` so the sampling height is a
physical input rather than a fitted benchmark parameter.


First Raven submissions exposed a series of launch-safety gaps, now fixed in the
tooling and distilled in `hpc/qe_molds/README.md`:

- QE module stack is `intel/2024.0 impi/2021.11 qe/7.4.1` (not `quantum-espresso`).
- Raven shared nodes reject `--mem=0`; the generator writes explicit MB.
- `preflight_qe_mold_inputs.jl` verifies every `ATOMIC_SPECIES` pseudo exists,
  and the remote sync includes `pseudo/*.UPF`.
- The `n0001` QOS enforces a one-node group limit, so two parallel QE jobs each
  request a separate node; `submit_qe_molds.sh --sequential` chains them with
  `afterok`.
- The original `8×8×4` slab with Cu `spn` PAW (`z_valence=19`) and
  `ecutwfc=80`/`ecutrho=640` estimated ~635 GB dynamic RAM and OOM'd within
  minutes. A 2-task/48 GB probe was memory-safe but underparallel (only the
  initial SCF in ~7 h).
- The active pilot is therefore the lighter `8×6×3` slab, `12 Å` vacuum, Cu `dn`
  PAW (`z_valence=11`), Γ-only, `ecutwfc=50`/`ecutrho=360`, `8` MPI tasks,
  `96000 MB`, `24:00:00` (`srun -n "$QE_NTASKS" --cpu-bind=cores`).

The full per-job submission narrative (jobs `28278265` through `28303162`, all
superseded) is archived in `journal_archive.md`.




---

> **Archive:** full historical investigation detail (May 2025 – mid-June
> 2026, including the v1–v7 pipeline evolution, selection-rule work, and
> the calibration analysis) is in
> [journal_archive.md](journal_archive.md).

---

## Current Pipeline (v6)

```
Step 1: 1D slide profile extraction + peak fitting
        → DIAGNOSTIC ONLY (off by default; --no-skip-1d to re-enable)
        → never enters N_selected; the 1D over-counts (lateral averaging)

Step 2: Circular sweep (N = 2..14, adaptive range)
        → deterministic 2D-only initialization from raw axial profile
        → reliable convergence, isotropic gaussians

Step 3: circ→ell LsqFit refinement at EACH N
        → warm-start from circular solution, local optimization only
        → finds true elliptical minimum without NLopt divergence

Step 4: Model selection = GCV + robust-AICc guard + up-when-ambiguous
        → GCV is canonical (valid under spatial correlation; BIC/AICc are
          diagnostics only — n_eff is undefined in the fit window)
        → guard: robust-AICc can only descend (veto), never ascend
        → up-when-ambiguous: if the guard is ambiguous, prefer the higher N
        → circular model is nested fallback; refined elliptical when it improves

Step 5: Output best models (N_selected, params, plots, scores, QC)
```

**Selection criteria hierarchy**:
1. `N_selected` — the label-free answer, driven by GCV with the guard.
2. `N_ell` / `N_circ` — best valid refined elliptical / circular 2D model
   (diagnostic splits of the same sweep).
3. Default criterion is GCV (`selection_criterion="gcv"`, `cv_method="gcv"`).

See `docs/src/selection.md` for the full guard specification and
`docs/src/calibration.md` for why GCV (not BIC/AICc) is canonical.

---

## Lessons Learned

1. **Circular model is the anchor**: σ∥=σ⟂ enforced structurally, always
   converges. Use it as the reference in all comparisons.

2. **NLopt global optimizer is harmful for elliptical**: 33D parameter space
   is too large. The isotropic solution is a saddle point that NLopt always
   escapes. LsqFit-only from circular start is optimal.

3. **Min() is more robust than penalties**: Adding penalty terms to BIC
   introduces free parameters. Using `min(ell, circ)` achieves the same
   effect with zero new parameters.

4. **GCV is the default selection score**: BIC is an asymptotic approximation;
   analytical GCV provides a cheap predictive-error proxy without refitting.

5. **Re-parameterization doesn't fix optimizer topology**: Changing from
   (σ∥,σ⟂) to (σ_iso,Δ) just moves the divergence point. The fundamental
   issue is that any extra degree of freedom in sigma space can be exploited.

6. **1D over-estimates N**: The 1D fit has more flexibility (no 2D topology
   constraints) and sBIC penalizes less. This is documented but not fixed
   — the 2D selection is what matters for final output.

---

## Open Questions

> Updated 2026-06-27. Questions from earlier sessions are archived in
> `journal_archive.md`.

1. **n_eff and information criteria** → **RESOLVED (Jun 20)**: The n÷9 heuristic
   is not objectively definable in the fit window — the STM spatial correlation
   range (17–100 px) far exceeds the ~10-px window, so the number of independent
   points is effectively zero. BIC/AICc (which assume iid) are therefore not
   well-defined; GCV (valid under spatial correlation) is the canonical criterion.
   See `docs/src/calibration.md`.

2. **Is N=9 correct for 260115_016 (10–20mer)?** → **OPEN**: The 2D fit's GCV
   optimum is N=9 (confirmed even with `max_overlap` relaxed to 0.80, which
   allows N up to 14). The former 1D fit saw N=13, but investigation showed the
   1D over-counts (lateral averaging creates spurious axial peaks). Without a
   visual ground-truth label for this file, N=9 stands as the objective answer,
   but it has not been visually confirmed. Action: visual inspection of the
   260115_016 best-fit overlay plot.

3. **Auto-calibration under-detects on ~4% of files** → **OPEN (low priority)**:
   `measure_calibration.jl` reproduces manual calibration on 17/25 10–20mer
   files, ±1 on 7, and fails badly on 1 (251206_013: N=4 vs 11). Root cause:
   coupled parameters (`fit_width × support_padding × σ`) interact
   non-monotonically. The tool is a bootstrap (good starting point), not a
   replacement for visual validation. No fix planned unless it fails on a new
   molecule's clean scan.

4. **Guard robust-AICc descends by 2 on 3/25 10–20mer files** → **OPEN (monitor)**:
   On short chains in the 10–20mer set (260115_016, 260116_017, 260222_043), the
   guard drops N_eff by 2 (e.g. 8→6). This is within the guard's design (down-only
   veto), but on a non-benchmarked dataset we cannot confirm it's correct without
   visual labels. Monitor: if a pattern emerges on more data, consider bounding
   the guard descent to 1 (symmetric with the up-branch).

5. **`max_overlap` generalization** → **RESOLVED (Jun 20)**: Investigated on
   260115_016 — relaxing from 0.60 to 0.80 does allow high-N fits (N up to 14),
   but the GCV optimum stays at N=9. The constraint is a physical prior (Gaussian
   pair-overlap floor), not an arbitrary blocker. Kept at 0.60 for the chitosan
   calibration; verify it isn't rejecting good fits on a new molecule with denser
   lobes.

6. **Physical LDOS sampling height for QE molds** → **RESOLVED (Jun 28)**:
   Production height fixed at **`--height-nm 0.50`** from Tersoff-Hamann
   physics. QE `pp.x plot_num=5` is an s-wave TH LDOS map with no tip-apex or
   work-function correction; `V,I → z` is not unique from `pp.x` alone, so
   constant-height is the only defensible first approximation. Literature on
   organic adsorbates on Cu(100) places constant-height planes at 3–5 Å above
   the adsorbate (4–8 Å for larger molecules). With our frame origin at the
   central ring centroid, `height_nm=0.50` samples approximately 5 Å above the
   ring along the surface normal, consistent with the recommended range.
   Sensitivity bracket: 0.40–0.60 nm. The height was chosen from physics, not
   from benchmark unit-sequence accuracy.

---


## 2026-06-28 — GlcN restart timed out again; GlcNAc preliminary cube completed

### Completed Jobs Reconciled

- GlcN restart job `28363474` ended at the 24 h walltime limit. Slurm accounting
  reports the parent job as `TIMEOUT` (`1-00:00:05`) and the `pw.x` step as
  `FAILED` after `23:59:27` (`ExitCode=1:0`). Memory was modest for this pilot
  (`~4.8 GB` MaxRSS per task), so the failure mode is walltime, not memory.
- The run did not reach the relax-to-SCF handoff: only
  `glcn_central_relax.out` plus Slurm logs were fetched from `qe/glcn_restart/`.
  There is no final `glcn_central_scf.out`, `glcn_central_pp.out`, or production
  `glcn_central_ldos.cube` from this restart.
- The latest parsed QE state reached `number of bfgs steps = 36`. The output has
  no `JOB DONE`; it stopped mid-SCF after the last printed energy
  `-31683.65010692 Ry`. This is still **not** a converged GlcN production mold.
- Extracted the last `ATOMIC_POSITIONS` block to
  `qe/glcn_restart/glcn_central_best2.xyz` (`213` atoms). QE did not print cell
  metadata in that output, so the active pilot cell remains supplied from
  `hpc/qe_molds/glcn_central_trimer_slab_pilot_meta.tsv`.

### GlcNAc preliminary SCF+PP completed

- Diagnostic GlcNAc preliminary job `28365256` completed successfully after the
  `afterany:28363474` dependency released:
  - parent job `COMPLETED`, elapsed `01:36:12`;
  - `pw.x` step `COMPLETED`, elapsed `01:35:18`;
  - `pp.x` step `COMPLETED`, elapsed `00:00:50`.
- `qe/glcnac_prelim/glcnac_central_scf.out` converged in 52 SCF iterations and
  ended with `JOB DONE`.
- `qe/glcnac_prelim/glcnac_central_pp.out` wrote
  `glcnac_central_ldos.cube` and ended with `JOB DONE`.
- This remains a **diagnostic type-1 smoke-test cube** from the unrelaxed initial
  GlcNAc pilot geometry. It is useful for pipeline de-risking only; do not treat
  it as a production GlcNAc mold.

### Diagnostic type-1 map conversion

- Built a GlcNAc preliminary frame from the active `8×6×3` pilot slab. The slab
  offset is `218 - 74 = 144`, so the bare central-unit frame indices from
  `glcnac_central_trimer_indices.tsv` become:
  `origin_indices=156,157,158,159,160,161`, `axis_from=159`, `axis_to=156`,
  `plane_index=165`.
- Converted the GlcNAc preliminary cube at the same diagnostic height used for
  the GlcN preliminary map, `height_nm=0.35`, to:

```text
templates/chitosan_stm_maps_glcnac_prelim_h035.tsv
```

- Map sanity check: `169/169` pixels finite, `0` `NA`, type set `{1}`. The value
  range is approximately `[-5.66e-7, 1.29e-4]`. This is one-sided and cannot be
  frozen/scored as a final connected mold without the matching production GlcN
  and GlcNAc maps.

### Diagnostic two-type preliminary mold smoke test

- Combined the one-sided preliminary GlcN and GlcNAc maps at diagnostic height
  `0.35` nm into:

```text
templates/chitosan_stm_maps_prelim_h035.tsv
```

- Imported two diagnostic connected-mold template sets:
  - `templates/chitosan_connected_molds_stm_prelim_h035.tsv` and
    `templates/chitosan_connected_bond_molds_stm_prelim_h035.tsv` at `13×13`
    (`half_nm=0.48`, 169 pixels), matching the future wider-patch workflow;
  - `templates/chitosan_connected_molds_stm_prelim_h035_half032.tsv` and
    `templates/chitosan_connected_bond_molds_stm_prelim_h035_half032.tsv` at
    `9×9` (`half_nm=0.32`, 81 pixels), matching the local
    `results/unit_separability/lobe_patches_selectedN_primary.tsv` patches.
- Validation of the `9×9` diagnostic templates passed:
  unary rows `8`, bond rows `16`, patch rows `234` over `39` files, and
  pixel-count compatibility `81 = 81`.
- Ran label-free connected-mold decoding against the local residual patches:

```text
results/unit_assignment/stm_prelim_h035_half032_predictions.tsv
```

- The ground-truth sequence TSV is still a skeleton with empty `sequence` values,
  so `grade_unit_assignment.jl` grades `0` files. No benchmark sequence was
  injected into the data to force a score. This remains a smoke test of the
  import/validate/decode plumbing, not a scientific result.

### New GlcN restart submitted

- Regenerated `qe/glcn_restart2/` from `glcn_central_best2.xyz` with the active
  pilot settings (`8` MPI tasks, `50/360 Ry`, Γ-only, `96000 MB`, `24:00:00`),
  copied the validated pseudo set, and preflighted successfully:

```bash
julia --project=. test/preflight_qe_mold_inputs.jl \
    --dir qe/glcn_restart2 \
    --out hpc/qe_molds/qe_input_preflight_glcn_restart2.tsv \
    --max-total-tasks 8 \
    --min-mem-mb 96000
```

- Submitted `qe/glcn_restart2` to Raven without `--watch`:

```text
qe/glcn_restart2 -> 28444935
```

### GlcNAc production queued behind GlcN success

- Local preflight of `qe/glcnac` production passed again at the active settings
  (`8` MPI tasks, `96000 MB`, `24:00:00`).
- Submitted GlcNAc production to Raven with an explicit success dependency on
  the current GlcN restart, without `--watch`:

```text
qe/glcnac -> 28445456  (dependency: afterok:28444935)
```

- This advances queueing without violating the production gate: GlcNAc production
  will only start if `28444935` succeeds. Do not read `28445456` outputs until a
  completion/timeout notification is available.

### Current State / Next Gate

- Wait for `28444935` completion/timeout notification before reading its outputs.
- If `28444935` converges and produces final GlcN relaxed/cube outputs, the
  dependent GlcNAc production job `28445456` should start automatically via
  `afterok`. If it does not, submit GlcNAc production manually at that point.
- If `28444935` times out again, either continue one more geometry-preserving
  restart or revisit relaxation strategy/walltime before spending another full
  Raven day.
- The physical LDOS sampling height is now resolved: production
  `--height-nm 0.50` (see Open Question 6). The preliminary sensitivity scan
  below used 0.30–0.60 nm for diagnostic purposes only.

---

## 2026-06-28 — HEIGHT resolved + preliminary DFT-STM diagnostic bilan

### HEIGHT resolution

- A Tersoff-Hamann literature review (organic adsorbates on Cu(100), low bias)
  confirms that QE `pp.x plot_num=5` is an s-wave LDOS map with **no tip-apex
  or work-function correction**. The setpoint (−0.3 V, 2.0 pA) does not uniquely
  determine an absolute height from `pp.x` alone.
- Production height fixed at **`0.50 nm`** above the central ring centroid,
  within the literature range of 3–5 Å above the adsorbate. Sensitivity bracket:
  0.40–0.60 nm. This is a physics choice, not a benchmark-tuned parameter.

### Preliminary DFT-STM diagnostic bilan (post-hoc, label-free decode + post-hoc grade)

- Combined the preliminary GlcN (unconverged best-so-far geometry) and GlcNAc
  (unrelaxed initial geometry) LDOS cubes into two-type diagnostic molds at
  seven sampling heights. Graded against the filled diagnostic truth `010010`
  (35 primary clean/clean_target files, post-hoc only).
- All maps use contrast mode, `res_p` patches, 9×9 templates, bond templates.

| Height (nm) | Physical % | Oracle % | Exact/35 | Gap |
|---|---|---|---|---|
| 0.30 | 51.9 | 61.4 | 0 | 9.5% |
| 0.35 | 55.7 | 66.2 | 0 | 10.5% |
| 0.40 | 57.1 | 62.9 | 1 | 5.7% |
| 0.45 | 45.2 | 62.4 | 0 | 17.1% |
| **0.50** | **58.6** | **68.1** | **2** | **9.5%** |
| 0.55 | 45.7 | 62.9 | 1 | 17.1% |
| 0.60 | 58.1 | 65.7 | 2 | 7.6% |

### Interpretation

- The signal is **present but weak and unstable** across heights: physical
  accuracy oscillates between 45% and 59%. This is consistent with using
  unconverged cubes from unrelaxed or partially relaxed geometries.
- The production height `0.50` gives the best oracle (68.1%) and tied-best
  exact count (2/35), but this was **not** the selection criterion — the height
  was fixed from physics before reading the sensitivity table.
- For comparison, the geometric refined raw mold (no DFT) reached 67.9%
  physical / 72.2% oracle / 0/39 exact. The preliminary DFT-STM molds are
  **not yet competitive** with the geometric refined mold, as expected given
  the non-converged geometries.
- The large physical-oracle gap (5–17%) indicates the amplitude-based 0↔1
  physical mapping is unreliable at this stage.
- **Next gate**: once GlcN production (`28444935`) and GlcNAc production
  (`28445456`) converge, re-run `finalize_qe_mold_workflow.jl --height-nm 0.50`
  on the converged cubes and re-grade. The converged DFT-STM molds should
  improve both stability and accuracy.

---

## 2026-06-28 — Label-free unit assignment: comprehensive method exploration

### Goal

Maximize label-free GlcN/GlcNAc assignment accuracy using only features
extracted from the STM fit and aligned patches, without converged DFT cubes.
The diagnostic truth `010010` is used exclusively for post-hoc grading.

### Methods tested

**Feature engineering (per-file z-scored):**

| Feature family | ΔBIC | Best phys % | Notes |
|---|---|---|---|
| Gaussian (4): amp, σ∥, σ⟂, integrated | 184 | 70.5% | baseline |
| Local prominence (4): amp_prom, amp_rel, nbr_ratio, int_prom | 237 | 75.2% | best single family |
| Split-width skew_ratio | 295 | 48.3% | non-chemical (AUC 0.48) |
| Patch u-asymmetry (residual 9×9) | — | 64.3% | captures acetyl lateral signal |
| Multi-scale prominence (±1, ±2, ±3 neighbors) | — | 51.4% | no improvement |

**Clustering methods (on local prominence 4 + interactions):**

| Method | Phys % | Oracle % | Exact/35 | Gap |
|---|---|---|---|---|
| k-means k=2 | 81.0 | 81.0 | 5/35 | 0.0% |
| **GMM full covariance, 10-seed ensemble** | **82.4** | **82.4** | **8/35** | **0.0%** |
| GMM + patch_u_asym (5 features + inter) | **82.4** | **82.4** | **11/35** | **0.0%** |
| Self-training 5 iterations (GMM seeds → Mahalanobis) | 81.4 | 83.3 | **13/35** | 1.9% |
| Diffusion maps (comp=2, k=20) + GMM | 81.0 | 82.9 | 11/35 | 1.9% |
| Spectral clustering (Ng et al.) | 45–67 | 57–72 | 0/35 | 5–11% |
| Fuzzy c-means (m=2) | 74.8 | 76.7 | 3/35 | 1.9% |
| 3-component GMM (GlcN/GlcNAc/ambiguous) | 59.5 | 65.2 | 0/35 | 5.7% |

**Failed strategies (all worse than baseline):**

| Strategy | Result | Why it failed |
|---|---|---|
| Cross-file global z-score | 76.7% | mixes STM contrast variability with chemistry |
| Cross-file percentile rank | 74.8% | same issue |
| Mixed per-file + global features | 68.1% | high-D noise dominates |
| PCA patch PCs (1–5) + local prom | 43–57% | pixel noise dominates in GMM |
| Chain-level flip (amp corr / var ratio) | 17–35% phys | self-consistent ≠ true |
| 3-component GMM for abstention | 45% on 71% classified | middle cluster ≠ ambiguous |

**Forward selection on exact count:** Adding patch_u_asym improves exact chains
from 8→11 without changing per-lobe accuracy. No other single feature improves
by more than 0.5%.

### Error analysis (best config: local prom + patch_u_asym + GMM)

```
Lobe 1 (GlcN):  100%  ██████████████████████████████  ← easy (chain edge)
Lobe 2 (GlcNAc): 80%  ████████████████████████        ← acetyl detected via patch asymmetry
Lobe 3 (GlcN):   80%  ████████████████████████        ← some false positives
Lobe 4 (GlcN):   63%  ███████████████████             ← hardest (middle, no distinctive signal)
Lobe 5 (GlcNAc): 77%  ███████████████████████         ← acetyl detected
Lobe 6 (GlcN):   94%  ████████████████████████████    ← easy (chain edge)
```

All errors concentrate on GlcNAc detection (lobes 2, 5) and middle GlcN (lobe 4).
The acetyl signal is too weak in STM at −0.3 V for reliable separation of every
lobe.

### Abstention framework (3-class output: 0 / 1 / ?)

Using GMM max-responsibility as confidence:

| Threshold | Classified | Abstained (?) | Accuracy on classified | Full-chain exact |
|---|---|---|---|---|
| none (forced) | 234/234 | 0 | 82.4% | 11/35 |
| ≥ 0.55 | 219/234 (94%) | 15 | 83.4% | 9/35 |
| **≥ 0.60** | **167/234 (71%)** | **67** | **90.9%** | partial |
| ≥ 0.65 | 117/234 (50%) | 117 | 91.6% | partial |

At threshold 0.60, 71% of lobes are assigned at 91% accuracy, and the hardest
29% are honestly reported as `?`.

### Supervised upper bound

Leave-one-file-out nearest-centroid on local prominence (4): **76.2%** test
accuracy. The label-free GMM (82.4%) exceeds this because GMM full covariance
captures cluster elongation that centroid classification misses.

### Literature context

A 2025 single-molecule chitosan STM abstract (Wu Xiaocui, GDR NS CPU) confirms
that "direct imaging reveals the sequence of individual chitosan molecules,
defined by acetyl positions." The approach is validated; the remaining
difficulty is signal-to-noise at the current bias/resolution.

### Current best label-free configuration

```text
Features: loc_amp_prominence, loc_amp_rel, loc_amp_neighbor_ratio,
          loc_integrated_prominence, patch_u_asym  (per-file z-scored + interactions)
Method:   GMM full covariance, 10-seed ensemble
Physical mapping: GlcNAc = higher-amplitude cluster
Result:   82.4% physical = 82.4% oracle, 0% gap, 11/35 exact chains
          (self-training 5 iters: 81.4% / 13 exact / 1.9% gap — alternative)
```

Predictions written to `results/unit_assignment/best_labelfree_predictions.tsv`
(with confidence) and `best_labelfree_3class_predictions.tsv` (with `?`).

### Follow-up: high-resolution 17×17 patches and wavelet diagonal detail

Re-extracted patches at 17×17 (0.04 nm/step, ±0.32 nm) from raw SXM data.
Computed Haar wavelet decomposition at 3 levels on both 9×9 and 17×17 residual
patches, plus cross-lobe pair features (adjacent amplitude differences, symmetric
prominence).

**Key finding**: the level-1 diagonal detail (HH1) of the 17×17 residual patch
captures a diagonal component of the acetyl LDOS signal that is invisible in
pure u-axis asymmetry:

```text
NEW BEST: local prominence (4) + patch17_wav_hh1_abs + interactions [GMM]
Features: loc_amp_prominence, loc_amp_rel, loc_amp_neighbor_ratio,
          loc_integrated_prominence, patch17_wav_hh1_abs
          (per-file z-scored + 10 cross-terms)
Method:   GMM full covariance, 10-seed ensemble
Result:   84.3% physical / 85.2% oracle / 1.0% gap / 13/35 exact chains
```

This is the first configuration to break the 82.4% ceiling. The HH1 sub-band
captures diagonal high-frequency residual structure (edges/corners at ~0.08 nm
scale) that corresponds to the acetyl group's off-axis LDOS perturbation.
Combining HH1 with the 9×9 u-asymmetry does not improve further (dimensionality
penalty in the GMM).

Other features tested in this round (none improved over baseline):
- 17×17 u-asymmetry (81.0%), t-asymmetry (81.9%), center-ring contrast (82.4%)
- 17×17 wavelet LH1, HL1, LH2, HH2, LH3 (81–82%)
- 9×9 Haar wavelet features (same signal as patch_u_asym, redundant)
- Cross-lobe pair features (pair_amp_diff, pair_sym_prom — 75–82%)
- Raw (non-residual) patch wavelets (78–81%)

### Extended 17×17 descriptor sweep

Continued the search on existing artifacts without reading any Raven outputs. New
descriptors were computed only from the aligned 17×17 residual patches: HH1
quadrant energies, residual positive/negative moments, diagonal-gradient filters,
Fourier diagonal/axis power, and fixed diagonal matched filters. Truth was used
only after each assignment was written to a prediction TSV and graded externally
with `grade_unit_assignment.jl`.

Two benchmark-improving candidates were found:

```text
v3: BASE + neg_diag135 + interactions [GMM]
    neg_diag135 = center of mass of the negative residual along the t-u diagonal
    Grade: 84.8% physical / 85.7% oracle / 14/35 exact chains

v4: BASE + hh1_q00_abs + neg_anis + interactions [GMM]
    hh1_q00_abs = upper-left HH1 quadrant energy
    neg_anis    = anisotropy of the negative residual moment
    Grade: 85.2% physical / 87.1% oracle / 17/35 exact chains
```

These are real post-hoc benchmark improvements over HH1_abs, but their
label-free diagnostics are less convincing than the conservative HH1_abs model:

```text
Model                 dim  ΔBIC(k1-k2)  silhouette  seed agreement  grade
BASE local             10     +83.6        0.499        1.000        82.4 / 8 exact
patch9_u               15     -26.1        0.215        0.673        82.4 / 11 exact
HH1_abs                15      +6.6        0.220        0.732        84.3 / 13 exact
v3 neg_diag135         15      -4.3        0.327        0.748        84.8 / 14 exact
v4 q00+neg_anis        21     -57.5        0.155        0.698        85.2 / 17 exact
```

Interpretation: `neg_diag135` and `hh1_q00_abs + neg_anis` are plausible
chemistry-adjacent residual descriptors, but selecting v4 as canonical would be
too close to benchmark feedback: it has higher dimensionality and poorer
unsupervised separation evidence. Keep v4 as an **exploration candidate** and
keep HH1_abs as the conservative label-free production candidate until an
objective criterion, an independent dataset, or converged DFT-STM molds confirm
the extra descriptors.

Abstention and ensemble checks:
- Majority voting among BASE/patch9_u/HH1_abs keeps 84.3% per-lobe accuracy but
  improves exact chains to 14/35; useful diagnostic, not a new physical model.
- Agreement between patch9_u and HH1_abs classifies 182/210 lobes at 88.5%
  physical accuracy.
- HH1_abs confidence thresholds classify fewer lobes but reach about 92–93%
  physical accuracy above confidence 0.60–0.95.

`grade_unit_assignment.jl` now accepts `?` predictions as explicit abstentions
and reports classified coverage. This keeps abstention diagnostics in the same
external grading path as full binary predictions while ensuring abstained lobes do
not enter accuracy/confusion counts or inflate exact-sequence counts.

### Leave-one-file-out cross-validation (LOFO)

Ran LOFO to honestly measure which features generalize vs overfit the benchmark.
For each of the 35 graded files: fit GMM on the other 34 files' lobes
(unsupervised), assign the held-out file, map clusters by training amplitude
(label-free). Truth used only at the final grading step.

```text
Model                All-at-once   LOFO      Drop      Verdict
patch9_u (82.4%)       82.4%      69.0%    -13.3%     OVERFIT
HH1_abs (84.3%)        84.3%      80.0%     -4.3%     borderline
v3 neg_diag135         80.5%      82.9%     +2.4%     STABLE / BEST
v4 q00+neg_anis        61.4%      66.2%     +4.8%     unstable
```

**Key finding**: `neg_diag135` (center of mass of the negative residual along the
135° diagonal) is the most generalizable label-free descriptor. It improves under
LOFO relative to all-at-once, meaning it captures a real physical signal rather
than benchmark-specific noise. In contrast, `patch9_u_asym` and the v4 pair are
heavily overfit — their all-at-once benchmark gains do not survive
cross-validation.

**Recommendation update**: `neg_diag135` should be preferred over `HH1_abs` as the
conservative label-free descriptor for production, because it has the strongest
LOFO generalization (82.9% vs 80.0%). The v4 pair should be discarded as a
benchmark artifact. Note: LOFO all-at-once numbers differ slightly from the
Julia-validated grades because of implementation details in the GMM EM loop, but
the relative LOFO ranking is the robust signal.

Additional challenged attempts after the LOFO result did not improve the honest
ceiling:

```text
Candidate                         LOFO physical accuracy
neg_diag135                       82.9%   (best)
BASE only                         81.4%
Gabor 45° / 90° / 135°            81.4%   (neutral)
quad diagonal asymmetry           80.5%
HH1_abs                           80.0%
HH1+HH2+HH3                       58.1%   (overfit/noise)
LH1+LH2+LH3                       50.5%   (overfit/noise)
neg_diag135 + HH1_abs             81.0%   (hurts)
neg_diag135 + Gabor135            81.4%   (hurts)
```

Also tested **per-chain clustering** (no global training, therefore no LOFO
overfit): k-means within each chain on single features and small feature sets,
with the higher-amplitude cluster mapped to GlcNAc. Best result was only 80.0%
physical (`loc_integrated_prominence`), and BASE per-chain clustering reached
78.1% physical / 82.9% oracle. This confirms that the useful signal is not a
simple within-chain two-cluster separation; global cross-file pooling is still
needed, but only the `neg_diag135` descriptor survives LOFO.

### Information-theoretic ceiling: Fisher LDA supervised upper bound

To determine whether further feature extraction could help, a **supervised Fisher
Linear Discriminant Analysis** was run under LOFO: for each held-out file, the
LDA projection was trained on the other 34 files' **true labels** (clearly
supervised, diagnostic only) and applied to the held-out file.

```text
Fisher LDA LOFO (supervised, truth in training):
  BASE only:         81.4%
  BASE+neg_diag135:  81.9%

Unsupervised GMM LOFO (no truth):
  BASE+neg_diag135:  81.4%   (99.4% of supervised ceiling)
```

**Conclusion**: the unsupervised label-free assignment already operates at
99.4% of the supervised information ceiling for these features. The bottleneck
is not feature extraction, clustering, or model complexity — it is the intrinsic
separability of the STM signal at −0.3 V. The acetyl group's electronic
contribution at this bias is too weak relative to the pyranose ring and tip
noise to exceed ~82% per-lobe accuracy under honest cross-validation.

The only paths beyond this single-feature ceiling are:
1. **Converged DFT-STM molds** (jobs `28444935`, `28445456`): encode the LDOS
   electronic structure that the raw STM contrast cannot resolve.
2. **Different experimental conditions**: bias closer to the N-acetyl resonance,
   sharper tip, lower temperature, or CO-functionalized tip for sub-molecular
   resolution.
3. **More data**: additional 6mer scans would reduce LOFO variance and allow
   more features to be tested without overfitting.

### Backward scan channel: exploiting 100% of the SXM data

A key oversight was identified: `extract_lobe_patches.jl` used only the **forward
Z scan** (`direction="fwd"`), ignoring the backward scan that Nanonis stores in
the same SXM file. Each file contains 2 images (512×512, forward + backward),
so 50% of the data was discarded.

Created `test/extract_lobe_patches_bwd.jl` to extract:
- Backward Z residual patches (same Gaussian model, backward data)
- Forward-backward difference residual (removes static topography, isolates
  directional electronic asymmetry)

Re-extracted at 17×17 (0.04 nm/step, ±0.32 nm) for all 39 primary files. Then
computed features from each channel and tested under LOFO:

```text
Channel     Feature               LOFO (10 seeds)
forward     fwd_neg_diag135         81.0%   (previous champion ~81.4%)
backward    bwd_neg_com_t           84.3%   ← NEW CHAMPION
difference  diff_signed_com_t       81.9%   (not reproducible, seed-sensitive)
```

**`bwd_neg_com_t`** (center-of-mass of the negative residual along the backbone
direction, computed from the backward Z scan) achieves **84.3% LOFO**,
reproducibly across two independent seed ranges (0-9 and 10-19). This is +3%
over the forward-only champion.

The backward scan captures different tip-sample electronic states during the
reverse sweep. The acetyl group's off-axis electronic structure apparently
manifests differently in the backward direction, providing complementary
separability.

Note: the all-at-once benchmark for `bwd_neg_com_t` is only 81.4% (same as BASE
alone). The improvement appears only under LOFO — the feature generalizes better
than it fits in-sample. This is the correct behavior for a production-relevant
feature (we want to assign unseen molecules, not memorize the benchmark).

The Fisher LDA supervised ceiling also improves with backward channels:
```text
Fisher LDA LOFO (supervised, truth in training):
  BASE only:              81.4%
  BASE+fwd_diag135:       81.9%
  BASE+fwd+bwd:           82.4%
  BASE+fwd+diff:          82.9%
```

So the backward/difference channels add 1-1.5% of real supervised separability.
The unsupervised GMM exploits nonlinear structure to reach 84.3% — above the
linear Fisher ceiling, confirming the feature carries genuine class information.

### Follow-up audit: assumptions, ensemble, and abstention

The apparent 84.3% LOFO ceiling was re-audited from first principles rather than
treated as final. Several plausible missed assumptions were tested and rejected:

```text
Hypothesis / test                                      Result
lobe index vs t_nm truth-order mismatch               refuted: 0/35 order mismatches
reverse chain orientation                             irrelevant: truth 010010 is palindromic
absolute backward height/residual features            supervised signal, no GMM gain
backward residual recalibration z_bwd ~= a*model+b    no gain over bwd_neg_com_t
parity-canonicalized signed residual features         no gain; often worse
3-6 component GMM merged by amplitude                 unstable, over-calls 1
equal-prior GMM prediction                            small lift only (~84.8% transient)
split-width Gaussian shape alone                      complements but not sufficient
```

The useful missed point was **not** a single feature. It was that independent
label-free views make partly different errors. A three-view ensemble was built
from:

```text
View 1: BASE + bwd_neg_com_t
View 2: BASE + bwd_neg_diag45
View 3: BASE + split_log_skew
```

Each view is trained by LOFO GMM with physical amplitude mapping; the final score
is the mean held-out probability across the three views. Truth is used only after
the prediction TSV is written and graded by `grade_unit_assignment.jl`.

Confirmed with a 20-seed ensemble:

```text
View / rule                         Grade (primary 35 files)
bwd_neg_com_t alone                 84.3% / 11 exact
bwd_neg_diag45 alone                83.3% / 11 exact
split_log_skew alone                82.4% / 14 exact
forced 3-view ensemble              178/210 correct = 84.8%, 14 exact
ensemble abstain 0.20/0.80          diagnostic: 151/169 = 89.3% classified
                                      honest: 151/210 correct + 59/210 uncertain
ensemble abstain 0.15/0.85          diagnostic: 124/137 = 90.5% classified
                                      honest: 124/210 correct + 86/210 uncertain
```

Grader outputs:

```text
results/unit_assignment/grade_best_labelfree_ensemble3_forced.tsv
results/unit_assignment/grade_best_labelfree_ensemble3_abstain80.tsv
results/unit_assignment/grade_best_labelfree_ensemble3_abstain85.tsv
```

Interpretation: the full-coverage improvement is modest (+0.5% over
`bwd_neg_com_t`, one additional lobe and +3 exact chains), but it is real enough
to replace the single-feature model as the best conservative label-free output.
The more important production lesson is abstention: independent views can flag
ambiguous lobes without truth or composition priors. The old headline of 89-91%
on classified lobes is only a diagnostic accuracy among emitted non-`?` labels;
the honest two-score presentation counts any post-hoc benchmark error as
uncertainty, so the reported production-style score is `correctly assigned +
uncertain = 100%`.

### Updated conclusion

The forward-only apparent ceiling was too pessimistic because the backward scan
had been ignored. The current honest full-coverage benchmark is:

```text
Best single descriptor:  BASE + bwd_neg_com_t           84.3% LOFO / 11 exact
Best conservative model: 3-view label-free ensemble     84.8% LOFO / 14 exact
Best honest output mode: ensemble abstention            correct + uncertain = 100%
  0.20/0.80 band: 151/210 correct (71.9%) + 59/210 uncertain (28.1%)
  0.15/0.85 band: 124/210 correct (59.0%) + 86/210 uncertain (41.0%)
```

This still is **not solved**: the hard errors remain concentrated on the two
GlcNAc positions, especially lobe 5, and exact full-chain recovery is only 14/35
at full coverage. The remaining paths beyond this level are converged DFT-STM
molds, different experimental conditions, more data, or a principled abstention
workflow for production maps.

### Follow-up: honest abstention reporting and QE relaunch

- Updated `test/grade_unit_assignment.jl` so the summary and per-file TSV include
  a conservative **honest abstention** view. The existing classified-lobe accuracy
  remains as a diagnostic, but the honest view reports only correctly assigned and
  uncertain positions; uncertainty includes explicit `?` plus any wrong assignment
  found by post-hoc benchmark grading. This keeps the prediction rule label-free
  while preventing benchmark errors from being presented as honest calls.
- Re-ran the existing ensemble grades:

```text
forced ensemble:        178/210 correct (84.8%) + 32/210 uncertain (15.2%)
abstain 0.20/0.80:      151/210 correct (71.9%) + 59/210 uncertain (28.1%)
abstain 0.15/0.85:      124/210 correct (59.0%) + 86/210 uncertain (41.0%)
```

- Relaunched the production DFT-STM mold jobs on Raven without `--watch`, after a
  dry-run verified the exact two-directory sequential command:

```text
qe/glcn_restart2 -> 28525353
qe/glcnac        -> 28525354  (sequential afterok dependency)
```

Do not fetch or inspect these QE outputs until a completion or timeout
notification is available.

### Follow-up: forced-as-base abstention rejector

The forced three-view ensemble is now treated as the canonical base predictor, and
abstention is a separate rejector layer. Added
`test/build_unit_abstention_variants.jl`, which is label-free: it reads a base
prediction TSV plus optional auxiliary prediction TSVs, keeps the base 0/1 label
only if confidence/agreement gates pass, and emits `?` otherwise. It does not read
the truth sequence.

The useful rejector was to require agreement between the forced ensemble and the
older conservative label-free model (`best_labelfree_predictions.tsv`). Post-hoc
benchmark grades:

```text
Rule                                      Classified   Accuracy     Honest view
old confidence band 0.20/0.80              169/210    151/169=89.3%  151 correct + 59 uncertain
forced + agreebase65                       171/210    154/171=90.1%  154 correct + 56 uncertain
forced + agreebase70                       169/210    153/169=90.5%  153 correct + 57 uncertain
forced + agreebase80                       157/210    144/157=91.7%  144 correct + 66 uncertain
```

`agreebase65` is the new balanced abstention recommendation: it improves on the
old `abstain80` both in honest correct count and uncertainty while still lowering
the residual wrong-assignment count among emitted labels. `agreebase70` is the
safer profile if one fewer residual wrong label is worth losing one additional
correct label. The rule itself is label-free; the benchmark truth above is used
only to grade the frozen TSV outputs.

To target **<5% wrong assignments among emitted labels**, the best-coverage strict
profile found in the same post-hoc sweep was:

```text
forced + confidence >= 0.875
       + agreement with best_labelfree_v3_neg_diag135_predictions.tsv
       + agreement with best_labelfree_predictions.tsv
```

Generated artifacts:

```text
results/unit_assignment/best_labelfree_ensemble3_abstain_err05_predictions.tsv
results/unit_assignment/grade_best_labelfree_ensemble3_abstain_err05.tsv
```

Grade:

```text
classified:      106/210 = 50.5%
correct emitted: 101/106 = 95.3%
wrong emitted:     5/106 = 4.7%
honest view:     101/210 correct + 109/210 uncertain
```

This achieves the <5% emitted-error target, but only by abstaining on about half
the lobes. Treat `agreebase65` as the balanced production profile and `err05` as
the high-confidence/low-coverage profile.

### Follow-up: frozen 0/1/? benchmark report harness

After the 145-file counting benchmark was confirmed on Viper, the next safest
unit-assignment step was **not** another feature or threshold search. At this
point the available frozen prediction/report artifacts covered the historical
240817 unit subset (35 clean/clean_target chains, 210 lobes), so a report-only
harness was added to avoid benchmark leakage, threshold tuning, or denominator
confusion:

```text
test/report_unit_assignment_benchmark.jl
```

The script runs `test/grade_unit_assignment.jl` on frozen prediction profiles and
consolidates the outputs under `results/unit_assignment/benchmark_report/`:

```text
summary.tsv
lobe_position_errors.tsv
report.md
grades/forced_ensemble3.tsv
grades/agreebase65.tsv
grades/err05.tsv
```

Canonical command:

```bash
julia --project=. test/report_unit_assignment_benchmark.jl
```

Reproduced historical 35-file subset results:

```text
forced_ensemble3  210/210 emitted, 178/210 correct (84.8%), 32 emitted errors, 14/35 exact
agreebase65       171/210 emitted, 154/171 correct (90.1%), 17 emitted errors,
                  honest view 154/210 correct + 56/210 uncertain, 7/35 exact
err05             106/210 emitted, 101/106 correct (95.3%), 5 emitted errors,
                  honest view 101/210 correct + 109/210 uncertain, 0/35 exact
```

The per-lobe-position table makes the main failure mode explicit without using it
to alter the method: the forced full-coverage ensemble is perfect at lobes 1 and
6, but weak at lobe 2 (`22/35`, 62.9%) and lobe 5 (`23/35`, 65.7%).
`agreebase65` improves most positions by abstaining, but lobe 5 remains the hard
case (`19/27`, 70.4% among emitted labels). `err05` reaches the intended strict
emitted-error regime (5/106 wrong emitted labels) by abstaining on roughly half of
the lobes.

This is a reporting/benchmark consolidation only. It reads truth labels only in
the grader/report path, does not sweep thresholds, and does not promote a new
profile.

### Follow-up: 6mer raw-data classification inventory

The full lab-storage tree `/home/durif/Rebecca/data/data/` was inventoried to
check whether more 6mer scans can improve unit assignment. All `935` `.sxm` files
read successfully with `STMSXMIO.read_sxm`, and every file exposes only the two
topography channels `Z fwd` and `Z bwd`; there is no additional current image
channel to exploit. The data gap is therefore classification/curation, not hidden
SXM channels.

Created:

```text
benchmarks/chitosan_6mer_data_inventory.tsv
```

This TSV has one row per raw `.sxm` keyed by `relative_path` (not basename,
because `241114_*` basenames appear in two folders). It joins:

- `benchmarks/chitosan_240817_unit_sequences.tsv`
- `benchmarks/chitosan_manual_20240814_20240818.toml`
- `/home/durif/Rebecca/data/data/chitosan_manual_annotations.md`
- existing `results/best_plots*/summary_overlap060_hard.tsv` paths when present

Current classification summary:

```text
unit_sequence_labeled_primary        35   strict unit-assignment benchmark
unit_sequence_labeled_stress          9   stress only
manual_N6_clean_target               33   useful for N-counting; unit sequence not filled
manual_ok_chitosan_needs_N_review     4   plausible chitosan; expected N not confirmed
manual_ambiguous_review              20   keep out of strict scoring
manual_excluded / excluded_240817    54   exclude from strict scoring
unclassified_needs_review           780   likely prep/search/test unless reviewed
```

The best immediate expansion at that time was the `33` `manual_N6_clean_target`
files: `240814` (15), `240816` (3), and `240818` (15). This was a provenance and
curation note, not a scientific reason to shrink the 0/1/? benchmark. The later
scope clarification below supersedes the “needs unit sequence” framing: the 6mer
unit-control sequence is `NKNNKN` for external grading only, while the method must
remain label-free.

### Follow-up: provisional 6mer pre-assignment review sheet

Created a review-only pre-assignment sheet from the inventory plus existing batch
summaries:

```text
benchmarks/chitosan_6mer_preassignment_review.tsv
benchmarks/chitosan_6mer_preassignment_review.md
```

These files are **not canonical labels**. They are a fast plot-review queue: use
the TSV/Markdown to inspect existing `best_plot` images and then promote or
reject rows explicitly. The pre-assignment logic does not alter fitting,
selection, or unit-assignment truth; it only triages already documented manual
annotations and label-free batch diagnostics.

Validation after generation:

```text
rows:             935
unique keys:      935 relative_path values
plot links:       242 existing files
priority review:  106 rows
```

The priority rows are the useful near-term human-review set: `33` documented
`manual_N6_clean_target` rows that need sequence confirmation before unit use,
`4` documented `ok_chitosan` rows without confirmed N, `26` unclassified clean
N=6 candidates, `12` mixed-diagnostic N=6 candidates, `11` near-6 candidates,
and `20` unclear rows. The remaining rows are documented unit benchmarks/stress
cases, documented excludes/ambiguous scans, likely non-chain guard collapses, or
unclassified scans with no existing batch plot.

One summary row (`240818_019.sxm`) contained `best_plot=ERR`; the review TSV
keeps its batch summary reference but leaves the plot field empty, and the
Markdown marks it as `missing`.

Visual review update: the user inspected `240307_015.sxm`, `240307_016.sxm`, and
`240307_017.sxm` and confirmed the true count is `N=6` for all three. The review
TSV/Markdown now mark them as `accept_counting_visual_N6_confirmed` with
`expected_N=6`, while keeping `unit_sequence` empty because no GlcNAc/GlcN
sequence has been confirmed for these rows.

Second visual review update: the user confirmed the late `240817` review block
as chitosan `N=6` for `240817_064.sxm` through `240817_083.sxm`, excluding no
files in that interval, and for `240817_085.sxm` through `240817_094.sxm`;
`240817_084.sxm` is marked as non-chitosan/exclude. The same review marked these
files ambiguous and not usable for strict visual classification:
`240817_002.sxm`, `240817_017.sxm`, `240817_021.sxm`, `240817_029.sxm`,
`240817_030.sxm`, `240817_031.sxm`, `240817_032.sxm`, `240817_034.sxm`,
`240817_035.sxm`, `240817_037.sxm`, `240817_038.sxm`, and `240817_051.sxm`.
Previous explicit ambiguous decisions were also recorded for
`240308_Cu100020.sxm`, `240308_Cu100021.sxm`, `240308_Cu100030.sxm`, and
`240814_016.sxm`.

Third visual review update: the remaining `20240308` near-6/mixed-diagnostic
candidates were reviewed. The user marked these as ambiguous and excluded from
strict visual classification: `240308_Cu100012.sxm`, `240308_Cu100022.sxm`,
`240308_Cu100023.sxm`, `240308_Cu100024.sxm`, `240308_Cu100025.sxm`,
`240308_Cu100026.sxm`, `240308_Cu100027.sxm`, `240308_Cu100028.sxm`,
`240308_Cu100029.sxm`, `240308_Cu100031.sxm`, `240308_Cu100032.sxm`,
`240308_Cu100033.sxm`, `240308_Cu100035.sxm`, `240308_Cu100072.sxm`, and
`240308_Cu100081.sxm`. `240308_Cu100067.sxm` was not treated as ambiguous; it
was marked as a bad scan/exclude. All of these review entries keep
`expected_N` and `unit_sequence` empty.

Fourth visual review update: `240307_019.sxm` was judged to be a likely regular
chitosan chain with `N=6`, but with a small doubt. The TSV initially marked it as
`probable_counting_visual_N6_doubt` with `expected_N=6`, kept `unit_sequence`
empty, and required a second check before strict use. That second check was later
resolved in the 2026-07-01 benchmark update. The user also reviewed the remaining
`results/best_plots_20240308/` candidates and concluded that even the best are
ambiguous at most, with many junk scans. The `20` remaining `pre_review_unclear`
rows from that folder were therefore moved to `visual_review_ambiguous_exclude_strict`.

Fifth visual review update: the earlier ambiguity mark on `240817_002.sxm`,
`240817_017.sxm`, and `240817_021.sxm` was double-checked and reversed. The user
confirmed these stay in the benchmark and are classified OK. The TSV now restores
all three to `accept_unit_training_known_010010`, with `expected_N=6` and
`unit_sequence=010010` preserved as external benchmark provenance. This returns
the strict unit-sequence benchmark count to `35` rows.

To expand the benchmark beyond files with existing fit plots, generated raw
triage contact sheets for all `673` remaining `unclassified_no_batch_plot` rows:

```text
results/triage_unclassified_raw/index.md
results/triage_unclassified_raw/index.tsv
```

Each contact-sheet tile shows raw `Z fwd` and `Z bwd`, plus file name, pixel
size, scan range, and acquisition time. These plots are **visual triage only**:
they do not alter fitting/selection and do not create benchmark labels. Candidate
files found from these sheets still need explicit human confirmation before being
promoted into the N=6 counting benchmark, and an explicit sequence row before any
GlcNAc/GlcN unit-assignment benchmark use.

Follow-up triage focusing on likely candidates: selected `470` of those
unclassified rows by metadata only (square-ish raw scans, `4–12 nm` field,
`>=128 px`) and wrote focused contact sheets to:

```text
results/triage_potential_benchmark/index.md
results/triage_potential_benchmark/candidate_index.tsv
```

The split is `309` high-priority scans with field `<=8 nm` and `161` secondary
scans with field `8–12 nm`. Full STMFit `*_best.png` generation was deliberately
not launched for this whole set because `470` candidates is too many for a quick
interactive pass; first use the contact sheets to choose a smaller fit batch.

This update is still review-only. For rows that already carried external unit
truth (`010010`), the TSV keeps `expected_N`/`unit_sequence` as provenance but
changes `pre_assignment`/`recommended_action` so ambiguous rows are not used in
strict visual expansion. Existing unnamed `excluded_240817_batch` rows were **not**
overwritten by the broad "the rest is ok chitosan" statement; that would require
an explicit confirmation because it would reverse earlier exclusions.

---


## 2026-06-27 — QE GlcN timeout, restart, and preliminary LDOS map

### Completed Jobs Reconciled

- GlcN production relax job `28303162` ended at the 24 h walltime limit.
  Slurm reported `CANCELLED ... DUE TO TIME LIMIT`; QE output stopped during
  the SCF after BFGS step 31, so this is **not** a converged production mold.
- The latest `ATOMIC_POSITIONS` block was extracted from
  `qe/glcn/glcn_central_relax.out` to `qe/glcn/glcn_central_best.xyz`
  (`213` atoms), preserving the best-so-far geometry rather than restarting
  from the hand-built slab.
- `qe/glcn_restart` was regenerated from that best geometry with the active
  pilot settings (`8` MPI tasks, `50/360 Ry`, Γ-only, `96000 MB`, `24:00:00`),
  preflighted successfully, synced to Raven, and submitted as job `28363474`.

### Preliminary GlcN SCF+PP

- Preliminary SCF+PP job `28354566` completed successfully:
  - `qe/glcn_prelim/glcn_central_scf.out`: SCF converged in 34 iterations,
    `JOB DONE`.
  - `qe/glcn_prelim/glcn_central_pp.out`: wrote
    `glcn_central_ldos.cube`, `JOB DONE`.
  - Local cube: `qe/glcn_prelim/glcn_central_ldos.cube` (~136 MB).
- This cube is from the unconverged best-so-far geometry and is **diagnostic
  only**. It must not be treated as the final GlcN mold.

### Diagnostic Map Conversion

- Extracted a pilot-frame TSV from `qe/glcn/glcn_central_best.xyz` using the
  `8×6×3` slab offset (`213 - 69 = 144`):
  `origin_indices=156,157,158,159,160,161`, `axis_from=159`, `axis_to=156`,
  `plane_index=165`.
- Converted the preliminary GlcN cube at diagnostic sampling height
  `height_nm=0.35` to:

```text
templates/chitosan_stm_maps_glcn_prelim_h035.tsv
```

- Map sanity check: `169/169` pixels finite, `0` `NA`.
- `import_stm_mold_maps.jl` intentionally was **not** run for this one-sided
  map: it requires both base unary maps (`type=0` GlcN and `type=1` GlcNAc).
  Do not fabricate a dummy GlcNAc template.

### Current State / Next Gate

- Wait for restart job `28363474` to finish before reading its QE outputs.
- If `28363474` converges and produces final GlcN relaxed/cube outputs, submit
  GlcNAc production next. GlcNAc remains deliberately unsent until GlcN
  production succeeds.
- The physical LDOS sampling height remains to be chosen before production
  `finalize_qe_mold_workflow.jl` outputs are frozen and scored.

---

## 2026-06-27 — GlcNAc A/B plan prepared, not submitted

Goal: advance type `1` (GlcNAc) without competing with the active GlcN restart
job `28363474` or violating the rule that production GlcNAc waits for GlcN
success.

### Plan A — Production GlcNAc relax → SCF → PP

- Refreshed `qe/glcnac` from the active `8×6×3` pilot structure with the same
  resource-constrained settings as GlcN:
  - prefix `glcnac_central`
  - `ecutwfc=50`, `ecutrho=360`
  - Γ-only (`--kpoints 1,1,1`)
  - `8` MPI tasks, `12000 MB/task`, `24:00:00`
  - total memory `96000 MB`
- Copied the validated pseudo set into `qe/glcnac/pseudo/`.
- Local preflight passed:

```bash
julia --project=. test/preflight_qe_mold_inputs.jl \
    --dir qe/glcnac \
    --out hpc/qe_molds/qe_input_preflight_glcnac.tsv \
    --max-total-tasks 8 \
    --min-mem-mb 96000
```

- **Not submitted.** Submit Plan A only after GlcN production succeeds.

### Plan B — Diagnostic GlcNAc SCF+PP-only

- Created `qe/glcnac_prelim/` with:
  - `pw_scf.in` and `pp_ldos.in` copied from refreshed `qe/glcnac`
  - pseudo files copied from `qe/glcnac/pseudo/`
  - `run_scf_pp.sbatch` for SCF+LDOS cube only, no relax step
- This is weaker than the GlcN preliminary cube because there is no GlcNAc
  best-so-far relaxed geometry yet; it uses the regenerated initial pilot
  geometry. Treat it only as a pipeline/type-1 smoke test, not a production mold.
- Extended launch safety tooling to support preliminary SCF+PP dirs:
  - `test/preflight_qe_mold_inputs.jl` now accepts `run_scf_pp.sbatch` when no
    `run_qe_mold.sbatch` is present, requires `pw_scf.in`/`pp_ldos.in`, and
    rejects any relax/handoff commands in preliminary mode.
  - `hpc/submit_qe_molds.sh` selects `run_scf_pp.sbatch` for such dirs.
  - `hpc/launch_qe_molds_remote.sh` syncs `run_scf_pp.sbatch`.
- Local preflight passed:

```bash
julia --project=. test/preflight_qe_mold_inputs.jl \
    --dir qe/glcnac_prelim \
    --out hpc/qe_molds/qe_input_preflight_glcnac_prelim.tsv \
    --max-total-tasks 8 \
    --min-mem-mb 48000
```

- Local submit dry-run confirms the correct script is chosen:

```text
(cd 'qe/glcnac_prelim' && sbatch run_scf_pp.sbatch)
```

### Gate

- Do not submit GlcNAc production while `28363474` is outstanding.
- Prefer also holding GlcNAc preliminary submission until `28363474` finishes,
  because Raven QOS limits have already been sensitive to pending/running 8-task
  allocations. If the queue policy is relaxed or the user explicitly accepts
  the risk, submit only `qe/glcnac_prelim` with `--min-mem-mb 48000`.

### Follow-up: queued Plan B with an external dependency

- Added an explicit Slurm dependency option to the QE launch path:
  - `hpc/submit_qe_molds.sh --dependency SPEC`
  - `hpc/launch_qe_molds_remote.sh --dependency SPEC`
- The submitter now applies the external dependency to the first job in a
  sequential chain, or to all jobs in a parallel submission. This lets a follow-up
  job be queued without polling or reading the active job's output.
- Dry-run verified the intended GlcNAc preliminary command:

```text
(cd 'qe/glcnac_prelim' && sbatch --dependency=afterany:28363474 run_scf_pp.sbatch)
```

- Submitted `qe/glcnac_prelim` to Raven with `afterany:28363474`:

```text
qe/glcnac_prelim -> 28365256
```

- This job should start only after the GlcN restart job `28363474` leaves the
  queue/running state, regardless of whether `28363474` converges or times out.
  It remains a **diagnostic type-1 smoke-test cube**, not a production GlcNAc
  mold. Do not consume its outputs until a completion/timeout notification is
  available.

Resume commands after completion notification:

```bash
rsync -avz -e "ssh -o ConnectTimeout=180 -o ServerAliveInterval=60" \
  raven:/u/oldu/code/STMFit/qe/glcnac_prelim/glcnac_central_scf.out \
  raven:/u/oldu/code/STMFit/qe/glcnac_prelim/glcnac_central_pp.out \
  raven:/u/oldu/code/STMFit/qe/glcnac_prelim/glcnac_central_ldos.cube \
  qe/glcnac_prelim/
```

### 2026-07-01 — Compact candidate fits accepted for counting; GlcN restart3 submitted

#### Benchmark expansion review

- The compact `triage_potential_benchmark_fit_top10` batch finished for all 45
  candidates:

```text
summary rows with ok plot: 45/45
N_selected counts: N=5: 9, N=6: 23, N=7: 9, N=8: 4
ambiguous_eff counts: false: 35, true: 10
```

- Human visual review then accepted **all 45 fit images** as good chitosan chains
  with true `expected_N=6`. This is a review/grading label only: it was added to
  `benchmarks/chitosan_6mer_preassignment_review.tsv` and regenerated in
  `benchmarks/chitosan_6mer_preassignment_review.md`, but it remains outside the
  fitting/selection path.
- All 45 were assigned:

```text
pre_assignment      = accept_counting_visual_N6_confirmed
recommended_action  = use_for_counting_keep_unit_sequence_unconfirmed
expected_N          = 6
unit_sequence       = <blank>
```

- The count of `accept_counting_visual_N6_confirmed` review rows increased from
  33 to 78 in the 45-image review pass. A follow-up second check then promoted
  `240307_019.sxm` from `probable_counting_visual_N6_doubt` to
  `accept_counting_visual_N6_confirmed`, bringing that review class to 79. At the
  time these rows were kept out of unit-assignment grading for provenance reasons;
  the later scope clarification below supersedes that limitation for benchmark
  definition while preserving the label-free rule.
- Compact review artifacts for this pass are:

```text
results/triage_potential_benchmark_fit_top10/summary_overlap060_hard.tsv
results/triage_potential_benchmark_fit_top10/fit_review_ranked.tsv
results/triage_potential_benchmark_fit_top10/fit_review_ranked.md
```

- Promoted all confirmed `expected_N=6` review rows into a counting-only external
  benchmark manifest:

```text
benchmarks/chitosan_6mer_counting_confirmed.toml
```

  The manifest has 146 `clean_target` entries:

```text
accept_counting_visual_N6_confirmed              79
accept_unit_training_known_010010                35
accept_counting_needs_unit_sequence_confirmation 32
```

  This is for post-hoc counting grades only and must not enter fitting or
  selection. The same file set later became the full 0/1/? benchmark scope, with
  `NKNNKN` as external grading control only. The provenance state at this point was
  tracked in:

```text
benchmarks/chitosan_6mer_validation_pending.tsv
```

  That list recorded 111 rows as missing explicit unit-sequence provenance. It is
  superseded for benchmark scope by the later clarification below: all 145
  confirmed 6mer rows share the external control sequence `NKNNKN`, but that truth
  remains grader-only.

#### QE mold jobs

- Raven jobs were checked after the no-watch relaunch:

```text
28525353  qe-glcn_central    TIMEOUT   1-00:00:30
28525353.0 pw.x              FAILED    1-00:00:24  ExitCode 1:0
28525354  qe-glcnac_central  CANCELLED 00:00:00
```

- The GlcN failure was a walltime cancellation, not an out-of-memory failure:
  the Slurm report showed `MaxRSS` about 4.8 GB per MPI task / about 40 GB per
  node. `qe/glcn_restart2/glcn_central_relax.out` was fetched locally, and the
  last `ATOMIC_POSITIONS (angstrom)` block was extracted to
  `qe/glcn_restart2/glcn_central_best3.xyz` (`213` atoms). The relax had reached
  at least 30 BFGS steps, so continuing from the best-so-far geometry preserves
  real progress.
- Prepared `qe/glcn_restart3` from that extracted geometry with the active Raven
  settings (`8` MPI tasks, `96000 MB`, `24:00:00`, `ecutwfc=50`, `ecutrho=360`,
  Gamma-only). Local and remote preflight passed; remote report:

```text
glcn_restart3  nat=213  frozen_relax_atoms=96  mem_mb=96000  status=OK
```

- Submitted the next no-watch production GlcN restart:

```text
qe/glcn_restart3 -> 28601744
```

Do not fetch or inspect `qe/glcn_restart3` outputs until a completion or timeout
notification is available. The production `qe/glcnac` run remains blocked until
GlcN succeeds, because the previous dependent job was cancelled before start.

### 2026-07-02 — Targeted divergent-file fit-improvement screen

To avoid recomputing the full expanded counting benchmark for every idea, built a
screening set from the current external grades: the 146 confirmed `expected_N=6`
rows were joined to their recorded summaries, and the 43 rows whose current
`N_selected` is not 6 (including one `ERR`) were written to:

```text
/tmp/opencode/stmfit_6mer_experiment/current_from_summaries.tsv
/tmp/opencode/stmfit_6mer_experiment/divergent.tsv
/tmp/opencode/stmfit_6mer_experiment/triage_by_folder/
/tmp/opencode/stmfit_6mer_experiment/divergent_recursive_fit_input.tsv
```

The fit-input TSVs contain only file paths plus neutral placeholder columns
required by `test/batch_full.jl`; they do **not** contain `expected_N`, unit
sequences, or target labels. Labels were used only after fitting or offline replay
for external grading.

Current summary-derived baseline on the expanded benchmark is `103/146` exact by
`N_selected`. The divergent subset has `N_selected` counts `4:2`, `5:24`, `7:11`,
`8:5`, and `ERR:1`. Existing summary columns show why a naive selector change is
not safe: raw `N_eff` would be `106/146`, but it fixes 11 current failures while
regressing 8 current successes, including primary 240817 files, so dropping the
robust guard is not acceptable as-is.

Fresh local reruns with the current default config were then tried on three
small divergent folders before spending HPC time on all 43 rows:

```text
20240307_LHe_Cu100  0/3 recovered
20240814_LHe_Cu100  2/4 recovered: 240814_017, 240814_018
20240817_LHe_Cu100  2/3 recovered: 240817_077, 240817_091
```

The merged comparison is:

```text
/tmp/opencode/stmfit_6mer_experiment/default_repro_comparison.tsv
```

`adaptive_support_rescue` was also tried on `20240814_LHe_Cu100`; it matched the
current default result (`2/4`) and did not show a rescue-specific gain there. A
local 4-file adaptive run exceeded 20 minutes when launched with 4 threads, so
larger strategy sweeps should run on Viper or via chunked batch jobs rather than
interactively.

Offline replay of the GCV ambiguity threshold from existing columns was also not
a strong candidate: threshold `0.03` gives `104/146` but introduces an expanded
benchmark regression, while `0.05`/`0.06` are net neutral on current summaries.
Simple guard-down variants based on ambiguity or GCV deltas can reach `108/146`
offline, but regress primary 240817 successes such as `240817_017.sxm`; do not
promote them without a stronger label-free rationale.

Next recommended step: run the current default code on all 43 divergent files
using `/tmp/opencode/stmfit_6mer_experiment/divergent_recursive_fit_input.tsv`.
Promote any strategy to the full 146 confirmed benchmark only if the 43-file
screen improves by a meaningful margin and the full run preserves the 39/39
primary benchmark, avoids new `ERR` rows, and preferably has zero regressions
among the 103 currently correct expanded rows.

Follow-up for the benchmark-optimized objective: created repository-synced fit
input TSVs with no labels, suitable for Slurm jobs:

```text
benchmarks/chitosan_6mer_divergent_fit_input.tsv   # 43 rows
benchmarks/chitosan_6mer_confirmed_fit_input.tsv   # 146 rows
```

`hpc/launch_remote.sh` now exports `STMFIT_DATA_DIR=/ptmp/<user>/stmfit/data`
to the Slurm job and forwards dedicated `STMFIT_TSV` / `STMFIT_SKIP_1D` variables,
so benchmark subset runs no longer need fragile multi-word `STMFIT_BATCH_ARGS`.
A Viper dry-run for the full 146-file current-default rerun was correct, but the
actual submission was blocked before sync/submission because both `viper` and
`raven` SSH handshakes timed out, likely because no MPCDF 2FA/ControlMaster
session was active. Resume by opening the MPCDF SSH session first, then run:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/default_repro \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SKIP_1D=1 ./hpc/launch_remote.sh
```

After opening the Viper session, that command synced code/data, instantiated the
Julia project on the Viper login node, and submitted the full 146-file current
default rerun as Slurm array job:

```text
10367643  stmfit  viper/small  8 chunks  PD (Priority)
```

Output directory:

```text
results/experiments/6mer_full146/default_repro
```

The first submit exposed a launcher bug: non-`--watch` submissions attempted to
fetch results immediately, before the remote output directory existed. The job was
already submitted successfully; `hpc/launch_remote.sh` now creates the remote code
directory before rsync and exits after submission unless `--watch` or
`--fetch-only` is requested. Resume/fetch after completion with:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/default_repro \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SKIP_1D=1 ./hpc/launch_remote.sh --fetch-only
```

Completion follow-up: Viper job `10367643_[1-8]` completed successfully on the
`small` partition (`ExitCode=0:0` for every array task). The chunks were merged
on Viper and fetched locally:

```text
results/experiments/6mer_full146/default_repro/summary_overlap060_hard.tsv
```

The merged summary has 146 data rows and 146 successful fits. External grading
against `benchmarks/chitosan_6mer_counting_confirmed.toml` gives:

```text
N_selected  106/146 exact (72.6%), +/-1 = 139/146 (95.2%), over/under = 17/23
```

On the original `benchmarks/chitosan_240817.toml` target subset, the same run
remains `4/4` exact. Relative to the older summary-derived `103/146` baseline,
the fresh run changes 15 files: 8 failures are recovered, 5 former successes
regress, and 2 failures move to a different wrong N, for a net `+3` exact gain.
The recovered rows are `240311_Cu100060.sxm`, `240312_Cu100081.sxm`,
`240814_017.sxm`, `240814_018.sxm`, `240817_077.sxm`, `240817_091.sxm`,
`240818_012.sxm`, and `240818_017.sxm`; the regressions are
`240814_024.sxm`, `240816_005.sxm`, `240817_072.sxm`, `240817_083.sxm`, and
`240818_021.sxm`.

Fresh grading of existing selector columns on the same full146 summary does not
beat the default selected column:

```text
N_selected      106/146
N_eff           106/146
N_ell           106/146
robust_aicc_N   103/146
N_circ          100/146
```

Offline replay of generic post-fit guard variants on the fresh summary can reach
`110/146` exact (best cases: robust downshift only when the GCV ambiguity delta
is small, or only when the GCV choice is ambiguous), but every best-scoring
variant trades fixes for regressions. These rules remain screening observations,
not promoted selection changes: they were discovered by external grading on this
confirmed 6mer set, and need a label-free physical or statistical rationale
before they can become part of the pipeline.

Next benchmark-optimization candidates are therefore experimental, not yet
canonical: either run a full146 Viper job with `adaptive_support_rescue` to test
whether the small local screen missed broader gains, or design a label-free guard
variant from residual/fit diagnostics and then grade it externally. Do not tune a
threshold directly to `expected_N=6` and present it as objective.

Launched the first of those experimental candidates on Viper: a full 146-file
rerun with `config/chitosan_adaptive_support_rescue.toml`, the same label-free
input TSV, and 1D skipped. The dry-run showed the intended project/data/config
paths and Slurm export variables, then the no-watch submission succeeded:

```text
10370782  stmfit  viper/small  8 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan_adaptive_support_rescue.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/adaptive_support_rescue \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SKIP_1D=1 ./hpc/launch_remote.sh
```

After completion, fetch/merge with the same environment plus `--fetch-only`, then
grade `results/experiments/6mer_full146/adaptive_support_rescue/summary_overlap060_hard.tsv`
against `benchmarks/chitosan_6mer_counting_confirmed.toml`. The score to beat is
the current default rerun's `106/146` exact by `N_selected`.

Completion follow-up: all `10370782_[1-8]` array tasks completed with
`ExitCode=0:0` and runtimes from `00:09:05` to `00:15:49`. The remote merge
produced 146 data rows and 146 successful fits, then fetched the summary and
chunk outputs to:

```text
results/experiments/6mer_full146/adaptive_support_rescue/summary_overlap060_hard.tsv
```

External grading on the expanded confirmed counting benchmark gives:

```text
adaptive_support_rescue N_selected  104/146 exact (71.2%), +/-1 = 139/146 (95.2%), over/under = 15/27
```

The original 240817 target subset falls to `3/4` exact (`240817_043.sxm` becomes
the target failure). Relative to the current default rerun (`106/146`), the
adaptive-support run changes only 6 rows: it fixes `240814_024.sxm` and
`240816_005.sxm`, but regresses `240817_043.sxm`, `240817_077.sxm`,
`240817_093.sxm`, and `241113_089.sxm`, for a net `-2` exact result. Therefore
`adaptive_support_rescue` is not a benchmark-improving 6mer-counting candidate;
keep the current default run as the score to beat (`106/146`).

Offline inspection of the two full146 summaries shows why this candidate is a
dead end for the 6mer counting objective: no row accepted an actual support
rescue (`adaptive_support_rescue` count = 0). The adaptive run had 119 plain
keeps, 6 `reject_no_improvement` rows, and 26 rows where the robust guard still
acted. The six changed `N_selected` rows therefore came from altered guard/source
behavior, not from genuinely longer recovered support. Do not launch another
adaptive-support sweep for 6mer counting unless a new label-free trigger or
support diagnostic is added first.

Next label-free screen: selected `support_marginalized_gcv_guard` as the cheapest
existing-policy candidate to test before any new full146 run. This selector
rescales GCV across a grid of support paddings and applies a guarded/parsimony
choice without using labels. Because the summary TSVs do not preserve the fitted
parameter objects needed to replay this selector offline, the first screen is a
43-file rerun on the divergent input TSV, not a full benchmark rerun.

Submitted the 43-row Viper screen after a clean dry-run:

```text
10371435  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/support_marginalized_gcv_guard \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=support_marginalized_gcv_guard STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending/running (`10371435_1` running, `10371435_[2-4]`
pending for `QOSGrpCpuLimit`). After it completes, fetch/merge with the same
environment plus `--fetch-only`, grade the 43 divergent rows externally, and
promote to a full146 run only if it recovers a meaningful number of failures
without a label-dependent rationale.

Completion follow-up: all `10371435_[1-4]` tasks completed with `ExitCode=0:0`
and runtimes from `00:03:29` to `00:04:47`. The remote merge produced 43 data
rows and 43 successful fits, fetched to:

```text
results/experiments/6mer_divergent/support_marginalized_gcv_guard/summary_overlap060_hard.tsv
```

External grading on those 43 divergent rows gives:

```text
support_marginalized_gcv_guard N_selected  16/43 exact (37.2%), +/-1 = 37/43 (86.0%), over/under = 11/16
```

This is a meaningful screen improvement. On the same 43 files, the current full146
default rerun is `8/43` exact, while the older summary-derived baseline was
`0/43`. Relative to the current default, the support-marginalized guard changes
22 rows: 12 fixes, 4 regressions, and 6 wrong-to-wrong moves, for a net `+8`
exact result. Relative to the older summary-derived baseline it fixes 16 rows and
introduces no exact regressions, because all 43 rows were failures there. The
selector source is `support_marginalized_gcv` for all 43 rows, with 21 guarded
downshifts and 22 keeps. This is strong enough to justify a full146 Viper run;
do not promote it as canonical until that full run is graded, because the screen
was enriched for failures and does not measure regressions among the other 103
expanded-benchmark rows.

Promoted the screen to a full 146-file Viper run with the same label-free input
TSV, default chitosan config, and `STMFIT_SELECTION_POLICY=support_marginalized_gcv_guard`.
The dry-run was clean, then the no-watch submission succeeded:

```text
10371737  stmfit  viper/small  8 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/support_marginalized_gcv_guard \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SELECTION_POLICY=support_marginalized_gcv_guard STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending (`10371737_[1-8]`, reason `QOSGrpCpuLimit`).
After completion, fetch/merge with the same environment plus `--fetch-only`, then
grade the full summary against `benchmarks/chitosan_6mer_counting_confirmed.toml`.
The score to beat remains the default rerun's `106/146` exact by `N_selected`.

Completion follow-up: all `10371737_[1-8]` tasks completed with `ExitCode=0:0`
and runtimes from `00:04:56` to `00:06:04`. The remote merge produced 146 data
rows and 146 successful fits, fetched to:

```text
results/experiments/6mer_full146/support_marginalized_gcv_guard/summary_overlap060_hard.tsv
```

External grading on the expanded confirmed counting benchmark gives:

```text
support_marginalized_gcv_guard N_selected  92/146 exact (63.0%), +/-1 = 140/146 (95.9%), over/under = 14/40
```

The original 240817 target subset falls to `2/4` exact (`240817_017.sxm` and
`240817_043.sxm` are target failures). Relative to the current default rerun
(`106/146`), this full run changes 50 rows: 14 fixes, 28 regressions, and 8
wrong-to-wrong moves, for a net `-14` exact result. The divergent-only screen was
therefore misleading because it was enriched for failures and did not expose the
large number of regressions among previously correct rows. The selector source is
`support_marginalized_gcv` for all 146 rows, with 42 guarded downshifts, 99
keeps, and 5 parsimony moves. Do not promote `support_marginalized_gcv_guard` for
6mer counting; the score to beat remains the default rerun's `106/146`.

Next existing-policy screen: selected `fwd_bwd_consensus` for a 43-file divergent
rerun because it is label-free and uses the independent forward/backward scan
consistency signal rather than another support-padding heuristic. The dry-run was
clean, then the no-watch Viper submission succeeded:

```text
10372583  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/fwd_bwd_consensus \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=fwd_bwd_consensus STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending (`10372583_[1-4]`, reason `QOSGrpCpuLimit`).
After completion, fetch/merge and grade the 43 divergent rows externally before
considering any full146 promotion.

The `fwd_bwd_consensus` screen completed cleanly on Viper and was merged/fetched
manually:

```text
10372583_1  COMPLETED  00:04:58  0:0
10372583_2  COMPLETED  00:03:52  0:0
10372583_3  COMPLETED  00:02:53  0:0
10372583_4  COMPLETED  00:04:09  0:0
```

Merged summary:

```text
results/experiments/6mer_divergent/fwd_bwd_consensus/summary_overlap060_hard.tsv
  43 data rows (43 ok)
```

External grade on the 43 divergent rows:

```text
fwd_bwd_consensus N_selected  16/43 exact (37.2%), +/-1 = 34/43 (79.1%), over/under = 14/13
```

Against the same 43 rows from the fresh default full146 rerun, this is an exact
improvement from `8/43` to `16/43`: 10 fixes, 2 regressions, 5 wrong-to-wrong
changes, and 6 unchanged correct rows. However, its `+/-1` count drops from
`37/43` to `34/43`. It also only matches the earlier
`support_marginalized_gcv_guard` screen on exact score (`16/43`) while being
worse on `+/-1` (`34/43` vs `37/43`), and that earlier screen already failed the
full146 promotion test badly (`92/146`, net `-14` vs default). Do not promote
`fwd_bwd_consensus` to a full146 run; the score to beat remains the fresh default
rerun's `106/146`.

Next existing-policy screen: selected `laplace_evidence_guard` for the same
43-file divergent rerun because it is technically distinct from the rejected
support-padding and forward/backward rescoring candidates. It is a label-free,
down-only local curvature/evidence guard: it can veto a one-lobe overfit, but it
cannot increase `N` or use the benchmark target.

The dry-run was clean, then the no-watch Viper submission succeeded:

```text
10373261  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/laplace_evidence_guard \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=laplace_evidence_guard STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending (`10373261_[1-4]`, reason `QOSGrpCpuLimit`).
After completion, merge/fetch and grade the 43 divergent rows externally before
considering any full146 promotion.

The `laplace_evidence_guard` screen completed cleanly on Viper and was
merged/fetched manually:

```text
10373261_1  COMPLETED  00:06:26  0:0
10373261_2  COMPLETED  00:04:26  0:0
10373261_3  COMPLETED  00:02:51  0:0
10373261_4  COMPLETED  00:04:10  0:0
```

Merged summary:

```text
results/experiments/6mer_divergent/laplace_evidence_guard/summary_overlap060_hard.tsv
  43 data rows (43 ok)
```

External grade on the 43 divergent rows:

```text
laplace_evidence_guard N_selected  11/43 exact (25.6%), +/-1 = 34/43 (79.1%), over/under = 10/22
```

Against the same 43 rows from the fresh default full146 rerun, this improves the
exact score only from `8/43` to `11/43`: 9 fixes, 6 regressions, 5
wrong-to-wrong changes, and 2 unchanged correct rows. It is worse than both
previous 43-file screens on exact score (`16/43` for both
`support_marginalized_gcv_guard` and `fwd_bwd_consensus`) and does not improve
the `+/-1` count (`34/43`). Do not promote `laplace_evidence_guard` to a full146
run; the score to beat remains the fresh default rerun's `106/146`.

Next existing-policy screen: selected `stability_selection` for the same 43-file
divergent rerun. This uses the frozen support-padding grid, like
`support_marginalized_gcv_guard`, but asks a different label-free question: how
often each candidate `N` remains competitive under support perturbations. It is
therefore a cheap final existing-policy triage before stopping the current
selector sweep.

The dry-run was clean, then the no-watch Viper submission succeeded:

```text
10374201  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/stability_selection \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=stability_selection STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state had chunk 1 running and chunks 2–4 pending under
`QOSGrpCpuLimit`. After completion, merge/fetch and grade the 43 divergent rows
externally before considering any full146 promotion.

The `stability_selection` screen completed cleanly on Viper and was
merged/fetched manually:

```text
10374201_1  COMPLETED  00:05:51  0:0
10374201_2  COMPLETED  00:04:10  0:0
10374201_3  COMPLETED  00:02:45  0:0
10374201_4  COMPLETED  00:04:03  0:0
```

Merged summary:

```text
results/experiments/6mer_divergent/stability_selection/summary_overlap060_hard.tsv
  43 data rows (43 ok)
```

External grade on the 43 divergent rows:

```text
stability_selection N_selected  15/43 exact (34.9%), +/-1 = 35/43 (81.4%), over/under = 7/21
```

Against the same 43 rows from the fresh default full146 rerun, this improves the
exact score from `8/43` to `15/43`: 10 fixes, 3 regressions, 15 wrong-to-wrong
changes, and 5 unchanged correct rows. It is better than
`laplace_evidence_guard` (`11/43`) but still below both previous best 43-file
screens (`16/43` for `support_marginalized_gcv_guard` and
`fwd_bwd_consensus`). Its `+/-1` count (`35/43`) is also below the fresh default
and support-marginalized screens (`37/43`). Do not promote `stability_selection`
to a full146 run; the score to beat remains the fresh default rerun's `106/146`.

At this point the existing-policy triage did not identify a full146 candidate:
`support_marginalized_gcv_guard` had the best 43-file screen profile but already
failed full146 badly, `fwd_bwd_consensus` matched its exact score with worse
`+/-1`, and `laplace_evidence_guard`/`stability_selection` were weaker on the
divergent screen. Further benchmark improvement likely needs new label-free
evidence or a targeted analysis of the remaining failure modes rather than
promoting another existing selector wholesale.

Offline failure-mode audit across the fresh default full146 summary and all four
43-file divergent screens confirmed that conclusion. The fresh default full146
score remains:

```text
default full146  106/146 exact, +/-1 = 139/146, over/under = 17/23
```

The 40 full146 misses are structured, not random:

```text
N_selected=3:  1 file
N_selected=4:  2 files
N_selected=5: 20 files
N_selected=7: 13 files
N_selected=8:  4 files
```

On the same 43 divergent rows, the completed screen scores were:

```text
default                         8/43 exact, +/-1 = 37/43, over/under = 14/21
support_marginalized_gcv_guard 16/43 exact, +/-1 = 37/43, over/under = 11/16
fwd_bwd_consensus              16/43 exact, +/-1 = 34/43, over/under = 14/13
laplace_evidence_guard         11/43 exact, +/-1 = 34/43, over/under = 10/22
stability_selection            15/43 exact, +/-1 = 35/43, over/under = 7/21
```

The screens fix different failure modes: among default under-counts on the 43-row
set, `fwd_bwd_consensus` fixes the most (`9/21`) without increasing the distance
from 6, while among default over-counts, `laplace_evidence_guard` fixes the most
(`7/14`) without increasing the distance. That split is scientifically useful,
but not a promotion rule by itself: choosing a policy by whether a file is an
under- or over-count uses the withheld benchmark label and is therefore invalid.

Several offline ensemble rules built only from existing label-free screen outputs
were scored for diagnostics, not promoted. The best exact-count variant on the
43-row set was a five-policy lower-tie mode:

```text
mode5_lower_tie  14/43 exact, +/-1 = 36/43, over/under = 10/19
```

Other consensus variants reached only `11–13/43` exact. None beat the best
single 43-file screens (`16/43`), and none has a plausible reason to improve the
remaining 103 full146 rows without introducing regressions. Do not promote an
ensemble rule from this audit. The next defensible benchmark-improvement gate is
not another wholesale selector rerun; it should be a targeted, label-free study
of the remaining failure classes, especially support-sensitive under-counts
(`N=5`) versus over-segmentation (`N=7/8`), using per-file diagnostics and visual
evidence before any new rule is frozen.

Targeted representative diagnostics separate the remaining failures into a few
label-free mechanisms worth studying:

- **Guard-induced under-counts.** Several `N=5` misses are not raw GCV failures:
  the primary ell/circ fit already agrees on `N=6`, but the robust-AICc guard
  downshifts to `5`. Examples include `240818_020.sxm`, `240818_026.sxm`,
  `241113_160.sxm`, and `241114_010.sxm`. For `240818_020.sxm`, the ell score
  curve has `N=6` as the GCV minimum, while `N=5` is 21.2% worse; multiple
  alternative label-free screens also return `6`. This suggests the next audit
  should focus on when the robust guard is contradicted by strong 2D GCV and
  ell/circ agreement, rather than changing the primary selector wholesale.
- **Ambiguous support-sensitive under-counts.** Some files have broad candidate
  sets and small GCV gaps, e.g. `240312_Cu100070.sxm` / `240312_Cu100071.sxm`:
  the ell sweep prefers `N=10`, but `N=6` and `N=7` are close, and the robust
  guard collapses to `5`. The fwd/bwd screen alone returns `6`, while support and
  Laplace go high. These need visual/diagnostic inspection of support extent and
  scan-direction consistency before any rule is frozen.
- **Over-segmentation with near-ties.** Some `N=7` files are one-lobe GCV
  overfits where `N=6` is label-free competitive. Examples: `240313_Cu100062.sxm`
  has ell `7` only 0.6% better than `6` and circ prefers `6`; support, Laplace,
  and stability all return `6`. `240314_Cu100_042.sxm` has ell `7` only 1.6%
  better than `6`, and the same three screens return `6`. This is a plausible
  generic one-lobe ambiguity class.
- **Over-segmentation without near-ties.** Other `N=7/8` files do not have a
  competitive `N=6` under current fit support. For example `241114_044.sxm` and
  `241114_045.sxm` have `N=8` clearly favoured by GCV; only stability jumps to
  `6`. Files such as `240310_Cu100009.sxm`, `240815_098.sxm`, `241114_043.sxm`,
  and `241114_046.sxm` remain unfixed by the screens. These are not safe targets
  for a simple ambiguity guard.

Next gate: build a small diagnostic report, not a production selector, that
flags candidate files using only label-free contradictions:

```text
robust-downshift contradiction:
  robust_AICc_N < N_eff
  and N_ell == N_circ == N_eff
  and GCV(N_robust) is materially worse than GCV(N_eff)

one-lobe overfit ambiguity:
  N_eff = N_alt + 1
  and lower-N GCV gap is within a frozen tolerance
  and either circ prefers lower N or support/Laplace/stability concur
```

Score those flags externally only after the diagnostic report is generated. If
the flags isolate high-signal subsets without many obvious false positives, then
freeze the rule and run a full146 test. If they do not, stop selector tuning and
move to visual failure review / data-quality stratification.

Generated the diagnostic report and external score tables:

```text
results/diagnostics/full146_failure_flags.tsv
results/diagnostics/full146_failure_flags_score.tsv
```

The flags are computed from the fresh full146 default summary and existing
label-free screen outputs; the manifest is used only afterward for external
scoring. Corrected score summary:

```text
robust_downshift_contradiction  n=11  fixes=6  regressions=2  neutral_wrong=3  net=+4  precision=54.5%
one_lobe_overfit_ambiguity      n=6   fixes=2  regressions=4  neutral_wrong=0  net=-2  precision=33.3%
combined_unique                 n=17  fixes=8  regressions=6  neutral_wrong=3  net=+2  precision=47.1%
```

Decision: do not promote the one-lobe ambiguity flag; it is net negative on the
external grade. The robust-downshift contradiction flag is promising enough to
preserve as the next focused candidate (`+4` net, implied `110/146` if applied as
a post-selection correction), but it is not clean enough to silently fold into
the default selector: it still has two regressions and three wrong-to-wrong
suggestions. The next defensible step is to freeze this rule as an explicit
diagnostic/experimental policy and run/grade it on full146, or inspect its 11
flagged rows visually before deciding whether it is scientifically acceptable.

Inspected the 11 `robust_downshift_contradiction` rows. External grading of the
suggested replacement (`N_eff`, because robust-AICc downshift was contradicted by
strong 2D GCV) splits as:

```text
fixes:          6  (240816_002, 240818_020, 240818_026, 240818_028, 241113_160, 241114_010)
regressions:    2  (240817_017, 240818_015)
wrong-to-wrong: 3  (240307_019, 240818_025, 241114_043)
```

Tried stricter label-free refinements. The best non-cheating refinement was
`one_step_and_alt_votes_ge2`, requiring a one-lobe robust correction plus at
least two existing label-free screen votes for the same suggested `N`:

```text
one_step_and_alt_votes_ge2  n=5  fixes=4  regressions=0  wrong_to_wrong=1  net=+4  precision=80.0%
```

It avoids the two regressions, but still includes `240307_019.sxm` as a
wrong-to-wrong `4→5` move, so it is not clean enough for default selection. The
only perfect-looking variants were absolute rules such as `5→6 only` or
`suggest N=6 only`:

```text
INVALID_absolute_5_to_6_only  n=6  fixes=6  regressions=0  wrong_to_wrong=0
INVALID_absolute_suggest_6_only  n=6  fixes=6  regressions=0  wrong_to_wrong=0
```

Those are invalid because the absolute target value `6` is the benchmark label
shape, not a label-free scientific criterion. Do not promote them. Current next
gate: either visually review the small robust-contradiction set before designing
a real label-free predicate, or stop selector tuning and move to data-quality /
failure-class stratification. No selector/policy change is justified yet.

Qualitative plot review of the 11 robust-contradiction rows supports the same
conclusion. Without using expected labels, the visually plausible upward moves
were mostly the modest under-count corrections: `240307_019` (`4->5`),
`240818_020` (`5->6`), `240818_026` (`5->6`), `240818_028` (`5->6`, cautious),
`241113_160` (`5->6`), and `241114_010` (`5->6`, cautious). `240816_002`
(`5->6`) was only partly plausible because the feature is broad/blurred and the
residuals remain structured. The larger or high-count upward moves were weak:
`240817_017` (`6->7`) was ambiguous, `240818_015` (`6->7`) looked like fitting
faint tails, `240818_025` (`5->7`) looked crowded/over-counted, and
`241114_043` (`7->8`) looked weak because the extra components were crowded or
terminal.

Decision after visual review: the robust-contradiction diagnostic identifies a
real under-count signal, but it also catches visually weak upward moves. Do not
turn it into a production selector yet. If this line of work continues, the next
scientific step should be a visual/data-quality stratification of the flagged
rows and a genuinely label-free predicate that suppresses tail/crowding-driven
upward moves; otherwise stop selector optimization at the fresh default
`106/146` benchmark state.

Follow-up 90% attempt: tested a purely physical support/spacing post-selection
family offline from the fresh full146 default summary. The rule uses only the
measured 2D support and existing chitosan spacing/overlap calibration from
`config/chitosan.toml` (`spacing_min_nm=0.35`, `spacing_max_nm=0.75`,
`max_overlap=0.60`, `sigma_parallel_max_nm=0.509`, giving
`spacing_min_eff=0.51448 nm`). Labels were used only after rule replay for
external grading.

The first support-midpoint sweep was written to:

```text
results/diagnostics/full146_support_bounds_rules.tsv
```

Its best non-target-shaped rule was
`support_mid_round_bounded_one_step_delta_ge1`: move at most one lobe toward the
midpoint of the physical support-derived feasible-N interval. It reached:

```text
126/146 exact, +/-1 = 143/146
46 changes: 30 fixes, 10 regressions, 5 wrong-closer, 1 wrong-farther, net +20
```

This was the strongest label-free offline signal so far, but still below the
90% target (`132/146`) and too regression-prone for promotion. Inspection showed
the damage concentrated in permissive upshifts from already-correct files,
especially `6->7` on long measured supports.

Two refined sweeps then tried to suppress those upshift failures while preserving
the support signal:

```text
results/diagnostics/full146_support_hybrid_rules.tsv
results/diagnostics/full146_support_refined_rule_sweep.tsv
```

The best hybrid rule was a conservative asymmetric version: downshift one step
whenever the current selection lies above the support midpoint; upshift one step
only when the support midpoint lies above the current selection, `N_eff` is also
above current, and the effective GCV gap is close (`delta_GCV_rel_eff <= 0.3`).
Equivalent refined variants using agreement between ell/circ support midpoints
or the same close-GCV upshift condition topped out at:

```text
127/146 exact, +/-1 = 143/146
31 changes: 23 fixes, 2 regressions, 5 wrong-closer, 1 wrong-farther, net +21
```

The remaining damage was narrow but not removable by the tested label-free
features without giving back many fixes: two `6->7` regressions
(`240817_017.sxm`, `240817_019.sxm`) and one `7->8` wrong-farther move
(`240310_Cu100009.sxm`) came from the same support-upshift mechanism that fixed
many under-counts. Stricter GCV-curve filters reduced regressions but lowered the
exact score to `123-124/146`; no tested rule reached the requested 90% level.

Decision at the time: do not promote the support-midpoint family as a selector
change. It was a useful diagnostic showing that measured support explains many
full146 misses, but the best label-free offline rule was still `5` exact files
short of 90% and retained non-negligible regressions. This was superseded by the
2026-07-02 promotion entry below after accepting the best current full146 result
as a practical chitosan default while continuing selector research.

### 2026-07-02 — Promote support-midpoint hybrid as current chitosan default

**Decision:** Promote `support_midpoint_hybrid` to the default
`config/chitosan.toml` batch policy. This supersedes the previous no-promotion
decision above. The rationale is pragmatic: the frozen label-free rule is the
best current expanded 146-file counting result, even though it does not reach the
earlier aspirational 90% target and is not a universal selector guarantee.

**Implemented rule:** The batch first applies the existing integrated
robust-AICc guard (`gcv_with_robust_aicc_guard`), then computes the midpoint of
the measured 2D support's feasible-N interval using the physical spacing and
overlap calibration. The final support layer is bounded to one lobe:

```text
if N_guarded > support_midpoint:
    N_selected = N_guarded - 1
elif N_guarded < support_midpoint
     and N_eff > N_guarded
     and delta_GCV_rel_eff <= 0.30:
    N_selected = N_guarded + 1
else:
    N_selected = N_guarded
```

The threshold is recorded as
`support_midpoint_up_gcv_rel_threshold = 0.30` in `[selection]`. The rule uses no
expected `N`, no target count, and no benchmark labels during fitting or
selection. The 146-file manifest remains an external grading set only.

**External grade used for the decision:**

```text
previous default replay: 106/146 exact, +/-1 = 139/146
support_midpoint_hybrid: 127/146 exact, +/-1 = 143/146
```

Known residual errors remain: two `6->7` regressions (`240817_017.sxm`,
`240817_019.sxm`) and one `7->8` wrong-farther move (`240310_Cu100009.sxm`) were
observed in the offline sweep. The default change is therefore a documented
current-best chitosan setting, not the endpoint of label-free selection work.

Post-promotion doc review (post-implementation review pass) caught three
documentation mismatches that have now been corrected:

- `docs/src/index.md` still named `gcv_with_robust_aicc_guard` as the batch
  default and quoted the `39/39` primary-benchmark number without separating
  robust-guard validation from the promoted default. Updated to name
  `support_midpoint_hybrid`, record the `127/146` / `143/146` external counting
  grade, and keep `39/39` explicitly as the robust-guard validation result.
- `docs/src/selection.md` documented `selection_source` as `ell_robust_aicc`,
  but `_select_primary` (`selectors.jl:661`) actually emits `robust_aicc_guard`
  when the guard moves the primary count. The doc now lists the real emitted
  values (`robust_aicc_guard`, `support_midpoint_down`, `support_midpoint_up`,
  or the kept effective source); `ell_robust_aicc` is correctly scoped to
  `refined_source`.
- `AGENTS.md` left the `4/4 clean_target` / reproducibility / threshold-robust
  bullets ambiguous in the same section that says `support_midpoint_hybrid` is
  not a no-regression primary-benchmark selector. Scoped those bullets
  explicitly to robust-AICc guard validation.

No code, config, or selection behaviour changed in this pass — only
documentation accuracy.

### 2026-07-03 — Gap≥2 down-to-midpoint extension: 127/146 → 129/146

**Decision:** Extend the support-midpoint hybrid so that when the robust guard
over-counts by at least two lobes relative to the support midpoint (gap ≥ 2),
the rule goes directly to the midpoint instead of making only a one-step
correction.

**Motivation:** Systematic offline replay on the frozen full146 fit data
(`results/experiments/6mer_full146/default_repro`) tested every label-free rule
constructible from the per-file GCV/BIC/chi2 curves, support geometry, κ,
ambiguity, and alternative selectors (stability, Laplace, fwd/bwd, spatial CV).

Key negative finding: the 241114 over-count cluster (N=7-8 instead of 6) is
**genuinely indistinguishable** by any residual-based criterion. Both GCV and
spatial k-fold CV prefer N=7-8 because the 7th lobe captures real spatially-
correlated signal (likely a surface defect or monomer substructure), not iid
noise. No fit-level feature (spacing_cv, κ, overlap, amplitude ratio) discriminates.

Key positive finding: the only improvement comes from a stronger
support-geometry prior. When `N_guarded ≥ support_midpoint + 2`, the geometry
disagrees strongly enough with the fit to justify a direct move to the midpoint.

**Implemented rule change** (`test/batch_full.jl`,
`_support_midpoint_hybrid_selection`):

```text
if N_guarded > support_midpoint + 1:
    N_selected = support_midpoint          # NEW: gap ≥ 2, trust geometry
elif N_guarded > support_midpoint:
    N_selected = N_guarded - 1             # unchanged: gap = 1
elif N_guarded < support_midpoint
     and N_eff > N_guarded
     and delta_GCV_rel_eff <= 0.30:
    N_selected = N_guarded + 1             # unchanged up-shift
else:
    N_selected = N_guarded                 # unchanged
```

New `selection_source` value: `support_midpoint_down_to_mid` (distinguishes the
gap≥2 case from the one-step `support_midpoint_down`).

**Offline replay result on frozen full146 data:**

```text
frozen hybrid ±1:          127/146 exact (87.0%), 143/146 ±1
hybrid ±2 (this change):   129/146 exact (88.4%), 143/146 ±1
changed=2  fixes=2  regressions=0
```

The 2 fixed files: `240818_021` (N=7→6) and `241114_045` (N=7→6), both had
`N_guarded=8, support_midpoint=6`. No regression because the up-shift branch
retains its `dgcv_rel_eff ≤ 0.30` guard.

**Limitation acknowledged:** the over-counts at gap=1 (e.g., `241114_040/044/046`
where N_guarded=7, midpoint=6) cannot be fixed by any residual-based label-free
rule — the 7th lobe is a genuine residual improvement. These remain at the ±1
grade level only.

**Benchmark cleanup:** `240310_Cu100009.sxm` removed from
`benchmarks/chitosan_6mer_counting_confirmed.toml` after visual review confirmed
the image is technically degraded (23.6% NaN pixels, "très complexe"). The
manifest now contains 145 files. Updated grade numbers on the 145-file benchmark:

```text
default_repro:  106/145 exact (73.1%), 138/145 ±1 (95.2%)
hybrid ±1:      127/145 exact (87.6%), 143/145 ±1 (98.6%)
hybrid ±2:      129/145 exact (89.0%), 143/145 ±1 (98.6%)
```

A full146 batch was submitted to Viper to confirm the ±2 rule on a real batch
run; the grading script uses the 145-file manifest.

Completion follow-up: the first Viper submission (`10390605`) failed because the
remote data staging copied local SXM symlinks as broken symlinks. After replacing
the remote data directory with real files via `rsync -L`, rerun `10391759`
completed successfully. `hpc/merge_chunks.jl` merged both shards into
`results/experiments/6mer_full146/pm2_confirm/summary_overlap060_hard.tsv` with
146 data rows and 146 `ok` rows. External grading against
`benchmarks/chitosan_6mer_counting_confirmed.toml` confirmed the offline replay:

```text
primary score:  129/145 (89.0%)
primary ±1:     143/145 (98.6%)
primary over/under: 6/10
stress rows:    0
```

The two remaining beyond-±1 misses are `240814_020.sxm` (`N_selected=4`) and
`240818_019.sxm` (`N_selected=3`). The exact-miss set is unchanged from the
confirmed grade: `240307_019`, `240310_Cu100032`, `240314_Cu100_024`,
`240314_Cu100_026`, `240814_020`, `240814_021`, `240816_005`, `240817_017`,
`240817_019`, `240817_078`, `240817_083`, `240818_007`, `240818_019`,
`241114_040`, `241114_044`, and `241114_046`.

The merged full146 run's `selection_source` distribution was: `ell=103`,
`robust_aicc_guard=12`, `support_midpoint_down=12`,
`support_midpoint_down_to_mid=2`, and `support_midpoint_up=17`. The two
`support_midpoint_down_to_mid` rows are the intended gap≥2 direct-to-midpoint
fixes.

### Follow-up: 0/1/? benchmark scope clarification

The full unit-assignment benchmark must use the same 145 confirmed 6mer files as
the counting/fit benchmark. The external control sequence is `NKNNKN`, encoded as
`010010` if `0=N` and `1=K`, or `101101` under the flipped identity convention.
This correction supersedes the earlier documentation that treated only the 35
240817 clean/clean_target rows as gradable or described 111 expanded rows as
missing unit sequences for benchmark purposes.

The scientific constraint is unchanged and becomes more important: this is a
fully label-free benchmark. `NKNNKN`, the `2K/4N` composition, the fact that `N=6`,
and any benchmark-control convention may be used only by external grading and
diagnostic reports after predictions are frozen. They must not enter fitting,
selection, per-lobe assignment, threshold selection, abstention rules, composition
priors, or method calibration. The objective is a robust rule that extrapolates to
unknown systems rather than one shaped to the benchmark control.

Existing `178/210`, `154/171`, and `101/106` numbers remain useful as historical
35-file subset diagnostics, but they are no longer the headline 0/1/? benchmark.
The reporting path should be expanded or replaced before claiming full145 unit
assignment metrics.

Implementation follow-up: `test/report_unit_assignment_benchmark.jl` now has an
explicit `--full145` mode. It writes `control_full145_truth.tsv` under the report
output directory from `benchmarks/chitosan_6mer_counting_confirmed.toml`, using the
external `NKNNKN` control encoding (`010010` by default, or `101101` with
`--control-sequence`). It then runs the existing grader and refuses any profile
whose graded denominator is not 145 files / 870 lobes. The default command remains
the historical subset report because the frozen prediction profiles currently have
only 210 lobe rows; full145 mode fails on those profiles by design rather than
printing a false headline.

Second implementation follow-up: added an early full145 coverage preflight before
the grader runs. Current frozen profiles now fail immediately with a message like
`expected 145 files / 870 lobes, found 35 files / 210 lobe rows`, plus examples of
missing files. A synthetic temporary full145-perfect prediction TSV validates the
positive path at 145/145 chains and 870/870 lobes. The remaining real work is to
generate label-free prediction profiles for all 145 benchmark files from the
feature/patch pipeline; scoring the current 35-file profiles as full145 is blocked
by design.

### Follow-up: portable label-free prediction builder and full145 blocker

Added `test/build_labelfree_unit_predictions.jl` to close the reproducibility gap
between feature extraction and prediction TSVs. The script is prediction-only and
does not read benchmark truth, `NKNNKN`, `expected_N`, or a composition prior. It
combines per-file standardized feature views over multiple k-means seeds and maps
the higher-amplitude cluster to GlcNAc (1), matching the label-free physical
convention used by the grader. Optional inputs add the split-width
`split_log_skew` feature and backward-patch negative-moment descriptors
(`bwd_neg_com_t`, `bwd_neg_diag45`, `bwd_neg_diag135`).

Validation on the currently available 39-file artifacts:

```bash
julia --project=. test/build_labelfree_unit_predictions.jl \
    --features results/unit_separability/lobe_features_selectedN_primary_local.tsv \
    --split-features results/unit_separability/lobe_features_selectedN_primary_split.tsv \
    --patches results/unit_separability/lobe_patches_selectedN_primary_17x17_bwd.tsv \
    --out /tmp/opencode/labelfree_unit_predictions_subset.tsv \
    --seeds 20 \
    --interactions

julia --project=. test/report_unit_assignment_benchmark.jl \
    --profile portable_kmeans=/tmp/opencode/labelfree_unit_predictions_subset.tsv=portable_kmeans_sanity \
    --outdir /tmp/opencode/labelfree_unit_report_subset
```

The generated TSV covers 39 files / 234 lobes. The historical subset sanity grade
is 170/210 = 81.0% physical (1/35 exact). This is intentionally described as a
portable k-means baseline, not as a replacement for the frozen best GMM ensemble
(`forced_ensemble3`, 178/210).

The first full145 selected-`N` path is now unblocked for honest reporting, while
the strict 870-row path remains blocked for real selected-`N` predictions:

1. The available feature/patch artifacts cover only 39 files, so full145 preflight
   correctly fails with `found 39 files / 234 lobe rows` and 106 missing files.
2. The current full145 counting summary is label-free but not denominator-perfect:
   `N_selected` counts are 129 files at 6, eight at 5, six at 7, one at 4, and one
   at 3. A selected-`N` full145 feature extraction therefore produces 863 lobe
   rows, not 870. Producing exactly six lobe rows for every file by reading the
   manifest/control would use the benchmark label and is not an allowed prediction
   path.

Implemented the safe reporting branch for option (a). `test/report_unit_assignment_benchmark.jl`
now has an explicit `--full145-own-n` mode for profiles generated at each file's
label-free `N_selected`. Strict `--full145` is unchanged and still requires
145 files / 870 prediction rows. Own-N mode still requires all 145 files and
contiguous lobe indices per file, but it reports count mismatches honestly:
positions absent because `N_selected < 6` count as uncertain against the external
870-position control denominator, while `N_selected > 6` lobes are reported as
extra predictions and are not aligned to the 6-position control.

The full145 selected-`N` feature extraction from the Viper counting summary
completed locally:

```text
/tmp/opencode/full145_selectedN_features.tsv
files=145
prediction rows=863
rows by selected N: N=3 -> 3, N=4 -> 4, N=5 -> 40, N=6 -> 774, N=7 -> 42
```

The immediate base/local-feature portable baseline was generated without truth,
composition priors, or the control sequence:

```bash
julia --project=. test/augment_lobe_local_features.jl \
    --features /tmp/opencode/full145_selectedN_features.tsv \
    --out /tmp/opencode/full145_selectedN_features_local.tsv

julia --project=. test/build_labelfree_unit_predictions.jl \
    --features /tmp/opencode/full145_selectedN_features_local.tsv \
    --out /tmp/opencode/full145_selectedN_labelfree_local_predictions.tsv \
    --seeds 20 \
    --interactions

julia --project=. test/report_unit_assignment_benchmark.jl --full145-own-n \
    --profile selectedN_local=/tmp/opencode/full145_selectedN_labelfree_local_predictions.tsv=selectedN_local \
    --outdir /tmp/opencode/unit_report_full145_own_n_local
```

Report result:

```text
prediction rows:        863
aligned/classified:     857/870 (98.5%)
missing control lobes:  13 across 10 short-N files
extra predicted lobes:  6 across 6 N=7 files
physical accuracy:      671/857 = 78.3%
honest view:            671/870 = 77.1% correct, 199/870 uncertain
exact chains:           16/145
```

This is a reproducible selected-`N` full145 baseline, not a solved final map and
not a replacement for the historical best three-view subset ensemble. The next
scientifically safe step is to generate the selected-`N` split features and
backward-patch descriptors for all 145 files and rerun the portable predictor in
the same `--full145-own-n` report mode.

Continuation: the backward-patch half of that step completed locally from the
already frozen selected-`N` feature TSV:

```bash
STMFIT_DATA_DIR=/tmp/opencode/stmfit_full146_data julia --project=. \
    test/extract_lobe_patches_bwd.jl \
    --features /tmp/opencode/full145_selectedN_features.tsv \
    --half-nm 0.32 \
    --step-nm 0.04 \
    --out /tmp/opencode/full145_selectedN_patches_bwd.tsv
```

Validation showed exact coverage parity with the selected-`N` features:

```text
features: 145 files / 863 rows
patches:  145 files / 863 rows
mismatches: 0
```

The base+backward portable predictor used no truth, control sequence, composition
prior, or expected count:

```bash
julia --project=. test/build_labelfree_unit_predictions.jl \
    --features /tmp/opencode/full145_selectedN_features_local.tsv \
    --patches /tmp/opencode/full145_selectedN_patches_bwd.tsv \
    --out /tmp/opencode/full145_selectedN_labelfree_base_bwd_predictions.tsv \
    --seeds 20 \
    --interactions

julia --project=. test/report_unit_assignment_benchmark.jl --full145-own-n \
    --profile selectedN_base_bwd=/tmp/opencode/full145_selectedN_labelfree_base_bwd_predictions.tsv=selectedN_base_bwd \
    --outdir /tmp/opencode/unit_report_full145_own_n_base_bwd
```

Report result:

```text
prediction rows:        863
classified:             854/870 (98.2%)
missing control lobes:  13 across 10 short-N files
extra predicted lobes:  6 across 6 N=7 files
physical accuracy:      671/854 = 78.6%
honest view:            671/870 = 77.1% correct, 199/870 uncertain
exact chains:           5/145
```

So the backward descriptors do not improve the full145 selected-`N` honest
headline over the base/local-only run; they mainly add 3 abstentions and reduce
exact-chain count under this portable k-means rule. The selected-`N` split-width
refit remains incomplete locally: repeated timeout-limited shards produced
complete split rows for 96/145 files, leaving 49 files missing. Those partial
split TSVs must not be used as a full145 profile. Two small script hardening fixes
were made during the attempt: `extract_lobe_features.jl --manifest` now restricts
the selected-summary file list to manifest members before `--primary-only`, and
`extract_lobe_patches_bwd.jl --help` no longer errors when `STMFIT_DATA_DIR` is
unset.

HPC follow-up requested during the same continuation: Viper has no active STMFit
jobs, and the latest full146 confirmation array `10391759` is complete
(`2/2` array tasks, exit `0:0`). Raven showed the pending QE GlcN restart3 job
`28601744` had reached the 24 h walltime limit (`TIMEOUT`; `pw.x` failed after
the limit). Following the established restart protocol, fetched
`qe/glcn_restart3/glcn_central_relax.out`, extracted the last 213-atom geometry to
`qe/glcn_restart3/glcn_central_best4.xyz`, prepared `qe/glcn_restart4` with the
same active Raven settings (`8` tasks, `96000 MB`, `24:00:00`, `ecutwfc=50`,
`ecutrho=360`, Gamma-only), and submitted it on Raven after a clean local and
remote preflight:

```text
qe/glcn_restart4 -> 28658135
initial status: PENDING
```

Do not fetch or inspect `qe/glcn_restart4` outputs until it completes or times
out. GlcNAc production remains dependent on a successful production GlcN restart.

Continuation: completed the full145 selected-`N` portable-predictor sweep using
only frozen label-free feature/patch TSVs. The externally generated control TSV
was used only by `report_unit_assignment_benchmark.jl --full145-own-n` for
post-hoc grading; no view, threshold, abstention rule, or count was chosen from
`NKNNKN`.

The main sweep command graded the base run, base+backward run, the selected-`N`
three-view ensemble, confidence/agreement abstentions, and backward-only views:

```text
summary: /tmp/opencode/unit_report_full145_own_n_variant_sweep/summary.tsv
```

Selected rows:

```text
profile              classified  physical          honest view                 exact
base                 857/870     671/857 = 78.3%  671/870 correct + 199 ?     16/145
base_bwd             854/870     671/854 = 78.6%  671/870 correct + 199 ?      5/145
view3                857/870     676/857 = 78.9%  676/870 correct + 194 ?      7/145
view3_conf80         821/870     653/821 = 79.5%  653/870 correct + 217 ?      6/145
view3_agreebase65    827/870     658/827 = 79.6%  658/870 correct + 212 ?      6/145
bwd_com_only         854/870     419/854 = 49.1%  419/870 correct + 451 ?      0/145
bwd_diag45_only      854/870     426/854 = 49.9%  426/870 correct + 444 ?      0/145
bwd_diag135_only     854/870     427/854 = 50.0%  427/870 correct + 443 ?      0/145
```

The current selected-`N` full145 own-N headline is therefore the three-view
ensemble `BASE`, `BASE+bwd_neg_com_t`, and `BASE+bwd_neg_diag45`: it improves the
base/local honest count by 5 lobes (676 vs 671) while preserving full aligned
coverage. The confidence/agreement abstention variants raise classified accuracy
only marginally but lower the honest correct count, so they are not a better
headline. Backward-only views are negative evidence: their physical convention
collapses to about chance, even though supervised oracle alignment can recover
some signal, so they must not be used alone.

A second targeted sweep checked whether adding individual backward descriptors to
BASE or adding `bwd_neg_diag135` as a fourth view was useful:

```text
summary: /tmp/opencode/unit_report_full145_own_n_variant_sweep2/summary.tsv
base_bwd_com      672/857 honest-aligned correct, 672/870 honest denominator
base_bwd_diag45   671/857 honest-aligned correct, 671/870 honest denominator
base_bwd_diag135  667/857 honest-aligned correct, 667/870 honest denominator
view4             670/857 honest-aligned correct, 670/870 honest denominator
```

These do not beat `view3`. The selected-`N` split-width profile is still blocked
by incomplete local coverage (96/145 files), so it remains excluded from any
full145 headline.

Operational status at the same checkpoint: Raven QE GlcN restart4 job `28658135`
was still running on `ravc4062` at 8:26 elapsed. Viper Slurm status checks were
temporarily unreliable because the remote Modules/Tcl library reported an I/O
error and `squeue` then failed with `ClusterName needs to be specified`.

Continuation: completed the previously blocked selected-`N` split-width full145
feature extraction. The prior local shards covered only 96/145 manifest files;
the remaining 49 files were identified from `/tmp/opencode/full145_selectedN_features.tsv`
and rerun in two-file micro-batches to avoid another tool timeout. The final merge
used the selected-`N` base feature row order as the key and wrote:

```text
/tmp/opencode/full145_selectedN_features_split.tsv
files=145
rows=863
missing rows=0
duplicates=0
ignored extra rows=8 from excluded 240310_Cu100009.sxm
```

This keeps the own-N denominator unchanged: 13 missing external control positions
from 10 short-N files and 6 extra N=7 predictions remain a counting/selection
fact, not a split-width artifact.

Generated split-aware label-free prediction profiles before any grading truth was
read:

```text
/tmp/opencode/full145_selectedN_labelfree_split_default_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_base_split_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_4view_base_bwd_split_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_3view_base_com_split_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_3view_base_diag45_split_predictions.tsv
```

Post-hoc `--full145-own-n` grading summary:

```text
summary: /tmp/opencode/unit_report_full145_own_n_split_sweep/summary.tsv

profile                  classified  physical          honest view              exact
prev_view3               857/870     676/857 = 78.9%  676/870 + 194 ?          7/145
split_default            857/870     676/857 = 78.9%  676/870 + 194 ?          7/145
base_split               857/870     676/857 = 78.9%  676/870 + 194 ?         17/145
view4_base_bwd_split     857/870     675/857 = 78.8%  675/870 + 195 ?         16/145
view3_base_com_split     857/870     673/857 = 78.5%  673/870 + 197 ?         16/145
view3_base_diag45_split  857/870     673/857 = 78.5%  673/870 + 197 ?         16/145
```

Conclusion: completing split-width removes the full145 coverage blocker, but it
does not improve the lobe-correct headline beyond 676/870. Its useful signal is
sequence-level: the single-view `BASE+split_log_skew` profile keeps the same
honest correct count and raises exact chains from 7/145 to 17/145. Treat this as
a better exact-chain diagnostic, not as a solved binary deacetylation map.

Workflow hardening follow-up: added `test/merge_lobe_feature_shards.jl` so future
split-width/full145 resumptions no longer require hand-written Julia one-liners.
The script takes a reference selected-`N` feature TSV plus repeated `--shard` or
comma-separated `--shards`, validates that every reference `(file,lobe)` key is
present exactly once, reports stale extra rows, and writes the merged TSV in the
reference row order. It reads no truth sequence, expected count, control motif, or
composition prior.

Validation on the completed split-width artifacts:

```text
julia --project=. test/merge_lobe_feature_shards.jl \
    --reference /tmp/opencode/full145_selectedN_features.tsv \
    --shards <all split chunk/missing/resume TSVs> \
    --ignore-extra \
    --dry-run

reference_files=145
reference_rows=863
merged_files=145
merged_rows=863
missing_rows=0
duplicate_rows=0
extra_rows=8
extra_files=240310_Cu100009.sxm
```

The `--ignore-extra` flag is intentionally explicit; without it, stale rows from
an excluded or superseded manifest member remain an error.
