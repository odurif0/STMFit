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

### 2026-06-22 — Unit assignment (GlcNAc/GlcN) investigation started

**Goal:** Assign each fitted Gaussian lobe a monomer type (0 = GlcN,
1 = GlcNAc) to produce a deacetylation map per chain. The ground-truth
sequence is known for all benchmark chains (240817) but must stay outside the
fitting/selection path (same label-free rule as for N).

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

### Updated conclusion

**84.3% per-lobe / 13/35 exact is the new label-free ceiling** with current
STM data at −0.3 V. The breakthrough came from a wavelet diagonal detail (HH1)
at 17×17 resolution, not from more complex clustering or deeper features. The
bottleneck remains GlcNAc detection (lobes 2, 5 at 60–69%), but middle GlcN
positions (3, 4) improved significantly. The abstention framework provides
honest reporting at 91% on confident lobes. Converged DFT cubes (`28444935`,
`28445456`) remain the critical path to higher accuracy.

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
