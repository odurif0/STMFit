# Unit Assignment (GlcNAc/GlcN per lobe)

## Motivation

Chitosan is a (1,4)-β-linked polysaccharide composed of two monomer units:
**GlcNAc** (N-acetyl-glucosamine) and **GlcN** (glucosamine). The fitted
Gaussian lobes correspond to individual monomer units, but the current pipeline
only counts them (N_selected) — it does not identify which lobe is which unit.

The goal is to assign each fitted lobe a type (0 = GlcN, 1 = GlcNAc) to produce
a **deacetylation map** per chain: the ordered sequence of GlcNAc/GlcN along
the molecular backbone.

## Label-free constraint (extended)

The same label-free rule that applies to N also applies to unit assignment:

- The **ground-truth/control sequence** for the 6mer benchmark is `NKNNKN`, encoded
  as `010010` or `101101` depending on the 0/1 identity convention, across the
  same 145 files used by the counting benchmark. It must **never enter the
  fitting, selection, assignment, thresholding, abstention, or calibration path**.
  It is used only by grading scripts (`grade_unit_assignment.jl`) and explicitly
  supervised diagnostics such as `--with-truth` cross-evaluation.
- Using the truth to choose the 0↔1 flip (the "oracle" convention) is
  **supervised** and must be clearly labeled as such. The **physical convention**
  (GlcNAc = highest-amplitude cluster, based on the larger acetyl group) is
  label-free.
- The assignment rule must not assume the number of GlcNAc/GlcN units in a
  chain. Even though the external control has two `K` units and four `N` units,
  rules such as "top-2 lobes are GlcNAc" are composition priors and are diagnostic
  only, not valid label-free assignment.

## Pipeline overview

The unit-assignment investigation proceeds in phases. **Phases 0–2a are
implemented as diagnostics; Phases 2–5 remain research directions** pending a
robust label-free unit-identity signal.

```
Phase 0: Grading framework (grade_unit_assignment.jl)
   │     Truth TSV + predictions → accuracy, confusion, edit distance
   │     4 alignments (identity/reverse/flip/reverse+flip)
   │     2 conventions (physical label-free + oracle supervised)
   │
Phase 1: Gaussian feature separability (analyze_unit_separability.jl)
   │     extract_lobe_features.jl → per-lobe (A, σ∥, σ⟂, integrated)
   │     Unimodal vs bimodal test (kmeans k=1 vs k=2, BIC)
   │     With --with-truth: AUC per feature, clustering accuracy
   │
Phase 1b: Non-Gaussian residual features — REMOVED (added noise: ΔBIC +184→+4.3)
   │
Phase 1c: Local/envelope-corrected features (augment_lobe_local_features.jl)
   │     Per-chain z-scores, local prominence, envelope residuals
   │     No truth and no composition constraint
   │
Phase 1d: Aligned patch diagnostics (extract_lobe_patches.jl)
   │     Raw and residual patches aligned to chain axis
   │     PCA/kmeans + optional supervised train/test diagnostic
   │
Phase 1e: Split-width Gaussian forward model (GaussianFit2D peak_profile=:split)
   │     σ∥ is split left/right around each lobe center
   │     skew_ratio = σright / σleft, fitted per lobe
   │     Tests whether STM resolves any lobe asymmetry before DFT-STM molds
   │
Phase 2a: Connected mold-template decoding (score_connected_mold_templates.jl)
    │     Apply GlcN/GlcNAc patch molds with global direction/phase/mirror states
    │     Enforces glycosidic connectivity/orientation, not composition
    │     refine_geometric_mold.jl searches local acetyl-site transforms label-free
    │
Phase 2: 1-type vs 2-type model selection [planned]
   │     shared_sigma_types ∈ {0,1,2} × spacing_model ∈ {free, alternating}
   │     GCV comparison → does a 2-type structure exist?
   │
Phase 3: Per-blob clustering assignment [planned]
   │     GMM 2-component on Phase 1+1b features → sequence 0/1 per chain
   │     Physical mapping {A,B} → {GlcN,GlcNAc}
   │
Phase 4: Template supervised validation [planned]
   │     Train/test split → nearest-centroid → generalization measure
   │
Phase 5: DFT-STM simulation [optional, planned]
         LDOS on Cu(100) → physical template
```

## Split-width asymmetry test (Phase 1e)

The current production model uses symmetric 2D Gaussians. This can erase the
GlcNAc acetyl shoulder by absorbing it into a symmetric width and leaving only a
noisy residual. Phase 1e adds an opt-in split-width profile:

```text
σleft  = σ∥ / sqrt(skew_ratio)
σright = σ∥ * sqrt(skew_ratio)
skew_ratio ∈ [1/skew_ratio_max, skew_ratio_max]
```

`skew_ratio = 1` is exactly the symmetric Gaussian. The split model is therefore
a nested test of whether the STM data justify a fitted left/right lobe
asymmetry. It remains label-free: the skew parameter is fitted independently per
lobe, with no truth sequence and no composition constraint.

Recommended experimental use:

1. Keep the production batch and N selection on `config/chitosan.toml`.
2. Refit features at the already selected N using `config/chitosan_split.toml`.
3. Compare GCV at fixed N between Gaussian and split profiles.
4. Only if split improves GCV and `skew_ratio` is stable/bimodal should we move
   to DFT-STM molecule molds.

`extract_lobe_features.jl` automatically fixes `n_min=n_max=N_selected` per file
when `--selected-summary` is provided. This avoids an unnecessary N sweep and is
important for split-width diagnostics, which add one fitted parameter per lobe.
For adaptive-support counts, retain the original summary's `refined_policy`:
accepted rescue reuses the configured rescue support in both circular and
elliptical fits, including split-profile extraction. N is never reselected.

Decision rule:

- If split does not improve GCV and `skew_ratio≈1`, the STM/tip conditions do
  not resolve the asymmetry; geometry/DFT molds are unlikely to help.
- If split improves GCV and `skew_ratio` separates lobes, asymmetry is a real STM
  signal and a physical mold is justified.

## Key design decisions

### Why connected molds rather than free asymmetric lobes?

The split-width test showed that the STM data contain real lobe asymmetry, but
that generic left/right width is dominated by overlap/envelope/tip effects rather
than GlcNAc/GlcN identity. A molecule mold must therefore encode a **specific
chemical geometry**, not just arbitrary skew. For chitosan, the β-(1→4)
connectivity restricts how pyranose rings and C2 substituents can be oriented
relative to the backbone. The connected-mold decoder tests only a small set of
global states:

```text
chain direction × pyranose parity phase × surface/mirror state
```

Within a global state, each lobe is scored against GlcN and GlcNAc templates at
the orientation allowed by connectivity. The sequence is chosen by template cost,
with no assumption about how many GlcNAc/GlcN units are present.

Pairwise context is represented by **sliding bonds**, not by disjoint dimers:

```text
(1,2), (2,3), ..., (N-1,N)
```

This supports both even and odd chains. The total decoded cost is:

```text
Σ unary_i(type_i, parity_i, mirror)
+ Σ bond_i(type_i, type_{i+1}, parity_i, mirror)
```

The transition motifs `00`, `01`, `10`, and `11` are adjacent-pair costs used by
Viterbi; they never impose a composition or tile the chain into fixed dimers.

The first implementation is deliberately patch-level (`score_connected_mold_templates.jl`):
it consumes already aligned lobe patches and externally generated templates. A
full forward-model fit with templates should only be added if patch-level mold
scoring carries a useful signal.

### Connected mold template format

The connected mold workflow has two input routes. The current source is a manual
geometric proxy-site TSV; the preferred future source is DFT-STM maps:

```text
manual geometric proxy sites (current)
  → templates/chitosan_geometric_sites.tsv
  → generate_connected_mold_templates.jl
  → unary + sliding-bond template TSVs
  → score_connected_mold_templates.jl --template-mode contrast

DFT-STM/LDOS maps for GlcN and GlcNAc in beta-(1->4) chain context (preferred)
  → templates/chitosan_stm_maps.tsv
  → import_stm_mold_maps.jl
  → unary + sliding-bond template TSVs
  → score_connected_mold_templates.jl --template-mode contrast
```

For the manual geometric path, `templates/chitosan_geometric_sites.tsv` directly
defines proxy sites in the aligned `(t,u)` patch frame:

```text
type    atom              t_nm    u_nm    weight  sigma_t_nm  sigma_u_nm
0       ring_center       0.000   0.000   ...
0       glcn_nh2          ...
1       ring_center       0.000   0.000   ...
1       acetyl_carbonyl   ...
```

These are not atoms and not benchmark-fitted parameters. They encode a physical
first guess: a shared pyranose backbone plus a short GlcN substituent or a longer
GlcNAc acetyl-side proxy. Identity scoring should use `--template-mode contrast`
so the shared backbone does not dominate the weak substituent signal.

Templates are wide TSV files with one row per `(type, parity, mirror)`:

```text
name    type    parity  mirror  p001    p002    ...
GlcN_p0_m0      0       0       0       ...
GlcNAc_p0_m0    1       0       0       ...
...
```

Required rows are all eight combinations:

```text
type ∈ {0,1}, parity ∈ {0,1}, mirror ∈ {0,1}
```

`type=0` means GlcN and `type=1` means GlcNAc. The `pNNN` columns must match the
patch grid size from `extract_lobe_patches.jl`, after stripping the patch prefix
(`raw_p001`/`res_p001` in the patch TSV corresponds to `p001` in the template
TSV). Template pixels should be generated from geometric proxy sites or DFT-STM
maps in the same aligned coordinate convention as the patches and then normalized;
the scorer standardizes both patches and templates before comparing them.

`score_connected_mold_templates.jl` supports two template modes:

| Mode | Meaning | Use |
|---|---|---|
| `full` | Score the full GlcN/GlcNAc mold. | Default; useful when the whole STM patch is type-specific. |
| `contrast` | Subtract the parity/mirror common mold before scoring. | Diagnostic for weak type-specific substituent signal, avoids common-backbone dominance. |

### DFT-STM / LDOS Map Format

The ideal mold source is not a molecular contour. For STM it should be a simulated
STM/LDOS image at the experimental bias, after adsorption geometry and tip/filter
effects are represented as well as possible. `test/import_stm_mold_maps.jl`
imports such maps once they are exported in the aligned lobe frame:

```text
type    t_nm    u_nm    value    parity  mirror
0       -0.48   -0.48   ...      0       0
0       -0.40   -0.48   ...      0       0
1       -0.48   -0.48   ...      0       0
```

Required columns are `type`, `t_nm`, `u_nm`, and `value`. `type=0` is GlcN and
`type=1` is GlcNAc. `parity` and `mirror` are optional; if they are absent, the
importer generates the beta-(1->4) orientation variants by flipping the base maps.
The maps must already use the same coordinate convention as extracted patches:
`t` along the fitted backbone, `u` transverse, and origin at the central lobe.

The recommended scientific input is not isolated monomers, but a central unit in
a short linked oligomer, for example:

```text
GlcN central unit in a beta-(1->4) trimer
GlcNAc central unit in a beta-(1->4) trimer
```

This keeps the glycosidic linkage, ring parity, and nearest-neighbor electronic
context while leaving the sequence free during decoding. If explicit dimer/pair
maps are available, the importer can also consume optional bond maps with columns:

```text
left_type  right_type  side  t_nm  u_nm  value  parity  mirror
```

where `side` is `left`/`right` or `l`/`r`. If bond maps are not provided, bond
templates are generated by concatenating the corresponding unary maps, matching
the current geometric-template behavior.

### Residual refinement rule

The manual geometric mold is allowed to be refined against STM patches/residuals
only if the objective is label-free. `test/refine_geometric_mold.jl` implements
the first version: it searches small global transforms of the GlcNAc acetyl proxy
sites (`t/u` shifts, transverse scale, weight scale, sigma scale), scores each
candidate against contrast templates, and ranks candidates by k=1 vs k=2 BIC of
the resulting per-lobe template-evidence margins. It never reads truth labels and
never imposes a GlcNAc/GlcN composition.

Valid future objectives include reconstruction error, mean residual-template
correlation, cross-validated stability across files, or a shared low-dimensional
correction to the proxy-site positions/weights. Invalid objectives include
maximizing the known `010010` benchmark accuracy, forcing two GlcNAc units per
chain, or choosing among candidate molds by exact-sequence score. The mold must
be frozen before running `grade_unit_assignment.jl`.

Optional sliding bond templates use one row per
`(left_type, right_type, parity, mirror)`:

```text
name    left_type   right_type  parity  mirror  l_p001  ...  r_p001  ...
00_p0_m0        0   0           0       0       ...
01_p0_m0        0   1           0       0       ...
10_p0_m0        1   0           0       0       ...
11_p0_m0        1   1           0       0       ...
```

Required rows are all 16 combinations:

```text
left_type,right_type ∈ {00,01,10,11}, parity ∈ {0,1}, mirror ∈ {0,1}
```

The `l_pNNN` and `r_pNNN` columns are scored against the left and right observed
lobe patches for each sliding edge. In reversed chain direction, the decoder
automatically reverses the chemical left/right order of the transition.

### Why not fit "cross sections" directly?

The STM image at constant height measures the LDOS convolved with the tip, not
a geometric van der Waals cross section. A direct fit of molecular templates to
the STM data would require DFT-simulated STM images (LDOS on Cu(100)), which is
molécule-specific and not generalist. The graduated approach (geometry →
empirical template → DFT) avoids this until necessary.

### Why per-blob clustering (Phase 3) rather than alternating model (Phase 2)?

Chitosan is a **random copolymer** (degree of deacetylation, DD). The sequence
GlcNAc/GlcN along a chain is not strictly alternating. The `alternating` spacing
model (`shared_sigma_types=2, chain_spacing_model="alternating"`) tests whether
a 2-type structure exists, but cannot assign types in a random sequence. The
per-blob clustering (Phase 3) handles arbitrary sequences.

### Orientation ambiguity

The PCA axis (see `core.jl:700-701`) has a deterministic orientation
("increasing y, then x"), but the ground truth may be encoded in a
molecule-relative convention (reducing → non-reducing end). The grading script
tests **4 alignments** (identity, reverse, flip, reverse+flip) to handle this.

### Two conventions for the 0↔1 flip

| Convention | Flip resolution | Label-free? | Usage |
|---|---|---|---|
| **Physical** | GlcNAc = highest-amplitude cluster (acetyl is larger) | Yes | Phase 3 (main), Phase 0 grading |
| **Oracle** | Best of 4 alignments (uses truth to choose flip) | No (supervised) | Phase 0 (upper bound), Phase 4 |

The gap between physical and oracle accuracy is itself a diagnostic: a small
gap validates the physical convention; a large gap indicates the amplitude
mapping is imperfect.

## Scripts

| Script | Phase | Role | Needs SXM? |
|---|---|---|---|
| `test/extract_lobe_features.jl` | 1 | Re-run fit, extract per-lobe Gaussian features + axis | Yes (`STMFIT_DATA_DIR`) |
| `test/analyze_unit_separability.jl` | 1 | Unimodal vs bimodal test, AUC, clustering accuracy | No (reads TSV) |
| `test/augment_lobe_local_features.jl` | 1c | Add local prominence and envelope-corrected features | No (reads TSV) |
| `test/extract_lobe_patches.jl` | 1d | Extract chain-axis-aligned raw/residual patches | Yes (`STMFIT_DATA_DIR`) |
| `test/analyze_lobe_patches.jl` | 1d | PCA/kmeans patch separability + optional supervised diagnostic | No (reads TSV) |
| `config/chitosan_split.toml` | 1e | Experimental split-width Gaussian profile config | No |
| `test/build_labelfree_unit_predictions.jl` | 1f | Build a reproducible label-free 0/1/? prediction TSV from feature/patch TSVs | No (reads TSVs) |
| `test/score_connected_mold_templates.jl` | 2a | Score connected GlcN/GlcNAc molds over global direction/phase/mirror states | No (reads TSVs) |
| `test/generate_connected_mold_templates.jl` | 2a | Generate unary and optional sliding-bond mold TSVs from aligned geometric/proxy sites | No (reads TSV) |
| `test/refine_geometric_mold.jl` | 2a | Label-free grid refinement of acetyl proxy-site geometry using template-evidence bimodality | No (reads TSVs) |
| `test/import_stm_mold_maps.jl` | 2a | Import DFT-STM/LDOS map TSVs into connected unary/bond mold templates | No (reads TSVs) |
| `test/cube_to_stm_maps.jl` | 2a | Sample QE/DFT cube files into aligned STM map TSVs | No (reads cubes) |
| `test/smoke_qe_mold_workflow.jl` | 2a | Synthetic no-QE smoke test for slab/QE-input/frame/cube/import handoffs | No SXM; writes temp files |
| `test/build_initial_chitosan_trimer_xyz.jl` | 2a | Generate deterministic initial X-GlcN-X and X-GlcNAc-X trimer XYZs plus frame-index TSVs | No SXM; writes XYZ/TSV |
| `test/validate_chitosan_trimer_structures.jl` | 2a | Validate generated trimer/slab atom counts, acetyl counts, distances, labels, and frame indices | No SXM; reads XYZ/TSV |
| `test/preflight_qe_mold_inputs.jl` | 2a | Validate prepared QE run directories before Slurm submission | No SXM; reads QE inputs |
| `test/finalize_qe_mold_workflow.jl` | 2a | Convert completed QE relaxed XYZ/cubes into STMFit connected mold templates | No SXM; reads QE outputs |
| `test/prepare_qe_mold_inputs.jl` | 2a | Generate QE relax/SCF/LDOS inputs from vetted slab+trimer XYZ files | No (reads XYZ) |
| `test/build_qe_slab_trimer_xyz.jl` | 2a | Assemble a reproducible Cu(100) slab below a supplied oriented trimer XYZ | No (reads XYZ) |
| `test/extract_qe_relaxed_xyz.jl` | 2a | Extract final relaxed coordinates and cell from QE `pw.x` output | No (reads QE output) |
| `test/update_qe_positions_from_xyz.jl` | 2a | Replace QE SCF `ATOMIC_POSITIONS` using relaxed XYZ coordinates | No (reads QE input + XYZ) |
| `test/extract_qe_mold_frame.jl` | 2a | Extract central-unit origin/t-axis/u-axis from a relaxed slab+trimer XYZ | No (reads XYZ) |
| `test/validate_connected_molds.jl` | 2a | Validate connected mold files, required combinations, and patch/template dimensions | No (reads TSVs) |
| `test/grade_unit_assignment.jl` | 0 | Grade predictions vs truth (4 alignments, 2 conventions) | No (reads TSVs) |

## Ground truth file

`benchmarks/chitosan_240817_unit_sequences.tsv` — one row per file, columns:
`file`, `sequence` (ordered 0/1 along t_nm increasing), `quality`, `target_N`,
`notes`. **Evaluation-only**: never read by the fitter.

## Commands

```bash
# Phase 1: extract per-lobe features (re-runs the fit, ~10-15 min/file)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. \
    test/extract_lobe_features.jl \
    --config config/chitosan.toml \
    --out results/unit_separability/lobe_features.tsv

# Phase 1: separability analysis (label-free)
julia --project=. test/analyze_unit_separability.jl \
    --features results/unit_separability/lobe_features.tsv \
    --out results/unit_separability

# Phase 1: separability analysis (with truth cross-evaluation)
julia --project=. test/analyze_unit_separability.jl \
    --features results/unit_separability/lobe_features.tsv \
    --truth benchmarks/chitosan_240817_unit_sequences.tsv \
    --out results/unit_separability

# Phase 1c: local/envelope-corrected features, no composition prior
julia --project=. test/augment_lobe_local_features.jl \
    --features results/unit_separability/lobe_features.tsv \
    --out results/unit_separability/lobe_features_local.tsv

# Phase 1d: aligned raw/residual patches
STMFIT_DATA_DIR=/path/to/data julia --project=. \
    test/extract_lobe_patches.jl \
    --features results/unit_separability/lobe_features.tsv \
    --out results/unit_separability/lobe_patches.tsv

# Phase 1d: patch PCA/kmeans, label-free unless --truth is supplied
julia --project=. test/analyze_lobe_patches.jl \
    --patches results/unit_separability/lobe_patches.tsv \
    --prefix res_p \
    --out results/unit_separability/patch_analysis_residual

# Phase 1e: split-width Gaussian features at fixed batch-selected N
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. \
    test/extract_lobe_features.jl \
    --config config/chitosan_split.toml \
    --selected-summary results/best_plots/summary_overlap060_hard.tsv \
    --manifest benchmarks/chitosan_240817.toml \
    --primary-only \
    --out results/unit_separability/lobe_features_selectedN_primary_split.tsv

julia --project=. test/analyze_unit_separability.jl \
    --features results/unit_separability/lobe_features_selectedN_primary_split.tsv \
    --features-list skew_ratio \
    --out results/unit_separability/selectedN_primary_split_skew

# Phase 1f: portable label-free predictor from frozen feature/patch TSVs.
# This does not read truth/control sequences; grade only after the TSV is written.
julia --project=. test/build_labelfree_unit_predictions.jl \
    --features results/unit_separability/lobe_features_selectedN_primary_local.tsv \
    --split-features results/unit_separability/lobe_features_selectedN_primary_split.tsv \
    --patches results/unit_separability/lobe_patches_selectedN_primary_17x17_bwd.tsv \
    --out results/unit_assignment/labelfree_unit_predictions.tsv \
    --seeds 20 \
    --interactions
```

## Challenger candidate manifest (T0 firewall + baseline contract)

The label-free challenger lane defined in
`.omo/plans/improve-unit-assignment-benchmark.md` is gated by a provenance-only
manifest at `config/unit_assignment_candidate.toml` and a dedicated checker at
`test/check_unit_assignment_candidate_manifest.jl`. The manifest is **not** a
model, **not** a selection config, and **not** a fit/calibration artifact. It
records the frozen provenance, exact feature lists, equal view weights, seeds,
bootstrap count, date parser rule, constant-current physical policy,
leave-date-out gate, common real gates, ranking rule, and promotion thresholds
so the challenger cannot silently retune any of them after the freeze.

### Firewall

The benchmark control motif/encoding, the benchmark truth/count column names,
and benchmark truth/grade paths may appear **only** inside the manifest's
dedicated `[grader_only]` section. They must never enter fitting, feature
construction, candidate selection, confidence, abstention, or calibration. The
checker conservatively rejects TOML multiline strings and requires exactly one
real top-level `[grader_only]` table. Its firewall has three complementary views:
it scans exact source bytes, it recursively scans parsed non-grader TOML
keys and string values case-insensitively after TOML decoding, and it scans
comments outside the `[grader_only]` table for the same forbidden token and path
vocabulary. The semantic walk includes strings nested in arrays, inline tables,
and ordinary nested tables, so four- or eight-digit Unicode escapes cannot hide a
forbidden token or path. Comments outside `[grader_only]` receive the same
treatment, so an inert comment cannot carry a forbidden reference either. Only
the sole root `[grader_only]` table is omitted from semantic inspection. The
distinct grader-only
denominator manifest `benchmarks/chitosan_6mer_counting_confirmed.toml` (fixing
the `145 files / 870 control positions` denominator) is referenced by field name
only from `[denominator]`; its path lives exclusively in `[grader_only]`.

### Lifecycle and hash binding

`grade_status` moves `locked` (T0) → `frozen_once` (T7) → `graded` (T8). While
`locked`, `provenance.status = "pending"` and no `frozen_hash` may be declared.
T7 hashes the exact non-grader source bytes, preserving comments, formatting,
and line endings. In this contract, “canonical hash” means precisely that exact
non-grader source-byte projection; it does not mean parsing and semantically
canonicalizing TOML. The projection excludes only bytes belonging to the one real
`[grader_only]` table plus the real `candidate.frozen_hash` and lifecycle-only
`candidate.grade_status` assignments. It writes that digest as
`candidate.frozen_hash` and switches to `frozen_once`; the status exclusion lets
the declared `frozen_once -> graded` transition retain the same digest. The
checker requires frozen provenance and a matching exact-source digest in both
states. In addition, every mandatory artifact provenance field—feature TSVs,
constant-current cubes, generated maps, molds, and configs—must hold a real
lowercase 64-hex SHA-256; pending, missing, uppercase, malformed non-hex, and
wrong-length values are rejected. This validates the bytes and state presented
on each invocation; it is not stateless historical proof of prior manifest
contents. A future scientific hypothesis requires a new versioned plan and
candidate manifest rather than an edit to frozen source.

### Verification

```bash
# Contract + firewall + hash-binding check (locked state at T0)
julia --project=. test/check_unit_assignment_candidate_manifest.jl \
    --config config/unit_assignment_candidate.toml --expect-locked

# Baseline firewall (existing unknown-production scripts) + T0 contract suite
julia --project=. test/test_unit_assignment_candidate_manifest.jl
```

The locked candidate declares the constant-current mean-height policy
(`0.50 nm`) and fixed sensitivity bracket (`0.40:0.05:0.60 nm`), the
leave-one-date-out rule (parse exactly one leading `YYYYMMDD` token; missing or
ambiguous fails), the 1-vs-2 identifiability gate (every held-out date fold
positive on per-lobe log-likelihood improvement AND scan-bootstrap 95% lower
confidence bound above zero), the common real no-truth gates, the no-truth
ranking order, and the final post-hoc promotion thresholds on the fixed
`145 / 870` denominator: `honest_correct >= 677`,
`physical_accuracy_classified >= 78.9%`, and `exact_chains >= 18`. Coverage,
short-`N` positions, and extra lobes remain separate report fields and must not
be hidden by the headline.

## Unknown chitosan sequence production

<!-- UNKNOWN-CHITOSAN-WORKFLOW:START -->

For an unlabeled chain, do not call the report/grading scripts while generating
assignments. Use the production wrapper to write both fixed label-free profiles,
then validate, plot, and queue chains for visual review:

```bash
julia --project=. test/run_unknown_unit_assignment.jl \
    --features results/unit_assignment/<sample>_features_local.tsv \
    --split-features results/unit_assignment/<sample>_features_split.tsv \
    --patches results/unit_assignment/<sample>_patches_bwd.tsv \
    --profile default \
    --outdir results/unit_assignment/<sample>_unknown

julia --project=. test/validate_unit_predictions.jl \
    --predictions results/unit_assignment/<sample>_unknown/predictions_base_split_log_skew.tsv \
    --features results/unit_assignment/<sample>_features_local.tsv

python3 test/plot_unit_assignment.py \
    --features results/unit_assignment/<sample>_features_local.tsv \
    --predictions results/unit_assignment/<sample>_unknown/predictions_base_split_log_skew.tsv \
    --out-dir results/unit_assignment/<sample>_unknown/plots \
    --mode all

julia --project=. test/summarize_unknown_unit_qc.jl \
    --predictions results/unit_assignment/<sample>_unknown/predictions_base_split_log_skew.tsv \
    --plots-dir results/unit_assignment/<sample>_unknown/plots \
    --out results/unit_assignment/<sample>_unknown/review_queue.tsv
```

Production artifacts are `predictions_*.tsv`, `summary.tsv`, `manifest.tsv`,
per-profile validation logs, plots, and `review_queue.tsv`. The review queue is
based only on prediction fields, confidence, lobe contiguity, optional view
coverage, N outliers within the run, and missing plot files.

<!-- UNKNOWN-CHITOSAN-WORKFLOW:END -->

### Challenger terminal status (T9 closure)

<!-- T9-TERMINAL-STATUS:START -->

The constant-current T3 lane is terminal `BLOCKED`: accepted GlcNAc cubes have
multiple nonunique isovalue branches, and no branch-selection policy was
predeclared. The hierarchical T5 identifiability gate passed all 13 held-out
dates and the 500-seed scan bootstrap, with a 95% lower bound of
`0.55621530226471905`. T6 integration and leakage checks were confirmed after
the URI-path and partial-view QC corrections.

T7 therefore ended as `NO_ELIGIBLE_CHALLENGER`. The hierarchical lane lacks
durable forward/backward evidence and the required common per-lobe audit table,
so no candidate was frozen. T8 is `SKIPPED_NO_ELIGIBLE_CHALLENGER`; the grader
invocation count was zero, and the one-shot grade budget remains unused.
`config/unit_assignment_candidate.toml` remains
`grade_status = "locked"` with `provenance.status = "pending"` and no frozen
hash. Its unchanged file SHA-256 is
`9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`.

Because no new grade ran, the existing benchmark headline remains historical
and current: `78.9%` classified physical accuracy and `17/145` exact chains.
The classified percentage uses only classified positions; it is not the
fixed-denominator honest view, which keeps missing/abstained benchmark positions
in the denominator. This closure does not establish promotion or new benchmark
validation.

## Exploration state of the art (2026-08-02, promoted)

The promotion bar (78.9% classified physical accuracy / 18 exact chains /
677 fixed-denominator honest) is MET by the label-free champion below,
validated by the cross-validated Fisher empirical mold (half-split: 66.3%
per-lobe, no overfitting).

Post-closure exploration (journal 2026-08-01/02, sections 8a-8k) produced a
label-free candidate set that improves the exact-chain record without reaching
the promotion bar. The best candidate is a **soft vote (mean of probabilities)
of the k-means 4-view and the GMM 1-view + per-channel constant-current
margin**:

| Metric | Soft candidate | Promotion bar |
|---|---|---|
| Classified physical accuracy | **79.3% (677/854)** | ≥ 78.9% ✓ |
| Fixed-denominator honest | **677/854** | ≥ 677 ✓ |
| Exact chains | **36/145** | ≥ 18 ✓ |
| Fisher mold CV (half-split) | 66.3% per-lobe, no overfit | — |

**Raw-input application prerequisite (2026-09-16):** the champion consumes the
precomputed `patch_u_asym` feature as well as the base, split, and backward
features. Its original producer is not present in the current repository or
local result scripts. The journal describes a first moment of a 9×9 residual
patch but does not specify the complete normalization/axis convention. Do not
replace it with a new descriptor or a 17×17 approximation. Recover the original
producer before claiming a raw-SXM-to-champion run on new 10–20mer data. The
builder now checks the required feature columns and matching 17×17 mold-patch
keys before starting computation. This does not change the accepted benchmark
or prevent regrading the frozen prediction below.

Frozen prediction: `results/unit_assignment/best_labelfree_cc_soft_20260802.tsv`
(145 files, 854 classified lobes, label-free construction; post-hoc grade only).
Key building blocks (all label-free):

- Adaptive-contour (constant-current isosurface) DFT-STM molds beat the
  constant-height molds as a GMM feature: per-channel forward/backward margins
  (not averaged) + self-training 2 iters give 78.4% / 38 exact alone.
- Soft voting between k-means and the GMM beats hard voting (always between
  components) and exceeds both components on the robust metric (672 > 666/668).
- The empirical mold is the **Fisher discriminant** (label-free GMM-on-PCA10
  cluster means + regularized noise covariance): the optimal linear score,
  66.3% per-lobe under half-split CV with no overfitting
  (`test/lib/empirical_fisher_mold.py`).
- `test/build_cc_soft_champion.py` + the mold builders rebuild templates,
  margins, both predictors, and the soft vote. The builder is prediction-only:
  it keeps every input lobe, accepts no control data, and never starts grading.
  Use Julia 1.13 (`--julia /path/to/julia` if needed), an explicit new `--workdir`,
  and a new `--out` path. Missing component predictions remain `?`; confidence
  is the existing soft-vote margin, not a calibrated probability.
  Benchmark reporting is a separate post-hoc command on the frozen benchmark
  prediction file:

  ```bash
  julia --project=. test/report_unit_assignment_benchmark.jl --full145-own-n \
      --profile champion=results/unit_assignment/best_labelfree_cc_soft_20260802.tsv \
      --outdir results/unit_assignment/champion_regression
  ```
- The 25x25 (0.48 nm half) context mold fails (54.4%): neighbor lobes dominate
  the NCC and dilute the central chemical signal; 17x17 (0.32 nm) is optimal.
- All abstention, height-bracket, registration, and seed-scaling levers are
  neutral or negative (journal sections 8d-8j).

The champion is promoted as the new unit-assignment state of the art:
label-free construction (empirical mold = unsupervised patch clustering with
the physical amplitude mapping; no truth, sequence, or composition prior),
post-hoc grading only. The empirical mold generalizes under half-split
cross-validation, confirming the bar is met without data re-use bias.

Final review hardened, but did not promote, both diagnostic lanes. Hierarchical
unstable/nonmonotone fits abstain, held-out evaluation rejects them before
scoring, and optional split/backward descriptors now reach fitting through a
records-based pipeline. Prediction provenance binds every consumed primary,
split, and backward artifact by stable role, normalized path, and byte SHA-256,
together with resolved views and model options. Constant-current calibration
scans the full declared interval under the validated
`--isovalue-scan-intervals` policy (default 1024), records
`isovalue_scan_intervals` in provenance, and accepts only one continuous
fixed-support root branch at that declared resolution. Multiple roots or support
discontinuities are rejected rather than resolved by a branch preference. The
accepted GlcNAc cube still violates that uniqueness contract, so this API
behavior does not unblock T3.

The manifest checker rejects multiline TOML strings, requires exactly one real
top-level `[grader_only]` table, and scans forbidden non-grader keys, values,
comments, and paths case-insensitively. Its digest binds every other non-grader
byte, including comments, formatting, and line endings, excluding only real
`candidate.frozen_hash` and lifecycle-only `candidate.grade_status` assignments.
Matching provenance and digest are mandatory in `frozen_once` and `graded`, and
the lifecycle transition does not rehash. This is not stateless historical
proof. The current candidate remains locked and byte-unchanged; these stronger
future-state checks did not freeze it or consume the grade budget.

Hierarchical merged-result publication is also recoverable as one bounded set,
not merely atomic one file at a time. The merged TSV and its gate are staged in
their destination directory under a durable `prepared` marker; the old gate is
backed up first, the merged TSV is installed next, and the new gate is committed
last before a durable `committed` marker. A rerun rolls a prepared transaction
back to the old complete generation or retains a committed new generation and
cleans its sidecars. The protocol replaces destination symlinks rather than
following them and refuses malformed or symlinked markers. It preserves public
paths, output bytes, return fields, and strict shard validation. Its guarantee
is limited to deterministic recovery through the tested protocol states; it
does not cover filesystem corruption, power loss between arbitrary syscalls, or
hostile concurrent mutation of transaction paths.

The constant-current diagnostic builder publishes its full output set—all
map/mask pairs together with the provenance writer—as one recoverable generation
under the same gate-last protocol, so provenance is visible only for a complete
old or complete new generation and a partial install never exposes a mixed set.

### Structured evaluator publication

The structured evaluator prefers Linux `renameat2(RENAME_NOREPLACE)` for its
already-fsynced private stage. Some filesystems reject directory publication
with `EINVAL`; the only fallback cases are numeric `EINVAL`, `ENOSYS`, and
`EOPNOTSUPP`/`ENOTSUP`. Other rename errors, including `EXDEV`, `EIO`, `EPERM`,
and `EACCES`, fail closed and do not enter the fallback.

The fallback exclusively reserves the destination with a mode-`0700` `mkdir`,
records and continuously revalidates its device/inode, and transfers sorted
non-receipt files with no-overwrite hard links. The staged `receipt.toml` is
hard-linked last and is the commit marker: visibility of that receipt makes the
publication ambiguous until destination/stage fsyncs, exact bytes and
inventory, file/directory modes, regular non-symlink single-link identities,
and production-context revalidation all pass. Only then is the result
`committed_verified`.

An interruption before reservation cleans only an identity-verified private
stage. A reservation failure before receipt leaves a no-receipt residue, which
is a collision on retry; failures at or after receipt never remove the
destination. Existing output is accepted only when it is the complete expected
set (nine report files or the receipt-only blocker set), with a mode-`0700`
directory, mode-`0644` regular non-symlink single-link files, and exact bytes.
Empty, extra, missing, symlinked, hard-linked, wrong-mode, or differing output
is never adopted or overwritten.

This receipt-gated fallback is the evaluator's ordinary output transaction; it
is distinct from the parent-owned immutable evidence publisher used for
authority roots and publication evidence.

The public hierarchical two-component EM fit validates all controls before data
or fit work: `n_starts` is a positive non-`Bool` integer representable as `Int`,
`first_seed` is a non-`Bool` integer representable as `Int` satisfying
`0 <= first_seed < n_starts`, and
`tol` and `cov_floor` are positive finite non-`Bool` reals representable as
`Float64`. Validation order is `n_starts`, `first_seed`, the
`first_seed < n_starts` relation, `tol`, then `cov_floor`, so malformed controls
produce parameter-specific `ArgumentError`s before data conversion, finite-row
checks, initialization, or EM fitting.
Defaults and valid custom controls remain deterministic, covariance floors are
enforced, and unstable/nonmonotone fits retain their explicit abstention paths.
The mixture priors remain fixed and equal at `[0.5, 0.5]`; validation adds no
composition prior and does not change assignment semantics.

Hierarchical identifiability is symmetric in the two component labels. If hard
assignments leave either component empty—whether all rows occupy component 1 or
all occupy component 2—the result is one-component evidence and the diagnostic
abstains. Separation, likelihood improvement, or amplitude spread cannot
override that empty-component rule.

`hierarchical_equalprior` remains executable for unknown-production diagnostics
but is diagnostic, not frozen or promoted. A reproducible label-free diagnostic
run is:

```bash
julia --project=. test/run_unknown_unit_assignment.jl \
    --features results/unit_assignment/<sample>_features_local.tsv \
    --profile hierarchical_equalprior \
    --outdir results/unit_assignment/<sample>_hierarchical_diagnostic

julia --project=. test/validate_unit_predictions.jl \
    --predictions results/unit_assignment/<sample>_hierarchical_diagnostic/predictions_hierarchical_equalprior.tsv \
    --features results/unit_assignment/<sample>_features_local.tsv
```

<!-- T9-TERMINAL-STATUS:END -->

```bash
# Merge resumable feature-extraction shards against a reference selected-N TSV.
# Use --ignore-extra only for known stale shard rows outside the current manifest.
julia --project=. test/merge_lobe_feature_shards.jl \
    --reference /tmp/opencode/full145_selectedN_features.tsv \
    --shards <comma-separated-shard-tsvs> \
    --out /tmp/opencode/full145_selectedN_features_split.tsv \
    --ignore-extra

# Phase 2a: connected geometric mold-template decoding, no composition prior
STMFIT_DATA_DIR=/path/to/data julia --project=. \
    test/extract_lobe_patches.jl \
    --features results/unit_separability/lobe_features_selectedN_primary.tsv \
    --out results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --half-nm 0.48 \
    --step-nm 0.08

julia --project=. test/generate_connected_mold_templates.jl \
    --atoms templates/chitosan_geometric_sites.tsv \
    --out templates/chitosan_connected_molds.tsv \
    --bond-out templates/chitosan_connected_bond_molds.tsv \
    --half-nm 0.48 \
    --step-nm 0.08

julia --project=. test/validate_connected_molds.jl \
    --atoms templates/chitosan_geometric_sites.tsv \
    --templates templates/chitosan_connected_molds.tsv \
    --bond-templates templates/chitosan_connected_bond_molds.tsv \
    --patches results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --prefix raw_p \
    --report results/unit_assignment/connected_mold_validation.txt

julia --project=. test/score_connected_mold_templates.jl \
    --patches results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --templates templates/chitosan_connected_molds.tsv \
    --prefix raw_p \
    --template-mode contrast \
    --out results/unit_assignment/connected_mold_predictions.tsv

# Optional Phase 2a refinement: label-free search around the geometric acetyl sites.
julia --project=. test/refine_geometric_mold.jl \
    --patches results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --sites templates/chitosan_geometric_sites.tsv \
    --out-sites templates/chitosan_geometric_sites_refined_raw.tsv \
    --report results/unit_assignment/geometric_mold_refinement_raw.tsv \
    --prefix raw_p

julia --project=. test/generate_connected_mold_templates.jl \
    --atoms templates/chitosan_geometric_sites_refined_raw.tsv \
    --out templates/chitosan_connected_molds_refined_raw.tsv \
    --bond-out templates/chitosan_connected_bond_molds_refined_raw.tsv \
    --half-nm 0.48 \
    --step-nm 0.08

julia --project=. test/score_connected_mold_templates.jl \
    --patches results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --templates templates/chitosan_connected_molds_refined_raw.tsv \
    --prefix raw_p \
    --template-mode contrast \
    --out results/unit_assignment/geometric_mold_predictions_refined_raw.tsv

# Future ideal path: import DFT-STM/LDOS maps instead of geometric proxy sites.
julia --project=. test/build_initial_chitosan_trimer_xyz.jl \
    --out-dir hpc/qe_molds

julia --project=. test/validate_chitosan_trimer_structures.jl \
    --dir hpc/qe_molds \
    --out hpc/qe_molds/structure_validation.tsv

julia --project=. test/build_qe_slab_trimer_xyz.jl \
    --molecule hpc/qe_molds/glcn_central_trimer.xyz \
    --out hpc/qe_molds/glcn_central_trimer_slab.xyz \
    --metadata hpc/qe_molds/glcn_central_trimer_slab_meta.tsv \
    --nx 8 --ny 8 --layers 4 \
    --center-indices 12,13,14,15,16,17 \
    --height-above-top 2.6 --vacuum 18.0

julia --project=. test/prepare_qe_mold_inputs.jl \
    --xyz hpc/qe_molds/glcn_central_trimer_slab.xyz \
    --cell-metadata hpc/qe_molds/glcn_central_trimer_slab_meta.tsv \
    --out-dir qe/glcn \
    --prefix glcn_central \
    --fix-below-z ZCUT \
    --sample-bias-ev -0.3

# Repeat for the GlcNAc-central trimer (required before preflighting both dirs):
julia --project=. test/build_qe_slab_trimer_xyz.jl \
    --molecule hpc/qe_molds/glcnac_central_trimer.xyz \
    --out hpc/qe_molds/glcnac_central_trimer_slab.xyz \
    --metadata hpc/qe_molds/glcnac_central_trimer_slab_meta.tsv \
    --nx 8 --ny 8 --layers 4 \
    --center-indices 12,13,14,15,16,17 \
    --height-above-top 2.6 --vacuum 18.0

julia --project=. test/prepare_qe_mold_inputs.jl \
    --xyz hpc/qe_molds/glcnac_central_trimer_slab.xyz \
    --cell-metadata hpc/qe_molds/glcnac_central_trimer_slab_meta.tsv \
    --out-dir qe/glcnac \
    --prefix glcnac_central \
    --fix-below-z ZCUT \
    --sample-bias-ev -0.3

julia --project=. test/preflight_qe_mold_inputs.jl \
    --dir qe/glcn \
    --dir qe/glcnac \
    --out hpc/qe_molds/qe_input_preflight.tsv \
    --max-total-tasks 8 \
    --sequential

bash hpc/submit_qe_molds.sh --watch --sequential

julia --project=. test/finalize_qe_mold_workflow.jl \
    --height-nm HEIGHT \
    --glcn-dir qe/glcn \
    --glcnac-dir qe/glcnac

# The generated run_qe_mold.sbatch performs these two handoff steps automatically.
julia --project=. test/extract_qe_relaxed_xyz.jl \
    --qe-out qe/glcn/glcn_central_relax.out \
    --out qe/glcn/glcn_central_relaxed.xyz \
    --metadata qe/glcn/glcn_central_relaxed_meta.tsv

julia --project=. test/update_qe_positions_from_xyz.jl \
    --input qe/glcn/pw_scf.in \
    --xyz qe/glcn/glcn_central_relaxed.xyz \
    --out qe/glcn/pw_scf_relaxed.in

julia --project=. test/extract_qe_mold_frame.jl \
    --xyz qe/glcn/glcn_central_relaxed.xyz \
    --origin-indices I,J \
    --axis-from I --axis-to J --plane-index K \
    --height-nm HEIGHT \
    --out qe/glcn/frame.tsv

julia --project=. test/cube_to_stm_maps.jl \
    --cube 0:qe/glcn/glcn_central_ldos.cube \
    --frame 0:qe/glcn/frame.tsv \
    --cube 1:qe/glcnac/glcnac_central_ldos.cube \
    --frame 1:qe/glcnac/frame.tsv \
    --cube-units bohr \
    --out templates/chitosan_stm_maps.tsv

julia --project=. test/import_stm_mold_maps.jl \
    --maps templates/chitosan_stm_maps.tsv \
    --out templates/chitosan_connected_molds_stm.tsv \
    --bond-out templates/chitosan_connected_bond_molds_stm.tsv \
    --half-nm 0.48 \
    --step-nm 0.08

julia --project=. test/score_connected_mold_templates.jl \
    --patches results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --templates templates/chitosan_connected_molds_stm.tsv \
    --bond-templates templates/chitosan_connected_bond_molds_stm.tsv \
    --prefix raw_p \
    --template-mode contrast \
    --out results/unit_assignment/stm_mold_predictions.tsv

# Phase 0: grading (once predictions exist)
julia --project=. test/grade_unit_assignment.jl \
    --predictions results/unit_assignment/assigned_sequences.tsv \
    --truth benchmarks/chitosan_240817_unit_sequences.tsv \
    --out results/benchmark_grades/unit_assignment.tsv

# Frozen-profile subset report (historical 35-file / 210-lobe subset)
julia --project=. test/report_unit_assignment_benchmark.jl

# Full145 denominator check/report mode once full145 predictions exist
julia --project=. test/report_unit_assignment_benchmark.jl --full145

# Full145 own-N report for predictions generated at each file's label-free N_selected
julia --project=. test/report_unit_assignment_benchmark.jl --full145-own-n \
    --profile selectedN_local=/tmp/opencode/full145_selectedN_labelfree_local_predictions.tsv=selectedN_local \
    --outdir /tmp/opencode/unit_report_full145_own_n_local
```

Prediction TSVs may use `0`, `1`, or `?` in the prediction column. `?` is an
explicit abstention. The grader reports two views:

- **Diagnostic classified accuracy** excludes `?` from accuracy/confusion counts
  and reports classified coverage as `classified/possible`.
- **Honest abstention view** reports only `correctly assigned` and `uncertain`,
  where `uncertain` includes explicit `?` plus any benchmark-detected wrong
  assignment. These two numbers sum to 100% over the graded positions.

The second view is post-hoc benchmark reporting only. It does not use truth to
choose thresholds or alter predictions; it prevents residual benchmark errors from
being presented as honest assignments.

The full 0/1/? benchmark scope is the same 145 confirmed 6mer files as
`benchmarks/chitosan_6mer_counting_confirmed.toml`. For external grading only,
the control sequence is `NKNNKN`: `010010` if `0=N` and `1=K`, or `101101` under
the flipped identity convention. This information must never enter fitting,
selection, assignment, threshold choice, abstention rules, composition priors, or
method calibration. The objective is a method robust enough to extrapolate to
unknown systems without a benchmark-specific control.

The current `test/report_unit_assignment_benchmark.jl` harness defaults to a
frozen-profile **subset** report, not the final full145 benchmark headline. It
runs the grader on the historical 240817 forced, balanced-abstention, and
strict-emitted-error profiles, writes per-profile grade TSVs under
`results/unit_assignment/benchmark_report/grades/`, and consolidates:

- `summary.tsv`: coverage, classified accuracy, honest correct/uncertain, emitted
  errors, and exact-chain counts for each frozen profile.
- `lobe_position_errors.tsv`: per-lobe-position concentration of errors and `?`
  calls, useful for diagnostics but not for choosing new thresholds.
- `report.md`: human-readable summary with the subset denominator (`35` primary
  files, `210` lobes) explicitly stated.

The script is deliberately report-only. It does not sweep thresholds, create a new
profile, or use truth to alter predictions. Its strict `--full145` mode generates
the grader-only `NKNNKN` control TSV from the 145-file manifest and enforces the
145-file / 870-lobe prediction denominator, so passing current 35-file frozen
profiles to that mode fails at a coverage preflight instead of producing a
misleading full-benchmark number. The explicit `--full145-own-n` mode is for
prediction profiles produced at each file's label-free `N_selected`: all 145 files
must be present, lobe rows must be contiguous per file, missing 6mer control
positions count uncertain in the honest view, and extra predicted lobes are
reported but not aligned to the 6-position external control.

## Joint proxy-mold diagnostic inference

The experimental joint posterior path is separate from production
`N_selected`. It preserves every valid candidate N, calibrates count and type
confidence on deterministic synthetic chains only, and permits explicit `?`
abstention. Chemistry evidence never changes the count posterior.

```bash
julia --project=. test/calibrate_joint_proxy_molds.jl \
  --config config/joint_proxy_molds.toml --seed 20260721 --cases 100 --fast \
  --out config/joint_proxy_calibration_dft_m030_h050_v1.toml

julia --project=. test/infer_joint_proxy_molds.jl \
  --config config/joint_proxy_molds.toml \
  --data-dir /path/to/20240817_LHe_Cu100 \
  --calibration config/joint_proxy_calibration_dft_m030_h050_v1.toml \
  --files 240817_017.sxm,240817_019.sxm \
  --outdir /tmp/opencode/joint_proxy_real_smoke

julia --project=. test/validate_joint_proxy_predictions.jl \
  --artifacts /tmp/opencode/joint_proxy_real_smoke \
  --calibration config/joint_proxy_calibration_dft_m030_h050_v1.toml
```

The output directory contains `candidate_n.tsv` (posterior over N and paired-view
diagnostics), `candidate_lobes.tsv` (per-candidate lobe geometry and type
probabilities), `predictions.tsv` (hard `0/1/?` calls for the selected candidate),
`chain_summary.tsv` (coverage, confidence, and abstention), and
`run_manifest.toml` (input list plus config/source/payload hashes). Input lists
accept plain `.sxm` basenames only; truth/benchmark metadata are rejected.

The active registry combines the tracked geometric provider with the versioned
`stm_dft_v1` source at equal weight. The DFT source is bound to corrected
`-0.300 V`, `0.50 nm` cube/map/template hashes; preliminary sources remain
disabled diagnostic history. Synthetic confidence calibration still cannot
establish real chemical correctness. Paired synthetic A/B checks establish only
that activation does not regress the existing recovery/abstention metrics. Never
choose or retune this provider from benchmark sequence accuracy.

## HPC parallelism

`extract_lobe_features.jl` supports `--chunk I/N` for job-array parallelism:

```bash
# On Raven (≤ 8 CPUs total, e.g. 4 chunks × 2 CPUs)
STMFIT_DATA_DIR=/data julia -t 2 --project=. \
    test/extract_lobe_features.jl \
    --config config/chitosan.toml \
    --chunk 1/4 \
    --out results/unit_separability/lobe_features_chunk01.tsv
# Then concatenate chunks manually (cat *_chunk*.tsv > lobe_features.tsv)
```

## Current status

- **Phase 0**: `grade_unit_assignment.jl` implemented, plus
  `report_unit_assignment_benchmark.jl` for the frozen-profile 0/1/? subset
  report. The full benchmark to process/score is the 145-file confirmed 6mer set;
  `NKNNKN` is the external control sequence for all benchmark chains. Existing
  `178/210`, `154/171`, and `101/106` metrics are historical 35-file subset
  grades and should not be presented as the full benchmark result. Use
  `report_unit_assignment_benchmark.jl --full145` only with strict six-lobe
  prediction profiles that actually contain all 145 benchmark chains. Use
  `--full145-own-n` for profiles generated at each file's label-free `N_selected`;
  current frozen subset profiles still fail the full145 preflight because they
  contain only 35 chains / 210 lobe rows.
- **Phase 1**: `extract_lobe_features.jl` + `analyze_unit_separability.jl`
  implemented and run on the corrected primary benchmark set using batch
  `N_selected` (39 files, 234 lobes). Gaussian features are **strongly
  bimodal**: `ΔBIC(k=1-k=2) = +184.1`.
- **Phase 1b**: residual-feature extraction with full baseline/tilt subtraction
  was tried and removed after it weakened the evidence to
  `ΔBIC(k=1-k=2) = +4.3` (weakly bimodal). Those residual features added more
  noise than signal and are no longer part of the maintained workflow.
- **Phase 1c**: `augment_lobe_local_features.jl` implemented and run. Local
  prominence/envelope features remain bimodal (`ΔBIC = +144.0`) but do not solve
  the chemistry: label-free assignment gives 69.2% physical accuracy and 1/39
  exact sequences against the withheld `NKNNKN` diagnostic truth.
- **Phase 1d**: aligned patch extraction and PCA/kmeans diagnostics implemented
  and run on the same 39 primary files (234 patches, 9×9 grid). Raw patches are
  bimodal (`ΔBIC = +66.2`) but supervised train/test accuracy is only 51.4%.
  Residual patches are weakly bimodal (`ΔBIC = +20.3`) with 61.1% supervised
  train/test accuracy. Current patches therefore do not yet carry enough
  generalizable unit-identity signal.
- **Phase 1e**: split-width Gaussian forward model implemented as
  `peak_profile = "split"` in `config/chitosan_split.toml` and run on the 39
  primary benchmark files at fixed batch-selected N. All files refit at N=6.
  Split improved fixed-N GCV on 36/39 files (median relative ΔGCV = -10.7%, mean
  = -12.6%) and `skew_ratio` is strongly bimodal (`ΔBIC = +294.8`). However,
  external diagnostic grading shows this asymmetry is not GlcNAc/GlcN identity:
  `skew_ratio` AUC is 0.477 (inverse 0.523), skew-only assignment is 48.3%
  physical / 63.7% oracle with 0/39 exact sequences, and Gaussian+skew degrades
  to 60.7% physical / 70.1% oracle. Conclusion: split captures real local shape
  asymmetry useful to the fit, but not the chemical unit label.
- **Phase 1f**: `build_labelfree_unit_predictions.jl` now provides a portable
  label-free baseline from existing feature/patch TSVs. It clusters per-file
  standardized feature views over multiple k-means seeds and maps the
  higher-amplitude cluster to GlcNAc (1). On the current 39-file feature
  artifacts, the three default views (`BASE+bwd_neg_com_t`,
  `BASE+bwd_neg_diag45`, `BASE+split_log_skew`) produce 234 prediction rows and
  grade at 170/210 = 81.0% on the historical 35-file subset. This is a
  sanity-checkable generator, not a replacement for the frozen best ensemble
  artifacts. Full145 selected-`N` profiles have now been generated from the
  current label-free counting summary without using the `NKNNKN` control: 145
  files / 863 prediction rows, with 13 missing control positions from 10 short-N
  files and 6 extra predicted lobes from 6 N=7 files. The base/local-feature
  baseline grades at 671/857 = 78.3% classified physical accuracy; the honest
  view is 671/870 = 77.1% correct plus 199/870 uncertain, with 16/145 exact
  chains. Adding all default backward descriptors directly does not improve the
  honest headline (671/870 correct, 5/145 exact), but a simple selected-`N`
  three-view ensemble over `BASE`, `BASE+bwd_neg_com_t`, and
  `BASE+bwd_neg_diag45` reaches 676/857 = 78.9% classified physical accuracy,
  honest 676/870 = 77.7% correct + 194/870 uncertain, and 7/145 exact chains.
  The selected-`N` split-width refit has since been completed for all 145 files
  (863 rows; 8 old rows from excluded `240310_Cu100009.sxm` ignored during the
  merge). Adding `split_log_skew` does not raise the lobe-correct headline:
  the default split+bwd views tie the same 676/870 honest correct result, and a
  four-view `BASE+bwd+split` variant falls slightly to 675/870. However,
  `BASE+split_log_skew` preserves 676/870 honest correct while increasing exact
  chains to 17/145, so it is the current best exact-chain selected-`N` full145
  profile. Confidence/agreement abstention variants raise classified accuracy
  only marginally (to about 79.5-79.6%) while reducing honest correct lobes, and
  backward-only views collapse to about 50% physical accuracy.
- **Workflow hardening**: `merge_lobe_feature_shards.jl` replaces the ad hoc
  Julia one-liners previously used to combine timeout-limited split-width shards.
  It orders rows by a reference selected-`N` TSV, rejects missing and duplicate
  `(file,lobe)` keys, reports extra stale rows, and writes the merged TSV only
  after coverage validation. This is the supported path for resumable split/full145
  feature extraction merges.
- **Phase 2a**: connected mold-template decoder implemented and switched to a
  manual geometric proxy-site source (`templates/chitosan_geometric_sites.tsv`).
  It tests global direction/phase/mirror states and applies oriented GlcN/GlcNAc
  patch templates plus optional sliding bond templates for transitions
  `00/01/10/11`, without truth labels or composition constraints. The first
  geometric diagnostic uses 13×13 raw patches (`±0.48 nm`) and contrast scoring.
  It validates technically but does not solve unit identity: raw contrast gives
  60.3% physical / 63.7% oracle with 0/39 exact sequences; adding bond templates
  does not improve the result. Full-template scoring reaches 67.9% oracle with
  4/39 exact sequences, but the physical 0↔1 mapping fails, so it is not a valid
  label-free assignment.
- **Phase 2a refinement**: `refine_geometric_mold.jl` implemented. Raw-patch
  refinement selected a shifted/narrower acetyl geometry by a fully label-free
  ΔBIC objective (`ΔBIC=417.3`, clusters 74/160). After amplitude-based physical
  remapping, the refined raw mold gives 67.9% physical / 72.2% oracle and 0/39
  physical exact sequences (1/39 oracle exact). Adding sliding bonds slightly
  worsens to 67.1% physical. Residual-patch refinement gives a stronger
  label-free ΔBIC (`625.8`) but poor physical grading (44.0%), so that residual
  mode is non-chemical. Conclusion: label-free refinement improves the raw mold
  but still does not reach a robust deacetylation map; DFT-STM or a better
  observable is still needed before forward-model integration.
- **DFT-STM mold path**: `import_stm_mold_maps.jl` implemented and validated on
  synthetic GlcN/GlcNAc maps. It converts long-form simulated STM/LDOS maps into
  the same connected template format as the geometric molds, with beta-(1->4)
  parity/mirror variants and optional bond-map support. Real DFT-STM maps are
  still needed before this path can be scientifically graded.
- **Validation**: `test/validate_connected_molds.jl` checks the connected-mold
  file chain before decoding: geometric/proxy site columns, required unary/bond
  combinations, and patch/template pixel-count compatibility. It does not read
  truth labels.
- **Current best conservative label-free output** (updated Jun 29): after the
  forward-only LOFO audit, the backward Z scan was added via
  `extract_lobe_patches_bwd.jl`. The best single descriptor is now
  `bwd_neg_com_t` at **84.3% LOFO / 11 exact chains**. A confirmed 20-seed
  three-view label-free ensemble (`BASE+bwd_neg_com_t`, `BASE+bwd_neg_diag45`,
  `BASE+split_log_skew`) gives the best historical 35-file subset result:
  **178/210 = 84.8% physical / 84.8% oracle / 14/35 exact**. Graded artifacts:
  `results/unit_assignment/best_labelfree_ensemble3_forced_predictions.tsv` and
  `results/unit_assignment/benchmark_report/grades/forced_ensemble3.tsv`.
- **Current honest high-confidence output**: the same three-view ensemble can emit
  `?` for ambiguous lobes without using truth or a composition prior. The current
  best abstention rule keeps the `forced` ensemble label only when (i) its
  confidence is at least 0.65 and (ii) an independent conservative label-free
  model (`best_labelfree_predictions.tsv`) agrees; otherwise it emits `?`. This
  `agreebase65` variant classifies 171/210 lobes at **154/171 = 90.1%**. In the
  honest two-score benchmark view it gives **154/210 = 73.3% correctly assigned**
  plus **56/210 = 26.7% uncertain**, improving on the old `0.20/0.80` confidence
  band (**151/210 correct + 59/210 uncertain**). More conservative variants are
  available: `agreebase70` gives **153/210 correct + 57/210 uncertain** at 90.5%
  classified accuracy, and `agreebase80` gives **144/210 correct + 66/210
  uncertain** at 91.7% classified accuracy. This is still not a solved binary map,
  but it is the most honest current mode when uncertain lobes are acceptable.
  A separate **strict <5% emitted-error profile** is also available: keep the
  forced label only when confidence is at least 0.875 and both
  `best_labelfree_v3_neg_diag135_predictions.tsv` and
  `best_labelfree_predictions.tsv` agree. This `err05` profile emits only 106/210
  labels, but grades at **101/106 = 95.3%** physical accuracy, i.e. **5/106 =
  4.7%** wrong assignments among emitted labels. In honest two-score form it is
  **101/210 = 48.1% correctly assigned** plus **109/210 = 51.9% uncertain**.
  The canonical consolidated report is
  `results/unit_assignment/benchmark_report/report.md`; generated summary tables
  are intentionally under `results/` and should be regenerated rather than
  committed.
- **Benchmark-exploration candidates** (not canonical): extended 17×17 residual
  descriptors (`hh1_q00_abs + neg_anis`, `patch9_u_asym`) improve post-hoc
  benchmark grades but collapse under LOFO cross-validation (69–66%). Absolute
  backward-height features, parity-canonicalized signed features, backward
  residual recalibration, equal-prior GMM prediction, and 3–6 component GMMs were
  also audited after the 84.3% result and did not beat the confirmed ensemble.
- **Phase 2–5**: planned. Decision point after DFT-STM molds or a stronger
  label-free observable is available.
- **Ground truth/control**: the benchmark control sequence is `NKNNKN`, encoded as
  `010010` or `101101` depending on the 0/1 identity convention. It is used
  exclusively for post-hoc grading, never in fit/selection/assignment-method
  calibration.

## Explicit Julia reconstruction (2026-09-16)

The original `patch_u_asym` producer is irrecoverable. The user approved a new,
explicitly defined Julia 1.13 implementation and complete pipeline, followed by
comparison with the frozen champion. This is `cc_soft_reconstructed_v1`, not a
claim to have recovered the old formula. The old champion and its reported
benchmark remain the historical reference. No variant is selected by benchmark
accuracy or agreement with that reference.

The descriptor `patch_u_asym_reconstructed` uses the backward residual patch
sampled on a 9×9 grid (half-width 0.32 nm, step 0.08 nm), normalized by the
existing extractor as `(p - median(p)) / std(p)` with the sample standard
deviation. It is not normalized again. With `u` transverse to the fitted chain,

```math
A_u = \frac{\sum_{u,t} \operatorname{sign}(u) p(u,t)}{\sum_{u,t}|p(u,t)|}.
```

The center row has weight zero. Pixel serialization follows the extractor's
`for u in coords, t in coords`: `t` varies fastest. This is deliberately tested
rather than inheriting `enrich_unit_features.py`'s conflicting array-axis comment.
Any nonfinite pixel or L1 mass ≤ `descriptor_zero_l1` yields `NA` with a reason;
no lobe is dropped. This is half-plane asymmetry, **not** the u-weighted first
moment mentioned elsewhere in the journal. The choice is explicit and made
before comparing reconstructed predictions.

Settings are in `config/unit_assignment_reconstructed.toml`. Synthetic tests
check the known u/t orientation, reflection sign, positive-scale invariance,
center-row behavior, and invalid-data handling. Complete binary agreement with
a frozen output would establish agreement on those inputs only; matching three
aggregate benchmark metrics is not a proof of formula identity.


### Native production command

Use Julia **1.13** (`julia --version` must show 1.13.x). The root Manifest is
not hand-edited. The production command launches only the same Julia executable:

```bash
GKSwstype=100 julia -t 4 --project=. test/run_reconstructed_chitosan.jl \
  --data-dir /path/to/unknown_sxm_directory \
  --count-config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --config config/unit_assignment_reconstructed.toml \
  --cube0 qe/glcn_restart5/glcn_central_ldos.cube \
  --frame0 qe/glcn_restart5/frame.tsv \
  --cube1 qe/glcnac/glcnac_central_ldos.cube \
  --frame1 qe/glcnac/frame.tsv \
  --outdir results/my_reconstructed_run
```

Add `--dry-run` first. The output directory must not exist. Raw files may be
nested, but their basenames must be unique. No benchmark manifest, expected N,
reference predictions, grading option, control sequence, or composition prior
is accepted by this command. The count config's `[preprocessing]` settings
(`stride`, `flatten`, `smooth_radius_px`) also drive all patch exports. Multi-file
base/split refits use the extractor's existing `--chunk` option in at most four
single-thread subprocesses, within the parent Julia thread budget. A strict
keyed merge runs before cohort-wide assignment. No clustering is split by chunk. The driver's historical triage-file default is
explicitly disabled. The 1D diagnostic stays disabled.

Without cached inputs, the stages are counting, fixed-selected-N base and split
refits, local amplitude features, normalized forward/backward 17×17 patches and
backward 9×9 patches, descriptor, native constant-current molds/Fisher, existing
Julia GMM/k-means predictors, soft vote, validation, QC and maps. `N_selected`
comes from the configured label-free counting policy. Fixed-N feature extraction
refits geometry; it does not restore the exact earlier fitted parameters.
For adaptive-support data, the original summary's `refined_policy` now replays
the already selected support settings from the same count config. The raw scan
and unchanged preprocessing reconstruct the support; a stored length alone is
not used to infer its endpoints. Refit values can still differ from saved fits.

Useful cache options are `--selected-summary` (label-free `filepath,N_selected`
plus `refined_policy` for adaptive-support refits),
`--features`, `--split-features`, `--patches-fwd`, `--patches-bwd`,
`--descriptor-patches`, and `--templates`. Cached geometry can skip counting;
its contiguous lobe keys define the previously selected N, or must match an
explicit selected summary. When base or split geometry must be refitted with
an adaptive-support config, keep the original selected summary even if another
geometry table is cached; N-only metadata is insufficient and fails explicitly.
Nonadaptive N-only summaries retain their existing behavior. Every downstream
table must cover exactly the selected keys. A subset is not silently substituted
for a failed chain. Only genuine split-fit geometry is valid for
`--split-features`; all-one skew ratios in a
base fit are not a split-fit cache. A supplied mold table replaces cube/frame
arguments. It must have been generated with the stated settings.

Outputs include:

- `predictions.tsv`: every selected lobe, `0/1/?`, uncalibrated soft-vote margin,
  `probability_1`, invalid reason and explicit reconstructed model name;
- `summary.tsv`: selected N and ordered assignment for each chain;
- `review_queue.tsv`: generic, label-free QC; this is not an accuracy estimate;
- `plots/standalone/*_chain.png`: fitted-coordinate maps with `0 = GlcN`,
  `1 = GlcNAc`, and gray `?`. The current plotter writes overview page 1 as
  `plots/summary_grid.png`, additional numbered pages as needed, and the complete
  file-to-page index in `plots/summary_pages.tsv` (at most eight chains per page);
- per-stage intermediates and `logs/`; a nonzero exit and `failures.tsv` identify
  incomplete runs. Such a run must not be presented as a complete application.

`?` means a required component is unavailable. Low soft-vote margins remain
visible in the TSV/QC; no benchmark-chosen abstention threshold is added.
Clustering and the Fisher fits use the supplied unlabeled cohort. Changing that
cohort can change predictions; a one-file timing run is not the 25-file result.
There is no known sequence for the 10–20mer application, so processing it cannot
establish chemical accuracy.

### First full-cohort execution: blocked at geometry export (2026-09-17)

Raven job `30271312` completed fresh counting for all 25 raw unknown scans:
all count rows have status `ok`, with 222 selected lobes and selected counts
from 5 to 16. This is label-free processing, not validation against a known
sequence. The base feature table contains 204 unique lobe keys across 23 files,
with no duplicate or extra keys. All nine selected lobes are missing for each
of `260215_022.sxm` and `260220_083.sxm`.

The two missing files are exactly those selected by adaptive support rescue.
The extractor used in that job imported only the selected N, then built a new
fit from the original, shorter support. It did not restore the rescued support
context. Both fixed-N refits reported `No chain model fit succeeded`.
The coverage check aborted, rather than silently discarding these 18 lobes.

No cohort predictions, assignment QC/maps, full146 reconstruction, champion
comparison or external grade were produced. `failures.tsv` lists all intended
25 files. The failed outputs and logs were fetched to
`results/reconstructed_cc_soft_v1/raven_20260916_native_v1/`; the Slurm log is
`results/reconstructed_cc_soft_v1/native-30271312.log`. Monitoring is cancelled;
no resubmission has occurred. The support handoff has since been corrected as
requested: `refined_policy` identifies accepted rescue, including its robust-
guard suffixes. The extractor copies the existing
`adaptive_rescue_support_noise_k` and `adaptive_rescue_support_padding_nm`
into both chain configs before a fixed-N fit. Keep/rejected/failed rescue modes
retain the original support; the image preprocessing and other settings do not
change. Unknown or missing adaptive metadata fails before any SXM read/output.
The implementation does not run selection again, infer endpoints from a saved
length, or add a retry with different parameters.

Focused tests check metadata, both profiles/configs and support feasibility
without optimization. Reading the saved 25-file summary preserves all 222
selected lobes and restores rescue for exactly the two affected files.

A subsequent real check on `260215_022.sxm` passes under Julia 1.13: both base and
split refits export exactly the selected nine lobes, with finite geometry and
GCV, using `refined_policy` from the unchanged saved counting summary. Both
circular and elliptical configs reconstruct support **5.552874796122559 nm**,
exactly the saved value, with 9,392 fit-mask pixels. The native geometry/support
verification passes 76 assertions. A fused, preprocessed STM overlay shows the
numbered fitted centers; these numbers are lobe indices, not chemical labels.

Base/split GCV values are 5.63285686e-6 / 1.51047396e-6; elapsed times are
74.49 / 84.94 s (one thread each, run in parallel). The base GCV differs from the
original successful counting fit (3.18253642e-6). This confirms the support
handoff, not recovery of identical fitted parameters. No refit was repeated or
retuned to reduce that difference. The split skew range, 0.500329–1.097235,
includes a value near the configured lower bound 0.5; it is not evidence of
chemical identification.

Artifacts are in `results/reconstructed_cc_soft_v1/support_handoff_260215_022_v1/`
(`features_base.tsv`, `features_split.tsv`, `support_check.tsv`, logs and
`geometry_overlay.png`). The old failed cohort outputs remain untouched.
The second check on **`260220_083.sxm`** also passes: nine base and nine split
lobes, support **5.90756511515228 nm**, **7,693 fit-mask pixels**, and another 76
native assertions. Its base/split wall times are 70.68 / 77.87 s, and GCV values
are 5.44278424e-6 / 2.55453904e-6 (original counting GCV 5.44586561e-6). Split skew
spans **0.500034–1.999874**, near both configured bounds; no bound was changed.
Its separate `support_handoff_260220_083_v1/` directory contains the same outputs,
and the geometry-only overlay was inspected. Neither check identifies chemical
classes, and neither recovered fit is claimed identical to the original fit.

The two base recoveries were merged with the original 204 valid rows into
complete **25-file / 222-lobe** base geometry. The two genuine split caches
provided **18 lobes / 2 files**, leaving **204 split lobes / 23 files** for the
continuation. Native preparation passes 29 assertions, including exact keys,
unchanged source rows and byte-identical original summary/template copies.
Prepared inputs remain in `results/reconstructed_cc_soft_v1/unknown25_resume_v1_inputs/`.
`hpc/resume_reconstructed_unknown25.sbatch` fits only the missing split files,
then passes complete geometry and the original selected summary to the unchanged
native application. The classifiers train on the full supplied 25-chain cohort,
not isolated files or the 23-file refit subset.

### Completed unknown25 application (2026-09-17)

The user-approved continuation, **Raven job 30278010**, completed from tested
source `aa422a3` under Julia 1.13: **exit 0:0, 6 min 16 s**, eight allocated CPUs
(four requested), MaxRSS about 3.70 GiB. It reused the complete base table and
18 cached split rows, fitted the other 204 split rows, and ran assignment/QC/maps
on all 25 chains. No counting/base rerun, parameter change, file omission,
full146 run or grader was included. The failed run and its artifacts remain
unchanged.

Outputs were fetched to
`results/reconstructed_cc_soft_v1/raven_20260917_unknown25_resume_v1/`; the final
application is in `unknown25/`. All 67 fetched hashes match Raven (66 output
files plus Slurm log), excluding raw-input symlink directories. Native local
verification passes **179 assertions**: all 12 intermediate/final lobe tables
have exactly the original 222 keys, the split merge preserves cached/new rows,
base fields remain unchanged through augmentation, and summary/QC/plot coverage
is complete. Reaggregating saved component votes and rerunning QC reproduce the
two saved TSVs byte-for-byte. No fit or classifier was rerun locally.

| Output diagnostic | Observed result |
|---|---|
| Selected scans / lobes | 25 / 222; per-scan N range 5–16 |
| Final class 0 / class 1 / `?` | 189 / 33 / 0 |
| GMM / k-means components | All 222 rows usable; 1 / 4 views each |
| Final confidence, min / median / mean / max | 0.0375 / 0.60 / 0.57016 / 0.7625 |
| QC review / no listed flag | 19 / 6 chains |
| Figures | 25 standalone maps and one grid |

These class counts are **model outputs, not the sample's chemical composition**.
`probability_1` is the arithmetic mean of the two component vote fractions;
`confidence = 2 * abs(probability_1 - 0.5)` is an uncalibrated vote margin, not a
probability of correctness. The fixed rule assigns class 1 when
`probability_1 >= 0.5` (including a tie; none occurred here). It emits `?` for an unavailable
component, not every low finite margin. Five margins are below 0.20 and the two
component binary assignments disagree at five keys; these are descriptive
counts, not new rejection rules. No abstention threshold was added.

All 19 QC review rows have mean confidence <0.60. `251206_038.sxm` (N=16) also
has `n_outlier`, relative to this cohort's IQR fences [3.5, 15.5], not an external
length label. Six `ok` rows mean only that the existing QC found no listed flag.
The weakest chain means are `260116_017.sxm` (0.43571) and `260221_044.sxm`
(0.4375). No missing-plot or noncontiguous-key flag occurs.

Read-only log review found no exception, Julia warning, missing-view message or
Fisher iteration-limit warning. The sole split-fit diagnostic is a coverage-floor
clamp from 6 to the already fixed N=5 for `260222_043.sxm`; all five selected rows
are retained. Absence of a log warning does not prove every optimizer start or
classifier seed converged; those diagnostics are not fully serialized.

Other saved-output diagnostics remain relevant to interpretation:

- Split skew spans 0.500031–1.999874, within the unchanged [0.5, 2] bounds. Sixteen
  rows across 14 scans are within 0.001 of a bound (five lower, eleven upper).
  At six-decimal export precision, parallel widths reach the effective optimizer
  endpoints in 1 lower / 17 upper rows, perpendicular widths in 1 lower / 78 upper
  rows. The nominal width interval is [0.191, 0.509] nm; the optimizer's [-5, 5]
  transform reaches about [0.193128, 0.506872] nm. These are constrained fits,
  not out-of-range values or a new validity test. No bound was retuned.
- The reconstructed descriptor is finite with reason `ok` for every row, but
  is concentrated toward +1 (median 0.964083; 128/222 absolute values >0.95).
  This does not by itself establish failure, identity with the lost descriptor,
  or chemical meaning. All Fisher scores and the three residual patch grids
  are finite; the patch rows have nonzero mass.
- Mold-only forward/backward binary decodings differ at **120/222 keys**. These
  are intermediate diagnostics: the predictor uses the nonnegative absolute
  `cost_margin` columns, not those decoded sequences. They therefore provide
  **no forward/backward chemical confirmation** of the final assignments.

The full grid and seven priority standalone maps were inspected: both rescued
scans, the two lowest-mean-confidence scans, the N=16 outlier, the N=5 clamp case,
and `260221_048.sxm` (minimum lobe margin). Markers and class labels are readable
in the standalone plots. These are **base fitted-coordinate maps, not raw STM
overlays**; colors do not encode margin magnitude. They do not establish correct
counts, separately resolved chemical units, or chemical identity.

Full146 reconstruction, external comparison and grading were subsequently
approved (see below), but have not run. The frozen champion's accepted metrics
are unchanged. This completes the unknown25 technical application and integrity
checks, not chemical validation or milestone 1 as a whole.

### Approved full146 comparison scope (2026-09-17)

After reviewing the unknown25 result, the user approved full146 reconstruction,
then external comparison and grading. The scope is one fixed configuration and
one cluster job, bounded at 24 hours; it is not a descriptor/threshold search.
`hpc/reconstructed_full146.sbatch` runs only this application. It reuses the
unchanged base geometry and native templates, computes genuine split fits at
fixed per-file N in at most four one-thread children, then runs all assignment
stages on the full 146-file cohort. No unknown25 rerun, base/count rerun or
benchmark label is part of the job. Comparison and grading happen locally only
after complete outputs are fetched and checked.

**The two count artifacts must not be confused.**

- The approved cache is
  `results/hierarchical_feature_extraction/9e3bcdeb83469fe0/features_merged.tsv`:
  **900 contiguous keys / 146 files**, all elliptical base geometry. Its SHA-256
  is `2cb1065cf5aa790dac1a5fec279d955c732cc91f7e60312103de32fddfb89484`.
  Cached N ranges from 4 to 10 and is not fixed at six.
- This cache comes from label-free GCV selection in the direct feature
  extractor, not the promoted batch `support_midpoint_hybrid` selector. Its
  archived config SHA-256 is byte-identical to current `config/chitosan.toml`.
  The producer reads dataset/file membership, not truth, expected N or reference
  predictions. No original selected-summary accompanies this nonadaptive cache.
- `results/experiments/6mer_full146/pm2_confirm/summary_overlap060_hard.tsv`
  covers the same 146 scans but has **871 selected lobes**, with **28 per-file
  count conflicts**. Native coverage checks reject pairing it with the cache.
  It stays separate; no count is lowered or geometry row dropped to force a
  match. The existing counting benchmark is not rerun or replaced here.
- Use the existing **features-only** native path. The runner writes
  `selected_from_features.tsv` from the 900 cached keys for fixed-N split fits.
  This is honestly derived metadata, not an alleged original counting summary.
  The source config has no adaptive rescue; no `refined_policy` is discarded.
  Do not generalize this to adaptive unknown25 inputs, where the original
  summary remains mandatory.
- All 900 cached `skew_ratio` values are one: this is not a split cache. An older
  39-file split table also has incompatible N on four scans. Neither is passed
  as `--split-features`; job 11786116 computed all 146 genuine split fits anew.

Prepared inputs are in
`results/reconstructed_cc_soft_v1/full146_v1_inputs/`: byte-identical
`base_geometry_full146.tsv`, `templates_cc.tsv` and 146 valid `full146_raw/`
symlinks. The raw names match the cache exactly. The wrapper's optional
`STMFIT_SELECTED_SUMMARY` is used only when an explicitly supplied original
summary matches the cache; it is unset for this run. It does not auto-discover
a summary. `STMFIT_OUTDIR` names the direct new full146 output directory.

**Completed September 18:** Viper job **11786116**, submitted once from `00528d2`,
finished `0:0` in **11m04s** (00:36:57–00:48:01 CEST). All 146/900 outputs have
been fetched and checked; external comparison and grading are complete below.
The eight-CPU allocation used Julia 1.13 without changing the Manifest or
scientific settings. Code and inputs were staged separately from the older Viper
tree. See `hpc/README.md` and the journal for paths, setup failures and recovery.
The monitor is stopped. **Do not resubmit.**

Completion requires all 146 files / 900 keys at each required stage, preserved
base geometry and settings, explicit invalid-component reasons where needed,
QC/maps, and a separate keyed comparison plus full145 own-N external grade.
Partial outputs cannot substitute for full coverage. The grade must report
missing control positions and extra predicted lobes as well as classified
accuracy and exact chains. Equality of scores or assignments would not prove
recovery of the lost descriptor. This fixed-cache comparison is **not a fresh
Julia 1.13 count benchmark or a rebuilt promoted counting-to-assignment chain**.

### Completed full146 reconstruction and grade (2026-09-18)

The fixed application is under
`results/reconstructed_cc_soft_v1/viper_20260917_full146_v1/`. All **184 output
files plus Slurm log** match remote hashes; the 358 remote job source/key-runtime
identities and 148 prepared inputs stayed unchanged during execution. Local Julia 1.13 saved-output verification
passes **73/73** checks, including full 146/900 key coverage in twelve tables,
genuine split chunks, base/config preservation, native component/final validators,
and byte-identical saved descriptor/join/vote/summary/QC replays. No fit or
classifier was rerun locally.

**Assignments and uncertainty.** The final table has **695 zero / 198 one /
7 `?`**, with 893 `ok` and seven `unavailable_gmm` reasons. Seven partial patches
remain in every table and invalidate the strict descriptor/Fisher input checks.
K-means retains two of its four views on those rows; GMM has no usable joint
view, so the final vote does not substitute the lone available component.
The abstaining keys are `240310_Cu100009:8`, `240818_015:1`,
`240818_019:1–4`, and `241113_088:6` (all filenames end `.sxm`). For
`240818_019`, two mold rows have infinite class costs, and the auxiliary scorer
prints `0000` with infinite total cost. That is not chemical evidence; the final
assignment correctly remains `????`. Only local cost margins, not decoded
sequences, enter the predictor.

Mean/median confidence is **0.79082 / 0.8** over all rows. There are **58 assigned
margins below 0.20**, including **twelve exact finite ties** that remain label 1
with confidence zero by the fixed `p >= 0.5` rule. These are not `?`. The two
binary components disagree on **146/893 jointly available lobes**, across 106
chains. Confidence remains an uncalibrated vote margin.

QC flags **43/146 chains**: 40 N outliers, six low mean margins and one high
uncertain fraction, with overlapping reasons. The N quartiles are both six, so
the existing cohort-relative IQR rule flags every non-six count; this is not an
expected-N prior or permission to change N. `241113_088` is QC `ok` despite one
`?`. Missing-view flags are not evaluated by this final-table auto-QC call.
Clean logs do not prove convergence of every optimizer/classifier seed.

**Diagnostics, not tuning criteria.** There are 94 split-skew rows at the rounded
effective optimizer bounds (45 lower / 49 upper), and 297 perpendicular widths
at their rounded upper bound. The reconstructed descriptor median is 0.958684;
634/893 finite values are at least 0.9. Fisher scores are positive on 875/893
valid rows, before the predictor negates the score. These distributions do not
establish chemical separation, calibrated probability or an original-descriptor
match. No bounds, descriptor or thresholds were retuned.

**Figures.** The original 146-panel overview had compressed/clipped panels and
overlapping labels. Plot-only fix `6fa3297` produces **19 indexed pages** in
`results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_figures_v1/`.
`summary_grid.png` is explicitly page 1; `summary_grid_002.png` through
`summary_grid_019.png` and `summary_pages.tsv` cover every file once. All 146
standalone PNGs are byte-identical to production. The layout regression passes
**267/267**, all 900 colored markers are present in the overview panels, and
representative pages plus sixteen standalone maps were visually inspected.
Original figures are preserved. These are fitted-coordinate maps, not STM
image overlays; class color does not encode confidence. A separate two-pixel
GR raster-rounding issue was reproduced before correcting the generated
integrity checker's exact-height assumption; it did not change production data.

**External evaluation.** Reports are under
`results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_external/`.
Comparison with `best_labelfree_cc_soft_20260802.tsv` gives **846/892 matching
labels (94.8%)**, zero missing reference keys and eight extra keys, all from
`240310_Cu100009.sxm`, outside the confirmed145 manifest. The filtered 145-file
table has exactly the same 892 keys as the reference. Differences are 39
reference-1 → reconstructed-0, one 0 → 1, and six binary → `?`; reported
confidences are not identical either.

| External full145 own-N quantity | Reconstructed result |
|---|---:|
| Control positions / predicted lobes | 870 / 892 |
| Classified control positions | 849/870 (97.6%) |
| Correct among classified | **666/849 (78.4%)** |
| Correct over all control positions | **666/870 (76.6%)** |
| Emitted errors / unclassified control positions | 183 / 21 |
| Missing control positions | 16, in 14 short-N files |
| Extra predicted lobes | 38, in 25 long-N files |
| Exact chains | **24/145** |

The 21 unclassified control positions are 16 missing plus five unavailable
under the existing grading alignment. The report's post-hoc `honest_uncertain`
value, **204/870**, also includes the 183 emitted errors; it is **not** a count
of production abstentions. Extra lobes are reported, not aligned to the six
external control positions or silently dropped from the production table.

The frozen reference remains **677/854 (79.3%) and 36 exact chains**. The new
result is lower and is **not promoted**. The authorized fixed reconstruction,
comparison and grade are complete, with no label feedback into production.
This does not recover the lost descriptor, rebuild the promoted counting chain,
or validate the chemistry of the unknown 10–20mer application.

### Saved native component endpoints (2026-09-19; no improvement)

Two fusion-removal comparisons were fixed before their external grade. They
copy the exact saved k-means or GMM decisions on the **same joint-validity
mask** as the reconstructed soft vote. All 900 cached keys and seven unavailable
keys remain. Neither head is retrained. Counts, features, seeds, class polarity,
thresholds and production defaults do not change. This is not the promoted
871-lobe counting lineage or an end-to-end counting-policy evaluation.

| Profile | Correct / all 870 controls | Correct / classified | Classified coverage | Emitted errors | Exact chains / 145 |
|---|---|---|---|---|---|
| Historical reference | 677/870 (77.8%) | 677/854 (79.3%) | 98.2% | 177 | 36 |
| Reconstructed soft vote | 666/870 (76.6%) | 666/849 (78.4%) | 97.6% | 183 | 24 |
| K-means endpoint, common mask | 662/870 (76.1%) | 662/849 (78.0%) | 97.6% | 187 | 16 |
| GMM endpoint, common mask | 622/870 (71.5%) | 622/849 (73.3%) | 97.6% | 227 | 19 |

Every profile has 16 missing control positions and 38 extra predicted lobes.
Native profiles have five unavailable aligned positions; the reference has
none. The sixth unavailable benchmark row is an extra; the seventh is outside
the 145-file manifest. Missing rows and emitted errors are not abstentions.
Neither endpoint improves the reference or reconstructed control. **No promotion,
new weight/threshold sweep, or restored historical labels follows this result.**
This repeatedly used benchmark supplies development evidence, not untouched
validation. Independent molecule/date/sequence confirmation remains necessary.

The uncalibrated vote fractions differ from calibrated chemical probabilities.
Head confidence is `max(p,1-p)`; fused confidence is `2abs(p-0.5)`. Their magnitudes
and QC cutoffs are not comparable. The exporter preserves the stored head label
even if its printed vote rounds to 0.5, and does not export amplitude for a
second physical remapping by the grader. It uses only saved predictions:

```bash
julia --startup-file=no --threads=1 --project=. test/export_saved_assignment_components.jl \
  results/reconstructed_cc_soft_v1/viper_20260917_full146_v1/predictions.tsv \
  results/reconstructed_cc_soft_v1/viper_20260917_full146_v1/pred_kmeans.tsv \
  results/reconstructed_cc_soft_v1/viper_20260917_full146_v1/pred_gmm.tsv \
  results/new_saved_component_endpoints
```

The destination must be new. Run manifest filtering only through the external
comparison script below, then pass the filtered endpoint TSVs to the existing
`--full145-own-n` report. The recorded four-profile result, paired failures and
read-only error audit are in `results/benchmark_recognition_20260919/`; the
French synthesis is `report.md`. Exporter tests are
`test/test_saved_assignment_components.jl` (Julia 1.13).

The source review separately verified asymmetric filtering in the fused fit
input: `(F + S(B))/2`, versus the symmetrically smoothed ROI map. This is an
old live implementation inconsistency, not a demonstrated cause of the grade
loss. A six-assertion synthetic native probe confirms direction-swap sensitivity;
no raw data or fit ran. Correction needs a separate, fixed comparison with
regenerated control and treatment geometry/features. Neither a correction nor a
third candidate is included in this saved-output experiment.

### Symmetric fused-fit comparison (2026-09-20; pending)

The user authorized one targeted correction and a matched benchmark comparison.
On `research/symmetric-fusion-20260920`, the fit uses `(F+B)/2` for the flattened,
unsmoothed directions; ROI detection still uses `(S(F)+S(B))/2`. This is the
existing single-view/alternate-channel filtering convention. The unchanged
control is `f2191f7`.

Both arms regenerate base and split geometry, patches and native assignments
on all 146 scans at the same 900 cached label-free counts. The original
`selected_from_features.tsv` supplies counts only; no geometry/patch cache is
passed to the paired runner. Physical and classifier configs, fit budgets,
seeds and CC templates are identical. `hpc/compare_fusion.sbatch` runs the two
arms sequentially in one four-CPU, 16 GB, two-hour Viper job, after a dry-run.
The external full145 own-N grade follows complete prediction outputs.

This is an assignment comparison at fixed counts, not revalidation of the
promoted 871-lobe counting result. The historical reference remains archival;
no recognition gain or method promotion is established at preparation time.

### Native numerical conventions

The constant-current port retains the old first-axis-fast cube index, sampling
grid increments, highest occupied z, and first-isovalue-below-target rule. Its
actual GlcN/GlcNAc template TSV is byte-identical to the Python reference on the
fixed inputs/settings. Target height is explicitly 0.50 nm, not the 0.55 nm
field found in the frame files. The Python implementation is used only in tests.

The Fisher port retains PCA10, full-covariance/free-weight GMM2, amplitude sign
mapping, sample latent covariance plus ridge, opposite even/odd-lobe folds,
raw-patch/centered-midpoint scoring, reflection maximum and six-decimal margins.
The old row-major last-axis reversal flips physical **t** for these serialized
patches, despite its historical `flip_u_disk` name. This is retained explicitly.
The native initialization/RNG, Lloyd stopping and responsibility arithmetic
differ from sklearn. Byte-identical sklearn predictions are not claimed.
Invalid or degenerate folds retain their rows with NA and an explicit reason.

### Comparison stays outside production

After prediction, compare with the frozen table without changing parameters:

```bash
julia --project=. test/compare_reconstructed_champion.jl \
  --predictions results/my_reconstructed_run/predictions.tsv \
  --reference results/unit_assignment/best_labelfree_cc_soft_20260802.tsv \
  --outdir results/my_reconstructed_comparison
```

This emits keyed label/confidence differences plus missing/extra lobes. It does
not fit, choose settings, align to a known sequence, or grade automatically.
Optional `--benchmark-manifest` is only an external reporting filter to create
`predictions_for_external_grade.tsv`; run the existing grader separately.
Compare like input cohorts and report refit/runtime/initialization differences.
The frozen champion's accepted 79.3%, 36 exact chains and 677/854 remain a
historical reference, not an accuracy claim for the reconstruction.

Focused tests (Python/NumPy is needed only for numerical reference tests):

```bash
julia --project=. test/test_reconstructed_unit_assignment.jl
julia --project=. test/test_selected_support_context.jl
julia --project=. test/test_cc_mold_native.jl
julia --project=. test/test_empirical_fisher_native.jl
GKSwstype=100 julia --project=. test/test_reconstructed_pipeline.jl --e2e
julia --project=. test/test_compare_reconstructed_champion.jl
julia --project=. test/test_reconstructed_unit_assignment_plots.jl
```

The synthetic pipeline test supplies extracted inputs. A separate real single-
file check is needed for SXM/count/refit wiring. Multi-file Gaussian refits belong
on the cluster after that timing check and a dry run, not on the local machine.
