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
- Fully reproducible: `test/build_cc_soft_champion.py` + the mold builders
  rebuild templates, margins, both predictors, the soft vote, and the grade.
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
