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

- The **ground-truth sequence** (ordered GlcNAc/GlcN per chain) is available
  for the benchmark (240817 dataset) but must **never enter the fitting or
  selection path**. It is used only by grading scripts (`grade_unit_assignment.jl`)
  and by the `--with-truth` cross-evaluation mode of `analyze_unit_separability.jl`.
- Using the truth to choose the 0↔1 flip (the "oracle" convention) is
  **supervised** and must be clearly labeled as such. The **physical convention**
  (GlcNAc = highest-amplitude cluster, based on the larger acetyl group) is
  label-free.
- The assignment rule must not assume the number of GlcNAc/GlcN units in a
  chain. Rules such as "top-2 lobes are GlcNAc" are composition priors and are
  diagnostic only, not valid label-free assignment.

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
    --emin-ev EMIN \
    --emax-ev EMAX

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
    --emin-ev EMIN \
    --emax-ev EMAX

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
```

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

- **Phase 0**: `grade_unit_assignment.jl` implemented. Awaiting predictions.
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
- **Current best label-free output**:
  `results/unit_assignment/assigned_sequences_selectedN_primary_gaussian.tsv`
  (Gaussian features only, primary files only, `N_selected` from the batch).
- **Phase 2–5**: planned. Decision point after DFT-STM molds or a stronger
  label-free observable is available.
- **Ground truth**: skeleton in `benchmarks/chitosan_240817_unit_sequences.tsv`
  (sequences to be filled in by the user).
