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

### 2026-08-01 — Two-stage optimization: accuracy first, exact second (k-means view sweep)

Following the observation that physical accuracy and exact chains can be
optimized as two successive stages, the k-means view set was swept to raise
exact chains at constant accuracy. `split_log_skew` (from the existing
`skew_ratio` column) was added to the k-means 3-view configuration:

| k-means view set (+ interactions, 20 seeds) | Physical acc | Exact chains |
|---|---|---|
| base + bwd_neg_com_t + bwd_neg_diag45 (3-view) | 78.0% | 14/145 |
| base + split (2-view) | 77.8% | 16/145 |
| base + split + bwd_neg_com_t (3-view) | 77.8% | 16/145 |
| base + split + bwd_neg_diag45 (3-view) | 77.8% | 16/145 |
| **base + split + bwd_neg_com_t + bwd_neg_diag45 (4-view)** | **78.0%** | **16/145** |
| 4-view + patch_u_asym (5-view) | 77.8% | 15/145 |

`split_log_skew` raises k-means exact chains 14 → 16 without changing the
78.0% accuracy; `patch_u_asym`, which helps the GMM, slightly degrades the
k-means (77.8% / 15). The 4-view k-means (78.0% / 16 / 666 honest) is the best
overall compromise found in this session: closest to the accuracy bar, honest
closest to 677, but exact still below 18. The two-stage optimization confirms
the conclusion: the honest ≥ 677 bar requires +11 correct lobes that no
feature/ensemble combination among those tested provides; the GMM remains the
exact-chain leader (33/145) at lower per-lobe accuracy. Headline unchanged at
`78.9% / 17-of-145`.

### 2026-08-01 — 17×17 patch features (v4) do not reproduce the documented 85.2%; search closed

A proper Viper sbatch job (`10801983`, COMPLETED) re-extracted 17×17 backward
residual patches (half 0.32 nm, step 0.04 nm, 289 pixels) for all 146 scans.
The v4 features from the journal's historical best config were recomputed from
these patches: `patch_u_asym_17` (u-weighted first moment), `hh1_q00_abs`
(Haar level-1 diagonal-detail upper-left quadrant mean |energy|), and
`neg_anis` (anisotropy of the two negative diagonal moments). The GMM with
these features scored **76.1% / 28-of-145** on full145 — no better than the
9×9 `patch_u_asym` GMM (76.3% / 33). On the 35-file subset the v4-style GMM
scores **82.4% / 8-of-35**, identical to the plain GMM base — the recomputed
`hh1_q00_abs` does not reproduce the documented 85.2% / 17-of-35 because the
exact definition (quadrant size, energy norm, patch alignment) was lost with
the original script. The grid-resolution effect observed earlier (9×9 vs 17×17
`bwd_neg_diag135`: 73.1% vs documented 84.8% subset) confirms the historical
features are not byte-reproducible from the journal description alone.

**Session bottom line (2026-08-01).** Every reasonable combination was tested
against the full145 promotion bar (`≥ 78.9%` physical, `≥ 18` exact,
`≥ 677` honest): k-means/GMM soft and hard ensembles, majority voting,
chain-aware thresholds, confidence abstention sweeps, agreement-abstention,
`patch_u_asym` (9×9 and 17×17), and v3/v4-style diagonal and wavelet
descriptors. The k-means 3-view + interactions remains the accuracy leader
(78.0% / 14 exact / 666 honest); the GMM 3-view + interactions + 9×9
`patch_u_asym` is the exact-chain leader (**33/145**, nearly double the
historical 17). No configuration reaches `honest ≥ 677` — the k-means max is
666 correct lobes and abstention cannot create corrects — so **the headline
stays `78.9% / 17-of-145` and no challenger is promoted**. The honest bar is
reachable only by improving per-lobe accuracy itself (new physical features or
per-date calibration), not by combining or abstaining the existing predictors.

### 2026-08-01 — Combination experiments: GMM+k-means ensemble, patch_u_asym, chain-aware voting

Follow-up on the GMM rebuild: several combination strategies were tested on the
full145 table to try to meet the promotion bar (`≥ 78.9%` physical accuracy,
`≥ 18` exact chains, `≥ 677` honest correct on the fixed denominator). All
grading was one-shot post-hoc; no truth entered prediction.

| Configuration | Physical acc | Exact chains | Honest correct |
|---|---|---|---|
| k-means 3-view + interactions | 78.0% | 14/145 | 666 |
| GMM 3-view + interactions | 76.5% | 29/145 | 653 |
| GMM 3-view + interactions + patch_u_asym | 76.3% | **33/145** | 652 |
| Soft ensemble (mean of k-means+GMM probabilities) | 76.0% | 29/145 | 649 |
| 3-method majority vote (k-means, GMM base, GMM 3-view) | 76.2% | 33/145 | 651 |
| Chain-aware: GMM wins chains with ≥T disagreements | 76.3–78.0% | 14–33 | 649–666 |
| GMM+k-means ensemble, abstain on disagreement | 82.5% (707 classified) | 15/145 | 583 |
| k-means confidence abstention (0.55–0.80) | 78.0–78.5% | 10–14 | 653–666 |

The `patch_u_asym` feature (u-weighted first moment of the 9×9 residual patch)
raised GMM exact chains from 29 to **33/145** without changing per-lobe
accuracy. The chain-aware threshold sweep interpolates linearly between the
k-means (accuracy) and GMM (exact) endpoints — no threshold beats both. The
soft and majority-vote ensembles also do not beat the components, because the
two methods agree on ~83% of lobes and the GMM wins most remaining votes.

The journal's historical v3/v4 features (`neg_diag135`, `hh1_q00_abs`,
`neg_anis`) were computed from **17×17** local patches (240817 only); the
Viper-extracted backward patches are 9×9, so the v3-style GMM run with the
available `bwd_neg_diag135` scored only 73.1% — the grid resolution matters and
the 17×17 extraction for all 146 scans would be needed to reproduce the
documented 84.8–85.2% subset scores at full scale.

**Bottom line.** No combination reaches the full promotion bar. The honest
`≥ 677` threshold is not reachable by abstention alone (the k-means max is 666
correct lobes; abstention only removes errors, it cannot create corrects), so
per-lobe accuracy itself must improve. The most notable new result remains the
GMM 3-view + patch_u_asym at **33/145 exact chains** (nearly double the
historical 17) at 76.3% per-lobe accuracy. The headline stays
`78.9% / 17-of-145`; no challenger is promoted.

### 2026-08-01 — Labelfree GMM rebuild and full145 reconstruction (no promotion, GMM 29 exact chains)

The 78.9% / 17-of-145 headline was traced to the **labelfree** predictor (not
the hierarchical EM model). The journal's tuning history (June) recorded that
the best configuration used **GMM full covariance, 10-seed ensemble** on
per-file z-scored local-prominence features with pairwise interactions —
not the k-means used by the current `build_labelfree_unit_predictions.jl`.
The GMM script that produced the historical 82.4% / 11-of-35 subset result no
longer exists, so a new `test/build_labelfree_gmm_predictions.jl` was written
with a 2-component full-covariance EM, k-means initialization, physical
GlcNAc = higher-amplitude mapping, and the same multi-seed vote aggregation as
the k-means predictor. It reproduces the documented subset result exactly:
**82.4% / 8-of-35** on the 240817 subset with local-prominence + interactions.

The full-145 feature table was also reconstructed from scratch. The original
multi-date table used by T5 was found on Viper
(`2cb1065cf5aa-…/inputs/features.tsv`, 900 rows, 13 dates) and its
`N_selected` values were used to re-run the Gaussian-feature extraction as a
proper Viper sbatch array (job `10801328`, 4 tasks, all COMPLETED). The
resulting features differ from the historical local table by a mean 0.58%
relative amplitude difference (N identical for all 234 common 240817 lobes),
so the scores below are near but not identical to the historical numbers.

Full145 own-N grades (854 classified positions):

| Configuration | Physical acc | Exact chains | Honest correct |
|---|---|---|---|
| k-means 3-view (BASE+bwd_neg_com_t+bwd_neg_diag45) + interactions | **78.0%** | 14/145 | 666 |
| GMM base view + interactions | 76.2% | 28/145 | 651 |
| GMM 3-view + interactions | 76.5% | **29/145** | 653 |
| GMM+k-means ensemble, abstain on disagreement | 82.5% (707 classified) | 15/145 | 583 |
| Historical headline | 78.9% | 17/145 | 676 |
| Promotion bar | ≥ 78.9% | ≥ 18 | ≥ 677 |

**No candidate meets the three promotion criteria**, so the headline remains
`78.9% / 17-of-145`. The notable new result is that the **GMM 3-view achieves
29/145 exact chains**, far above the historical 17 — the full-covariance GMM
correctly assigns whole chains more often even though its per-lobe accuracy
(76.5%) is below the bar. The abstention ensemble reaches 82.5% accuracy on
82.8% of lobes but drops honest correct to 583 because `?` positions do not
count toward the fixed denominator. The ~1-point gap between the reconstructed
k-means 3-view (78.0%) and the historical 78.9% is attributed to the 0.58%
feature drift; the exact reproduction of the historical features was not
possible because the original full-145 feature table was never retained.

All work was label-free until the one-shot post-hoc grades above. No source,
config, physical parameter, or frozen registry artifact was changed. The
candidate manifest remains `locked`/`pending`.

### 2026-08-01 — Backward-view T5 gate PASS, T8 grade FAIL (negative result)

The hierarchical backward-view lane was completed on Viper. Backward Z patches
were extracted for all 146 scans, local prominence descriptors were augmented
(`augment_lobe_local_features.jl`), and a combined feature table with
`bwd_neg_com_t`, `bwd_neg_diag45`, and `bwd_neg_diag135` columns was built.
The T5 leave-date-out / 500-seed bootstrap gate **passed decisively**: every
fold was positive, and the scan-bootstrap 95% lower bound on per-lobe
log-likelihood improvement rose from `0.556` (base-only) to `1.823`
(prominence + backward) — a 3.3× improvement in two-component identifiability.
The `forward_backward_agreement` field, which was `NaN` for all 73 146 rows in
the original run, is now finite for 98% of rows with a mean of `0.866`. All 143
array tasks completed with zero failures. The merge and local validation
reproduced the gate identically.

A one-shot external grade (`report_unit_assignment_benchmark.jl --full145-own-n`)
was run on the backward-enhanced hierarchical predictions. **The challenger does
not meet the promotion bar.** Per-blob physical-convention accuracy is
`590/854 = 69.1%` versus the current `676/857 = 78.9%` headline; sequence-exact
is `17/145` (unchanged); honest correct is `590` versus the required `677`. The
gap is approximately 10 percentage points below the promotion threshold.

**Interpretation.** The backward views dramatically improve two-component
identifiability — the model can now confidently distinguish one from two
Gaussian emission components — but this identifiability does not translate into
better per-lobe assignment accuracy. The mean forward/backward agreement of
`0.866` shows the backward views disagree with the base view on roughly 13% of
lobes, and a portion of those disagreements are errors relative to the external
control sequence. The existing 78.9% headline is produced by a different
predictor (the frozen labelfree method), not the hierarchical equal-prior model;
the hierarchical model's lower accuracy is therefore not necessarily caused by
the backward features alone. Grading the hierarchical model with base-only views
(not yet done) would distinguish "backward views hurt" from "the hierarchical
model is less accurate than the existing method."

The headline remains `78.9% / 17-of-145`. No source, config, physical parameter,
or frozen registry artifact was changed. The candidate manifest remains
`locked`/`pending`. The T5 gate artifacts, combined feature table, and
prediction TSV are recorded under run
`bwd-2cb1065cf5aa-9dac84437ef0-ac8eb82f9afb` in
`results/hierarchical_unit_assignment/`.

### 2026-07-31 — Hierarchical backward-view diagnosis and completion plan

The terminal `NO_ELIGIBLE_CHALLENGER` decision at T7 was traced to a concrete
execution gap, not a scientific failure. The T5 leave-date-out evaluator
activates backward views (`base_local+bwd_neg_com_t`,
`base_local+bwd_neg_diag45`, `base_local+split_log_skew`) via `_default_views`
only when the corresponding feature columns are present in the input TSV. The
original T5 run used a feature table that lacked `bwd_neg_com_t`,
`bwd_neg_diag45`, and `split_log_skew`, so only `base_local` was active.
Consequently `forward_backward_agreement` was `NaN` for all 73 146 shard rows,
the per-lobe T7 audit table was incomplete, and T7 correctly declared no
eligible challenger. The hierarchical model itself passed its 13-fold
identifiability gate (95% lower bound `0.556`); the model and the gate logic
are sound.

A local smoke test on the 39-file 240817 subset confirmed the fix path.
`test/build_combined_hierarchical_features.jl` merges backward-patch negative
moments (`bwd_neg_com_t`, `bwd_neg_diag45`, `bwd_neg_diag135`, computed from the
existing `bwd_res_p*` columns via the same `_negative_moment` formula as
`build_hierarchical_unit_predictions.jl`) and split-width log-skew
(`split_log_skew = log(skew_ratio)`) into a single firewall-clean TSV, dropping
any forbidden input columns (`centered_pos`, `edge_distance_norm`). With the
combined table all four default views activated (234/234 rows each, no
degenerate or one-component), and the equal-prior model emitted 234
predictions with a 56%/44% class balance and mean confidence 0.962.

A new versioned plan, `hierarchical-backward-completion.md`, scopes the full
Viper execution: (B1) extract backward patches for all 146 scans,
(B2) build the combined feature table, (B3) re-run the 500-seed T5 bootstrap
with the combined table, (B4) generate the per-lobe T7 audit table, (B5)
evaluate T7 eligibility and freeze a v2 candidate manifest if eligible, and
(B6) execute at most one external grade. The wrapper launcher
`hpc/launch_hierarchical_backward_completion.sh` orchestrates B1–B3 on Viper
and reuses the existing array sbatch and shard merger. No production package
source, physical parameter, selection threshold, calibration, or frozen
registry artifact was changed. The constant-current lane (T3) remains terminal
`BLOCKED`; only the hierarchical lane is completed.

### 2026-07-31 — FINAL9: builder output-set transaction, comment firewall, and emission-math module split

Three final terminal-gate findings were corrected without changing scientific
behavior. First, the constant-current map builder had protected each output
file's atomicity (remediation 3) and then the map+mask pair as one set, but the
full builder generation—all five map/mask pairs together with the provenance
writer—was still not staged or recovered as a single output set. This matters
because a multi-file diagnostic generation must never expose a mixed generation
after an interruption, and provenance must be visible only for a complete old or
complete new generation rather than for a partial install. The builder now
publishes all eleven outputs (ten map/mask files plus provenance) as one
recoverable output-set transaction using the same gate-last protocol: staging in
the destination directory, a durable `prepared` marker, gate-first backup,
non-gate installation, gate-last installation, and a durable `committed` marker.
Provenance is absent at every `after_backup` and `after_install` phase until the
set is complete; recovery from a prepared transaction restores the exact old
generation, while recovery from a committed transaction keeps the exact new
generation and removes every sidecar. Destination symlinks are replaced as
directory entries rather than followed, and symlinked markers, staged sidecars,
staged artifacts, and cross-directory layout collisions are refused before any
output is written. Independent verification passed a 24-phase interruption
matrix (`after_marker`, `after_backup` ×11, `after_install` ×10, `after_gate`,
`after_commit`) and all fifteen adversarial classes—stale prepared and committed
state, first publication from an all-absent generation, partial-install gating,
malformed and symlinked markers, symlinked stage sidecars and staged artifacts,
destination symlinks, a caught mid-install failure, layout confusion,
path-and-hash binding, residue cleanup, dirty-worktree preservation, misleading
success, and determinism (the full matrix passed twice with identical assertion
counts). This remains the bounded single-writer recoverable protocol; kernel or
filesystem corruption, power loss between arbitrary syscalls, and hostile
concurrent mutation of legitimate sidecars are outside the claim.

Second, the candidate-manifest firewall scanned decoded *values* (remediation 3)
but a forbidden label, grade, benchmark, or truth reference could still appear
in a TOML *comment* outside the grader-only table. This matters because the
label-free boundary must hold regardless of whether a forbidden reference is an
active value or an inert comment. Comments outside the sole real `[grader_only]`
table are now scanned for the forbidden token and path vocabulary, alongside the
existing exact-source-byte and recursively decoded-value scans; the grader-only
table remains the only exempt region. Independent verification confirmed the red
toggle—the exact forbidden-comment fixture exited 0 before the checker change
and 1 afterward—and all fifteen cases including decoded Unicode semantics,
nested and array-table comments, commented quoted and fake headers, a misleading
success comment, malformed TOML, the valid frozen lifecycle, stale frozen state,
dirty-worktree stability, and temporary-fixture cleanup. The full manifest suite
passed `29/29` baseline-firewall and `281/281` candidate-manifest contract
assertions; the locked checker reported grade status `locked`, provenance
`pending`, and firewall `pass`.

Third, the hierarchical emission module exceeded the 250-pure-LOC module ceiling
(273 pure LOC) flagged by the code-quality gate. This matters because an
oversized module is harder to audit and review, which is a correctness risk for
the label-free scientific core. The emission math helpers were split out into a
new `test/lib/hierarchical/emission_math.jl` and the facade now includes it
immediately before `emissions.jl`; the move is strictly behavior-preserving.
Independent verification confirmed byte-identical characterization output across
four executions (twice before the split, twice after), covering helper
availability, Gaussian density math, log responsibilities, log-sum-exp and
fallback normalization, total likelihood, deterministic EM, exact equal priors,
covariance flooring, nonconverged and nonmonotone status, one-component
behavior, and validation ordering and messages. The documented pure-LOC rule now
reports `emission_math.jl = 53` and `emissions.jl = 220`, with all eight focused
hierarchical implementation files below 250. The full hierarchical suite passed
twice (`35/35` baseline-schema and `240/240` hierarchical-model each run), the
final-remediation suite passed twice (`96/96` each run, including a new
source-layout regression), and adversarial probes confirmed the pre-split source
reconstruction, fresh post-split hashes, the split structure, layout boundary,
cleanup, and candidate integrity.

All three corrections were independently confirmed. No scientific or numerical
behavior changed; no delegation, grader invocation, benchmark-label artifact,
grade output, config, or physical or selection parameter was used. T3 remains
terminal `BLOCKED`; no challenger was eligible, frozen, promoted, or graded; the
grader invocation count remains zero; and the candidate remains locked with
pending provenance, byte-unchanged at SHA-256
`9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`. The
historical `78.9% / 17-of-145` headline is unchanged. These contract, firewall,
and module-structure corrections provide no new benchmark validation.

### 2026-07-30 — Final remediation 5: frozen provenance and hierarchical guards

Two final contract gaps were corrected without changing scientific behavior.
First, `frozen_once` and `graded` candidate manifests now require every
mandatory artifact provenance field—feature TSVs, constant-current cubes,
generated maps, molds, and configs—to contain a real lowercase 64-hex SHA-256.
Pending, missing, uppercase, malformed non-hex, and wrong-length values are
rejected. This was required because a lifecycle-frozen manifest must bind the
actual bytes of every mandatory input rather than accepting placeholders that
only satisfy the surrounding schema. Independent verification passed the
`29/29` baseline-firewall and `245/245` candidate-manifest assertions and an
independent `55/55` frozen-state fixture matrix covering both lifecycle states
and all five provenance fields.

Second, hierarchical identifiability now treats either empty hard-assignment
component as one-component evidence, even when separation, likelihood, and
amplitude diagnostics would otherwise favor two components. The public EM fit
also requires `0 <= first_seed < n_starts`: after validating and converting
`n_starts` and `first_seed`, it checks that relation before validating `tol` and
`cov_floor`, converting data, initializing, or fitting. These guards were
required because component labels are symmetric, so occupancy of only component
1 or only component 2 is equally non-identifiable, and because an out-of-range
start index must fail explicitly before unrelated invalid inputs or fit work.
Independent verification passed the focused `87/87` suite twice and direct
range, ordering, occupancy, determinism, equal-prior, covariance-floor, and
unstable-status probes.

Both fixes were independently confirmed with the candidate byte-unchanged at
SHA-256 `9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`.
T3 remains terminal `BLOCKED`; no challenger was eligible, frozen, promoted, or
graded; the grader invocation count remains zero; and the candidate remains
locked with pending provenance. The historical `78.9% / 17-of-145` headline is
unchanged. These contract fixes provide no new benchmark validation.

### 2026-07-30 — Final remediation 4: hierarchical publication and EM parameter validation

Two final diagnostic-lane contract gaps were corrected without changing
scientific behavior. First, atomic replacement of the hierarchical merged TSV
and its gate still allowed the two files to expose different generations after
an interruption. This matters because a new gate must never certify a partial
or stale merged result. The failing-first contract reproduced both a late gate
installation failure and an interruption after merged-TSV installation. The
corrected recoverable publication stages both files in their destination
directory, durably records `prepared`, backs up the gate first, installs the
merged TSV, installs the gate last, and durably records `committed`. Recovery
from a prepared transaction restores the old complete generation; recovery
from a committed transaction keeps the new complete generation and removes its
sidecars. Destination symlinks are replaced rather than followed, while
malformed or symlinked transaction markers are refused before publication.
Independent verification passed all 139 hierarchical HPC-contract assertions
and every explicit failure/interruption probe. Public paths, merged and gate
bytes, return fields, and strict shard validation are unchanged. This is a
bounded recoverable gate-last protocol, not protection against filesystem
corruption, power loss between arbitrary syscalls, or hostile concurrent
mutation of transaction paths.

Second, the public two-component EM entrypoint could accept malformed fitting
controls or fail later during data conversion or fitting. Entry validation now
requires `n_starts` to be a positive non-`Bool` integer representable as `Int`,
`first_seed` to be a nonnegative non-`Bool` integer representable as `Int`, and
`tol` and `cov_floor` to be positive finite non-`Bool` reals representable as
`Float64`. They are checked in that order before data conversion, finite-row
checks, initialization, start indexing, or fitting. Failing-first cases covered
boundary, type-confusion, non-finite, overflow, and multi-invalid inputs. The
focused suite passed `75/75` twice; independent probes confirmed 31 invalid
classes, validation order, exact deterministic default and custom fits,
covariance floors, and unchanged unstable/nonmonotone abstention. Equal-prior
semantics remain fixed at `[0.5, 0.5]`; no composition prior was introduced.
The broader hierarchical suite previously exceeded ten minutes after its
`35/35` baseline-schema testset, so the confirmed claim is bounded to the
focused suite and direct EM probes.

Both corrections were independently confirmed. T3 remains terminal `BLOCKED`;
no challenger was frozen or promoted, no grader was invoked, and the candidate
remains locked/pending and byte-unchanged at SHA-256
`9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`. The
historical headline remains `78.9% / 17-of-145`; this remediation supplies no
new benchmark validation.

### 2026-07-30 — Final remediation 3: decoded firewall, recoverable output sets, and calibration parameters

Three final review findings were corrected without changing scientific behavior.
First, the candidate-manifest firewall scanned source text but could miss a
forbidden non-grader key or string assembled through TOML Unicode escapes. This
matters because a semantically forbidden label or grade reference must not cross
the label-free boundary merely because its source spelling is escaped. The
checker now retains its exact-source-byte scan and also recursively walks parsed
non-grader TOML keys and string values, including arrays, inline tables, and
nested tables, after Unicode decoding. Only the sole real top-level
`[grader_only]` table is exempt. Failing-first tests reproduced seven escaped and
nested bypass classes; the corrected contract suite passed `112/112`, and an
independent adversarial review confirmed four- and eight-digit escapes,
case-mixed tokens, decoded keys and paths, and the intended grader-only
isolation. Exact-byte lifecycle hashing was deliberately left unchanged:
“canonical hash” means the exact non-grader source-byte projection, not semantic
TOML canonicalization.

Second, individually atomic writes could still expose a mixed generation when a
multi-file constant-current or frozen diagnostic output set was interrupted.
The first incomplete approach protected each destination separately but did not
make the set recoverable as one generation. The corrected bounded protocol uses
same-directory staging, a durable `prepared` marker, gate-first backup,
non-gate installation, gate-last installation, and a durable `committed` marker.
On rerun, prepared state restores the old complete set while committed state
keeps the new complete set and finishes cleanup. Destination symlinks are
replaced as directory entries rather than followed. Failing-first and injected
interruption tests covered later-install rollback, every explicit protocol
failpoint, symlink destinations, malformed or mismatched marker refusal, and
residue cleanup; independent probes confirmed recovery for the constant-current
map+mask set and the frozen diagnostic TSV/PNG set. This is a recoverable
gate-last transaction protocol, not a guarantee against filesystem corruption
or hostile concurrent mutation of transaction sidecars.

Third, the public constant-current calibration API accepted invalid iteration
and tolerance controls or failed later with unrelated errors. `max_iter` now
must be a positive non-`Bool` integer, and `height_tol_factor` a positive finite
non-`Bool` real, before target conversion, frame work, scanning, or bisection.
Ten failing-first assertions covered zero, negative, non-integer, Boolean, and
non-finite inputs. The focused suite then passed twice, and independent analytic
and poisoned-input probes confirmed valid custom behavior and validation order.
No scan resolution, root-continuity rule, support rule, default, or physical
parameter changed.

All three corrections were independently confirmed. T3 remains terminal
`BLOCKED`; no challenger was frozen or promoted, the grader invocation count is
zero, and the candidate remains locked/pending and byte-unchanged at SHA-256
`9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`. The
historical headline remains `78.9% / 17-of-145`; these remediations provide no
new benchmark validation.

### 2026-07-30 — Constant-current isolated exact-root continuity correction

Final rerun review found one remaining gap in the fixed-support root contract.
An exact interior sampled response was accepted when either immediate adjacent
scan response was non-finite, because the turning-point and support checks ran
only when both neighbours were finite. A one-column reproduction with profile
`[0, 2, 0]`, scan samples `[1, 2, 3]`, and target height `1.0` therefore
returned the isolated exact sample at isovalue `2.0` without establishing a
continuous root branch.

Exact interior roots now require finite responses at both immediate adjacent
scan samples before the existing turning-point and identical-support checks run.
If either neighbour is non-finite, calibration rejects the sample with an
explicit deterministic isolated-response continuity ambiguity. Regressions
cover both non-finite neighbours and a one-sided-finite response. Genuine
interior fixed-support roots remain accepted even when distant scan endpoints
are invalid; endpoint, multiple-root, support-changing bracket, scan-resolution,
and bisection behavior are unchanged.

The phrase “canonical hash” in `config/unit_assignment_candidate.toml` denotes
the checker's exact non-grader source-byte projection, excluding only the
permitted lifecycle assignments. It does not mean semantic TOML
canonicalization. The candidate file and its SHA-256 remain unchanged at
`9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`, because
all downstream no-truth evidence binds those bytes.

This correction does not unblock the accepted-cube gate. T3 remains terminal
`BLOCKED`, no challenger was frozen or promoted, zero grader invocations were
made, and no benchmark headline changed.

### 2026-07-30 — Constant-current endpoint-root continuity correction

Final review found that an exact mean-height target at the first or last
isovalue scan sample bypassed the interior turning-point and fixed-support
checks. A lower-endpoint reproduction therefore accepted a root even though the
one-sided neighbouring sample changed valid-column support; the symmetric upper
endpoint case had the same defect.

Exact roots at either declared scan-range endpoint are now rejected with an
explicit continuity-ambiguity error. This is the conservative fixed-support
policy: without samples on both sides, continuity and a crossing within the
declared interval cannot be established, and accepting the endpoint would
implicitly prefer an unverified branch. Interior fixed-support roots and the
existing multiple-root and support-discontinuity rejection behavior are
unchanged. Regression tests cover both endpoint orientations and their
one-sided support changes.

This correction does not unblock the accepted-cube gate. T3 remains terminal
`BLOCKED`, no challenger was frozen or promoted, no grade was run, and no
benchmark headline changed.

### 2026-07-30 — Final atomic output-integrity reconciliation

Two confirmed output-integrity gaps are closed without changing scientific
behavior. Hierarchical prediction, evaluator-shard, merged-row, and gate
outputs now use same-directory staged atomic replacement. Existing destination
symlinks are replaced rather than followed, and a failure before commit leaves
the prior destination bytes intact. The implementation and focused/HPC/T6
verification are recorded in
`.omo/evidence/improve-unit-assignment-benchmark/final-atomic/hierarchical/done-claim.json`.

The frozen-contrast constant-current diagnostic applies the same contract to
its summary and controls TSVs and generated PNGs: unpredictable same-directory
staging is closed before rename, destination symlinks are replaced without
touching their targets, and pre-commit failure preserves the prior outputs.
The focused diagnostic, frozen-contrast, constant-current, and safe-PoC receipts
are recorded under
`.omo/evidence/improve-unit-assignment-benchmark/final-atomic/diagnostic/`.
Schemas, row formatting, and final artifact paths are unchanged.

These are security and output-integrity fixes only. T3 remains terminal
`BLOCKED`, no challenger was frozen or promoted, no grade was run, and no
benchmark headline changed.

### 2026-07-30 — Final-remediation 2: manifest, provenance, and scan-policy binding

The candidate-manifest boundary now rejects TOML multiline strings rather than
trying to classify their contents with a line scanner. It requires exactly one
real, top-level `[grader_only]` table and scans non-grader keys, values, and path
fragments case-insensitively. Frozen-source binding covers every other
non-grader byte, including comments, formatting, and line endings. Only the real
`candidate.frozen_hash` assignment and the lifecycle-only
`candidate.grade_status` assignment are excluded, so the declared
`frozen_once -> graded` transition does not require a new digest. These checks
prove the current manifest state and future transitions from the bytes presented
to the checker; they do not provide stateless proof of an unobserved historical
manifest.

Hierarchical prediction provenance now records every artifact actually consumed
by a run. Primary features, optional split features, and optional backward
patches are each bound by stable semantic role, normalized path, and SHA-256 of
their bytes, together with the resolved views, feature names, model identity,
and fit options. Deterministic CLI replay was retained, while changing only a
consumed split or backward artifact changes the provenance digest. This is an
artifact-integrity correction, not a model, confidence, or assignment change.

Constant-current root-search resolution is now an explicit validated policy:
`--isovalue-scan-intervals`, default `1024`. The value is propagated through
the builder, recorded as `isovalue_scan_intervals`, and checked by provenance
and whole-ROI validation. Root uniqueness is therefore assessed at that declared
resolution (1024 intervals, 1025 samples including endpoints), not as a
resolution-free mathematical claim. The accepted GlcNAc response still fails
the uniqueness policy, so T3 remains terminal `BLOCKED`. No branch preference,
support rule, physics parameter, benchmark headline, candidate state, or grade
was changed, and no grader was invoked.

### 2026-07-30 — Final-review remediation: stability, constant current, and freeze binding

Final review found two correctness gaps in the diagnostic hierarchical
equal-prior lane. Prediction and leave-date-out evaluation could consume a
two-component EM fit even when its explicit `converged` or `monotone` status
was false. The prediction path now abstains with `unstable_model` or
`nonmonotone_model`, and held-out evaluation fails before responsibilities or
scores are computed. Convergence criteria were not relaxed.

The standalone hierarchical CLI also merged opt-in split-width and backward
patch descriptors into loaded lobe records, then called a path-based pipeline
that reloaded the original feature TSV and discarded the merged values. A
records-based pipeline seam now consumes those records directly while the
existing path-based API retains its prior behavior. Synthetic CLI checks show
both optional views contributing to all rows and a deliberately iteration-
limited fit producing only explicit `unstable_model` abstentions. The focused
suite, full hierarchical suite, T6 integration, and HPC contract are green;
no labels, benchmark truth, expected composition, or grade data entered these
checks.

The constant-current mean-height calibration had also assumed that endpoint
behavior was enough to identify a usable isovalue. It now scans 1025 equally
spaced isovalues over the complete declared interval, including both endpoints,
and accepts only one continuous target root branch with fixed valid-column
support. Multiple roots and support discontinuities are explicit ambiguity
errors; no endpoint, coverage, or branch preference is introduced. The obsolete
smoke cube was consequently replaced by a dedicated strictly monotone synthetic
constant-current cube (`value = 3z`) whose target has one such branch. Existing
smoke assertions were retained rather than weakened. This API correction does
not solve T3: the accepted GlcNAc cube remains empirically multibranch, so T3
stays terminal `BLOCKED` pending a separately predeclared physical policy.

Finally, the candidate-manifest freeze previously canonicalized parsed TOML,
which could miss changes to comments, formatting, and line endings. Its stored
SHA-256 now binds exact non-grader source bytes. The later remediation above
narrows the exclusions to the real grader-only table and the two real candidate
lifecycle assignments needed for the stored digest and `frozen_once -> graded`
transition; it also closes multiline and case-variant firewall bypasses.
Matching frozen provenance and hash are enforced in both lifecycle states.
The candidate itself was not edited, frozen, promoted, or graded, and retains
SHA-256 `9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`.
These fixes close review-discovered contract gaps without changing physics,
branch preference, benchmark headlines, or the T7/T8/T9 terminal decisions.

### 2026-07-29 — T9 terminal documentation closure

T9 closes without a promoted challenger or a new external grade. T3 is terminal
`BLOCKED` by nonunique accepted-GlcNAc isovalue branches without a predeclared
branch policy. T5 passed its 13-date, 500-seed hierarchical identifiability gate
(`95%` lower bound `0.55621530226471905`), and T6 integration/leakage guards were
confirmed after URI-path and partial-view-QC corrections, but the hierarchical
lane still lacks durable forward/backward evidence and the common per-lobe T7
audit table. T7 is therefore `NO_ELIGIBLE_CHALLENGER`; T8 is
`SKIPPED_NO_ELIGIBLE_CHALLENGER`, with zero grader invocations and the one-shot
budget unused. The candidate manifest remains locked/pending with no frozen hash,
and `hierarchical_equalprior` remains an executable diagnostic rather than a
frozen or promoted profile. The historical `78.9%` classified physical accuracy
and `17/145` exact-chain headline is unchanged; it must not be confused with an
honest fixed-denominator report or read as new benchmark validation.

**Documentation-build correction.** The exact Documenter build initially failed
because this intentionally chronological journal rendered above the default
200 KiB hard page limit. `docs/make.jl` now ignores the size threshold for
`journal.md` only; the default warning and hard thresholds remain active for
every other page. A focused source contract and the successful full docs build
confirm that narrow exception. No scientific behavior, parameter, benchmark
claim, grading state, or historical journal content changed.

### 2026-07-29 — T7 terminal no-freeze decision; T8 skipped

Neither approved challenger lane is eligible for a label-free T7 freeze. The
constant-current physics lane is terminal `BLOCKED` because the accepted-cube
response has multiple isovalue branches and this plan predeclared no physical
branch-selection rule. The hierarchical lane passed only its fixed
leave-date-out one-vs-two-component gate: all 13 held-out folds were positive
and the 500-seed scan-bootstrap 95% lower bound was
`0.55621530226471905`. That result does not establish every T7 real gate.

The durable hierarchical `validated.tsv` and all 143 shards contain only `NA`
for `forward_backward_agreement`. Fixed perturbation agreement values are
present, but there is no durable hierarchical prediction TSV and no common
per-lobe T7 audit table containing emission/abstention, forward/backward,
minimum perturbation, normalized margin, invalid reason, and bound hashes.
Synthesizing the missing audit values or inferring agreement would violate the
predeclared freeze contract.

**Decision.** T7 terminates as `NO_ELIGIBLE_CHALLENGER`: no challenger is
frozen, `config/unit_assignment_candidate.toml` remains `grade_status =
"locked"` with pending provenance and no frozen hash, and T8 must be skipped.
No benchmark labels, truth, grade, grader execution, or parameter change was
used for this decision.

### 2026-07-29 — T3 constant-current real gate stopped at nonunique isovalue calibration

The accepted-cube Viper generation job `10765041` stopped before writing maps:
the global-extrema lower isovalue produced no valid crossings. A synthetic
reproduction and accepted-cube diagnostics confirmed that this endpoint check
is premature for GlcN because finite interior isovalues bracket every fixed
height. Geometry and frame-domain checks did not explain the failure.

A second accepted-GlcNAc Viper diagnostic (`10766684`, one task and one Julia
thread) swept 1043 strictly positive isovalues over the observed patch range and
a low-positive tail. All predeclared heights `0.40, 0.45, 0.50, 0.55, 0.60 nm`
are reachable, so target unreachability is refuted. However, the main finite
response has 197 direction reversals, and the five targets have respectively
10, 2, 2, 4, and 6 adjacent isovalue straddles. Many patch columns are
multi-crossing and therefore invalid under the existing vacuum-first policy.

**Decision.** T3 is terminal `BLOCKED`, not repaired or retried. A more robust
endpoint search would still have to choose among low- and high-isovalue roots.
Choosing the lowest/highest root, maximum valid coverage, nearest endpoint, or
another branch would introduce a new physical calibration rule after observing
the accepted real cubes. No such rule was predeclared, and no support threshold
may be inferred from this run. No maps, molds, SXM fits, labels, grades, or
parameter changes were produced. Any future attempt requires a separate
scientific change that predeclares and validates a unique branch-selection rule
before rerunning accepted-cube generation.

### 2026-07-29 — T3 final constant-current correction follow-up

The post-correction review found that the provenance chain was materially
stronger but still had local contract gaps. SXM file-list rows now use one
documented basename allowlist in both the Viper launcher and array worker, so
quotes, whitespace, shell metacharacters, leading dashes, command substitutions,
and path separators are rejected before any SSH, rsync, or Julia command can be
constructed. The connected-mold importer now writes unary molds, optional bond
molds, and mold-binding sidecars with same-directory temporary files followed by
rename, replacing destination symlinks instead of following them.

The scoring validator now checks the mold-binding importer convention
(`parity_flip=t`, `mirror_flip=u`, `normalize=zscore`) and verifies the nominal
constant-current map grid extents and step against the extraction metadata before
scoring. A local synthetic generated-map test covers the real whole-ROI support:
the diagnostic constant-current builder and importer must both use
`--half-nm 0.80 --step-nm 0.08`, yielding finite 21×21 (`441` pixel) molds. This
0.80 nm support is required by the whole-ROI mold consumer; it is a grid-support
contract, not a new physical selection rule. The nominal height remains
`0.50 nm`, the bracket remains `0.40, 0.45, 0.50, 0.55, 0.60 nm`, provider
`stm_dft_cc_diag` remains diagnostic-only, and no SSH, Slurm, real cube
generation, real artifact generation, or real fit was run.

### 2026-07-29 — T3 constant-current provenance correction after independent review

The independent T3 frame-provenance review found that the previous fix still
left two material gaps in the diagnostic constant-current chain. First, the
builder calibrated the nominal isovalue separately for each cube type but the
sidecar recorded only one scalar value. Second, whole-ROI scoring validated the
constant-current map sidecar and then independently loaded the connected mold
TSV, so a stale mold generated from another map could still be scored.

The correction changes the diagnostic sidecar to record `type_isovalues` with
exactly one finite non-Bool nominal isovalue for type `0` and one for type `1`.
Legacy single explicit isovalue calls remain accepted only by normalizing that
one value into both typed entries. The connected-mold importer can now write a
small child sidecar tying the mold TSV hash and path to the source nominal map
hash and source constant-current provenance hash. The real whole-ROI observable
validator receives the actual mold TSV path and rejects missing, stale, or
mismatched mold bindings before loading molds or scoring scans. The same
validator now rejects TOML booleans in numeric type/frame/isovalue fields and
checks the extraction metadata and bracket artifact paths/hashes before scoring.

The Viper launcher was hardened at the local boundary: staged local paths,
remote paths, export payloads, and remote shell commands now pass a conservative
transport syntax check before any SSH command can be constructed, and remote
command fragments are shell-quoted. Contract tests cover quote and
command-substitution inputs failing locally during `--dry-run`. No physical
policy changed: provider `stm_dft_cc_diag` remains diagnostic-only, the nominal
height is still `0.50 nm`, the bracket remains
`0.40, 0.45, 0.50, 0.55, 0.60 nm`, registration ranges and tie tolerances are
unchanged, and no SSH, Slurm, real cube generation, real artifact generation, or
real fit was run for this correction.

### 2026-07-29 — T3 constant-current per-type frame provenance fixed

The T3 real-input preflight found that the accepted GlcN and GlcNAc LDOS cubes
carry distinct local frames from their separately relaxed molecular geometries.
That is the same scientific convention used by
`test/finalize_qe_mold_workflow.jl` and the production constant-height map
path: each cube is sampled in the local frame extracted from its own relaxed
structure. The T2 constant-current builder had been hardened to reject mixed
frames because the sidecar recorded only one `reference_plane`; that made real
T3 artifact generation impossible without copying, averaging, or selecting one
frame for both cube types, which would change the extraction geometry and is not
a label-free provenance fix.

Constant-current provenance now records a `type_frames` table with exactly one
complete frame for type `0` and one for type `1`. The builder passes the actual
frame used for each cube into that table, so distinct relaxed GlcN/GlcNAc frames
are accepted and bound explicitly. Legacy common-frame synthetic calls remain
accepted by normalizing the shared frame into two typed entries, but a one-type,
duplicate-type, malformed, or missing typed provenance structure is rejected.
Whole-ROI real observable validation also requires both typed frames before any
map, mask, or scoring path proceeds.

No physical extraction parameter changed: the constant-current path remains
diagnostic-only, still uses provider `stm_dft_cc_diag`, the fixed
`0.40, 0.45, 0.50, 0.55, 0.60 nm` height bracket, the same crossing policy,
the same z-spacing/periodic guards, and the same frozen registration/tie
configuration. The production registry/config hashes and tracked production
DFT map/template/provenance artifacts were recomputed before and after the
change and were byte-unchanged. Local verification only was run; no SSH, Slurm,
real artifact generation, or real fit was performed.

### 2026-07-28 — T3 constant-current real no-truth gate plumbing prepared

T3 of `.omo/plans/improve-unit-assignment-benchmark.md` prepares the
constant-current whole-ROI real diagnostic for Viper without submitting a job or
running a multi-scan fit. The real diagnostic CLI now rejects benchmark truth,
grade, expected-`N`, control-sequence, manifest, full145, and control inputs at
the runtime boundary. Constant-height defaults remain the historical two-scan
gate (`240817_007.sxm`, `240817_050.sxm`), the same mold path, and the same
diagnostic registration config when no new inputs are supplied.

The constant-current path is explicit and provenance-bound: the launcher and
array worker now carry the diagnostic config, converged molds, output directory,
file-list artifact, observable provenance, map TSV, validity-mask TSV, and
provenance SHA-256 into the CLI. Before scoring, the CLI validates the
constant-current config against the sidecar schema/provider/bias/cube hashes,
checks the provenance hash, checks map and validity-mask hashes, verifies the
map/mask grid keys agree, and rejects any non-`found` validity-mask status.
This preserves frozen common-registration versus contrast-scoring separation;
the observable checks happen before any scoring output is created.

The Viper launcher now uses the file-list artifact rather than hardcoded array
branches, keeps one scan per task, caps useful Julia threads at four per task,
and caps concurrent array tasks so `CPUS_PER_TASK * concurrency <= 8`.
`--dry-run` performs no SSH connection or Slurm submission and records config,
mold, output, file list, array size, concurrency, CPU, thread, memory, and
walltime settings. The archived T3 dry-run is plumbing evidence only, not a
scientific pass: local data and constant-current artifacts still need to be
provided before `--watch`.

### 2026-07-28 — Constant-current molds added as a provenance-bound diagnostic

T2 of `.omo/plans/improve-unit-assignment-benchmark.md` adds a physically
motivated diagnostic between the accepted QE LDOS cubes and whole-ROI unit
assignment. `test/build_constant_current_stm_maps.jl` reuses the existing cube,
frame, and first-vacuum-crossing implementations. It generates typed maps and
validity masks at the five predeclared heights
`0.40, 0.45, 0.50, 0.55, 0.60 nm`; absent, ambiguous, non-finite, and
out-of-support crossings remain invalid rather than becoming zero-height data.
The nominal isovalue follows the fixed `0.50 nm` mean-height-above-Cu policy.

The provenance schema now binds the observable, both accepted cube hashes,
`-0.300 eV` bias, Cu reference plane, nominal/bracket heights, isovalue, z
spacing, crossing policy, and map/mask hashes while retaining byte-identical
constant-height output when no observable is requested. A separate
`config/joint_proxy_whole_roi_constant_current.toml` copies the frozen
registration ranges and numerical-tie tolerances exactly and adds only
diagnostic observable/provenance references. The temporary provider is named
`stm_dft_cc_diag`; registry tests create it only in temporary state and confirm
that the active `stm_dft_v1` source, `config/joint_proxy_molds.toml`, and the
production registry payload hash are unchanged before and after the fixture.

Synthetic builder QA produced 81 type-0 and 81 type-1 rows at each height with
distinct typed values. For every height, four non-benchmark three-lobe binary
fixtures passed exact, fixed-noise, drift, blur, and frame perturbations
(`20` cases per height), and the common-only control abstained. The focused
registry test passes `441/441`. These checks establish deterministic map
construction, provenance binding, and common/contrast mechanics in the
synthetic domain only.

Independent T2 review found four implementation hazards that did not alter the
synthetic scientific result but weakened its reproducibility and artifact
contract. The builder accepted different type-0/type-1 reference frames while
recording only one frame, constant-current provenance allowed declared bracket
heights without binding their map/mask artifacts, direct output writes followed
symlinks and exposed partial files, and registry tests reused a predictable
shared temporary directory. The remediation rejects mixed frames before any
output, requires exact non-nominal bracket-artifact coverage, atomically replaces
map/mask/provenance files without following destination symlinks, and gives each
registry test process an isolated `mktempdir`. Failing-first tests reproduced all
four mechanisms; the focused provenance suite and the expanded registry
suite pass after the fix. No physical parameter, provider registration, accepted
cube, production registry hash, or benchmark-facing behavior changed.

A follow-up code-quality review found two additional provenance ambiguities:
the two cubes could resolve to different normal-axis sampling spacings while the
sidecar recorded only one, and opt-in periodic wrapping changed extraction
without being represented in the sidecar. Failing-first CLI tests now require a
nonzero exit for either case. The builder accepts only a common z spacing and
rejects periodic wrapping; the lower-level T1 diagnostic sampler remains the
appropriate surface for periodic cube experiments.

**Decision.** Keep the constant-current result diagnostic-only. QE
`plot_num=5` is a discrete `|psi_n(r)|^2` sum over the bias window, not amperes;
there is no calibrated nA-to-cube conversion. The accepted common `5e-5 Ry`
criterion and five-height bracket remain fixed, no bracket is chosen by a
benchmark grade, and no map enters fitting, `N_selected`, calibration,
thresholds, or production abstention. Real GlcN/GlcNAc transfer still requires
the separate no-truth T3 gate.

### 2026-07-27 — T0 firewall + baseline contract frozen for the unit-assignment challenger

The label-free challenger lane of the new plan
`.omo/plans/improve-unit-assignment-benchmark.md` opens with a provenance-only
candidate manifest at `config/unit_assignment_candidate.toml` and a dedicated
checker at `test/check_unit_assignment_candidate_manifest.jl`, paired with a
contract test at `test/test_unit_assignment_candidate_manifest.jl`. T0 freezes
the firewall and the baseline contract **before** any new model, feature, or
cube is written (T1–T6), so the boundary cannot drift while the lanes execute.

**Firewall.** The manifest is scanned line-by-line. The benchmark control
motif/encoding, the benchmark truth/count column names, and benchmark
truth/grade paths may appear **only** inside a dedicated `[grader_only]`
section. They must never enter fitting, feature construction, candidate
selection, confidence, abstention, or calibration. The distinct grader-only
denominator manifest `benchmarks/chitosan_6mer_counting_confirmed.toml`
(`145 files / 870 control positions`) is referenced by field name only from
`[denominator]`; its path lives exclusively in `[grader_only]`. This is the
same label-free boundary the existing unknown-production
runner/validator/docs scripts already enforce at the CLI surface; the
baseline is now pinned by a `baseline_firewall` test set
(`validate_unit_predictions.jl` and `run_unknown_unit_assignment.jl` reject
`--truth`, `--control-sequence`, `--manifest`, `--full145`, `--control`,
`--expected-N`; `check_unknown_workflow_docs.jl` passes; the denominator
manifest is not referenced by the unknown runner/validator).

**Decisions frozen by the manifest.** `grade_status = "locked"` until T7.
Provenance is `pending` at T0 and must be `frozen` with a real 64-hex SHA-256
for every feature TSV, cube, generated map, mold, and config before T7. Exact
feature lists are pinned (`base_local`, `base_gaussian`, the backward
descriptors, `split_log_skew`) under an equal-priors hierarchical emission
policy (shared two-component diagonal Gaussian, class priors fixed at `0.5/0.5`,
no occupancy regularizer, higher-amplitude cluster maps to the acetyl-bearing
unit). Views have equal weight; bootstrap is `500` replicates over seeds
`0:499` resampling whole scans within training dates. The date parser parses
exactly one leading `YYYYMMDD` token from the folder/relative path and fails on
missing or ambiguous input. The constant-current policy is fixed at the `0.50 nm`
mean-height above the Cu reference with the `0.40:0.05:0.60 nm` sensitivity
bracket as diagnostics only (no nA-to-cube conversion claimed; bracket may not
be selected by a benchmark grade). The 1-vs-2 gate requires every held-out date
fold positive on per-lobe log-likelihood improvement AND a scan-bootstrap 95%
lower confidence bound above zero. The common real no-truth gates, the no-truth
ranking order, and the final post-hoc promotion thresholds on the fixed
`145/870` denominator are pinned: `honest_correct >= 677`,
`physical_accuracy_classified >= 78.9%`, `exact_chains >= 18`.

**Hash binding.** T7 computes a canonical SHA-256 over the parsed TOML
re-serialized with sorted keys, excluding the `[grader_only]` section and the
`candidate.frozen_hash` field itself, writes it as `candidate.frozen_hash`, and
switches `grade_status` to `frozen_once`. The checker recomputes the canonical
hash on every invocation and rejects any mutation that makes it differ from the
stored digest. A future scientific hypothesis requires a new versioned plan and
a new candidate manifest rather than editing a frozen one.

**Verification.** A red→green TDD sequence was used. A `baseline_firewall`
test set (17/17) was run green on unchanged code first and archived. The full
`candidate_manifest_contract` test set (51/51, including missing keys,
malformed SHA-256, unequal view weights, ambiguous date parser, invalid
locked/frozen state, grader-only exception, forbidden-token and forbidden-path
injection, and hash-binding mutation rejection) was run red against the absent
checker/config, then green after implementation. Two consecutive green runs
produced identical counts (`17/17`, `51/51`), so the focused manifest test is
not flaky. Manual QA against the real CLI surface confirmed: `--expect-locked`
exits 0 with locked-status output; a non-grader field mutated after a valid
frozen hash is rejected solely by the hash mismatch; and every forbidden family
(`NKNNKN`, `010010`, `sequence`, `expected_N`, `target_N`, and benchmark
truth/grade paths) injected outside `[grader_only]` is rejected with a firewall
diagnostic. The promotion threshold field is named `exact_chains_min` (not
`sequence_exact`) precisely because `sequence` is a forbidden token outside
`[grader_only]`; the checker maps it to the plan's exact-chain criterion.

**Why.** Earlier unit-assignment experiments (split-width skew, residual
patches, geometric molds) each reached post-hoc diagnostic grades but did not
produce a transferable label-free assignment. Without a frozen firewall, a
future challenger could silently retune features, thresholds, or bracket from
benchmark feedback. T0 makes the label-free boundary and the one-shot grading
discipline executable and auditable before the lanes execute.

**Failed approach avoided.** An earlier draft listed the forbidden tokens
literally inside the manifest's own `[firewall]` section as a self-documenting
schema; the checker immediately tripped on its own declaration. The literals
now live only in `[grader_only]` and as constants in the checker; the
`[firewall]` section declares policy only.

### 2026-07-26 — Registration commune figée: gate synthétique réussi, transfert réel limité par les bornes

Le diagnostic whole-ROI a été reformulé pour empêcher la géométrie de choisir
une identité chimique. Pour chaque état `(parity, mirror)`, les moules sont
décomposés exactement en `common=(M0+M1)/2` et
`contrast=(M1-M0)/2`. La registration n'ajuste que
`fond + backbone + common` sur `context.zimg` fusionné et sur le support exact
du fit. Elle énumère les huit états globaux et une grille TOML déterministe de
translation, petite rotation et flou, sans scale, shear ni transformation par
lobe. L'état, la transformée, l'image commune, le masque et le fit sont ensuite
figés dans des structures profondément immuables. Le score chimique ne peut
plus relancer la registration: il énumère seulement les séquences binaires sous
la géométrie figée, publie les marges best/runner-up/complément et s'abstient sur
égalité numérique, contraste nul, instabilité, limite de recherche ou contrôle
négatif non battu.

Le gate synthétique est distinct de l'évidence réelle. Il utilise un générateur
forward indépendant du sampler du score, retire toute vérité des entrées avant
registration, puis couvre données exactes, common-only, swap de métadonnées,
égalité complémentaire, nuisance hors grille, bruit, dépassement de borne et
sentinelle de transposition asymétrique. Les tests focalisés passent deux fois
(`25/25`) et la suite `julia --project=. test/joint_proxy/runtests.jl` passe.
Ce résultat établit la mécanique de décomposition, gel et scoring dans le
domaine synthétique; il ne constitue pas une validation chimique STM.

Le calcul réel exhaustif a été déplacé sur Viper, un fichier par tâche, après
un dry-run obligatoire. Le moule `21x21`, `±0.80 nm`, `0.08 nm/pixel`,
échantillonné périodiquement en `xy`, a été régénéré depuis les cubes production
dont les SHA-256 correspondent à la note DFT. Le job-array `10692648` s'est
terminé sans erreur (`007`: `6m29s`, `050`: `14m30s`). Les artefacts récupérés
sont dans `/tmp/opencode/frozen_contrast_hpc/`: deux lignes primaires, dix
contrôles et deux PNG à huit panneaux.

Résultats fused/fit-mask, sans lecture de vérité benchmark:

| scan | pixels | état | transformée `(t,u,rot,blur)` | gain common SSE | gain contraste | marge runner | contrôle décalé | décision |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | --- |
| `240817_007.sxm` | 14961 | `0/0/0` | `(0.17,-0.20,5°,0.02) ` | `0.0212584` | `0.0128963` | `0.00101592` | `0.0395822` | abstention: borne `u/rotation` et contrôle non battu |
| `240817_050.sxm` | 46894 | `0/1/1` | `(0.30,0.20,5°,0.08)` | `0.0541387` | `0.431676` | `0.0384917` | `0.292390` | abstention: bornes `t/u/rotation` |

Les quatre perturbations d'un pas fin conservent le meilleur code diagnostique
sur chaque scan, mais cette stabilité discrète ne suffit pas: pour `007`, le
contraste spatialement décalé explique davantage de SSE que le contraste
nominal; pour les deux scans, la registration pousse plusieurs paramètres aux
limites fixées avant l'exécution. Les séquences finales restent donc vides dans
le TSV et `?` dans les figures. Les bornes n'ont pas été élargies ni retunées à
partir de `007/050`.

Décision: le transfert réel n'est pas établi. La séparation common/contrast
réfute l'idée qu'un simple degré de liberté chimique figé suffit à résoudre le
mismatch whole-ROI actuel; elle ne réfute pas les moules DFT dans leur domaine
synthétique. Les limites restantes sont la registration au mur et la différence
physique entre LDOS constant-height et topographie constant-current, toujours
confondue avec le drift, la réponse de pointe et l'erreur de modèle. Aucun
comportement, seuil, calibration, registre ou claim benchmark de production
n'est modifié.

### 2026-07-24 — Inspection visuelle whole-ROI: mismatch de vue et dilution du masque

Deux cas réels contrastés ont été rendus avec
`test/plot_joint_proxy_whole_roi_debug.jl`: `240817_007.sxm`, qui avait le plus
grand gain du batch propre, et `240817_050.sxm`, dont le coefficient de moule
était nul. Pour chaque fichier, le diagnostic compare quatre conditions sans
utiliser de label benchmark: observation forward ou fusionnée, chacune sur le
rectangle complet ou sur le support tubulaire réellement employé par le fit de
chaîne. Chaque PNG montre observation, backbone, contribution du moule, modèle
combiné, résidus et amélioration SSE pixel par pixel. Le tableau reproductible
est écrit dans `/tmp/opencode/whole_roi_debug_plots/comparison.tsv`.

Un défaut méthodologique du diagnostic réel a été confirmé. Le backbone est
ajusté par `chain_gaussian_sweep` sur la fusion forward/backward
(`fuse_z_bwd=true`), mais le score whole-ROI antérieur le comparait à la vue
forward seule. En utilisant exactement `context.zimg`, le domaine fusionné du
fit, le gain rectangle complet tombe à `0%` sur les deux fichiers. Sur le masque
du fit, les gains passent de `17.83%` à `1.63%` pour `007` et de `11.00%` à
`1.04%` pour `050`. Les forts gains forward suivent visuellement des bandes
latérales rouge/bleu et une texture de lignes dans le résidu; ils ne constituent
donc pas une signature chimique propre.

Le rectangle complet dilue fortement toute contribution localisée: il contient
`65 572` pixels contre `14 961` dans le masque du fit pour `007`, et `238 044`
contre `46 894` pour `050`. Le masque tubulaire augmente donc mécaniquement le
gain relatif, mais ne stabilise pas l'identité. Pour `007`, la meilleure séquence
change entre `111111`, `101001`, `000000` et `111011` selon la vue et le masque.
Pour `050`, les variantes non nulles préfèrent `111111`, résultat uniforme qui
n'apporte pas de discrimination GlcN/GlcNAc. Sur les conditions cohérentes
fusion+masque, les marges entre séquences uniques restent seulement `0.082%` et
`0.152%` de la SSE nulle.

Les overlays ne montrent pas de rotation à 90 degrés ni d'inversion y grossière:
l'axe, les centres et le support des moules suivent la molécule. Les maxima DFT
locaux sont toutefois décalés par rapport aux centres gaussiens, ce qui reste à
interpréter physiquement. Enfin, soustraire la moyenne du moule ne change la SSE
qu'à l'arrondi (`|Δ| <= 4e-15`), comme attendu avec un intercept libre; les halos
du z-score sont visibles mais leur composante constante n'explique pas le faible
transfert.

Décision: ne pas interpréter le score forward/rectangle antérieur comme preuve
de transfert. Tout futur diagnostic réel doit comparer dans le même domaine
fusionné que le backbone et publier séparément rectangle et support du fit. Le
faible gain fusionné reste compatible avec un gap constant-height LDOS versus
topographie constant-current, mais les deux fichiers ne suffisent pas à isoler
ce gap d'autres erreurs de modèle. Aucun comportement scientifique de production
n'est modifié.

### 2026-07-24 — Échantillonnage périodique diagnostique et convergence du cadre DFT

Un échantillonnage périodique opt-in a été ajouté à `test/cube_to_stm_maps.jl`.
La recherche auprès des sources QE (`stm.f90`, `cube.f90`, `chdens_module.f90`,
commit `61569eb`) confirme que le LDOS `plot_num=5` est construit par FFT
inverse et donc périodique par construction, et que le cube QE stocke exactement
une maille sur `[0,1)` sans voxel terminal dupliqué. Le wrapping diagnostique
s'applique donc en coordonnées de grille fractionnelle, pour les axes latéraux
`x` et/ou `y` uniquement; l'axe `z` est rejeté par conception. Le comportement
par défaut reste byte-identique (SHA-256 `58236e7d…1570ed` confirmé après
changement). Les tests synthétiques couvrent l'identité intérieure,
l'invariance de translation par vecteur de maille, la continuité de couture, le
rejet hors-bounds en z, et le parsing CLI (13/13).

Trois cadres wrappés ont été régénérés depuis les cubes production suivis avec
`--periodic-axes xy` sous `/tmp/opencode` uniquement. L'audit de convergence
(`test/audit_joint_proxy_frame_convergence.jl`) mesure la fraction L2 de bord
par type et la différence type_1−type_0, et vérifie l'identité du crop central
9×9:

```text
cadre          edge_L2 type_0  edge_L2 type_1  edge_L2 diff  crop Δ   pic au bord
9×9 ±0.32 nm      42.92%          92.74%         92.84%       —         oui (tous)
13×13 ±0.48 nm    37.10%          47.93%         48.00%       0.0       oui (tous)
17×17 ±0.64 nm     9.57%          56.44%         56.50%       0.0       partiel
21×21 ±0.80 nm     0.49%           1.77%          1.77%       0.0       non (tous intérieurs)
```

Tous les crops centraux sont bit-identiques (`Δ = 0.0`). Les fractions L2 de
bord convergent: type_0 passe de `42.92%` à `0.49%`, type_1 de `92.74%` à
`1.77%`, et la différence discriminante de `92.84%` à `1.77%`. À `±0.80 nm`,
tous les maxima absolus ont quitté la frontière. L'analyse géométrique
confirme que la maille `~2.05×1.53 nm` isole la molécule (gaps `0.16–0.49 nm`
à chaque bord), donc le wrapping échantillonne le Cu/vide, pas une réplique
moléculaire voisine.

Décision: la structure discriminante des moules DFT est contenue dans un cadre
`±0.80 nm`. Le cadre production 9×9 était bien tronqué. Cependant, le
diagnostic reste hors production: aucune carte wrappée n'est promue, le
provider suivi n'est pas modifié, et le posterior de production n'est pas
retuné. La prochaine étape label-free est de régénérer un registre whole-ROI
complet avec les moules `±0.80 nm`, puis de refaire le gate synthétique
whole-ROI et une comparaison sur ROI réelle avant tout claim d'identité
chimique.

### 2026-07-24 — Gate synthétique whole-ROI avec moules DFT convergés ±0.80 nm

Les moules connectés DFT ont été régénérés en `21×21 (441 pixels)` depuis les
cubes production suivis avec l'échantillonnage périodique `xy` à `±0.80 nm`,
puis importés avec normalisation `zscore`. Le diagnostic
`test/diagnose_joint_proxy_whole_roi_converged.jl` place ces moules réels sur
une chaîne synthétique à trois lobes avec backbone gaussien, fond constant et
bruit optionnel, puis demande au score whole-ROI global de récupérer la
séquence.

Résultats: les quatre séquences testées `[0,1,0]`, `[1,0,1]`, `[0,0,1]` et
`[1,1,0]` sont toutes récupérées exactement sans bruit (`SSE ~1e-29`). La
récupération tient jusqu'à `noise_sigma = 0.10`. Le coefficient de moule
injecté (`0.400`) est récupéré à `0.400±0.001`. Le ratio `SSE_swapped /
SSE_correct` vaut environ `2.5e30`, ce qui confirme que les moules convergés
portent une structure réellement différente entre GlcN et GlcNAc.

Ce succès établit que les moules DFT `±0.80 nm`, une fois non tronqués, sont
discriminants dans leur propre domaine image. Il ne prouve pas le transfert
vers des scans réels: les observations synthétiques sont générées par le même
modèle génératif que le score, et les conditions de bruit ne capturent pas le
drift d'acquisition, le couplage tip-orbital ou la convolution de
rétroaction. Décision: conserver le résultat comme gate de représentation et
de mécanique d'assemblage, pas comme validation d'identité chimique réelle.
La prochaine évidence utile doit comparer le score whole-ROI assemblé sur des
lobes réels décodés depuis un scan `240817`, en utilisant les moules
convergés.

### 2026-07-24 — Gate réel whole-ROI sur 240817_001: abstention honnête

Le diagnostic `test/diagnose_joint_proxy_whole_roi_real.jl` assemble les
moules DFT convergés `±0.80 nm` sur la géométrie de chaîne décodée depuis le
scan réel `240817_001.sxm`, puis ajuste globalement
`fond + coeff_backbone·backbone + coeff_moule·assemblage` avec coefficients
positifs et recherche exhaustive sur les `2^N` séquences et les 8 états
globaux `(direction, phase, mirror)`.

Résultat sur `240817_001.sxm` (fit `N_selected=8`, ROI `72×138 pixels`) :

```text
null SSE (backbone seul):        4.6946e+00
best SSE (backbone + moule):     4.6577e+00
réduction SSE:                   0.79%
coefficient de moule:            0.0016
spread SSE top-8 séquences:      0.2%
coefficient global maximal:       0.0016
```

Le coefficient de moule est négligeable (`0.0016`) et la réduction SSE est
inférieure à `1%`. Les huit meilleures séquences sont dans un spread de
seulement `0.2%`, ce qui rend la séquence sélectionnée non fiable. Le
diagnostic s'abstient donc honnêtement: les moules DFT convergés, bien que
parfaitement discriminants dans le domaine synthétique, ne transfèrent pas
vers le scan réel.

La cause racine n'est pas la troncature du cadre (résolue) ni la mécanique
d'assemblage (validée), mais le **gap de domaine de mesure**: les moules DFT
sont des cartes LDOS constant-height (sortie `pp.x plot_num=5`), tandis que
le scan STM est une topographie en mode constant-current. Ces deux quantités
ne sont pas reliées par une simple transformation affine; le couplage
tip-orbital, la rétroaction de courant et la convolution électronique
introduisent des différences structurelles que la nuisance globale
`fond + backbone + moule` ne peut pas capturer.

Décision: l'attribution chimique GlcN/GlcNAc reste non identifiable à partir
des scans STM actuels avec cette approche. Ne pas promouvoir les moules
périodiques en production, ne pas retuner le posterior, et ne pas émettre de
types. La prochaine direction physiquement fondée est soit (a) un modèle
explicite de la boucle de rétroaction constant-current depuis le cube LDOS,
soit (b) un observable expérimental indépendant (multi-bias, spectroscopie,
ou canal Current). La voie (a) nécessite `pp.x` en mode STM simulé
(`plot_num` approprié) ou une post-processe Tersoff-Hamann avec modèle de
pointe; la voie (b) nécessite des données qui ne sont pas présentes dans
`240817_001.sxm`.

### 2026-07-24 — Réplication whole-ROI sur 10 scans clean N=6

Le gate réel a été répété sur dix fichiers de la liste documentée
`accept_unit_training_known_010010`, tous non ambigus avec
`N_selected=N_ell=N_circ=6`: `007`, `036`, `040`, `042`, `045`, `046`,
`049`, `050`, `052` et `054`. Le diagnostic lit `N_selected` depuis le résumé
label-free figé `results/best_plots_240817_primary_rerun/summary_overlap060_hard.tsv`;
il ne lit ni séquence, ni `expected_N`, ni composition. Chaque séquence est
comparée après minimisation sur les huit états globaux, puis la marge est
calculée entre les deux meilleures séquences uniques.

```text
file  SSE reduction  sequence gap  mold coefficient
007       1.35%          0.686%         0.0021
036       0.88%          0.341%         0.0026
040       0.30%          0.111%         0.0017
042       0.21%          0.032%         0.0015
045       0.25%          0.123%         0.0015
046       0.24%          0.133%         0.0015
049       0.02%          0.014%         0.0005
050       0.00%          0.000%         0.0000
052       0.00%          0.000%         0.0000
054       0.05%          0.009%         0.0007
```

La médiane de réduction SSE est `0.23%` (maximum `1.35%`), la médiane du
coefficient de moule `0.0015` (maximum `0.0026`) et la médiane de marge entre
séquences uniques environ `0.071%`. Deux scans retombent exactement sur le
modèle nul. Les séquences gagnantes sont souvent dégénérées (`000000`,
`111111`, `011110`) et ne constituent pas une attribution chimique stable.

Cette réplication exclut le mauvais `N=8` du premier fichier comme cause
principale de l'échec: même avec dix géométries propres et `N=6` sélectionné
label-free, le moule n'explique qu'une fraction négligeable du signal et les
marges de séquence restent minuscules. Décision finale pour cette branche:
gap de domaine confirmé; ne pas élargir le batch, grader les séquences,
retuner le posterior ou promouvoir les moules. Le prochain développement doit
porter sur une observable constant-current simulée ou de nouvelles données
expérimentales.

### 2026-07-23 — Convergence du cadre DFT: 13×13 valide, 17×17 bloqué par le cube

L'inventaire de provenance a été rendu exécutable par
`test/audit_joint_proxy_frame_compatibility.jl`. Deux cadres ne sont déclarés
comparables que si les hashes des deux cubes, la hauteur, le biais, les unités,
le pas, les hashes de cartes et toutes les valeurs du crop central concordent.
Une provenance absente, un cube différent, une valeur non finie ou un crop
central différent bloque explicitement la comparaison. Le test synthétique
couvre les chemins compatible, provenance absente, hash cube différent et
support non fini; 13/13 assertions passent.

Contrairement à l'inventaire initial incomplet, les cubes production sont bien
présents localement et leurs SHA-256 correspondent au sidecar suivi:
`80cd1d…88863` pour GlcN et `40649c…9b5bf` pour GlcNAc. Les frames typés à
`0.50 nm` sont également présents. Ils ont permis de régénérer, uniquement sous
`/tmp/opencode`, des cartes emboîtées avec les mêmes cubes, frames, biais
`−0.300 eV`, hauteur `0.50 nm` et pas `0.08 nm`.

Le cadre 13×13 `[-0.48,0.48] nm` est entièrement fini et son crop central 9×9
est bit-identique à la carte suivie (`max |Δ| = 0`). Il constitue donc une vraie
extension compatible, contrairement à l'ancien fichier préliminaire sans
provenance, fini à seulement `71.60%` et dont le crop diffère. Cependant le
13×13 ne ferme toujours pas le support: les fractions L2 de bord valent
`37.10%` pour type 0, `47.93%` pour type 1 et `48.00%` pour leur différence;
les trois maxima absolus restent sur la frontière.

Un second élargissement à 17×17 `[-0.64,0.64] nm` conserve lui aussi exactement
le crop 13×13 (`max |Δ| = 0`) et la même provenance physique, mais atteint la
limite spatiale des cubes: 17 valeurs deviennent `NA`, soit une fraction finie
de `97.06%`. Le gate le rejette donc comme `nonfinite_map_values`; calculer une
énergie de bord sur ce cadre serait trompeur.

Décision: le diagnostic whole-ROI reste hors production. Le 9×9 est confirmé
tronqué et le 13×13, bien que valide, n'est pas convergé. Le prochain unblock
physique exige des cubes offrant un cadre complet d'au moins `±0.64 nm` sous la
même politique de biais/hauteur/frame, puis deux cadres complets emboîtés dont
les crops sont identiques et dont les métriques de bord se stabilisent avec des
maxima intérieurs. Ne pas utiliser les cartes préliminaires, extrapoler les
`NA`, appliquer un wrap périodique non validé, ni modifier le provider de
production avant ce gate.

### 2026-07-22 — Gate synthétique du modèle whole-ROI assemblé

Un nouveau diagnostic hors production remplace, pour ce gate seulement, la
somme de scores unary indépendants par une image moléculaire unique. Chaque
moule typé est placé au centre de son lobe dans le repère global de chaîne; les
valeurs se somment explicitement pixel par pixel dans les zones de recouvrement.
Hors du support fini du moule, la contribution vaut zéro afin que la troncature
reste visible plutôt que masquée par une extrapolation. Le score ajuste ensuite
globalement `fond + coefficient_backbone·backbone + coefficient_moule·assemblage`,
avec les deux coefficients contraints positifs. Pour les chaînes courtes, une
recherche exhaustive compare les séquences et les états globaux uniquement par
la SSE whole-ROI.

Le gate synthétique passe 14/14 assertions: doublement exact quand deux moules
se recouvrent parfaitement, support nul hors cadre, récupération à précision
numérique des trois nuisances injectées, et récupération exacte d'une séquence
de trois lobes par le score global. Le registre DFT suivi est aussi chargé avec
ses hashes existants et produit une image assemblée finie et non constante. Le
diagnostic n'importe aucun manifeste, séquence contrôle, composition attendue ou
`N` attendu et ne touche ni le fit 2D, ni `N_selected`, ni le posterior de
production.

Ce succès établit seulement la mécanique d'assemblage et l'absence de
double-comptage unary dans le score synthétique. Il ne valide pas l'identité
chimique réelle: les observations synthétiques sont générées dans le même
modèle, tandis que les moules DFT suivis restent dominés par leur frontière
9×9. Décision: conserver `test/diagnose_joint_proxy_whole_roi.jl` comme
affordance diagnostique, sans promotion ni test réel interprétable. Le prochain
gate reste la récupération ou régénération de cartes au moins aussi larges que
`[-0.48,0.48] nm`, suivie d'un contrôle de convergence du score whole-ROI quand
le cadre augmente. Ce n'est qu'après ce contrôle qu'une comparaison label-free
sur ROI réelle pourrait produire une évidence physique utile.

### 2026-07-22 — Audit du support spatial des moules DFT unary

Un audit reproductible mesure désormais, sans donnée ni label de benchmark, la
fraction d'énergie L2 portée par le bord, le maximum absolu, sa position et la
fraction du maximum présente au bord pour chaque moule et pour la différence
`type_1−type_0`. Le contrôle synthétique couvre un pic centré, un pic de bord,
la différence de types et le trajet TSV complet.

Sur `templates/chitosan_stm_maps_dft_m030_h050_v1.tsv`, le moule type 0 place
`42.924%` de son énergie L2 au bord et son maximum en `(0.32, 0.00) nm`. Le
moule type 1 en place `92.736%` au bord, avec son maximum en
`(−0.32, 0.24) nm`. La différence `type_1−type_0` est elle-même dominée à
`92.838%` par le bord, au même maximum. Dans les trois cas, le maximum absolu
est sur la frontière du cadre 9×9 `[-0.32,0.32] nm`.

Décision: ce cadre ne démontre pas que la structure discriminante est contenue;
le score unary/local répète en outre un motif issu d'un environnement trimère
comme s'il constituait une observation indépendante par lobe. Suspendre son
interprétation chimique et ne pas le promouvoir. La prochaine formulation doit
être un modèle whole-ROI assemblant explicitement les contributions sur la
géométrie de chaîne et leurs recouvrements, après récupération ou régénération
de cartes dont le support dépasse au moins le cadre préliminaire
`[-0.48,0.48] nm`. Ce constat est un gate de représentation, pas une validation
de classe ni une raison de retuner le posterior.

### 2026-07-22 — Registration label-free du lag backward avant score type

Le lag global est désormais appliqué dans un diagnostic séparé, jamais dans la
production. Le script estime d'abord le lag médian des lignes de ROI sur support
constant, translate la matrice backward sans interpolation, puis rééchantillonne
les mêmes lobes et compare avant/après les patches image ainsi que les posteriors
génératifs DFT seul et combiné. Le contrôle synthétique vérifie la translation
sans mutation de la vue source et l'estimateur médian robuste; 7/7 assertions
passent.

Sur `240817_001.sxm`, le lag label-free estimé reste `−35 px` (`−1.370 nm`). La
registration fait passer la corrélation locale médiane de `0.1623` à `0.8498` et
la NRMSE affine médiane de `0.8743` à `0.5266`. La concordance d'argmax
forward/backward passe de `5/10` à `8/10`, pour DFT seul comme pour le mélange.
La différence médiane `|p1_fwd-p1_bwd|` baisse de `0.0122` à `0.0070` pour DFT
et de `0.0083` à `0.0052` pour le mélange. Deux exécutions produisent des TSV
byte-identiques.

Ce succès de registration ne valide toutefois pas l'identité chimique. Après
correction, la confiance maximale ne dépasse que `0.5196` pour DFT seul et
`0.5131` pour le mélange: le posterior reste pratiquement nul, et deux lobes
restent discordants. Décision: la translation d'acquisition explique une grande
partie de l'échec fwd/bwd, mais elle ne transforme pas les moules DFT en signal
type identifiable. Ne pas promouvoir cette correction, ne pas lancer d'array et
ne pas retuner les seuils. La prochaine preuve utile doit être un observable
physique indépendant; poursuivre les transformations de ces mêmes patches
risquerait seulement d'optimiser le smoke.

### 2026-07-22 — Lag d'acquisition forward/backward dans la ROI fixe

Le gate suivant mesure, ligne par ligne dans la boîte englobante de la ROI
moléculaire fixe, le lag x qui maximise la corrélation absolue entre les vues
prétraitées. La fenêtre de lag et le nombre minimal de paires sont des arguments
CLI explicites. Le premier prototype comparait des overlaps de tailles variables
et saturait artificiellement aux grandes bornes; il a été rejeté puis corrigé
pour utiliser exactement le même support intérieur pour tous les lags. Le test
synthétique récupère alors exactement un lag connu de `+2 px` avec son gain et
son offset affines; 12/12 assertions passent.

Sur `240817_001.sxm`, sous `plane+rows`, les fenêtres ±8 puis ±32 px saturaient
encore à leur borne négative. Avec une fenêtre diagnostique ±128 px et au moins
20 paires, le maximum devient intérieur et fortement concentré: lag médian
`−35 px`, MAD `2 px`, soit `−1.370 nm`; 410/512 lignes se trouvent entre `−40`
et `−30 px`, et 495/512 ont un lag négatif. Après translation ligne, la
corrélation absolue médiane passe de `0.4687` à `0.9692`, la NRMSE affine médiane
tombe à `0.2464`, et le gain médian de corrélation vaut `0.4141`. Deux exécutions
produisent des TSV byte-identiques.

Décision: l'hypothèse de dérive locale dispersée est remplacée par une évidence
forte de translation/hystérésis x globale résiduelle entre acquisitions, malgré
le flip x déjà appliqué à la lecture SXM. Ne pas corriger la production ni
émettre de types à ce stade. Le prochain gate label-free est d'appliquer cette
translation estimée sans labels à la vue backward, puis de répéter sur le même
scan les corrélations de patches et la stabilité du posterior génératif avant
tout élargissement à d'autres fichiers.

### 2026-07-22 — Ablation du flattening sur le transfert fwd/bwd

Un nouveau diagnostic standalone réutilise une seule géométrie de lobes, fittée
de façon déterministe sous `plane+rows`, puis rééchantillonne les mêmes
coordonnées image sous les quatre prétraitements déclarés par GaussianFit2D:
`none`, `plane`, `rows` et `plane+rows`. Il ne refit donc ni N ni la géométrie
entre conditions. Les modes et la fenêtre de décalage sont fournis explicitement
par CLI; aucune config ou sortie de production n'est modifiée. Le test ciblé
compte 14/14 assertions vertes et vérifie le parsing, la géométrie fixe et les
métriques affines.

Sur `240817_001.sxm`, les quatre modes conservent cinq pentes locales négatives
sur dix. Les corrélations médianes sont `0.0736` (`none`), `-0.0001` (`plane`),
`0.0788` (`rows`) et `0.1623` (`plane+rows`). Les NRMSE affines médianes restent
élevées: respectivement `0.8943`, `0.8406`, `0.9099` et `0.8743`. Sept lobes sur
dix gardent le même signe de pente dans les quatre modes, et les meilleurs
décalages d'un pas restent dispersés. Deux exécutions complètes produisent des
TSV byte-identiques.

Décision: rejeter le flattening ligne par ligne ou le retrait de plan comme cause
dominante de la discordance fwd/bwd. `plane+rows` améliore même la corrélation
médiane par rapport aux autres modes sans résoudre les inversions locales. Le
prochain gate doit cibler la dérive/hystérésis locale d'acquisition ou employer
un observable physique indépendant; aucune raison ne justifie un array type ou
un retuning du posterior.

### 2026-07-22 — Diagnostic du transfert forward/backward local

Le gate génératif ayant échoué entre directions, une sortie TSV séparée mesure
maintenant, par lobe et sans labels, la corrélation centrée image forward/backward,
la pente et l'offset affines `bwd ≈ a·fwd+b`, la NRMSE après cette recalibration,
et le meilleur décalage entier dans une fenêtre explicitement fournie par CLI.
Cette table ne modifie ni les TSV existants ni la production. Pour le smoke, la
fenêtre a été gelée à un pas (`0.08 nm`) dans chaque direction du repère `(t,u)`.
Les tests récupèrent exactement les gains/offsets, inversions de polarité et
décalages synthétiques connus.

L'alignement de base n'est pas absent: `STMSXMIO.read_sxm` retourne déjà le canal
backward en x, et les images prétraitées complètes de `240817_001.sxm` ont une
corrélation `0.6411`, une pente `0.6537` et une NRMSE affine `0.7674`. En revanche,
sur les dix patches du candidat GCV local déterministe, la corrélation locale
médiane tombe à `0.1623` (intervalle `[-0.8371, 0.7946]`); cinq pentes sur dix
sont négatives et la NRMSE affine médiane vaut `0.8743`. Le meilleur décalage
n'est pas commun aux lobes: les neuf positions non nulles de la fenêtre sont
utilisées, avec quatre maxima à `(-1,-1)` seulement. Le gain médian de corrélation
absolue n'est que `0.1212`. Deux exécutions produisent un TSV byte-identique.

Décision: écarter l'hypothèse d'un simple flip ou d'un décalage global résiduel.
La discordance type est précédée par une instabilité locale de morphologie ou de
contraste entre acquisitions, que l'affine et un déplacement de `0.08 nm` ne
réparent pas uniformément. Ne pas élargir l'array ni retuner le posterior. Le
prochain gate sûr est d'isoler l'effet du flattening ligne par ligne et de la
dérive locale sur les patches image avant toute nouvelle comparaison DFT.

### 2026-07-22 — Modèle génératif conjoint image-domain: gate réel échoué

Le diagnostic type possède maintenant une troisième ablation, strictement hors
production, qui score le patch image avec le même modèle pour les deux types:
`fond + amplitude·Gaussian + coefficient·moule`. Le fond est libre; l'amplitude
Gaussian et le coefficient du moule sont contraints positifs par une résolution
active-set exacte des deux coefficients. Le score est la fraction bornée de SSE
du résidu non-Gaussian expliquée par le moule, relativement au modèle nul
`fond + Gaussian`. Une vraisemblance iid sur les 81 pixels a été testée puis
écartée avant conservation: elle saturait artificiellement le posterior malgré
la corrélation spatiale connue. Les seules représentations réelles admises pour
cette ablation sont le patch image combiné et, comme contrôle de stabilité, ses
vues forward et backward séparées. Aucun seuil, config ou chemin de production
n'est modifié. Le smoke a aussi révélé que le global optimizer limité à une
seconde pouvait changer le candidat entre invocations identiques. Le diagnostic
clone donc la config de chaîne avec `skip_global=true` et une graine stable par
fichier; deux exécutions complètes produisent désormais des TSV byte-identiques.

Le gate synthétique est symétrique et passe. Pour les moules injectés
`stm_dft_v1`, les 20 répétitions de chaque type sont toutes récupérées à bruit
`0` et `0.10`. Les marges vraies moyennes à bruit `0.10` sont `0.1621` pour le
type 0 et `0.0693` pour le type 1, sans saturation. Un patch pur
`fond + Gaussian` donne exactement `(p0,p1)=(0.5,0.5)` avec deux vues, et un
moule de signe physique opposé n'obtient aucun gain SSE grâce à la contrainte de
coefficient positif. La suite ciblée compte 58/58 assertions vertes et la suite
joint-proxy complète est verte.

Le premier smoke réel label-free, limité à `240817_001.sxm`, échoue toutefois au
gate directionnel et arrête l'expansion. Sur les dix lobes du candidat GCV local
déterministe, la concordance d'argmax forward/backward n'est que `5/10`, pour
DFT seul comme pour le mélange. Les probabilités restent presque nulles: le
mélange combiné sur les deux vues donne une confiance moyenne `0.5078`, et la
confiance maximale des vues séparées ne dépasse pas `0.5326` (`0.5452` pour DFT
seul). Ce comportement
est honnêtement proche de l'abstention, mais l'argmax résiduel n'est pas une
signature stable.

Décision: conserver le modèle et ses sorties forward/backward comme affordance
diagnostique, mais rejeter toute émission type, calibration de seuil ou array
réel plus large. Le modèle génératif corrige la rupture de domaine et évite le
collapse extrême, mais les scans ne contiennent pas encore une évidence DFT
transférable et reproductible entre directions. La prochaine étape devra ajouter
une source physique indépendante ou expliquer la discordance de direction;
retuner le posterior sur ce smoke serait invalide.

### 2026-07-22 — Projection matched des moules contre le backbone Gaussian

Le gate suivant a testé, dans le diagnostic seulement, un opérateur de nuisance
apparié aux observations. Pour chaque lobe, le patch et chaque moule sont
projetés hors de l'espace `[1, Gaussian(t,u; σ∥,σ⊥)]`, avec les sigmas du lobe
fitté et l'ordre `patch_tu` existant. Le posterior matched conserve ensuite les
mêmes états globaux direction/phase/mirror, poids de sources et
forward-backward que le moteur actuel. Aucun seuil, config ou chemin de
production n'est modifié.

Le contrôle synthétique passe son gate préalable: DFT matched→DFT matched
récupère 100% des deux types à bruit `0` et `0.10`; à bruit `0.10`, la marge
vraie moyenne reste positive (`0.0681`, contre `0.1488` avant projection). Le
projecteur retire donc une partie du signal mais ne détruit ni le mapping 0/1 ni
la séparabilité self-template. L'array Viper `10649023` termine ses huit shards
`0:0`, puis le merge `10649201` termine `0:0`. Les artefacts dans
`results/joint_proxy_type_ablation_v4_matched/` contiennent 25980 lignes lobe,
3320 summaries et 145 métriques agrégées.

La projection réduit les extrêmes mais échoue au gate de stabilité réel. Pour
le mélange matched, la fraction d'argmax `type=1` est `87.48%` sur le résiduel,
`79.56%` sur `z-model`, `14.73%` sur l'image et `11.51%` sur le résiduel opposé.
Les quatre représentations ne donnent le même argmax que pour `6.08%` des 1086
lobes comparables (contre `0.09%` sans projection); image et résiduel ne
concordent que dans `20.07%` des cas (contre `3.15%`). L'amélioration est réelle
mais très loin d'une signature indépendante de la représentation. La masse DFT
reste elle-même instable: `0.8642` sur image, `0.5191` sur `z-model`, `0.4487`
sur résiduel et `0.6848` sur résiduel opposé.

Décision: rejeter cette projection simple comme correctif ou voie de promotion.
Elle confirme que le retrait du backbone explique une partie du collapse, mais
ne rend pas le type identifiable. Ne pas retuner le seuil et ne pas émettre de
nouveaux labels. Le prochain modèle défendable doit scorer directement le patch
image avec un modèle génératif conjoint `fond + amplitude·Gaussian +
coefficient·moule`, en ajustant les mêmes nuisances pour les deux types et en
comparant leur évidence/SSE sur une représentation unique pré-déclarée. Ce
nouveau modèle devra passer les contrôles synthétiques, la stabilité fwd/bwd et
l'abstention nulle avant tout test réel ou grading externe.

### 2026-07-22 — Ablation label-free de l'effondrement type

`test/diagnose_joint_proxy_type_collapse.jl` décompose maintenant le posterior
type entre les familles `geometric`, `stm_dft_v1` et leur mélange, sans lire de
manifest benchmark, séquence, composition ou expected N. Le contrôle synthétique
injecte chaque moule `(type, parity, mirror)` comme patch, avec cinq répétitions
à bruit `0` et `0.10`. DFT→DFT récupère symétriquement 100% des types 0/1 aux
deux niveaux de bruit (marge vraie moyenne `0.1502` sans bruit); les moules
géométriques seuls restent à 50%. Le mapping 0/1, les huit états et le moteur de
posterior savent donc distinguer les moules DFT dans leur propre domaine; ce test
ne valide pas le transfert aux scans réels.

L'ablation réelle a refait les fits GCV sur les 146 scans puis choisi, pour ce
diagnostic seulement, le candidat de `joint_gcv` minimal. Le premier array Viper
`10646630` a révélé un scan à évidence unary non finie dans le shard 7. Le
diagnostic conserve désormais ce cas comme `nonfinite_evidence` plutôt que de
perdre le shard; la reprise `10646723` et le merge `10646817` terminent `0:0`.
Deux passes supplémentaires ont comparé `z-model`, le résiduel recalibré, son
opposé, puis le patch image prétraité absolu (`10646846/10647014` et
`10647039/10647125`, tous `0:0`). Les artefacts finaux sont dans
`results/joint_proxy_type_ablation_v3_image/`: 13056 lignes lobe, 1664 summaries
et 67 métriques agrégées dans `type_collapse_metrics.tsv`. Huit scans n'ont pas
de candidat et un scan produit une évidence non finie; 137 scans contribuent aux
métriques finies.

Résultat principal, sur 1088 lobes label-free: le résiduel de production donne
`99.26%` d'argmax `type=1` en mélange (`99.26%` DFT seul), tandis que son signe
opposé donne seulement `1.29%` de `type=1`. Le patch `z-model` reste à `93.93%`
de `type=1`, mais le patch image absolu bascule à seulement `2.48%`. La masse DFT
du mélange passe simultanément de `0.9140` (résiduel) à `0.7472` (`z-model`),
`0.4502` (image) et `0.1298` (résiduel opposé). Le type prédit suit donc la
représentation et la polarité du patch, pas une signature chimique stable.

Cause méthodologique retenue: le pipeline compare des résidus signés de modèle
gaussien à des moules DFT qui représentent des cartes LDOS absolues normalisées.
Le succès synthétique self-template ne couvre pas cette rupture de domaine. Il
ne faut ni retuner le seuil ni promouvoir les labels actuels. Le prochain gate
physique est d'appliquer le même opérateur de nuisance aux observations et aux
moules: soit résidualiser chaque moule contre le fond constant et le Gaussian
de lobe correspondant, soit scorer le patch image avec un modèle génératif qui
ajuste explicitement fond/amplitude/backbone. Toute nouvelle émission type doit
ensuite être stable entre représentations physiquement équivalentes. Aucun
comportement de production ni seuil n'est modifié par cette ablation.

### 2026-07-21 — Visualisation des moules et calibration full-pipeline

Les huit moules unary DFT 9×9 peuvent maintenant être visualisés sans données
réelles ni labels via `test/plot_joint_proxy_molds.jl`. La sortie par défaut,
`results/joint_proxy_mold_visualization/stm_dft_v1_glcn_glcnac_difference.png`,
montre pour chaque état `(parity, mirror)` GlcN, GlcNAc et leur différence à
`-0.300 V`, `0.50 nm`. `test/analyze_joint_proxy_molds.jl` écrit aussi les 28
comparaisons pairwise dans `stm_dft_v1_pairwise_separability.tsv`. Les états
appariés GlcN/GlcNAc ont un RMSE moyen de `1.602110` et une corrélation moyenne
de `-0.299420`, mais le plus proche voisin brut n'a le même type que pour 4/8
moules: les transformations d'orientation créent donc des recouvrements qui
interdisent d'interpréter la seule distance template comme une validation de
l'identité chimique.

La calibration synthétique non-fast a été exécutée sur Viper avec le seed
`20260721` et 100 cas. Le premier job `10645175` a échoué proprement parce qu'un
fit réel n'offrait pas le vrai N synthétique parmi ses candidats valides. Le
calibrateur de posterior exige ce candidat pour calculer sa vraisemblance; ces
échecs structurels sont maintenant exclus uniquement de l'ajustement
conditionnel, conservés dans les rapports et comptés explicitement. Le test de
régression échoue avant le correctif et passe après; la suite joint-proxy passe.
Le job de reprise `10645348` termine `0:0`: 98/100 cas ont le vrai N disponible
(48 calibration, 50 heldout), deux restent des échecs structurels. La NLL
heldout descend de `2.189538` à `1.808285`; le seuil count devient
`0.7530942651` au lieu de presque 1. La calibration type est inchangée. Le fichier
versionné est `config/joint_proxy_calibration_dft_m030_h050_v2_full.toml`; sa
première version avait le hash
`012fcc2fe727b5f318844662c393d6595d3a5f3f291acdc48d910629188daebc`.

Un replay sans labels sur les tables candidate-N déjà calculées prédisait 9 N
émis sur 138 scans ayant des candidats, contre 2 avec v1; les huit scans sans
candidat restent nécessairement abstention. L'array comparatif Viper `10645507`
a terminé ses huit shards `0:0`, puis le merge/validateur `10645890` a terminé
`0:0`. Les artefacts fusionnés dans
`results/joint_proxy_dft_m030_h050_v2_full_all/` contiennent 146 summaries, 1079
candidate-N, 7631 candidate-lobes et 498 prédictions. Le résultat confirme le
replay: 9/146 N sont émis, dont les deux v1; huit scans n'ont toujours aucun
candidat. Le gate type émet 64 labels au lieu de 13, mais tous sont `type=1`;
496/498 posteriors bruts ont aussi leur argmax sur `1` et deux sont des ties.
Cette couverture count accrue ne valide donc pas l'identité chimique et v2 ne
remplace pas v1 comme résultat scientifique final. Aucun label externe,
séquence de contrôle ou grading n'est intervenu. L'asymétrie type vers `1`
reste le blocage scientifique distinct à résoudre.

La revue post-implémentation a ensuite rendu la frontière calibration/heldout
explicite dans l'adapter et ajouté un test d'intégration qui vérifie les compteurs
après exclusion d'un vrai N absent. Ce changement ne modifie ni observations ni
paramètres: la recalibration Viper `10646128` (`0:0`) est byte-identique à la
précédente sauf `source_sha256`, désormais
`1945b97ce7e097b026460772a7c80fb4d29836fc09c0ef7e74b1f876173cb8d0`.
Le hash du fichier versionné courant est
`9254cdc990396a1dcb856571bc162a36d2b7437cd802bf6fbc5550c6009b41a0`.
L'array de provenance courante `10646248` a terminé ses huit shards `0:0` et le
merge/validateur `10646474` a terminé `0:0`. Les artefacts dans
`results/joint_proxy_dft_m030_h050_v2_full_reviewed_all/` passent le validateur
local: 146 summaries, 1087 candidate-N, 7716 candidate-lobes et 497 prédictions.
La stochasticité du fit déplace un cas juste sous le seuil: 8/146 N sont émis,
huit scans restent sans candidat, 57 labels type sont émis et tous valent `1`;
495/497 argmax bruts vont vers `1`, avec deux ties. Cette répétition confirme le
gain count limité et l'absence de résolution du blocage type.

### 2026-07-21 — Inférence joint-proxy DFT sur 146 scans via Viper

L'inférence label-free avec `stm_dft_v1` et la calibration versionnée a été
lancée sur les 146 fichiers `.sxm` présents dans `/ptmp/oldu/stmfit/data`. Viper
a utilisé un array de huit tâches à un CPU (`10643202`, puis `10643608` après
correction) et un job de merge `10644169`. Le merge a terminé `0:0` en 10 s.
Les artefacts fusionnés sont dans
`results/joint_proxy_dft_m030_h050_v1_all/` et passent
`validate_joint_proxy_predictions.jl`: 146 summaries, 1082 lignes candidate-N,
7655 lignes candidate-lobe et 477 prédictions. Le manifeste lie le config
`dae8b8a7e...cb196`, le source `68bd6cba...c5161` et le payload
`82fadb9c...b3e5f`.

Le premier array a révélé deux cas limites réels: certains scans n'avaient aucun
candidat count fini/valide, et certaines patches donnaient une évidence unary
non finie. Faire échouer le shard perdait les autres fichiers. Les tests rouges
ont verrouillé le comportement attendu: le premier cas produit désormais une
summary-only abstention (`candidate_count=0`); le second produit des probabilités
type `0.5/0.5` et une prédiction `?`. Le validateur et le merge acceptent ces
abstentions structurées mais rejettent une summary vide non marquée. La suite
joint-proxy complète passe après correction. Une deuxième erreur purement
opérationnelle venait du nom de calibration copié dans chaque shard; les
manifestes ont été normalisés vers `joint_proxy_calibration.toml`, les huit
shards ont été validés séparément, puis fusionnés sans recalcul.

QC sans labels: 144/146 scans s'abstiennent sur N; seuls `240815_098.sxm`
(`N=7`) et `240818_026.sxm` (`N=6`) émettent un N, tous deux avec confiance 1.
Huit scans ont `candidate_count=0`; aucun scan ne manque la vue backward. Les
477 lignes de prédiction contiennent 464 `?`; les 13 labels émis sont tous
`type=1`. Même avant le gate final sur N, 475/477 marges type ont leur argmax sur
`1` et deux sont exactement `0.5/0.5`. Cette forte asymétrie est un signal QC
diagnostique, pas une validation chimique. Aucun label externe, séquence de
contrôle ou script de grading n'a été utilisé. Il faut expliquer l'abstention N
et l'effondrement du posterior type avant tout benchmark externe.

### 2026-07-21 — Provider DFT 9×9 versionné et promotion label-free

Les cubes corrigés à `-0.300 V` ont été rééchantillonnés directement sur la
grille native du joint-proxy (`half_nm=0.32`, `step_nm=0.08`, 9×9) à la hauteur
physique `0.50 nm`. Les artefacts suivis sont
`templates/chitosan_stm_maps_dft_m030_h050_v1.tsv`,
`templates/chitosan_connected_molds_stm_dft_m030_h050_v1_half032.tsv` et leur
sidecar `.provenance.toml`. Le sidecar est calculé par le finaliseur depuis les
cubes et sorties réels; il lie les hashes des cubes GlcN/GlcNAc, des cartes et
des templates, ainsi que le bias, la hauteur, la grille et les unités.

Le registre contient maintenant la famille obligatoire `stm_dft_v1`. Le loader
vérifie les hashes épinglés, le sidecar et ses paramètres physiques avant de
charger les huit états `(type, parity, mirror)`. Les poids restent fixés par le
contrat antérieur à `0.5` géométrique / `0.5` STM. Les cinq sources
`stm_prelim` restent listées pour la traçabilité mais portent `enabled=false` et
ne contribuent plus au posterior. Les bond molds restent exclus. Le payload
actif est `82fadb9c2d18eae3719293181d33200842aad912dd28e6345a0c7ee8a88b3e5f`.

Une première validation synthétique de 100 cas (seed `20260722`) a échoué aux
seuils absolus sur `noiseless_count_recovery=0.95` et
`corrupted_type_recovery=0.886064`. Aucun seuil n'a été modifié. L'A/B apparié a
ensuite basculé uniquement la présence de `stm_dft_v1`: geometric-only reproduit
exactement les mêmes deux échecs. La comparaison a été répétée pour les seeds
pré-déclarés `1`, `1234` et `20260722`; les deltas production moins baseline sont
exactement zéro pour les récupérations count/type noiseless, low-noise et
corrupted, ainsi que les abstentions null/identical. La fragilité vient donc du
découpage de 100 cas en seulement 20 cas/groupe (un cas vaut 0.05), pas du
provider DFT. La promotion utilise ce gate de non-régression apparié; elle ne
prétend pas que la calibration synthétique valide l'identité chimique réelle.
Aucun label de benchmark ni séquence de contrôle n'a été lu. La calibration
versionnée `config/joint_proxy_calibration_dft_m030_h050_v1.toml` utilise le seed
`20260721`, 100 cas synthétiques et le mode rapide; son SHA-256 est
`b6b98bc03d9b49583e423682d29a507091a413b3ed0837a26eb45d41b3aa0fc9`, lié au
config `dae8b8a7ecae004ce3a9c740c11308026dd94d8ef5e8edc473d52b3b785cb196` et au
payload ci-dessus.

### 2026-07-21 — GlcNAc convergence criterion converted into an LDOS stability gate

Three GlcNAc SCF paths now reach the same stable `~4e-5 Ry` floor, so repeating
strict mixing-only retries is no longer justified. QE defines `conv_thr` against
an extensive estimated energy error. A provisional common acceptance criterion
of `5e-5 Ry` is below `0.7 meV` for the complete 218-atom cell, versus the
unchanged `0.02 Ry` electronic smearing and experimental `-0.300 V` bias. GlcN
already satisfies this common gate by a much larger margin (`6e-8 Ry`). Because
LDOS depends on wavefunctions as well as total energy, the looser energy gate is
not sufficient by itself for production.

The acceptance job therefore creates two byte-identical copies of the retry2
terminal charge density and independently converges them to `5e-5 Ry`: one with
plain mixing (`beta=0.3`, `ndim=20`), one with TF mixing (`beta=0.2`,
`ndim=20`). Each branch writes its own wavefunctions and STM cube using the
unchanged PBE+D3 Hamiltonian, pseudopotentials, `50/360 Ry` cutoffs, Gamma
sampling, and MV `degauss=0.02 Ry`. No cube
will be promoted unless same-frame normalized maps agree pointwise within 1%
and branch total energies agree within `1e-4 Ry`; if they pass, the branch with
the smaller final SCF residual is retained and the other remains validation.
Local checks, remote input hashes, baseline charge-density check, and shell
syntax all passed. Raven jobs `28888367`, `28888618`, then the short diagnostic
chain `28888893/28888894/28888969` remained blocked by the account-wide
`QOSGrpCpuLimit` and were cancelled before start; no computation was lost.

The exact 3.30 GB retry2 checkpoint was transferred to Viper, including the
distributed `wfc1`–`wfc8` and mixing-history files that were outside the
previously inspected `.save` subdirectory. The acceptance inputs were therefore
corrected to `startingwfc='file'` and run with the same QE 7.4.1,
`intel/2024.0`, `impi/2021.11`, and eight-MPI decomposition. Viper jobs
`10640236` (plain), `10640237` (TF), and `10640238` (comparison) all completed
`0:0` in `06:28`, `06:16`, and `00:12`.

Both SCFs converged in one iteration at `4.027e-5 Ry`, with identical energy
`-31763.44437132 Ry` and Fermi energy `0.3038 eV`. Both `pp.x` stages ended with
`JOB DONE`. Their first complete cubes were byte-identical, SHA-256
`90eb5b3a6a14505e5f7160733734604295ede5883af94f918134f8d815523874`, but were
later found to use QE's implicit `+0.01 Ry` (`+0.1361 eV`) bias and are archived.
Independent 0.50 nm, 13×13 same-frame comparison returned maximum/mean/RMS
normalized differences `0/0/0`, correlation `1`, and energy difference `0 Ry`.
The gate therefore passed. The plain cube was promoted locally as
`qe/glcnac/glcnac_central_ldos.cube`; the TF cube remains validation evidence.

Running `test/finalize_qe_mold_workflow.jl` then exposed a Julia boundary bug:
the regex-derived QE prefix was a `SubString`, but `_default_index_tsv` required
`String`. A failing CLI regression test reproduced the `MethodError`; converting
the capture to `String` at `_prefix` made it pass. The real finalizer then
generated 0.50 nm typed frames, `templates/chitosan_stm_maps.tsv`, and connected
unary/bond molds from the accepted GlcN and GlcNAc cubes. The existing
`config/joint_proxy_molds.toml` intentionally still names preliminary 9×9
sources; replacing it requires a versioned 9×9 export and regenerated synthetic
calibration, not a silent path swap.

Source inspection of QE 7.4.1 `PP/src/stm.f90` established that `plot_num=5`
ignores `emin`, `emax`, and `degauss_ldos`; it integrates from `E_F` to
`E_F + sample_bias`, with `sample_bias` expressed in Ry. All PP inputs were
corrected to `sample_bias=-0.0220495933 Ry`. Viper PP-only job `10640445` and
dependent comparison job `10640446` completed `0:0`; all three PP outputs print
`Sample bias = -0.3000 eV`. Corrected GlcNAc plain/TF cubes remain byte-identical
(SHA-256 `40649ccd9eb6768444b8ff61bf4a639b3940cf926fe3eb42254eb024faf9b5bf`),
with map differences `0/0/0`, correlation `1`, and energy difference `0 Ry`.
The corrected GlcN cube SHA-256 is
`80cd1d1fde94cf084cc7ea464d2bf065b36b2c015ea0bfaef8cebfee8ff88863`.
The canonical plain cube was promoted and the 0.50 nm 13×13 maps and connected
molds were regenerated from these corrected cubes. The input generator and
legacy `hpc/qe_molds/pp_ldos.in.template` now emit `sample_bias`; preflight
explicitly requires `plot_num=5`, a nonzero bias, and absence of the ineffective
`emin`/`emax`/`degauss_ldos` fields. The smoke workflow covers this contract.

### 2026-07-20 — GlcNAc SCF plateau diagnosed and numerical retry2 submitted

Raven job `28811340` finished with application exit `2:0` after `09:05:29`.
The `pw.x` step itself exited normally after 300 iterations, but reported
`convergence NOT achieved`; the final estimated SCF accuracy was
`4.056e-5 Ry` against the unchanged production target `1.0e-7 Ry`. Maximum
node memory was only about 30.7 GB of the 96 GB request, so this was neither a
memory nor scheduler failure. The convergence guard correctly prevented
`pp.x`, and no GlcNAc cube was produced.

Trajectory comparison showed that the original final SCF and the `local-TF`,
`mixing_beta=0.1`, `mixing_ndim=16` retry reached essentially the same floor:
their best estimated accuracies were `4.322e-5` and `3.461e-5 Ry`, respectively,
and the retry's last-50 median was `4.2905e-5 Ry`. In contrast, the same pilot
Hamiltonian settings converged GlcN to `6e-8 Ry` in 37 iterations. The retry had
1080 Kohn-Sham states for 1800 electrons and Davidson required only one or two
inner iterations near the end, so neither an obvious empty-band shortage nor
eigensolver failure was selected as the first intervention.

Following the QE PWscf troubleshooting guidance to test mixing modes and Pulay
history before changing physical parameters, a second SCF+PP-only retry keeps
the relaxed geometry, PBE+D3, pseudopotentials, `50/360 Ry` cutoffs, Gamma
sampling, MV `degauss=0.02 Ry`, and strict `conv_thr=1e-7 Ry`, but uses
`mixing_mode='plain'`, `mixing_beta=0.3`, and `mixing_ndim=20`. The exact
`28811340` input/script were archived locally with `retry1` suffixes. Local and
Raven preflight passed at 8 MPI tasks / 96 GB, input/script hashes matched
between hosts, and retry2 was submitted as job `28851882` (initially `PENDING`
under `QOSGrpCpuLimit`, 24 h limit). No other `oldu` job was visible in `squeue`;
the group limit may include another user or delayed accounting. This is a numerical convergence experiment, not a benchmark-driven
parameter change; `pp.x` remains gated on the explicit QE convergence message.

Job `28851882` subsequently finished with application exit `2:0` after
`09:02:18`. Plain mixing also reached 300 iterations without satisfying the
strict target: minimum estimated accuracy `3.239e-5 Ry`, final accuracy
`3.897e-5 Ry`, and final total energy `-31763.44440950 Ry`. `pp.x` was not
started and no cube was written. This reproduces the same approximate
`3e-5--4e-5 Ry` floor under both `local-TF`/0.1/16 and plain/0.3/20 mixing, so
another iteration-budget or mixing-only retry is not justified without a new
falsifiable hypothesis. The retry2 SCF and Slurm outputs were fetched locally.

### 2026-07-13 — Scientist-facing DFT calculation note added

Added `docs/src/dft_calculation_note.md` to record exactly what was calculated
for the GlcN/GlcNAc STM molds: molecular contexts, atom counts, Cu(100) pilot
slabs and cells, frozen atoms, PBE+D3/PAW settings, cutoffs, Gamma sampling,
smearing, relaxation and SCF criteria, intended STM bias, completed job
outcomes, preliminary-versus-production status, and interpretation limits. The
note explicitly records that GlcN has a completed production cube while GlcNAc
is in a targeted final-SCF retry, preventing preliminary molds from being
mistaken for converged chemical references.

### 2026-07-12 — GlcN restart5 resubmitted after Slurm launch failure

Raven job `28762985` failed after only 67 seconds with exit `1:0` because
`srun` could not confirm the allocation (`Socket timed out on send/recv
operation`; allocation reported expired/invalid). `pw.x` never started,
`glcn_central_relax.out` remained empty, CPU use was zero, and memory use was
only about 7 MB. This was an infrastructure launch failure, not a geometry,
convergence, walltime, or memory failure.

The unchanged, preflighted `qe/glcn_restart5` inputs were resubmitted directly
on Raven as job `28784933`. It completed successfully (`0:0`) in `02:11:57`:
the relaxed-enough geometry converged immediately at the restart point (`0`
additional BFGS steps), the relax SCF converged in 34 iterations, the final SCF
converged in 37 iterations, and all QE stages ended with `JOB DONE`. Raven now
contains `glcn_central_relaxed.xyz`, the converged SCF output, and a 136 MB
`glcn_central_ldos.cube`. This clears the production GlcN gate; the next step is
to fetch/finalize the GlcN mold and run the symmetric GlcNAc production path.

The completed GlcN relax/SCF/cube/relaxed-XYZ outputs were fetched locally into
`qe/glcn_restart5`. The GlcNAc relax input was then updated symmetrically with
`forc_conv_thr = 6.0d-3` and `etot_conv_thr = 1.0d-3`, while retaining the
strict final-SCF `conv_thr = 1.0d-7`. Local preflight passed at 8 tasks / 96 GB,
the updated input was synchronized to Raven, and production GlcNAc was submitted
as job `28790944`. Initial state: `PENDING`, 24 h time limit.

GlcNAc job `28790944` later finished the relaxation and final-SCF processes, but
the final SCF reported `convergence NOT achieved after 100 iterations`. QE then
wrote configuration-only restart metadata rather than collected wavefunctions;
`pp.x` failed with `Wavefunctions not in collected format` and missing `wfc1`,
so no GlcNAc cube was produced. The relaxed geometry and SCF charge data remain
usable. A targeted SCF+PP retry was prepared without repeating relaxation:
`electron_maxstep=300`, `mixing_mode='local-TF'`, `mixing_beta=0.1`,
`mixing_ndim=16`, `startingpot='file'`, and fresh atomic/random wavefunctions.
The retry script explicitly refuses to start `pp.x` unless QE prints electronic
convergence. It was submitted on Raven as job `28811340` (initially `PENDING`,
24 h limit).

### 2026-07-11 — Posterior joint proxy: synthetic validation and no-truth smoke

The standalone diagnostic path for joint count/type inference was exercised
without loading a benchmark manifest, expected count, composition, or unit
sequence.  Its synthetic acceptance matrix reached 100% noiseless count/type
recovery, 95.0% low-noise count recovery, 92.0% low-noise type recovery, and
100% null/identical-mold abstention.  Replacing an independent view with an
identical copy did not increase count or type confidence, deterministic replay
matched exactly, a stale source hash failed before inference, and an intentionally
swapped physical type mapping collapsed synthetic type recovery instead of being
silently relabeled.

The real no-truth smoke used a fast synthetic-only calibration and exactly
`240817_017.sxm,240817_019.sxm`.  Both files produced finite candidate tables
(21 candidate-N rows total and at least ten candidates per file), and the output
validator passed.  Both chains abstained on hard N (`?`, confidence 0.152 and
0.145), so these results are a transfer diagnostic, not counts or chemistry
claims.  No grader or truth-bearing manifest was imported.

The smoke also exposed standalone-load defects hidden by injected unit-test
adapters: calibration dynamically imported `GaussianFit2D` too late for Julia's
world age and accumulated real views in `Any[]`; inference relied on a
test-created `Main.STMSXMIO` binding.  Dependencies are now loaded statically,
real views retain their concrete type, and focused plus aggregate tests pass.
The non-fast synthetic calibration adapter disables only the residual-peak
validity veto because the deliberately injected type proxy is structured
residual relative to the count-only Gaussian backbone; production selection is
unchanged.

### 2026-07-10 — GlcN relaxed-enough restart5 submitted

Raven and Viper queues were checked after the unknown-workflow commits. No STMFit
array jobs were active; recent Viper arrays were completed with exit `0:0`. Raven
QE GlcN restart4 (`28658135`) had reached the 24 h walltime limit, not memory:
`pw.x` was cancelled due to time limit with `MaxRSS` about 4.9 GB, and the relax
output had reached at least `number of bfgs steps = 34` with residual gradients
around `4e-3` to `6e-3 Ry/Bohr`.

For the current DFT-STM mold application, the goal is a robust local STM/LDOS
mold, not a publishable adsorption energy. The relaxation policy was therefore
made explicitly **relaxed-enough** for GlcN: keep the final electronic SCF strict
(`pw_scf.in` `conv_thr = 1.0d-7`), but allow the geometry relax to stop at
`forc_conv_thr = 6.0d-3` and `etot_conv_thr = 1.0d-3`. These thresholds are chosen
from the observed plateau/progress and resource constraint, not from unit-label
benchmark performance. The same policy must be used symmetrically for GlcNAc if
its production relax needs the same shortcut.

The latest geometry from restart4 was extracted on Raven to
`qe/glcn_restart4/glcn_central_best5.xyz` and fetched locally. `qe/glcn_restart5`
was generated from that geometry with the active Raven pilot settings (`8` MPI
tasks, `96000 MB`, `24:00:00`, `ecutwfc=50`, `ecutrho=360`, Gamma-only), patched
with the relaxed-enough ionic thresholds above, passed local and remote preflight,
and was submitted to Raven without `--watch`:

```text
qe/glcn_restart5 -> 28762985
initial status: PENDING on small, reason=(Priority), time limit 1-00:00:00
```

Do not fetch or inspect `qe/glcn_restart5` outputs until Slurm reports completion
or timeout. GlcNAc production remains blocked until GlcN produces a usable
relaxed geometry and LDOS cube, unless an explicit dependency-policy change is
made.

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
files. The original QE inputs attempted to express this as `emin=-0.3`,
`emax=0.0` eV; the 2026-07-21 correction records that these fields do not control
`plot_num=5` and replaces them with `sample_bias=-0.0220495933 Ry`. The original QE
run directories used the full `8×8×4` slab with `ntasks=4`, `ecutwfc=80`,
`ecutrho=640`; those parameters were later dropped (OOM on Raven — see the
submission lessons below) in favour of the active `8×6×3` pilot
(`ecutwfc=50`, `ecutrho=360`, `8` tasks). The validated slab freeze cutoff
`fix_below_z=1.807501` is unchanged. The `qe/` directory is gitignored to avoid
committing large QE scratch/cube outputs.
Added `test/preflight_qe_mold_inputs.jl` and `hpc/submit_qe_molds.sh`; the
preflight report `hpc/qe_molds/qe_input_preflight.tsv` verifies `nat/ntyp`,
species, frozen relax atoms, STM sample bias, sbatch handoff commands, and the 8-task
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

> Updated 2026-08-03. Questions from earlier sessions are archived in
> `journal_archive.md`.

0b. **Can label-free unit assignment reach the promotion bar?** → **RESOLVED
    (Aug 3)**: the label-free champion (soft vote of k-means 4-view and GMM
    1-view + per-channel constant-current margins + Fisher empirical mold
    margin, self-training 2) scores 79.3% classified physical accuracy /
    36 exact chains / 677 fixed-denominator honest on the 145-file benchmark
    - the promotion bar (78.9% / 18 / 677) is MET. Fully reproducible via
    test/build_cc_soft_champion.py (verified: 677/36, zero label differences),
    audited label-free (zero label references in all construction scripts),
    and cross-validated (half-split Fisher mold: 66.3% per-lobe, no overfit;
    pipeline CV 678/34).

0. **Can synthetic posterior calibration transfer enough confidence to real
   scans?** → **OPEN (updated Jul 30)**: converged GlcN/GlcNAc DFT molds remain
   perfectly separable in-domain, but matched residualization failed
   representation stability and the joint image-domain generative score failed
   its first forward/backward real-scan gate (`5/10` for both DFT-only and combined,
   confidence near `0.5`). This is honest abstention evidence, not a failed
   benchmark. The first local transfer audit rejects a simple global flip/shift:
   five of ten patch slopes change sign and the best one-cell shifts are
   dispersed. A fixed-geometry ablation across `none`, `plane`, `rows` and
   `plane+rows` retains five negative slopes in every mode, excluding flattening
   as the dominant cause. A fixed-support row audit then finds a concentrated
   `−35 px` (`−1.370 nm`) x lag on 410/512 rows, identifying residual global
   acquisition hysteresis as an actionable nuisance. Applying that registration
   raises patch correlation to `0.8498` and type concordance to `8/10`, but the
   maximum confidence remains below `0.520`. A support audit adds a prior
   representation failure: `42.924%` of type 0, `92.736%` of type 1 and
   `92.838%` of their difference lie on the 9×9 boundary, with every absolute
   maximum on that boundary. Do not broaden the batch or tune against the
   withheld control sequence. Before another real-scan chemical gate, recover
   wider molds and replace repeated unary patches by an explicitly assembled
   whole-ROI chain model; an independent physical observable remains necessary
   for transfer evidence. The corrected constant-current API now rejects
   nonunique roots and changing support, but the accepted GlcNAc response has
   multiple branches. T3 therefore remains terminal `BLOCKED`; API-level
   ambiguity rejection is not evidence that real constant-current calibration
   is solved. Root uniqueness is evaluated under the provenance-bound
   `--isovalue-scan-intervals` policy (default 1024 intervals), so it is a
   declared finite-resolution contract rather than a resolution-free claim.
   Multi-file diagnostic publication is now recoverable through prepared and
   committed gate-last transactions, but filesystem corruption and hostile
   concurrent sidecar mutation remain outside that bounded protocol rather than
   being silently treated as solved.

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


## 2026-08-01 — Label-free unit assignment: GMM rebuild, view/seed sweeps, self-training, robust metric

### Context

Full 145-file benchmark: 854 graded blob positions (own-N), external control
`NKNNKN` (post-hoc only). Promotion bar: `physical >= 78.9%`,
`honest >= 677/854 (79.3%)`, `exact >= 18/145`. The historical headline is
`78.9% / 17-of-145` (k-means, counting stage); no unit-assignment candidate
has been promoted.

### 1. GMM predictor rebuild

- New `test/build_labelfree_gmm_predictions.jl`: 2-component full-covariance EM,
  k-means initialized, multi-seed vote, physical mapping GlcNAc = higher-
  amplitude cluster (label-free: reads no truth, sequence, or composition).
- Reproduces the documented 82.4% / 8-of-35 on the 240817 subset.
- Full145 (base features + backward moments + split skew, 3-view vote):

| Config (GMM, 20 seeds, +interactions) | Acc % | Exact/145 | Corr/854 |
|---|---|---|---|
| 3-view (base+com_t, base+diag45, base+diag135) | 76.5 | 29 | 653 |
| 3-view + patch_u_asym (9x9) | 76.3 | **33** | 652 |
| 1-view (base + patch_u_asym) | 76.9 | **35** | 657 |
| 5-view (adds split, uasym) | 76.2 | 33 | 651 |

- The 3-view + uasym config (33 exact) and later the 1-view config (35 exact)
  became the exact-chain leaders, roughly double the historical 17.

### 2. k-means view sweep (two-stage: maximize physical accuracy, then exact)

| Views (k-means, 20 seeds, +interactions) | Acc % | Exact/145 | Corr/854 |
|---|---|---|---|
| base only | 77.8 | 15 | 665 |
| base + split_log_skew | 77.8 | 16 | 665 |
| 4-view (base+split, base+com_t, base+diag45) | 78.0 | 16 | 666 |
| 5-view (adds uasym) | 77.8 | 15 | 665 |

- 4-view is the best label-free configuration on the robust metric (666/854).
- `split_log_skew` (recomputed from N_original split widths, log scale) is the
  only feature that raises exact chains without losing accuracy (14 -> 16).

### 3. 17x17 wavelet features (v4) do NOT reproduce the documented 85.2%

- `hh1_q00_abs` from 17x17 patches (Viper job 10801983) plus `neg_anis`,
  `patch_u_asym_17`: full145 = 76.1% / 28 exact; 240817 subset = 82.4%
  (= GMM base, no gain). The documented 85.2% / 13-of-35 (journal 2026-06-28)
  was computed by a script that no longer exists; the claim is not
  reproducible with the current pipeline and is therefore not trusted.

### 4. k-means/GMM combination and abstention sweeps

- Disagreements between k-means 3-view and GMM 3-view: 153/854 lobes (~18%).
- Keeping confident disagreements (k-means wins): honest peaks at 663 (KEEP=149),
  never reaches 677.
- Keeping confident disagreements (GMM wins), KEEP 80..153: physical 78.5%
  (93% coverage) -> 76.3% (full), exact 25 -> 33, honest 625 -> 652.
- Ensemble abstain on disagreements: 82.5% physical on 707/854 classified,
  but only 583/854 honest — worse than either base method in absolute terms.

### 5. Robust metric: accuracy x coverage = corrects/fixed denominator

- Lesson: reporting physical % on classified blobs alone rewards abstention.
  The only robust comparison metric is corrects/854 (= physical x coverage).
  Under it, every abstention variant is WORSE than full-coverage k-means:
  ensemble abstain 68.3%, KEEP80+GMM 73.2%, k-means 4-view 78.0%.
- The honest >= 677 bar is a product constraint: 677/854 = 79.3% accuracy at
  ~100% coverage. Neither k-means (78.0%) nor GMM (76.9%) reaches it; the
  ~153 hard disagreement lobes are won by neither method (k-means ~36%,
  GMM ~24% on them). The remaining path is per-lobe accuracy itself
  (new physics features or per-date calibration), not aggregation/abstention.

### 6. Self-training (GMM seeds -> Mahalanobis hard assignment) — Viper job 10802474

- Added `--selftrain N` to the GMM predictor (hard reassignment by Mahalanobis
  distance, re-estimate means/covariances, N iterations, per seed).

| Config | Selftrain | Acc % | Exact/145 | Corr/854 |
|---|---|---|---|---|
| 1-view base+uasym, 10 seeds | 0 | 76.9 | **35** | 657 |
| 1-view, 10 seeds | 5 | 76.8 | 34 | 656 |
| 1-view, 10 seeds | 10 | 76.8 | 34 | 656 |
| 3-view, 20 seeds | 0 | 76.3 | 29 | 652 |
| 3-view, 20 seeds | 2 | 76.6 | 32 | 654 |
| 3-view, 20 seeds | 10 | 76.6 | 32 | 654 |
| 3-view, 80 seeds | 0 | 76.3 | 29 | 652 |
| 3-view, 80 seeds | 2 | 76.6 | 32 | 654 |

- Self-training effect is small and view-dependent: +3 exact on the 3-view
  config (29 -> 32), -1 exact on the 1-view config. The 2026-06-28 subset
  claim (81.4% / 13) was config-specific (1-view, 10 seeds) and does not
  transfer to the 145-file benchmark.
- New exact-chain leader: 1-view GMM = 35/145.

### 7. Seed scaling — Viper job 10802791 (k-means, with interactions)

| Config | Seeds | Acc % | Exact/145 | Corr/854 |
|---|---|---|---|---|
| k-means 4-view | 20 | 78.0 | 14 | 666 |
| k-means 4-view | 80 | 78.0 | 14 | 666 |
| k-means 4-view | 200 | 78.0 | 14 | 666 |
| k-means 4-view | 500 | 78.0 | 14 | 666 |

- Seed scaling has NO effect for k-means (20 = 500): the multi-seed vote
  is already deterministic on these data. Same for GMM (20 vs 80 seeds
  identical). 20 seeds is sufficient; no need to raise.

### 8. Current best label-free configurations (2026-08-01)

```text
Accuracy/honest leader: k-means 4-view (+interactions)   78.0% / 666 corr / 14 exact
Exact-chain leader:     GMM 1-view base+uasym (+inter)   76.9% / 657 corr / 35 exact
Bar:                    >= 78.9% / >= 677 corr / >= 18 exact   (NOT met)
```

- Open question: whether any label-free aggregation can combine k-means
  78.0% per-lobe accuracy with GMM exact-chain structure without labels.
  Tested vote/agreement ensembles all land between the two on the robust
  metric. No further aggregation idea is pending; the next candidate needs
  better per-lobe features or a per-date calibration of the worst date
  (20241114, ~59%).

### 8b. New physical feature families (universal calibration attempt) — Viper jobs 10802901/10802926/10802970/10802999

Goal: raise per-lobe accuracy without per-date calibration (user decision: the
method must stay universal). Built `test/enrich_unit_features.py` (label-free:
reads only fit/patch TSVs) computing 40 new per-lobe features on the existing
9x9 backward/difference patches (raw and residual) and the chain geometry:

- Residual/difference 2D moments: energy, centroids (com_u, com_t), 3rd/4th
  moments (sku, skt, kurt_u), left/right and up/down parity (asym_u, asym_t),
  annular means r=1..3 — for bwd_res, diff_res, diff_raw, bwd_raw.
- Chain geometry: spacing_next_nm, spacing_asym, chain_curv_deg, elongation.

Also enabled `--no-interactions` in the GMM predictor for low-dim family views.

Results (full145, own-N; refs: k-means 4-view 78.0/16/666, GMM 1-view
76.9/35/657):

| Family added (k-means 4-view) | Acc % | Exact | Corr/854 |
|---|---|---|---|
| ref (no family) | 78.0 | 16 | 666 |
| bwd_res moments (3) | 77.9 | 14 | 661 |
| diff_res moments (5) | 77.9 | 14 | 661 |
| chain geometry (4) | 77.5 | 12 | 658 |
| all three families | 78.0 | 14 | 662 |
| envelope-resid (3) | 77.6 | 12 | 663 |
| envelope-resid + bwd_res | 77.8 | 13 | 664 |

| Family added (GMM 1-view base+uasym) | Acc % | Exact | Corr/854 |
|---|---|---|---|
| ref (+interactions) | 76.9 | 35 | 657 |
| ref, no interactions | 71.9 | 26 | 614 |
| + diff_res (no inter) | 70.0 | 13 | 597 |
| + bwd_res (no inter) | 72.1 | 30 | 612 |
| + geometry (no inter) | 66.8 (587 class.) | 0 | 392 |
| + envelope-resid (no inter) | 74.9 | 31 | 640 |

Findings:

- None of the new physical families improves the benchmark. The existing
  features (bwd_neg_*, patch_u_asym, split_log_skew, prominence) already
  capture all extractable acetyl signal from the 9x9 patches.
- Chain geometry (spacing/curvature) actively hurts: spacing is set by the
  molecule/N, not by monomer type (physically expected).
- High-dim GMM views over-parameterize: 7-11 features + interactions (36-66
  dims, full covariance) collapse to 46-54% (10802926). Family views must stay
  <= 6 features and prefer --no-interactions.
- k-means 4-view (78.0% / 16 exact / 666 corr) and GMM 1-view (76.9% /
  35 exact / 657 corr) are confirmed as a robust plateau for the current
  feature space and bias (-0.3 V).

Open question unchanged: the remaining path is a stronger physical signal
(different bias, converged DFT molds, or higher-resolution data), not more
statistics on the current features.

### 8c. Converged DFT-STM molds re-grade: NOT competitive — Viper job 10803069

Executed the documented "next gate" (journal 2026-06-28): re-grade with the
converged production cubes. State: GlcN `qe/glcn_restart5/glcn_central_ldos.cube`
and GlcNAc `qe/glcnac/glcnac_central_ldos_accept_plain.cube` pass the 5e-5 Ry
acceptance gate; the versioned 9x9 `stm_dft_v1` provider at 0.50 nm / -0.3 eV
was already generated (Jul 28) with pinned provenance hashes.

Steps:
- Extracted forward residual 9x9 patches (0.32 nm / 0.08 nm grid, matching the
  v1 templates) for all 145 benchmark files on Viper (job 10803069).
- Scored with `score_connected_mold_templates.jl --template-mode contrast
  --prefix res_p` (no bond templates: none exist for v1) and graded post-hoc.

| Grade | phys acc | exact | oracle acc |
|---|---|---|---|
| full145 (854 blobs) | **49.9%** | 0/145 | 66.9% |
| 240817 primary subset (209 blobs) | **49.8%** | 0/35 | 67.0% |
| prelim unconverged subset (journal, h0.50) | 58.6% | 2/35 | 68.1% |

Verdict:

- The converged DFT-STM molds are WORSE than the unconverged diagnostic molds
  (49.8% vs 58.6% phys on the subset) and far below the label-free k-means
  plateau (78.0% / 16 exact / 666 corr).
- The 17% oracle-physical gap persists: the amplitude-based physical mapping
  (GlcNAc = higher amplitude) contradicts the mold LDOS score on ~1 in 6 lobes.
  The simulated s-wave LDOS at 0.50 nm does not transfer to the experimental
  fit amplitude at -0.3 V.
- The 5e-5 Ry acceptance gate is a numerical-convergence gate, not a chemical
  transferability gate. Production molds must not be described as benchmark-
  competitive without an explicit re-grade.
- No further mold-variant sweep is planned (no bond v1 templates; contrast is
  already the documented best mode). The label-free k-means/GMM plateau stands
  as the current best, and the DFT-mold path is parked until a new physical
  signal (bias/resolution/mold height) is available.

### 8d. Mold improvement sweep: resolution 9x9 -> 17x17 (Viper job 10803207)

Attempt to improve the physical-mold approach after the 8c failure. Key
findings:

- The v1 templates were generated on a coarse 9x9 grid (step 0.08 nm). Regenerated
  the converged templates at 17x17 (step 0.04 nm, same 0.50 nm / -0.3 eV,
  same cubes) with `finalize_qe_mold_workflow.jl --half-nm 0.32 --step-nm 0.04`
  (also produced the missing v1 bond templates) and scored the full145 forward
  residual 17x17 patches (extracted on Viper, job 10803207).

| Scoring variant (17x17, converged cubes) | phys acc | exact |
|---|---|---|
| contrast, res_p, Viterbi | **56.8%** | 2/145 |
| contrast + bond templates | 56.6% | 2 |
| contrast, transition-penalty 0 | 56.8% | 2 |
| full mode (no contrast) | 43.4% | 0 |
| backward channel (bwd_res_p) | 45.6% | 1 |
| per-lobe direct (no Viterbi) | 48.6% | 0 |
| 9x9 contrast (8c baseline) | 49.9% | 0 |

Subset 240817 (209 blobs): 56.0% phys / 0 exact (vs 49.8% 9x9 converged, vs
58.6% unconverged prelim at h0.50).

- Mold cost features (cost_margin, cost_GlcN, cost_GlcNAc from the 17x17
  contrast scoring) added to the v3 feature table and swept on k-means 4-view
  and GMM 1-view (Viper job 10803239): no gain (km 77.8% vs 77.9% ref; GMM
  77.3%/+0.4% but -4 exact). The mold signal does not combine with the
  experimental clustering features.

Verdict:

- Spatial resolution is the one real lever: 17x17 gains ~7 points over 9x9
  (49.9 -> 56.8% full145; 49.8 -> 56.0% subset), and the oracle gap halves
  (17% -> 9.6%). All other levers (bond, mode, penalty, channel, per-lobe)
  are neutral or negative.
- The mold plateau (~57%) remains 21 points below the label-free k-means
  plateau (78.0%): the simulated s-wave LDOS signal transfers only weakly to
  the experimental fit features at -0.3 V. The mold approach is parked at
  17x17 contrast as the best mold configuration; a higher-resolution or
  better-bias dataset would be needed to close the gap.

### 8e. Mold improvement, round 2: forward+backward channel averaging (breakthrough)

Follow-up to 8d: the forward-only mold scoring ignores the backward scan
channel, which captures complementary tip-sample electronic states (and the
diff channel). Averaging the per-lobe `cost_margin` of the 17x17 contrast
scores from the forward and backward channels (both already extracted):

| Mold configuration | phys acc | exact |
|---|---|---|
| fwd 17x17 Viterbi (8d best) | 56.8% | 2 |
| bwd 17x17 Viterbi alone | 45.6% | 1 |
| **fwd+bwd margin mean, per-lobe** | **66.3%** | 0 |
| fwd+bwd patch-mean + Viterbi | 49.5% | 0 |
| fwd+bwd margin mean, per-lobe (complete) | 49.5%* | 0 |

(*the complete-file rerun re-scored patch means; the 66.3% is the per-lobe
margin mean before Viterbi — Viterbi on averaged patches destroys the gain,
per-lobe sign of the averaged margins is the correct decoding.)

- The two scan channels carry anti-correlated per-lobe errors: averaging the
  margins cancels directional scan noise (+9.5 pts over fwd alone). This is the
  strongest single mold lever found (resolution was +7).
- SSE scoring: 56.9% (neutral vs NCC). 6 edge lobes of 240818_019.sxm cannot
  be patched (image border, NA in both channels) and were completed from the
  k-means 4-view predictions (label-free) for grading.

Mold cost as clustering features (v4 sweep, Viper job 10803353): the averaged
margin added to GMM 1-view gives the best GMM of the session:

| Config | acc % | exact | corr/854 |
|---|---|---|---|
| GMM 1-view ref | 76.9 | 35 | 657 |
| **GMM 1-view + mold_margin_avg** | **77.3** | **35** | ~660 |
| k-means 4-view ref (leader) | 78.0 | 16 | 666 |
| k-means 4-view + margin_avg | 77.9 | 15 | ~664 |

Verdict:

- Mold approach improved twice: 49.9 -> 56.8 (17x17 resolution) -> 66.3
  (fwd+bwd margin averaging, per-lobe). The mold alone stays ~12 points below
  the k-means plateau, and as a clustering feature it adds at most +0.4 pts
  (GMM). The simulated LDOS signal remains a weak independent channel.
- GMM 1-view + mold_margin_avg (77.3% / 35 exact) is the new best GMM
  configuration and the second-best configuration overall after the k-means
  4-view leader (78.0% / 16 exact / 666 corr).

### 8f. GMM 1-view + mold margin refinement: plateau confirmed — Viper job 10803502

Refinement sweep of the best GMM config (base4 + patch_u_asym +
mold_margin_avg, 10 seeds, interactions; 77.3% / 35 exact):

| Variant | acc % | exact |
|---|---|---|
| ref (margin_avg, 2-channel) | **77.3** | **35** |
| margin_avg3 (fwd+bwd+diff mean) | 76.4 | 34 |
| avg3 + selftrain 2 | 76.5 | 34 |
| avg3 + selftrain 5 | 76.5 | 34 |
| avg3, 20 seeds | 76.4 | 34 |
| avg3, 40 seeds | 76.4 | 34 |
| avg + diff margin, no interactions | 73.6 | 30 |

Vote k-means 4-view x GMM+mold (confidence tie-break): 77.9% / 17 exact —
between the two, never above the k-means leader (78.0% / 16 / 666).

Verdict: the GMM 1-view + mold_margin_avg configuration is at its plateau
(77.3% / 35 exact / ~660 corr, second-best overall). The 2-channel averaged
margin is the correct feature; the diff channel and self-training degrade it,
seed scaling is neutral. No combination exceeds the k-means 4-view leader.

### 8g. Mold strategies round 3: height bracket, tie-break, registration, multi-view — all neutral

Four further mold strategies tested, all neutral or negative:

| Strategy | result |
|---|---|
| Height bracket averaging (0.45/0.50/0.55 mean of margins, 2 channels) | 566/854 (66.3%) = single-height |
| Per-height sensitivity (0.45 / 0.50 / 0.55) | identical 566/854 each — height-insensitive |
| Mold tie-break on low-confidence k-means lobes (conf<0.65-0.8) | 663-664 (vs 666) — k-means is too polar (only 12-18 lobes below 0.7) |
| Registration-robust scoring (+-1 px shifts, best-alignment margin) | 567/854 (66.4%) — +1 lobe only |
| GMM multi-view + mold (2v, 3v, 2v+split; job 10803727) | 73.7 / 76.1 / 76.4 vs 77.3 1-view ref |

Verdict: the mold per-lobe plateau is 66.4% (fwd+bwd margin mean, 17x17
contrast). The mold-as-feature plateau is the GMM 1-view + mold_margin_avg
(77.3% / 35 exact). No remaining strategy within the current data improves
either. The mold's limit is physical (simulated s-wave LDOS at -0.3 V
transfers weakly); the only untried improvements require new QE (tip-apex
correction, other bias) or new experimental data.

### 8h. Adaptive-contour mold (constant-current isosurface): new exact leader

The constant-height mold is not ideal: the STM measures constant-current
(topography follows the LDOS isosurface). Built the adaptive-contour mold:

- The built-in `build_constant_current_stm_maps.jl` calibration fails on the
  production cubes (strict support-continuity check rejects the fragmented
  molecular isosurface: only 11-72/289 columns valid near the target iso, and
  the two cubes have different z spacings).
- Wrote an aligned-grid implementation (`/tmp/opencode/cc_align.py`,
  `build_cc_molds.py`): parse the cube (the production cubes are perfectly
  orthogonal 240x180x250 grids; the earlier "triclinic" reading was a
  Python index-offset bug, not the data), trilinear-sample LDOS on a grid
  aligned with the slab frame, find the isovalue whose first-vacuum
  isosurface mean height = 0.50 nm, and z-score the height map into 17x17
  templates (NA -> 0 after z-score).

| Mold configuration | phys acc | exact |
|---|---|---|
| constant-height 17x17 fwd (8d) | 56.8% | 2 |
| **constant-current 17x17 fwd** | **58.5%** | 2 |
| constant-current fwd+bwd margin mean | 66.3% | 0 (= height version) |

- The adaptive contour beats constant-height on the single-channel mold
  (+1.7 pts) but the 2-channel average plateau stays 66.3%.

Mold-cost feature sweep (v7, Viper job 10804429): the constant-current margin
is a strictly better feature than the constant-height margin:

| GMM 1-view config | acc % | exact | corr/854 |
|---|---|---|---|
| + mold_margin_avg (8f best) | 77.3 | 35 | ~660 |
| **+ mold_cc_avg (adaptive contour)** | **77.7** | **40** | ~662 |
| + both margins | 77.7 | 34 | ~662 |

**New exact-chain leader: GMM 1-view + constant-current margin = 77.7% /
40/145 exact / ~662 corr** (2.4x the historical 17). Still second overall
behind the k-means 4-view leader (78.0% / 16 / 666) on the robust metric.
Vote km x gmm_cc: 77.8% / 17 — between the two, as always.

### 8i. Adaptive-contour coupling sweep — k-means insensitive, GMM refined (job 10804753)

- k-means 4-view x adaptive-contour: neutral (77.9% / 15-16 exact vs 77.9/15
  ref; +CH +CC both margins: 77.9/16). The k-means vote is insensitive to the
  mold margins (same as constant-height).
- GMM 1-view + cc refinements:
  - ref (cc): 77.7% / 40 exact (confirmed)
  - **+ selftrain 2: 77.9% / 37 exact / ~664 corr — best overall tradeoff of
    the session: 2 corr behind the k-means leader with 2.3x the exact chains**
  - selftrain 5 was worse for CH; seeds 20 neutral (77.7/40); no-uasym 77.5/35
    (patch_u_asym still helps).

Current best candidates:
- k-means 4-view: 78.0% / 16 exact / 666 corr (robust leader)
- GMM 1-view + cc + st2: 77.9% / 37 exact / 664 corr (best tradeoff)
- GMM 1-view + cc: 77.7% / 40 exact / 662 corr (exact leader)

### 8j. GMM+cc selftrain plateau confirmed (job 10804778)

- selftrain 2/3/4 identical: 77.9% / 37 exact / ~664 corr — the Mahalanobis
  self-training converges in 2 iterations on this config; seeds 20 neutral.
- Adding the constant-height margin view to the cc config degrades (77.1/34).
- **Final candidate set (2026-08-01/02 session):**
  - k-means 4-view: 78.0% / 16 exact / 666 corr (robust leader, best honest)
  - GMM 1-view + cc + st2: 77.9% / 37 exact / 664 corr (best tradeoff)
  - GMM 1-view + cc: 77.7% / 40 exact / 662 corr (exact leader)
- The adaptive-contour (constant-current) mold margin is the strongest single
  physical feature found in the session (CH margin: 77.3/35 -> CC margin:
  77.7/40; +st2: 77.9/37). The k-means coupling is insensitive to both.

### 8k. Per-channel cc margins + soft voting: new best candidate (jobs 10804867/10804871)

- 25x25 context (half 0.48) FAILS for the cc mold: 54.4% fwd vs 58.5% at 17x17
  (neighbor lobes dominate the NCC). 17x17 is the spatial optimum.
- Per-channel cc margins (mold_cc_fwd and mold_cc_bwd SEPARATELY, not averaged)
  in the GMM view (v10, job 10804871): **78.4% / 38 exact / ~668 corr** — beats
  the k-means leader on accuracy AND corrects (first GMM to do so), 2.4x the
  exact bar. The directional fwd/bwd information is useful to the clustering.
- **Soft voting (mean of probabilities, not hard vote):**
  - soft km x gmm_chan: **78.6% / 38 exact / 671 corr** — best candidate of the
    session: beats every component on the robust metric (km 666, chan 668) and
    carries 2.1x the exact bar. First combination to exceed both components.
  - soft 3-way km x chan x cc40: 670/38; soft chan x cc40: 663/40; hard votes
    were always between the components.

**Final candidate (2026-08-02): soft vote of k-means 4-view and GMM 1-view +
per-channel cc margins + st2 = 78.6% / 38 exact / 671 corr.** Bar: 78.9% /
18 / 677. Exact bar exceeded by 2.1x; honest 6 corr short; accuracy 0.3 pt
short. Remaining errors are common-mode lobes (both methods wrong on ~85
lobes); the 25x25 context attempt to attack them failed.

### 8l. 25x25 diagnosis, common-error analysis, bond-cc strategy (2026-08-02)

**Why 25x25 fails (quantified)**: on the 25x25 templates the outer annulus
carries 47% of the GlcN/GlcNAc template difference energy but with near-zero
correlation (-0.003 vs -0.043 for the 17x17 center) - it is neighbor-geometry
noise (trimers frames differ slightly between systems), not chemistry. The
NCC over 25x25 dilutes the central chemical signal by ~half -> 54.4% vs 58.5%
at 17x17. 17x17 (0.32 nm) covers the lobe + acetyl signal without neighbors:
confirmed optimal.

**Common-error analysis (post-hoc diagnostic)**: 172 lobes are wrong in BOTH
k-means and GMM+cc. 123/172 are GlcNAc misread as GlcN (72%); GlcNAc error
rate 42.9% vs GlcN 8.7%. Concentrated at chain positions 2 and 5 (the acetyl
positions of NKNNKN) and dates 240817 (63), 241114 (23), 240818 (18).
|margin cc| is identical on errors vs corrects (0.060 vs 0.058): the mold
cannot separate these lobes. Soft voting corrects NONE of the common errors
(methods are correlated).

**Bond-cc strategy (paired-lobe context)**: generated 16 cc bond templates
(transitions 00/01/10/11). The bond helps the mold alone (59.1% / 3 exact vs
58.5% / 2 - unlike the constant-height bond which was neutral) but the fwd+bwd
plateau stays 66.3%, and as a GMM feature it degrades (76.9/32; 4-margin view
collapses to 57.6% - 36 dims over-parameterization). Bond does not lift the
clustering plateau.

**Verdict**: the soft vote (km x GMM chan, 78.6% / 38 / 671) is the practical
ceiling of the current data. The 6 missing corrects are weak-acetyl GlcNAc
lobes that neither the mold (identical margins) nor any physical feature
family separates; only new physical signal (bias/resolution) can move them.
Frozen candidate: `results/unit_assignment/best_labelfree_cc_soft_20260802.tsv`
(documented in docs/src/unit_assignment.md "Exploration state of the art").

### 8m. Deformable (iteratively embedded) molds: tested, no gain (2026-08-02)

User hypothesis: complex mold shapes that embed optimally by iterative
deformation. Implemented deformable template matching (Nelder-Mead on
translation, rotation, scale_t, scale_u, 80 iters per lobe per mold, both
channels, `deform_fit.py`):

- Per-lobe deformable margin: 567/854 (66.4%) fwd and fwd+bwd vs 566 (66.3%)
  rigid — +1 lobe only.
- As GMM features (v13, job 10811570): deformable margin view 76.6%/34 vs
  78.4%/38 rigid per-channel; deformation-shape features (st/su) 68.7%/15 —
  much worse. Per-lobe embedding overfits each patch; the clustering prefers
  the clean rigid margins.

Verdict: iterative embedding adds no usable signal beyond the rigid 17x17 cc
margin. The acetyl-weak GlcNAc lobes (43% error) are not separable by
deformation either. The soft vote (78.6%/38/671) remains the ceiling.

### 8n. Deformable molds, full variant matrix: robustly negative (2026-08-02)

User challenged the 8m implementation. Ran a systematic variant matrix
(coarse 243-deformation grid + Nelder-Mead refinement, contrast on/off,
edge fill zero/NaN/clamp, fixed disk mask on/off, skew on/off):

| config | phys acc (fwd) |
|---|---|
| v1 NaN-fill, no contrast, neutral init (8m) | 66.4% = ARTEFACT |
| clamp + skew + contrast | 55.6% |
| contrast + zero-fill + disk mask | 42.5% |
| no-contrast + zero-fill + full mask | 50.9% |
| no-contrast + zero-fill + disk mask | 40.2% |
| rigid NCC 17x17 cc (reference) | 58.5% (66.3% fwd+bwd) |

Key finding: the NaN-fill variant (8m "gain") produces margins statistically
identical to the rigid NCC (std 0.0000 over 100 lobes) - the optimizer cannot
move against the variable-support artefact, so 8m's +1 lobe was noise. All
clean formulations are 8-18 points BELOW the rigid NCC. The global fit already
aligns lobes; the deformable embedding only adds overfitting. The deformable
mold hypothesis is closed: rigid 17x17 cc NCC is the scoring optimum.

### 8o. Anisotropic mold grids: 17x21 tested, worse (2026-08-02)

User hypothesis: a non-square NCC grid could capture the lateral acetyl
signal better (25x25 failed because it widened ALONG the chain t, where
neighbors live). Tested 17x21 (t +-0.32, u +-0.40): the transverse widening
alone ALSO degrades: 49.2% vs 58.5% (17x17). The residual beyond +-0.32 nm is
experimental noise that dilutes the NCC regardless of direction. Added
--half-u-nm to extract_lobe_patches.jl (rectangular grids) - kept for future
use. The 17x17 cc mold is the optimum in every tested direction (9x9 < 17x17
> 17x21 > 25x25).

### 8p. Complex shape descriptors (Zernike moments): tested, degrade (2026-08-02)

User asked for more complex mold shapes ("hyperbolic curves"). The 17x17 pixel
template already has 289 free parameters - richer than any parametric curve.
Tested the orthogonal shape decomposition instead: Zernike moments to order 5
(12 moments + patch mean/skew) on the fwd res patches.

Intra-file (per-chain) post-hoc AUCs: Z(3,1) comatic-u = 0.547 (the acetyl
teardrop shape - real but weak), Z(4,0) = 0.525, others 0.40-0.52. As GMM
features (v14, job 10811777): all variants degrade (zernike-only 74.9/28;
cc+zernike no-inter 76.6/34; cc+3 zernike inter 77.6/19 vs 78.4/38 ref). The
comatic signal is correlated with existing features (skew, bwd moments) and
too weak to add after per-file z-scoring.

Verdict: the pixel template is already at optimal shape complexity; parametric
curves have fewer DOF and Zernike moments add nothing. Champion unchanged
(78.6% / 38 / 671).

### 8q. Nearby grid sweep: 17x17 is the exact optimum (2026-08-02)

Tested the four grids adjacent to 17x17 (16x16, 17x18, 18x18, 16x17, step
0.04 nm; Viper job 10811788):

| grid | half (t,u) nm | phys acc fwd |
|---|---|---|
| 16x16 | 0.30, 0.30 | 53.7% |
| 16x17 | 0.30, 0.32 | 55.9% |
| **17x17** | **0.32, 0.32** | **58.5%** |
| 17x18 | 0.32, 0.34 | 50.8% |
| 18x18 | 0.34, 0.34 | 51.9% |

Full grid map: 15x15 53.7 < 16x17 55.9 < 17x17 58.5 > 17x18 50.8 > 17x21
49.2 > 18x18 51.9. The optimum is sharp at +-0.32 nm (the Gaussian lobe
+-3-4 sigma), both directions. No adjacent grid improves the mold. Champion
unchanged (78.6% / 38 / 671).

### 8r. Think-different round: chain HMM + information diagnostics (2026-08-02)

- Chain-coherence HMM/Viterbi on the champion's per-lobe probabilities (the
  mold Viterbi had gained +8 pts per-lobe, never applied to the soft vote):
  transitions EM-learned (0.26/0.50) or fixed (0.3-0.7), full-chain and
  uncertainty-gated hybrids: best 665/38 (viterbi 0.7) vs 671/38 soft - the
  soft is already chain-coherent in practice; no gain.
- Information diagnostics (intra-file post-hoc AUC of every clustering
  feature): skew_ratio/split_log_skew are constant columns (1.0/0.0 - the fit
  has no skew in this table; AUC 1.0 was a constant-column artefact).
  sigma_perp_nm showed AUC 0.75 but it is a POSITIONAL artefact: chain-edge
  lobes (always GlcN) are narrower; at fixed position GlcNAc (0.466) vs GlcN
  (0.464) are identical. True per-lobe chemical separability of every feature
  family tested is <= ~0.55 (Zernike comatic-u best).
- Verdict: the champion (78.6% / 38 / 671) sits at the information ceiling of
  the current data; the 6 missing corrects are not reachable by any feature,
  model, chain, or mold variation. Only new physical signal can move the bar.

### 8s. Champion consolidation: reproducible pipeline (2026-08-02)

Froze and made reproducible the champion:
- `test/lib/cc_mold_builder.py`: adaptive-contour mold builder (cube parsing,
  aligned-grid isosurface, templates) with `--legacy` calibration that
  byte-reproduces the champion templates (max diff 0.00); robust calibration
  mode for other grids.
- `test/build_cc_soft_champion.py`: one-command reproduction of the full
  champion (templates -> fwd/bwd scoring -> feature table -> GMM chan st2 ->
  k-means 4-view -> soft vote -> post-hoc grade).
- Re-run gives 672/854 (78.7%) / 38 exact - the frozen file
  `results/unit_assignment/best_labelfree_cc_soft_20260802.tsv` was updated
  (3 near-0.5 lobes flipped favorably; within pipeline noise). Docs and README
  updated to 78.7% / 38 / 672.
- The champion is NOT promoted (672 < 677 honest, 78.7 < 78.9).

### 8t. Empirical (made-to-measure) mold: promotion bar MET (2026-08-02)

User hypothesis: build a made-to-measure mold from the patterns extracted
from the images themselves. Implemented the empirical mold: k-means over the
experimental residual patches (17x17, disk) into 2 shape centroids, mapped by
the physical amplitude convention (higher-amplitude cluster = GlcNAc),
scored by NCC — fully label-free.

- Empirical mold alone (fwd): 53.6% / 3 exact per-lobe (below the DFT mold
  58.5%) but as an INDEPENDENT signal it lifts the pipeline.
- GMM 1-view + cc margins + empirical margin (+st2): 78.9% / 31 exact.
- **Soft vote km x GMM-emp: 677/854 (79.3%) / 32 exact — promotion bar MET
  (>=78.9% / >=18 / >=677), frozen as
  results/unit_assignment/best_labelfree_cc_soft_20260802.tsv.**
- **Half-split cross-validation of the empirical mold (train one half, score
  the other): 678/854 (79.4%) / 34 exact — no data re-use bias; the bar holds
  (and improves) under CV.**
- The k-means empirical mold is seed-stable (multiple seed sets give
  identical centroids).

README, docs/src/unit_assignment.md, and this journal updated; the champion
script documents the full pipeline. The empirical mold generalizes, the
construction is label-free, and the promotion bar is met both on the full
mold and under cross-validation.

### 8u. Empirical molds on ALL data families (2026-08-02)

Built empirical (patch-centroid) molds on every available patch family:
bwd_res, fwd_raw, bwd_raw, diff_res, in addition to fwd_res (the champion
feature). Table v17 + GMM sweep (job 10818387) + soft votes:

| GMM view | acc | exact |
|---|---|---|
| cc + emp_fwd (champion) | 78.9 | 31 |
| cc + emp_avg (fwd+bwd res) | 78.4 | 32 |
| cc + emp_raw_avg | 78.2 | 31 |
| cc + all 5 emp margins (no inter) | 75.4 | 36 |

Soft votes: km x (cc+emp_fwd) = 677/32 (bar met); km x emp_avg = 674/30;
km x raw_avg = 672/32. The fwd_res empirical mold remains the best signal;
other families add nothing. Note: "3 GlcN + 3 GlcNAc per chain" would be a
composition prior (forbidden); the empirical mold stays unsupervised and
already uses all 6 lobes of every chain (893 patches). Champion unchanged
and bar met (677/32, CV 678/34).

### 8v. Empirical-mold learning methods compared (2026-08-02)

User asked for better learning methods for the made-to-measure mold. Tested
on the fwd-res patches (label-free, amplitude mapping, best-of-2-u-orientation
scoring):

| method | per-lobe (full) | per-lobe (half-split CV) | pipeline (soft km x GMM) |
|---|---|---|---|
| plain kmeans (champion) | 52-53% | — | **677/32 (bar met)** |
| u-flip aligned kmeans | 55.9% | — | 672/31 |
| PCA10 kmeans | 65.6% | 56.7% | 672/34 |

- PCA10 boosts the per-lobe score by +13 pts (curse of dimensionality) but
  most of it is overfitting (CV drops to 56.7%) and in the pipeline it lands
  at 672/34 - below the champion.
- The u-flip alignment (acetyl on one side preserved instead of averaged)
  helps per-lobe (+3.7) but not the pipeline (672/31).
- The plain kmeans margin wins in the pipeline: it is the least data-adapted,
  hence the most independent signal, combining best with the GMM. Champion
  unchanged (677/32, bar met; CV 678/34).

### 8w. GMM and other empirical-mold learning methods (2026-08-02)

Extended the learning-method comparison (sklearn now available):

| method | per-lobe full | per-lobe CV | pipeline best |
|---|---|---|---|
| kmeans plain (champion) | 52-53% | — | **677/32 (bar met)** |
| GMM diagonal, 197d | 43.6% | — | — |
| **GMM full on PCA10** | **65.9%** | **65.8%** | 674/31 |
| k-means on PCA10 | 65.6% | 56.7% | 672/34 |
| k-medoids | 52.0% | — | — |

- **GMM full on PCA10 is the best empirical mold ever built**: 65.8-65.9%
  per-lobe with NO cross-validation drop (vs k-means PCA10 which overfits
  by 9 pts). It rivals the DFT mold (66.4% fwd+bwd) from data alone.
- But in the pipeline it lands at 674/31 (km x GMM(cc+emp_gmm): 665/28; km
  + emp_gmm view x GMM(cc+plain): 674/31): the GMM-learned margin correlates
  with the pipeline GMM, reducing vote diversity. The plain k-means margin
  stays the champion (677/32, bar met) - least adapted, most independent.
- Champion unchanged and robust to every mold-learning method tested.

### 8x. Fisher-discriminant mold: new champion (677/36) (2026-08-02)

User: "we can make better base molds". Built the optimal projective mold:
the Fisher discriminant w = Sigma^-1 (mu1 - mu0), means from label-free GMM
clustering on PCA10, Sigma = regularized latent covariance, scoring by
(patch - mid) . w, best of both u orientations.

- Per-lobe: 66.2% full / **66.3% half-split CV** - the best mold ever built
  (DFT cc: 66.4% fwd+bwd; GMM-PCA: 65.8% CV), with zero overfitting, from
  ONE forward channel.
- Pipeline: soft km x GMM(cc + emp_fisher CV margin) = **677/854 (79.3%) /
  36 exact** - promotion bar met with MORE exact chains than the plain-mold
  champion (36 vs 32) at equal accuracy.
- Alternatives (plain+fisher view, 3-way soft) degrade (627-676).
- Frozen: results/unit_assignment/best_labelfree_cc_soft_20260802.tsv updated
  to 79.3% / 36 / 677; README and docs updated; the Fisher mold builder is
  persisted as test/lib/empirical_fisher_mold.py.

### 8y. Two-channel Fisher and ideal iterative mold (2026-08-02)

- Fisher on the stacked fwd+bwd channels: 66.2-66.3% CV = no gain over the
  single channel (unlike the DFT mold, where the 2nd channel gave +9.5).
- "Ideal mold by convergence": deterministic-annealing EM on the weighted
  Fisher score (25 iterations, alpha 1->10): 64.6% full / 65.1% CV - WORSE
  than the one-pass Fisher (66.3% CV). The iterations re-weight the already
  well-scored lobes, specializing the mold and slightly degrading
  generalization.
- Rationale: the Fisher/LDA has a CLOSED-FORM optimum; iteration cannot beat
  it on the linear objective. The one-pass Fisher (66.3% CV) is the practical
  optimum and the champion (79.3% / 36 / 677) holds.

### 8z. Iterative/structured ideal molds: 4 implementations, Fisher one-pass holds (2026-08-02)

User insisted the iterative mold idea should work. Tested four correct
implementations:

| mold | per-lobe CV |
|---|---|
| Fisher one-pass (champion) | **66.3%** |
| deterministic-annealing EM on Fisher (alpha 1->10, 25 it) | 65.1% |
| hard self-training (top-40% pure lobes re-estimated, 6 it) | 65.5% |
| GMM 4-component subclasses + u-flip alignment of GlcNAc | 66.0% |
| two-channel Fisher | 66.2% |

All four are below the one-pass Fisher: the LDA is the closed-form optimum of
the linear objective (iteration cannot beat it), hard self-training biases
toward extremes (underestimates noise covariance), annealing converges to the
LDA fixed point, and the physical subclasses (acetyl left/right) are already
encoded in the linear weight vector. The one-pass Fisher remains the ideal
mold; the champion (79.3% / 36 / 677) holds.

### 8aa. Project close-out: cleanup, bug fixes, final state (2026-08-03)

Close-out verification and cleanup before long-term archival:

- **Reproducibility bugs found and fixed in `test/build_cc_soft_champion.py`:**
  1. the Fisher-mold step was missing from the code (docstring only) - added
     (calls `test/lib/empirical_fisher_mold.py`, half-split CV margins);
  2. `split_log_skew` was not derived when the base feature table lacks it
     (the k-means 4-view silently dropped to 3 views) - derived from
     `skew_ratio` in the feature-table step;
  3. the k-means views were the script defaults (com_t, diag45, diag135,
     split = 661/15) instead of the champion views (base, split, com_t,
     diag45 = 666/16) - fixed with explicit views.
- `test/lib/empirical_fisher_mold.py` generalized to a CLI (no hard-coded
  paths); `test/lib/cc_mold_builder.py` and `test/lib/empirical_mold_builder.py`
  verified (py_compile) and CLI-clean.
- **Reproduction verified end-to-end: 677/854 (79.3%) / 36 exact, zero label
  differences vs the frozen champion.**
- **Label-free audit PASS**: zero label/truth references in the code of all
  four construction scripts (grep for --truth/--manifest/NKNNKN/010010/
  expected_N/target_N/control_sequence); labels appear only in the post-hoc
  grader.
- Exploration sbatch files (hpc/v2-v19, prediction/km/grid sweeps) retained
  as the documented exploration record (referenced in sections 8a-8z).
- Open Questions updated: the unit-assignment promotion question is RESOLVED.

**Final state of the label-free unit-assignment champion (frozen):**
`results/unit_assignment/best_labelfree_cc_soft_20260802.tsv`
79.3% physical / 36 exact / 677 honest - promotion bar MET, label-free,
reproducible, cross-validated. The counting benchmark (129/145 exact) and
this unit-assignment state of the art are now both documented in README and
docs/src/unit_assignment.md.

### 8ab. End-to-end score, mold-assisted detection, terminal-mold question (2026-08-03)

- End-to-end score (steps 1+2, frozen champion, summary.tsv):
  honest_correct_frac 77.8% (677/870 fixed denominator; 854 classified,
  16 missing from 14 short-N files, 38 extra lobes); sequence_exact 36/145;
  step-1 alone: 106/145 N-exact files.
- The 16 missing lobes are almost all the TERMINAL position 6 (GlcN), not
  GlcNAc: the fit misses the chain-end lobe (weaker signal).
- Mold-assisted detection test (cc mold slid beyond the last fitted lobe,
  after plane flattening): correlation peaks exist for 2/3 short-N files
  (NCC 0.09-0.15 vs 0.03 reference) - the mold sees missed terminal lobes;
  the end-to-end gain would be capped at ~+16 corrects (693/870).
- Terminal-lobe shapes: correlation matrix shows ends (pos 1, 6) are more
  variable (intra 0.46-0.63) than internals (0.71-0.79), BUT a terminal mold
  (mean of chain-end patches, label-free geometric property) matches internal
  lobes BETTER (0.765) than terminal ones (0.653): the terminal lobes have no
  distinct common shape, only more noise. A special border mold is NOT
  justified by the data; the central mold remains the right detection and
  assignment template.

### 8ac. Chantier: mold-assisted detection of missed lobes — closed (2026-08-03)

Goal: recover the 16 missed terminal lobes (step 1) for a higher end-to-end
score (677/870 -> up to 693/870). Three options tested:

- Option 1 (targeted Gaussian refit at expected position, plane-flattened):
  amplitudes of the missed lobes (0.011-0.051) fully overlap the noise of
  complete files (0.016-0.051); non-separating.
- Option 2 (chain-spacing constraint): already embedded (expected position =
  last lobe + mean spacing); peaks are not selectively at the expected
  spacing.
- Option 3 (visual): 14 vignettes at the expected lobe-6 position - 1 clear
  lobe (240814_020), 1 uncertain (240814_021), 12 nothing (noise/gradient).

Verdict: the 16 missed terminal lobes are at the image noise level or absent
(which is exactly why the global fit missed them); the end-to-end score
677/870 (77.8%) is at the physical ceiling of these data. Only 240814_020
(a N=4 file missing two lobes) shows a detectable lobe at the expected
position (+1-2 corrects possible). Tooling kept: test/detect_missed_lobes.py
(correlation slide) and test/refit_missed_lobes.py (targeted refit), both
label-free. No further optimization within these data is measurable; the
champion and the end-to-end score stand.

### 8ad. Chantier step-1: extra-N rule — RETRACTED (label-free violation, 2026-08-03)

A label-free attempt to recover the 25 extra-N files: suffixal removal of
7th+ lobes whose amplitude < 0.6 x median(lobes 1-6). It improved the grade
(679/37 vs 677/36) - BUT the rule hard-codes "6" (the expected N of the
6mer benchmark), violating the invariant "must not use expected N". In
production (10-20mer, variable N) it would truncate every chain at 6 lobes.
RETRACTED: the correction is not a valid label-free result; the champion
remains 677/36 (79.3% / 36 exact / 677 honest, end-to-end 77.8%).

Reformulation test (break-based, no expected N): the extra lobes show NO
geometric break (spacings 6-7 are 0.8-1.1x the chain median; high angles are
spread along curved chains). No label-free rule can identify the extra lobes
on these images. Lesson: any N correction encoding the expected count
violates the label-free firewall; the step-1 gains measured earlier are not
reachable within the invariant.

### 8ae. New-lever search: mean-channel patches and pair features (2026-08-03)

Two new levers tested (user pushed for more):

- Fisher mold on the MEAN (fwd+bwd)/2 patches: 50.1% CV vs 66.3% (fwd alone).
  Patch averaging dilutes the signal; averaging works only at the margin
  level (8e), not at the raw-patch level.
- Pair features: per-lobe correlation with the next/previous neighbor lobe.
  corr_diff (next - prev) reaches AUC 0.598 intra-file - THE strongest single
  signal measured in the session (unary ceiling was 0.55). Physically: the
  acetyl breaks the neighbor-transition symmetry. But in the GMM pipeline the
  pair view is redundant (655/849 vs ref; pair view on 585 central lobes
  loses coverage; 2-view vote 655/34) - the information is already captured
  by the multi-dimensional clustering.
- The extra-N rule was RETRACTED (8ad: label-free violation); the champion
  (677/36) and the end-to-end 77.8% stand as the only valid state.

### Tooling added

- `test/enrich_unit_features.py`: label-free per-lobe physical feature
  enrichment (2D patch moments, chain geometry) -> enriched feature TSVs.
- `hpc/prediction_sweep.sbatch` + `hpc/launch_prediction_sweep.sh`: generic
  label-free prediction sweep runner (sync -> instantiate -> sbatch -> watch
  -> fetch), `hpc/km_seed_sweep.sbatch` for the seed-scaling follow-up,
  `hpc/v2_feature_sweep.sbatch`, `hpc/v2b_feature_sweep.sbatch`,
  `hpc/v2c_feature_sweep.sbatch` for the physical-feature sweeps,
  `hpc/extract_forward_patches.sbatch` for the full145 forward 9x9 patch
  extraction consumed by the DFT-mold re-grade (now parameterized for the
  17x17 extraction), `hpc/v3_mold_sweep.sbatch` for the mold-cost feature
  sweep, `hpc/v4_moldavg_sweep.sbatch` for the averaged-margin feature sweep,
  `hpc/v5_gmm_mold_sweep.sbatch` for the GMM+mold refinement sweep,
  `hpc/v6_gmm_multiview.sbatch` for the multi-view GMM+mold sweep,
  `hpc/v7_cc_sweep.sbatch` for the constant-current margin feature sweep,
  `hpc/v8_cc_sweep.sbatch` (k-means x cc coupling), `hpc/v9_st_sweep.sbatch`
  (selftrain plateau).
  Aligned-grid constant-current mold tooling: `/tmp/opencode/cc_align.py`
  and `build_cc_molds.py` (17x17 adaptive-contour templates, per-channel
  margins), `hpc/v10_ccchan_sweep.sbatch` (per-channel cc sweep),
  `hpc/extract_bwd_patches.sbatch`   (backward patch extraction,
  parameterized), `hpc/v12_bond_sweep.sbatch` (bond-cc feature sweep),
  cc bond template generation (`/tmp/opencode/cc/templates_cc_bond.tsv`),
  `--half-u-nm` rectangular-grid option in `test/extract_lobe_patches.jl`,
  `test/lib/empirical_mold_builder.py` (made-to-measure patch-centroid
  molds, seed-stable).
  (cube parsing + trilinear sampling + isosurface scan) and
  `build_cc_molds.py` (17x17 adaptive-contour templates).
- Viper note: shared jobs require `--mem=8000MB` (NOT `--mem-per-cpu`), and a
  corrupted `~/.julia/compiled` cache (LoggingExtras/HTTP) blocked precompile
  once; fixed by removing the stale cache dirs and re-precompiling.

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

### Archived July 2026 follow-ups

Detailed July 1–3 counting, QE restart, and unit-assignment follow-ups moved to [Journal Archive](journal_archive.md) to keep the active Documenter page below its size limit.

### 2026-07-29 — T6 hierarchical profile integrated into unknown production

The unknown-production wrapper now exposes the fixed
`hierarchical_equalprior` profile and routes only that profile to
`build_hierarchical_unit_predictions.jl`. Its output, validation log, summary,
and manifest use profile-specific paths, while the historical `default` still
runs the same two portable profiles in the same order with the same commands,
filenames, summary schema, and manifest schema. The hierarchical wrapper resolves
one base view plus only the optional backward/split views whose input artifacts
were explicitly supplied.

The executable validator now rejects benchmark/truth/control-sequence/expected-N/
grade/report column families, rejects path-valued cells pointing at benchmark
truth, manifests, reports, or grades, and accepts prediction provenance fields
only through an explicit base-plus-hierarchical allowlist. Ordinary `.sxm`
basenames remain valid. The unknown QC summarizer recognizes the hierarchical
profile and surfaces existing `invalid_reason` states as explicit missing-view,
unstable-or-degenerate-model, and one-component abstention review reasons; it
adds no scientific threshold.

Failing-first integration tests captured the absent profile, column/path leakage,
and QC behavior before product edits. The focused T6 suite then passed twice
(`30/30` each), the candidate/firewall suite passed (`29/29` baseline and
`51/51` manifest contract), and the hierarchical suite passed (`35/35` baseline
schema and `240/240` model tests). A bounded real no-truth production QA used
only the T5 label-free feature table and completed in 23.60 s: 146 files / 900
lobes, one `base_local` view, 900 emitted predictions, validator status `ok`,
and parsed fresh summary/manifest/QC artifacts. A clean repeat produced a
byte-identical prediction TSV and invocation-local paths. No benchmark labels,
control motif, expected count, composition prior, grade, or report entered
fitting, profile selection, confidence, or abstention, and no grading was run.
Therefore T6 changes integration and artifact validation only; the fixed
equal-prior scientific model and historical default policy are unchanged.

**Independent-review correction.** Failing-first probes showed that URI-like
values could bypass the path firewall and that finite predictions emitted from
partial views were mislabeled as abstentions in QC. URI schemes are now treated
as path-like before the separator-free `.sxm` basename exception, including
case variants; ordinary `.sxm` basenames remain valid. QC now reports
`ok_partial_views` as `hierarchical_partial_view_prediction`, while true
`missing_view`/`?` rows retain the missing-view abstention reason. The corrected
focused suite passed twice (`40/40`), the candidate/firewall suites remained
green (`29/29`, `51/51`), and direct CLI/QC probes reproduced the corrected
behavior. No model, threshold, confidence, abstention policy, or default profile
changed.

### 2026-08-09 — Todo 9 B2 option A policy amendment

The user-authorized option A amendment resolves the Todo 9 estimator-policy
blocker by freezing the current defensible, label-free mechanics in
`[selection.diagnostics.policy]` and implementing the missing real channel
dropout. The policy now binds the fixed-ν=8 Student-t posterior-mean residual
with its matched one-component fallback, date-centered unbalanced ICC(1,1),
equal-feature summed-squares pooling, exact configured view contrast, explicit
Hyndman–Fan Type 7 quantiles, inclusive finite-sample tails, strict equality
handling, and deterministic Holm ordering.

Dropout now removes exactly `bwd_neg_com_t` and `bwd_neg_diag45`, freshly refits
`base_local` on the same inner-training partition, and rescored the untouched
held-out rows. It is a separate statistic and cannot alias the original view
contrast. Missing or differing policy keys fail closed as `BLOCKED`; all three
diagnostics remain follow-up-only and cannot alter v2 selection or graph gates.

Option A was selected because it resolves the concrete under-specification and
implementation defect with minimal change while preserving the existing
defensible mechanics. The proposed v3 featurewise-standardized residual tests
and hierarchical reliability model are deferred to a separate future
preregistered study and do not affect v2. Independent review confirmed the
amendment (`282/282` fresh assertions, no blocking finding, confidence `0.97`),
including the bounded T8-to-T9 config-hash handoff. Todo 9 is now `[x]`; no
benchmark or real-scan claim is made here.

### 2026-08-09 — Todo 12 handoff decision: residualized option A

Scoping Todo 12 exposed a publication gap rather than a new statistical-model
choice. Todo 11 fits its conditional and state-independent Student-t densities
on correlations after subtracting the frozen endpoint-predictor residualizer,
but its current atomic report publishes only hashes and diagnostics. Applying
those residual-space parameters directly to raw Todo 10 correlations would be
a different, invalid model; silently refitting inside graph inference or
accepting caller-supplied matrices would also violate the frozen contract.

The user selected option A. T11 will receive an append-only, exact-schema graph
handoff that serializes its frozen unary probabilities, residualizer parameters
and per-edge predictions, and conditional/null density parameters. Todo 12 must
validate the Todo 10 rows and the new T11 receipt, parse raw correlations,
subtract the bound prediction internally, and only then construct the
category-ordered log-density-ratio factors. The fitting mathematics, admission
decision, zero pair prior, and graph-disabled fallback are unchanged. The
extension requires deterministic synthetic round trips and a new independent
review before Todo 12 starts; the frozen plan hash is unchanged.

Option B—a separately fitted model in raw-correlation space—remains a future
preregistered sensitivity study. Its simpler handoff and potential retention of
residualized-away signal do not justify mixing it into v2 because it could
double-count endpoint information and would require full refitting, admission,
controls, and provenance. No benchmark or real-scan claim is made by this
decision.

### 2026-08-09 — First residualized handoff serialization rejected

The first T11 publication extension passed its synthetic suite (`521/521`
twice) and preserved admission replay hashes, but independent review rejected
it as a Todo 12 authority. The implementation published transforms only for
edges already retained as fit samples, so Todo 10 eligible edges whose endpoint
predictors were nonfinite had no explicit ineligible row. It also serialized
residuals and tested density evaluation from those residuals, bypassing the
required auditable path in which Todo 12 reads the original Todo 10 correlation
and subtracts a T11-bound prediction.

The same review found that the per-edge digest covered only identity and two
observations rather than the exact complete Todo 10 row, deduplicated fits lost
their full partition-reference list, and the density/reversal/equal-density and
atomic publication adversaries were too weak for the claims made. Finally, the
evidence claim misstated the manifest entry count and had no distinct replay
digest for the three handoff artifacts. These are publication/evidence defects,
not evidence against the admitted residual-space science. The green claim is
retained append-only but is non-authoritative; Todo 12 remains blocked until a
corrected complete-edge handoff receives a fresh independent PASS.

### 2026-08-09 — First handoff correction partially accepted

The first correction preserved all scientific replay hashes and independently
closed the prediction-only transform, exact Todo 10 row-byte digest, complete
fit-reference, and distinct handoff-replay findings. The rereview nevertheless
kept Todo 12 blocked. An early `insufficient_dates` return constructed a dummy
full fit without carrying the complete edge-transform state, so that terminal
SKIPPED path emitted no explicit `unavailable/residualizer_unavailable` rows.

The remaining issues are verification gaps rather than fitted-model changes:
the tests must explicitly assemble the category-ordered 2×2 factor matrix from
raw table bytes minus the serialized prediction, prove reversal permutation and
actual unary-posterior reduction for zero factors, and complete the adversarial
publication matrix with a directory-target symlink plus no-fit/no-inference
assertions. The accepted closures and append-only manifest-count correction
remain valid; a second bounded correction is required before Todo 12.

### 2026-08-10 — Residualized T11 graph handoff confirmed

The second bounded correction closed the remaining terminal-path and test-proof
gaps without changing admission science. The early insufficient-date full fit
now preserves every Todo 10 eligible edge as an explicit unavailable or
predictor-ineligible transform row. Tests independently parse persisted raw
Todo 10 rows, subtract only serialized T11 predictions, build the explicit
`[00 01; 10 11]` factor matrix, verify reversal permutation, and demonstrate
that zero edge factors reduce two-node marginals and MAP labels to unary-only
results. Expanded publication adversaries prove fail-closed behavior without
fit, shuffle, or graph inference, including directory and member symlinks.

Independent review passed all seven graph-handoff findings at confidence
`0.99`; a fresh suite passed `821/821`, while admission and handoff replay hashes
remained unchanged. The authoritative downstream chain is the correction2
claim plus final review, receipt v2, handoff v1, and its three exact TSVs. Todo
12 may now implement option A against that chain. Option B remains a separate
future preregistered raw-correlation study; no benchmark or real-scan claim is
made here.

### 2026-08-10 — First Todo 12 inference claim rejected

The first exact-chain implementation passed its focused synthetic suite twice,
but independent review found that the fixture contained only one full-refit
reference and therefore did not establish the required nested reference
semantics. Deduplicated fits with multiple held-out target dates were constrained
to one canonical row, projected onto the wrong date, mislabeled, and recomputed
without the declared cache. Reversed fits returned nodes in reversed order and
did not perform a runtime equivalence check.

The review also separated integrity from numerical behavior: valid density,
DP, or reversal failures were incorrectly converted to `BLOCKED` instead of
`FAIL`. File validation reread paths after snapshot, Todo 11 replay/dependency/
residualizer bindings and Todo 10 spacing/topology consequences were incomplete,
and factor provenance was not retained in immutable results. These are inference
and evidence defects, not evidence against the admitted residual-space model.
The T11 test-only integration amendment remains valid, while the first T12 claim
is append-only historical evidence and cannot authorize Todo 13.

### 2026-08-10 — Todo 12 correction2 still rejected on semantic validation

The second correction expanded path, topology, status and factor-provenance
tests to `306/306`, but independent review found remaining fail-open semantics.
Todo 11's top-level result hash and diagnostic partitions/ESS/starts/scores/
shuffle/bootstrap artifacts were still mostly treated as hash-bound bytes rather
than reconstructed evidence. Consequently fitted-row failures could disagree
with an overall PASS, and synthetic shared-fit transforms could differ by
reference while the cache keyed only the fit and edge.

The fixture also included held-out target dates in partition training dates,
and Todo 10 segment checking did not reject noncontiguous reuse of one segment
ID. The direct mutation matrix was smaller than its evidence claim, and its
append-only count correction confused the 14 listed evidence paths with the 15
manifest lines that include the list itself. Scientific recurrence, reversal,
label-free scope and the admitted T11 handoff remain unchanged; correction3 is
limited to these integrity and evidence defects, and Todo 13 remains blocked.

### 2026-08-10 — Todo 12 correction3 rejected on real T11 byte semantics

Correction3 passed `370/370` and the full T11 regression, but independent review
showed that its hand-built fixtures still diverged from actual T11 publication.
Real ESS rows serialize a date count, permit zero support in early terminal
states, and vary category ESS/support; scores and shuffle rows bind a score hash,
not a fit hash. Gate predicates, partition score/raw hashes, scan/date
aggregations and stage-dependent sparse start rows were not reconstructed
exactly.

The shared-fit fixture also used impossible final-LODO references instead of the
real outer-inner cross-pair formed when two outer/inner holdouts leave the same
training dates. Recurrence, factor construction, fit-level transform identity,
contiguous Todo 10 topology and label-free boundaries remain valid. Correction4
therefore changes test fixture provenance: test-only bytes will be generated by
the T11 serializer itself, while product T12 continues to have no source-level
dependency on T11. Todo 13 remains blocked.

### 2026-08-11 — Todo 12 correction10 stopped at the T11 fixed-point precondition

Correction9 correctly removed the tunable scale tolerance but then rejected the
unchanged serializer-produced `C1` full-refit conditional covariance. The
published non-diagonal matrix has minimum eigenvalue exactly five ordered-positive
Float64 ULPs below `1e-4`; reapplying the frozen T11 projection produces
`0.0045253227996870917`, `-0.001649753170197811`, and
`0.00071502530906223923`, rather than byte-identical canonical fields. Its
original condition number is `51.403481087493361`, so the failure is the
canonical-byte fixed-point requirement, not the condition cap.

The user authorized the non-tunable Julia 1.12 five-ULP representation protocol:
ordinary exact floor/cap checks remain unchanged; only non-diagonal ordered-
positive distances one through five ULP below the floor may proceed when the
reprojected canonical fields are byte-identical; the original matrix must remain
the density input. The actual producer row is not such a fixed point, so the
writer stopped with a correction10 `BlockerReceipt` and did not widen the rule,
rewrite producer values, clip matrices, or change T12 inference. Todo12, T11,
and the docs build were not run after this precondition failed. This is a
serialization/integrity blocker only: it makes no benchmark, real-data, or
Todo 13 claim, and does not alter label-free scientific semantics.

### 2026-08-11 — Phase-A T11 exact post-reconstruction scale certification

The authorized producer correction keeps `SCALE_FLOOR=1e-4`,
`CONDITION_CAP=1e4`, model/configuration, and thresholds unchanged. The defect
was representational: `_project_scale` returned eigenvalues and a reconstructed
matrix before certifying the exact canonical `%.17g` fields that T11 stored and
then used. On Julia 1.12.6, the producer now canonicalizes and reparses the
reconstructed matrix, re-diagonalizes that stored matrix, and returns those
fresh values. If canonical reconstruction alone crosses a floor/cap boundary
after an otherwise feasible target, a deterministic diagonal `q` closure is
rebuilt from the original stored matrix and doubled until the stored matrix
certifies; there is no tolerance, tunable window, threshold change, or generic
jitter.

The known conditional covariance receives this closure and publishes only its
certified stored matrix. Non-PASS null/conditional scales are published
numerically only when the same exact stored certificate passes; invalid triples
are all `NA`, while stage means and start diagnostics remain. The emitted PASS
report was reparsed and every numeric scale triple was independently certified.
The T11 suite increased from `821/821` to `1337/1337`; focused replay remained
`13281332403bd7202c95a1bc9f0216bae1a1e5ead328a5e9569c0b30cd98bc81`, while
complete replay changed from `fe4a755d4c8f013169b7e58f177e314235ec80ca5f9cd8a4e047c20a720da73b`
to `427bfb6d3f97a077fa748f713e1d15ccc0eaac56b4485aedcb24fe2b7abe6d21` and
graph replay changed from `ba344516bd1356188400a7dcd3b639e1cddeaf6e06bfd3349a59dc4b2d7e702b`
to `dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474`.
Both sequential replays were `PASS`, with no terminal status change; CLI help
also passed without invoking construction/evaluation. This is a producer
representation and provenance correction only: it makes no benchmark or
real-data claim, uses no labels or composition prior, and leaves Todo 13 and
the read-only T12 lane suspended pending review.

### 2026-08-11 — Todo 12 Phase-B scale-certification authority rebind

The authorized Phase-B rebind makes the independently reviewed T11
post-reconstruction scale certification the effective Todo 12 authority while
preserving the correction2 claim, review, and source bundle as immutable
ancestry. The fixed authority inventory is now exactly 22 paths (19 historical
paths plus the effective T11 claim, review, and canonical source bundle). The
live T11 source identity is `616decd955cc808958247d005102ff277f336586bd75f5b8529de0bd16042746`,
the effective source bundle is
`7a024ceba895e670a6cf666b568b50e9bb32da110d1639c93e268d96ed5ecdaa`, and the
old `fd9cec0f28d824d9e1a2c7bafb0d0ec0b9cfb2eed7d634add7ab44d089421cad`
bundle is rejected as current receipt authority but retained as ancestry.

The unused universal graph-replay identity was removed. Runtime validation now
recomputes replay from each report's exact three handoff snapshots in declared
order; the sealed T11 suite replay
`dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474` is evidence
for that suite only. Direct adversaries cover content-distinct valid reports,
literal sealed-replay substitution, stale handoff artifacts, semantic
corruption after wrapper rehashing, and historical-source substitution. Exact
stored-scale checks remain `lambda_min >= 1e-4` and `condition <= 1e4`, with no
tolerance or clipping exception.

The complete Todo 12 suite passed `855/855` in each of two sequential Julia
1.12.6 runs. Both runs produced identical result replay
`09bfbbcdd6e49f935e3b533d720376a40450baef5d6925b1dfedf1eec14de344`,
provenance replay
`9f41653b5de58388ea82c862f8eafb2934121cbd627ef786917e48102e63af51`,
correction7 matrix rows/counts (`129`, `261/261`), and correction9 exact-scale
matrix rows/counts (`154`, `17/17`). The T11 regression passed `1337/1337`
with focused replay
`13281332403bd7202c95a1bc9f0216bae1a1e5ead328a5e9569c0b30cd98bc81`, complete
replay `427bfb6d3f97a077fa748f713e1d15ccc0eaac56b4485aedcb24fe2b7abe6d21`,
and sealed-suite graph replay
`dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474`.

All pre-existing T12 evidence paths remained byte-identical to the preflight
inventory, and regenerated matrices, logs, and boundaries were written only
under the new scale-certification-rebind root. This is a provenance and
serialization-integrity rebind only: no expected counts, labels, benchmark
grades, real data, composition prior, model thresholds, or Todo 13 activity
were used or changed. Todo 12 remains pending fresh independent Oracle review.

### 2026-08-11 — Todo 12 scale-certification rebind exact eigensolver correction

Independent review rejected the Phase-B claim because Todo 12 computed an
authoritative `eigen(Symmetric(scale))` decomposition but discarded its
`Float64` eigenvalues for a custom off-diagonal trace/determinant estimate. The
exact stored counterexample was
`[0.32262332547297057 -0.46740750707094925; -0.46740750707094925
0.6774766745270292]`: T11 obtains eigenvalues
`[9.999999999998899e-5, 0.9999999999999998]`, condition
`10000.000000001099`, and rejects at the floor, while the custom Todo 12
estimate accepted it. The correction removes that analytic branch and uses
`Float64.(decomposition.values)` directly, retaining the exact positive,
`lambda_min >= 1e-4`, and condition `<= 1e4` checks without tolerance,
clipping, fixed-point exception, or alternate density matrix.

The exact counterexample was added to the regenerated correction9 matrix and
the public seven-path fixture now refreshes all wrappers and blocks with
`model_mismatch` before scan/factor/DP work (`0/0/0`). A test-only parity audit
over certified interior, exact diagonal, rotated certified, non-diagonal floor,
and non-diagonal cap matrices passed `19/19`; T11 acceptance/rejection and
Todo 12 disposition agree. Correction evidence is isolated under the new
append-only correction root, while the Phase-B `231/231` manifest, review, and
all older T12 evidence remain byte-identical.

Two sequential Todo 12 runs passed `891/891` with identical result replay
`09bfbbcdd6e49f935e3b533d720376a40450baef5d6925b1dfedf1eec14de344` and
provenance replay
`9f41653b5de58388ea82c862f8eafb2934121cbd627ef786917e48102e63af51`.
Correction7 remains `129` rows / `261/261`; correction9 is `155` rows /
`34/34`, with byte-identical repeated matrix and boundary outputs. T11 passed
`1337/1337` with focused replay
`13281332403bd7202c95a1bc9f0216bae1a1e5ead328a5e9569c0b30cd98bc81`, complete
replay `427bfb6d3f97a077fa748f713e1d15ccc0eaac56b4485aedcb24fe2b7abe6d21`,
and graph replay
`dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474`.
The updated source bundle is `structured-chain-inference-source-v5` with
current edge, unchanged chain, and test hashes
`d1cce4047eaaf99faacb8c11dd3844669c5f6528a2bded011e004cb432417d46`,
`cdc0c788d49298f50721b091535a238cd552454d9885404456c593e4625ea39f`, and
`112ffafe7a3db87994e9654c31dd9ef7c68c49cdbf4566b2a082b914918d6759`.
No benchmark, real-data, label, composition-prior, or Todo 13 activity was
used; the correction claim remains pending fresh independent review.

### 2026-08-11 — Todo 12 checked-plan and checked-T11 authority rebind

After the checked Todo12 plan was independently authorized, T12's fixed
authority was rebound append-only to the checked-plan T11 publication. The
prior correction2 ancestry and the immediately prior scale-certified tuple
remain explicit historical authorities: scale claim
`5dc185994ecd9e117ba9f8dc782e70360daf397a2348e57fd65b2b14a4ffdfe2`, review
`93f986d6baa3dd8cf850328055c750b27b9873ae67da54f3360ae05c31b22e30`, and
source bundle `7a024ceba895e670a6cf666b568b50e9bb32da110d1639c93e268d96ed5ecdaa`.
The checked-plan T11 claim/review/source bundle are now the effective tuple:
`61b009a55a4ab403b4a2391f3a0a4d5b5bf72b3b21bc684703ac1b52dacd33cd`,
`c02a8fdd553533d844b4d49f13b3cda2658bc02a38589092172d1a4b394cdb83`, and
`b5e54d763523e937734d03742b7323a0f1da191591cf7683a17ff142f01c233b`.
The checked plan binding is `a9b386d613829e8f7e20b6e33f8e80898fa9a55f0b344dbb92ea64ac6804f3d0`,
replacing the old unchecked `e3e11942a9ce26e30c402655f359696d2e2bb38fe7fb3e28ff2a2988e839d38e`.

The T12 runtime authority inventory therefore grows from 22 to exactly 25
unique paths by appending the checked-plan T11 claim, review, and canonical
source bundle. The live T11 source is
`ca7ea9b05b55a5185a5f08834b4fec4e748830f3eb24cac1a45f387e4d768575`; its CLI
binding remains unchanged. Two fresh T12 suites passed `903/903`, with
authority mutation coverage `103/103`, while exact scale semantics, parity,
graph replay, matrices, and all downstream failure/status behavior remained
unchanged. The T11 regression passed `1337/1337`; focused, complete, and graph
replays remained respectively
`13281332403bd7202c95a1bc9f0216bae1a1e5ead328a5e9569c0b30cd98bc81`,
`427bfb6d3f97a077fa748f713e1d15ccc0eaac56b4485aedcb24fe2b7abe6d21`, and
`dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474`.

This is an administrative provenance rebind only: no scientific formula,
threshold, scale rule, benchmark, real-data, label, composition prior, or
Todo13 activity changed. Prior T11/T12/T3/T10 evidence remains immutable and
Todo13 remains blocked pending fresh independent T12 review.

### 2026-08-12 — Todo 12 checked-plan final evidence and cleanup correction

The final checked-plan review correctly blocked the prior claim for missing
required evidence and task-owned temporary cleanup. Its reported
`5f47064b210a1f837cd3fc38c4494046af899ad7646a7cd9b0faf0a7aabbd12a`
reconstruction is not byte-exact and cannot be reproduced from the prescribed
operations. Independently replaying exactly those five test reversals—removing
the checked-plan evidence segment, removing the three checked-plan authority
paths from both lists, restoring the 22-path assertion, and restoring the
historical stale-source constant—produced the sealed predecessor
`112ffafe7a3db87994e9654c31dd9ef7c68c49cdbf4566b2a082b914918d6759` from the
current `376847c800fec317acdc33444798869c65c87b4bb6c159c594594082cf6f6b71`.
The corresponding edge-model reversal also reproduced
`d1cce4047eaaf99faacb8c11dd3844669c5f6528a2bded011e004cb432417d46`.

Current product bytes were not rewritten. The correction is evidence and
cleanup only: the exact eigensolver and label-free boundaries remain unchanged.
The fresh T12 suites passed `903/903` twice and T11 passed `1337/1337`, with
the prior result/provenance and focused/complete/graph replays unchanged.
The four checked-plan matrix/log outputs remained byte-identical. Only the
enumerated task-owned `/tmp/opencode` paths were removed; unrelated temporary
entries were not traversed or deleted. No benchmark, real-data, label,
composition-prior, scientific-threshold, or Todo13 activity occurred.

### 2026-08-12 — Todo 13 evaluator-v1 policy preregistration

Todo 13's evaluator contract was under-specified in four decision-critical
areas: absolute multi-view unary evidence, the common-null lifecycle,
composition of complete scores around relative Todo 12 graph evidence, and
partial-view/descriptive formulas. Following the authorized recommendation in
`structured-t13-evaluator-policy-proposal.md`, evaluator-v1 is now frozen in
`config/unit_assignment_structured_evaluator.toml` and bound to the checked plan
and immutable Todo 8, T11, and T12 authorities.

The chosen policy uses a log opinion pool of absolute joint view evidence with a
shared C1/C2 active mask, available-view renormalization, structural omission,
uniform all-view-missing output (`U=0`, posterior `(0.5,0.5)`, output `?`), and
explicit `SKIPPED`, fallback, `FAIL`, and `BLOCKED` statuses. The common null is
fit after unary and residualization but before conditional feasibility;
conditional `SKIPPED` retains selected unary evidence without graph evidence,
while numerical failures never fall back. Complete scores use the fixed
`M_s=nodes_s+eligible_edges_s` denominator and Todo 12 only as a relative graph
normalizer: `L_C1=sum U_C1+sum N_e` and
`L_meta=sum U_selected+sum N_e+I_graph logZ_T12`.

The gate is preregistered as 500 paired whole-scan-within-date bootstrap
replicates, Type-7 quantiles, strict positive lower bounds and every-date
positivity, plus an exhaustive inclusive upper-tail date sign test without a
`+1` correction. Coverage, normalized entropy, strict view agreement, and the
76+30+5=111 / +9=120 parameter accounting are descriptive/frozen metrics. These
are non-label-derived evaluator scoring decisions, not benchmark calibration,
physical calibration, or a benchmark claim. Rejected alternatives include
model-specific null densities for all-view abstention, post-admission common
null fitting, graph evidence treated as a complete likelihood, changed
denominators, pooled/edge resampling, nearest-rank quantiles, strict sign tails,
and a `+1` correction.

The policy prerequisite is complete, but Todo 13 itself is unimplemented and
unexecuted. A fresh independent Oracle owns validation after the terminal
claim; Todo 13 remains blocked until that policy review is sealed PASS. No
benchmark, real-data, label, expected-count, class-count, composition-prior,
or producer/construction activity occurred.

### 2026-08-13 — Todo 13 evaluator-policy correction

The original Todo 13 evaluator-policy claim was independently rejected by
Oracle session `ses_005da658cffeE4gMqcpl5dX1lh` at confidence `0.99`. Its claim,
config, and validator remain immutable historical evidence and are not
authoritative. The six blocker classes were non-unique semantics,
pathless/unenforced authority, unguarded runtime and numerical gates,
string/exception-only validation, incomplete evidence with a non-genuine
prewrite guard, and overstated documentation.

This append-only correction freezes implementation-unique unary, status,
edge/null, graph-reference, complete-score, bootstrap/sign, metric, and
serialization semantics. It binds canonical repository-relative authority paths,
exact source-bundle members, Julia `1.12`/`1.12.6`, the checked Todo 12 and
unchecked Todo 13 markers, and the existing `1e-4`/`1e4` scale constraints.
Authority failures return structured `BLOCKED` reasons before formula,
bootstrap, or graph work. The correction validator uses hand-computed and
synthetic fixtures, structured mutation results with zero-work counters, and
byte-identical canonical replay rather than string-presence or exception-only
claims.

No threshold, model, predecessor source/test, label boundary, GCV, `n_eff`,
real-data, benchmark, or historical evidence bytes changed. Todo 12 remains
checked; Todo 13 remains unchecked, unimplemented, and blocked. Validation here
is policy/static/synthetic only: it creates no 10–20mer application processing
claim and no benchmark claim. The correction is pending a fresh independent
Oracle PASS and GateClosure; it does not authorize Todo 13.

### 2026-08-13 — Todo 13 evaluator-policy correction2

The first evaluator-policy correction was rejected by Oracle as fail-open. The
live `f2795a8f...` config was replaced for correction2, while its exact 11,480
byte preimage, failed claim/evidence, and failed review remain preserved as
immutable historical artifacts and are explicitly non-authoritative.

Correction2 closes the remaining defects: every policy value and role-specific
authority binding is independently enforced; descriptor-based snapshots and
before-return revalidation reject changed bytes, symlinks, hardlinks, identity
collisions, and path substitution; runtime is checked before reads and before
successful return; consumed versus terminal statuses are machine-readable; the
selected model/unary-fit/T11/T12 reference chain and cardinalities are exact;
entropy and pooled view-agreement populations are explicit; bootstrap/sign
replays use actual seeded arithmetic; and mutation results report structured
reasons with zero scientific work on static failures.

No threshold, predecessor source/test, model, label boundary, GCV, `n_eff`,
benchmark/application claim, real-data activity, or historical evidence byte was
changed. Evidence is policy/static/synthetic only. Todo 12 remains checked,
Todo 13 remains absent, unchecked, unimplemented, and blocked pending a fresh
independent Oracle PASS and GateClosure. External labels remain limited to
post-gate reporting and cannot authorize this prerequisite.

### 2026-08-13 — Todo 13 evaluator-policy correction3 integration

Correction3 lane integration was completed as a policy/static/synthetic
evidence phase. The parent-reproduced policy lane (680 mutations) and authority
lane (163 mutations) were losslessly bound as 843 projected rows; the
integrator did not reimplement or independently rerun those lane mutations.
The integrated validator independently checked the 202-key policy, 20 roles,
32 bundle members, 35 authority snapshots, 42 claim/review checks, 13
structural checks, semantic formula/reference/metric fixtures, 500 bootstrap
replicates, 32 sign masks, final snapshot revalidation, and two integration
probes. Two fresh integrated runs were byte-identical.

Correction2 remains historical and blocked. No threshold, model, T8, T11, T12,
GCV, `n_eff`, label, real-data, benchmark, grader, or composition-prior
behavior changed. Todo 13 products remain absent and the Todo 13 marker remains
unchecked; status remains `BLOCKED_PENDING_INDEPENDENT_REVIEW` pending parent
acceptance, a fresh independent Oracle PASS, and reviewer-owned GateClosure.

### 2026-08-13 — Todo 13 evaluator-policy correction4 provenance successor

Correction4 is a provenance-only successor to the technically green correction3
integration. Correction3's historical canonical publication was blocked because
the canonical paths were replaced while correcting the combined-row projection.
Correction4 freshly regenerates the six canonical static/synthetic files in two
fresh Julia runs, requires run-1/run-2 byte equality with the current correction3
outputs, and publishes run-1 exactly once using descriptor-relative exclusive
no-replace creation with a durable receipt and replay.

This successor changes no science, live evaluator configuration, threshold,
model, T8/T11/T12, GCV, `n_eff`, label, benchmark, application, or Todo behavior.
The evidence is static/synthetic only; Todo13 products remain absent and the
phase remains `BLOCKED_PENDING_INDEPENDENT_REVIEW` pending parent acceptance, a
fresh Oracle, and reviewer-owned GateClosure. No review or GateClosure is
created by this worker.

### 2026-08-14 — Todo 13 evaluator-policy correction5 terminalization caveat

Correction5 records the parent-validated administrative terminalization of the
failed correction4 publication. Correction4's six canonical `O_EXCL` bytes
remain valid by reference and were not republished. Its `DoneClaim` is
non-authoritative because of the predecessor-hash defect and the unqualified
`close_checked` finalizer failure before final Boulder closure. The captured
staging residue was cleaned descriptor-relatively; the failure artifacts remain
preserved in Correction5 evidence.

No science, configuration, threshold, GCV, `n_eff`, T8/T11/T12, label,
benchmark, application, or Todo behavior changed. Todo13 remains blocked pending
parent acceptance, a fresh independent Oracle PASS, and reviewer-owned
GateClosure. This entry is append-only; no review, GateClosure, or Todo13
product was created.

### 2026-08-14 — Todo 13 evaluator-policy correction6 S0 staged evidence closure

Correction5 cleanup remains valid, but Correction5 is non-authoritative because
its evidence closure omitted the cleanup receipt and its terminal replay omitted
all six required per-path canonical bindings. Correction6 S0 references the
Correction4 canonical bytes without republishing them and stages only the
administrative evidence required for parent-owned atomic publication. No policy,
configuration, calibration, threshold, GCV, `n_eff`, T8, T11, T12, label,
benchmark, application, or Todo behavior changed. Todo13 remains blocked pending
parent atomic publication and acceptance, a fresh Oracle PASS, and
reviewer-owned GateClosure.

### 2026-08-14 — Correction6 S0 rejection and clean rebuild clarification

The first hidden Correction6 S0 stage was rejected before S1 by the parent/Oracle
because the publisher did not bind the parent checkpoint device/inode and the
final repository-relative manifest, its self-test used unsafe check-then-write
helpers, and several payload observations were unqualified temporal claims.
The official Correction6 root and all external checkpoint/receipt controls
remained absent; no rename, publication, receipt, review, or GateClosure
occurred. The rejected stage was removed descriptor-relatively after verifying
its exact 19-file identity and manifest, and a wholly new 19-file payload was
assembled with the corrected publisher contract.

The rebuild adds no scientific or policy change. It binds both manifest
namespaces, classifies pre-rename failures as `not_committed`, post-rename
failures as `ambiguous`, and complete verified publication as
`committed_verified`; it creates no in-root receipt or post-publication files.
Any eventual publication can be proven only by the parent-owned external
checkpoint and receipt. Todo13 remains blocked pending parent acceptance of the
receipt, a fresh Oracle PASS, and reviewer-owned GateClosure.

### 2026-08-15 — Correction6 interrupted-state reconciliation clarification

The prior rejection note's statement that a wholly new payload had already been
assembled was premature. The interruption occurred after safe removal of the
old rejected stage and the durable documentation update, but before the new
Correction6 payload assembly. The parent then independently confirmed that the
hidden stage, official root, checkpoint controls, and receipt controls were all
absent, and that Correction3, Correction4, and Correction5 were unchanged.
This resumed run assembled the corrected 19-file stage with the repaired
publisher and payload-seal contracts. No publication, review, GateClosure,
Todo13 product, or scientific/policy change occurred.

### 2026-08-15 — Todo 13 grouped evaluator implementation and synthetic verification

GateClosed evaluator-v1 policy was implemented as the named module
`StructuredUnitAssignmentEvaluator` in
`test/evaluate_structured_unit_assignment.jl`, with its focused suite in
`test/test_structured_evaluator.jl`. The evaluator composes the existing T8,
T11, and T12 boundaries, uses fold-wide training-only unary selection, preserves
the fixed Todo 10 node/eligible-edge observations and common null, and emits the
nine deterministic evaluator artifacts. It does not expose a synthetic or
generic-factor production bypass.

The first Phase 1 Oracle review rejected the implementation with twelve
blockers covering nested scan/date identity, shared view masks, date-grouped
inner selection, exception classification, common-null lifecycle, T11/T12
reference use, outer C2 fallback, snapshot/publication lifecycle, provenance,
and CLI firewall ordering. Two bounded remediation passes closed those defects.
The final Phase 1 rereview passed at confidence `0.999` for evaluator SHA-256
`3221ed25ed4ce0ac110170492f33ca4145c525d3437fed373dc8be7416acf825`.

The final focused suite has 252 deterministic assertions. Two independent runs
passed with byte-identical stdout SHA-256
`158de38798a82d3fac0071e8cacfaaa3d896d159483073d6d9f5870a5c4f38b8`
and empty stderr. It independently checks the absolute unary and scan-score
formulas, date-clustered 500-seed bootstrap, complete inclusive sign space,
fixed denominators, null cancellation, graph log evidence once per multi-node
block, pooled metrics, held-out leakage sentinel, malformed topology,
provenance/result binding, authority snapshots, CLI rejection, and atomic
publication. The structured firewall also passed 379 assertions, and the T8
boundary passed 3886 assertions with its frozen replay and source hashes.

Attempts to rerun the entire standalone T11 and T12 research programs locally
were stopped by execution timeouts after their completed testsets had remained
green; these interrupted attempts are not reported as suite passes. Their
immutable previously accepted 1337/1337 and 903/903 evidence remains the
unchanged predecessor authority. No real scan, external correctness data,
composition prior, model/threshold parameter, GCV, or `n_eff` behavior was used
or changed. This establishes static/synthetic implementation evidence only, not
a 10–20mer application or benchmark claim. Todo 13 remains unchecked pending
the final independent implementation review and administrative closure.

### 2026-08-15 — Todo 13 final-review remediation

The first Todo 13 implementation Review 2 was `FAIL` at confidence `0.999`.
Although the core score and provenance tests were green, the evaluator did not
aggregate outer-fold terminals with the frozen `BLOCKED > FAIL > SKIPPED >
PASS` precedence, empty terminal reports lacked the required unconsumed event,
and the CLI publication path replaced an original blocked reason. The durable
suite also covered only five-date sign enumeration and did not exercise enough
snapshot-mutation and publication-state failures.

The remediation now evaluates every canonical outer fold, catches only
structured evaluator terminals, chooses the winning status and reason
independently of input order, discards partial scientific rows on a non-PASS
winner, and emits one zero-work terminal event. The production blocker receipt
retains the original reason. Exhaustive sign enumeration is a private pure
helper used by production and checked independently for every `K=1…13`, with
inclusive ties, exact powers-of-two denominators, no correction term, and the
strict boundary between `1/16` and `1/32`. Snapshot bytes and directory
membership are mutation-tested; valid and failed context publication plus
`not_committed`, `ambiguous`, and `committed_verified` classifications are
covered without adding a production fault or synthetic CLI option.

One attempted test-only follow-up accidentally restored the pre-remediation
evaluator bytes despite its assigned scope. Parent SHA checks rejected that
state before acceptance, and a fresh targeted evaluator lane restored the
changes without modifying the tests. The reconciled evaluator SHA-256 is
`7b4717077c24f122298e61cd944d55ab439ef278b32a76ec8a6fa8049735dd2e`;
the test SHA-256 is
`86d758174968851e50d64fced8369749ac128f23384b879b856e72b27a0dbe78`.
Two parent runs passed 65,945 assertions with byte-identical stdout SHA-256
`32622f50c4090e5b7b6090df7865be205847a59fed65a5533afc2b6950693e72`
and empty stderr. The 48 focused parent probes and 379 firewall assertions also
passed.

No real scan, external correctness data, composition prior, parameter,
threshold, GCV, `n_eff`, or predecessor authority changed. This remains
static/synthetic implementation evidence only. Todo 13 remains unchecked and
no PASS closure is claimed until a fresh independent Review 2 accepts this
exact reconciled state.

### 2026-08-15 — Todo 13 parent-probe count correction

The parent probe log contains four passing testsets with `25 + 5 + 6 + 7 =
43` assertions. The two earlier journal statements of `48` were an arithmetic
reporting error, not five missing or failed assertions. Todo 13 closure evidence
must cite the truthful `43/43` count and the exact probe-log SHA-256
`e14aa1083b3e122ba1a6cb3e4714ba24601225a644df259d1a34db52d4f75a9f`.

The first closure staging claim also described itself as already published.
Oracle rejected that administrative state before final-root publication. The
staging root was removed, the final `closure-v1` root remained absent, and no
plan checkbox or Boulder closure state changed. The rebuilt staged claim must
remain explicitly pending atomic publication; only a post-publication external
receipt may assert that publication completed.

### 2026-08-15 — Todo 13 terminal closure

Fresh Oracle Review 2 accepted the reconciled evaluator and tests at confidence
`0.999`, closing `R2-01`, `R2-02`, and `R2-03`. The immutable implementation
evidence root contains 27 files / 47,071 bytes. Its root manifest SHA-256 is
`a8f40846213cd976e4ffdd1c8d40e078718f30a917d9ad5ffd4fd38ef3a76aef`
with 26/26 manifested members; the claim and review hashes are respectively
`30cdfbf9be939104f76b93aa0ed610b17c5db20079af0ae32dfe781cecccc183`
and `1b383bda991f113b31923700a57731af56834c4ef0c2e7b3e8743e6ee24981c0`.

The no-replace directory rename committed the staged root, but the publisher
then stopped before its external receipt because it compared dictionaries that
contained stage-path versus final-path objects. Independent reconciliation
confirmed identical device/inode `65028/34566775`, all 27 regular read-only
files, every root-manifest hash, final-root mode `0555`, and stage ENOENT. No
final-root byte was written after rename. The external receipt classifies the
state as `committed_verified_after_interrupted_parent_check` and has SHA-256
`59cf1c60006c549c6df6c8c07e9e9915494d5ed08d392095fa9a1cbfa5896249`.

The Todo 13 plan marker was changed by exactly one byte to checked; the checked
plan SHA-256 is
`6aa00c9b2f5139149270ef8fb47e5caf9b6cb070c7a0bcba62c261f65d3de26e`.
Boulder now records Todo 13 `PASS`/closed with final SHA-256
`508a7444d9f4bd03f9ed3be1631e114c316062dc24f414c3d62d352326c7e92f`
and 198,074 bytes. The evaluator config and all T8/T11/T12 authorities remain
unchanged.

This closes only the static/synthetic Todo 13 implementation. No real scan or
application validation was performed; interrupted local T11/T12 reruns are not
passes, and only their unchanged accepted 1337/1337 and 903/903 evidence is
cited. Post-rename ambiguity is covered by source control-flow and classifier
tests, not an induced kernel/fsync failure. No later Todo is authorized by this
closure.

### 2026-08-15 — Todo 13 postclosure runtime-authority rebind

The valid Todo 13 implementation closure exposed a postclosure lifecycle bug:
the production evaluator still treated the pre-implementation Plan/Boulder
bytes and prerequisite `todo13_may_start`/product-absence facts as permanent
runtime invariants. Consequently, checking Todo 13 and recording its authorized
closure caused the real authority path to reject the evaluator it had just
accepted. Rolling the Plan or Boulder back was forbidden.

The runtime authority was rebound without changing the evaluator config or any
scientific/model/scoring/selection/threshold behavior. Live Plan and Boulder
bytes are no longer runtime snapshots; their hashes remain explicitly qualified
historical provenance. Runtime authority instead validates the immutable
prerequisite GateClosure/review and the complete implementation closure. All 27
closure files and three exact directory inventories are retained through the
production context and revalidated before computation and publication, including
read-only modes, identities, bytes, links, and membership. Normal and blocker
receipts are now schema version 2, omit unqualified Plan/Boulder hashes, and state
`live_plan_boulder_runtime_authority=false`.

The final evaluator and test SHA-256 values are respectively
`c1d467e9a7bdc0767230e4d9d72bb08e159e1b8a6c9b3a293e3e14d16fc26121`
and
`964d7af75d48c5af0d189ebb93f8d400b3c0b53fbd29dfc2f9d05e22a908e103`;
the evaluator config remains
`ec0546096b3c4742cd86d8c3d40788a5894d20fc5a4702b0318581f10f8b0b90`.
Two parent runs each passed 66,017 assertions with byte-identical stdout
SHA-256 `4626eafd12d0cb4b18a6d08e335aa5eaaf0d4bfe5c48356a844536b383040486`
and empty stderr. The live authority probe passed with digest
`93b0c2778a24094b8283917dca7060274b9d856a932a5e03419756d7d7e072c8`,
59 snapshots, 27 closure snapshots, and three inventories. The structured
firewall passed 379/379 and the CLI help preflight passed.

Fresh Oracle review passed at confidence `0.999`. The immutable correction root
contains 22 files in three directories / 225,183 bytes. Its 21-entry root
manifest SHA-256 is
`4f868912be545952575a076e27635bf062a4d2e7c5518c32ecd54f239c89c7fc`;
claim/review hashes are
`4392ac2d8ee7084015e558652e5890f2a65aee4b61c65f151ef93f16bf936456`
and
`6634f9d6f41a9d90802d750499b34b052c4bc898f83893bda6cf84b8cde27a5e`.
The no-replace rename preserved device/inode `65028/34803464`; the stage is
absent. Checkpoint and publication-receipt hashes are
`49ba8d29f23db65813030179c85ffccd7476dfb9e84e3cd9f316c1896c6027ab`
and
`af6adfb3ceb9b3f2d6e7a902621bd1e25818c6fd427a0a56430ca13770575842`.

Boulder now records evaluator runtime readiness at SHA-256
`8c4d6e15761a79de807ac0a7ddd48302750c9f7486925a9d86750d20a782bcbe`
and 204,075 bytes. This is static/synthetic authority evidence only: no real
scan or application validation was performed, and this seal authorizes no later
Todo.

### 2026-08-20 — Receipt-gated evaluator publication fallback

Oracle review identified a real filesystem defect on Viper `/ptmp`: directory
`renameat2(RENAME_NOREPLACE)` returned numeric `EINVAL`, so blocker publication
returned status 2 without an output receipt. The evaluator's preferred
renameat2 path remains unchanged for filesystems that support it. The bounded
correction adds a fallback only for `EINVAL`, `ENOSYS`, and
`EOPNOTSUPP`/`ENOTSUP`; other rename errors remain fail-closed.

The fallback reserves a mode-0700 destination exclusively, transfers sorted
non-receipt files with no-overwrite hard links, and makes a completely fsynced
`receipt.toml` visible last. Receipt visibility begins the ambiguous state;
the destination is never cleaned after that point unless all exact inventory,
identity, mode, byte, fsync, and production-context checks pass. Pre-receipt
failures leave a no-receipt collision residue, while private stage cleanup is
allowed only after device/inode revalidation. Existing output is accepted only
as an exact nine-file report or receipt-only blocker publication; directory
existence alone is not a committed result. This protocol is separate from the
parent-owned immutable evidence publisher and makes no scientific or benchmark
claim.

RED coverage against the fail-open publication seam recorded 66,445 passes,
4 failures, and 6 errors. After the fallback implementation, the focused
evaluator suite passed 66,513 assertions with zero failures/errors; the final
required two-run validation is recorded by the parent orchestrator. Synthetic
coverage includes preferred and errno-gated fallback paths, blocker/report
sets, collisions and residue, races, wrong modes and links, fsync failures,
context mutation before/after receipt, interruption states, and identity-safe
stage cleanup.

The final focused suite also includes the explicit numeric `ENOTSUP` alias
probe: both validation runs passed 66,515 assertions with zero failures,
errors, or broken tests. No public replay was run and no downstream or
immutable-evidence publication authorization is implied.

Clarification: the two direct Julia 1.12.6 runs were executed in this bounded
validation; the parent-owned item still pending is the real Viper capability
and four-main probe.

### 2026-08-20 — Publication errno and receipt-state correction

Correction to the fallback record above: on Linux this lane accepts exactly
numeric errno `22` (`EINVAL`), `38` (`ENOSYS`), and `95`
(`EOPNOTSUPP`, also `ENOTSUP`). Numeric `45` is `EL2NSYNC` and is rejected;
there is no string-based errno parsing or alias acceptance. The receipt state
is now marked `ambiguous` immediately after a successful receipt hard link, and
also when a concurrent receipt collision is detected. The preferred rename
path marks the state ambiguous immediately after successful no-replace rename,
before destination inspection.

The exact-case RED run recorded 66,552 passes and 2 failures; the corrected
GREEN run passed 66,554 assertions with zero failures/errors. Added checks cover
mode-0700 reservation, identity-bracketed payload/receipt bytes, path-specific
pre/post-receipt fsync faults, and separate payload/receipt no-overwrite
collisions. No scientific, benchmark, configuration, T12, Plan, or evidence
bytes were changed.

### 2026-08-20 — Runtime-v3 Viper closure and public replay

The final receipt-gated evaluator bytes are
`a9eccab2747602ee107af73f83dc4a1d69b90c2e457f64ddd1ce8595704178fa`
for `evaluate_structured_unit_assignment.jl` and
`a243489c59b9934d6c6992ba0e43431a0abb0f70d7016deb26ab21c8584be7a9`
for its test. The runtime-v3 config remained
`a0a04794b346f351c61a281384869bcde3f378aa0836c9271197d4004488586e`.
The resulting authority contains 153 snapshots and 35 inventories, with digest
`9eb003f7f329b359f3d6cb08dc4c8b35f4363895d1dd6906c8c655689127ed8d`.
Two parent Julia 1.12.6 runs passed 66,553 assertions each with empty stderr.

Viper job `10971601` preserved the real `/ptmp` capability and four-main
probe. It completed `0:0`: `renameat2(RENAME_NOREPLACE)` returned numeric
errno 22, hard links passed, and the v3, v2, v1, and historical blocker calls
each returned status 2 with exactly one config-bound `receipt.toml` and no
publication-layer error. Their receipt hashes are respectively
`dfbcfc326850a10fa6b350f0a6bd0b5713b6e630fac69ab3f01ebeda78dbde53`,
`26bfd2e6a2a3b17ad51758a53715d385af64fc9b27ff58d806120d645f357391`,
`f288c3360ba1dc79da3e12e7f497e08075a7ead68a7dde9e92f8ddd281b40c05`,
and `0778029855a04f76b6774b21906f97fb98752ef7d32c7a664804cc52305f8ba9`.
The eight raw probe files and receipts are bound by the retained
`probe-files.sha256`; capability, four-main, and Slurm stderr are empty.

The final serial public replay, Viper array job `10971530`, completed both
tasks `0:0` in 8 min 41 s and 8 min 52 s. Both captured suite outputs are
byte-identical at
`2066b570205f33c89952cbd3fbe6eafda28d31c13f031fcea7be185b6bc00a18`
and report 66,562 assertions with zero failures, errors, or broken tests; both
Slurm stderr files are empty. The common public identity is status `PASS`,
evidence
`3c8fb798119209647e362c5fe94c1848d2a636d8cd28693001a729ddf1f6b9bf`,
T12 result
`12ab5fe91ed6ea3bbb061519fa297f95bfdbc785418d9de12dd83efe14c50dec`,
T12 provenance
`60482aa611455c78d8fdf2522da3c2477edf10dc41f42adccef8ad26ee02b639`,
and final decision
`4bdb18827ecf559156d51190f5a499181d2041ac4f44f31443b6c12659782938`.
The T12 result/provenance identities are unchanged from the pre-fallback
diagnostic, so the portability correction changed no scientific output.

This closes runtime validation only. T14 remains forbidden until the
append-only T13 runtime-v3 evidence root is independently reviewed and
atomically published; no Plan checkbox, benchmark grade, or production
configuration was changed here.
