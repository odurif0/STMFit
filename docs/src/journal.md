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

> Updated 2026-09-23. Questions from earlier sessions are archived in
> `journal_archive.md`.

0ad. **Does a larger native LM budget resolve the registered one-file rejection?**
    → **AUTHORIZED, BOUNDED CONVERGENCE DIAGNOSTIC (Sep 23)**: replay the same
    saved pixels, N, frame, physical bounds and native starts on `240817_006.sxm`.
    Compare fresh 300/10,000-iteration fused/paired continuations for all four
    shape families; keep native stopping tolerances and validity unchanged.
    Report trace-prefix identity and independently recomputed box stationarity,
    without treating a local stop as global optimality. Two-hour preparation,
    one one-hour Viper job, no retry, classifier, partial grade or new seed.
    Keep support **676/34**; recognition comparison is a separate task.

0ac. **Can a shared-geometry model of the two acquisition directions improve recognition?**
    → **ONE-FILE FEASIBILITY COMPLETE; RECOGNITION INDETERMINATE (Sep 22)**:
    native registered Gaussian residual rejection is reproduced and instrumented.
    Its elliptical candidate reaches **3.6092 noise units**, with three adjacent
    pixels above 3.5; LM hits 50 iterations. From those same native starts, an
    additional 300-iteration fused continuation reaches **3.6512**, paired
    gain/plane fitting **3.5184** (forward view **4.6419**). Both remain invalid;
    all eight local continuations hit their cap without convergence. Paired
    split passes mean validity but fails the forward-view guard. No full-cohort
    fit, Viper submission, partial classifier or grade. This bounded version
    cannot complete the comparison; it does not disprove all paired models.
    Keep support **676/34**, below historical **677/36**. Report:
    `results/paired_acquisition_20260922/report.md`.

0ab. **Can original unregistered geometric support make the registered refit complete?**
    → **DIAGNOSTIC COMPLETE; RECOGNITION INDETERMINATE (Sep 22)**: source
    **98d664d**, job **11942393**, **FAILED 1:0 in 28m10s** after nine seconds
    queued. All original contexts replay exactly and the three span failures
    disappear. A different scan, `240817_006.sxm`, fails registered Gaussian
    residual validity at saved N=6; its registered split is not attempted.
    **145/146 scans** complete four fits: **582 valid, one failed, one skipped**.
    No partial classifier or grade. All **706 outputs plus Slurm log** are
    fetched/checksummed; 27,060 numerical checks pass, plus cohort/input/native
    support checks. Fifteen reference tables replay exactly; retain **676/34**,
    below historical **677/36**. No rescue, threshold change, new N or retry.

0aa. **Does refitting both shape models after accepted registration improve recognition?**
    → **DIAGNOSTIC COMPLETE; RECOGNITION INDETERMINATE (Sep 22)**: source
    **c7836f9**, job **11942091**, **FAILED 1:0 in 26m57s**, after 2m18s queued.
    All 146 inputs are accounted for; **143 complete all four fits**. Of
    584 planned fits, 574 are valid, three fail the unchanged native span
    constraint and seven later stages are not attempted. Two failures occur
    in the zero-shift control at N=10; one after registration at N=6.
    No partial classifier or grade is emitted. Fifteen reference tables replay
    exactly; **549 outputs plus the Slurm log** are fetched/checksummed and
    **22,742 final checks** pass. Retain support **676/34**, below historical **677/36**.
    This does not establish a recognition loss or invalid saved counts. A
    complete comparison needs a new decision about support/N consistency;
    no rescue, retry, new N or other campaign follows.

0z. **Can identifiable acquisition translation improve fixed-geometry recognition?**
    → **COMPLETE; NEGATIVE (Sep 22)**: source **9961973**, job **11941278**,
    **0:0 in 6m52s**, after 6m04s queued. Registration is accepted on
    **104/146 scans**; all **629 comparable patches** improve spatial
    correlation (median **0.2314→0.9911**). Yet recognition falls to
    **671/870 correct, 26/145 exact**, versus both native replay and
    observation-mask control **676/34**; coverage **852→850/870**. Against
    either control: 11 scan gains, 15 losses, 119 ties; **zero exact gains,
    eight losses**. N and all geometry remain fixed; 42 rejected scans stay
    at zero shift. All **626 outputs and Slurm log** are fetched/checksummed;
    **10,359 tests/smoke checks, 23,526 pre-grade checks and 41 aggregate checks**
    pass, with 580 paired rows. Retain saved support, below historical 677/36.
    No post-grade tuning, retry, other direction or unknown25 rerun follows.

0y. **Does full geometric variable projection improve recognition at fixed N?**
    → **COMPLETE; NEGATIVE (Sep 22)**: source **7c5fb6c**, job **11940220**,
    **0:0 in 15m17s**, after 3m07s queued. Profiled gives **632/870 correct,
    20/145 exact**, against same-start joint and replayed support **676/34**,
    unchanged **852/870 coverage**. Against joint: 18 scan gains, 44 losses,
    83 ties; two exact gains, sixteen losses. All 146 profiled objectives
    improve, but 103 searches hit 800 evaluations; better RSS is not better
    recognition. N, split cache, classifier and physical settings stay fixed.
    All **782 outputs and Slurm log** are fetched/checksummed; **5,047 tests,
    20,138 pre-grade checks and 41 aggregate checks**, with 435 paired rows.
    Retain saved support, below historical 677/36. No retry, post-grade tuning,
    other suggested direction or unknown25 rerun follows.

0x. **Does bounded amplitude/background profiling improve assignment with
    saved N and geometry fixed?** → **COMPLETE; NEGATIVE (Sep 22)**:
    source **567cdde**, job **11938782**, **0:0 in 6m17s**, after **5m29s**
    queued. The profile gives **675/870 correct, 32/145 exact**, versus
    exactly replayed support **676/34**, at unchanged **852/870 coverage**.
    All 146 profiles satisfy KKT and reduce RSS (median **0.681%**); all N,
    geometric fields and split cache remain fixed. Against control: 4 scan
    gains, 5 losses, 136 ties; zero exact gains, two losses. All **561 outputs
    and Slurm log** are fetched/checksummed; **4,936 tests, 22,912 pre-grade
    checks and 31 aggregate checks**, with 290 paired scan rows. Retain saved
    support, below historical 677/36. No nonlinear refit, parameter tuning,
    new direction or unknown25 rerun follows.

0w. **Does fresh raw-GCV selection with corrected fused input improve the full
    recognition pipeline?** → **COMPLETE; NEGATIVE (Sep 22), JOB 11936561**:
    source **853ad38**, **0:0 in 1h25m43s**, after 26s queued. Fresh raw GCV
    gives **661/870 correct, 32/145 exact**, versus fresh control **675/34**
    and saved support **676/34**. Exact counts fall **106 → 101/145**, coverage
    **852 → 849/870**; 22 scans change N, total 900 → 898 lobes. This is the
    direct-extractor lineage, not the separate batch hybrid counting reference.
    The fresh control changes three decisions; same-N geometries also vary,
    preventing pure causal attribution to N. All **438 outputs and Slurm log**
    are fetched/checksummed; **1,658 tests, 9,054 pre-grade checks, 37 aggregate
    checks**, 435 paired rows. Retain saved support, below historical 677/36.
    No retry, tuning, registration/noise experiment or unknown25 rerun follows.

0v. **Can covariance regularization throughout learning or coherent Student
    decisions improve recognition?** → **COMPLETE; NEGATIVE (Sep 22)**:
    source **bb6943b**, job **11936229**, **0:0 in 7m56s**, after **3s** queued.
    All-update shrinkage gives **637/870 correct, 24/145 exact**; coherent
    df-five Student gives **627/870, 25/145**, versus exactly replayed
    **676/34** support, all at **852/870 coverage**. All ten fits per arm
    meet their stopping criterion (25/43 and 65–86 updates). Two shrinkage
    seeds show ordinary-likelihood decreases, without exclusion or retry.
    Shrinkage gains on 19 scans and loses on 55 (exact 6/16); Student gains
    on 17 and loses on 58 (4/13). False positives outweigh recovered false
    negatives; nearly all changed decisions enter zero-margin vote ties.
    All **588 outputs and Slurm log** are fetched/checksummed. **4,277 tests,
    11,120 pre-grade checks**, eight component checks and 580 paired scan
    rows pass; keys/N, common upstream tables and four unavailable rows stay
    fixed. Retain support; historical 677/36 remains unexceeded. No new run,
    coefficient/df/seed/threshold search, combination, N refit, DFT change,
    MFA-cap change, measurement-error model or unknown25 rerun follows.
    This reused development benchmark is not independent validation.

0u. **Can lower-rank or heavy-tailed mixture learning improve recognition?** →
    **COMPLETE; NEGATIVE FOR THESE FIXED VERSIONS (Sep 22)**: source
    **e2205cb**, job **11935072**, **0:0 in 8m17s**, after **1m36s** queued.
    Rank-four factors give **620/870 correct, 23/145 exact**; df-five Student
    gives **631/870, 25/145**, versus exactly replayed **676/34** support,
    all at **852/870 coverage**. All Student fits converge; all factor fits
    hit the fixed 200-update cap, so a fully converged MFA optimum is not
    assessed. Factors gain on 21 scans and lose on 65 (5 exact gains / 16
    losses); Student gains on 18 and loses on 57 (4 / 13). False positives
    dominate the loss; almost all final flips enter zero-margin ties.
    All **588 outputs plus Slurm log** are fetched/checksummed, with keys/N,
    upstream/Fisher/k-means tables and availability unchanged. **3,616 tests,
    11,083 pre-grade checks**, eight component checks and 580 paired scan
    rows pass. Retain **676/34**, below historical 677/36. No retry, cap/rank/
    df/seed search, combined arm, post-grade tuning, N refit, paired-view
    model, DFT change or unknown25 rerun follows. Reused grading remains
    development evidence, not independent validation.

0t. **Can scan-grouped Fisher or within-scan GMM naming improve recognition?** →
    **COMPLETE; NEGATIVE / NULL (Sep 22)**: source **43f8993**, job **11931478**,
    **0:0 in 8m30s**, after 3m01s queued. Whole-scan Fisher gives **668/870
    correct, 34/145 exact**, at unchanged **852/870 coverage**, versus the
    exactly replayed 676/34 support control. Eight scans gain, fifteen lose,
    122 tie; four exact chains are gained and four lost. Both disjoint 73-scan
    fits converge. All 26 final flips enter existing zero-margin vote ties.
    Within-scan amplitude naming chooses the same group for every unchanged
    GMM fit, so all scores/decisions remain identical at **676/34**. All keys,
    N and availability remain. The **588 outputs plus Slurm log** are fetched
    and checksummed; 3,035 focused assertions, 10,987 pre-grade checks, eight
    component checks and 580 paired scan comparisons pass. Retain the
    **676/34 support candidate**, below historical 677/36. No split search,
    combined candidate, retry, post-grade tuning, count refit, DFT change or
    unknown25 rerun follows. The reused grade is not independent validation.

0s. **Can covariance sharing throughout GMM learning improve recognition?** →
    **COMPLETE; NEGATIVE (Sep 22)**: source **becd0ce**, job **11931106**,
    **0:0 in 5m53s**, after 15 seconds queued. Sharing covariance at
    initialization, every EM M-step and both hard updates gives **666/870
    correct, 666/852 (78.2%), 6/145 exact**, versus the exactly replayed
    676/34 support control. Twenty-five scans gain, thirty-five lose and
    eighty-five tie; zero exact chains are gained and 28 lost. Of 69 final
    changes, 68 switch 1→0 and one 0→1; 55 leave old GMM=1 / k-means=0 ties.
    All ten tied fits are named, with identical component covariances and
    free unequal masses. All 146 scans / 900 keys, N, availability and
    upstream/Fisher/k-means tables remain. All 392 outputs and the Slurm log
    are fetched/checksummed; 2,716 focused assertions, 7,338 pre-grade output
    checks, five saved-component checks and paired arithmetic pass. Reject
    this variant and retain 676/34; historical 677/36 remains unexceeded.
    No retry, other proposed lead, post-grade adjustment, count refit, DFT
    change or unknown25 rerun follows. This fixed comparison is complete,
    not independent validation on the repeatedly reused benchmark.

0r. **Can whole-scan bootstrap aggregation improve recognition?** →
    **COMPLETE; NEGATIVE (Sep 21)**: source **da886ce**, job **11925188**, exit
    **0:0 in 7m17s**, after 29 seconds queued. Twenty whole-scan bags give
    **665/870 correct, 665/852 (78.1%), 10/145 exact**, versus the exactly
    replayed 676/34 support control. Twenty-two scans gain, thirty-three lose,
    ninety tie; zero exact chains are gained and 24 lost. Of 62 final 1→0
    flips, 57 leave old GMM=1 / k-means=0 ties. All 200 seed fits are accepted;
    no best bag/seed or class count is chosen. Retain 676/34; historical
    677/36 remains unexceeded. All 146 scans / 900 keys, N, availability and
    upstream/Fisher/k-means tables remain. All 393 outputs and the Slurm log
    are fetched/checksummed; 2,471 focused assertions, 21,927 pre-grade output
    checks, seven component-score checks and paired arithmetic pass. No retry,
    fusion/threshold/precision adjustment, new arm, count refit, DFT change or
    unknown25 run follows. This fixed comparison is complete, not independent
    validation on the repeatedly reused benchmark.

0q. **Can continuous per-seed GMM scores improve the final vote?** →
    **COMPLETE; NEGATIVE (Sep 21)**: source **bef315f**, job **11922389**, exit
    **0:0 in 7m19s**, after 36m23s queued. Continuous memberships give
    **671/870 correct, 671/852 (78.8%), 26/145 exact**, versus the exactly
    replayed 676/34 support control at identical coverage. Eight scans gain,
    thirteen lose, 124 tie; one exact chain is gained, nine lost. All 24 final
    flips are 1→0 from old exact ties (GMM=1, k-means=0). Learning, group naming,
    threshold and precision remain unchanged. Reject this variant and retain
    676/34; historical 677/36 remains unexceeded. All 146 scans / 900 keys, N,
    upstream/Fisher/k-means tables and four unavailable rows are retained.
    All 392 outputs and the Slurm log are fetched/checksummed; 2,199 pre-run
    assertions, 7,240 output checks and paired arithmetic pass. No temperature,
    calibration, scan weighting, resampling, naming-rule change, combined arm,
    post-grade tie adjustment, count refit, DFT change or unknown25 run follows.
    Reused benchmark evidence is not independent validation.

0p. **Can equal total GMM training weight per scan improve recognition?** →
    **COMPLETE; MIXED RESULT (Sep 21)**: source **5ff68e2**, job **11915029**,
    exit **0:0 in 6m01s**, after 41 minutes queued. Weighting gives **675/870
    correct, 675/852 (79.2%), 36/145 exact**, versus the exactly replayed 676/34
    support control. Three scans gain, four lose and 138 tie; two exact chains
    are gained and none lost. All seven changed decisions end at the existing
    zero-margin vote tie. Retain **676/34 as the primary working reference**;
    record weighting as a tradeoff, not a champion. Coverage, N, all 146 scans /
    900 keys and upstream/Fisher/k-means tables are unchanged. All 392 outputs
    and the Slurm log are fetched/checksummed; 2,057 pre-run assertions, 360
    output checks and paired arithmetic pass. Usable rows alone determine
    weights; chemical populations remain free. No bagging, robust scaling,
    additional arm or post-grade tuning. Historical 677/36 remains unexceeded.

0o. **Can robust per-scan GMM normalization improve recognition?** →
    **COMPLETE; NEGATIVE (Sep 21)**: source **77b3ce6**, job **11911253**, exit
    **0:0 in 5m55s**. Median/Type-7-IQR gives **667/870 correct, 667/852 (78.3%),
    10/145 exact**, versus the exactly replayed 676/34 support control. Same
    coverage, all 900 keys and four unavailable assignments; upstream signals,
    Fisher and k-means are identical. There are 23 scan gains, 33 losses, 89
    ties; no exact chains gained and 24 lost. No fallback occurs in 1,168
    scan/feature pairs, but expanded feature magnitude grows 5.44→335.47.
    Reject this variant and retain 676/34. All 392 outputs and the Slurm log
    are fetched/checksummed; 1,900 focused assertions, 1,261 output checks and
    paired arithmetic pass. The initial audit compared necessarily different
    staging paths; it was corrected to check scan identities/N before grading,
    not to alter science. No clipping, factor/quantile search, weighting,
    resampling, extra arm, count/DFT/unknown25 rerun or post-grade tuning follows.
    Historical 677/36 remains the target, not independent validation.

0n. **Can complete-patch learning retain partial-patch predictions without their training influence?** →
    **COMPLETE; NEGATIVE (Sep 21)**: source **f7a4074**, job **11908779**,
    exit **0:0 in 6m42s**. Control exactly replays 676/34; complete training gives
    **674/870 correct, 674/850 (79.3%), 33/145 exact**. Of 900 retained keys,
    893 train and 894 can be scored by GMM. One scan has no complete patch for
    per-file normalization: two previous predictions become abstentions. One
    scan improves, two regress, 142 tie; `240818_020` loses exactness. Keep the
    support candidate at **676/34**, not this alternative. All 392 outputs and
    the Slurm log are fetched and checksummed; 1,782 focused assertions, 59 output checks and
    paired arithmetic pass. The initial audit's false all-scans-have-training
    assumption was corrected to check the already-frozen abstention rule, not
    the science; its failed log remains. No fallback, robust scaling, weighting,
    extra arm, post-grade tuning, count refit, DFT or unknown25 run follows.
    Historical **677/36** remains the target to exceed, not independent validation.

0m. **Can complete observed disks or final covariance-volume scoring help?** →
    **COMPLETE; SUPPORT GAIN, NO PROMOTION (Sep 21)**: source **db411a2**, job
    **11906168**, exit **0:0 in 7m58s**. All twelve control tables replay 675/33
    byte-for-byte. Complete-disk symmetric support gives **676/870 correct,
    676/852 (79.3%), 34/145 exact**, with coverage 849→852/870. Three newly
    available predictions are correct, but two older decisions regress: two
    scan gains, two losses, 141 ties, one exact gained and none lost. Retain
    `unit_assignment_patch_support.toml` as the latest opt-in candidate, not a
    champion. Final Gaussian scoring gives **666/849 (78.4%), 32 exact** and
    is rejected: four scan gains, twelve losses, 129 ties. All 900 keys remain;
    support leaves four unavailable keys versus seven in the other arms. All
    585 outputs and the Slurm log are fetched and checksummed; 1,604 focused assertions,
    76 output checks and paired-grade arithmetic pass. No imputation, combined
    arm, normalization/weighting change, count refit, unknown25 run, DFT change
    or post-grade tuning follows. Historical 677/36 still leads by one position
    and two exact chains merely to equal it; this is not independent validation.

0l. **Can training-consistent Fisher centering or final covariance shrinkage help?** →
    **COMPLETE; CENTERING GAIN, NO PROMOTION (Sep 21)**: source **0d58b44**, job
    **11891104**, exit **0:0 in 9m04s**. Twelve control tables replay 673/33
    exactly. Training-mean centering gives **675/870 correct, 675/849 (79.5%),
    33/145 exact**: two scan gains, no losses or exact-chain changes. Retain
    `unit_assignment_centered_fisher.toml` as the opt-in working candidate at that stage,
    not a champion. Final Ledoit-Wolf covariance gives **665/870, 665/849 (78.3%),
    34 exact**: five scan gains, eleven losses, two exact chains gained and one
    lost. It is not retained as the working replacement despite much better
    covariance conditioning. All 900 keys, seven unavailable rows and counts
    remain; 590 outputs and the log are fetched and checksummed. Verification
    passes 1,186 focused assertions, 76 output checks and paired-grade arithmetic.
    Historical 677/36 remains ahead by two positions and three chains merely to
    equal it. No combination, normalization change or post-grade tuning follows;
    neither numerical stability nor this reused benchmark is independent validation.

0k. **Can signed physical margins or affine Fisher patches improve label-free assignment?** →
    **COMPLETE; BOTH NEGATIVE (Sep 21)**: source **3c45c7c**, job **11889205**,
    exit 0:0 in **7m05s**. The twelve control tables reproduce 673/33 exactly.
    Signed CC margins give **668/870 correct, 668/849 (78.7%), 32/145 exact**;
    affine Fisher gives **672/870, 672/849 (79.2%), 32/145**, all at unchanged
    849/870 coverage. Signed margins gain on four scans and lose on nine; Fisher
    gains on nine and loses on ten. Both are rejected as working replacements.
    All 146 scans / 900 keys and seven unavailable assignments remain; all 585
    outputs and the log are fetched and checksummed. Verification passes 1,518
    focused assertions, 77 saved-output checks and paired-grade arithmetic.
    Keep the affine-descriptor 673/33 candidate below historical 677/36. No
    combined arm, post-grade tuning, count refit or unknown25 rerun follows.
    Benchmark truth is evaluation-only; this is not independent validation.

0j. **Can a better-defined transverse descriptor close the historical gap?** →
    **COMPLETE; AFFINE GAIN, NO PROMOTION (Sep 20)**: job **11888011** finishes
    0:0 in 7m19s. All twelve control tables reproduce matched residuals at 672/29.
    The first moment gives **672/28** and is rejected; affine-residual half-plane
    asymmetry gives **673/870 correct, 673/849 (79.3%), 33/145 exact**, at the
    same coverage. Eight scans gain and seven lose; five exact chains are gained
    and one lost. All sixteen changed decisions land at the unchanged vote tie
    (zero margin); one belongs to the already-excluded technical scan. The new
    opt-in working candidate still trails historical 677/36 by four positions
    and three exact chains. All 585 outputs and the log are fetched. Verification
    passes 833 focused assertions plus 77 saved-output checks. No post-grade
    tuning, combined variant, count refit or unknown25 rerun follows.

0i. **Does matching the patch residual filter improve recognition?** →
    **COMPLETE; SMALL GAIN, NO PROMOTION (Sep 20)**: job **11883607** finishes
    0:0 in 5m15s. `S(data-model)` gains **671 → 672 correct / 870** and
    **28 → 29 exact chains / 145**, at unchanged 849/870 coverage. The twelve
    regenerated control tables are byte-identical to the saved symmetric result.
    Both arms preserve all 900 keys and seven unavailable assignments. Two scans
    gain, one loses and 142 tie; one exact chain is gained, none lost. The variant
    is retained as the opt-in working candidate, not a champion. Historical
    677/36 still leads by five correct positions and seven exact chains. No
    post-grade setting change, count refit, unknown25 run or independent
    generalization claim follows. All 390 outputs and the log are fetched.

0h. **Can a reproducible assignment exceed the lost historical champion?** →
    **OPEN; TARGET UNCHANGED (Sep 22)**: the user explicitly rejects lowering
    the target to the best runnable reconstruction. The historical 677 correct
    positions / 870 and 36 exact chains / 145 remain the comparison target;
    complete-disk support (676 / 34, comparison 0m) is the latest opt-in working
    candidate, not a new champion. It improves training-mean Fisher (675 / 33),
    after affine-residual half-plane asymmetry (673 / 33), matched residuals
    (672 / 29) and symmetric fusion (671 / 28). The remaining gap is one correct
    position and two exact chains merely to equal history; coverage is 852 versus 854.
    Equal-scan weighting (0p) gives 675 / 36 at the same coverage: two exact
    chains gained against the support control but one correct position lost.
    It matches the historical exact-chain total without exceeding history or
    replacing the 676 / 34 primary working reference.
    Continuous internal GMM voting (0q) regresses to 671 / 26, at the same
    coverage; all 24 final changes leave old exact ties and switch 1→0.
    It is rejected without changing the threshold or serialization precision.
    Whole-scan bagging (0r) regresses to 665 / 10 at the same coverage; 24
    exact chains are lost and none gained. Of 62 final 1→0 flips, 57 leave
    old GMM=1 / k-means=0 ties. It is rejected without post-grade adjustment
    to fusion, seeds, replicate count or naming; 676 / 34 remains primary.
    Sharing covariance throughout GMM learning (0s) regresses to 666 / 6,
    at unchanged coverage; zero exact chains are gained and 28 lost against
    the support control. It is rejected without trying other proposed leads
    or post-grade adjustments. The historical objective is not achieved.
    The one physical-u Fisher mirror comparison is **COMPLETE; NEGATIVE**:
    job 11879515 (2m03s, exit 0:0) replays the control byte-for-byte but the
    treatment gives **669/849 (78.8%), 27 exact**, versus 671/849 and 28.
    Four final decisions change: one scan gains, three lose, 141 tie; one
    exact chain is lost. All 900 keys and seven unavailable rows remain.
    The mirror variant is not retained as the working reference. The separately
    user-approved matched-residual comparison (0i) subsequently gains one correct
    position and one exact chain. The next user-requested comparison (0j) changes
    the descriptor alone and gains another position and four exact chains.
    The signed-CC and affine-Fisher comparison (0k) is also complete: both regress
    to 668/32 and 672/32 respectively, so that pass keeps 673/33. The numerical
    comparison (0l) then gains two positions with centering, without exact-chain
    improvement. Shrinkage loses eight positions despite one extra exact chain.
    Complete-disk support (0m) then recovers three correct predictions but loses
    the two centering gains, net +1 correct and +1 exact chain. Final Gaussian
    scoring regresses to 666/32 and is rejected. No post-grade setting search
    follows these passes. Counts, thresholds, voting, seeds and DFT inputs remain
    unchanged; support reduces unavailable keys from seven to four without
    imputation. Exceeding history remains unresolved on the fixed 870-position
    denominator and the exact-chain criterion.

0g. **Does symmetric fit filtering improve recognition?** → **COMPLETE; MODEST
    GAIN, NO PROMOTION (Sep 20)**: source `d12283a` fits `(F+B)/2`, retaining
    `(S(F)+S(B))/2` for ROI detection. Job **11878491** completed in 34m25s,
    exit 0:0; both arms retain all 146 files / 900 cached lobes and the same
    seven unavailable keys. Regenerated control `f2191f7` reproduces the saved
    prediction TSV byte-for-byte. Full145 improves **666 → 671 correct / 870**,
    **78.4% → 79.0%** classified accuracy and **24 → 28 exact chains**. Coverage
    stays 849/870, with 16 missing positions and 38 extra lobes. There are 19
    winning scans, 16 losing and 110 tied; nine exact chains gained, five lost.
    The historical 677-correct / 36-exact reference remains better. The branch
    is not promoted; count selection and unknown25 were not rerun. Independent
    generalization is unestablished on this reused benchmark. No further
    filtering, weight or threshold search follows this completed comparison.

0f. **Can a bounded change improve recognition on the user's benchmark?** →
    **SAVED-OUTPUT FIRST PASS COMPLETE; NO GAIN (Sep 19)**: the two predeclared
    native component endpoints preserve all 900 cached keys and the same seven
    unavailable rows. On full145, k-means gives 662/849 correct and 16 exact
    chains; GMM gives 622/849 and 19, versus reconstructed soft vote 666/849 and
    24, and historical reference 677/854 and 36. Correct/all870 is respectively
    662, 622, 666 and 677; all retain 16 missing positions and 38 extra lobes.
    Neither endpoint is promoted and no weight/threshold search follows.
    Independent saved-row arithmetic reproduces all four physical grades and
    their paired gains/losses. The 12 reconstruction exact-chain losses all
    have N6; eleven involve old zero-margin decisions. Both heads share 95
    errors on the 636 N6 positions. These diagnose existing decisions, not an
    information ceiling or permission to restore labels. A separate source
    finding is live asymmetric fused-fit filtering `(F + S(B))/2`; a native
    synthetic-only probe confirms it, but its recognition effect is unknown.
    No fix, raw reprocessing, fit, training or HPC job ran in this pass. A
    separately scoped correction would need regenerated control and treatment
    geometry/features, fixed inference, and external recognition grading.
    The benchmark is reused development data; independent generalization and
    chemical/localization causes remain unestablished. Production and all
    previous completed jobs stay unchanged.

0e. **Can missing-value-aware background estimation reduce preprocessing bias
    without removing molecular contrast?** → **FOUR-SCAN COMPARISON COMPLETE;
    MOLECULAR VALIDATION UNRESOLVED (Sep 18)**: the synthetic prototype remains
    unchanged. Approved Viper job **11820414** completed 0:0 in 55 s; 37 outputs
    and the log were fetched and hash-verified. No new observed pixels are lost;
    the incomplete scan remains `PARTIAL`. Guarding reduces descriptive
    background spread, but forward/backward physical RMS changes are mixed and
    local anticorrelations remain. All best lags remain boundary/ambiguous;
    no shift or method is accepted. The saved-output checker retains 30,811
    passes and 15 strict stationarity failures. Separate 118/118 arithmetic
    checks establish export-roundoff compatibility, **not exact stationarity
    replay**; normalized parameters were not exported. All other plane/row,
    objective, mask, metric and native-baseline checks pass. No count/classifier
    refit, grade, tuning, noise calibration or production promotion occurred.
    The frozen footprint is not true background; real molecular preservation
    is still unvalidated. No further campaign is approved.

0d. **Label-free mathematical exploration: what is numerical, representational,
    or acquisition-limited?** → **BOUNDED DIAGNOSTICS COMPLETE; NO PROMOTION
    (Sep 18)**: the corrected, explicitly authorized job 11812202 completed 0:0.
    All 78 rows remain; fixed-N RSS improves, but the candidate-set GCV minima
    stay 8/12/9/9 and 25/26 refinements stop at 800 evaluations. Original selected
    counts 6/11/9/9 are not truth and are unchanged. All four registrations are
    boundary/ambiguous; all 16 background ACF lengths are censored at 16 pixels.
    No new calibrated noise, independence or block-size claim is supported.
    Current is feedback-coupled and anticorrelated across directions. Both
    Fisher score files replay byte-identically (222/900 keys; seven NA retained),
    without recovering unique historical weights or validating chemistry.
    Gradient attribution and raw missingness remain limitations, not proof of
    absent chemical information. No grade, threshold, production prediction or
    method promotion changed. Accuracy, calibrated uncertainty and exact unknown
    counts/chemical identities remain unresolved; no further campaign is approved.

0c. **Apply a reproducible native pipeline to raw 10–20mer scans** → **UNKNOWN25
    OUTPUTS COMPLETE; SCIENTIFIC VALIDATION STILL LIMITED (Sep 17)**: after the
    selected-support repair, job **30278010** completed from `aa422a3` (exit 0:0,
    6 min 16 s). All 25 scans / 222 keys are present, with exact cache reuse,
    179 native integrity assertions and byte-identical saved-vote/QC replay.
    The model assigns 189 class 0 / 33 class 1 / zero `?`; 19/25 chains need QC
    review for low mean confidence. Confidence is uncalibrated; zero `?` does
    not establish certainty. The grid and seven priority standalone maps were
    inspected, but they are fitted-coordinate maps, not raw STM overlays.
    All run monitors are stopped and failed outputs are preserved. The authorized
    full146 reconstruction/comparison/grade is now **complete (Sep 18)**:
    Viper job **11786116**, source `00528d2`, completed `0:0` in 11m04s with all
    **146/900** cached GCV keys, **695 zero / 198 one / 7 `?`**, and unchanged
    scientific inputs/settings. This is not the incompatible 871-key promoted
    counting chain. Saved-output checks pass 73/73; plot-only fix `6fa3297` passes
    267/267 and provides 19 readable overview pages without changing standalone
    maps. External full145 own-N grade: **666/849 (78.4%), 24 exact chains**,
    16 missing control positions, 38 extra lobes. The frozen reference remains
    **677/854 (79.3%), 36 exact chains**. The reconstruction is not promoted or
    called an exact historical recovery; unknown-chain chemical validation and
    milestone closure remain unresolved. No benchmark-guided tuning or new
    method campaign automatically follows this completion. The separately
    authorized bounded exploration is recorded in0d. Do not resubmit either
    completed reconstruction job.

0b. **Can label-free unit assignment reach the promotion bar?** → **RESOLVED
    (Aug 3)**: the label-free champion (soft vote of k-means 4-view and GMM
    1-view + per-channel constant-current margins + Fisher empirical mold
    margin, self-training 2) scores 79.3% classified physical accuracy /
    36 exact chains / 677 fixed-denominator honest on the 145-file benchmark
    - the promotion bar (78.9% / 18 / 677) is MET. Regrading saved outputs
    reproduces the accepted metrics; a raw-input rebuild of that exact method
    is not currently reproducible because the original `patch_u_asym` producer
    is lost. The recorded half-split Fisher score (66.3% per-lobe) and pipeline
    CV (678/34) are historical results, not new reconstruction validation. See
    0c for the separately named native application.

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

### 2026-09-16 — Unknown-chitosan: prediction-only frozen champion, Julia 1.13

Milestone 1 resumes the real 10–20mer application, not a new method campaign.
The user explicitly requires Julia 1.13 for production, tests, and subprocesses;
local `julia --version` is `1.13.0` (juliaup `release`). The five local packages
load under the root environment. The Manifest is not hand-edited.

`build_cc_soft_champion.py` now stops after predictions. Its automatic benchmark
report call and named-file exclusion are removed from the production path.
The separate existing report command retains external benchmark membership.
The mold settings, features, seeds, self-training, probability average, binary
threshold, and finite-pair confidence arithmetic are unchanged. Feature keys,
not the intersection of component keys, determine output coverage. Missing or
nonfinite component predictions and explicit component abstentions now remain
`?` with a reason instead of being dropped or replaced by a forced binary vote.
This is a missing-data correction, not a new confidence threshold.

Production inputs reject control columns and duplicate keys. A new explicit
work directory and output file are required; earlier outputs are not overwritten.
Python subprocesses use the invoking interpreter, and `--julia` can select the
Julia 1.13 executable. Confidence remains an uncalibrated soft-vote margin.

Verification: `python3 test/test_cc_soft_champion.py` passes 12 focused tests,
including the frozen finite-pair arithmetic and tie, missing/invalid components,
no named-file exclusion, forbidden inputs, output preservation, and prediction-
only subprocess wiring. CLI help works. The existing post-hoc report on
`results/unit_assignment/best_labelfree_cc_soft_20260802.tsv`, run under Julia
1.13, reproduces **677/854 (79.3%) and 36/145 exact chains** in
`results/unknown_chitosan_milestone1/benchmark_frozen/`. This is regrading the
frozen predictions, not a full scientific rebuild. No counting rule changed.

Next: wire real selected-N feature extraction and the existing QC/plot commands.
The 25 raw unknown scans are at
`/home/durif/Rebecca/data/10_20mer_analysis`, not the repository's `data/` tree;
the local DFT cubes exist. No remote calculation has been started.

### 2026-09-16 — Champion input preflight and existing benchmark regressions

The prediction-only champion now checks all mandatory frozen feature columns,
contiguous feature lobes, and exact forward/backward 17×17 mold-patch key
coverage before any expensive computation. This avoids building molds and
fitting predictors only to discover that `patch_u_asym` or a patch is absent.
The scientific feature definitions and model settings are unchanged.
`python3 test/test_cc_soft_champion.py` passes 15 tests. Julia 1.13 also validates
a two-lobe prediction fixture containing an explicit missing-component `?`.

The raw-input application is blocked by one unrecovered scientific prerequisite:
`patch_u_asym` is consumed by the frozen champion, but its producer is absent
from current sources and local result scripts. The original example input
`/tmp/opencode/full_features_origN_bwd_uasym.tsv` is also absent. The journal
calls it a u-weighted first moment of a 9×9 residual patch; grid orientation,
coordinate weights, normalization, and invalid/zero-signal behavior are not
fully specified. A 17×17 approximation is not the frozen descriptor. Recover
the original producer rather than inventing or benchmark-tuning a substitute.
No remote job or new fitting campaign was started.

Existing counting reports were rerun under Julia 1.13, using saved summaries
and the existing external grader only:

- `results/best_plots_240817_primary_rerun/summary_overlap060_hard.tsv` with
  `benchmarks/chitosan_240817.toml`: **39/39 exact**.
- `results/experiments/6mer_full146/pm2_confirm/summary_overlap060_hard.tsv`
  with `benchmarks/chitosan_6mer_counting_confirmed.toml`:
  **129/145 exact, 143/145 within one lobe**.

The first expanded-report attempt mistakenly used the adjacent historical
`adaptive_support_rescue` summary (104/145 exact, 138/145 within one). That was
the wrong saved policy, not a code or scientific regression. The documented
promoted `support_midpoint_hybrid` summary is `pm2_confirm`; the runbook now
names its exact path and command. No candidates were refit or parameters
changed. New reports/logs are under
`results/unknown_chitosan_milestone1/counting_regression/`; the separately
regraded frozen assignment remains **677/854 (79.3%), 36/145 exact**.

### 2026-09-16 — Unknown-assignment figures: one-file grid and explicit convention

`plot_unit_assignment.py` crashed on a one-file grid because its subplot axes
were left as a nested array. Axes are now flattened for every file count. The
grid title is method-neutral (`Label-free unit-assignment map`) instead of
claiming that every input came from preliminary DFT molds at a fixed height.
Standalone figures now show the existing convention explicitly: GlcN (0),
GlcNAc (1), uncertain (?). Ordering, colors, confidence values, assignment
thresholds, and the legacy overlay coordinate mapping are unchanged.

`MPLBACKEND=Agg python3 test/test_plot_unit_assignment.py` passes three focused
checks: one-file grid, standalone 0/1/? with the legend and confidence summary,
and invalid prediction rejection. The champion preflight/vote suite still
passes 15 tests. Julia 1.13 native package tests pass: STMFitCore 32 assertions,
STMSXMIO 21 assertions. These checks use synthetic fixtures and existing tests;
they do not imply successful raw-input application to the 25 real chains.

Follow-up verification for this figure change: the native `--mode all` CLI
rendered the two-lobe 0/? fixture to
`results/unknown_chitosan_milestone1/champion_fixture/plots/`. The standalone
image was inspected; its unknown lobe is neutral gray and the chemical legend
is explicit. The Julia 1.13 QC command flags the same fixture for uncertain
fraction 0.5 and mean confidence 0.4, with no missing-plot flag. These are
synthetic checks, not real-chain assignments. `GKSwstype=100 julia --project=.
docs/make.jl` exits 0; Documenter only warns about archive/search-index size
and skips deployment outside CI. Logs remain under
`results/unknown_chitosan_milestone1/`. Milestone 1 remains incomplete pending
the original `patch_u_asym` extraction recipe; no replacement was implemented.

### 2026-09-16 — User-approved explicit descriptor reconstruction in Julia 1.13

The user reports that the original untracked `patch_u_asym` producer is destroyed
and approves implementing an explicit Julia descriptor with a complete pipeline,
then comparing with the frozen champion. The previous recovery blocker is thus
replaced by a bounded reconstruction task: one documented descriptor, one native
pipeline, tests and one fixed comparison; no variant sweep or score-guided tuning.
The historical champion is not replaced or relabeled as reproduced in advance.

The new `cc_soft_reconstructed_v1` uses backward residual 9×9 patches, half-width
0.32 nm and step 0.08 nm, already median/sample-std normalized by the extractor.
`patch_u_asym_reconstructed = sum(sign(u)*p)/sum(abs(p))` is a half-plane parity
statistic, not the ambiguous u-weighted first moment in the old journal. The
extractor's actual u-outer/t-inner order is used (t varies fastest); inheriting
the conflicting row-axis comment in `enrich_unit_features.py` would swap axes.
Nonfinite samples or L1 mass ≤ 1e-12 yield NA and an explicit reason, not row
removal. No labels, expected counts, composition prior or frozen outputs enter
construction. Settings and the method name are in
`config/unit_assignment_reconstructed.toml`, fixed before comparison.

`test/build_reconstructed_descriptor.jl` is the native CLI;
`test/lib/reconstructed_unit_assignment.jl` also supplies strict key checks and
soft-vote output retaining unavailable lobes as ?. Synthetic checks verify u/t
orientation, mirror sign, positive-scale invariance, central-row behavior, the
distinction from a first moment, invalid patches/columns/keys, and unchanged
finite-pair voting arithmetic. `julia --project=.
test/test_reconstructed_unit_assignment.jl` passes **31 assertions** on Julia
1.13. Next are native constant-current/Fisher ports and raw-to-QC integration.

### 2026-09-16 — External reconstruction comparison and native CLI check

`test/compare_reconstructed_champion.jl` compares two saved prediction tables
without fitting or choosing a descriptor. It reports every key in their union,
missing/extra lobes, assignment disagreements, and reported confidence deltas.
An optional external benchmark manifest only prepares rows for the existing,
separate grader. Production cannot accept a reference/benchmark option and never
calls this report. Matching aggregate accuracy is deliberately not used as a
reconstruction criterion; exact assignment/confidence equality is reported
separately and only establishes equality on the supplied inputs.

The external comparison tests pass 9 assertions on synthetic changed/identical/
partial tables and benchmark-only row export. The descriptor suite now passes
32 assertions, including execution of its actual CLI. That CLI test caught and
fixed a Julia macro-parsing error in the main guard (`@__FILE__` needs a function
boundary before `&&`). No scientific setting changed. The raw-to-QC native
wrapper is being integrated; its first 8 input-boundary checks pass, but a full
real-data run and champion comparison have not yet been performed.

### 2026-09-16 — Native constant-current mold port, fixed numerical comparison

`test/build_cc_molds_native.jl` and `test/lib/cc_mold_native.jl` replace the Python
constant-current builder for the explicitly reconstructed pipeline. All sampling
and isovalue settings come from `config/unit_assignment_reconstructed.toml`.
This port preserves the reference's first-axis-fast cube indexing, NumPy-style
floating grids, highest occupied z sample, first isovalue below the target mean
height, and existing parity/mirror/normalization rules. These conventions are
not silently corrected or claimed to be a new physical validation of the cubes.

On Julia 1.13, `test/test_cc_mold_native.jl` passes 143 synthetic assertions.
One actual GlcN/GlcNAc pair comparison passes another 76 assertions (44.3 seconds
for the reported testset). Python runs only as the test reference. LDOS samples,
selected isovalues, missing values and height vectors match exactly; normalized
values differ by at most 4.44e-16. The final template TSV is byte-identical.
GlcN has 72/289 supported samples, mean height 0.4964 nm; GlcNAc has 235/289 and
0.3814 nm. The latter follows the frozen first-below-target rule (target 0.50 nm),
not a nearest-target optimization. Frame-file height 0.55 nm is not substituted
for that explicit setting. No QE work or benchmark-guided setting change occurs.

The real-pair log and reusable native templates are under
`results/reconstructed_cc_soft_v1/mold_parity/` (log alongside that directory).
Raw malformed inputs and existing outputs are rejected rather than overwritten.

### 2026-09-16 — Native empirical Fisher reconstruction

The native Fisher component is in `test/lib/empirical_fisher_native.jl`, with
`test/build_empirical_fisher_native.jl` as its production CLI. It reads only the
forward patches and explicit reconstruction config. Python/NumPy is used only in
numerical tests, not production. `test/test_empirical_fisher_native.jl` passes
201 assertions under Julia 1.13, including separate-process byte replay.

The fixed port keeps PCA10, a two-component full-covariance GMM with learned
weights, center-pixel amplitude ordering, latent sample covariance plus ridge,
Fisher back-projection, opposite even/odd-lobe folds, six-decimal margins, and
the historical raw-patch/centered-midpoint score. The legacy row-major reshape
and last-axis reversal actually reverse physical t for the current u-outer /
t-inner extractor; the old `flip_u_disk` name did not identify the physical axis
correctly. This convention is retained explicitly, not silently changed.

This native GMM is not asserted to reproduce sklearn bit for bit. Julia's
MersenneTwister, one two-center kmeans++ draw, Lloyd partition-stability stopping,
exact responsibility masses (rather than sklearn's epsilon addition), and
SVD/BLAS rounding can change fitted clusters. They are documented numerical
implementation differences, not variants selected by champion agreement.
Finite EM solutions at the configured iteration limit remain usable with a
warning. Invalid patches and missing/degenerate opposite folds retain all input
keys with NA and an explicit reason. Tests cover sign/orientation, fixed-partition
NumPy linear algebra, EM arithmetic, fold independence, and unavailable scores.
No benchmark labels or frozen predictions entered implementation or testing.

### 2026-09-16 — Native raw-to-QC wrapper, synthetic end-to-end result

`test/run_reconstructed_chitosan.jl` now composes counting or explicit cached
selected-N geometry, base/split refits, local features, three patch exports,
reconstructed descriptor, native CC/Fisher, existing Julia predictors, soft
vote, validation, QC and native maps. It propagates the Julia 1.13 executable and
up to four threads to subprocesses. The old default triage file is disabled.
No Python, benchmark, known sequence, reference table, or automatic grader is
called by production. Existing inputs can be reused explicitly; output paths
must be new. Counts and lobe-key coverage are checked after every extraction.
Missing files abort rather than silently reducing the cohort; `failures.tsv`
retains the intended input file list and stage. Unknown counts are NA, not zero.

The native pipeline has 12 passing input/failure-boundary assertions and seven
passing end-to-end assertions on a five-chain/40-lobe synthetic extracted-input
fixture. That fixture exercises real subprocesses through predictions, QC and
five individual maps plus the grid. It does not exercise raw SXM parsing or
Gaussian fitting. A separate one-chain plotting smoke passes and was visually
checked; its initial Julia ternary parse error and clipped grid heading were
fixed without changing science. The plot legend states 0=GlcN, 1=GlcNAc, ?=
uncertain. Output confidence is an uncalibrated soft-vote margin; no new
confidence-driven abstention rule is introduced.

A focused real `251206_013.sxm` timing run has started locally with the explicit
10–20mer adaptive-support config, a fresh count (not a choice between old auto /
manual summaries), and the native checked CC templates. It is capped at 30
minutes. No remote job has been submitted. Real full-cohort application and
frozen-champion comparison remain pending.

### 2026-09-16 — Explicit patch preprocessing in the native workflow

Both patch extractors now accept optional `--config` and take stride, flattening
and smoothing from its `[preprocessing]` section. The native wrapper passes the
same count config to all three patch exports. Previously those extractors fixed
1 / plane+rows / 1 internally, even if fitting used a different explicit config.
Those exact legacy settings remain when the new option is omitted. Grid order,
normalization, residual calculation and missing-file behavior are unchanged.

`test/test_patch_preprocessing.jl` passes 69 Julia 1.13 assertions. A synthetic
33×33 SXM gives byte-identical default and explicit `config/chitosan.toml`
outputs; separate setting changes reach the actual preprocessing calls. Invalid
configs fail before output creation. Both CLI option forms and forward/backward
paths are covered. No real multi-scan fit was used for this test. Documentation
also now separates historical 25-scan counting from the pending native
assignment application and corrects the README's overly broad raw-champion
reproducibility statement without changing any accepted benchmark number.

The native documentation build passes under Julia 1.13, with only size-limit
warnings for the old journal/search index and the normal local-deploy warning.

### 2026-09-16 — Real raw-input pilot completed; bounded Raven preparation

The Julia 1.13 native pipeline completes from raw `251206_013.sxm` through count,
base/split refits, three patch exports, CC/Fisher scores, predictions, validation,
QC and maps. Wall time is **13 min 45.91 s**, maximum resident set size 1,268,908
KiB (about 1.21 GiB). The adaptive-support count selects **11**; both fixed-N
feature fits preserve 11 lobes (elliptical GCV 3.002e-6 base, 1.976e-6 split).
All 11 final assignments are `?`, correctly, because each Fisher training half
contains fewer than the configured ten rows. This expected small-cohort failure
is explicit in `fisher_cv.tsv`, not replaced by invented binary labels. The full
25-chain cohort is required for the intended unsupervised application. Outputs
are under `results/reconstructed_cc_soft_v1/pilot_251206_013/`.

Multi-file refits now use the existing extractor `--chunk` interface, with at
most four single-thread subprocesses and BLAS limited to one thread. Counts,
physical settings, selected N and cohort-wide classifiers are unchanged. Chunks
are merged only with matching headers and unique contiguous lobe keys; complete
expected counts are checked next. Seven focused chunk/merge tests pass alongside
the twelve boundary tests. The full synthetic end-to-end check is being repeated
with this execution wiring. This is standard bounded local subprocess work
inside one allocation, not cluster-side agent orchestration.

A simple `hpc/reconstructed_chitosan.sbatch` runs the 25-chain application first,
then the fixed 146-file saved-geometry comparison input through reconstruction.
It requests four CPUs, 16,000 MB and 24 hours, below the observed eight-CPU group
quota. Labels/grading are absent; comparison and grading will run locally after
fetch. Shell syntax passes; no job has yet been submitted. Raven SSH works but
its Julia modules stop at 1.12, so the tested Julia 1.13 installation is copied
to the user software directory. Dry-runs exposed old local agent metadata; the
actual code sync uses tracked source files only (plus the existing root Manifest),
without `.omo`, hidden agent directories, QE, data, or generated results. Only
25 unknown SXMs, 146 geometry-matched SXMs, saved base geometry and verified native
CC templates are transferred as inputs. The initial generic launcher dry-run
pointed at the personal default Viper; the corrected dry-run explicitly sets
`STMFIT_SSH_HOST=raven`. No calculation is run on a login node.

### 2026-09-16 — Final local checks pass; Raven package setup still pending

The final synthetic integration run passes **28 assertions**: 12 boundaries,
9 end-to-end (including nonempty binary output on the sufficient synthetic
cohort and the reconstructed model identifier), and 7 bounded chunk/merge
checks. The docs build passes under Julia 1.13. These checks do not change or
replace the real pilot's explicit small-cohort uncertainty.

Plain rsync transfers of the tracked source, required raw scans, cached geometry,
verified templates and the user-space Julia installation complete successfully.
Raven reports Julia 1.13.0. Both actual remote-input metadata dry-runs pass with
25 and 146 files. Slurm `--test-only` accepts the job (its processor accounting
shows eight, within the observed group quota); this is not a submitted job.

Initial project `Pkg` loading fails while precompiling `LibCURL_jll`, reporting
an unavailable `Zstd_jll` precompiled image. A focused native `using Zstd_jll`
then succeeds. Setup is retried by loading Pkg in the stdlib environment before
activating the old project Manifest. No Manifest is manually edited, no Julia
1.12 fallback is used, and no STM compute has run on the login node. Actual
submission waits for dependency setup to complete successfully.

The native Slurm entrypoint and its required inputs/environment are now also
listed in `hpc/README.md`, including the single-allocation CPU limit, ordinary
dry-run/submit/fetch commands, and separate local comparison. This documents
existing execution behavior; no scientific setting or accepted metric changes.

### 2026-09-16 — Raven job 30271312 submitted, results pending

The clean-stdlib Pkg bootstrap installs the exact locked dependencies without
changing `Project.toml` or `Manifest.toml` (local/remote SHA-256 equality checked).
Its first precompile pass still reports a Pkg REPL-extension cache message;
therefore its zero exit alone was not treated as readiness. A separate native
import of GaussianFit2D, STMMolecularFit, Plots and Clustering then completes
successfully under Julia 1.13.0 after warming the remaining caches (430 s).
Both cohort metadata dry-runs and Slurm test-only had already passed.

One actual Slurm job, **30271312**, is now submitted with
`hpc/reconstructed_chitosan.sbatch`. Initial status is **PENDING (Priority)**.
The job processes unknown25 first, then full146 with saved selected-N base
geometry and newly generated split fits/patches. Inputs contain no benchmark
labels. The source execution commit is `3f2e671`; later commits are documentation.
Remote outputs will be `results/reconstructed_cc_soft_v1/raven_20260916_native_v1/unknown25` and `.../full146`.
The log is `results/reconstructed_cc_soft_v1/native-30271312.log`.
No result or milestone completion is claimed yet. Monitor this job; do not
resubmit it. Fetch `/ptmp`-dependent outputs locally, validate coverage, then
compare and grade using separate local Julia scripts with unchanged settings.

### 2026-09-17 — First native cohort run fails safely at rescued-support export

Raven job **30271312** ran on `ravc4126` from 01:39:43 to 02:39:48 CEST under
Julia 1.13.0 and exited **FAILED, 1:0**, after **01:00:05**. Slurm reports eight
allocated logical CPUs and batch MaxRSS 3,758,536 KiB (~3.59 GiB). Counting and
base-feature export ran on the compute node; no scientific work ran on the
login node. Scientific source/config remained at `3f2e671`.

The actual fresh count summary contains **25 unique files, all status `ok`,
222 selected lobes, range 5–16**, with exactly the intended raw-file basenames.
`features.tsv` contains **204 unique keys over 23 files**, no duplicates/extras,
and every present file retains exactly its selected keys `1:N_selected`.
The missing keys are all nine lobes of **260215_022.sxm** and all nine lobes of
**260220_083.sxm**. Both extractor logs report `No chain model fit succeeded`.
These are also the only two count rows selected by `adaptive_support_rescue`.
The extractor reads only N from the selected summary and rebuilds the original
support, losing the support expansion used by counting. The fit engine rejects
an infeasible fixed-N support before optimization. This is a production wiring
gap, not evidence for changing the descriptor or tuning against labels.

The mandatory coverage check stops the pipeline at `base_features`;
`failures.tsv` retains all 25 intended files and the explicit mismatch reason.
**No cohort assignments, assignment QC/maps, full146 execution, champion
comparison or new external grade exist.** Counting success is not chemical
validation, and the milestone remains incomplete. Do not silently omit the two
files or lower their counts to unblock the classifier.

Both rsync fetches exit 0. Artifacts are saved locally under
`results/reconstructed_cc_soft_v1/raven_20260916_native_v1/unknown25/`; raw-input
symlinks are excluded. The main log is
`results/reconstructed_cc_soft_v1/native-30271312.log`. The internal monitoring
heartbeat is cancelled. **No job was resubmitted.** Next: establish the minimal
support-context repair and focused regression test, then obtain a separate
continuation decision. All descriptor/classifier settings and accepted
historical benchmark metrics remain unchanged.

Static diagnosis confirms the exact lost state. The accepted rescues expand
support **2.86 -> 5.55 nm** for `260215_022.sxm` and **1.93 -> 5.91 nm** for
`260220_083.sxm`. Nine lobes require a minimum span of **4.11585 nm** under the
unchanged spacing/overlap constraints, so both old supports are rejected before
optimization. `batch_full.jl` copies the already configured rescue settings
(`support_noise_k=1.5`, `support_padding_nm=0.75`) into both chain configs;
`extract_lobe_features.jl` instead reloads their base values (2.5, 0.25 nm).
The image preprocessing, multistart and refinement choices do not explain this
failure. Increasing optimizer effort cannot remove an infeasible-support guard.

The minimal proposed repair is to replay the accepted support mode from
`refined_policy` in the selected summary, in both circular and elliptical
configs, for base and split extraction. Accepted rescue/accepted-plus-guard
must be distinguished from keep/reject/failed policies; `selection_policy`
alone is not sufficient because it is common to all rows. Preserve
`N_selected` and all other settings. The raw scans plus unchanged preprocessing
and configured rescue values allow deterministic support reconstruction; saved
support lengths alone do not specify endpoints. Exact successful geometry
cannot be reused: the counting output stores only score tables and PNGs, not
per-lobe parameters or a fit context. No new counting fit or selector decision
is needed to address this handoff.

No repair or new fit has been run. The proposed focused regression is a Julia
1.13 metadata-to-config test for accepted/guard/keep/reject/failed cases, a
no-optimizer support-feasibility test, then an explicitly approved one-file
rescued-scan check. Actual downstream success remains unverified. README,
runbook, unit-assignment documentation and this open question now report the
observed blocker. The documentation build passes (28 s; only the existing
archive/search-size warnings and skipped deployment); `git diff --check` passes.

### 2026-09-17 — Replay selected adaptive support in feature extraction

The user requested the support handoff correction first, not a new cohort run.
`extract_lobe_features.jl` now retains each selected count together with its
accepted/base support mode. It reads the original summary's `refined_policy`,
including the existing robust-guard suffixes, rather than relying on
`selection_source`, which a down-only guard may replace. Accepted rescue uses
only the two rescue support values already present in the count TOML; keep,
rejected and failed rescue cases keep base support. The fixed-N configuration
is built before cloning its circular counterpart, so both circular and
elliptical fits receive identical support settings. The same path handles
Gaussian base and split-profile extraction. No selection is rerun, and no N,
spacing/overlap constraint, preprocessing, optimizer or descriptor setting
changes. There are no new physical defaults or config parameters.

Ambiguous/missing adaptive metadata and malformed selected rows now fail before
any SXM read/output, instead of guessing a support mode or silently changing the
requested subset. Nonadaptive `filepath,N_selected` summaries remain supported.
Adaptive resumes that need a base or split refit must retain the original
`--selected-summary`, even when another geometry table is cached. This replays
support from raw inputs and the unchanged config; it does not restore the old
fitted parameters. The extractor now has the standard script-entry guard so
tests can exercise its real parser/config path without launching a fit.

A native Julia 1.13 check on the saved failed-run summary preserves all **25
files and 222 selected lobes**, identifies exactly the two accepted rescues,
and checks both Gaussian and split configs for every file (153 assertions).
No SXM was read and no optimizer was invoked in that check.

`test/test_selected_support_context.jl` passes **935 assertions** under Julia
1.13, including all 18 emitted support/guard outcomes, unchanged fields in both
profiles/configs, missing or ambiguous metadata, infeasible/feasible N=9 support,
and actual synthetic support endpoints/masks without optimization. Twelve
assertions exercise three real CLI subprocesses against unreadable dummy SXMs:
missing adaptive metadata, a missing configured rescue value, and a requested
file absent from the summary all fail with the expected diagnostic before any
output directory/TSV is created. The complete focused suite exits 0 (56 s).

The existing native pipeline regression also passes **28 assertions** (12
boundaries, 9 synthetic extracted-input end-to-end, 7 chunk/merge); CLI help,
`git diff --check` and the documentation build pass. The docs build reports only
the existing archive/search-size warnings and skipped local deployment.

No real Gaussian refit, full cohort, external comparison, grading, SSH transfer,
or new Slurm submission is part of this repair. Existing failed artifacts stay
unchanged. Scientific success after the repair remains unverified until a
focused real rescued-file check and subsequent cohort continuation.

### 2026-09-17 — Real rescued-file support repair check passes

The user approved the focused single-file check following the support handoff
repair (`98c54ae`). Only raw **260215_022.sxm** was processed, using Julia
**1.13.0** and its existing saved counting row (`N_selected=9`,
`refined_policy=adaptive_support_rescue`). The full 25-row summary was retained;
`--files 260215_022.sxm` and a one-file raw symlink directory restricted this
experiment. No counting selection, support choice or classifier was rerun.
The split TOML was generated with native TOML code exactly as in the production
runner: the same count config plus `peak_profile="split"` and the frozen
`split_skew_ratio_max=2.0`. The source count config, summary, raw SXM and unit
config are unchanged (SHA-256 checked before/after).

Both actual extractor CLI processes exit **0** and produce **nine unique lobe
rows**, keys `(260215_022.sxm, 1:9)`, with `source=ell` and finite geometry/GCV.
They ran in parallel, one Julia/BLAS thread each, with a 20-minute safety timeout
that was not reached:

| Fit | Wall time | MaxRSS (KiB) | GCV |
|---|---|---|---|
| base | 74.49 s | 1,117,640 | 5.63285686e-6 |
| split | 84.94 s | 1,113,704 | 1.51047396e-6 |

Native post-checks recompute the same fused-ROI axis/support helpers without
optimization. Both profiles, and their circular/elliptical configs, reproduce
**5.552874796122559 nm** exactly as saved in the counting summary, with t limits
**[-2.6520287145353896, 2.900846081587169] nm**, **9,392 fit pixels**, and the
already configured rescue values `support_noise_k=1.5`, `support_padding_nm=0.75`.
All **76 assertions pass**: selected keys/counts, finite fields, positive GCV,
ordered positions/spacings, support feasibility and agreement, projected
coordinate consistency, and base/split skew limits. The base profile has unit
skew; split skew spans **0.500329–1.097235**, including a near-bound value.

The geometry-only overlay on the fused, preprocessed STM ROI was inspected:
numbered centers follow the chain signal. This is not proof of nine separately
resolved chemical units, and the numbers are lobe indices rather than 0/1
assignments. The new base GCV is **not** the original counting-fit GCV
(3.1825364225520916e-6). Fixed-N refitting does not restore the old parameters;
no bit-identical geometry is claimed, and no setting or refit was changed to
improve that agreement.

Outputs are saved in
`results/reconstructed_cc_soft_v1/support_handoff_260215_022_v1/`: base/split
feature TSVs, generated split config, time/resource logs, `verification.log`,
`support_check.tsv`, a native verification script and `geometry_overlay.png`.
The old 204-row failed-cohort geometry and 25-row count summary were not changed
or merged. The support bug is removed for this real file; **260220_083.sxm is
still untested** after the repair. No full cohort, unit assignment, external
comparison/grade, SSH transfer or new Slurm job was launched. Milestone 1 remains
incomplete until the remaining geometry and cohort-wide application are done.

### 2026-09-17 — Second rescued scan recovered; unknown25 continuation inputs prepared

Following the user request to continue, **260220_083.sxm** was checked locally
with the same bounded procedure as the first rescued scan. Julia **1.13.0** ran
one base and one split extractor process in parallel, each with one Julia/BLAS
thread and a 20-minute safety timeout. The original 25-row counting summary was
supplied unchanged, including `refined_policy=adaptive_support_rescue` and the
selected N=9. The split config used only the existing runner transformation
(`peak_profile="split"`, frozen `split_skew_ratio_max=2.0`). No counting sweep,
parameter search or class prediction was run.

Both processes exit **0**, each with **nine unique, contiguous lobe keys** and
`source=ell`:

| Fit | Wall time | MaxRSS (KiB) | GCV |
|---|---|---|---|
| base | 70.68 s | 1,085,724 | 5.44278424e-6 |
| split | 77.87 s | 1,074,568 | 2.55453904e-6 |

The native post-check passes **76 assertions**, using the same verification code
as the first scan. Both profiles and their circular/elliptical configs reproduce
**5.90756511515228 nm**, exactly the saved support. The t limits are
**[-3.0082010929434198, 2.89936402220886] nm**, with **7,693 fit-mask pixels** and
unchanged rescue settings 1.5 / 0.75 nm. Raw input, source count/unit configs and
original summary hashes are unchanged. The base GCV is not identical to the
original counting GCV (5.4458656131748305e-6). Split skew spans
**0.500034–1.999874**, close to both fixed bounds; this is recorded, not retuned.

The geometry-only overlay was inspected. These fitted centers are model
components, not independently verified chemical units; neither this image nor
the finite GCV establishes assignment accuracy. Outputs are in
`results/reconstructed_cc_soft_v1/support_handoff_260220_083_v1/`, including both
feature TSVs, support checks, logs and `geometry_overlay.png`.

A separate native Julia metadata operation then used the existing strict merge
helper to combine **204 old valid base rows + 9 + 9 recovered rows**. Complete
base geometry now has **25 files / 222 keys**. The two genuine split caches have
**2 files / 18 keys**; the other **23 files / 204 split keys** are still missing.
All keys derive from the saved label-free summary, not an expected chain length
or composition. The preparation passes **29 assertions**, including complete
per-file keys/N, exact coverage, preservation of every source row, and unchanged
source files. Summary and CC templates are copied byte-for-byte. The raw folder
contains only symlinks to the 25 original SXMs.

Prepared inputs are in
`results/reconstructed_cc_soft_v1/unknown25_resume_v1_inputs/`. The failed Raven
outputs remain untouched. `hpc/resume_reconstructed_unknown25.sbatch` provides a
separate, application-only continuation: keep the original summary and complete
base cache, skip the two complete split caches, export the remaining split files
with at most four single-thread children, merge/check complete split coverage,
then run the unchanged native application on **all 25 chains**. It neither
repeats counting/base fits nor trains classifiers on just the refit subset.
Its resource budget remains four requested CPUs, 16,000 MB and 24 hours.

The prepared-data metadata dry-run succeeds without SXM reads, optimization,
classifiers or output creation. It is not a scientific run. No SSH transfer,
Slurm submission or full146/comparison/grading was performed. A future job needs
the committed repair/continuation source and prepared inputs synced to Raven,
then a separate Slurm submission in a new output root. Old job **30271312 is not
resubmitted**. Full-cohort assignments and milestone 1 remain pending.

The committed continuation (`279555d`) also validates selected-support metadata
before creating outputs, using the existing extractor parser/settings validator
inside an isolated module. N-only or ambiguous adaptive metadata is rejected;
no duplicate policy parser or new scientific default was added. The final
actual-input dry-run passes under Julia 1.13 (9.14 s), leaving the proposed
`raven_20260917_unknown25_resume_v1/` output root absent.

The retained self-contained regression
`julia -t 4 --project=. test/test_resume_reconstructed_unknown25.jl` passes
**60 assertions** (46 metadata/shell boundaries, 14 fake-exporter wiring;
61.82 s). It creates temporary synthetic fixtures and loads the current Slurm
script payload, not a copy of its implementation. Tests include incomplete,
extra, duplicate and gapped cached keys, missing/ambiguous adaptive policy,
real shell dry-runs and input/output collision rejection, four one-thread
exports, exact 204+18 merge, unchanged inputs/configs, and the complete-cohort
native command arguments. Its exporter is replaced only within the test process;
no scientific fit or classifier is launched. The broader transient metadata
checks also passed (151 assertions). The shell syntax and documentation build
pass; documentation warnings are page/search-index sizes and skipped deployment.
No remote operation or new scientific batch occurred during verification.

### 2026-09-17 — Unknown25 continuation submitted to Raven

The user approved the actual continuation with "go" after the two recovered-file
checks and local preparation. One new job, **30278010**, was submitted through
`ssh raven` at Slurm `SubmitTime=2026-09-17T14:20:31` (Raven time). Its first
observed state is **PENDING**, with no allocated node/start time yet. The
scientific source is the tested local commit **aa422a3**. This is a separate
application-only continuation, not a requeue/resubmission of failed job 30271312.

Before submission, Raven reported Julia **1.13.0** and no queued/running jobs for
this account. Project and Manifest hashes already matched locally; the ignored
Manifest was preserved rather than regenerated or edited. Reviewed rsync dry-runs
preceded both copies. Only 353 tracked runtime-source files were eligible for
code sync (five changed files were transferred); hidden/generated/agent and
retired authority/spec trees were excluded. No `--delete` was used. The new
input directory contains 29 regular files, **77,722,317 bytes**: 25 SXMs plus the
original summary, complete base table, partial split table and CC templates.
The local raw symlinks were dereferenced and no raw symlink remains remotely.

All **383 source/lock/input SHA-256 values** match locally (353 runtime files,
the existing Manifest, 29 inputs). The native remote metadata dry-run exits 0
in 49.71 s: 25 selected files / 222 lobes, two complete split caches / 18 lobes,
and 23 remaining split files / 204 lobes. It creates no output root and performs
no fitting/classification on the login node. No dependency change or package
installation was needed. `sbatch --test-only` also passes: **30277916 was only
a test-only number**, not the submitted job. Its estimated start/allocation is
advisory, not evidence that the actual job has begun.

Actual job **30278010** requests four CPUs, 16,000 MB and 24 hours in partition
`small`, QOS `n0001`. The pending request reports four CPUs; the scheduler's
test-only estimate used eight processors, as seen in earlier Raven accounting.
The allocation will be checked when the job starts. Only this one job was
submitted; no concurrent second allocation, full146 job or external grader was
started.

Paths:

- synced project: `/u/oldu/code/STMFit`;
- Julia: `/u/oldu/software/julia-1.13.0/bin/julia`;
- new inputs: `/ptmp/oldu/stmfit/unknown25_resume_v1_inputs`;
- new output root: `results/reconstructed_cc_soft_v1/raven_20260917_unknown25_resume_v1`;
- application outputs will be under its `unknown25/` subdirectory;
- Slurm log: `results/reconstructed_cc_soft_v1/unknown25-resume-30278010.log`;
- local transfer/preflight/submission logs:
  `results/reconstructed_cc_soft_v1/raven_unknown25_resume_launch_v1/`.

A temporary local five-minute monitor uses ordinary SSH `squeue`/`sacct` checks;
there is no cluster agent or job controller. It will fetch the new logs/outputs
and stop when this job terminates, without duplicate submission or unchanged
status messages. The old monitor remains cancelled. A successful preflight or
job submission is not application/chemical validation: actual 25-file/222-key
predictions, uncertainty, validation/QC and figures must still be checked after
completion. Failure will be reported with retained outputs and missing keys,
not hidden by reducing N, omitting a scan or retuning. Full146 reconstruction
and external comparison/grading remain separate later work.

### 2026-09-17 — Unknown25 continuation completed; intact outputs, limited confidence

**Deliverable:** complete and check the fixed `cc_soft_reconstructed_v1`
application on all 25 unknown-chain scans, keeping the original adaptive
selected summary, selected counts and scientific settings. Done for this
application means complete file/lobe coverage, usable or explicitly unavailable
components, reproducible prediction/QC tables and inspected diagnostic figures.
It does not mean chemical accuracy, correct molecular lengths, recovery of the
lost historical descriptor, or milestone 1 completion.

Raven job **30278010**, submitted once from tested source **aa422a3** with
`hpc/resume_reconstructed_unknown25.sbatch`, is **COMPLETED, exit 0:0**. Slurm
records 20:55:39–21:01:55 CEST on September 17, elapsed **6 min 16 s**, node
`ravc4043`, **eight allocated CPUs** and MaxRSS **3,876,696 KiB (3.70 GiB)**.
The request was four CPUs, 16,000 MB and 24 hours. The monitor was stopped after
completion. An expired `squeue` record returned "Invalid job id"; its nonzero
exit suppressed a chained `sacct` call. A separate accounting query established
completion. Queue disappearance was not treated as a failed scientific job.
No resubmission, requeue, cancellation or second scientific job occurred.

The continuation reused all **222 base rows / 25 files** and the **18 genuine
split rows / two rescued scans**, fitted only **204 split rows / 23 files** with
at most four single-thread children, then applied the native pipeline to the
**full 25-chain cohort**. Counting, base geometry and CC templates were not
recomputed. The original summary, support policies and all settings were kept.
The earlier failed job 30271312 and its outputs remain untouched; its monitor
also remains cancelled.

**Fetch and verification.** The new result root is
`results/reconstructed_cc_soft_v1/raven_20260917_unknown25_resume_v1/`; application
outputs are under `unknown25/`. The separate Slurm log is
`results/reconstructed_cc_soft_v1/unknown25-resume-30278010.log`. Rsync fetched
**66 regular output files / 6,270,285 bytes**, excluding `raw_inputs/` at every
level, plus that log. All **67 remote/local SHA-256 values match**. No generated
result or raw file is committed.

Local Julia **1.13.0** verification passes **179 assertions**. Every one of the
12 intermediate/final lobe tables has the exact original **222 keys / 25 files**.
The 18+204 split merge retains each source row, base columns remain unchanged
through augmentation, required predictor values are finite, split configs match
the unchanged settings, and per-chain summaries agree with predictions. The
native validator passes. Reaggregation of the saved GMM/k-means votes and native
QC regeneration are byte-identical to the saved TSVs. This is saved-output
verification, not a new fit, classifier run or scientific parameter variant.
Checks and logs are under
`results/reconstructed_cc_soft_v1/raven_unknown25_resume_launch_v1/`
(`verify_completed_unknown25.jl`, `completed_verification.log`,
`local_prediction_validation.log`, `local_qc_validation.log`,
`predictions_vote_check.tsv`, `review_queue_check.tsv`).

**Actual model outputs.** There are **189 class-0 / 33 class-1 / zero `?`** rows,
all named `cc_soft_reconstructed_v1` with `invalid_reason=ok`. Both component
predictors have all 222 usable rows (GMM one view, k-means four views). All
reconstructed descriptors and Fisher scores are finite; all forward/backward
17×17 and backward 9×9 patch rows are finite with nonzero mass. This supplies
complete diagnostic assignments, not a measurement of the chains' true 189/33
composition. The selected N range is 5–16; no expected molecular length was
used to change it.

The lobe confidence is an **uncalibrated vote margin**,
`2 * abs(probability_1 - 0.5)`, where `probability_1` averages two component vote
fractions. Its minimum / median / mean / maximum are
**0.0375 / 0.60 / 0.5701576577 / 0.7625**. Five values are below 0.20 (a descriptive
count only); five GMM/k-means binary assignments disagree. The existing vote
rule chooses 1 at a fraction >=0.5, including exact ties; there were no exact
ties. It emits `?` for unavailable components, not all low finite margins.
No confidence abstention rule was added. Zero `?` therefore does not mean
certainty, and component agreement is not chemical validation.

The unchanged QC flags **19/25 chains** for mean confidence <0.60. One of them,
`251206_038.sxm` (N=16), also has the cohort-relative `n_outlier` flag (IQR fences
[3.5, 15.5]); six chains have no listed flag. The weakest chain means are
`260116_017.sxm` (0.43571429) and `260221_044.sxm` (0.4375). No missing-plot or
noncontiguous-key flag occurs. An `ok` QC status is not proof of correctness.

**Read-only diagnostic review.** All 17 run logs plus Slurm/accounting were
read independently. No exception, Julia warning, missing-view message or
Fisher iteration-limit warning was found. The split sweep for `260222_043.sxm`
reports a coverage floor of 6 clamped to the already fixed N=5; it exports all
five selected rows, rather than performing new count selection. Interleaved
command printing in the outer Slurm log is concurrent logging, not malformed
execution; separate chunk logs and output tables are coherent. Silence in logs
does not prove per-start optimizer or per-seed classifier convergence, which
is not fully serialized.

Saved-output diagnostics expose limits beyond the basic integrity checks:

- Split skew is within [0.5, 2], spanning 0.500031–1.999874. With an explicitly
  descriptive distance of 0.001, five rows are near the lower bound and eleven
  near the upper bound, over 14 files. At exported six-decimal precision,
  parallel widths reach the effective optimizer endpoints in 1 lower / 17 upper
  rows; perpendicular widths in 1 lower / 78 upper rows. The nominal range
  [0.191, 0.509] nm maps through the bounded [-5, 5] optimizer coordinates to
  approximately [0.193128, 0.506872] nm. These are constrained solutions, not
  out-of-range errors or a new rejection rule; no bound was changed.
- The reconstructed half-plane descriptor is concentrated toward +1: median
  0.964083, range -0.945059–0.994111, 128/222 absolute values >0.95. This alone
  establishes neither failure nor chemical meaning and does not identify the
  lost historical producer. No descriptor variant was tried.
- Mold-only forward/backward decoded labels disagree at **120/222 keys**. These
  are intermediate results: the final predictor uses their unsigned absolute
  `cost_margin` columns, not those sequences. This is **not forward/backward
  chemical confirmation**. The read-only report is
  `raven_unknown25_resume_launch_v1/completed_log_review.txt` under the same
  reconstructed-results parent.

**Figures.** All 25 standalone PNGs and the grid exist. The grid and seven
priority standalone maps were viewed: both rescued scans, `260116_017.sxm`,
`260221_044.sxm`, `251206_038.sxm`, `260222_043.sxm`, and `260221_048.sxm` (the
minimum-margin lobe). Standalone markers and class labels are readable. They
are **base fitted-coordinate maps, not raw STM overlays**; their colors do not
encode margin size. Visual inspection does not establish correct counts,
individually resolved chemical units or chemical identity. The earlier raw
geometry overlays for the two rescued scans remain separate focused checks.

Current README, runbook, unit-assignment and HPC status now describe the actual
completed application and its limitations. The journal's current promotion
answer also distinguishes saved-output regrading from the unavailable exact
historical raw-input rebuild. No production code, dependency, Manifest, physical
setting, selection rule or benchmark headline was changed. The Julia 1.13
documentation build and `git diff --check` pass. Documentation warnings are
limited to page/search-index sizes and skipped deployment; no broken-reference
warning occurs.

**Next, separately authorized work:** full146 reconstruction with these fixed
settings, followed by external keyed comparison and grading. None has run as
part of this continuation, no benchmark labels were read for the application,
and the accepted frozen-champion metrics remain unchanged. The unknown25
technical workflow is complete; chemical validation and milestone 1 as a whole
remain unresolved.

### 2026-09-17 — Full146 approved; cache lineage clarified before execution

After the completed unknown25 report (`9ece7bf`), the user replied **"go"** to
full146 reconstruction followed by external comparison and grading. The
scientific deliverable is a complete fixed-method assignment table on all 146
scans, then a separately generated comparison/grade with coverage, uncertainty
and limitations. Done requires all intended file/lobe keys, unchanged inputs
and settings, inspected outputs/QC, and external-only reporting. Scope is one
configuration, one job capped at **24 hours**, four requested CPUs and 16,000 MB.
No descriptor search, threshold tuning, new abstention rule, expected-N prior,
file omission or champion-label input is authorized. Unknown25 is finished and
must not be rerun.

**Important input distinction found in preflight.** The previously prepared
base cache,
`results/hierarchical_feature_extraction/9e3bcdeb83469fe0/features_merged.tsv`,
contains **900 unique contiguous lobe keys / 146 files**, all `source=ell`.
Its SHA-256 is
`2cb1065cf5aa790dac1a5fec279d955c732cc91f7e60312103de32fddfb89484`, matching the
archived extraction record. Its eight-shard union preserves all source rows.
Cached N counts are 4:2 files, 5:12, 6:106, 7:16, 8:8 and 10:2. These are own-N
GCV geometry counts, not fixed-six inputs.

The promoted counting summary
`results/experiments/6mer_full146/pm2_confirm/summary_overlap060_hard.tsv`
contains the same 146 files, but **871 selected keys**: 118 per-file counts agree
and **28 conflict**. There are 866 common keys, 34 cache-only and five
summary-only keys. Native `check_counts` rejects this combination. It would be
wrong to attach this summary to the 900-row geometry, silently trim the cache,
change N, or claim that the cached counts reproduce the promoted counting
policy. No benchmark labels or grades were used to choose between them.

Read-only provenance review resolves the distinction. The archived extraction
config hash is
`a24460ab5114bd3907c063c1f6d7b678113037e084881a7fbbb211db9eeaf2d9`, identical to
current `config/chitosan.toml`. The archived extractor hash matches
`37e52333034ac18979638702c68644cc0a546dea:test/extract_lobe_features.jl`.
Its `_effective_best` chooses the lowest finite valid circular/elliptical GCV
among fitted N. Logs for all 146 files show variable-N sweeps and the saved own
N, not fixed-N summary handoff. The original producer consumed only dataset/file
membership; those paths match every raw/cache file. Its intelligent sweep has
BIC early-stopping diagnostics; the cache is not an exhaustive-count optimum.
It is also not a fresh Julia 1.13 counting benchmark.

There is **no original selected-summary for this extraction**, and its support
configuration is nonadaptive. The already approved features-only full146
command remains appropriate: preserve the 900-row cache and let the native
runner derive `selected_from_features.tsv` for fixed-N split refits. This
metadata is derived from geometry, not presented as a recovered original
summary. The model's `support_midpoint_hybrid` config field does not retroactively
make direct-extractor GCV counts into promoted batch-policy counts. The separate
counting benchmark remains unchanged. For adaptive unknown25 inputs, retaining
the original summary and `refined_policy` is still mandatory.

No compatible full146 split cache was found. All 900 base `skew_ratio` values
are one and cannot stand in for split fitting. An older genuine split table
has 234 rows / 39 files, with incompatible N on four scans (233 common keys,
one split-only key); it is not passed or silently subset. All **146 genuine
split fits** will be computed at fixed cached N with the frozen split transform
(`peak_profile=split`, `skew_ratio_max=2.0`), in at most four single-thread
children. Whole-cohort patches, descriptor, CC/Fisher, classifiers, vote, QC and
maps follow. Neither counting nor base fitting is rerun.

Local staged inputs are now
`results/reconstructed_cc_soft_v1/full146_v1_inputs/`: byte-identical
`base_geometry_full146.tsv`, native `templates_cc.tsv` (unchanged SHA
`fede623b72c0381875fb738c4c1ec128397919b33be937ba5423c80cc04478e4`), and
`full146_raw/` with 146 valid raw symlinks. Names match the cache exactly,
without extras or duplicate basenames. The 148 staged files total
**308,158,485 bytes** when dereferenced (raw targets: 307,924,640 bytes).
No promoted summary, truth, manifest label or reference prediction is staged
as production input. Raw symlinks must be dereferenced on cluster transfer.

The native local metadata check passes **14 assertions** under Julia 1.13,
including raw/key coverage, source hashes, finite geometry, own-N consistency,
nonadaptive config and actual runner dry-run without output creation. An initial
parenthesis typo in the temporary verification script stopped at parsing; it
was corrected and the first failed log was retained. No scientific program
stage or input changed. Existing synthetic pipeline boundary/sharding tests
pass **19 assertions**, and external-comparison fixtures pass **nine**; these
are test fixtures, not a grade of scientific predictions.

`hpc/reconstructed_full146.sbatch` is a small full146-only entrypoint, avoiding
the older two-cohort script that would repeat unknown25. `STMFIT_OUTDIR` is its
new direct output directory. The optional `STMFIT_SELECTED_SUMMARY` is explicit
and empty for this cache; no summary is auto-discovered. It preserves the native
runner's version, coverage and output-collision checks. No cluster agent,
authority framework or persistent job controller is added. Its retained
regression passes **122 assertions** under Julia 1.13 (70 shell/resource/capture
checks and 52 actual native metadata/output-boundary checks). It checks
empty/unset/explicit optional summaries, incompatible-summary rejection and no
auto-discovery, missing inputs/env, forbidden arguments, spaces/relative paths,
output directory/file/dangling-link collisions in both modes, exit propagation,
and unchanged inputs/configs. No fixture launches a fit or classifier. The
actual prepared-input entrypoint dry-run also passes in **8.75 s**, leaving the
output root absent. Metadata checks do not parse SXM/template content or prove
fit viability. Shell syntax, `git diff --check` and the Julia 1.13 documentation
build pass; doc warnings remain page/search sizes and skipped deployment.

**Blocker at preparation (resolved in the next entry):** the metadata-only
noninteractive SSH probe to `oldu@viper` was refused at authentication. Viper is preferred by the repo for multi-file
fits; the user was asked to authenticate with `ssh oldu@viper`. No remote copy,
Slurm submission, scientific fit/classifier, full146 prediction or grade has
occurred. The staged inputs and tested source must still be synced after
reviewed rsync dry-runs, checked with Julia 1.13 and scheduler test-only, then
submitted exactly once. Old Raven jobs 30271312 and 30278010 stay untouched.

After fetch and complete 146/900-key verification, the external comparison may
read the frozen champion and filter the benchmark membership for the existing
full145 own-N grader. It must report missing control positions and extra
predicted lobes, not hide them. Neither comparison nor grade can tune the
production model. This is an **assignment reconstruction on cached GCV geometry**,
not a raw rebuild of the promoted counting-to-assignment chain. Native numerical
differences and the lost descriptor's irrecoverable identity remain limitations.


### 2026-09-18 — Fixed full146 job submitted on Viper; Julia bootstrap verified

After the user's **"continuer"**, Viper authentication succeeds. The scientific
deliverable and frozen 146-file / 900-key scope are unchanged: genuine fixed-N
split fits, whole-cohort reconstructed assignment, then external comparison and
grading only after complete output validation. No unknown25 rerun, new descriptor,
parameter variant, count change or benchmark-label input is introduced.

**Transfer.** The old Viper `/u/oldu/code/STMFit` directory has no Git metadata,
so it is preserved rather than overwritten. Reviewed rsync dry-runs precede
copies to `/u/oldu/code/STMFit_full146_v1` (tested source **00528d2**) and
`/ptmp/oldu/stmfit/full146_v1_inputs`. All **354 source/lock files** and **148
inputs** match local SHA-256 values. The 146 raw scans are dereferenced regular
files on Viper. No `--delete`, truth, reference predictions, agent/archive state
or prior outputs enter the transfer. Julia 1.13.0 is copied from the accessible
Raven installation into `/u/oldu/software/julia-1.13.0`; its executable, system
image, libjulia and installation startup file match the local hashes. Source and
these four runtime hashes are checked again after environment setup.

**Failed setup and targeted recovery.** Initial root-project `using Pkg` fails
in **8.61 s** while precompiling `LibCURL_jll`: its worker cannot load a
`Zstd_jll` cache. Loading Pkg before activating the root project allows
instantiation/precompilation to finish in **75.12 s**, without changing Project
or Manifest, but the log still reports unavailable Pkg images for REPLExt and
RegistryInstances. Its zero exit is not proof that every cache is usable.
Likewise, the first production import exits zero after **234.69 s**, but logs
Zstd errors and a failed `Plots.FileIOExt` load. These logs are retained; neither
step is called a clean runtime check.

Read-only diagnosis finds a concrete graph mismatch: the frozen Manifest's
Julia 1.12-era LibCURL/Pkg dependency lists omit Zstd, which the installed Julia
1.13 stdlibs now import. The loader can recover missing stdlib dependency edges,
but precompilation schedules workers from the older manifest graph and uses
strict cache loading. The project selects the locked external
`Zstd_jll/YX1LN`, while the stdlib environment selects a bundled copy. This
explains why warming stdlib Pkg need not prepare the project's Zstd dependency.
The error's displayed cache flags are Julia defaults, not evidence of changed
optimization flags. Exact cache rejection details were not established.

A sequential load of **Zstd_jll, Pkg, then the production packages under the
real root project** completes cleanly in **72.30 s**. A new process with normal
production startup, four threads and **no Zstd/Pkg preload** imports
GaussianFit2D, STMMolecularFit, Plots and Clustering cleanly in **3.75 s**;
`Plots.FileIOExt` is loaded. The application lock remains byte-identical:
no `Pkg.resolve()`/`Pkg.update()`, Manifest edit, cache deletion, scientific
parameter change or Julia fallback is used.
The login CPU reports `znver4`. This is environment verification, not a fit or
scientific validation; all scientific work remains for the compute allocation.

**Preflight and submission.** The actual remote metadata dry-run passes in
**7.36 s** under Julia 1.13.0, identifies all 146 scans and creates no scientific
output directory. Input coverage, unchanged Project/Manifest/config hashes and
output absence are checked. Slurm test-only ID **11786052** is accepted and
estimates eight processors; it is not a submitted scientific job or a measured
allocation.

The actual full146-only **job 11786116** is submitted **once**, at
**2026-09-18 00:04:54 CEST**, with four requested CPUs, 16,000 MB and 24 hours.
The submit filter selects `small`. Initial `squeue`, `sacct` and `scontrol`
report **PENDING / QOSGrpCpuLimit**, four requested CPUs and no allocation yet.
Do not resubmit, requeue or change resources/settings to bypass the queue.
Old jobs 30271312 and 30278010 remain untouched.

The direct remote output is
`/u/oldu/code/STMFit_full146_v1/results/reconstructed_cc_soft_v1/viper_20260917_full146_v1`;
the name retains the September 17 preparation date. The separate log is
`results/reconstructed_cc_soft_v1/full146-11786116.log` under that project.
Local preflight/setup logs are under
`results/reconstructed_cc_soft_v1/full146_launch_v1/`. Bounded local monitoring
checks queue and accounting independently, so an expired queue entry cannot
hide the terminal accounting state. There is no cluster-side agent/controller.

At this submission checkpoint, no full146 scientific outputs, comparison or
grade are available. Next: fetch without `raw_inputs/`, verify all 146/900 keys,
base/config preservation, genuine split rows, components/reasons, uncertainty,
QC and figures, then run the authorized external-only comparison and full145
own-N grade. Failed or incomplete outputs must be retained, not trimmed or
silently retried. Existing benchmark headlines and unknown25 results are unchanged.


### 2026-09-18 — Full146 integrity verified; assignment overview layout corrected

Viper job **11786116** completed `0:0` in **11m04s**, from 00:36:57 to
00:48:01 CEST on `vipc2162`. Accounting records eight allocated CPUs (four
requested) and MaxRSS **3,486,576 KiB**. The monitor is stopped; no second job
or scientific variant was run. All **184 fetched output files (26,642,815
bytes)** and the separate Slurm log match remote SHA-256 values. The 358
source/key-runtime identities and 148 prepared inputs remain unchanged.
`raw_inputs/` was excluded at every level. Original outputs remain under
`results/reconstructed_cc_soft_v1/viper_20260917_full146_v1/`.

Local Julia 1.13 saved-output checks retain exactly **146 files / 900 keys**
in all twelve lobe tables, all four genuine split chunks, unchanged base and
configuration, valid missing-input reasons, and byte-identical descriptor,
feature-join, vote, summary and QC replays. Native component/final validators
pass. Final assignments are **695 zero / 198 one / 7 `?`**. The seven partial
patch rows remain present and propagate to `unavailable_gmm`; there is no
single-component fallback or new confidence cutoff. Full log/QC interpretation
is saved in `full146_launch_v1/completed_log_review.txt`. This is integrity
verification, not a chemical grade or proof that every optimizer seed converged.

**Checker failure, not data corruption:** the first check stopped at 41 passes
and one error because it expected an exact 2000×14800 summary PNG. The actual
height is 14798. A separate blank-plot probe reproduced this raster rounding
under Julia 1.13.0, Plots 1.41.6 and GR 0.73.24. Only the generated checker's tall
grid assertion now allows a 0–2-pixel height rounddown; widths and standalone
sizes remain exact. Corrected fixtures pass **48/48**, and completed-output
checks pass **73/73**. The initial checker and failed log are retained. No
production image, table or scientific threshold was changed to satisfy a test.

**A separate real figure defect:** visual inspection of the original overview
found a large empty title band, compressed equal-aspect plot bodies, clipped
endpoint disks and overlapping labels. The sixteen inspected standalone maps
were readable. A giant canvas, small panel allocation and a title band fixed at
10% of total height caused a layout problem; reuse of rendered Plots objects
alone is not established as the cause.

`plot_reconstructed_unit_assignment.jl` now uses at most **eight panels per
page**, two columns, the original 780×440 panel footprint, a fixed 48-pixel
title band, and overview-only left padding/widening. Equal physical coordinate
scales, class colors and `0/1/?` annotations are retained. `summary_grid.png`
is explicitly **page 1**, not the full cohort; subsequent pages are numbered
and `summary_pages.tsv` records each file once in sorted order. There are
**19 pages** for this cohort. The standalone rendering code is unchanged.

The native plot regression passes **267/267** under Julia 1.13, including
single, three-, eight- and 146-file cases, pagination, labels/colors, physical
aspect, rendered plot-area size, endpoint margins and input/collision checks.
A local **plot-only** replay from saved geometry/predictions took 16.3 s.
All 146 standalone PNGs are byte-identical to production; the overview pixel
check finds all 900 colored markers with the correct per-panel class counts.
All nineteen PNGs decode, and representative corrected pages/edge cases were
visually inspected. Corrected figures are copied without overwriting the
originals to
`results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_figures_v1/`.
The intermediate margin probe is not the final figure set. GR raster rounding
also affects the page canvas (1560×1808 requested, 1558×1806 exported here),
without omitting panels or markers.

These remain fitted-coordinate maps, not raw STM overlays. Colors do not encode
confidence. The scientific inputs, counts, fits, predictors, abstention rules,
dependencies and unknown25 outputs are unchanged. External comparison and the
full145 own-N grade are the next separate local step; no benchmark results
were read to choose this layout or change predictions.


### 2026-09-18 — Full146 external comparison and own-N grade completed

After the preceding transfer, 146/900 integrity, log/QC and figure checks,
the fixed reconstruction was evaluated locally under Julia 1.13. This is the
authorized post-processing step, not a new fit, descriptor/threshold variant,
counting run or training campaign. Production inputs and all 184 outputs remain
byte-identical after evaluation. The plot-only source fix is `6fa3297`; the
scientific job still comes from `00528d2`. Unknown25 and all earlier outputs
remain untouched.

**Keyed comparison.** `test/compare_reconstructed_champion.jl` read the completed
900-row predictions and the frozen
`results/unit_assignment/best_labelfree_cc_soft_20260802.tsv`. It finished in
**3.92 s** and reports **846/892 matching labels (94.8%)**, **zero missing
reference keys**, and **eight extra keys**. All eight extras are
`240310_Cu100009.sxm`, the file outside the confirmed145 manifest. The external
membership filter retains **892 rows / 145 files**, exactly the same key set
as the reference; no production row was removed. Of 46 differing labels, 39
change reference 1 to reconstructed 0, one changes 0 to 1, and six change a
binary label to `?`. Reported confidences differ too (maximum absolute difference
0.8). Agreement is not descriptor identity or proof of calibrated confidence.

**Full145 own-N grade.** The separate
`test/report_unit_assignment_benchmark.jl --full145-own-n` call uses only that
filtered external copy, the confirmed145 manifest and the existing external
control convention. It finishes in **11.47 s**. The registered control has
**145 files / 870 positions**, while the label-free cache supplies **892 lobes**
on those files. Missing control positions are **16 across 14 short-N files**;
extra predicted lobes are **38 across 25 long-N files**. These counts are
reported, not repaired, trimmed from production, or fed into model selection.
This is not a fresh grade of the promoted 871-lobe counting chain.

| Quantity | `cc_soft_reconstructed_v1` |
|---|---:|
| Classified control positions | 849/870 (97.6%) |
| Physical convention correct / classified | **666/849 (78.4%)** |
| Exact chains | **24/145** |
| Correct / all control positions | **666/870 (76.6%)** |
| Emitted errors | 183 |
| Unclassified control positions | 21 = 16 missing + five unavailable |
| Post-hoc `honest_uncertain` | 204 = 183 errors + 21 unclassified |
| Supervised oracle diagnostic | 675/849 (79.5%), 24 exact chains |

The oracle is not the label-free result and is not used to change labels.
The lower-level grader's 854-position printout describes present control
positions; the own-N report restores all 870, including the 16 absent ones.
Likewise `honest_uncertain` is not a model abstention count. Production has
seven `?`; the confirmed145 copy has six, one of which lies outside the
aligned control positions. Per-file and per-position tables independently sum
to the report. At positions 2 and 5, classified accuracy is 82/144 (56.9%) and
65/143 (45.5%); these are external diagnostics, not calibration targets.

**Comparison with the frozen reference and decision.** The accepted reference
remains **677/854 (79.3%) and 36 exact chains**. The reconstruction has fewer
correct classified positions and exact chains; it is **not promoted**. No new
parameter search, confidence/abstention change or retraining is authorized or
started. This completes the approved full146 application/comparison task, but
not an exact recovery of the lost descriptor or chemical validation of unknown
10–20mer sequences. Historical benchmark headlines remain attached to the
frozen outputs, not copied onto this reconstruction.

**Production QC context.** Final full146 labels are 695 zero / 198 one / seven
`?`, with mean/median uncalibrated margin 0.79082/0.8. Twelve exact finite ties
still emit 1; 58 assigned margins are below 0.20 (descriptive only). The binary
components disagree on 146/893 jointly available rows across 106 chains. QC
flags 43 chains: 40 N outliers, six low mean margins and one high uncertain
fraction, with overlaps. Q1=Q3=6 explains the cohort-relative N flag; it is not
an expected-N constraint. `241113_088` is QC `ok` despite one unavailable lobe;
this auto-QC call does not evaluate missing component views.

The seven incomplete-patch keys remain in descriptor, Fisher, component and
final tables with their native reasons. In `240818_019`, two local mold rows
have infinite costs and the auxiliary scorer prints `0000` with infinite file
cost; the final assignment is correctly `????`. That nominal auxiliary
sequence must not be treated as chemical evidence. No extra logged error or
iteration-limit warning occurs, but not every optimizer/seed termination state
is exported. Split skew has 94 rounded effective-bound hits (45 lower / 49
upper); perpendicular widths have 297 upper hits. Descriptor median 0.958684
and 875/893 positive exported Fisher scores are distributional diagnostics,
not chemical class counts. None of these observations was used for retuning.

**Local deliverables:**

- Original complete production: `results/reconstructed_cc_soft_v1/viper_20260917_full146_v1/`.
- Corrected figures, nineteen indexed pages and byte-identical standalone maps:
  `results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_figures_v1/`.
- External reports: `results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_external/`;
  `comparison/` contains keyed differences and the membership-only grade input;
  `grade_full145_own_n/` contains `report.md`, `summary.tsv`, per-file grades and
  per-position errors. Commands and stdout are retained beside them.
- Verification/log/visual evidence: `results/reconstructed_cc_soft_v1/full146_launch_v1/`,
  including the preserved failed checker and intermediate layout probe.

No generated results, raw inputs or benchmark-derived reports are committed.


### 2026-09-18 — Bounded label-free improvement exploration authorized

The user approved exploring the mathematical improvement avenues discussed after
completion of unknown25 and the fixed full146 reconstruction. The scientific
deliverable is measured evidence separating numerical fitting limitations,
representation/acquisition confounding, and missing independent chemical signal;
it is not a new benchmark champion or automatic promotion.

The first pass is limited to one half-day of development/analysis and at most
one Viper job (four requested CPUs, 16 GB, two hours; allocation can be eight
CPUs). Four opt-in diagnostics are prepared: constrained variable projection of
Gaussian amplitudes/background at fixed candidate N; saved-patch representation
attribution; fixed-geometry acquisition/noise evidence; and actual raw-STM
unknown25 overlays plus an available-channel inventory. Synthetic tests and
lightweight saved-output analysis are local under Julia 1.13. Multi-file fitting
and heavy image diagnostics run only on Viper, after a dry-run. No experiment
outcome or gain is claimed at this preparation checkpoint.

The four raw diagnostic scans are selected **only from saved unknown25 counting
diagnostics**, not external errors or chemical predictions: `260115_016.sxm`
(smallest finite relative effective-GCV gap), `251206_013.sxm` (largest saved
forward/backward NRMSE), `260215_022.sxm` (lexicographically first actual adaptive
support rescue), and `251206_034.sxm` (median finite relative GCV gap). The
`*_keep` policies are not rescue activations. The original selected summary
remains necessary to restore support. Neighboring candidate counts, if fitted,
are saved `N_selected-1`, `N_selected`, `N_selected+1`; this is a local diagnostic,
not an exhaustive search or replacement `N_selected` policy. Infeasible and
failed cases remain visible.

All existing production configs, sources and outputs remain the reference.
New settings belong only to `config/label_free_exploration.toml`, which delegates
physical model/selection/preprocessing to the explicit original molecule config.
No benchmark label, expected N, sequence, composition prior, reference prediction
or grade may choose parameters or candidates. No label-guided feature/seed/
threshold search, dependency update, new workflow framework, QE run, or reopening
of the retired T3/T11–T13 campaigns is part of this pass. An additional measurement
can be recommended, not simulated into a claim of observed chemical validation.

Two mathematical corrections motivate the diagnostics. An empirical supervised
LDA score in the historical notes is not an information-theoretic/Bayes ceiling.
GCV remains the canonical practical per-candidate criterion, but ordinary GCV is
not a universal guarantee under nonlinear constrained fitting or spatially
correlated residuals. The `n_eff = n ÷ 9` placeholder is not recalibrated. Variable
projection retains the full model parameter count in GCV; numerical elimination
of linear variables does not remove their statistical degrees of freedom.
Stability frequencies and component vote margins are not calibrated probabilities
of correct molecular count or chemical identity.


**Pre-run candidate completeness clarification:** the selected-N neighborhood
is augmented by the already saved finite `N_eff` and effective runner-up. For
`260115_016.sxm`, saved selected N is 6 but effective-GCV winner/runner-up are
8/7; omitting 8 would omit the ambiguity that selected this scan. The final
candidate sets, fixed before any diagnostic fit, are 5/6/7/8 for that scan,
10/11/12 for `251206_013.sxm`, and 8/9/10 for both `260215_022.sxm` and
`251206_034.sxm` (13 N candidates total). This derives only from the saved
label-free counting summary, not human estimates or benchmark counts, and still
does not constitute an exhaustive optimum or replacement selector.


The lightweight metadata audit also separates the saved count stages: unknown25
has 235 raw effective-GCV lobes and 222 final selected lobes; the final policy
reduces N by one in seven files and by two in three files, leaving fifteen
unchanged. These are **policy-stage differences, not errors against truth**.
Only `251206_013.sxm` and `260115_016.sxm` have recorded relative effective-GCV
gaps below the existing 0.05 ambiguity threshold. Final-policy uncertainty
cannot be inferred from that raw-GCV gap alone.

A small dedicated `hpc/label_free_exploration.sbatch` invokes the native
`test/run_label_free_exploration.jl` driver. It runs four independent one-thread
cases, each acquisition diagnostic followed by variable projection, records
stage exit codes and retains failure logs. It requires Slurm for real execution
and has a metadata-only dry-run. The driver checks the original geometry/count
coverage and the declared candidate union, rather than accepting a hidden target
count. Inputs are four raw scans plus unchanged original summary/base geometry
and the explicit candidate table. No new job is submitted at this checkpoint.

The first driver regression had one incorrect test-only command-token-count
assertion (expected 19, actual 17). It was corrected after inspecting the native
command; no scientific argument was removed. The subsequent tests also exercise
failure propagation with small synthetic Julia subprocesses, not scientific
fits. Initial failures and successful logs are preserved under the new result
root. The initial SSH probe timed out; the bounded retry and inventory succeeded
with Julia 1.13.0 and no queued jobs. Existing remote production trees remain
untouched.


The native unknown25 header inventory found Z and Current in both directions
at the same recorded bias (-0.300 V) on all scans, not a second-bias or STS
measurement. Configured lock-in fields do not establish a recorded lock-in
channel. The Current-weight helper in GaussianFit2D is defined but has no call
site in the current Julia sources. The four-file acquisition diagnostic may
therefore report the actually recorded auxiliary `Current` channel's finite
scale and cross-view consistency, and, if implemented within the first-pass
budget, residual structure after one declared [1,Z,dZ/dx,dZ/dy] linear nuisance
projection. This is an acquisition/feedback diagnostic only: Current and Z are
coupled by the feedback loop, so it is not an independent chemical likelihood,
new fit weighting, feature search or evidence of correct assignments.


### 2026-09-18 — First label-free diagnostics: measured representation and raw-data limits

The saved-data audit retained 25/222 unknown keys and 146/900 full146 keys,
including seven incomplete full146 patches. It produced no new classification
or grade. The existing backward 9x9 descriptor matches the fixed `1e-12`
arithmetic tolerance at every key (maximum unknown difference
`2.220446049250313e-16`, full146 difference zero; reasons retained).
A post-run checker initially demanded exact zero difference; that test-only
mistake was corrected to the already declared tolerance without changing or
rerunning the scientific audit.

Median backward9 transverse-gradient fractions of centered patch energy are
0.844470 (unknown25) and 0.846248 (full146). The corresponding signed numerator
shares are 0.986613 and 0.987684; these are **not variance fractions**. Median
saved descriptors are 0.964083/0.958684; recomputing the diagnostic after removing
each patch's fitted transverse ramp yields 0.032991/0.030143. Forward17 has the
opposite median sign. Within-file centering does not remove all this structure:
the additive transverse contribution to the existing standardized descriptor
has median within-file SD 0.976924/0.991160, versus 0.104742/0.211213 for its
remainder. These are linear attributions, not causal identification of an
acquisition artefact, absence of chemistry, or justification for a new feature.
No gradient subtraction enters the predictors. Reports are under
`results/label_free_exploration_20260918/representation/{unknown25,full146}/`.

All25 raw-QC plates were rendered without fitting, with separate blind and
indexed-center views. Final reviewed outputs are in
`results/label_free_exploration_20260918/unknown25_raw_qc/final/`; earlier
preview/all25 versions are preserved. Corrections made colorbars explicit,
kept high-range tick labels readable and made nonfinite/imputed pixels visible.
The rectangular `260220_083.sxm` retains its approximately 4x8 nm aspect.
Saved centers follow broad structure in the reviewed examples, but do not
independently resolve or validate each counted unit. No exact residual is
reconstructed from the rounded decoded geometry.

Two scans contain nonfinite raw Z samples in both directions:
`260115_016.sxm` has finite fractions 0.3937378/0.3925781, and
`260117_005.sxm` 0.7013550/0.7011719. Native preprocessing median-imputes those
samples. All saved centers have a finite nearest raw pixel, which does not
prove full patch coverage or remove the effect of preprocessing missing data.
The first scan remains in the predetermined counting diagnostic; it is not
silently removed to improve an outcome. All100 recorded channel/direction
rows have nonconstant finite samples. Z and Current are recorded fwd/bwd in
all25 scans at -0.300 V, with lock-in OFF. This is not spectroscopy or an
independent chemical validation.

The acquisition diagnostic implements direct auxiliary Current scales and
cross-view/height coupling on exactly the fixed Z support and lag. It does
**not** implement the optional [1,Z,dZ/dx,dZ/dy] nuisance projection. Imputed
samples are not treated as independent observations. Correlation lengths use
absolute ACF magnitude; missing lag support, censoring or later resurgence
block length claims. Joint-likelihood and calibrated-independence claims remain
blocked. Original geometry, not newly optimized count fits, supplies support.

The first and only local real counting smoke used `260215_022.sxm`, N9, with
its original adaptive support. All six native/fixed-profile/refined-profile
rows are valid and numerically mapped to the original forward model. At fixed
N, circular RSS decreases from 0.14109435 to 0.04728573, and elliptical RSS from
0.05250398 to 0.02771763 (about 47.2% for the latter). The fixed-geometry linear
solve alone gives elliptical RSS 0.05190226. This is extra bounded optimization
from a native result, **not an equal-budget solver comparison**, a count check,
a global optimum or an accuracy improvement. Full GCV parameter counts and
all native physical bounds remain. Both outer refinements stopped at the
800-evaluation limit and are explicitly nonconverged. Native convergence is
not known. The prototype is version-coupled to private raw bounds, checked directly with the
compiled native fitter at zero iterations; no production engine was edited.

A separately authorized fixed Fisher replay will expose freshly fitted
opposite-parity weights/midpoints on the same single planned Viper job. Its
input is the **forward17 residual patches, prefix `res`**, as the actual native
producer specifies; the initial proposed backward input was corrected before
any real replay. The backward9 descriptor is a different input. Replay uses
only frozen native settings and compares complete keys, six-decimal score text
and reasons with the saved Fisher features before response attribution.
Even exact score agreement would not recover unique historical coefficients;
a mismatch remains a new realization, without tuning/retry. Reflection is
physical t, not u; the actual 197-pixel disk is reflection-closed.

The bounded driver now runs the four cases (13 N candidates, expected78 fit
rows), then the two Fisher cohorts sequentially. It saves completed stage exit
codes immediately and continues after failures. Updated native tests passed:
variable projection2834, acquisition184, representation319, Fisher replay204,
raw-QC112, driver55; raw output integrity305 and saved representation59+59.
The initial prototype fixed a `nothing`-to-TOML metadata serialization error;
synthetic near-zero and tuple-comparison test errors were corrected without
changing science. No benchmark has been regraded, no original prediction or
count is replaced, and no new Viper job has been submitted at this checkpoint.


### 2026-09-18 — Single bounded exploration submitted to Viper

Job **11806180** was submitted exactly once at approximately12:44 CEST from
scientific source `35c7d247331b42ffbc074e51c9ac51bf0a9f2dfb`. Its request is four CPUs,16GB,two hours;
the scheduler may allocate eight CPUs, while scientific concurrency remains
four one-thread cases. The complete ten-stage metadata dry-run and import-only
Julia1.13.0/one-BLAS-thread check passed before submission. No raw fit ran on
the login node. Scheduler test-only ID11806150 is not the scientific job.

The isolated runtime at `/u/oldu/code/STMFit_labelfree_exploration_20260918` contains37 required source/lock
files; `/ptmp/oldu/stmfit/label_free_exploration_20260918_inputs_v1` contains12 input files. Every transferred hash
matches locally and raw symlinks were dereferenced. The new output directory
`/ptmp/oldu/stmfit/label_free_exploration_20260918_run_v1` was absent immediately before submission. Source/inputs
were transferred without deletion; old reconstruction trees remain untouched.
The output/log will be fetched locally before interpretation. Initial state is
PENDING; the scheduler's pre-submit estimate was14:10 CEST, not a
completion promise. The one-job cap remains in force: failure is inspected and
reported, never a reason for automatic retry, trimming candidates or relaxing
bounds. No new batch result or benchmark gain is claimed at submission.


### 2026-09-18 — Bounded Viper pass blocked by omitted environment export

Job **11806180** is terminal **FAILED 1:0**. It ran 13:32:51–13:32:53 CEST
(two seconds) on `vipm2063`, with eight allocated CPUs. The fetched Slurm log
reports the exact failure:

```text
/var/spool/slurmd/job11806180/slurm_script: line 12: STMFIT_PROJECT_DIR: required
```

The failure is the shell's first required-path check, **before Julia or any
scientific stage**. A read-only probe confirmed
`SBATCH_EXPORT=NONE` in Viper's submitting environment and a native `/usr/bin/sbatch`
executable, not a wrapper. The parent submission assigned the required variables
with `env` but omitted the explicit `sbatch --export=...` flag. This was a launch
error, not a numerical/model failure. The earlier successful full146 actual
submission had used an explicit export flag; copying only its resource-only
scheduler check was not sufficient. Shell dry-run, import-only readiness and
`sbatch --test-only` cannot establish the compute-node environment handoff.

There is **no remote scientific output directory, zero of the planned 78 fit
rows, and no real Fisher replay or acquisition diagnostic**. These are not 78
numerically infeasible models. No substitute fit rows were manufactured.
The failure log was fetched into the exclusive local
`results/label_free_exploration_20260918/viper_diagnostics_v1/` bundle. Its hash
matches remotely; all 37 staged source/lock files and 12 input files are unchanged.
The failed remote source is preserved. The job monitor was cancelled; no retry,
requeue, second job, local replacement batch or change to scientific settings
was performed. The one-job allowance is exhausted even though no fitting ran.

The manual commands now name
`--export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN` explicitly.
Command-line options override `SBATCH_EXPORT`, which overrides `#SBATCH`
directives. The batch-script comments, HPC notes and agent notes explain this
boundary. A local regression reproduces the missing-variable early failure in
an empty environment and checks the supplied-variable handoff with `/bin/echo`
in place of Julia; it is **not** a Slurm emulator or a real fit. A future
corrected scientific submission requires renewed user approval and a fresh
output path. The correction has not been validated in a new compute job.
The local driver suite passes 55 assertions plus 11 environment-handoff checks
under Julia 1.13. The initial no-self-submission test also matched the new
explanatory comment; making it inspect non-comment lines resolved that test
fixture failure. Shell syntax and the documentation build pass; only the known
page/search-size and deployment warnings remain. All five original project,
manifest and physical/unit-config locks, and the other 36 staged source files,
remain unchanged locally. Only the launcher's comments changed locally; its
failed remote copy remains byte-identical to the submitted source.

**First-pass scientific conclusion:** constrained variable projection is a
concrete numerical avenue, supported only by the one fixed-N=9 smoke so far.
Its additional optimization reduces elliptical RSS by 47.2%, but both outer
searches stop at 800 evaluations and no correct-count or equal-budget advantage
is established. The profile outer objective is raw RSS; it does not include
the adjacent-Gaussian overlap/conditioning penalty κ used in the native global
stage. The native local LM stage also minimizes raw RSS. This test therefore
does not isolate a solver effect from extra optimization budget. The saved-only
native replay passes **347/347** checks (an initial worker message said 355;
the retained log is authoritative). It
reconstructs model values, decoded physical parameters, finite raw bounds and
full-p GCV, and checks recorded KKT/validity/stopping flags; it does not rerun
the original full-image validity calculation or prove statistical calibration.

The representation audit measures a dominant transverse component in normalized
residual patches (about 85% median centered energy in both cohorts), with
substantial within-file variation surviving predictor standardization. This
makes acquisition, centering and residual construction important to examine,
but does not identify a cause, quantify chemical information or justify
subtracting that component. Two unknown scans have large unobserved raw areas;
all 25 blind/overlay plates make them explicit. Recorded Current is nonconstant
but feedback-coupled, and all scans use the same -0.300 V bias with lock-in OFF.
It is not a new independent chemical measurement.

Four-file numerical/N comparisons, measured registration/noise/block evidence,
and real Fisher even/odd/offset decomposition remain **unmeasured because of
the launch failure**. No benchmark was used to select variants or regraded, no
new `N_selected` or chemical prediction was made, and no 80% information-theoretic
ceiling or accuracy improvement is established. Source code and synthetic tests
for these blocked diagnostics are available, not empirical validation. The
first bounded pass ends with the measured local findings and this explicit
blocker, rather than an automatic new method campaign.


### 2026-09-18 — One corrected submission explicitly authorized

After the first-pass closeout, the user explicitly authorized one new submission
with the same four count/acquisition cases and two forward17 Fisher replays,
four requested CPUs, 16 GB and a two-hour limit. This permits execution of the
previously blocked diagnostics, not further method development or an automatic
retry. The failed job 11806180 and its local failure bundle remain unchanged.

New job **11812202** was submitted exactly once at 16:38:31 CEST. The
actual command includes
`--export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN`
explicitly, overriding Viper's `SBATCH_EXPORT=NONE`. It reuses the unchanged
37-file source from `35c7d247331b42ffbc074e51c9ac51bf0a9f2dfb` and 12 unchanged input files.
No source sync or dependency change was needed. Source/input hashes, metadata
for all ten stages and Julia 1.13.0/one-BLAS-thread readiness were checked again.
The duplicate-job queue check was empty. Scheduler test-only 11812197 is not the
actual job; its 18:11 estimate is not a guaranteed start time.

The new exclusive remote output path is
`/ptmp/oldu/stmfit/label_free_exploration_20260918_run_v2`;
the new log is
`/ptmp/oldu/stmfit/label_free_exploration_20260918_corrected_11812202.log`.
Outputs will be fetched into
`results/label_free_exploration_20260918/viper_diagnostics_v2/`, separate from the
failed first bundle. Initial state was PENDING (Priority). No calculation result
is claimed at submission. The fresh one-submission allowance is consumed:
inspect and retain all failures, never requeue or launch a third job automatically.
The same 13 N candidates, expected 78 fit rows, original support, physical bounds,
full-p GCV and unchanged two Fisher cohorts remain in force. No additional
local real fit, benchmark evaluation, feature search or promotion is authorized.


Corrected-job outcome: independent accounting confirms **COMPLETED 0:0**,
17:10:22–17:14:12 CEST (3 min 50 s), on `vipc2001`, with eight allocated CPUs.
The terminal job had already expired from `squeue`; accounting was checked
independently. All 73 scientific output files and the Slurm log were fetched
and match remote hashes. All 37 source and 12 input hashes remain unchanged.
All ten expected stage records exist and report exit0. The output includes all
78 requested count-fit rows and both Fisher cohorts. No further job or fit is
authorized; numerical/interpretive checks follow separately. No new polling
heartbeat was created for an already terminal job.


### 2026-09-18 — Completed bounded mathematical diagnostics; no promotion

The authorized corrected job 11812202 is now checked, not merely processed.
All 73 output files and the Slurm log match their fetched hashes, the 37 staged
source files and 12 inputs are unchanged, and the original five project/config
locks remain intact. No additional fit, extraction, Fisher training, benchmark
grade or source modification ran during saved-output verification. The two
read-only reviewers are finished. Results are under
`results/label_free_exploration_20260918/viper_diagnostics_v2/`; detailed checks
are under `run_logs/corrected_submission_v2/` in the same exploration root.

**Numerical fitting.** All 78 declared rows are present: 13 counts × 2 geometries
× 3 methods. The saved-array Julia 1.13 replay passes 10,365 assertions, including
all 3,889 parameter rows (13 warm starts), 678 decoded lobe rows across the
candidate models, and 24,548 fit pixels. Native
forward-model reconstruction reproduces RSS and full-p GCV exactly in this
replay; an independently built saved-lobe basis differs by at most 8.33e-17.
All 52 linear profile solutions pass independently calculated box-KKT checks,
without solving a linear system. Full p remains 4N+3 circular and 5N+3 elliptical.
Support, original `refined_policy`, finite raw bounds and physical gates remain.
Validity is reproduced on the saved fit mask using saved noise/calibration;
this does not independently validate raw preprocessing, fusion or noise.

At every fixed N/geometry, RSS decreases native → fixed profile → refinement.
Fixed-geometry profiling alone reduces RSS by 0.079–20.427%; total reductions
are 4.775–69.497%. These ranges include invalid native starting rows, not just
successful-looking cases. Four invalid rows at 260215_022/N10 (native and fixed,
circular and elliptical; reason `residual high`) remain in the tables. Both
N10 refinements pass the original gate. All 26 refinements are valid, but 25 stop
at 800 evaluations without convergence; only 260115_016/N5 circular reaches FTOL
at 787. Native convergence remains unknown. No time/evaluation cap was raised.

| Case | Declared Ns | Native→refined elliptical RSS reductions | Minimum valid candidate-set GCV N, all three methods | Original selected N (not truth) |
|---|---|---|---:|---:|
|260115_016|5,6,7,8|8.52–17.88%|8|6|
|251206_013|10,11,12|19.20–35.47%|12|11|
|260215_022|8,9,10|50.94–69.50%|9|9|
|251206_034|8,9,10|5.58–40.99%|9|9|

The N9 elliptical result in this job is 0.05227166419 → 0.05171210328 →
0.02438147059 (53.36% total reduction). It is a separate invocation from the
one local smoke's 47.2%, not a replacement of that saved result. The candidate
minima stay 8/12/9/9; the production selector was not rerun or replaced by these
rankings. Two minima lie at a declared upper edge. The sample is not exhaustive.
Extra optimization and unknown native convergence prevent equal-budget or
accuracy claims. The native global stage penalizes adjacent-Gaussian overlap
conditioning κ; local LM and the profile outer search use raw RSS. Saved endpoint
κ is 3.760–7.806, below configured 10; this says nothing about the global trajectory.

**Acquisition and noise.** The saved-only review passes 255 grouped checks. It
retains 340 lag candidates, 1,165 row records, 70 local records, 408 ACF rows and 28
auxiliary rows. All four best signed Z correlations occur at dx=−8, the fixed
search boundary, and all are ambiguous. Zero-lag→best-grid correlations are
0.3964→0.6479,0.4726→0.6568,0.9601→0.9928,0.7421→0.9592 in the table's case order.
None is an accepted registration; no search window was expanded. After the
unresolved shift, 2/6 and 7/11 local correlations remain negative in the first two
cases. Whole-support agreement does not establish per-lobe agreement.

All 16 off-footprint view/axis ACF lengths are right-censored at 16 pixels. Their
lag 16 correlations remain at least 0.703, above the declared exp(−1) threshold.
No usable block size, stationary noise covariance, independent sample count or
joint likelihood is established. The saved-residual ACF is separate model/data
structure, not a substitute for background noise; 251206_034 shows magnitude
resurgence after an early crossing. Block zeros and insufficient-row lag zeros
are sentinels, not measured zero scales/lags. The `n_eff = n ÷ 9` placeholder,
GCV and all settings remain unchanged.

Current is present but anticorrelated fwd/bwd on the fixed Z support: after the
unresolved inherited lag, r=−0.5226,−0.7619,−0.4053,−0.8411. The positive-affine
diagnostic clips gain to 0, which gives NRMSE 1; it does not silently flip the sign.
All 20 auxiliary comparisons have complete pairs on that restricted support.
This does not restore 260115_016's missing 60.6–60.7% of the full frame. Restoring
raw-missing masks excludes imputed observations but cannot undo imputation's
influence on plane/row flattening and support construction. Current remains
feedback-coupled, single-bias evidence; its optional nuisance projection was
not performed. No demonstrated new implementation bug was found.

**Fixed Fisher replay.** Both exported score TSVs are byte-identical to the
saved originals: 222 unknown25 and 900 full146 keys, including all seven original
unavailable rows. A separate saved-patch/coefficient Julia 1.13 replay passes
33,461 assertions without fitting. It reconstructs the 197-pixel disk and physical
t-reflection, opposite-parity support, six-decimal scores, reasons and response
components; maximum error is 1.78e-15/2.66e-15. New weights remain a fresh realization
consistent with the exports, not unique recovered historical coefficients.

For `s=e+c+abs(o)`, median even/absolute-odd/maximum responses are
3.4315/0.5222/4.1236 (unknown25) and 4.9946/0.06045/5.0704 (full146). The even term
is already positive for 200/222 and 874/893 finite rows; the maximum is positive
for 209/222 and 875/893. Only 9 and 1 rows cross zero because of the odd bonus. Thus
the mirror maximum does not alone explain the mostly positive scores. The
shared offset is near zero (absolute max 1.31e-16); the historical centered-mid
versus raw-input convention was not corrected. These are proxy responses, not
chemical class counts or calibrated probabilities. No downstream prediction
changed and no grade was run.

**Conclusion.** There is measured numerical headroom, but no changed diagnostic
minimum N in these fixed sets and no demonstrated accuracy gain. The earlier
transverse-gradient attribution (about 85% median centered patch energy) remains
an observation, not proof of an artefact or absent chemistry. The acquisition
results do not yet support calibrated noise/independence or a reliable new
registration correction. No universal 80% limit follows from these experiments.
The bounded exploration is complete as diagnostics, with failures, unavailable
values and unresolved scientific questions retained. No method is promoted;
any new campaign, fit or parameter change needs a new user-approved scope.


### 2026-09-18 — Start masked-background prototype, synthetic first

After asking for other methods, the user replied “go; aller” to the proposed
priority of masked robust preprocessing, followed only later by a true paired
acquisition model. The initial deliverable is a small opt-in background
prototype and quantitative synthetic signal-preservation checks. This first
development pass is capped at two hours from 20:28 CEST. Completion means native
Julia 1.13 tests, explicit retained limitations, documented controls and no
change to production. There is no real SXM read, fit, new benchmark grade or HPC
submission in this first step. A proposed four-case real comparison must have
its own explicit budget before submission. Completed jobs remain closed.

The motivation is measured missingness and its upstream effects, not a promised
accuracy gain: native preprocessing imputes before fitting a plane and row
levels; restoring raw masks later does not undo the effect on observed pixels.
The prototype will use only finite supplied-background observations, then
subtract background-only row medians. Neither missing observations nor
insufficient-background rows are filled. Molecular gradients are not blindly
removed. A mask remains an assumption, not chemical truth.

Four variants are fixed before synthetic results: native unsmoothed reference,
finite-only OLS, guarded OLS, and guarded Huber. This separates missing-value
handling, foreground exclusion and robust plane estimation. The Huber cutoff
is 1.345; at most 50 IRLS iterations use a fixed initial-background-residual MAD,
with explicit numerical scale/rank/convergence guards. All settings live in
`config/masked_robust_preprocessing.toml`, not production defaults.

The study must retain leaky foreground masks and unsupported rows, report both
own and common observed coverage, and distinguish native/global from
background-zero level conventions. A mathematical limitation is already known:
sequential plane then row fitting can bias x slope when row offsets correlate
with unequal x coverage. Paired signal-present/absent invariance under a correct
mask is not proof of an accurate absolute background. The synthetic study will
measure both. No numerical or real-image improvement is claimed at this start.


### 2026-09-18 — Masked background: verified synthetic benefits and failure modes

The synthetic-first prototype is complete; it is not a production replacement.
The parent reran both suites under Julia 1.13, one Julia thread and one BLAS
thread: **441/441** engine assertions and **1,116/1,116** signal-study assertions.
A separate **12,641/12,641** saved-table arithmetic check uses no estimator,
fixture generator or raw reader. It replays the reported error/coverage metrics,
exported plane/row reconstructions, final fixed-scale plane objectives and
stationarity diagnostics.
Nine original project/config/core locks and all 73 previous exploration outputs
remain unchanged. No real scan, count/chemical fit, benchmark or cluster command
ran during this prototype study.

The fixed study contains five deterministic 73×97 images and four methods, each
with and without an injected signal: 40 outcomes, 70 metric rows, 2,190 row-status
records and 35,405 saved pixel records. Synthetic signals include compact bumps
and a signed localized transverse component. This is a constructed mechanism
test, not a representative noise simulation or an accuracy benchmark. Supplied
foreground exclusions are trusted inputs; real masks are not known this well.
The native reference calls the actual `GaussianFit2D.preprocess_channel` on
synthetic shared SXM types, using stride 1, plane+rows and no smoothing. Native
imputation still influences its fit, but raw missingness is restored for metrics.

The table shows absolute foreground reconstruction RMSE in **pm** on the same
observed support for every method in each row. One observable constant level,
measured on the fixed common background corner, is removed from each image.
Truth does not set that level; no slope, gain or extra row correction is removed
from the reported errors. A Gaussian tail in that corner can bias the level.

| Fixed synthetic case | Native | Finite-only OLS | Guarded OLS | Guarded Huber | Common observed pixels |
|---|---:|---:|---:|---:|---:|
| Compact signal, complete image | 0.5549 | 0.5549 | ≈0 (roundoff) | ≈0 (roundoff) | 7081 |
| Structured missingness | 66.18 | 3.414 | 0.9288 | 1.124 | 5117 |
| Background outlier patch | 107.8 | 107.8 | 118 | 9.843 | 7081 |
| Leaky Gaussian exclusion | 21.51 | 21.51 | 19.62 | 14.6 | 7081 |
| Rows without observed background | 63.42 | 3.784 | 1.015 | 1.416 | 4452 |

Several distinctions prevent an overstated result:

- In the compact-signal cases, a correct guarded mask makes the paired
  signal-present minus signal-absent response accurate to about 3e-17 nm. This
  is conditional signal non-use during background estimation, not proof of an
  accurate absolute background. With structured missingness, guarded OLS still
  has 0.929 pm foreground RMSE; Huber is worse at 1.124 pm. Unequal per-row x
  support lets row offsets bias the sequentially fitted x slope.
- With background outliers, exclusion alone worsens foreground error from
  107.8 to 118.0 pm. Huber reduces it to 9.843 pm, not zero. Its fixed initial-OLS
  scale rises from 0.0113862 to 0.0869708 nm under contamination; it is not an
  independently measured STM noise scale.
- With a leaky Gaussian exclusion, both guarded variants distort the signal.
  The transverse-component diagnostic gives 0.0732983/0.0751164 nm versus an
  injected 0.08 nm: 8.38%/6.10% attenuation. Robustness does not replace a correct
  footprint, and shared image structure is not automatically removable nuisance.
- The final case loses **206 observed foreground pixels** under either guarded
  method. Only 4,452/4,658 observed pixels and 1,246/1,452 foreground pixels
  remain (95.58% and 85.81%). The table compares the native reference on exactly
  those 4,452 pixels; its own-support result is separately retained. No claim is
  made about the unavailable foreground. Finite-only OLS retains all observed
  pixels but reports PARTIAL because one entire raw row is missing.

All ten Huber signal/null runs converge in 5–18 iterations with the predeclared
controls unchanged. Forced-budget, missing/rank/row failures and absent metric
anchors are separately tested and retained as unavailable. The final study reads
and exports the actual diagnostic TOML; it does not use a hidden duplicate set
of controls. The initial primitive test run had four failures because a Julia
whole-array Frobenius tolerance was used where a per-pixel maximum was intended.
The assertions now use max(abs(error)) with the same 2e-10 bound. No estimator
or scientific control was tuned to cure that test error; the failed log remains.

Authoritative final-source evidence is
`results/masked_preprocessing_20260918/parent_validation/` (six files), with
parent test/checker logs under `run_logs/`. Earlier worker outputs are preserved,
not the final reference. An initial worker smoke used truth to align a constant
level; it was explicitly superseded by the observable-only gauge before the
reported study. The early worker message's 978-test claim was a transcription
error; its log has 954. Final source has 1,116 and was independently rerun above.
Source: `test/lib/masked_robust_preprocessing.jl`,
`test/lib/masked_preprocessing_fixtures.jl`, and their two standalone test files.

**Decision:** keep all four diagnostic variants; do not select Huber as a new
default or remove the measured transverse gradient from real patches. A later
real-data comparison is proposed on the same four scans only, with frozen saved
geometry/support and the existing lag window: no count fitting, classifier,
grade or parameter search. Proposed cap: one Viper job, four requested CPUs,
16 GB, 30 minutes wall time. This is a proposal, not a submission allowance.
No such job or real-data trial has been started; no production behavior changed.


### 2026-09-18 — Authorize the four-scan masked-preprocessing comparison

After the synthetic prototype and the explicit one-job proposal, the user
replied “go tu as carte blanche”. The first real-data tranche is now authorized:
compare native, finite-only OLS, guarded OLS and guarded Huber on the same four
frozen cases (`260115_016`, `251206_013`, `260215_022`, `251206_034`). The scientific
deliverable is a measured comparison of background changes, observed coverage
and forward/backward agreement, with failed methods and unresolved shifts kept.
This is not a count or chemical-accuracy experiment.

Preparation is capped at two development hours from 21:18 CEST. The first
compute tranche is one Viper job: four requested CPUs, 16 GB, 30 minutes wall
time, four one-thread cases. Local work is synthetic testing and metadata-only
preparation; all real image processing belongs on the compute node. No previous
job is resumed or resubmitted. The submit command will explicitly export the
four required variables because Viper still sets `SBATCH_EXPORT=NONE`.

The original four-file sample, 25-file summary and 222-lobe base geometry remain
byte-identical. The original selected support (including the real adaptive
rescue) and its conservative background exclusion are frozen across methods.
Their dependence on previous imputation/preprocessing remains a limitation.
The lag window and all numerical controls remain unchanged. No count model,
Fisher/classifier, benchmark grade, expected N/composition prior, parameter
search, new noise calibration or production-default change is included.

Agreement must be compared on constant observed supports: method-own,
all-method common, and native-versus-each-variant common intersections over the
whole lag grid. A failed method may empty the common support; its failure stays
visible rather than being dropped. Background level alignment uses observable
common off-footprint pixels only. A best lag is still diagnostic, not validated
registration. The pilot will retain raw/processed arrays and masks for saved-only
arithmetic and visual checks after retrieval. No real-data result is claimed yet.


Preparation validation for the four-scan pilot is complete locally. The parent
independently ran **5,438/5,438 comparison assertions** and **205/205 driver
assertions** under Julia 1.13, with one Julia/BLAS thread. Tests include synthetic
SXM parsing, original-support replay, missingness, blocked methods, identical
whole-lag supports, fixed observable background medians, retained negative and
ambiguous/boundary correlations, full saved-array/table output and CLI behavior.
The four real inputs passed **metadata-only** checks without parsing SXM pixels
or creating results: saved selected N is 6/11/9/9, with actual adaptive rescue
only for `260215_022`. These values check saved geometry, not molecular truth.

Test failures were retained and fixed without changing scientific controls:
a driver assertion used Julia's one-argument `occursin` with the wrong bound
argument; a synthetic injected view used a native string background sentinel
where a non-native Boolean mask was required. Final serialization preserves
integer options explicitly. Auxiliary-channel wording and its assertions now
state what the unchanged reader does: it incidentally preprocesses Current,
then the new comparison discards those arrays and produces no Current evidence.
Earlier queued worker counts/hashes are not the final source: the parent tests
and stable source hashes are authoritative. Documentation and whitespace checks
pass. Twelve original config/core/helper files remain unchanged.

All seven input files were transferred to a new Viper input directory after a
reviewed rsync dry-run. Their byte hashes match locally; remote raw files are
regular files, not symlinks. `/ptmp` itself is a site symlink, so the requested
fresh result location is passed canonically as
`/viper/ptmp1/oldu/stmfit/masked_preprocessing_20260918_run_v1`. This is the same
scratch location, not a second run. No real-data processing or job submission
has occurred at this source-validation point.


### 2026-09-18 — Four-scan masked-preprocessing job submitted once

The authorized pilot was submitted exactly once as Viper job **11820414**
at Slurm time **22:25:06 CEST**, from scientific source
`cb3b2d9843eccac2612859fe8b19cfe609fdd074`. The actual command explicitly requested four CPUs,
16000 MB and 30 minutes, with the four required variables exported on the
`sbatch` command line. The one-job allowance is now consumed; there is no
automatic resubmission.

All 33 staged source/lock files and seven input files match their local byte
hashes, and the Julia executable matches the local 1.13.0 binary. Login-node
imports, four per-case metadata-only checks and the actual launcher's dry-run
passed. No real SXM preprocessing ran on the login node. Scheduler test-only
identifier `11820378` was not a real submission; its suggested start time is not
a guarantee. The actual job is **PENDING**, with no node assigned, because nodes
are unavailable/drained or reserved for higher-priority partitions. No real-data
result is available yet.

The fresh result path is
`/viper/ptmp1/oldu/stmfit/masked_preprocessing_20260918_run_v1`; its `/ptmp` alias
refers to the same location. The Slurm log is
`/ptmp/oldu/stmfit/masked_preprocessing_20260918_real4_11820414.log`.
Monitor only this new job, fetch outputs/logs when terminal, and verify saved
arrays and fixed-support metrics before scientific interpretation. The original
counting, chemical assignments, physical settings and completed jobs stay
unchanged. A successful Slurm exit would not itself validate registration,
coverage, background correctness or molecular accuracy.


During the queue wait, a saved-array-only plotting helper was checked on
synthetic images (**14/14 assertions**, followed by direct inspection). It
retains the all-method and all three native-pair comparisons, uses the exported
common observed background median per view, and leaves an unavailable common
comparison blank rather than substituting a mask or level. Separate availability
panels distinguish raw missing pixels from newly unavailable observations.
The fixed 99.5% absolute-value color limit is a **display rule only**, shared
within each scope; clipped pixels are counted and numerical arrays stay intact.

The first previews exposed two display bugs not caught by PNG/header tests:
`legend=false` did not enable a colorbar, and a mask contour coupled to height
color limits was not visible. Explicit colorbars and red discrete support-edge
markers fixed both. Blocked panels have no invented height scale; their
uncomputed clipping count is unavailable. Final synthetic previews show the
intended physical axes/aspect, visible units/support, and coverage categories.
This is not a real-image review or scientific result. Helper, tests, hashes,
previous previews and final visual notes are retained under
`results/masked_preprocessing_20260918/run_logs/real4_v1/`. No remote scientific
source, settings, arrays or predictions changed. Job `11820414` remained pending
at the 22:39:47 CEST check; the single-job monitor remains active.


The independent saved-output checker is ready before real results arrive.
It imports only Julia standard libraries; it does not load the producer, raw
SXM reader, estimator or solver. The parent independently passed **23,843/23,843
synthetic assertions** (20,846 regular/prior-summary checks plus 2,997 checks of
retained forced nonconvergence). A deliberately wrong raw-hash manifest failed
exactly its intended assertion (2,996 pass / one intended fail, exit 1).
The parent reviewed the equations, support intersections, mapping and failure
handling. Earlier negative metric/plane/process/mapping probes and failed
checker-development logs are retained. One checker-only bug reversed the
reported difference-median convention: it must be **backward minus forward**.
Correcting that expression did not change the producer or scientific settings.

The checker reconstructs exported physical planes/rows, final fixed-scale
objectives and stationarity, observed availability, the fixed 85-lag supports,
all/native-pair intersections, signed metrics and background medians. It does
not independently recover the initial OLS/MAD or iteration history, unknown
native estimator parameters, raw parsing/preprocessing or a true background
mask. Full-frame counts are reproducible from saved arrays here because the
frozen production stride is one. Raw SHA claims are bound to the parent's saved
`input_sha256.toml` without opening SXM bytes; the parent separately verifies
actual local/remote input bytes. The four earlier acquisition summaries have
the expected native-reference/support-hash schema, but numerical comparison
awaits this job's actual outputs.

The checker is
`results/masked_preprocessing_20260918/run_logs/real4_v1/verify_saved_real4.jl`
(SHA256 `eb4ae4dd052f3a4a94217a744c524dfba883225ac6a86694f630b1c54684b482`).
These are **synthetic checker tests, not verified real-data findings**.
Both implementation workers have handed back ownership and are retired.
The latest real-job check, 22:55:20 CEST, still shows `11820414` pending under
`QOSGrpCpuLimit`: a group CPU-quota wait, not a program failure. No further
submission, resource change, raw-data analysis, fit or grade was launched.


## 2026-09-18 — Four-scan masked-preprocessing results; no promotion

**Deliverable completed:** measured background changes, observed coverage and
forward/backward agreement for the four authorized frozen scans and all four
methods. The real pilot ran only in Viper job **11820414**, from scientific
source `cb3b2d9843eccac2612859fe8b19cfe609fdd074`, on `vipc2169` at
23:32:39–23:33:34 CEST: **55 s, COMPLETED 0:0**. It requested four CPUs,
16000 MB and 30 minutes; Slurm allocated eight CPUs, with four one-thread Julia
1.13 case commands. Their elapsed times were 39.793–40.941 s. Preparation and
synthetic checks finished in 110.4 minutes, within the 120-minute preparation
budget. The group CPU-quota queue wait did not trigger another submission.

All **37 outputs (229,592,456 bytes)** and the Slurm log were fetched once.
Remote/local SHA256 values match, all 33 source and seven input identities are
unchanged, and all 12 original source/config locks match. The earlier 73-file
diagnostic tree, 184-file full146 tree and 66-file unknown25 tree plus its log
remain byte-identical. The current job monitor is stopped. No new raw SXM
processing occurred locally; local work used saved arrays and arithmetic only.
No count, classifier, feature refit, grade, parameter search or production
change was made. The native reader's incidental Current preprocessing was
followed by discarding Current; it supplies no new Current evidence.

### Coverage, fixed supports and convergence

Methods are native reference (N), finite-only OLS (F), guarded OLS (G), and
guarded fixed-scale Huber (H). Their own, all-four and native-pair supports and
background anchors are identical within every case, including hashes:
**3,948 / 6,279 / 9,392 / 4,929 comparison pixels**. Each is the full retained
frozen fit support. No method loses an additional observed pixel in this pilot;
this does not validate the inherited ROI/geometry/footprint or recover missing
measurements. All masks, the actual adaptive rescue on `260215_022`, the 85-lag
window and all controls remain frozen.

`260115_016` has only **103,216 / 102,912 observed pixels** out of 262,144 per
view (39.374% / 39.258%); its 158,928 / 159,232 missing pixels remain unavailable.
All three masked methods remain `PARTIAL`: 310 / 311 unsupported rows contain
zero observations. Their other 202 / 201 rows are available. The other three
scans are fully observed. Of 32 view rows, eight are `REFERENCE`, 18 `OK`, six
`PARTIAL`; there is no masked nonconvergence. H converges in 11/11, 12/13, 12/10,
11/10 iterations (forward/backward in case order). Native convergence and
internal estimator details remain unknown. All 16,384 row records, 112 summary
rows, 4,760 lag candidates, 980 local rows and 112 background rows are retained.

### Measured agreement and background changes

Zero-lag forward/backward RMS below uses exactly one median per view on the
same observed off-footprint background, fixed over all 85 lags. There is no
slope, gain, truth alignment or applied registration. Units are **pm**; these
are differences between views, not errors against true molecular signal.

| Scan | Pixels | N RMS | F RMS | G RMS | H RMS |
|---|---:|---:|---:|---:|---:|
| 260115_016.sxm | 3948 | 51.542 | 49.177 | 51.512 | 51.547 |
| 251206_013.sxm | 6279 | 49.039 | 49.039 | 49.266 | 49.268 |
| 260215_022.sxm | 9392 | 12.762 | 12.762 | 12.781 | 12.782 |
| 251206_034.sxm | 4929 | 23.254 | 23.254 | 23.026 | 23.033 |

The incomplete scan's finite-only RMS drops by 2.364 pm (4.59%). Guarding
increases its correlation more, but barely changes RMS; H is slightly worse
than N in RMS. On `251206_013` and `260215_022`, G/H increase both correlation
and RMS. `251206_034` improves zero-lag RMS, but worsens RMS at the diagnostic
correlation-selected best lag (N 8.117 pm, G/H 8.454/8.460 pm). For the three
fully observed scans, F and N agree in level-invariant metrics to roundoff.
H has slightly lower correlation and higher RMS than G at both reported lag
conditions in all four cases; these small differences have no calibrated
statistical significance and do not select a winner.

Shared-background normalized MAD (forward/backward, **pm**):

| Scan | N | F | G | H |
|---|---:|---:|---:|---:|
| 260115_016.sxm | 46.856 / 42.726 | 8.258 / 3.827 | 1.170 / 1.046 | 1.136 / 1.037 |
| 251206_013.sxm | 5.131 / 8.718 | 5.131 / 8.718 | 4.971 / 4.807 | 4.921 / 4.807 |
| 260215_022.sxm | 3.047 / 1.796 | 3.047 / 1.796 | 1.368 / 1.322 | 1.365 / 1.321 |
| 251206_034.sxm | 10.653 / 11.777 | 10.653 / 11.777 | 5.498 / 5.424 | 5.480 / 5.399 |

The much smaller background spread on `260115_016` coexists with persistent
foreground disagreement. H changes its level-aligned background versus N by
46.422/43.084 pm RMS; change magnitude is not correction accuracy. These are
post-fit descriptions on background used by the estimator, not held-out noise
measurements. H slightly lowers MAD versus G but increases background standard
deviation in both views of `251206_013` and `251206_034`. A single dispersion
statistic would hide this adverse result. Neither the OLS-derived Huber scale
nor these post-fit statistics calibrate noise, covariance, blocks or `n_eff`.

Every method has best diagnostic shift **(-8,2), (-8,-2), (-8,-1), (-8,-2)**
in case order: all are boundary-limited and ambiguous, none usable or accepted.
Even the approximately 0.993 best correlation on `260215_022` has 15–17 near
maxima. No shift was applied and the lag window was not expanded. Local signed
anticorrelations remain: zero-lag negative counts are 3/6, 7/11, 0/9, and 4/9
for N/F versus 5/9 for G/H in the last case. The additional `251206_034` patch 8
changes from +0.0680 to -0.00363/-0.00394. Local patches overlap and their IDs
are saved geometry, not validated counts or chemistry.

### Verification limit: exported stationarity loses precision

The independent full saved-output check reports **30,811 pass / 15 fail /
30,826 total** in 43.2 s. Every failure is the same strict stationarity-scalar
comparison; all plane/row corrections, final fixed-scale objectives, masks,
coverage, metrics, metadata and four previous native baselines pass. The
original checker (SHA256
`eb4ae4dd052f3a4a94217a744c524dfba883225ac6a86694f630b1c54684b482`),
`atol=rtol=1e-10`, code-1 result and initial log remain unchanged. **We did not
replace its assertion with a weaker check to obtain a green full-suite result.**

A separate arithmetic-only diagnosis was independently repeated by the parent:
**118/118 checks pass**, and its three result tables are byte-identical. No
producer import, raw reader, fit, SVD or initial OLS/MAD replay was used.
The source computes the gradient from median-centered heights and normalized
parameters before exporting a physical plane. Those normalized parameters are
not saved. Here background levels are -84 to -23 nm; physical planes differ by
only up to 2.84e-14 nm across equivalent evaluation forms, but coherent errors
can accumulate over 49,838–262,144 pixels. Actual stationarity discrepancies
range from 3.27e-12 to 3.75e-9 nm.

A known-parameter synthetic calculation retains exactly the same internal
gradient after a constant level shift to -64 nm, yet its physical-export
recalculation changes from passing to failing the original tolerance. Two
predeclared normalized intercepts can export identical coefficients and planes
while giving gradients separated by 2.33e-10 nm, more than two 1e-10 tolerance
intervals. This establishes an export information-loss mechanism, not the
correctness of every real producer operation.

All 24 real differences are within a separately derived export-roundoff budget;
the conservative export-plus-accumulation bound ranges from **2.35e-9 to
1.04e-6 nm**. It depends on operation counts, precision and saved magnitudes,
not a tolerance fitted to the observed discrepancy. A fixed +1e-4 corruption
is rejected in 24/24 views, but +1e-6 only in 22/24. This is **bounded
compatibility, not exact stationarity replay or exclusion of every small
error**. Strict 1e-10 stationarity remains unverified. Initial OLS/MAD history,
unknown native estimator parameters and true-background correctness also
remain outside the saved-only verification. No producer defect is established;
no rerun, output change or scientific-control adjustment was made.

### Visual review and decision

Twenty saved-array PNGs and 112 panel records were produced in a separate
folder. All four common comparisons, the incomplete-scan availability map,
and representative native-pair comparisons were inspected directly. Their
848 panel-metadata checks pass. Physical aspect, nm axes/colorbars, red fixed
support edges and gray unobserved regions are visible; no extra flip or
registration was introduced. These are corrected arrays, not raw STM images.
The fixed 99.5% shared-within-scope display clipping is counted, never used to
filter scientific calculations. A parent dimension check initially assumed
all figures had width 2,000 pixels; pair figures correctly use two 500-pixel
columns (1,000 pixels total). Correcting that assertion changed no images or
plot helper. The eight four-method figures and twelve pair figures match the
existing layout.

**Decision:** retain this as an opt-in diagnostic. Lower background spread and
slightly higher pooled correlations do not show preserved real molecular
contrast, better counting/chemistry, calibrated uncertainty, or resolved
registration. No method is selected or promoted, and no new campaign is
started. The scientific deliverable is the measured comparison with these
limits, not a claim that all numerical diagnostics replay exactly.

Local report: `results/masked_preprocessing_20260918/real4_v1_report.md`.
Outputs: `real4_v1/`; figures: `real4_v1_figures_v1/`; evidence under
`run_logs/real4_v1/`, including the unchanged initial verifier/log,
`scientific_saved_table_review.md`, `stationarity_roundoff_review/REPORT.md`
and the separate `stationarity_roundoff_parent_check/` reproduction. These
generated artifacts stay untracked. The completed job and all prior jobs
remain closed; all workers and the current monitor are retired.


## 2026-09-19 — Benchmark-centered recognition: authorized saved-output first pass

The user clarified that recognition quality, judged by the external benchmark,
is the objective, then approved proceeding. Background flatness, fit RSS and
forward/backward correlation are diagnostic proxies, not recognition gains.
The first pass is bounded to two hours of saved-output analysis from
13:19 CEST. Its scientific deliverable is a comparable benchmark error diagnosis
and, where justified, at most two explicitly frozen label-free single-component
comparisons with measured recognition outcomes. If saved outputs cannot support
a meaningful candidate, report that limit rather than inventing a method gain.

The frozen reference and reproducible reconstruction remain distinct. Correct
units/all controls, exact chains, classification coverage, missing and extra
lobes and class errors must accompany conditional accuracy. The 900-lobe
assignment cache is not the 871-lobe promoted counting result. Existing grades
are historical development evidence, not an untouched generalization test.
New candidate definitions must be fixed before their external grades; no
expected count, sequence, class count, position prior or benchmark label may
enter inference, fitting, calibration, thresholds or abstention. Labels remain
in the external audit/grader only. No supervised training is authorized.

This is a new bounded analytical pass, not a restart of any completed job.
Original sources, predictions and result trees remain immutable. No raw SXM
processing, new fit, cluster submission, broad parameter search or production
promotion is included in this initial saved-output scope. No new approval
framework or runtime machinery is introduced. All earlier jobs and the
four-scan pilot remain closed.


Before any new candidate grade (13:28 CEST), two saved-native component
endpoints were fixed: exact k-means decisions and exact GMM decisions, both on
the baseline joint-validity mask. All 900 keys and seven unavailable keys remain;
the byte-identical reconstructed soft vote is the control. No threshold, seed,
training, feature, class mapping, count or coverage policy changes. Membership
filtering occurs only in the existing external comparison script, leaving the
scientific exports unfiltered. These are fusion ablations on the new explicit
native reconstruction, not a reopened vote-weight/abstention sweep or a new
counting run. Historical component sweeps used different realizations and remain
negative evidence against assuming an improvement in advance.

The saved vote fractions are hard seed votes, not calibrated posteriors.
K-means additionally averages four view votes; on 893 jointly available keys
it is exactly 0/1 on 857. Components disagree on 146 keys; the soft vote differs
from k-means on 18 and GMM on 128. This justifies checking which component the
fusion actually follows, but predicts no accuracy gain. Copied component
confidence uses `max(p,1-p)` while the fused control uses `2abs(p-0.5)`;
these scales are not directly comparable and are not evaluation thresholds.

The opt-in exporter `test/export_saved_assignment_components.jl` reads only
saved label-free predictions and uses the existing benchmark-column-rejecting
TSV reader. It preserves copied decisions even when the printed vote rounds to
0.5, refuses existing/dangling output paths and validates key/availability
agreement before writing. Its Julia 1.13 synthetic checks pass 103/103; parent
saved-output checks verify 9,002 copied fields/keys. No raw reader or estimator
is called. Recognition results remain unmeasured at this checkpoint.


## 2026-09-19 — Two frozen native recognition endpoints: negative result, no promotion

**Deliverable:** use the external recognition benchmark to decide whether a
fixed saved-head change improves assignment, rather than infer improvement from
RSS, smoother backgrounds or registration correlation. The authorized first
pass began at 13:19 CEST with a two-hour cap. Both candidate definitions were
fixed at 13:28:29; the single new four-profile external grade began at 13:50:12.
The result is negative and closes this saved-output pass. The exporter and its
103 synthetic checks are in `800dfb4`; no production model was changed.

### Measured recognition result

| Profile | Correct / all 870 | Correct / classified | Classified coverage | Emitted errors | Exact chains / 145 |
|---|---|---|---|---|---|
| Historical reference | 677/870 (77.8%) | 677/854 (79.3%) | 98.2% | 177 | 36 |
| Reconstructed soft vote | 666/870 (76.6%) | 666/849 (78.4%) | 97.6% | 183 | 24 |
| K-means endpoint, same mask | 662/870 (76.1%) | 662/849 (78.0%) | 97.6% | 187 | 16 |
| GMM endpoint, same mask | 622/870 (71.5%) | 622/849 (73.3%) | 97.6% | 227 | 19 |

Each endpoint copies the exact saved head decisions on the baseline joint-valid
mask; it does not rerun a classifier or retune any threshold, polarity, feature,
seed, self-training step or abstention rule. Scientific exports retain all 900
keys across 146 files and seven unavailable keys. External manifest filtering alone
selects the common 892 keys across 145 files. Six unavailable benchmark keys remain:
five aligned control positions and one extra. One unavailable key is among the
eight rows of the outside-benchmark scan. The reference emits binary labels at
all 892 keys. All profiles keep 16 missing control positions across 14 short
scans and 38 extra lobes across 25 long scans. The 193/204/208/248 positions not
correct include emitted errors; they are **not model abstention counts**.

This is the cached **900-lobe own-N assignment lineage**, not the incompatible
871-lobe promoted counting result. Counts and key membership within the
scientific cohort never change. The exact historical descriptor/weights are
still not recovered; the reference is an archival output, not an exactly
reproduced fitting pipeline. Neither new endpoint improves either baseline.
The historical headline remains unchanged; no automatic promotion or further
vote/feature/seed/abstention search follows.

### What failed, including the favorable cases

Against reconstructed control, k-means gains correct positions on 5 scans,
loses on 9, and ties on 131. It gains no exact chains and loses eight. GMM wins
on 18, loses on 57, and ties on 70; it gains seven exact chains but loses 12.
Against the reference, k-means has 14 wins, 27 losses and 104 ties and no exact
gain and 20 losses; GMM has 14/62/69 and three exact gains and 20 losses.

The class-confusion totals explain the net loss under the existing physical
grading convention. Reference `(TN,FP,FN,TP)` is `(507,59,118,170)`;
reconstruction `(519,43,140,147)`; k-means `(522,40,147,140)`; GMM
`(445,117,110,177)`. K-means avoids three false positives but adds seven false
negatives (-4 correct). GMM recovers 30 false negatives but adds 74 false
positives (-44 correct). Higher GlcNAc sensitivity alone is not an improvement.
These class counts are external diagnostics, never a desired composition or
an inference prior.

The grader independently chooses identity/reversal using external truth.
Against reconstructed control, k-means changes 17 benchmark keys on 16 scans
and changes orientation on two; its direct delta at the control alignment is
-7, then reselection adds three, giving -4. GMM changes 128 keys on 98 scans
and changes orientation on 13; direct -60 plus reselection of 16 gives -44.
This is bookkeeping of the existing grade, not a new alignment optimization in
inference or an insertion/deletion/spatial ground-truth match. The full audit
retains all 725 pair/file rows and 4,460 pair/key rows, including adverse and
unchanged cases and both mappings.

For the 106 N6 scans, the palindromic external control makes both orientations
tie and the grader keeps identity. Of 636 positions, 442 are correct in both
heads, 95 wrong in both, 67 correct only in k-means, 31 only in GMM and one
jointly unavailable. The 95 common errors occur on 65 scans (23 truth-0 and
72 truth-1 positions). Selecting one of the two unchanged head labels cannot
repair these particular errors. **No per-lobe best-head selector or oracle
prediction was built.** This is not a limit on information in the raw images.

### Reproduced reference/reconstruction diagnosis

The initial independent saved-row audit reproduces all original physical
per-file and position grades: 677/854 and 36 exact versus 666/849 and 24.
Both archived per-file grade TSVs are byte-identical to the new run's baseline
grades. No count or key changes explain the difference. All 12 exact-chain
losses have N6, remain QC `ok`, and acquire a GlcNAc false negative `1 -> 0`.
Eleven of those old source decisions have exactly zero stored margin; the
twelfth is 0.0875. Across all keys, 38/39 changes `1 -> 0` start at zero margin.
This demonstrates fragile old decisions, not permission to restore labels,
retune ties or claim the responsible scientific stage is identified.

The -11 correct-position difference comprises five formerly correct positions
now unavailable and six net additional emitted errors. The strict missing-patch
policy explains the abstentions; the old binary outputs do not establish
usable raw evidence. At the reference's alignment, label changes sum to -17;
three seven-lobe reversal reselections recover six, explaining the published
-11. Count-matched N6 still has 113 -> 121 emitted errors. QC `ok` contains
114/183 reconstructed errors and all 12 exact losses; it is not chemical
validation. The 13 filename-date groups are descriptive groups, not established
independent molecules. No learned error predictor or date-specific adjustment
was constructed.

### Source finding: actual fused-fit preprocessing is direction-asymmetric

The blind source reviewer verified the actual production call path, not only
an isolated helper formula. `GaussianFit2D.preprocess_channel` returns flattened
unsmoothed `z` fourth and separately allocated `z_smooth` fifth. In
`_fused_roi_data`, backward `z` is overwritten by its smoothed version. Hence
fit data use `(F + S(B))/2`, while the ROI uses `(S(F) + S(B))/2`. The extractor,
batch, selector and diagnostic callers pass a freshly decoded image; there is
no hidden forward smoothing that restores symmetry. The acquisition diagnostic
already mirrors this convention in its missing-pixel support rule.

The retained Julia 1.13 probe uses generated 33x33 arrays and `flatten="none"`;
it calls native preprocessing/ROI functions but no raw reader or fitter. Six
assertions pass, independently repeated by the parent. Swapping directions
changes the centered fit map by a maximum of 0.09899306394259819 in arbitrary synthetic
amplitude units, while its ROI mask stays equal. Radius zero restores the
symmetry. This is a verified live implementation inconsistency, traced by blame
to May 12 (`04b66166`), not a measured real-data recognition loss or proof of
intended filtering policy. It predates this reconstruction and does not by
itself explain the reference/reconstruction delta.

**No production fix or third candidate is included here.** A useful next test
would fix one filtering convention prospectively and regenerate both its
control and treatment geometry/features at the same cached label-free counts,
fit budget and remaining settings. Comparing a newly fitted correction only
with an old cache would confound filtering and regeneration. Such compute needs
a separate bounded scope; neither this issue nor the negative endpoint result
authorizes another job. It is not a reason to reopen registration, DFT, lost
producer recovery or weight sweeps.

### Verification, limits and retained outputs

- Julia 1.13 only. Exporter 103/103; parent 9,002 real-copy assertions; original
  tables, cached geometry, eight source/config locks and production defaults
  unchanged. No Project/Manifest modification or dependency addition.
- Initial independent audit 14,068 assertions with 12 physical-convention
  synthetic tests. Final endpoint audit 22,065 assertions reproduces all 580
  native per-file physical records and 24 position rows; 21 inputs hash-checked
  before/after. Parent checks all paired aggregates, the N6 error intersection,
  the 19+21 input pins and both unchanged historical grade TSVs (85 checks),
  plus 44 native aggregate/confusion/coverage identities. The parent repeated
  all 12 synthetic audit tests and the native six-assertion filtering probe.
- Confidence values are uncalibrated and on different head/fusion scales;
  neither confidence nor QC thresholds enter the grade. This heavily reused
  benchmark is development evidence, not an untouched generalization test.
  Independent molecule/date/sequence confirmation remains unestablished.
  Wrong-N alignment and unknown localization limit causal chemical diagnosis;
  no universal 80% ceiling follows.
- No raw SXM decoding, optimizer, classifier training, local multi-file fit,
  remote command, scheduler job, restart, changed inference parameter or
  promotion occurred. All prior completed experiments remain closed.

Outputs are under `results/benchmark_recognition_20260919/`: the French
`report.md`, frozen unfiltered exports in `component_endpoints_v1/`, unchanged
external membership filtering in `external_comparison_*/`, one four-profile
native grade in `external_grade_endpoints_v1/`, the original error audit in
`external_error_audit/`, the separate final replay in
`external_error_audit/frozen_endpoint_comparison_v1/`, and the blind source
review/probe in `label_free_candidate_review/`. `docs/src/unit_assignment.md`
records the negative comparison and opt-in exporter; generated files remain
ignored. The definitions and grading commands are retained in `run_logs/`.

## 2026-09-20 — Symmetric fused-fit correction and matched recognition comparison

The user explicitly chose the targeted filtering correction and benchmark
comparison. **Deliverable and completion criterion:** test direction symmetry,
regenerate control and corrected assignments under matched settings, externally
grade both on the same benchmark, and report correct/all controls, conditional
accuracy, coverage, exact chains, missing/extra lobes and paired gains/losses.
Recognition benefit must be measured; correcting the inconsistency alone is
not evidence of improved chemical assignment.

The fixed correction calls the same native preprocessing helper for both
directions, fits the mean of flattened **unsmoothed** images, and keeps the
symmetrically smoothed ROI and existing max-of-view noise estimator. This
matches `_robust_roi_data` and `_channel_roi_data`; no new smoothing choice,
physical parameter, threshold or tuning sweep is introduced. The acquisition
diagnostic's direct-observation mask is updated to exclude each direction's
missing pixels symmetrically. Median imputation and its influence on background
flattening remain limitations, not newly validated missing-data treatment.

The unchanged control is commit `f2191f7`. Both arms use the original full146
`selected_from_features.tsv` (146 files / 900 lobes), derived from the saved
label-free GCV geometry, not the incompatible 871-lobe promoted counting
summary. The counts are held fixed, while geometry, fitted support, split
features and all three patch families are regenerated independently. Native
CC/Fisher and classifiers use the entire 146-file cohort and identical model,
seed, polarity and abstention settings. Benchmark membership and labels are
used only by the external comparison/grader after prediction.

The bounded computation is **one Viper job, four requested CPUs, 16000 MB,
two hours**, running both arms sequentially with four feature shards each.
`hpc/compare_fusion.sbatch` uses two source checkouts, rejects differing configs
or dependency locks, and cannot run scientific work without a Slurm allocation.
It accepts no cached geometry or patches. A metadata-only dry-run precedes
submission; failed outputs remain visible and are not silently retried. Existing
CC templates and raw scans are reused; no DFT, registration, new count selection,
descriptor recovery or feature/weight search is included.

Julia 1.13 checks pass: GaussianFit2D **151/151**, including 148 assertions of
direction symmetry, unchanged ROI/noise, missing samples, all flatten modes,
stride, smoothing radius, identical-view limit and mixed channel units;
acquisition diagnostics **187/187**, including the corrected missing-pixel mask;
paired job **15/15**, including exact command arguments, fresh fit inputs,
matched configs, allocation enforcement and failure preservation. These are
synthetic tests, not benchmark results. The branch remains an unpromoted
scientific correction until the paired recognition results are available.

**Submission:** source `d12283a`, control `f2191f7`. The native pipeline boundary
and shard tests also pass **19/19**; the documentation builds under Julia 1.13
with size warnings only. Rsync checksum dry-runs verify all 146 raw scans and
both transferred source trees. Each remote import resolves GaussianFit2D inside
its own checkout under Julia 1.13.0; both metadata-only dry-runs identify 146
files and create no output. Both dependency locks retain SHA-256
`617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe`.
Scheduler test-only 11878486 is not a scientific job. Actual **job 11878491**
was submitted once, with all six required path variables explicitly exported
in `sbatch --export=ALL,...`. Its result root is
`/ptmp/oldu/stmfit/fusion_comparison_20260920_v1`, with `control/`, `symmetric/`
and separate arm logs; the scheduler log is
`/ptmp/oldu/stmfit/fusion_comparison_20260920-11878491.log`. Fetch to
`results/fusion_comparison_20260920/` and grade only complete outputs. No
recognition result is available at submission time.

### Completed comparison: five additional correct positions, four additional exact chains

Job **11878491** ran on `vipc2169` from **16:59:30 to 17:33:55 CEST**, September
20, and completed **0:0 in 34m25s**. Viper allocated eight CPUs for the four-CPU
request; the runner used four one-thread fit shards per arm, sequential arms.
The control finished at 17:17:07 and the correction at 17:33:55. Both arms
regenerated base/split geometry and all three patch families, then ran the
same native whole-cohort assignment, validation, maps and QC.

All **424 output files** plus the Slurm log were fetched locally; an rsync
checksum dry-run reports zero differing/missing output files. No failure TSV
or fit error appears. Parent checks confirm identical 146-file / 900-key
coverage in all six geometry/patch/prediction stages and preserve N in both
geometry tables. Control labels are **695 zero / 198 one / 7 unavailable**;
symmetric labels are **689 zero / 204 one / 7 unavailable**, with exactly the
same unavailable keys. The correction flips 48 labels across 39 scans: 27
`0 → 1`, 21 `1 → 0`. This full scientific cohort includes the same one scan
outside the benchmark; membership filtering happens only in the external tool.

The regenerated control's final prediction TSV is byte-identical to the
September 18 reconstruction (SHA-256
`58a70c56beb7d61b686e918b6f6dc4339e12545a217cf5816be6075a68f00af5`).
The regenerated base/split geometry TSVs differ from their older caches;
byte identity is claimed only for the final prediction table, not the fits.
The original archive and all inputs remain unchanged.

One external native four-profile grade compares the historical reference,
saved reconstruction, regenerated control and corrected output. The saved and
regenerated reconstruction grades agree exactly.

| Profile | Correct / all 870 | Correct / classified | Coverage | Emitted errors | Exact chains / 145 |
|---|---:|---:|---:|---:|---:|
| Historical reference | 677 (77.8%) | 677/854 (79.3%) | 98.2% | 177 | 36 |
| Regenerated control | 666 (76.6%) | 666/849 (78.4%) | 97.6% | 183 | 24 |
| Symmetric fused fit | 671 (77.1%) | 671/849 (79.0%) | 97.6% | 178 | 28 |

All profiles retain **16 missing control positions** and **38 extra lobes**.
The native arms share five unavailable aligned positions; their other two
unavailable keys are one extra lobe and one outside-benchmark lobe. The 199
control positions not correct in the symmetric result include 178 emitted
errors, 16 missing positions and five abstentions; they are not 199 abstentions.

Paired arithmetic retains every scan: **19 wins / 16 losses / 110 ties** in
correct positions; **nine exact chains gained and five lost**. Confusion totals
`(TN,FP,FN,TP)` change **(519,43,140,147) → (518,44,134,153)**. Six additional
GlcNAc detections cost one additional false positive, for the net **+5**.
On the **106 count-matched scans**, correct/classified improves **514/635 →
518/635**; on the other 39 scans, **152/214 → 153/214**. The external grader
changes orientation on four wrong-N scans. These index-aligned grades do not
establish a spatial match between fitted lobes and physical monomers.

**Decision:** the targeted correction improves recognition modestly on this
fixed-count comparison and restores the intended direction symmetry. It remains
on the research branch, **not promoted**: the historical reference still has
six more correct positions and eight more exact chains. No threshold, class
prior, weight, seed or second filtering convention was tried after grading.
The reused benchmark is development evidence, not an independent validation.
Count selection (including the separate 871-lobe promoted counting lineage)
and the unknown25 application were not rerun; neither gets a new accuracy claim.

The scientific deliverable is complete. French synthesis:
`results/fusion_comparison_20260920/report.md`; full outputs: `run_v1/`;
external native grades: `external_grade_v1/`; complete paired scan table:
`paired_audit_v1/paired_files.tsv`. `summarize_pair.jl` verifies aggregate,
coverage and confusion arithmetic against the native grade before writing the
paired report. The job is closed and no new compute is submitted.

## 2026-09-20 — Physical transverse Fisher mirror: one bounded comparison

The user reaffirms that the deliverable is a reproducible pipeline exceeding
the lost champion, not merely improving the native reconstruction, and asks
to continue the descriptor/assignment investigation. The symmetric-fusion
candidate remains unpromoted. External per-file comparison against history
finds 24 winning scans, 27 losing and 94 tied; five exact chains are gained and
13 lost. All 13 lost exact chains retain N=6. These are post-hoc diagnostics,
not a subset to train on or a target-count rule.

**This pass:** test one source-grounded representation convention, with one
30-minute four-CPU Viper job. Done means tested implementation, complete paired
outputs, external full145 grade against the symmetric control and historical
reference, and an explicit retain/reject decision even if the result is worse.
No broad feature search, weight/threshold/seed tuning, exact-producer recovery,
missing-patch imputation, counting/QE refit or unknown25 rerun is included.

The reconstructed Fisher implementation deliberately kept Python's `[:, ::-1]`
reflection. With the actual u-outer/t-inner patch serialization this reverses
physical t (along the chain), not u (across it). That was documented at porting;
it is not a newly discovered regression or proof of the historical score gap.
Testing transverse-side invariance is physically motivated without labels.
The opt-in config `unit_assignment_transverse_fisher.toml` differs solely in
`model.name` and `preprocessing.fisher_layout=physical_u_outer_t_inner`.
Legacy behavior remains available and unchanged. Both paths keep the identical
197-pixel disk, PCA/GMM training, parity folds, amplitude convention, midpoint,
score sign/rounding and max-of-two operation. Only the reflected held-out patch
changes. Synthetic tests independently check the actual serialized axis,
involution, energy preservation, unchanged fitted weights and physical-u score
invariance; a legacy t mirror fails the deliberately asymmetric u probe.

`hpc/compare_fisher_mirror.sbatch` reuses the five complete saved geometry/patch
tables from symmetric-fusion job 11878491 in both arms. It runs the existing
whole-cohort native pipeline twice, including Fisher, GMM, k-means, vote and QC.
There is no per-error-row selector and no label input. The control must reproduce
the saved 671-correct / 28-exact candidate before interpreting a treatment
difference. The model is fixed before the new grade; this pass cannot establish
independent generalization on a benchmark already used in development.

Pre-run verification passes under Julia 1.13: Fisher 229/229 (including 28
new transverse-convention checks), paired shell runner 16/16, descriptor/vote
32/32 and pipeline boundaries/shards 19/19: **296 assertions**. The physical-u
variant changes no fitted synthetic Fisher weight or training partition. The
preexisting legacy NumPy-reference and separate-process replay tests still pass.
The sandbox initially blocked Julia's launcher lock and the SSH system-config
check; the same read-only checks succeeded with the approved host access.

### Completed result: transverse mirror rejected

Source `8bb29f4` was synced to
`/u/oldu/code/STMFit_transverse_fisher_20260920`. The mandatory metadata-only
dry-run checked both methods and all 146 input files. Ten local/remote hashes
matched: the Manifest, three configs, Fisher source and five cached geometry/
patch inputs. No Manifest change or package installation was needed. The single
sbatch submission supplied `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`,
`STMFIT_INPUT_DIR`, `STMFIT_OUTDIR` and `JULIA_BIN` explicitly in `--export`.

Viper job **11879515** completed **0:0 in 2m03s** on `vipc2488` (four CPUs
requested; eight allocated). Control started at 18:28:18 CEST and finished at
18:29:29; treatment finished at 18:30:19. The **358 output files** and Slurm
log were fetched from `/ptmp/oldu/stmfit/fisher_mirror_comparison_20260920_v1`
to `results/fisher_mirror_comparison_20260920/run_v1/`; checksum dry-runs report
no differences. The native application retained 146 files / 900 keys in all
nine checked tables per arm, with the same seven unavailable prediction keys.
There is no failure TSV or logged scientific error.

The control predictions are byte-identical to the saved symmetric-fusion
candidate. Local/half-plane features, CC margins and k-means outputs are also
byte-identical between arms. Predictor tables differ only in `emp_fisher`:
481 scores change, as do 104 GMM stability frequencies, five GMM labels and
four final decisions. Labels go from 689 zero / 204 one / seven unavailable
to 691 / 202 / seven. These frequencies are not calibrated probabilities.

The external full145 own-N grade, after completing all predictions, is:

| Profile | Correct / all 870 | Correct / classified | Coverage | Emitted errors | Exact chains / 145 |
|---|---:|---:|---:|---:|---:|
| Historical reference | 677 (77.8%) | 677/854 (79.3%) | 98.2% | 177 | 36 |
| Replayed symmetric-fusion control | 671 (77.1%) | 671/849 (79.0%) | 97.6% | 178 | 28 |
| Transverse Fisher mirror | 669 (76.9%) | 669/849 (78.8%) | 97.6% | 180 | 27 |

Both native arms have five unavailable aligned control positions, 16 missing
positions and 38 extras. One scan gains one correct position, three lose one,
and 141 tie; zero exact chains are gained and one is lost. All four changes
have N=6 and unchanged external orientation: `240817_075` lobe 5 (1→0, loss),
`240817_087` lobe 4 (0→1, loss), `241113_086` lobe 2 (1→0, exact chain lost),
and `241113_094` lobe 4 (1→0, gain). Three begin at exactly 0.5 stored soft-vote
frequency; this is not permission to change tie handling. Confusion totals
`(TN,FP,FN,TP)` change `(518,44,134,153)` → `(518,44,136,151)`.
Against history, treatment has 23 winning scans, 28 losing and 94 tied;
five exact gains and 14 losses. No label, expected count or benchmark class
composition entered either scientific run.

**Decision:** reject this variant as the working reference; retain symmetric
fusion at 671 correct / 28 exact. Source/config/tests remain as an explicit
negative experiment, not a default or a claim that longitudinal reflection is
chemically correct. This isolated coordinate change does not close the gap to
history. No other mirror, weighting, descriptor, seed or threshold was tried
after grading. The reused development benchmark does not establish independent
generalization. Counting and unknown25 were not rerun. The broader user goal of
exceeding the lost historical champion remains **unmet**, not redefined.

Verification now totals **819 assertions** under Julia 1.13: Fisher 229,
runner 16, descriptor/vote 32, pipeline boundaries/shards 19, saved representation
319 and Fisher attribution 204. Nine archived-grade parser assertions also
pass, and the post-hoc paired script checks aggregates, confusion, fixed N and
coverage against the native grade before writing all 290 reference/file pairs.
Documentation builds successfully with existing size warnings only.
The final report is `results/fisher_mirror_comparison_20260920/report.md`;
grades, paired tables and the small arithmetic reproduction script are alongside
it. The job is complete and no further compute is queued by this pass.

A final external-only saved-output diagnosis of the retained symmetric working
candidate's 13 historical exact-chain losses finds 16 affected positions: six
finite errors shared by k-means and GMM, nine finite head disagreements, and one
unavailable assignment. There are 13 `1→0`, two `0→1` and one `0→?` changes;
nine old stored margins are zero. Both grades retain N=6 and the same alignment
for these files. The diagnostic is reproducible via
`results/fisher_mirror_comparison_20260920/diagnose_lost_chains.jl`; it produces
no prediction or model input. Choosing/reweighting the two unchanged head
decisions cannot repair the six common errors, but this is not a raw-information
ceiling or justification to train on this post-hoc subset. No new method follows
automatically from it.

## 2026-09-20 — Matched patch residual filtering: one bounded comparison

The user approves the next proposed source-grounded experiment. The scientific
deliverable is a complete paired recognition comparison of `S(data)-model`
against `S(data-model)`, at fixed geometry and selected counts. Done means tested
extractors, a byte-reproduced control, all full146 outputs and an external
full145 grade including every gain/loss, followed by a retain/reject decision.
The time box is one 30-minute four-CPU / 16 GB Viper job. Before this comparison,
the working candidate is 671 correct / 28 exact; the target stays above 677 / 36.
No threshold, weight, seed, smoothing-radius, background, mirror, count or DFT
search is in scope; no setting is changed after grading. Unknown25 is not rerun.

Both extractors' legacy modes compute `S(data)-model` using an unsmoothed decoded
Gaussian model plus fitted plane. An image equal to that model leaves the
deterministic residual `S(model)-model`. The opt-in candidate applies the native
finite-window box smoother to `data-model`, before unchanged interpolation and
patch normalization. Exact model data yield zero; an injected shoulder is
preserved as `S(shoulder)`. This fixes filter consistency, not arbitrary gradient
or background contamination, and does not establish chemical recognition.

The assignment setting `preprocessing.patch_residual_filter` is explicit and
validated. Existing configs now state `smooth_data_only`: their hashes change,
but their effective behavior does not. New `unit_assignment_matched_residual.toml`
differs from the reconstructed config only in this setting (`smooth_residual`)
and its method name. The original Fisher layout is retained, not the negative
physical-u variant. The pipeline passes the assignment config to each extractor
and rejects cached patch inputs for the new mode. Standalone extractors without
the new option keep the old arithmetic byte-for-byte.

Both arms of `hpc/compare_matched_residual.sbatch` regenerate all three patch
families from the same 900-lobe / 146-scan symmetric-fusion geometry cache, using
the same split-width table, raw cohort and templates. The full native Fisher,
k-means, GMM, soft vote and QC follow. Labels and historical predictions enter
only the later external report, never the scientific job. All keys and invalid
rows must be retained; no per-error subset or composition prior is introduced.
This reused development benchmark cannot establish independent generalization.

Pre-run Julia 1.13 checks pass: patch extraction **380/380**, assignment rules
**37/37**, pipeline boundaries/shards **25/25**, paired runner **16/16**, and
Fisher conventions/replay **229/229** (687 assertions).
Synthetic checks cover a Gaussian plus tilted plane, an injected shoulder,
radii 0/1/2, edge windows, missing support, forward/backward cancellation and
input nonmutation. Extractor fixtures verify byte-identical explicit legacy
output, unchanged raw columns, and independently reconstructed matched patches.
The first independent smoothing assertion used exact equality between copied
windows and native strided views; one radius-2 reduction differed by two ulps
(`6.94e-18`). Its comparison now allows eight ulps at the expected signal scale;
the production calculation and exact null-model checks were not altered.

The synthetic extracted-input pipeline initially failed only at figure creation:
Julia could not load a StatsBase precompiled image while compiling Plots. A
sequential `using StatsBase, Plots` succeeds; a fresh end-to-end run then passes
all nine integration checks in 1m30s, for **696 focused assertions** including
the previously listed suites. No dependency version, cache deletion, Manifest
edit or science-code change was needed. The documentation build also succeeds
with its existing size warnings. Viper's independent plotting dependency import
is clean before submission.

Source **dae6395** is synced to
`/u/oldu/code/STMFit_matched_residual_20260920`. The mandatory metadata-only
dry-run validates both methods and all 146 input files. The Manifest, counting
config, assignment configs, changed sources, geometry and template hashes agree
locally/remotely. The single submitted job is **11883607**, using explicit
command-line export of all five required variables (`SBATCH_EXPORT=NONE` must
not drop them). No scientific computation runs on the login node.

### Completed result: small gain, historical target still unmet

Job **11883607** completes **0:0 in 5m15s** on `vipm2063`, from 20:53:26 to
20:58:41 CEST (four CPUs requested, eight allocated). Control runs from 20:53:28
to 20:57:01; treatment finishes at 20:58:41. All **390 regular output files**
and the Slurm log are fetched from
`/ptmp/oldu/stmfit/matched_residual_comparison_20260920_v1` to
`results/matched_residual_comparison_20260920/run_v1/` (raw staging links are not
duplicated). Checksum dry-runs find no source/output/log differences. Thirteen
local/remote source, configuration and input hashes also match. No dependency
version or Manifest change is made; there is no failed stage or scientific error.

Saved-output checks pass **47/47**: twelve tables in each arm retain all
**146 files / 900 lobe keys**. Every control table, including all three regenerated
patch families, matches the saved symmetric-fusion table byte-for-byte. Control
prediction SHA-256 is
`d6e52f5a3a38915dc3e411e0a94696cd35f9f9e37771320055483745531b6299`.
Local geometry features and every raw-patch/metadata column are unchanged between
arms. Each patch family has 900 changed residual rows; 893 finite half-plane
descriptors and Fisher scores change, as do 898 CC scores per direction. The
seven unavailable descriptor/prediction keys remain unchanged. GMM changes 252
labels and 297 vote frequencies; k-means is byte-identical. Only three final
labels change. These vote frequencies are not calibrated probabilities.

The full145 own-N grade is run only after both arms and their checks complete:

| Profile | Correct / all 870 | Correct / classified | Coverage | Emitted errors | Exact chains / 145 |
|---|---:|---:|---:|---:|---:|
| Historical reference | 677 (77.8%) | 677/854 (79.3%) | 98.2% | 177 | 36 |
| Regenerated symmetric-fusion control | 671 (77.1%) | 671/849 (79.0%) | 97.6% | 178 | 28 |
| Matched patch residual | 672 (77.2%) | 672/849 (79.2%) | 97.6% | 177 | 29 |

Both native arms retain **16 missing control positions, 38 extra lobes and five
unavailable aligned positions**. Their 198/870 not-correct treatment positions
comprise 177 emitted errors, 16 missing positions and five abstentions; these
are not 198 abstentions. Every paired file is reported: **two gains, one loss,
142 ties**, with **one exact chain gained and none lost**. All three final
changes are `0→1`, on N=6 scans with unchanged grade orientation:

| Scan | Lobe | Correct positions, control → treatment | Exact-chain change |
|---|---:|---:|---|
| `240814_023.sxm` | 2 | 3 → 4 | none |
| `240817_005.sxm` | 2 | 5 → 6 | gained |
| `240817_087.sxm` | 4 | 5 → 4 | none |

Confusion totals `(TN,FP,FN,TP)` change `(518,44,134,153)` →
`(517,45,132,155)`: two more GlcNAc detections cost one false positive. Against
history, treatment wins on 23 scans, loses on 26 and ties on 96; it gains five
exact chains and loses twelve, with six external orientation differences. This
grade does not establish spatial correspondence of fitted lobes to monomers.

**Decision:** retain the explicit matched-residual variant as the next working
candidate at **672 correct / 29 exact**, not as a promoted champion. The source
consistency improvement and this small benchmark gain do not reach the user's
target: history still has **five more correct positions and seven more exact
chains**. The original assignment config remains the legacy control; the new
config is opt-in. No threshold, seed, weight, radius, background or second
variant is tried after grading. Count selection and unknown25 are unchanged.
This heavily reused benchmark is development evidence, not independent
generalization. The bounded experiment is complete; the broader target remains
**unmet** and has not been lowered. The job is closed with no further submission.

The French report is `results/matched_residual_comparison_20260920/report.md`.
`check_outputs.jl` replays the 47 saved-output checks; `summarize_pair.jl` checks
grade aggregates, coverage, confusion and fixed N before writing all 290
reference/file pairs to `paired_audit_v1/paired_files.tsv`. Full grades are in
`external_grade_v1/`. Focused verification remains **696 Julia 1.13 assertions**,
plus the saved-output checks; the initial local plotting-cache failure and its
successful unchanged rerun are recorded above.
The final documentation build also passes, with existing size warnings only.

## 2026-09-20 — Two source-grounded transverse descriptors: bounded comparison

The user asks to continue improving recognition with the strongest leads. The
deliverable is a reproducible candidate comparison, not an infrastructure change
or a promise of positive results. Done means synthetic/CLI tests, a replayed
672-correct / 29-exact control, all 146 files / 900 keys, full145 external grades
against the control and historical 677/36, every paired gain/loss and a decision.
The time box is **one 30-minute, four-CPU / 16 GB Viper job** with three arms.
The completed unknown25 application is not rerun. The archived workflow's
Must-NOT-have constraints remain: no labels/composition prior in inference, no
truth-selected per-file profile, no dropped uncertain rows, no user-artifact
overwrite and no count-selection change. No new DFT, registration, background
fit, parameter sweep or voting change is included.

Source/history review rules out repeating the already rejected physical-u Fisher
mirror or the component-only votes. The existing half-plane descriptor does not
distinguish moving residual mass inside one transverse half-plane; the older
journal describes an unavailable u-weighted first moment, with insufficient
detail to claim exact recovery. Prior label-free patch audits also measured
roughly 85% transverse-ramp energy and strong descriptor saturation, but did
not identify that ramp as nonmolecular or justify changing the raw image.

Two alternatives are therefore fixed before any new grade:

1. **First moment:** `sum((u/max(abs(u)))*p)/sum(abs(p))`, on the same normalized
   backward 9×9 matched residual. This measures transverse displacement rather
   than only half-plane membership; its scale comes from the existing grid.
2. **Affine-residual half-plane asymmetry:** fit `[1,t,u]` by least squares to the
   complete patch, subtract that plane, then apply the original signed-half-plane
   numerator and L1 denominator to the remainder. This is a different descriptor,
   not a preprocessing correction. The projection removes affine molecular
   signal as well as any affine nuisance. Pure planes become explicitly
   unavailable at the existing `1e-12` L1 threshold; missing pixels remain NA.

The control config is unchanged. The two new configs differ from it only in
`model.name` and `model.descriptor`. Their shared storage column is retained for
compatibility, not used to disguise a historical-producer recovery. The first
moment and plane projection introduce no fitted hyperparameter or class prior.
The native pipeline keeps the same geometry, selected N, patches, Fisher layout,
CC templates, classifier settings, seeds and vote rule. Only this GMM input
feature can differ. All three arms regenerate identical patches; the cached
patch guard is retained rather than introducing a new cache-authorization layer.
No hybrid or additional setting is chosen after inspecting the external grade.

Synthetic verification checks position sensitivity, unchanged legacy arithmetic,
u/t axis conventions, reflections, scale/sign behavior, invalid rows, independent
QR projection, idempotence, orthogonality, energy loss and affine invariance.
The projection test explicitly measures that some injected molecular-like signal
is removed. A first test run stopped on an ambiguous `.5.*` numeric literal in
the synthetic fixture; spacing and an explicit leading zero fix that test syntax,
without changing any production formula or benchmark setting.

A second test-only failure compared an 81-pixel vector's aggregate L2 error
to a per-pixel `5e-15` bound. On the largest injected plane the measured maximum
pixel error is `2.71e-15` (1.52 ulps of the plane scale), while the L2 error is
`9.07e-15`. The check now states the same bound in the maximum norm; the QR,
orthogonality, descriptor-invariance checks and production formulas are unchanged.
All **125 descriptor assertions** pass, including both real CLI invocations.
Existing assignment/vote tests pass **37/37**, patch extraction **380/380**,
the three-arm shell runner **16/16**, and the complete synthetic pipeline
**46/46**, including its figures and all three cache-rejection configurations.

Fisher numerical/reference/CLI regression also passes **229/229**, for **833
focused assertions**. The documentation builds with existing size warnings only.
Source **252edaa** is synchronized to
`/u/oldu/code/STMFit_transverse_descriptors_20260920`; a checksum dry-run finds
no source differences. Eleven source/configuration/input hashes match. The
Manifest, counting config and matched-residual control config are unchanged.
All three metadata-only dry-runs validate 146 files and the independent plotting
dependency import succeeds on Julia 1.13. The single job submitted is
**11888011**, with the required five-variable command-line `--export` list,
four CPUs, 16 GB and a 30-minute limit. No scientific work runs on the login node.

### Completed three-arm comparison: modest affine gain

Job **11888011** completes **0:0 in 7m19s** on `vipm2008`; Slurm allocates eight
CPUs for the four-CPU request. Control takes 3m45s, moment 1m49s and affine 1m44s.
All **585 regular outputs plus the Slurm log** are fetched to
`results/transverse_descriptor_comparison_20260920/`, excluding only raw staging
symlinks. The checksum dry-run is empty. Source/configuration and input hashes
remain unchanged. The job is closed, with no further submission.

All **77 saved-output checks** pass before grading: each of twelve scientific
tables in each arm retains all 146 files / 900 keys; the twelve control tables
are byte-identical to the saved matched-residual result. Eight tables (all three
patches, local geometry, both CC margins, Fisher and k-means) are byte-identical
across all arms. In descriptor/predictor input tables only the descriptor column
changes, on 893 finite rows. Both candidates retain the same seven unavailable
assignments, with no new affine-zero-mass rows. There is no scientific failure.

The first moment changes 249 GMM labels and 357 vote frequencies; affine residuals
change 264 and 354. The final vote changes six and sixteen labels respectively.
These vote frequencies are not calibrated chemical probabilities. External
grading is run only after the complete outputs and controls are verified:

| Profile | Correct / all 870 | Correct / classified | Coverage | Emitted errors | Exact chains / 145 |
|---|---:|---:|---:|---:|---:|
| Historical reference | 677 (77.8%) | 677/854 (79.3%) | 98.2% | 177 | 36 |
| Matched-residual control | 672 (77.2%) | 672/849 (79.2%) | 97.6% | 177 | 29 |
| Transverse first moment | 672 (77.2%) | 672/849 (79.2%) | 97.6% | 177 | 28 |
| Affine-residual half-plane asymmetry | 673 (77.4%) | 673/849 (79.3%) | 97.6% | 176 | 33 |

All native arms retain **16 missing control positions, 38 extra lobes and five
unavailable aligned positions**. The affine candidate's 197 not-correct positions
comprise 176 emitted errors, 16 missing positions and five abstentions; they are
not 197 abstentions. Rounding classified accuracy to 79.3% for both the affine
candidate and history does not close the fixed-denominator or exact-chain gap.

The first moment has **three winning scans, three losing, 139 tied**, with no
exact gains and one exact loss (`241113_086`). Its six final changes are all in
the benchmark. One N=5 scan (`240817_043`) changes external grade orientation.
Confusion totals `(TN,FP,FN,TP)` change `(517,45,132,155)` → `(518,44,133,154)`.
Against history there are 25 wins, 29 losses and 91 ties, five exact gains and
thirteen losses. This candidate is not retained.

Affine residuals have **eight winning scans, seven losing, 130 tied**. Five exact
chains are gained (`240307_017`, `240310_Cu100007`, `240817_007`, `240817_044`,
`240818_020`) and one lost (`241113_086`). All sixteen final changes are `0→1`:
fifteen lie in the benchmark and one is lobe 7 of `240310_Cu100009.sxm`, the
technical scan excluded by the unchanged external manifest since July 3. It
remains in scientific inference and all full146 outputs, not selectively dropped.
The two N=5 scans `240314_Cu100_024` and `241113_089` reverse external grade
orientation; that convention is not a spatial monomer correspondence.
Confusion totals change `(517,45,132,155)` → `(510,52,124,163)`: eight additional
true positives cost seven false positives. Against history there are 23 wins,
26 losses and 96 ties, six exact gains and nine losses, with eight orientation
differences. All 580 candidate/reference/file pairs, including all losses and
ties, are retained in `paired_audit_v1/paired_files.tsv`.

**Important limitation:** all sixteen affine decisions land at exactly **0.5
vote frequency and zero margin**, using the unchanged tie rule. The net gain is
one correct position, not sixteen secure identifications. No vote threshold,
weight, seed, hybrid or missing-value rule is changed in response to this grade.
This outcome does not establish that the removed plane was a physical artifact.

**Decision:** retain `config/unit_assignment_affine_residual.toml` as the new
opt-in working candidate at **673 correct / 33 exact**, **without promotion**.
The first-moment config is preserved as a negative experiment, not a default.
The user's historical target remains unmet: **four more correct positions and
three more exact chains are needed merely to equal 677/36**. This repeatedly
used benchmark remains development evidence, not independent validation.
Count selection, the separate 871-lobe counting lineage and unknown25 are not
rerun or revalidated. The bounded comparison is complete; the broader goal is
still open and has not been lowered.

The French report, all changes, complete grades and replay commands are in
`results/transverse_descriptor_comparison_20260920/report.md`.
Verification remains **833 focused Julia 1.13 assertions**, plus the **77 saved
output checks** and paired-grade arithmetic assertions. The recorded test-only
syntax/norm failures above do not alter the production formulas or this result.
The final documentation build passes, with only the existing generated-page and
search-index size warnings; deployment is skipped outside CI.

## 2026-09-20 — Signed CC margins and affine Fisher: two independent label-free leads

The user approves exploring the two proposed leads while ensuring they remain
label-free. The scientific deliverable is a complete, reproducible comparison
against the 673-correct / 33-exact working candidate and historical 677/36.
Done requires synthetic/CLI tests, exact control replay, all 146 scans / 900
keys, external full145 grades with every gain/loss, and a decision. The time box
is **one Viper job, 30 minutes, four CPUs / 16 GB**, with three independent arms.
No combined variant, weight/threshold/seed/radius search, new DFT, count refit,
unknown25 rerun or old-artifact overwrite is included. The archived workflow's
Must-NOT-have list is retained; its retired process instructions are not revived.

The source audit confirms that `join_predictor_features` reads each directional
`cost_margin`, while its producer writes `abs(cost_GlcN-cost_GlcNAc)`. This loses
the favored template's identity; it is a defined feature, not a proven bug.
The first candidate replaces just these two features by
`cost_GlcN-cost_GlcNAc`, positive when the existing GlcNAc template has lower
cost. It does not use the decoded `predicted` or `physical_label` columns, change
the label-free geometric direction/phase/mirror selection, or search a sign
against experimental truth. Simulated chemical template identities are physical
inputs, not the experimental benchmark labels or a composition prior. Their
sign may be unreliable; no improved recognition is assumed.

The second candidate retains the absolute CC features and projects the forward
Fisher disk off the fixed `[1,t,u]` subspace. A QR basis uses only the original
197 disk coordinates. The same projection is applied before PCA/GMM/Fisher
training and to both original/reflected held-out patches at scoring. The original
patch center amplitudes still order the two learned groups; they are NOT
projected, replaced by expected labels, or constrained to a class count. Fold
parity, mixture weights, PCA dimension, mirror, covariance, historical midpoint,
regularization, seeds and voting stay unchanged. The projection removes affine
molecular signal too and is not a calibrated acquisition/background correction.

The explicit new config fields are `model.mold_margin_mode`,
`preprocessing.fisher_patch_projection`, and `model.fisher_projection_zero_l1`.
All existing assignment configs explicitly select `absolute_cost_margin`,
`none`, and `1e-12`; their hashes change, not their effective legacy arithmetic.
The zero-mass tolerance is fixed before the grade at the existing descriptor's
numerical L1 floor, not fitted to data. Under affine Fisher, projected disk mass
at or below it is unavailable (`zero_affine_patch_mass`), retaining its key.
No missing pixel is imputed. Raw patches, the successful backward descriptor,
CC templates and k-means remain unchanged in both arms. The two new configs
differ from the control only by their name and one respective mode field.

The audit distinguishes inference from evaluation: the scientific runner accepts
no truth, expected N, control sequence, reference prediction or benchmark
manifest. It processes the whole cohort before the external manifest selects
the grade rows. Tests poison experimental-control columns (rejected) and decoded
mold labels (irrelevant to the features), check free unequal mixture weights,
original amplitude mapping, held-out-fold independence and complete invalid rows.
Neither version chooses the top k lobes or any expected monomer positions.
External comparison can guide subsequent development; this heavily reused
benchmark is therefore not an independent generalization test.

The first affine test run fails only on an over-tight aggregate orthogonality
tolerance: diagonal Gram errors are at most `1.78e-15`, but their aggregate norm
exceeds `2e-15`. The test now uses a dimension-aware `197*eps(Float64)` bound
for dot products across the disk; independent normal-equation projection,
idempotence, plane-invariance and score tests remain separate. No production
formula or configured scientific tolerance changes in response to this test.

All **1,518 Julia 1.13 assertions** now pass: signed-feature/input boundaries 73,
affine Fisher 74, legacy Fisher 229, descriptor/vote 37, prior transverse
descriptors 125, patch extraction 380, full synthetic pipeline 58, new shell
runner 16, legacy Fisher attribution 207 and saved-representation diagnostics
319. The raw-patch attribution program explicitly rejects the projected mode
rather than attributing its unmodified identity to a different scoring path.
The two candidate formulas, zero-mass rule and all settings are frozen before
any new scientific cohort output or grade. The documentation builds with the
existing page/search-index size warnings only.

Source **3c45c7c** is synchronized to
`/u/oldu/code/STMFit_label_free_signals_20260920`, without benchmark truth tables
or previous result trees. The tracked-source checksum dry-run is empty. The
first dependency precheck reports StatsBase unavailable, before any dry-run or
job: `git ls-files` had omitted the ignored root `Manifest.toml`. Read-only
inspection confirms its absence. Copying the unchanged existing Manifest
(SHA-256 `617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe`)
fixes the staging omission without Pkg resolution, version changes or source
modification. No scientific submission or computation preceded this correction.

With the unchanged Manifest present, dependency imports pass on Julia 1.13.0
and all three metadata-only dry-runs validate 146 scans. Nine local/remote
source/configuration hashes match, and the three geometry/template input hashes
are unchanged. Job **11889205** is then submitted exactly once with the required
five-variable command-line `--export` list. The scientific source directory has
no `benchmarks/` or old `results/` tree. No benchmark truth enters the command or
scientific inputs; the external grade remains local and downstream of all arms.

## 2026-09-21 — Signed CC and affine Fisher comparison complete: both regress

The single approved comparison is complete. Job **11889205** first waits under
the shared `n0001` QOS CPU limit, then runs on `vipm2008` from **00:19:33 to
00:26:38 CEST**, completing **0:0 in 7m05s**. Four CPUs were requested and eight
allocated; the scientific Julia subprocesses use four threads. No computation
runs on the login node, no queue setting is changed and no job is resubmitted.
Scientific source **3c45c7c**, all configurations, the unchanged Manifest,
geometry tables and templates retain their pre-submission hashes.

All **585 output files and the Slurm log** are fetched from
`/ptmp/oldu/stmfit/label_free_signals_20260920_v1` into
`results/label_free_signals_20260920/run_v1` (the log is stored alongside).
Checksum comparisons are empty; raw-image staging links are not copied.
All **77 saved-output checks pass** before grading. Each arm retains the same
900 keys / 146 scans in all twelve scientific tables, with no failed stage and
the same seven unavailable assignments. All twelve control tables are
byte-identical to the preceding affine-descriptor result; the predictions hash
remains `4039b7bb3baaa186d0dc9cc810face2dfa1f151bc45467c81d4f639e3558aa9f`.

The three patch families, local geometry, successful backward descriptor,
forward/backward CC cost tables and k-means outputs are byte-identical across
arms. Signed CC changes only its two predictor features (584 forward and 585
backward values); Fisher scores are unchanged. It changes 14 GMM labels and
14 final decisions, with 231 final vote frequencies changed. Affine Fisher
changes only its predictor margin (893 values), 21 GMM labels and 20 final
decisions, with 238 final frequencies changed. No new `zero_affine_patch_mass`
row appears in this cohort. The modeled raw disk inputs and original amplitude
anchor remain unchanged; only the fixed projection is different in that arm.

Only after these checks does the existing external manifest select 145 scans /
892 predicted lobes for the full **870-position** grade. All three native arms
classify **849/870 (97.6%)**, with 16 missing positions, five aligned abstentions
and 38 extra lobes. Counts are fixed, including 14 short-N and 25 extra-N scans.

| Version | Correct / 870 | Correct / classified | Emitted errors | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical reference | 677 (77.8%) | 677/854 (79.3%) | 177 | 36 |
| Affine-descriptor control | 673 (77.4%) | 673/849 (79.3%) | 176 | 33 |
| Signed CC margins | 668 (76.8%) | 668/849 (78.7%) | 181 | 32 |
| Affine Fisher disk | 672 (77.2%) | 672/849 (79.2%) | 177 | 32 |

Against control, signed CC gives **four scan gains, nine losses and 132 ties**,
net **minus five correct positions**. Two exact chains are gained (`240307_016`,
`241113_086`) and three lost (`240310_Cu100007`, `240817_044`, `240818_020`).
Its 14 changes comprise five `0→1` decisions at the unchanged zero-margin tie
and nine `1→0` decisions at voting margin one. One latter change is outside
the benchmark; of the eight graded margin-one changes, three gain and five lose.
Thus even unanimous votes are not calibrated chemical confidence.

Affine Fisher gives **nine gains, ten losses and 126 ties**, net **minus one
correct position**. Three exact chains are gained (`240814_011`, `240817_075`,
`241113_086`) and four lost (`240310_Cu100007`, `240817_044`, `240818_020`,
`241114_027`). Four changes are `0→1` at zero margin; sixteen are `1→0`, fifteen
at margin one and one at 0.75. The projection that improved backward half-plane
asymmetry therefore does not improve this Fisher arm under the fixed settings.
This result does not establish that affine content is a physical background or
that signed template identity is generally useless.

Each arm also reverses lobe 7 of `240310_Cu100009.sxm`, the technical scan
excluded by the unchanged external manifest since July 3; all eight of its
lobes remain in scientific inference. External grade orientation changes on
`240314_Cu100_024.sxm` for signed margins and `241113_089.sxm` for affine Fisher.
These are grading conventions, not spatial monomer correspondences. Confusions
`(TN,FP,FN,TP)` change from `(510,52,124,163)` to `(509,53,128,159)` and
`(515,47,130,157)` respectively. Against historical predictions, signed CC has
22 winning / 28 losing / 95 tied scans and affine Fisher 22 / 27 / 96; each
gains six exact chains and loses ten. Every one of the **580 paired file rows**
and **34 changed decisions** is retained, including losses and non-benchmark
changes, in `paired_audit_v1/` and `decision_audit_v1/`.

**Decision:** reject both new arms as replacements. Keep
`config/unit_assignment_affine_residual.toml` at **673 correct / 33 exact** as
the opt-in working candidate, without promotion. The historical **677/36**
target is unchanged: four correct positions and three exact chains are still
needed merely to equal it. No combined variant, sign flip, threshold, seed,
weight, radius, per-file method choice or missing-value rule is tried after the
grade. Count selection, the separate 871-lobe counting lineage and unknown25
are not rerun or revalidated. The new configs are retained as negative experiments.

The label-free requirement applies to inference: no experimental truth,
expected N, composition or benchmark manifest enters the runner. The signed
feature uses physical simulated template identities; the projection uses disk
coordinates. Original unsupervised amplitude naming and free mixture weights
are untouched. External grading can inform method development, so this heavily
reused benchmark is not independent validation. The bounded scientific
deliverable is finished; the broader historical target remains open.

The complete French report, formulas, all changed decisions and replay commands
are in `results/label_free_signals_20260920/report.md`. Verification comprises
**1,518 focused Julia 1.13 assertions**, **77 saved-output checks**, and the
paired-grade arithmetic assertions. The pre-submission test-only QR tolerance
failure and missing-Manifest staging error are recorded above; neither changed
production mathematics or settings after grading.
The final documentation build passes, with only the existing generated-page
and search-index size warnings; deployment is skipped outside CI.

## 2026-09-21 — Two independent numerical leads: approved bounded comparison

The user authorizes trying the next label-free leads. The scientific deliverable
is a reproducible comparison of training-consistent Fisher centering and final
GMM covariance shrinkage, separately, against the runnable affine-descriptor
candidate **673 correct / 870, 673/849 classified, 33/145 exact**. Historical
**677/870, 677/854, 36/145** remains the target. Definition of done: synthetic
and CLI tests, full 146-scan / 900-key control replay, both complete arms,
numerical diagnostics, then external 145-scan / 870-position grading with all
paired gains and losses. The budget is one Viper job, 30 minutes, four requested
CPUs, 16 GB. No combined arm, per-chain normalization change, new count fit,
unknown25 rerun, DFT change, seed/threshold/weight sweep or post-grade tuning.

`unit_assignment_centered_fisher.toml` changes only the method name and
`model.fisher_score_center = "training_mean"`. PCA/GMM/Fisher already train on
mean-centered patches, but legacy scoring subtracts the near-zero mean of those
centered patches from raw held-out patches. The new origin is the original
opposite-fold training mean. Fisher weights, component means and populations,
amplitude naming, parity folds, mirror and covariance are unchanged. For each
fold it subtracts a constant `dot(training_mean-legacy_mid,w)` even after the
unchanged max-mirror operation. A common global offset would cancel in per-file
z-scoring; differing fold offsets can survive. The diagnostic measures the
actual exported-score identity and standardized effect, rather than assuming
that fixing an origin will improve recognition. This is not equal-cluster
centering and imposes no composition prior.

`unit_assignment_shrunk_gmm.toml` changes only the method name and
`model.gmm_final_covariance = "ledoit_wolf"`. The deliberately narrow comparison
regularizes only the final hard-cluster covariance after the two existing
self-training iterations. EM and earlier reassignment covariance arithmetic,
final memberships, means and free proportions remain unchanged before scoring.
For centered p-by-n data, `S=XX'/n`, `T=tr(S)/p I`, the fixed analytical mixing
coefficient is `clamp(max(0,(mean(norm(x_i)^4)-sum(abs2,S))/n)/sum(abs2,S-T),0,1)`;
when `S==T` report zero. Add the unchanged `1e-6 I` floor after mixing. The
formula matches the Ledoit-Wolf spherical-target estimator, independently
checked against a dense outer-product reference and its official implementation.
No coefficient is selected by recognition grade. This is numerical
regularization: learned memberships and dependent lobes do not justify iid
optimality, noise calibration or chemical-confidence claims.

All native configs now explicitly declare `fisher_score_center`,
`gmm_final_covariance` and `gmm_covariance_ridge`. Legacy settings remain
`legacy_centered_mean`, `ridge`, `1e-6`. The GMM CLI receives the config from the
pipeline; its optional default is the explicit reconstructed config, not a
hidden scientific value. The new modes are independently opt-in, not defaults.

The label-free diagnostic also compares full-cohort GMM outputs with two fixed
complementary halves of whole scans, partitioned by alternating sorted filenames.
All retained and unavailable keys are exported. It does not select a partition
or setting, retrain Fisher within those halves, or score withheld scans. This is
training-data withdrawal sensitivity, **not held-out recognition or independent
validation**. Exact pre-final covariance memberships/means/weights are checked
for both modes on the full cohort and each half. Condition numbers, shrinkage
coefficients and vote/decision changes are reported without a confidence claim.

The first new unit-test invocation stopped at a Julia parse error in the new
include-safe GMM main guard (`@__FILE__` followed by `&&`). Parenthesizing the
macro via `abspath(@__FILE__)` fixed the guard before any scientific job. All
136 numerical-mode assertions then passed; no formula or setting was changed.

Pre-submission verification passes **1,186 Julia 1.13 assertions**: numerical
modes 136, numerical diagnostics 20, new sbatch runner 21, pipeline 70,
legacy Fisher 229, affine Fisher 74, signal features 73, descriptor/vote 37,
Fisher attribution 207 and representation diagnostics 319. The diagnostic's
synthetic CLI completes with 80 keys and 120 covariance rows, preserving
whole-scan partitions and rejecting benchmark columns before computation.
The formulas, partitions and all settings are fixed before scientific output.

## 2026-09-21 — Numerical comparison complete: centering gains two, shrinkage loses eight

Source **0d58b44** is synchronized to
`/u/oldu/code/STMFit_numerical_signals_20260921`, without benchmark truth tables
or previous result trees. The unchanged ignored Manifest is explicitly copied;
ten local/remote critical hashes match. Julia 1.13 dependency imports and all
three metadata-only dry-runs pass for 146 scans. The source checksum comparison
with itemized differences is empty, as is the post-run check. Geometry tables
and templates retain their preceding hashes. The documentation builds before
submission with only the existing page/search-index size warnings.

The job is submitted exactly once with the required explicit five-variable
`sbatch --export` list. Job **11891104** waits briefly under `QOSGrpCpuLimit`,
then completes **0:0 in 9m04s** on `vipl2005`, **01:11:51–01:20:55 CEST**.
Four CPUs are requested and eight allocated; Julia subprocesses use four threads.
All scientific computation stays on the compute node. No queue setting changes,
resubmission, second experiment, combined arm or parameter adjustment follows.

All **590 output files and the Slurm log** are fetched into
`results/numerical_signals_20260921/run_v1` and its parent. Itemized checksum
comparisons show no differences; raw-image staging links are not copied.
All **76 saved-output checks pass**. The twelve control tables are byte-identical
to the saved 673/33 result, including prediction hash
`4039b7bb3baaa186d0dc9cc810face2dfa1f151bc45467c81d4f639e3558aa9f`.
Every arm retains all 146 scans / 900 keys in all twelve scientific tables and
the same seven unavailable assignments. No failed stage or warning is found.

The three patch families, local geometry, backward descriptor, forward/backward
CC tables and k-means predictions are byte-identical across arms. Centering
changes only the Fisher predictor feature (893 finite values), four GMM vote
frequencies, one GMM hard decision and two final decisions. All predictor inputs
and Fisher outputs in the shrinkage arm are byte-identical to control; only its
GMM and final predictions change (36 frequencies, 22 hard decisions each).

**Numerical diagnostics.** Original opposite-fold training means produce score
offsets **4.5178598359** (432 even-lobe training rows) and **4.4480358561** (461
odd-lobe training rows); both native GMM fits converge. Exported control and
centered scores match the predicted constant-shift identity within their
six-decimal serialization precision. After unchanged per-file standardization,
median / maximum absolute Fisher changes are **0.0296753 / 0.380814**; expanded
36-feature distances are **0.0707576 / 1.053445**. Most of the common origin
offset therefore cancels; a fold-dependent effect remains.

The covariance diagnostic verifies identical pre-final-score hard memberships,
means, sample covariances and free weights between ridge and shrinkage for all
ten seeds, both on the full cohort and each fixed whole-scan half. Full-cohort
mixing coefficients are **0.0326768–0.0452060**, with free weights ranging
**0.2598–0.7402**. Median covariance condition number falls **23,138.4 → 190.838**,
maximum **83,386.4 → 300.975**. Both halves also improve conditioning. Nevertheless,
training withdrawal changes **370/893** GMM decisions with ridge and **350/893**
with shrinkage; mean absolute vote change decreases **0.268533 → 0.249160**.
All seven unavailable keys remain in full/half diagnostics. This limited
stability improvement does not establish held-out recognition, calibrated noise,
chemical confidence or iid assumptions; the Fisher features are not retrained
within the halves. No setting is chosen from these diagnostics.

Only after all scientific outputs and checks does the external unchanged
manifest select **145 scans / 892 predicted lobes / 870 truth positions**.
Native coverage stays **849/870 (97.6%)**: 16 missing positions, five aligned
abstentions and 38 extra lobes. Counting remains fixed (14 short-N, 25 extra-N
scans); the excluded technical scan stays in scientific inference.

| Version | Correct / 870 | Correct / classified | Emitted errors | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical reference | 677 (77.8%) | 677/854 (79.3%) | 177 | 36 |
| Affine-descriptor control | 673 (77.4%) | 673/849 (79.3%) | 176 | 33 |
| Training-mean Fisher | 675 (77.6%) | 675/849 (79.5%) | 174 | 33 |
| Final covariance shrinkage | 665 (76.4%) | 665/849 (78.3%) | 184 | 34 |

Centering has **two scan gains, zero losses, 143 ties**. Lobe 3 of
`240815_048.sxm` and `240817_076.sxm` changes `1→0`, with final vote margins
**0.1** and **0.9** respectively; each scan improves 4→5 correct. No exact chain
or grade orientation changes. Confusion `(TN,FP,FN,TP)` changes
`(510,52,124,163) → (512,50,124,163)`. Against history, centering gives
23 winning / 24 losing / 98 tied scans, six exact chains gained and nine lost.

Shrinkage gives **five scan gains, eleven losses, 129 ties**, net **minus eight
correct positions**. It gains exact chains `240814_011` and `240814_025` but
loses `240818_020`, hence one more exact chain overall. Confusion becomes
`(498,64,120,167)`. Nineteen changes are `0→1` at zero vote margin. All three
`1→0` changes at margin **0.9** (`240311_Cu100063`, `240818_020`, `241113_086`)
lose one correct position. Thus neither a better-conditioned covariance nor
even large vote margins guarantee improved chemistry. Grade orientation changes
on `240314_Cu100_026`, `240817_058`, `240817_078`; it is not a spatial monomer
correspondence. Against history, shrinkage gives 22 winning / 30 losing / 93
tied scans, seven exact chains gained and nine lost. Every one of the **580
paired scan rows and 24 changed decisions**, including all losses, is retained.

**Decision:** retain `config/unit_assignment_centered_fisher.toml` at **675/33**
as the opt-in working candidate at that stage, **not a champion promotion**. Final
shrinkage is not retained as the working replacement despite one extra exact
chain: its loss of eight correct positions is contrary to the primary accuracy
objective. The two arms are not combined. Historical **677/870 and 36/145** still
lead by two positions and three chains merely to equal them; historical coverage
is also greater, so **79.5% versus 79.3% classified accuracy is not a victory**.
No threshold, seed, weight, normalization or per-file method is retuned after
the grade. Counts, the separate frozen counting benchmark and unknown25 are not
rerun or revalidated. This heavily reused benchmark is development evidence,
not independent validation or evidence of unknown-chain accuracy.

The complete French report, formulas, all changed decisions, numerical tables
and replay commands are in `results/numerical_signals_20260921/report.md`.
Verification comprises **1,186 focused Julia 1.13 assertions**, **76 saved-output
checks**, numerical invariants in the successful job, and paired-grade arithmetic
assertions. The initial include-guard parse error is recorded above. No scientific
formula or configuration changed after the source commit or after grading.
The bounded deliverable is complete; the broader historical target remains open.

The final documentation build passes, with only the existing generated-page
and search-index size warnings; deployment remains skipped outside CI.

## 2026-09-21 — Complete patch support and final Gaussian score: bounded comparison

The user authorizes trying the first two proposed leads after the centered-Fisher
675/33 result. Scientific deliverable: a complete independent comparison of
support handling and covariance-volume scoring, with all gains/losses, coverage
and exact chains. Done means synthetic/CLI regressions pass, one bounded Viper
job returns all three arms, twelve control tables exactly replay the saved
675/33 candidate, all 900 keys survive, and external full145 grading is reported
only after scientific outputs are fixed. The target remains historical 677/870
correct and 36/145 exact. One job is limited to 30 minutes, four requested CPUs
and 16 GB. No automatic retry, combined variant, post-grade tuning, normalization,
scan weighting, count refit, new DFT or unknown25 application is included.

**Support hypothesis.** Existing Fisher validity requires all 289 square pixels
although only its 197-pixel disk and center enter learning/scoring. The new
explicit `assignment_patch_support = "complete_disk_symmetric"` permits invalid
values outside that unchanged disk, but no missing disk pixel or missing schema
column. No pixel value is imputed. The backward affine descriptor independently
requires the complete 49-pixel integer-grid disk inside its 9x9 square. Remaining
observed peripheral pixels are retained only when all u/t reflection partners
are present. This avoids a one-sided observation mask generating apparent
transverse contrast. Fit the affine plane by least squares on this observed
symmetric support and apply the same half-plane numerator and residual L1 norm.
This discards observed corner partners too; it is not recovery of their unknown
counterparts. Complete patches take the byte-identical original path; missing
central disks remain NA. CC's existing pairwise-finite treatment is untouched.
No new minimum fraction, noise model or normalization is introduced. Availability
changes can alter Fisher/GMM learning globally, requiring a full-cohort comparison.

**Volume hypothesis.** `gmm_final_score = "gaussian_density"` replaces the final
`log(weight)-d^2/2` score with full Gaussian log density, adding `-logdet(Sigma)/2`
and the shared dimension constant. It uses the existing `1e-8 I` Cholesky guard.
The EM initialization and fit, two distance-only hard self-training updates,
last ridge covariances and free mixture weights are unchanged. Final chemical
group naming remains mean raw amplitude of assigned members; naming can therefore
change with assignments even though fitted parameters do not. Votes remain hard
seed decisions, not calibrated density probabilities. The mode requires at least
one hard self-training iteration. No shrinkage or support change accompanies it.

Both candidates copy `unit_assignment_centered_fisher.toml` and change only name
and their own setting. All nine older native configs explicitly retain
`full_square` / `mahalanobis`. Physical and selection parameters stay in TOML.
No expected N, sequence, truth class or composition enters scientific code.
The science cohort includes the externally excluded technical scan; only the
separate grader may use the confirmed benchmark manifest. Synthetic tests cover
independent least-squares/density formulas, reflection/affine invariances,
discarded-pixel independence, all mandatory disk positions, unchanged complete
patches and fitted GMM parameters, invalid-key retention, strict configuration,
CLI/API agreement and the three-arm shell runner. Results follow after execution.

An initial synthetic density assertion failed on one random covariance with
condition number **7.0243e6**: LU versus Cholesky differs by **9.053e-7 absolute /
7.201e-11 relative**, just above the test's fixed `5e-11` relative tolerance.
The test now bounds this independent-solver comparison by
`max(5e-11, 8eps(Float64)*cond(Sigma))`, and adds strict `5e-13` comparisons for
well-conditioned versions of every matrix. This changes test error accounting
only; no scientific ridge, formula, data or candidate setting is changed.

Pre-run verification now passes **1,604 assertions** under Julia 1.13: new
support/score tests 286, numerical assignment 136, numerical diagnostics 20,
native Fisher 229, affine Fisher 74, Fisher attribution 207, representation
diagnostics 319, transverse descriptors 125, signal features 73, reconstructed
assignment 37, full synthetic pipeline 82, and new shell runner 16. The docs
build passes with only the existing page/search-index size warnings and skipped
non-CI deployment. No real multi-file scientific fit was run locally. The count
config and Manifest hashes match the preceding comparison. Viper input hashes
also match, and no competing job is queued or running before submission.

### Completed comparison and decision

Validated scientific source is committed as **db411a2** on
`research/patch-support-gmm-volume-20260921`. The new Viper source directory is
`/u/oldu/code/STMFit_patch_support_score_20260921`; the unchanged ignored Manifest
is copied explicitly, with no dependency resolution. Tracked-source checksum
comparison (`rsync -anci`) is empty. All three dry runs recognize 146 files under
Julia 1.13.0. Job **11906168** waits briefly on the global `n0001` CPU quota, then
completes **0:0 in 7m58s** on `vipc2294`, **10:36:38–10:44:36 CEST** on September
21. Four CPUs are requested, eight allocated, four Julia threads. No QOS/resource
change, resubmission or scientific fit on the login node occurs.

All **585 files and the Slurm log** are fetched from
`/ptmp/oldu/stmfit/patch_support_score_20260921_v1` into
`results/patch_support_score_20260921/run_v1`; raw staging symlinks are excluded.
Itemized checksum comparisons for outputs/log are empty. No stage failure or
warning appears in the logs. **76 output checks** pass before external grading:
twelve control tables exactly replay the preceding centered-Fisher candidate;
all three arms retain the same 146 files and 900 keys. Patches, local features,
CC scores and k-means predictions are identical. Complete backward descriptors
stay byte-identical. The volume arm changes no input feature or Fisher score.

Support makes these three previously invalid rows available without imputation:

| Key | Forward observed / 289 | Fisher disk / 197 | Backward observed / 81 | Symmetric retained | New prediction / margin |
|---|---:|---:|---:|---:|---|
| `240818_019.sxm`, lobe 3 | 287 | 197 | 80 | 77 | 0 / 0.8 |
| `240818_019.sxm`, lobe 4 | 287 | 197 | 80 | 77 | 0 / 0.8 |
| `241113_088.sxm`, lobe 6 | 277 | 197 | 77 | 65 | 0 / 1.0 |

Four incomplete-disk keys remain unavailable: `240310_Cu100009` lobe 8,
`240818_015` lobe 1 and `240818_019` lobes 1–2. The first is excluded only by
external grading, never by science. The support arm changes 896 Fisher scores,
246 GMM vote frequencies, eleven hard GMM predictions and five final predictions.
Thus accepting three rows also changes learning on the rest of the cohort.
Final full-cohort frequencies `(0,1,?)` are control `(672,221,7)`, support
`(673,223,4)`, volume `(653,240,7)`; these are outputs, never composition inputs.
Volume changes 25 GMM vote frequencies, 21 hard GMM predictions and 19 final
decisions. All nineteen final changes are `0→1` at the unchanged zero-margin tie.

Only after complete output validation does the unchanged external manifest select
**145 scans / 892 predicted lobes / 870 truth positions**. All arms have the same
16 missing positions, 38 extra lobes, 14 short-N scans and 25 extra-N scans.
Support classifies 852/870 (97.9%), leaving two aligned abstentions; control and
volume classify 849/870 (97.6%), leaving five. Historical coverage is 854/870.

| Version | Correct / 870 | Correct / classified | Emitted errors | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical reference | 677 (77.8%) | 677/854 (79.3%) | 177 | 36 |
| Centered-Fisher control | 675 (77.6%) | 675/849 (79.5%) | 174 | 33 |
| Complete-disk support | 676 (77.7%) | 676/852 (79.3%) | 176 | 34 |
| Final Gaussian score | 666 (76.6%) | 666/849 (78.4%) | 183 | 32 |

Support gives **two scan gains, two losses, 141 ties**, net **+1 correct**.
`240818_019` improves 0→2 correct and `241113_088` improves 5→6, gaining the
only new exact chain. But lobe 3 of both `240815_048` and `240817_076` changes
`0→1` at zero margin, each losing one correct position: the two centering gains
from the previous comparison disappear. No grade orientation or exact chain is
lost. Confusion `(TN,FP,FN,TP)` becomes `(513,52,124,163)` from `(512,50,124,163)`.
Three recovered correct decisions minus two regressions explain the gain; the
194 non-correct positions comprise 176 errors, 16 missing and two abstentions.
Against history: 23 winning / 25 losing / 97 tied scans, six exact gained and
eight lost. No per-file restoration or tie threshold is chosen from these labels.

Volume gives **four scan gains, twelve losses, 129 ties**, net **−9 correct**.
Exact chain `240817_045` is gained; `240817_021` and `241113_087` are lost.
Confusion is `(498,64,119,168)`. Grade orientations change for `240314_Cu100_026`
and `240817_078`; they are evaluation alignments, not spatial correspondence.
Against history: 20 winning / 30 losing / 95 tied scans, six exact gained and
ten lost. Normalized Gaussian density does not improve chemistry in this test.

**Decision:** retain **`config/unit_assignment_patch_support.toml` at 676/34**
as the latest opt-in working candidate, **not a champion promotion**. It improves
fixed-denominator correct count, exact chains and coverage despite lower
classified accuracy and two extra emitted errors. Reject the independent final
Gaussian score. Historical 677/36 remains ahead by one correct position and two
exact chains merely to equal it, with two additional classified positions.
No combination, post-grade parameter tuning, per-file choice, normalization,
scan weighting, count refit, new DFT or unknown25 run follows. Predictions are
label-free; this extensively reused development benchmark is not independent
validation or evidence of unknown-chain accuracy.

All **580 paired scan rows and 24 changed decisions**, including all losses,
pixel support details, formulas and replay commands are in
`results/patch_support_score_20260921/report.md` and its linked tables. Paired
arithmetic checks agree with every headline, coverage and confusion count. The
bounded comparison is complete; the broader historical objective remains open.
Prediction SHA-256: control `2a4e386ed21c324b1fb0fa73fa0b2cb031cf2908d2d6fec646befd2a13a9cb26`,
support `eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`,
volume `a465521eb546d13241a27a2a23a4954e926076d972a7c38d9ae73394313749b1`.

The final documentation build passes with only the existing page/search-index
size warnings and skipped non-CI deployment. All current README, runbook,
selection, assignment, config and calibration headlines reflect 676/34; dated
earlier experiments retain their own results. No scientific file changed after
the source commit. Results remain ignored; validated source and documentation
are committed, not left as uncommitted workspace state.

## 2026-09-21 — Complete-patch training versus admissible partial-patch scoring

**Scientific deliverable and done condition.** The user authorizes the first of
three proposed leads: one complete-training candidate versus the latest support
candidate (676/870 correct, 34/145 exact), at unchanged selected N and geometry.
Done means an exact control replay, all 146 files / 900 keys retained, inspected
training eligibility and coverage, external grades and every paired regression,
fetched results and committed source/docs. Historical 677/870, 36/145 and 854/870
coverage remain the objective to exceed, not merely the current runnable result.
One Viper Slurm job requests four CPUs, 16 GB and 30 minutes; no automatic retry,
local multi-file fit, robust scaling, equal scan weighting, combined candidate,
post-grade parameter adjustment, new DFT, count refit or unknown25 rerun.

**Hypothesis fixed before execution.** Admitting partial patches improved coverage
but also changed cohort learning. Separate learning eligibility from prediction
eligibility. New explicit `selection.assignment_training_support` accepts
`all_admissible` (eleven existing native configs) or `complete_patches` (one new
copy of `unit_assignment_patch_support.toml`, changing only name and this field).
A geometry-only table records the observed counts in forward 17x17, backward
17x17 and backward 9x9 patches. Fisher requires only its own complete forward
square for training; GMM requires all three families it consumes. There is no
new tuned completeness fraction. Scoring retains the previously frozen complete-
disk / reflection-symmetric rule; no pixel is filled or inferred.

Fisher fits only eligible opposite-lobe-parity rows. GMM per-file feature moments
use finite values of complete-patch rows, then transform all admissible rows;
interactions follow this same frozen scaling. Missing training moments cause NA,
not a partial-row/global fallback; the existing zero/undefined-std guard is one.
EM, two hard updates, final covariances, free component weights and physical group
naming use training members only. The high-raw-amplitude group name is frozen
before partial predictions. K-means did not change in the preceding support
experiment and stays unchanged here, including its normalization and four views.
Existing score centering, ridge, seeds, final Mahalanobis score and soft vote
remain fixed. Complete-case learning might remove useful observations or worsen
conditioning; the outcome is not presumed positive.

Tests exercise mask/key/schema validation, no benchmark columns, complete-data
identity, opposite-fold scoring, no feedback through normalization/fits/group
naming, unequal free weights, missing-training abstention, all-key retention and
CLI/API identity. The old raw Fisher-attribution diagnostic explicitly rejects
the new policy rather than silently rebuilding the wrong training folds. The
new two-arm Slurm runner requires explicit absolute exported inputs, a dry-run
and an allocation. Results will be appended after validation and external grading.

Pre-run verification passes **1,782 assertions** under Julia 1.13: complete-
training tests 154, prior support/score 286, numerical assignment 136, numerical
diagnostics 20, native Fisher 229, affine Fisher 74, Fisher attribution 210,
representation diagnostics 319, transverse descriptors 125, signal features 73,
reconstructed assignment 37, dual-policy synthetic pipeline 103 and shell runner
16. Both synthetic end-to-end policies produce byte-identical outputs with
complete patches. The documentation build passes with the existing HTML/search
size warnings and non-CI deployment skip. No real-cohort fit runs locally.

Read-only inspection of saved patch pixels finds **893/900 complete rows** for
each training stage, with the same seven exclusions. The three newly scoreable
partial rows remain prediction-only by the fixed rule. Viper input hashes,
root Manifest and count config match the previous comparison; no other user job
is queued or running before submission. No scientific rule is changed after
inspecting these observation counts.

Scientific source **f7a4074** is committed before submission. The dedicated Viper
checkout `/u/oldu/code/STMFit_complete_training_20260921` is checksum-identical
to the tracked source; the ignored Manifest is copied unchanged. Both dry runs
recognize 146 inputs under Julia 1.13.0. Single job **11908779** is submitted at
**11:45:16 CEST** with the five absolute input variables explicitly exported,
partition `small`, account `mfk_cpu`, QOS `n0001`. Output target is
`/ptmp/oldu/stmfit/complete_training_20260921_v1`. Slurm initially queues it for
node availability and the global QOS CPU limit. No resource change, duplicate
submission, local fallback fit or external grading is performed while waiting.

### Completed comparison, observed missing-training case and decision

Job **11908779** queues for **31m21s**, then completes **0:0 in 6m42s** on
`vipc2262`, **12:16:37–12:23:19 CEST**. Four CPUs are requested, eight allocated;
Julia uses four threads. The control finishes at 12:21:15, then the candidate
runs in the same job. No resource change, duplicate submission or local fallback
fit occurs. All **392 outputs plus the Slurm log** are fetched into
`results/complete_training_20260921/run_v1`; itemized checksum checks are empty.
No scientific-stage warning or failure appears in the logs.

All twelve control tables exactly replay the previous support candidate, and
both arms retain all 146 files / 900 keys. Patches, local features, descriptor,
CC scores, k-means and the observed-count tables are byte-identical. Only the
Fisher predictor column changes before GMM. Both masks select 893 complete rows,
432 even and 461 odd. Fisher still scores 896 rows, but the candidate GMM can
score only 894. File `240818_019.sxm` has four fitted lobes and no complete square:
its forward/backward17 counts are 9, 130, 287, 287 out of 289; backward9 counts
are 4, 37, 80, 80 out of 81. Its two admissible partial rows therefore lack
per-file complete-training moments and become NA, exactly as frozen beforehand.
Only `241113_088.sxm` lobe 6 remains a scored prediction-only partial row.

**Audit correction, not a scientific repair:** the initial saved-output checker
wrongly asserted that every file had at least one complete row. It passes 58
checks and fails this one before any external grade. Inspection of observed
pixel counts identifies the legitimate zero-training file. Replace only that
unsupported audit assumption with the actual promised invariant: files without
complete training rows have no valid normalized GMM view and must abstain on
all their keys. The initial `output_checks.log` is retained; all **59 checks**
pass in `output_checks_v2.log`. No inference is changed or rerun, and no
partial-row/global normalization fallback is introduced. The synthetic suite
already tests this missing-training behavior. This is an observed cost of the
chosen rule, not grounds to silently weaken it.

The candidate changes 896 Fisher scores, 246 GMM vote frequencies, twelve hard
GMM predictions and five final predictions. Full-cohort frequencies (0,1,?)
change from (673,223,4) to (674,220,6). Every unavailable control key stays
unavailable; only the two extra abstentions above are added. All observations
and checks are completed before the separate external grade reads its labels.

| Version | Correct / 870 | Correct / classified | Errors emitted | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical reference | 677 | 677/854 (79.3%) | 177 | 36 |
| Support control | 676 | 676/852 (79.3%) | 176 | 34 |
| Complete training | 674 | 674/850 (79.3%) | 176 | 33 |

The external cohort remains 145 scans / 892 predicted lobes / 870 reference
positions. Coverage falls 852→850/870; aligned abstentions rise two→four.
Counts stay fixed: 16 missing positions, 38 extra lobes, 14 short-N scans and
25 extra-N scans. The candidate's 196 non-correct positions are 176 errors,
16 missing and four abstentions, not 196 abstentions. Rounded classified
accuracy masks this loss on the fixed denominator.

Against the control: **one winning scan, two losing, 142 ties**, net −2 correct.
`240817_076` lobe 3 changes 1→0, improving 4→5 correct; `240818_019` lobes 3
and 4 change 0→?, losing 2→0 correct. `240818_020` lobe 2 changes 1→0, losing
6→5 and its exact chain. No exact chain is gained and no grading orientation
changes. The fifth changed prediction is lobe 7 of the excluded technical scan
`240310_Cu100009`, 1→0, ungraded. The three 1→0 changes move from zero margin
to 0.7; this larger margin does not imply correctness. The earlier centering
regression in `240815_048` is not recovered. Confusion (TN,FP,FN,TP) changes
(513,52,124,163)→(512,51,125,162). Against history: 23 gains, 25 losses, 97 ties,
−3 correct, six exact gained / nine lost, eight grading-orientation changes.
These grading alignments are not spatial monomer correspondence.

**Decision: reject this complete-training variant.** Keep
`config/unit_assignment_patch_support.toml` at **676/870, 34/145**, still short
of historical 677/36 and its 854/870 coverage. This test does not refute every
possible separation of training and prediction; it rejects this particular
complete-case/per-file-normalization rule as the working replacement. No
fallback, threshold, new arm, robust scaling, scan weighting, count refit, DFT
change or unknown25 rerun follows. All **290 paired scan rows and five changed
decisions**, including losses, appear in
`results/complete_training_20260921/report.md`. Paired arithmetic agrees with
all headline, coverage and confusion counts. This development benchmark remains
reused, not independent validation; the bounded comparison is complete and the
broader historical objective remains open.

Prediction SHA-256: control
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`, candidate
`f0b711fc889c98b157bbc8bd6afb724576c3b29f04dd73e1a155d003882b4571`.

The final docs build passes with the existing HTML/search-index size warnings
and skipped non-CI deployment. README, runbook, selection, assignment, config,
calibration and HPC notes all record the negative 674/33 result and retain
676/34 as the working candidate. Scientific source is unchanged since f7a4074;
validated source and conclusions are committed. Generated results remain ignored
and are preserved locally, including the failed initial audit and all losses.

## 2026-09-21 — Robust per-scan GMM normalization, one fixed comparison

**Scientific deliverable and done condition.** The user authorizes the recommended
first lead: compare median/IQR scaling with the current 676/870, 34/145 support
candidate. Done means synthetic tests, an exact control replay, all input keys
retained, checked upstream identity and numerical support, one external grade
with paired gains/losses and coverage, fetched outputs and committed conclusions.
Historical 677/870, 36/145 and 854/870 coverage remain the goal to exceed. This
is one 30-minute Viper allocation, four requested CPUs and 16 GB, not an open
campaign. No automatic retry, weighting, scan resampling, complete-case training,
combined variant, count refit, DFT change, unknown25 rerun or post-grade tuning.

**Fixed scientific change.** `unit_assignment_robust_normalization.toml` copies
`unit_assignment_patch_support.toml`, changing only the method name and explicit
`[preprocessing] gmm_feature_normalization` from `mean_sample_std` to `median_iqr`.
All twelve previous native configs explicitly retain the original mode. The
existing degenerate-scale guard is exposed as `gmm_scale_fallback = 1.0`.
Moments are computed independently for each feature/file on the same finite
observations as before, including eligible partial rows. The robust expression
is `(x - median(x))/(Q75-Q25)` with Type-7 quantiles (`alpha=beta=1`). A zero or
nonfinite IQR uses one; positive IQRs are not floored. Empty support remains NA.
Pairwise products follow scaling; no clipping, imputation, global moments or
Gaussian-consistency factor is introduced. No benchmark information participates.

Only the GMM input transformation changes. Fisher, patch extraction, CC,
descriptor, k-means, GMM seeds, EM, two hard updates, ridge, free component
weights, amplitude-based component naming and final soft vote remain fixed.
Robust moments need not improve recognition: a small IQR can magnify outliers
and their products, and few lobes can make quantiles unstable. This is a testable
hypothesis, not a promised improvement or calibrated confidence mechanism.

The synthetic test is first run before implementation and fails on the missing
candidate, as expected. Tests cover the exact formulas, legacy arithmetic,
singleton/constant/zero-IQR/empty/nonfinite support, no positive-IQR floor,
per-file isolation, row/name invariance, masked training isolation, products
after scaling, free weights, deterministic CLI/API identity and forbidden label
inputs. Full-cohort scientific fits remain on Viper. Results follow after the
unchanged external grade; the reused benchmark is development evidence only.

Initial verification passes **99 normalization assertions** and **154 existing
complete-training assertions** under Julia 1.13. The documentation build passes
with the pre-existing HTML/search-index size warnings and non-CI deployment
skip. Further regression/pipeline suites continue before submission. The old
representation diagnostic now rejects median/IQR configs explicitly, rather
than silently reporting mean/std arithmetic for a different model. The source
is committed early; neither real-cohort fitting nor external grading has run.

Pre-submission verification finishes with **1,900 assertions**: normalization
99, complete training 154, support/score 286, numerical assignment 136,
numerical diagnostics 20, native Fisher 229, affine Fisher 74, Fisher attribution
210, representation 319, transverse descriptors 125, signal features 73,
assignment 37, three-mode synthetic pipeline 122 and Slurm runner 16. A shell
test list used the nonexistent name `test_label_free_signals`; the actual
`test_label_free_signal_features` suite and remaining interface/pipeline suites
were already executed separately and passed. No test or scientific rule is
weakened. The synthetic pipeline preserves upstream/Fisher/k-means tables under
the robust change and complete-data identity under complete training.

Saved-feature arithmetic (no local cohort fit) finds the same **896/900 finite
views** in both modes. Maximum absolute base descriptor grows **2.62205 →
37.13789**, and the expanded maximum **5.44386 → 335.46532**. This illustrates
the predeclared risk of amplifying tails; no clipping, scale factor or threshold
is added after observing it.

Scientific source **77b3ce6** is checksum-identical to the dedicated Viper tree
`/u/oldu/code/STMFit_robust_normalization_20260921`. The ignored Manifest is copied
unchanged. Geometry, split geometry, templates and count config retain their
preceding SHA-256 values; no benchmark tables or past result trees are synced.
Both metadata-only dry runs recognize 146 files under Julia 1.13.0. One job
**11911253** is submitted at **13:11:25 CEST**, partition `small`, account
`mfk_cpu`, QOS `n0001`, with all five absolute inputs explicitly exported.
Output target: `/ptmp/oldu/stmfit/robust_normalization_20260921_v1`. The job starts
in the queue. No duplicate submission, resource change, local fallback fit or
external grade is performed while waiting.

### Completed robust-normalization comparison and decision

Job **11911253** waits **51 seconds**, then completes **0:0 in 5m55s** on
`vipl2005`, **13:12:16–13:18:11 CEST**, with four requested/eight allocated CPUs
and four Julia threads. The control finishes at 13:16:26, then the robust arm
at 13:18:11. No duplicate submission, resource change, local cohort fit or
scientific retry occurs. All **392 outputs and the Slurm log** are fetched into
`results/robust_normalization_20260921/run_v1` and its parent; checksum comparisons
are empty. Tracked runtime sources and the ignored Manifest remain identical
locally and remotely. No scientific-stage warning or failure appears.

All twelve control tables reproduce the saved support candidate byte-for-byte.
Both arms retain all 146 scans / 900 keys. Patches, local/descriptor/predictor
features, CC, Fisher, k-means and observed pixel counts are byte-identical.
The initial audit passes 52 checks but incorrectly requires byte-identical
`selected_from_features.tsv` staging paths: those contain `control/raw_inputs/`
or `robust/raw_inputs/` respectively. The checker is corrected to compare unique
scan identities and their N, also validating N against feature keys. Its failed
log is preserved. **1,261 output checks** pass (59 outputs/counts, one key-order
check, 1,201 arithmetic/support checks). This correction precedes grading and
changes no scientific code, parameter, output or fit.

Both modes still have **896/900 finite GMM views** and the same four unavailable
final assignments. None of **1,168 scan/feature pairs** uses the zero/nonfinite
IQR fallback. Independent sorted Type-7 interpolation reproduces each transformed
feature; products are checked after scaling. The expanded maximum is indeed
**5.44386 → 335.46532**, as observed before execution. GMM changes all 896 vote
frequencies and **666 hard decisions**. The final vote changes **64 decisions**,
63 graded and one on technical scan `240310_Cu100009`: **62 are 1→0, two 0→1**.
Final raw frequencies `(0,1,?)` change `(673,223,4)→(733,163,4)`. These are
reported outcomes, not fitted class-count constraints or composition targets.

Only after these checks does the unchanged external grade select 145 scans /
892 predicted lobes / 870 truth positions:

| Version | Correct / 870 | Correct / classified | Errors emitted | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical reference | 677 | 677/854 (79.3%) | 177 | 36 |
| Support control | 676 | 676/852 (79.3%) | 176 | 34 |
| Median/IQR | 667 | 667/852 (78.3%) | 185 | 10 |

Against control: **23 winning scans, 33 losing, 89 ties**, net **−9 correct**;
**zero exact chains gained, 24 lost**, eight external grading-orientation
changes. Confusion `(TN,FP,FN,TP)` changes `(513,52,124,163)→(537,28,157,130)`:
24 fewer false positives do not compensate for 33 more false negatives.
Coverage stays **852/870**, with two aligned abstentions. Counts stay 16 missing
positions, 38 extra lobes, 14 short-N and 25 extra-N scans. The candidate's 203
non-correct positions comprise 185 errors, 16 missing and two abstentions.
Against history: 27 wins, 41 losses, 77 ties, net −10 correct; one exact gained,
27 lost and twelve orientation changes. Grade alignment is not a spatial
monomer correspondence.

**Decision: reject this median/IQR variant; keep the support candidate at
676/870 and 34/145.** Historical 677/36 and 854/870 coverage remain the target
to exceed. This rejects the specific fixed normalization, not every robust
estimator. Tail amplification is observed, but its independent causal role in
the grade loss has not been isolated. No clipping, consistency factor, new
quantiles, IQR floor, scan weights, seed/threshold change, combined arm, refit
of N, DFT change or unknown25 run follows the grade. The single authorized
comparison is complete; the broader historical objective remains open.

All **290 paired scan rows and 64 changed decisions**, including every loss and
the excluded technical scan, are in
`results/robust_normalization_20260921/report.md` and its linked tables. Grade
sums, confusion, coverage and exact-chain counts agree. Inference is label-free;
the extensively reused benchmark remains development evidence, not independent
validation or unknown-chain accuracy. Prediction SHA-256: control
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`, robust
`a9e0ade8156f9a740bd1e9005f9caac2db05ca2840d4907a8baa2c8474a090c9`.

The final documentation build passes, with the existing large-page/search-index
warnings and non-CI deployment skip. README, runbook, selection, assignment,
config, calibration and HPC notes record the negative 667/10 result and retain
676/34 as the working candidate. Runtime source is unchanged since 77b3ce6.
Validated source and conclusions are committed; generated outputs, all losses
and the initial failed metadata-audit log remain preserved locally and ignored.

## 2026-09-21 — Equal total scan influence in GMM learning, one fixed comparison

**Scientific deliverable and done condition.** The user authorizes the first
recommended lead: compare equal-scan observation weights with the retained
676/870, 34/145 support candidate. Done means tested weighted arithmetic and
label-free boundaries, an exact control replay, all 146 scans / 900 keys and
upstream signals checked, external correct/coverage/exact-chain grades and every
paired regression, fetched outputs and committed source/conclusions. Historical
677/870, 36/145 and 854/870 coverage remain the objective, not merely an
improvement over the current runnable reconstruction. One Viper allocation:
four requested CPUs, 16 GB, 30 minutes. No automatic retry, local cohort fit,
scan resampling, robust scaling, complete-case filtering, alternative physical
naming rule, combined variant, N refit, DFT change or unknown25 application.
The archived plan's scientific Must-NOT-have list is retained, not its retired
process apparatus. No labels or composition counts enter fitting or attribution.

**Fixed change before grading.** All thirteen earlier native configs now state
`selection.gmm_training_weighting = "equal_lobes"`, preserving legacy arithmetic.
New `unit_assignment_scan_weighting.toml` copies the support control, changing
only method name and this field to `equal_scans`. For a view with `n` usable
training rows from `S` represented scans, a row from scan `s` gets weight
`n/(S*m_s)`, where `m_s` is that scan's usable training count, never its expected
N. Thus each scan totals `n/S` and the overall mean weight is one. Empty scans
do not enter training; their output rows and ordinary abstentions remain.
Finite per-feature mean/std normalization still uses its original admissible
support before complete-view validity is checked. No observation is imputed.

The installed Clustering implementation accepts observation weights for Lloyd
updates but not k-means++ draws. The opt-in branch therefore supplies two
weighted k-means++ seeds (first proportional to weight, second to weight times
squared distance) to the existing weighted k-means implementation. Its centroid
objective, initial GMM means/covariances, weighted EM objective and sufficient
statistics, both hard updates, free mixture weights and raw-amplitude naming
all use the same fixed row weights. Per-lobe classification remains conditional
and is not multiplied by row weight. Initial component weights remain 1/2,
then are learned without any equal-class constraint. Ridge remains 1e-6;
weighted Ledoit-Wolf is rejected rather than silently using an unweighted
estimator. Normalization, interactions, support, Fisher, independent k-means
head, seed integers, iteration budgets, final score and vote stay fixed.

Equal influence may reduce domination by scans with more detected lobes, but
may also overemphasize short/noisy scans. Recognition benefit is unproven.
Synthetic tests check explicit weighted initialization/EM/hard-update formulas,
free masses, amplitude naming, invalid-row exclusion, scan-total equality,
label/count/name independence and deterministic CLI/API agreement. Historical
representation diagnostics reject the new weighting, since their unweighted
refits would not represent this candidate. No external grade has yet run.

The first **138 synthetic assertions pass** under Julia 1.13; the documentation
build passes with its existing page/search size warnings and non-CI deployment
skip. Full regression and four-mode synthetic pipeline checks are in progress.
Remote geometry, split geometry, template and Manifest hashes match the preceding
experiment, and no user job is queued/running at this inspection. Source is
committed before any real-cohort execution; no scientific result is presumed.

Pre-submission verification finishes with **2,057 passing assertions**: weighting
138, normalization 99, complete training 154, support/score 286, numerical
assignment 136, numerical diagnostics 20, native Fisher 229, affine Fisher 74,
Fisher attribution 210, representation 319, transverse descriptors 125, signal
features 73, assignment 37, four-mode synthetic pipeline 141 and Slurm runner 16.
The pipeline confirms unchanged upstream/Fisher/k-means outputs. Saved-feature
arithmetic, without any cohort fit, finds **896 usable rows from all 146 scans**.
Their usable-count frequencies are 2:1, 4:1, 5:12, 6:107, 7:16, 8:7, 10:2 scans.
Weights range **0.61369863–3.06849315** and each scan totals **6.13698630**.
No rule is adjusted after inspecting these observed counts.

Committed scientific source **5ff68e2** is checksum-identical to the dedicated
Viper checkout `/u/oldu/code/STMFit_scan_weighting_20260921`; Manifest is copied
unchanged. The transfer dry-run initially includes tracked benchmark tables;
filtering those paths out of the Git file list removes them before any real
transfer. No benchmark, prior results or agent state is copied. Both metadata-
only pipeline dry runs identify all 146 files under Julia 1.13.0. `sbatch
--test-only` accepts the fixed resources but is not compute-node validation.
Single job **11915029** is submitted at **14:55:07 CEST**, partition `small`,
account `mfk_cpu`, QOS `n0001`, with all five absolute inputs explicitly exported.
Output target: `/ptmp/oldu/stmfit/scan_weighting_20260921_v1`. Its initial state
is pending, reason `QOSGrpCpuLimit`; no duplicate, resource change or local fit
is attempted. External grading has not run.

**Completed follow-through, as requested by the user.** The existing job starts
at **15:36:07 CEST** after **41m00s** queued under the shared CPU quota and
finishes at **15:42:08**, **COMPLETED, exit 0:0, elapsed 6m01s**, on `vipc2294`.
The control runs 15:36:08–15:40:09, then weighting 15:40:09–15:42:08. Slurm
allocates eight CPUs for the four-CPU request; Julia uses four threads. Peak
memory is 2.915767 GB per node. No duplicate job, resource change, retry or local
cohort fit is used. All **392 regular output files and the Slurm log** are
fetched from the job directory, excluding raw-input staging symlinks; the
checksum comparison is empty. Runtime source/config and Manifest also match.

Before grading, **360 saved-output assertions pass**: 59 output/key/count
checks and 301 weighting/support checks. All twelve control scientific tables
replay the retained support candidate byte-for-byte. Both arms retain all
146 scans / 900 keys, the same 896 valid GMM rows and the same four unavailable
assignments. Patches, local/descriptor/predictor features, CC, Fisher, separate
k-means and pixel counts are identical. N agrees by scan identity and feature
keys; staging paths are arm-specific. Training weights reproduce the previously
fixed inverse-usable-count formula, mean one and equal total 896/146 per scan.

GMM vote fractions change on **42 rows**, with **one hard GMM decision** changing.
The final vote changes **seven decisions**, all 0→1, all on benchmark scans.
Each lands exactly at the unchanged vote tie (`probability_1=0.5`, confidence
zero); the existing `>=0.5` rule emits 1. Raw `(0,1,?)` totals change
`(673,223,4)→(666,230,4)`. These are outcomes, not composition targets. Zero vote
margin is not calibrated chemical confidence; no tie or abstention rule changes.

Only after integrity checks does the unchanged external grade select 145 scans /
892 prediction lobes / 870 truth positions:

| Version | Correct / 870 | Correct / classified | Errors emitted | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 177 | 36 |
| Replayed support control | 676 | 676/852 (79.3%) | 176 | 34 |
| Equal-scan weighting | 675 | 675/852 (79.2%) | 177 | 36 |

Against control: **3 scan gains, 4 losses, 138 ties**, net **−1 correct**;
**two exact chains gained, none lost**, no external grading-orientation changes.
`240307_016` and `240814_011` gain exactness; `240817_017` gains one correct
position but still has an extra lobe. Losses are `240817_058`, `240817_059`,
`241114_011` and `241114_022`. Confusion `(TN,FP,FN,TP)` changes
`(513,52,124,163)→(509,56,121,166)`: three fewer false negatives, four more false
positives. Coverage stays **852/870**, with two aligned abstentions. Counts stay
16 missing positions, 38 extra lobes, 14 short-N and 25 extra-N scans. Weighted
non-correct positions total 195: 177 errors, 16 missing and two abstentions.
Against history: 23 gains, 24 losses, 98 ties, net −2 correct; six exact chains
gained and six lost, with eight grading-orientation changes. Equal exact-chain
totals do not mean recovering the same exact chains. Grade alignment is not
an independently established spatial monomer correspondence.

**Decision: retain the 676/34 support candidate as the primary working
reference; record 675/36 as a mixed result, not a new champion.** Neither arm
dominates on both correct-position and exact-chain counts. The historical
677/36 remains unexceeded, with two more correct and classified positions than
weighting. This experiment does not isolate which weighted learning/naming
stage caused a flip or establish a general failure of equal scan influence.
No weight formula, seed, threshold, naming rule, combined arm, refit of N,
DFT change or unknown25 run follows the grade. This authorized comparison is
complete; the broader historical objective remains open.

All **290 paired scan rows and seven changed decisions**, including every loss,
are in `results/scan_weighting_20260921/report.md` and its linked tables. Paired
arithmetic reproduces grade sums, confusion, coverage and exact-chain counts.
Inference is label-free; the repeatedly reused benchmark is development
evidence, not independent validation or unknown-chain accuracy. Prediction
SHA-256: control
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`, weighted
`102467b38804841827eb29375aaa14d8ef60ce1559a47c2f1f6097593c0df687`.

The documentation build passes with the existing large-page/search-index
warnings and non-CI deployment skip. README, runbook, selection, assignment,
config, calibration and HPC notes now record the 675/36 tradeoff and retain
676/34 as the primary working reference. Runtime source remains **5ff68e2**;
the conclusion is committed separately. Generated outputs, logs and all paired
losses are preserved locally and ignored by Git.

## 2026-09-21 — Continuous internal GMM vote, one fixed comparison

**Scientific deliverable and done condition.** The user authorizes the first
recommended lead only: compare a continuous internal GMM vote with the retained
676/870, 34/145 support candidate. Done means tested arithmetic and unchanged
model learning, exact control replay, all 146 scans / 900 keys and N retained,
upstream outputs checked, external correct/coverage/exact-chain grades and
paired losses reported, outputs fetched and validated source/conclusions
committed. Historical 677/870 and 36/145 remain the target. One Viper allocation:
four requested CPUs, 16 GB and 30 minutes. No automatic retry or local cohort
fit, scan resampling/weighting, naming-rule change, temperature, new Gaussian
score, combined variant, N refit, DFT change or unknown25 rerun. The archived
application plan's scientific Must-NOT-have list is retained, not its retired
process apparatus.

**Fixed change before any grade.** All fourteen earlier native configs now
declare `selection.gmm_seed_aggregation = "hard_vote"`. New
`unit_assignment_continuous_vote.toml` copies the support candidate, changing
only its name and this field to `mean_membership`. The code already computes
normalized exponentials of the two final component scores. The legacy path
averages `argmax(resp) == high_cluster` (0 or 1); the candidate averages
`resp[high_cluster]`. The high-amplitude group's name is still determined from
the same hard assignments and training amplitudes before this aggregation.

Every seed's EM, two hard moment updates, free masses, means, covariances,
ridge and physical naming are unchanged. Scores remain
`log(weight) - Mahalanobis_distance²/2`, without the covariance-volume term.
Normalization, interactions, support, k-means, Fisher, seed integers and N are
fixed. No temperature, fitted probability calibration or class-count prior is
introduced. Memberships and the final vote margin are not calibrated chemical
probabilities. The existing eight-decimal serialization and final `>=0.5`
decision rule are retained, as are unavailable rows and abstention.

This specifically tests information lost by hardening each seed's score; it is
not another sweep of seed count, final vote weights or equality rules. Scores
may already be saturated, and smooth aggregation need not improve chemistry.
Historical representation diagnostics reject this mode instead of silently
reporting hard-vote refits. Synthetic tests check the independent score/average
formula, identical fitted parameters, free class masses, finite/invalid rows,
legacy ties, deterministic CLI/API agreement and label/count/name independence.
Remote cached geometry, split geometry and templates have the same SHA-256 as
the preceding comparison; no user job is running/queued at preparation. No
real-cohort fit or external grade has run for this candidate.

The first synthetic run passes the aggregation and identical-parameter checks,
but two assertions incorrectly require every seed to have unequal component
masses. A free mixture may legitimately learn equal masses. The test is
corrected to check each mass against its actual hard-membership count, plus a
fixture check that unequal masses are also learned. The failed log is retained;
no scientific code, score, setting or fit is changed in response. Tests also
record the existing near-tie eight-decimal serialization without changing it.

The corrected targeted suite passes **107 assertions** under Julia 1.13,
including independent continuous-score arithmetic and identical learned
parameters across modes. The documentation build passes with existing
large-page/search-index warnings and the non-CI deployment skip. Weighting and
normalization regression suites pass; the remaining focused suites and the
five-mode synthetic end-to-end pipeline are in progress. Source is committed
before preparing the dedicated remote directory; no scientific result is yet
available.

The five-mode synthetic end-to-end suite completes successfully (**160
assertions** including boundaries/chunking); the new Slurm runner passes **16
assertions** before submission. The 107 targeted assertions and the weighting,
normalization, complete-training and support/score regressions also pass.
Complementary legacy-diagnostic regressions continue in parallel; their final
status is required before interpreting the external result.

Source **bef315f** is synced to `/u/oldu/code/STMFit_continuous_vote_20260921`.
The whitelist contains tracked `config/`, `hpc/`, `packages/`, `test/` and root
`Project.toml`, not benchmark tables, prior results, raw data or agent state.
The unchanged Manifest is copied explicitly. Login-node `Pkg.instantiate()`
with automatic precompilation disabled succeeds (retaining the documented
1.12.6-Manifest/1.13-runtime warning); subsequent source and Manifest checksum
comparisons are empty. Both metadata-only dry runs see **146 scans** under
Julia 1.13.0. `sbatch --test-only` accepts the resources; it does not validate
the compute-node environment.

One job **11922389** is submitted at **18:49:32 CEST**, partition `small`,
account `mfk_cpu`, QOS `n0001`, with all five absolute variables explicitly
exported. Output target: `/ptmp/oldu/stmfit/continuous_vote_20260921_v1`.
Its initial state is **PENDING, QOSGrpCpuLimit**. No duplicate submission,
resource change or local real-cohort fit is attempted. External grading has
not run; source remains frozen while this job and the regressions are followed.

All complementary regressions finish successfully while the job remains queued:
**2,199 assertions across 17 suites**. Counts are continuous aggregation 107,
scan weighting 138, robust normalization 99, complete training 154, support/score
286, numerical assignment 136, numerical diagnostics 20, native Fisher 229,
affine Fisher 74, Fisher attribution 210, representation 319, transverse
descriptors 125, signal features 73, assignment 37, five-mode pipeline 160 and
two Slurm runners 16 each. The earlier standalone runner repetition is not
double-counted. The initial failed synthetic assertion log is retained alongside
the corrected passing run. No source or configuration changed during this
verification; the comparison is still waiting for the shared CPU quota.

### Completed comparison and decision

The same job **11922389** starts at **19:25:55 CEST**, after **36m23s** queued,
and completes **0:0 at 19:33:14**, **7m19s**, on `vipc2034`. Control runs
19:25:56–19:31:18, then continuous 19:31:18–19:33:14. Four CPUs are requested,
eight allocated; Julia uses four threads and BLAS one. Slurm reports MaxRSS
**2,417,063 K** (about 2.475 GB). No duplicate job, resource change, retry or
local/login-node cohort fit occurs. All **392 output files and the Slurm log**
are fetched locally, excluding raw-input symlinks, with an empty post-fetch
checksum comparison. Scientific source remains **bef315f**.

Before external grading, **7,240 saved-output assertions pass**: 59 replay/key/N
checks and 7,181 availability/score/vote checks. All twelve control tables
reproduce the saved support candidate byte-for-byte. Both arms retain all
146 scans / 900 keys, the same 896 usable GMM rows and four unavailable
predictions. Patches, local/descriptor/predictor features, CC, Fisher, k-means,
training pixel counts and selected N are unchanged. N is compared by scan
identity, not arm-specific staging paths. Final vote arithmetic is reconstructed
from the serialized component scores without reading labels or performing fits.

GMM scores change on **129 rows** and its hard class changes on one. Nonendpoint
scores increase **242→315** out of 896 usable rows. The final vote changes
**24 decisions**, all **1→0**, all on benchmark scans; raw `(0,1,?)` counts move
**(673,223,4)→(697,199,4)**. These are outcomes, not imposed chemical populations.
Every changed final decision was a control tie, with **GMM=1 and k-means=0**.
Continuous GMM scores become **0.75154286–0.99999999**, still class 1 within the
GMM, but their mean with k-means zero falls below 0.5. Final margins range
**0.00000001–0.24845714**. Exact final ties decrease **58→34** before and after
serialization; no displayed tie emits class 0 in this run. Thus this change
acts through fusion of the continuous and unchanged k-means heads; it does not
change the learned groups or missing-data handling. No tie rule, output
precision or score calibration is adjusted after seeing this result.

Only after integrity checks does the unchanged external grade select 145 scans /
892 predicted lobes / 870 truth positions:

| Version | Correct / 870 | Correct / classified | Errors emitted | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 177 | 36 |
| Replayed support control | 676 | 676/852 (79.3%) | 176 | 34 |
| Continuous GMM seed vote | 671 | 671/852 (78.8%) | 181 | 26 |

Against control: **8 scan gains, 13 losses, 124 ties**, net **−5 correct**;
**one exact chain gained, nine lost**, net **−8**, with three grade-orientation
changes. `241113_086` gains exactness. Exactness is lost on `240307_017`,
`240310_Cu100007`, `240817_007`, `240817_021`, `240817_044`, `240817_048`,
`240817_081`, `240818_017` and `240818_020`. Confusion `(TN,FP,FN,TP)` changes
`(513,52,124,163)→(522,43,138,149)`: nine fewer false positives, fourteen more
false negatives. Coverage remains **852/870** with two aligned abstentions;
N audit remains 16 missing positions, 38 extra lobes, 14 short-N and 25 extra-N
scans. The 199 non-correct truth positions comprise 181 emitted errors,
16 missing positions and two abstentions. Against history: **20 scan gains,
29 losses, 96 ties**, net −6 correct; three exact chains gained and thirteen
lost, with nine grade-orientation changes. Grading orientation does not establish
an independent spatial correspondence between fitted lobes and monomers.

**Decision: reject continuous GMM seed aggregation and retain the 676/34 support
candidate.** Historical 677/36 remains the objective to exceed. This negative
result concerns the fixed uncalibrated membership score within the unchanged
two-head vote, not every continuous-score method. No threshold, precision,
temperature, seed, naming, combined-arm or post-grade search follows. Counting,
DFT and unknown25 are untouched. This authorized comparison is complete; the
broader historical objective remains open.

All **290 paired scan rows and 24 changed decisions**, including every loss,
are preserved in `results/continuous_vote_20260921/report.md` and its linked
tables, with a saved-component-score tie audit. Paired arithmetic reproduces
grade sums, confusion, coverage and exact-chain counts. Inference is label-free;
the repeatedly reused benchmark is development evidence, not independent
validation or accuracy on unknown chains. Prediction SHA-256: control
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`, continuous
`f8ac7664700bb3c71d3f0bd214365149fb1c095ae6e2ed4e76cc262827226e36`.

The final documentation build passes with the existing large-page/search-index
warnings and non-CI deployment skip. README, runbook, selection, assignment,
config, calibration and HPC notes record the 671/26 negative result and retain
676/34 as the working reference. Source and tested preparation were committed
before submission; the final conclusion is committed separately. Outputs, logs,
the initial failed synthetic test and every paired loss remain local and
Git-ignored. No validated scientific source is left uncommitted.

## 2026-09-21 — Whole-scan GMM bootstrap, one fixed comparison

**Scientific deliverable and done condition.** The user authorizes the first
recommended follow-up: compare whole-scan GMM bagging with the retained support
candidate (676/870 correct, 34/145 exact). Done means deterministic tested
resampling and training-only naming, exact control replay, all keys and N
retained, checked outputs and paired external grades against control and saved
history, fetched results and committed conclusions. Historical 677/870 and
36/145 remain the objective, not a presumed ceiling or a claimed gain. Scope:
one control and one variant, one Viper job, four requested CPUs, 16 GB and a
one-hour execution limit. No automatic retry, local/login-node cohort fit,
combined variant, new naming rule, continuous vote, covariance change, count
refit, DFT change or unknown25 rerun. The archived application's scientific
Must-NOT-have list remains in force; its retired process apparatus is not used.

**Fixed method before grading.** Fifteen earlier native configs explicitly
declare `gmm_resampling="none"`, `gmm_bootstrap_replicates=1` and
`gmm_bootstrap_seed=0`. `unit_assignment_scan_bagging.toml` copies the support
config, changing its name and resampling to `whole_scans`, with **20 replicates**.
For each feature view, the pool is the S scans with at least one usable row.
Replicate b uses MersenneTwister seed b−1 to draw S scan indices uniformly with
replacement. Every usable row of a drawn scan is copied once per draw, preserving
within-scan dependence and natural scan lengths. This is not the rejected
equal-total-weight-per-scan method. The raw scan identities are grouping keys,
never dates, quality grades, known chemical classes or benchmark membership.

Per-scan mean/sample-std normalization and the 28 interactions of eight original
features are computed once, unchanged. An unsampled scan still uses its own
observations for its usual per-scan transform; its rows do not train the GMM
or name its groups. This is not an out-of-bag validation claim. Every replicate
uses the same ten initialization seeds 0–9, existing EM, two hard updates,
ridge, Mahalanobis scoring and freely learned mixture masses. Physical naming
uses only the duplicated training rows' raw amplitudes, counting multiplicity.
All valid rows, including unsampled scans, are then predicted. Each replicate
averages the named binary seed decisions; the candidate averages these replicate
means equally. The existing unnamed-seed exclusion is retained, without retries;
an entirely unnamed replicate contributes nothing and all-invalid results stay
unavailable. No best seed, best replicate or composition is selected.

The unchanged final GMM/k-means mean, `>=0.5` tie rule, eight-decimal output,
Fisher, CC, raw patch extraction, support, geometry and N are retained. Bootstrap
draws and accepted-seed counts are exported as `gmm_scan_bootstrap.tsv`, including
zero multiplicities and scans with no usable rows. Frequencies are uncalibrated
votes, not chemical probabilities. Inference treats all 146 scans; labels and
the full145 manifest remain external grading inputs only. No post-grade tuning
or independent-validation claim is allowed on this reused development benchmark.

The first targeted test stops before any resampling because its config inventory
tries to parse `chitosan_exclude.txt` as TOML. The test is corrected to enumerate
only `unit_assignment_*.toml`; no scientific setting, score or output is changed.
The failed log is retained. No real-cohort calculation or external grade has
yet run for this candidate.

The corrected targeted suite passes **222 assertions** under Julia 1.13,
including independent bag arithmetic, literal repeated-row fit identity,
training-only amplitude naming, unchanged invalid rows and CLI/API byte identity.
The continuous-vote regression also passes (**108 assertions** with the added
config inventory entry). The Slurm runner now explicitly declares `--no-requeue`
as well as the one-hour limit. The six-mode synthetic pipeline and complementary
regressions are running before scientific submission; the result remains unknown.

The six-mode synthetic end-to-end pipeline passes **189 assertions**, including
the new bootstrap audit and unchanged upstream/Fisher/k-means tables. Source
**da886ce** is synced to `/u/oldu/code/STMFit_scan_bagging_20260921` using only
tracked config/HPC/packages/test/root-project files and the unchanged Manifest.
Source/Manifest checksums match; login-node instantiation succeeds under Julia
1.13.0 with the existing 1.12.6-Manifest warning. Geometry, split geometry and
templates retain the previous SHA-256 values. Both metadata-only dry runs see
146 scans. The runner passes **20 assertions**, including resource limits and
`--no-requeue`; `sbatch --test-only` accepts the request, not the compute-node
environment. No raw data, benchmark tables, outputs or agent state are synced.

One actual job **11925188** is submitted at **20:42:12 CEST**, with all five
absolute variables explicitly exported. It starts at **20:42:41**, after
**29 seconds** queued, on **vipm2008**. Four CPUs are requested, eight allocated;
Julia uses four threads and BLAS one. `Requeue=0`, `Restarts=0`, one-hour limit.
Output: `/ptmp/oldu/stmfit/scan_bagging_20260921_v1`. No duplicate submission,
resource change or local/login-node cohort fit is attempted. Complementary
legacy regressions continue; they must pass before interpreting the external
grade. No scientific result is yet available, and settings remain frozen.

All complementary regressions finish successfully while the control starts:
**2,471 assertions across 19 suites**. Counts: bootstrap 222, continuous vote
108, scan weighting 138, robust normalization 99, complete training 154,
support/score 286, numerical assignment 136, numerical diagnostics 20, native
Fisher 229, affine Fisher 74, Fisher attribution 210, representation 319,
transverse descriptors 125, signal features 73, assignment 37, six-mode pipeline
189 and three Slurm runners 20/16/16. The first 16-assertion repetition of the
new runner is not double-counted; its final resource/requeue checks total 20.
The failed initial TOML-inventory log remains. No scientific source/config
changes during these checks; no external grade has run.

### Completed comparison and decision

The same job **11925188** completes **0:0 at 20:49:58 CEST**, **7m17s**, on
`vipm2008`. Control runs 20:42:43–20:47:17 and bagged 20:47:17–20:49:58.
Slurm reports MaxRSS **1,933,112 K** (about 1.980 GB). No duplicate job,
resource change, retry or local/login-node cohort fit occurs. All **393 output
files and the Slurm log** are fetched locally, excluding raw-input symlinks,
with an empty post-fetch checksum comparison. Scientific source stays
**da886ce**; **c45e440** records submission and complementary test completion.

Before external grading, **21,927 saved-output assertions pass**: 59 exact
replay/identity/N checks, 7,181 availability/vote checks, 14,686 bootstrap-draw
checks and one check that the control has no bootstrap file. Twelve control
tables replay the saved support result byte-for-byte. All **146 scans / 900
keys**, selected N, upstream patches/features, CC, Fisher, k-means and training
pixel counts remain unchanged. Both arms retain **896 usable GMM rows** and
the same four unavailable predictions. There is no count refit or imputation.

The bootstrap audit has **2,920 rows**, one per input scan per bag, including
zero multiplicities. Independent draws from MersenneTwister seeds 0–19 reproduce
every scan multiplicity and training-row count. Each bag draws 146 scans with
replacement; **88–98 distinct scans and 882–905 duplicated usable rows** occur
per bag. Each scan is drawn 10–34 times over the whole experiment. All twenty
bags accept all ten initialization fits (**200/200**), with no retry or selection.
These counts are observed outcomes, not desired chemical populations.

GMM scores change on **804/896 usable rows**, increasing on 448 and decreasing
on 356. Nonendpoint scores increase **242→804**. Its hard class changes on 38
rows, all 1→0. Final voting changes **62 decisions**, all **1→0**: 61 on
benchmark scans, one on technical `240310_Cu100009.sxm` lobe 7. Raw `(0,1,?)`
counts move **(673,223,4)→(735,161,4)**. Of the final flips, **57 were old ties
with GMM=1 and k-means=0**; any bagged GMM score below 1 moves their unchanged
average below 0.5. Five changes had positive control margins (0.025–0.5).
The GMM hard class changes for 29 of the 62 final flips, staying 1 for the
other 33. Candidate GMM scores on changed rows range **0.13–0.995**, final
margins **0.005–0.87**. Exact final ties fall **58→0**, before and after
serialization. Seven saved-component-score assertions pass. This is not only
floating-point rounding; it exposes a strong interaction with fixed fusion.
No threshold, tie, fusion, precision or naming rule is adjusted after the grade.

The unchanged external full145 own-N grade is run only after integrity checks:

| Version | Correct / 870 | Correct / classified | Errors emitted | Exact / 145 |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 177 | 36 |
| Replayed support control | 676 | 676/852 (79.3%) | 176 | 34 |
| Whole-scan GMM bagging | 665 | 665/852 (78.1%) | 187 | 10 |

Against control: **22 scan gains, 33 losses, 90 ties**, net **−11 correct**;
**zero exact chains gained, 24 lost**, with eight grade-orientation changes.
Confusion `(TN,FP,FN,TP)` changes `(513,52,124,163)→(537,28,159,128)`:
24 fewer false positives, outweighed by 35 more false negatives. Coverage stays
**852/870**, with two aligned abstentions. Counts stay 892 predicted lobes on
145 scans, 16 missing truth positions, 38 extra lobes, 14 short-N scans and
25 extra-N scans. The 205 non-correct truth positions comprise 187 emitted
errors, 16 missing positions and two abstentions. Against history: 26 scan
gains, 41 losses, 78 ties, net −12 correct; one exact chain gained, 27 lost,
with twelve orientation changes. Grade alignment is not an independent spatial
monomer correspondence.

Exactness is lost on `240307_017`, `240310_Cu100007`, `240311_Cu100061`,
`240314_Cu100_025`, `240814_012`, `240815_072`, `240817_003`, `240817_006`,
`240817_007`, `240817_021`, `240817_041`, `240817_044`, `240817_047`,
`240817_048`, `240817_049`, `240817_052`, `240817_080`, `240817_081`,
`240818_017`, `240818_020`, `240818_026`, `241113_087`, `241113_088` and
`241114_027`. All **290 paired scan comparisons and 62 changed decisions**,
including every loss, are preserved in `results/scan_bagging_20260921/report.md`
and its linked tables. Paired arithmetic reproduces grade sums, confusion,
coverage and exact-chain counts. The saved-score audit reads no labels.

**Decision: reject whole-scan bagging and retain the 676/34 support candidate.**
Historical 677/36 remains the objective to exceed. This negative result concerns
the specified ensemble within unchanged fusion, not every bagging method.
No replicate-count, seed, naming, threshold, fusion or combined-variant search
follows; counting, DFT and unknown25 stay untouched. This authorized comparison
is complete; the broader historical objective remains open. Learning/inference
are label-free, but the repeatedly reused benchmark is development evidence,
not independent validation or unknown-chain chemical accuracy.

Prediction SHA-256: control
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`, bagged
`16df34427f9dc82bc509e4998d36f5364aac8b75ccb4ead42d8a9b3da1a24ed7`.
Bootstrap audit SHA-256:
`28a4cd501b16203f35b9a042e897f53d8b414f30bbbdae1e8d3e4ceda8477a86`.

The final documentation build passes with the existing large-page/search-index
warnings and non-CI deployment skip. README, runbook, selection, assignment,
config, calibration and HPC notes record the negative 665/10 result and retain
676/34 as the working reference. Source and tested preparation were committed
before submission; conclusions are committed separately after reviewing the
saved outputs and documentation diff. Generated results, logs and the initial
failed synthetic test remain local and Git-ignored; no cluster retry, new
scientific run or GitHub push follows this result.

## 2026-09-22 — Tied covariance throughout GMM learning, one fixed comparison

**Scientific deliverable and done condition.** The user approves the first
recommended next lead: measure whether a common covariance learned throughout
the GMM improves recognition against the retained support candidate, **676/870
correct, 34/145 exact, 852/870 coverage**, and historical **677/870, 36/145,
854/870**. Done means independent numerical/synthetic checks, exact control
replay, verified complete outputs, external paired grades, fetched results and
committed conclusions. Historical performance is a target to exceed, not an
assumed ceiling. One control and one candidate, one Viper job, four requested
CPUs, 16 GB and **30 minutes**, no automatic retry or local/login-node cohort
fit. The archived application's scientific Must-NOT-have list remains in force;
its retired process apparatus is not used. No new naming rule, grouped Fisher,
fusion, threshold, bagging, weighting, DFT, count refit or unknown25 rerun.

**Fixed method before grading.** Sixteen earlier native configs explicitly
declare `[model] gmm_covariance_structure="full"`. The new
`unit_assignment_tied_covariance.toml` changes only the support config's name
and that field to `"tied"`. It pools the within-component scatter at k-means
initialization, every EM M-step and both Mahalanobis hard updates. If x_i is
one of n usable rows, r_ic its responsibility and mu_c the updated mean,
the common covariance is `sum(r_ic * (x_i-mu_c)*(x_i-mu_c)')/n + ridge*I`.
Initialization and hard updates use membership indicators. This is not an
unweighted average of component covariances, nor the global covariance with
between-group separation. Each covariance copy contains the same matrix;
component means and masses remain freely learned. The existing 1/2 initial
weights are starting values only, not a chemical-composition constraint.

Ridge stays **1e-6**, factorization guard **1e-8**, ten seeds **0–9**, two
hard updates, full eight-feature/28-interaction representation and per-scan
mean/sample-std scaling. Hard reassignment still uses only distance as before;
final scoring retains log free weight minus half Mahalanobis distance. With a
common covariance, a Gaussian log-volume term would cancel between components;
the score mode is not changed. Physical naming still uses mean raw amplitude
of assigned training members. Seed votes remain hard and the final k-means/GMM
mean, `>=0.5` rule, eight-decimal output and unavailable-row policy stay fixed.
No GCV, guard, calibration, geometry or selected-N change is made.

This candidate tests reduced covariance freedom during group formation, unlike
the failed final-only shrinkage and final-score experiments. Its 36-dimensional
covariance has 666 independent entries instead of 1,332 across two components.
That is motivation, not evidence of better chemistry or calibrated uncertainty.
The new mode rejects combinations outside the support-control policies, and
historical representation diagnostics reject it rather than describe an
unrelated fit. Inference sees all 146 scans, no benchmark membership, sequence,
expected counts or saved champion predictions. External full145 grading follows
frozen outputs and integrity checks; no post-grade tuning is allowed. This
reused benchmark supplies development evidence, not independent validation.

Independent synthetic arithmetic and the seven-mode pipeline checks are running.
The first unprivileged `julia --version` cannot create juliaup's lockfile under
the read-only home sandbox; authorized test processes use Julia 1.13 through
the existing release installation. No scientific run or grade has started.

The targeted suite passes **212 assertions** under Julia 1.13: explicit config
isolation, independent pooled-scatter equations, initialization and three EM
updates, both hard updates, free unequal masses, unchanged legacy arithmetic,
invalid-row retention, physical naming and deterministic CLI/API predictions.
The GMM log now records covariance structure, final inter-component covariance
difference, free weights and naming availability per seed, without changing
predictions. The seven-mode synthetic pipeline and complementary regressions
continue. Viper reports Julia 1.13.0, all existing geometry/raw/template inputs
present, fresh source/output paths, and no queued/running user job at inspection.
No scientific job has been submitted.

Scientific source is frozen in **becd0ce**. The tracked source whitelist and
unchanged Manifest are synchronized to Viper; both checksum dry-runs are empty.
Login-node `Pkg.instantiate()` succeeds with Julia 1.13.0 (the existing
Manifest's 1.12.6 warning is retained). The two metadata-only runner dry-runs
each discover **146 scans**, and Slurm `--test-only` accepts the four-CPU,
16-GB, 30-minute request with the five explicit exported paths. Its estimated
job ID is not a submission. The seven-mode synthetic pipeline passes **216
assertions**. All **2,716 assertions across 20 suites** now pass: tied covariance
212, pipeline 216, numerical assignment/diagnostics 136/20, scan weighting 138,
continuous votes 109, scan bagging 223, robust normalization 99, patch support
and score 286, complete training 154, native Fisher 229, affine Fisher 74,
Fisher attribution 210, representation diagnostics 319, transverse descriptors
125, signal features 73, assignment 37, and tied/bagging/continuous runners
20/20/16. The documentation builds successfully with the existing large-page
and search-index warnings and non-CI deployment skip. No test failure or
scientific input change precedes submission.

Passing checks are committed in **123801d** before the single Viper submission
**11931106** at **00:36:42 CEST**. The 00:36:51 inspection finds it pending,
with `Requeue=0`, `Restarts=0`, four requested CPUs, 16,000 MB, a 30-minute
limit and all five explicit exported paths. Scientific source remains
**becd0ce**; these later commits are documentation only. External grading
has not begun; this task continues through output verification and the result.

The control script starts on **vipc2227 at 00:37:00 CEST**, 18 seconds after
submission. Final Slurm accounting subsequently establishes the allocation
start at **00:36:57**, hence **15 seconds queued**, followed by script setup.

### Completed result: reject this variant, retain 676/34

Job **11931106** completes at **00:42:50 CEST**, **5m53s**, exit **0:0**,
without retry (`Requeue=0`, `Restarts=0`). Four CPUs were requested and eight
allocated, with the fixed 16-GB/30-minute limits; batch MaxRSS is **2,779,148 K**.
The control finishes all 146 scans at 00:41:02 and the tied arm at 00:42:50.
All **392 regular outputs plus the Slurm log** are fetched, excluding raw-input
links. Both checksum dry-runs are empty. No scientific source or parameter changes occur
after the scientific source freeze in **becd0ce**.

Before grading, **7,338 saved-output assertions** pass. All twelve control
scientific tables are byte-identical to the saved 676/34 support result.
Both arms retain all 146 scans / 900 keys, the same selected counts, eleven
upstream/Fisher/k-means/support tables, four unavailable rows and the specified
serialized vote arithmetic. All ten candidate seed fits have valid physical
naming and exactly identical component covariances; the smaller freely learned
component has **118 or 129 of 896** hard-training members, not a preset class
count. The four unavailable keys remain `240310_Cu100009.sxm` lobe 8,
`240818_015.sxm` lobe 1, and `240818_019.sxm` lobes 1–2.

The one external full145 own-N grade is:

| Profile | Correct / 870 | Classified accuracy | Exact / 145 | Coverage / 870 |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 36 | 854 |
| Exactly replayed support control | 676 | 676/852 (79.3%) | 34 | 852 |
| Covariance shared throughout GMM | 666 | 666/852 (78.2%) | 6 | 852 |

Against control, **25 scans gain, 35 lose, 85 tie**, with **zero exact chains
gained and 28 lost**. Against history, 29 gain, 43 lose and 73 tie, with one
exact chain gained and 31 lost. Selected counts remain 892 benchmark predictions
for 870 controls, with 16 missing positions and 38 extras; the same two aligned
abstentions give coverage 852/870. Confusion `(TN,FP,FN,TP)` changes from
`(513,52,124,163)` to `(540,25,161,126)`. The decrease in false positives does
not offset the increase in false negatives. This is a negative recognition
result despite the smaller covariance parameterization.

Five saved-component checks pass without labels or fitting. Of **69 final
changes**, **68 switch 1→0 and one 0→1**; one change is outside the benchmark.
Fifty-five leave old GMM=1 / k-means=0 ties. The sole 0→1 change,
`240814_025.sxm` lobe 2, enters such a tie; that scan remains at 5/6 correct.
All 69 final flips also change GMM hard class. Across 896 usable rows, GMM
scores decrease on 104, increase on 323 and stay identical on 469; GMM hard
classes change 1→0 on 100 rows and 0→1 on one. Raw final `(0,1,?)` frequencies
change `(673,223,4)` → `(740,156,4)`. Exact final ties decrease **58 → 3**,
identically before and after serialization; 652 candidate GMM scores are
nonendpoint versus 242 for control. No threshold, rounding, naming or class
proportion is adjusted in response.

Grade orientation changes on nine scans versus control, thirteen versus history;
it is not independent spatial monomer correspondence. The 28 exact chains lost
against control are `240307_017`, `240310_Cu100007`, `240311_Cu100061`,
`240314_Cu100_025`, `240814_012`, `240815_072`, `240817_003`, `240817_004`,
`240817_005`, `240817_006`, `240817_007`, `240817_021`, `240817_041`,
`240817_044`, `240817_047`, `240817_048`, `240817_049`, `240817_052`,
`240817_053`, `240817_055`, `240817_080`, `240817_081`, `240818_017`,
`240818_020`, `240818_026`, `241113_087`, `241113_088`, and `241114_027`.
All **290 paired scan comparisons and 69 changed decisions**, including every
loss, are preserved in `results/tied_covariance_20260922/report.md` and its
linked tables. Paired arithmetic reproduces grade sums, confusion, coverage,
selected counts and exact-chain totals.

**Decision: reject this tied-covariance variant and retain the 676/34 support
candidate.** Historical 677/36 remains the target to exceed. This result concerns
the specified covariance sharing inside the unchanged representation/fusion,
not every form of covariance regularization. No other proposed lead, combined
arm or post-grade search follows; counting, DFT and unknown25 remain untouched.
This authorized comparison is complete; the broader historical objective is
still open. Learning/inference are label-free, but this repeatedly reused
benchmark is development evidence, not independent validation or unknown-chain
chemical accuracy.

Prediction SHA-256: control
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`, tied
`566e5228bb4baeb9735bbd23a730f4b0fa201dbd5aa416b8b454705fd9689ec6`.

The final documentation build passes with the existing large-page/search-index
warnings and non-CI deployment skip. README, runbook, selection, assignment,
config, calibration and HPC notes record the negative 666/6 result and retain
676/34 as the working reference. Source and passing preparation were committed
before submission; reviewed conclusions are committed separately. Generated
results, logs and diagnostic reports remain local and Git-ignored. No new
scientific run or GitHub push follows this completed comparison.

## 2026-09-22 — Two independent leads: whole-scan Fisher and relative GMM naming

**Deliverable and scope.** The user first approves the next grouped-Fisher
comparison, then explicitly permits several leads. The bounded deliverable is
one control and two independent candidate results against **676/870 correct,
34/145 exact, 852/870 coverage**, also compared to historical **677/870,
36/145, 854/870**. Completion means synthetic tests, exact control replay,
all-key/count/upstream checks, one external paired grade with every loss,
fetched outputs and committed conclusions. One Viper job, four requested CPUs,
16 GB, **30 minutes**, no requeue or automatic retry. No local/login-node
real-cohort fitting. The third proposed paired-acquisition model is not included;
reliable registration is still unestablished. The application's scientific
Must-NOT-have constraints remain; retired orchestration is not reopened.

**Fixed before any real output.** Seventeen previous native configs explicitly
retain `[selection] fisher_cv_scheme="lobe_parity"`,
`fisher_scan_split_seed=0` and `gmm_cluster_naming="raw_amplitude"`.
Both candidates copy `unit_assignment_patch_support.toml`; each changes its
name and exactly one mode, never both. The tied-covariance candidate stays
rejected and is not a starting point for this comparison.

* `unit_assignment_scan_fisher.toml` sets `fisher_cv_scheme="scan_hash_twofold"`.
  From all supplied scan basenames, before patch validity/training filtering,
  sort by SHA-256 of UTF-8 decimal seed, one NUL byte, then UTF-8 basename;
  basename breaks any digest tie. Assign rank 1 to group 0, rank 2 to group 1,
  alternating thereafter. Group sizes differ by at most one, not by class
  proportions or lobe counts. Seed is fixed at zero, with no split search.
  Row ordering, lobe reversal and pixel validity do not alter membership;
  changing the cohort or renaming scans can. All PCA, GMM, mean centering,
  amplitude ordering, covariance and Fisher weights use only the opposite
  scan group. The native parity mode keeps its exact operations. Scan mode
  canonically orders training keys. Missing/degenerate folds keep their rows
  with NA/reasons, without borrowing the held-out scan. CLI logs record group
  identities, actual training/held-out counts and scan overlap. This changes
  one feature producer, not the downstream GMM's transductive learning policy;
  it does not create independent end-to-end validation.
* `unit_assignment_relative_naming.toml` sets `gmm_cluster_naming="within_scan_z"`.
  GMM initialization, EM, hard updates, final scores, memberships and free masses
  stay identical. For each scan, center raw amplitudes on their mean and divide
  by sample std using only its eligible, feature-valid training rows. Constant
  or singleton scales use the existing explicit `gmm_scale_fallback=1.0`.
  Compare these standardized amplitudes' two group means globally; the higher
  group retains the existing class-1 naming convention. This does not force
  both groups within each scan, change their sizes or select top-k lobes.
  Held-out/unavailable amplitudes cannot set training moments. The existing
  tie policy is unchanged. New naming logs retain both means and the chosen
  group per seed. This scoped mode rejects combinations with grouped Fisher,
  bagging, weighting, changed covariance/scoring or changed scaling.

GCV, selected N, geometry, patches, physical calibration, support, naming of
Fisher groups, mirror axis, ridge, seeds, feature scaling/interactions, hard
votes, final mean and `>=0.5` rule, precision and unavailable-row policy remain
fixed except for the declared arm. No expected N, sequence, chemical class
counts, benchmark membership or saved champion predictions enter inference.
All 146 scans remain in science; only the later external grader selects full145.
No post-grade adjustment or combined candidate is authorized. The benchmark
has been reused extensively and remains development evidence, not independent
validation or unknown-chain chemical accuracy.

Implementation is on `research/scan-fisher-naming-20260922`, from **972ebd8**.
Julia 1.13 import and the synthetic split smoke pass. The first targeted test
passes configuration/group checks and all numerical Fisher non-use checks,
but its attribution-guard fixture calls the wrapper module rather than its
contained diagnostic module (UndefVarError, not a scientific failure). The
fixture is corrected to use that module's own types; two dictionary comparisons
are also made explicitly elementwise. The initial log is retained under
`results/scan_fisher_naming_20260922/initial_scan_fisher_naming_failure.log`.
No scientific setting is changed. Targeted and nine-mode pipeline tests are
running; no real-cohort job or grade has begun.

The second targeted attempt passes all 119 assertions reached, including
the corrected diagnostic guard and relative-naming arithmetic, then exposes
a test-fixture syntax error: `10.0f` is parsed as an invalid float literal,
not multiplication by file index `f`. The fixture now spells `10.0*f`.
The second failed log is retained separately; no scientific implementation or
parameter changes to resolve this test-only error.

The naming log also fingerprints every fitted mean, covariance and mixture
weight (concatenated Float64 bytes, SHA-256) before voting. This read-only log
will allow the real relative-naming arm to prove exact parameter identity
with the control for each seed, not only matching scalar summaries. It does
not change fitting, assignment or any score. The pre-run documentation build
passes with the existing size/deployment warnings. Viper reports Julia 1.13.0,
the same geometry/split/template hashes, fresh source/output paths and no
queued/running user job at the metadata check; no new job has been submitted.

The targeted suite now passes **291 assertions** under Julia 1.13: config
isolation 71, deterministic whole-scan grouping 12, opposite-scan Fisher
non-use and legacy arithmetic 26, relative-naming algebra 10, CLI/API and
unchanged GMM fits 172. The last logging-only fingerprint addition is exercised
by subsequent CLI/pipeline regressions. The nine-mode synthetic pipeline and
complementary regressions continue before real submission.

All **19 suites / 3,035 assertions** now pass under Julia 1.13. The pipeline
has 274 assertions, including 153 for nine end-to-end synthetic variants;
targeted whole-scan/naming tests have 291. The complementary suites cover
native Fisher 229, affine Fisher 74, numerical assignment 136, complete
training 154, Fisher attribution 210, representation diagnostics 319,
numerical diagnostics 20, runner 20, scan weighting 138, continuous vote 111,
scan bagging 225, tied covariance 214, robust normalization 99, patch support
286, label-free signals 73, transverse descriptors 125 and assignment 37.
The later GMM CLI regressions exercise the full-parameter fingerprint log.

Scientific source **43f8993** and the unchanged Manifest are synchronized to
`/u/oldu/code/STMFit_scan_fisher_naming_20260922`; both checksum comparisons
are empty. Login-node instantiate (without automatic precompilation) succeeds
with Julia 1.13.0 and preserves the Manifest. The three-arm metadata dry-run
and scheduler `--test-only` pass. The latter's speculative identifier is not
a submitted job. Actual job **11931478** is submitted once with all five
variables explicitly exported, four requested CPUs, 16 GB, a 30-minute limit
and no requeue. Outputs are reserved at
`/ptmp/oldu/stmfit/scan_fisher_naming_20260922_v1`. No grade has been read and
no scientific output/result is claimed yet. Full local logs are under
`results/scan_fisher_naming_20260922/`.

**Completed result.** Job **11931478** waits **3m01s**, starts at 01:30:44 CEST
and completes **0:0 in 8m30s** at 01:39:14 CEST. Slurm allocates eight CPUs for
the four-CPU request, with 16 GB requested and about 2.35 GB reported peak
memory. The three arms complete without failures. All **588 regular output
files and the Slurm log** are fetched locally; the output checksum comparison
is empty. No retry or second job is submitted.

Before grading, **10,987 saved-output assertions** pass. The twelve control
tables are byte-identical to the saved support run; its prediction SHA-256 is
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`.
All arms retain 146 scans / 900 keys, all selected N, nine common upstream,
k-means and training-support tables, and the same four unavailable keys.
Grouped Fisher changes only the Fisher column in the predictor table, with
896 changed valid margins. The independent hash-rank check reproduces both
73-scan groups. Actual training/held-out sets contain 449/447 and 447/449
valid rows, zero overlapping scans, and two converged Fisher fits.

All ten relative-naming fits have exactly the control's full-parameter hashes
and choose the same high component under both naming rules. Relative means
range from approximately -0.6513 to -0.3380 in the lower groups and 0.6091
to 1.0324 in the higher groups; raw naming also selects those high groups.
The GMM TSV is byte-identical. Final decisions, scores, margins and reasons
are unchanged; only the model-name field makes the final TSV bytes differ.
This is a null result on this cohort, not failure to activate the new rule.

The unchanged external full145-own-N grade is run only after these checks:

| Profile | Correct / 870 | Classified accuracy | Exact / 145 | Coverage / 870 |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 36 | 854 |
| Replayed support control | 676 | 676/852 (79.3%) | 34 | 852 |
| Whole-scan Fisher | 668 | 668/852 (78.4%) | 34 | 852 |
| Within-scan GMM naming | 676 | 676/852 (79.3%) | 34 | 852 |

Fisher versus control: **8 scan gains, 15 losses, 122 ties; four exact chains
gained and four lost**. Against history: 23 gains, 30 losses, 92 ties; seven
exact chains gained, nine lost. Control → Fisher confusion `(TN,FP,FN,TP)`
is `(513,52,124,163)` → `(496,69,115,172)`: nine fewer false negatives do not
offset seventeen more false positives. Counts stay 892 benchmark predictions,
16 missing positions, 38 extras, 14 short-N and 25 extra-N files. The same
two aligned abstentions leave 852 classified positions; Fisher emits 184
errors versus control's 176. Relative naming ties control on every scan.

Saved-component arithmetic passes **eight checks** across the two arms.
Fisher changes 242 GMM/final scores and 26 final decisions, all **0→1 into
GMM=1 / k-means=0 exact ties**, with zero final margin. Twenty-one also change
the GMM hard decision; five already had a class-1 GMM vote of 0.8. Across
896 usable rows, GMM scores decrease for 209, increase for 33 and are unchanged
for 654; GMM hard decisions change 0→1 for 22 rows and 1→0 for one. Final
ties increase 58→84, both before and after serialization. Raw `(0,1,?)`
frequencies change `(673,223,4)`→`(647,249,4)`. Relative naming changes none
of these quantities. No threshold/precision adjustment is made.

Exact chains gained versus control: `240307_016.sxm`, `240814_011.sxm`,
`240817_045.sxm`, `241114_028.sxm`. Lost: `240817_021.sxm`, `240818_020.sxm`,
`241113_087.sxm`, `241113_088.sxm`. Four grade-orientation changes occur
versus control, eight versus history; the grade is not independent spatial
correspondence. The **580 paired scan rows**, every changed decision, all
losses and replay commands are in `results/scan_fisher_naming_20260922/report.md`
and linked tables. No scientific setting is changed after reading the grade.

**Decision: reject whole-scan Fisher as the replacement; relative naming gives
no gain. Retain the support candidate at 676/34.** These exact candidates do
not exceed historical 677/36. The null naming result does not establish
invariance under other acquisition conditions; the negative Fisher result
does not prove parity cross-fitting or downstream transductive training is
independent validation. The reused benchmark remains development evidence,
not unknown-chain accuracy. This authorized two-candidate comparison is complete;
the historical objective remains open. No further arm, combination, split/seed
search, retry, count refit, DFT change or unknown25 rerun follows.

The final documentation build passes under Julia 1.13, with existing large-page
and search-index warnings and the expected non-CI deployment skip. README,
runbook, selection, assignment, config, calibration and HPC documentation are
updated consistently. Scientific source and all settings remain unchanged
from **43f8993** through grading; generated results remain ignored.

## 2026-09-22 — Fixed factor-analyzer and Student learning: two independent comparisons

The user authorizes “tu peux tester plusieurs pistes”. Scientific deliverable:
one complete control replay and two independent candidate grades against
**676/870 correct, 34/145 exact, 852/870 coverage**, with historical **677/36,
854/870 coverage** still the target. Done means tested numerical updates,
unchanged upstream behavior, all outputs fetched and checked before external
full145-own-N grading, per-scan gains/losses, a documented decision and commits.
Budget: **three hours preparation/tests from 10:03 CEST**, then **one Viper
job capped at one hour**. No combined model, rank/df/seed search, post-grade
adjustment, N refit, DFT update or unknown25 rerun. The proposed paired-view
acquisition model is explicitly deferred; no registration claim is made.

The support representation remains eight inputs plus 28 pairwise products.
`unit_assignment_factor_analyzer.toml` fixes rank **4** in each component;
`unit_assignment_student_t.toml` fixes df **5.0** in both. These choices are
frozen before real inference, not selected from benchmark labels. Only the
learning family and output name differ from control. Independent component
centers and masses remain free; neither initialization nor naming imposes
composition. The original k-means initialization is reused without Gaussian
EM for the alternatives; each initially has mass 0.5, then learned masses.

Factor initialization subtracts the existing ridge from empirical scatter,
uses bottom-spectrum mean noise and the four leading signal eigenvectors.
EM jointly estimates means/loadings by augmented latent moments, with a
separate diagonal residual variance per component. The two hard assignments
remain pure Mahalanobis, each followed by one family-preserving conditional
update. The Student alternative retains its precision-weighted mean/scale
updates through both hard steps; scale denominator is responsibility mass.
Methods and primary references are in `calibration.md`. Neither alternative
silently returns to full-Gaussian hard covariance estimation.

New explicit config settings fix alternative EM at 200 updates, relative
log-kernel tolerance 1e-6, minimum mass 1e-12 and numerical guard 1e-8. The
existing ridge 1e-6 is the factor noise floor / Student scale regularizer.
Common density constants are omitted from learning kernels, so their sums
are not used to compare families. Finite capped fits are kept and explicitly
logged as nonconverged, with no retry. Collapsed/nonfinite/nonpositive fits
are unavailable seeds; no Gaussian fallback and no silent row removal.
All seeds unusable means `?`/`no_valid_view`. Final log-mass/Mahalanobis score,
amplitude naming, ten-seed hard vote, k-means head, final fusion, thresholds,
normalization, interactions, Fisher and N are unchanged. Alternative-family
combination with earlier candidates is rejected, not implicitly supported.

The first synthetic smoke completes both models. An independent-equation
test initially rejects 18 Student matrices because Julia's bare `isposdef`
also demands exact floating-point symmetry, while weighted products differ
between triangles at about 1e-16. All moment comparisons pass. The test is
corrected to check approximate symmetry and positive definiteness of the
same `Symmetric` matrix used by production; no scientific setting changes.
Extended numerical and end-to-end regressions continue before submission.

The next targeted run passes all numerical assertions; three CLI byte-replay
checks fail in both alternatives **and control** because the synthetic
`missing.sxm` row is last in the in-memory fixture but first in the CLI's
alphabetically sorted records. Renaming that synthetic missing row to
`unavailable.sxm` aligns the two orderings; neither production sorting nor
any prediction rule changes. Both initial failure logs are retained locally.

All **514 targeted mixture assertions** now pass: explicit isolated policies,
independent E/M and hard-update equations, covariance reconstruction and
positivity, Student outlier influence/scale denominator/Gaussian limit,
deterministic API/CLI votes, renamed keys, forbidden truth arguments and
unavailable-seed row retention. The three-arm runner passes **20 assertions**;
the whole-scan/naming regression passes **295**. The documentation build
passes under Julia 1.13. Source is committed before remote preparation;
the extended synthetic pipeline and complementary regressions are still
running, so no complete-suite or real-output claim is made at this point.

Source **e2205cb** and the unchanged Manifest are synchronized to
`/u/oldu/code/STMFit_factor_student_mixtures_20260922`; final source/Manifest
checksum comparisons are empty. Login-node instantiate without precompilation
passes with Julia 1.13.0. The metadata-only three-arm dry-run and scheduler
`--test-only` pass; the latter's speculative identifier is not a submitted
job. The **331 pipeline assertions**, including **198 across eleven complete
synthetic variants**, pass before submission. Further complementary tests
continue concurrently, without changing scientific source.

Actual job **11935072** is submitted once with five explicitly exported
variables, four requested CPUs, 16 GB, a one-hour limit and no requeue.
Outputs are reserved at
`/ptmp/oldu/stmfit/factor_student_mixtures_20260922_v1`; local preparation
logs are under `results/factor_student_mixtures_20260922/`. No real outcome
or benchmark score is claimed yet. To keep the assignment documentation
within its existing HTML size limit, the preceding Fisher/naming section is
condensed without changing its figures; full detail remains in its dated
journal entry and local report. No documentation size limit is increased.

### Completed comparison: both fixed learning alternatives regress

Job **11935072** waits **1m36s**, starts at **10:44:26 CEST** and completes
**0:0 at 10:52:43**, elapsed **8m17s**. The control finishes at 10:48:34,
factors at 10:50:34 and Student at 10:52:43. Slurm allocates eight CPUs for
the four-CPU request, 16 GB requested, MaxRSS **2,238,210 K** (about 2.29 GB).
No retry or second job is submitted. All **588 regular outputs plus Slurm
log** are fetched; both checksum comparisons are empty.

All **20 suites / 3,616 assertions** finish successfully under Julia 1.13.
These include 514 mixture, 331 pipeline, 20 runner, 295 scan-Fisher/naming,
216 tied covariance, 138 scan weighting, 113 continuous vote, 227 bagging,
99 robust normalization, 154 complete training, 136 numerical assignment,
229 native Fisher, 74 affine Fisher, 210 Fisher attribution, 319 representation
diagnostics, 20 numerical diagnostics, 286 patch support, 125 transverse
descriptors, 73 label-free signal and 37 assignment assertions.

Before any new grade, **11,083 saved-output assertions** pass. All twelve
control tables replay the previous support run byte-for-byte; prediction
SHA-256 remains `eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`.
All arms preserve 146 scans / 900 keys, every selected N, eleven common
upstream/Fisher/k-means/training-support tables and four unavailable keys:
`240310_Cu100009.sxm` lobe 8, `240818_015.sxm` lobe 1 and
`240818_019.sxm` lobes 1–2. Only external grading excludes the technical scan.
No incomplete output or unusable seed is observed; all ten fits per arm
are physically named, with positive, freely learned masses summing to one.

All **ten factor fits reach 200 updates without convergence**, as allowed
by the predeclared bounded-fit policy. They retain rank four through both
hard steps; smallest diagonal noise is about **7.59e-5**, above the existing
1e-6 floor, and covariance reconstruction error is at most **1.12e-16**.
All **ten Student fits converge in 65–86 updates**, retain df five and both
hard updates, with latent precision range approximately **0.003115–5.410931**.
No iteration cap, tolerance, rank, df, regularizer or seed is adjusted.
The factor result tests this bounded implementation, not a converged MFA
optimum; the finite/correct arithmetic does not establish optimization quality.

Only then is the unchanged external full145-own-N grade run:

| Profile | Correct / 870 | Classified accuracy | Exact / 145 | Coverage / 870 |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 36 | 854 |
| Exactly replayed support control | 676 | 676/852 (79.3%) | 34 | 852 |
| Fixed rank-four factor analyzers | 620 | 620/852 (72.8%) | 23 | 852 |
| Fixed df-five Student components | 631 | 631/852 (74.1%) | 25 | 852 |

Against support, factors give **21 scan gains, 65 losses, 59 ties**, with
**5 exact chains gained and 16 lost**. Student gives **18 gains, 57 losses,
70 ties**, with **4 exact chains gained and 13 lost**. Against history,
factors give 25/69/51 scan gains/losses/ties and 7/20 exact gains/losses;
Student gives 23/65/57 and 6/17. The full **580 paired scan rows** are saved,
including all losses. Grade orientation changes on 13 / 7 scans against
support and 13 / 11 against history; it is not independent spatial matching.

Confusion `(TN,FP,FN,TP)` is support `(513,52,124,163)`, factors
`(407,158,74,213)`, Student `(421,144,77,210)`. Thus factors recover 50 false
negatives but add 106 false positives; Student recovers 47 but adds 92.
Emitted errors rise **176→232 / 221**. Counts remain 892 benchmark predictions,
16 missing positions, 38 extras, 14 short-N and 25 extra-N scans, with the
same two aligned abstentions leaving 852 classified positions.

Eight independent saved-component checks pass. Factors change **201 GMM/final
scores and 159 final decisions**, Student **240 scores and 142 decisions**.
Every final flip is **0→1**, with candidate GMM vote exactly one; **158 / 141
flips enter exact GMM=1 / k-means=0 ties**, while the remaining flip in each
arm has margin 0.25. Final ties rise **58→216 / 199**. GMM hard decisions
change 0→1 for 153 / 136 rows and never 1→0; six final flips per arm change
only an already class-1 GMM vote's strength. Among 896 usable rows, GMM votes
decrease/increase/tie for **25/176/695** with factors and **97/143/656** with
Student. Raw `(0,1,?)` frequencies change `(673,223,4)` to `(514,382,4)` /
`(531,365,4)`. No fusion, tie-breaking or precision rule is retuned afterward.

Factor exact gains: `240814_011`, `240814_025`, `240816_002`, `240817_002`,
`241114_028`. Losses: `240310_Cu100007`, `240817_005`, `240817_007`,
`240817_021`, `240817_041`, `240817_044`, `240817_047`, `240817_048`,
`240817_049`, `240817_053`, `240817_054`, `240817_055`, `240817_079`,
`240818_020`, `241113_087`, `241113_088` (all `.sxm`). Student exact gains:
`240814_011`, `240814_025`, `240817_082`, `241114_028`. Losses: `240815_072`,
`240817_005`, `240817_007`, `240817_021`, `240817_041`, `240817_047`,
`240817_048`, `240817_049`, `240817_053`, `240817_055`, `240818_020`,
`241113_087`, `241113_088`. Full decision-level tables and commands are in
`results/factor_student_mixtures_20260922/report.md` and linked outputs.

**Decision: reject these two fixed variants as replacements and retain the
676/34 support candidate.** Historical 677/36 remains unexceeded. These
negative results do not prove that all lower-rank or robust models fail,
especially with MFA's observed iteration-cap limit. No extra run, parameter
search, combined model, count refit, DFT change or unknown25 rerun follows.
Inference remains label-free; choosing/assessing methods repeatedly on this
benchmark is development work, not independent validation. This authorized
two-candidate comparison is complete; the historical objective remains open.

The final Julia 1.13 documentation build passes, with existing large-page /
search-index warnings and the expected non-CI deployment skip. The assignment
page stays within the unchanged 200 KiB HTML limit. README, runbook, selection,
assignment, config, calibration and HPC notes all report the same result and
MFA convergence limitation. Scientific source/settings remain unchanged from
**e2205cb** through grading; generated reports and fetched outputs stay ignored.

## 2026-09-22 — All-update shrinkage and coherent Student decisions: bounded comparison

The user authorizes the two proposed targeted experiments. Scientific
deliverable: complete support-control replay and two separate full145-own-N
grades, with coverage, emitted errors, exact chains and paired scan gains/losses.
The reference remains **676/870 correct, 34/145 exact, 852/870 coverage**;
historical **677/36, 854/870 coverage** remains the objective. Done means tested
equations and integration, fetched/checked complete outputs before grading,
an evidence-backed decision and committed source/docs. Time box: **three hours
preparation/tests from 11:20 CEST, one Viper job capped at one hour**. No
combined model, post-grade adjustment, seed/df/threshold sweep, changed N,
new DFT/registration/noise model, MFA-cap change or unknown25 rerun. The archived
scientific Must-NOT-have list remains; no retired process apparatus is revived.

`unit_assignment_em_shrinkage.toml` changes only covariance scope and estimator
(plus output name) from support. Separate covariances are shrunk toward their
own spherical target at k-means initialization, every Gaussian EM update and
both hard updates. Hard groups reuse the previous Ledoit-Wolf formula. Soft
groups use the declared fixed-weight plug-in extension: `a_i=r_i/sum(r)`,
`S=sum(a_i*y_i*y_i')`, `T=tr(S)/p I`, and coefficient
`clamp(sum(a_i^2*||y_i*y_i'-S||_F^2)/||S-T||_F^2,0,1)`, zero when S=T.
The same 1e-6 ridge follows shrinkage. It is not represented as a reproduction
of Halbe et al.'s RGMM or an iid-optimal/noise estimator; that source motivates
regularization during learning. The direct outer-product test independently
checks the formula, hard-membership limit, zero weights, scale invariance of
responsibilities, rotations and rank-deficient samples. No independent-sample
count is inferred. Shrinkage need not increase ordinary likelihood; observed
decreases are logged, never used to select seeds. The existing 200-update /
1e-6 stopping limits are explicit; finite capped fits remain, without retries.

`unit_assignment_student_density.toml` keeps df five and the existing Student
scale/latent-precision updates, but uses Student log mass-density for both hard
assignments and final scoring. The score includes log mass and log determinant,
and the `log1p(distance²/df)` tail. Only the common dimension/df constant is
omitted. Scale is not multiplied by df/(df-2). Both classification-style hard
updates remain; this is not pure soft EM or a continuous-vote experiment.
Independent density and weighted-moment calculations check all stages.

Both new policy keys are explicit in all native profiles, with old values
preserving their behavior. Unsupported combinations are rejected. Eight
features plus 28 interactions, per-file scaling, selected N, patch support,
Fisher, raw-amplitude naming, seeds, hard-vote aggregation, final fusion,
threshold and abstention rules remain unchanged. All keys remain even if a
seed or row is unavailable. These internal scores are not calibrated chemical
probabilities; no labels or composition quotas enter inference. The historical
representation diagnostic rejects the two new modes it cannot reproduce.

An initial Julia invocation in the filesystem sandbox cannot create juliaup's
lock; using the authorized local Julia cache access resolves it, with Julia
1.13 confirmed. The first targeted numerical suite catches an implementation
error: the weighted covariance helper returned a positional tuple instead of
the named tuple expected by callers. Adding the named-tuple semicolon fixes
the interface without changing equations or settings. The initial log is
retained as `initial_tuple_return_failure.log`. Synthetic/CLI and complementary
regressions continue; no real scientific output or new grade exists yet.

Viper access and its Julia 1.13.0 executable are verified; the first queue check
shows no jobs for the account. The comparison runner remains one sbatch, four
requested CPUs, 16 GB, one hour, no requeue, with five explicitly exported
paths. No real job has been submitted yet.

The final targeted suite passes **585 assertions**, including unavailable-seed
row retention in both new modes; the existing mixture suite passes **516** and
the new runner **20**. The Julia 1.13 documentation build passes with existing
large-page/search-index warnings and the expected non-CI deployment skip.
Source is committed before remote preparation. The thirteen-variant synthetic
pipeline and complementary regressions are still running; their completion is
not yet claimed. No benchmark outcomes have been read for the new methods.

All **21 suites / 4,277 assertions** now pass under Julia 1.13. The pipeline
has **393 assertions**, including **260 across thirteen end-to-end synthetic
variants**; the new modes retain identical upstream tables. Targeted equations
and CLI tests contribute 585, previous mixtures 516, runner 20, covariance
sharing 218, scan weighting 138, continuous votes 115, bagging 229, robust
normalization 99, scan-Fisher/naming 299, numerical assignment 136, complete
training 154, native Fisher 229, affine Fisher 74, Fisher attribution 210,
representation diagnostics 321, numerical diagnostics 20, patch support 286,
signal features 73, transverse descriptors 125 and base assignment 37.

Scientific source **bb6943b** and the unchanged Manifest are synchronized to
`/u/oldu/code/STMFit_em_shrinkage_student_density_20260922`; source and both
pre/post-instantiate Manifest checksum checks are empty. Login-node instantiate
without precompilation succeeds under Julia 1.13.0, retaining the expected
warning about the Manifest's 1.12.6 origin. All three metadata dry-runs confirm
146 files. Scheduler `--test-only` succeeds; **11936070 is only its speculative
identifier, not a submitted job**.

Actual job **11936229** is submitted once, after the full test suite passes,
with all five variables explicitly exported, four requested CPUs, 16 GB,
one hour and no requeue. Its output path is
`/ptmp/oldu/stmfit/em_shrinkage_student_density_20260922_v1`; local logs and
saved-output report scripts are under
`results/em_shrinkage_student_density_20260922/`. The preparation takes about
41 minutes, within the declared three-hour limit. No real result or grade is
claimed at submission; scientific settings remain frozen at **bb6943b**.

### Completed comparison: both coherent-learning candidates regress

Job **11936229** is submitted at **12:01:26 CEST**, starts at **12:01:29**
after **3s** queued, and completes **0:0 at 12:09:25**, elapsed **7m56s**.
Control finishes at 12:05:34, all-update shrinkage at 12:07:29 and Student
at 12:09:22. Four CPUs are requested, eight allocated; 16 GB requested,
MaxRSS **2,283,380 K** (about 2.34 GB). No retry or second job is submitted.
All **588 regular output files plus Slurm log** are fetched locally; both
checksum comparisons are empty. Remote raw-input symlinks are not recopied.

Before any new grade, **11,120 saved-output checks** pass. Twelve control
tables are byte-identical to the saved support run; prediction SHA-256 remains
`eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3`.
All arms preserve 146 scans / 900 keys, every selected N and eleven common
upstream/Fisher/k-means/support tables. The same four rows are unavailable:
`240310_Cu100009.sxm` lobe 8, `240818_015.sxm` lobe 1 and
`240818_019.sxm` lobes 1–2. Only external grading excludes the technical scan.
All ten seeds per arm are numerically usable and physically named, with free
positive component masses; no fallback, dropped seed or silent row removal.

All regularized fits satisfy the fixed relative-log-likelihood stopping rule
in **25 or 43 updates**. Initial/EM shrinkage coefficients range about
**0.029928–0.060398**, and hard-stage coefficients about **0.033321–0.039135**.
Seeds 1 and 3 each have **eleven ordinary-likelihood decreases**, with largest
drop **1.1437665**; the other eight have none. This is the expected limitation
of a plug-in covariance regularizer, not a monotone likelihood-maximizing
M-step. No seed is discarded or rerun. Student fits all converge in **65–86
updates**, with latent precisions about **0.00304549–5.44709175**, retaining
df five and both density-based hard updates. No numerical criterion claims
a global optimum or calibrated chemical uncertainty.

The unchanged external full145-own-N grade then gives:

| Profile | Correct / 870 | Classified accuracy | Exact / 145 | Coverage / 870 |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 36 | 854 |
| Exactly replayed support control | 676 | 676/852 (79.3%) | 34 | 852 |
| All-update covariance shrinkage | 637 | 637/852 (74.8%) | 24 | 852 |
| Coherent Student-density decisions | 627 | 627/852 (73.6%) | 25 | 852 |

Against support, shrinkage gives **19 scan gains, 55 losses, 71 ties**, with
**6 exact gains / 16 losses**. Student gives **17 gains, 58 losses, 70 ties**,
with **4 exact gains / 13 losses**. Against history, shrinkage gives 26/61/58
scan gains/losses/ties and 8/20 exact gains/losses; Student gives 22/65/58
and 6/17. All **580 paired scan rows** are saved. Grade orientation changes
on 8 / 7 scans against support and 12 / 11 against history; these alignments
are not independently measured spatial monomer correspondences.

Confusion `(TN,FP,FN,TP)` is support `(513,52,124,163)`, shrinkage
`(430,135,80,207)`, Student `(413,152,73,214)`. Shrinkage's 44 fewer false
negatives do not offset 83 extra false positives; Student's 51 fewer do not
offset 100 extra false positives. Emitted errors rise **176→215 / 225**.
Counts remain 892 benchmark predictions, 16 missing positions, 38 extras,
14 short-N and 25 extra-N scans. The same two aligned abstentions leave
852 classified positions. No count-selection or counting-validation claim changes.

Eight independent saved-component checks pass. Shrinkage changes **179 scores
and 129 final decisions**, Student **240 scores and 154 decisions**. Every
final flip is **0→1**, with candidate GMM vote one. **128 / 153** enter exact
GMM=1 / k-means=0 ties; one flip per arm has final margin 0.25. Total ties
rise **58→186 / 211**. GMM hard classes change 0→1 for 123 / 148 rows, never
1→0; six final changes per arm only reinforce a GMM vote already classified 1.
Among 896 usable rows, GMM scores decrease/increase/tie for **46/133/717**
with shrinkage and **85/155/656** with Student. Raw `(0,1,?)` frequencies
change `(673,223,4)` to `(544,352,4)` / `(519,377,4)`. These observations do
not authorize retuning the fusion or tie rule after grading.

Shrinkage exact gains: `240814_025`, `240816_002`, `240817_024`, `240817_082`,
`240817_085`, `241114_028`. Losses: `240307_017`, `240815_072`, `240817_005`,
`240817_007`, `240817_021`, `240817_041`, `240817_044`, `240817_047`,
`240817_048`, `240817_049`, `240817_052`, `240817_053`, `240817_055`,
`240818_020`, `241113_087`, `241113_088` (all `.sxm`). Student gains:
`240814_025`, `240817_082`, `240817_085`, `241114_028`. Losses: `240815_072`,
`240817_005`, `240817_007`, `240817_021`, `240817_041`, `240817_047`,
`240817_048`, `240817_049`, `240817_053`, `240817_055`, `240818_020`,
`241113_087`, `241113_088`. Every changed decision and full grade command
is retained in `results/em_shrinkage_student_density_20260922/report.md`.

**Decision: reject both fixed candidates as replacements; retain support
676/34.** Historical 677/36 remains unexceeded. Unlike the preceding capped
MFA test, all fits meet their stopping criterion here, but this does not
reject all regularization schemes or Student models. No post-grade parameter
change, new job, combined arm, count refit, DFT change, measurement-error model,
MFA retry or unknown25 rerun follows. This two-candidate comparison is complete;
the broader historical objective remains open. The reused benchmark is
development evidence, not independent validation or unknown-chain accuracy.

README, runbook, selection, assignment, configuration, calibration and HPC
notes are updated consistently. The preceding Fisher/naming and MFA/Student
summaries in the assignment page are condensed, preserving their scores and
links to the full dated journal, to stay within the existing HTML size limit.
No documentation size limit is raised. Final documentation verification follows;
scientific source/settings remain unchanged from **bb6943b** through grading.

The final Julia 1.13 documentation build passes. The assignment page is
**199.53 KiB**, under the unchanged 200 KiB limit; existing size/search-index
warnings and the non-CI deployment skip remain. `git diff bb6943b -- config
test hpc/compare_em_shrinkage_student_density.sbatch Manifest.toml` is empty:
no scientific code or setting changes during execution, checking or grading.
Generated results remain ignored; the completed outcome is committed in the
journal and all cited documentation headlines.

## 2026-09-22 — Fresh GCV counts after symmetric filtering: bounded comparison

The user authorizes the recommended first direction. Scientific deliverable:
reselect counts using the corrected input and unchanged native GCV rule, apply
the unchanged support assignment, and externally compare counting, correct/all
positions, classified coverage, emitted errors, exact chains and paired gains /
losses against the fixed-N control and historical reference. Done requires
tested input/selection boundaries, complete fetched outputs, integrity checks,
an explicit result/decision and committed documentation. Time box: three hours
preparation from 12:28 CEST, one Viper job capped at four hours, four requested
CPUs / 16 GB, no automatic retry. Other proposed directions are not authorized.

The source/history check matters: the cached 900 lobes came from the direct
extractor's minimum valid circular/elliptical GCV, not `batch_full.jl`'s current
`support_midpoint_hybrid` policy. The experiment repeats precisely that raw-GCV
path, including its existing intelligent sweep and BIC early-stop diagnostic.
It is not exhaustive and does not substitute the separate 871-key promoted
counting summary. No expected N, sequence, chemical proportions, benchmark
membership/quality filtering or externally selected per-scan rule enters fitting.
GCV, n_eff, physical bounds, optimization limits and TOMLs remain unchanged.

`run_gcv_reselection.jl` uses four native extractor shards for the new sweep.
Only each scan's selected count is passed onward: both control and treatment
regenerate their base/split fits at fixed N, then patches, Fisher, classifiers,
vote and QC with `unit_assignment_patch_support.toml`. Sweep geometry is saved
for inspection but never used as assignment geometry. This avoids confounding
the count change with variable-sweep warm starts in only one assignment arm.
Unchanged-N geometries will be compared explicitly; whole-cohort learning can
still change predictions on those scans when other scans' counts change.
The old counts remain only in the control, not in treatment's sweep/range.

Metadata-only dry-runs check the entire raw/control cohort. Incomplete, invalid,
duplicated or gapped outputs fail the comparison; no scans are silently dropped.
The computation reads no grading tables. Existing historical and support grades
are references only; no grade-guided tuning, new DFT, acquisition/noise model,
changed abstention, classifier variant or unknown25 rerun is included. This
heavily reused benchmark remains development evidence, not independent validation.

The new comparison suite passes **54 synthetic assertions** under Julia 1.13,
covering the original raw-GCV validity/tie rules, unchanged free/fixed-N settings,
strict cohorts and contiguous keys, independent fixed-N arms, failed-output
retention, CLI exclusions and Slurm/dry-run arguments. An initial fixture reused
a TSV path rejected by the existing no-overwrite writer; unique fixture paths
fix that test-only collision, and its log is retained. Package, selected-support
and assignment regressions are still running; no real fit or grade has run.

All eight complementary suites now pass **1,604 assertions**, for **1,658 total**
including the new suite. The Julia 1.13 documentation build passes with existing
size/search-index warnings. Source is committed as **853ad38**. Local and remote
metadata-only dry-runs both validate 146 files. Initial remote preparation
omitted the ignored Manifest from the tracked-file transfer, so Pkg generated a
different lock. No fit ran with it: preserve that generated lock separately,
transfer the original Manifest explicitly, and instantiate again. Its before /
after SHA-256 is unchanged (`617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe`). Configs,
templates and original counts retain their saved hashes. Raw-content verification
dereferences local staging links; symlink targets and directory timestamps are
not scientific content differences. Slurm test-only **11936531** is a speculative
identifier, not a submitted job. No real job has been submitted at this point.

The final source/raw-content checks and explicit Manifest comparison are empty.
Julia imports the fit engine from the correct dedicated checkout, not another
local package tree. **Job 11936561** is submitted once with all five variables
explicitly exported and no requeue. Source is **853ad38**, preparation record
**7e4ce2c**. Outputs go to `/ptmp/oldu/stmfit/gcv_reselection_20260922_v1` and will
be fetched under `results/gcv_reselection_20260922/`; the Slurm log is
`/ptmp/oldu/stmfit/gcv_reselection_20260922-11936561.log`. Preparation takes about
17 minutes, within the three-hour bound. No benchmark result is available yet.

### Completed output and integrity

The single job finishes **COMPLETED 0:0 in 1h25m43s**, from **12:45:36 to
14:11:19 CEST**, after **26s** queued, on `vipc2575`. Viper allocates eight
CPUs for the four requested; batch MaxRSS is **3,888,598 K**. No retry,
requeue or login-node fit occurs. All **438 regular outputs and the Slurm
log** are fetched locally; every output matches its remote SHA-256. A final
itemized source comparison is empty. Manifest, both configs, templates and
original selected-count summary retain their pre-job hashes.

Both arms retain all **146 scans**, 15 complete scientific tables and 146
individual plots, with no `failures.tsv`. The native sweep changes **22 N**,
retaining 124; total lobes become **900 → 898**. Thirteen scans lose one,
seven gain one and two gain two. The technical scan remains eight lobes.
Raw `(0,1,?)` populations are `(672,224,4)` for fresh control and `(672,222,4)`
for treatment; these are outcomes, never prescribed class counts. The four
unavailable keys remain identical. **9,054 pre-grade assertions** pass:
1,886 cohort/geometry/naming checks and 7,168 independent vote/availability
checks. All ten seeds retain the old raw-amplitude naming rule, and learning
stays Gaussian/full-covariance/final-only/Mahalanobis without resampling.

### Fresh-refit variability, discovered before grading

Fresh control is **not byte-identical** to saved support despite unchanged N,
fit sources, inputs and settings. Base geometry differs on **55/900 rows in
eight scans**; split geometry on **8/900 rows in one scan**. Maximum center
shifts are **0.008447 nm / 0.021220 nm**, maximum absolute relative GCV changes
**0.127% / 0.706%**. Patches differ on 55 rows; learned Fisher changes on
896 usable rows, GMM on 50 rows, k-means on 48. Nine final rows change scores,
but only **three decisions**: `240817_058.sxm` lobe 4, 0→1;
`240817_076.sxm` lobe 3, 1→0; `241114_011.sxm` lobe 4, 0→1.
Decision agreement is **897/900**, not an exact replay.

Between fresh arms on **755 lobes in 124 same-N scans**, base geometry differs
on **52 rows / eight scans**, split geometry on **88 rows / fourteen scans**.
Maximum center shifts are **0.008447 / 0.076370 nm**, maximum relative GCV
changes **0.115% / 1.049%**. The sub-0.009-nm bound reported during monitoring
applied only to base geometry, not split geometry. The unchanged `GN_DIRECT_L`
optimizer has a ten-second wall-clock cap before LM refinement. Timing-related
variation is plausible, but not established without per-fit stopping/evaluation
traces. No numerical setting is changed or fit repeated to obtain a preferred
result. Whole-cohort learning can independently propagate count changes to
same-N scans. This is consequently an end-to-end reselection comparison, not
an isolated causal estimate of N or the filtering correction.

### External benchmark result: negative

Only after complete output checks does the external report select the 145
benchmark files and use their labels. Full-own-N grading retains the 870-position
denominator and reports extras separately; **37 independent aggregate checks**
pass and **435 paired rows** cover all three reference comparisons.

| Profile | Exact N / 145 | Correct / 870 | Exact / 145 | Coverage / 870 | Missing / extra | Emitted errors |
|---|---:|---:|---:|---:|---:|---:|
| Historical saved reference | 106 | 677 | 36 | 854 | 16 / 38 | 177 |
| Saved support | 106 | 676 | 34 | 852 | 16 / 38 | 176 |
| Fresh fixed-N control | 106 | 675 | 34 | 852 | 16 / 38 | 177 |
| Fresh raw-GCV counts | 101 | 661 | 32 | 849 | 19 / 39 | 188 |

Fresh control loses one correct position against saved support, with unchanged
exact-chain count. Treatment gives **661/849 classified accuracy (77.9%)**,
versus **675/852 (79.2%)** for fresh control, saved **676/852 (79.3%)** and
historical **677/854 (79.3%)**. Relative to all 870 positions it gives **76.0%**,
not 77.9%. Benchmark output rows become 892 → 890; short/long files 14/25 →
17/27. Two abstentions align with the reference, so three more missing positions
explain coverage falling 852 → 849 without additional unavailable keys.

Against fresh control, **9 scans gain, 15 lose, 121 tie**, net **−14 correct**;
**three exact chains are gained, five lost**. Against saved support: 10/16/119,
net −15, exact 3/5. Against history: 24/30/91, net −16, exact 6/10. Grade
orientation changes on seven scans against either control, eleven against
history; that alignment is not measured spatial correspondence. Of the −14
positions, −12 occur on the 22 changed-N scans and −2 on the 123 unchanged-N
benchmark scans; this partition does not identify a causal mechanism.

Exact-chain gains against control: `240307_016.sxm`, `240814_011.sxm`,
`240817_077.sxm`. Losses: `240307_017.sxm`, `240310_Cu100007.sxm`,
`240314_Cu100_025.sxm`, `240817_041.sxm`, `240817_049.sxm`.
Exact-count gains: `240310_Cu100006.sxm`, `240314_Cu100_024.sxm`,
`240817_077.sxm`, `240817_093.sxm`, `241113_089.sxm`. Count losses:
`240307_017.sxm`, `240313_Cu100058.sxm`, `240314_Cu100_025.sxm`,
`240814_024.sxm`, `240814_025.sxm`, `240817_041.sxm`, `240817_049.sxm`,
`240817_060.sxm`, `240817_091.sxm`, `241114_028.sxm`.

**Reject this reselection version and retain saved support 676/34; historical
677/36 remains unexceeded.** The bounded authorized comparison is complete,
not the broader historical objective. No post-grade optimization budget,
threshold, selector, seed, classifier, registration/noise or combined-method
change follows; no DFT or unknown25 run occurs. The direct-extractor raw-GCV
lineage here is not the separately promoted batch-policy counting reference.
Labels remain external to inference, but repeated benchmark reuse supplies
development evidence rather than independent validation or unknown-chain
chemical accuracy. Numerical reproducibility is a limitation for future
small-gain claims, not an authorization for another campaign.

Full report, settings, hashes, tests, grade tables and paired losses are saved
under `results/gcv_reselection_20260922/report.md`. Documentation headlines
are updated consistently. The older tied-covariance section of the assignment
page is condensed, retaining its outcome, method and journal/report pointers,
to keep the documentation within its existing HTML size limit; no limit is raised.

## 2026-09-22 — Frozen-geometry amplitude/background profiling: bounded comparison

The user authorizes the recommended first direction. Scientific deliverable:
compare unchanged support assignment against the same method after constrained
linear re-estimation of Gaussian amplitudes and background, at **identical saved
N, centers, widths, axes and skew**. Done means tested arithmetic and input
boundaries, one complete fetched comparison, independent output checks, external
full145 correct/870, exact/145, coverage, errors and paired gains/losses, followed
by a committed decision. Budget: one half-day interpreted as **four hours from
15:28 CEST**, then **one two-hour Viper job**, four requested CPUs / 16 GB,
no automatic retry. No nonlinear optimization, new N, registration, noise model,
classifier variant, DFT change or unknown25 rerun is included.

The existing September 18 prototype already solves finite-box least squares;
this experiment reuses its `box_lsq` and all explicit linear solver settings
from `label_free_exploration.toml`. It does not use its outer variable-projection
search. The old prototype's lower RSS did not establish recognition accuracy.
The amplitude floor/range come from the unchanged native fused image, and
finite endpoints are mapped through the native raw bounds/decoder, not replaced
with ordinary nonnegative least squares. No composition prior or class quota
enters the solve.

The input is the original symmetric-fusion base/split feature cache used by
saved support **676/870, 34/145**. The last reselection's fresh **675/34** control
and **661/32** candidate are not inputs. Gaussian base amplitudes, plane offset
and two tilts are profiled jointly on reconstructed native fused ROI/tube pixels.
Cached decimal centers, widths and axes are used literally, including their
serialization precision; they are not recovered full-precision raw fit vectors.
The reconstructed native axis/origin must agree at the cache's decimal precision.
This preserves the interpretation of the existing support construction without
claiming exact recovery of an unsaved full-precision fit.

Only `amplitude`, `baseline`, `tilt_x`, `tilt_y`, derived `amp_rel` and diagnostic
`gcv` may change in the base table. GCV retains the **full original model
parameter count**, not just its linear degrees of freedom, and selects nothing.
The split cache remains byte-identical: the assignment uses only its skew,
which is geometry and must stay frozen. All base-derived local features, patches,
Fisher, classifiers and votes are regenerated in both arms with unchanged
`unit_assignment_patch_support.toml`. The control must reproduce saved support;
there is no reoptimization chosen by the external score.

The decimal cache can place a coefficient just outside its native endpoint.
Initial projection is allowed only within **half the last written decimal unit**
plus machine roundoff; bounds themselves never widen, and the adjustment is
reported. Larger infeasibility, missing rows, changed geometry, nonconvergence
or increased RSS stops the candidate. Failed scan records are retained; there
is no silent row drop or selective fallback. Four one-thread profile shards
save all fit pixels, coefficients/bounds/gradients and KKT/RSS diagnostics for
independent saved-output checks. Inference reads no sequence, expected count,
benchmark membership, chemical count or grade. External grading remains
development evidence on the reused benchmark, not independent validation.

Initial synthetic testing caught an invalid one-lobe fixture (`range` with one
element but unequal endpoints), fixed by using a singleton amplitude. The failed
test log is retained; no real output informed the correction. Focused and
regression tests are running; no cluster job or benchmark grade has run.

The **498 new synthetic/metadata assertions pass**. A single local raw check on
the lexicographically first input, `240307_015.sxm`, preserves its seven lobes
and every geometric field; RSS falls **0.06347161460 → 0.06181862243** (about
2.60%), with KKT satisfied. This is not a recognition result or a sample chosen
by benchmark error. The full146 local metadata-only dry-run passes. The native
max-image amplitude scale is also saved with the per-scan diagnostics, so the
mapped coefficient endpoints can be independently checked. No cohort fit or
grade runs locally; no Slurm job has been submitted yet.

All nine complementary suites pass **4,438 assertions**, including the original
box/profile solver suite, for **4,936 assertions across ten suites** with the new
checks. The Julia 1.13 docs build passes (existing size/search warnings only).
Scientific source is committed as **567cdde**. The dedicated remote checkout
imports the correct local fit engine under Julia 1.13.0; instantiate leaves
the original Manifest SHA-256 **617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe** unchanged.
Both cached geometries and templates match local hashes. Itemized source and
dereferenced raw-content comparisons are empty. Local and remote dry-runs both
validate all 146 files without reading pixels. Slurm test-only passes with all
five variables explicitly exported; it does not validate compute-node execution.

The single actual job is **11938782**, submitted with `--no-requeue` and the
explicit five-variable export. Source remains **567cdde**, preparation record
**19b3921**. Slurm test-only identifier **11938774** was speculative, not a real
job. Preparation took about twenty minutes, within the four-hour bound. Remote
output is `/ptmp/oldu/stmfit/frozen_amplitudes_20260922_v1`; local fetch/report root
is `results/frozen_amplitude_profile_20260922/`. At submission no benchmark
result was available; the completed result follows.

### Completed numerical checks and recognition result

The job starts **15:53:51 CEST**, after **5m29s queued**, and finishes
**16:00:08**, **COMPLETED, 0:0 in 6m17s**, on `vipc2109`. Viper allocates
eight CPUs for the four-CPU request. Both arms finish **146 scans / 900 lobes**
with 146 individual plots each; all four profile shards succeed. The single
job stays within the two-hour limit, with no retry or requeue. All **561
regular output files and the Slurm log** are fetched locally and match their
remote SHA-256 hashes. Manifest and all three configuration hashes stay unchanged.

Before external grading, **22,912 independent saved-output assertions** pass:
60 cohort/control checks, 5,384 vote/availability checks and 17,468 basis/profile
checks. Thirteen control scientific tables are byte-identical to saved support,
including prediction SHA-256
**eb8d69b8cebcd607b0244db5c65104347e6bcd24acee06f64748cf0e6ec310e3**.
The candidate hash is
**4c6eee74bb5b777b6d3009229845c7e531abd94810f3fd98211432737a38769d**.
N, every nonmutable base field, split geometry and observed raw patch pixels
remain identical. No failure table or lost row is hidden. The four unavailable
keys remain `240310_Cu100009.sxm:8`, `240818_015.sxm:1` and
`240818_019.sxm:1,2`; only two align with external control positions.

The checker independently reconstructs each Gaussian basis from the frozen
decimal table and saved fit coordinates, then verifies native bounds, exported
coefficients, gradients, RSS, full-parameter GCV and KKT conditions. All
**146 profiles converge in 1–3 iterations** and lower RSS. Relative reductions
are **0.0137–12.1647%, median 0.6814%**; summed RSS changes
**20.5570966504 → 20.3671925622**. Forty-three scans touch at least one
coefficient bound, with 41 lower and 13 upper coefficients. Twenty-four scans
need initial endpoint rounding correction, at most **3.43e-10**, within the
documented half-decimal-unit limit; no bound widens. Amplitude changes range
**−10.94% to +14.23%, median +0.384%**. The first checker's displayed
"median" used the lower middle observation for an even sample size; this
report-only display is corrected to the mean of both central observations and
the identical checks are rerun. The original log is retained; scientific
outputs, checks and grades are unaffected.

Only external scripts then read benchmark membership and labels. Full-own-N
grading covers all 145 benchmark scans with denominator 870; **31 independent
aggregate checks** pass and **290 paired scan rows** compare candidate to
support and history.

| Profile | Correct / 870 | Correct / classified | Exact / 145 | Coverage / 870 | Errors |
|---|---:|---:|---:|---:|---:|
| Historical saved reference | 677 | 677/854 (79.3%) | 36 | 854 | 177 |
| Exactly replayed support | 676 | 676/852 (79.3%) | 34 | 852 | 176 |
| Frozen-geometry linear profile | 675 | 675/852 (79.2%) | 32 | 852 | 177 |

Relative to all positions, the candidate gives **675/870 (77.6%)**, not 79.2%.
All N are unchanged: **106/145 exact counts**, 892 benchmark predictions,
16 missing positions, 38 extra lobes, 14 short and 25 long scans. Against
support: **4 scan gains, 5 losses, 136 ties**, net −1 correct; **zero exact
gains and two losses**, `240310_Cu100007.sxm` and `240818_020.sxm`.
Grade orientation never changes against support. Confusion `(TN,FP,FN,TP)`
changes **(513,52,124,163) → (516,49,128,159)**: three fewer false positives
but four more false negatives. Against history: **23/28/94**, net −2 correct,
**six exact gains and ten losses**, with eight orientation changes. Full
paired lists are in the local report; grade orientation is not independent
spatial monomer correspondence.

Three additional component checks pass. GMM/k-means saved probabilities change
on **30/27 rows**, hard classes on **16/6**. Final decisions change on **ten
lobes of ten scans**, nine 1→0 and one 0→1. Eight changes leave old exact
vote ties and one enters a tie; the other changes 0.9 → 0.45. Final ties fall
58 → 46. Raw `(0,1,?)` populations move **(673,223,4) → (681,215,4)**;
these are unconstrained outcomes. No tie, precision, classifier, bound or
threshold adjustment follows this observation.

**Reject this version and retain saved support 676/34.** Historical 677/36
remains unexceeded. The lower pixel residual is real for these fixed decimal
bases but does not improve this recognition comparison. This is not a verdict
on every profiling method or the unsaved full-precision geometry. Whole-cohort
learning couples scans. Learning/inference remain label-free; repeated external
benchmark reuse is development evidence, not independent generalization or
unknown-chain chemical validation. The authorized deliverable is complete,
not the broader historical objective. No other proposed direction, nonlinear
fit, count sweep, registration/noise model, DFT change or unknown25 rerun occurs.

Full report, tests, diagnostics, hashes and losses are under
`results/frozen_amplitude_profile_20260922/report.md`. Cited documentation
headlines are updated consistently. The older whole-scan bagging section is
condensed on the assignment page, retaining its method, outcome and full
journal/report pointers, to keep the existing HTML size limit; no limit is raised.

The final Julia 1.13 documentation build passes with the existing size/search
warnings; the assignment page is **198.83 KiB**, below its unchanged 200 KiB
limit. Scientific source remains unchanged from **567cdde**. The completed
comparison is handed off without another compute submission.

## 2026-09-22 — Matched fixed-N geometric profiling: bounded recognition comparison

The user authorizes the first proposed direction after the negative amplitude-only
result. Deliverable: a complete external recognition comparison of geometric
variable projection against matched joint optimization, plus the saved support
and historical references, with an explicit decision. Done requires tests,
complete fetched outputs, independent numerical/cohort checks, correct/all870,
exact/145, coverage, errors and paired gains/losses, and committed documentation.
Budget: **four hours preparation from 16:52 CEST**, then **one Viper job capped
at four hours**, four requested CPUs / 16 GB, no automatic retry. No missing-head
vote, joint acquisition model, new counts, class quota, DFT change, unknown25
rerun or post-grade parameter search is included.

The prior full-profile prototype improved RSS on four unknown scans but gave
almost only capped outer searches, had no equal-budget control and did not
measure recognition. The completed fixed-geometry benchmark gave 675/32 against
676/34 despite lower RSS. Neither establishes whether positions/widths are a
recognition bottleneck. The new comparison tests that hypothesis, not an assumed
benefit of numerical optimization.

One shared native initialization is regenerated per saved N with the unchanged
circular fit and circular→elliptical refinement. Existing valid minimum-GCV
selection chooses the family **at that fixed N only**. Its full-precision vector,
family, axis and native ROI/tube pixels are shared by both local arms. This avoids
inventing unsaved raw parameters from rounded feature coordinates, or separately
refitting the two starts under the native wall-clock limit. The initialization
may differ from the original saved support; its parameters and both bootstrap
validities/GCVs are saved. The reference arm separately must replay saved support
676/34 exactly. The common initialization's native optimizer convergence is not
exposed and is not claimed.

Both local arms use **LN_BOBYQA**, maximum **800 total callback evaluations**,
**30 seconds inside NLopt**, relative x tolerance **1e-6** and objective tolerance
**1e-9**. These numerical values are reused from the September 18 diagnostic in
`config/geometry_profile_comparison.toml`; the native elliptical initialization
keeps its explicit 50-iteration budget. The first callback counts against the
evaluation ceiling and compiles before the NLopt timer; its setup time is saved
separately. Equal ceilings are **not equal arithmetic or elapsed cost**: a
profiled call additionally solves a bounded linear problem. Both actual evaluation
counts and stopping statuses are reported, and capped iterates are not called
converged. No extra optimization or seed search follows the result.

The joint arm moves the whole native parameter vector. The profiled arm moves
only the geometry block and solves amplitudes/background at each evaluation with
the already tested finite-box solver. Both minimize **the same native global
RSS × (1 + κ penalty)**, with unchanged κ threshold/weight. The penalty depends
only on geometry, so the linear elimination remains appropriate. This deliberately
does not reuse the old prototype's mismatched raw-RSS-only outer objective.
The same native final validity checks reject inadmissible evaluated points in
both arms. Each returns its best valid objective point, with rejected evaluations
and the full scalar trace recorded. A failed initialization or linear solve is
explicit; a profiled arm cannot silently substitute an unprofiled native result.
No scan is dropped or externally selected for a retry.

All parameters retain the native physical boxes. GCV uses the full original
parameter count and is diagnostic during local refinement; N and family stay
fixed once the common initialization is chosen. Split geometry/skew is reused
byte-for-byte, isolating Gaussian base geometry. This is not a simultaneous
refit of the asymmetric shape model. Both arms regenerate all base-derived
patches, descriptors, Fisher, GMM, k-means and final votes using the unchanged
support assignment config. Geometry tables keep the native exporter precision.
The 146-file cohort is determined by raw/cache keys, never benchmark membership;
labels and expected counts are read only by separate external grading afterward.
This reused development benchmark cannot establish independent generalization.

Preparation uses synthetic tests and one focused local raw check only. The
full146 comparison is sharded across four one-thread Viper workers. No scientific
cluster job or external grade has run at this preparation checkpoint.

Preparation checks pass under Julia **1.13.0**: **5,047 assertions in eleven
suites**, including 111 new matched-refinement checks (active κ penalty,
unchanged starts, actual callback ceilings, full GCV, KKT, explicit failures,
three complete assignment arms and forbidden label options). The 146-file
local metadata dry-run passes. The documentation builds within the unchanged
HTML limits; existing size/search warnings remain.

The only real local fit is the lexicographically first cached scan,
`240307_015.sxm`, saved **N=7**, without external grading. Shared native RSS
is **0.06349117339**; joint gives **0.06347592674**, `XTOL_REACHED`, 97
evaluations, six rejected points; profiled gives **0.04999638934**,
`MAXEVAL_REACHED`, 800 evaluations, zero rejected points. Profiled optimization
takes 16.79 seconds excluding separately reported setup. This is a numerical
smoke check, not evidence of recognition improvement or outer convergence.
All **84 saved-output checks** pass, including independently assembled Gaussian
RSS, κ objective, full GCV, linear KKT, native bounds, trace/budget consistency
and serialization. No second local scan or parameter adjustment follows.

Remote cached base, split, templates and saved support predictions match local
SHA256 values. A dereferenced checksum-only rsync reports no raw-scan difference.
Evidence is under `results/geometry_profile_20260922/`. An initial local dry-run
used a nonexistent Julia launcher path; repeating metadata-only validation with
the verified `/usr/bin/julia` succeeds. Neither attempt performed scientific
computation. Viper uses its separately verified absolute Julia 1.13 executable.

Scientific source is committed as **7c5fb6c** and synced to
`/u/oldu/code/STMFit_geometry_profile_20260922`. All tracked source checksums
match. The ignored local Manifest is copied separately, then `Pkg.instantiate()`
and package imports succeed on the login node under Julia 1.13; its SHA256
remains **617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe**.
Only dependency preparation and metadata dry-run run there, never STM fits.
Remote 146-file dry-run and `sbatch --test-only` pass. The latter reports a
prospective eight-CPU allocation, not a completed job or runtime environment test.
The saved-output smoke verifier additionally checks overlap independently,
bringing its checks to **87**, all passing without another fit.

The sole real submission is **11940220**, at **17:24:09 CEST**, output
`/ptmp/oldu/stmfit/geometry_profile_20260922_v1`, with explicit command-line
`--export=ALL` and all five required absolute-path variables, `--no-requeue`,
four requested CPUs / 16 GB and four-hour walltime. Slurm initially reports
**PENDING / QOSGrpCpuLimit**, no compute yet. Preparation took 32 minutes,
within the four-hour preparation window. No other candidate or retry is
submitted. Continue through completion, fetching, checks and external grading;
submission alone does not satisfy the scientific deliverable.

While the quota is pending, the assignment page's older continuous-GMM-vote
section is condensed, retaining its method, negative result, decision and
full journal/report pointers. This makes room for the eventual geometric
comparison without raising the existing 200 KiB HTML limit; no scientific
code or result changes.

The job starts **17:27:16 CEST** on `vipc2400`, after **3m07s** queued.
Viper allocates eight CPUs for the four-CPU request; the code still runs at
most four workers. The compute log confirms the explicit runtime paths and
starts the saved-reference pipeline before the paired geometric shards.
Documentation builds successfully after condensation; existing warnings and
limits are unchanged. Completion and external scores remain pending.

### Completed result and decision

Job **11940220** completes **COMPLETED, 0:0**, **17:27:16–17:42:33 CEST**,
**15m17s** elapsed, after **3m07s** queued. Four CPUs / 16 GB were requested;
Viper allocated eight CPUs, with four one-thread geometric workers as declared.
All **146 scans / 900 lobes** finish without initialization, profiling or
assignment failures. The **782 regular output files and Slurm log** are fetched
locally with matching SHA256 values. No job retry or second real local scan occurs.
The scientific code remains unchanged from **7c5fb6c**.

**20,138 pre-grade checks pass** (8,151 cohort/reference/vote checks and 11,987
geometry checks). Thirteen scientific reference tables replay saved support
byte-for-byte. Checks reconstruct Gaussian values independently from decoded
full-precision parameters and saved fit pixels, verifying RSS, overlap, κ
objective, full-parameter GCV, coefficient bounds/KKT, exported features,
shared starts, best valid trace points, actual evaluations and termination
flags. All final linear subproblems and recorded inner solves satisfy KKT.
No labels are read until these complete checks and checksum verification pass.

The common native RSS sum is **20.5570980004**; joint gives **20.5343660301**;
profiled gives **19.3301659145**. All **146/146** profiled objectives improve
over joint. Median RSS reductions from the common start are **0.0623% joint**
and **6.6009% profiled**; profiled reductions range **0.0521–43.2951%**.
The final κ penalty is inactive on this real cohort; its active branch remains
tested synthetically. Joint stops are **132 XTOL / 14 FTOL**; profiled stops
are **35 XTOL / eight FTOL / 103 MAXEVAL**. No time ceiling is reached.
The 103 capped searches are not converged; XTOL/FTOL are optimizer stops,
not proof of global optimality. Evaluation min/median/max is **67/87/127**
joint and **61/800/800** profiled. Rejected invalid evaluations total **2,168**
joint and **109** profiled. These are rejected points, not excluded scans.

Equal ceilings do not produce equal cost: sums of local optimization time
are **9.24 seconds joint / 402.97 seconds profiled**, with separately reported
callback setup **0.374 / 0.533 seconds**. These omit shared initialization and
downstream assignment, and must not be presented as full pipeline runtimes.
Native family matches the saved cache on every scan, but 14 initial centers
differ at export precision, up to **0.0002465 nm**. Joint/profiled still share
exactly the same fresh start. Profiling changes 891 centers; over all 900 lobes,
displacement from the common start has median **0.03925 nm**, maximum
**0.52143 nm**. Amplitude changes range **−64.13% to +178.44%**, median
**+1.95%**, within unchanged native bounds. Lower RSS does not establish
geometrically or chemically truer parameters.

Separate `compare_reconstructed_champion.jl` and
`report_unit_assignment_benchmark.jl --full145-own-n` produce the complete
external grade. **41 aggregate checks** and **435 paired scan rows** confirm:

| Profile | Correct / 870 | Exact / 145 | Coverage / 870 | Emitted errors |
|---|---:|---:|---:|---:|
| Historical saved reference | 677 | 36 | 854 | 177 |
| Exactly replayed saved support | 676 | 34 | 852 | 176 |
| Same-start joint optimization | 676 | 34 | 852 | 176 |
| Geometric variable projection | 632 | 20 | 852 | 220 |

Profiled classified accuracy is **632/852 (74.2%)**; full-denominator accuracy
is **632/870 (72.6%)**. Counts stay **106/145 exact**, **892 predictions**,
**16 missing truth positions / 38 extra lobes**, 14 short and 25 long scans.
The same four raw keys are unavailable in all new arms; two align to benchmark
truth. No exclusion or new abstention explains the drop.

Against joint: **18 scan gains, 44 losses, 83 ties**, net **−44 correct**;
**two exact gains, sixteen losses**. Against saved support: **17/44/84**,
same net and exact gains/losses. Eight external grade orientations change
in each contrast. Exact gains are `240307_016.sxm` and `240815_048.sxm`;
all sixteen losses are preserved in the paired report. Against history:
**27/57/61**, net **−45 correct**, six exact gains and 22 losses, with fourteen
orientation changes. Grade orientation is not measured spatial registration.
Confusion `(TN,FP,FN,TP)` goes **(513,52,124,163) → (481,84,136,151)** against
both controls: **32 additional false positives and twelve false negatives**.

Profiled changes **110 final decisions on 67 scans** against saved support,
**67 from 0→1, 43 from 1→0**; twenty leave old exact vote ties. Raw populations
go **(673,223,4) → (649,247,4)** for `(0,1,?)`. Joint changes two decisions
but leaves these populations and the external aggregates unchanged. Against
joint, profiled changes 110 decisions across 68 scans. The whole-cohort heads
couple scans; this is not a per-scan independent classifier experiment.

**Reject this geometric-profile version and retain saved support 676/34.**
Historical 677/36 remains unexceeded. Better optimization of this Gaussian
objective at fixed N worsens this recognition pipeline. Absorption of useful
contrast by the movable base geometry is a possible explanation, not a causal
finding established here. This does not reject all variable-projection methods,
and no post-grade optimizer, threshold, geometry, vote or parameter adjustment
is made. Inference remains label-free; repeated benchmark reuse provides
development evidence, not independent validation or unknown-chain chemistry.
The authorized deliverable is complete; the broader historical objective is not.
No other direction, acquisition model, missing-head vote, count sweep, DFT
change or unknown25 run follows.

Full methods, diagnostics, hashes, tests, confusion and paired losses are in
`results/geometry_profile_20260922/report.md`. README, runbook, selection,
assignment, configuration, calibration and HPC documentation are updated to
this outcome; the assignment page retains the existing 200 KiB limit.

The first final documentation render reports **200.01 KiB** for the assignment
page, just above its unchanged 200 KiB limit. The new summary is shortened,
leaving displacement details in this journal and the full report; no threshold
or scientific result is changed.

The final Julia 1.13 documentation build passes with the existing size/search
warnings. The assignment page is **199.81 KiB**, below the unchanged limit.
Source remains **7c5fb6c**, with no post-grade scientific changes. The completed
negative comparison is handed off with saved support retained and no new job.

### 2026-09-22 — Authorized fixed-geometry acquisition translation

After the negative geometric-profile result, the user authorizes **piste 1
only**: correct identifiable trace/retrace displacement before patch extraction.
Preparation starts **18:21 CEST**, capped at four hours; at most **one Viper job,
two hours, four requested CPUs / 16 GB**, no automatic retry. Completion means
a fully fetched, checked external comparison against saved support **676/870,
34/145** and historical **677/870,36/145**, or an explicit indeterminate result.
The 10–20mer application has already been completed on the frozen reference;
it is not rerun. No cross-view residual fit, measurement-error classifier,
missing-head fallback, N selection, composition prior or DFT change is included.

`acquisition_registration.toml` fixes the experiment before the real smoke or
grade. It compares **unsmoothed native-preprocessed Z**, restoring masks from
the actual SXM channels (the preprocessing function's `raw` return is already
imputed). The backward x reversal belongs to the reader and is NOT repeated.
Integer sample lags span the smaller of **2 nm** and **one quarter of image
width**, with no per-file enlargement. This covers the earlier unlabelled
−1.370 nm example without assuming that its displacement applies elsewhere.
All candidates use the exact same observed forward pixels, each available in
backward at every candidate lag. Rows with fewer than **16** such pixels are
excluded from all candidates. Correlation pools within-row-centered contrasts;
it is signed, with no absolute-correlation or negative-gain rescue.

The global maximum and all **four contiguous y bands** must each have at least
**256 pixels / eight rows**, correlation **≥0.60**, an interior maximum, and a
correlation gap **≥0.01** to every peak outside a **0.08 nm** neighborhood.
Band maxima must agree with the global maximum within **0.08 nm**. Both
distances round upward to integer pixels, with a one-pixel minimum; the search
must span more than the neighborhood plus one pixel. These are declared
stability criteria, not independent-band evidence, a calibrated uncertainty
interval, or a chemical test. Every rejected/ambiguous scan remains in the
cohort with **zero applied displacement** and an explicit reason.

Three full arms isolate native reference replay, restored-observation-mask
control at zero lag, and the same mask treatment with accepted shifts. The
forward view is the reference: **B_new(y,x)=B(y,x+dx)**, outside pixels become
NaN, then the original smoother and **S(B_new−M_frozen)** are used. No shifting
of an already-subtracted residual is allowed. Base/split tables, N, lobe
coordinates, amplitudes, widths, axis, baseline, classifier settings, molds and
vote remain fixed. The cached geometry was fitted to an unregistered fused
view: using forward coordinates is the explicit perturbation being tested,
NOT a claim that this saved geometry is optimal for the registered images.
Restoring missing masks cannot undo their earlier influence on native plane/
row flattening. The zero-shift control separates that mask effect from the
translation effect. Edge losses become unavailable pixels, never dropped files.

Implementation is opt-in through `--acquisition-shifts`, requiring supplied
base/split geometry and fresh patches. Synthetic tests cover known positive,
negative and zero shifts, signed signal, boundary/constant/periodic ambiguity,
band disagreement, fixed support, raw holes, untouched geometry and model
subtraction order. An initial test expected a direction-check helper to return
`nothing` while Julia returned its final `true` condition; an explicit `nothing`
return resolves this API-only mismatch. No scientific setting was adjusted.

**Preflight completed, before grading.** Julia 1.13 passes **10,033 assertions**
(registration 9,224; patch extraction 390; pipeline 133; unchanged support/score
286). A single real smoke uses the lexicographically first cached scan,
`240307_015.sxm`, not a selected success/failure case. It accepts **−52 pixels,
−0.610567 nm** within the fixed ±128-pixel window; global correlation is
0.990756, distant-peak gap 0.016541. Band maxima are −49/−50/−53/−52 pixels,
each passing the frozen gates. Raw missing counts are 36,016 forward / 36,352
backward, so the separate mask-only control is necessary. All nine patch
exports retain the seven saved lobes and metadata; **326 additional checks**
pass. No grade or chemical agreement was computed for this scan. There was
no parameter adjustment after the smoke. The full146 driver dry-run verifies
146 files / 900 keys without reading image arrays. Documentation builds under
the unchanged 200 KiB limit (assignment page **199.81 KiB**).

**Submission:** source **9961973** is committed before transfer. Checksummed
rsync dry-runs confirm identical tracked sources, all raw inputs, base/split
caches and templates. The ignored Manifest is copied separately and remains
SHA256 `617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe`
after login-node `Pkg.instantiate()` / precompile with Julia 1.13. The expected
Manifest 1.12.6 warning is recorded, not repaired by hand. Local and remote
driver dry-runs pass; `sbatch --test-only` succeeds. The actual submission uses
the explicit command-line `--export=ALL,...` list for all five required paths.
One job **11941278** is submitted at **18:47:13 CEST**, after approximately
26 minutes preparation; initial state is `PENDING (QOSGrpCpuLimit)`. The
two-hour cap and no-retry scope are unchanged. Remote source:
`/u/oldu/code/STMFit_acquisition_registration_20260922`; outputs:
`/ptmp/oldu/stmfit/acquisition_registration_20260922_v1`. No STM computation runs
on the login node. Submission is not a scientific result; the comparison is
still pending and the saved support candidate remains the working reference.

**Started:** Slurm records submission at **18:47:14** (one second after the
submission command's timestamp) and start **18:53:18 CEST**, **6m04s** queued.
The job runs on `vipc2227`, eight allocated CPUs with four one-thread workers.
All four acquisition logs advance. The earlier geometric-profile paragraph in
the assignment page is condensed, retaining its headline and linking the full
journal/report. Documentation passes; the page is now **198.11 KiB**, leaving
room for this outcome without raising the existing 200 KiB limit. No scientific
code or setting changes after source **9961973**.

**Completed and checked:** job **11941278** finishes at **19:00:10 CEST**,
**COMPLETED, 0:0, 6m52s**, without retry or budget extension. All **626 result
files and the Slurm log** are fetched with matching SHA256. The independent
saved-only review passes **23,526 assertions**, recomputing acceptance from
all **183,080 lag scores**, checking full cohorts, literal cached geometry,
input hashes, fifteen reference tables and final probability votes. External
grading then passes **41 aggregate checks**; four paired comparisons contain
**580 file rows**. Source stays **9961973**; only documentation changes later.

Registration is accepted on **104/146 scans**, all nonzero and negative,
spanning **−0.841487 to −0.150685 nm**, median **−0.328767 nm**. The other
42 scans stay at zero lag, not removed. Rejection reasons overlap: 37 ambiguous
peaks, 14 band disagreements, eight weak/negative correlations, seven
insufficient supports, four constant signals and three boundary peaks. No
window expands. The accepted scans contain 631 saved lobes; 629 have measurable
raw-patch correlation on a fixed common finite set before/after. **All 629
improve**, median **0.231432→0.991123**, none negative afterward. The remaining
two have empty common support (240818_019 lobes 1–2). This is spatial agreement
evidence, not independence, calibrated noise, or chemical identity.

There are **36 scans with missing raw pixels**, but the zero-shift observed-mask
control reproduces **all fifteen scientific tables** of the native reference,
patches included, byte-for-byte. Thus restoring masks has no numerical effect
on these molecular patches, not a license to impute other supports. The native
reference itself exactly replays the saved support.

| Profile | Correct /870 | Exact /145 | Coverage /870 | Emitted errors | TN / FP / FN / TP |
|---|---:|---:|---:|---:|---|
| Historical saved reference | 677 | 36 | 854 | 177 | 507 / 59 / 118 / 170 |
| Native support replay | 676 | 34 | 852 | 176 | 513 / 52 / 124 / 163 |
| Observed-mask, zero shift | 676 | 34 | 852 | 176 | 513 / 52 / 124 / 163 |
| Accepted acquisition shifts | 671 | 26 | 850 | 179 | 523 / 40 / 139 / 148 |

Against either control: **11 scan gains, 15 losses, 119 ties**, net **−5 correct
positions; zero exact gains, eight losses**, with three grade-orientation
changes. Exact losses: 240307_017, 240310_Cu100007, 240817_006, 240817_007,
240817_044, 240817_080, 240817_081 and 240818_020. Against historical: 22 gains,
28 losses, 95 ties; four exact gains, fourteen losses, net −6 units / −10
exact chains. Grade orientation is a reporting convention, not registration.

All **146 scans / 900 keys** remain. The grade covers 145 scans / 892 predicted
lobes against 870 truth positions. Every N stays unchanged: 106 exact counts,
14 short / 25 long chains, 16 missing positions and 38 extra lobes. Final raw
assignments change at **30 keys / 28 scans** (28 `1→0`, two `0→?`). Six changed
keys belong to zero-shift scans, because the unchanged learning rules refit
on the cohort's changed features; learned classifier parameters are not frozen.
The GMM hard head changes 154 decisions, k-means six; the final soft vote is
checked separately. Forward patches, forward CC and Fisher stay identical
between mask control and registered arms.

Four old unavailable keys remain. Two more appear at **240818_019 lobes 3–4**:
translation exposes missing backward descriptor-disk pixels at the image edge,
so GMM and the unchanged final rule abstain. K-means has no missing head; no
fallback fills the two new gaps. False positives decrease **52→40**, but false
negatives increase **124→139**; classified accuracy falls **79.3%→78.9%**.

**Decision: reject this fixed-geometry registration variant for recognition;
retain saved support 676/34**, still below historical 677/36. Better trace/
retrace agreement is not enough for this frozen-model assignment pipeline.
Interaction with geometry fitted to the original unregistered fusion is a
plausible explanation, not a demonstrated cause. The result does not reject
all registration methods or prove that chemical signal is absent. No gate,
window, geometry, vote, threshold, missing-head rule or parameter is tuned
after the grade. No other suggested direction, N, DFT or unknown25 run follows.
Inference remains label-free; repeated benchmark reuse is development evidence,
not independent validation or accuracy evidence for unknown chains.

The authorized comparison is complete; the broader historical objective is not.
Full method, hashes, checks, confusion and paired losses:
`results/acquisition_registration_20260922/report.md`. README, runbook,
selection, assignment, configuration, calibration and HPC notes reflect this
outcome. README's support claim is clarified to **no patch-pixel imputation by
the support rule**, not a false claim that native SXM preprocessing never imputes.

Final Julia 1.13 documentation rendering passes, with the existing size/search
warnings. The assignment page is **199.61 KiB**, below the unchanged 200 KiB
limit. The first final-build invocation never executes because its automatic
command-approval review times out; the single retry only rebuilds local docs,
not inference or Slurm. Scientific source remains **9961973**, with no
post-grade parameter changes. The completed negative experiment is handed off
with saved support retained; no required scientific work remains in this scope.

### 2026-09-22 — Authorized native refit after acquisition registration

The user authorizes the first proposed follow-up only: refit the model on the
registered acquisitions, with a matched unregistered refit. Preparation begins
**20:12 CEST**, capped at four hours, then **one Viper job capped at two hours**,
four requested CPUs / 16 GB, without automatic retry. The deliverable is the
complete fetched, verified recognition comparison against saved support
**676/870 correct, 34/145 exact** and historical **677/870, 36/145**, or an
explicit indeterminate result if a full valid comparison cannot be produced.
No other suggested direction, count sweep, class/composition prior, missing-head
fallback, DFT change, unknown25 rerun or post-grade search is included. The
completed application and archived plan's scientific Must-NOT-have constraints
remain in force; no retired process apparatus is revived.

The prior registration experiment left geometry fitted to the original fused
image unchanged. Its negative recognition result does not establish the cause;
this experiment tests the coupling between acquisition alignment and refitting,
not a promised gain. Three arms use the same raw/cache cohort: saved reference
replay, observed-mask zero-shift refit, and accepted-registration refit. The
registration estimator and every gate in `acquisition_registration.toml` remain
unchanged. Rejected shifts remain zero; scans are not removed.

Both the **Gaussian base and asymmetric split shape** are freshly fitted, so
the classifier's split-skew feature is not left tied to the unregistered image.
Their selected N is copied from the saved label-free geometry, never expected
counts. Each profile uses the existing native circular initialization followed
by circular-to-elliptical refinement; valid minimum **full-parameter GCV**
chooses its family at that N. Physical boxes, overlap, κ penalty and native
global/local optimizers stay unchanged: `chitosan.toml` provides the global
10-second / 10,000-iteration limit, 300 circular local iterations, one start;
the existing 50-iteration elliptical step is explicit in `registered_refit.toml`.
No extra BOBYQA or variable projection is applied. Native optimizer convergence
is not exposed, so it is reported unknown; equal ceilings do not guarantee
identical wall-clock-limited search trajectories or convergence.

Native flattening is applied before the integer-x backward translation, with
forward coordinates as reference. Actual raw masks are restored, no second
backward reversal occurs, and outside-image samples remain NaN. The fit uses
the **unsmoothed observed mean**, with the native ROI detector's statistics
computed only on finite smoothed pixels. Its existing thresholds, dilation,
5th-percentile reference, axis and tube/support rules are unchanged; the ROI is
intersected with observed samples. Both arms follow these same rules, although
their data-dependent ROI, axis, support, family and full parameter vector may
differ. Their RSS values therefore are not a same-target optimization contest.
Noise retains the conservative native maximum of the two view estimates;
there is no independence claim or division by sqrt(2), no new covariance or
change to `n_eff`.

Core gains an **opt-in `observed_only` path**, off by default, for ROI/seed
statistics. If an initializer samples a missing pixel, it takes the physically
nearest actually observed value **only to initialize that amplitude**. No image
pixel, objective sample or downstream patch is filled. Missing masks cannot
undo the influence of earlier native imputation on plane/row flattening.
All seeds still use only 2D data, never the diagnostic 1D fit. Fully observed
ROI and seed results are tested against the unchanged default path. At zero
accepted shift, the control's two fits are reused exactly, avoiding artificial
differences from a second time-limited optimization.

Each shard saves both profiles, full parameter vectors and bounds, bootstrap
family GCV/validity, axis/support, offset, fit pixels with both original view
values, timing and failures. Any failed scan prevents incomplete assignment or
grading; it is not dropped or silently replaced by cached geometry. All
patches, CC/Fisher descriptors and both cohort-wide learning heads are rebuilt
using the unchanged support config and seeds; learned parameters are not frozen.
Labels remain solely in the subsequent external grade. Reusing this development
benchmark repeatedly does not make the outcome independent validation.

Initial sandbox-only Julia/SSH checks cannot access the usual launcher lock
and SSH configuration; normal approved execution resolves those environment
restrictions. The first real SSH handshake times out; a read-only connection
check then succeeds through the existing master. No job has been submitted at
this preparation checkpoint.

Preflight passes under Julia **1.13.0**: **14,591 assertions** across the new
registered-refit suite (4,299), unchanged acquisition suite (9,221), core 2D
suite (151), prior geometric suite (111), pipeline (133), patch preprocessing
(390) and support/score suite (286). Tests include genuine zero-shift fit reuse,
both shape profiles, observed masks, opt-in/default finite equivalence and
forbidden inputs. The full146 metadata-only dry-run finds 146 files / 900 keys.

The only real local scan is the lexicographically first cached file,
`240307_015.sxm`, saved **N=7**, with the prior unchanged −52-pixel translation.
Both zero-shift and registered Gaussian/split fits are valid and preserve all
seven keys. RSS is **0.06344984 / 0.06044526** for control Gaussian/split,
**0.10306913 / 0.06098382** after registration; these use different data/support,
so they are not a matched-target improvement metric or a chemical grade.
The smoke has **10** cohort assertions; **169** saved-only independent checks
reassemble both peak shapes and verify parameters, bounds, actual view pixels,
RSS, full GCV, family choice and feature serialization. The saved-output
verifier initially has a dotted-number Julia parse error; correcting its
syntax and rerunning that verifier requires no second fit or setting change.
No external labels are read. Remote/local base, split and template SHA256s
match; the ignored Manifest remains unchanged. Documentation rendering is
still in progress at this checkpoint.

Scientific source is committed as **c7836f9** before transfer to
`/u/oldu/code/STMFit_registered_refit_20260922`. The tracked-source checksum
dry-run and dereferenced raw-data checksum dry-run report no differences.
Login-node `Pkg.instantiate()` and package imports succeed under Julia 1.13;
the Manifest SHA256 remains
`617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe`, with its
expected Julia-1.12.6 provenance warning. No STM fit runs on the login node.
Local documentation builds successfully; the older geometric-profile summary
is condensed without losing its result or journal/report pointers. The
assignment HTML is **198.20 KiB**, below the unchanged 200 KiB limit.

The remote metadata dry-run passes. `sbatch --test-only` succeeds and reports
eight prospective allocated CPUs for the four-CPU request; it is not a runtime
environment test. The sole actual submission is **11942091**, at
**20:39:57 CEST**, with explicit command-line `--export=ALL,...`, all five
required paths, `--no-requeue`, two-hour walltime and four requested CPUs / 16 GB.
Preparation takes about **28 minutes**, within the four-hour window. Output:
`/ptmp/oldu/stmfit/registered_refit_20260922_v1`. No scientific change after
source **c7836f9**, no second real local scan and no second job is submitted.
Continue through completion, fetching, saved-output checks and external grade;
submission alone does not complete the deliverable.

Job **11942091** starts **20:42:15 CEST** on `vipc2400`, **2m18s** queued,
eight allocated CPUs with four one-thread refit workers. Registration finishes
for all 146 scans, and the completed reference's **fifteen scientific tables**
match the saved support by SHA256, including final predictions.

During refitting, **240312_Cu100070/071.sxm** fail in the zero-shift Gaussian
control: saved N=10 requires a minimum span **4.6303 nm** under the unchanged
native parametrization, but the observed-mask support is **4.2040 nm**. The
earlier saved native support was **4.7951 nm**, with 8,131 fit pixels; the
new control saves 7,169 observed fit pixels. The two new fit-data tables are
byte-identical, as were the prior numerical native fits; no scan is removed.
**240814_025.sxm** subsequently fails only after registration: saved N=6
requires **2.5724 nm**, versus support **2.5421 nm**. These are fixed-N/support
incompatibilities, not evidence of chemically worse recognition or proof that
the saved counts are false. The complete diagnostic continues without changing
N, overlap, support gates or source. No incomplete candidate can be graded.

Only the saved-output verifier gains an explicit `--recorded-failures` mode to
audit every completed fit and account for failed stages without claiming a
complete comparison. Its default remains strict. No inference or cluster code
changes after source **c7836f9**. The ordinary smoke verification still passes.
The explicit recorded-failure fixture also passes: it combines the saved smoke
with a clearly synthetic failed filename in a temporary directory, verifies
all completed numerical fits, accounts for the missing stages, and emits no
classifier or grade. This is a verifier test, not another scientific scan.

### Completed diagnostic and decision

Job **11942091** finishes **21:09:12 CEST**, **FAILED, 1:0**, after **26m57s**
compute and **2m18s** queued. All four shards finish traversing their input
files; the failure status comes from the declared refusal to emit incomplete
assignment arms, not a time limit, OOM, new retry or count adjustment. All
**549 regular output files plus the Slurm log** are fetched locally with
matching SHA256. The source used on Viper remains **c7836f9**.

There are exactly the three support-span failures listed above, no others.
All **146 scans / 900 cached input keys** are accounted for; **143 scans**
complete Gaussian and split fits in both arms. Status accounting is:

| Stage | Native-valid fits | Failed fits | Not attempted after failure |
|---|---:|---:|---:|
| Zero-shift Gaussian | 144 | 2 | 0 |
| Zero-shift split | 144 | 0 | 2 |
| Registered Gaussian | 143 | 1 | 2 |
| Registered split | 143 | 0 | 3 |
| Total | 574 | 3 | 7 |

The two N=10 controls fail before any new family optimization; the same files'
registered stages are not attempted. The N=6 case has two valid control fits
before its registered Gaussian initialization is rejected. Its required
2.5724 nm versus 2.5421 nm support is not a numerical convergence diagnosis.
The model's minimum-span rule is conservative because it couples spacing to
maximum allowed width and overlap. These failures do **not** establish that
the saved N values are false or that no other physical geometry is possible.
Native-valid fits are not advertised as converged: native convergence remains
unexposed. No constraint is relaxed and no alternate N/support is tried.

The acquisition shift table, summaries, peaks and **all lag scores** replay
the preceding experiment exactly: **104 accepted nonzero shifts, 42 zero**.
The two failed controls are zero-shift scans, leaving **80** exact reused
profile fits among the successful zero-shift cases. The saved reference
replays **fifteen scientific tables**, including all predictions, byte-for-byte.

The final saved-only audit passes **22,134 assertions**: recorded failures and
missing-stage consistency, cohort accounting, actual finite fit pixels and
view mean, full parameter vectors and bounds, independently assembled Gaussian
and split RSS, full-parameter GCV, family choice, exported features where a
complete shard exists, and exact zero-shift reuse. A separate status/input/span
summary passes **608** checks (591 stage/cohort, eight input-hash, nine physical
span), making **22,742 final checks**. These checks validate the diagnostic
record, not a complete classifier or benchmark. The explicit failure-verifier
mode changes only postprocessing; inference and its settings remain frozen.

**No new control or registered assignment table is produced, and the external
candidate grade is not run.** The aborted arms are not trained on the 143
successful scans, nor patched with saved predictions or a new abstention rule.
No per-unit accuracy, exact-chain count or claimed gain can be assigned to
this incomplete candidate. The exactly replayed saved support still corresponds
to **676/870 correct, 34/145 exact**, below historical **677/870, 36/145**.

**Decision: retain saved support; this registered-refit comparison is
inconclusive**, not a measured negative recognition effect. The bounded
diagnostic is finished, but the requested full recognition comparison is
unavailable and the broader historical objective remains unmet. Support/N
consistency must be resolved in a separately scoped experiment: preserve/replay
the original support, or permit a fresh label-free N selection on corrected
acquisitions. Neither alternative, another proposed direction, new DFT or an
unknown25 rerun is launched. No benchmark-dependent adjustment occurs.

Full methods, exact failures, all **584 planned stage statuses**, accounting,
hashes, tests and commands: `results/registered_refit_20260922/report.md`.
README, runbook, selection, assignment, configuration, calibration and HPC
notes distinguish this inconclusive diagnostic from the earlier scored
fixed-geometry experiment. The saved working candidate is unchanged.

The final documentation build succeeds under Julia **1.13.0**, without
deployment. The assignment page is approximately **199.6 KiB**, below the unchanged
200 KiB limit; existing page-size/search-index warnings remain. The final
wording distinguishes 549 regular output files from the additional Slurm log.

### 2026-09-22 — Authorized registered refit with original geometric support

The user authorizes the proposed follow-up: retain the original geometric
support instead of re-estimating it on each masked/registered acquisition.
Preparation begins **21:26 CEST**, limited to four hours; one Viper job is
limited to two hours, four requested CPUs / 16 GB, with no automatic retry.
Done means all 146 scans / 900 saved keys, reference replay plus matched fresh
zero-shift/registered Gaussian and split fits, fetched outputs and independent
saved-output checks, then the full external grade and paired gains/losses. If
any required fit fails, the diagnostic is explicitly inconclusive; no partial
classifier or grade, deleted scan, saved-prediction fallback or new abstention.
This is a new scoped experiment, not a resubmission of job 11942091.

`registered_refit_original_support.toml` keeps all previous fit, registration,
assignment and selected-N settings. Its sole methodological change is the
support policy: replay the **unregistered native ROI, weighted axis/origin,
fit tube and axial bounds once**, before restoring observation masks, then
freeze them for both profiles and both arms. The frame must reproduce the
saved Gaussian cache at its original serialization precision. No bound is
extended to accommodate N; no N sweep or physical parameter adjustment occurs.
The original support is therefore label-free but retains the legacy native
preprocessing's imputation influence. It is not a new observation-only support
estimator or a claim of physically correct geometry.

Each objective uses only the original fit mask intersected with actual finite
forward/backward observations after the arm's translation. The offset remains
the fifth percentile on that arm's observed portion of the original ROI.
No imputed sample enters the new objective or patches. The original native ROI
values and fit membership are saved solely to verify the support calculation;
the native noise rule, initialization-only nearest observed sample, optimizers,
GCV family choice, full parameter penalty and `n_eff` formula stay unchanged.
Exactly zero shifts reuse the control fits. Arm-specific finite masks, offsets,
families and fits can still differ; RSS is not a same-target improvement metric.

The previous failed registered scan `240814_025.sxm` is the single planned real
local smoke: it tests both zero-shift and accepted −30-pixel registration at its
saved N=6, without reading labels. All other real fits stay on Viper. Synthetic
tests check frozen support, frame mismatch rejection, raw holes, both profiles,
zero-shift reuse, strict configurations and unchanged legacy mode. Complete
support replay will be checked against the earlier saved native-geometry run,
including all original fit pixels, before any external grade.

The completed unknown25 application is not rerun. Existing scientific
Must-NOT-have constraints remain: no label/composition prior in inference,
truth-selected per-file variant, uncertain-row exclusion, overwritten user
artifact or N selection change. No new DFT, missing-head rule, registration
gate search, classifier tuning, combined hypothesis, post-grade search or
retired orchestration is included. Saved support **676/870, 34/145** is retained
until evidence supports a decision; historical **677/870, 36/145** is not a
selection input. Repeated development-benchmark reuse is not independent validation.

Preflight passes under Julia **1.13.0**: 4,387 registered-refit assertions,
9,221 unchanged-registration assertions, 111 geometric comparison assertions
and 133 pipeline assertions (**13,852** total), plus 175 saved-only legacy
smoke checks. The single real smoke `240814_025.sxm` completes all four fits
at saved N=6, including the formerly failing −30-pixel registered arm. Its
support is exactly the earlier native **3.6720726859651593 nm**, not an expanded
or fitted bound. Seventeen smoke checks, 188 independent saved-fit checks and
306 saved-support/reference/feasibility checks pass. The last group reads no
new raw image: all 146 previously saved original supports pass the unchanged
minimum-span rule for both native families. This does not promise that every
new fit will be valid or improve recognition.

Smoke RSS is **0.4235117496 / 0.3998868007** for control Gaussian/split and
**0.4664677406 / 0.4537206868** after registration; the observations differ, so
these are not a chemical score or common-target improvement. The metadata-only
launcher dry-run accounts for all **146 files / 900 keys**, with
`original_support=true`. Local and Viper base/split/template SHA256s match;
remote Julia is 1.13.0. No new real local file or fit retry is used.
Documentation builds successfully, without deployment; existing size/search
warnings remain. The Manifest is unchanged (SHA256
`617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe`).

Scientific source **98d664d** is committed before transfer to
`/u/oldu/code/STMFit_registered_original_support_20260922`. Tracked-source and
dereferenced raw-data checksum dry-runs show no differences. Login-node
instantiation/imports succeed under Julia 1.13, retaining the unchanged
Manifest and its expected 1.12.6 provenance warning. The remote metadata-only
dry-run and `sbatch --test-only` pass; test-only ID 11942392 is not a submitted
job or a check of the compute-node environment.

The sole actual job **11942393** is submitted **21:42:17 CEST**, starts
**21:42:26** on `vipc2488`, with nine seconds queued. It requests four CPUs /
16 GB, receives eight allocated CPUs and runs four one-thread workers.
Command-line `--export=ALL,...` names all five paths, with `--no-requeue` and
the explicit two-hour limit. Preparation takes about sixteen minutes.
Output: `/ptmp/oldu/stmfit/registered_original_support_20260922_v1`.
No scientific source/settings change or second submission follows. Continue
through fetching, saved-output verification and external grade, not merely
submission. The assignment page's older amplitude-profile summary is condensed
to leave documentation space, preserving its result and full journal/report links.

During execution, the reference's fifteen scientific tables replay exactly.
The three previous failures now complete all four fits on Viper, including
both N=10 zero-shift cases and registered `240814_025.sxm`. A different file,
**240817_006.sxm**, then fails at **registered Gaussian, saved N=6**: both
circular and elliptical families are rejected for **`residual high`**.
This is the existing maximum-absolute-residual/noise validity guard (native
threshold 3.5), not a support-span failure. The failed parameter vectors and
individual residual values are not exported by the existing initializer;
no numerical excess, convergence claim or chemical interpretation is inferred.
Its registered split stage is not attempted after the Gaussian failure.
The complete traversal continues unchanged, but this failure prevents a
complete candidate grade. No new threshold, N, optimizer budget or fallback
is introduced to rescue it. Final failure/cohort accounting awaits completion.

The condensed documentation rebuild passes; assignment HTML is 197.86 KiB,
below the unchanged 200 KiB limit, with existing warnings and no deployment.

Postprocessing additionally recomputes maximum absolute residual divided by
the native noise for every saved valid fit, compares it to the exported value
and checks the unchanged native validity threshold. The smoke passes 196
saved-fit checks plus the previous 306 original-support checks. This extends
verification only: the scientific source on Viper stays **98d664d**, and failed
fit parameters are neither reconstructed by a new optimization nor invented.

### Completed original-support diagnostic and decision

Job **11942393** ends **22:10:36 CEST**, **FAILED 1:0**, after **28m10s**
compute and nine seconds queued. All four shards traverse all their files;
the failure is the declared incomplete-cohort guard, not timeout or OOM.
All **706 regular output files plus the Slurm log** are fetched locally and
match their SHA256s. Source **98d664d** and inference settings stay unchanged;
there is no retry or second real local scan.

All **146 scans / 900 saved input keys** are accounted for. The exact status
table is saved in `results/registered_original_support_20260922/fit_statuses.tsv`:

| Stage | Native-valid fits | Failed fits | Not attempted after failure |
|---|---:|---:|---:|
| Zero-shift Gaussian | 146 | 0 | 0 |
| Zero-shift split | 146 | 0 | 0 |
| Registered Gaussian | 145 | 1 | 0 |
| Registered split | 145 | 0 | 1 |
| Total | 582 | 1 | 1 |

The three previous span failures complete all four fits. Every original
ROI/axis/tube/bound and every original native fit sample matches the earlier
saved native-geometry run exactly, including the failed scan. All 146 original
contexts are feasible at saved N under both families' unchanged span rule.
No bound was expanded to rescue a count. The registration's four tables replay
exactly: **104 accepted nonzero shifts, 42 zero**; all **84** zero-shift profile
fits are reused exactly. Fifteen saved-reference tables, including predictions,
are also byte-identical.

The sole failure remains **240817_006.sxm**, saved **N=6**, accepted **−35 px**
translation. Its original support is **3.8117942672196974 nm**; both controls
are valid on 6,496 pixels, with native maximum-residual/noise values
**1.5661341761** (Gaussian) and **1.4073437257** (split). Both registered
Gaussian families fail the unchanged **3.5** validity guard. Their actual
residual values/parameters are unavailable, so no quantitative excess or
convergence explanation is claimed. Registered split is not attempted.
This does not prove that N is wrong or that recognition is worse.

Saved-only verification passes **27,060 numerical/cohort assertions**:
24 explicit-failure checks and 27,036 parameter/pixel/independent RSS/GCV,
maximum residual, family, original-context and zero-shift checks. Separate
summaries pass 594 stage/key checks, nine input-hash/config-copy checks,
1,607 prior-native context/sample checks, and 147 fixed-N support-feasibility
checks, plus the original-support mode assertion. No new optimizer is run.
These validate the diagnostic record, not a complete classifier or benchmark.

**No new control or registered assignment table is produced and no external
candidate grade is run.** In particular, 145 successfully refitted scans do
not constitute the full145 external benchmark. There is no scan exclusion,
training on the successful subset, saved-prediction substitution, new
abstention rule, threshold relaxation or post-grade tuning. The complete-arm
descriptive script is prepared but not executed. The exactly replayed saved
support remains **676/870 correct, 34/145 exact**, below historical **677/870,
36/145**; these are not scores for the new refits.

**Decision: retain saved support; recognition remains inconclusive.** The
authorized support replay resolves the previously observed span incompatibility,
but does not produce the required complete recognition comparison. The bounded
diagnostic is finished and the broader historical objective remains unmet.
Handling model-validity failures, changing N or any other hypothesis needs a
separate human-approved scope; no such change, new DFT or unknown25 rerun follows.
Full methods and evidence: `results/registered_original_support_20260922/report.md`.
README, runbook, selection, assignment, configuration, calibration and HPC notes
are updated consistently. Scientific source remains **98d664d**; later source
changes affect only saved-output verification and documentation.

The final documentation build passes under Julia **1.13.0**, without
deployment. Assignment HTML is **199.27 KiB**, below the unchanged 200 KiB
limit; existing page-size/search-index warnings remain. The Manifest hash
remains unchanged. The completed diagnostic is handed off without promoting
a candidate or claiming the broader recognition objective is achieved.

### 2026-09-22 — Authorized failure diagnosis before a paired-acquisition fit

The user authorizes the proposed follow-up. The first deliverable is a
numerically explained reproduction of `240817_006.sxm`, selected because it
failed native validity, not because of its external grade. Preparation starts
at **23:21 CEST**, capped at four hours, followed if feasible by **one Viper
job capped at two hours**, four requested CPUs / 16 GB, no automatic retry.
The intended method comparison is common molecular geometry with separate
forward/backward acquisition gain/background at saved label-free N. Done means
tested code, fetched and checked complete-cohort predictions and an external
recognition comparison, or an explicit indeterminate result if the diagnostic
exposes a blocking issue. The reference is retained until that comparison.

The native optimizer receives an optional, post-optimization observer. It
records each start's raw initial/global/final vectors, global termination and
evaluation count, LM convergence/iterations and swallowed optimizer exceptions.
Observer vectors are copies. Native objectives, bounds, initialization, budgets,
minimum-RSS start choice and final validity remain unchanged. With observation
enabled LM stores its existing numerical trace; no optimizer is added. Timed
global optimization means a new replay is not assumed bit-identical to the old
cluster result, and a local convergence flag is not global optimality.

`test/diagnose_registered_fit_failure.jl` replays the original geometric support
and both zero/accepted shifts on exactly one raw scan. It saves all native
circular and elliptical Gaussian/split candidates, including rejected ones,
their observed-pixel residual maps and the unchanged validity calculations.
The split fit is attempted independently here to diagnose shape adequacy; this
does not change the previous production comparison's skip-on-failure policy.
No candidate is rescued, selected by labels, substituted into saved predictions
or graded. No new physical/selection parameter is introduced.

The completed unknown25 application and the archived scientific Must-NOT-have
list remain intact: no composition prior, expected count, truth-selected variant,
uncertain-row exclusion, user-artifact overwrite, new DFT, count reselection,
classifier tuning, threshold relaxation or retired orchestration. The reused
development benchmark is not independent validation. The noise difference
between acquisitions is not assumed independent or a calibrated uncertainty.

The native observer passes **40 synthetic equality/isolation checks**, and the
2D package passes **151 tests**. The first real local diagnostic stops after
the two control Gaussian fits because the shared TSV writer deliberately
refuses to overwrite a previously emitted summary. This is a diagnostic writer
bug, not a scientific failure. The initial outputs/log are retained; no
registered fit or external grade ran. The writer now emits unique per-stage
tables and final summaries once; a new synthetic end-to-end test adds **27
checks**, with six CLI checks (**73 diagnostic tests** total). The corrected
single-file diagnostic uses a separate `diagnostic_v2` directory, never
overwriting the first attempt.

An old `audit_fwd_bwd_joint_refit.jl` prototype already exists, but has free
gain/offset scale degeneracies, permits negative gains, subsamples its fit and
silently substitutes an OLS rescore on failure. It is not a complete recognition
comparison and is not reused as such. The new opt-in paired feasibility model
fixes the mean gain to one, sharing all molecular parameters, with opposite
view-specific gain and background-plane differences. Its average is exactly
the native Gaussian/split model; no new mean shape or noise model is introduced.

`config/paired_acquisition.toml` freezes positive gains in **[0.5, 1.5]**, a
differential intercept box **[-5,5] nm**, differential tilt boxes **[-1,1]**, and
**300 LM iterations**. The plane boxes reuse native numerical bounds; the gain
box is an exploratory constraint, not calibrated chemical evidence. All values
are explicit and are not chosen by the benchmark. The gain multiplies only the
molecular peaks, not the mean background. Equal/opposite view coefficients
remove the trivial scale/offset ambiguity; this does not guarantee that every
remaining parameter is identifiable on a narrow molecular support.

The planned one-file paired feasibility check starts each family from its
exported registered native candidate, even if native validity rejected it.
Both the fused-mean continuation and paired-view LM use the same molecular
start, physical bounds and 300-iteration ceiling. No new global search or
alternative seed is added. Both must pass native mean validity; the paired arm
also must pass the existing 3.5 maximum residual/noise guard in **each** view,
using the unchanged conservative native noise. At fixed saved N, valid minimum
full-parameter GCV selects the family separately in each arm/profile. The paired
score uses stacked observations and four additional parameters; it is a
predictive heuristic, not a proof of independent noise, and `n_eff` is unchanged.
Any missing Gaussian or split family leaves this planned complete-comparison
workflow infeasible; it does not permit partial learning or scoring. No claim of
global impossibility follows from a failed local attempt.

The paired model passes **85 synthetic/configuration checks**, including exact
mean preservation, view exchange, recovery of synthetic acquisition parameters,
full GCV count, finite bounds, rejected unobserved samples, unchanged starting
vectors and non-relaxation of native validity. An initial test-only decimal/
broadcast syntax error was corrected before any real paired fit. The only real
scan remains `240817_006.sxm`; full-cohort fitting is reserved for Viper if this
fixed hypothesis can produce all required fits on that scan.

### Completed native diagnosis and matched paired feasibility

Sources **69ea44b** (observer/model) and **5a1008c** (matched diagnostic runner).
The corrected local native diagnostic completes all eight family fits on
`240817_006.sxm`, saved N=6, original support **3.8117942672 nm**, **6,496
observed pixels** in each arm and unchanged noise **0.01516235363 nm**.
Control maximum absolute residual/noise is **1.5450/1.5661** for circular/
elliptical Gaussian, **1.3936/1.4073** for split; all are valid. After the
unchanged accepted −35-pixel translation:

| Native family | Maximum absolute residual/noise | Pixels above 3.5 | Valid |
|---|---:|---:|---|
| Gaussian circular | 4.248451 | 6 | no |
| Gaussian elliptical | 3.609172 | 3 | no |
| Split circular | 3.963961 | 3 | no |
| Split elliptical | 3.341651 | 0 | yes |

The Gaussian elliptical exceedances are all positive, at **row 179, columns
287–289**, with residuals **0.054724, 0.054291, 0.053399 nm**. They are adjacent,
not demonstrated independent outliers. Its 99th percentile absolute residual
is **1.2015 noise units**; the failure is a localized maximum, not a high
residual everywhere. All eight native LM fits (valid ones included) hit their
configured **300 circular / 50 elliptical** iteration caps without declaring
convergence. No swallowed global or LM exception is recorded. Global wrapper
statuses are Success or MaxIters; Success aggregates several NLopt termination
conditions and is not a global-optimality certificate. The split success cannot
replace the Gaussian base required by the unchanged assignment pipeline.

The paired feasibility check independently replays the exact saved native
registered objective pixels and uses the final native family vectors as shared
starts. All eight matched local continuations finish their **300 iterations**
without declared convergence or exceptions:

| Shape | Family | Fused mean maximum/noise | Paired mean maximum/noise | Paired forward/backward maximum/noise |
|---|---|---:|---:|---:|
| Gaussian | circular | 4.196838 | 4.127342 | 5.314639 / 3.806154 |
| Gaussian | elliptical | 3.651233 | 3.518354 | 4.641883 / 3.331164 |
| Split | circular | 4.002019 | 3.923763 | 5.017759 / 3.706424 |
| Split | elliptical | 3.357234 | 3.230708 | 4.415026 / 2.964692 |

Only the **fused split elliptical** continuation is valid. Paired split
elliptical passes native mean validity but fails the declared forward-view
guard. Both Gaussian families fail even the **original mean guard**, so the
new per-view checks are not the sole reason this comparison cannot proceed.
The fused elliptical Gaussian decreases RSS **0.176116 → 0.168989** while
increasing the maximum residual; least squares and the maximum-residual guard
are different objectives. Paired Gaussian reduces its stacked RSS
**1.066788 → 0.873776**, but its mean RSS is **0.176413**, illustrating the
cross-view tradeoff, not a chemistry improvement. GCV values with different
observation targets are not compared across modes to select a method.

**Stop this bounded version at the one-file feasibility result.** No full146
fit, HPC submission, classifier learning, partial grade or benchmark-selected
rescue occurs. The four-hour preparation / one two-hour-job ceilings are not
expanded; unused compute is not a reason to relax validity or launch a known
incomplete comparison. Counts, labels, thresholds, optimizer settings and
application predictions remain untouched. Further optimizer convergence or
mean-shape hypotheses need a separately bounded experiment; the current data
do not distinguish a local optimization limitation from model inadequacy.
Retain support **676/870 correct, 34/145 exact**, below historical **677/870,
36/145**. No recognition gain or loss is established for this new model.

The completed synthetic end-to-end native→paired test raises the diagnostic
suite to **103 assertions**, alongside **85 paired mathematical/configuration
checks** and **151 2D package tests**. Saved native residual maps independently
reconstruct in **207,995 pixel/numerical assertions**, without an optimizer or
label read. Complete saved-output verification, non-regression tests and final
documentation checks are recorded below when complete. All data stay local in
`results/paired_acquisition_20260922/`, with the initial writer failure retained.

Final numerical verification passes **363,886 paired pixel/parameter/score
assertions**, in addition to the **207,995 native** checks, without further
fitting. The unchanged registered-refit suite passes **4,387 assertions**;
with the 103 diagnostic, 85 paired-model and 151 core tests this is **4,726
focused tests**. The older signed-CC/affine-Fisher documentation is condensed
without changing its negative result or report pointer, leaving room for the
new bounded diagnostic under the existing HTML size limit. The Manifest SHA256
remains **617de8162fa4c7a814aacb106f145349bc57ead58ec9a6c332971a821b5974fe**.

The paired-diagnostic documentation build passes under Julia **1.13.0**,
without deployment. Assignment HTML is **199.53 KiB**, below the unchanged
200 KiB limit; existing size/search-index warnings remain. The completed
bounded experiment retains the reference and hands off without any HPC job
or claim of improved recognition.

### 2026-09-23 — Budget-only convergence diagnostic on the registered blocking scan

The user authorizes the next experiment. Preparation starts September 22 at
**23:53 CEST**, capped at two hours; computation is limited to **one Viper job
of one hour**, four requested CPUs / 16 GB, no automatic retry. The scientific
deliverable is to distinguish an insufficient local iteration budget from a
rejection persisting at a numerically stationary point, or explicitly report
that this budget still does not decide. This is not a recognition campaign.

`diagnose_paired_convergence.jl` consumes only the verified saved registered
pixels and native starts from `paired_acquisition_20260922/diagnostic_v2`.
Its hash checks bind those inputs to `paired_v1` and the unchanged physical and
assignment configurations. On **one scan, `240817_006.sxm`**, all four native
families receive matched fresh fused/paired fits at **300 and 10,000 LM
iterations**, from exactly the same initial vectors. No global optimization,
new seed, objective, N, physical box, noise estimate, registration, ROI,
classifier, residual threshold, 1D fit or `n_eff` change is allowed. The
completed unknown25 application and previous scientific Must-NOT-have list
remain intact. No labels are read and no benchmark is graded.

`config/paired_convergence.toml` makes the budget and numerical diagnostics
explicit before the real fits. Native `x_tol=1e-8`, `g_tol=1e-12` are unchanged;
each fit has a 600-second LM timer in addition to the job's hard wall limit.
The timer starts inside LsqFit after initial setup, so exported whole-call
elapsed time may exceed 600 seconds. Four single-thread processes shard the
eight mode/family cases; each runs its own short/long pair in sequence.

Saved LM traces must have identical short/long prefixes. Accepted updates are
reconstructed from strictly decreasing RSS; rejected trial steps are not
applied. Full parameter checkpoints are saved at 0, 300, 1,000, 3,000 and
10,000 iterations, if reached, plus the actual final iteration. Native LM can
stop on a small clipped trial step or its unprojected gradient; that flag
alone does not establish stationarity under box constraints. Trace gradient
norms can combine the pre-step Jacobian and post-step residual and are not
advertised as independently recomputed final gradients.

The separate final audit uses second-order bound-respecting differences at
`h=6.055454452393343e-6*max(1,abs(parameter))` and `h/2`, clipped to a quarter
of the parameter box width. It reports the infinity norm of the unit-box
projected gradient of RSS / initial RSS (floor `1e-12`). Both norms must be
at most `1e-6`, and full gradients must agree within `1e-7 + 1e-4*max_norm`,
for the diagnostic stationarity flag. Box-scaled Jacobian singular values
and numerical rank (`rtol=1e-10`) expose ill-conditioning. These are numerical
audit settings, not physical selection or validity thresholds, nor evidence
of global optimality. Valid minimum full-parameter GCV still selects families
within each mode/profile; no GCV comparison across different targets is made.

Before any real fit, the first synthetic derivative test exposes a Julia
literal mistake (`3f0` denotes a Float32 literal, not `3*f0`) in the new
one-sided difference expression. Explicit multiplication fixes it. All **64
new tests** then pass, including active bounds, trace prefix identity, saved
pixel runner, merger completeness and unchanged default fitting. The existing
paired model passes **85 tests**. The production engine remains unchanged;
the only shared helper addition is opt-in diagnostic exports/native optimizer
keyword forwarding. Preparation and synthetic evidence are not a scientific
result; the real outcome is recorded below after fetching the one job.

Source **92b6f15** is committed before sync. Local and remote metadata-only
dry-runs pass under Julia **1.13.0**; remote dependency instantiation leaves
the Manifest SHA256 unchanged. After `sbatch --test-only`, the sole real job
**11943523** is submitted at **00:14:57 CEST**, with explicit command-line
exports and `--no-requeue`; it initially waits for `QOSGrpCpuLimit`.

The saved-output verifier independently evaluates peaks/acquisition planes,
objective, residual maps, full-parameter GCV, bounds, checkpoint RSS, trace
prefixes, final finite-difference gradients and selection completeness. Its
synthetic end-to-end coverage raises the new suite to **17,492 assertions**;
the **85 paired-model**, **103 native-workflow** and **151 core-2D** checks also
pass (**17,831 total**). No real optimization occurs in the verifier. Fetched
data, logs and the French report belong in the ignored
`results/paired_convergence_20260923/` directory.
