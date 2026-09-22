# Chitosan Benchmark and 10–20mer Runbook

Latest September 22 comparison: **fixed-geometry acquisition registration is
rejected**. It gives **671/870 correct, 26/145 exact**, versus both native
replay and observation-mask control **676/34**; coverage falls **852→850/870**.
Registration passes on 104/146 scans and improves all 629 comparable local
image correlations (median 0.23→0.99), yet loses eight exact chains and gains
none. All N and geometric fields stay fixed. Keep saved support, below
historical 677/36; no promotion or post-grade tuning. Source **9961973**, job
**11941278**, **0:0 in 6m52s**; all **626 outputs and Slurm log** are fetched
and checksummed. Full report: `results/acquisition_registration_20260922/report.md`.
No new Gaussian fit, count selection, missing-head fallback, DFT or unknown25
run follows this completed experiment.

Earlier September 22 comparison: **geometric variable projection is rejected**.
At fixed cached N and unchanged support-assignment settings, it gives
**632/870 correct, 20/145 exact**, versus **676/34** for both the same-start
joint optimizer and byte-identical saved support replay. Coverage remains
**852/870**. Against joint: 18 scan gains, 44 losses, 83 ties; two exact
chains gained, sixteen lost. Better pixel fit is not better recognition here.
Keep the **saved 676/34 support candidate**, below historical 677/36; the
new geometric outputs are not a replacement cache. Source **7c5fb6c**, job
**11940220**, **0:0 in 15m17s**; all **782 outputs and Slurm log** are fetched
and checksummed. Full comparison: `results/geometry_profile_20260922/report.md`.
The bounded experiment is complete. No retry, new method, production promotion,
counting or unknown25 rerun is authorized by this result.

September 20 research branch: the completed symmetric-fusion comparison improves
the regenerated control from **666/870 correct and 24/145 exact chains** to
**671/870 and 28/145**, at the same cached label-free counts and coverage
(666/849 → 671/849 classified; 78.4% → 79.0%). The historical reference remains
677/870 and 36/145; the correction is not promoted. Both geometries were
regenerated in Viper job **11878491** (34m25s, exit 0:0), with control source
`f2191f7` and corrected source `d12283a`. This comparison does not rerun count
selection or unknown25. The counting results below remain pre-correction
results. See the September 20 journal entry and `hpc/compare_fusion.sbatch`.

The follow-up physical-u Fisher mirror comparison (job 11879515, 2m03s) is
negative: **669/870 correct, 27/145 exact**, with unchanged 849/870 coverage,
against the exactly reproduced symmetric-fusion 671/28. Keep the original
Fisher layout in the working candidate; the opt-in transverse config is a
recorded negative experiment. The historical 677/36 target remains unmet.

The later matched-residual comparison (job **11883607**, 5m15s, exit 0:0)
improves the reproduced control **671/28 → 672/29** (correct positions / exact
chains), with unchanged **849/870 coverage** and **79.2% classified accuracy**.
`config/unit_assignment_matched_residual.toml` became the opt-in working
candidate at that stage, not a production promotion: historical 677/36 led by five
positions and seven exact chains. Use this explicit assignment config to
reproduce the candidate; it rejects cached patches. Counting, unknown25,
thresholds and all other settings remain unchanged. The full paired report is
`results/matched_residual_comparison_20260920/report.md`; see the journal and
`hpc/compare_matched_residual.sbatch`. Do not resubmit the completed comparison.

The subsequent two-descriptor job **11888011** completes in **7m19s**, exit 0:0.
Its control reproduces all twelve saved matched-residual tables byte-for-byte.
The first moment gives **672/28** and is rejected; affine-residual half-plane
asymmetry gives **673/870 correct, 673/849 (79.3%), 33/145 exact**, with unchanged
counts and coverage. Use **`config/unit_assignment_affine_residual.toml`** for
this new opt-in working candidate; it also rejects cached patches. It is not
promoted: historical 677/36 remains ahead by four positions and three chains.
Eight scans gain, seven lose; five exact chains are gained and one lost. All
sixteen changed decisions have zero vote margin (one is outside the benchmark),
so do not interpret the gain as calibrated confidence or tune the tie rule.
All 585 outputs and the Slurm log are fetched. The complete report is
`results/transverse_descriptor_comparison_20260920/report.md`; the source is
`252edaa` and the runner is `hpc/compare_transverse_descriptors.sbatch`.
No counting or unknown25 rerun follows; do not resubmit this completed job.

The September 21 comparison of signed CC margins and affine Fisher (job
**11889205**, **7m05s**, exit 0:0; source **3c45c7c**) is negative for both
independent arms: **668/870 correct, 668/849 (78.7%), 32/145 exact** for signed
margins; **672/870, 672/849 (79.2%), 32/145** for affine Fisher. All twelve
control tables exactly reproduce the 673/33 affine-descriptor candidate. The
three arms retain all 146 scans / 900 keys and the same seven unavailable rows;
external grading alone selects the unchanged full145 manifest. Signed margins
gain on four scans and lose on nine; Fisher gains on nine and loses on ten.
Keep `config/unit_assignment_affine_residual.toml` at **673/33**, without
promotion or post-grade tuning. Historical 677/36 remains ahead. All 585 outputs
and the Slurm log are fetched and checksummed; the full report is
`results/label_free_signals_20260920/report.md`. The runner is
`hpc/compare_label_free_signals.sbatch`; this completed job must not be resubmitted
automatically. Counting and unknown25 remain unchanged.

The bounded numerical comparison (job **11891104**, **9m04s**, exit 0:0;
source **0d58b44**) improves **673 → 675 correct / 870** with training-mean Fisher
centering. Classified accuracy is **675/849 (79.5%)**, with **33/145 exact chains
unchanged**. Both changed decisions improve a scan; none loses correctness.
Use **`config/unit_assignment_centered_fisher.toml`** to replay that opt-in
working candidate, not a promoted champion. Historical **677/870, 36/145** still
leads by two correct positions and three exact chains merely to equal it.

The independent final-covariance Ledoit-Wolf arm gives **665/849 (78.3%),
34/145 exact**: eight fewer correct positions despite one extra exact chain,
so it is not retained as the working replacement. All twelve control tables
replay exactly; all arms keep 146 scans / 900 keys, seven unavailable rows and
849/870 classified coverage. All 590 outputs and the Slurm log are fetched and
checksummed. The label-free numerical diagnostics improve covariance conditioning
but do not establish recognition or calibrated confidence. The full report is
`results/numerical_signals_20260921/report.md`; the runner is
`hpc/compare_numerical_signals.sbatch`. This comparison is complete: no automatic
resubmission, combination, post-grade tuning, counting or unknown25 rerun.
The extensively reused benchmark is development evidence, not independent validation.

The preceding independent support/score comparison (source **db411a2**, job
**11906168**, **7m58s**, exit **0:0**) improves **675/33 → 676/34** with complete
observed disks and symmetric affine-descriptor support. Use
**`config/unit_assignment_patch_support.toml`** for the latest opt-in candidate,
not a promoted champion. Classified coverage rises **849 → 852/870**; accuracy
is **676/852 (79.3%)**. Three newly usable predictions are correct, but two older
decisions regress: two scan gains, two losses, 141 ties, one exact chain gained
and none lost. Historical 677/870 and 36/145 still lead by one position and two
chains merely to equal them; historical coverage is 854/870.

The separate final Gaussian-density score gives **666/849 (78.4%), 32/145 exact**
and is rejected. All 146 scans / 900 keys remain; no imputation, normalization,
scan weighting, count refit or unknown25 rerun is included. The twelve control
tables replay 675/33 byte-for-byte; all 585 outputs and the Slurm log are fetched
and checksummed. Verification passes 1,604 focused assertions, 76 output checks
and paired-grade arithmetic. Full losses and replay commands are in
`results/patch_support_score_20260921/report.md`; the runner is
`hpc/compare_patch_support_score.sbatch`. This comparison is complete. Do not
resubmit it, combine its variants or tune settings after the grade. Predictions
remain label-free; reused benchmark results are not independent validation.

The subsequent complete-training comparison (source **f7a4074**, job
**11908779**, **6m42s**, exit **0:0**) is **negative**: **674/870 correct,
674/850 (79.3%), 33/145 exact**, against the exactly replayed **676/852,
34 exact** support control. Keep `config/unit_assignment_patch_support.toml`;
do not promote `unit_assignment_complete_training.toml`. The fixed training
rule uses 893 complete rows but the GMM can score only 894 rows: `240818_019`
has no complete patch for its per-file moments, so its two partial predictions
return to abstention. One scan improves, two regress, 142 tie; `240818_020`
loses exactness. All 146 files / 900 keys remain, including the excluded technical
scan in science. All 392 outputs and the Slurm log are fetched and checksummed.
Verification passes 1,782 focused assertions, 59 saved-output checks and paired
grade arithmetic. The initial audit's assumption that every scan has a complete
patch was false; the audit now checks the already-specified abstention behavior,
with the failed log retained. No scientific output or setting was changed.
Full losses: `results/complete_training_20260921/report.md`. This single-job
comparison is complete; no retry, normalization fallback, extra method, count
refit, DFT change or unknown25 rerun follows. Historical **677/36** remains ahead.

The next authorized comparison changes only per-scan GMM feature normalization
to median/Type-7-IQR (source **77b3ce6**, job **11911253**, **5m55s**, exit
**0:0**). It is **negative: 667/870 correct, 667/852 (78.3%), 10/145 exact**,
versus the exactly reproduced **676/852, 34 exact**. Coverage, selected N and
all upstream/Fisher/k-means tables remain identical. Keep
`config/unit_assignment_patch_support.toml`, not
`unit_assignment_robust_normalization.toml`. Against the control, 23 scans gain,
33 lose and 89 tie; no exact chain is gained and 24 are lost. All 392 outputs
and the Slurm log are fetched and checksummed. Verification passes 1,900 focused
assertions, 1,261 saved-output checks and paired arithmetic. An initial audit
compared arm-specific staging paths byte-for-byte; it was corrected to compare
scan identities and N before grading, without touching inference. See
`results/robust_normalization_20260921/report.md` for all losses. This comparison
is complete: no clipping, quantile search, scale-factor adjustment, weighting,
extra method, retry or post-grade tuning. Historical 677/36 remains the target;
the reused benchmark is development evidence, not independent validation.

The equal-scan GMM-weighting comparison (source **5ff68e2**, job **11915029**,
**6m01s**, exit **0:0**, after 41 minutes queued) gives a **tradeoff: 675/870
correct, 675/852 (79.2%), 36/145 exact**, versus the exactly replayed **676/852,
34 exact** support control. Retain `config/unit_assignment_patch_support.toml`
as the primary working reference; the weighting config is not promoted.
There are three scan gains, four losses and 138 ties; two exact chains are
gained and none lost. All seven changed decisions land at the existing vote
tie, with zero margin. All 146 scans / 900 keys, N and upstream/Fisher/k-means
tables remain unchanged. All 392 outputs and the Slurm log are fetched and
checksummed; 2,057 pre-run assertions, 360 output checks and paired arithmetic
pass. See `results/scan_weighting_20260921/report.md` for every gain and loss.
This fixed comparison is complete: no retry, additional arm, threshold change,
count refit, DFT change or unknown25 rerun follows. Historical **677/36** remains
unexceeded; the reused benchmark is not independent validation.

The continuous GMM seed-vote comparison (source **bef315f**, job **11922389**,
**7m19s**, exit **0:0**, after **36m23s** queued) is **negative: 671/870 correct,
671/852 (78.8%), 26/145 exact**, versus the exactly replayed **676/852,
34 exact** support control. Keep `config/unit_assignment_patch_support.toml`;
do not promote `unit_assignment_continuous_vote.toml`. Eight scans gain,
thirteen lose and 124 tie; one exact chain is gained and nine lost. All 24
changed decisions are 1→0 from old exact ties (GMM=1, k-means=0); only GMM
aggregation changes, not learning, naming, the final threshold or precision.
All 146 scans / 900 keys, N, availability and upstream/Fisher/k-means tables
are retained. All 392 outputs and the Slurm log are fetched/checksummed;
2,199 pre-run assertions, 7,240 output checks and paired arithmetic pass.
See `results/continuous_vote_20260921/report.md` for all gains and losses.
This comparison is complete: no retry, post-grade tie adjustment, new arm,
count refit, DFT change or unknown25 rerun follows. Historical **677/36** remains
unexceeded; inference is label-free, not independently validated by this reused
development benchmark.

The whole-scan GMM-bagging comparison (source **da886ce**, job **11925188**,
**7m17s**, exit **0:0**, after **29 seconds** queued) is **negative: 665/870
correct, 665/852 (78.1%), 10/145 exact**, versus the exactly replayed **676/852,
34 exact** support control. Keep `config/unit_assignment_patch_support.toml`;
the twenty-replicate bootstrap is not promoted. Twenty-two scans gain,
thirty-three lose and ninety tie; zero exact chains are gained and 24 lost.
Of 62 final 1→0 changes, 57 leave old GMM=1 / k-means=0 ties. All 200 seed
fits are accepted; each bag contains 88–98 distinct scans and 882–905 duplicated
usable rows. All 146 scans / 900 keys, selected N, availability and upstream/
Fisher/k-means tables are retained. All **393 outputs and the Slurm log** are
fetched/checksummed; **2,471 focused assertions, 21,927 pre-grade output checks**,
seven component-score checks and paired arithmetic pass. See
`results/scan_bagging_20260921/report.md` for every gain and loss. No retry,
post-grade fusion/threshold adjustment, new arm, count refit, DFT change or
unknown25 rerun follows. Historical **677/36** remains unexceeded; inference
is label-free, but this reused benchmark is not independent validation.

The tied-covariance comparison (source **becd0ce**, job **11931106**, **5m53s**,
exit **0:0**, after 15 seconds queued) is **negative: 666/870 correct,
666/852 (78.2%), 6/145 exact**, versus the exactly replayed **676/852,
34 exact** support control. Twenty-five scans gain, thirty-five lose and
eighty-five tie; zero exact chains are gained and 28 lost. Of 69 final changes,
68 switch 1→0 and one switches 0→1; 55 leave old GMM=1 / k-means=0 ties.
All ten candidate seed fits are named and have identical component covariances,
with freely learned masses. All 146 scans / 900 keys, selected N, availability
and upstream/Fisher/k-means tables remain. All **392 outputs and the Slurm log**
are fetched/checksummed; **2,716 focused assertions, 7,338 pre-grade output
checks**, five component-score checks and paired arithmetic pass. Keep
`config/unit_assignment_patch_support.toml`; historical **677/36** remains
unexceeded. Full gains/losses: `results/tied_covariance_20260922/report.md`.
No retry, other proposed lead, post-grade setting change, count refit, DFT
change or unknown25 rerun follows. This reused benchmark is development
evidence, not independent validation.

The next independent whole-scan Fisher / relative GMM-naming comparison
(source **43f8993**, job **11931478**, **8m30s**, exit **0:0**, after **3m01s**
queued) brings **no improvement**. Fisher gives **668/870 correct, 668/852
(78.4%), 34/145 exact**; relative naming gives **676/870, 676/852 (79.3%),
34/145**, identical to the exactly replayed support control. Fisher gains
on eight scans, loses on fifteen and ties on 122; four exact chains are
gained and four lost. All 26 final flips enter zero-margin GMM=1 / k-means=0
ties. Its two 73-scan groups are disjoint and both fits converge. Relative
naming preserves all ten fitted-parameter hashes and chooses the same group
names for every seed, leaving scores and decisions unchanged.

All **588 outputs and the Slurm log** are fetched/checksummed. Verification
passes **3,035 focused assertions**, **10,987 pre-grade output checks**, eight
saved-component checks and all 580 paired scan comparisons. All 146 scans /
900 keys, selected N, availability and common upstream/k-means outputs stay
fixed. Keep `config/unit_assignment_patch_support.toml`; historical **677/36**
remains unexceeded. Full gains/losses: `results/scan_fisher_naming_20260922/report.md`.
This two-candidate comparison is complete: no retry, combined candidate,
post-grade setting change, count refit, DFT change or unknown25 rerun follows.
Inference remains label-free; the reused grade is not independent validation.

The fixed factor-analyzer / Student learning comparison (source **e2205cb**,
job **11935072**, **0:0 in 8m17s**, after **1m36s** queued) is negative:
**620/870 correct, 23/145 exact** and **631/870, 25/145**, respectively, at
unchanged **852/870 coverage**. The support control replays **676/34** exactly.
Factors gain on 21 scans and lose on 65 (5 exact gains / 16 losses);
Student gains on 18 and loses on 57 (4 exact gains / 13 losses). All Student
fits converge in 65–86 updates; all factor fits reach the fixed 200-update
cap. That limits the MFA conclusion to this bounded implementation.

All **588 outputs plus Slurm log** are fetched/checksummed. Verification:
**3,616 assertions in 20 suites**, **11,083 pre-grade checks**, eight component
checks and 580 paired scan rows. All keys, N, upstream/Fisher/k-means tables
and availability remain unchanged. False positives rise 52→158 / 144, while
false negatives fall 124→74 / 77. Almost all final flips enter zero-margin
vote ties. Retain `config/unit_assignment_patch_support.toml`, still below
historical 677/36. No cap/rank/df/seed search, fusion change, combined arm,
count refit, DFT change or unknown25 rerun follows. Complete gains/losses:
`results/factor_student_mixtures_20260922/report.md`. The reused benchmark
remains development evidence, not independent generalization.

The subsequent all-update shrinkage / coherent Student-density comparison
(source **bb6943b**, job **11936229**, **0:0 in 7m56s**, after **3s** queued)
is also negative: **637/870 correct, 24/145 exact** and **627/870, 25/145**,
at unchanged **852/870 coverage**, versus the exactly replayed **676/34**
support control. Shrinkage gains on 19 scans and loses on 55 (6 exact gains /
16 losses); Student gains on 17 and loses on 58 (4 / 13). All ten fits per
candidate meet their stopping criterion (25/43 and 65–86 updates). Two
regularized fits have nonmonotone ordinary likelihood, as allowed and reported;
none is discarded or retried. False positives rise 52→135 / 152, exceeding
the recovered false negatives. Almost all changed decisions enter vote ties.

All **588 outputs and Slurm log** are fetched/checksummed; **4,277 assertions
in 21 suites**, **11,120 pre-grade checks**, eight component checks and 580
paired scan rows pass. All 146 scans / 900 keys, N, common upstream tables
and four unavailable rows remain. Keep `config/unit_assignment_patch_support.toml`;
historical 677/36 remains unexceeded. No post-grade coefficient/df/threshold
search, combination, count refit, DFT change or unknown25 rerun follows.
Details: `results/em_shrinkage_student_density_20260922/report.md`. This reused
benchmark is development evidence, not independent validation.

The raw-GCV reselection (source **853ad38**, job **11936561**, **0:0 in
1h25m43s**, after **26s** queued) is **negative**: **661/870 correct, 32/145
exact**, versus **675/34 for the fresh fixed-N control** and saved support
**676/34**. Exact counts fall **106 → 101/145**, coverage **852 → 849/870**;
emitted errors rise **177 → 188** against fresh control. Of 146 scans, 22
change N and the total falls 900 → 898 lobes. Counting gains/losses are 5/10;
recognition scan gains/losses/ties are 9/15/121, exact-chain gains/losses 3/5.
This repeats the original direct-extractor GCV rule, not the separately
promoted batch hybrid policy; the existing intelligent sweep is not exhaustive.

Both arms regenerate fixed-N geometry and all assignment stages. The fresh
control differs from saved support on three decisions, giving one fewer correct
position. Even same-N base/split geometry varies between fresh arms (maximum
center shifts 0.00845/0.07637 nm). Therefore the comparison is not a pure causal
N effect. All **438 outputs and Slurm log** are fetched/checksummed; **1,658
tests, 9,054 pre-grade checks and 37 external aggregate checks** pass, with 435
paired scan rows. No cohort loss or extra abstention is hidden. Retain saved
support, below historical **677/36**; no tuning, retry or unknown25 rerun follows.
Full evidence: `results/gcv_reselection_20260922/report.md`. This reused grade
does not establish independent generalization or unknown-chain chemical accuracy.

The fixed-geometry amplitude/background comparison (source **567cdde**, job
**11938782**, **0:0 in 6m17s**, after **5m29s** queued) is **negative**:
**675/870 correct, 675/852 (79.2%), 32/145 exact**, versus exactly replayed
support **676/34**. Coverage stays **852/870**; all N, centers, widths, axes
and split geometry remain identical. Against control: **4 scan gains,
5 losses, 136 ties; zero exact gains, two losses**. Errors rise **176→177**.

All **146 bounded linear profiles** satisfy KKT in one to three iterations;
RSS decreases by a median **0.681%**. This numerical improvement does not
translate into recognition gain. All **561 outputs and Slurm log** are
fetched/checksummed; **4,936 tests, 22,912 pre-grade checks and 31 external
aggregate checks** pass, with 290 paired scan rows. Thirteen control tables
replay byte-for-byte. Retain saved support, below historical **677/36**;
no nonlinear refit, selection, threshold change, retry or unknown25 rerun.
Full evidence: `results/frozen_amplitude_profile_20260922/report.md`.
The reused benchmark remains development evidence, not independent validation.

This page is the hand-off document for reproducing the current chitosan workflow
without relying on prior conversation context.  The goal is a label-free fitting
and model-selection pipeline that is externally graded on the benchmark and then
reused unchanged for curated 10–20mer images, except for the allowed `N` range.

## Scientific constraints

- Do not use an expected `N`, target `N`, or preferred `N` inside fitting or
  model selection.
- Benchmark labels are evaluation-only.  They may be used by grading scripts, but
  not by `test/batch_full.jl`, TOML calibration, or model-selection logic.
- Keep the circular-to-elliptical 2D pipeline: circular sweep first, then local
  elliptical refinement from circular solutions.
- Do not globally widen support for all files as a shortcut; global support
  changes were found to regress clean benchmark files.
- Ambiguous/suspicious 10–20mer images should be retained and annotated by QC or
  confidence fields, not silently excluded.

## Current configs

| Purpose | Config | Notes |
|---|---|---|
| Production/default short-chain chitosan | `config/chitosan.toml` | Current default is `support_midpoint_hybrid`, promoted as the best label-free full145 counting rule so far. |
| Benchmark-aligned adaptive workflow | `config/chitosan_adaptive_support_rescue.toml` | Experimental generic support-rescue policy. |
| 10–20mer adaptive workflow | `config/chitosan_10_20mer_adaptive_support_rescue.toml` | Same policy as benchmark adaptive workflow, with `n_max = 24`. |
| Legacy 10–20mer baseline | `config/chitosan_10_20mer.toml` | Raw GCV baseline/comparison. |
| Legacy 10–20mer support-rescue diagnostics | `config/chitosan_10_20mer_rescue*.toml` | Comparison/audit only, not ground truth. |

## Generic adaptive-support policy

`selection_policy = "adaptive_support_rescue"` performs the following steps for
each file:

1. Fit the standard support using the configured support detector.
2. Select `N_eff` from the circular/elliptical 2D sweeps using the configured
   criterion, currently GCV.
3. Trigger a permissive support-rescue pass only if `N_eff` is at the objective
   support-feasibility ceiling.  The ceiling uses the detected support length and
   physical spacing/sigma bounds; it does not use labels.
4. Accept the rescue only if support length increases, selected `N` increases,
   circular/elliptical counts remain coherent, and the selected count is feasible
   on the rescued support.
5. Apply the robust-AICc guard down-only on the active support.  If the guard
   fails, keep the active `N_eff` rather than failing the file.

The optional `adaptive_robust_guard_max_drop` parameter exists only for diagnostic
experiments.  It is not used by the benchmark-aligned workflow.

## Reproduce the short-chain benchmark workflow

Use the benchmark data folder and the robust-AICc guard override when the goal is
to reproduce the original 240817 primary-benchmark validation (`39/39`). The
plain `config/chitosan.toml` default is now optimized for the expanded full145
external counting grade instead.

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 39 \
  --data-dir /home/durif/Rebecca/data/data/20240817_LHe_Cu100 \
  --outdir results/best_plots_240817_robust_guard \
  --tsv /tmp/opencode/chitosan_240817_primary_files.tsv \
  --config config/chitosan.toml \
  --selection-policy gcv_with_robust_aicc_guard
```

Then grade externally against the benchmark manifest.  Labels stay outside the
fit/selection run:

```bash
julia --project=. test/grade_chitosan_benchmark.jl \
  --manifest benchmarks/chitosan_240817.toml \
  --results results/best_plots_240817_robust_guard/summary_overlap060_hard.tsv \
  --out results/benchmark_grades/chitosan_240817_robust_guard_N_selected.tsv \
  --column N_selected
```

Known validation from the current development pass: the robust-AICc guard alone
reaches `N_selected = 39/39` exact on the 240817 primary benchmark, with all four
`clean_target` files (`017`, `019`, `043`, `058`) reporting `N_selected = 6`.
`043` is recovered by the up-when-ambiguous guard branch (its `N_eff = 5`, but
`robust_AICc_N = 6` on an ambiguous file; see the Research Journal §2026-06-17).
The promoted `config/chitosan.toml` default is now `support_midpoint_hybrid`
because it improves the expanded 145-file external counting grade from `106/145`
exact to `129/145` exact (`143/145` within one lobe), confirmed by a full Viper
batch run with 146 `ok` rows and grading against the 145-file manifest.
Re-measure any time with the grade script below.

Reproducibility note: the batch is deterministic run-to-run on a given machine
(verified identical `N_selected` across 3 consecutive runs on 2026-06-17).
Divergences between a past recorded number and a fresh run indicate intervening
code changes, not run-to-run noise.

### Expanded counting benchmark

The canonical 240817 benchmark above remains the reproducible validation set.
For broader external counting grades, use:

```text
benchmarks/chitosan_6mer_counting_confirmed.toml
```

The saved production-policy counting summary is
`results/experiments/6mer_full146/pm2_confirm/summary_overlap060_hard.tsv`
(`support_midpoint_hybrid`). Reproduce its external report without refitting:

```bash
julia --project=. test/grade_chitosan_benchmark.jl \
    --manifest benchmarks/chitosan_6mer_counting_confirmed.toml \
    --results results/experiments/6mer_full146/pm2_confirm/summary_overlap060_hard.tsv \
    --file-column filepath --column N_selected \
    --out results/benchmark_grades/confirmed145_promoted_N_selected.tsv
```

Expected report: 129/145 exact and 143/145 within one lobe. The nearby
`adaptive_support_rescue` directory contains a historical policy, not this
production regression input.

This manifest contains 145 visually/manually confirmed `expected_N=6` files
derived from `benchmarks/chitosan_6mer_preassignment_review.tsv`. It is still an
external grading manifest only and must not be read by fitting or selection code.
For unit assignment, the same 145 files are the benchmark scope and the external
control sequence is `NKNNKN` (`010010` or `101101`, depending on the 0/1 identity
convention). That sequence is grading/control information only: it must not enter
the assignment method, threshold choice, abstention rule, composition prior, or
any calibration intended to extrapolate to unknown systems. Older pending lists
track provenance/curation state and should not be read as a license to shrink the
0/1/? benchmark to the 240817 subset.

## Run the 10–20mer adaptive workflow

Use the same adaptive policy with the 10–20mer config.  `--skip-1d` is
recommended because the 1D panels are not required for this long-chain workflow
and can be expensive.

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_adaptive_support_rescue \
  --tsv results/10_20mer_analysis_adaptive_support_rescue/triage_unused.tsv \
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --skip-1d
```

For a targeted smoke test of the support-rescue behavior:

```bash
julia --project=. test/batch_full.jl 1 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir /tmp/opencode/stmfit_10_20_adaptive_target \
  --tsv /tmp/opencode/one_260220_083.tsv \
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --skip-1d
```

Expected targeted result from the current pass: `260220_083.sxm` accepts support
rescue and reports `N_selected = 9`.

## Restarting and HPC usage

`test/batch_full.jl` appends to `summary_overlap060_hard.tsv` and skips files
already present in the selected `--outdir`.  If a run times out, rerun the same
command with the same output directory to continue from remaining files.

For HPC jobs, keep the same command-line arguments and only change scheduler
details, thread count, and output paths.  Prefer writing to a fresh output
directory per experiment.  Keep the final `summary_overlap060_hard.tsv`, the
per-file folders, and the exact config used.

## Outputs to inspect

The main table is:

```text
<outdir>/summary_overlap060_hard.tsv
```

Important columns:

- `N_eff`: raw effective 2D selection before refined policy/guard.
- `N_selected`: primary reported count after the configured batch policy.
- `selection_policy`: policy requested by config/CLI.
- `selection_source`: source of the primary selection, e.g. `ell` or
  `ell_robust_aicc`, plus `support_midpoint_down` or `support_midpoint_up` when
  the default support-midpoint layer makes the final one-step adjustment.
- `refined_policy`: audit trail for whether the robust-AICc/adaptive support
  stage was kept, accepted, rejected, or guarded before the final selection.
- `robust_aicc_N`: robust-AICc diagnostic/guard count.
- `support_2D_ell_nm`, `support_2D_circ_nm`: active support length after any
  accepted rescue.
- `best_plot`, `file_dir`: locations of visual outputs and per-file artifacts.

The 1D slide fit is **off by default** (it never enters `N_selected`; it is a
diagnostic only). `N_1D` and 1D comparison columns are therefore expected to be
`NA`, and 1D panels are not drawn. To re-enable the 1D diagnostic (e.g. to cross-
check a suspected 2D under-detection), pass `--no-skip-1d`.

## Minimal validation checklist after code/config changes

Run the default config smoke checks:

```bash
julia --project=. test/batch_full.jl 0 --config config/chitosan.toml
julia --project=. test/batch_full.jl 0 --config config/chitosan_10_20mer.toml --skip-1d
git diff --check
```

Run adaptive targeted checks when touching support rescue or selection logic:

```bash
julia --project=. test/batch_full.jl 1 \
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir /tmp/opencode/stmfit_10_20_adaptive_target \
  --tsv /tmp/opencode/one_260220_083.tsv \
  --skip-1d

julia --project=. test/batch_full.jl 2 \
  --config config/chitosan_adaptive_support_rescue.toml \
  --data-dir /home/durif/Rebecca/data/data/20240817_LHe_Cu100 \
  --outdir /tmp/opencode/stmfit_240817_adaptive_guard_targets \
  --tsv /tmp/opencode/two_240817_guard_regressions.tsv
```

Expected targeted results:

- `260220_083.sxm`: `N_selected = 9`.
- `240817_058.sxm`: `N_selected = 6`.
- `240817_019.sxm`: `N_selected = 6`.

## Known caveats

- The support-midpoint hybrid default is the best current label-free full145
  counting rule, not a universal solution; re-validate before using it as a new
  molecule's default.
- The adaptive configs are still experimental and are not the default production
  config until explicitly promoted.
- The legacy 10–20mer finalizer output is a comparison artifact, not ground
  truth.
- Full 10–20mer runs can be slow locally because the robust guard refits a broad
  range.  Use restartable output directories or HPC for complete reruns.
- If a file looks chemically or visually suspect, do not silently exclude it from
  10–20mer analysis; record the QC concern and keep the result available.

## Calibrating a new molecule

The pipeline is molecule-agnostic in its core (chain-of-Gaussians model,
label-free selection); only the calibration constants differ. To analyse a new
chain-like molecule under similar STM conditions:

1. **Copy the template**: `cp config/template.toml config/<molecule>.toml`.
2. **Re-derive the `[model]` values** from a few representative scans — the
   template comments explain each (FWHM → sigma, observed pitch → spacing,
   support length). These are the load-bearing changes.
3. **Treat the `[selection]` defaults** (`gcv_ambiguity_rel_threshold = 0.05`,
   `robust_guard_nu = 8.0`, `support_midpoint_up_gcv_rel_threshold = 0.30`) as
   chitosan-calibrated starting points. If the new molecule's lobe statistics
   differ markedly from chitosan, run
   `test/sensitivity_thresholds.jl` to check whether `N_selected` is robust to
   the threshold; re-calibrate only if it is sensitive.
4. **Exclude non-target files** via `--exclude-from results/<molecule>_exclude.txt`
   (one filename per line) rather than editing the batch code.

The selection path never uses an expected `N` or benchmark label, so the same
logic is label-free. What may need attention is *how often* the
up-when-ambiguous and support-midpoint branches fire for a molecule whose support
geometry or GCV curve has a different shape — hence the sensitivity check.

## Unit assignment workflow (GlcNAc/GlcN)

See [Unit Assignment](unit_assignment.md) for the full pipeline description.

### Unknown chitosan sequence production

<!-- UNKNOWN-CHITOSAN-WORKFLOW:START -->

For a chain with no unit-identity labels, start from the selected-N feature TSVs
produced by the label-free counting pipeline and write production artifacts only:

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

The production directory contains fixed-profile prediction TSVs, `summary.tsv`,
`manifest.tsv`, validation logs, plots, and `review_queue.tsv`. Treat `?` as an
explicit abstention requiring review, not as a dropped lobe.

<!-- UNKNOWN-CHITOSAN-WORKFLOW:END -->

#### Approved native reconstruction (September 16)

The lost `patch_u_asym` producer cannot be recovered. The approved replacement
is explicitly named `cc_soft_reconstructed_v1`; it is not presented as the frozen
champion under another name. Its fixed half-plane descriptor and native Julia
1.13 raw-to-QC command are documented in
[Explicit Julia reconstruction](unit_assignment.md#Explicit-Julia-reconstruction-(2026-09-16)).
Use `test/run_reconstructed_chitosan.jl --help` and `--dry-run` before a run.
The command has no benchmark/reference/grading input. The separate
`test/compare_reconstructed_champion.jl` reports keyed output differences after
prediction; external grading stays separate. Do not search descriptor variants
to recover the old aggregate score.

Native CC/Fisher tests and the extracted-input synthetic pipeline pass. One real
raw-input pilot (`251206_013.sxm`) also completes in 13 min 46 s: N_selected=11,
all 11 rows retained, QC and maps produced. Its isolated cohort is too small for
PCA10 in each Fisher half, so all assignments are explicitly `?`. This is a
technical pilot, not chemical validation; the pilot is not milestone completion.

**Initial September 17 attempt:** Raven job `30271312` failed at the base
geometry coverage check after 1 h 00 min 05 s. Fresh counting covers all 25 scans
with 222 selected lobes (range 5–16); geometry export retains 204 lobes over 23
files. `260215_022.sxm` and `260220_083.sxm` each lack all nine selected lobes.
These are the two scans whose counts used adaptive support rescue; the fixed-N
extractor in that run read only `N_selected` and lost that support context.
The wrapper correctly stopped before assignment. `failures.tsv` preserves all
25 intended files; no full146 run, comparison or new external grade occurred.
Fetched artifacts are under
`results/reconstructed_cc_soft_v1/raven_20260916_native_v1/unknown25/`.
The internal monitor is cancelled and no job has been resubmitted. The requested
handoff repair now reads `refined_policy` and reuses the existing rescue support
settings for both base and split refits, without changing N or any other setting.
Retain the original `--selected-summary` when resuming adaptive-support data;
N-only cache metadata cannot identify that support and is rejected when a
geometry refit is needed.
Both real rescued-file checks now pass under Julia 1.13. Each scan exports nine
base and nine split lobes, and its native geometry/support check passes 76
assertions. The geometry-only overlays were inspected; they do not validate
chemical assignments or prove that all fitted units are separately resolved.

| Scan | Reconstructed support (nm) | Base / split time (s) | Base / split GCV |
|---|---|---|---|
| `260215_022.sxm` | 5.552874796122559 | 74.49 / 84.94 | 5.63285686e-6 / 1.51047396e-6 |
| `260220_083.sxm` | 5.90756511515228 | 70.68 / 77.87 | 5.44278424e-6 / 2.55453904e-6 |

Outputs are under `results/reconstructed_cc_soft_v1/support_handoff_<scan>_v1/`
(without `.sxm` in `<scan>`). The original counting-fit GCV values were
3.18253642e-6 and 5.44586561e-6 respectively. This is fixed-N refitting on the
selected support, not exact parameter replay; no settings were retuned.

Native preparation merged the old 204 valid base rows with the two nine-row
recoveries, preserving every source row. Prepared inputs are in
`results/reconstructed_cc_soft_v1/unknown25_resume_v1_inputs/`: original counting
summary, complete base geometry (**25 files / 222 lobes**), genuine partial split
geometry (**2 files / 18 lobes**), fixed CC templates and 25 raw-input symlinks.
The preparation passes 29 assertions; the old failed run is unchanged.

The new `hpc/resume_reconstructed_unknown25.sbatch` is a separate application-only
continuation. It keeps the original summary, skips the two cached split files,
fits only the remaining **23 files / 204 lobes**, checks and merges the complete
split table, then runs cohort-wide assignment/QC/maps. It does not repeat
counting/base fits or run full146/grading. Its `--dry-run` checks metadata without
reading SXM pixels, fitting or creating outputs. See `hpc/README.md` for staging
and submission. **Job 30278010 completed on Raven on September 17**, from tested
source `aa422a3`: exit `0:0`, 6 min 16 s, eight allocated CPUs (four requested),
MaxRSS about 3.70 GiB. The original count/base tables and 18 cached split rows
were reused; the other 204 split rows were fitted. All downstream stages used
the full 25-chain cohort. No settings, selected counts or files were changed.
The failed run remains untouched and both monitors are stopped. Do not repeat
this submission.

The fetched application is in
`results/reconstructed_cc_soft_v1/raven_20260917_unknown25_resume_v1/unknown25/`:

- `predictions.tsv`: **222 keys / 25 files**, with **189 class-0, 33 class-1 and
  zero `?`** model outputs; all component rows are available.
- `summary.tsv`: the same selected counts (range 5–16), sequences and per-chain
  assignment counts. These are predictions, not measured chemical composition.
- `review_queue.tsv`: **19 chains flagged for low mean confidence**; one of
  those, `251206_038.sxm` (N=16), also has the cohort-relative `n_outlier` flag.
  Six chains have no listed review flag; that is not scientific validation.
- `plots/`: 25 standalone fitted-coordinate maps and `summary_grid.png`. The
  grid and seven priority standalone maps were inspected. These are not raw STM
  overlays, and their class colors do not display the confidence margin.

All 67 fetched file hashes (66 outputs plus Slurm log) match Raven. Local Julia
1.13 verification passes **179 assertions**, covering all 12 intermediate/final
lobe tables, exact cache reuse, saved-component vote replay, summary/QC and plot
coverage. The native prediction validator passes; regenerated QC and the final
vote table are byte-identical to their saved versions. No classifier or fit was
rerun locally.

`probability_1` is the mean of two component vote fractions, not a calibrated
probability; `confidence = 2 * abs(probability_1 - 0.5)`. The lobe-wise mean is
0.57016 and the minimum is 0.0375. Five margins are below 0.20 (a descriptive
count, not a new threshold). `?` denotes unavailable components, not all low
finite margins. The existing QC threshold is mean confidence <0.60. Read these
fields alongside the maps rather than interpreting zero `?` as certainty.
Further constrained-shape and channel-disagreement diagnostics are recorded in
[Unit assignment](unit_assignment.md#Completed-unknown25-application-(2026-09-17)).

No logged exception or Julia warning occurs. The split log for
`260222_043.sxm` reports a coverage-floor clamp from 6 to its already selected
fixed N=5; the five rows are retained, with no new count selection. This unknown25
application does not establish chemical accuracy. Its outputs are unchanged.

#### Completed full146 reconstruction and external evaluation (2026-09-18)

The single authorized **Viper job 11786116**, from `00528d2`, completed `0:0` in
**11m04s**. All **146 files / 900 cached GCV lobe keys** are retained through
assignment; this is not the distinct 871-lobe promoted counting summary or a
fresh count benchmark. See the [fixed-cache scope](unit_assignment.md#Approved-full146-comparison-scope-(2026-09-17)).

- Production: `results/reconstructed_cc_soft_v1/viper_20260917_full146_v1/`.
  Final labels are **695 zero / 198 one / 7 `?`**. Seven incomplete patch rows
  propagate to unavailable GMM input; no row is dropped or label forced.
- Integrity: all 184 outputs plus Slurm log match remote hashes. Julia 1.13
  saved-value checks pass **73/73**, including native validators and byte-identical
  descriptor/join/vote/summary/QC replays. All inputs and scientific settings
  remain unchanged. Clean logs do not prove every optimizer seed converged.
- QC: **43/146 chains flagged**, including 40 cohort-relative N outliers,
  six low mean margins and one high uncertain fraction (overlapping reasons).
  There are twelve finite zero-margin votes, still label 1 by the frozen rule.
  Mean confidence is 0.79082, not a calibrated chemical probability.
- Figures: the original giant overview was poorly laid out. Plot-only fix
  `6fa3297` leaves the 146 standalone PNGs byte-identical and creates **19 readable
  overview pages** plus `summary_pages.tsv`, under
  `results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_figures_v1/`.
  The original figures are preserved; these are fitted-coordinate maps, not raw
  STM overlays. The layout regression passes **267/267**.
- External comparison/grade:
  `results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_external/`.
  Agreement is **846/892 overlapping labels**, with zero missing reference keys
  and eight extra keys from the one file outside the confirmed145 manifest.
  Full145 own-N grading gives **666/849 correct (78.4%)**, **24/145 exact chains**,
  **16 missing control positions** and **38 extra predicted lobes**. Correctness
  over all 870 control positions is **666/870 (76.6%)**. There are 183 emitted
  errors and 21 unclassified control positions (16 missing + five unavailable).

The reference remains **677/854 (79.3%) / 36 exact chains**. This reconstruction
is not promoted, not an exact recovery of the lost descriptor, and not evidence
of unknown-chain chemical accuracy. No benchmark result fed back into fitting,
assignment, thresholds or abstention. The authorized full146 task is complete;
no old terminal challenger status below is reopened.


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
top-level `[grader_only]` table, and scans forbidden non-grader keys, values, and
paths case-insensitively. Its digest binds every other non-grader byte, including
comments, formatting, and line endings, excluding only real
`candidate.frozen_hash` and lifecycle-only `candidate.grade_status` assignments.
Matching provenance and digest are mandatory in `frozen_once` and `graded`, and
the lifecycle transition does not rehash. This is not stateless historical
proof. The current candidate remains locked and byte-unchanged; these stronger
future-state checks did not freeze it or consume the grade budget.

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

### Benchmark and post-hoc grading

The workflow for benchmark diagnostics is:

1. **Fill the ground truth** in
   `benchmarks/chitosan_240817_unit_sequences.tsv` (column `sequence`:
   ordered 0/1 along t_nm increasing; 0 = GlcN, 1 = GlcNAc).

2. **Extract per-lobe features** (re-runs the fit, ~10–15 min/file):
   ```bash
   STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. \
       test/extract_lobe_features.jl \
       --config config/chitosan.toml \
       --out results/unit_separability/lobe_features.tsv
   ```

3. **Run separability analysis** (label-free + with-truth):
   ```bash
   julia --project=. test/analyze_unit_separability.jl \
       --features results/unit_separability/lobe_features.tsv \
       --truth benchmarks/chitosan_240817_unit_sequences.tsv \
       --out results/unit_separability
   ```
   Check `results/unit_separability/separability_report.txt`:
   - ΔBIC (k=1−k=2) > 10 → bimodal (two distinct populations)
   - Physical accuracy > 75% → good separability
   - Physical accuracy < 60% → poor; use patch/mold diagnostics or DFT-STM maps

**Decision point:** if the separability report shows physical accuracy > 60%,
continue with label-free assignment diagnostics. If < 60%, prefer the connected
mold / DFT-STM path or a different STM bias/condition.

The ground truth is **never read by the fitter** — only by the separability
analysis (`--with-truth`) and the grading script.
