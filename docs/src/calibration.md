# Calibration: what each parameter is based on

A parameter is label-free at inference when no label is read while fitting or
predicting. It is *independently calibrated* only when its value was also
chosen without benchmark grades. This page states which parameters satisfy
which condition.

## Counting (`config/chitosan.toml`)

| Parameter | Value | Basis |
|---|---|---|
| `sigma_parallel_min/max_nm` | 0.191 / 0.509 | Apparent lobe FWHM 0.45–1.20 nm divided by 2.355; a physical range, not an isolated-monomer measurement |
| `spacing_min/max_nm` | 0.35 / 0.75 | ±~30% around the chitosan repeat pitch (~0.52 nm); an assumption |
| `max_overlap` | 0.60 | Physical constraint on adjacent-lobe overlap. If high-N fits disappear, inspect it before changing selection logic |
| `fit_width_nm` | 0.16 | Tube half-width. Chosen 0.15 → 0.16 on benchmark counting sweeps |
| `support_noise_k`, `support_padding_nm` | 2.5, 0.25 | Support threshold and padding. Padding chosen 0.20 → 0.25 on benchmark sweeps |
| `kappa_max` | 10 | Condition-number penalty. Chosen 8 → 10 on benchmark sweeps |
| `selection_criterion` | `gcv` | Methodological choice ([Model selection](selection.md)) |
| `robust_guard_nu`, `gcv_ambiguity_rel_threshold`, `support_midpoint_up_gcv_rel_threshold` | 8, 0.05, 0.30 | Hybrid rule thresholds, adopted after benchmark grades |

**Provenance limit.** The hybrid rule, support padding, fit width and κ
threshold were selected while looking at known-count benchmark grades. The
promoted method reads no label, but its counting stage is not independently
calibrated, and its benchmark counts are development evidence, not an
independent validation. Do not tune these values further to recover known cases.

## Molecule consensus (`config/molecule_consensus.toml`)

These settings add no fitted physical parameter and were declared before the
consensus stage was graded:

- `grid_step_nm = 0.04`, `highpass_sigma_nm = 0.32`, `frame_margin_nm = 0.3`,
  `roi_margin_nm = 1.2` follow the chitosan lobe scale (~0.6 nm spacing).
- `max_shift_nm = 1.6`, `max_center_distance_nm = 6.0` bound thermal drift
  between consecutive scans.
- `ncc_min = 0.5` separates unrelated structure (NCC near 0) from an unchanged
  molecule (near 1). The label-free distribution of the 111 registered
  benchmark pairs has an empty gap between 0.39 and 0.55, so no tested link sits
  near the threshold.
- `min_track_scans = 3` with a strict majority is the smallest evidence that
  can outvote one deviating scan.
- Fusion EM starts from a neutral (π, θ0, θ1) = (0.5, 0.1, 0.9); the fit is
  start-independent on the benchmark cohort. `fusion_min_scans = 3`.

The consensus needs repeated scans of one molecule. Cohorts of single scans are
unchanged by it.

## Unit assignment (`config/unit_assignment_*.toml`)

Assignment settings were fixed before each graded comparison. The kept configs
mark a lineage of declared single changes (intermediate configs are in the
archive tag):

`unit_assignment_reconstructed.toml` (native reconstruction) →
`unit_assignment_matched_residual.toml` (`S(data − model)` residuals) →
`unit_assignment_patch_support.toml` (after three accepted steps: affine
descriptor, training-mean Fisher centring, complete-disk patch support) →
`unit_assignment_corroborated_training.toml` (learning restricted to corroborated
scans; promoted). `unit_assignment_transverse_fisher.toml` is a test fixture
for the physical Fisher layout.

The mold construction height (0.50 nm target) and isovalue grid come from the
DFT workflow, not from grades. The history of these choices is in the
[journal](journal.md): many single-change variants were graded and rejected,
and none was retuned after grading.

## Measurement audit (`test/measure_calibration.jl`)

The tool reports apparent widths and spacings in both scan directions,
missing measurements, and where the former bootstrap substituted fallback
values. It does **not** write a production config. On the full146 cohort
(September 24), the old bootstrap needed fallbacks for widths on 245/292 views
and for spacings on 261/292. The observed-only measurements are more available
but not more consistent between forward and backward views. No automatic
physical calibration is justified by them. Settings: `config/calibration_measurements.toml`;
cohort comparison: `test/summarize_calibration_measurements.jl`.

## Effective sample size

The `n ÷ 9` effective sample size is a placeholder used only by BIC/AICc
diagnostics and guards. The residual correlation range exceeds the fit tube,
so an effective sample size is not identifiable there. GCV needs none and drives
selection ([Model selection](selection.md#Why-GCV)).

## Calibrating a new molecule

1. Copy `config/template.toml`. Inspect raw scans in both directions with
   `test/measure_calibration.jl` and `test/inspect_one_file.jl`.
2. Separate unavailable values, apparent measurements and independent physical
   evidence. Do not take a median of silently substituted defaults.
3. Set widths, spacing, overlap and `n_max` from independent evidence (FWHM,
   known pitch, chain length), and record their basis.
4. Freeze the config before fitting and before any external grading. Visual QC
   and repeated benchmark scores do not establish unknown-chain chemistry.
5. Noise and correlation depend on tip, session and preprocessing. Recheck them
   when transferring a calibration, even with the same instrument.
