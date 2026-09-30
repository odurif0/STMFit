# Configuration

Every physical, selection and assignment parameter lives in a TOML file under
`config/`, never in a hidden code default. Each file has `[model]`,
`[selection]` and `[preprocessing]` sections. Configs used by a promoted run
are frozen: a change is a new experiment with its own journal entry, and it
must update this page and [Calibration](calibration.md).

| File | Used by | Status |
|---|---|---|
| `chitosan.toml` | Per-scan counting (6-mer benchmark, default) | Promoted |
| `chitosan_10_20mer_adaptive_support_rescue.toml` | Counting for 10–20mer application | Application config |
| `molecule_consensus.toml` | Registration, consensus, repeat-scan collection, fusion | Promoted |
| `unit_assignment_corroborated_training.toml` | Per-scan assignment | Promoted |
| `unit_assignment_patch_support.toml` | Same, learning on all scans | Previous reference (control) |
| `unit_assignment_reconstructed.toml`, `_matched_residual.toml`, `_transverse_fisher.toml` | Lineage and test fixtures | Not for production |
| `calibration_measurements.toml` | `test/measure_calibration.jl` audit | Diagnostic |
| `template.toml` | Starting point for a new molecule | Template |

## Counting configs

`[model]` (read into `GaussianFit2D.ChainSweepConfig`):

| Key | chitosan | Meaning |
|---|---|---|
| `n_min`, `n_max` | 2, 14 (struct default) | Sweep bounds; the 10–20mer config sets `n_max = 24` |
| `spacing_min_nm`, `spacing_max_nm` | 0.35, 0.75 | Inter-lobe spacing bounds |
| `max_overlap` | 0.60 | Maximum adjacent-lobe overlap fraction |
| `sigma_parallel_min_nm`, `sigma_parallel_max_nm` | 0.191, 0.509 | Axial σ bounds (FWHM 0.45–1.20 nm) |
| `fit_width_nm` | 0.16 | Tube half-width around the axis |
| `support_noise_k`, `support_padding_nm`, `support_min_length_nm`, `support_baseline_quantile` | 2.5, 0.25, 1.0, 0.10 | Axial support detection |
| `global_maxtime`, `global_maxiter`, `max_iter`, `multistart` | 0, 10000, 300, 1 | NLopt and LsqFit budgets. `global_maxtime = 0` disables the time limit, so the global search is bounded by `global_maxiter` and fits are deterministic (2026-09-30) |
| `support_threshold_rule`, `support_threshold_fraction` | `profile_quantile_noise`, 0.5 | Axial support threshold. `half_maximum_cap` caps the legacy threshold at off-ROI background + fraction × (peak − background); available, not selected ([Model selection](selection.md)) |
| `selection_criterion`, `cv_method` | `gcv`, `gcv` | Per-candidate score (`bic`, `aicc`, `cv` are diagnostics) |
| `selection_policy` | `support_midpoint_hybrid` | Batch policy ([Model selection](selection.md)); 10–20mer: `adaptive_support_rescue` |
| `cv_folds`, `bic_cv_margin`, `cv_ratio_threshold` | 5, 100, 2.0 | Used only by k-fold CV diagnostics |
| `kappa_max`, `kappa_weight` | 10, 1 | Condition-number penalty |
| `min_amplitude_fraction` | 0.3 | Minimum lobe amplitude (fraction of the data maximum) |
| `chain_tilted_baseline`, `chain_circular_sigmas` | true, false | Model variants |
| `adaptive_rescue_*` | 10–20mer only | Rescue-pass support settings (`noise_k = 1.5`, `padding_nm = 0.75`, `min_support_gain_nm = 0.75`, `feasible_margin = 0`, `max_circ_ell_delta = 1`) |

`[selection]`: `gcv_ambiguity_rel_threshold = 0.05` (two N are
indistinguishable below this relative GCV gap), `robust_guard_nu = 8.0`,
`support_midpoint_up_gcv_rel_threshold = 0.30`.

`[preprocessing]`: `stride = 1`, `flatten = "plane+rows"`,
`smooth_radius_px = 1`. The same settings drive the assignment patch export.

Other `ChainSweepConfig` fields keep their struct defaults
(`sigma_perp_min/max_nm = 0.10/0.55`, `intelligent_sweep = true`,
`early_stop_patience = 3`, `student_nu = 4.0`); see
`packages/GaussianFit2D.jl/src/core.jl`.

### `test/batch_full.jl` options

`[N_files]` (first positional argument), `--config`, `--data-dir`, `--outdir`,
`--chunk i/n`, `--tsv` (optional triage input), `--skip-1d` (default) /
`--no-skip-1d`, `--selection-policy`, `--gcv-ambiguity-rel-threshold`,
`--robust-guard-nu`, `--exclude-from FILE`, and `--plot-manifest TOML` with
`--skip-plot-quality` (plots only; never selection). The batch appends to `summary_overlap060_hard.tsv` and
skips files already done in `--outdir`, so it can be resumed.

## `molecule_consensus.toml`

| Key | Value | Meaning |
|---|---|---|
| `[model] window_half_nm` | 3.0 | Registration window around the earlier scan's chain centroid |
| `grid_step_nm` | 0.04 | Absolute-frame resampling step |
| `highpass_sigma_nm` | 0.32 | Removes structure broader than ~half a lobe spacing |
| `max_shift_nm`, `coarse_step_px` | 1.6, 4 | Largest drift between consecutive scans; coarse search step |
| `roi_margin_nm` | 1.2 | Margin of the registered-ROI refit |
| `fusion_em_max_iter`, `fusion_em_tol` | 500, 1e-10 | Latent-class EM stopping |
| `fusion_init_pi`, `fusion_init_theta0`, `fusion_init_theta1` | 0.5, 0.1, 0.9 | Neutral EM start |
| `[selection] ncc_min` | 0.5 | Link threshold |
| `min_overlap_px`, `max_center_distance_nm` | 500, 6.0 | Registration overlap and centroid distance limits |
| `min_track_scans`, `majority` | 3, `strict` | Consensus evidence |
| `frame_margin_nm` | 0.3 | Registered footprint must lie this far inside the frame |
| `collect_max_gap_min`, `collect_max_offset_nm`, `collect_max_range_nm` | 30, 5, 20 | `test/collect_repeat_scans.jl` candidate rules |
| `fusion`, `fusion_min_scans` | `latent_class`, 3 | Fusion model and minimum scans per physical lobe |
| `[preprocessing] channel`, `flatten` | `Z`, `plane+rows` | Registration image |

## Assignment configs (`unit_assignment_*.toml`)

Values of the promoted `unit_assignment_corroborated_training.toml`:

| Group | Keys |
|---|---|
| Identity | `model.name = "cc_soft_patch_support_corroborated_v1"` |
| Descriptor | `descriptor = "affine_residual_half_plane_asymmetry"`, `descriptor_column = "patch_u_asym_reconstructed"`, `descriptor_channel = "bwd_res"`, `descriptor_half_nm = 0.32`, `descriptor_step_nm = 0.08`, `descriptor_zero_l1 = 1e-12`, `preprocessing.descriptor_normalization = "median_sample_std"` |
| Split fit | `split_skew_ratio_max = 2.0` |
| Molds | `mold_half_nm = 0.32`, `mold_step_nm = 0.04`, `mold_target_height_nm = 0.50`, `mold_z_min/max/step_nm = -0.5/2.6/0.005`, `mold_isovalue_min_log10 = -5.5`, `mold_isovalue_max_fraction = 0.8`, `mold_isovalue_count = 80`, `mold_calibration = "first_below_target"`, `mold_margin_mode = "absolute_cost_margin"` |
| Fisher | `fisher_pca_components = 10`, `fisher_noise_regularization = 0.01`, `fisher_gmm_regularization = 0.001`, `fisher_gmm_maxiter = 300`, `fisher_gmm_tolerance = 0.001`, `fisher_seed = 0`, `fisher_score_center = "training_mean"`, `fisher_projection_zero_l1 = 1e-12`, `selection.fisher_cv_scheme = "lobe_parity"`, `fisher_scan_split_seed = 0`, `preprocessing.fisher_patch_projection = "none"`, `preprocessing.fisher_layout = "legacy_row_major"` |
| GMM | `gmm_learning_family = "gaussian"`, `gmm_covariance_scope = "final_only"`, `gmm_covariance_structure = "full"`, `gmm_hard_assignment = "mahalanobis"`, `gmm_final_score = "mahalanobis"`, `gmm_final_covariance = "ridge"`, `gmm_covariance_ridge = 1e-6`, `gmm_learning_maxiter = 200`, `gmm_learning_tolerance = 1e-6`, `gmm_learning_min_mass = 1e-12`, `selection.gmm_cluster_naming = "raw_amplitude"`, `gmm_resampling = "none"`, `gmm_seed_aggregation = "hard_vote"`, `gmm_training_weighting = "equal_lobes"`, `gmm_seeds = 10`, `gmm_selftrain = 2`, `preprocessing.gmm_feature_normalization = "mean_sample_std"`, `gmm_scale_fallback = 1.0` |
| k-means and vote | `kmeans_seeds = 20`, `first_seed = 0`, `interactions = true`, `vote_threshold = 0.5` |
| Training | `assignment_training_support = "all_admissible"`, `assignment_training_scans = "corroborated_counts"` |
| Patches | `assignment_patch_support = "complete_disk_symmetric"`, `patch_residual_filter = "smooth_residual"`, `pixel_order = "u_outer_t_inner"` |

Several keys are policy switches whose other historical values were rejected.
`test/build_labelfree_gmm_predictions.jl` now accepts only the promoted value and
fails otherwise. The consensus driver requires `patch_residual_filter =
"smooth_residual"`.

## `calibration_measurements.toml`

Descriptive settings of the measurement audit: bright-axis quantile, strip and
bin widths, width quantiles, prominence filter, and the former bootstrap's
fallbacks (reported, never substituted). It is not a fitting calibration.

## `template.toml`

Annotated starting point for a new molecule. Only `[model]` is mandatory;
see [Calibration](calibration.md#Calibrating-a-new-molecule).
