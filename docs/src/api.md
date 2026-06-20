# API Reference

## GaussianFit2D — Core Fitting

Main public entry points:

| Object | Purpose |
|--------|---------|
| `STMSXMIO.read_sxm(path)` | Read Nanonis STM `.sxm` files. Owned by `STMSXMIO`; re-exported by `GaussianFit2D` and `STMMolecularFit`. |
| `STMSXMIO.SXMImage` | Parsed SXM image (header, channels, geometry). |
| `GaussianFit2D.PatternConfig` | Image preprocessing, ROI, channel and output configuration. |
| `GaussianFit2D.ChainSweepConfig` | 2D chain sweep, support, physical constraints, optimization and selection configuration. |
| `GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)` | Adaptive 2D chain sweep over N. Circular runs use deterministic 2D-only initialization. |
| `GaussianFit2D.chain_direct_fit(img, pcfg, ccfg)` | Single-N direct 2D fit, useful for diagnostics or externally supplied initialization. |
| `GaussianFit2D.fit_chain_consensus(img, pcfg, ccfg)` | Fit chain models on multiple channels/views for consistency checks. |
| `GaussianFit2D.fit_chain_batch(files, pcfg, ccfg)` | Batch helper for package-level chain sweeps. |
| `GaussianFit2D.MolecularFeature` | Fitted Gaussian lobe: amplitude, position, widths. |
| `GaussianFit2D.ChainModelResult` | Fit result for one N, including scores, validity, residual diagnostics and metrics. |

Important internal pipeline functions documented conceptually in
[Pipeline Architecture](pipeline.md): `_active_t_support`,
`_deterministic_chain_seed`, `_fit_chain_n`, and the batch-level
`_refine_circ_to_ell` helper.

## GaussianFit1D — 1D Profile Fitting

| Object | Purpose |
|--------|---------|
| `GaussianFit1D.FitConfig` | 1D peak fitting configuration. |
| `GaussianFit1D.build_config(dict)` | Build a `FitConfig` from key/value settings. |
| `GaussianFit1D.run_fit(cfg)` | Sweep and fit 1D multi-Gaussian profiles. |
| `GaussianFit1D.FitResult` | 1D fit output with scores, residual diagnostics and covariance metrics. |

## STMMolecularFit — Orchestration

| Object | Purpose |
|--------|---------|
| `STMMolecularFit.SlideConfig` | 1D slide extraction configuration. Support is noise-based via `support_noise_k` and `support_padding_nm`. |
| `STMMolecularFit.FitSlideConfig` | Configuration forwarded to the 1D Gaussian fitter. |
| `STMMolecularFit.extract_slide(img, cfg)` | Extract an axial slide/profile from an STM molecular ROI. |
| `STMMolecularFit.fit_slide(slide, cfg)` | Fit the extracted slide with `GaussianFit1D`. |
| `STMMolecularFit.extract_and_fit_slide(path; ...)` | Convenience end-to-end 1D extraction and fit. |

## STMFitCore — Shared Utilities

| Object | Purpose |
|--------|---------|
| `STMFitCore.effective_spacing_min` | Converts spacing and overlap constraints into an effective minimum spacing. |
| `STMFitCore.kappa_penalty` | Soft penalty for ill-conditioned adjacent lobes. |
| `STMFitCore.adjacent_kappa_max` | Maximum adjacent condition number across a chain. |
| `STMFitCore.endpoint_overrun` | Measures lobe span exceeding detected support. |
| `STMFitCore.overlap_condition_number` | Pairwise overlap conditioning utility. |
| `STMFitCore.durbin_watson` | Residual autocorrelation diagnostic. |
| `STMFitCore.runs_test` | Residual sign-run diagnostic. |
| `STMFitCore.ResidualDiagnostics` | Container for residual QC metrics. |
| `STMFitCore.compute_residual_diagnostics` | Compute residual QC metrics for a fitted model. |
