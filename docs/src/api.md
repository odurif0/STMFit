# API Reference

## GaussianFit2D — Core Fitting

```@docs
GaussianFit2D.read_sxm
GaussianFit2D.ChainSweepConfig
GaussianFit2D.PatternConfig
GaussianFit2D.chain_gaussian_sweep
GaussianFit2D.chain_direct_fit
GaussianFit2D.fit_chain_consensus
GaussianFit2D.fit_chain_batch
GaussianFit2D.MolecularFeature
GaussianFit2D.ChainModelResult
```

## GaussianFit1D — 1D Profile Fitting

```@docs
GaussianFit1D.FitConfig
GaussianFit1D.run_fit
GaussianFit1D.build_config
GaussianFit1D.FitResult
```

## STMMolecularFit — Orchestration

```@docs
STMMolecularFit.SlideConfig
STMMolecularFit.FitSlideConfig
STMMolecularFit.extract_slide
STMMolecularFit.fit_slide
STMMolecularFit.extract_and_fit_slide
```

## STMFitCore — Shared Utilities

```@docs
STMFitCore.effective_spacing_min
STMFitCore.kappa_penalty
STMFitCore.adjacent_kappa_max
STMFitCore.endpoint_overrun
STMFitCore.overlap_condition_number
STMFitCore.durbin_watson
STMFitCore.runs_test
STMFitCore.ResidualDiagnostics
STMFitCore.compute_residual_diagnostics
```
