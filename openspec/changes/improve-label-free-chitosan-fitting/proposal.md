## Why

STMFit has a strong production baseline, but deciding whether a new preprocessing, fit model, or selector is genuinely better is still spread across research scripts and benchmark-facing reports. We need one reproducible, label-free challenger protocol that identifies the best physically defensible chitosan fit from image evidence alone, while keeping benchmark labels strictly external.

## What Changes

- Add a common candidate interface for chitosan preprocessing, 2D fitting, and `N_selected` selection, with every physical and selection parameter supplied through versioned TOML configuration.
- Add a label-free comparison protocol driven primarily by GCV, with residual structure, fit stability, physical feasibility, runtime, and abstention reported as diagnostics or guards.
- Evaluate candidates across representative short-chain and 10–20mer scans without using expected `N`, sequence labels, class counts, or benchmark outcomes during fitting, calibration, ranking, thresholding, or promotion.
- Require perturbation and held-out-scan stability evidence before a challenger can replace the current production fit path.
- Keep the 1D fit diagnostic-only and off by default; it may cross-check but never determine `N_selected`.
- Keep benchmark grading in a separate, locked post-selection step used only to report the externally measured effect of a frozen candidate.
- Record every experiment, failed approach, parameter decision, and terminal verdict in the scientific journal, with heavyweight comparisons executed on Viper.

## Capabilities

### New Capabilities
- `label-free-fit-selection`: Defines the candidate contract and the leakage-free rules for fitting, model selection, `N_selected`, physical guards, and abstention.
- `fit-robustness-evidence`: Defines reproducible candidate comparison across scans and perturbations using GCV-led, label-free evidence and a frozen external-grade boundary.

### Modified Capabilities

None.

## Impact

- Affects `packages/GaussianFit2D.jl/src/core.jl`, `packages/STMMolecularFit.jl/src/selectors.jl`, shared preprocessing owned by `packages/STMSXMIO.jl`, and orchestration in `test/batch_full.jl`.
- Adds focused contract tests, candidate configs under `config/`, comparison/report tooling under `test/`, Viper launch evidence under `hpc/`, and dated entries in `docs/src/journal.md`.
- Does not promote `STMMolecularFitGUI` as a production entrypoint and does not change the production fit until a frozen challenger passes the label-free evidence contract and a separate review.
- Does not use unit-assignment labels or the current constant-current diagnostic to select lobe count.
