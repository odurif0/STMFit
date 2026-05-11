# Multi-Gaussian Fit

Core Julia package for multi-Gaussian fitting of STM/AFM line profiles.

This repository contains **only the fitting core**. The web GUI lives in a separate repository:
[`STMMolecularFitGUI.jl`](https://github.com/od/STMMolecularFitGUI.jl).

## Quick Start

### Installation

```julia
using Pkg
Pkg.develop(path="/home/durif/Git/GaussianFit1D.jl")
```

### CLI usage

```bash
cd /home/durif/Git/GaussianFit1D.jl
julia --project=. run.jl data.txt --min-spacing 0.4
```

### As a library

```julia
using GaussianFit1D

x, y = load_data("data.txt")
cfg = build_config(Dict("min_spacing" => 0.4))
fr = run_fit(x, y, cfg)
best = best_result(fr)
println("Best: n=$(best.n_peaks), BIC=$(best.bic)")
```

### GUI

The web GUI is in a separate repository: `STMMolecularFitGUI.jl`.

## Model

The profile intensity is modeled as a sum of N Gaussian peaks
plus a baseline y₀ fixed at zero (data are auto-offset):

```
I(x) = Σᵢ Aᵢ · exp( -(x - μᵢ)² / (2 σᵢ²) )
```

Peak positions are parameterized as μ_{i+1} = μ_i + δ_i to
enforce ordered centers (3 parameters per peak: A, δ, σ).
An optional **asymmetric edge** mode allows the outermost peak on
each side to use a different width (σ_in ≠ σ_out):

```
x < μ₁ :  σ = σ₁^outer    x ≥ μ₁ :  σ = σ₁
x < μ_N :  σ = σ_N        x ≥ μ_N :  σ = σ_N^outer
```

| Peaks | Free parameters (symmetric) | Free parameters (asymmetric) |
|-------|--------------------------|---------------------------|
| $N$ | $3N+1$ | $3N+3$ |

Where 3N+1 = baseline y₀ + N×(amplitude A, spacing δ, width σ). Asymmetric mode
adds 2 extra free sigma parameters for the outer edges of the first and last peaks,
yielding 3N+3 total.

## Fitting Procedure

1. **Global optimization** — NLopt `GN_DIRECT_L` (derivative-free,
   bound-constrained) finds the basin of the global minimum.
2. **Local refinement** — Levenberg-Marquardt (`LsqFit.jl`, finite differences)
   polishes the solution and estimates the parameter covariance matrix.
3. **Bidirectional sweep** — Starting from a physically-motivated
   `center_n`, the sweep expands outward in both directions,
   fitting $N \in [N_\text{min}, N_\text{max}]$.
   Every model uses an equidistant midpoint initialization.
4. **Model selection** — The best model is chosen by **BIC**
   (Bayesian Information Criterion). Competing models within
   `bic_competition_threshold` are flagged.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_spacing` | 0.4 nm | Minimum distance between adjacent peaks |
| `max_spacing` | 0.675 nm | Maximum peak spacing |
| `fwhm_min` | 0.45 nm | Minimum allowed FWHM |
| `fwhm_max` | 1.2 nm | Maximum allowed FWHM |
| `offset_to_zero` | true | Subtract the global minimum from Y-values |
| `amplitude_min_fraction` | 0.3 | Minimum peak amplitude as fraction of max(Y) |
| `asymmetric_edges` | true | Allow different inner/outer sigma for edge peaks |
| `edge_sigma_min/max` | 1.0 / 4.0 | Edge margin factor (× FWHM→sigma) for edge peaks |
| `global_maxiter` | 5000 | NLopt global optimizer max iterations |
| `global_maxtime` | 15 s | NLopt time limit per fit |
| `global_tol` | 1e-5 | NLopt convergence tolerance |
| `early_stop_patience` | 3 | Consecutive failed fits before stopping sweep |
| `early_stop_dbic` | 100 | BIC degradation threshold for early stop |
| `bic_competition_threshold` | 20 | ΔBIC cutoff for marking models as competitive |

All parameters are adjustable from the GUI or via `build_config()` in the CLI.
