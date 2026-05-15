# Mathematical Background

## Condition Number κ

When two Gaussian peaks overlap, their parameters become correlated. The
**condition number** κ of the 2×2 Gram matrix quantifies this:

```
ρ = exp(-d²/(4σ²))       normalized correlation
κ = (1 + ρ)/(1 - ρ)      condition number
```

where `d` is the peak separation and `σ = max(σ_i, σ_{i+1})` is the larger
width of the pair.

| Regime | κ | Meaning |
|--------|---|---------|
| Well-separated (d >> σ) | ≈ 1 | Independent parameters |
| Typical chitosan | 4–8 | Moderate correlation |
| Strongly overlapping | > 25 | Parameter uncertainty amplified 25× |
| Coincident (d = 0) | ∞ | Parameters unresolvable |

### Progressive Penalty

The penalty is zero below `kappa_max` (default 8), then ramps quadratically:

```
penalty(κ) = 0                                    if κ ≤ κ_max
           = weight × ((κ - κ_max)/κ_max)²        if κ > κ_max

objective = RSS × (1 + penalty(κ))
```

Applied **only in the NLopt global stage** — the LsqFit LM refinement is
unpenalized to keep the final fit unbiased.

## Pseudo-Voigt Profile

The pseudo-Voigt is a weighted sum of Gaussian and Lorentzian:

```
V(x; η) = (1 - η) · G(x) + η · L(x)

G(x) = exp(-½(x/σ)²)           Gaussian
L(x) = 1 / (1 + (x/σ)²)        Lorentzian
```

- `η = 0` → pure Gaussian (default, STM tip-limited broadening)
- `η = 1` → pure Lorentzian (homogeneous broadening)
- `η ∈ (0, 1)` → mixed profile

The mixing parameter η is **global** (shared across all peaks in a fit),
justified by the fact that peak shape in STM is determined by tip condition
and imaging parameters, not by individual molecules.

The pseudo-Voigt adds exactly **one parameter** to the fit regardless of the
number of peaks. Model selection (BIC/AICc) automatically determines whether
the extra parameter is warranted by the data.

## Student-t BIC

Standard BIC assumes Gaussian residuals. STM data often has heavier tails
(outlier pixels, tip changes). The Student-t NLL is more robust:

```
NLL_student = Σ_i  ½(ν+1) · log(1 + r_i²/(ν·σ²))
BIC_student = 2 · NLL_student + k · log(n_eff)
```

where `ν = 4` (degrees of freedom, heavier tails than Gaussian), `σ` is the
noise estimate (MAD-based), and `n_eff` is the effective sample size.

For 2D images, `n_eff = max(10, length(zfit) ÷ 9)` — each 3×3 pixel block
counts as one independent observation, accounting for spatial correlation.

## Residual Diagnostics

### Durbin-Watson Statistic

Tests for autocorrelation in fit residuals:

```
DW = Σ(r_i - r_{i-1})² / Σ r_i²
```

| DW value | Interpretation |
|----------|---------------|
| ≈ 2 | No autocorrelation — good fit |
| < 1.5 | Positive autocorrelation — missed structure |
| > 2.5 | Negative autocorrelation — overfitting |

A low DW indicates systematic structure remaining in the residuals, suggesting
the model is missing a peak or using an incorrect profile shape.

### Wald-Wolfowitz Runs Test

Tests whether the signs of residuals are random:

```
A "run" is a maximal sequence of consecutive residuals with the same sign.
Expected runs = 1 + 2·n_pos·n_neg / n
```

- Too few runs → systematic bias (fit consistently above or below data)
- Too many runs → oscillation (model is overfitting noise)

Both tests return a p-value; p < 0.05 suggests a problem with the fit.

## Uncertainty Quantification

### 1D: Laplace Approximation

After LsqFit's Levenberg-Marquardt refinement, the covariance matrix is:

```
Σ = σ̂² · (J'J)⁻¹
```

where J is the Jacobian at the optimum. The **correlation matrix** is:

```
corr(i,j) = Σ(i,j) / √(Σ(i,i) · Σ(j,j))
```

### Center-Center Correlation

Peak centers in the delta-parameterization are cumulative sums:

```
c_0 = μ_0
c_i = μ_0 + Σ_{j=1}^{i} Δ_j
```

Their covariance is propagated via the Jacobian of this transformation:

```
Σ_centers = J · Σ_params · J'
```

Strong center-center correlations (> 0.9) indicate that peak positions are
poorly constrained relative to each other — the κ penalty helps avoid this.

### 2D: Jacobian-based Errors

The 2D fit also uses LsqFit for local refinement. Parameter errors are
extracted from `estimate_covar(fit)` and stored per result. These are
errors in the *encoded* (sigmoid-transformed) parameter space.

## Generalized Cross-Validation (GCV)

The default CV method (`cv_method="gcv"`) replaces expensive k-fold
cross-validation with an analytical approximation:

```
GCV = (1/n) · Σ ρ_ν(r_i · n/(n-p))
```

where `n` is the number of data points, `p` the number of model parameters,
`r_i` the residuals, and `ρ_ν` the Student-t loss function.

GCV is derived from leave-one-out CV under the linear approximation of the
model at the optimum. For converged nonlinear fits it provides an excellent
approximation of LOOCV at zero computational cost.

**Performance**: k-fold CV requires `k` full refits per N value in the sweep,
typically adding 80–100s per N for 5-fold CV on 2D images. GCV replaces this
with a single O(n) pass over the residuals, reducing the total 2D sweep time
by roughly 5×.

The k-fold method remains available via `cv_method="kfold"` for validation
or when the linear approximation is suspected to be inadequate.
