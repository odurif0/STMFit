# diagnose_neff.jl — compare 3 estimators of the effective sample size (n_eff)
# on real chitosan fits, to decide which to use instead of the n÷9 heuristic.
#
# Usage:
#   julia -t 2 --project=. test/diagnose_neff.jl [file1.sxm file2.sxm ...]
#   (default: a representative subset of the chitosan benchmark)
#
# For each file, fits a chain at the batch-configured N and reports:
#   n_pixels | n_heuristic (n÷9) | n_dw (Durbin-Watson AR(1)) | n_vario (2D variogram)
#            | DW | rho | vario_range_px
#
# Decision rule: if n_dw ≈ n_vario across files, DW suffices (cheap, already
# computed). If they diverge, 2D correlation matters and the variogram is needed.

using STMSXMIO, GaussianFit2D, STMFitCore
using Statistics, Printf, LinearAlgebra, TOML

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const CONFIG = joinpath(@__DIR__, "..", "config", "chitosan.toml")

"Build the same ChainSweepConfig the batch uses (from chitosan.toml)."
function make_ccfg()
    m = TOML.parsefile(CONFIG)["model"]
    GaussianFit2D.ChainSweepConfig(
        n_min=get(m,"n_min",2), n_max=get(m,"n_max",14),
        spacing_min_nm=m["spacing_min_nm"], spacing_max_nm=m["spacing_max_nm"],
        fit_width_nm=m["fit_width_nm"], support_noise_k=m["support_noise_k"],
        support_padding_nm=m["support_padding_nm"], max_overlap=m["max_overlap"],
        global_maxtime=m["global_maxtime"], global_maxiter=m["global_maxiter"],
        sigma_parallel_min_nm=m["sigma_parallel_min_nm"], sigma_parallel_max_nm=m["sigma_parallel_max_nm"],
        sigma_perp_min_nm=m["sigma_parallel_min_nm"], sigma_perp_max_nm=m["sigma_parallel_max_nm"],
        chain_tilted_baseline=true, intelligent_sweep=true, fuse_z_bwd=true)
end

"Isotropic empirical variogram of residuals; returns (range_px, sill)."
function variogram_residuals(x, y, r; max_lag_px=10)
    # Pixel resolution (assume uniform grid): infer from nearest distinct x spacing.
    ux = sort(unique(x)); uy = sort(unique(y))
    length(ux) < 2 && return (1.0, var(r))
    dx = minimum(diff(ux))
    area = (ux[end]-ux[1]+dx) * (uy[end]-uy[1]+dx)
    n = length(r)
    # Bin pairwise |Δr|² by lag distance (in pixel units).
    gamma = zeros(max_lag_px); counts = zeros(Int, max_lag_px)
    for i in 1:n
        xi, yi, ri = x[i], y[i], r[i]
        for j in (i+1):n
            h = sqrt(((x[j]-xi)/dx)^2 + ((y[j]-yi)/dx)^2)  # lag in px
            lag = round(Int, h)
            (1 <= lag <= max_lag_px) || continue
            gamma[lag] += (r[j]-ri)^2
            counts[lag] += 1
        end
    end
    for k in 1:max_lag_px
        counts[k] > 0 && (gamma[k] /= 2*counts[k])
    end
    sill = var(r)
    # Range = first lag where gamma reaches 95% of sill (clamped to max_lag).
    a = Float64(max_lag_px)
    for k in 1:max_lag_px
        if counts[k] > 0 && gamma[k] >= 0.95*sill
            a = Float64(k); break
        end
    end
    return (a, sill, area, dx)
end

"n_eff from variogram: n * (range²) / domain_area (in consistent px² units)."
n_eff_variogram(x, y, r) = let
    n = length(r)
    a, sill, area, dx = variogram_residuals(x, y, r)
    area_px = area / dx^2
    max(10, round(Int, n * a^2 / area_px))
end

"n_eff from Durbin-Watson: n * (1-ρ)/(1+ρ) with ρ = 1 - DW/2."
function n_eff_dw(resid)
    dw, _ = STMFitCore.durbin_watson(resid)
    isnan(dw) && return (max(10, length(resid)÷9), dw, NaN)
    rho = clamp(1 - dw/2, -0.999, 0.999)
    return (max(10, round(Int, length(resid)*(1-rho)/(1+rho))), dw, rho)
end

function diagnose_file(fp, ccfg)
    img = read_sxm(fp)
    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp/neff_diag", no_plot=true)
    res, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, deepcopy(ccfg))
    # Use the best fit's residuals at the selected N.
    r = best
    # Recover the fit-mask residuals via the same path as _finalize_chain_result!.
    xs, ys, zimg, mask, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    xfit, yfit, zfit, ac, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ctx.axisctx_full, ccfg)
    pred = GaussianFit2D._chain_model_values(xfit, yfit, r.params, r.n, ac, ccfg; amp_min=r.amp_min, amp_range=r.amp_range)
    resid = zfit .- pred
    n = length(zfit)
    n_heur = max(10, n ÷ 9)
    n_dw, dw, rho = n_eff_dw(resid)
    a, sill, area, dx = variogram_residuals(xfit, yfit, resid)
    n_var = n_eff_variogram(xfit, yfit, resid)
    @printf("%-16s N=%-2d npix=%-5d | heur=%-4d dw=%-4d (DW=%.2f ρ=%.2f) | vario=%-4d (range=%dpx) | dw/heur=%.2f var/heur=%.2f\n",
        basename(fp), r.n, n, n_heur, n_dw, dw, rho, n_var, round(Int,a), n_dw/n_heur, n_var/n_heur)
    return (basename(fp), n, n_heur, n_dw, n_var, dw, rho, a)
end

const DEFAULT_FILES = ["240817_002.sxm","240817_017.sxm","240817_019.sxm","240817_043.sxm",
                       "240817_058.sxm","240817_032.sxm","240817_044.sxm","240817_066.sxm"]

files = isempty(ARGS) ? DEFAULT_FILES : ARGS
ccfg = make_ccfg()
println("=== n_eff diagnostic (chitosan) ===")
println("file             N  npix    | heuristic  DW-derived      variogram        ratios")
rows = map(f -> diagnose_file(joinpath(DATA_DIR, f), ccfg), files)
println()
# Summary stats
heurs = [r[3] for r in rows]; dws = [r[4] for r in rows]; vars = [r[5] for r in rows]
@printf("median n_heuristic = %.0f, n_dw = %.0f, n_vario = %.0f\n", median(heurs), median(dws), median(vars))
@printf("median ratio dw/heur = %.2f, vario/heur = %.2f\n", median(dws./heurs), median(vars./heurs))
@printf("median DW = %.2f, median ρ = %.3f, median vario range = %.0f px\n", median([r[6] for r in rows]), median([r[7] for r in rows]), median([r[8] for r in rows]))
