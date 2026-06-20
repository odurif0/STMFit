# diagnose_fullimg_autocorr.jl — check whether the spatial correlation range
# is OBSERVABLE in the full preprocessed image (where it wasn't in the small
# fit window). This determines whether the geostatistical n_eff method is
# feasible: estimate range from the full image, apply to the fit window.
#
# Usage: julia -t 1 --project=. test/diagnose_fullimg_autocorr.jl [file.sxm ...]

using STMSXMIO, GaussianFit2D, STMFitCore
using Statistics, Printf, TOML

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const CONFIG = joinpath(@__DIR__, "..", "config", "chitosan.toml")

"2D isotropic autocorrelation ρ(h) by integer lag, computed directly (no FFT).
ρ(h) = mean over pairs (i,i+h) of zc[i]·zc[i+h], normalized by the variance
mean(zc²). O(n_pixels · max_lag²) — fine for max_lag ≤ 40."
function autocorr_by_lag(z::Matrix{Float64}; max_lag=40)
    zc = z .- mean(z)
    ny, nx = size(z)
    sums = zeros(max_lag); counts = zeros(Int, max_lag)
    for iy in 1:ny, ix in 1:nx
        v = zc[iy, ix]
        for dy in 0:min(max_lag, ny-iy), dx in 0:min(max_lag, nx-ix)
            (dy == 0 && dx == 0) && continue
            h = round(Int, sqrt(dy^2 + dx^2))
            (1 <= h <= max_lag) || continue
            sums[h] += v * zc[iy+dy, ix+dx]
            counts[h] += 1
        end
    end
    # Normalize each lag by its pair count; then normalize so ρ(0)=1 by dividing
    # by the per-pixel variance estimate (mean of zc²).
    var_per_pixel = sum(abs2, zc) / length(zc)
    var_per_pixel > 0 || return fill(NaN, max_lag)
    return [counts[k] > 0 ? (sums[k]/counts[k]) / var_per_pixel : NaN for k in 1:max_lag]
end

"Range = first lag where ρ drops below 1/e ≈ 0.368 (exponential model)."
function range_1e(rho::Vector{Float64})
    for k in eachindex(rho)
        isnan(rho[k]) && continue
        rho[k] <= exp(-1) && return Float64(k)
    end
    return Float64(length(rho))  # never reaches 1/e within window
end

function diagnose(fp)
    img = read_sxm(joinpath(DATA_DIR, fp))
    ch = STMSXMIO.get_channel(img, "Z"; direction="fwd")
    m = TOML.parsefile(CONFIG)
    pp = m["preprocessing"]
    pcfg = GaussianFit2D.PatternConfig(filepath=joinpath(DATA_DIR,fp), channel="Z", direction="fwd",
        stride=pp["stride"], flatten=pp["flatten"], smooth_radius_px=pp["smooth_radius_px"],
        output_dir="/tmp/neff_full", no_plot=true)
    xs, ys, raw, z, z_smooth, unit, noise = GaussianFit2D.preprocess_channel(img, ch, pcfg)
    # Subsample every 4th pixel — enough to estimate the range, 16× faster.
    zsub = z_smooth[1:4:end, 1:4:end]
    rho = autocorr_by_lag(zsub; max_lag=min(40, size(zsub)...))
    a = range_1e(rho)
    # ρ at first few lags
    r1 = length(rho) >= 1 ? rho[1] : NaN
    r3 = length(rho) >= 3 ? rho[3] : NaN
    r10 = length(rho) >= 10 ? rho[10] : NaN
    # n_eff in a fit window of ~10×10 px: independent points = window_area / (pi·range²)
    win_px = 10.0  # approximate fit-window side in px
    n_fit = win_px^2
    a_indep = pi * a^2
    n_eff_win = a_indep > 0 ? n_fit / a_indep : NaN
    @printf("%-16s img=%dx%d | ρ(1)=%.2f ρ(3)=%.2f ρ(10)=%.2f | range(1/e)=%5.1f px | n_eff(10px win)=%.1f\n",
        fp, size(z_smooth)..., r1, r3, r10, a, n_eff_win)
    return a, r1, r3, r10
end

const DEFAULT = ["240817_002.sxm","240817_043.sxm","240817_058.sxm","240817_017.sxm"]
files = isempty(ARGS) ? DEFAULT : ARGS
println("=== Full-image 2D autocorrelation diagnostic ===")
println("(range observable if ρ drops below 1/e≈0.37 within ~40px)")
for f in files
    try
        diagnose(f)
    catch e
        println("$f: ERROR $(typeof(e))")
    end
end
