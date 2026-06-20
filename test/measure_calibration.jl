# measure_calibration.jl — objectively derive all objectivable calibration
# parameters from a single clean scan (no hand-tuning). Outputs a TOML-ready
# calibration, generalizable to any chain-like molecule on the same STM.
#
# Usage:
#   julia -t 2 --project=. test/measure_calibration.jl <clean_scan.sxm> [--n-lobe N]
#
# Each measured parameter is annotated with its measurement method:
#   [measured]   = derived from the scan itself (objective)
#   [principled] = derived from a physical/numerical principle with one fixed choice
#   [free]       = left to default (not objectively measurable)

using STMSXMIO, GaussianFit2D, STMFitCore
using Statistics, Printf, TOML, LinearAlgebra

const FP = ARGS[1]
const TARGET_N = findfirst(a -> a == "--n-lobe", ARGS) !== nothing ?
    parse(Int, ARGS[findfirst(a -> a == "--n-lobe", ARGS) + 1]) : nothing

img = read_sxm(FP)
ch = STMSXMIO.get_channel(img, "Z"; direction="fwd")

# ───────────────────────────────────────────────────────────────────
# 1. Noise level (1.4826·MAD of high-frequency component)  [measured]
#    Uses the standard preprocessing pipeline (flatten + smooth) to get the
#    denoised field, then MAD of (raw - smoothed) is the noise.
# ───────────────────────────────────────────────────────────────────
function measure_noise(img, ch)
    pcfg = GaussianFit2D.PatternConfig(filepath=img.filepath, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp/calib", no_plot=true)
    xs, ys, raw, z, z_smooth, unit, noise = GaussianFit2D.preprocess_channel(img, ch, pcfg)
    return noise  # already computed robustly by the pipeline (1.4826·MAD, in nm)
end

# ───────────────────────────────────────────────────────────────────
# 2. Pixel resolution  [measured]
# ───────────────────────────────────────────────────────────────────
px_nm = img.range_nm[1] / img.width  # nm/pixel along x

# ───────────────────────────────────────────────────────────────────
# Shared: preprocess + extract chain-axis profile (smoothed, in nm).
# Returns (profile_x_nm, profile_z_nm, axis, px_nm).
# ───────────────────────────────────────────────────────────────────
function chain_profile(img, ch)
    pcfg = GaussianFit2D.PatternConfig(filepath=img.filepath, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp/calib", no_plot=true)
    xs, ys, raw, z, z_smooth, unit, noise = GaussianFit2D.preprocess_channel(img, ch, pcfg)
    px = xs[2] - xs[1]
    # Weighted PCA axis on bright pixels.
    finite = isfinite.(z_smooth)
    thr = quantile(z_smooth[finite], 0.70)
    mask = finite .& (z_smooth .>= thr)
    ny, nx = size(z_smooth)
    xb = Float64[]; yb = Float64[]; zb = Float64[]
    for iy in 1:ny, ix in 1:nx
        mask[iy,ix] || continue
        push!(xb, xs[ix]); push!(yb, ys[iy]); push!(zb, z_smooth[iy,ix])
    end
    w = zb .- minimum(zb) .+ 1e-12; sw = sum(w)
    ox = sum(xb .* w)/sw; oy = sum(yb .* w)/sw
    X = hcat(xb .- ox, yb .- oy)
    _, _, V = svd(Diagonal(w ./ maximum(w)) * X)
    ax = collect(V[:,1]); ax ./= max(norm(ax), 1e-12); ax[2] < 0 && (ax .*= -1)
    t = (xb .- ox).*ax[1] + (yb .- oy).*ax[2]
    u = -(xb .- ox).*ax[2] + (yb .- oy).*ax[1]
    strip = abs.(u) .<= 3px  # narrow perpendicular strip
    ts = t[strip]; zs = zb[strip]
    # Bin at 3× pixel resolution (smoother profile).
    binw = 3px
    tlo, thi = minimum(ts), maximum(ts)
    edges = collect(tlo:binw:thi); nbin = length(edges)-1
    prof = zeros(nbin); cnt = zeros(Int,nbin)
    for k in eachindex(ts)
        b = min(nbin, max(1, Int(fld(ts[k]-tlo, binw))+1)); prof[b]+=zs[k]; cnt[b]+=1
    end
    prof ./= max.(cnt,1)
    cents = edges[1:nbin] .+ binw/2
    return cents, prof, ax, px
end

# ───────────────────────────────────────────────────────────────────
# 3. FWHM range of lobes  [measured]
#    Detects peaks in the smoothed axis profile, fits a half-max width.
# ───────────────────────────────────────────────────────────────────
function measure_fwhm_range(img, ch)
    cents, prof, ax, px = chain_profile(img, ch)
    isempty(prof) && return (0.3, 1.0)
    thr_peak = median(prof) + 2*std(prof)
    fwhms = Float64[]
    nbin = length(prof)
    for k in 2:(nbin-1)
        prof[k] > thr_peak && prof[k] >= prof[k-1] && prof[k] >= prof[k+1] || continue
        lo, hi = max(1, k-3), min(nbin, k+3)
        win = lo:hi; p = prof[win]
        amp = prof[k] - minimum(p); amp <= 0 && continue
        half = minimum(p) + amp/2
        above = win[p .>= half]; isempty(above) && continue
        fwhm = cents[min(last(above)+1,nbin)] - cents[first(above)]
        fwhm > 2px && fwhm < 5.0 && push!(fwhms, fwhm)
    end
    isempty(fwhms) && return (0.3, 1.0)
    return (quantile(fwhms, 0.05), quantile(fwhms, 0.95))
end

# ───────────────────────────────────────────────────────────────────
# 4. Repeat spacing (median distance between adjacent peaks)  [measured]
# ───────────────────────────────────────────────────────────────────
function measure_spacing(img, ch)
    cents, prof, ax, px = chain_profile(img, ch)
    isempty(prof) && return 0.5
    thr_peak = median(prof) + 2*std(prof)
    peaks = Float64[]
    nbin = length(prof)
    for k in 2:(nbin-1)
        if prof[k] > thr_peak && prof[k] >= prof[k-1] && prof[k] >= prof[k+1]
            push!(peaks, cents[k])
        end
    end
    length(peaks) >= 2 || return 0.5
    return median(diff(sort(peaks)))
end

# ───────────────────────────────────────────────────────────────────
# 5. Spatial correlation range (autocorrelation 1/e on full image)  [measured]
# ───────────────────────────────────────────────────────────────────
function autocorr_range(img, ch; max_lag=40)
    pcfg = GaussianFit2D.PatternConfig(filepath=img.filepath, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp/calib", no_plot=true)
    xs, ys, raw, z, z_smooth, unit, noise = GaussianFit2D.preprocess_channel(img, ch, pcfg)
    px = xs[2]-xs[1]
    zsub = z_smooth[1:4:end, 1:4:end]
    zc = zsub .- mean(zsub)
    var_px = sum(abs2, zc) / length(zc)
    var_px > 0 || return NaN
    ny, nx = size(zc)
    sums = zeros(max_lag); counts = zeros(Int, max_lag)
    for iy in 1:ny, ix in 1:nx
        v = zc[iy,ix]
        for dy in 0:min(max_lag, ny-iy), dx in 0:min(max_lag, nx-ix)
            (dy==0 && dx==0) && continue
            h = round(Int, sqrt(dy^2+dx^2))
            (1<=h<=max_lag) || continue
            sums[h] += v*zc[iy+dy,ix+dx]; counts[h]+=1
        end
    end
    rho = [counts[k]>0 ? (sums[k]/counts[k])/var_px : NaN for k in 1:max_lag]
    for k in eachindex(rho)
        isnan(rho[k]) && continue
        rho[k] <= exp(-1) && return Float64(k) * 4  # ×4 (subsample)
    end
    return Float64(max_lag * 4)
end

# ───────────────────────────────────────────────────────────────────
# Run all measurements
# ───────────────────────────────────────────────────────────────────
println("┌─ Objective calibration from: $(basename(FP))")
println("│ image: $(img.width)×$(img.height), range $(round(img.range_nm[1],digits=2))×$(round(img.range_nm[2],digits=2)) nm")

noise_nm = measure_noise(img, ch)
px = px_nm
fwhm_lo, fwhm_hi = measure_fwhm_range(img, ch)
spacing = measure_spacing(img, ch)
range_corr_px = autocorr_range(img, ch)

# ───────────────────────────────────────────────────────────────────
# Derive calibration from measurements  [principled]
# ───────────────────────────────────────────────────────────────────
sigma_lo = fwhm_lo / 2.355
sigma_hi = fwhm_hi / 2.355
spacing_min = 0.7 * spacing
spacing_max = 1.3 * spacing
fit_width = sigma_lo  # tube half-width ≈ min lobe half-width
support_min = 3 * spacing  # at least 3 repeats to call it a chain
n_max = round(Int, img.range_nm[1] / spacing_min + 2)  # generous, from chain axis length

println("│")
println("│ ── measured ──")
@printf("│ noise σ           = %.4f nm  [1.4826·MAD of HF band]\n", noise_nm)
@printf("│ pixel resolution  = %.4f nm/px\n", px)
@printf("│ FWHM range        = [%.3f, %.3f] nm  [5/95 pct of peak fits]\n", fwhm_lo, fwhm_hi)
@printf("│ repeat spacing    = %.3f nm  [median peak-to-peak]\n", spacing)
@printf("│ corr. range (1/e) = %.1f px (%.2f nm)\n", range_corr_px, range_corr_px*px)
println("│")
println("│ ── derived (principled) ──")
@printf("│ sigma_parallel    = [%.3f, %.3f] nm  [FWHM / 2.355]\n", sigma_lo, sigma_hi)
@printf("│ spacing           = [%.3f, %.3f] nm  [±30%% around measured]\n", spacing_min, spacing_max)
@printf("│ fit_width_nm      = %.3f nm  [= σ_min (tube half-width)]\n", fit_width)
@printf("│ support_min_length= %.3f nm  [3× spacing]\n", support_min)
@printf("│ n_max             = %d  [axis_length / spacing_min + 2]\n", n_max)
println("└─")

# ───────────────────────────────────────────────────────────────────
# Emit TOML
# ───────────────────────────────────────────────────────────────────
calib = Dict(
    "model" => Dict(
        "sigma_parallel_min_nm" => sigma_lo,
        "sigma_parallel_max_nm" => sigma_hi,
        "spacing_min_nm" => spacing_min,
        "spacing_max_nm" => spacing_max,
        "fit_width_nm" => fit_width,
        "support_min_length_nm" => support_min,
        "n_max" => n_max,
        "max_overlap" => 0.60,                 # [principled: Gaussian pair-overlap]
        "support_noise_k" => 2.5,              # [principled: SNR threshold k·σ]
        "support_padding_nm" => round(fit_width, digits=3),  # [derived: ≈ tube half-width]
        "global_maxtime" => 10.0,              # [free: optimizer budget]
        "global_maxiter" => 10000,
        "chain_tilted_baseline" => true,       # [free: model form]
        "chain_circular_sigmas" => false,
        "selection_criterion" => "gcv",        # [principled: valid under spatial corr.]
        "selection_policy" => "gcv_with_robust_aicc_guard",
    ),
    "selection" => Dict(
        "gcv_ambiguity_rel_threshold" => 0.05, # [measured robust on benchmark]
        "robust_guard_nu" => 8.0,
    ),
    "preprocessing" => Dict(
        "channel" => "Z",                      # [free: depends on acquisition]
        "direction" => "fwd",
        "stride" => 1,
        "flatten" => "plane+rows",             # [principled: STM scan-line correction]
        "smooth_radius_px" => 1,
    ),
)
out = replace(basename(FP), ".sxm" => "_calibration.toml")
open(out, "w") do io
    println(io, "# Auto-derived calibration from $(basename(FP))")
    println(io, "# Measured: noise=$(round(noise_nm,digits=4))nm, px=$(round(px,digits=4))nm/px,")
    println(io, "#           FWHM=[$(round(fwhm_lo,digits=3)), $(round(fwhm_hi,digits=3))]nm, spacing=$(round(spacing,digits=3))nm")
    println(io, "#           corr_range=$(round(range_corr_px,digits=1))px ($(round(range_corr_px*px,digits=2))nm)")
    TOML.print(io, calib)
end
println("\n→ wrote $out")
