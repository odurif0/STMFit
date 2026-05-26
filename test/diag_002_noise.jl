#!/usr/bin/env julia --project=.
# Test: fixed global noise for Student-t BIC (fix for 1D)
using GaussianFit1D, STMMolecularFit, Statistics
using Printf

filepath = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"

slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30,
    support_noise_k=2.5, support_padding_nm=0.25, output_dir="/tmp/stmfit_test", no_plot=true)
img = STMMolecularFit.read_sxm(filepath)
slide = STMMolecularFit.extract_slide(img, slide_cfg)
x, y = slide.x, slide.y

# Fit N=4,5,6
results = GaussianFit1D.FitResult[]
for n in [4,5,6]
    cfg = GaussianFit1D.build_config(Dict{String,Any}(
        "filepath" => filepath, "min_spacing" => 0.35, "max_spacing" => 0.75,
        "fwhm_min" => 0.45, "fwhm_max" => 1.20, "max_overlap" => 0.60,
        "kappa_max" => 10.0, "kappa_weight" => 1.0,
        "global_maxtime" => 10.0, "global_maxiter" => 5000,
        "use_student_bic" => true, "no_show" => true,
        "peak_profile" => :gaussian, "output_dir" => "/tmp/stmfit_test",
    ))
    r = GaussianFit1D._fit_one(n, x, y, cfg)
    r !== nothing && push!(results, r)
end

# --- CURRENT (per-fit) noise ---
println("\n=== CURRENT: per-fit noise estimate ===")
println(rpad("N",4), rpad("sBIC",10), rpad("RSS",12), rpad("noise",10), rpad("R²",8))
for r in results
    res = y .- r.y_fit
    noise_cur = max(std(res)*0.1, median(abs.(res))*1.4826, 1e-12)
    nu = 4.0
    nll = sum(0.5 * (nu+1) .* log1p.((res ./ noise_cur) .^ 2 ./ nu))
    n_eff = length(y)
    sbic = 2*nll + r.n_params * log(n_eff)
    @printf("%4d %10.1f %12.6f %10.6f %8.4f\n", r.n_peaks, sbic, r.rss, noise_cur, r.r_squared)
end

# --- FIXED noise (MAD of N_max residuals) ---
println("\n=== FIXED noise: MAD of N=$(results[end].n_peaks) residuals ===")
r_ref = results[argmin([r.rss for r in results])]
res_ref = y .- r_ref.y_fit
noise_fixed = 1.4826 * median(abs.(res_ref .- median(res_ref)))
println("  Reference: N=$(r_ref.n_peaks), RSS=$(round(r_ref.rss,digits=6))")
println("  Fixed noise = $(round(noise_fixed, digits=6))")

println(rpad("N",4), rpad("sBIC(fixed)",12), rpad("RSS",12), rpad("resid/noise",12), rpad("R²",8))
for r in results
    res = y .- r.y_fit
    nu = 4.0
    nll = sum(0.5 * (nu+1) .* log1p.((res ./ noise_fixed) .^ 2 ./ nu))
    n_eff = length(y)
    sbic = 2*nll + r.n_params * log(n_eff)
    @printf("%4d %12.1f %12.6f %12.3f %8.4f\n", r.n_peaks, sbic, r.rss, sqrt(mean(res .^ 2))/noise_fixed, r.r_squared)
end

# --- ALTERNATIVE: MAD of differences (data-based, fit-independent) ---
println("\n=== FIXED noise: MAD of diff(y) ===")
noise_diff = 1.4826 * median(abs.(diff(y) .- median(diff(y)))) / sqrt(2)
println("  noise_diff = $(round(noise_diff, digits=6))")
for r in results
    res = y .- r.y_fit
    nu = 4.0
    nll = sum(0.5 * (nu+1) .* log1p.((res ./ noise_diff) .^ 2 ./ nu))
    n_eff = length(y)
    sbic = 2*nll + r.n_params * log(n_eff)
    @printf("%4d %12.1f %12.6f\n", r.n_peaks, sbic, r.rss)
end
