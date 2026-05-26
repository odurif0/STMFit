#!/usr/bin/env julia --project=.
using GaussianFit1D, STMMolecularFit, Statistics, Printf

filepath = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"
slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30,
    support_noise_k=2.5, support_padding_nm=0.25, output_dir="/tmp/stmfit_test", no_plot=true)
img = STMMolecularFit.read_sxm(filepath)
slide = STMMolecularFit.extract_slide(img, slide_cfg)
x, y = slide.x, slide.y

common = Dict{String,Any}(
    "filepath" => filepath, "min_spacing" => 0.35, "max_spacing" => 0.75,
    "fwhm_min" => 0.45, "fwhm_max" => 1.20, "max_overlap" => 0.60,
    "kappa_max" => 10.0, "global_maxtime" => 10.0, "global_maxiter" => 5000,
    "use_student_bic" => true, "no_show" => true, "output_dir" => "/tmp/stmfit_test",
)

# Fit N=6 with both profiles
for profile in [:gaussian, :pseudo_voigt]
    cfg = GaussianFit1D.build_config(Dict{String,Any}(common..., "peak_profile" => profile))
    r = GaussianFit1D._fit_one(6, x, y, cfg)
    
    res = y .- r.y_fit
    noise_perfit = max(std(res)*0.1, 1.4826*median(abs.(res)), 1e-12)
    noise_mad = 1.4826 * median(abs.(res .- median(res)))
    
    println("=== $(profile) ===")
    @printf("  pops: %d params, RSS=%.6f, R²=%.6f\n", length(r.popt), r.rss, r.r_squared)
    if profile == :pseudo_voigt
        @printf("  η=%.6f, η_err=%.6f\n", r.popt[end], r.perr[end])
    end
    @printf("  noise(per-fit)=%.6f, noise(MAD)=%.6f\n", noise_perfit, noise_mad)
    
    # Student NLL with per-fit noise
    nu = 4.0
    nll_perfit = sum(0.5*(nu+1) .* log1p.((res ./ noise_perfit) .^ 2 ./ nu))
    # Student NLL with MAD noise
    nll_mad = sum(0.5*(nu+1) .* log1p.((res ./ noise_mad) .^ 2 ./ nu))
    
    n_eff = length(y)
    @printf("  sBIC(per-fit)=%.1f  sBIC(MAD)=%.1f  bic(Gaussian)=%.1f\n",
            2*nll_perfit + r.n_params*log(n_eff),
            2*nll_mad + r.n_params*log(n_eff),
            r.bic)
    
    # Also compute what happens if we apply Gaussian's noise to PV
    if profile == :gaussian
        println("  Centers: $(join([@sprintf("%.3f",c) for c in GaussianFit1D._params_to_centers(r.popt,6)], " "))")
        println("  Amps: $(join([@sprintf("%.3f",GaussianFit1D._get_amplitude(r.popt,i)) for i in 0:5], " "))")
    end
end

# Cross-comparison: use Gaussian noise for both
println("\n=== Cross-comparison (shared noise) ===")
r_g = GaussianFit1D._fit_one(6, x, y, GaussianFit1D.build_config(Dict{String,Any}(common..., "peak_profile"=>:gaussian)))
r_pv = GaussianFit1D._fit_one(7, x, y, GaussianFit1D.build_config(Dict{String,Any}(common..., "peak_profile"=>:pseudo_voigt)))
# Note: pseudo_voigt with n_peaks=6 has 1+3*6+1=20 params (η added)

noise_shared = 1.4826 * median(abs.((y .- r_g.y_fit) .- median(y .- r_g.y_fit)))
println("Shared noise (from Gaussian): $(round(noise_shared, digits=6))")

nu = 4.0; n_eff = length(y)
for (label, r, k) in [("Gaussian", r_g, r_g.n_params), ("Pseudo-Voigt", r_pv, r_pv.n_params)]
    res = y .- r.y_fit
    nll = sum(0.5*(nu+1) .* log1p.((res ./ noise_shared) .^ 2 ./ nu))
    sbic = 2*nll + k*log(n_eff)
    @printf("  %-14s: k=%d, RSS=%.6f, NLL=%.3f, sBIC=%.1f\n", label, k, sum(abs2, res), nll, sbic)
end
