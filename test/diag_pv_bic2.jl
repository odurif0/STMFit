#!/usr/bin/env julia --project=.
using GaussianFit1D, STMMolecularFit, Statistics, Printf

filepath = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"
slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30, support_threshold_fraction=0.20,
    support_noise_k=2.5, support_padding_nm=0.20, output_dir="/tmp/stmfit_test", no_plot=true)
img = STMMolecularFit.read_sxm(filepath)
slide = STMMolecularFit.extract_slide(img, slide_cfg)
x, y = slide.x, slide.y

common = Dict{String,Any}("filepath"=>filepath,"min_spacing"=>0.35,"max_spacing"=>0.75,
    "fwhm_min"=>0.45,"fwhm_max"=>1.20,"max_overlap"=>0.60,"kappa_max"=>8.0,
    "global_maxtime"=>10.0,"global_maxiter"=>5000,"use_student_bic"=>true,
    "no_show"=>true,"output_dir"=>"/tmp/stmfit_test")

# Gaussian baseline (reference for noise)
r_g = GaussianFit1D._fit_one(6, x, y,
    GaussianFit1D.build_config(Dict{String,Any}(common..., "peak_profile"=>:gaussian)))
noise_shared = 1.4826 * median(abs.((y .- r_g.y_fit) .- median(y .- r_g.y_fit)))
println("Gaussian N=6: RSS=$(round(r_g.rss, digits=6)), noise=$(round(noise_shared, digits=6))")

# Pseudo-Voigt N=6
r_pv = GaussianFit1D._fit_one(6, x, y,
    GaussianFit1D.build_config(Dict{String,Any}(common..., "peak_profile"=>:pseudo_voigt)))
println("PV N=6: RSS=$(round(r_pv.rss, digits=6)), eta=$(round(r_pv.popt[end], digits=6)), n_params=$(r_pv.n_params)")

# Compare with SHARED noise
nu = 4.0; neff = length(y)
for (label, r) in [("Gaussian", r_g), ("PV", r_pv)]
    res = y .- r.y_fit
    nll = sum(0.5 * (nu + 1) .* log1p.((res ./ noise_shared) .^ 2 ./ nu))
    sbic = 2 * nll + r.n_params * log(neff)
    gbic = neff * log(r.rss / neff) + r.n_params * log(neff)
    @printf("%-10s: k=%d  RSS=%.6f  NLL=%.2f  sBIC=%.1f  gBIC=%.1f\n",
            label, r.n_params, r.rss, nll, sbic, gbic)
end

println("\n-> PV loses on BOTH criteria (higher RSS + 1 extra param).")
