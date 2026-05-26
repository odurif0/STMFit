#!/usr/bin/env julia --project=.
# Test the fix: run_model_comparison sweep on 002
using GaussianFit1D, STMMolecularFit, Statistics
using Printf

filepath = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"
slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30,
    support_noise_k=2.5, support_padding_nm=0.25, output_dir="/tmp/stmfit_test", no_plot=true)
img = STMMolecularFit.read_sxm(filepath)
slide = STMMolecularFit.extract_slide(img, slide_cfg)
x, y = slide.x, slide.y
println("Support: $(round(slide.support_length_nm,digits=2)) nm, $(length(x)) pts\n")

cfg = GaussianFit1D.build_config(Dict{String,Any}(
    "filepath" => filepath, "min_spacing" => 0.35, "max_spacing" => 0.75,
    "fwhm_min" => 0.45, "fwhm_max" => 1.20, "max_overlap" => 0.60,
    "kappa_max" => 10.0, "kappa_weight" => 1.0,
    "global_maxtime" => 10.0, "global_maxiter" => 5000,
    "use_student_bic" => true, "no_show" => true,
    "peak_profile" => :gaussian, "output_dir" => "/tmp/stmfit_test",
))

results = GaussianFit1D.run_model_comparison(x, y, cfg)
println()

best = GaussianFit1D.best_result(GaussianFit1D.FitRunResult(x, y, results, cfg))
println("BEST: N=$(best.n_peaks), sBIC=$(round(best.student_bic,digits=1)), R²=$(round(best.r_squared,digits=4))")
if best.residual_diagnostics !== nothing
    rd = best.residual_diagnostics
    println("  DW=$(round(rd.durbin_watson,digits=3)), runs=$(rd.runs_n), RMS=$(round(rd.residual_rms,digits=6))")
end
