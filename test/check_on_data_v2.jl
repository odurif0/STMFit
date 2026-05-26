#!/usr/bin/env julia --project=.
# Run 1D fit on chitosan data: compare profiles with SHARED noise
using GaussianFit1D, STMMolecularFit, Statistics, Printf

FILES = [
    "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm",
    "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_007.sxm",
    "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_014.sxm",
]
RESULTS_DIR = "/tmp/stmfit_test/results"

for filepath in FILES
    fname = basename(filepath)
    println("\n" * "="^90)
    println("  $fname")
    println("="^90)

    slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30,
        support_noise_k=2.5, support_padding_nm=0.25, output_dir=RESULTS_DIR, no_plot=true)
    slide = nothing
    try
        img = STMMolecularFit.read_sxm(filepath)
        slide = STMMolecularFit.extract_slide(img, slide_cfg)
    catch e
        println("  SKIP: slide extraction failed ($e)")
        continue
    end
    slide === nothing && continue
    x, y = slide.x, slide.y
    println("  Support: $(round(slide.support_length_nm, digits=2)) nm, $(length(x)) points")

    # Run Gaussian FIRST (establishes the noise reference)
    noise_ref = NaN  # will be populated by Gaussian sweep
    fr_g = GaussianFit1D.run_fit(x, y, GaussianFit1D.build_config(Dict{String,Any}(
        "filepath" => filepath, "min_spacing" => 0.35, "max_spacing" => 0.75,
        "fwhm_min" => 0.45, "fwhm_max" => 1.20, "max_overlap" => 0.60,
        "kappa_max" => 10.0, "kappa_weight" => 1.0,
        "global_maxtime" => 8.0, "global_maxiter" => 5000,
        "use_student_bic" => true, "no_show" => true,
        "peak_profile" => :gaussian, "output_dir" => RESULTS_DIR,
    )); save_cache=false, verbose=false)
    # Gaussian sweep auto-populates noise_estimate in the config
    noise_ref = fr_g.cfg.noise_estimate
    best_g = !isempty(fr_g.all_results) ? GaussianFit1D.best_result(fr_g) : nothing

    # Run PV with Gaussian's noise
    fr_pv = GaussianFit1D.run_fit(x, y, GaussianFit1D.build_config(Dict{String,Any}(
        "filepath" => filepath, "min_spacing" => 0.35, "max_spacing" => 0.75,
        "fwhm_min" => 0.45, "fwhm_max" => 1.20, "max_overlap" => 0.60,
        "kappa_max" => 10.0, "kappa_weight" => 1.0,
        "global_maxtime" => 8.0, "global_maxiter" => 5000,
        "use_student_bic" => true, "no_show" => true,
        "peak_profile" => :pseudo_voigt, "output_dir" => RESULTS_DIR,
        "noise_estimate" => noise_ref,
    )); save_cache=false, verbose=false)
    best_pv = !isempty(fr_pv.all_results) ? GaussianFit1D.best_result(fr_pv) : nothing

    # Run Lorentzian with Gaussian's noise
    fr_lor = GaussianFit1D.run_fit(x, y, GaussianFit1D.build_config(Dict{String,Any}(
        "filepath" => filepath, "min_spacing" => 0.35, "max_spacing" => 0.75,
        "fwhm_min" => 0.45, "fwhm_max" => 1.20, "max_overlap" => 0.60,
        "kappa_max" => 10.0, "kappa_weight" => 1.0,
        "global_maxtime" => 8.0, "global_maxiter" => 5000,
        "use_student_bic" => true, "no_show" => true,
        "peak_profile" => :lorentzian, "output_dir" => RESULTS_DIR,
        "noise_estimate" => noise_ref,
    )); save_cache=false, verbose=false)
    best_lor = !isempty(fr_lor.all_results) ? GaussianFit1D.best_result(fr_lor) : nothing

    println("\n  Shared noise: $(round(noise_ref, digits=6)) (from Gaussian sweep)")
    println("  " * rpad("Profile", 16) * rpad("N", 5) * rpad("sBIC", 10) *
            rpad("R²", 9) * rpad("η", 8) * rpad("DW", 8) * rpad("Runs", 6) *
            rpad("RMS", 10) * rpad("κ", 6))
    println("  " * repeat("-", 78))

    for (profile, best) in [("gaussian", best_g), ("pseudo_voigt", best_pv), ("lorentzian", best_lor)]
        if best === nothing
            println("  " * rpad(profile, 16) * "FAILED")
            continue
        end
        eta_str = profile == "pseudo_voigt" ? @sprintf("%.4f", best.popt[end]) : "---"
        rd = best.residual_diagnostics
        dw = rd !== nothing && isfinite(rd.durbin_watson) ? @sprintf("%.3f", rd.durbin_watson) : "NaN"
        runs = rd !== nothing ? string(rd.runs_n) : "?"
        rms = rd !== nothing ? @sprintf("%.6f", rd.residual_rms) : "?"
        println("  " * rpad(profile, 16) * rpad(string(best.n_peaks), 5) *
                rpad(@sprintf("%.1f", best.student_bic), 10) *
                rpad(@sprintf("%.4f", best.r_squared), 9) *
                rpad(eta_str, 8) * rpad(dw, 8) * rpad(runs, 6) *
                rpad(rms, 10) * rpad(@sprintf("%.1f", best.kappa_max_adj), 6))
    end
end
println("\nDone.")
