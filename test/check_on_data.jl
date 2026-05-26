#!/usr/bin/env julia --project=.
# Run 1D fit on chitosan data: compare Gaussian vs Lorentzian vs Pseudo-Voigt
# with residual diagnostics and covariance
using GaussianFit1D, STMMolecularFit
using STMFitCore: ResidualDiagnostics
using Printf

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

    # Load via STMMolecularFit (extract 1D slide from SXM)
    pcfg = STMMolecularFit.PreprocessConfig(channel="Z", direction="fwd", stride=1,
                            flatten="plane+rows", smooth_radius_px=1)
    slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30,
                            support_noise_k=2.5, support_padding_nm=0.25,
                            output_dir=RESULTS_DIR, no_plot=true)
    
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

    # Run all 3 profiles
    profiles = [:gaussian, :lorentzian, :pseudo_voigt]
    results = Dict{Symbol, Any}()

    for profile in profiles
        cfg = GaussianFit1D.build_config(Dict{String,Any}(
            "filepath" => filepath,
            "min_spacing" => 0.35,
            "max_spacing" => 0.75,
            "fwhm_min" => 0.45,
            "fwhm_max" => 1.20,
            "max_overlap" => 0.60,
            "kappa_max" => 10.0,
            "kappa_weight" => 1.0,
            "global_maxtime" => 8.0,
            "global_maxiter" => 5000,
            "use_student_bic" => true,
            "no_show" => true,
            "peak_profile" => profile,
            "output_dir" => RESULTS_DIR,
        ))
        fr = GaussianFit1D.run_fit(x, y, cfg; save_cache=false, verbose=false)
        results[profile] = fr
    end

    # Print comparison table
    println("\n  " * rpad("Profile", 16) * rpad("N", 4) * rpad("sBIC", 10) *
            rpad("R²", 10) * rpad("η", 8) * rpad("DW", 8) * rpad("DW_p", 8) *
            rpad("Runs", 6) * rpad("RMS", 10) * rpad("κ_max", 8))
    println("  " * repeat("-", 88))

    for profile in profiles
        fr = results[profile]
        isempty(fr.all_results) && continue
        best = GaussianFit1D.best_result(fr)
        best === nothing && continue

        eta_str = profile == :pseudo_voigt ? @sprintf("%.3f", best.popt[end]) : "---"
        
        rd = best.residual_diagnostics
        dw_str = rd !== nothing && isfinite(rd.durbin_watson) ? @sprintf("%.3f", rd.durbin_watson) : "NaN"
        dwp_str = rd !== nothing && isfinite(rd.durbin_watson_p) ? @sprintf("%.3f", rd.durbin_watson_p) : "NaN"
        runs_str = rd !== nothing ? string(rd.runs_n) : "?"
        rms_str = rd !== nothing ? @sprintf("%.6f", rd.residual_rms) : "?"

        println("  " * rpad(string(profile), 16) * rpad(string(best.n_peaks), 4) *
                rpad(@sprintf("%.1f", best.student_bic), 10) *
                rpad(@sprintf("%.4f", best.r_squared), 10) *
                rpad(eta_str, 8) *
                rpad(dw_str, 8) * rpad(dwp_str, 8) *
                rpad(runs_str, 6) * rpad(rms_str, 10) *
                rpad(@sprintf("%.1f", best.kappa_max_adj), 8))
    end

    # Gaussian best: show center-center correlation
    best_g = GaussianFit1D.best_result(results[:gaussian])
    if best_g !== nothing && best_g.center_center_corr !== nothing && best_g.n_peaks > 1
        cc = best_g.center_center_corr
        max_corr = 0.0
        for i in 1:size(cc,1), j in (i+1):size(cc,2)
            abs(cc[i,j]) > abs(max_corr) && (max_corr = cc[i,j])
        end
        println("  Max center-center corr (Gaussian): $(round(max_corr, digits=3))")
    end
end

println("\nDone.")
