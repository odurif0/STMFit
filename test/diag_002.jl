#!/usr/bin/env julia --project=.
# Diagnostic: why does 002 prefer N=4 over N=6 despite R²=0.997?
using GaussianFit1D, STMMolecularFit, Statistics
using Printf

filepath = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"

# Extract 1D
slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30,
    support_threshold_fraction=0.20, support_noise_k=2.5,
    support_padding_nm=0.20, output_dir="/tmp/stmfit_test", no_plot=true)
img = STMMolecularFit.read_sxm(filepath)
slide = STMMolecularFit.extract_slide(img, slide_cfg)
x, y = slide.x, slide.y
println("Support: $(round(slide.support_length_nm,digits=2)) nm, $(length(x)) pts")
println("1D noise estimate: $(round(slide.noise_1d,digits=6))")

# Force-fit n=4 and n=6 with Gaussian, verbose mode to see details
for n in [4, 5, 6]
    cfg = GaussianFit1D.build_config(Dict{String,Any}(
        "filepath" => filepath,
        "min_spacing" => 0.35,
        "max_spacing" => 0.75,
        "fwhm_min" => 0.45,
        "fwhm_max" => 1.20,
        "max_overlap" => 0.60,
        "kappa_max" => 8.0,
        "kappa_weight" => 1.0,
        "global_maxtime" => 10.0,
        "global_maxiter" => 5000,
        "use_student_bic" => true,
        "no_show" => true,
        "peak_profile" => :gaussian,
        "output_dir" => "/tmp/stmfit_test",
        "noise_estimate" => slide.noise_1d,  # use 1D noise
    ))
    println("\n--- Force-fitting N=$n ---")
    r = GaussianFit1D._fit_one(n, x, y, cfg)
    if r !== nothing
        nll_g = r.rss / length(x) * log(r.rss / length(x))  # approx Gaussian NLL per point
        actual_rss = r.rss
        res = y .- r.y_fit
        noise_est = max(std(res)*0.1, median(abs.(res))*1.4826, 1e-12)
        nu = 4.0
        student_nll = sum(0.5 * (nu + 1) .* log1p.((res ./ noise_est).^2 ./ nu))
        
        @printf("  sBIC=%.1f  R²=%.4f  RSS=%.6f  n_params=%d\n",
                r.student_bic, r.r_squared, r.rss, r.n_params)
        @printf("  Student NLL=%.3f  Gaussian BIC=%.1f  noise=%.6f\n",
                student_nll, r.bic, noise_est)
        
        # Check spacings
        centers = GaussianFit1D._params_to_centers(r.popt, n)
        if n > 1
            spac = diff(centers)
            @printf("  Centers: %s\n", join([@sprintf("%.3f",c) for c in centers], " "))
            @printf("  Spacings: %s\n", join([@sprintf("%.3f",s) for s in spac], " "))
            @printf("  Mean spacing: %.3f, min: %.3f, max: %.3f\n",
                    mean(spac), minimum(spac), maximum(spac))
        end
        
        # Residual diagnostics
        if r.residual_diagnostics !== nothing
            rd = r.residual_diagnostics
            @printf("  DW=%.3f  RMS=%.6f  runs=%d\n",
                    rd.durbin_watson, rd.residual_rms, rd.runs_n)
        end
        
        # Check amplitudes
        amps = [GaussianFit1D._get_amplitude(r.popt, i) for i in 0:(n-1)]
        @printf("  Amps: %s\n", join([@sprintf("%.3f",a) for a in amps], " "))
        @printf("  Min amp fraction: %.3f of data max=%.3f\n",
                minimum(amps)/maximum(y), maximum(y))
    else
        println("  FAILED")
    end
end
