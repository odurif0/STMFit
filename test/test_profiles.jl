#!/usr/bin/env julia
# Compare Gaussian vs Lorentzian vs Pseudo-Voigt profiles on 1D slide data.
using STMMolecularFit, GaussianFit1D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_003.sxm"
const OUT = "results/profile_compare"
mkpath(OUT)

scfg = STMMolecularFit.SlideConfig(
    width_nm=0.30, support_threshold_fraction=0.20,
    support_noise_k=2.5, support_padding_nm=0.20,
    output_dir=OUT, no_plot=true)

img = STMMolecularFit.read_sxm(FILE)
slide = STMMolecularFit.extract_slide(img, scfg)

for prof in (:gaussian, :lorentzian, :pseudo_voigt)
    println("\n═══ Profile: $prof ═══")
    fcfg = STMMolecularFit.FitSlideConfig(
        min_spacing=0.35, max_spacing=0.75, max_overlap=0.6,
        peak_profile=prof, output_dir=OUT)
    fit = STMMolecularFit.fit_slide(slide, fcfg)
    best = GaussianFit1D.best_result(fit.fit_run)

    rd = best.residual_diagnostics
    @printf("  N=%d  sBIC=%.1f  R²=%.4f  κ=%.1f\n",
        best.n_peaks, best.student_bic, best.r_squared, best.kappa_max_adj)
    if rd !== nothing
        @printf("  DW=%.3f (p=%.4f)  Runs=%d/%.0f (p=%.4f)  RMS=%.5f\n",
            rd.durbin_watson, rd.durbin_watson_p,
            rd.runs_n, rd.runs_expected, rd.runs_p,
            rd.residual_rms)
    end
    if prof == :pseudo_voigt && best.popt !== nothing
        np1 = length(best.popt)
        eta = best.popt[end]
        eta_err = (best.perr !== nothing && length(best.perr) == np1) ? best.perr[end] : NaN
        @printf("  η=%.4f ± %.4f\n", eta, eta_err)
    end

    # Decode centers from popt: [y0, A0, mu0, s0, A1, d1, s1, ...]
    p = best.popt
    n = best.n_peaks
    cents = Float64[p[3]]  # mu0
    for i in 1:(n-1)
        push!(cents, cents[end] + p[3 + 3*i])  # mu_{i} = mu_{i-1} + delta_i
    end
    if length(cents) > 1
        sp = diff(cents)
        @printf("  Centers: %s nm\n", join([@sprintf("%.3f", c) for c in cents], ", "))
        @printf("  Spacings: %s nm\n", join([@sprintf("%.3f", s) for s in sp], ", "))
    end
end
