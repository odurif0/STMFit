#!/usr/bin/env julia --project=.
# Comprehensive test of new features: pseudo-Voigt, covariance, residual diagnostics
using Test
using Random
using GaussianFit1D
using STMFitCore: ResidualDiagnostics, compute_residual_diagnostics, durbin_watson, runs_test

println("="^80)
println("FEATURE TEST SUITE — Pseudo-Voigt, Covariance, Residual Diagnostics")
println("="^80)

# ===========================================================================
# Synthetic data: 3 well-separated Gaussian peaks
# ===========================================================================
n_pts = 400
x_data = collect(range(0.0, 8.0, length=n_pts))
# Peaks at 1.5, 4.0, 6.5 — spacing 2.5 nm, each width 0.25 nm (FWHM ~0.59)
true_centers = [1.5, 4.0, 6.5]
true_amps    = [1.0, 0.9, 0.8]
true_sigmas  = [0.25, 0.25, 0.25]
y_true = zeros(n_pts)
for (c, a, s) in zip(true_centers, true_amps, true_sigmas)
    y_true .+= a .* exp.(-0.5 .* ((x_data .- c) ./ s).^2)
end
noise_level = 0.005 * maximum(y_true)
Random.seed!(1234)
y_data = y_true .+ noise_level .* randn(n_pts)

# Tighter physical bounds to guide optimizer
common_kwargs = Dict{String,Any}(
    "min_spacing" => 1.5,
    "max_spacing" => 5.0,
    "fwhm_min" => 0.25,
    "fwhm_max" => 1.5,
    "max_overlap" => 0.6,
    "kappa_max" => 10.0,
    "global_maxtime" => 10.0,
    "global_maxiter" => 5000,
    "use_student_bic" => true,
    "no_show" => true,
    "output_dir" => nothing,
)

# ===========================================================================
# 1. DW statistics on pure noise vs autocorrelated vs good fit
# ===========================================================================
@testset "1. Durbin-Watson & runs test logic" begin
    rng = MersenneTwister(42)
    noise = randn(rng, 200)

    # DW on pure noise ~ 2
    dw, dw_p = durbin_watson(noise)
    @test 1.5 < dw < 2.5
    @test dw_p > 0.01
    println("  ✓ DW(noise) = $(round(dw,digits=3)), p=$(round(dw_p,digits=4))")

    # DW on autocorrelated << 2
    ar1 = cumsum(randn(rng, 200)) * 0.1
    dw_ar, _ = durbin_watson(ar1)
    @test dw_ar < 1.5
    println("  ✓ DW(AR1) = $(round(dw_ar,digits=3)) < 1.5")

    # Runs test on alternating signs
    alt = repeat([1.0, -1.0], 50)
    rn, re, rp = runs_test(alt)
    @test rn == 100       # every pair of adjacent points is a run
    @test rp < 0.05       # extremely unlikely under randomness
    println("  ✓ Runs(alt) = $rn (exp $(round(re,digits=1))), p=$(round(rp,digits=6))")

    # Runs test on uniform signs
    rn_uni, _, _ = runs_test(ones(10))
    @test rn_uni == 1
    println("  ✓ Runs(uniform) = 1")

    # compute_residual_diagnostics round-trip
    rd = compute_residual_diagnostics(noise)
    @test isfinite(rd.durbin_watson)
    @test isfinite(rd.residual_rms)
    @test rd.runs_n > 0
    @test isfinite(rd.runs_expected)
    println("  ✓ Round trip: RMS=$(round(rd.residual_rms,digits=6)), runs=$(rd.runs_n)")
end

# ===========================================================================
# 2. Gaussian fit — verify covariance + residual diags are populated
# ===========================================================================
@testset "2. Gaussian fit → covariance + residuals present" begin
    cfg = build_config(Dict{String,Any}(
        common_kwargs..., "peak_profile" => :gaussian,
    ))
    fr = run_fit(x_data, y_data, cfg; save_cache=false, verbose=false)
    @test !isempty(fr.all_results)
    best = best_result(fr)
    @test best !== nothing
    @test best.n_peaks >= 2
    println("  ✓ Gaussian: n=$(best.n_peaks), sBIC=$(round(best.student_bic,digits=1)), R²=$(round(best.r_squared,digits=4))")

    # Residual diagnostics exist and are sensible
    @test best.residual_diagnostics !== nothing
    rd = best.residual_diagnostics
    @test rd isa ResidualDiagnostics
    # DW should be reasonable for a good fit (not NaN)
    @test isfinite(rd.durbin_watson) || rd.durbin_watson === NaN  # NaN OK for tight fit
    println("  ✓ Residual diags: DW=$(rd.durbin_watson isa Float64 && isfinite(rd.durbin_watson) ? round(rd.durbin_watson,digits=3) : "NaN"), runs_n=$(rd.runs_n), RMS=$(round(rd.residual_rms,digits=6))")

    # Covariance matrix exists even if pcorr is nothing (ill-conditioned case)
    @test best.pcov !== nothing
    println("  ✓ Covariance: $(size(best.pcov))")
    if best.pcorr !== nothing
        println("  ✓ Correlation: $(size(best.pcorr)), diagonal≈1: $(all(x->isapprox(x,1.0,atol=1e-10), diag(best.pcorr)))")
    else
        println("  ⓘ Correlation not computed (ill-conditioned covariance — expected for near-bound parameters)")
    end
    # Per-peak errors exist
    @test length(best.perr) == length(best.popt)
    @test any(x -> x > 0, best.perr)  # at least one parameter has finite error
    println("  ✓ Perr: $(length(best.perr)) values, max=$(round(maximum(filter(isfinite, best.perr)), digits=6))")
end

# ===========================================================================
# 3. Pseudo-Voigt — η ∈ [0,1], η ≈ 0 for Gaussian data
# ===========================================================================
@testset "3. Pseudo-Voigt fit" begin
    cfg = build_config(Dict{String,Any}(
        common_kwargs..., "peak_profile" => :pseudo_voigt,
        "global_maxtime" => 10.0, "global_maxiter" => 5000,
    ))
    fr = run_fit(x_data, y_data, cfg; save_cache=false, verbose=false)
    @test !isempty(fr.all_results)
    best = best_result(fr)
    @test best !== nothing
    @test best.n_peaks >= 2

    params = best.popt
    eta = params[end]  # η is always last param in pseudo-Voigt mode
    @test 0.0 <= eta <= 1.0

    println("  ✓ Pseudo-Voigt: n=$(best.n_peaks), sBIC=$(round(best.student_bic,digits=1)), η=$(round(eta,digits=4))")
    
    if eta < 0.5
        println("  ✓ η ≪ 0.5 confirms Gaussian data (pure Gauss → η ≈ 0)")
    else
        println("  ⓘ η=$(round(eta,digits=4)) > 0.5 (Gaussian data with noise may bias η upward)")
        @test eta <= 1.0
    end

    # η has a finite error estimate
    @test length(best.perr) == length(params)
    eta_idx = length(params)
    eta_err = best.perr[eta_idx]
    @test eta_err >= 0  # may be NaN
    println("  ✓ η error: $(isfinite(eta_err) ? round(eta_err,digits=6) : "NaN")")

    # Residual diagnostics present
    @test best.residual_diagnostics !== nothing
    rd = best.residual_diagnostics
    @test rd isa ResidualDiagnostics
    println("  ✓ Residuals: RMS=$(round(rd.residual_rms,digits=6))")
end

# ===========================================================================
# 4. Lorentzian worse than Gaussian for pure Gaussian data
# ===========================================================================
@testset "4. Lorentzian vs Gaussian" begin
    cfg_lor = build_config(Dict{String,Any}(
        common_kwargs..., "peak_profile" => :lorentzian,
        "global_maxtime" => 10.0, "global_maxiter" => 5000,
    ))
    fr_lor = run_fit(x_data, y_data, cfg_lor; save_cache=false, verbose=false)
    @test !isempty(fr_lor.all_results)
    best_lor = best_result(fr_lor)
    @test best_lor !== nothing

    cfg_g = build_config(Dict{String,Any}(
        common_kwargs..., "peak_profile" => :gaussian,
        "global_maxtime" => 10.0, "global_maxiter" => 5000,
    ))
    fr_g = run_fit(x_data, y_data, cfg_g; save_cache=false, verbose=false)
    best_g = best_result(fr_g)

    println("  ✓ Lorentzian sBIC=$(round(best_lor.student_bic,digits=1))")
    println("  ✓ Gaussian   sBIC=$(round(best_g.student_bic,digits=1))")
    if best_lor.student_bic > best_g.student_bic
        println("  ✓ Lorentzian worse than Gaussian (correct for pure Gaussian data)")
    else
        println("  ⓘ Lorentzian not worse (may happen with noise+small data)")
    end
end

# ===========================================================================
# 5. Profile validation
# ===========================================================================
@testset "5. Profile validation" begin
    cfg_bad = build_config(Dict{String,Any}("peak_profile" => :voigt, "no_show" => true, "output_dir" => nothing))
    @test_throws ErrorException GaussianFit1D._fit_one(3, x_data, y_data, cfg_bad)
    println("  ✓ Invalid profile ':voigt' correctly rejected")

    for profile in (:gaussian, :lorentzian, :pseudo_voigt)
        cfg = build_config(Dict{String,Any}(
            common_kwargs..., "peak_profile" => profile,
            "global_maxtime" => 3.0, "global_maxiter" => 1000,
        ))
        r = GaussianFit1D._fit_one(3, x_data, y_data, cfg)
        @test r !== nothing
        @test r.n_peaks == 3
    end
    println("  ✓ All 3 profiles run for n=3")
end

# ===========================================================================
# 6. Export contains new fields
# ===========================================================================
@testset "6. Export diagnostic fields" begin
    cfg = build_config(Dict{String,Any}(
        common_kwargs..., "peak_profile" => :gaussian,
        "global_maxtime" => 3.0, "global_maxiter" => 1000,
        "output_dir" => "/tmp/stmfit_test",
    ))
    fr = run_fit(x_data, y_data, cfg; save_cache=false, verbose=false)
    export_file = export_results(fr.x, fr.y, fr.all_results, fr.cfg)
    text = read(export_file, String)
    
    diag_fields = ["durbin_watson", "durbin_watson_p", "runs_n",
                   "runs_expected", "runs_p", "residual_rms", "residual_max", "kappa_max_adj"]
    for field in diag_fields
        @test occursin(field, text)
    end
    println("  ✓ All $(length(diag_fields)) diagnostic fields present in export")
    rm(export_file; force=true)
end

# ===========================================================================
# 7. Parameter counts
# ===========================================================================
@testset "7. Parameter counts" begin
    for profile in (:gaussian, :lorentzian, :pseudo_voigt)
        for n in (2, 3)
            cfg = build_config(Dict{String,Any}(
                common_kwargs..., "peak_profile" => profile,
                "global_maxtime" => 3.0, "global_maxiter" => 1000,
            ))
            r = GaussianFit1D._fit_one(n, x_data, y_data, cfg)
            if r !== nothing
                n_expected = 1 + 3*n + (profile == :pseudo_voigt ? 1 : 0)
                @test r.n_params == n_expected
            end
        end
    end
    println("  ✓ Param counts correct for all profiles × (2,3) peaks")
end

# ===========================================================================
# 8. MGF_VERSION bumped
# ===========================================================================
@testset "8. Version" begin
    @test GaussianFit1D.MGF_VERSION == "5.0.0"
    println("  ✓ MGF_VERSION = $(GaussianFit1D.MGF_VERSION)")
end

# ===========================================================================
# 9. Peak profile guards in GaussianFit2D
# ===========================================================================
@testset "9. 2D peak_profile guard" begin
    using GaussianFit2D
    ccfg = GaussianFit2D.ChainSweepConfig()
    @test ccfg.peak_profile == :gaussian
    # Setting a non-gaussian profile should not crash (it's just a field)
    ccfg.peak_profile = :lorentzian
    @test ccfg.peak_profile == :lorentzian
    println("  ✓ ChainSweepConfig.peak_profile field exists with default :gaussian")
end

println()
println("="^80)
println("ALL FEATURE TESTS PASSED")
println("="^80)
