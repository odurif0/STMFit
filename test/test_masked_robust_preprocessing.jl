#!/usr/bin/env julia
# Synthetic arithmetic only. No real scans, fits, jobs, grades, or production imports.
# julia --startup-file=no --threads=1 --project=. test/test_masked_robust_preprocessing.jl
using Test
using LinearAlgebra
using Statistics
using TOML

include(joinpath(@__DIR__, "lib", "masked_robust_preprocessing.jl"))
using .MaskedRobustPreprocessing

# Every parameter is explicit, including small support counts for tiny fixtures.
function option_dict(; changes...)
    section = Dict{String,Any}(
        "huber_delta" => 1.345, "maxiter" => 50,
        "coefficient_rtol" => 1e-8, "weight_atol" => 1e-8,
        "scale_floor_nm" => 1e-9, "rank_rtol" => 1e-12,
        "min_background_pixels" => 3, "min_background_fraction" => 0.0,
        "min_row_background_pixels" => 2)
    for (key, value) in changes
        section[String(key)] = value
    end
    return Dict("masked_robust_preprocessing" => section)
end
options(; changes...) = read_options(option_dict(; changes...))
affine(xs, ys; a=0.7, bx=0.2, by=-0.3) = [a + bx*x + by*y for y in ys, x in xs]
coeffvec(result) = collect(values(result.coefficients))

@testset "masked robust preprocessing (synthetic only)" begin
    @test VERSION >= v"1.13.0" && VERSION < v"1.14.0"
    @test Threads.nthreads() == 1

    @testset "explicit options and input contracts" begin
        cfg = option_dict()
        @test options().huber_delta == 1.345
        @test options().maxiter == 50
        io = IOBuffer(); TOML.print(io, cfg)
        @test read_options(TOML.parse(String(take!(io)))).rank_rtol == 1e-12
        @test_throws ArgumentError read_options(Dict())
        @test_throws ArgumentError read_options(Dict("masked_robust_preprocessing" => 3))
        for key in keys(cfg["masked_robust_preprocessing"])
            incomplete = deepcopy(cfg)
            delete!(incomplete["masked_robust_preprocessing"], key)
            @test_throws ArgumentError read_options(incomplete)
        end
        extra = deepcopy(cfg); extra["masked_robust_preprocessing"]["impute"] = false
        @test_throws ArgumentError read_options(extra)
        for (key, invalid) in (
            (:huber_delta, 0), (:huber_delta, -1), (:huber_delta, NaN),
            (:huber_delta, Inf), (:huber_delta, true), (:huber_delta, "1.345"),
            (:maxiter, 0), (:maxiter, 2.0), (:maxiter, false),
            (:coefficient_rtol, -1), (:coefficient_rtol, Inf), (:coefficient_rtol, 1),
            (:weight_atol, -1), (:weight_atol, NaN), (:weight_atol, 1),
            (:scale_floor_nm, 0), (:scale_floor_nm, Inf), (:scale_floor_nm, false),
            (:rank_rtol, 0), (:rank_rtol, 1), (:rank_rtol, Inf),
            (:min_background_pixels, 0), (:min_background_pixels, 2.5),
            (:min_background_fraction, -0.1), (:min_background_fraction, 1.1),
            (:min_background_fraction, NaN), (:min_background_fraction, true),
            (:min_row_background_pixels, 0), (:min_row_background_pixels, 2.0))
            @test_throws ArgumentError options(; (key => invalid,)...)
        end
        xs, ys = collect(0.0:4.0), collect(0.0:2.0)
        z = affine(xs, ys); mask = trues(size(z)); opts = options()
        @test_throws ArgumentError masked_background(xs, ys, z, mask, opts; loss=:soft)
        @test_throws ArgumentError masked_background(reverse(xs), ys, z, mask, opts)
        @test_throws ArgumentError masked_background(xs, reverse(ys), z, mask, opts)
        @test_throws ArgumentError masked_background([0,1,1,3,4], ys, z, mask, opts)
        @test_throws ArgumentError masked_background([0,1,NaN,3,4], ys, z, mask, opts)
        @test_throws ArgumentError masked_background(xs, [0,1,Inf], z, mask, opts)
        @test_throws ArgumentError masked_background(Float64[], ys, z, mask, opts)
        @test_throws ArgumentError masked_background(xs, ys, z, ones(Int, size(z)), opts)
        @test_throws DimensionMismatch masked_background(xs, ys, transpose(z), mask, opts)
        @test_throws DimensionMismatch masked_background(xs, ys, z, trues(2,2), opts)
    end

    @testset "affine recovery, supplied footprint, finite support, immutability" begin
        xs = collect(range(-2, 3; length=13)); ys = collect(range(-1, 2; length=7))
        truth = affine(xs, ys)
        mask = trues(size(truth)); mask[2:6, 6:8] .= false
        signal = zeros(size(truth)); signal[2:6, 6:8] .= 0.24
        z = truth + signal
        z[1,1] = NaN; z[3,2] = Inf; z[6,12] = -Inf
        saved = deepcopy((xs, ys, z, mask))
        for loss in (:ols, :huber)
            result = masked_background(xs, ys, z, mask, options(); loss)
            @test result.status == :OK
            @test result.reason == :ok
            @test result.converged
            @test result.design_rank == 3
            @test coeffvec(result) ≈ [0.7,0.2,-0.3] atol=2e-14
            @test result.plane ≈ truth atol=2e-14
            @test maximum(abs.(result.row_offsets)) < 2e-14
            @test result.observed == isfinite.(z)
            @test result.background_mask == (mask .& isfinite.(z))
            @test result.valid == isfinite.(z)
            @test isnan.(result.corrected) == .!isfinite.(z)
            @test isnan.(result.background) == .!isfinite.(z)
            @test all(isfinite, result.plane) # predictions, not observations
            @test result.corrected[result.valid] ≈ signal[result.valid] atol=2e-14
            @test result.background[result.valid] ≈ truth[result.valid] atol=2e-14
            @test all(==(:ok), result.row_status)
            @test result.row_background_counts == vec(sum(mask .& isfinite.(z); dims=2))
            for j in eachindex(ys)
                @test abs(median(result.corrected[j, result.background_mask[j,:]])) < 2e-14
            end
            @test isequal((xs, ys, z, mask), saved)
        end

        clean = copy(truth); extreme = copy(truth)
        extreme[.!mask] .= [isodd(i) ? 1e300 : -1e300 for i in 1:count(.!mask)]
        for loss in (:ols, :huber)
            base = masked_background(xs, ys, clean, mask, options(); loss)
            changed = masked_background(xs, ys, extreme, mask, options(); loss)
            @test changed.status == :OK
            @test changed.coefficients == base.coefficients
            @test changed.robust_scale == base.robust_scale
            @test changed.plane == base.plane
            @test changed.row_offsets == base.row_offsets
            @test changed.background == base.background
            @test changed.background_mask == mask
        end
    end

    @testset "row medians and nonidentifiable decomposition" begin
        xs = collect(-3.0:3.0); ys = collect(-3.0:3.0)
        truth = affine(xs, ys)
        shifts = [0.06,-0.04,-0.02,0.0,-0.02,-0.04,0.06]
        mask = trues(size(truth)); mask[:,4] .= false # balanced x support
        signal = zeros(size(truth)); signal[:,4] .= 0.2
        z = truth .+ shifts .+ signal
        ols = masked_background(xs, ys, z, mask, options(); loss=:ols)
        @test ols.coefficients.intercept_nm ≈ 0.7 atol=2e-14
        @test ols.coefficients.slope_x ≈ 0.2 atol=2e-14
        @test ols.coefficients.slope_y ≈ -0.3 atol=2e-14
        @test ols.row_offsets ≈ shifts atol=2e-14 # zero mean and zero y projection
        for loss in (:ols, :huber)
            result = masked_background(xs, ys, z, mask, options(); loss)
            @test result.status == :OK
            @test result.background ≈ truth .+ shifts atol=3e-14
            @test result.corrected ≈ signal atol=3e-14
        end
        # A linear row shift is absorbed by the first-stage y tilt. The total
        # background is identified here; its generating decomposition is not.
        row_linear = truth .+ (0.4 .* ys)
        result = masked_background(xs, ys, row_linear, trues(size(truth)), options(); loss=:ols)
        @test result.coefficients.slope_y ≈ 0.1 atol=2e-14
        @test maximum(abs.(result.row_offsets)) < 2e-14
        @test result.background ≈ row_linear atol=2e-14
        @test maximum(abs.(result.corrected)) < 2e-14

        # Asymmetric support correlates row shifts with x: this sequential
        # method cannot promise recovery of the generating x slope either.
        asymmetric = falses(size(truth))
        for j in eachindex(ys)
            asymmetric[j, (j <= 3 ? (1:3) : (5:7))] .= true
        end
        asymmetric_z = truth .+ [0.2,0.2,0.2,-0.15,-0.15,-0.15,-0.15]
        biased = masked_background(xs, ys, asymmetric_z, asymmetric, options(); loss=:ols)
        @test biased.status == :OK
        @test abs(biased.coefficients.slope_x - 0.2) > 0.01
        for j in eachindex(ys)
            @test abs(median(biased.corrected[j, asymmetric[j,:]])) < 2e-14
        end
    end

    @testset "zero, constant, missing, rank and global support failures" begin
        xs = collect(-2.0:2.0); ys = collect(-2.0:2.0); mask = trues(5,5)
        for value in (0.0, -7.25), loss in (:ols, :huber)
            z = fill(value, 5,5)
            result = masked_background(xs, ys, z, mask, options(); loss)
            @test result.status == :OK
            @test result.coefficients == (intercept_nm=value, slope_x=0.0, slope_y=0.0)
            @test all(==(0.0), result.corrected)
            @test result.background == z
            @test result.robust_scale == 1e-9
            @test result.scale_floor_active
            @test result.converged
        end
        for value in (NaN, Inf, -Inf), loss in (:ols, :huber)
            z = fill(value, 5,5); before = copy(z)
            result = masked_background(xs, ys, z, mask, options(); loss)
            @test result.status == :UNAVAILABLE
            @test result.reason == :insufficient_background_pixels
            @test !any(result.observed)
            @test !any(result.background_mask)
            @test !any(result.valid)
            @test all(isnan, result.corrected)
            @test all(isnan, result.background)
            @test all(isnan, result.plane)
            @test !result.converged
            @test result.iterations == 0
            @test isequal(z, before)
        end
        z = affine(xs, ys)
        for support in (:single_row, :single_column, :diagonal)
            selected = falses(5,5)
            support == :single_row && (selected[3,:] .= true)
            support == :single_column && (selected[:,3] .= true)
            support == :diagonal && (selected[diagind(selected)] .= true)
            for loss in (:ols, :huber)
                result = masked_background(xs, ys, z, selected, options(); loss)
                @test result.status == :UNAVAILABLE
                @test result.reason == :rank_deficient
                @test result.design_rank == 2
                @test all(isnan, result.corrected)
                @test all(isnan, result.plane)
                @test !any(result.valid)
                @test result.background_mask == selected
            end
        end
        too_few = masked_background(xs, ys, z, mask, options(min_background_pixels=26))
        @test too_few.reason == :insufficient_background_pixels
        selected = copy(mask); selected[:,1:3] .= false
        fraction = masked_background(xs, ys, z, selected, options(min_background_fraction=0.5))
        @test fraction.reason == :insufficient_background_fraction
        missing = copy(z); missing[:,1:3] .= NaN
        fraction_missing = masked_background(xs, ys, missing, mask, options(min_background_fraction=0.5))
        @test fraction_missing.reason == :insufficient_background_fraction # denominator = whole grid
        empty = masked_background(xs, ys, z, falses(5,5), options())
        @test empty.reason == :insufficient_background_pixels
        severe_rank = masked_background(xs, ys, z, mask, options(rank_rtol=0.99))
        @test severe_rank.reason == :rank_deficient
        @test severe_rank.design_rank == 1
        single_x = masked_background([0.0], ys, z[:,1:1], trues(5,1), options(min_row_background_pixels=1))
        @test single_x.reason == :rank_deficient
    end

    @testset "row support stays unavailable without extrapolation" begin
        xs = collect(0.0:7.0); ys = collect(0.0:5.0)
        z = affine(xs, ys); mask = trues(size(z))
        mask[2,:] .= false # finite foreground, no background
        mask[3,:] .= false; mask[3,1] = true # below minimum of two
        mask[4,:] .= false; mask[4,1:2] .= true # exactly the row threshold
        z[5,:] .= NaN # no raw observations
        saved = deepcopy((xs, ys, z, mask))
        for loss in (:ols, :huber)
            result = masked_background(xs, ys, z, mask, options(); loss)
            @test result.status == :PARTIAL
            @test result.reason == :insufficient_row_background
            @test result.row_background_counts == [8,0,1,2,0,8]
            @test result.row_status == [:ok,:insufficient_background,:insufficient_background,
                                        :ok,:insufficient_background,:ok]
            @test result.converged # plane convergence does not imply row availability
            @test all(isfinite, result.plane)
            @test all(isnan, result.corrected[[2,3,5],:])
            @test all(isnan, result.background[[2,3,5],:])
            @test all(isnan, result.row_offsets[[2,3,5]])
            @test !any(result.valid[[2,3,5],:])
            @test all(result.valid[[1,4,6],:])
            @test result.background_mask[3,1] # still used by global plane
            @test isequal((xs, ys, z, mask), saved)
        end
        none = masked_background(xs, ys, z, mask, options(min_row_background_pixels=9); loss=:ols)
        @test none.status == :UNAVAILABLE
        @test none.reason == :no_supported_rows
        @test none.converged
        @test all(isnan, none.corrected)
        @test !any(none.valid)
    end

    @testset "controlled outliers, fixed scale, weighted solve and stopping" begin
        xs = collect(range(-2, 3; length=31)); ys = collect(range(-1, 2; length=21))
        truth = affine(xs, ys)
        z = truth + [0.006*(sin(3x+2y) + 0.6cos(4x-3y)) for y in ys, x in xs]
        z[2:4:20, 23:2:31] .+= 2.0
        mask = trues(size(z)); opts = options()
        before = deepcopy((xs, ys, z, mask))
        ols = masked_background(xs, ys, z, mask, opts; loss=:ols)
        huber = masked_background(xs, ys, z, mask, opts; loss=:huber)
        @test ols.status == :OK
        @test huber.status == :OK
        @test huber.converged
        @test huber.iterations > 1
        @test norm(huber.plane - truth) < 0.5 * norm(ols.plane - truth)
        @test norm(coeffvec(huber) - [0.7,0.2,-0.3]) < 0.5 * norm(coeffvec(ols) - [0.7,0.2,-0.3])
        initial_residuals = vec(z - ols.plane)
        expected_scale = max(median(abs.(initial_residuals .- median(initial_residuals))) /
                             0.6744897501960817, opts.scale_floor_nm)
        @test huber.robust_scale ≈ expected_scale rtol=1e-13
        @test huber.robust_scale == ols.robust_scale
        @test ols.objective_nm2 ≈ sum(abs2, initial_residuals) / 2 rtol=1e-13
        @test ols.initial_objective_nm2 == ols.objective_nm2
        cutoff = opts.huber_delta * huber.robust_scale
        rho(r) = abs(r) <= cutoff ? r^2 / 2 : cutoff * (abs(r) - cutoff / 2)
        @test huber.initial_objective_nm2 ≈ sum(rho, initial_residuals) rtol=1e-13
        @test huber.objective_nm2 ≈ sum(rho, z - huber.plane) rtol=1e-13
        @test huber.objective_nm2 < huber.initial_objective_nm2
        xx_centered = [2*(x - (first(xs)+last(xs))/2)/(last(xs)-first(xs)) for y in ys, x in xs]
        yy_centered = [2*(y - (first(ys)+last(ys))/2)/(last(ys)-first(ys)) for y in ys, x in xs]
        centered_design = hcat(ones(length(z)), vec(xx_centered), vec(yy_centered))
        expected_gradient = maximum(abs.(centered_design' * vec(clamp.(z - huber.plane, -cutoff, cutoff))))
        @test huber.stationarity_inf_nm ≈ expected_gradient atol=1e-12 # sum of 651 roundoff terms
        @test ols.stationarity_inf_nm < 1e-12
        @test huber.coefficient_change <= opts.coefficient_rtol
        @test huber.weight_change <= opts.weight_atol
        @test isequal((xs, ys, z, mask), before)

        limited = masked_background(xs, ys, z, mask, options(maxiter=1); loss=:huber)
        @test limited.status == :NONCONVERGED
        @test limited.reason == :maxiter
        @test !limited.converged
        @test limited.iterations == 1
        @test limited.coefficient_change > opts.coefficient_rtol || limited.weight_change > opts.weight_atol
        @test all(isfinite, limited.plane) # retained final estimate, not a success
        @test all(isfinite, coeffvec(limited))
        @test all(isnan, limited.corrected)
        @test all(isnan, limited.background)
        @test all(isnan, limited.row_offsets)
        @test !any(limited.valid)
        @test limited.observed == trues(size(z))
        @test limited.background_mask == mask
        @test all(==(:plane_unavailable), limited.row_status)
        @test limited.objective_nm2 ≈ sum(rho, z - limited.plane) rtol=1e-13
        @test limited.objective_nm2 <= limited.initial_objective_nm2
        @test isfinite(limited.stationarity_inf_nm)
        @test limited.stationarity_inf_nm > huber.stationarity_inf_nm

        # Independent uncentered physical-coordinate QR for the first IRLS step.
        xx = [x for y in ys, x in xs]; yy = [y for y in ys, x in xs]
        design = hcat(ones(length(z)), vec(xx), vec(yy))
        w = min.(1.0, (opts.huber_delta * expected_scale) ./ abs.(initial_residuals))
        qr_beta = (design .* sqrt.(w)) \ (vec(z) .* sqrt.(w))
        @test coeffvec(limited) ≈ qr_beta rtol=1e-12 atol=1e-13
        @test vec(limited.plane) ≈ design * qr_beta rtol=1e-12 atol=1e-13

        coefficient_only = masked_background(xs, ys, z, mask,
            options(maxiter=1, coefficient_rtol=0.9, weight_atol=0.0); loss=:huber)
        @test coefficient_only.coefficient_change <= 0.9
        @test coefficient_only.weight_change > 0.0
        @test coefficient_only.status == :NONCONVERGED
        weight_only = masked_background(xs, ys, z, mask,
            options(maxiter=1, coefficient_rtol=0.0, weight_atol=0.99); loss=:huber)
        @test weight_only.coefficient_change > 0.0
        @test weight_only.weight_change <= 0.99
        @test weight_only.status == :NONCONVERGED
        exact = masked_background(xs, ys, truth, mask, options(maxiter=1); loss=:huber)
        @test exact.status == :OK # a limit of one is not automatically failure
        @test exact.iterations == 1

        # z scaling changes the numerical floor in the same units, not the model.
        for factor in (0.01, 1e3, -2.0), loss in (:ols, :huber)
            base = loss == :ols ? ols : huber
            scaled = masked_background(xs, ys, factor .* z, mask,
                options(scale_floor_nm=abs(factor)*opts.scale_floor_nm); loss)
            @test scaled.status == base.status
            @test scaled.plane ≈ factor .* base.plane rtol=2e-12 atol=1e-12
            @test scaled.corrected ≈ factor .* base.corrected rtol=2e-10 atol=1e-12
            @test coeffvec(scaled) ≈ factor .* coeffvec(base) rtol=2e-12 atol=1e-12
            @test scaled.robust_scale ≈ abs(factor)*base.robust_scale rtol=2e-12
        end
        for loss in (:ols, :huber)
            base = loss == :ols ? ols : huber
            shifted = masked_background(xs, ys, z .+ 13.7, mask, opts; loss)
            @test shifted.status == base.status
            @test shifted.corrected ≈ base.corrected atol=1e-13
            @test shifted.coefficients.intercept_nm ≈ base.coefficients.intercept_nm + 13.7 atol=1e-13
            @test shifted.robust_scale ≈ base.robust_scale rtol=1e-12
            tx, ty, sx, sy = 1e4, -5e3, 0.01, 300.0
            transformed = masked_background(tx .+ sx .* xs, ty .+ sy .* ys, z, mask, opts; loss)
            @test transformed.status == base.status
            # This tolerance is per pixel, not the Frobenius norm of 651 errors.
            @test maximum(abs.(transformed.corrected - base.corrected)) <= 2e-10
            @test maximum(abs.(transformed.plane - base.plane)) <= 2e-10
            @test transformed.coefficients.slope_x ≈ base.coefficients.slope_x/sx rtol=2e-8
            @test transformed.coefficients.slope_y ≈ base.coefficients.slope_y/sy rtol=2e-8
            @test transformed.coefficients.intercept_nm ≈ base.coefficients.intercept_nm -
                base.coefficients.slope_x*tx/sx - base.coefficients.slope_y*ty/sy rtol=2e-8
        end
    end

    @testset "large coordinate origins retain centered evaluation" begin
        xs = 1e12 .+ collect(0.0:0.125:1.5)
        ys = -1e10 .+ collect(0.0:0.25:1.5)
        z = [0.7 + 0.2*(x-1e12) - 0.3*(y+1e10) for y in ys, x in xs]
        for loss in (:ols, :huber)
            result = masked_background(xs, ys, z, trues(size(z)), options(); loss)
            @test result.status == :OK
            @test maximum(abs.(result.plane - z)) < 2e-14
            @test maximum(abs.(result.corrected)) < 2e-14
            @test result.coefficients.slope_x ≈ 0.2 atol=2e-14
            @test result.coefficients.slope_y ≈ -0.3 atol=2e-14
        end
    end

    @testset "native y,x order and explicit transpose" begin
        xs = collect(range(-4, 7; length=11)); ys = collect(range(2, 5; length=6))
        z = affine(xs, ys; a=-1.2, bx=0.43, by=-0.17)
        mask = trues(size(z)); mask[2,3] = false; z[5,9] = NaN
        for loss in (:ols, :huber)
            native = masked_background(xs, ys, z, mask, options(); loss)
            @test size(native.corrected) == (6,11)
            @test coeffvec(native) ≈ [-1.2,0.43,-0.17] atol=1e-14
            transposed = masked_background(ys, xs, permutedims(z), permutedims(mask), options(); loss)
            @test transposed.status == :OK
            @test coeffvec(transposed) ≈ [-1.2,-0.17,0.43] atol=1e-14
            @test transposed.observed == permutedims(native.observed)
            @test transposed.background_mask == permutedims(native.background_mask)
            @test isnan.(transposed.corrected) == permutedims(isnan.(native.corrected))
            @test transposed.plane ≈ permutedims(native.plane) atol=1e-14
            @test_throws DimensionMismatch masked_background(xs, ys, permutedims(z), mask, options(); loss)
        end
    end
end
