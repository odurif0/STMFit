#!/usr/bin/env julia

using Test
using Random
using SHA
using Printf
using Statistics
using TOML
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "lib", "structured_assignment", "robust_emissions.jl"))
using .StructuredRobustEmissions
const SRE = StructuredRobustEmissions

include(joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl"))
const HUA = HierarchicalUnitAssignment

function independent_logsumexp(values)
    m = maximum(values)
    return m + log(sum(exp(v - m) for v in values))
end

function observed_loglik(X, means, scales)
    total = 0.0
    for i in axes(X, 1)
        terms = [log(0.5) + log_student_t_diag(view(X, i, :),
                                               view(means, k, :),
                                               view(scales, k, :)) for k in 1:2]
        total += independent_logsumexp(terms)
    end
    return total
end

function canonical_result_hash(fit::StudentTTwoFit)
    io = IOBuffer()
    println(io, Int(fit.status), '|', Int(fit.invalid_reason), '|', fit.best_start,
            '|', fit.high_amplitude_component)
    for value in fit.means
        println(io, @sprintf("%.17g", value))
    end
    for value in fit.scales
        println(io, @sprintf("%.17g", value))
    end
    println(io, @sprintf("%.17g", fit.loglik))
    for value in fit.responsibilities
        println(io, @sprintf("%.17g", value))
    end
    return bytes2hex(sha256(take!(io)))
end

function canonical_trace_hash(fit::StudentTTwoFit)
    io = IOBuffer()
    for trace in fit.traces
        println(io, trace.start_index, '|', @sprintf("%.17g", trace.quantile_delta),
                '|', trace.converged, '|', trace.monotone, '|', trace.valid,
                '|', Int(trace.invalid_reason), '|', trace.rejected_iteration)
        println(io, @sprintf("%.17g", trace.rejected_loglik))
        for value in trace.accepted_loglik
            println(io, @sprintf("%.17g", value))
        end
        for value in trace.final_means
            println(io, @sprintf("%.17g", value))
        end
        for value in trace.final_scales
            println(io, @sprintf("%.17g", value))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function two_cluster_fixture(; contaminated=false)
    left = [-2.0 + 0.18 * sin(0.71 * i) for i in 1:40]
    right = [2.0 + 0.16 * cos(0.63 * i) for i in 1:40]
    values = vcat(left, right)
    amplitudes = vcat([1.0 + 0.03 * sin(i) for i in 1:40],
                      [3.0 + 0.03 * cos(i) for i in 1:40])
    if contaminated
        append!(values, [-30.0, -25.0, 25.0, 30.0])
        append!(amplitudes, [1.0, 1.0, 3.0, 3.0])
    end
    return reshape(values, :, 1), amplitudes
end

@testset "fixed-nu analytic density contract" begin
    @test FIXED_NU === 8
    @test FIXED_NU isa Int
    @test EQUAL_MIXTURE_WEIGHTS === (0.5, 0.5)
    @test SCALE_FLOOR == 1.0e-4
    @test MAX_ITER === 200
    @test CONVERGENCE_TOL == 1.0e-6
    @test DECREASE_TOL == 1.0e-10
    @test START_QUANTILE_DELTAS === (0.25, 0.40, 0.45, 0.475, 0.49)

    model = TOML.parsefile(joinpath(ROOT, "config", "unit_assignment_structured_model.toml"))["model"]
    @test model["student_t_nu"] == FIXED_NU
    @test Tuple(model["mixture_weights"]) == EQUAL_MIXTURE_WEIGHTS
    @test model["mixture_weight_policy"] == "fixed_not_estimated"
    @test model["covariance_floor"] == SCALE_FLOOR
    @test model["max_iter"] == MAX_ITER
    @test model["tol"] == CONVERGENCE_TOL
    @test model["objective_decrease_tolerance"] == DECREASE_TOL
    @test Tuple(model["unary_start_quantile_deltas"]) == START_QUANTILE_DELTAS
    @test model["physical_orientation_feature"] == "amplitude"

    @test SRE._loggamma_integer_or_half_twice(8) ≈ log(6.0) atol=2e-15 rtol=0
    @test SRE._loggamma_integer_or_half_twice(9) ≈ 2.4537365708424423 atol=2e-15 rtol=0
    @test SRE._loggamma_integer_or_half_twice(10) ≈ log(24.0) atol=2e-15 rtol=0
    @test SRE._loggamma_integer_or_half_twice(11) ≈ 3.9578139676187165 atol=3e-15 rtol=0
    @test_throws ArgumentError SRE._loggamma_integer_or_half_twice(0)

    fixtures = [
        ([1.25], [-0.5], [2.0], -2.0847866970055247),
        ([1.0, -2.0], [0.25, -1.5], [0.5, 3.0], -2.7439462679132856),
        ([0.5, -1.25, 2.0], [0.1, -0.75, 1.0], [0.8, 1.5, 2.5], -3.7228418858340011),
    ]
    for (x, mu, scale, expected) in fixtures
        @test log_student_t_diag(x, mu, scale) ≈ expected atol=2e-14 rtol=2e-14
    end
    @test_throws DimensionMismatch log_student_t_diag([1.0], [0.0, 1.0], [1.0])
    @test_throws ArgumentError log_student_t_diag([1.0], [0.0], [0.0])
    @test_throws ArgumentError log_student_t_diag([NaN], [0.0], [1.0])
end

@testset "fixed weights, exact ECM denominators, and projection" begin
    X = [0.0 1.0; 2.0 3.0; 4.0 5.0]
    means = [-0.5 0.5; 3.5 4.5]
    scales = [1.2 0.8; 0.7 1.6]
    resp = [0.8 0.2; 0.4 0.6; 0.1 0.9]
    proposal = SRE._ecm_proposal(X, means, scales, resp)
    @test proposal.reason == ROBUST_OK

    expected_means = zeros(2, 2)
    expected_scales = zeros(2, 2)
    for k in 1:2
        u = Float64[]
        for i in axes(X, 1)
            delta = sum((X[i, j] - means[k, j])^2 / scales[k, j] for j in axes(X, 2))
            push!(u, (FIXED_NU + size(X, 2)) / (FIXED_NU + delta))
        end
        ru = resp[:, k] .* u
        expected_means[k, :] .= vec(sum(ru .* X; dims=1)) ./ sum(ru)
        for j in axes(X, 2)
            numerator = sum(resp[i, k] * u[i] * (X[i, j] - expected_means[k, j])^2
                            for i in axes(X, 1))
            expected_scales[k, j] = max(numerator / sum(resp[:, k]), SCALE_FLOOR)
        end
    end
    @test proposal.means ≈ expected_means atol=2e-15 rtol=2e-15
    @test proposal.scales ≈ expected_scales atol=2e-15 rtol=2e-15

    one_row = reshape([0.7, -0.4], 1, :)
    two_means = [0.0 0.0; 1.0 -1.0]
    two_scales = [1.0 2.0; 0.5 3.0]
    r = student_t_responsibilities(one_row, two_means, two_scales)
    raw = [log(0.5) + log_student_t_diag(vec(one_row), vec(two_means[k, :]),
                                         vec(two_scales[k, :])) for k in 1:2]
    expected_r1 = exp(raw[1] - independent_logsumexp(raw))
    @test r[1, 1] ≈ expected_r1 atol=2e-15 rtol=2e-15
    @test r[1, 2] ≈ 1.0 - expected_r1 atol=2e-15 rtol=2e-15
    @test sum(r; dims=2) ≈ ones(1, 1) atol=2e-15 rtol=0

    old_means = [1.0 2.0; 3.0 4.0]
    old_scales = fill(1.0, 2, 2)
    new_means = fill(99.0, 2, 2)
    new_scales = fill(SCALE_FLOOR, 2, 2)
    rejected = SRE._accept_projected_state(old_means, old_scales, 10.0,
                                            new_means, new_scales, 9.9)
    @test !rejected.accepted
    @test rejected.means == old_means
    @test rejected.scales == old_scales
    @test rejected.loglik == 10.0
    accepted = SRE._accept_projected_state(old_means, old_scales, 10.0,
                                            new_means, new_scales, 10.0 - 0.5e-10)
    @test accepted.accepted
    @test accepted.means == new_means
    @test accepted.scales == new_scales
    @test accepted.loglik == 10.0 - 0.5e-10

    collapsed_resp = [1.0 0.0; 1.0 0.0; 1.0 0.0]
    collapsed = SRE._ecm_proposal(X, means, scales, collapsed_resp)
    @test collapsed.reason == COLLAPSED_COMPONENT
    @test SRE._best_valid_start_index([1.0, 1.0, 0.5], [true, true, true]) == 1
    @test SRE._best_valid_start_index([1.0, 1.1, 1.1], [true, true, true]) == 2
    @test SRE._best_valid_start_index([1.0, 2.0], [false, false]) == 0
end

@testset "deterministic robust mixture, starts, traces, and amplitude orientation" begin
    clean_X, clean_amp = two_cluster_fixture()
    contaminated_X, contaminated_amp = two_cluster_fixture(contaminated=true)
    clean = fit_student_t_two_component(clean_X, clean_amp)
    contaminated = fit_student_t_two_component(contaminated_X, contaminated_amp)

    @test clean.status == ROBUST_VALID
    @test contaminated.status == ROBUST_VALID
    @test clean.invalid_reason == ROBUST_OK
    @test contaminated.invalid_reason == ROBUST_OK
    @test clean.weights === EQUAL_MIXTURE_WEIGHTS
    @test contaminated.weights === EQUAL_MIXTURE_WEIGHTS
    @test length(clean.traces) == 5
    @test length(contaminated.traces) == 5
    @test [trace.start_index for trace in clean.traces] == collect(1:5)
    @test Tuple(trace.quantile_delta for trace in clean.traces) == START_QUANTILE_DELTAS
    @test all(trace -> all(trace.accepted_loglik[i] >= trace.accepted_loglik[i - 1] - DECREASE_TOL
                           for i in 2:length(trace.accepted_loglik)), clean.traces)
    @test all(trace -> !trace.valid || (trace.converged && trace.monotone), clean.traces)
    @test clean.traces[clean.best_start].valid
    @test clean.loglik == clean.traces[clean.best_start].accepted_loglik[end]
    @test clean.means == clean.traces[clean.best_start].final_means
    @test clean.scales == clean.traces[clean.best_start].final_scales
    @test clean.loglik ≈ observed_loglik(clean_X, clean.means, clean.scales) atol=2e-12 rtol=2e-12
    @test all(clean.scales .>= SCALE_FLOOR)
    @test all(contaminated.scales .>= SCALE_FLOOR)

    replay = fit_student_t_two_component(clean_X, clean_amp)
    @test canonical_result_hash(replay) == canonical_result_hash(clean)
    @test canonical_trace_hash(replay) == canonical_trace_hash(clean)

    reversed_orientation = fit_student_t_two_component(clean_X, 4.0 .- clean_amp)
    @test reversed_orientation.status == ROBUST_VALID
    @test reversed_orientation.means == clean.means
    @test reversed_orientation.scales == clean.scales
    @test reversed_orientation.loglik == clean.loglik
    @test reversed_orientation.responsibilities == clean.responsibilities
    @test reversed_orientation.high_amplitude_component == 3 - clean.high_amplitude_component

    gaussian_clean = HUA.fit_em_two_component(clean_X, 0; n_starts=5,
                                               cov_floor=SCALE_FLOOR,
                                               max_iter=MAX_ITER,
                                               tol=CONVERGENCE_TOL)
    gaussian_contaminated = HUA.fit_em_two_component(contaminated_X, 0; n_starts=5,
                                                      cov_floor=SCALE_FLOOR,
                                                      max_iter=MAX_ITER,
                                                      tol=CONVERGENCE_TOL)
    robust_mean_move = norm(sort(vec(contaminated.means)) - sort(vec(clean.means)))
    gaussian_mean_move = norm(sort(vec(gaussian_contaminated.means)) -
                              sort(vec(gaussian_clean.means)))
    robust_scale_move = norm(sort(vec(contaminated.scales)) - sort(vec(clean.scales)))
    gaussian_scale_move = norm(sort(vec(gaussian_contaminated.vars)) -
                               sort(vec(gaussian_clean.vars)))
    @test robust_mean_move < 0.4 * gaussian_mean_move
    @test robust_scale_move < 0.4 * gaussian_scale_move

    rng = MersenneTwister(7007)
    states = rand(rng, 0:1, 91)
    random_X = Matrix{Float64}(undef, length(states), 2)
    random_amp = Matrix{Float64}(undef, length(states), 1)
    for i in eachindex(states)
        if states[i] == 0
            random_X[i, :] .= (-1.8, -0.9) .+ (0.18 * randn(rng), 0.15 * randn(rng))
            random_amp[i] = 1.0 + 0.04 * randn(rng)
        else
            random_X[i, :] .= (2.0, 1.1) .+ (0.16 * randn(rng), 0.17 * randn(rng))
            random_amp[i] = 3.0 + 0.04 * randn(rng)
        end
        if i % 23 == 0
            random_X[i, :] .+= (9.0, -8.0)
        end
    end
    random_fit = fit_student_t_two_component(random_X, vec(random_amp))
    @test random_fit.status == ROBUST_VALID
    predicted = [argmax(view(random_fit.responsibilities, i, :)) ==
                 random_fit.high_amplitude_component ? 1 : 0 for i in eachindex(states)]
    @test mean(predicted .== states) > 0.90
end

@testset "exact amplitude ties abstain under component relabeling" begin
    separated = reshape(vcat(
        [-3.20, -3.10, -3.00, -2.90, -2.80, -3.15, -2.95, -2.85],
        [2.80, 2.90, 3.00, 3.10, 3.20, 2.85, 2.95, 3.15],
    ), :, 1)
    tied_amplitudes = fill(2.0, size(separated, 1))
    tied_original = fit_student_t_two_component(separated, tied_amplitudes)
    tied_relabelled = fit_student_t_two_component(-separated, tied_amplitudes)

    for fit in (tied_original, tied_relabelled)
        @test fit.status == ROBUST_ABSTAINED
        @test fit.invalid_reason == AMPLITUDE_ORIENTATION_UNDEFINED
        @test fit.high_amplitude_component == 0
        @test fit.status != ROBUST_VALID
        hard = [fit.responsibilities[i, 1] >= fit.responsibilities[i, 2] ? 1 : 2
                for i in axes(fit.responsibilities, 1)]
        class_means = [mean(tied_amplitudes[hard .== component]) for component in 1:2]
        @test class_means == [2.0, 2.0]
    end
    @test tied_original.weights === EQUAL_MIXTURE_WEIGHTS
    @test tied_relabelled.weights === EQUAL_MIXTURE_WEIGHTS

    distinct_X, distinct_amplitudes = two_cluster_fixture()
    distinct_original = fit_student_t_two_component(distinct_X, distinct_amplitudes)
    distinct_relabelled = fit_student_t_two_component(-distinct_X, distinct_amplitudes)
    @test distinct_original.status == ROBUST_VALID
    @test distinct_relabelled.status == ROBUST_VALID
    original_physical = distinct_original.responsibilities[:,
                                                          distinct_original.high_amplitude_component]
    relabelled_physical = distinct_relabelled.responsibilities[:,
                                                              distinct_relabelled.high_amplitude_component]
    @test original_physical ≈ relabelled_physical atol=1e-12 rtol=1e-12
    @test distinct_original.high_amplitude_component ==
          3 - distinct_relabelled.high_amplitude_component
end

@testset "matched one-component null and fail-closed outcomes" begin
    rng = MersenneTwister(808)
    gaussian_like = reshape(randn(rng, 180), :, 1)
    unrelated_amplitude = 2.0 .+ 0.15 .* randn(rng, size(gaussian_like, 1))
    null_fit = fit_student_t_one_component(gaussian_like)
    split_fit = fit_student_t_two_component(gaussian_like, unrelated_amplitude)
    @test null_fit.status == ROBUST_VALID
    @test null_fit.invalid_reason == ROBUST_OK
    @test null_fit.converged
    @test null_fit.monotone
    @test null_fit.loglik == null_fit.loglik_trace[end]
    @test split_fit.status == ROBUST_ABSTAINED
    @test split_fit.invalid_reason in (ONE_COMPONENT_EVIDENCE,
                                       AMPLITUDE_ORIENTATION_UNDEFINED)
    @test split_fit.one_component.status == ROBUST_VALID

    empty_fit = fit_student_t_two_component(zeros(0, 1), Float64[])
    @test empty_fit.status == ROBUST_ABSTAINED
    @test empty_fit.invalid_reason == EMPTY_INPUT
    @test empty_fit.best_start == 0

    tiny_fit = fit_student_t_two_component(reshape([0.0], :, 1), [1.0])
    @test tiny_fit.status == ROBUST_ABSTAINED
    @test tiny_fit.invalid_reason == TINY_INPUT

    nonfinite_fit = fit_student_t_two_component(reshape([0.0, NaN, 1.0], :, 1),
                                                [1.0, 2.0, 3.0])
    @test nonfinite_fit.status == ROBUST_FAILED
    @test nonfinite_fit.invalid_reason == NONFINITE_INPUT

    nonfinite_amplitude = fit_student_t_two_component(reshape([-1.0, 1.0], :, 1),
                                                       [1.0, Inf])
    @test nonfinite_amplitude.status == ROBUST_FAILED
    @test nonfinite_amplitude.invalid_reason == NONFINITE_INPUT

    collapsed_fit = fit_student_t_two_component(fill(1.0, 20, 2),
                                                 collect(1.0:20.0))
    @test collapsed_fit.status == ROBUST_ABSTAINED
    @test collapsed_fit.invalid_reason == VARIANCE_COLLAPSE

    empty_null = fit_student_t_one_component(zeros(0, 2))
    @test empty_null.status == ROBUST_ABSTAINED
    @test empty_null.invalid_reason == EMPTY_INPUT

    @test_throws DimensionMismatch fit_student_t_two_component(zeros(3, 2), [1.0, 2.0])
    @test_throws ArgumentError fit_student_t_two_component(zeros(3, 0), [1.0, 2.0, 3.0])
end

const HASH_FIXTURE = two_cluster_fixture(contaminated=true)
const HASH_X = HASH_FIXTURE[1]
const HASH_AMP = HASH_FIXTURE[2]
const HASH_FIT = fit_student_t_two_component(HASH_X, HASH_AMP)
println("T7_RESULT_SHA256=", canonical_result_hash(HASH_FIT))
println("T7_TRACE_SHA256=", canonical_trace_hash(HASH_FIT))
println("T7_FIXED_NU=", FIXED_NU)
