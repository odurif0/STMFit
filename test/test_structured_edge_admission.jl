#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using SHA
using Statistics
using Test
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const EDGE_ADMISSION_SOURCE =
    joinpath(@__DIR__, "lib", "structured_assignment", "edge_admission.jl")
const EDGE_ADMISSION_CLI = joinpath(@__DIR__, "evaluate_structured_edge_admission.jl")

include(joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl"))
const HUA_BASELINE = HierarchicalUnitAssignment
include(joinpath(@__DIR__, "lib", "structured_assignment", "robust_emissions.jl"))
const SRE_BASELINE = StructuredRobustEmissions
const COMPLETE_REPLAY_SHA256 = Ref("")
const GRAPH_HANDOFF_REPLAY_SHA256 = Ref("")

@testset "reused class-free normalization and posterior characterization" begin
    names = ["amp_prominence", "amp_neighbor_ratio"]
    records = HUA_BASELINE.LobeRecord[
        HUA_BASELINE.LobeRecord("20260101_a.sxm", 1, 1.0, Dict(names[1] => 1.0, names[2] => 2.0)),
        HUA_BASELINE.LobeRecord("20260101_a.sxm", 2, 2.0, Dict(names[1] => 3.0, names[2] => 4.0)),
        HUA_BASELINE.LobeRecord("20260101_a.sxm", 3, 3.0, Dict(names[1] => 5.0, names[2] => 8.0)),
        HUA_BASELINE.LobeRecord("20260102_b.sxm", 1, 1.0, Dict(names[1] => 100.0, names[2] => 7.0)),
        HUA_BASELINE.LobeRecord("20260102_b.sxm", 2, 2.0, Dict(names[1] => 110.0, names[2] => 7.0)),
        HUA_BASELINE.LobeRecord("20260102_b.sxm", 3, 3.0, Dict(names[1] => 130.0, names[2] => 7.0)),
    ]
    X, valid = HUA_BASELINE.feature_matrix(records, names)
    @test all(valid)
    Z = HUA_BASELINE.normalize_per_scan(records, X, names)
    @test Z[1:3, 1] == [-1.0, 0.0, 1.0]
    @test Z[1:3, 2] == [-1.0, 0.0, 2.0]
    @test Z[4:6, 1] == [-1.0, 0.0, 2.0]
    @test all(isnan, Z[4:6, 2])

    means = [-1.0 0.5; 1.0 -0.5]
    scales = [0.8 1.2; 1.1 0.7]
    posterior = SRE_BASELINE.student_t_responsibilities([0.2 -0.1], means, scales)
    @test size(posterior) == (1, 2)
    @test sum(posterior[1, :]) ≈ 1.0 atol=2e-15 rtol=0
    @test SRE_BASELINE.FIXED_NU == 8
    @test SRE_BASELINE.EQUAL_MIXTURE_WEIGHTS == (0.5, 0.5)
end

if !isfile(EDGE_ADMISSION_SOURCE) || !isfile(EDGE_ADMISSION_CLI)
    @testset "structured edge admission RED: declared products are required" begin
        @test isfile(EDGE_ADMISSION_SOURCE)
        @test isfile(EDGE_ADMISSION_CLI)
    end
    exit(1)
end

include(EDGE_ADMISSION_SOURCE)
using .StructuredEdgeAdmission
const SEA = StructuredEdgeAdmission

function independent_full_log_t(x, mu, scale)
    delta = dot(x - mu, scale \ (x - mu))
    return log(24.0) - log(6.0) - log(8.0 * π) -
           0.5 * log(det(scale)) - 5.0 * log1p(delta / 8.0)
end

function independent_edge_log_ratios(residual, conditional_means, conditional_scale,
                                     null_mean, null_scale)
    null_density = independent_full_log_t(residual, null_mean, null_scale)
    return [independent_full_log_t(residual, vec(conditional_means[row, :]),
                                   conditional_scale) - null_density
            for row in 1:4]
end

function independent_edge_gain(residual, weights, conditional_means,
                               conditional_scale, null_mean, null_scale)
    ratios = independent_edge_log_ratios(
        residual, conditional_means, conditional_scale, null_mean, null_scale)
    return sum(Float64(weights[row]) * ratios[row] for row in 1:4), ratios
end

function independent_synthetic_edge_row_hash(edge)
    fields = [
        edge.file, string(edge.left_lobe), string(edge.right_lobe),
        @sprintf("%.17g", edge.left_t_nm), @sprintf("%.17g", edge.right_t_nm),
        @sprintf("%.17g", edge.right_t_nm - edge.left_t_nm), edge.segment_id,
        edge.segment_id, "eligible", "none", @sprintf("%.17g", edge.observation[1]),
        @sprintf("%.17g", edge.observation[2]), repeat("0", 64), repeat("0", 64),
        repeat("0", 64), repeat("0", 64), repeat("0", 64), repeat("0", 64),
    ]
    return bytes2hex(sha256(codeunits(join(fields, '\t') * "\n")))
end

function independent_two_node_joint(unary_left, unary_right, log_factors)
    weights = zeros(2, 2)
    for left_state in 0:1, right_state in 0:1
        left_prior = left_state == 1 ? unary_left : 1.0 - unary_left
        right_prior = right_state == 1 ? unary_right : 1.0 - unary_right
        weights[left_state + 1, right_state + 1] =
            left_prior * right_prior * exp(log_factors[left_state + 1, right_state + 1])
    end
    weights ./= sum(weights)
    marginals = [sum(weights[row, :]) for row in 1:2],
                [sum(weights[:, column]) for column in 1:2]
    return weights, marginals
end

function independent_null_proposal(residuals, mean0, scale0)
    n = size(residuals, 1)
    latent = [10.0 / (8.0 + dot(view(residuals, i, :) - mean0,
                                  scale0 \ (view(residuals, i, :) - mean0)))
              for i in 1:n]
    mean1 = vec(sum(latent .* residuals; dims=1)) / sum(latent)
    raw = zeros(2, 2)
    for i in 1:n
        difference = vec(residuals[i, :]) - mean1
        raw .+= latent[i] .* (difference * difference')
    end
    return mean1, raw / n
end

function independent_conditional_proposal(residuals, weights, means, scale, null_scale)
    n = size(residuals, 1)
    latent = zeros(n, 4)
    for i in 1:n, category in 1:4
        difference = vec(residuals[i, :]) - vec(means[category, :])
        latent[i, category] = 10.0 / (8.0 + dot(difference, scale \ difference))
    end
    updated = zeros(4, 2)
    for category in (1, 4)
        weighted = weights[:, category] .* latent[:, category]
        updated[category, :] .= vec(sum(weighted .* residuals; dims=1)) ./ sum(weighted)
    end
    mixed_weight = weights[:, 2] .* latent[:, 2] .+ weights[:, 3] .* latent[:, 3]
    mixed = vec(sum(mixed_weight .* residuals; dims=1)) ./ sum(mixed_weight)
    updated[2, :] .= mixed
    updated[3, :] .= mixed
    raw = zeros(2, 2)
    denominator = sum(weights)
    for i in 1:n, category in 1:4
        difference = vec(residuals[i, :]) - vec(updated[category, :])
        raw .+= weights[i, category] * latent[i, category] .* (difference * difference')
    end
    raw ./= denominator
    lambda = min(0.5, 3.0 / n)
    return updated, (1.0 - lambda) .* raw .+ lambda .* null_scale, lambda
end

function old_project_scale_for_test(matrix)
    symmetric = (Matrix{Float64}(matrix) + Matrix{Float64}(matrix)') / 2.0
    decomposition = eigen(Symmetric(symmetric))
    eigenvalues = max.(Float64.(decomposition.values), SEA.SCALE_FLOOR)
    projected = decomposition.vectors * Diagonal(eigenvalues) * decomposition.vectors'
    projected = (projected + projected') / 2.0
    condition_number = maximum(eigenvalues) / minimum(eigenvalues)
    return projected, eigenvalues, condition_number
end

function published_scale(row, columns, prefix)
    fields = [row[columns["$(prefix)_scale_$(suffix)"]] for suffix in ("ff", "fb", "bb")]
    if all(==("NA"), fields)
        return nothing
    end
    all(!=("NA"), fields) || error("partial numeric $prefix scale")
    values = parse.(Float64, fields)
    return [values[1] values[2]; values[2] values[3]], fields
end

function assert_certified_published_scale(row, columns, prefix)
    parsed = published_scale(row, columns, prefix)
    parsed === nothing && return
    scale, fields = parsed
    certificate = SEA._stored_scale_certificate(scale)
    @test certificate.status == :ok
    @test certificate.eigenvalues == eigvals(Symmetric(scale))
    @test certificate.condition_number ==
          maximum(certificate.eigenvalues) / minimum(certificate.eigenvalues)
    @test fields == [@sprintf("%.17g", scale[1, 1]),
                     @sprintf("%.17g", scale[1, 2]),
                     @sprintf("%.17g", scale[2, 2])]
end

@testset "fixed-nu full-scale density and projection" begin
    @test SEA.PLAN_SHA256 ==
          "a9b386d613829e8f7e20b6e33f8e80898fa9a55f0b344dbb92ea64ac6804f3d0"
    @test SEA.FIXED_NU == 8
    @test SEA.SCALE_FLOOR == 1.0e-4
    @test SEA.CONDITION_CAP == 1.0e4
    @test SEA.MAX_ITER == 200
    @test SEA.CONVERGENCE_TOL == 1.0e-6
    @test SEA.DECREASE_TOL == 1.0e-10
    @test SEA.START_ALPHAS == (0.0, 0.25, 0.5, 0.75, 1.0)

    x = [1.0, -2.0]
    mu = [0.25, -1.5]
    scale = [0.5 0.1; 0.1 3.0]
    @test SEA.log_student_t_full(x, mu, scale) ≈ -2.7722364321680395 atol=2e-15 rtol=0
    @test SEA.log_student_t_full(x, mu, scale) ≈ independent_full_log_t(x, mu, scale) atol=2e-15 rtol=0

    projected = SEA._project_scale([0.5 0.03; 0.01 5.0e-5])
    @test projected.status == :ok
    @test projected.scale == projected.scale'
    @test minimum(eigvals(Symmetric(projected.scale))) >= SEA.SCALE_FLOOR
    @test cond(projected.scale, 2) <= SEA.CONDITION_CAP
    @test SEA._project_scale([NaN 0.0; 0.0 1.0]).status == :nonfinite
    @test SEA._project_scale([2.0 0.0; 0.0 SEA.SCALE_FLOOR]).status == :ill_conditioned
end

@testset "canonical stored scale certification and deterministic closure" begin
    offending = [
        0.0045253227996870935 -0.0016497531701978115
        -0.0016497531701978115 0.00071502530906223933
    ]
    projected = SEA._project_scale(offending)
    @test projected.status == :ok
    @test projected.scale == projected.scale'
    @test minimum(projected.eigenvalues) >= SEA.SCALE_FLOOR
    @test projected.condition_number <= SEA.CONDITION_CAP
    @test projected.scale == [
        parse(Float64, "0.0045253227996870926") parse(Float64, "-0.001649753170197811")
        parse(Float64, "-0.001649753170197811") parse(Float64, "0.00071502530906224009")
    ]
    @test projected.eigenvalues == eigvals(Symmetric(projected.scale))
    @test projected.condition_number ==
          maximum(projected.eigenvalues) / minimum(projected.eigenvalues)

    exact_diagonal = SEA._project_scale([SEA.SCALE_FLOOR 0.0; 0.0 1.0])
    @test exact_diagonal.status == :ok
    @test exact_diagonal.scale == [SEA.SCALE_FLOOR 0.0; 0.0 1.0]
    @test exact_diagonal.eigenvalues == [SEA.SCALE_FLOOR, 1.0]
    @test exact_diagonal.condition_number == SEA.CONDITION_CAP

    previous_floor = SEA._project_scale([prevfloat(SEA.SCALE_FLOOR) 0.0;
                                         0.0 prevfloat(SEA.SCALE_FLOOR)])
    @test previous_floor.status == :ok
    @test minimum(previous_floor.eigenvalues) >= SEA.SCALE_FLOOR

    rotated_floor = SEA._project_scale([
        0.00500005 0.00490005
        0.00490005 0.00500005
    ])
    @test rotated_floor.status == :ok
    @test minimum(rotated_floor.eigenvalues) >= SEA.SCALE_FLOOR

    exact_cap_rotated = SEA._project_scale([
        0.50005 0.49995
        0.49995 0.50005
    ])
    @test exact_cap_rotated.status == :ok
    @test exact_cap_rotated.condition_number <= SEA.CONDITION_CAP

    over_cap = SEA._project_scale([2.0 0.0; 0.0 SEA.SCALE_FLOOR])
    @test over_cap.status == :ill_conditioned
    @test SEA._stored_scale_certificate([0.0 0.0; 0.0 1.0]).status == :nonpositive
    @test SEA._stored_scale_certificate([NaN 0.0; 0.0 1.0]).status == :nonfinite
    @test SEA._closure_projection([floatmax(Float64) 0.0; 0.0 floatmax(Float64)]).status == :no_progress

    for matrix in ([0.5 0.03; 0.03 3.0], [1.2 -0.2; -0.2 0.8],
                   [0.25 0.01; 0.01 0.7])
        old_scale, old_eigenvalues, old_condition = old_project_scale_for_test(matrix)
        current = SEA._project_scale(matrix)
        @test current.status == :ok
        @test current.scale == old_scale
        @test current.eigenvalues == eigvals(Symmetric(current.scale))
        @test current.condition_number ==
              maximum(current.eigenvalues) / minimum(current.eigenvalues)
        @test old_condition <= SEA.CONDITION_CAP
        @test all(old_eigenvalues .> 0.0)
    end
end

@testset "null and conditional exact update order" begin
    residuals = [
        -0.8 0.2
        -0.2 0.9
         0.4 -0.5
         1.1 0.3
    ]
    mean0 = [0.15, 0.1]
    scale0 = [0.9 0.12; 0.12 0.7]
    expected_mean, expected_raw = independent_null_proposal(residuals, mean0, scale0)
    proposal = SEA._null_proposal(residuals, mean0, scale0)
    @test proposal.status == :ok
    @test proposal.mean ≈ expected_mean atol=2e-15 rtol=2e-15
    @test proposal.unprojected_scale ≈ expected_raw atol=2e-15 rtol=2e-15
    @test proposal.scale ≈ SEA._project_scale(expected_raw).scale atol=2e-15 rtol=2e-15

    weights = [
        0.54 0.06 0.36 0.04
        0.12 0.48 0.08 0.32
        0.05 0.05 0.45 0.45
        0.02 0.18 0.08 0.72
    ]
    @test all(sum(weights; dims=2) .== 1.0)
    means = [-0.6 0.1; 0.0 0.3; 0.0 0.3; 0.8 -0.1]
    expected_means, expected_shrunk, expected_lambda =
        independent_conditional_proposal(residuals, weights, means, scale0, scale0)
    conditional = SEA._conditional_proposal(residuals, weights, means, scale0, scale0)
    @test conditional.status == :ok
    @test conditional.means ≈ expected_means atol=3e-15 rtol=3e-15
    @test conditional.means[2, :] == conditional.means[3, :]
    @test conditional.unprojected_scale ≈
          (expected_shrunk .- expected_lambda .* scale0) ./ (1.0 - expected_lambda) atol=3e-15 rtol=3e-15
    @test conditional.shrinkage == expected_lambda
    @test conditional.shrunk_scale ≈ expected_shrunk atol=3e-15 rtol=3e-15
    @test conditional.scale ≈ SEA._project_scale(expected_shrunk).scale atol=3e-15 rtol=3e-15
end

@testset "five starts, accepted objective, and failure semantics" begin
    residuals = reduce(vcat, ([sin(0.31i) cos(0.23i)] for i in 1:120))
    probabilities = [0.08 + 0.84 * (isodd(i) ? 1.0 : 0.0) for i in 1:121]
    weights = reduce(vcat, (reshape(collect(SEA.pair_weights(probabilities[i], probabilities[i + 1])), 1, :)
                            for i in 1:120))
    null_fit = SEA.fit_state_independent(residuals)
    @test null_fit.status == :ok
    @test null_fit.converged
    @test all(null_fit.objective_trace[i] >= null_fit.objective_trace[i - 1] - SEA.DECREASE_TOL
              for i in 2:length(null_fit.objective_trace))

    initial = SEA._conditional_initial_states(residuals, weights, null_fit)
    @test Tuple(state.alpha for state in initial) == SEA.START_ALPHAS
    @test all(state.means[2, :] == state.means[3, :] for state in initial)
    @test initial[1].means == repeat(null_fit.mean', 4, 1)

    conditional = SEA.fit_conditional(residuals, weights, null_fit)
    @test length(conditional.starts) == 5
    @test Tuple(start.alpha for start in conditional.starts) == SEA.START_ALPHAS
    @test conditional.best_start in 1:5
    @test conditional.means[2, :] == conditional.means[3, :]
    @test all(start -> all(start.objective_trace[i] >= start.objective_trace[i - 1] - SEA.DECREASE_TOL
                           for i in 2:length(start.objective_trace)), conditional.starts)

    @test SEA._select_best_start([1.0, 1.0 + 0.5e-12, 0.0, -1.0, -2.0], trues(5)) == 1
    @test SEA._select_best_start([1.0, 1.0 + 2.0e-12, 0.0, -1.0, -2.0], trues(5)) == 2
    accepted = SEA._accept_objective(10.0, 10.0 - 0.5e-10)
    rejected = SEA._accept_objective(10.0, 10.0 - 2.0e-10)
    @test accepted == :accepted
    @test rejected == :nonmonotone
    @test SEA.fit_state_independent(residuals; iteration_limit=0).reason == :nonconverged
    @test SEA.fit_conditional(residuals, weights, null_fit; iteration_limit=0).reason == :no_valid_start
    @test SEA.fit_state_independent([0.0 0.0; NaN 1.0]).reason == :nonfinite
end

@testset "ridge residualizer and untouched held-out gain" begin
    endpoint_rows = [
        collect(range(-1.4, 1.4; length=14)) .+ 0.07i for i in 1:32
    ]
    design = hcat(ones(32), reduce(vcat, (row' for row in endpoint_rows)))
    coefficients = reshape([0.12 + 0.003i for i in 1:30], 15, 2)
    outputs = design * coefficients
    fit = SEA.fit_residualizer(reduce(vcat, (row' for row in endpoint_rows)), outputs)
    ridge = max(1.0e-12, 1.0e-6 * tr(design[:, 2:end]' * design[:, 2:end]) / 14)
    penalty = Diagonal(vcat(0.0, fill(ridge, 14)))
    expected = (design' * design + penalty) \ (design' * outputs)
    @test fit.ridge == ridge
    @test fit.coefficients ≈ expected atol=5e-14 rtol=5e-14

    held_x = collect(range(-0.9, 0.9; length=14))
    held_y = [0.33, -0.21]
    untouched = copy(held_y)
    residual = SEA.residualize(fit, held_x, held_y)
    @test residual ≈ held_y - vec(vcat(1.0, held_x)' * expected) atol=3e-15 rtol=3e-15
    @test held_y == untouched

    conditional_means = [-0.4 0.1; 0.2 0.5; 0.2 0.5; 0.6 -0.2]
    conditional_scale = [0.7 0.08; 0.08 0.9]
    null_mean = [0.05, 0.0]
    null_scale = [1.1 0.02; 0.02 0.8]
    edge_weights = collect(SEA.pair_weights(0.3, 0.8))
    expected_gain = sum(edge_weights[category] *
                        (independent_full_log_t(residual, vec(conditional_means[category, :]), conditional_scale) -
                         independent_full_log_t(residual, null_mean, null_scale))
                        for category in 1:4)
    @test SEA.edge_gain(residual, edge_weights, conditional_means, conditional_scale,
                        null_mean, null_scale) ≈ expected_gain atol=3e-15 rtol=3e-15
end

@testset "complete-vector cyclic shift and grouped aggregation" begin
    file = "20260101_phase.sxm"
    segment = "0123456789abcdef"
    seed = 17
    digest = sha256(codeunits(string(seed, '\0', file, '\0', segment)))
    integer = foldl((value, byte) -> (value << 8) | UInt64(byte), digest[1:8]; init=UInt64(0))
    expected_offset = 1 + Int(integer % UInt64(4))
    @test SEA.cyclic_offset(seed, file, segment, 5) == expected_offset
    @test 1 <= expected_offset <= 4

    vectors = [
        (0.7, 0.1, 0.15, 0.05),
        (0.1, 0.2, 0.3, 0.4),
        (0.25, 0.25, 0.25, 0.25),
        (0.05, 0.15, 0.2, 0.6),
        (0.4, 0.3, 0.2, 0.1),
    ]
    rotated = SEA.rotate_complete_weights(vectors, expected_offset)
    @test sort(collect(rotated)) == sort(collect(vectors))
    @test all(abs(sum(vector) - 1.0) <= SEA.WEIGHT_SUM_ROUNDOFF for vector in rotated)
    @test SEA.rotate_complete_weights(vectors, 5) == vectors

    scan_scores = Dict(
        "20260101" => Dict("a" => 1.0, "b" => 3.0),
        "20260102" => Dict("c" => 4.0, "d" => 6.0),
    )
    aggregate = SEA.aggregate_scan_date(scan_scores)
    @test aggregate.date_means == Dict("20260101" => 2.0, "20260102" => 5.0)
    @test aggregate.overall == 3.5
    constant_scores = Dict(
        "20260101" => Dict("a" => 2.0, "b" => 2.0),
        "20260102" => Dict("c" => 5.0, "d" => 5.0),
    )
    @test SEA.whole_scan_bootstrap(constant_scores, 0) == 3.5
    @test length(SEA.bootstrap_distribution(constant_scores)) == 500
    @test all(==(3.5), SEA.bootstrap_distribution(constant_scores))
end

@testset "ESS, reversal, terminal states, and pre-shuffle blocking" begin
    weights = [
        0.5 0.2 0.2 0.1
        0.4 0.3 0.2 0.1
        0.3 0.3 0.3 0.1
    ]
    expected = [(sum(weights[:, category])^2) / sum(abs2, weights[:, category])
                for category in 1:4]
    @test SEA.kish_ess(weights) ≈ expected atol=2e-15 rtol=2e-15
    @test SEA.reversal_equivalent([0.2, -0.1], [0.2 + 0.5e-12, -0.1])
    @test !SEA.reversal_equivalent([0.2, -0.1], [0.2 + 2.0e-12, -0.1])

    @test SEA.combine_terminal_states(["PASS", "PASS"]) == "PASS"
    @test SEA.combine_terminal_states(["PASS", "SKIPPED"]) == "SKIPPED"
    @test SEA.combine_terminal_states(["PASS", "FAIL"]) == "FAIL"
    @test SEA.combine_terminal_states(["PASS", "BLOCKED"]) == "BLOCKED"

    for term in (:slot, :terminal, :gap, :position, :transition)
        before = SEA._SHUFFLE_CALL_COUNT[]
        error = try
            SEA.check_shuffle_terms((term,))
            nothing
        catch caught
            caught
        end
        @test error isa SEA.EdgeAdmissionError
        @test error.code == :forbidden_shuffle_term
        @test SEA._SHUFFLE_CALL_COUNT[] == before
    end
    @test SEA.check_shuffle_terms(()) === nothing
end

function synthetic_admission_fixture(; mechanism::Symbol=:informative)
    patterns = (
        (0, 0, 1, 0, 1, 1),
        (1, 1, 0, 1, 0, 0),
        (0, 1, 1, 0, 0, 1),
        (1, 0, 0, 1, 1, 0),
    )
    nodes = SEA.AdmissionNode[]
    hidden = Dict{Tuple{String,Int},Int}()
    for date_index in 1:7
        date = @sprintf("202601%02d", date_index)
        for scan_index in 1:4
            file = @sprintf("%s_scan%02d.sxm", date, scan_index)
            pattern = patterns[scan_index]
            scan_shift = 0.17 * sin(0.4 * date_index + 0.8 * scan_index)
            for lobe in 1:6
                state = pattern[lobe]
                hidden[(file, lobe)] = state
                signal = state == 1 ? 2.3 : -2.3
                predictors = ntuple(7) do feature
                    signal * (1.0 + 0.025 * feature) + scan_shift +
                    0.11 * sin(0.31 * lobe + 0.17 * feature + 0.13 * scan_index)
                end
                amplitude = 1.0 + 2.0 * state + 0.03 * cos(0.5 * lobe + scan_index)
                push!(nodes, SEA.AdmissionNode(file, date, lobe,
                                               0.5 * (lobe - 1), amplitude,
                                               predictors))
            end
        end
    end
    sort!(nodes; by=node -> (node.file, node.t_nm, node.lobe))
    hashes = Dict(
        "candidate_config_sha256" => repeat("1", 64),
        "model_config_sha256" => repeat("2", 64),
        "universe_receipt_sha256" => repeat("3", 64),
        "edge_receipt_sha256" => repeat("4", 64),
        "input_sha256" => repeat("5", 64),
    )
    shell = SEA.AdmissionData(nodes, SEA.AdmissionEdge[], hashes)
    normalized = SEA.normalize_admission_data(shell)
    true_coefficients = zeros(15, 2)
    true_coefficients[1, :] .= (0.035, -0.025)
    for row in 2:15
        true_coefficients[row, 1] = 0.0025 * sin(0.37 * row)
        true_coefficients[row, 2] = 0.0020 * cos(0.29 * row)
    end
    category_means = Dict(
        (0, 0) => [-0.31, 0.12],
        (0, 1) => [0.34, -0.10],
        (1, 0) => [0.34, -0.10],
        (1, 1) => [-0.25, 0.15],
    )
    edges = SEA.AdmissionEdge[]
    for file in sort!(unique(node.file for node in nodes))
        chain = filter(node -> node.file == file, nodes)
        segment = bytes2hex(sha256(codeunits("synthetic-segment\nfile=$file\n")))[1:16]
        scan_index = parse(Int, match(r"scan(\d+)", file).captures[1])
        for ordinal in 1:5
            left = chain[ordinal]
            right = chain[ordinal + 1]
            left_index = normalized.node_index[(file, left.lobe)]
            right_index = normalized.node_index[(file, right.lobe)]
            endpoint = vcat(collect(normalized.normalized_predictors[left_index, :]),
                            collect(normalized.normalized_predictors[right_index, :]))
            baseline = vec(vcat(1.0, endpoint)' * true_coefficients)
            category = (hidden[(file, left.lobe)], hidden[(file, right.lobe)])
            date_index = parse(Int, file[7:8])
            conditional = if mechanism == :informative
                category_means[category]
            elseif mechanism == :endpoint_explained
                [0.0, 0.0]
            elseif mechanism == :state_independent
                [0.075 * sin(0.43 * ordinal + 0.61 * date_index + 0.29 * scan_index),
                 0.070 * cos(0.37 * ordinal + 0.47 * date_index + 0.31 * scan_index)]
            elseif mechanism == :scan_artifact
                [0.18 * sin(0.73 * date_index + 1.17 * scan_index),
                 -0.16 * cos(0.91 * date_index + 1.31 * scan_index)]
            elseif mechanism == :phase_shuffled
                shifted_left = hidden[(file, mod1(left.lobe + scan_index, 6))]
                shifted_right = hidden[(file, mod1(right.lobe + scan_index, 6))]
                category_means[(shifted_left, shifted_right)]
            else
                error("unknown synthetic mechanism")
            end
            noise = [0.018 * sin(0.71 * ordinal + 0.19 * scan_index),
                     0.017 * cos(0.53 * ordinal + 0.23 * scan_index)]
            observation = baseline .+ conditional .+ noise
            @assert all(abs.(observation) .< 1.0)
            push!(edges, SEA.AdmissionEdge(
                file,
                left.date,
                left.lobe,
                right.lobe,
                left.t_nm,
                right.t_nm,
                segment,
                ordinal,
                (observation[1], observation[2]),
            ))
        end
    end
    sort!(edges; by=edge -> (edge.file, edge.left_t_nm, edge.right_t_nm,
                             edge.left_lobe, edge.right_lobe))
    return SEA.AdmissionData(nodes, edges, hashes)
end

function final_gate_for_mechanism(mechanism::Symbol)
    data = synthetic_admission_fixture(mechanism=mechanism)
    normalized = SEA.normalize_admission_data(data)
    dates = sort!(unique(node.date for node in data.nodes))
    cache = Dict{String,SEA.ConditionalFit}()
    evaluations = SEA.PartitionEvaluation[]
    for heldout in dates
        training = [date for date in dates if date != heldout]
        fit = SEA.fit_partition_model(normalized, "C1", training)
        reversed_fit = SEA.fit_partition_model(normalized, "C1", training;
                                               reversed=true)
        observed = SEA.evaluate_partition(fit, heldout; run_shuffle=true,
                                          shuffle_fit_cache=cache)
        reversed = SEA.evaluate_partition(reversed_fit, heldout;
                                          run_shuffle=false)
        push!(evaluations, SEA._with_reversal(observed, reversed))
    end
    return SEA.evaluate_gate(evaluations)
end

function write_fixture_file(root::String, relative::String, bytes)
    path = joinpath(root, split(relative, '/')...)
    mkpath(dirname(path))
    write(path, bytes)
    return path
end


function copy_cli_dependencies(root::String)
    relatives = String[
        ".omo/plans/structured-label-free-unit-assignment.md",
        ".omo/evidence/structured-label-free-unit-assignment/t2/correction2/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase3-t3-universe/correction/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/t7/correction/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/review/AdversarialVerify.json",
        "config/chitosan.toml",
        "config/unit_assignment_structured_candidate.toml",
        "config/unit_assignment_structured_model.toml",
        "test/build_label_free_edge_features.jl",
        "test/evaluate_structured_edge_admission.jl",
        "test/extract_lobe_patches.jl",
        "test/extract_lobe_patches_bwd.jl",
        "test/test_structured_edge_features.jl",
        "test/test_structured_robust_emissions.jl",
        "test/lib/hierarchical_unit_assignment.jl",
        "test/lib/structured_assignment/edge_admission.jl",
        "test/lib/structured_assignment/edge_features.jl",
        "test/lib/structured_assignment/firewall.jl",
        "test/lib/structured_assignment/robust_emissions.jl",
        "test/lib/structured_assignment/universe.jl",
    ]
    append!(relatives, [relpath(path, ROOT) for path in
                        readdir(joinpath(ROOT, "test", "lib", "hierarchical");
                                join=true) if isfile(path)])
    for relative in relatives
        destination = joinpath(root, split(relative, '/')...)
        mkpath(dirname(destination))
        cp(joinpath(ROOT, split(relative, '/')...), destination; force=true)
    end
    return nothing
end

function toml_line(key::String, value)
    return value isa Integer ? "$key = $value" : "$key = $(repr(String(value)))"
end

file_sha256(path::String) = bytes2hex(sha256(read(path)))

function synthetic_original_feature_bytes(data::SEA.AdmissionData)::Vector{UInt8}
    header = ["file", "N", "lobe", "amplitude", "t_nm", "u_nm"]
    counts = Dict(file => count(node -> node.file == file, data.nodes)
                  for file in unique(node.file for node in data.nodes))
    rows = [join([
        node.file,
        string(counts[node.file]),
        string(node.lobe),
        @sprintf("%.17g", node.amplitude),
        @sprintf("%.17g", node.t_nm),
        "0",
    ], '\t') for node in data.nodes]
    return Vector{UInt8}(codeunits(join(header, '\t') * "\n" *
                                   join(rows, "\n") * "\n"))
end

function synthetic_patch_bytes(
    data::SEA.AdmissionData,
    channel::Symbol,
)::Vector{UInt8}
    header = channel == :forward ?
        SEA.StructuredEdgeFeatures.FORWARD_PATCH_HEADER :
        SEA.StructuredEdgeFeatures.BACKWARD_PATCH_HEADER
    payload = fill("0", length(header) - 5)
    rows = [join(vcat([
        node.file,
        string(node.lobe),
        @sprintf("%.17g", node.t_nm),
        "0",
        @sprintf("%.17g", node.amplitude),
    ], payload), '\t') for node in data.nodes]
    return Vector{UInt8}(codeunits(join(header, '\t') * "\n" *
                                   join(rows, "\n") * "\n"))
end

function synthetic_patch_receipt_text(
    root::String,
    channel::Symbol,
    features::String,
    data_dir::String,
    output::String,
    feature_sha256::String,
    keys_sha256::String,
    row_count::Int,
    command::String,
)::String
    source = channel == :forward ?
        joinpath(root, "test", "extract_lobe_patches.jl") :
        joinpath(root, "test", "extract_lobe_patches_bwd.jl")
    header = channel == :forward ?
        SEA.StructuredEdgeFeatures.FORWARD_PATCH_HEADER :
        SEA.StructuredEdgeFeatures.BACKWARD_PATCH_HEADER
    fields = [
        toml_line("schema", "stmfit-structured-patch-producer-receipt-v1"),
        "schema_version = 1",
        toml_line("status", "PASS"),
        toml_line("channel", String(channel)),
        toml_line("execution_site", "Viper_Slurm_compute_node"),
        "julia_threads = 1",
        toml_line("root", realpath(root)),
        toml_line("cwd", realpath(root)),
        toml_line("features_path", realpath(features)),
        toml_line("data_path", realpath(data_dir)),
        toml_line("output_path", realpath(output)),
        toml_line("source_path", realpath(source)),
        toml_line("source_sha256", file_sha256(source)),
        toml_line("feature_sha256", feature_sha256),
        toml_line("output_sha256", file_sha256(output)),
        toml_line("candidate_config_sha256", SEA.CANDIDATE_CONFIG_SHA256),
        toml_line("model_config_sha256", SEA.MODEL_CONFIG_SHA256),
        toml_line("grid_sha256", SEA.StructuredEdgeFeatures.GRID_SHA256),
        toml_line("keys_sha256", keys_sha256),
        toml_line("table_header_sha256",
                  bytes2hex(sha256(codeunits(join(header, '\t') * "\n")))),
        "row_count = $row_count",
        toml_line("command", command),
    ]
    return join(fields, "\n") * "\n"
end

function replace_toml_value!(path::String, key::String, value)
    lines = split(chomp(read(path, String)), '\n')
    prefix = "$key = "
    indices = findall(line -> startswith(line, prefix), lines)
    length(indices) == 1 || error("expected one TOML field $key")
    lines[only(indices)] = toml_line(key, value)
    write(path, join(lines, "\n") * "\n")
    return nothing
end

function replace_edge_binding!(fixture, column::Int, value::String)
    edge_path = joinpath(fixture.edges, "edge_observations.tsv")
    lines = split(chomp(read(edge_path, String)), '\n')
    fields = String.(split(lines[2], '\t'; keepempty=true))
    fields[column] = value
    lines[2] = join(fields, '\t')
    write(edge_path, join(lines, "\n") * "\n")
    replace_toml_value!(joinpath(fixture.edges, "receipt.toml"),
                        "edge_observations_sha256", file_sha256(edge_path))
    return nothing
end

function load_cli_fixture(fixture)
    return SEA.load_admission_data(
        fixture.root;
        features=relpath(fixture.features, fixture.root),
        candidate_config="config/unit_assignment_structured_candidate.toml",
        model_config="config/unit_assignment_structured_model.toml",
        universe_dir=relpath(fixture.universe, fixture.root),
        edge_dir=relpath(fixture.edges, fixture.root),
        forward_receipt=relpath(fixture.forward_receipt, fixture.root),
        backward_receipt=relpath(fixture.backward_receipt, fixture.root),
    )
end

function expect_pre_shuffle_rejection(mutation::Function, fixture)
    edge_path = joinpath(fixture.edges, "edge_observations.tsv")
    receipt_path = joinpath(fixture.edges, "receipt.toml")
    original_edge = read(edge_path)
    original_receipt = read(receipt_path)
    before = SEA._SHUFFLE_CALL_COUNT[]
    error = try
        mutation()
        load_cli_fixture(fixture)
        nothing
    catch caught
        caught
    finally
        write(edge_path, original_edge)
        write(receipt_path, original_receipt)
    end
    @test error isa SEA.EdgeAdmissionError
    @test SEA._SHUFFLE_CALL_COUNT[] == before
    return error
end

function expect_cli_pass(fixture)
    output = joinpath(fixture.root, "synthetic", "valid-cli-output")
    try
        command = `$(Base.julia_cmd()) --project=$(ROOT) $(EDGE_ADMISSION_CLI) --root $(fixture.root) --features $(relpath(fixture.features, fixture.root)) --candidate-config config/unit_assignment_structured_candidate.toml --model-config config/unit_assignment_structured_model.toml --universe-dir $(relpath(fixture.universe, fixture.root)) --edge-dir $(relpath(fixture.edges, fixture.root)) --forward-receipt $(relpath(fixture.forward_receipt, fixture.root)) --backward-receipt $(relpath(fixture.backward_receipt, fixture.root)) --out-dir $(relpath(output, fixture.root))`
        stdout = IOBuffer()
        stderr = IOBuffer()
        process = run(pipeline(command; stdout=stdout, stderr=stderr); wait=false)
        wait(process)
        stdout_text = String(take!(stdout))
        stderr_text = String(take!(stderr))
        @test process.exitcode == 0
        @test occursin("status=PASS", stdout_text)
        @test occursin(r"result_sha256=[0-9a-f]{64}", stdout_text)
        @test isempty(stderr_text)
        @test sort(readdir(output)) == [
            "bootstrap.tsv", "ess.tsv", "fitted_edge_models.tsv",
            "fitted_edge_transforms.tsv", "fitted_unary_nodes.tsv", "models.tsv",
            "partitions.tsv", "receipt.toml", "scores.tsv", "shuffle.tsv", "starts.tsv",
        ]
        receipt = TOML.parsefile(joinpath(output, "receipt.toml"))
        @test receipt["status"] == "PASS"
        @test receipt["schema"] == "stmfit-structured-edge-admission-receipt-v2"
        @test receipt["graph_handoff_schema"] == SEA.GRAPH_HANDOFF_SCHEMA
        @test receipt["hashes"]["plan_sha256"] == SEA.PLAN_SHA256
        @test occursin(r"^[0-9a-f]{64}$", receipt["result_sha256"])
        @test Set(keys(receipt["artifacts"])) ==
              Set(setdiff(readdir(output), ["receipt.toml"]))
        for (name, digest) in receipt["artifacts"]
            @test file_sha256(joinpath(output, name)) == digest
        end
        @test isempty(filter(name -> startswith(name, ".structured-admission-stage-"),
                             readdir(dirname(output))))
    finally
        rm(output; recursive=true, force=true)
    end
    return nothing
end

function expect_cli_blocked(
    mutation::Function,
    fixture,
    name::String,
    expected_reason::String,
)
    edge_path = joinpath(fixture.edges, "edge_observations.tsv")
    receipt_path = joinpath(fixture.edges, "receipt.toml")
    original_edge = read(edge_path)
    original_receipt = read(receipt_path)
    output = joinpath(fixture.root, "synthetic", "blocked-$name")
    try
        mutation()
        command = `$(Base.julia_cmd()) --compile=min --project=$(ROOT) $(EDGE_ADMISSION_CLI) --root $(fixture.root) --features $(relpath(fixture.features, fixture.root)) --candidate-config config/unit_assignment_structured_candidate.toml --model-config config/unit_assignment_structured_model.toml --universe-dir $(relpath(fixture.universe, fixture.root)) --edge-dir $(relpath(fixture.edges, fixture.root)) --forward-receipt $(relpath(fixture.forward_receipt, fixture.root)) --backward-receipt $(relpath(fixture.backward_receipt, fixture.root)) --out-dir $(relpath(output, fixture.root))`
        stdout = IOBuffer()
        stderr = IOBuffer()
        process = run(pipeline(command; stdout=stdout, stderr=stderr); wait=false)
        wait(process)
        @test process.exitcode == 2
        @test isempty(String(take!(stdout)))
        @test occursin("BLOCKED", String(take!(stderr)))
        @test readdir(output) == ["receipt.toml"]
        blocked = TOML.parsefile(joinpath(output, "receipt.toml"))
        @test blocked["status"] == "BLOCKED"
        @test blocked["reason"] == expected_reason
        @test isempty(filter(name -> startswith(name, ".structured-admission-stage-"),
                             readdir(dirname(output))))
    finally
        write(edge_path, original_edge)
        write(receipt_path, original_receipt)
        rm(output; recursive=true, force=true)
    end
    return nothing
end

different_sha256(value::String) =
    value == repeat("0", 64) ? repeat("1", 64) : repeat("0", 64)

function make_cli_fixture(; forbidden_column::Bool=false, cleanup::Bool=true)
    root = mktempdir(; prefix="stmfit-edge-admission-cli-", cleanup=cleanup)
    copy_cli_dependencies(root)
    data = synthetic_admission_fixture()
    original_features = write_fixture_file(
        root,
        "synthetic/original-features.tsv",
        synthetic_original_feature_bytes(data),
    )
    feature_sha = file_sha256(original_features)
    source_sha = bytes2hex(sha256(codeunits("synthetic-universe-source-v1\n")))
    structured_keys = Tuple(
        SEA.StructuredEdgeFeatures.StructuredUniverse.LobeKey(node.file, node.lobe)
        for node in data.nodes
    )
    keys_sha = bytes2hex(sha256(
        SEA.StructuredEdgeFeatures.StructuredUniverse._key_identity_bytes(
            structured_keys,
        ),
    ))

    scan_header = [
        "schema_version", "file", "date", "lobe_count", "feature_sha256",
        "candidate_config_sha256", "model_config_sha256", "source_sha256",
        "forward_producer_sha256", "backward_producer_sha256", "grid_sha256",
    ]
    scan_rows = String[]
    for file in sort!(unique(node.file for node in data.nodes))
        date = only(unique(node.date for node in data.nodes if node.file == file))
        push!(scan_rows, join([
            "1", file, date, "6", feature_sha, SEA.CANDIDATE_CONFIG_SHA256,
            SEA.MODEL_CONFIG_SHA256, source_sha,
            SEA.StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256,
            SEA.StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256,
            SEA.StructuredEdgeFeatures.GRID_SHA256,
        ], '\t'))
    end
    scan_bytes = Vector{UInt8}(codeunits(join(scan_header, '\t') * "\n" *
                                         join(scan_rows, "\n") * "\n"))
    universe_sha = bytes2hex(sha256(scan_bytes))
    key_header = "schema_version\tfile\tlobe\tfeature_sha256\tuniverse_sha256\tkeys_sha256\n"
    key_rows = [join(["1", node.file, string(node.lobe), feature_sha,
                      universe_sha, keys_sha], '\t') for node in data.nodes]
    patch_key_bytes = Vector{UInt8}(codeunits(key_header * join(key_rows, "\n") * "\n"))
    universe_artifacts = Dict(
        "scan_universe.tsv" => scan_bytes,
        "patch_keys.tsv" => patch_key_bytes,
        "folds.tsv" => Vector{UInt8}(codeunits("schema_version\n1\n")),
        "seeds.tsv" => Vector{UInt8}(codeunits("schema_version\n1\n")),
        "perturbations.tsv" => Vector{UInt8}(codeunits("schema_version\n1\n")),
        "shards.tsv" => Vector{UInt8}(codeunits("schema_version\n1\n")),
    )
    universe_dir = joinpath(root, "synthetic", "universe")
    mkpath(universe_dir)
    for (name, bytes) in universe_artifacts
        write(joinpath(universe_dir, name), bytes)
    end
    universe_receipt_lines = [
        toml_line("schema", "stmfit-structured-scan-universe-receipt-v1"),
        "schema_version = 1",
        toml_line("feature_sha256", feature_sha),
        toml_line("candidate_config_sha256", SEA.CANDIDATE_CONFIG_SHA256),
        toml_line("model_config_sha256", SEA.MODEL_CONFIG_SHA256),
        toml_line("source_sha256", source_sha),
        toml_line("forward_producer_sha256", SEA.StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256),
        toml_line("backward_producer_sha256", SEA.StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256),
        toml_line("grid_sha256", SEA.StructuredEdgeFeatures.GRID_SHA256),
        toml_line("universe_sha256", universe_sha),
        toml_line("keys_sha256", keys_sha),
        "",
        "[artifacts]",
    ]
    for name in sort!(collect(keys(universe_artifacts)))
        push!(universe_receipt_lines,
              "$(repr(name)) = $(repr(bytes2hex(sha256(universe_artifacts[name]))))")
    end
    universe_receipt = join(universe_receipt_lines, "\n") * "\n"
    write(joinpath(universe_dir, "receipt.toml"), universe_receipt)
    universe_receipt_sha = bytes2hex(sha256(codeunits(universe_receipt)))

    producer_data = joinpath(root, "synthetic", "producer-data")
    mkpath(producer_data)
    forward_patches = write_fixture_file(
        root,
        "synthetic/patches/forward.tsv",
        synthetic_patch_bytes(data, :forward),
    )
    backward_patches = write_fixture_file(
        root,
        "synthetic/patches/backward.tsv",
        synthetic_patch_bytes(data, :backward),
    )
    commands = SEA.StructuredEdgeFeatures.expected_producer_commands(
        root;
        features=original_features,
        data_dir=producer_data,
        forward_patches=forward_patches,
        backward_patches=backward_patches,
    )
    forward_receipt = write_fixture_file(
        root,
        "synthetic/patches/forward-receipt.toml",
        synthetic_patch_receipt_text(
            root, :forward, original_features, producer_data, forward_patches,
            feature_sha, keys_sha, length(data.nodes), commands.forward,
        ),
    )
    backward_receipt = write_fixture_file(
        root,
        "synthetic/patches/backward-receipt.toml",
        synthetic_patch_receipt_text(
            root, :backward, original_features, producer_data, backward_patches,
            feature_sha, keys_sha, length(data.nodes), commands.backward,
        ),
    )

    edge_source_sha = first(SEA.StructuredEdgeFeatures._source_snapshots())
    forward_hash = file_sha256(forward_patches)
    backward_hash = file_sha256(backward_patches)
    node_rows = String[]
    edge_rows = String[]
    by_file = Dict(file => filter(node -> node.file == file, data.nodes)
                   for file in sort!(unique(node.file for node in data.nodes)))
    edge_by_key = Dict((edge.file, edge.left_lobe, edge.right_lobe) => edge
                       for edge in data.edges)
    for file in sort!(collect(keys(by_file)))
        chain = by_file[file]
        segment = SEA.StructuredEdgeFeatures._segment_id(file,
                                                         [node.lobe for node in chain])
        for (index, node) in enumerate(chain)
            left_reason = index == 1 ? "start" : "eligible"
            right_reason = index == length(chain) ? "end" : "eligible"
            push!(node_rows, join([
                file, string(node.lobe), @sprintf("%.17g", node.t_nm), segment,
                "connected", left_reason, right_reason, feature_sha,
                SEA.MODEL_CONFIG_SHA256, edge_source_sha,
            ], '\t'))
        end
        for index in 1:(length(chain) - 1)
            left = chain[index]
            right = chain[index + 1]
            edge = edge_by_key[(file, left.lobe, right.lobe)]
            push!(edge_rows, join([
                file, string(left.lobe), string(right.lobe),
                @sprintf("%.17g", left.t_nm), @sprintf("%.17g", right.t_nm),
                @sprintf("%.17g", right.t_nm - left.t_nm), segment, segment,
                "eligible", "none", @sprintf("%.17g", edge.observation[1]),
                @sprintf("%.17g", edge.observation[2]),
                SEA.StructuredEdgeFeatures.GRID_SHA256, forward_hash,
                backward_hash, feature_sha, SEA.MODEL_CONFIG_SHA256,
                edge_source_sha,
            ], '\t'))
        end
    end
    node_bytes = Vector{UInt8}(codeunits(
        join(SEA.StructuredEdgeFeatures.NODE_HEADER, '\t') * "\n" *
        join(node_rows, "\n") * "\n"))
    edge_bytes = Vector{UInt8}(codeunits(
        join(SEA.StructuredEdgeFeatures.EDGE_HEADER, '\t') * "\n" *
        join(edge_rows, "\n") * "\n"))
    edge_dir = joinpath(root, "synthetic", "edges")
    mkpath(edge_dir)
    write(joinpath(edge_dir, "node_segments.tsv"), node_bytes)
    write(joinpath(edge_dir, "edge_observations.tsv"), edge_bytes)
    receipt_fields = [
        "schema" => "stmfit-structured-edge-feature-receipt-v1",
        "schema_version" => 1,
        "status" => "PASS",
        "edge_observations_file" => "edge_observations.tsv",
        "node_segments_file" => "node_segments.tsv",
        "edge_header" => join(SEA.StructuredEdgeFeatures.EDGE_HEADER, '\t'),
        "node_header" => join(SEA.StructuredEdgeFeatures.NODE_HEADER, '\t'),
        "edge_header_sha256" => bytes2hex(sha256(codeunits(join(SEA.StructuredEdgeFeatures.EDGE_HEADER, '\t') * "\n"))),
        "node_header_sha256" => bytes2hex(sha256(codeunits(join(SEA.StructuredEdgeFeatures.NODE_HEADER, '\t') * "\n"))),
        "edge_observations_sha256" => bytes2hex(sha256(edge_bytes)),
        "node_segments_sha256" => bytes2hex(sha256(node_bytes)),
        "edge_row_count" => length(edge_rows),
        "node_row_count" => length(node_rows),
        "feature_sha256" => feature_sha,
        "candidate_config_sha256" => SEA.CANDIDATE_CONFIG_SHA256,
        "model_config_sha256" => SEA.MODEL_CONFIG_SHA256,
        "universe_receipt_sha256" => universe_receipt_sha,
        "keys_sha256" => keys_sha,
        "grid_sha256" => SEA.StructuredEdgeFeatures.GRID_SHA256,
        "forward_patch_sha256" => forward_hash,
        "backward_patch_sha256" => backward_hash,
        "forward_patch_receipt_sha256" => file_sha256(forward_receipt),
        "backward_patch_receipt_sha256" => file_sha256(backward_receipt),
        "forward_producer_sha256" => SEA.StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256,
        "backward_producer_sha256" => SEA.StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256,
        "forward_command" => commands.forward,
        "backward_command" => commands.backward,
        "source_sha256" => edge_source_sha,
        "t2_review_sha256" => SEA.StructuredEdgeFeatures.T2_REVIEW_SHA256,
        "t3_review_sha256" => SEA.StructuredEdgeFeatures.T3_REVIEW_SHA256,
    ]
    write(joinpath(edge_dir, "receipt.toml"),
          join((toml_line(key, value) for (key, value) in receipt_fields), "\n") * "\n")

    feature_header = copy(SEA.ADMISSION_FEATURE_HEADER)
    forbidden_column && push!(feature_header, "slot_index")
    feature_rows = String[]
    for node in data.nodes
        row = [node.file, string(node.lobe), @sprintf("%.17g", node.t_nm),
               @sprintf("%.17g", node.amplitude)]
        append!(row, @sprintf("%.17g", value) for value in node.predictors)
        forbidden_column && push!(row, string(node.lobe))
        push!(feature_rows, join(row, '\t'))
    end
    features = write_fixture_file(
        root,
        forbidden_column ? "synthetic/admission-forbidden.tsv" :
                           "synthetic/admission.tsv",
        join(feature_header, '\t') * "\n" * join(feature_rows, "\n") * "\n",
    )
    return (
        root=root,
        data=data,
        original_features=original_features,
        producer_data=producer_data,
        forward_patches=forward_patches,
        backward_patches=backward_patches,
        forward_receipt=forward_receipt,
        backward_receipt=backward_receipt,
        commands=commands,
        features=features,
        universe=universe_dir,
        edges=edge_dir,
    )
end

@testset "Todo 10 provenance is validated before shuffle" begin
    fixture = make_cli_fixture()
    try
        expect_cli_pass(fixture)
        @test length(load_cli_fixture(fixture).edges) == 140
        receipt_path = joinpath(fixture.edges, "receipt.toml")
        receipt = TOML.parsefile(receipt_path)
        sha_fields = [
            "edge_header_sha256", "node_header_sha256",
            "edge_observations_sha256", "node_segments_sha256",
            "feature_sha256", "candidate_config_sha256", "model_config_sha256",
            "universe_receipt_sha256", "keys_sha256", "grid_sha256",
            "forward_patch_sha256", "backward_patch_sha256",
            "forward_patch_receipt_sha256", "backward_patch_receipt_sha256",
            "forward_producer_sha256", "backward_producer_sha256",
            "source_sha256", "t2_review_sha256", "t3_review_sha256",
        ]
        for key in sha_fields
            error = expect_pre_shuffle_rejection(fixture) do
                replace_toml_value!(receipt_path, key, "not-a-sha")
            end
            error isa SEA.EdgeAdmissionError && @test error.code == :invalid_hash
        end

        uppercase_error = expect_pre_shuffle_rejection(fixture) do
            replace_toml_value!(receipt_path, "forward_patch_sha256", repeat("A", 64))
        end
        uppercase_error isa SEA.EdgeAdmissionError &&
            @test uppercase_error.code == :invalid_hash

        forward_row_error = expect_pre_shuffle_rejection(fixture) do
            replace_edge_binding!(fixture, 14,
                                  different_sha256(receipt["forward_patch_sha256"]))
        end
        forward_row_error isa SEA.EdgeAdmissionError &&
            @test forward_row_error.code == :edge_hash_mismatch

        backward_row_error = expect_pre_shuffle_rejection(fixture) do
            replace_edge_binding!(fixture, 15,
                                  different_sha256(receipt["backward_patch_sha256"]))
        end
        backward_row_error isa SEA.EdgeAdmissionError &&
            @test backward_row_error.code == :edge_hash_mismatch

        patch_receipt_error = expect_pre_shuffle_rejection(fixture) do
            replace_toml_value!(receipt_path, "forward_patch_receipt_sha256",
                                different_sha256(receipt["forward_patch_receipt_sha256"]))
        end
        patch_receipt_error isa SEA.EdgeAdmissionError &&
            @test patch_receipt_error.code == :patch_receipt_hash_mismatch

        t2_error = expect_pre_shuffle_rejection(fixture) do
            replace_toml_value!(receipt_path, "t2_review_sha256",
                                different_sha256(receipt["t2_review_sha256"]))
        end
        t2_error isa SEA.EdgeAdmissionError &&
            @test t2_error.code == :dependency_hash_mismatch

        t3_error = expect_pre_shuffle_rejection(fixture) do
            replace_toml_value!(receipt_path, "t3_review_sha256",
                                different_sha256(receipt["t3_review_sha256"]))
        end
        t3_error isa SEA.EdgeAdmissionError &&
            @test t3_error.code == :dependency_hash_mismatch

        command_error = expect_pre_shuffle_rejection(fixture) do
            replace_toml_value!(receipt_path, "forward_command",
                                receipt["forward_command"] * " --noncanonical")
        end
        command_error isa SEA.EdgeAdmissionError &&
            @test command_error.code == :producer_command_mismatch

        cli_cases = [
            ("not-a-sha", "invalid_hash",
             () -> replace_toml_value!(receipt_path,
                                       "forward_patch_receipt_sha256",
                                       "not-a-sha")),
            ("uppercase", "invalid_hash",
             () -> replace_toml_value!(receipt_path,
                                       "forward_patch_sha256",
                                       repeat("A", 64))),
            ("forward-row-mismatch", "edge_hash_mismatch",
             () -> replace_edge_binding!(
                 fixture,
                 14,
                 different_sha256(receipt["forward_patch_sha256"]),
             )),
            ("backward-row-mismatch", "edge_hash_mismatch",
             () -> replace_edge_binding!(
                 fixture,
                 15,
                 different_sha256(receipt["backward_patch_sha256"]),
             )),
            ("patch-receipt-mismatch", "patch_receipt_hash_mismatch",
             () -> replace_toml_value!(
                 receipt_path,
                 "forward_patch_receipt_sha256",
                 different_sha256(receipt["forward_patch_receipt_sha256"]),
             )),
            ("t2-review-mismatch", "dependency_hash_mismatch",
             () -> replace_toml_value!(
                 receipt_path,
                 "t2_review_sha256",
                 different_sha256(receipt["t2_review_sha256"]),
             )),
            ("t3-review-mismatch", "dependency_hash_mismatch",
             () -> replace_toml_value!(
                 receipt_path,
                 "t3_review_sha256",
                 different_sha256(receipt["t3_review_sha256"]),
             )),
            ("noncanonical-command", "producer_command_mismatch",
             () -> replace_toml_value!(
                 receipt_path,
                 "forward_command",
                 receipt["forward_command"] * " --noncanonical",
             )),
        ]
        for (name, reason, mutation) in cli_cases
            expect_cli_blocked(mutation, fixture, name, reason)
        end
    finally
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "training-only partition sentinels and informative controls" begin
    data = synthetic_admission_fixture()
    normalized = SEA.normalize_admission_data(data)
    @test length(normalized.audits) == 28
    @test all(audit.node_count == 6 for audit in normalized.audits)
    @test all(audit.file == only(unique(node.file for node in data.nodes
                                       if node.file == audit.file))
              for audit in normalized.audits)
    dates = sort!(unique(node.date for node in data.nodes))
    training_dates = dates[1:5]

    c1 = SEA.fit_partition_model(normalized, "C1", training_dates)
    c2 = SEA.fit_partition_model(normalized, "C2", training_dates)
    @test c1.status == "PASS"
    @test c2.status == "PASS"
    @test c1.support.dates == 5
    @test c1.support.scans == 20
    @test c1.support.edges == 100
    @test all(c1.support.category_ess .>= 10.0)
    @test all(c1.support.minimum_support_scans .>= 3)
    @test c1.residualizer.training_row_count == 100
    @test c1.conditional_fit.means[2, :] == c1.conditional_fit.means[3, :]
    @test length(c1.conditional_fit.starts) == 5
    nonfinite_nodes = copy(data.nodes)
    nonfinite_index = findfirst(node -> node.date == dates[end] && node.lobe == 1,
                                nonfinite_nodes)
    nonfinite_node = nonfinite_nodes[nonfinite_index]
    nonfinite_nodes[nonfinite_index] = SEA.AdmissionNode(
        nonfinite_node.file, nonfinite_node.date, nonfinite_node.lobe,
        nonfinite_node.t_nm, nonfinite_node.amplitude,
        ntuple(index -> index == 1 ? NaN : nonfinite_node.predictors[index], 7),
    )
    nonfinite_data = SEA.AdmissionData(nonfinite_nodes, copy(data.edges), copy(data.hashes))
    nonfinite_fit = SEA.fit_partition_model(
        SEA.normalize_admission_data(nonfinite_data), "C1", training_dates)
    @test length(nonfinite_fit.transform_edges) == length(data.edges)
    @test any(state -> state.status == "ineligible" &&
                        state.reason == :nonfinite_endpoint_predictor,
              nonfinite_fit.transform_edges)
    @test all(state -> state.status == "eligible" ||
                       state.reason == :nonfinite_endpoint_predictor,
              nonfinite_fit.transform_edges)

    changed_nodes = copy(data.nodes)
    heldout_index = findfirst(node -> node.date == dates[end] && node.lobe == 1,
                              changed_nodes)
    changed = changed_nodes[heldout_index]
    changed_nodes[heldout_index] = SEA.AdmissionNode(
        changed.file,
        changed.date,
        changed.lobe,
        changed.t_nm,
        changed.amplitude,
        ntuple(index -> changed.predictors[index] + 50.0 * index, 7),
    )
    changed_edges = copy(data.edges)
    heldout_edge_index = findfirst(edge -> edge.date == dates[end], changed_edges)
    changed_edge = changed_edges[heldout_edge_index]
    changed_edges[heldout_edge_index] = SEA.AdmissionEdge(
        changed_edge.file,
        changed_edge.date,
        changed_edge.left_lobe,
        changed_edge.right_lobe,
        changed_edge.left_t_nm,
        changed_edge.right_t_nm,
        changed_edge.segment_id,
        changed_edge.ordinal,
        (changed_edge.observation[1] + 0.4,
         changed_edge.observation[2] - 0.3),
    )
    mutated = SEA.AdmissionData(changed_nodes, changed_edges, copy(data.hashes))
    mutated_normalized = SEA.normalize_admission_data(mutated)
    mutated_fit = SEA.fit_partition_model(mutated_normalized, "C1", training_dates)
    @test mutated_fit.status == "PASS"
    @test mutated_fit.training_sha256 == c1.training_sha256
    @test mutated_fit.fit_sha256 == c1.fit_sha256
    training_files = Set(sample.edge.file for sample in c1.samples
                         if sample.edge.date in training_dates)
    @test training_files == Set(c1.unary.training_scans)

    untouched = Dict(SEA._edge_identity(sample.edge) => copy(sample.residual)
                     for sample in c1.samples if sample.edge.date == dates[6])
    evaluation = SEA.evaluate_partition(c1, dates[6]; run_shuffle=true)
    @test evaluation.status == "PASS"
    @test length(evaluation.shuffle_values) == 500
    @test length(evaluation.shuffle_fit_sha256) == 500
    @test evaluation.date_mean > 0.0
    @test evaluation.control_date_mean > evaluation.shuffle_quantile
    @test evaluation.shuffle_pass
    @test all(untouched[SEA._edge_identity(sample.edge)] == sample.residual
              for sample in c1.samples if sample.edge.date == dates[6])

    reversed_fit = SEA.fit_partition_model(normalized, "C1", training_dates;
                                           reversed=true)
    reversed_evaluation = SEA.evaluate_partition(reversed_fit, dates[6];
                                                  run_shuffle=false)
    checked = SEA._with_reversal(evaluation, reversed_evaluation)
    original_residuals = Dict(SEA._edge_identity(sample.edge) => sample.residual
                              for sample in c1.samples if sample.edge.date == dates[6])
    reversed_residuals = Dict(SEA._edge_identity(sample.edge) => sample.residual
                              for sample in reversed_fit.samples if sample.edge.date == dates[6])
    reversal_residual_max_error = maximum(norm(original_residuals[key] -
                                               reversed_residuals[key], Inf)
                                          for key in keys(original_residuals))
    @test reversal_residual_max_error <= 1.0e-12
    reversal_keys = sort!(collect(keys(evaluation.gains)))
    reversal_max_error = maximum(abs(evaluation.gains[key] -
                                     reversed_evaluation.gains[key])
                                 for key in reversal_keys)
    @test reversal_max_error <= 1.0e-12
    @test checked.status == "PASS"
    @test checked.reversal_pass
    sample = first(c1.samples)
    reversed_sample = only(filter(candidate ->
        SEA._edge_identity(candidate.edge) == SEA._edge_identity(sample.edge),
        reversed_fit.samples))
    canonical_prediction = SEA._residualizer_prediction(
        c1.residualizer, sample.endpoint_predictors)
    reversed_prediction = SEA._residualizer_prediction(
        reversed_fit.residualizer, reversed_sample.endpoint_predictors)
    @test canonical_prediction ≈ reversed_prediction atol=2e-15 rtol=2e-15
end

@testset "fully nested outer and final admission partitions" begin
    data = synthetic_admission_fixture()
    normalized = SEA.normalize_admission_data(data)
    report = SEA.evaluate_admission(data)
    COMPLETE_REPLAY_SHA256[] = report.result_sha256
    @test report.status == "PASS"
    result = report.models["C1"]
    robust_result = report.models["C2"]
    @test result.status == "PASS"
    @test result.reason == :ok
    @test length(result.outer_folds) == 7
    @test all(fold.inner_gate.status == "PASS" for fold in result.outer_folds)
    @test all(length(fold.inner_gate.partitions) == 6 for fold in result.outer_folds)
    @test all(fold.outer_score.status == "PASS" for fold in result.outer_folds)
    @test all(fold.outer_score.reversal_pass for fold in result.outer_folds)
    @test result.final_gate.status == "PASS"
    @test length(result.final_gate.partitions) == 7
    @test all(length(partition.shuffle_values) == 500
              for fold in result.outer_folds
              for partition in fold.inner_gate.partitions)
    @test all(length(partition.shuffle_values) == 500
              for partition in result.final_gate.partitions)
    @test result.full_fit.status == "PASS"
    @test length(result.result_sha256) == 64
    @test robust_result.status == "PASS"
    @test length(robust_result.outer_folds) == 7
    @test robust_result.final_gate.status == "PASS"
    @test robust_result.full_fit.status == "PASS"
    @test !isdefined(SEA, :chain_inference)
    @test !isdefined(SEA, :infer_chain)
    @test !any(occursin("chain_inference", String(name))
               for name in names(SEA; all=true))
    @test begin
        t11_source = read(EDGE_ADMISSION_SOURCE, String)
        !occursin("edge_model.jl", t11_source) &&
        !occursin("chain_inference.jl", t11_source)
    end

    six_dates = sort!(unique(node.date for node in data.nodes))[1:6]
    six_data = SEA.AdmissionData(
        [node for node in data.nodes if node.date in six_dates],
        [edge for edge in data.edges if edge.date in six_dates],
        copy(data.hashes),
    )
    six_result = SEA.evaluate_model_admission(
        SEA.normalize_admission_data(six_data), "C1")
    @test six_result.status == "SKIPPED"
    @test six_result.reason == :insufficient_dates
    @test length(six_result.full_fit.transform_edges) == length(six_data.edges)
    @test all(state -> state.status == "unavailable" &&
                       state.reason == :residualizer_unavailable,
              six_result.full_fit.transform_edges)
    six_report = SEA.AdmissionReport(
        six_result.status,
        six_result.reason,
        Dict("C1" => six_result),
        copy(six_data.hashes),
        six_result.result_sha256,
    )
    six_files = SEA.report_files(six_report)
    _, six_edge_rows = SEA.StructuredEdgeFeatures._parse_table(
        six_files["fitted_edge_transforms.tsv"], "six-date fitted edge transforms")
    @test length(six_edge_rows) == length(six_data.edges)
    @test all(row -> row[21] == "unavailable" && row[22] == "residualizer_unavailable" &&
                       row[19] == "NA" && row[20] == "NA", six_edge_rows)
    @test Set(row[16] for row in six_edge_rows) ==
          Set(edge.raw_row_sha256 for edge in six_data.edges)
    six_publish_root = mktempdir(; prefix="stmfit-edge-admission-six-date-")
    try
        SEA.publish_report(six_publish_root, "six-date", six_report)
        six_published = read(joinpath(
            six_publish_root, "six-date", "fitted_edge_transforms.tsv"))
        _, six_published_rows = SEA.StructuredEdgeFeatures._parse_table(
            six_published, "six-date published fitted edge transforms")
        @test length(six_published_rows) == length(six_data.edges)
        @test all(row -> row[21] == "unavailable" &&
                           row[22] == "residualizer_unavailable" &&
                           row[19] == "NA" && row[20] == "NA",
                  six_published_rows)
    finally
        rm(six_publish_root; recursive=true, force=true)
    end

    nonfinite_six_nodes = copy(six_data.nodes)
    nonfinite_index = findfirst(node -> node.date == six_dates[end] && node.lobe == 1,
                                nonfinite_six_nodes)
    nonfinite_six_node = nonfinite_six_nodes[nonfinite_index]
    nonfinite_six_nodes[nonfinite_index] = SEA.AdmissionNode(
        nonfinite_six_node.file, nonfinite_six_node.date, nonfinite_six_node.lobe,
        nonfinite_six_node.t_nm, nonfinite_six_node.amplitude,
        ntuple(index -> index == 1 ? NaN : nonfinite_six_node.predictors[index], 7),
    )
    nonfinite_six = SEA.evaluate_model_admission(
        SEA.normalize_admission_data(SEA.AdmissionData(
            nonfinite_six_nodes, copy(six_data.edges), copy(six_data.hashes))), "C1")
    @test any(state -> state.status == "ineligible" &&
                        state.reason == :nonfinite_endpoint_predictor,
              nonfinite_six.full_fit.transform_edges)

    handoff_files = SEA.report_files(report)
    repeated_handoff_files = SEA.report_files(report)
    @test handoff_files == repeated_handoff_files
    GRAPH_HANDOFF_REPLAY_SHA256[] = ""
    @test sort!(collect(keys(handoff_files))) == [
        "bootstrap.tsv", "ess.tsv", "fitted_edge_models.tsv",
        "fitted_edge_transforms.tsv", "fitted_unary_nodes.tsv", "models.tsv",
        "partitions.tsv", "receipt.toml", "scores.tsv", "shuffle.tsv", "starts.tsv",
    ]
    let admission_receipt = TOML.parse(String(handoff_files["receipt.toml"]))
        GRAPH_HANDOFF_REPLAY_SHA256[] = admission_receipt["graph_handoff_replay_sha256"]
        @test admission_receipt["status"] == "PASS"
        @test admission_receipt["schema"] == "stmfit-structured-edge-admission-receipt-v2"
        @test admission_receipt["schema_version"] == 2
        @test admission_receipt["graph_handoff_schema"] == SEA.GRAPH_HANDOFF_SCHEMA
        @test admission_receipt["category_order"] == collect(SEA.GRAPH_CATEGORY_ORDER)
        @test admission_receipt["state_order"] == collect(SEA.GRAPH_STATE_ORDER)
        @test admission_receipt["bootstrap_count"] == 500
        @test admission_receipt["shuffle_count"] == 500
        @test admission_receipt["fixed_nu"] == 8
        @test admission_receipt["node_priors"] == [0.5, 0.5]
        @test admission_receipt["pair_prior"] == [0.0, 0.0, 0.0, 0.0]
        @test Set(keys(admission_receipt["artifacts"])) ==
              Set(setdiff(collect(keys(handoff_files)), ["receipt.toml"]))
    end

    model_header, model_rows = SEA.StructuredEdgeFeatures._parse_table(
        handoff_files["fitted_edge_models.tsv"], "synthetic fitted edge models")
    unary_header, unary_rows = SEA.StructuredEdgeFeatures._parse_table(
        handoff_files["fitted_unary_nodes.tsv"], "synthetic fitted unary nodes")
    edge_header, edge_rows = SEA.StructuredEdgeFeatures._parse_table(
        handoff_files["fitted_edge_transforms.tsv"], "synthetic fitted edge transforms")
    @test model_header == SEA.GRAPH_MODEL_HEADER
    @test unary_header == SEA.GRAPH_UNARY_HEADER
    @test edge_header == SEA.GRAPH_EDGE_HEADER
    references = SEA._graph_fit_references(report)
    fit_hashes = [reference.fit.fit_sha256 for reference in references]
    @test length(model_rows) == length(unique(fit_hashes))
    partition_header, partition_rows = SEA.StructuredEdgeFeatures._parse_table(
        handoff_files["partitions.tsv"], "synthetic partitions")
    partition_fit_hashes = [row[12] for row in partition_rows]
    model_fit_hashes = [row[3] for row in model_rows]
    model_columns = Dict(name => index for (index, name) in enumerate(model_header))
    unary_columns = Dict(name => index for (index, name) in enumerate(unary_header))
    edge_columns = Dict(name => index for (index, name) in enumerate(edge_header))
    for row in model_rows
        assert_certified_published_scale(row, model_columns, "null")
        assert_certified_published_scale(row, model_columns, "conditional")
    end
    @test all(row[model_columns["status"]] == "PASS" for row in model_rows)
    @test all(count(==(fit_hash), model_fit_hashes) == 1
              for fit_hash in unique(vcat(partition_fit_hashes,
                                          [result.full_fit.fit_sha256
                                           for result in values(report.models)])))
    for row in model_rows
        reference_list = split(row[model_columns["references"]], ';')
        @test length(reference_list) == parse(Int, row[model_columns["reference_count"]])
        @test SEA._hash_lines([join(split(reference, '|'), '\t')
                               for reference in reference_list]) ==
              row[model_columns["references_sha256"]]
    end
    expected_reference_keys = Set{String}()
    for model in sort!(collect(keys(report.models)))
        result_model = report.models[model]
        for item in SEA._logical_evaluations(result_model)
            evaluation = item.evaluation
            push!(expected_reference_keys, join((
                "partition", model, item.scope, item.outer_date,
                evaluation.target_date, evaluation.fit.fit_sha256,
            ), '|'))
        end
        push!(expected_reference_keys, join((
            "full_refit", model, "full_refit", "NA", "NA",
            result_model.full_fit.fit_sha256,
        ), '|'))
    end
    observed_reference_keys = Set{String}()
    for row in model_rows
        union!(observed_reference_keys,
               split(row[model_columns["references"]], ';'))
    end
    @test observed_reference_keys == expected_reference_keys
    @test length(unary_rows) ==
          sum(length(reference.fit.unary.node_identities) for reference in references)
    @test length(edge_rows) == sum(length(reference.fit.transform_edges)
                                   for reference in references)
    @test GRAPH_HANDOFF_REPLAY_SHA256[] == SEA.graph_handoff_replay_hash(handoff_files)

    first_edge = first(edge_rows)
    first_fit = only(filter(reference -> reference.fit.fit_sha256 ==
                             first_edge[edge_columns["fit_sha256"]], references))
    first_sample = only(filter(sample -> SEA._graph_edge_identity(sample.edge) ==
                               first_edge[edge_columns["edge_identity"]],
                               first_fit.fit.samples))
    raw = collect(first_sample.edge.observation)
    prediction = [parse(Float64, first_edge[edge_columns["pred_fwd"]]),
                  parse(Float64, first_edge[edge_columns["pred_bwd"]])]
    @test raw - prediction ≈ first_sample.residual atol=2e-15 rtol=2e-15
    @test first_edge[edge_columns["edge_row_raw_sha256"]] ==
          first_sample.edge.raw_row_sha256
    @test first_sample.edge.raw_row_sha256 ==
          independent_synthetic_edge_row_hash(first_sample.edge)
    @test !any(occursin(name, join(edge_header, '\t'))
               for name in ("raw_fwd", "raw_bwd", "residual_fwd", "residual_bwd"))
    @test first_edge[edge_columns["residualizer_sha256"]] ==
          SEA._residualizer_sha256(first_fit.fit.residualizer)

    unary_by_key = Dict(
        (row[unary_columns["fit_sha256"]], row[unary_columns["file"]],
         parse(Int, row[unary_columns["lobe"]])) => row
        for row in unary_rows
    )
    left_unary = unary_by_key[(first_fit.fit.fit_sha256,
                               first_sample.edge.file, first_sample.edge.left_lobe)]
    right_unary = unary_by_key[(first_fit.fit.fit_sha256,
                                first_sample.edge.file, first_sample.edge.right_lobe)]
    left_probability = parse(Float64, left_unary[unary_columns["p1"]])
    right_probability = parse(Float64, right_unary[unary_columns["p1"]])
    parsed_model = only(filter(row -> row[model_columns["fit_sha256"]] ==
                               first_fit.fit.fit_sha256, model_rows))
    parsed_means = [
        [parse(Float64, parsed_model[model_columns["conditional_$(category)_mean_$(output)"]])
         for output in SEA.GRAPH_OUTPUT_NAMES]
        for category in SEA.GRAPH_CATEGORY_ORDER
    ]
    parsed_scale = reshape(Float64[
        parse(Float64, parsed_model[model_columns["conditional_scale_ff"]]),
        parse(Float64, parsed_model[model_columns["conditional_scale_fb"]]),
        parse(Float64, parsed_model[model_columns["conditional_scale_fb"]]),
        parse(Float64, parsed_model[model_columns["conditional_scale_bb"]]),
    ], 2, 2)
    parsed_null_mean = [
        parse(Float64, parsed_model[model_columns["null_mean_fwd"]]),
        parse(Float64, parsed_model[model_columns["null_mean_bwd"]]),
    ]
    parsed_null_scale = reshape(Float64[
        parse(Float64, parsed_model[model_columns["null_scale_ff"]]),
        parse(Float64, parsed_model[model_columns["null_scale_fb"]]),
        parse(Float64, parsed_model[model_columns["null_scale_fb"]]),
        parse(Float64, parsed_model[model_columns["null_scale_bb"]]),
    ], 2, 2)
    parsed_residual = raw - prediction
    parsed_gain, parsed_ratios = independent_edge_gain(
        parsed_residual,
        collect(SEA.pair_weights(left_probability, right_probability)),
        reduce(vcat, (reshape(mean, 1, :) for mean in parsed_means)),
        parsed_scale,
        parsed_null_mean,
        parsed_null_scale,
    )
    @test length(parsed_ratios) == 4
    @test SEA.edge_gain(
        parsed_residual,
        collect(SEA.pair_weights(left_probability, right_probability)),
        reduce(vcat, (reshape(mean, 1, :) for mean in parsed_means)),
        parsed_scale,
        parsed_null_mean,
        parsed_null_scale,
    ) ≈ parsed_gain atol=1e-14 rtol=1e-14
    @test parsed_gain ≈ SEA._sample_gain(first_sample, first_fit.fit) atol=1e-14 rtol=1e-14
    @test parsed_model[model_columns["conditional_01_mean_corr_fwd"]] ==
          parsed_model[model_columns["conditional_10_mean_corr_fwd"]]
    @test parsed_model[model_columns["conditional_01_mean_corr_bwd"]] ==
          parsed_model[model_columns["conditional_10_mean_corr_bwd"]]
    @test parsed_model[model_columns["category_order"]] == "00,01,10,11"
    @test parsed_model[model_columns["state_order"]] == "0,1"
    zero_means = repeat(parsed_null_mean', 4, 1)
    zero_gain, zero_ratios = independent_edge_gain(
        parsed_residual,
        collect(SEA.pair_weights(left_probability, right_probability)),
        zero_means,
        parsed_null_scale,
        parsed_null_mean,
        parsed_null_scale,
    )
    @test zero_ratios == zeros(4)
    @test zero_gain == 0.0
    @test SEA.edge_gain(
        parsed_residual,
        collect(SEA.pair_weights(left_probability, right_probability)),
        zero_means,
        parsed_null_scale,
        parsed_null_mean,
        parsed_null_scale,
    ) == 0.0
    joint_zero, marginal_zero = independent_two_node_joint(0.2, 0.7, zeros(2, 2))
    @test marginal_zero[1] ≈ [0.8, 0.2] atol=1e-12 rtol=0
    @test marginal_zero[2] ≈ [0.3, 0.7] atol=1e-12 rtol=0
    @test argmax(marginal_zero[1]) == argmax([0.8, 0.2])
    @test argmax(marginal_zero[2]) == argmax([0.3, 0.7])
    @test joint_zero ≈ [0.24 0.56; 0.06 0.14] atol=1e-15 rtol=0

    persisted_fixture = make_cli_fixture()
    try
        persisted_data = load_cli_fixture(persisted_fixture)
        persisted_edge = only(filter(edge ->
            edge.file == first_sample.edge.file &&
            edge.left_lobe == first_sample.edge.left_lobe &&
            edge.right_lobe == first_sample.edge.right_lobe,
            persisted_data.edges))
        persisted_lines = split(chomp(read(joinpath(
            persisted_fixture.edges, "edge_observations.tsv"), String)), '\n')[2:end]
        persisted_line = only(filter(line -> begin
            fields = split(line, '\t'; keepempty=true)
            fields[1] == first_sample.edge.file &&
            parse(Int, fields[2]) == first_sample.edge.left_lobe &&
            parse(Int, fields[3]) == first_sample.edge.right_lobe
        end, persisted_lines))
        persisted_fields = split(persisted_line, '\t'; keepempty=true)
        persisted_raw = [parse(Float64, persisted_fields[11]),
                         parse(Float64, persisted_fields[12])]
        @test bytes2hex(sha256(codeunits(persisted_line * "\n"))) ==
              persisted_edge.raw_row_sha256
        @test persisted_edge.observation == Tuple(persisted_raw)
        @test persisted_raw - prediction ≈ first_sample.residual atol=1e-14 rtol=1e-14
        persisted_ratios = independent_edge_log_ratios(
            persisted_raw - prediction,
            reduce(vcat, (reshape(mean, 1, :) for mean in parsed_means)),
            parsed_scale,
            parsed_null_mean,
            parsed_null_scale,
        )
        persisted_matrix = [persisted_ratios[1] persisted_ratios[2];
                            persisted_ratios[3] persisted_ratios[4]]
        in_memory_ratios = independent_edge_log_ratios(
            persisted_raw - prediction,
            first_fit.fit.conditional_fit.means,
            first_fit.fit.conditional_fit.scale,
            first_fit.fit.null_fit.mean,
            first_fit.fit.null_fit.scale,
        )
        in_memory_matrix = [in_memory_ratios[1] in_memory_ratios[2];
                            in_memory_ratios[3] in_memory_ratios[4]]
        @test persisted_matrix ≈ in_memory_matrix atol=1e-14 rtol=1e-14

        reversed_fit = SEA.fit_partition_model(
            normalized, first_fit.fit.model_id, first_fit.fit.training_dates;
            reversed=true)
        reversed_sample = only(filter(sample ->
            SEA._graph_edge_identity(sample.edge) ==
            SEA._graph_edge_identity(first_sample.edge), reversed_fit.samples))
        reversed_prediction = SEA._residualizer_prediction(
            reversed_fit.residualizer, reversed_sample.endpoint_predictors)
        reversed_ratios = independent_edge_log_ratios(
            persisted_raw - reversed_prediction,
            reversed_fit.conditional_fit.means,
            reversed_fit.conditional_fit.scale,
            reversed_fit.null_fit.mean,
            reversed_fit.null_fit.scale,
        )
        reversed_matrix = [reversed_ratios[1] reversed_ratios[2];
                           reversed_ratios[3] reversed_ratios[4]]
        @test persisted_raw - reversed_prediction ≈ reversed_sample.residual atol=1e-14 rtol=1e-14
        @test reversed_matrix ≈ persisted_matrix' atol=1e-12 rtol=1e-12
    finally
        rm(persisted_fixture.root; recursive=true, force=true)
    end

    isolated_nodes = copy(data.nodes)
    push!(isolated_nodes, SEA.AdmissionNode(
        "20260101_singleton.sxm", "20260101", 99, 99.0, 1.0,
        ntuple(index -> 0.2 * index, 7),
    ))
    sort!(isolated_nodes; by=node -> (node.file, node.t_nm, node.lobe))
    isolated_data = SEA.AdmissionData(isolated_nodes, copy(data.edges), copy(data.hashes))
    isolated_normalized = SEA.normalize_admission_data(isolated_data)
    isolated_fit = SEA.fit_partition_model(
        isolated_normalized,
        "C1",
        sort!(unique(node.date for node in isolated_nodes)),
    )
    @test isolated_fit.status == "PASS"
    @test ("20260101_singleton.sxm", "20260101", 99, 99.0) in
          isolated_fit.unary.node_identities
    isolated_result = SEA.ModelAdmissionResult(
        result.status,
        result.reason,
        result.model_id,
        result.outer_folds,
        result.final_gate,
        isolated_fit,
        result.result_sha256,
    )
    isolated_report = SEA.AdmissionReport(
        report.status,
        report.reason,
        Dict("C1" => isolated_result, "C2" => robust_result),
        copy(report.hashes),
        report.result_sha256,
    )
    isolated_files = SEA.report_files(isolated_report)
    _, isolated_unary_rows = SEA.StructuredEdgeFeatures._parse_table(
        isolated_files["fitted_unary_nodes.tsv"], "isolated fitted unary nodes")
    @test any(row -> row[4] == "20260101_singleton.sxm" && row[6] == "99",
              isolated_unary_rows)

    publication_error = function(mutation::Function)
        mutated = copy(handoff_files)
        mutation(mutated)
        fit_before = SEA._FIT_CALL_COUNT[]
        before = SEA._SHUFFLE_CALL_COUNT[]
        error = try
            SEA._validate_graph_handoff_artifacts(report, mutated)
            nothing
        catch caught
            caught
        end
        @test SEA._SHUFFLE_CALL_COUNT[] == before
        @test SEA._FIT_CALL_COUNT[] == fit_before
        @test error isa SEA.EdgeAdmissionError
        return error
    end
    @test publication_error(files -> delete!(files, "fitted_edge_models.tsv")).code ==
          :graph_handoff_schema
    @test publication_error(files -> files["extra.tsv"] = Vector{UInt8}(codeunits("x\n"))).code ==
          :graph_handoff_schema
    @test publication_error(files -> files["fitted_unary_nodes.tsv"] =
                            replace(files["fitted_unary_nodes.tsv"],
                                    UInt8('\t') => UInt8('\n'); count=1)).code ==
          :graph_handoff_schema
    tampered_receipt = replace(handoff_files["receipt.toml"], UInt8('2') => UInt8('3'); count=1)
    receipt_error = try
        SEA._validate_receipt_bytes(report,
                                    Dict(name => bytes for (name, bytes) in handoff_files
                                         if name != "receipt.toml"),
                                    tampered_receipt)
        nothing
    catch caught
        caught
    end
    @test receipt_error isa SEA.EdgeAdmissionError
    @test receipt_error.code in (:graph_handoff_schema, :schema_mismatch)

    mutate_tsv = function(bytes::Vector{UInt8}, field::String, value::String)
        lines = split(chomp(String(copy(bytes))), '\n')
        header = split(lines[1], '\t'; keepempty=true)
        index = only(findall(==(field), header))
        fields = split(lines[2], '\t'; keepempty=true)
        fields[index] = value
        lines[2] = join(fields, '\t')
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    mutate_header = function(bytes::Vector{UInt8})
        lines = split(chomp(String(copy(bytes))), '\n')
        fields = split(lines[1], '\t'; keepempty=true)
        lines[1] = join(vcat(fields[2:end], fields[1]), '\t')
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    append_row = function(bytes::Vector{UInt8})
        text = chomp(String(copy(bytes)))
        lines = split(text, '\n')
        push!(lines, lines[end])
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    remove_row = function(bytes::Vector{UInt8})
        lines = split(chomp(String(copy(bytes))), '\n')
        deleteat!(lines, 2)
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    append_field = function(bytes::Vector{UInt8})
        lines = split(chomp(String(copy(bytes))), '\n')
        lines[1] *= "\textra_field"
        lines[2] *= "\textra_value"
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    remove_field = function(bytes::Vector{UInt8})
        lines = split(chomp(String(copy(bytes))), '\n')
        for row in eachindex(lines)
            fields = split(lines[row], '\t'; keepempty=true)
            deleteat!(fields, length(fields))
            lines[row] = join(fields, '\t')
        end
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    for file in SEA.GRAPH_HANDOFF_FILES
        @test publication_error(files -> delete!(files, file)).code == :graph_handoff_schema
        @test publication_error(files -> files[file] = append_row(files[file])).code ==
              :graph_handoff_schema
        @test publication_error(files -> files[file] = remove_row(files[file])).code ==
              :graph_handoff_schema
        @test publication_error(files -> files[file] = append_field(files[file])).code ==
              :graph_handoff_schema
        @test publication_error(files -> files[file] = remove_field(files[file])).code ==
              :graph_handoff_schema
        @test publication_error(files -> files[file] = mutate_header(files[file])).code ==
              :graph_handoff_schema
    end
    for (file, field, value) in (
        ("fitted_edge_models.tsv", "reference_count", "999"),
        ("fitted_edge_models.tsv", "references_sha256", repeat("0", 64)),
        ("fitted_edge_models.tsv", "references", "tampered"),
        ("fitted_unary_nodes.tsv", "state_order", "1,0"),
        ("fitted_unary_nodes.tsv", "p0", "0.25"),
        ("fitted_unary_nodes.tsv", "node_order_sha256", repeat("0", 64)),
        ("fitted_edge_transforms.tsv", "edge_row_raw_sha256", repeat("0", 64)),
        ("fitted_edge_transforms.tsv", "pred_fwd", "0.25"),
        ("fitted_edge_transforms.tsv", "status", "tampered"),
        ("fitted_edge_transforms.tsv", "reason", "tampered"),
    )
        @test publication_error(files -> files[file] =
                                 mutate_tsv(files[file], field, value)).code ==
              :graph_handoff_schema
    end
    @test publication_error(files -> files["fitted_edge_models.tsv"] =
                            mutate_header(files["fitted_edge_models.tsv"])).code ==
          :graph_handoff_schema

    receipt_without_self = Dict(name => bytes for (name, bytes) in handoff_files
                                if name != "receipt.toml")
    mutate_receipt_line = function(bytes::Vector{UInt8}, key::String, value::String)
        lines = split(chomp(String(copy(bytes))), '\n')
        indices = findall(line -> startswith(line, key * " = "), lines)
        isempty(indices) && error("receipt field missing: $key")
        index = only(indices)
        lines[index] = key * " = " * value
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    for (key, value) in (
        ("graph_handoff_replay_sha256", repr(repeat("0", 64))),
        ("category_order", "[\"00\", \"10\", \"01\", \"11\"]"),
        ("state_order", "[\"1\", \"0\"]"),
        ("coefficient_order", "[\"tampered\"]"),
    )
        bad_receipt = mutate_receipt_line(
            SEA.report_files(report)["receipt.toml"], key, value)
        @test begin
            fit_before = SEA._FIT_CALL_COUNT[]
            before = SEA._SHUFFLE_CALL_COUNT[]
            error = try
                SEA._validate_receipt_bytes(report, receipt_without_self, bad_receipt)
                nothing
            catch caught
                caught
            end
            error isa SEA.EdgeAdmissionError &&
            SEA._SHUFFLE_CALL_COUNT[] == before &&
            SEA._FIT_CALL_COUNT[] == fit_before
        end
    end
    remove_receipt_line = function(bytes::Vector{UInt8}, key::String)
        lines = split(chomp(String(copy(bytes))), '\n')
        index = only(findall(line -> startswith(line, key * " = "), lines))
        deleteat!(lines, index)
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    extra_receipt_key = function(bytes::Vector{UInt8})
        text = String(copy(bytes))
        Vector{UInt8}(codeunits("extra_required = \"tampered\"\n" * text))
    end
    reorder_receipt = function(bytes::Vector{UInt8})
        lines = split(chomp(String(copy(bytes))), '\n')
        lines[1], lines[2] = lines[2], lines[1]
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    remove_receipt_section = function(bytes::Vector{UInt8}, section::String)
        text = String(copy(bytes))
        marker = "\n[" * section * "]\n"
        occursin(marker, text) || error("receipt section missing: $section")
        Vector{UInt8}(codeunits(replace(text, marker => "\n"; count=1)))
    end
    reorder_receipt_sections = function(bytes::Vector{UInt8})
        lines = split(chomp(String(copy(bytes))), '\n')
        hashes = only(findall(==("[hashes]"), lines))
        models = only(findall(==("[models]"), lines))
        lines[hashes], lines[models] = lines[models], lines[hashes]
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    mutate_artifact_receipt = function(bytes::Vector{UInt8})
        lines = split(chomp(String(copy(bytes))), '\n')
        index = first(findall(line -> startswith(line, "\"fitted_edge_models.tsv\" = "),
                             lines))
        lines[index] = "\"fitted_edge_models.tsv\" = \"$(repeat("0", 64))\""
        Vector{UInt8}(codeunits(join(lines, '\n') * "\n"))
    end
    receipt_mutations = Vector{Vector{UInt8}}([
        remove_receipt_line(SEA.report_files(report)["receipt.toml"], "schema"),
        remove_receipt_section(SEA.report_files(report)["receipt.toml"], "graph_handoff"),
        extra_receipt_key(SEA.report_files(report)["receipt.toml"]),
        reorder_receipt(SEA.report_files(report)["receipt.toml"]),
        reorder_receipt_sections(SEA.report_files(report)["receipt.toml"]),
        mutate_receipt_line(SEA.report_files(report)["receipt.toml"],
                            "graph_handoff_replay_sha256", repr(repeat("0", 64))),
        mutate_receipt_line(SEA.report_files(report)["receipt.toml"],
                            "edge_transform_row_count", "0"),
        mutate_receipt_line(SEA.report_files(report)["receipt.toml"],
                            "fit_hashes", "[\"tampered\"]"),
        mutate_artifact_receipt(SEA.report_files(report)["receipt.toml"]),
    ])
    for bad_receipt in receipt_mutations
        @test begin
            fit_before = SEA._FIT_CALL_COUNT[]
            before = SEA._SHUFFLE_CALL_COUNT[]
            error = try
                SEA._validate_receipt_bytes(report, receipt_without_self, bad_receipt)
                nothing
            catch caught
                caught
            end
            error isa SEA.EdgeAdmissionError &&
            SEA._SHUFFLE_CALL_COUNT[] == before &&
            SEA._FIT_CALL_COUNT[] == fit_before
        end
    end

    publication_root = mktempdir(; prefix="stmfit-edge-admission-publication-")
    try
        first_hashes = SEA.publish_report(publication_root, "output", report)
        second_hashes = SEA.publish_report(publication_root, "output", report)
        @test first_hashes == second_hashes
        @test SEA.graph_handoff_replay_hash(handoff_files) ==
              SEA.graph_handoff_replay_hash(handoff_files)
        file_victim = joinpath(publication_root, "file-victim")
        write(file_victim, "victim\n")
        file_symlink_output = joinpath(publication_root, "file-symlink-output")
        SEA.publish_report(publication_root, "file-symlink-output", report)
        rm(joinpath(file_symlink_output, "fitted_edge_models.tsv"))
        symlink(file_victim, joinpath(file_symlink_output, "fitted_edge_models.tsv"))
        file_symlink_error = try
            SEA.publish_report(publication_root, "file-symlink-output", report)
            nothing
        catch caught
            caught
        end
        @test file_symlink_error isa SEA.EdgeAdmissionError
        @test file_symlink_error.code == :publication_collision
        write(joinpath(publication_root, "output", "fitted_edge_models.tsv"), "changed\n")
        collision = try
            SEA.publish_report(publication_root, "output", report)
            nothing
        catch caught
            caught
        end
        @test collision isa SEA.EdgeAdmissionError
        @test collision.code == :publication_collision

        victim = joinpath(publication_root, "victim")
        write(victim, "victim\n")
        linked = joinpath(publication_root, "linked-output")
        symlink(victim, linked)
        symlink_error = try
            SEA.publish_report(publication_root, "linked-output", report)
            nothing
        catch caught
            caught
        end
        @test hasproperty(symlink_error, :code)
        @test symlink_error.code == :symlink_rejected

        directory_victim = joinpath(publication_root, "directory-victim")
        mkpath(directory_victim)
        directory_sentinel = joinpath(directory_victim, "sentinel")
        write(directory_sentinel, "untouched\n")
        linked_directory = joinpath(publication_root, "linked-directory-output")
        symlink(directory_victim, linked_directory)
        directory_symlink_error = try
            SEA.publish_report(publication_root, "linked-directory-output", report)
            nothing
        catch caught
            caught
        end
        @test hasproperty(directory_symlink_error, :code)
        @test directory_symlink_error.code == :symlink_rejected
        @test read(directory_sentinel, String) == "untouched\n"

        failed_output = joinpath(publication_root, "failed-output")
        SEA.StructuredEdgeFeatures._PUBLICATION_INSTALL_HOOK[] =
            (stage, destination) -> error("synthetic publication failure")
        failure = try
            SEA.publish_report(publication_root, "failed-output", report)
            nothing
        catch caught
            caught
        finally
            SEA.StructuredEdgeFeatures._PUBLICATION_INSTALL_HOOK[] = nothing
        end
        @test failure !== nothing
        @test !ispath(failed_output)
        @test isempty(filter(name -> startswith(name, ".structured-admission-stage-"),
                             readdir(publication_root)))
    finally
        rm(publication_root; recursive=true, force=true)
    end

    for bytes in values(handoff_files)
        text = lowercase(String(copy(bytes)))
        for token in ("benchmark", "truth", "grader", "expected_n", "nknnkn",
                      "composition", "top_k", "class_count", "run_length",
                      "transition", "terminal_position")
            @test !occursin(token, text)
        end
    end
end

@testset "non-PASS scale publication is certified or NA" begin
    data = synthetic_admission_fixture()
    report = SEA.evaluate_admission(data)
    base_result = report.models["C1"]
    base_fit = base_result.full_fit
    base_reference = only(filter(reference ->
        reference.fit.fit_sha256 == base_fit.fit_sha256 &&
        reference.fit_role == "full_refit",
        SEA._graph_fit_references(report)))
    columns = Dict(name => index for (index, name) in enumerate(SEA.GRAPH_MODEL_HEADER))

    invalid_conditional = SEA.ConditionalFit(
        :fail,
        :conditional_fit_failed,
        copy(base_fit.conditional_fit.means),
        fill(NaN, 2, 2),
        -Inf,
        copy(base_fit.conditional_fit.starts),
        0,
        false,
    )
    invalid_conditional_fit = SEA._invalid_partition(
        "FAIL", :conditional_fit_failed, base_fit.model_id,
        copy(base_fit.training_dates), base_fit.reversed, base_fit.unary;
        residualizer=base_fit.residualizer,
        samples=copy(base_fit.samples),
        transform_edges=copy(base_fit.transform_edges),
        support=base_fit.support,
        null_fit=base_fit.null_fit,
        conditional_fit=invalid_conditional,
    )
    invalid_conditional_row = SEA._graph_model_row(
        merge(base_reference, (fit=invalid_conditional_fit,)))
    @test published_scale(invalid_conditional_row, columns, "conditional") === nothing
    @test published_scale(invalid_conditional_row, columns, "null") !== nothing
    @test all(invalid_conditional_row[columns["conditional_scale_$(suffix)"]] == "NA"
              for suffix in ("ff", "fb", "bb"))

    invalid_null = SEA.NullFit(
        :fail, :nonfinite, copy(base_fit.null_fit.mean), fill(NaN, 2, 2),
        -Inf, Float64[], false, 0, NaN)
    invalid_null_fit = SEA._invalid_partition(
        "FAIL", :nonfinite, base_fit.model_id,
        copy(base_fit.training_dates), base_fit.reversed, base_fit.unary;
        residualizer=base_fit.residualizer,
        samples=copy(base_fit.samples),
        transform_edges=copy(base_fit.transform_edges),
        support=base_fit.support,
        null_fit=invalid_null,
        conditional_fit=base_fit.conditional_fit,
    )
    invalid_null_row = SEA._graph_model_row(
        merge(base_reference, (fit=invalid_null_fit,)))
    @test published_scale(invalid_null_row, columns, "null") === nothing
    @test published_scale(invalid_null_row, columns, "conditional") !== nothing
    @test all(invalid_null_row[columns["null_scale_$(suffix)"]] == "NA"
              for suffix in ("ff", "fb", "bb"))

    valid_failed_conditional = SEA.ConditionalFit(
        :fail,
        :conditional_fit_failed,
        copy(base_fit.conditional_fit.means),
        copy(base_fit.conditional_fit.scale),
        base_fit.conditional_fit.objective,
        copy(base_fit.conditional_fit.starts),
        0,
        false,
    )
    valid_failed_fit = SEA._invalid_partition(
        "FAIL", :conditional_fit_failed, base_fit.model_id,
        copy(base_fit.training_dates), base_fit.reversed, base_fit.unary;
        residualizer=base_fit.residualizer,
        samples=copy(base_fit.samples),
        transform_edges=copy(base_fit.transform_edges),
        support=base_fit.support,
        null_fit=base_fit.null_fit,
        conditional_fit=valid_failed_conditional,
    )
    valid_failed_row = SEA._graph_model_row(
        merge(base_reference, (fit=valid_failed_fit,)))
    @test published_scale(valid_failed_row, columns, "conditional") !== nothing
    assert_certified_published_scale(valid_failed_row, columns, "conditional")
end

@testset "state-independent, phase, scan, and endpoint controls abstain" begin
    controls = Dict{Symbol,SEA.GateEvaluation}()
    for mechanism in (:state_independent, :phase_shuffled, :scan_artifact,
                      :endpoint_explained)
        controls[mechanism] = final_gate_for_mechanism(mechanism)
        @test controls[mechanism].status == "SKIPPED"
        @test controls[mechanism].reason in (
            :nonpositive_date,
            :bootstrap_not_positive,
            :shuffle_control_failed,
        )
        @test length(controls[mechanism].partitions) == 7
        @test all(length(partition.shuffle_values) == 500
                  for partition in controls[mechanism].partitions)
    end
    null_gate = controls[:state_independent]
    for partition in null_gate.partitions
        rank = count(value -> value < partition.control_date_mean,
                     partition.shuffle_values) / 500
        @test 0.0 <= rank <= 1.0
    end
end

@testset "strict Todo 10/universe loading and BLOCKED CLI receipt" begin
    fixture = make_cli_fixture()
    try
        loaded = load_cli_fixture(fixture)
        @test length(loaded.nodes) == 168
        @test length(loaded.edges) == 140
        @test [node.file for node in loaded.nodes] ==
              [node.file for node in fixture.data.nodes]
        @test all(length(value) == 64 for value in values(loaded.hashes))
        @test loaded.hashes["plan_sha256"] == SEA.PLAN_SHA256
        @test loaded.hashes["t7_review_sha256"] == SEA.T7_REVIEW_SHA256
        @test loaded.hashes["t10_review_sha256"] == SEA.T10_REVIEW_SHA256
        @test loaded.hashes["candidate_config_sha256"] ==
              SEA.CANDIDATE_CONFIG_SHA256
        @test loaded.hashes["model_config_sha256"] == SEA.MODEL_CONFIG_SHA256
        edge_lines = split(chomp(read(joinpath(fixture.edges,
                                               "edge_observations.tsv"), String)), '\n')[2:end]
        @test all(edge.raw_row_sha256 ==
                  bytes2hex(sha256(codeunits(line * "\n")))
                  for (edge, line) in zip(loaded.edges, edge_lines))
    finally
        rm(fixture.root; recursive=true, force=true)
    end

    forbidden = make_cli_fixture(forbidden_column=true)
    try
        error = try
            SEA.load_admission_data(
                forbidden.root;
                features=relpath(forbidden.features, forbidden.root),
                candidate_config="config/unit_assignment_structured_candidate.toml",
                model_config="config/unit_assignment_structured_model.toml",
                universe_dir=relpath(forbidden.universe, forbidden.root),
                edge_dir=relpath(forbidden.edges, forbidden.root),
                forward_receipt=relpath(forbidden.forward_receipt, forbidden.root),
                backward_receipt=relpath(forbidden.backward_receipt, forbidden.root),
            )
            nothing
        catch caught
            caught
        end
        @test error isa SEA.EdgeAdmissionError
        @test error.code == :admission_feature_schema

        output = joinpath(forbidden.root, "synthetic", "blocked-output")
        hashes = SEA.publish_blocker_receipt(
            forbidden.root,
            relpath(output, forbidden.root),
            error;
            hashes=Dict("input_sha256" => bytes2hex(sha256(read(forbidden.features)))),
        )
        @test Set(keys(hashes)) == Set(["receipt.toml"])
        @test readdir(output) == ["receipt.toml"]
        receipt_text = read(joinpath(output, "receipt.toml"), String)
        receipt = TOML.parse(receipt_text)
        @test receipt["status"] == "BLOCKED"
        @test receipt["reason"] == "admission_feature_schema"
        @test receipt["plan_sha256"] == SEA.PLAN_SHA256
        @test !occursin("slot_index", receipt_text)

        cli_output = joinpath(forbidden.root, "synthetic", "blocked-cli-output")
        command = `$(Base.julia_cmd()) --project=$(ROOT) $(EDGE_ADMISSION_CLI) --root $(forbidden.root) --features $(relpath(forbidden.features, forbidden.root)) --candidate-config config/unit_assignment_structured_candidate.toml --model-config config/unit_assignment_structured_model.toml --universe-dir $(relpath(forbidden.universe, forbidden.root)) --edge-dir $(relpath(forbidden.edges, forbidden.root)) --forward-receipt $(relpath(forbidden.forward_receipt, forbidden.root)) --backward-receipt $(relpath(forbidden.backward_receipt, forbidden.root)) --out-dir $(relpath(cli_output, forbidden.root))`
        stdout = IOBuffer()
        stderr = IOBuffer()
        process = run(pipeline(command; stdout=stdout, stderr=stderr); wait=false)
        wait(process)
        @test process.exitcode == 2
        @test isempty(String(take!(stdout)))
        @test occursin("BLOCKED", String(take!(stderr)))
        @test isfile(joinpath(cli_output, "receipt.toml"))
        @test TOML.parsefile(joinpath(cli_output, "receipt.toml"))["status"] ==
              "BLOCKED"
    finally
        rm(forbidden.root; recursive=true, force=true)
    end
end

@testset "support guards and deterministic SKIPPED/FAIL receipts" begin
    full = synthetic_admission_fixture()
    kept_dates = sort!(unique(node.date for node in full.nodes))[1:6]
    small = SEA.AdmissionData(
        [node for node in full.nodes if node.date in kept_dates],
        [edge for edge in full.edges if edge.date in kept_dates],
        copy(full.hashes),
    )
    normalized = SEA.normalize_admission_data(small)
    skipped = SEA.evaluate_model_admission(normalized, "C1")
    @test skipped.status == "SKIPPED"
    @test skipped.reason == :insufficient_dates
    skipped_report = SEA.AdmissionReport(
        "SKIPPED",
        :insufficient_dates,
        Dict("C1" => skipped),
        copy(small.hashes),
        skipped.result_sha256,
    )
    skipped_files = SEA.report_files(skipped_report)
    @test TOML.parse(String(skipped_files["receipt.toml"]))["status"] ==
          "SKIPPED"
    @test length(split(chomp(String(skipped_files["bootstrap.tsv"])), '\n')) ==
          501

    failed_model = SEA.ModelAdmissionResult(
        "FAIL",
        :nonfinite,
        "C1",
        skipped.outer_folds,
        skipped.final_gate,
        skipped.full_fit,
        bytes2hex(sha256(codeunits("synthetic-complete-numerical-failure\n"))),
    )
    failed_report = SEA.AdmissionReport(
        "FAIL",
        :numerical_failure,
        Dict("C1" => failed_model),
        copy(small.hashes),
        failed_model.result_sha256,
    )
    failed_receipt = TOML.parse(String(SEA.report_files(failed_report)["receipt.toml"]))
    @test failed_receipt["status"] == "FAIL"

    base_fit = SEA.fit_partition_model(SEA.normalize_admission_data(full), "C1",
                                       sort!(unique(node.date for node in full.nodes))[1:5])
    concentrated = SEA.EdgeSample[
        SEA.EdgeSample(sample.edge, sample.endpoint_predictors, sample.residual,
                       (1.0, 0.0, 0.0, 0.0), sample.raw_sha256)
        for sample in base_fit.samples
        if sample.edge.date in base_fit.training_dates
    ]
    support = SEA._support_evidence(concentrated, base_fit.training_dates)
    @test !support.sufficient
    @test support.reason == :insufficient_category_ess
    @test support.category_ess[2:4] == zeros(3)

    training_samples = [sample for sample in base_fit.samples
                        if sample.edge.date in base_fit.training_dates]
    four_dates = base_fit.training_dates[1:4]
    date_limited = [sample for sample in training_samples
                    if sample.edge.date in four_dates]
    @test SEA._support_evidence(date_limited, four_dates).reason ==
          :insufficient_dates

    removed_scan = first(sort!(unique(sample.edge.file for sample in training_samples)))
    scan_limited = [sample for sample in training_samples
                    if sample.edge.file != removed_scan]
    @test length(scan_limited) == 95
    @test SEA._support_evidence(scan_limited, base_fit.training_dates).reason ==
          :insufficient_scans

    edge_limited = SEA.EdgeSample[]
    by_scan = Dict{String,Vector{SEA.EdgeSample}}()
    for sample in training_samples
        push!(get!(by_scan, sample.edge.file, SEA.EdgeSample[]), sample)
    end
    for file in sort!(collect(keys(by_scan)))
        append!(edge_limited, sort!(by_scan[file]; by=sample -> sample.edge.ordinal)[1:4])
    end
    @test length(edge_limited) == 80
    @test SEA._support_evidence(edge_limited, base_fit.training_dates).reason ==
          :insufficient_edges

    sparse_category = SEA.EdgeSample[]
    for sample in training_samples
        scan_number = parse(Int, match(r"scan(\d+)", sample.edge.file).captures[1])
        weights = scan_number <= 2 ? (0.25, 0.25, 0.25, 0.25) :
                                     (1 / 3, 1 / 3, 1 / 3, 0.0)
        push!(sparse_category,
              SEA.EdgeSample(sample.edge, sample.endpoint_predictors,
                             sample.residual, weights, sample.raw_sha256))
    end
    sparse_support = SEA._support_evidence(sparse_category,
                                           base_fit.training_dates)
    @test all(sparse_support.category_ess .>= 10.0)
    @test sparse_support.minimum_support_scans[4] == 2
    @test sparse_support.reason == :insufficient_scan_support
end

println("structured_edge_admission_result_sha256=", SEA.focused_replay_hash())
println("structured_edge_admission_complete_sha256=", COMPLETE_REPLAY_SHA256[])
println("structured_edge_admission_graph_handoff_sha256=", GRAPH_HANDOFF_REPLAY_SHA256[])
