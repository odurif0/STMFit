module StructuredEdgeAdmission

using LinearAlgebra
using Printf
using Random
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "..", "hierarchical_unit_assignment.jl"))
include(joinpath(@__DIR__, "robust_emissions.jl"))
include(joinpath(@__DIR__, "edge_features.jl"))

const FIXED_NU::Int = 8
const SCALE_FLOOR::Float64 = 1.0e-4
const CONDITION_CAP::Float64 = 1.0e4
const MAX_ITER::Int = 200
const CONVERGENCE_TOL::Float64 = 1.0e-6
const DECREASE_TOL::Float64 = 1.0e-10
const START_TIE_TOL::Float64 = 1.0e-12
const START_ALPHAS = (0.0, 0.25, 0.5, 0.75, 1.0)
const DENOMINATOR_FLOOR::Float64 = 1.0e-12
const WEIGHT_SUM_ROUNDOFF::Float64 = 8.0 * eps(Float64)
const MINIMUM_CATEGORY_ESS::Float64 = 10.0
const MINIMUM_EDGES::Int = 90
const MINIMUM_DATES::Int = 5
const MINIMUM_SCANS::Int = 20
const MINIMUM_SUPPORT_SCANS_PER_DATE::Int = 3
const BOOTSTRAP_SEEDS = 0:499
const SHUFFLE_SEEDS = 0:499
const PREDICTOR_NAMES = (
    "amp_prominence",
    "amp_neighbor_ratio",
    "integrated_prominence",
    "amp_rel",
    "bwd_neg_com_t",
    "bwd_neg_diag45",
    "split_log_skew",
)
const GRAPH_HANDOFF_SCHEMA = "stmfit-structured-edge-model-handoff-v1"
const GRAPH_CATEGORY_ORDER = ("00", "01", "10", "11")
const GRAPH_STATE_ORDER = ("0", "1")
const GRAPH_COEFFICIENT_NAMES = (
    "intercept",
    "left_amp_prominence",
    "left_amp_neighbor_ratio",
    "left_integrated_prominence",
    "left_amp_rel",
    "left_bwd_neg_com_t",
    "left_bwd_neg_diag45",
    "left_split_log_skew",
    "right_amp_prominence",
    "right_amp_neighbor_ratio",
    "right_integrated_prominence",
    "right_amp_rel",
    "right_bwd_neg_com_t",
    "right_bwd_neg_diag45",
    "right_split_log_skew",
)
const GRAPH_OUTPUT_NAMES = ("corr_fwd", "corr_bwd")
const START_REASON_OK = :ok
const _SHUFFLE_CALL_COUNT = Ref(0)
const _FIT_CALL_COUNT = Ref(0)
const PLAN_SHA256 =
    "a9b386d613829e8f7e20b6e33f8e80898fa9a55f0b344dbb92ea64ac6804f3d0"
const T7_REVIEW_SHA256 =
    "8b0fb2c1a324a345de19d4e718278a52a6d568c0e474939cfb6bfe643e061635"
const T10_REVIEW_SHA256 =
    "02d421bb51dcdb4d658da791f14cc7a6cf5f6ced26f3f25216f0b6e6cb88cb77"
const T7_SOURCE_SHA256 =
    "0ca4a3863a4b9ed13aee5cd3a0229df4f2764a8524b764eeb0fd27238f14f241"
const T7_TEST_SHA256 =
    "c3deb8e82c2ba2a8720474d91f2294898b92372cf0bcbf53e1470e92a4260824"
const T10_CLI_SHA256 =
    "c2133b0bc45e5ae68ab49f8ae39470c7bd2ee4fe85a52da548fbc3f84096022b"
const T10_SOURCE_SHA256 =
    "c29fe06c7fc0ba33c4e168169fec5ced5a7907bbfcec626c183e3ced044da0f9"
const T10_TEST_SHA256 =
    "21f4896cc1db508ef2dd9d5acc21196a51c7ef066e475468725db73ece98bb22"
const CANDIDATE_CONFIG_SHA256 = StructuredEdgeFeatures.CANDIDATE_CONFIG_SHA256
const MODEL_CONFIG_SHA256 = StructuredEdgeFeatures.MODEL_CONFIG_SHA256
const ADMISSION_FEATURE_HEADER = vcat(
    ["file", "lobe", "t_nm", "amplitude"],
    collect(PREDICTOR_NAMES),
)
const EDGE_RECEIPT_SHA_FIELDS = (
    "edge_header_sha256",
    "node_header_sha256",
    "edge_observations_sha256",
    "node_segments_sha256",
    "feature_sha256",
    "candidate_config_sha256",
    "model_config_sha256",
    "universe_receipt_sha256",
    "keys_sha256",
    "grid_sha256",
    "forward_patch_sha256",
    "backward_patch_sha256",
    "forward_patch_receipt_sha256",
    "backward_patch_receipt_sha256",
    "forward_producer_sha256",
    "backward_producer_sha256",
    "source_sha256",
    "t2_review_sha256",
    "t3_review_sha256",
)
const VIEW_SPECS = (
    ("base_local", (1, 2, 3, 4)),
    ("base_local+bwd_neg_com_t", (1, 2, 3, 4, 5)),
    ("base_local+bwd_neg_diag45", (1, 2, 3, 4, 6)),
    ("base_local+split_log_skew", (1, 2, 3, 4, 7)),
)

export FIXED_NU,
       SCALE_FLOOR,
       CONDITION_CAP,
       MAX_ITER,
       CONVERGENCE_TOL,
       DECREASE_TOL,
       START_TIE_TOL,
       START_ALPHAS,
       MINIMUM_CATEGORY_ESS,
       MINIMUM_EDGES,
       MINIMUM_DATES,
       MINIMUM_SCANS,
       MINIMUM_SUPPORT_SCANS_PER_DATE,
       PREDICTOR_NAMES,
       EdgeAdmissionError,
       ScaleProjection,
       NullFit,
       ConditionalStartTrace,
       ConditionalFit,
       ResidualizerFit,
       log_student_t_full,
       pair_weights,
       kish_ess,
       fit_state_independent,
       fit_conditional,
       fit_residualizer,
       residualize,
       edge_gain,
       cyclic_offset,
       rotate_complete_weights,
       aggregate_scan_date,
       whole_scan_bootstrap,
       bootstrap_distribution,
       reversal_equivalent,
       combine_terminal_states,
       check_shuffle_terms,
       focused_replay_hash,
       AdmissionNode,
       AdmissionEdge,
       AdmissionData,
       NormalizedAdmissionData,
       UnaryFitResult,
       PartitionFitResult,
       PartitionEvaluation,
       GateEvaluation,
       ModelAdmissionResult,
       AdmissionReport,
       normalize_admission_data,
       fit_unary_responsibilities,
       fit_partition_model,
       evaluate_partition,
       evaluate_gate,
       evaluate_model_admission,
       evaluate_admission,
       PLAN_SHA256,
       T7_REVIEW_SHA256,
       T10_REVIEW_SHA256,
       CANDIDATE_CONFIG_SHA256,
       MODEL_CONFIG_SHA256,
       ADMISSION_FEATURE_HEADER,
       GRAPH_HANDOFF_SCHEMA,
       GRAPH_CATEGORY_ORDER,
       GRAPH_STATE_ORDER,
       GRAPH_COEFFICIENT_NAMES,
       GRAPH_OUTPUT_NAMES,
       graph_handoff_replay_hash,
       load_admission_data,
       report_files,
       publish_report,
       publish_blocker_receipt,
       source_bundle_sha256

struct EdgeAdmissionError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, error::EdgeAdmissionError) =
    print(io, "BLOCKED [", error.code, "]: ", error.message)

_fail(code::Symbol, message::AbstractString) =
    throw(EdgeAdmissionError(code, String(message)))

struct ScaleProjection
    status::Symbol
    scale::Matrix{Float64}
    eigenvalues::Vector{Float64}
    condition_number::Float64
end

struct NullProposal
    status::Symbol
    reason::Symbol
    mean::Vector{Float64}
    unprojected_scale::Matrix{Float64}
    scale::Matrix{Float64}
end

struct NullFit
    status::Symbol
    reason::Symbol
    mean::Vector{Float64}
    scale::Matrix{Float64}
    objective::Float64
    objective_trace::Vector{Float64}
    converged::Bool
    iterations::Int
    rejected_objective::Float64
end

struct ConditionalInitialState
    alpha::Float64
    means::Matrix{Float64}
    scale::Matrix{Float64}
    objective::Float64
end

struct ConditionalProposal
    status::Symbol
    reason::Symbol
    means::Matrix{Float64}
    unprojected_scale::Matrix{Float64}
    shrinkage::Float64
    shrunk_scale::Matrix{Float64}
    scale::Matrix{Float64}
end

struct ConditionalStartTrace
    start_index::Int
    alpha::Float64
    objective_trace::Vector{Float64}
    converged::Bool
    monotone::Bool
    valid::Bool
    reason::Symbol
    rejected_iteration::Int
    rejected_objective::Float64
    final_means::Matrix{Float64}
    final_scale::Matrix{Float64}
end

struct ConditionalFit
    status::Symbol
    reason::Symbol
    means::Matrix{Float64}
    scale::Matrix{Float64}
    objective::Float64
    starts::Vector{ConditionalStartTrace}
    best_start::Int
    converged::Bool
end

struct ResidualizerFit
    ridge::Float64
    coefficients::Matrix{Float64}
    training_row_count::Int
    training_sha256::String
    evaluation_order::Vector{Int}
end

_sha256(bytes) = bytes2hex(sha256(bytes))
_format_float(value::Float64) = @sprintf("%.17g", value)

function _as_residual_matrix(values)::Matrix{Float64}
    matrix = Matrix{Float64}(values)
    size(matrix, 2) == 2 || throw(DimensionMismatch("edge residuals must have two columns"))
    return matrix
end

function _as_weight_matrix(values, rows::Int)::Matrix{Float64}
    matrix = Matrix{Float64}(values)
    size(matrix) == (rows, 4) ||
        throw(DimensionMismatch("soft endpoint-pair weights must be N by 4"))
    return matrix
end

function _invalid_projection(reason::Symbol)::ScaleProjection
    return ScaleProjection(reason, fill(NaN, 2, 2), fill(NaN, 2), Inf)
end

function _canonical_scale_matrix(matrix::AbstractMatrix)
    size(matrix) == (2, 2) || return nothing
    fields = (
        _format_float(Float64(matrix[1, 1])),
        _format_float(Float64(matrix[1, 2])),
        _format_float(Float64(matrix[2, 2])),
    )
    values = try
        parse.(Float64, collect(fields))
    catch
        return nothing
    end
    all(isfinite, values) || return nothing
    return [values[1] values[2]; values[2] values[3]], fields
end

function _stored_scale_certificate(matrix::AbstractMatrix)
    size(matrix) == (2, 2) ||
        return (status=:shape, eigenvalues=Float64[], condition_number=Inf)
    matrix[1, 2] == matrix[2, 1] ||
        return (status=:asymmetric, eigenvalues=Float64[], condition_number=Inf)
    all(isfinite, matrix) ||
        return (status=:nonfinite, eigenvalues=Float64[], condition_number=Inf)
    decomposition = try
        eigen(Symmetric(Matrix{Float64}(matrix)))
    catch
        return (status=:projection_failed, eigenvalues=Float64[], condition_number=Inf)
    end
    all(isfinite, decomposition.values) && all(isfinite, decomposition.vectors) ||
        return (status=:nonfinite, eigenvalues=Float64[], condition_number=Inf)
    eigenvalues = Float64.(decomposition.values)
    all(value -> value > 0.0, eigenvalues) ||
        return (status=:nonpositive, eigenvalues=eigenvalues, condition_number=Inf)
    minimum(eigenvalues) >= SCALE_FLOOR ||
        return (status=:floor, eigenvalues=eigenvalues,
                condition_number=maximum(eigenvalues) / minimum(eigenvalues))
    condition_number = maximum(eigenvalues) / minimum(eigenvalues)
    isfinite(condition_number) ||
        return (status=:nonfinite_condition, eigenvalues=eigenvalues,
                condition_number=condition_number)
    condition_number <= CONDITION_CAP ||
        return (status=:condition, eigenvalues=eigenvalues,
                condition_number=condition_number)
    return (status=:ok, eigenvalues=eigenvalues, condition_number=condition_number)
end

function _closure_projection(stored::Matrix{Float64})::ScaleProjection
    q = max(
        nextfloat(stored[1, 1]) - stored[1, 1],
        nextfloat(stored[2, 2]) - stored[2, 2],
    )
    isfinite(q) && q > 0.0 || return _invalid_projection(:no_progress)
    while true
        candidate_input = stored + q .* Matrix{Float64}(I, 2, 2)
        all(isfinite, candidate_input) ||
            return _invalid_projection(:closure_nonfinite)
        canonical = _canonical_scale_matrix(candidate_input)
        canonical === nothing && return _invalid_projection(:closure_parse)
        candidate, _ = canonical
        certificate = _stored_scale_certificate(candidate)
        if certificate.status == :ok
            return ScaleProjection(:ok, candidate, certificate.eigenvalues,
                                   certificate.condition_number)
        elseif certificate.status ∉ (:floor, :condition)
            return _invalid_projection(Symbol("closure_", certificate.status))
        end
        next_q = q * 2.0
        isfinite(next_q) && next_q > q ||
            return _invalid_projection(:range_exhausted)
        q = next_q
    end
end

function _project_scale(matrix)::ScaleProjection
    input = Matrix{Float64}(matrix)
    size(input) == (2, 2) || throw(DimensionMismatch("edge scale must be 2 by 2"))
    all(isfinite, input) || return _invalid_projection(:nonfinite)
    symmetric = (input + input') / 2.0
    all(isfinite, symmetric) || return _invalid_projection(:nonfinite)
    decomposition = try
        eigen(Symmetric(symmetric))
    catch
        return _invalid_projection(:projection_failed)
    end
    all(isfinite, decomposition.values) && all(isfinite, decomposition.vectors) ||
        return _invalid_projection(:nonfinite)
    eigenvalues = max.(Float64.(decomposition.values), SCALE_FLOOR)
    projected = decomposition.vectors * Diagonal(eigenvalues) * decomposition.vectors'
    projected = (projected + projected') / 2.0
    all(isfinite, projected) || return _invalid_projection(:nonfinite)
    condition_number = maximum(eigenvalues) / minimum(eigenvalues)
    isfinite(condition_number) || return _invalid_projection(:nonfinite)
    condition_number <= CONDITION_CAP ||
        return ScaleProjection(:ill_conditioned, projected, eigenvalues, condition_number)
    canonical = _canonical_scale_matrix(projected)
    canonical === nothing && return _invalid_projection(:canonicalization_failed)
    stored, _ = canonical
    certificate = _stored_scale_certificate(stored)
    certificate.status == :ok &&
        return ScaleProjection(:ok, stored, certificate.eigenvalues,
                               certificate.condition_number)
    certificate.status in (:floor, :condition) ||
        return _invalid_projection(Symbol("stored_", certificate.status))
    return _closure_projection(stored)
end

function _scale_terms(scale::AbstractMatrix)
    size(scale) == (2, 2) || throw(DimensionMismatch("edge scale must be 2 by 2"))
    all(isfinite, scale) || throw(ArgumentError("edge scale must be finite"))
    scale[1, 2] == scale[2, 1] || throw(ArgumentError("edge scale must be symmetric"))
    first_diagonal = Float64(scale[1, 1])
    off_diagonal = Float64(scale[1, 2])
    second_diagonal = Float64(scale[2, 2])
    determinant = first_diagonal * second_diagonal - off_diagonal * off_diagonal
    isfinite(determinant) && determinant > 0.0 ||
        throw(ArgumentError("edge scale must be positive definite"))
    return first_diagonal, off_diagonal, second_diagonal, determinant
end

function _mahalanobis_two(
    residual::AbstractVector,
    mean::AbstractVector,
    scale::AbstractMatrix,
)::Float64
    length(residual) == 2 && length(mean) == 2 ||
        throw(DimensionMismatch("edge density vectors must have dimension two"))
    all(isfinite, residual) && all(isfinite, mean) ||
        throw(ArgumentError("edge density vectors must be finite"))
    first_diagonal, off_diagonal, second_diagonal, determinant = _scale_terms(scale)
    first = Float64(residual[1]) - Float64(mean[1])
    second = Float64(residual[2]) - Float64(mean[2])
    delta = (second_diagonal * first * first -
             2.0 * off_diagonal * first * second +
             first_diagonal * second * second) / determinant
    isfinite(delta) && delta >= -eps(Float64) ||
        throw(ArgumentError("edge Mahalanobis distance is invalid"))
    return max(delta, 0.0)
end

function log_student_t_full(
    residual::AbstractVector,
    mean::AbstractVector,
    scale::AbstractMatrix,
)::Float64
    _, _, _, determinant = _scale_terms(scale)
    delta = _mahalanobis_two(residual, mean, scale)
    value = log(24.0) - log(6.0) - log(FIXED_NU * π) -
            0.5 * log(determinant) -
            0.5 * (FIXED_NU + 2) * log1p(delta / FIXED_NU)
    isfinite(value) || throw(ArgumentError("edge log density is nonfinite"))
    return value
end

function pair_weights(probability_left::Real, probability_right::Real)
    left = Float64(probability_left)
    right = Float64(probability_right)
    isfinite(left) && 0.0 <= left <= 1.0 ||
        throw(ArgumentError("left unary probability must be finite in [0,1]"))
    isfinite(right) && 0.0 <= right <= 1.0 ||
        throw(ArgumentError("right unary probability must be finite in [0,1]"))
    result = (
        (1.0 - left) * (1.0 - right),
        (1.0 - left) * right,
        left * (1.0 - right),
        left * right,
    )
    total = sum(result)
    abs(total - 1.0) <= WEIGHT_SUM_ROUNDOFF ||
        throw(ArgumentError("soft endpoint-pair weights do not sum to one"))
    return result
end

function _validate_weights(weights::Matrix{Float64})
    all(isfinite, weights) || return :nonfinite
    all(value -> value >= 0.0, weights) || return :negative_weight
    for row in axes(weights, 1)
        abs(sum(view(weights, row, :)) - 1.0) <= WEIGHT_SUM_ROUNDOFF ||
            return :weight_sum
    end
    return :ok
end

function kish_ess(weights)::Vector{Float64}
    matrix = Matrix{Float64}(weights)
    size(matrix, 2) == 4 || throw(DimensionMismatch("Kish ESS needs four categories"))
    _validate_weights(matrix) == :ok || throw(ArgumentError("invalid pair weights"))
    result = Vector{Float64}(undef, 4)
    for category in 1:4
        values = view(matrix, :, category)
        denominator = sum(abs2, values)
        result[category] = denominator > 0.0 ? sum(values)^2 / denominator : 0.0
    end
    return result
end

function _null_objective(
    residuals::Matrix{Float64},
    mean::Vector{Float64},
    scale::Matrix{Float64},
)::Float64
    objective = 0.0
    for edge in axes(residuals, 1)
        objective += log_student_t_full(view(residuals, edge, :), mean, scale)
    end
    return objective
end

function _conditional_objective(
    residuals::Matrix{Float64},
    weights::Matrix{Float64},
    means::Matrix{Float64},
    scale::Matrix{Float64},
)::Float64
    objective = 0.0
    for edge in axes(residuals, 1)
        objective += weights[edge, 1] *
                     log_student_t_full(view(residuals, edge, :),
                                        view(means, 1, :), scale)
        objective += (weights[edge, 2] + weights[edge, 3]) *
                     log_student_t_full(view(residuals, edge, :),
                                        view(means, 2, :), scale)
        objective += weights[edge, 4] *
                     log_student_t_full(view(residuals, edge, :),
                                        view(means, 4, :), scale)
    end
    return objective
end

function _null_proposal(
    residuals_input,
    mean_input::AbstractVector,
    scale_input::AbstractMatrix,
)::NullProposal
    residuals = _as_residual_matrix(residuals_input)
    mean = Float64.(collect(mean_input))
    scale = Matrix{Float64}(scale_input)
    n = size(residuals, 1)
    n > 0 || return NullProposal(:fail, :empty, mean, fill(NaN, 2, 2), scale)
    all(isfinite, residuals) && all(isfinite, mean) && all(isfinite, scale) ||
        return NullProposal(:fail, :nonfinite, mean, fill(NaN, 2, 2), scale)
    latent = Vector{Float64}(undef, n)
    denominator = 0.0
    for edge in 1:n
        delta = _mahalanobis_two(view(residuals, edge, :), mean, scale)
        latent[edge] = (FIXED_NU + 2.0) / (FIXED_NU + delta)
        denominator += latent[edge]
    end
    isfinite(denominator) && denominator > DENOMINATOR_FLOOR ||
        return NullProposal(:fail, :denominator, mean, fill(NaN, 2, 2), scale)
    proposed_mean = zeros(Float64, 2)
    for edge in 1:n
        proposed_mean .+= latent[edge] .* view(residuals, edge, :)
    end
    proposed_mean ./= denominator
    all(isfinite, proposed_mean) ||
        return NullProposal(:fail, :nonfinite, mean, fill(NaN, 2, 2), scale)
    unprojected = zeros(Float64, 2, 2)
    for edge in 1:n
        difference = vec(residuals[edge, :]) - proposed_mean
        unprojected .+= latent[edge] .* (difference * difference')
    end
    unprojected ./= n
    projection = _project_scale(unprojected)
    projection.status == :ok ||
        return NullProposal(:fail, projection.status, proposed_mean,
                            unprojected, projection.scale)
    return NullProposal(:ok, :ok, proposed_mean, unprojected, projection.scale)
end

_accept_objective(old::Real, new::Real) =
    isfinite(new) && Float64(new) >= Float64(old) - DECREASE_TOL ?
    :accepted : :nonmonotone

function _invalid_null(reason::Symbol; mean=zeros(2), scale=fill(SCALE_FLOOR, 2, 2),
                       objective=-Inf, trace=Float64[], iterations=0,
                       rejected_objective=NaN)
    return NullFit(:fail, reason, Float64.(mean), Matrix{Float64}(scale),
                   Float64(objective), Float64.(trace), false, Int(iterations),
                   Float64(rejected_objective))
end

function fit_state_independent(residuals_input; iteration_limit::Int=MAX_ITER)::NullFit
    0 <= iteration_limit <= MAX_ITER ||
        throw(ArgumentError("iteration_limit must be in 0:$MAX_ITER"))
    residuals = _as_residual_matrix(residuals_input)
    n = size(residuals, 1)
    n > 0 || return _invalid_null(:empty)
    all(isfinite, residuals) || return _invalid_null(:nonfinite)
    mean = vec(sum(residuals; dims=1)) ./ n
    initial_scale = zeros(Float64, 2, 2)
    for edge in 1:n
        difference = vec(residuals[edge, :]) - mean
        initial_scale .+= difference * difference'
    end
    initial_scale ./= n
    projection = _project_scale(initial_scale)
    projection.status == :ok ||
        return _invalid_null(projection.status; mean=mean, scale=projection.scale)
    scale = projection.scale
    objective = try
        _null_objective(residuals, mean, scale)
    catch
        -Inf
    end
    isfinite(objective) ||
        return _invalid_null(:nonfinite_objective; mean=mean, scale=scale,
                             objective=objective, trace=[objective])
    trace = Float64[objective]
    for iteration in 1:iteration_limit
        proposal = try
            _null_proposal(residuals, mean, scale)
        catch
            NullProposal(:fail, :nonfinite, mean, fill(NaN, 2, 2), scale)
        end
        proposal.status == :ok ||
            return _invalid_null(proposal.reason; mean=mean, scale=scale,
                                 objective=objective, trace=trace,
                                 iterations=iteration - 1)
        proposed_objective = try
            _null_objective(residuals, proposal.mean, proposal.scale)
        catch
            -Inf
        end
        _accept_objective(objective, proposed_objective) == :accepted ||
            return _invalid_null(:nonmonotone; mean=mean, scale=scale,
                                 objective=objective, trace=trace,
                                 iterations=iteration - 1,
                                 rejected_objective=proposed_objective)
        previous = objective
        mean = proposal.mean
        scale = proposal.scale
        objective = proposed_objective
        push!(trace, objective)
        if abs(objective - previous) <= CONVERGENCE_TOL * (1.0 + abs(previous))
            return NullFit(:ok, :ok, mean, scale, objective, trace, true,
                           iteration, NaN)
        end
    end
    return _invalid_null(:nonconverged; mean=mean, scale=scale,
                         objective=objective, trace=trace,
                         iterations=iteration_limit)
end

function _raw_category_centers(
    residuals::Matrix{Float64},
    weights::Matrix{Float64},
)
    centers = zeros(Float64, 4, 2)
    for category in (1, 4)
        denominator = sum(view(weights, :, category))
        isfinite(denominator) && denominator > DENOMINATOR_FLOOR ||
            return nothing
        for edge in axes(residuals, 1)
            centers[category, :] .+= weights[edge, category] .* view(residuals, edge, :)
        end
        centers[category, :] ./= denominator
    end
    mixed_denominator = sum(weights[:, 2]) + sum(weights[:, 3])
    isfinite(mixed_denominator) && mixed_denominator > DENOMINATOR_FLOOR ||
        return nothing
    for edge in axes(residuals, 1)
        centers[2, :] .+= (weights[edge, 2] + weights[edge, 3]) .* view(residuals, edge, :)
    end
    centers[2, :] ./= mixed_denominator
    centers[3, :] .= centers[2, :]
    return centers
end

function _conditional_initial_states(
    residuals_input,
    weights_input,
    null_fit::NullFit,
)::Vector{ConditionalInitialState}
    residuals = _as_residual_matrix(residuals_input)
    weights = _as_weight_matrix(weights_input, size(residuals, 1))
    _validate_weights(weights) == :ok || return ConditionalInitialState[]
    null_fit.status == :ok || return ConditionalInitialState[]
    raw_centers = _raw_category_centers(residuals, weights)
    raw_centers === nothing && return ConditionalInitialState[]
    states = ConditionalInitialState[]
    repeated_null = repeat(null_fit.mean', 4, 1)
    for alpha in START_ALPHAS
        means = (1.0 - alpha) .* repeated_null .+ alpha .* raw_centers
        means[3, :] .= means[2, :]
        objective = try
            _conditional_objective(residuals, weights, means, null_fit.scale)
        catch
            -Inf
        end
        push!(states, ConditionalInitialState(alpha, means, copy(null_fit.scale),
                                              objective))
    end
    return states
end

function _conditional_proposal(
    residuals_input,
    weights_input,
    means_input::AbstractMatrix,
    scale_input::AbstractMatrix,
    null_scale_input::AbstractMatrix,
)::ConditionalProposal
    residuals = _as_residual_matrix(residuals_input)
    weights = _as_weight_matrix(weights_input, size(residuals, 1))
    means = Matrix{Float64}(means_input)
    scale = Matrix{Float64}(scale_input)
    null_scale = Matrix{Float64}(null_scale_input)
    size(means) == (4, 2) || throw(DimensionMismatch("conditional means must be 4 by 2"))
    _validate_weights(weights) == :ok ||
        return ConditionalProposal(:fail, :invalid_weights, means,
                                   fill(NaN, 2, 2), NaN,
                                   fill(NaN, 2, 2), scale)
    all(isfinite, residuals) && all(isfinite, means) && all(isfinite, scale) &&
    all(isfinite, null_scale) ||
        return ConditionalProposal(:fail, :nonfinite, means,
                                   fill(NaN, 2, 2), NaN,
                                   fill(NaN, 2, 2), scale)
    n = size(residuals, 1)
    latent = Matrix{Float64}(undef, n, 4)
    for edge in 1:n, category in 1:4
        delta = _mahalanobis_two(
            view(residuals, edge, :),
            view(means, category, :),
            scale,
        )
        latent[edge, category] = (FIXED_NU + 2.0) / (FIXED_NU + delta)
    end
    proposed_means = zeros(Float64, 4, 2)
    for category in (1, 4)
        denominator = 0.0
        for edge in 1:n
            product = weights[edge, category] * latent[edge, category]
            denominator += product
            proposed_means[category, :] .+= product .* view(residuals, edge, :)
        end
        isfinite(denominator) && denominator > DENOMINATOR_FLOOR ||
            return ConditionalProposal(:fail, :denominator, means,
                                       fill(NaN, 2, 2), NaN,
                                       fill(NaN, 2, 2), scale)
        proposed_means[category, :] ./= denominator
    end
    mixed_denominator = 0.0
    for edge in 1:n
        product = (weights[edge, 2] + weights[edge, 3]) * latent[edge, 2]
        mixed_denominator += product
        proposed_means[2, :] .+= product .* view(residuals, edge, :)
    end
    isfinite(mixed_denominator) && mixed_denominator > DENOMINATOR_FLOOR ||
        return ConditionalProposal(:fail, :denominator, means,
                                   fill(NaN, 2, 2), NaN,
                                   fill(NaN, 2, 2), scale)
    proposed_means[2, :] ./= mixed_denominator
    proposed_means[3, :] .= proposed_means[2, :]
    all(isfinite, proposed_means) ||
        return ConditionalProposal(:fail, :nonfinite, means,
                                   fill(NaN, 2, 2), NaN,
                                   fill(NaN, 2, 2), scale)

    denominator = sum(weights)
    isfinite(denominator) && denominator > DENOMINATOR_FLOOR ||
        return ConditionalProposal(:fail, :denominator, means,
                                   fill(NaN, 2, 2), NaN,
                                   fill(NaN, 2, 2), scale)
    unprojected = zeros(Float64, 2, 2)
    for edge in 1:n
        difference00 = vec(residuals[edge, :]) - vec(proposed_means[1, :])
        unprojected .+= weights[edge, 1] * latent[edge, 1] .*
                       (difference00 * difference00')
        difference_mixed = vec(residuals[edge, :]) - vec(proposed_means[2, :])
        unprojected .+= (weights[edge, 2] + weights[edge, 3]) * latent[edge, 2] .*
                       (difference_mixed * difference_mixed')
        difference11 = vec(residuals[edge, :]) - vec(proposed_means[4, :])
        unprojected .+= weights[edge, 4] * latent[edge, 4] .*
                       (difference11 * difference11')
    end
    unprojected ./= denominator
    shrinkage = min(0.5, 3.0 / n)
    shrunk = (1.0 - shrinkage) .* unprojected .+ shrinkage .* null_scale
    projection = _project_scale(shrunk)
    projection.status == :ok ||
        return ConditionalProposal(:fail, projection.status, proposed_means,
                                   unprojected, shrinkage, shrunk,
                                   projection.scale)
    return ConditionalProposal(:ok, :ok, proposed_means, unprojected,
                               shrinkage, shrunk, projection.scale)
end

function _invalid_conditional(reason::Symbol; starts=ConditionalStartTrace[])
    return ConditionalFit(:fail, reason, zeros(4, 2), fill(SCALE_FLOOR, 2, 2),
                          -Inf, starts, 0, false)
end

function _run_conditional_start(
    residuals::Matrix{Float64},
    weights::Matrix{Float64},
    null_fit::NullFit,
    state::ConditionalInitialState,
    start_index::Int,
    iteration_limit::Int,
)::ConditionalStartTrace
    means = copy(state.means)
    scale = copy(state.scale)
    objective = state.objective
    trace = Float64[objective]
    if !isfinite(objective)
        return ConditionalStartTrace(start_index, state.alpha, trace, false, true,
                                     false, :nonfinite_objective, 0, NaN,
                                     means, scale)
    end
    for iteration in 1:iteration_limit
        proposal = try
            _conditional_proposal(residuals, weights, means, scale,
                                  null_fit.scale)
        catch
            ConditionalProposal(:fail, :nonfinite, means, fill(NaN, 2, 2),
                                NaN, fill(NaN, 2, 2), scale)
        end
        if proposal.status != :ok
            return ConditionalStartTrace(start_index, state.alpha, trace, false,
                                         true, false, proposal.reason, iteration,
                                         NaN, means, scale)
        end
        proposed_objective = try
            _conditional_objective(residuals, weights, proposal.means,
                                   proposal.scale)
        catch
            -Inf
        end
        if _accept_objective(objective, proposed_objective) != :accepted
            return ConditionalStartTrace(start_index, state.alpha, trace, false,
                                         false, false, :nonmonotone, iteration,
                                         proposed_objective, means, scale)
        end
        previous = objective
        means = proposal.means
        scale = proposal.scale
        objective = proposed_objective
        push!(trace, objective)
        if abs(objective - previous) <= CONVERGENCE_TOL * (1.0 + abs(previous))
            return ConditionalStartTrace(start_index, state.alpha, trace, true,
                                         true, true, :ok, 0, NaN, means, scale)
        end
    end
    return ConditionalStartTrace(start_index, state.alpha, trace, false, true,
                                 false, :nonconverged, 0, NaN, means, scale)
end

function _select_best_start(objectives::AbstractVector, valid::AbstractVector{Bool})::Int
    length(objectives) == length(valid) ||
        throw(DimensionMismatch("start objective and validity vectors differ"))
    best = 0
    for index in eachindex(objectives)
        valid[index] && isfinite(objectives[index]) || continue
        if best == 0 || Float64(objectives[index]) >
                        Float64(objectives[best]) + START_TIE_TOL
            best = index
        end
    end
    return best
end

function fit_conditional(
    residuals_input,
    weights_input,
    null_fit::NullFit;
    iteration_limit::Int=MAX_ITER,
)::ConditionalFit
    0 <= iteration_limit <= MAX_ITER ||
        throw(ArgumentError("iteration_limit must be in 0:$MAX_ITER"))
    residuals = _as_residual_matrix(residuals_input)
    weights = _as_weight_matrix(weights_input, size(residuals, 1))
    null_fit.status == :ok || return _invalid_conditional(:null_failed)
    all(isfinite, residuals) || return _invalid_conditional(:nonfinite)
    weight_status = _validate_weights(weights)
    weight_status == :ok || return _invalid_conditional(weight_status)
    states = _conditional_initial_states(residuals, weights, null_fit)
    length(states) == length(START_ALPHAS) ||
        return _invalid_conditional(:denominator)
    starts = ConditionalStartTrace[]
    for (start_index, state) in enumerate(states)
        push!(starts, _run_conditional_start(
            residuals,
            weights,
            null_fit,
            state,
            start_index,
            iteration_limit,
        ))
    end
    final_objectives = [start.valid ? last(start.objective_trace) : -Inf
                        for start in starts]
    best = _select_best_start(final_objectives, [start.valid for start in starts])
    best > 0 || return _invalid_conditional(:no_valid_start; starts=starts)
    selected = starts[best]
    return ConditionalFit(:ok, :ok, selected.final_means,
                          selected.final_scale, last(selected.objective_trace),
                          starts, best, true)
end

function _matrix_hash(matrix::AbstractMatrix)::String
    io = IOBuffer()
    println(io, size(matrix, 1), 'x', size(matrix, 2))
    for value in matrix
        println(io, _format_float(Float64(value)))
    end
    return _sha256(take!(io))
end

function fit_residualizer(endpoint_predictors_input, raw_outputs_input)::ResidualizerFit
    endpoint_predictors = Matrix{Float64}(endpoint_predictors_input)
    raw_outputs = Matrix{Float64}(raw_outputs_input)
    size(endpoint_predictors, 2) == 14 ||
        throw(DimensionMismatch("residualizer needs fourteen endpoint predictors"))
    size(raw_outputs) == (size(endpoint_predictors, 1), 2) ||
        throw(DimensionMismatch("residualizer outputs must be N by 2"))
    size(endpoint_predictors, 1) > 0 ||
        throw(ArgumentError("residualizer training support is empty"))
    all(isfinite, endpoint_predictors) && all(isfinite, raw_outputs) ||
        throw(ArgumentError("residualizer training values must be finite"))
    design = hcat(ones(size(endpoint_predictors, 1)), endpoint_predictors)
    ridge = max(
        1.0e-12,
        1.0e-6 * tr(design[:, 2:end]' * design[:, 2:end]) / 14.0,
    )
    penalty = Diagonal(vcat(0.0, fill(ridge, 14)))
    coefficients = try
        (design' * design + penalty) \ (design' * raw_outputs)
    catch error
        throw(EdgeAdmissionError(:residualizer_failed, sprint(showerror, error)))
    end
    all(isfinite, coefficients) ||
        _fail(:residualizer_failed, "ridge coefficients are nonfinite")
    material = hcat(endpoint_predictors, raw_outputs)
    return ResidualizerFit(ridge, coefficients, size(endpoint_predictors, 1),
                           _matrix_hash(material), collect(1:15))
end

function _reverse_residualizer(
    canonical::ResidualizerFit,
    reversed_predictors::Matrix{Float64},
    raw_outputs::Matrix{Float64},
)::ResidualizerFit
    permutation = vcat(1, collect(9:15), collect(2:8))
    coefficients = canonical.coefficients[permutation, :]
    material = hcat(reversed_predictors, raw_outputs)
    return ResidualizerFit(canonical.ridge, coefficients,
                           canonical.training_row_count,
                           _matrix_hash(material), permutation)
end

function residualize(
    fit::ResidualizerFit,
    endpoint_predictors::AbstractVector,
    raw_output::AbstractVector,
)::Vector{Float64}
    length(endpoint_predictors) == 14 ||
        throw(DimensionMismatch("residualization needs fourteen endpoint predictors"))
    length(raw_output) == 2 ||
        throw(DimensionMismatch("edge observation must have dimension two"))
    all(isfinite, endpoint_predictors) && all(isfinite, raw_output) ||
        throw(ArgumentError("held-out residualization values must be finite"))
    prediction = _residualizer_prediction(fit, endpoint_predictors)
    result = Float64.(raw_output) - prediction
    all(isfinite, result) || _fail(:residualizer_failed, "held-out residual is nonfinite")
    return result
end

function _residualizer_prediction(
    fit::ResidualizerFit,
    endpoint_predictors::AbstractVector,
)::Vector{Float64}
    length(endpoint_predictors) == 14 ||
        throw(DimensionMismatch("prediction needs fourteen endpoint predictors"))
    all(isfinite, endpoint_predictors) ||
        throw(ArgumentError("prediction predictors must be finite"))
    design = vcat(1.0, Float64.(endpoint_predictors))
    prediction = zeros(Float64, 2)
    for index in fit.evaluation_order
        prediction[1] += design[index] * fit.coefficients[index, 1]
        prediction[2] += design[index] * fit.coefficients[index, 2]
    end
    all(isfinite, prediction) || _fail(:residualizer_failed, "prediction is nonfinite")
    return prediction
end

function edge_gain(
    residual::AbstractVector,
    weights::AbstractVector,
    conditional_means::AbstractMatrix,
    conditional_scale::AbstractMatrix,
    null_mean::AbstractVector,
    null_scale::AbstractMatrix,
)::Float64
    length(weights) == 4 || throw(DimensionMismatch("edge gain needs four weights"))
    size(conditional_means) == (4, 2) ||
        throw(DimensionMismatch("edge gain needs four two-dimensional means"))
    weight_values = Float64.(weights)
    all(isfinite, weight_values) && all(value -> value >= 0.0, weight_values) ||
        throw(ArgumentError("edge gain weights are invalid"))
    abs(sum(weight_values) - 1.0) <= WEIGHT_SUM_ROUNDOFF ||
        throw(ArgumentError("edge gain weights do not sum to one"))
    total = 0.0
    null_density = log_student_t_full(residual, null_mean, null_scale)
    total += weight_values[1] *
             (log_student_t_full(residual, view(conditional_means, 1, :),
                                 conditional_scale) - null_density)
    total += (weight_values[2] + weight_values[3]) *
             (log_student_t_full(residual, view(conditional_means, 2, :),
                                 conditional_scale) - null_density)
    total += weight_values[4] *
             (log_student_t_full(residual, view(conditional_means, 4, :),
                                 conditional_scale) - null_density)
    isfinite(total) || throw(ArgumentError("edge gain is nonfinite"))
    return total
end

function cyclic_offset(
    seed::Integer,
    file::AbstractString,
    segment_id::AbstractString,
    edge_count::Integer,
)::Int
    edge_count >= 2 || throw(ArgumentError("cyclic offset needs at least two edges"))
    0 <= seed <= typemax(Int) || throw(ArgumentError("shuffle seed is invalid"))
    file_text = String(file)
    segment_text = String(segment_id)
    isempty(file_text) && throw(ArgumentError("shuffle file is empty"))
    isempty(segment_text) && throw(ArgumentError("shuffle segment is empty"))
    occursin('\0', file_text) && throw(ArgumentError("shuffle file contains NUL"))
    occursin('\0', segment_text) && throw(ArgumentError("shuffle segment contains NUL"))
    digest = sha256(codeunits(string(seed, '\0', file_text, '\0', segment_text)))
    integer = UInt64(0)
    for byte in digest[1:8]
        integer = (integer << 8) | UInt64(byte)
    end
    return 1 + Int(integer % UInt64(edge_count - 1))
end

function rotate_complete_weights(vectors, offset::Integer)
    sequence = collect(vectors)
    count = length(sequence)
    count > 0 || return sequence
    0 <= offset <= count || throw(ArgumentError("cyclic offset is outside the segment"))
    offset == 0 && return sequence
    normalized = mod(offset, count)
    normalized == 0 && return sequence
    _SHUFFLE_CALL_COUNT[] += 1
    return [sequence[mod1(index - normalized, count)] for index in 1:count]
end

function check_shuffle_terms(enabled_terms)
    isempty(collect(enabled_terms)) ||
        _fail(:forbidden_shuffle_term,
              "a data-dependent edge term invalidates phase exchangeability")
    return nothing
end

function aggregate_scan_date(scan_scores::AbstractDict)
    isempty(scan_scores) && throw(ArgumentError("date score map is empty"))
    date_means = Dict{String,Float64}()
    for date in sort!(String.(collect(keys(scan_scores))))
        scans = scan_scores[date]
        isempty(scans) && throw(ArgumentError("date $date has no scan scores"))
        values = Float64[Float64(scans[scan]) for scan in sort!(collect(keys(scans)))]
        all(isfinite, values) || throw(ArgumentError("scan scores must be finite"))
        date_means[date] = mean(values)
    end
    overall = mean(date_means[date] for date in sort!(collect(keys(date_means))))
    return (date_means=date_means, overall=overall)
end

function whole_scan_bootstrap(scan_scores::AbstractDict, seed::Integer)::Float64
    seed in BOOTSTRAP_SEEDS || throw(ArgumentError("bootstrap seed must be in 0:499"))
    isempty(scan_scores) && throw(ArgumentError("date score map is empty"))
    random = MersenneTwister(Int(seed))
    date_values = Float64[]
    for date in sort!(collect(keys(scan_scores)))
        scans = scan_scores[date]
        names = sort!(collect(keys(scans)))
        isempty(names) && throw(ArgumentError("bootstrap date has no scans"))
        sampled = [Float64(scans[names[rand(random, eachindex(names))]]) for _ in eachindex(names)]
        all(isfinite, sampled) || throw(ArgumentError("bootstrap scan score is nonfinite"))
        push!(date_values, mean(sampled))
    end
    return mean(date_values)
end

bootstrap_distribution(scan_scores::AbstractDict)::Vector{Float64} =
    [whole_scan_bootstrap(scan_scores, seed) for seed in BOOTSTRAP_SEEDS]

function reversal_equivalent(
    observed::AbstractVector,
    reversed::AbstractVector;
    tolerance::Float64=1.0e-12,
)::Bool
    length(observed) == length(reversed) || return false
    all(isfinite, observed) && all(isfinite, reversed) || return false
    return all(abs(Float64(left) - Float64(right)) <= tolerance
               for (left, right) in zip(observed, reversed))
end

function combine_terminal_states(states)::String
    values = String.(collect(states))
    isempty(values) && throw(ArgumentError("terminal state list is empty"))
    allowed = Set(["PASS", "FAIL", "SKIPPED", "BLOCKED"])
    all(in(allowed), values) || _fail(:invalid_terminal_state, "unknown terminal state")
    "BLOCKED" in values && return "BLOCKED"
    "FAIL" in values && return "FAIL"
    "SKIPPED" in values && return "SKIPPED"
    return "PASS"
end

function focused_replay_hash()::String
    density = log_student_t_full([1.0, -2.0], [0.25, -1.5],
                                 [0.5 0.1; 0.1 3.0])
    weights = pair_weights(0.3, 0.8)
    offset = cyclic_offset(17, "20260101_phase.sxm", "0123456789abcdef", 5)
    io = IOBuffer()
    println(io, "schema=structured-edge-admission-focused-replay-v1")
    println(io, "density=", _format_float(density))
    println(io, "weights=", join(_format_float.(collect(weights)), ','))
    println(io, "offset=", offset)
    println(io, "bootstrap_count=", length(BOOTSTRAP_SEEDS))
    println(io, "shuffle_count=", length(SHUFFLE_SEEDS))
    return _sha256(take!(io))
end

struct AdmissionNode
    file::String
    date::String
    lobe::Int
    t_nm::Float64
    amplitude::Float64
    predictors::NTuple{7,Float64}
end

struct AdmissionEdge
    file::String
    date::String
    left_lobe::Int
    right_lobe::Int
    left_t_nm::Float64
    right_t_nm::Float64
    segment_id::String
    ordinal::Int
    observation::NTuple{2,Float64}
    raw_row_sha256::String
end

struct AdmissionData
    nodes::Vector{AdmissionNode}
    edges::Vector{AdmissionEdge}
    hashes::Dict{String,String}
end

struct NormalizationAudit
    file::String
    node_count::Int
    raw_sha256::String
    normalized_sha256::String
end

struct NormalizedAdmissionData
    data::AdmissionData
    normalized_predictors::Matrix{Float64}
    node_index::Dict{Tuple{String,Int},Int}
    audits::Vector{NormalizationAudit}
    normalization_sha256::String
end

struct UnaryViewModel
    name::String
    feature_indices::Vector{Int}
    family::Symbol
    status::Symbol
    reason::Symbol
    high_component::Int
    means::Matrix{Float64}
    scales::Matrix{Float64}
end

struct UnaryFitResult
    status::String
    reason::Symbol
    model_id::String
    probabilities_1::Vector{Float64}
    views::Vector{UnaryViewModel}
    training_dates::Vector{String}
    training_scans::Vector{String}
    node_identities::Vector{Tuple{String,String,Int,Float64}}
    node_order_sha256::String
    fit_sha256::String
    probability_sha256::String
end

struct EdgeSample
    edge::AdmissionEdge
    endpoint_predictors::Vector{Float64}
    residual::Vector{Float64}
    weights::NTuple{4,Float64}
    raw_sha256::String
end

struct EdgeTransformState
    edge::AdmissionEdge
    endpoint_predictors::Vector{Float64}
    status::String
    reason::Symbol
end

struct SupportEvidence
    dates::Int
    scans::Int
    edges::Int
    category_ess::Vector{Float64}
    minimum_support_scans::Vector{Int}
    free_parameters::Int
    effective_observation_ratio::Float64
    sufficient::Bool
    reason::Symbol
end

struct PartitionFitResult
    status::String
    reason::Symbol
    model_id::String
    training_dates::Vector{String}
    reversed::Bool
    unary::UnaryFitResult
    residualizer::Union{Nothing,ResidualizerFit}
    null_fit::Union{Nothing,NullFit}
    conditional_fit::Union{Nothing,ConditionalFit}
    samples::Vector{EdgeSample}
    transform_edges::Vector{EdgeTransformState}
    support::SupportEvidence
    training_sha256::String
    fit_sha256::String
end

struct PartitionEvaluation
    status::String
    reason::Symbol
    model_id::String
    target_date::String
    training_dates::Vector{String}
    fit::PartitionFitResult
    scan_scores::Dict{String,Float64}
    date_mean::Float64
    control_scan_scores::Dict{String,Float64}
    control_date_mean::Float64
    shuffle_values::Vector{Float64}
    shuffle_fit_sha256::Vector{String}
    shuffle_quantile::Float64
    shuffle_pass::Bool
    reversal_pass::Bool
    gains::Dict{String,Float64}
    score_sha256::String
end

struct GateEvaluation
    status::String
    reason::Symbol
    dates::Vector{String}
    date_means::Dict{String,Float64}
    bootstrap_values::Vector{Float64}
    bootstrap_lower::Float64
    every_date_positive::Bool
    every_shuffle_pass::Bool
    reversal_pass::Bool
    partitions::Vector{PartitionEvaluation}
    gate_sha256::String
end

struct OuterFoldEvidence
    outer_date::String
    inner_gate::GateEvaluation
    outer_score::PartitionEvaluation
    fit_sha256::String
end

struct ModelAdmissionResult
    status::String
    reason::Symbol
    model_id::String
    outer_folds::Vector{OuterFoldEvidence}
    final_gate::GateEvaluation
    full_fit::PartitionFitResult
    result_sha256::String
end

struct AdmissionReport
    status::String
    reason::Symbol
    models::Dict{String,ModelAdmissionResult}
    hashes::Dict{String,String}
    result_sha256::String
end

function _hash_lines(lines)::String
    io = IOBuffer()
    for line in lines
        println(io, line)
    end
    return _sha256(take!(io))
end

function _node_identity(node::AdmissionNode)::String
    return string(node.file, '\t', node.lobe, '\t', _format_float(node.t_nm))
end

function _node_identity_tuple(
    identity::Tuple{String,String,Int,Float64},
)::String
    file, date, lobe, t_nm = identity
    return string(file, '\t', date, '\t', lobe, '\t', _format_float(t_nm))
end

function _graph_node_identity(
    identity::Tuple{String,String,Int,Float64},
)::String
    file, date, lobe, t_nm = identity
    return join((file, date, string(lobe), _format_float(t_nm)), '|')
end

function _graph_edge_identity(edge::AdmissionEdge)::String
    return join((edge.file, string(edge.left_lobe), string(edge.right_lobe),
                 edge.segment_id), '|')
end

function _node_order_sha256(
    identities::Vector{Tuple{String,String,Int,Float64}},
)::String
    lines = String["schema=structured-edge-todo10-node-order-v1"]
    append!(lines, _node_identity_tuple(identity) for identity in identities)
    return _hash_lines(lines)
end

function _edge_identity(edge::AdmissionEdge)::String
    return string(
        edge.file,
        '\t',
        edge.left_lobe,
        '\t',
        edge.right_lobe,
        '\t',
        edge.segment_id,
    )
end

function _synthetic_edge_row_bytes(
    file::String,
    left_lobe::Int,
    right_lobe::Int,
    left_t_nm::Float64,
    right_t_nm::Float64,
    segment_id::String,
    observation::NTuple{2,Float64},
)::Vector{UInt8}
    fields = [
        file, string(left_lobe), string(right_lobe), _format_float(left_t_nm),
        _format_float(right_t_nm), _format_float(right_t_nm - left_t_nm),
        segment_id, segment_id, "eligible", "none",
        _format_float(observation[1]), _format_float(observation[2]),
        repeat("0", 64), repeat("0", 64), repeat("0", 64),
        repeat("0", 64), repeat("0", 64), repeat("0", 64),
    ]
    return Vector{UInt8}(codeunits(join(fields, '\t') * "\n"))
end

function AdmissionEdge(
    file,
    date,
    left_lobe,
    right_lobe,
    left_t_nm,
    right_t_nm,
    segment_id,
    ordinal,
    observation,
)
    bytes = _synthetic_edge_row_bytes(
        String(file), Int(left_lobe), Int(right_lobe), Float64(left_t_nm),
        Float64(right_t_nm), String(segment_id),
        (Float64(observation[1]), Float64(observation[2])),
    )
    return AdmissionEdge(
        String(file), String(date), Int(left_lobe), Int(right_lobe),
        Float64(left_t_nm), Float64(right_t_nm), String(segment_id),
        Int(ordinal), (Float64(observation[1]), Float64(observation[2])),
        _sha256(bytes),
    )
end

function _edge_endpoint_predictors(
    normalized::NormalizedAdmissionData,
    edge::AdmissionEdge;
    reversed::Bool=false,
)::Vector{Float64}
    left_index = normalized.node_index[(edge.file, edge.left_lobe)]
    right_index = normalized.node_index[(edge.file, edge.right_lobe)]
    left = collect(view(normalized.normalized_predictors, left_index, :))
    right = collect(view(normalized.normalized_predictors, right_index, :))
    return reversed ? vcat(right, left) : vcat(left, right)
end

function _edge_transform_states(
    normalized::NormalizedAdmissionData;
    reversed::Bool=false,
    residualizer::Union{Nothing,ResidualizerFit}=nothing,
)::Vector{EdgeTransformState}
    states = EdgeTransformState[]
    for edge in normalized.data.edges
        endpoint = _edge_endpoint_predictors(normalized, edge; reversed=reversed)
        finite = all(isfinite, endpoint)
        status, reason = if !finite
            "ineligible", :nonfinite_endpoint_predictor
        elseif residualizer === nothing
            "unavailable", :residualizer_unavailable
        else
            "eligible", :none
        end
        push!(states, EdgeTransformState(edge, endpoint, status, reason))
    end
    return states
end

function _validate_admission_data(data::AdmissionData)
    isempty(data.nodes) && _fail(:empty_input, "admission input has no nodes")
    node_keys = [(_node.file, _node.lobe) for _node in data.nodes]
    length(node_keys) == length(unique(node_keys)) ||
        _fail(:duplicate_node, "admission input repeats a node key")
    sorted_nodes = sort(data.nodes; by=node -> (node.file, node.t_nm, node.lobe))
    data.nodes == sorted_nodes ||
        _fail(:node_order, "admission nodes are not in observed-t order")
    index = Dict(key => row_index for (row_index, key) in enumerate(node_keys))
    edge_keys = Tuple{String,Int,Int}[]
    for edge in data.edges
        edge.file == basename(edge.file) ||
            _fail(:invalid_file, "scan identity must be a basename")
        edge.date == StructuredEdgeFeatures.StructuredUniverse.parse_scan_date(edge.file) ||
            _fail(:date_mismatch, "edge date differs from its scan identity")
        left_key = (edge.file, edge.left_lobe)
        right_key = (edge.file, edge.right_lobe)
        haskey(index, left_key) && haskey(index, right_key) ||
            _fail(:edge_key_mismatch, "edge endpoint is absent from nodes")
        left = data.nodes[index[left_key]]
        right = data.nodes[index[right_key]]
        left.t_nm == edge.left_t_nm && right.t_nm == edge.right_t_nm ||
            _fail(:edge_metadata_mismatch, "edge endpoint t differs from its node")
        left.t_nm < right.t_nm ||
            _fail(:edge_order, "edge endpoints are not in increasing observed t")
        all(isfinite, edge.observation) ||
            _fail(:nonfinite_edge, "edge observation is nonfinite")
        _require_sha(edge.raw_row_sha256, "edge raw row SHA-256")
        push!(edge_keys, (edge.file, edge.left_lobe, edge.right_lobe))
    end
    length(edge_keys) == length(unique(edge_keys)) ||
        _fail(:duplicate_edge, "admission input repeats an edge")
    data.edges == sort(data.edges; by=edge -> (
        edge.file,
        edge.left_t_nm,
        edge.right_t_nm,
        edge.left_lobe,
        edge.right_lobe,
    )) || _fail(:edge_order, "admission edges are not canonically ordered")
    return index
end

function normalize_admission_data(data::AdmissionData)::NormalizedAdmissionData
    node_index = _validate_admission_data(data)
    normalized = fill(NaN, length(data.nodes), 7)
    audits = NormalizationAudit[]
    for file in sort!(unique(node.file for node in data.nodes))
        indices = findall(node -> node.file == file, data.nodes)
        records = HierarchicalUnitAssignment.LobeRecord[
            HierarchicalUnitAssignment.LobeRecord(
                node.file,
                node.lobe,
                node.amplitude,
                Dict(PREDICTOR_NAMES[column] => node.predictors[column]
                     for column in 1:7),
            ) for node in data.nodes[indices]
        ]
        raw = reduce(vcat, (reshape(collect(node.predictors), 1, :)
                            for node in data.nodes[indices]))
        transformed = HierarchicalUnitAssignment.normalize_per_scan(
            records,
            raw,
            collect(PREDICTOR_NAMES),
        )
        normalized[indices, :] .= transformed
        push!(audits, NormalizationAudit(
            file,
            length(indices),
            _matrix_hash(raw),
            _matrix_hash(transformed),
        ))
    end
    material = String[]
    for audit in audits
        push!(material, join((audit.file, audit.node_count, audit.raw_sha256,
                             audit.normalized_sha256), '\t'))
    end
    return NormalizedAdmissionData(data, normalized, node_index, audits,
                                   _hash_lines(material))
end

function _strict_amplitude_orientation(
    probabilities::Matrix{Float64},
    amplitudes::Vector{Float64},
)
    assignments = [probabilities[row, 1] >= probabilities[row, 2] ? 1 : 2
                   for row in axes(probabilities, 1)]
    counts = [count(==(component), assignments) for component in 1:2]
    all(>(0), counts) || return 0, :collapsed_component
    means = [mean(amplitudes[assignments .== component]) for component in 1:2]
    all(isfinite, means) || return 0, :orientation_nonfinite
    means[1] == means[2] && return 0, :orientation_tie
    return means[1] > means[2] ? 1 : 2, :ok
end

function _invalid_unary(
    model_id::String,
    reason::Symbol,
    node_count::Int,
    training_dates::Vector{String},
    training_scans::Vector{String},
    views::Vector{UnaryViewModel};
    node_identities::Vector{Tuple{String,String,Int,Float64}}=
        Tuple{String,String,Int,Float64}[],
    status::String="FAIL",
)
    fit_hash = _hash_lines([
        "schema=structured-edge-unary-fit-v1",
        "model=$model_id",
        "status=$status",
        "reason=$(String(reason))",
        "training_dates=$(join(training_dates, ','))",
        "training_scans=$(join(training_scans, ','))",
    ])
    probabilities = fill(0.5, node_count)
    length(node_identities) == node_count ||
        _fail(:internal_contract, "invalid unary node identity count differs")
    return UnaryFitResult(status, reason, model_id, probabilities, views,
                          training_dates, training_scans, node_identities,
                          _node_order_sha256(node_identities), fit_hash,
                          _matrix_hash(reshape(probabilities, :, 1)))
end

function _fit_unary_view(
    normalized::NormalizedAdmissionData,
    model_id::String,
    name::String,
    feature_indices::Vector{Int},
    training_indices::Vector{Int},
)::UnaryViewModel
    finite_training = [index for index in training_indices
                       if all(isfinite, view(normalized.normalized_predictors,
                                             index, feature_indices))]
    length(finite_training) >= 2 ||
        return UnaryViewModel(name, feature_indices, Symbol(model_id), :skipped,
                              :insufficient_view_support, 0, zeros(2, length(feature_indices)),
                              fill(SCALE_FLOOR, 2, length(feature_indices)))
    training_matrix = normalized.normalized_predictors[finite_training, feature_indices]
    amplitudes = [normalized.data.nodes[index].amplitude for index in finite_training]
    if model_id == "C1"
        fit = try
            HierarchicalUnitAssignment.fit_em_two_component(
                training_matrix,
                0;
                n_starts=5,
                cov_floor=SCALE_FLOOR,
                max_iter=MAX_ITER,
                tol=CONVERGENCE_TOL,
            )
        catch
            return UnaryViewModel(name, feature_indices, :C1, :failed,
                                  :unary_fit_failed, 0,
                                  zeros(2, length(feature_indices)),
                                  fill(SCALE_FLOOR, 2, length(feature_indices)))
        end
        fit.converged && fit.monotone ||
            return UnaryViewModel(name, feature_indices, :C1, :failed,
                                  fit.monotone ? :unary_nonconverged : :unary_nonmonotone,
                                  0, fit.means, fit.vars)
        responsibilities = HierarchicalUnitAssignment.responsibilities(fit, training_matrix)
        high, reason = _strict_amplitude_orientation(responsibilities, amplitudes)
        reason == :ok ||
            return UnaryViewModel(name, feature_indices, :C1, :skipped, reason,
                                  0, fit.means, fit.vars)
        return UnaryViewModel(name, feature_indices, :C1, :ok, :ok, high,
                              fit.means, fit.vars)
    elseif model_id == "C2"
        fit = StructuredRobustEmissions.fit_student_t_two_component(
            training_matrix,
            amplitudes,
        )
        if fit.status == StructuredRobustEmissions.ROBUST_VALID
            return UnaryViewModel(name, feature_indices, :C2, :ok, :ok,
                                  fit.high_amplitude_component, fit.means,
                                  fit.scales)
        elseif fit.status == StructuredRobustEmissions.ROBUST_ABSTAINED
            return UnaryViewModel(name, feature_indices, :C2, :skipped,
                                  Symbol(string(fit.invalid_reason)), 0,
                                  fit.means, fit.scales)
        end
        return UnaryViewModel(name, feature_indices, :C2, :failed,
                              Symbol(string(fit.invalid_reason)), 0,
                              fit.means, fit.scales)
    end
    throw(ArgumentError("unary model must be C1 or C2"))
end

function _view_probabilities(
    fitted_view::UnaryViewModel,
    normalized::NormalizedAdmissionData,
)::Vector{Float64}
    probabilities = fill(NaN, length(normalized.data.nodes))
    fitted_view.status == :ok || return probabilities
    for index in eachindex(normalized.data.nodes)
        values = Base.view(normalized.normalized_predictors, index,
                           fitted_view.feature_indices)
        all(isfinite, values) || continue
        component_probabilities = if fitted_view.family == :C1
            fit = HierarchicalUnitAssignment.TwoComponentFit(
                fitted_view.means,
                fitted_view.scales,
                [0.5, 0.5],
                0.0,
                Float64[],
                true,
                true,
            )
            HierarchicalUnitAssignment.responsibilities(fit, reshape(collect(values), 1, :))
        else
            StructuredRobustEmissions.student_t_responsibilities(
                reshape(collect(values), 1, :),
                fitted_view.means,
                fitted_view.scales,
            )
        end
        probabilities[index] =
            component_probabilities[1, fitted_view.high_component]
    end
    return probabilities
end

function _unary_fit_hash(
    model_id::String,
    dates::Vector{String},
    scans::Vector{String},
    views::Vector{UnaryViewModel},
    training_nodes::Vector{AdmissionNode},
)::String
    lines = [
        "schema=structured-edge-unary-fit-v1",
        "model=$model_id",
        "training_dates=$(join(dates, ','))",
        "training_scans=$(join(scans, ','))",
    ]
    append!(lines, "node=$(_node_identity(node))" for node in training_nodes)
    for fitted in views
        push!(lines, join(("view", fitted.name, fitted.family, fitted.status,
                          fitted.reason, fitted.high_component), '\t'))
        for value in fitted.means
            push!(lines, "mean=$(_format_float(value))")
        end
        for value in fitted.scales
            push!(lines, "scale=$(_format_float(value))")
        end
    end
    return _hash_lines(lines)
end

function fit_unary_responsibilities(
    normalized::NormalizedAdmissionData,
    model_id::AbstractString,
    training_dates_input,
)::UnaryFitResult
    model = String(model_id)
    model in ("C1", "C2") || throw(ArgumentError("unary model must be C1 or C2"))
    training_dates = sort!(unique(String.(collect(training_dates_input))))
    isempty(training_dates) && throw(ArgumentError("unary training dates are empty"))
    training_indices = findall(node -> node.date in training_dates,
                               normalized.data.nodes)
    isempty(training_indices) && throw(ArgumentError("unary training nodes are empty"))
    training_scans = sort!(unique(normalized.data.nodes[index].file
                                  for index in training_indices))
    node_identities = Tuple{String,String,Int,Float64}[
        (node.file, node.date, node.lobe, node.t_nm)
        for node in normalized.data.nodes
    ]
    views = UnaryViewModel[]
    for (name, indices) in VIEW_SPECS
        push!(views, _fit_unary_view(normalized, model, name, collect(indices),
                                    training_indices))
    end
    any(view -> view.status == :failed, views) &&
        return _invalid_unary(model, :unary_numerical_failure,
                              length(normalized.data.nodes), training_dates,
                              training_scans, views;
                              node_identities=node_identities, status="FAIL")
    valid_views = filter(view -> view.status == :ok, views)
    isempty(valid_views) &&
        return _invalid_unary(model, :no_valid_unary_view,
                              length(normalized.data.nodes), training_dates,
                              training_scans, views;
                              node_identities=node_identities, status="SKIPPED")
    per_view = [_view_probabilities(view, normalized) for view in valid_views]
    probabilities = fill(0.5, length(normalized.data.nodes))
    for index in eachindex(probabilities)
        logits = Float64[]
        for values in per_view
            probability = values[index]
            isfinite(probability) || continue
            clipped = clamp(probability, 1.0e-12, 1.0 - 1.0e-12)
            push!(logits, log(clipped / (1.0 - clipped)))
        end
        if !isempty(logits)
            average = mean(logits)
            probabilities[index] = 1.0 / (1.0 + exp(-average))
        end
    end
    all(probability -> isfinite(probability) && 0.0 <= probability <= 1.0,
        probabilities) || _fail(:unary_probability, "unary probability is invalid")
    fit_hash = _unary_fit_hash(
        model,
        training_dates,
        training_scans,
        views,
        normalized.data.nodes[training_indices],
    )
    return UnaryFitResult(
        "PASS",
        :ok,
        model,
        probabilities,
        views,
        training_dates,
        training_scans,
        node_identities,
        _node_order_sha256(node_identities),
        fit_hash,
        _matrix_hash(reshape(probabilities, :, 1)),
    )
end

function _empty_support(reason::Symbol)::SupportEvidence
    return SupportEvidence(0, 0, 0, zeros(4), zeros(Int, 4), 9, 0.0,
                           false, reason)
end

function _support_evidence(samples::Vector{EdgeSample}, training_dates::Vector{String})
    weights = reduce(vcat, (reshape(collect(sample.weights), 1, :) for sample in samples);
                     init=zeros(0, 4))
    category_ess = isempty(samples) ? zeros(4) : kish_ess(weights)
    scans = sort!(unique(sample.edge.file for sample in samples))
    support_counts = fill(typemax(Int), 4)
    for category in 1:4
        for date in training_dates
            supported = Set{String}()
            for sample in samples
                sample.edge.date == date || continue
                sample.weights[category] > 0.0 && push!(supported, sample.edge.file)
            end
            support_counts[category] = min(support_counts[category], length(supported))
        end
    end
    isempty(training_dates) && (support_counts .= 0)
    ratio = length(samples) / 9.0
    conditions = (
        length(training_dates) >= MINIMUM_DATES,
        length(scans) >= MINIMUM_SCANS,
        length(samples) >= MINIMUM_EDGES,
        all(ess -> ess >= MINIMUM_CATEGORY_ESS, category_ess),
        all(count -> count >= MINIMUM_SUPPORT_SCANS_PER_DATE, support_counts),
        ratio >= 10.0,
    )
    reason = if !conditions[1]
        :insufficient_dates
    elseif !conditions[2]
        :insufficient_scans
    elseif !conditions[3]
        :insufficient_edges
    elseif !conditions[4]
        :insufficient_category_ess
    elseif !conditions[5]
        :insufficient_scan_support
    elseif !conditions[6]
        :insufficient_parameter_ratio
    else
        :ok
    end
    return SupportEvidence(length(training_dates), length(scans), length(samples),
                           category_ess, support_counts, 9, ratio,
                           all(conditions), reason)
end

function _invalid_partition(
    status::String,
    reason::Symbol,
    model_id::String,
    training_dates::Vector{String},
    reversed::Bool,
    unary::UnaryFitResult;
    residualizer=nothing,
    samples=EdgeSample[],
    transform_edges=EdgeTransformState[],
    support=_empty_support(reason),
    null_fit=nothing,
    conditional_fit=nothing,
)
    training_hash = _hash_lines([
        "schema=structured-edge-partition-training-v1",
        "model=$model_id",
        "dates=$(join(training_dates, ','))",
        "reversed=$reversed",
        "unary=$(unary.fit_sha256)",
    ])
    fit_hash = _hash_lines([
        "schema=structured-edge-partition-fit-v1",
        "status=$status",
        "reason=$(String(reason))",
        "training=$training_hash",
    ])
    return PartitionFitResult(status, reason, model_id, training_dates,
                              reversed, unary, residualizer, null_fit,
                              conditional_fit, samples, transform_edges,
                              support, training_hash, fit_hash)
end

function _partition_fit_hash(
    model_id::String,
    training_dates::Vector{String},
    reversed::Bool,
    unary::UnaryFitResult,
    residualizer::ResidualizerFit,
    null_fit::NullFit,
    conditional_fit::ConditionalFit,
    samples::Vector{EdgeSample},
)::String
    lines = [
        "schema=structured-edge-partition-fit-v1",
        "model=$model_id",
        "dates=$(join(training_dates, ','))",
        "reversed=$reversed",
        "unary=$(unary.fit_sha256)",
        "residualizer=$(residualizer.training_sha256)",
        "ridge=$(_format_float(residualizer.ridge))",
        "null=$(_format_float(null_fit.objective))",
        "conditional=$(_format_float(conditional_fit.objective))",
        "best_start=$(conditional_fit.best_start)",
    ]
    for value in residualizer.coefficients
        push!(lines, "B=$(_format_float(value))")
    end
    for value in null_fit.mean
        push!(lines, "null_mean=$(_format_float(value))")
    end
    for value in null_fit.scale
        push!(lines, "null_scale=$(_format_float(value))")
    end
    for value in conditional_fit.means
        push!(lines, "conditional_mean=$(_format_float(value))")
    end
    for value in conditional_fit.scale
        push!(lines, "conditional_scale=$(_format_float(value))")
    end
    append!(lines, "edge=$(_edge_identity(sample.edge))" for sample in samples)
    return _hash_lines(lines)
end

function fit_partition_model(
    normalized::NormalizedAdmissionData,
    model_id::AbstractString,
    training_dates_input;
    reversed::Bool=false,
)::PartitionFitResult
    _FIT_CALL_COUNT[] += 1
    model = String(model_id)
    training_dates = sort!(unique(String.(collect(training_dates_input))))
    transform_edges = _edge_transform_states(normalized; reversed=reversed)
    unary = fit_unary_responsibilities(normalized, model, training_dates)
    unary.status == "FAIL" &&
        return _invalid_partition("FAIL", unary.reason, model, training_dates,
                                  reversed, unary;
                                  transform_edges=transform_edges)
    unary.status == "SKIPPED" &&
        return _invalid_partition("SKIPPED", unary.reason, model,
                                  training_dates, reversed, unary;
                                  transform_edges=transform_edges)
    training_design = Vector{Vector{Float64}}()
    canonical_training_design = Vector{Vector{Float64}}()
    training_outputs = Vector{Vector{Float64}}()
    eligible_edges = Tuple{AdmissionEdge,Vector{Float64}}[]
    for state in transform_edges
        state.reason == :nonfinite_endpoint_predictor && continue
        edge = state.edge
        endpoint = state.endpoint_predictors
        canonical_endpoint = _edge_endpoint_predictors(normalized, edge)
        push!(eligible_edges, (edge, endpoint))
        if edge.date in training_dates
            push!(training_design, endpoint)
            push!(canonical_training_design, canonical_endpoint)
            push!(training_outputs, collect(edge.observation))
        end
    end
    isempty(training_design) &&
        return _invalid_partition("SKIPPED", :insufficient_edges, model,
                                  training_dates, reversed, unary;
                                  transform_edges=transform_edges)
    design_matrix = reduce(vcat, (reshape(row, 1, :) for row in training_design))
    canonical_design_matrix = reduce(
        vcat,
        (reshape(row, 1, :) for row in canonical_training_design),
    )
    output_matrix = reduce(vcat, (reshape(row, 1, :) for row in training_outputs))
    residualizer = try
        canonical = fit_residualizer(canonical_design_matrix, output_matrix)
        reversed ? _reverse_residualizer(canonical, design_matrix,
                                         output_matrix) : canonical
    catch
        return _invalid_partition("FAIL", :residualizer_failed, model,
                                  training_dates, reversed, unary;
                                  transform_edges=transform_edges)
    end
    transform_edges = _edge_transform_states(
        normalized;
        reversed=reversed,
        residualizer=residualizer,
    )
    samples = EdgeSample[]
    for (edge, endpoint) in eligible_edges
        left_index = normalized.node_index[(edge.file, edge.left_lobe)]
        right_index = normalized.node_index[(edge.file, edge.right_lobe)]
        left_probability = unary.probabilities_1[left_index]
        right_probability = unary.probabilities_1[right_index]
        weights = reversed ? pair_weights(right_probability, left_probability) :
                             pair_weights(left_probability, right_probability)
        raw = collect(edge.observation)
        residual = residualize(residualizer, endpoint, raw)
        raw_hash = edge.raw_row_sha256
        push!(samples, EdgeSample(edge, endpoint, residual, weights, raw_hash))
    end
    training_samples = [sample for sample in samples
                        if sample.edge.date in training_dates]
    support = _support_evidence(training_samples, training_dates)
    support.sufficient ||
        return _invalid_partition("SKIPPED", support.reason, model,
                                  training_dates, reversed, unary;
                                  residualizer=residualizer, samples=samples,
                                  transform_edges=transform_edges,
                                  support=support)
    residuals = reduce(vcat, (reshape(sample.residual, 1, :)
                              for sample in training_samples))
    weights = reduce(vcat, (reshape(collect(sample.weights), 1, :)
                            for sample in training_samples))
    null_fit = fit_state_independent(residuals)
    null_fit.status == :ok ||
        return _invalid_partition("FAIL", null_fit.reason, model,
                                  training_dates, reversed, unary;
                                  residualizer=residualizer, samples=samples,
                                  transform_edges=transform_edges,
                                  support=support, null_fit=null_fit)
    conditional = fit_conditional(residuals, weights, null_fit)
    conditional.status == :ok ||
        return _invalid_partition("FAIL", conditional.reason, model,
                                  training_dates, reversed, unary;
                                  residualizer=residualizer, samples=samples,
                                  transform_edges=transform_edges,
                                  support=support, null_fit=null_fit,
                                  conditional_fit=conditional)
    training_hash = _hash_lines([
        "schema=structured-edge-partition-training-v1",
        "model=$model",
        "dates=$(join(training_dates, ','))",
        "reversed=$reversed",
        "normalization=$(_hash_lines([
            join((audit.file, audit.raw_sha256, audit.normalized_sha256), '\t')
            for audit in normalized.audits
            if StructuredEdgeFeatures.StructuredUniverse.parse_scan_date(audit.file) in training_dates
        ]))",
        "unary=$(unary.fit_sha256)",
        "residualizer=$(residualizer.training_sha256)",
        ["edge=$(_edge_identity(sample.edge))" for sample in training_samples]...,
    ])
    fit_hash = _partition_fit_hash(model, training_dates, reversed, unary,
                                   residualizer, null_fit, conditional,
                                   training_samples)
    return PartitionFitResult("PASS", :ok, model, training_dates, reversed,
                              unary, residualizer, null_fit, conditional,
                              samples, transform_edges, support,
                              training_hash, fit_hash)
end

function _sample_gain(sample::EdgeSample, fit::PartitionFitResult)::Float64
    conditional = fit.conditional_fit::ConditionalFit
    null_fit = fit.null_fit::NullFit
    return edge_gain(sample.residual, collect(sample.weights), conditional.means,
                     conditional.scale, null_fit.mean, null_fit.scale)
end

function _contiguous_blocks(samples::Vector{EdgeSample})::Vector{Vector{EdgeSample}}
    grouped = Dict{Tuple{String,String},Vector{EdgeSample}}()
    for sample in samples
        push!(get!(grouped, (sample.edge.file, sample.edge.segment_id), EdgeSample[]),
              sample)
    end
    blocks = Vector{Vector{EdgeSample}}()
    for key in sort!(collect(keys(grouped)))
        ordered = sort!(grouped[key]; by=sample -> sample.edge.ordinal)
        current = EdgeSample[]
        previous = 0
        for sample in ordered
            if !isempty(current) && sample.edge.ordinal != previous + 1
                push!(blocks, current)
                current = EdgeSample[]
            end
            push!(current, sample)
            previous = sample.edge.ordinal
        end
        isempty(current) || push!(blocks, current)
    end
    return blocks
end

function _rotated_weight_matrix(samples::Vector{EdgeSample}, seed::Int)
    weights = [sample.weights for sample in samples]
    sample_index = Dict(_edge_identity(sample.edge) => index
                        for (index, sample) in enumerate(samples))
    permutable = Set{String}()
    for block in _contiguous_blocks(samples)
        length(block) >= 2 || continue
        offset = cyclic_offset(seed, first(block).edge.file,
                               first(block).edge.segment_id, length(block))
        rotated = rotate_complete_weights([sample.weights for sample in block], offset)
        for (sample, vector) in zip(block, rotated)
            identity = _edge_identity(sample.edge)
            weights[sample_index[identity]] = vector
            push!(permutable, identity)
        end
    end
    matrix = reduce(vcat, (reshape(collect(weight), 1, :) for weight in weights);
                    init=zeros(0, 4))
    return matrix, permutable
end

function _permutable_ids(samples::Vector{EdgeSample})::Set{String}
    result = Set{String}()
    for block in _contiguous_blocks(samples)
        length(block) >= 2 || continue
        for sample in block
            push!(result, _edge_identity(sample.edge))
        end
    end
    return result
end

function _scan_means(samples::Vector{EdgeSample}, gains::Dict{String,Float64};
                     restrict_to::Union{Nothing,Set{String}}=nothing)
    by_scan = Dict{String,Vector{Float64}}()
    for sample in samples
        identity = _edge_identity(sample.edge)
        restrict_to === nothing || identity in restrict_to || continue
        push!(get!(by_scan, sample.edge.file, Float64[]), gains[identity])
    end
    result = Dict{String,Float64}()
    for scan in sort!(collect(keys(by_scan)))
        isempty(by_scan[scan]) || (result[scan] = mean(by_scan[scan]))
    end
    return result
end

function _conditional_fit_hash(fit::ConditionalFit)::String
    lines = [
        "schema=structured-edge-conditional-fit-v1",
        "status=$(fit.status)",
        "reason=$(fit.reason)",
        "best_start=$(fit.best_start)",
        "objective=$(_format_float(fit.objective))",
    ]
    for start in fit.starts
        push!(lines, join(("start", start.start_index,
                          _format_float(start.alpha), start.valid,
                          start.reason,
                          isempty(start.objective_trace) ? "NA" :
                          _format_float(last(start.objective_trace))), '\t'))
    end
    return _hash_lines(lines)
end

function _invalid_evaluation(
    status::String,
    reason::Symbol,
    fit::PartitionFitResult,
    target_date::String,
)
    score_hash = _hash_lines([
        "schema=structured-edge-partition-score-v1",
        "status=$status",
        "reason=$(String(reason))",
        "fit=$(fit.fit_sha256)",
        "target=$target_date",
    ])
    return PartitionEvaluation(status, reason, fit.model_id, target_date,
                               fit.training_dates, fit,
                               Dict{String,Float64}(), NaN,
                               Dict{String,Float64}(), NaN, Float64[], String[],
                               NaN, false, false, Dict{String,Float64}(),
                               score_hash)
end

function evaluate_partition(
    fit::PartitionFitResult,
    target_date::AbstractString;
    run_shuffle::Bool=true,
    shuffle_fit_cache::Union{Nothing,Dict{String,ConditionalFit}}=nothing,
)::PartitionEvaluation
    target = String(target_date)
    target in fit.training_dates &&
        _fail(:partition_leakage, "score date appears in training dates")
    fit.status == "FAIL" && return _invalid_evaluation("FAIL", fit.reason, fit, target)
    fit.status == "SKIPPED" &&
        return _invalid_evaluation("SKIPPED", fit.reason, fit, target)
    heldout = [sample for sample in fit.samples if sample.edge.date == target]
    isempty(heldout) &&
        return _invalid_evaluation("BLOCKED", :missing_heldout_edges, fit, target)
    gains = Dict(_edge_identity(sample.edge) => _sample_gain(sample, fit)
                 for sample in heldout)
    scan_scores = _scan_means(heldout, gains)
    isempty(scan_scores) &&
        return _invalid_evaluation("BLOCKED", :missing_heldout_scores, fit, target)
    date_mean = mean(scan_scores[scan] for scan in sort!(collect(keys(scan_scores))))
    permutable = _permutable_ids(heldout)
    control_scan_scores = _scan_means(heldout, gains; restrict_to=permutable)
    isempty(control_scan_scores) &&
        return _invalid_evaluation("SKIPPED", :no_permutable_segment, fit, target)
    control_date_mean = mean(control_scan_scores[scan]
                             for scan in sort!(collect(keys(control_scan_scores))))
    shuffle_values = Float64[]
    shuffle_hashes = String[]
    if run_shuffle
        training_samples = [sample for sample in fit.samples
                            if sample.edge.date in fit.training_dates]
        residuals = reduce(vcat, (reshape(sample.residual, 1, :)
                                  for sample in training_samples))
        null_fit = fit.null_fit::NullFit
        for seed in SHUFFLE_SEEDS
            training_weights, _ = _rotated_weight_matrix(training_samples, seed)
            cache_key = string(fit.fit_sha256, '|', seed)
            shuffled_fit = if shuffle_fit_cache === nothing
                fit_conditional(residuals, training_weights, null_fit)
            else
                get!(shuffle_fit_cache, cache_key) do
                    fit_conditional(residuals, training_weights, null_fit)
                end
            end
            shuffled_fit.status == :ok ||
                return _invalid_evaluation("FAIL", :shuffle_fit_failed, fit,
                                           target)
            heldout_weights, shuffled_permutable =
                _rotated_weight_matrix(heldout, seed)
            shuffled_gains = Dict{String,Float64}()
            for (index, sample) in enumerate(heldout)
                identity = _edge_identity(sample.edge)
                identity in shuffled_permutable || continue
                shuffled_gains[identity] = edge_gain(
                    sample.residual,
                    view(heldout_weights, index, :),
                    shuffled_fit.means,
                    shuffled_fit.scale,
                    null_fit.mean,
                    null_fit.scale,
                )
            end
            scan_values = _scan_means(heldout, shuffled_gains;
                                      restrict_to=shuffled_permutable)
            isempty(scan_values) &&
                return _invalid_evaluation("FAIL", :shuffle_score_failed, fit,
                                           target)
            push!(shuffle_values,
                  mean(scan_values[scan]
                       for scan in sort!(collect(keys(scan_values)))))
            push!(shuffle_hashes, _conditional_fit_hash(shuffled_fit))
        end
    else
        shuffle_values = fill(-Inf, length(SHUFFLE_SEEDS))
        shuffle_hashes = fill("not_run", length(SHUFFLE_SEEDS))
    end
    shuffle_quantile = run_shuffle ? quantile(shuffle_values, 0.975) : -Inf
    shuffle_pass = !run_shuffle || control_date_mean > shuffle_quantile
    lines = [
        "schema=structured-edge-partition-score-v1",
        "fit=$(fit.fit_sha256)",
        "target=$target",
        "date_mean=$(_format_float(date_mean))",
        "control_mean=$(_format_float(control_date_mean))",
        "shuffle_quantile=$(_format_float(shuffle_quantile))",
        "shuffle_pass=$shuffle_pass",
    ]
    for scan in sort!(collect(keys(scan_scores)))
        push!(lines, "scan=$scan:$(_format_float(scan_scores[scan]))")
    end
    for identity in sort!(collect(keys(gains)))
        push!(lines, "gain=$identity:$(_format_float(gains[identity]))")
    end
    append!(lines, "shuffle=$(_format_float(value))" for value in shuffle_values)
    return PartitionEvaluation("PASS", :ok, fit.model_id, target,
                               fit.training_dates, fit, scan_scores, date_mean,
                               control_scan_scores, control_date_mean,
                               shuffle_values, shuffle_hashes,
                               shuffle_quantile, shuffle_pass, true, gains,
                               _hash_lines(lines))
end

function _with_reversal(
    evaluation::PartitionEvaluation,
    reversed::PartitionEvaluation,
)::PartitionEvaluation
    if evaluation.status != "PASS" || reversed.status != "PASS"
        status = combine_terminal_states([evaluation.status, reversed.status])
        reason = status == "FAIL" ? :reversal_fit_failed :
                 status == "BLOCKED" ? :reversal_blocked :
                 :reversal_unavailable
        return PartitionEvaluation(status, reason, evaluation.model_id,
                                   evaluation.target_date,
                                   evaluation.training_dates, evaluation.fit,
                                   evaluation.scan_scores,
                                   evaluation.date_mean,
                                   evaluation.control_scan_scores,
                                   evaluation.control_date_mean,
                                   evaluation.shuffle_values,
                                   evaluation.shuffle_fit_sha256,
                                   evaluation.shuffle_quantile,
                                   evaluation.shuffle_pass, false,
                                   evaluation.gains,
                                   evaluation.score_sha256)
    end
    keys_original = sort!(collect(keys(evaluation.gains)))
    keys_reversed = sort!(collect(keys(reversed.gains)))
    equivalent = keys_original == keys_reversed && reversal_equivalent(
        [evaluation.gains[key] for key in keys_original],
        [reversed.gains[key] for key in keys_reversed],
    )
    equivalent ||
        return PartitionEvaluation("FAIL", :reversal_inequivalent,
                                   evaluation.model_id, evaluation.target_date,
                                   evaluation.training_dates, evaluation.fit,
                                   evaluation.scan_scores,
                                   evaluation.date_mean,
                                   evaluation.control_scan_scores,
                                   evaluation.control_date_mean,
                                   evaluation.shuffle_values,
                                   evaluation.shuffle_fit_sha256,
                                   evaluation.shuffle_quantile,
                                   evaluation.shuffle_pass, false,
                                   evaluation.gains,
                                   evaluation.score_sha256)
    return PartitionEvaluation(evaluation.status, evaluation.reason,
                               evaluation.model_id, evaluation.target_date,
                               evaluation.training_dates, evaluation.fit,
                               evaluation.scan_scores, evaluation.date_mean,
                               evaluation.control_scan_scores,
                               evaluation.control_date_mean,
                               evaluation.shuffle_values,
                               evaluation.shuffle_fit_sha256,
                               evaluation.shuffle_quantile,
                               evaluation.shuffle_pass, true,
                               evaluation.gains, evaluation.score_sha256)
end

function _invalid_gate(
    status::String,
    reason::Symbol,
    dates::Vector{String};
    partitions::Vector{PartitionEvaluation}=PartitionEvaluation[],
)
    hash = _hash_lines([
        "schema=structured-edge-gate-v1",
        "status=$status",
        "reason=$(String(reason))",
        "dates=$(join(dates, ','))",
    ])
    return GateEvaluation(status, reason, dates, Dict{String,Float64}(),
                          Float64[], NaN, false, false, false,
                          partitions, hash)
end

function evaluate_gate(evaluations::Vector{PartitionEvaluation})::GateEvaluation
    isempty(evaluations) && return _invalid_gate("BLOCKED", :missing_partitions,
                                                String[])
    dates = sort!([evaluation.target_date for evaluation in evaluations])
    length(dates) == length(unique(dates)) ||
        return _invalid_gate("BLOCKED", :duplicate_partition, dates)
    combined = combine_terminal_states(evaluation.status for evaluation in evaluations)
    combined == "BLOCKED" &&
        return _invalid_gate("BLOCKED", :partition_blocked, dates;
                             partitions=evaluations)
    combined == "FAIL" &&
        return _invalid_gate("FAIL", :partition_failed, dates;
                             partitions=evaluations)
    combined == "SKIPPED" &&
        return _invalid_gate("SKIPPED", :partition_skipped, dates;
                             partitions=evaluations)
    date_means = Dict(evaluation.target_date => evaluation.date_mean
                      for evaluation in evaluations)
    every_positive = all(value -> value > 0.0, values(date_means))
    every_shuffle = all(evaluation -> evaluation.shuffle_pass, evaluations)
    reversal = all(evaluation -> evaluation.reversal_pass, evaluations)
    scan_scores = Dict(
        evaluation.target_date => copy(evaluation.scan_scores)
        for evaluation in evaluations
    )
    bootstrap = bootstrap_distribution(scan_scores)
    lower = quantile(bootstrap, 0.025)
    status = every_positive && lower > 0.0 && every_shuffle && reversal ?
             "PASS" : "SKIPPED"
    reason = !every_positive ? :nonpositive_date :
             !(lower > 0.0) ? :bootstrap_not_positive :
             !every_shuffle ? :shuffle_control_failed :
             !reversal ? :reversal_inequivalent : :ok
    lines = [
        "schema=structured-edge-gate-v1",
        "status=$status",
        "reason=$(String(reason))",
        "dates=$(join(dates, ','))",
        "bootstrap_lower=$(_format_float(lower))",
        "every_positive=$every_positive",
        "every_shuffle=$every_shuffle",
        "reversal=$reversal",
    ]
    for date in dates
        push!(lines, "date=$date:$(_format_float(date_means[date]))")
    end
    append!(lines, "bootstrap=$(_format_float(value))" for value in bootstrap)
    return GateEvaluation(status, reason, dates, date_means, bootstrap, lower,
                          every_positive, every_shuffle, reversal, evaluations,
                          _hash_lines(lines))
end

function _partition_cache_key(model::String, dates::Vector{String}, reversed::Bool)
    return string(model, '|', join(dates, ','), '|', reversed)
end

function evaluate_model_admission(
    normalized::NormalizedAdmissionData,
    model_id::AbstractString,
)::ModelAdmissionResult
    model = String(model_id)
    dates = sort!(unique(node.date for node in normalized.data.nodes))
    length(dates) >= MINIMUM_DATES + 2 || begin
        transform_edges = _edge_transform_states(normalized; reversed=false)
        dummy_unary = _invalid_unary(model, :insufficient_dates,
                                     length(normalized.data.nodes), dates,
                                     sort!(unique(node.file for node in normalized.data.nodes)),
                                     UnaryViewModel[];
                                     node_identities=Tuple{String,String,Int,Float64}[
                                         (node.file, node.date, node.lobe, node.t_nm)
                                         for node in normalized.data.nodes
                                     ],
                                     status="SKIPPED")
        dummy_fit = _invalid_partition("SKIPPED", :insufficient_dates,
                                       model, dates, false, dummy_unary;
                                       transform_edges=transform_edges)
        gate = _invalid_gate("SKIPPED", :insufficient_dates, dates)
        hash = _hash_lines(["model=$model", "status=SKIPPED",
                            "reason=insufficient_dates"])
        return ModelAdmissionResult("SKIPPED", :insufficient_dates, model,
                                    OuterFoldEvidence[], gate, dummy_fit, hash)
    end
    cache = Dict{String,PartitionFitResult}()
    shuffle_fit_cache = Dict{String,ConditionalFit}()
    get_fit = function(training_dates::Vector{String}, reversed::Bool)
        ordered = sort!(copy(training_dates))
        key = _partition_cache_key(model, ordered, reversed)
        return get!(cache, key) do
            fit_partition_model(normalized, model, ordered; reversed=reversed)
        end
    end
    outer_folds = OuterFoldEvidence[]
    for outer_date in dates
        outer_training = [date for date in dates if date != outer_date]
        inner_evaluations = PartitionEvaluation[]
        for inner_date in outer_training
            inner_training = [date for date in outer_training if date != inner_date]
            fit = get_fit(inner_training, false)
            reversed_fit = get_fit(inner_training, true)
            observed = evaluate_partition(
                fit,
                inner_date;
                run_shuffle=true,
                shuffle_fit_cache=shuffle_fit_cache,
            )
            reversed_score = evaluate_partition(reversed_fit, inner_date;
                                                  run_shuffle=false)
            push!(inner_evaluations, _with_reversal(observed, reversed_score))
        end
        inner_gate = evaluate_gate(inner_evaluations)
        outer_fit = get_fit(outer_training, false)
        reversed_outer_fit = get_fit(outer_training, true)
        outer_score = evaluate_partition(outer_fit, outer_date;
                                         run_shuffle=false)
        reversed_outer_score = evaluate_partition(reversed_outer_fit, outer_date;
                                                   run_shuffle=false)
        outer_score = _with_reversal(outer_score, reversed_outer_score)
        push!(outer_folds, OuterFoldEvidence(outer_date, inner_gate,
                                             outer_score,
                                             outer_fit.fit_sha256))
    end
    final_evaluations = PartitionEvaluation[]
    for heldout_date in dates
        training_dates = [date for date in dates if date != heldout_date]
        fit = get_fit(training_dates, false)
        reversed_fit = get_fit(training_dates, true)
        observed = evaluate_partition(
            fit,
            heldout_date;
            run_shuffle=true,
            shuffle_fit_cache=shuffle_fit_cache,
        )
        reversed_score = evaluate_partition(reversed_fit, heldout_date;
                                              run_shuffle=false)
        push!(final_evaluations, _with_reversal(observed, reversed_score))
    end
    final_gate = evaluate_gate(final_evaluations)
    full_fit = get_fit(dates, false)
    states = String[fold.inner_gate.status for fold in outer_folds]
    append!(states, [fold.outer_score.status for fold in outer_folds])
    push!(states, final_gate.status)
    push!(states, full_fit.status)
    status = combine_terminal_states(states)
    reason = status == "PASS" ? :ok :
             status == "FAIL" ? :numerical_failure :
             status == "BLOCKED" ? :integrity_failure :
             :conditional_mechanism_not_admitted
    lines = [
        "schema=structured-edge-model-admission-v1",
        "model=$model",
        "status=$status",
        "reason=$(String(reason))",
        "final_gate=$(final_gate.gate_sha256)",
        "full_fit=$(full_fit.fit_sha256)",
    ]
    for fold in outer_folds
        push!(lines, join(("outer", fold.outer_date,
                          fold.inner_gate.gate_sha256,
                          fold.outer_score.score_sha256,
                          fold.fit_sha256), '\t'))
    end
    return ModelAdmissionResult(status, reason, model, outer_folds,
                                final_gate, full_fit, _hash_lines(lines))
end

function evaluate_admission(
    data::AdmissionData;
    enabled_terms=(),
)::AdmissionReport
    check_shuffle_terms(enabled_terms)
    normalized = normalize_admission_data(data)
    models = Dict{String,ModelAdmissionResult}()
    for model in ("C1", "C2")
        models[model] = evaluate_model_admission(normalized, model)
    end
    status = combine_terminal_states(result.status for result in values(models))
    reason = status == "PASS" ? :ok :
             status == "FAIL" ? :numerical_failure :
             status == "BLOCKED" ? :integrity_failure :
             :conditional_mechanism_not_admitted
    lines = [
        "schema=structured-edge-admission-report-v1",
        "status=$status",
        "reason=$(String(reason))",
        "normalization=$(normalized.normalization_sha256)",
    ]
    for model in sort!(collect(keys(models)))
        push!(lines, "model=$model:$(models[model].result_sha256)")
    end
    for key in sort!(collect(keys(data.hashes)))
        push!(lines, "input=$key:$(data.hashes[key])")
    end
    result_hash = _hash_lines(lines)
    hashes = copy(data.hashes)
    hashes["normalization_sha256"] = normalized.normalization_sha256
    return AdmissionReport(status, reason, models, hashes, result_hash)
end

struct UniverseInput
    directory::String
    receipt_sha256::String
    feature_sha256::String
    universe_sha256::String
    keys_sha256::String
    scans::Dict{String,String}
    keys::Vector{Tuple{String,Int}}
end

struct EdgeBundleInput
    directory::String
    receipt_sha256::String
    edge_sha256::String
    node_sha256::String
    feature_sha256::String
    source_sha256::String
    t3_review_sha256::String
    nodes::Dict{Tuple{String,Int},Tuple{Float64,String}}
    edges::Vector{AdmissionEdge}
end

function _repository_root()::String
    return normpath(joinpath(@__DIR__, "..", "..", ".."))
end

function _sha256_file(path::AbstractString)::String
    return _sha256(read(String(path)))
end

function _require_sha(value, context::AbstractString)::String
    value isa String || _fail(:invalid_hash, "$context is not a string")
    text = String(value)
    occursin(r"^[0-9a-f]{64}$", text) ||
        _fail(:invalid_hash, "$context is not a lowercase SHA-256")
    return text
end

function _validated_edge_receipt_hashes(receipt::AbstractDict)::Dict{String,String}
    hashes = Dict{String,String}()
    for field in EDGE_RECEIPT_SHA_FIELDS
        hashes[field] = _require_sha(
            _required_string(receipt, field, "edge receipt"),
            "edge receipt.$field",
        )
    end
    return hashes
end

function _exact_keys(table::AbstractDict, expected::Set{String}, context::String)
    observed = Set(String(key) for key in keys(table))
    observed == expected || _fail(:schema_mismatch, "$context keys differ")
    return nothing
end

function _required_string(table::AbstractDict, key::String, context::String)::String
    value = get(table, key, nothing)
    value isa String || _fail(:schema_mismatch, "$context.$key is not a string")
    return String(value)
end

function _required_integer(table::AbstractDict, key::String, context::String)::Int
    value = get(table, key, nothing)
    value isa Integer && !(value isa Bool) ||
        _fail(:schema_mismatch, "$context.$key is not an integer")
    return Int(value)
end

function _finite_text(text::String, context::String)::Float64
    value = tryparse(Float64, text)
    value === nothing && _fail(:invalid_number, "$context is not numerical")
    isfinite(value) || _fail(:invalid_number, "$context is nonfinite")
    return value
end

function _optional_text(text::String, context::String)::Float64
    text == "NA" && return NaN
    return _finite_text(text, context)
end

function _positive_integer_text(text::String, context::String)::Int
    occursin(r"^[1-9][0-9]*$", text) ||
        _fail(:invalid_key, "$context is not a canonical positive integer")
    value = tryparse(Int, text)
    value === nothing && _fail(:invalid_key, "$context is outside Int")
    return value
end

function _verify_fixed_file(path::String, expected::String, context::String)
    isfile(path) && !islink(path) || _fail(:missing_dependency, "$context is absent")
    observed = _sha256_file(path)
    observed == expected || _fail(:dependency_hash_mismatch,
                                  "$context SHA-256 differs")
    return observed
end

function _verify_dependencies(root::String)
    dependencies = Dict(
        joinpath(root, ".omo", "plans", "structured-label-free-unit-assignment.md") =>
            PLAN_SHA256,
        joinpath(root, ".omo", "evidence", "structured-label-free-unit-assignment",
                 "t7", "correction", "review", "AdversarialVerify.json") =>
            T7_REVIEW_SHA256,
        joinpath(root, ".omo", "evidence", "structured-label-free-unit-assignment",
                 "provenance-rebind", "phase4-t10-edge-features", "review", "AdversarialVerify.json") =>
            T10_REVIEW_SHA256,
        joinpath(root, "test", "lib", "structured_assignment", "robust_emissions.jl") =>
            T7_SOURCE_SHA256,
        joinpath(root, "test", "test_structured_robust_emissions.jl") =>
            T7_TEST_SHA256,
        joinpath(root, "test", "build_label_free_edge_features.jl") =>
            T10_CLI_SHA256,
        joinpath(root, "test", "lib", "structured_assignment", "edge_features.jl") =>
            T10_SOURCE_SHA256,
        joinpath(root, "test", "test_structured_edge_features.jl") =>
            T10_TEST_SHA256,
    )
    for path in sort!(collect(keys(dependencies)))
        _verify_fixed_file(path, dependencies[path], relpath(path, root))
    end
    for path in (
        joinpath(root, ".omo", "evidence", "structured-label-free-unit-assignment",
                 "t7", "correction", "review", "AdversarialVerify.json"),
        joinpath(root, ".omo", "evidence", "structured-label-free-unit-assignment",
                 "provenance-rebind", "phase4-t10-edge-features", "review", "AdversarialVerify.json"),
    )
        occursin(r"\"verdict\"\s*:\s*\"confirmed\"", read(path, String)) ||
            _fail(:dependency_verdict_mismatch,
                  "confirmed dependency review is not confirmed")
    end
    return dependencies
end

function _todo10_receipt_artifacts(
    root::String,
    supplied::AbstractString,
    channel::Symbol,
)
    context = "$(String(channel)) producer receipt"
    path = StructuredEdgeFeatures._resolve_existing(root, supplied, context)
    bytes, snapshot = StructuredEdgeFeatures._snapshot(path, context)
    receipt = StructuredEdgeFeatures._parse_toml(bytes, context)
    StructuredEdgeFeatures._exact_keys(
        receipt,
        StructuredEdgeFeatures.PATCH_RECEIPT_KEYS,
        context,
    )
    features_path = StructuredEdgeFeatures._resolve_existing(
        root,
        StructuredEdgeFeatures._required_string(receipt, "features_path", context),
        "$context feature table",
    )
    output_path = StructuredEdgeFeatures._resolve_existing(
        root,
        StructuredEdgeFeatures._required_string(receipt, "output_path", context),
        "$context patch table",
    )
    return (
        path=path,
        snapshot=snapshot,
        receipt=receipt,
        features_path=features_path,
        output_path=output_path,
    )
end

function _exact_table_row_bytes(bytes::Vector{UInt8}, context::String)
    text = String(copy(bytes))
    isvalid(text) || _fail(:invalid_encoding, "$context is not valid UTF-8")
    occursin('\r', text) && _fail(:carriage_return, "$context must use LF line endings")
    endswith(text, '\n') || _fail(:missing_final_lf, "$context lacks final LF")
    lines = split(chomp(text), '\n'; keepempty=true)
    length(lines) >= 2 || _fail(:empty_input, "$context has no data rows")
    return [Vector{UInt8}(codeunits(line * "\n")) for line in lines[2:end]]
end

function _todo10_frozen_adapter(
    root::String,
    feature_path::String,
    feature_bytes::Vector{UInt8},
    universe::UniverseInput,
)
    structured_universe = StructuredEdgeFeatures.StructuredUniverse
    scans, keys = structured_universe._parse_features(feature_bytes)
    observed_keys = [(key.file, key.lobe) for key in keys]
    observed_keys == universe.keys ||
        _fail(:universe_key_mismatch,
              "Todo 10 feature keys differ from the supplied universe")
    observed_scans = Dict(scan.file => scan.date for scan in scans)
    observed_scans == universe.scans ||
        _fail(:universe_key_mismatch,
              "Todo 10 feature scans differ from the supplied universe")
    keys_sha256 = _sha256(structured_universe._key_identity_bytes(keys))
    keys_sha256 == universe.keys_sha256 ||
        _fail(:universe_key_mismatch,
              "Todo 10 canonical key digest differs from the supplied universe")
    return structured_universe.FrozenUniverse(
        root,
        feature_path,
        joinpath(root, "config", "unit_assignment_structured_candidate.toml"),
        joinpath(root, "config", "unit_assignment_structured_model.toml"),
        universe.feature_sha256,
        CANDIDATE_CONFIG_SHA256,
        MODEL_CONFIG_SHA256,
        repeat("0", 64),
        StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256,
        StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256,
        StructuredEdgeFeatures.GRID_SHA256,
        scans,
        keys,
        keys_sha256,
        (),
        (),
        (),
        (),
        (),
    )
end

function _validate_todo10_provenance(
    root::String,
    forward_receipt::AbstractString,
    backward_receipt::AbstractString,
    universe::UniverseInput,
    hashes::Dict{String,String},
    forward_command::String,
    backward_command::String,
)
    for (field, parts, expected, context) in (
        ("t2_review_sha256",
         ("t2", "correction2", "review", "AdversarialVerify.json"),
         StructuredEdgeFeatures.T2_REVIEW_SHA256,
         "Todo 2 correction review"),
        ("t3_review_sha256",
         ("provenance-rebind", "phase3-t3-universe", "correction", "review", "AdversarialVerify.json"),
         StructuredEdgeFeatures.T3_REVIEW_SHA256,
         "Todo 3 correction review"),
    )
        hashes[field] == expected ||
            _fail(:dependency_hash_mismatch, "$context receipt binding differs")
        path = joinpath(root, ".omo", "evidence",
                        "structured-label-free-unit-assignment", parts...)
        StructuredEdgeFeatures._validate_fixed_review(path, expected, context)
    end

    forward_info = _todo10_receipt_artifacts(root, forward_receipt, :forward)
    backward_info = _todo10_receipt_artifacts(root, backward_receipt, :backward)
    forward_info.snapshot.sha256 == hashes["forward_patch_receipt_sha256"] ||
        _fail(:patch_receipt_hash_mismatch,
              "forward producer receipt bytes differ from the edge receipt")
    backward_info.snapshot.sha256 == hashes["backward_patch_receipt_sha256"] ||
        _fail(:patch_receipt_hash_mismatch,
              "backward producer receipt bytes differ from the edge receipt")
    forward_info.features_path == backward_info.features_path ||
        _fail(:producer_path_mismatch,
              "forward and backward producer receipts bind different feature tables")

    feature_bytes = read(forward_info.features_path)
    _sha256(feature_bytes) == universe.feature_sha256 ||
        _fail(:feature_hash_mismatch,
              "Todo 10 feature bytes differ from the edge receipt")
    _sha256_file(forward_info.output_path) == hashes["forward_patch_sha256"] ||
        _fail(:patch_hash_mismatch,
              "forward patch bytes differ from the edge receipt")
    _sha256_file(backward_info.output_path) == hashes["backward_patch_sha256"] ||
        _fail(:patch_hash_mismatch,
              "backward patch bytes differ from the edge receipt")

    frozen = _todo10_frozen_adapter(
        root,
        forward_info.features_path,
        feature_bytes,
        universe,
    )
    nodes = StructuredEdgeFeatures._parse_feature_nodes(feature_bytes, frozen)
    forward = StructuredEdgeFeatures._parse_patch_table(
        forward_info.output_path,
        :forward,
        nodes,
        frozen,
    )
    backward = StructuredEdgeFeatures._parse_patch_table(
        backward_info.output_path,
        :backward,
        nodes,
        frozen,
    )
    forward_provenance = StructuredEdgeFeatures._validate_patch_receipt(
        forward_info.path,
        :forward,
        root,
        forward_info.features_path,
        forward,
        backward_info.output_path,
        frozen,
    )
    backward_provenance = StructuredEdgeFeatures._validate_patch_receipt(
        backward_info.path,
        :backward,
        root,
        forward_info.features_path,
        backward,
        forward_info.output_path,
        frozen,
    )
    forward_provenance.data_path == backward_provenance.data_path ||
        _fail(:producer_path_mismatch,
              "forward and backward producer receipts bind different data directories")
    forward_command == forward_provenance.command ||
        _fail(:producer_command_mismatch,
              "edge forward command differs from the canonical Todo 10 command")
    backward_command == backward_provenance.command ||
        _fail(:producer_command_mismatch,
              "edge backward command differs from the canonical Todo 10 command")
    return nothing
end

function _validate_runtime_contract(candidate_path::String, model_path::String)
    _verify_fixed_file(candidate_path, CANDIDATE_CONFIG_SHA256,
                       "structured candidate config")
    _verify_fixed_file(model_path, MODEL_CONFIG_SHA256,
                       "structured model config")
    model = TOML.parsefile(model_path)
    edge = model["model"]["edge"]
    edge["basis"] == ["corr_fwd", "corr_bwd"] ||
        _fail(:config_contract_mismatch, "edge basis differs")
    edge["student_t_nu"] == FIXED_NU ||
        _fail(:config_contract_mismatch, "edge nu differs")
    Tuple(edge["conditional_start_alphas"]) == START_ALPHAS ||
        _fail(:config_contract_mismatch, "edge starts differ")
    edge["mixed_state_mean"] == "tied_01_10" ||
        _fail(:config_contract_mismatch, "mixed-state mean differs")
    edge["scale"] == "shared_full_2x2" ||
        _fail(:config_contract_mismatch, "edge scale differs")
    edge["scale_shrinkage_formula"] == "lambda=min(0.5,3/N_edges)" ||
        _fail(:config_contract_mismatch, "shrinkage formula differs")
    edge["eigenvalue_floor"] == SCALE_FLOOR ||
        _fail(:config_contract_mismatch, "edge eigenvalue floor differs")
    edge["condition_number_cap"] == CONDITION_CAP ||
        _fail(:config_contract_mismatch, "edge condition cap differs")
    edge["pair_prior"] == [0.0, 0.0, 0.0, 0.0] ||
        _fail(:config_contract_mismatch, "edge pair prior is nonzero")
    edge["pair_prior_policy"] == "zero_data_independent" ||
        _fail(:config_contract_mismatch, "edge pair-prior policy differs")
    edge["responsibility_policy"] == "frozen_unary_soft_pairs_only" ||
        _fail(:config_contract_mismatch, "edge responsibilities are not frozen")
    edge["max_iter"] == MAX_ITER && edge["tol"] == CONVERGENCE_TOL &&
    edge["objective_decrease_tolerance"] == DECREASE_TOL &&
    edge["start_tie_tolerance"] == START_TIE_TOL ||
        _fail(:config_contract_mismatch, "edge iteration contract differs")
    preprocessing = model["preprocessing"]["residualization"]
    preprocessing["features"] == collect(PREDICTOR_NAMES) ||
        _fail(:config_contract_mismatch, "residualization features differ")
    preprocessing["endpoint_order"] ==
        "left_then_right_in_increasing_observed_t_nm" ||
        _fail(:config_contract_mismatch, "endpoint order differs")
    preprocessing["intercept_penalized"] === false ||
        _fail(:config_contract_mismatch, "residualizer intercept is penalized")
    preprocessing["predictor_count"] == 14 ||
        _fail(:config_contract_mismatch, "residualizer predictor count differs")
    selection = model["selection"]
    selection["bootstrap"]["count"] == 500 &&
    selection["bootstrap"]["edge_shuffle_count"] == 500 ||
        _fail(:config_contract_mismatch, "resampling count differs")
    feasibility = selection["feasibility"]
    feasibility["minimum_dates"] == MINIMUM_DATES &&
    feasibility["minimum_scans"] == MINIMUM_SCANS &&
    feasibility["minimum_category_kish_ess"] == MINIMUM_CATEGORY_ESS &&
    feasibility["minimum_distinct_edges"] == MINIMUM_EDGES &&
    feasibility["minimum_support_scans_per_training_date"] ==
        MINIMUM_SUPPORT_SCANS_PER_DATE ||
        _fail(:config_contract_mismatch, "edge feasibility differs")
    return nothing
end

function _load_universe(root::String, supplied::AbstractString)::UniverseInput
    directory = StructuredEdgeFeatures._resolve_existing_directory(
        root,
        supplied,
        "universe bundle",
    )
    expected_names = [
        "folds.tsv",
        "patch_keys.tsv",
        "perturbations.tsv",
        "receipt.toml",
        "scan_universe.tsv",
        "seeds.tsv",
        "shards.tsv",
    ]
    sort(readdir(directory)) == expected_names ||
        _fail(:universe_schema, "universe bundle file set is not exact")
    for name in expected_names
        path = joinpath(directory, name)
        isfile(path) && !islink(path) ||
            _fail(:universe_schema, "universe bundle member is not a regular file")
    end
    receipt_path = joinpath(directory, "receipt.toml")
    receipt_bytes = read(receipt_path)
    receipt = TOML.parse(String(copy(receipt_bytes)))
    _exact_keys(receipt, Set([
        "schema", "schema_version", "feature_sha256",
        "candidate_config_sha256", "model_config_sha256", "source_sha256",
        "forward_producer_sha256", "backward_producer_sha256", "grid_sha256",
        "universe_sha256", "keys_sha256", "artifacts",
    ]), "universe receipt")
    _required_string(receipt, "schema", "universe receipt") ==
        "stmfit-structured-scan-universe-receipt-v1" ||
        _fail(:universe_schema, "universe receipt schema differs")
    _required_integer(receipt, "schema_version", "universe receipt") == 1 ||
        _fail(:universe_schema, "universe receipt version differs")
    feature_sha = _require_sha(
        _required_string(receipt, "feature_sha256", "universe receipt"),
        "universe feature SHA-256",
    )
    _required_string(receipt, "candidate_config_sha256", "universe receipt") ==
        CANDIDATE_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "universe candidate hash differs")
    _required_string(receipt, "model_config_sha256", "universe receipt") ==
        MODEL_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "universe model hash differs")
    _required_string(receipt, "forward_producer_sha256", "universe receipt") ==
        StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256 ||
        _fail(:dependency_hash_mismatch, "universe forward producer differs")
    _required_string(receipt, "backward_producer_sha256", "universe receipt") ==
        StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256 ||
        _fail(:dependency_hash_mismatch, "universe backward producer differs")
    _required_string(receipt, "grid_sha256", "universe receipt") ==
        StructuredEdgeFeatures.GRID_SHA256 ||
        _fail(:grid_hash_mismatch, "universe grid differs")
    universe_sha = _require_sha(
        _required_string(receipt, "universe_sha256", "universe receipt"),
        "universe table SHA-256",
    )
    keys_sha = _require_sha(
        _required_string(receipt, "keys_sha256", "universe receipt"),
        "universe key SHA-256",
    )
    artifacts = receipt["artifacts"]
    artifacts isa AbstractDict ||
        _fail(:universe_schema, "universe artifacts table is absent")
    expected_artifacts = setdiff(expected_names, ["receipt.toml"])
    _exact_keys(artifacts, Set(expected_artifacts), "universe artifacts")
    for name in expected_artifacts
        expected = _require_sha(_required_string(artifacts, name,
                                                 "universe artifacts"),
                                "universe artifact hash")
        _sha256_file(joinpath(directory, name)) == expected ||
            _fail(:universe_hash_mismatch, "universe artifact bytes differ")
    end
    _sha256_file(joinpath(directory, "scan_universe.tsv")) == universe_sha ||
        _fail(:universe_hash_mismatch, "scan universe digest differs")

    scan_header, scan_rows = StructuredEdgeFeatures._parse_table(
        read(joinpath(directory, "scan_universe.tsv")),
        "scan universe",
    )
    expected_scan_header = [
        "schema_version", "file", "date", "lobe_count", "feature_sha256",
        "candidate_config_sha256", "model_config_sha256", "source_sha256",
        "forward_producer_sha256", "backward_producer_sha256", "grid_sha256",
    ]
    scan_header == expected_scan_header ||
        _fail(:universe_schema, "scan-universe header differs")
    scans = Dict{String,String}()
    scan_counts = Dict{String,Int}()
    for (row_index, fields) in enumerate(scan_rows)
        fields[1] == "1" || _fail(:universe_schema, "scan row version differs")
        file = fields[2]
        file == basename(file) || _fail(:invalid_file, "universe file is not a basename")
        date = StructuredEdgeFeatures.StructuredUniverse.parse_scan_date(file)
        fields[3] == date || _fail(:date_mismatch, "universe date differs")
        haskey(scans, file) && _fail(:duplicate_scan, "universe repeats a scan")
        count = _positive_integer_text(fields[4], "scan row $row_index lobe count")
        fields[5] == feature_sha && fields[6] == CANDIDATE_CONFIG_SHA256 &&
        fields[7] == MODEL_CONFIG_SHA256 &&
        fields[9] == StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256 &&
        fields[10] == StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256 &&
        fields[11] == StructuredEdgeFeatures.GRID_SHA256 ||
            _fail(:universe_hash_mismatch, "scan row binding differs")
        scans[file] = date
        scan_counts[file] = count
    end
    key_header, key_rows = StructuredEdgeFeatures._parse_table(
        read(joinpath(directory, "patch_keys.tsv")),
        "universe patch keys",
    )
    key_header == ["schema_version", "file", "lobe", "feature_sha256",
                   "universe_sha256", "keys_sha256"] ||
        _fail(:universe_schema, "patch-key header differs")
    keys = Tuple{String,Int}[]
    for (row_index, fields) in enumerate(key_rows)
        fields[1] == "1" || _fail(:universe_schema, "key row version differs")
        file = fields[2]
        haskey(scans, file) || _fail(:universe_key_mismatch, "key scan is absent")
        lobe = _positive_integer_text(fields[3], "key row $row_index lobe")
        fields[4] == feature_sha && fields[5] == universe_sha &&
        fields[6] == keys_sha ||
            _fail(:universe_hash_mismatch, "key row binding differs")
        push!(keys, (file, lobe))
    end
    length(keys) == length(unique(keys)) ||
        _fail(:duplicate_key, "universe repeats a node key")
    keys == sort(keys) || _fail(:universe_key_order, "universe keys are not sorted")
    for file in Base.keys(scans)
        count(key -> key[1] == file, keys) == scan_counts[file] ||
            _fail(:universe_key_mismatch, "scan lobe count differs from keys")
    end
    return UniverseInput(directory, _sha256(receipt_bytes), feature_sha,
                         universe_sha, keys_sha, scans, keys)
end

function _load_edge_bundle(
    root::String,
    supplied::AbstractString,
    universe::UniverseInput,
    forward_receipt::AbstractString,
    backward_receipt::AbstractString,
)::EdgeBundleInput
    directory = StructuredEdgeFeatures._resolve_existing_directory(
        root,
        supplied,
        "edge bundle",
    )
    sort(readdir(directory)) ==
        ["edge_observations.tsv", "node_segments.tsv", "receipt.toml"] ||
        _fail(:edge_bundle_schema, "edge bundle file set is not exact")
    for name in ("edge_observations.tsv", "node_segments.tsv", "receipt.toml")
        path = joinpath(directory, name)
        isfile(path) && !islink(path) ||
            _fail(:edge_bundle_schema, "edge bundle member is not a regular file")
    end
    receipt_bytes = read(joinpath(directory, "receipt.toml"))
    receipt = TOML.parse(String(copy(receipt_bytes)))
    receipt_keys = Set([
        "schema", "schema_version", "status", "edge_observations_file",
        "node_segments_file", "edge_header", "node_header",
        "edge_header_sha256", "node_header_sha256", "edge_observations_sha256",
        "node_segments_sha256", "edge_row_count", "node_row_count",
        "feature_sha256", "candidate_config_sha256", "model_config_sha256",
        "universe_receipt_sha256", "keys_sha256", "grid_sha256",
        "forward_patch_sha256", "backward_patch_sha256",
        "forward_patch_receipt_sha256", "backward_patch_receipt_sha256",
        "forward_producer_sha256", "backward_producer_sha256",
        "forward_command", "backward_command", "source_sha256",
        "t2_review_sha256", "t3_review_sha256",
    ])
    _exact_keys(receipt, receipt_keys, "edge receipt")
    _required_string(receipt, "schema", "edge receipt") ==
        "stmfit-structured-edge-feature-receipt-v1" ||
        _fail(:edge_bundle_schema, "edge receipt schema differs")
    _required_integer(receipt, "schema_version", "edge receipt") == 1 ||
        _fail(:edge_bundle_schema, "edge receipt version differs")
    _required_string(receipt, "status", "edge receipt") == "PASS" ||
        _fail(:edge_bundle_status, "edge receipt is not PASS")
    _required_string(receipt, "edge_observations_file", "edge receipt") ==
        "edge_observations.tsv" ||
        _fail(:edge_bundle_schema, "edge table filename differs")
    _required_string(receipt, "node_segments_file", "edge receipt") ==
        "node_segments.tsv" ||
        _fail(:edge_bundle_schema, "node table filename differs")
    hashes = _validated_edge_receipt_hashes(receipt)
    forward_command = _required_string(receipt, "forward_command", "edge receipt")
    backward_command = _required_string(receipt, "backward_command", "edge receipt")
    edge_path = joinpath(directory, "edge_observations.tsv")
    node_path = joinpath(directory, "node_segments.tsv")
    edge_sha = _sha256_file(edge_path)
    node_sha = _sha256_file(node_path)
    edge_sha == hashes["edge_observations_sha256"] ||
        _fail(:edge_hash_mismatch, "edge table bytes differ")
    node_sha == hashes["node_segments_sha256"] ||
        _fail(:edge_hash_mismatch, "node table bytes differ")
    feature_sha = hashes["feature_sha256"]
    feature_sha == universe.feature_sha256 ||
        _fail(:feature_hash_mismatch, "edge and universe feature hashes differ")
    hashes["candidate_config_sha256"] == CANDIDATE_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "edge candidate hash differs")
    hashes["model_config_sha256"] == MODEL_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "edge model hash differs")
    hashes["universe_receipt_sha256"] == universe.receipt_sha256 ||
        _fail(:universe_hash_mismatch, "edge universe receipt hash differs")
    hashes["keys_sha256"] == universe.keys_sha256 ||
        _fail(:universe_key_mismatch, "edge universe-key hash differs")
    hashes["grid_sha256"] == StructuredEdgeFeatures.GRID_SHA256 ||
        _fail(:grid_hash_mismatch, "edge grid hash differs")
    hashes["forward_producer_sha256"] ==
        StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256 ||
        _fail(:dependency_hash_mismatch, "edge forward producer differs")
    hashes["backward_producer_sha256"] ==
        StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256 ||
        _fail(:dependency_hash_mismatch, "edge backward producer differs")
    source_sha = hashes["source_sha256"]
    current_source_sha, _ = StructuredEdgeFeatures._source_snapshots()
    source_sha == current_source_sha ||
        _fail(:dependency_hash_mismatch, "edge source bundle differs")
    _validate_todo10_provenance(
        root,
        forward_receipt,
        backward_receipt,
        universe,
        hashes,
        forward_command,
        backward_command,
    )

    node_header, node_rows = StructuredEdgeFeatures._parse_table(
        read(node_path),
        "edge node table",
    )
    node_header == StructuredEdgeFeatures.NODE_HEADER ||
        _fail(:edge_bundle_schema, "node table header differs")
    _required_string(receipt, "node_header", "edge receipt") ==
        join(node_header, '\t') ||
        _fail(:edge_bundle_schema, "node receipt header differs")
    hashes["node_header_sha256"] ==
        _sha256(codeunits(join(node_header, '\t') * "\n")) ||
        _fail(:edge_bundle_schema, "node header digest differs")
    _required_integer(receipt, "node_row_count", "edge receipt") ==
        length(node_rows) || _fail(:edge_bundle_schema, "node row count differs")
    nodes = Dict{Tuple{String,Int},Tuple{Float64,String}}()
    ordered_node_keys = Tuple{String,Int}[]
    for (row_index, fields) in enumerate(node_rows)
        file = fields[1]
        lobe = _positive_integer_text(fields[2], "node row $row_index lobe")
        t_nm = _finite_text(fields[3], "node row $row_index t")
        segment = fields[4]
        occursin(r"^[0-9a-f]{16}$", segment) ||
            _fail(:segment_mismatch, "node segment identifier is invalid")
        fields[5] in ("connected", "isolated") ||
            _fail(:edge_bundle_schema, "node status is invalid")
        node_feature_sha = _require_sha(fields[8],
                                        "node row $row_index feature hash")
        node_model_sha = _require_sha(fields[9],
                                      "node row $row_index model hash")
        node_source_sha = _require_sha(fields[10],
                                       "node row $row_index source hash")
        node_feature_sha == feature_sha && node_model_sha == MODEL_CONFIG_SHA256 &&
        node_source_sha == source_sha ||
            _fail(:edge_hash_mismatch, "node row binding differs")
        key = (file, lobe)
        haskey(nodes, key) && _fail(:duplicate_key, "node table repeats a key")
        nodes[key] = (t_nm, segment)
        push!(ordered_node_keys, key)
    end
    Set(ordered_node_keys) == Set(universe.keys) ||
        _fail(:universe_key_mismatch, "node keys differ from universe")
    ordered_node_keys == sort(ordered_node_keys; by=key ->
        (key[1], nodes[key][1], key[2])) ||
        _fail(:edge_order, "node rows are not in observed-t order")

    edge_bytes = read(edge_path)
    edge_row_bytes = _exact_table_row_bytes(edge_bytes, "edge observations")
    edge_header, edge_rows = StructuredEdgeFeatures._parse_table(
        edge_bytes,
        "edge observations",
    )
    length(edge_row_bytes) == length(edge_rows) ||
        _fail(:edge_bundle_schema, "edge row byte count differs")
    edge_header == StructuredEdgeFeatures.EDGE_HEADER ||
        _fail(:edge_bundle_schema, "edge table header differs")
    _required_string(receipt, "edge_header", "edge receipt") ==
        join(edge_header, '\t') ||
        _fail(:edge_bundle_schema, "edge receipt header differs")
    hashes["edge_header_sha256"] ==
        _sha256(codeunits(join(edge_header, '\t') * "\n")) ||
        _fail(:edge_bundle_schema, "edge header digest differs")
    _required_integer(receipt, "edge_row_count", "edge receipt") ==
        length(edge_rows) || _fail(:edge_bundle_schema, "edge row count differs")
    expected_pairs = Tuple{String,Int,Int}[]
    for file in sort!(collect(keys(universe.scans)))
        chain = sort!([key for key in universe.keys if key[1] == file];
                      by=key -> (nodes[key][1], key[2]))
        for index in 1:max(0, length(chain) - 1)
            push!(expected_pairs, (file, chain[index][2], chain[index + 1][2]))
        end
    end
    observed_pairs = Tuple{String,Int,Int}[]
    eligible_edges = AdmissionEdge[]
    ordinal_by_segment = Dict{Tuple{String,String},Int}()
    for (row_index, fields) in enumerate(edge_rows)
        file = fields[1]
        left = _positive_integer_text(fields[2], "edge row $row_index left lobe")
        right = _positive_integer_text(fields[3], "edge row $row_index right lobe")
        left_t = _finite_text(fields[4], "edge row $row_index left t")
        right_t = _finite_text(fields[5], "edge row $row_index right t")
        _finite_text(fields[6], "edge row $row_index gap") == right_t - left_t ||
            _fail(:edge_metadata_mismatch, "edge gap differs from endpoints")
        nodes[(file, left)][1] == left_t && nodes[(file, right)][1] == right_t ||
            _fail(:edge_metadata_mismatch, "edge t differs from node table")
        status = fields[9]
        reason = fields[10]
        status in ("eligible", "split") ||
            _fail(:edge_bundle_schema, "edge status is invalid")
        row_hashes = [
            _require_sha(fields[column],
                         "edge row $row_index $(edge_header[column])")
            for column in 13:18
        ]
        row_hashes[1] == StructuredEdgeFeatures.GRID_SHA256 &&
        row_hashes[2] == hashes["forward_patch_sha256"] &&
        row_hashes[3] == hashes["backward_patch_sha256"] &&
        row_hashes[4] == feature_sha && row_hashes[5] == MODEL_CONFIG_SHA256 &&
        row_hashes[6] == source_sha ||
            _fail(:edge_hash_mismatch, "edge row binding differs")
        push!(observed_pairs, (file, left, right))
        if status == "eligible"
            reason == "none" ||
                _fail(:edge_bundle_schema, "eligible edge has a split reason")
            fields[7] == fields[8] ||
                _fail(:segment_mismatch, "eligible edge crosses a segment")
            forward = _finite_text(fields[11], "edge row $row_index forward")
            backward = _finite_text(fields[12], "edge row $row_index backward")
            -1.0 <= forward <= 1.0 && -1.0 <= backward <= 1.0 ||
                _fail(:edge_observation_range, "edge observation is outside [-1,1]")
            segment_key = (file, fields[7])
            ordinal = get(ordinal_by_segment, segment_key, 0) + 1
            ordinal_by_segment[segment_key] = ordinal
            date = universe.scans[file]
            push!(eligible_edges, AdmissionEdge(file, date, left, right, left_t,
                                                right_t, fields[7], ordinal,
                                                (forward, backward),
                                                _sha256(edge_row_bytes[row_index])))
        else
            reason != "none" ||
                _fail(:edge_bundle_schema, "split edge lacks a reason")
            fields[11] == "NA" && fields[12] == "NA" ||
                _fail(:edge_bundle_schema, "split edge exposes an observation")
            fields[7] != fields[8] ||
                _fail(:segment_mismatch, "split edge does not split segments")
        end
    end
    observed_pairs == expected_pairs ||
        _fail(:edge_key_mismatch, "edge rows do not cover consecutive nodes")
    return EdgeBundleInput(directory, _sha256(receipt_bytes), edge_sha,
                           node_sha, feature_sha, source_sha,
                           hashes["t3_review_sha256"], nodes,
                           eligible_edges)
end

function _load_admission_nodes(
    root::String,
    supplied::AbstractString,
    universe::UniverseInput,
    edge::EdgeBundleInput,
)
    path = StructuredEdgeFeatures._resolve_existing(root, supplied,
                                                    "admission feature table")
    bytes, _ = StructuredEdgeFeatures._snapshot(path, "admission feature table")
    header, rows = StructuredEdgeFeatures._parse_table(bytes,
                                                       "admission feature table")
    header == ADMISSION_FEATURE_HEADER ||
        _fail(:admission_feature_schema,
              "admission feature header is not the exact allowlist")
    firewall = StructuredEdgeFeatures.StructuredUniverse.InputBoundary
    firewall.validate_declared_features(vcat(["amplitude"],
                                             collect(PREDICTOR_NAMES)))
    nodes = AdmissionNode[]
    keys = Tuple{String,Int}[]
    for (row_index, fields) in enumerate(rows)
        file = fields[1]
        file == basename(file) || _fail(:invalid_file, "admission file is not a basename")
        date = StructuredEdgeFeatures.StructuredUniverse.parse_scan_date(file)
        lobe = _positive_integer_text(fields[2], "admission row $row_index lobe")
        t_nm = _finite_text(fields[3], "admission row $row_index t")
        amplitude = _finite_text(fields[4], "admission row $row_index amplitude")
        predictors = ntuple(7) do index
            _optional_text(fields[4 + index],
                           "admission row $row_index $(PREDICTOR_NAMES[index])")
        end
        key = (file, lobe)
        haskey(edge.nodes, key) ||
            _fail(:admission_key_mismatch, "admission key is absent from edge nodes")
        edge.nodes[key][1] == t_nm ||
            _fail(:edge_metadata_mismatch, "admission t differs from edge node")
        universe.scans[file] == date ||
            _fail(:date_mismatch, "admission date differs from universe")
        push!(keys, key)
        push!(nodes, AdmissionNode(file, date, lobe, t_nm, amplitude,
                                  predictors))
    end
    length(keys) == length(unique(keys)) ||
        _fail(:duplicate_key, "admission feature table repeats a key")
    Set(keys) == Set(universe.keys) ||
        _fail(:admission_key_mismatch, "admission keys differ from universe")
    nodes == sort(nodes; by=node -> (node.file, node.t_nm, node.lobe)) ||
        _fail(:node_order, "admission feature rows are not in observed-t order")
    return nodes, _sha256(bytes), path
end

function source_bundle_sha256(root::AbstractString=_repository_root())::String
    canonical_root = StructuredEdgeFeatures._canonical_root(root)
    paths = String[
        joinpath(canonical_root, "test", "evaluate_structured_edge_admission.jl"),
        joinpath(canonical_root, "test", "lib", "structured_assignment",
                 "edge_admission.jl"),
        joinpath(canonical_root, "test", "lib", "hierarchical_unit_assignment.jl"),
        joinpath(canonical_root, "test", "lib", "structured_assignment",
                 "robust_emissions.jl"),
        joinpath(canonical_root, "test", "lib", "structured_assignment",
                 "edge_features.jl"),
        joinpath(canonical_root, "test", "lib", "structured_assignment",
                 "universe.jl"),
        joinpath(canonical_root, "test", "lib", "structured_assignment",
                 "firewall.jl"),
    ]
    append!(paths, sort!(filter(isfile, readdir(
        joinpath(canonical_root, "test", "lib", "hierarchical"); join=true))))
    length(paths) == length(unique(paths)) ||
        _fail(:source_collision, "admission source path is repeated")
    lines = String["structured-edge-admission-source-v1"]
    for path in sort!(paths; by=path -> relpath(path, canonical_root))
        isfile(path) && !islink(path) ||
            _fail(:missing_dependency, "admission source is absent")
        push!(lines, string(relpath(path, canonical_root), '=', _sha256_file(path)))
    end
    return _hash_lines(lines)
end

function load_admission_data(
    root::AbstractString;
    features::AbstractString,
    candidate_config::AbstractString,
    model_config::AbstractString,
    universe_dir::AbstractString,
    edge_dir::AbstractString,
    forward_receipt::AbstractString,
    backward_receipt::AbstractString,
)::AdmissionData
    canonical_root = StructuredEdgeFeatures._canonical_root(root)
    _verify_dependencies(canonical_root)
    candidate_path = StructuredEdgeFeatures._resolve_existing(
        canonical_root,
        candidate_config,
        "candidate config",
    )
    model_path = StructuredEdgeFeatures._resolve_existing(
        canonical_root,
        model_config,
        "model config",
    )
    _validate_runtime_contract(candidate_path, model_path)
    universe = _load_universe(canonical_root, universe_dir)
    edge = _load_edge_bundle(
        canonical_root,
        edge_dir,
        universe,
        forward_receipt,
        backward_receipt,
    )
    nodes, input_sha, input_path = _load_admission_nodes(
        canonical_root,
        features,
        universe,
        edge,
    )
    hashes = Dict(
        "plan_sha256" => PLAN_SHA256,
        "t7_review_sha256" => T7_REVIEW_SHA256,
        "t10_review_sha256" => T10_REVIEW_SHA256,
        "candidate_config_sha256" => CANDIDATE_CONFIG_SHA256,
        "model_config_sha256" => MODEL_CONFIG_SHA256,
        "universe_receipt_sha256" => universe.receipt_sha256,
        "universe_sha256" => universe.universe_sha256,
        "universe_keys_sha256" => universe.keys_sha256,
        "edge_receipt_sha256" => edge.receipt_sha256,
        "edge_observations_sha256" => edge.edge_sha256,
        "node_segments_sha256" => edge.node_sha256,
        "edge_source_sha256" => edge.source_sha256,
        "edge_feature_sha256" => edge.feature_sha256,
        "t3_review_sha256" => edge.t3_review_sha256,
        "input_sha256" => input_sha,
        "source_sha256" => source_bundle_sha256(canonical_root),
    )
    _ = input_path
    return AdmissionData(nodes, edge.edges, hashes)
end

const MODEL_RESULT_HEADER = [
    "schema_version", "model", "status", "reason", "outer_fold_count",
    "final_gate_sha256", "full_fit_sha256", "result_sha256",
]
const PARTITION_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date",
    "status", "reason", "training_dates", "training_scan_count",
    "heldout_scan_count", "training_sha256", "fit_sha256", "score_sha256",
    "raw_heldout_sha256", "raw_at_score_sha256", "reversal_pass",
    "shuffle_pass",
]
const ESS_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date",
    "category", "dates", "scans", "edges", "kish_ess",
    "minimum_support_scans", "free_parameters", "effective_observation_ratio",
    "sufficient", "reason", "fit_sha256",
]
const START_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date",
    "fit_role", "start_index", "alpha", "status", "reason", "iterations",
    "objective_initial", "objective_final", "rejected_iteration",
    "rejected_objective", "trace_sha256", "fit_sha256",
]
const SCORE_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date",
    "level", "file", "left_lobe", "right_lobe", "segment_id", "gain",
    "edge_count", "score_sha256",
]
const SHUFFLE_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date",
    "seed", "date_gain", "conditional_fit_sha256", "observed_control_gain",
    "upper_quantile", "passed", "score_sha256",
]
const BOOTSTRAP_HEADER = [
    "schema_version", "model", "scope", "outer_date", "seed", "mean_gain",
    "lower_quantile", "every_date_positive", "every_shuffle_pass",
    "reversal_pass", "status", "reason", "gate_sha256",
]
const GRAPH_MODEL_HEADER = vcat([
    "schema_version", "fit_role", "fit_sha256", "model", "scope",
    "outer_date", "target_date", "training_dates", "reversed", "status",
    "reason", "training_sha256", "unary_fit_sha256",
    "reference_count", "references_sha256", "references",
    "unary_probability_sha256", "node_order_sha256", "residualizer_sha256",
    "residualizer_training_sha256", "residualizer_ridge", "student_t_nu",
    "category_order", "state_order", "null_mean_fwd", "null_mean_bwd",
    "null_scale_ff", "null_scale_fb", "null_scale_bb",
], [
    "conditional_$(category)_mean_$(output)"
    for category in GRAPH_CATEGORY_ORDER for output in GRAPH_OUTPUT_NAMES
], [
    "conditional_scale_ff", "conditional_scale_fb", "conditional_scale_bb",
    "selected_start", "selected_alpha", "selected_objective",
    "selected_converged",
], [
    "$(name)_$(output)"
    for name in GRAPH_COEFFICIENT_NAMES for output in GRAPH_OUTPUT_NAMES
])
const GRAPH_UNARY_HEADER = [
    "schema_version", "fit_sha256", "model", "file", "date", "lobe",
    "t_nm", "node_identity", "state_order", "p0", "p1",
    "unary_fit_sha256", "unary_probability_sha256", "node_order_sha256",
]
const GRAPH_EDGE_HEADER = [
    "schema_version", "fit_sha256", "model", "scope", "outer_date",
    "target_date", "file", "date", "left_lobe", "right_lobe",
    "left_t_nm", "right_t_nm", "segment_id", "ordinal", "edge_identity",
    "edge_row_raw_sha256", "residualizer_sha256", "unary_fit_sha256",
    "pred_fwd", "pred_bwd", "status", "reason",
]
const GRAPH_HANDOFF_FILES = [
    "fitted_edge_models.tsv",
    "fitted_unary_nodes.tsv",
    "fitted_edge_transforms.tsv",
]
const REPORT_ARTIFACT_FILES = [
    "models.tsv", "partitions.tsv", "ess.tsv", "starts.tsv", "scores.tsv",
    "shuffle.tsv", "bootstrap.tsv",
]

function _table_bytes(header::Vector{String}, rows)
    io = IOBuffer()
    println(io, join(header, '\t'))
    for row in rows
        length(row) == length(header) ||
            _fail(:internal_contract, "report row width differs")
        any(value -> occursin('\t', value) || occursin('\n', value) ||
                     occursin('\r', value), row) &&
            _fail(:internal_contract, "report field contains a delimiter")
        println(io, join(row, '\t'))
    end
    return take!(io)
end

_out_float(value::Float64) = isfinite(value) ? _format_float(value) : "NA"
_out_bool(value::Bool) = value ? "true" : "false"

function _public_scale_fields(scale::AbstractMatrix)
    size(scale) == (2, 2) || return nothing
    scale[1, 2] == scale[2, 1] || return nothing
    all(isfinite, scale) || return nothing
    canonical = _canonical_scale_matrix(scale)
    canonical === nothing && return nothing
    stored, fields = canonical
    _stored_scale_certificate(stored).status == :ok || return nothing
    return fields
end

function _logical_evaluations(result::ModelAdmissionResult)
    logical = NamedTuple[]
    for outer in result.outer_folds
        for partition in outer.inner_gate.partitions
            push!(logical, (scope="outer_inner", outer_date=outer.outer_date,
                            evaluation=partition))
        end
        push!(logical, (scope="outer_score", outer_date=outer.outer_date,
                        evaluation=outer.outer_score))
    end
    for partition in result.final_gate.partitions
        push!(logical, (scope="final_lodo", outer_date="NA",
                        evaluation=partition))
    end
    return logical
end

function _raw_heldout_hash(evaluation::PartitionEvaluation)::String
    lines = String[]
    for sample in evaluation.fit.samples
        sample.edge.date == evaluation.target_date || continue
        push!(lines, string(_edge_identity(sample.edge), '=', sample.raw_sha256))
    end
    return _hash_lines(lines)
end

function _trace_hash(trace::Vector{Float64})::String
    return _hash_lines(_format_float.(trace))
end

function _append_fit_evidence!(
    ess_rows::Vector{Vector{String}},
    start_rows::Vector{Vector{String}},
    model::String,
    scope::String,
    outer_date::String,
    heldout_date::String,
    fit::PartitionFitResult,
)
    for category in 1:4
        push!(ess_rows, [
            "1", model, scope, outer_date, heldout_date, string(category - 1),
            string(fit.support.dates), string(fit.support.scans),
            string(fit.support.edges), _out_float(fit.support.category_ess[category]),
            string(fit.support.minimum_support_scans[category]),
            string(fit.support.free_parameters),
            _out_float(fit.support.effective_observation_ratio),
            _out_bool(fit.support.sufficient), String(fit.support.reason),
            fit.fit_sha256,
        ])
    end
    null_fit = fit.null_fit
    if null_fit !== nothing
        fitted_null = null_fit::NullFit
        push!(start_rows, [
            "1", model, scope, outer_date, heldout_date, "null", "0", "NA",
            fitted_null.status == :ok ? "PASS" : "FAIL",
            String(fitted_null.reason), string(fitted_null.iterations),
            isempty(fitted_null.objective_trace) ? "NA" :
                _out_float(first(fitted_null.objective_trace)),
            isempty(fitted_null.objective_trace) ? "NA" :
                _out_float(last(fitted_null.objective_trace)),
            "0", _out_float(fitted_null.rejected_objective),
            _trace_hash(fitted_null.objective_trace), fit.fit_sha256,
        ])
    end
    conditional = fit.conditional_fit
    if conditional !== nothing
        fitted = conditional::ConditionalFit
        for start in fitted.starts
            push!(start_rows, [
                "1", model, scope, outer_date, heldout_date, "conditional",
                string(start.start_index), _out_float(start.alpha),
                start.valid ? "PASS" : "FAIL", String(start.reason),
                string(max(0, length(start.objective_trace) - 1)),
                isempty(start.objective_trace) ? "NA" :
                    _out_float(first(start.objective_trace)),
                isempty(start.objective_trace) ? "NA" :
                    _out_float(last(start.objective_trace)),
                string(start.rejected_iteration),
                _out_float(start.rejected_objective),
                _trace_hash(start.objective_trace), fit.fit_sha256,
            ])
        end
    end
    return nothing
end

function _graph_fit_references(report::AdmissionReport)
    references = NamedTuple[]
    for model in sort!(collect(keys(report.models)))
        result = report.models[model]
        for item in _logical_evaluations(result)
            evaluation = item.evaluation
            push!(references, (
                fit=evaluation.fit,
                fit_role="partition",
                model=model,
                scope=item.scope,
                outer_date=item.outer_date,
                target_date=evaluation.target_date,
            ))
        end
        push!(references, (
            fit=result.full_fit,
            fit_role="full_refit",
            model=model,
            scope="full_refit",
            outer_date="NA",
            target_date="NA",
        ))
    end
    grouped = Dict{String,Vector{NamedTuple}}()
    for reference in references
        push!(get!(grouped, reference.fit.fit_sha256, NamedTuple[]), reference)
    end
    result = NamedTuple[]
    for fit_sha256 in sort!(collect(keys(grouped)))
        entries = sort!(grouped[fit_sha256]; by=reference -> (
            reference.model,
            reference.fit_role == "partition" ? 0 : 1,
            reference.scope,
            reference.outer_date,
            reference.target_date,
        ))
        canonical = first(entries)
        reference_lines = [join((
            entry.fit_role,
            entry.model,
            entry.scope,
            entry.outer_date,
            entry.target_date,
            entry.fit.fit_sha256,
        ), '\t') for entry in entries]
        references = [join((
            entry.fit_role,
            entry.model,
            entry.scope,
            entry.outer_date,
            entry.target_date,
            entry.fit.fit_sha256,
        ), '|') for entry in entries]
        push!(result, merge(canonical, (
            reference_count=length(entries),
            references_sha256=_hash_lines(reference_lines),
            references_serialized=join(references, ';'),
        )))
    end
    return result
end

function _residualizer_sha256(fit::ResidualizerFit)::String
    lines = String[
        "schema=structured-edge-residualizer-v1",
        "training_sha256=$(fit.training_sha256)",
        "training_row_count=$(fit.training_row_count)",
        "ridge=$(_format_float(fit.ridge))",
        "evaluation_order=$(join(fit.evaluation_order, ','))",
    ]
    for (name, row) in zip(GRAPH_COEFFICIENT_NAMES, axes(fit.coefficients, 1))
        for (output, column) in zip(GRAPH_OUTPUT_NAMES, axes(fit.coefficients, 2))
            push!(lines, "$(name)_$(output)=$(_format_float(fit.coefficients[row, column]))")
        end
    end
    return _hash_lines(lines)
end

function _graph_model_row(reference)::Vector{String}
    fit = reference.fit
    row = fill("NA", length(GRAPH_MODEL_HEADER))
    positions = Dict(name => index for (index, name) in enumerate(GRAPH_MODEL_HEADER))
    put(name, value) = (row[positions[name]] = string(value))
    put("schema_version", "1")
    put("fit_role", reference.fit_role)
    put("fit_sha256", fit.fit_sha256)
    put("model", reference.model)
    put("scope", reference.scope)
    put("outer_date", reference.outer_date)
    put("target_date", reference.target_date)
    put("training_dates", join(fit.training_dates, ','))
    put("reversed", _out_bool(fit.reversed))
    put("status", fit.status)
    put("reason", String(fit.reason))
    put("training_sha256", fit.training_sha256)
    put("unary_fit_sha256", fit.unary.fit_sha256)
    put("reference_count", reference.reference_count)
    put("references_sha256", reference.references_sha256)
    put("references", reference.references_serialized)
    put("unary_probability_sha256", fit.unary.probability_sha256)
    put("node_order_sha256", fit.unary.node_order_sha256)
    put("student_t_nu", "8")
    put("category_order", join(GRAPH_CATEGORY_ORDER, ','))
    put("state_order", join(GRAPH_STATE_ORDER, ','))

    residualizer = fit.residualizer
    if residualizer !== nothing
        fitted_residualizer = residualizer::ResidualizerFit
        put("residualizer_sha256", _residualizer_sha256(fitted_residualizer))
        put("residualizer_training_sha256", fitted_residualizer.training_sha256)
        put("residualizer_ridge", _format_float(fitted_residualizer.ridge))
        for (name, row_index) in zip(GRAPH_COEFFICIENT_NAMES,
                                     axes(fitted_residualizer.coefficients, 1))
            for (output, column_index) in zip(GRAPH_OUTPUT_NAMES,
                                              axes(fitted_residualizer.coefficients, 2))
                put("$(name)_$(output)",
                    _format_float(fitted_residualizer.coefficients[row_index, column_index]))
            end
        end
    end

    null_fit = fit.null_fit
    if null_fit !== nothing
        fitted_null = null_fit::NullFit
        put("null_mean_fwd", _out_float(fitted_null.mean[1]))
        put("null_mean_bwd", _out_float(fitted_null.mean[2]))
        null_fields = _public_scale_fields(fitted_null.scale)
        if fitted_null.status == :ok && null_fields === nothing
            _fail(:graph_handoff_schema, "PASS null scale is not certified")
        elseif null_fields !== nothing
            put("null_scale_ff", null_fields[1])
            put("null_scale_fb", null_fields[2])
            put("null_scale_bb", null_fields[3])
        end
    end

    conditional = fit.conditional_fit
    if conditional !== nothing
        fitted_conditional = conditional::ConditionalFit
        fitted_conditional.means[2, :] == fitted_conditional.means[3, :] ||
            _fail(:graph_handoff_schema, "conditional 01/10 means are not exactly tied")
        for (category_index, category) in enumerate(GRAPH_CATEGORY_ORDER)
            for (output_index, output) in enumerate(GRAPH_OUTPUT_NAMES)
                put("conditional_$(category)_mean_$(output)",
                    _out_float(fitted_conditional.means[category_index, output_index]))
            end
        end
        conditional_fields = _public_scale_fields(fitted_conditional.scale)
        if fitted_conditional.status == :ok && conditional_fields === nothing
            _fail(:graph_handoff_schema, "PASS conditional scale is not certified")
        elseif conditional_fields !== nothing
            put("conditional_scale_ff", conditional_fields[1])
            put("conditional_scale_fb", conditional_fields[2])
            put("conditional_scale_bb", conditional_fields[3])
        end
        if fitted_conditional.status == :ok
            selected = fitted_conditional.starts[fitted_conditional.best_start]
            put("selected_start", fitted_conditional.best_start)
            put("selected_alpha", _format_float(selected.alpha))
            put("selected_objective", _format_float(fitted_conditional.objective))
            put("selected_converged", _out_bool(fitted_conditional.converged))
        end
    end
    return row
end

function _graph_handoff_artifacts(report::AdmissionReport)
    references = _graph_fit_references(report)
    model_rows = [_graph_model_row(reference) for reference in references]
    unary_rows = NamedTuple[]
    edge_rows = NamedTuple[]
    for reference in references
        fit = reference.fit
        unary = fit.unary
        length(unary.node_identities) == length(unary.probabilities_1) ||
            _fail(:graph_handoff_schema, "unary node identity count differs")
        unary.node_order_sha256 == _node_order_sha256(unary.node_identities) ||
            _fail(:graph_handoff_schema, "unary node order binding differs")
        for index in eachindex(unary.probabilities_1)
            identity = unary.node_identities[index]
            p1 = unary.probabilities_1[index]
            isfinite(p1) && 0.0 <= p1 <= 1.0 ||
                _fail(:graph_handoff_schema, "unary probability is invalid")
            push!(unary_rows, (
                fit_sha256=fit.fit_sha256,
                model=reference.model,
                file=identity[1],
                date=identity[2],
                lobe=identity[3],
                t_nm=identity[4],
                node_identity=_graph_node_identity(identity),
                p0=1.0 - p1,
                p1=p1,
                unary_fit_sha256=unary.fit_sha256,
                unary_probability_sha256=unary.probability_sha256,
                node_order_sha256=unary.node_order_sha256,
            ))
        end
        residualizer = fit.residualizer
        residualizer_sha = residualizer === nothing ? "NA" :
                           _residualizer_sha256(residualizer::ResidualizerFit)
        transform_ids = [_graph_edge_identity(state.edge) for state in fit.transform_edges]
        length(transform_ids) == length(unique(transform_ids)) ||
            _fail(:graph_handoff_schema, "transform edge coverage is duplicated")
        for state in fit.transform_edges
            state.status in ("eligible", "ineligible", "unavailable") ||
                _fail(:graph_handoff_schema, "transform edge status is invalid")
            expected_reason = state.status == "eligible" ? :none :
                              state.status == "ineligible" ?
                              :nonfinite_endpoint_predictor :
                              :residualizer_unavailable
            state.reason == expected_reason ||
                _fail(:graph_handoff_schema, "transform edge reason is invalid")
            prediction = if state.status == "eligible"
                residualizer === nothing &&
                    _fail(:graph_handoff_schema, "eligible transform lacks residualizer")
                _residualizer_prediction(residualizer::ResidualizerFit,
                                         state.endpoint_predictors)
            else
                nothing
            end
            push!(edge_rows, (
                fit_sha256=fit.fit_sha256,
                model=reference.model,
                scope=reference.scope,
                outer_date=reference.outer_date,
                target_date=reference.target_date,
                file=state.edge.file,
                date=state.edge.date,
                left_lobe=state.edge.left_lobe,
                right_lobe=state.edge.right_lobe,
                left_t_nm=state.edge.left_t_nm,
                right_t_nm=state.edge.right_t_nm,
                segment_id=state.edge.segment_id,
                ordinal=state.edge.ordinal,
                edge_identity=_graph_edge_identity(state.edge),
                edge_row_raw_sha256=state.edge.raw_row_sha256,
                residualizer_sha256=residualizer_sha,
                unary_fit_sha256=unary.fit_sha256,
                pred_fwd=prediction === nothing ? "NA" : _format_float(prediction[1]),
                pred_bwd=prediction === nothing ? "NA" : _format_float(prediction[2]),
                status=state.status,
                reason=String(state.reason),
            ))
        end
    end
    sort!(unary_rows; by=row -> (row.fit_sha256, row.file, row.t_nm, row.lobe))
    sort!(edge_rows; by=row -> (
        row.fit_sha256, row.file, row.left_t_nm, row.right_t_nm,
        row.left_lobe, row.right_lobe,
    ))
    return Dict(
        "fitted_edge_models.tsv" => _table_bytes(GRAPH_MODEL_HEADER, model_rows),
        "fitted_unary_nodes.tsv" => _table_bytes(
            GRAPH_UNARY_HEADER,
            [[
                "1", row.fit_sha256, row.model, row.file, row.date,
                string(row.lobe), _format_float(row.t_nm), row.node_identity,
                join(GRAPH_STATE_ORDER, ','), _format_float(row.p0),
                _format_float(row.p1), row.unary_fit_sha256,
                row.unary_probability_sha256, row.node_order_sha256,
            ] for row in unary_rows],
        ),
        "fitted_edge_transforms.tsv" => _table_bytes(
            GRAPH_EDGE_HEADER,
            [[
                "1", row.fit_sha256, row.model, row.scope, row.outer_date,
                row.target_date, row.file, row.date, string(row.left_lobe),
                string(row.right_lobe), _format_float(row.left_t_nm),
                _format_float(row.right_t_nm), row.segment_id, string(row.ordinal),
                row.edge_identity, row.edge_row_raw_sha256,
                row.residualizer_sha256, row.unary_fit_sha256,
                row.pred_fwd, row.pred_bwd,
                row.status, row.reason,
            ] for row in edge_rows],
        ),
    )
end

function graph_handoff_replay_hash(files::AbstractDict)::String
    Set(GRAPH_HANDOFF_FILES) ⊆ Set(String.(collect(keys(files)))) ||
        _fail(:graph_handoff_schema, "graph-handoff replay inputs are incomplete")
    io = IOBuffer()
    write(io, codeunits("schema=structured-edge-graph-handoff-replay-v1\n"))
    for name in GRAPH_HANDOFF_FILES
        bytes = Vector{UInt8}(files[name])
        write(io, codeunits("file=$(name)\nlength=$(length(bytes))\n"))
        write(io, bytes)
        write(io, UInt8('\n'))
    end
    return _sha256(take!(io))
end

function _parse_report_table(bytes::Vector{UInt8}, header::Vector{String}, name::String)
    text = String(copy(bytes))
    isvalid(text) || _fail(:invalid_encoding, "$name is not valid UTF-8")
    occursin('\r', text) && _fail(:carriage_return, "$name must use LF line endings")
    endswith(text, '\n') || _fail(:missing_final_lf, "$name lacks final LF")
    lines = split(chomp(text), '\n'; keepempty=true)
    observed, rows = if length(lines) == 1
        String.(split(first(lines), '\t'; keepempty=true)), Vector{Vector{String}}()
    else
        StructuredEdgeFeatures._parse_table(bytes, name)
    end
    observed == header || _fail(:graph_handoff_schema, "$name header differs")
    return rows
end

function _validate_graph_handoff_artifacts(
    report::AdmissionReport,
    artifacts::Dict{String,Vector{UInt8}},
)
    expected = _graph_handoff_artifacts(report)
    Set(keys(artifacts)) == Set(vcat(REPORT_ARTIFACT_FILES, GRAPH_HANDOFF_FILES)) ||
        _fail(:graph_handoff_schema, "publication file set differs")
    Set(name for name in keys(artifacts) if name in GRAPH_HANDOFF_FILES) ==
        Set(GRAPH_HANDOFF_FILES) ||
        _fail(:graph_handoff_schema, "graph-handoff file set differs")
    for name in GRAPH_HANDOFF_FILES
        haskey(artifacts, name) ||
            _fail(:graph_handoff_schema, "graph-handoff artifact is absent")
        artifacts[name] == expected[name] ||
            _fail(:graph_handoff_schema, "$name bytes differ from in-memory values")
    end
    _parse_report_table(artifacts["fitted_edge_models.tsv"], GRAPH_MODEL_HEADER,
                        "fitted edge models")
    _parse_report_table(artifacts["fitted_unary_nodes.tsv"], GRAPH_UNARY_HEADER,
                        "fitted unary nodes")
    _parse_report_table(artifacts["fitted_edge_transforms.tsv"], GRAPH_EDGE_HEADER,
                        "fitted edge transforms")
    return nothing
end

function _table_data_row_count(bytes::Vector{UInt8}, name::String)::Int
    text = String(copy(bytes))
    endswith(text, '\n') || _fail(:graph_handoff_schema, "$name lacks final LF")
    return max(0, length(split(chomp(text), '\n')) - 1)
end

function _graph_fit_hashes(report::AdmissionReport)
    references = _graph_fit_references(report)
    return sort!(unique(reference.fit.fit_sha256 for reference in references)),
           sort!(unique(result.full_fit.fit_sha256 for result in values(report.models)))
end

function _receipt_bytes(
    report::AdmissionReport,
    artifacts::Dict{String,Vector{UInt8}},
)::Vector{UInt8}
    fit_hashes, full_fit_hashes = _graph_fit_hashes(report)
    all_tsv = sort!(collect(keys(artifacts)))
    row_counts = Dict(name => _table_data_row_count(artifacts[name], name)
                      for name in all_tsv)
    io = IOBuffer()
    println(io, "schema = ", repr("stmfit-structured-edge-admission-receipt-v2"))
    println(io, "schema_version = 2")
    println(io, "status = ", repr(report.status))
    println(io, "reason = ", repr(String(report.reason)))
    println(io, "result_sha256 = ", repr(report.result_sha256))
    println(io, "model_count = ", length(report.models))
    println(io, "bootstrap_count = 500")
    println(io, "shuffle_count = 500")
    println(io, "fixed_nu = 8")
    println(io, "node_priors = [0.5, 0.5]")
    println(io, "pair_prior = [0.0, 0.0, 0.0, 0.0]")
    println(io, "graph_handoff_schema = ", repr(GRAPH_HANDOFF_SCHEMA))
    println(io, "graph_handoff_replay_sha256 = ",
            repr(graph_handoff_replay_hash(artifacts)))
    println(io, "graph_handoff_files = [",
            join(repr.(GRAPH_HANDOFF_FILES), ", "), "]")
    println(io, "category_order = [", join(repr.(GRAPH_CATEGORY_ORDER), ", "), "]")
    println(io, "state_order = [", join(repr.(GRAPH_STATE_ORDER), ", "), "]")
    coefficient_order = ["$(name)_$(output)"
                         for name in GRAPH_COEFFICIENT_NAMES
                         for output in GRAPH_OUTPUT_NAMES]
    println(io, "coefficient_order = [", join(repr.(coefficient_order), ", "), "]")
    println(io, "\n[hashes]")
    for key in sort!(collect(keys(report.hashes)))
        println(io, repr(key), " = ", repr(report.hashes[key]))
    end
    println(io, "\n[models]")
    for model in sort!(collect(keys(report.models)))
        result = report.models[model]
        value = join((result.status, String(result.reason), result.result_sha256), '|')
        println(io, repr(model), " = ", repr(value))
    end
    println(io, "\n[row_counts]")
    for name in sort!(collect(keys(row_counts)))
        println(io, repr(name), " = ", row_counts[name])
    end
    println(io, "\n[graph_handoff]")
    println(io, "fit_hashes = [", join(repr.(fit_hashes), ", "), "]")
    println(io, "full_refit_hashes = [",
            join(repr.(full_fit_hashes), ", "), "]")
    println(io, "model_row_count = ", row_counts["fitted_edge_models.tsv"])
    println(io, "unary_node_row_count = ", row_counts["fitted_unary_nodes.tsv"])
    println(io, "edge_transform_row_count = ", row_counts["fitted_edge_transforms.tsv"])
    println(io, "\n[graph_handoff.schemas]")
    for name in GRAPH_HANDOFF_FILES
        println(io, repr(name), " = ", repr(GRAPH_HANDOFF_SCHEMA))
    end
    println(io, "\n[artifacts]")
    for name in sort!(collect(keys(artifacts)))
        println(io, repr(name), " = ", repr(_sha256(artifacts[name])))
    end
    return take!(io)
end

function _validate_receipt_bytes(
    report::AdmissionReport,
    artifacts::Dict{String,Vector{UInt8}},
    receipt_bytes::Vector{UInt8},
)
    receipt_bytes == _receipt_bytes(report, artifacts) ||
        _fail(:graph_handoff_schema, "receipt bytes are not canonical")
    receipt = try
        TOML.parse(String(copy(receipt_bytes)))
    catch error
        _fail(:graph_handoff_schema, "receipt TOML is invalid: $(sprint(showerror, error))")
    end
    _exact_keys(receipt, Set([
        "schema", "schema_version", "status", "reason", "result_sha256",
        "model_count", "bootstrap_count", "shuffle_count", "fixed_nu",
        "node_priors", "pair_prior", "graph_handoff_schema",
        "graph_handoff_replay_sha256",
        "graph_handoff_files", "category_order", "state_order",
        "coefficient_order", "hashes", "models", "row_counts",
        "graph_handoff", "artifacts",
    ]), "admission receipt")
    receipt["schema"] == "stmfit-structured-edge-admission-receipt-v2" &&
    receipt["schema_version"] == 2 ||
        _fail(:graph_handoff_schema, "admission receipt schema differs")
    receipt["graph_handoff_schema"] == GRAPH_HANDOFF_SCHEMA ||
        _fail(:graph_handoff_schema, "graph-handoff schema differs")
    receipt["graph_handoff_replay_sha256"] ==
        graph_handoff_replay_hash(artifacts) ||
        _fail(:graph_handoff_schema, "graph-handoff replay digest differs")
    receipt["graph_handoff_files"] == GRAPH_HANDOFF_FILES ||
        _fail(:graph_handoff_schema, "graph-handoff file declaration differs")
    receipt["category_order"] == collect(GRAPH_CATEGORY_ORDER) ||
        _fail(:graph_handoff_schema, "category order differs")
    receipt["state_order"] == collect(GRAPH_STATE_ORDER) ||
        _fail(:graph_handoff_schema, "state order differs")
    expected_coefficient_order = ["$(name)_$(output)"
                                  for name in GRAPH_COEFFICIENT_NAMES
                                  for output in GRAPH_OUTPUT_NAMES]
    receipt["coefficient_order"] == expected_coefficient_order ||
        _fail(:graph_handoff_schema, "coefficient order differs")
    receipt["hashes"] == report.hashes ||
        _fail(:graph_handoff_schema, "receipt dependency hashes differ")
    expected_models = Dict(
        model => join((report.models[model].status,
                       String(report.models[model].reason),
                       report.models[model].result_sha256), '|')
        for model in keys(report.models)
    )
    receipt["models"] == expected_models ||
        _fail(:graph_handoff_schema, "receipt model bindings differ")
    expected_artifacts = Dict(name => _sha256(bytes) for (name, bytes) in artifacts)
    receipt["artifacts"] == expected_artifacts ||
        _fail(:graph_handoff_schema, "receipt artifact hashes differ")
    expected_row_counts = Dict(name => _table_data_row_count(bytes, name)
                               for (name, bytes) in artifacts)
    receipt["row_counts"] == expected_row_counts ||
        _fail(:graph_handoff_schema, "receipt row counts differ")
    fit_hashes, full_fit_hashes = _graph_fit_hashes(report)
    graph = receipt["graph_handoff"]
    _exact_keys(graph, Set([
        "fit_hashes", "full_refit_hashes", "model_row_count",
        "unary_node_row_count", "edge_transform_row_count", "schemas",
    ]), "receipt graph handoff")
    graph["fit_hashes"] == fit_hashes &&
    graph["full_refit_hashes"] == full_fit_hashes ||
        _fail(:graph_handoff_schema, "receipt fit bindings differ")
    graph["model_row_count"] == expected_row_counts["fitted_edge_models.tsv"] &&
    graph["unary_node_row_count"] == expected_row_counts["fitted_unary_nodes.tsv"] &&
    graph["edge_transform_row_count"] == expected_row_counts["fitted_edge_transforms.tsv"] ||
        _fail(:graph_handoff_schema, "receipt graph row counts differ")
    graph["schemas"] == Dict(name => GRAPH_HANDOFF_SCHEMA for name in GRAPH_HANDOFF_FILES) ||
        _fail(:graph_handoff_schema, "receipt graph schemas differ")
    return nothing
end

function report_files(report::AdmissionReport)::Dict{String,Vector{UInt8}}
    model_rows = Vector{Vector{String}}()
    partition_rows = Vector{Vector{String}}()
    ess_rows = Vector{Vector{String}}()
    start_rows = Vector{Vector{String}}()
    score_rows = Vector{Vector{String}}()
    shuffle_rows = Vector{Vector{String}}()
    bootstrap_rows = Vector{Vector{String}}()
    for model in sort!(collect(keys(report.models)))
        result = report.models[model]
        push!(model_rows, [
            "1", model, result.status, String(result.reason),
            string(length(result.outer_folds)), result.final_gate.gate_sha256,
            result.full_fit.fit_sha256, result.result_sha256,
        ])
        logical = _logical_evaluations(result)
        for item in logical
            evaluation = item.evaluation
            raw_hash = _raw_heldout_hash(evaluation)
            heldout_scans = sort!(collect(keys(evaluation.scan_scores)))
            push!(partition_rows, [
                "1", model, item.scope, item.outer_date,
                evaluation.target_date, evaluation.status,
                String(evaluation.reason), join(evaluation.training_dates, ','),
                string(length(evaluation.fit.unary.training_scans)),
                string(length(heldout_scans)), evaluation.fit.training_sha256,
                evaluation.fit.fit_sha256, evaluation.score_sha256, raw_hash,
                raw_hash, _out_bool(evaluation.reversal_pass),
                _out_bool(evaluation.shuffle_pass),
            ])
            _append_fit_evidence!(ess_rows, start_rows, model, item.scope,
                                  item.outer_date, evaluation.target_date,
                                  evaluation.fit)
            sample_by_id = Dict(_edge_identity(sample.edge) => sample
                                for sample in evaluation.fit.samples
                                if sample.edge.date == evaluation.target_date)
            for identity in sort!(collect(keys(evaluation.gains)))
                sample = sample_by_id[identity]
                push!(score_rows, [
                    "1", model, item.scope, item.outer_date,
                    evaluation.target_date, "edge", sample.edge.file,
                    string(sample.edge.left_lobe), string(sample.edge.right_lobe),
                    sample.edge.segment_id,
                    _out_float(evaluation.gains[identity]), "1",
                    evaluation.score_sha256,
                ])
            end
            for scan in sort!(collect(keys(evaluation.scan_scores)))
                edge_count = count(sample -> sample.edge.file == scan,
                                   values(sample_by_id))
                push!(score_rows, [
                    "1", model, item.scope, item.outer_date,
                    evaluation.target_date, "scan", scan, "NA", "NA", "NA",
                    _out_float(evaluation.scan_scores[scan]),
                    string(edge_count), evaluation.score_sha256,
                ])
            end
            push!(score_rows, [
                "1", model, item.scope, item.outer_date,
                evaluation.target_date, "date", "NA", "NA", "NA", "NA",
                _out_float(evaluation.date_mean),
                string(length(evaluation.gains)), evaluation.score_sha256,
            ])
            for seed in SHUFFLE_SEEDS
                shuffle_value = length(evaluation.shuffle_values) == 500 ?
                                evaluation.shuffle_values[seed + 1] : NaN
                shuffle_fit_hash =
                    length(evaluation.shuffle_fit_sha256) == 500 ?
                    evaluation.shuffle_fit_sha256[seed + 1] : "not_run"
                push!(shuffle_rows, [
                    "1", model, item.scope, item.outer_date,
                    evaluation.target_date, string(seed),
                    _out_float(shuffle_value), shuffle_fit_hash,
                    _out_float(evaluation.control_date_mean),
                    _out_float(evaluation.shuffle_quantile),
                    _out_bool(evaluation.shuffle_pass), evaluation.score_sha256,
                ])
            end
        end
        _append_fit_evidence!(ess_rows, start_rows, model, "full_refit", "NA",
                              "NA", result.full_fit)
        gates = [("outer_inner", fold.outer_date, fold.inner_gate)
                 for fold in result.outer_folds]
        push!(gates, ("final_lodo", "NA", result.final_gate))
        for (scope, outer_date, gate) in gates
            for seed in BOOTSTRAP_SEEDS
                value = length(gate.bootstrap_values) == 500 ?
                        gate.bootstrap_values[seed + 1] : NaN
                push!(bootstrap_rows, [
                    "1", model, scope, outer_date, string(seed),
                    _out_float(value), _out_float(gate.bootstrap_lower),
                    _out_bool(gate.every_date_positive),
                    _out_bool(gate.every_shuffle_pass),
                    _out_bool(gate.reversal_pass), gate.status,
                    String(gate.reason), gate.gate_sha256,
                ])
            end
        end
    end
    artifacts = Dict(
        "models.tsv" => _table_bytes(MODEL_RESULT_HEADER, model_rows),
        "partitions.tsv" => _table_bytes(PARTITION_HEADER, partition_rows),
        "ess.tsv" => _table_bytes(ESS_HEADER, ess_rows),
        "starts.tsv" => _table_bytes(START_HEADER, start_rows),
        "scores.tsv" => _table_bytes(SCORE_HEADER, score_rows),
        "shuffle.tsv" => _table_bytes(SHUFFLE_HEADER, shuffle_rows),
        "bootstrap.tsv" => _table_bytes(BOOTSTRAP_HEADER, bootstrap_rows),
    )
    graph_artifacts = _graph_handoff_artifacts(report)
    merge!(artifacts, graph_artifacts)
    _validate_graph_handoff_artifacts(report, artifacts)
    files = copy(artifacts)
    files["receipt.toml"] = _receipt_bytes(report, artifacts)
    _validate_receipt_bytes(report, artifacts, files["receipt.toml"])
    return files
end

function _compare_files(destination::String, files::Dict{String,Vector{UInt8}})
    isdir(destination) && !islink(destination) ||
        _fail(:publication_collision, "output destination is not a regular directory")
    sort(readdir(destination)) == sort!(collect(keys(files))) ||
        _fail(:publication_collision, "output file set differs")
    for name in keys(files)
        path = joinpath(destination, name)
        isfile(path) && !islink(path) && read(path) == files[name] ||
            _fail(:publication_collision, "output bytes differ")
    end
    return Dict(name => _sha256(bytes) for (name, bytes) in files)
end

function _publish_files(
    root::String,
    output_directory::AbstractString,
    files::Dict{String,Vector{UInt8}},
)
    destination = StructuredEdgeFeatures._resolve_output(
        root,
        output_directory,
        "edge admission output",
    )
    ispath(destination) && return _compare_files(destination, files)
    parent = dirname(destination)
    stage = mktempdir(parent; prefix=".structured-admission-stage-", cleanup=false)
    installed = false
    try
        for name in sort!(collect(keys(files)))
            StructuredEdgeFeatures._write_synced(joinpath(stage, name), files[name])
        end
        StructuredEdgeFeatures._fsync_directory(stage)
        StructuredEdgeFeatures._run_publication_install_hook(stage, destination)
        StructuredEdgeFeatures._rename_noreplace(stage, destination)
        installed = true
        StructuredEdgeFeatures._fsync_directory(parent)
        return _compare_files(destination, files)
    finally
        !installed && ispath(stage) && rm(stage; recursive=true, force=true)
    end
end

function publish_report(
    root::AbstractString,
    output_directory::AbstractString,
    report::AdmissionReport,
)
    canonical_root = StructuredEdgeFeatures._canonical_root(root)
    return _publish_files(canonical_root, output_directory, report_files(report))
end

function publish_blocker_receipt(
    root::AbstractString,
    output_directory::AbstractString,
    error::EdgeAdmissionError;
    hashes::AbstractDict=Dict{String,String}(),
)
    canonical_root = StructuredEdgeFeatures._canonical_root(root)
    io = IOBuffer()
    println(io, "schema = ", repr("stmfit-structured-edge-admission-receipt-v2"))
    println(io, "schema_version = 2")
    println(io, "status = ", repr("BLOCKED"))
    println(io, "reason = ", repr(String(error.code)))
    println(io, "message_sha256 = ", repr(_sha256(codeunits(error.message))))
    println(io, "plan_sha256 = ", repr(PLAN_SHA256))
    println(io, "t7_review_sha256 = ", repr(T7_REVIEW_SHA256))
    println(io, "t10_review_sha256 = ", repr(T10_REVIEW_SHA256))
    println(io, "\n[hashes]")
    for key in sort!(String.(collect(keys(hashes))))
        value = String(hashes[key])
        occursin(r"^[0-9a-f]{64}$", value) || continue
        println(io, repr(key), " = ", repr(value))
    end
    files = Dict("receipt.toml" => take!(io))
    return _publish_files(canonical_root, output_directory, files)
end

end
