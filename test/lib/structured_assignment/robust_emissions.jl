module StructuredRobustEmissions

using Statistics

export FIXED_NU, EQUAL_MIXTURE_WEIGHTS, SCALE_FLOOR, MAX_ITER,
       CONVERGENCE_TOL, DECREASE_TOL, START_QUANTILE_DELTAS,
       RobustStatus, ROBUST_VALID, ROBUST_ABSTAINED, ROBUST_FAILED,
       RobustInvalidReason, ROBUST_OK, EMPTY_INPUT, TINY_INPUT,
       INVALID_DIMENSION, NONFINITE_INPUT, VARIANCE_COLLAPSE,
       COLLAPSED_COMPONENT, NONFINITE_LIKELIHOOD, NONMONOTONE_PROPOSAL,
       NONCONVERGED, NO_VALID_START, ONE_COMPONENT_EVIDENCE,
       AMPLITUDE_ORIENTATION_UNDEFINED,
       StudentTStartTrace, StudentTOneFit, StudentTTwoFit,
       log_student_t_diag, student_t_responsibilities,
       fit_student_t_one_component, fit_student_t_two_component

const FIXED_NU::Int = 8
const EQUAL_MIXTURE_WEIGHTS = (0.5, 0.5)
const SCALE_FLOOR::Float64 = 1.0e-4
const MAX_ITER::Int = 200
const CONVERGENCE_TOL::Float64 = 1.0e-6
const DECREASE_TOL::Float64 = 1.0e-10
const START_QUANTILE_DELTAS = (0.25, 0.40, 0.45, 0.475, 0.49)
const CLASS_MASS_FLOOR::Float64 = 1.0e-12

@enum RobustStatus::UInt8 begin
    ROBUST_VALID = 0
    ROBUST_ABSTAINED = 1
    ROBUST_FAILED = 2
end

@enum RobustInvalidReason::UInt8 begin
    ROBUST_OK = 0
    EMPTY_INPUT = 1
    TINY_INPUT = 2
    INVALID_DIMENSION = 3
    NONFINITE_INPUT = 4
    VARIANCE_COLLAPSE = 5
    COLLAPSED_COMPONENT = 6
    NONFINITE_LIKELIHOOD = 7
    NONMONOTONE_PROPOSAL = 8
    NONCONVERGED = 9
    NO_VALID_START = 10
    ONE_COMPONENT_EVIDENCE = 11
    AMPLITUDE_ORIENTATION_UNDEFINED = 12
end

struct StudentTStartTrace
    start_index::Int
    quantile_delta::Float64
    accepted_loglik::Vector{Float64}
    converged::Bool
    monotone::Bool
    valid::Bool
    invalid_reason::RobustInvalidReason
    rejected_iteration::Int
    rejected_loglik::Float64
    final_means::Matrix{Float64}
    final_scales::Matrix{Float64}
end

struct StudentTOneFit
    status::RobustStatus
    invalid_reason::RobustInvalidReason
    mean::Vector{Float64}
    scale::Vector{Float64}
    loglik::Float64
    loglik_trace::Vector{Float64}
    converged::Bool
    monotone::Bool
end

struct StudentTTwoFit
    status::RobustStatus
    invalid_reason::RobustInvalidReason
    means::Matrix{Float64}
    scales::Matrix{Float64}
    weights::NTuple{2,Float64}
    loglik::Float64
    responsibilities::Matrix{Float64}
    traces::Vector{StudentTStartTrace}
    best_start::Int
    one_component::StudentTOneFit
    high_amplitude_component::Int
end

function _as_matrix(X)
    Z = if X isa AbstractVector
        reshape(Float64.(X), :, 1)
    elseif X isa AbstractMatrix
        Matrix{Float64}(X)
    else
        throw(ArgumentError("input must be a real vector or matrix"))
    end
    size(Z, 2) > 0 || throw(ArgumentError("input must have positive dimension"))
    return Z
end

function _loggamma_integer_or_half_twice(twice_argument::Int)
    twice_argument > 0 || throw(ArgumentError("gamma argument must be positive"))
    if iseven(twice_argument)
        integer_argument = twice_argument ÷ 2
        result = 0.0
        for value in 1:(integer_argument - 1)
            result += log(Float64(value))
        end
        return result
    end
    half_steps = (twice_argument - 1) ÷ 2
    result = 0.5 * log(π)
    for value in 0:(half_steps - 1)
        result += log(Float64(value) + 0.5)
    end
    return result
end

function log_student_t_diag(x::AbstractVector, mean::AbstractVector,
                            scale::AbstractVector)
    p = length(x)
    p > 0 || throw(ArgumentError("density dimension must be positive"))
    length(mean) == p && length(scale) == p ||
        throw(DimensionMismatch("x, mean, and scale dimensions must match"))
    delta = 0.0
    log_scale_term = 0.0
    @inbounds for j in 1:p
        xj = Float64(x[j])
        meanj = Float64(mean[j])
        scalej = Float64(scale[j])
        isfinite(xj) && isfinite(meanj) ||
            throw(ArgumentError("x and mean must be finite"))
        isfinite(scalej) && scalej > 0.0 ||
            throw(ArgumentError("scale must be positive and finite"))
        difference = xj - meanj
        delta += difference * difference / scalej
        log_scale_term += log(FIXED_NU * π * scalej)
    end
    isfinite(delta) || return -Inf
    return _loggamma_integer_or_half_twice(FIXED_NU + p) -
           _loggamma_integer_or_half_twice(FIXED_NU) -
           0.5 * log_scale_term -
           0.5 * (FIXED_NU + p) * log1p(delta / FIXED_NU)
end

function _logsumexp_pair(first::Float64, second::Float64)
    maximum_value = max(first, second)
    isfinite(maximum_value) || return -Inf
    return maximum_value + log(exp(first - maximum_value) +
                               exp(second - maximum_value))
end

function student_t_responsibilities(X, means::AbstractMatrix,
                                    scales::AbstractMatrix)
    Z = _as_matrix(X)
    n, p = size(Z)
    size(means) == (2, p) || throw(DimensionMismatch("means must be 2 by p"))
    size(scales) == (2, p) || throw(DimensionMismatch("scales must be 2 by p"))
    all(isfinite, Z) || throw(ArgumentError("input must be finite"))
    responsibilities = Matrix{Float64}(undef, n, 2)
    fixed_log_weight = log(0.5)
    for i in 1:n
        first = fixed_log_weight +
                log_student_t_diag(view(Z, i, :), view(means, 1, :),
                                   view(scales, 1, :))
        second = fixed_log_weight +
                 log_student_t_diag(view(Z, i, :), view(means, 2, :),
                                    view(scales, 2, :))
        normalizer = _logsumexp_pair(first, second)
        isfinite(normalizer) || throw(ArgumentError("non-finite responsibility"))
        responsibilities[i, 1] = exp(first - normalizer)
        responsibilities[i, 2] = exp(second - normalizer)
    end
    return responsibilities
end

function _observed_loglik(Z::Matrix{Float64}, means::Matrix{Float64},
                          scales::Matrix{Float64})
    total = 0.0
    fixed_log_weight = log(0.5)
    for i in axes(Z, 1)
        first = fixed_log_weight +
                log_student_t_diag(view(Z, i, :), view(means, 1, :),
                                   view(scales, 1, :))
        second = fixed_log_weight +
                 log_student_t_diag(view(Z, i, :), view(means, 2, :),
                                    view(scales, 2, :))
        total += _logsumexp_pair(first, second)
    end
    return total
end

function _one_observed_loglik(Z::Matrix{Float64}, mean::Vector{Float64},
                              scale::Vector{Float64})
    total = 0.0
    for i in axes(Z, 1)
        total += log_student_t_diag(view(Z, i, :), mean, scale)
    end
    return total
end

function _ecm_proposal(Z::AbstractMatrix, means::AbstractMatrix,
                       scales::AbstractMatrix, responsibilities::AbstractMatrix)
    n, p = size(Z)
    components = size(means, 1)
    size(means, 2) == p && size(scales) == size(means) ||
        throw(DimensionMismatch("parameter dimensions do not match input"))
    size(responsibilities) == (n, components) ||
        throw(DimensionMismatch("responsibility dimensions do not match"))
    new_means = zeros(Float64, components, p)
    new_scales = zeros(Float64, components, p)
    latent_weights = zeros(Float64, n)
    for component in 1:components
        mass = sum(Float64(responsibilities[i, component]) for i in 1:n)
        isfinite(mass) && mass > CLASS_MASS_FLOOR ||
            return (means=Matrix{Float64}(means), scales=Matrix{Float64}(scales),
                    reason=COLLAPSED_COMPONENT)
        weighted_mass = 0.0
        for i in 1:n
            delta = 0.0
            for j in 1:p
                difference = Float64(Z[i, j]) - Float64(means[component, j])
                delta += difference * difference / Float64(scales[component, j])
            end
            latent_weights[i] = (FIXED_NU + p) / (FIXED_NU + delta)
            weighted_mass += Float64(responsibilities[i, component]) *
                             latent_weights[i]
        end
        isfinite(weighted_mass) && weighted_mass > CLASS_MASS_FLOOR ||
            return (means=Matrix{Float64}(means), scales=Matrix{Float64}(scales),
                    reason=COLLAPSED_COMPONENT)
        for j in 1:p
            numerator = 0.0
            for i in 1:n
                numerator += Float64(responsibilities[i, component]) *
                             latent_weights[i] * Float64(Z[i, j])
            end
            new_means[component, j] = numerator / weighted_mass
        end
        for j in 1:p
            numerator = 0.0
            for i in 1:n
                difference = Float64(Z[i, j]) - new_means[component, j]
                numerator += Float64(responsibilities[i, component]) *
                             latent_weights[i] * difference * difference
            end
            new_scales[component, j] = max(numerator / mass, SCALE_FLOOR)
        end
    end
    all(isfinite, new_means) && all(isfinite, new_scales) ||
        return (means=Matrix{Float64}(means), scales=Matrix{Float64}(scales),
                reason=NONFINITE_LIKELIHOOD)
    return (means=new_means, scales=new_scales, reason=ROBUST_OK)
end

function _accept_projected_state(old_means::AbstractMatrix,
                                 old_scales::AbstractMatrix, old_loglik::Real,
                                 new_means::AbstractMatrix,
                                 new_scales::AbstractMatrix, new_loglik::Real)
    accepted = isfinite(new_loglik) &&
               Float64(new_loglik) >= Float64(old_loglik) - DECREASE_TOL
    if accepted
        return (accepted=true, means=Matrix{Float64}(new_means),
                scales=Matrix{Float64}(new_scales), loglik=Float64(new_loglik))
    end
    return (accepted=false, means=Matrix{Float64}(old_means),
            scales=Matrix{Float64}(old_scales), loglik=Float64(old_loglik))
end

function _best_valid_start_index(logliks::AbstractVector,
                                 valid::AbstractVector{Bool})
    length(logliks) == length(valid) || throw(DimensionMismatch("start vectors differ"))
    best = 0
    for index in eachindex(logliks)
        valid[index] && isfinite(logliks[index]) || continue
        if best == 0 || Float64(logliks[index]) > Float64(logliks[best])
            best = index
        end
    end
    return best
end

function _start_parameters(Z::Matrix{Float64}, delta::Float64)
    _, p = size(Z)
    means = zeros(Float64, 2, p)
    for j in 1:p
        values = collect(view(Z, :, j))
        means[1, j] = quantile(values, 0.5 - delta)
        means[2, j] = quantile(values, 0.5 + delta)
    end
    return means, ones(Float64, 2, p)
end

function _run_two_start(Z::Matrix{Float64}, start_index::Int,
                        delta::Float64)
    means, scales = _start_parameters(Z, delta)
    loglik = _observed_loglik(Z, means, scales)
    trace = Float64[loglik]
    if !isfinite(loglik)
        return StudentTStartTrace(start_index, delta, trace, false, true, false,
                                  NONFINITE_LIKELIHOOD, 0, NaN, means, scales)
    end
    for iteration in 1:MAX_ITER
        responsibilities = student_t_responsibilities(Z, means, scales)
        proposal = _ecm_proposal(Z, means, scales, responsibilities)
        if proposal.reason != ROBUST_OK
            return StudentTStartTrace(start_index, delta, trace, false, true, false,
                                      proposal.reason, 0, NaN, means, scales)
        end
        proposed_loglik = _observed_loglik(Z, proposal.means, proposal.scales)
        decision = _accept_projected_state(means, scales, loglik,
                                           proposal.means, proposal.scales,
                                           proposed_loglik)
        if !decision.accepted
            return StudentTStartTrace(start_index, delta, trace, false, false, false,
                                      NONMONOTONE_PROPOSAL, iteration,
                                      proposed_loglik, decision.means,
                                      decision.scales)
        end
        previous_loglik = loglik
        means = decision.means
        scales = decision.scales
        loglik = decision.loglik
        push!(trace, loglik)
        if abs(loglik - previous_loglik) <=
           CONVERGENCE_TOL * (1.0 + abs(previous_loglik))
            return StudentTStartTrace(start_index, delta, trace, true, true, true,
                                      ROBUST_OK, 0, NaN, means, scales)
        end
    end
    return StudentTStartTrace(start_index, delta, trace, false, true, false,
                              NONCONVERGED, 0, NaN, means, scales)
end

function _invalid_one(p::Int, status::RobustStatus,
                      reason::RobustInvalidReason)
    return StudentTOneFit(status, reason, zeros(Float64, p),
                          fill(SCALE_FLOOR, p), -Inf, Float64[], false, true)
end

function _invalid_two(n::Int, p::Int, status::RobustStatus,
                      reason::RobustInvalidReason;
                      traces=StudentTStartTrace[], one_component=_invalid_one(p, status, reason),
                      means=zeros(Float64, 2, p), scales=fill(SCALE_FLOOR, 2, p),
                      loglik=-Inf, responsibilities=fill(NaN, n, 2), best_start=0)
    return StudentTTwoFit(status, reason, means, scales, EQUAL_MIXTURE_WEIGHTS,
                          loglik, responsibilities, traces, best_start,
                          one_component, 0)
end

function _has_variance_collapse(Z::Matrix{Float64})
    n = size(Z, 1)
    for j in axes(Z, 2)
        center = sum(view(Z, :, j)) / n
        variance = sum((Z[i, j] - center)^2 for i in 1:n) / n
        isfinite(variance) && variance > SCALE_FLOOR || return true
    end
    return false
end

function fit_student_t_one_component(X)
    Z = _as_matrix(X)
    n, p = size(Z)
    n == 0 && return _invalid_one(p, ROBUST_ABSTAINED, EMPTY_INPUT)
    all(isfinite, Z) || return _invalid_one(p, ROBUST_FAILED, NONFINITE_INPUT)
    _has_variance_collapse(Z) &&
        return _invalid_one(p, ROBUST_ABSTAINED, VARIANCE_COLLAPSE)
    mean = [median(view(Z, :, j)) for j in 1:p]
    scale = [max(sum((Z[i, j] - mean[j])^2 for i in 1:n) / n,
                 SCALE_FLOOR) for j in 1:p]
    loglik = _one_observed_loglik(Z, mean, scale)
    trace = Float64[loglik]
    isfinite(loglik) ||
        return StudentTOneFit(ROBUST_FAILED, NONFINITE_LIKELIHOOD, mean, scale,
                              loglik, trace, false, true)
    for _ in 1:MAX_ITER
        proposal = _ecm_proposal(Z, reshape(mean, 1, p), reshape(scale, 1, p),
                                 ones(Float64, n, 1))
        proposal.reason == ROBUST_OK ||
            return StudentTOneFit(ROBUST_FAILED, proposal.reason, mean, scale,
                                  loglik, trace, false, true)
        proposed_mean = vec(proposal.means)
        proposed_scale = vec(proposal.scales)
        proposed_loglik = _one_observed_loglik(Z, proposed_mean, proposed_scale)
        decision = _accept_projected_state(reshape(mean, 1, p),
                                           reshape(scale, 1, p), loglik,
                                           proposal.means, proposal.scales,
                                           proposed_loglik)
        if !decision.accepted
            return StudentTOneFit(ROBUST_FAILED, NONMONOTONE_PROPOSAL,
                                  mean, scale, loglik, trace, false, false)
        end
        previous_loglik = loglik
        mean = vec(decision.means)
        scale = vec(decision.scales)
        loglik = decision.loglik
        push!(trace, loglik)
        if abs(loglik - previous_loglik) <=
           CONVERGENCE_TOL * (1.0 + abs(previous_loglik))
            return StudentTOneFit(ROBUST_VALID, ROBUST_OK, mean, scale,
                                  loglik, trace, true, true)
        end
    end
    return StudentTOneFit(ROBUST_FAILED, NONCONVERGED, mean, scale,
                          loglik, trace, false, true)
end

function _one_component_evidence(n::Int, p::Int, two_loglik::Float64,
                                 one_fit::StudentTOneFit)
    one_fit.status == ROBUST_VALID || return false
    improvement = 2.0 * (two_loglik - one_fit.loglik)
    threshold = 2.0 * p * log(max(n, 2))
    return improvement < threshold
end

function _amplitude_orientation(responsibilities::Matrix{Float64},
                                amplitudes::Vector{Float64})
    sums = zeros(Float64, 2)
    counts = zeros(Int, 2)
    for i in axes(responsibilities, 1)
        component = responsibilities[i, 1] >= responsibilities[i, 2] ? 1 : 2
        sums[component] += amplitudes[i]
        counts[component] += 1
    end
    all(>(0), counts) || return 0, COLLAPSED_COMPONENT
    first_mean = sums[1] / counts[1]
    second_mean = sums[2] / counts[2]
    isfinite(first_mean) && isfinite(second_mean) ||
        return 0, AMPLITUDE_ORIENTATION_UNDEFINED
    first_mean == second_mean && return 0, AMPLITUDE_ORIENTATION_UNDEFINED
    return first_mean > second_mean ? 1 : 2, ROBUST_OK
end

function fit_student_t_two_component(X, amplitudes)
    Z = _as_matrix(X)
    n, p = size(Z)
    amplitude_values = Float64.(collect(amplitudes))
    length(amplitude_values) == n ||
        throw(DimensionMismatch("amplitudes must have one value per row"))
    n == 0 && return _invalid_two(n, p, ROBUST_ABSTAINED, EMPTY_INPUT)
    n >= 2 || return _invalid_two(n, p, ROBUST_ABSTAINED, TINY_INPUT)
    all(isfinite, Z) && all(isfinite, amplitude_values) ||
        return _invalid_two(n, p, ROBUST_FAILED, NONFINITE_INPUT)
    _has_variance_collapse(Z) &&
        return _invalid_two(n, p, ROBUST_ABSTAINED, VARIANCE_COLLAPSE)

    one_fit = fit_student_t_one_component(Z)
    one_fit.status == ROBUST_VALID ||
        return _invalid_two(n, p, ROBUST_FAILED, one_fit.invalid_reason;
                            one_component=one_fit)
    traces = StudentTStartTrace[]
    for (start_index, delta) in enumerate(START_QUANTILE_DELTAS)
        push!(traces, _run_two_start(Z, start_index, delta))
    end
    logliks = [isempty(trace.accepted_loglik) ? -Inf : trace.accepted_loglik[end]
               for trace in traces]
    valid = [trace.valid for trace in traces]
    best_start = _best_valid_start_index(logliks, valid)
    best_start > 0 ||
        return _invalid_two(n, p, ROBUST_FAILED, NO_VALID_START;
                            traces=traces, one_component=one_fit)

    best = traces[best_start]
    means = copy(best.final_means)
    scales = copy(best.final_scales)
    responsibilities = student_t_responsibilities(Z, means, scales)
    occupied = falses(2)
    for i in axes(responsibilities, 1)
        occupied[responsibilities[i, 1] >= responsibilities[i, 2] ? 1 : 2] = true
    end
    if count(occupied) < 2
        return _invalid_two(n, p, ROBUST_ABSTAINED, COLLAPSED_COMPONENT;
                            traces=traces, one_component=one_fit, means=means,
                            scales=scales, loglik=best.accepted_loglik[end],
                            responsibilities=responsibilities, best_start=best_start)
    end
    if _one_component_evidence(n, p, best.accepted_loglik[end], one_fit)
        return _invalid_two(n, p, ROBUST_ABSTAINED, ONE_COMPONENT_EVIDENCE;
                            traces=traces, one_component=one_fit, means=means,
                            scales=scales, loglik=best.accepted_loglik[end],
                            responsibilities=responsibilities, best_start=best_start)
    end
    high_component, orientation_reason =
        _amplitude_orientation(responsibilities, amplitude_values)
    if orientation_reason != ROBUST_OK
        return _invalid_two(n, p, ROBUST_ABSTAINED, orientation_reason;
                            traces=traces, one_component=one_fit, means=means,
                            scales=scales, loglik=best.accepted_loglik[end],
                            responsibilities=responsibilities, best_start=best_start)
    end
    return StudentTTwoFit(ROBUST_VALID, ROBUST_OK, means, scales,
                          EQUAL_MIXTURE_WEIGHTS, best.accepted_loglik[end],
                          responsibilities, traces, best_start, one_fit,
                          high_component)
end

end
