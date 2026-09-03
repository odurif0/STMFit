# Physical component mapping (higher-amplitude cluster -> GlcNAc (1)), the
# per-view fit, the in-sample one-component identifiability check, and the
# equal-weight view log-likelihood combination.

function physical_high_amplitude_cluster(records::Vector{LobeRecord},
                                         assignments::Vector{Int})
    K = maximum(assignments; init=2)
    amps = [r.amplitude for r in records]
    mean_per_cluster = zeros(K)
    count_per_cluster = zeros(Int, K)
    for (i, c) in enumerate(assignments)
        1 <= c <= K || continue
        mean_per_cluster[c] += amps[i]
        count_per_cluster[c] += 1
    end
    best_c = 1
    best_amp = -Inf
    for c in 1:K
        if count_per_cluster[c] > 0
            m = mean_per_cluster[c] / count_per_cluster[c]
            if m > best_amp
                best_amp = m
                best_c = c
            end
        end
    end
    return best_c
end

mutable struct ViewFit
    name::String
    features::Vector{String}
    fit2::TwoComponentFit          # 2-component model
    fit1::OneComponentFit          # 1-component companion
    valid_lobe_idx::Vector{Int}    # indices of records with finite features in this view
    high_cluster::Int              # 1 or 2 — cluster mapped to GlcNAc (1)
    one_component_evidence::Bool   # means collapsed / empty cluster / LL near 1-comp
    degenerate::Bool               # zero-MAD across all scans in every feature
    responsibility_for_label1::Vector{Float64}  # per-record (NaN if invalid)
    log_odds_for_label1::Vector{Float64}
end

# Heuristic identifiability check (in-sample; the held-out gate is T5). The
# 2-component model has one-component evidence when (a) the two cluster means
# are essentially indistinguishable in normalized space, (b) the BIC
# improvement of the 2-component model over the 1-component companion is
# negligible given the extra parameters, or (c) the raw amplitudes do not
# separate the two clusters (the physical higher-amplitude mapping is then
# undefined). Weights are fixed at (0.5, 0.5), so only means/variances count
# as free parameters.
function _one_component_evidence(fit2::TwoComponentFit, fit1::OneComponentFit,
                                 Z::AbstractMatrix,
                                 valid_records::Vector{LobeRecord},
                                 hard_assignments::Vector{Int};
                                 sep_abs::Float64=0.25)
    n, p = size(Z)
    occupied = falses(2)
    for c in hard_assignments
        1 <= c <= 2 && (occupied[c] = true)
    end
    count(occupied) < 2 && return true
    sep = maximum(abs(fit2.means[1, j] - fit2.means[2, j]) for j in 1:p)
    sep < sep_abs && return true
    extra_params = 2p
    ll_improvement = 2 * (fit2.loglik - fit1.loglik)
    bic_threshold = extra_params * log(max(n, 2))
    ll_improvement < bic_threshold && return true
    amps = Float64[r.amplitude for r in valid_records]
    amp_med = _median_finite(amps)
    amp_mad = _mad_finite(amps, amp_med)
    if isfinite(amp_mad) && amp_mad > 0.0
        K = maximum(hard_assignments; init = 2)
        sum_per = zeros(K)
        cnt_per = zeros(Int, K)
        for (a, c) in zip(amps, hard_assignments)
            1 <= c <= K || continue
            sum_per[c] += a
            cnt_per[c] += 1
        end
        means_amp = [cnt_per[c] > 0 ? sum_per[c] / cnt_per[c] : NaN for c in 1:K]
        valid_means = filter(isfinite, means_amp)
        if length(valid_means) >= 2
            amp_gap = maximum(valid_means) - minimum(valid_means)
            if amp_gap < 0.5 * amp_mad
                return true
            end
        end
    end
    return false
end

function fit_view(records::Vector{LobeRecord}, features::Vector{String},
                  name::AbstractString;
                  first_seed::Integer=0, n_starts::Int=default_n_starts,
                  cov_floor::Float64=cov_floor_default,
                  max_iter::Int=200, tol::Float64=1e-6)
    for fn in features
        is_forbidden_feature(fn) &&
            error("View $name uses a forbidden feature: $fn")
    end
    X, row_valid = feature_matrix(records, features)
    Z = normalize_per_scan(records, X, features)
    valid_idx = Int[i for i in 1:length(records) if all(isfinite, view(Z, i, :))]
    n_valid = length(valid_idx)
    degenerate = n_valid == 0
    one_comp = false
    fit2 = nothing
    fit1 = nothing
    high_c = 1
    resp_for1 = fill(NaN, length(records))
    logit_for1 = fill(NaN, length(records))
    if degenerate
        p = length(features)
        fit2 = TwoComponentFit(zeros(2, p), fill(cov_floor, 2, p),
                               copy(EQUAL_PRIOR_WEIGHTS), NaN, Float64[],
                               false, true)
        fit1 = OneComponentFit(zeros(p), fill(cov_floor, p), NaN)
    else
        Zv = Z[valid_idx, :]
        fit2 = fit_em_two_component(Zv, first_seed; n_starts=n_starts,
                                    cov_floor=cov_floor, max_iter=max_iter, tol=tol)
        fit1 = fit_one_component(Zv; cov_floor=cov_floor)
        if fit2.converged && fit2.monotone
            resp = responsibilities(fit2, Zv)
            hard = [argmax(view(resp, i, :)) for i in 1:size(resp, 1)]
            valid_records = records[valid_idx]
            one_comp = _one_component_evidence(fit2, fit1, Zv, valid_records, hard)
            if !one_comp
                high_c = physical_high_amplitude_cluster(valid_records, hard)
                for (k, i) in enumerate(valid_idx)
                    r = resp[k, high_c]
                    if isfinite(r) && 0.0 < r < 1.0
                        resp_for1[i] = r
                        rc = clamp(r, 1e-12, 1.0 - 1e-12)
                        logit_for1[i] = log(rc / (1.0 - rc))
                    elseif isfinite(r)
                        resp_for1[i] = r
                        logit_for1[i] = r >= 0.5 ? 30.0 : -30.0
                    end
                end
            end
        end
    end
    return ViewFit(String(name), String[String(f) for f in features],
                   fit2, fit1, valid_idx, high_c, one_comp, degenerate,
                   resp_for1, logit_for1)
end

struct CombinedPrediction
    probability_1::Vector{Float64}
    views_used::Vector{Int}
    invalid_reason::Vector{String}
end

# Equal-weight view combination in log-odds space: missing views are omitted
# and the count of used views is reported.
function combine_views(fits::Vector{ViewFit}, n_records::Int)
    p1 = fill(NaN, n_records)
    used = zeros(Int, n_records)
    reason = fill("ok", n_records)
    logit_sum = zeros(n_records)
    for vf in fits
        (!vf.fit2.converged || !vf.fit2.monotone) && continue
        vf.one_component_evidence && continue
        vf.degenerate && continue
        for rec_idx in 1:n_records
            lo = vf.log_odds_for_label1[rec_idx]
            if isfinite(lo)
                logit_sum[rec_idx] += lo
                used[rec_idx] += 1
            end
        end
    end
    for i in 1:n_records
        if used[i] == 0
            p1[i] = NaN
            if all(vf.degenerate for vf in fits)
                reason[i] = "degenerate_view"
            elseif any(!vf.fit2.converged for vf in fits)
                reason[i] = "unstable_model"
            elseif any(!vf.fit2.monotone for vf in fits)
                reason[i] = "nonmonotone_model"
            elseif all(vf.one_component_evidence for vf in fits)
                reason[i] = "one_component_evidence"
            else
                reason[i] = "missing_view"
            end
        else
            mean_logit = logit_sum[i] / used[i]
            p1[i] = 1.0 / (1.0 + exp(-mean_logit))
            if used[i] < length(fits) && length(fits) > 0
                reason[i] = "ok_partial_views"
            else
                reason[i] = "ok"
            end
        end
    end
    return CombinedPrediction(p1, used, reason)
end
