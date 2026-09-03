module StructuredChainInference

using LinearAlgebra
using Printf
using SHA
using TOML

const _DP_CALL_COUNT = Ref(0)
const _FACTOR_BUILD_CALL_COUNT = Ref(0)
const _SCAN_DECODE_CALL_COUNT = Ref(0)
const _SNAPSHOT_RETURN_HOOK = Ref{Union{Nothing,Function}}(nothing)
const _REVERSAL_CHECK_HOOK = Ref{Union{Nothing,Function}}(nothing)
const _POST_VALIDATION_HOOK = Ref{Union{Nothing,Function}}(nothing)

struct _InferenceFailed <: Exception
    code::Symbol
    message::String
end
_infer_fail(code::Symbol, message) = throw(_InferenceFailed(code, String(message)))

struct ChainFactorResult
    file::String
    date::String
    left_lobe::Int
    right_lobe::Int
    left_t_nm::Float64
    right_t_nm::Float64
    segment_id::String
    factor::NTuple{4,Float64}
    factor_sha256::String
    raw_row_sha256::String
    fit_sha256::String
    model::String
    unary_fit_sha256::String
    residualizer_sha256::String
end

struct ChainNodeResult
    file::String
    date::String
    lobe::Int
    t_nm::Float64
    log_marginals::NTuple{2,Float64}
    marginals::NTuple{2,Float64}
    label::Int
    unary_fit_sha256::String
end

struct ChainBlockResult
    segment_id::String
    nodes::Tuple{Vararg{ChainNodeResult}}
    factors::Tuple{Vararg{ChainFactorResult}}
    log_evidence::Float64
    viterbi_score::Float64
    labels::Tuple{Vararg{Int}}
end

struct ChainReferenceResult
    fit_sha256::String
    model::String
    fit_role::String
    scope::String
    outer_date::String
    target_date::String
    status::Symbol
    reason::Symbol
    blocks::Tuple{Vararg{ChainBlockResult}}
    log_evidence::Float64
    viterbi_score::Float64
end

struct ChainInferenceReport
    status::Symbol
    reason::Symbol
    references::Tuple{Vararg{ChainReferenceResult}}
    blocks::Tuple{Vararg{ChainBlockResult}}
    log_evidence::Float64
    viterbi_score::Float64
    result_sha256::String
    provenance_sha256::String
end

export ChainNodeResult, ChainFactorResult, ChainBlockResult, ChainReferenceResult,
       ChainInferenceReport, infer_structured_chains

include(joinpath(@__DIR__, "edge_model.jl"))

function _empty_report(status::Symbol, reason::Symbol; provenance="")
    result = _hash_lines(["schema=structured-chain-inference-v1", "status=$status", "reason=$reason", "provenance=$provenance"])
    return ChainInferenceReport(status, reason, (), (), 0.0, 0.0, result, provenance)
end

function _reference_result_from_failure(reference::_ReferenceData)
    return ChainReferenceResult(reference.fit_sha256, reference.model, reference.fit_role,
                                reference.scope, reference.outer_date, reference.target_date,
                                :FAIL, :todo11_fail,
                                (), 0.0, 0.0)
end

function _report_hash(status, reason, references, provenance)
    lines = String["schema=structured-chain-inference-v1", "status=$status", "reason=$reason"]
    for reference in references
        push!(lines, join((reference.fit_sha256, reference.model, reference.fit_role,
                           reference.scope, reference.outer_date, reference.target_date,
                           reference.status, reference.reason), '|'))
        for block in reference.blocks
            push!(lines, "block=$(block.segment_id):$(block.log_evidence):$(block.viterbi_score)")
            for factor in block.factors
                push!(lines, join((factor.file, factor.date, factor.left_lobe,
                                   factor.right_lobe, factor.segment_id,
                                   factor.factor_sha256, factor.raw_row_sha256,
                                   factor.fit_sha256, factor.model,
                                   factor.unary_fit_sha256,
                                   factor.residualizer_sha256), '|'))
                push!(lines, join(factor.factor, '|'))
            end
            for node in block.nodes
                push!(lines, join((node.file, node.date, node.lobe, @sprintf("%.17g", node.t_nm),
                                   node.log_marginals[1], node.log_marginals[2],
                                   node.marginals[1], node.marginals[2], node.label), '|'))
            end
        end
    end
    push!(lines, "provenance=$provenance")
    return _hash_lines(lines)
end

function _build_report(references::Tuple, provenance::String)
    statuses = [reference.status for reference in references]
    status = any(==( :FAIL), statuses) ? :FAIL : any(==( :PASS), statuses) ? :PASS : :SKIPPED
    reason = status == :FAIL ? :todo11_fail : status == :PASS ? :graph_inference : :unary_only
    blocks = Tuple(vcat((collect(reference.blocks) for reference in references)...))
    log_evidence = sum(reference.log_evidence for reference in references)
    viterbi_score = sum(reference.viterbi_score for reference in references)
    result = _report_hash(status, reason, references, provenance)
    return ChainInferenceReport(status, reason, references, blocks, log_evidence,
                                viterbi_score, result, provenance)
end

function _transpose_factor(factor::NTuple{4,Float64})
    return (factor[1], factor[3], factor[2], factor[4])
end

function _factor_sha(edge::_Todo10Edge, model::_ModelData, factor::NTuple{4,Float64})
    return _hash_lines([
        "schema=structured-chain-factor-v1",
        "fit_sha256=$(model.fit_sha256)",
        "model=$(model.model)",
        "file=$(edge.file)",
        "date=$(edge.date)",
        "left_lobe=$(edge.left_lobe)",
        "right_lobe=$(edge.right_lobe)",
        "left_t_nm=$(_fmt(edge.left_t_nm))",
        "right_t_nm=$(_fmt(edge.right_t_nm))",
        "segment_id=$(edge.segment_id)",
        "raw_row_sha256=$(edge.raw_sha256)",
        "unary_fit_sha256=$(model.unary_fit_sha256)",
        "residualizer_sha256=$(model.residualizer_sha256)",
        "factor=$(join(_fmt.(collect(factor)), '|'))",
    ])
end

function _assert_reversal_equivalence(natural, reversed)
    hook = _REVERSAL_CHECK_HOOK[]
    _REVERSAL_CHECK_HOOK[] = nothing
    hook !== nothing && hook(natural, reversed)
    isapprox(natural[1], reversed[1]; atol=1e-12, rtol=1e-12) || _infer_fail(:reversal_mismatch, "reversal log evidence differs")
    isapprox(natural[5], reversed[5]; atol=1e-12, rtol=1e-12) || _infer_fail(:reversal_mismatch, "reversal Viterbi score differs")
    isapprox(natural[3], reversed[3][end:-1:1, :]; atol=1e-12, rtol=1e-12) || _infer_fail(:reversal_mismatch, "reversal marginals differ")
    isapprox(natural[2], reversed[2][end:-1:1, :]; atol=1e-12, rtol=1e-12) || _infer_fail(:reversal_mismatch, "reversal log marginals differ")
    natural[4] == reversed[4][end:-1:1] || _infer_fail(:reversal_mismatch, "reversal labels differ")
    return nothing
end

function _unary_only_block(segment_id::String, nodes::Vector{_Todo10Node}, unary_map)
    node_results = ChainNodeResult[]
    labels = Int[]
    log_evidence = 0.0
    viterbi = 0.0
    for node in nodes
        unary = unary_map[(node.file, node.lobe)]
        log_unary = (unary.p0 == 0.0 ? -Inf : log(unary.p0), unary.p1 == 0.0 ? -Inf : log(unary.p1))
        isfinite(log_unary[1]) || isfinite(log_unary[2]) || _infer_fail(:unary_failure, "unary has no feasible state")
        normalizer = _logsumexp(collect(log_unary))
        marginals = _binary_softmax(log_unary[1], log_unary[2])
        label = _argmax_lower(collect(log_unary)) - 1
        push!(node_results, ChainNodeResult(node.file, node.date, node.lobe, node.t_nm,
                                            log_unary, marginals, label,
                                            unary.unary_fit_sha256))
        push!(labels, label)
        log_evidence += normalizer
        viterbi += log_unary[label + 1]
    end
    return ChainBlockResult(segment_id, Tuple(node_results), (), log_evidence,
                            viterbi, Tuple(labels))
end

function _assert_distinct_paths(root::String, authority::Vector{_Snapshot}, paths)
    seen = Set(snapshot.path for snapshot in authority)
    supplied = Set{String}()
    for (path, context) in paths
        absolute = _resolve_under_root(root, path, context)
        absolute in supplied && _blocked(:source_collision, "duplicate supplied path: $context")
        absolute in seen && _blocked(:source_collision, "supplied path duplicates fixed authority: $context")
        push!(supplied, absolute)
    end
end

function infer_structured_chains(
    root::AbstractString;
    edge_observations::AbstractString,
    node_segments::AbstractString,
    todo10_receipt::AbstractString,
    fitted_edge_models::AbstractString,
    fitted_unary_nodes::AbstractString,
    fitted_edge_transforms::AbstractString,
    todo11_receipt::AbstractString,
)::ChainInferenceReport
    _DP_CALL_COUNT[] = 0
    _FACTOR_BUILD_CALL_COUNT[] = 0
    _SCAN_DECODE_CALL_COUNT[] = 0
    validation_complete = false
    try
        canonical_root = _canonical_root(root)
        authority = _authority_snapshots(canonical_root)
        registry = Dict(snapshot.relative => snapshot for snapshot in authority)
        _assert_distinct_paths(canonical_root, authority, [
            (edge_observations, "Todo 10 edge observations"),
            (node_segments, "Todo 10 node segments"),
            (todo10_receipt, "Todo 10 receipt"),
            (fitted_edge_models, "Todo 11 fitted edge models"),
            (fitted_unary_nodes, "Todo 11 fitted unary nodes"),
            (fitted_edge_transforms, "Todo 11 fitted edge transforms"),
            (todo11_receipt, "Todo 11 receipt"),
        ])
        for (path, context) in [
            (edge_observations, "Todo 10 edge observations"),
            (node_segments, "Todo 10 node segments"),
            (todo10_receipt, "Todo 10 receipt"),
            (fitted_edge_models, "Todo 11 fitted edge models"),
            (fitted_unary_nodes, "Todo 11 fitted unary nodes"),
            (fitted_edge_transforms, "Todo 11 fitted edge transforms"),
            (todo11_receipt, "Todo 11 receipt"),
        ]
            _register_snapshot!(canonical_root, path, context, registry)
        end
        todo10 = _validate_todo10(canonical_root, edge_observations, node_segments,
                                  todo10_receipt, registry)
        todo11 = _validate_todo11(canonical_root, fitted_edge_models, fitted_unary_nodes,
                                   fitted_edge_transforms, todo11_receipt, registry,
                                   todo10)
        validation_complete = true
        hook = _POST_VALIDATION_HOOK[]
        _POST_VALIDATION_HOOK[] = nothing
        if hook !== nothing
            hook()
            for snapshot in values(registry)
                _resnapshot(snapshot)
            end
        end
        references = ChainReferenceResult[]
        factor_cache = Dict{Any,Any}()
        scan_cache = Dict{Any,Any}()
        ordered_models = sort!(collect(todo11.models); by=model -> model.fit_sha256)
        for model in ordered_models
            for serialized_reference in model.references
                if todo11.status == "FAIL"
                    push!(references, _reference_result_from_failure(serialized_reference))
                else
                    push!(references, _reference_result(serialized_reference, model, todo10, todo11,
                                                        factor_cache, scan_cache;
                                                        unary_only=todo11.status == "SKIPPED"))
                end
            end
        end
        isempty(references) && _blocked(:reference_mismatch, "Todo 11 contains no references")
        hook = _SNAPSHOT_RETURN_HOOK[]
        _SNAPSHOT_RETURN_HOOK[] = nothing
        hook !== nothing && hook()
        for snapshot in values(registry)
            _resnapshot(snapshot)
        end
        authority_lines = sort([snapshot.relative * "=" * snapshot.sha256 for snapshot in authority])
        provenance = _hash_lines(["authority=$(join(authority_lines, '\n'))", "todo10=$(todo10.provenance_sha256)", "todo11=$(todo11.provenance_sha256)"])
        return _build_report(Tuple(references), provenance)
    catch error
        if error isa _ChainBlocked
            return _empty_report(:BLOCKED, error.code; provenance=_hash_lines(["blocked=$(error.code)", error.message]))
        end
        if error isa _InferenceFailed || validation_complete
            code = error isa _InferenceFailed ? error.code : :internal_inference
            return _empty_report(:FAIL, code; provenance=_hash_lines(["failed=$(code)", sprint(showerror, error)]))
        end
        return _empty_report(:BLOCKED, :internal_validation; provenance=_hash_lines([sprint(showerror, error)]))
    end
end

end
