module StructuredUnitAssignment

using Printf
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "hierarchical_unit_assignment.jl"))

module StructuredFirewall
include(joinpath(@__DIR__, "structured_assignment", "firewall.jl"))
end

module FrozenChampionAdapter
include(joinpath(@__DIR__, "structured_assignment", "champion_adapter.jl"))
end

include(joinpath(@__DIR__, "structured_assignment", "robust_emissions.jl"))

const HUA = HierarchicalUnitAssignment
const SAF = StructuredFirewall
const FCA = FrozenChampionAdapter
const SRE = StructuredRobustEmissions
const ROOT = realpath(joinpath(@__DIR__, "..", ".."))
const MODEL_CONFIG_PATH = joinpath(ROOT, "config", "unit_assignment_structured_model.toml")
const MODEL_CONFIG_SHA256 =
    "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const CHAMPION_ADAPTER_SHA256 =
    "470a4c6676eaa3b0a00a51a64f9fce96ac7282fd1adb84c921dc3b967f6aa83c"
const ROBUST_EMISSIONS_SHA256 =
    "0ca4a3863a4b9ed13aee5cd3a0229df4f2764a8524b764eeb0fd27238f14f241"
const STRUCTURED_FIREWALL_SHA256 =
    "a189898d31759352d9dac0de8ff281d6001e1dc2db94df2cb725ff597537a9bb"
const SUPPORTED_MODELS = ["gaussian_one", "gaussian_two", "student_t_two"]
const OUTPUT_COLUMNS = [
    "file", "lobe", "predicted", "probability_1", "confidence", "model",
    "status", "invalid_reason", "views_used", "segment_id", "edge_status",
    "edge_admission_hash", "trigger_hash", "source_hash", "config_hash",
    "input_hash",
]
const BASE_FEATURES = [
    "amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel",
]
const VIEW_SPECS = [
    "base_local" => copy(BASE_FEATURES),
    "base_local+bwd_neg_com_t" => vcat(BASE_FEATURES, ["bwd_neg_com_t"]),
    "base_local+bwd_neg_diag45" => vcat(BASE_FEATURES, ["bwd_neg_diag45"]),
    "base_local+split_log_skew" => vcat(BASE_FEATURES, ["split_log_skew"]),
]

export SUPPORTED_MODELS, OUTPUT_COLUMNS, ModelContract, UnaryPrediction,
       PredictionBundle, load_contract, load_champion_unaries,
       build_predictions, prediction_tsv_bytes, write_prediction_tsv,
       source_hash, trigger_hash, HierarchicalUnitAssignment,
       StructuredRobustEmissions

struct StructuredUnaryError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, error::StructuredUnaryError) =
    print(io, error.code, ": ", error.message)

struct FileSnapshot
    path::String
    bytes::Vector{UInt8}
    sha256::String
end

struct ModelContract
    config_path::String
    config_hash::String
    views::Vector{Pair{String,Vector{String}}}
    n_starts::Int
    covariance_floor::Float64
    max_iter::Int
    tolerance::Float64
    snapshot::FileSnapshot
end

struct UnaryPrediction
    file::String
    lobe::Int
    predicted::String
    probability_1::Union{Nothing,Float64}
    confidence::Float64
    model::String
    status::String
    invalid_reason::String
    views_used::Int
    segment_id::String
    edge_status::String
    edge_admission_hash::String
    trigger_hash::String
    source_hash::String
    config_hash::String
    input_hash::String
end

struct PredictionBundle
    rows::Vector{UnaryPrediction}
    snapshots::Vector{FileSnapshot}
end

struct InputBundle
    records::Vector{HUA.LobeRecord}
    snapshot::FileSnapshot
end

struct RequiredViewFit
    fit::HUA.ViewFit
    reason::String
    failed::Bool
    observed::BitVector
    normalized_valid::BitVector
end

struct RequiredViewBlocker
    reason::String
    failed::Bool
end

const _before_noreplace_hook = Ref{Union{Nothing,Function}}(nothing)

_fail(code::Symbol, message::String) = throw(StructuredUnaryError(code, message))
_sha256(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))

function _has_parent(path::String)
    return any(==(".."), split(replace(path, '\\' => '/'), '/'; keepempty=true))
end

function _absolute_path(path::AbstractString)
    supplied = String(path)
    isempty(supplied) && _fail(:invalid_path, "path is empty")
    occursin('\0', supplied) && _fail(:invalid_path, "path contains NUL")
    _has_parent(supplied) && _fail(:path_escape, "parent traversal is forbidden")
    return normpath(isabspath(supplied) ? supplied : joinpath(ROOT, supplied))
end

function _reject_symlinks(path::String)
    cursor = first(splitpath(abspath(path)))
    islink(cursor) && _fail(:symlink_rejected, "path contains a symbolic link")
    for component in splitpath(abspath(path))[2:end]
        cursor = joinpath(cursor, component)
        islink(cursor) && _fail(:symlink_rejected, "path contains a symbolic link")
    end
    return nothing
end

function _snapshot(path::AbstractString; expected_path::Union{Nothing,String}=nothing,
                   expected_hash::Union{Nothing,String}=nothing)
    absolute = _absolute_path(path)
    _reject_symlinks(absolute)
    isfile(absolute) || _fail(:missing_input, "required input is absent")
    if expected_path !== nothing
        absolute == expected_path || _fail(:path_substitution, "immutable path differs")
    end
    bytes = read(absolute)
    digest = _sha256(bytes)
    if expected_hash !== nothing
        digest == expected_hash || _fail(:hash_mismatch, "immutable bytes differ")
    end
    return FileSnapshot(absolute, bytes, digest)
end

function _verify(snapshot::FileSnapshot)
    _reject_symlinks(snapshot.path)
    isfile(snapshot.path) || _fail(:stale_input, "input disappeared")
    read(snapshot.path) == snapshot.bytes || _fail(:stale_input, "input bytes changed")
    return nothing
end

function _expect(document::AbstractDict, key::String, expected, context::String)
    get(document, key, nothing) == expected ||
        _fail(:config_contract, "$context.$key differs from the frozen contract")
end

function _table(document::AbstractDict, key::String, context::String)
    value = get(document, key, nothing)
    value isa AbstractDict || _fail(:config_contract, "$context.$key is absent")
    return value
end

function load_contract(path::AbstractString)
    snapshot = _snapshot(path; expected_path=MODEL_CONFIG_PATH,
                         expected_hash=MODEL_CONFIG_SHA256)
    document = try
        TOML.parse(String(copy(snapshot.bytes)))
    catch error
        _fail(:config_parse, sprint(showerror, error))
    end
    model = _table(document, "model", "config")
    _expect(model, "schema", "structured_label_free_unit_assignment_model_v2", "model")
    _expect(model, "version", 2, "model")
    _expect(model, "model_ids", ["C1", "C2", "edge_admission", "admitted_graph"],
            "model")
    _expect(model, "reference_family", "two_component_diagonal_gaussian", "model")
    _expect(model, "challenger_family", "two_component_diagonal_student_t", "model")
    _expect(model, "node_priors", [0.5, 0.5], "model")
    _expect(model, "mixture_weights", [0.5, 0.5], "model")
    _expect(model, "mixture_weight_policy", "fixed_not_estimated", "model")
    _expect(model, "covariance", "diagonal", "model")
    _expect(model, "covariance_floor", SRE.SCALE_FLOOR, "model")
    _expect(model, "student_t_nu", SRE.FIXED_NU, "model")
    _expect(model, "n_starts", length(SRE.START_QUANTILE_DELTAS), "model")
    _expect(model, "unary_start_quantile_deltas", collect(SRE.START_QUANTILE_DELTAS),
            "model")
    _expect(model, "max_iter", SRE.MAX_ITER, "model")
    _expect(model, "tol", SRE.CONVERGENCE_TOL, "model")
    _expect(model, "objective_decrease_tolerance", SRE.DECREASE_TOL, "model")
    _expect(model, "physical_orientation_feature", "amplitude", "model")

    views = _table(model, "views", "model")
    _expect(views, "names", first.(VIEW_SPECS), "model.views")
    _expect(views, "base_local", BASE_FEATURES, "model.views")
    _expect(views, "backward_descriptors", ["bwd_neg_com_t", "bwd_neg_diag45"],
            "model.views")
    _expect(views, "split_descriptor", ["split_log_skew"], "model.views")
    _expect(views, "fusion", "equal_logit_mean", "model.views")
    _expect(views, "fixed_equal_weights", [0.25, 0.25, 0.25, 0.25], "model.views")
    _expect(views, "missing_view_policy", "marginalize_and_report_count", "model.views")
    _expect(views, "all_views_missing_posterior", [0.5, 0.5], "model.views")

    preprocessing = _table(document, "preprocessing", "config")
    _expect(preprocessing, "normalization", "per_scan_median_mad", "preprocessing")
    _expect(preprocessing, "normalization_scope", "each_scan_uses_only_its_own_lobes",
            "preprocessing")
    _verify(snapshot)
    return ModelContract(
        snapshot.path,
        snapshot.sha256,
        [name => copy(features) for (name, features) in VIEW_SPECS],
        Int(model["n_starts"]),
        Float64(model["covariance_floor"]),
        Int(model["max_iter"]),
        Float64(model["tol"]),
        snapshot,
    )
end

load_champion_unaries(arguments...) = FCA.load_frozen_champion_unaries(arguments...)

function _source_entries()
    hierarchical = joinpath(@__DIR__, "hierarchical")
    structured = joinpath(@__DIR__, "structured_assignment")
    return [
        "test/build_structured_unit_predictions.jl" =>
            joinpath(ROOT, "test", "build_structured_unit_predictions.jl"),
        "test/lib/hierarchical/emission_math.jl" => joinpath(hierarchical, "emission_math.jl"),
        "test/lib/hierarchical/emissions.jl" => joinpath(hierarchical, "emissions.jl"),
        "test/lib/hierarchical/firewall.jl" => joinpath(hierarchical, "firewall.jl"),
        "test/lib/hierarchical/io_helpers.jl" => joinpath(hierarchical, "io_helpers.jl"),
        "test/lib/hierarchical/loading.jl" => joinpath(hierarchical, "loading.jl"),
        "test/lib/hierarchical/nuisance.jl" => joinpath(hierarchical, "nuisance.jl"),
        "test/lib/hierarchical/pipeline.jl" => joinpath(hierarchical, "pipeline.jl"),
        "test/lib/hierarchical/views.jl" => joinpath(hierarchical, "views.jl"),
        "test/lib/hierarchical_unit_assignment.jl" =>
            joinpath(@__DIR__, "hierarchical_unit_assignment.jl"),
        "test/lib/structured_assignment/champion_adapter.jl" =>
            joinpath(structured, "champion_adapter.jl"),
        "test/lib/structured_assignment/firewall.jl" => joinpath(structured, "firewall.jl"),
        "test/lib/structured_assignment/robust_emissions.jl" =>
            joinpath(structured, "robust_emissions.jl"),
        "test/lib/structured_unit_assignment.jl" => abspath(@__FILE__),
    ]
end

function _source_bundle()
    snapshots = FileSnapshot[]
    material = IOBuffer()
    println(material, "structured-unary-source-v1")
    for (relative, path) in sort(_source_entries(); by=first)
        snapshot = _snapshot(path)
        if endswith(relative, "champion_adapter.jl")
            snapshot.sha256 == CHAMPION_ADAPTER_SHA256 ||
                _fail(:dependency_hash, "Todo 6 adapter bytes differ")
        elseif endswith(relative, "robust_emissions.jl")
            snapshot.sha256 == ROBUST_EMISSIONS_SHA256 ||
                _fail(:dependency_hash, "Todo 7 emission bytes differ")
        elseif relative == "test/lib/structured_assignment/firewall.jl"
            snapshot.sha256 == STRUCTURED_FIREWALL_SHA256 ||
                _fail(:dependency_hash, "structured firewall bytes differ")
        end
        push!(snapshots, snapshot)
        println(material, relative, '=', snapshot.sha256)
    end
    return _sha256(take!(material)), snapshots
end

function source_hash()
    digest, snapshots = _source_bundle()
    digest == LOADED_SOURCE_HASH ||
        _fail(:stale_source, "loaded source differs from current bytes")
    foreach(_verify, LOADED_SOURCE_SNAPSHOTS)
    foreach(_verify, snapshots)
    return LOADED_SOURCE_HASH
end

function trigger_hash(model::AbstractString, config_hash::AbstractString)
    String(model) in SUPPORTED_MODELS || _fail(:unsupported_model, "model is not configured")
    occursin(r"^[0-9a-f]{64}$", config_hash) ||
        _fail(:invalid_hash, "config hash is malformed")
    material = "structured-unary-trigger-v1\nmodel=$(String(model))\n" *
               "config_hash=$(String(config_hash))\n"
    return _sha256(collect(codeunits(material)))
end

function _parse_feature_bytes(snapshot::FileSnapshot, contract::ModelContract)
    bytes = snapshot.bytes
    isvalid(String, bytes) || _fail(:invalid_encoding, "feature table is not UTF-8")
    isempty(bytes) && _fail(:empty_input, "feature table is empty")
    last(bytes) == UInt8('\n') || _fail(:malformed_input, "feature table lacks final LF")
    UInt8('\r') in bytes && _fail(:malformed_input, "feature table contains CR")
    UInt8('\0') in bytes && _fail(:malformed_input, "feature table contains NUL")
    lines = split(String(copy(bytes)), '\n'; keepempty=true)
    pop!(lines) == "" || _fail(:malformed_input, "feature table termination differs")
    length(lines) > 1 || _fail(:empty_input, "feature table has no rows")
    any(isempty, lines) && _fail(:malformed_input, "feature table has a blank row")
    header = String.(split(first(lines), '\t'; keepempty=true))
    any(isempty, header) && _fail(:invalid_columns, "feature table has an empty column")
    length(header) == length(unique(header)) ||
        _fail(:duplicate_column, "feature table repeats a column")
    declared = unique(vcat(["amplitude"], reduce(vcat, last.(contract.views))))
    required = vcat(["file", "lobe"], declared)
    Set(header) == Set(required) && length(header) == length(required) ||
        _fail(:invalid_columns, "feature columns differ from the configured allowlist")
    rows = Dict{String,String}[]
    for (row_index, line) in enumerate(lines[2:end])
        values = String.(split(line, '\t'; keepempty=true))
        length(values) == length(header) ||
            _fail(:invalid_row_columns, "feature row $row_index has the wrong width")
        push!(rows, Dict(header[index] => values[index] for index in eachindex(header)))
    end
    SAF.build_feature_batch(rows, declared, declared;
                            required_metadata=["file", "lobe"])
    return rows
end

function _load_input(path::AbstractString, contract::ModelContract)
    snapshot = _snapshot(path)
    rows = _parse_feature_bytes(snapshot, contract)
    records = try
        HUA.load_records(snapshot.path)
    catch error
        _fail(:feature_load, sprint(showerror, error))
    end
    length(records) == length(rows) || _fail(:key_mismatch, "feature loader changed row count")
    row_keys = Set((row["file"], parse(Int, row["lobe"])) for row in rows)
    record_keys = Set((record.file, record.lobe) for record in records)
    row_keys == record_keys || _fail(:key_mismatch, "feature loader changed keys")
    _verify(snapshot)
    return InputBundle(records, snapshot)
end

function _placeholder_view(name::String, features::Vector{String}, n::Int,
                           covariance_floor::Float64; converged::Bool=false,
                           monotone::Bool=true, one_component::Bool=false,
                           degenerate::Bool=false)
    p = length(features)
    fit2 = HUA.TwoComponentFit(zeros(2, p), fill(covariance_floor, 2, p),
                               [0.5, 0.5], -Inf, Float64[], converged, monotone)
    fit1 = HUA.OneComponentFit(zeros(p), fill(covariance_floor, p), -Inf)
    return HUA.ViewFit(name, copy(features), fit2, fit1, Int[], 0,
                       one_component, degenerate, fill(NaN, n), fill(NaN, n))
end

function _strict_gaussian_orientation!(fit, records, features)
    fit.degenerate && return "degenerate_view", false
    !fit.fit2.converged && return "unstable_model", true
    !fit.fit2.monotone && return "nonmonotone_model", true
    fit.one_component_evidence && return "one_component_evidence", false
    X, _ = HUA.feature_matrix(records, features)
    normalized = HUA.normalize_per_scan(records, X, features)
    valid = fit.valid_lobe_idx
    isempty(valid) && return "degenerate_view", false
    responsibilities = HUA.responsibilities(fit.fit2, normalized[valid, :])
    hard = [argmax(view(responsibilities, index, :)) for index in axes(responsibilities, 1)]
    counts = [count(==(component), hard) for component in 1:2]
    if any(iszero, counts)
        fit.one_component_evidence = true
        fill!(fit.responsibility_for_label1, NaN)
        fill!(fit.log_odds_for_label1, NaN)
        return "collapsed_component", false
    end
    amplitudes = [records[valid[index]].amplitude for index in eachindex(valid)]
    means = [mean(amplitudes[hard .== component]) for component in 1:2]
    if !all(isfinite, means) || means[1] == means[2]
        fit.one_component_evidence = true
        fit.high_cluster = 0
        fill!(fit.responsibility_for_label1, NaN)
        fill!(fit.log_odds_for_label1, NaN)
        return "amplitude_orientation_undefined", false
    end
    return "ok", false
end

function _view_evidence(records, features::Vector{String})
    X, _ = HUA.feature_matrix(records, features)
    normalized = HUA.normalize_per_scan(records, X, features)
    observed = BitVector([
        all(isfinite, view(X, index, :)) for index in eachindex(records)
    ])
    normalized_valid = BitVector([
        all(isfinite, view(normalized, index, :)) for index in eachindex(records)
    ])
    return observed, normalized_valid, normalized
end

function _gaussian_views(records, contract::ModelContract)
    assessments = RequiredViewFit[]
    for (name, features) in contract.views
        observed, normalized_valid, _ = _view_evidence(records, features)
        fit, reason, failed = if !any(observed)
            (_placeholder_view(name, features, length(records),
                               contract.covariance_floor;
                               converged=true, degenerate=true),
             "missing_view", false)
        elseif !any(normalized_valid)
            (_placeholder_view(name, features, length(records),
                               contract.covariance_floor;
                               converged=true, degenerate=true),
             "invalid_required_view_normalization", false)
        else
            candidate = try
                HUA.fit_view(records, features, name;
                             first_seed=0, n_starts=contract.n_starts,
                             cov_floor=contract.covariance_floor,
                             max_iter=contract.max_iter, tol=contract.tolerance)
            catch
                nothing
            end
            if candidate === nothing
                (_placeholder_view(name, features, length(records),
                                   contract.covariance_floor),
                 "gaussian_fit_failure", true)
            else
                candidate_reason, candidate_failed =
                    _strict_gaussian_orientation!(candidate, records, features)
                (candidate, candidate_reason, candidate_failed)
            end
        end
        push!(assessments, RequiredViewFit(
            fit, reason, failed, observed, normalized_valid))
    end
    return assessments
end

function _student_view(records, features::Vector{String}, name::String,
                       contract::ModelContract)
    observed, normalized_valid, normalized = _view_evidence(records, features)
    if !any(observed)
        fit = _placeholder_view(
            name, features, length(records), contract.covariance_floor;
            converged=true, degenerate=true)
        return RequiredViewFit(fit, "missing_view", false,
                               observed, normalized_valid)
    elseif !any(normalized_valid)
        fit = _placeholder_view(
            name, features, length(records), contract.covariance_floor;
            converged=true, degenerate=true)
        return RequiredViewFit(fit, "invalid_required_view_normalization", false,
                               observed, normalized_valid)
    end
    valid = findall(normalized_valid)
    fit = try
        SRE.fit_student_t_two_component(
            normalized[valid, :], [records[index].amplitude for index in valid])
    catch
        placeholder = _placeholder_view(
            name, features, length(records), contract.covariance_floor)
        return RequiredViewFit(placeholder, "student_t_fit_failure", true,
                               observed, normalized_valid)
    end
    trace = fit.best_start > 0 ? copy(fit.traces[fit.best_start].accepted_loglik) : Float64[]
    valid_fit = fit.status == SRE.ROBUST_VALID
    failed = fit.status == SRE.ROBUST_FAILED
    selected_trace = fit.best_start > 0 ? fit.traces[fit.best_start] : nothing
    monotone = selected_trace === nothing ? !failed : selected_trace.monotone
    fit2 = HUA.TwoComponentFit(copy(fit.means), copy(fit.scales), [0.5, 0.5],
                               fit.loglik, trace, !failed, monotone)
    fit1 = HUA.OneComponentFit(copy(fit.one_component.mean),
                               copy(fit.one_component.scale),
                               fit.one_component.loglik)
    probability = fill(NaN, length(records))
    log_odds = fill(NaN, length(records))
    if valid_fit
        for (fit_index, record_index) in enumerate(valid)
            value = fit.responsibilities[fit_index, fit.high_amplitude_component]
            probability[record_index] = value
            bounded = clamp(value, 1.0e-12, 1.0 - 1.0e-12)
            log_odds[record_index] = log(bounded / (1.0 - bounded))
        end
    end
    reason = valid_fit ? "ok" : lowercase(string(fit.invalid_reason))
    view_fit = HUA.ViewFit(name, copy(features), fit2, fit1, valid,
                           fit.high_amplitude_component, !valid_fit, false,
                           probability, log_odds)
    return RequiredViewFit(view_fit, reason, failed, observed, normalized_valid)
end

function _student_views(records, contract::ModelContract)
    assessments = RequiredViewFit[]
    for (name, features) in contract.views
        push!(assessments, _student_view(records, features, name, contract))
    end
    return assessments
end

function _specific_reason(default::String, reasons::Vector{String})
    all(==("amplitude_orientation_undefined"), reasons) &&
        return "amplitude_orientation_undefined"
    all(reason -> occursin("one_component", reason), reasons) &&
        return "one_component_evidence"
    all(==("degenerate_view"), reasons) && return "degenerate_view"
    return default
end

function _required_view_blockers(assessments::Vector{RequiredViewFit}, n_records::Int)
    blockers = Union{Nothing,RequiredViewBlocker}[nothing for _ in 1:n_records]
    invalid = RequiredViewFit[
        assessment for assessment in assessments
        if any(assessment.observed) && assessment.reason != "ok"
    ]
    if !isempty(invalid)
        reasons = [assessment.reason for assessment in invalid]
        failed = any(assessment.failed for assessment in invalid)
        reason = if any(==("invalid_required_view_normalization"), reasons)
            "invalid_required_view_normalization"
        elseif any(==("amplitude_orientation_undefined"), reasons)
            "amplitude_orientation_undefined"
        elseif any(value -> occursin("one_component", value) ||
                            value == "collapsed_component", reasons)
            "one_component_evidence"
        elseif failed
            first(assessment.reason for assessment in invalid if assessment.failed)
        else
            "invalid_required_view_fit"
        end
        fill!(blockers, RequiredViewBlocker(reason, failed))
        return blockers
    end

    for index in 1:n_records
        for assessment in assessments
            if assessment.observed[index] && !assessment.normalized_valid[index]
                blockers[index] = RequiredViewBlocker(
                    "invalid_required_view_normalization", false)
                break
            elseif assessment.normalized_valid[index] &&
                   !isfinite(assessment.fit.log_odds_for_label1[index])
                blockers[index] = RequiredViewBlocker("invalid_required_view_fit", false)
                break
            end
        end
    end
    return blockers
end

function _prediction_rows(records, combined, model::String, reasons::Vector{String},
                           failures::Vector{Bool}, source::String, config::String,
                           input::String, trigger::String,
                           blockers::Vector{Union{Nothing,RequiredViewBlocker}})
    rows = UnaryPrediction[]
    for (index, record) in enumerate(records)
        value = combined.probability_1[index]
        views_used = combined.views_used[index]
        reason = _specific_reason(combined.invalid_reason[index], reasons)
        blocker = blockers[index]
        if blocker !== nothing && blocker.failed
            push!(rows, UnaryPrediction(
                record.file, record.lobe, "?", nothing, 0.0, model, "failed",
                blocker.reason, 0, "pending", "pending", "pending", trigger,
                source, config, input,
            ))
        elseif blocker !== nothing
            push!(rows, UnaryPrediction(
                record.file, record.lobe, "?", 0.5, 0.0, model, "abstained",
                blocker.reason, 0, "pending", "pending", "pending", trigger,
                source, config, input,
            ))
        elseif isfinite(value) && value != 0.5
            predicted = value > 0.5 ? "1" : "0"
            confidence = abs(2.0 * value - 1.0)
            push!(rows, UnaryPrediction(
                record.file, record.lobe, predicted, value, confidence, model,
                "ok", combined.invalid_reason[index], views_used,
                "pending", "pending", "pending", trigger, source, config, input,
            ))
        elseif isfinite(value)
            push!(rows, UnaryPrediction(
                record.file, record.lobe, "?", 0.5, 0.0, model, "abstained",
                "probability_tie", views_used, "pending", "pending", "pending",
                trigger, source, config, input,
            ))
        elseif any(failures)
            push!(rows, UnaryPrediction(
                record.file, record.lobe, "?", nothing, 0.0, model, "failed",
                reason, 0, "pending", "pending", "pending", trigger, source,
                config, input,
            ))
        else
            push!(rows, UnaryPrediction(
                record.file, record.lobe, "?", 0.5, 0.0, model, "abstained",
                reason, 0, "pending", "pending", "pending", trigger, source,
                config, input,
            ))
        end
    end
    return rows
end

function _one_component_rows(records, contract::ModelContract, model::String,
                             source::String, input::String, trigger::String)
    availability = zeros(Int, length(records))
    for (_, features) in contract.views
        X, _ = HUA.feature_matrix(records, features)
        normalized = HUA.normalize_per_scan(records, X, features)
        valid = Int[index for index in eachindex(records)
                    if all(isfinite, view(normalized, index, :))]
        isempty(valid) || HUA.fit_one_component(
            normalized[valid, :]; cov_floor=contract.covariance_floor)
        for index in valid
            availability[index] += 1
        end
    end
    return [UnaryPrediction(
        record.file, record.lobe, "?", 0.5, 0.0, model, "abstained",
        availability[index] == 0 ? "degenerate_view" : "one_component_evidence",
        0, "pending", "pending", "pending", trigger, source,
        contract.config_hash, input,
    ) for (index, record) in enumerate(records)]
end

function build_predictions(config_path::AbstractString, feature_path::AbstractString,
                           model::AbstractString)
    contract = load_contract(config_path)
    model_name = String(model)
    model_name in SUPPORTED_MODELS || _fail(:unsupported_model, "model is not configured")
    input = _load_input(feature_path, contract)
    source, source_snapshots = _source_bundle()
    source == LOADED_SOURCE_HASH ||
        _fail(:stale_source, "loaded source differs from current bytes")
    foreach(_verify, LOADED_SOURCE_SNAPSHOTS)
    trigger = trigger_hash(model_name, contract.config_hash)
    rows = if model_name == "gaussian_one"
        _one_component_rows(input.records, contract, model_name, source,
                            input.snapshot.sha256, trigger)
    else
        assessments = model_name == "gaussian_two" ?
            _gaussian_views(input.records, contract) :
            _student_views(input.records, contract)
        fits = [assessment.fit for assessment in assessments]
        reasons = [assessment.reason for assessment in assessments]
        failures = [assessment.failed for assessment in assessments]
        combined = HUA.combine_views(fits, length(input.records))
        blockers = _required_view_blockers(assessments, length(input.records))
        _prediction_rows(input.records, combined, model_name, reasons, failures,
                         source, contract.config_hash, input.snapshot.sha256, trigger,
                         blockers)
    end
    snapshots = vcat(source_snapshots, [contract.snapshot, input.snapshot])
    foreach(_verify, snapshots)
    return PredictionBundle(rows, snapshots)
end

function _field_text(value::Union{Nothing,Float64})
    value === nothing && return "NA"
    return @sprintf("%.17g", value)
end

function _validate_row(row::UnaryPrediction)
    row.predicted in ("0", "1", "?") || _fail(:invalid_output, "prediction is invalid")
    row.status in ("ok", "abstained", "failed") ||
        _fail(:invalid_output, "status is invalid")
    row.model in SUPPORTED_MODELS || _fail(:invalid_output, "model is invalid")
    0 <= row.views_used <= 4 || _fail(:invalid_output, "view count is invalid")
    row.segment_id == row.edge_status == row.edge_admission_hash == "pending" ||
        _fail(:invalid_output, "pre-admission metadata is invalid")
    for digest in (row.trigger_hash, row.source_hash, row.config_hash, row.input_hash)
        occursin(r"^[0-9a-f]{64}$", digest) ||
            _fail(:invalid_output, "output hash is malformed")
    end
    return nothing
end

function prediction_tsv_bytes(bundle::PredictionBundle)
    io = IOBuffer()
    println(io, join(OUTPUT_COLUMNS, '\t'))
    for row in bundle.rows
        _validate_row(row)
        fields = [
            row.file,
            string(row.lobe),
            row.predicted,
            _field_text(row.probability_1),
            @sprintf("%.17g", row.confidence),
            row.model,
            row.status,
            row.invalid_reason,
            string(row.views_used),
            row.segment_id,
            row.edge_status,
            row.edge_admission_hash,
            row.trigger_hash,
            row.source_hash,
            row.config_hash,
            row.input_hash,
        ]
        any(field -> occursin('\t', field) || occursin('\n', field), fields) &&
            _fail(:invalid_output, "output field contains a delimiter")
        println(io, join(fields, '\t'))
    end
    return take!(io)
end

function _rename_noreplace(source::String, destination::String)
    Sys.islinux() || _fail(
        :atomic_noreplace_unsupported,
        "atomic no-replace publication requires Linux renameat2",
    )
    result = try
        ccall(
            :renameat2,
            Cint,
            (Cint, Cstring, Cint, Cstring, Cuint),
            Cint(-100), source, Cint(-100), destination, Cuint(1),
        )
    catch error
        _fail(:atomic_noreplace_unsupported,
              "Linux renameat2 is unavailable: $(sprint(showerror, error))")
    end
    result == 0 && return nothing
    error_number = Base.Libc.errno()
    if error_number in (17, 21, 39)
        _fail(:output_collision, "output destination already has an owner")
    elseif error_number in (22, 38, 95)
        _fail(:atomic_noreplace_unsupported,
              "renameat2 RENAME_NOREPLACE is unsupported (errno=$error_number)")
    else
        _fail(:atomic_publish_failure,
              "renameat2 RENAME_NOREPLACE failed (errno=$error_number)")
    end
end

function write_prediction_tsv(path::AbstractString, bundle::PredictionBundle)
    destination = _absolute_path(path)
    parent = dirname(destination)
    isdir(parent) || _fail(:missing_output_parent, "output parent is absent")
    _reject_symlinks(parent)
    bytes = prediction_tsv_bytes(bundle)
    temporary, io = mktemp(parent; cleanup=false)
    committed = false
    try
        write(io, bytes)
        flush(io)
        close(io)
        foreach(_verify, bundle.snapshots)
        hook = _before_noreplace_hook[]
        hook === nothing || hook(temporary, destination)
        _rename_noreplace(temporary, destination)
        committed = true
    finally
        isopen(io) && close(io)
        !committed && ispath(temporary) && rm(temporary; force=true)
    end
    return destination
end

const LOADED_SOURCE_HASH, LOADED_SOURCE_SNAPSHOTS = _source_bundle()

end
