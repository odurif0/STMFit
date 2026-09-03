#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using SHA
using Test
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SCALE_REBIND_CORRECTION_EVIDENCE_DIR = joinpath(
    ROOT, ".omo", "evidence", "structured-label-free-unit-assignment",
    "t12", "scale-certification-rebind", "correction", "checked-plan-rebind",
)
include(joinpath(@__DIR__, "lib", "structured_assignment", "chain_inference.jl"))
using .StructuredChainInference
const SCI = StructuredChainInference

# The producer serializer is test-only fixture authority.  Product code remains
# independent of T10/T11 producer modules; only canonical report bytes are
# borrowed here.
module SerializerFixtureAuthority
include(joinpath(@__DIR__, "lib", "structured_assignment", "edge_admission.jl"))
end
const SEA = SerializerFixtureAuthority.StructuredEdgeAdmission

const T10_SOURCE_SHA = "285809e3711706d7a059ba2a0dd4139b44c1874dec79308e2762a63799aa954e"
const CANONICAL_T10_SOURCE_SHA = "285809e3711706d7a059ba2a0dd4139b44c1874dec79308e2762a63799aa954e"
const OBSOLETE_T10_SOURCE_SHA = "c29fe06c7fc0ba33c4e168169fec5ced5a7907bbfcec626c183e3ced044da0f9"
const T10_CANDIDATE_SHA = "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const T10_MODEL_SHA = "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const T10_GRID_SHA = "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e"
const T10_FORWARD_SHA = "f79258197cd123f833e0541647c4e7107044149deb558ba765fbeee0545737ad"
const T10_BACKWARD_SHA = "4f866a1e27289a0cf010a6ed4b5e0b92454dc55f97807d33cc95edb9e87182ba"

sha(bytes) = bytes2hex(sha256(bytes))
fmt(value) = @sprintf("%.17g", Float64(value))
hash_lines(lines) = sha(codeunits(join(String.(collect(lines)), "\n") * "\n"))

function write_text(path, text)
    mkpath(dirname(path))
    write(path, text)
    return path
end

function copy_authority(root)
    relatives = [
        ".omo/plans/structured-label-free-unit-assignment.md",
        "config/unit_assignment_structured_candidate.toml",
        "config/unit_assignment_structured_model.toml",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase3-t3-universe/correction/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase3-t3-universe/correction/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/correction/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/correction/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/claim-correction/ClaimCorrection.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/canonical-t11-source-bundle.bytes",
        "test/lib/structured_assignment/edge_admission.jl",
        "test/evaluate_structured_edge_admission.jl",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/canonical-t11-source-bundle.bytes",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/DoneClaim.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/review/AdversarialVerify.json",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/canonical-t11-source-bundle.bytes",
    ]
    for relative in relatives
        destination = joinpath(root, relative)
        mkpath(dirname(destination))
        cp(joinpath(ROOT, relative), destination; force=true)
    end
end

function key_hash(files)
    return hash_lines(vcat(["schema_version\tfile\tlobe"], ["1\t$file\t$lobe" for (file, lobe, _) in files]))
end

function segment_id(file, first_lobe, last_lobe)
    lo, hi = minmax(first_lobe, last_lobe)
    sha(codeunits("segment-v1\nfile=$(basename(file))\nfirst_lobe=$lo\nlast_lobe=$hi\ngrid_hash=$T10_GRID_SHA\n"))[1:16]
end

function make_todo10(root; dates=("20260101",))
    input = joinpath(root, "todo10")
    mkpath(input)
    files = ["$(date)_chain.sxm" for date in dates]
    nodes = [(file, lobe, t) for file in files for (lobe, t) in enumerate((1.0, 1.4, 1.8))]
    segments = Dict(file => segment_id(file, 1, 3) for file in files)
    feature = repeat("1", 64)
    forward = repeat("2", 64)
    backward = repeat("3", 64)
    node_rows = [join([file, string(lobe), fmt(t), segments[file], "connected",
                       lobe == 1 ? "start" : "eligible",
                       lobe == 3 ? "end" : "eligible",
                       feature, T10_MODEL_SHA, T10_SOURCE_SHA], '\t')
                 for (file, lobe, t) in nodes]
    edge_values = [(file, left, right, fwd, bwd)
                   for file in files for (left, right, fwd, bwd) in
                   [(1, 2, 0.40, -0.10), (2, 3, 0.20, 0.30)]]
    t_by_lobe = Dict(lobe => t for (_, lobe, t) in nodes if lobe <= 3)
    edge_rows = [join([file, string(left), string(right), fmt(t_by_lobe[left]),
                       fmt(t_by_lobe[right]), fmt(t_by_lobe[right] - t_by_lobe[left]), segments[file], segments[file],
                       "eligible", "none", fmt(fwd), fmt(bwd), T10_GRID_SHA,
                       forward, backward, feature, T10_MODEL_SHA, T10_SOURCE_SHA], '\t')
                 for (file, left, right, fwd, bwd) in edge_values]
    edge_bytes = Vector{UInt8}(codeunits(join(SCI._T10_EDGE_HEADER, '\t') * "\n" * join(edge_rows, "\n") * "\n"))
    node_bytes = Vector{UInt8}(codeunits(join(SCI._T10_NODE_HEADER, '\t') * "\n" * join(node_rows, "\n") * "\n"))
    edge_path = joinpath(input, "edge_observations.tsv")
    node_path = joinpath(input, "node_segments.tsv")
    write(edge_path, edge_bytes)
    write(node_path, node_bytes)
    universe_hash = repeat("a", 64)
    keys = key_hash(nodes)
    receipt_pairs = [
        "schema" => "stmfit-structured-edge-feature-receipt-v1",
        "schema_version" => 1,
        "status" => "PASS",
        "edge_observations_file" => "edge_observations.tsv",
        "node_segments_file" => "node_segments.tsv",
        "edge_header" => join(SCI._T10_EDGE_HEADER, '\t'),
        "node_header" => join(SCI._T10_NODE_HEADER, '\t'),
        "edge_header_sha256" => sha(codeunits(join(SCI._T10_EDGE_HEADER, '\t') * "\n")),
        "node_header_sha256" => sha(codeunits(join(SCI._T10_NODE_HEADER, '\t') * "\n")),
        "edge_observations_sha256" => sha(edge_bytes),
        "node_segments_sha256" => sha(node_bytes),
        "edge_row_count" => length(edge_rows),
        "node_row_count" => length(node_rows),
        "feature_sha256" => feature,
        "candidate_config_sha256" => T10_CANDIDATE_SHA,
        "model_config_sha256" => T10_MODEL_SHA,
        "universe_receipt_sha256" => universe_hash,
        "keys_sha256" => keys,
        "grid_sha256" => T10_GRID_SHA,
        "forward_patch_sha256" => forward,
        "backward_patch_sha256" => backward,
        "forward_patch_receipt_sha256" => repeat("4", 64),
        "backward_patch_receipt_sha256" => repeat("5", 64),
        "forward_producer_sha256" => T10_FORWARD_SHA,
        "backward_producer_sha256" => T10_BACKWARD_SHA,
        "forward_command" => "synthetic-forward",
        "backward_command" => "synthetic-backward",
        "source_sha256" => T10_SOURCE_SHA,
        "t2_review_sha256" => SCI._T2_REVIEW_SHA256,
        "t3_review_sha256" => SCI._T3_REVIEW_SHA256,
    ]
    receipt_text = join([key * " = " * (value isa Integer ? string(value) : repr(String(value))) for (key, value) in receipt_pairs], "\n") * "\n"
    receipt_path = write_text(joinpath(input, "receipt.toml"), receipt_text)
    return (edge=edge_path, node=node_path, receipt=receipt_path, file=first(files),
            files=files, dates=String.(dates), nodes=nodes, edges=edge_values,
            segment=segments[first(files)])
end

function matrix_hash(values)
    io = IOBuffer()
    println(io, length(values), 'x', 1)
    for value in values
        println(io, fmt(value))
    end
    return sha(take!(io))
end

function make_todo11_legacy(root, t10; shared_references=false, reversed=false, terminal_status="PASS")
    output = joinpath(root, "todo11")
    mkpath(output)
    fit_specs = shared_references ? [
        (repeat("d", 64), repeat("a", 64), "partition", "final_lodo", "NA", t10.dates[1],
         "partition|C1|final_lodo|NA|$(t10.dates[1])|$(repeat("d", 64));partition|C1|final_lodo|NA|$(t10.dates[2])|$(repeat("d", 64))"),
        (repeat("f", 64), repeat("b", 64), "full_refit", "full_refit", "NA", "NA",
         "full_refit|C1|full_refit|NA|NA|$(repeat("f", 64))"),
    ] : [
        (repeat("f", 64), repeat("a", 64), "full_refit", "full_refit", "NA", "NA",
         "full_refit|C1|full_refit|NA|NA|$(repeat("f", 64))"),
    ]
    fit_rows = String[]
    unary_rows = String[]
    transform_rows = String[]
    partition_rows = Vector{Vector{String}}()
    ess_rows = Vector{Vector{String}}()
    start_rows = Vector{Vector{String}}()
    score_rows = Vector{Vector{String}}()
    shuffle_rows = Vector{Vector{String}}()
    bootstrap_rows = Vector{Vector{String}}()
    p1_base = [0.2, 0.7, 0.4]
    node_order = hash_lines(vcat(["schema=structured-edge-todo10-node-order-v1"],
        [join((file, file[1:8], string(lobe), fmt(t)), '\t') for (file, lobe, t) in t10.nodes]))
    segment_by_file = Dict(file => segment_id(file, 1, 3) for file in t10.files)
    t_by_lobe = Dict(lobe => t for (_, lobe, t) in t10.nodes if lobe <= 3)
    for (fit_index, (fit, unary_fit, role, scope, outer_date, target_date, references)) in enumerate(fit_specs)
        model_dates = role == "full_refit" ? t10.dates : t10.dates[3:end]
        training_dates = join(model_dates, ',')
        p1s = [p1_base[mod1(index + fit_index - 1, 3)] for index in 1:length(t10.nodes)]
        probability = matrix_hash(p1s)
        residualizer_training = repeat(fit_index == 1 ? "e" : "c", 64)
        training_edge_count = length(t10.edges)
        residualizer = terminal_status == "PASS" ? hash_lines(vcat([
            "schema=structured-edge-residualizer-v1", "training_sha256=$residualizer_training",
            "training_row_count=$training_edge_count", "ridge=0",
            "evaluation_order=$(join(reversed ? vcat(1, collect(9:15), collect(2:8)) : collect(1:15), ','))",
        ], ["$(name)_$(output)=0" for name in SCI._GRAPH_COEFFICIENT_NAMES for output in SCI._GRAPH_OUTPUT_NAMES])) : "NA"
        model_row = fill("NA", length(SCI._GRAPH_MODEL_HEADER))
        positions = Dict(name => i for (i, name) in enumerate(SCI._GRAPH_MODEL_HEADER))
        put(name, value) = (model_row[positions[name]] = string(value))
        model_status = terminal_status
        model_reason = terminal_status == "PASS" ? "ok" : lowercase(terminal_status)
        put_values = (("schema_version", "1"), ("fit_role", role), ("fit_sha256", fit),
            ("model", "C1"), ("scope", scope), ("outer_date", outer_date),
            ("target_date", target_date), ("training_dates", training_dates),
            ("reversed", string(reversed)), ("status", model_status), ("reason", model_reason),
            ("training_sha256", repeat(Char('b' + fit_index - 1), 64)), ("unary_fit_sha256", unary_fit),
            ("reference_count", string(length(split(references, ';')))),
            ("references_sha256", hash_lines([join(split(reference, '|'), '\t') for reference in split(references, ';')])),
            ("references", references), ("unary_probability_sha256", probability),
            ("node_order_sha256", node_order), ("residualizer_sha256", residualizer),
            ("residualizer_training_sha256", terminal_status == "PASS" ? residualizer_training : "NA"),
            ("residualizer_ridge", terminal_status == "PASS" ? "0" : "NA"),
            ("student_t_nu", "8"), ("category_order", "00,01,10,11"), ("state_order", "0,1"),
            ("null_mean_fwd", terminal_status == "PASS" ? "0" : "NA"),
            ("null_mean_bwd", terminal_status == "PASS" ? "0" : "NA"),
            ("null_scale_ff", terminal_status == "PASS" ? "1" : "NA"),
            ("null_scale_fb", terminal_status == "PASS" ? "0" : "NA"),
            ("null_scale_bb", terminal_status == "PASS" ? "1" : "NA"),
            ("conditional_scale_ff", terminal_status == "PASS" ? "1" : "NA"),
            ("conditional_scale_fb", terminal_status == "PASS" ? "0" : "NA"),
            ("conditional_scale_bb", terminal_status == "PASS" ? "1" : "NA"),
            ("selected_start", "1"), ("selected_alpha", "0"), ("selected_objective", "0"),
            ("selected_converged", "true"))
        for (name, value) in put_values
            put(name, value)
        end
        means = [("00", 0.10, 0.00), ("01", 0.00, 0.00), ("10", 0.00, 0.00), ("11", -0.10, 0.00)]
        for (category, fwd, bwd) in means
            put("conditional_$(category)_mean_corr_fwd", terminal_status == "PASS" ? fmt(fwd) : "NA")
            put("conditional_$(category)_mean_corr_bwd", terminal_status == "PASS" ? fmt(bwd) : "NA")
        end
        for name in SCI._GRAPH_COEFFICIENT_NAMES, output_name in SCI._GRAPH_OUTPUT_NAMES
            put("$(name)_$(output_name)", terminal_status == "PASS" ? "0" : "NA")
        end
        push!(fit_rows, join(model_row, '\t'))
        for (index, (file, lobe, t)) in enumerate(t10.nodes)
            p1 = p1s[index]
            date = file[1:8]
            push!(unary_rows, join(["1", fit, "C1", file, date, string(lobe), fmt(t),
                join((file, date, string(lobe), fmt(t)), '|'), "0,1", fmt(1.0 - p1), fmt(p1), unary_fit, probability, node_order], '\t'))
        end
        reference_list = split(references, ';')
        canonical_reference = first(reference_list)
        reference = split(canonical_reference, '|')
        reference_scope, reference_outer, reference_target = reference[3:5]
        for (index, (file, left, right, _, _)) in enumerate(t10.edges)
            date = file[1:8]
            raw_hash = sha(codeunits(readlines(t10.edge)[index + 1] * "\n"))
            status = terminal_status == "PASS" ? "eligible" : "unavailable"
            reason = terminal_status == "PASS" ? "none" : "residualizer_unavailable"
            pred_fwd = terminal_status == "PASS" ? fmt(index % 2 == 0 ? 0.00 : 0.01) : "NA"
            pred_bwd = terminal_status == "PASS" ? fmt(index % 2 == 0 ? 0.03 : -0.02) : "NA"
            segment = segment_by_file[file]
            push!(transform_rows, join(["1", fit, "C1", reference_scope, reference_outer, reference_target, file,
                date, string(left), string(right), fmt(t_by_lobe[left]), fmt(t_by_lobe[right]), segment,
                string(index), join((file, string(left), string(right), segment), '|'), raw_hash,
                terminal_status == "PASS" ? residualizer : "NA", unary_fit, pred_fwd, pred_bwd, status, reason], '\t'))
        end
        for reference_text in reference_list
            reference = split(reference_text, '|')
            reference_scope, reference_outer, reference_target = reference[3:5]
            if role == "partition"
                push!(partition_rows, ["1", "C1", reference_scope, reference_outer, reference_target,
                    terminal_status, model_reason, training_dates, string(length(model_dates)), "1",
                     repeat(Char('b' + fit_index - 1), 64), fit, repeat(Char('d' + fit_index - 1), 64),
                     repeat(Char('e' + fit_index - 1), 64), repeat(Char('e' + fit_index - 1), 64), "true", "true"])
            end
        end
        for reference_text in reference_list
            reference = split(reference_text, '|')
            reference_scope, reference_outer, reference_target = reference[3:5]
            for category in 0:3
                push!(ess_rows, ["1", "C1", reference_scope, reference_outer, reference_target, string(category), string(length(model_dates)),
                    string(length(model_dates)), string(length(t10.edges)), "1", "1", "1", "1", "true", "ok", fit])
            end
            for start_role in ("null", "conditional"), start_index in (start_role == "null" ? (0,) : 1:5)
                alpha = start_role == "null" ? "NA" : ("0", "0.25", "0.5", "0.75", "1")[start_index]
                start_status = terminal_status
                numeric = terminal_status == "PASS" ? ("1", "0", "0", "1", "NA") : ("0", "NA", "NA", "0", "NA")
                push!(start_rows, ["1", "C1", reference_scope, reference_outer, reference_target, start_role,
                    string(start_index), alpha, start_status, terminal_status == "PASS" ? "ok" : lowercase(terminal_status),
                    numeric[1], numeric[2], numeric[3], numeric[4], numeric[5], repeat(Char('e' + fit_index - 1), 64), fit])
            end
            if role == "partition"
                score_hash = repeat(Char('d' + fit_index - 1), 64)
                target_file = "$(reference_target)_chain.sxm"
                target_edges = [(file, left, right) for (file, left, right, _, _) in t10.edges if file == target_file]
                for (file, left, right) in target_edges
                    push!(score_rows, ["1", "C1", reference_scope, reference_outer, reference_target, "edge",
                        file, string(left), string(right), segment_by_file[file], "1", "1", score_hash])
                end
                push!(score_rows, ["1", "C1", reference_scope, reference_outer, reference_target, "scan",
                    target_file, "NA", "NA", "NA", "1", string(length(target_edges)), score_hash])
                push!(score_rows, ["1", "C1", reference_scope, reference_outer, reference_target, "date",
                    "NA", "NA", "NA", "NA", "1", string(length(target_edges)), score_hash])
                for seed in 0:499
                    push!(shuffle_rows, ["1", "C1", reference_scope, reference_outer, reference_target, string(seed), "1",
                        terminal_status == "PASS" ? unary_fit : "not_run", "1", "1", terminal_status == "PASS" ? "true" : "false", score_hash])
                end
            end
        end
    end
    gate_hash = SCI._diag_gate_hash("PASS", "ok", String[], Dict{String,Float64}(), fill(1.0, 500), 1.0, true, true, true)
    for seed in 0:499
        push!(bootstrap_rows, ["1", "C1", "final_lodo", "NA", string(seed), "1", "1", "true", "true",
            "true", terminal_status, terminal_status == "PASS" ? "ok" : lowercase(terminal_status), gate_hash])
    end
    model_bytes = Vector{UInt8}(codeunits(join(SCI._GRAPH_MODEL_HEADER, '\t') * "\n" * join(fit_rows, "\n") * "\n"))
    unary_bytes = Vector{UInt8}(codeunits(join(SCI._GRAPH_UNARY_HEADER, '\t') * "\n" * join(unary_rows, "\n") * "\n"))
    transform_bytes = Vector{UInt8}(codeunits(join(SCI._GRAPH_EDGE_HEADER, '\t') * "\n" * join(transform_rows, "\n") * "\n"))
    final_gate = gate_hash
    full_fit = shared_references ? repeat("f", 64) : repeat("f", 64)
    overall_status = terminal_status
    overall_reason = terminal_status == "PASS" ? "ok" : terminal_status == "FAIL" ? "numerical_failure" : "conditional_mechanism_not_admitted"
    model_result = hash_lines([
        "schema=structured-edge-model-admission-v1", "model=C1", "status=$(overall_status)",
        "reason=$(overall_reason)", "final_gate=$(final_gate)", "full_fit=$(full_fit)",
    ])
    normalization = repeat("d", 64)
    hashes = Dict("plan_sha256" => SCI._PLAN_SHA256, "t7_review_sha256" => SCI._T7_REVIEW_SHA256,
        "t10_review_sha256" => SCI._T10_REVIEW_FILE_SHA256, "candidate_config_sha256" => T10_CANDIDATE_SHA,
        "model_config_sha256" => T10_MODEL_SHA, "universe_receipt_sha256" => repeat("a", 64),
        "universe_sha256" => repeat("b", 64), "universe_keys_sha256" => key_hash(t10.nodes),
        "edge_receipt_sha256" => sha(read(t10.receipt)), "edge_observations_sha256" => sha(read(t10.edge)),
        "node_segments_sha256" => sha(read(t10.node)), "edge_source_sha256" => T10_SOURCE_SHA,
        "edge_feature_sha256" => repeat("1", 64), "t3_review_sha256" => SCI._T3_REVIEW_SHA256,
        "input_sha256" => repeat("1", 64), "source_sha256" => SCI._T11_SOURCE_BUNDLE_SHA256,
        "normalization_sha256" => normalization)
    result_lines = ["schema=structured-edge-admission-report-v1", "status=$(overall_status)",
        "reason=$(overall_reason)", "normalization=$(normalization)", "model=C1:$(model_result)"]
    for key in sort!(collect(keys(hashes)))
        key == "normalization_sha256" && continue
        push!(result_lines, "input=$(key):$(hashes[key])")
    end
    overall_result = hash_lines(result_lines)
    artifacts = Dict{String,Vector{UInt8}}(
        "fitted_edge_models.tsv" => model_bytes, "fitted_unary_nodes.tsv" => unary_bytes,
        "fitted_edge_transforms.tsv" => transform_bytes,
        "models.tsv" => Vector{UInt8}(codeunits(join(SCI._MODEL_RESULT_HEADER, '\t') * "\n" *
            join(["1", "C1", overall_status, overall_reason, shared_references ? "0" : "0", final_gate, full_fit, model_result], '\t') * "\n")),
        "partitions.tsv" => Vector{UInt8}(codeunits(join(SCI._PARTITION_HEADER, '\t') * "\n" * (isempty(partition_rows) ? "" : join([join(row, '\t') for row in partition_rows], "\n") * "\n"))),
        "ess.tsv" => Vector{UInt8}(codeunits(join(SCI._ESS_HEADER, '\t') * "\n" * (isempty(ess_rows) ? "" : join([join(row, '\t') for row in ess_rows], "\n") * "\n"))),
        "starts.tsv" => Vector{UInt8}(codeunits(join(SCI._START_HEADER, '\t') * "\n" * join([join(row, '\t') for row in start_rows], "\n") * "\n")),
        "scores.tsv" => Vector{UInt8}(codeunits(join(SCI._SCORE_HEADER, '\t') * "\n" * (isempty(score_rows) ? "" : join([join(row, '\t') for row in score_rows], "\n") * "\n"))),
        "shuffle.tsv" => Vector{UInt8}(codeunits(join(SCI._SHUFFLE_HEADER, '\t') * "\n" * (isempty(shuffle_rows) ? "" : join([join(row, '\t') for row in shuffle_rows], "\n") * "\n"))),
        "bootstrap.tsv" => Vector{UInt8}(codeunits(join(SCI._BOOTSTRAP_HEADER, '\t') * "\n" * join([join(row, '\t') for row in bootstrap_rows], "\n") * "\n")),
    )
    for (name, bytes) in artifacts
        write(joinpath(output, name), bytes)
    end
    row_counts = Dict(name => length(split(chomp(String(copy(bytes))), '\n')) - 1 for (name, bytes) in artifacts)
    artifact_hashes = Dict(name => sha(bytes) for (name, bytes) in artifacts)
    fit_hashes = sort([spec[1] for spec in fit_specs])
    graph = Dict("fit_hashes" => fit_hashes, "full_refit_hashes" => [full_fit],
        "model_row_count" => row_counts["fitted_edge_models.tsv"],
        "unary_node_row_count" => row_counts["fitted_unary_nodes.tsv"],
        "edge_transform_row_count" => row_counts["fitted_edge_transforms.tsv"],
        "schemas" => Dict(name => SCI._GRAPH_HANDOFF_SCHEMA for name in SCI._GRAPH_HANDOFF_FILES))
    replay_io = IOBuffer()
    write(replay_io, codeunits("schema=structured-edge-graph-handoff-replay-v1\n"))
    for name in SCI._GRAPH_HANDOFF_FILES
        bytes = artifacts[name]
        write(replay_io, codeunits("file=$(name)\nlength=$(length(bytes))\n")); write(replay_io, bytes); write(replay_io, UInt8('\n'))
    end
    receipt = Dict("schema" => "stmfit-structured-edge-admission-receipt-v2", "schema_version" => 2,
        "status" => overall_status, "reason" => overall_reason, "result_sha256" => overall_result,
        "model_count" => 1, "bootstrap_count" => 500, "shuffle_count" => 500, "fixed_nu" => 8,
        "node_priors" => [0.5, 0.5], "pair_prior" => [0.0, 0.0, 0.0, 0.0],
        "graph_handoff_schema" => SCI._GRAPH_HANDOFF_SCHEMA,
        "graph_handoff_replay_sha256" => sha(take!(replay_io)), "graph_handoff_files" => SCI._GRAPH_HANDOFF_FILES,
        "category_order" => collect(SCI._GRAPH_CATEGORY_ORDER), "state_order" => collect(SCI._GRAPH_STATE_ORDER),
        "coefficient_order" => ["$(name)_$(output)" for name in SCI._GRAPH_COEFFICIENT_NAMES for output in SCI._GRAPH_OUTPUT_NAMES],
        "hashes" => hashes,
        "models" => Dict("C1" => "$(overall_status)|$(overall_reason)|$model_result"),
        "row_counts" => row_counts, "graph_handoff" => graph, "artifacts" => artifact_hashes)
    receipt_path = write_text(joinpath(output, "receipt.toml"), sprint(TOML.print, receipt))
    return (models=joinpath(output, "fitted_edge_models.tsv"), unary=joinpath(output, "fitted_unary_nodes.tsv"),
        transforms=joinpath(output, "fitted_edge_transforms.tsv"), receipt=receipt_path, output=output,
        fit=full_fit, partition_fit=shared_references ? repeat("d", 64) : full_fit)
end

function _serializer_admission_data(t10)
    nodes = SEA.AdmissionNode[]
    for (file, lobe, t) in t10.nodes
        push!(nodes, SEA.AdmissionNode(file, file[1:8], lobe, Float64(t), 1.0,
                                       ntuple(index -> 0.1 * index, 7)))
    end
    edges = SEA.AdmissionEdge[]
    raw_lines = readlines(t10.edge)
    for (index, (file, left, right, fwd, bwd)) in enumerate(t10.edges)
        left_t = only(t for (edge_file, lobe, t) in t10.nodes if edge_file == file && lobe == left)
        right_t = only(t for (edge_file, lobe, t) in t10.nodes if edge_file == file && lobe == right)
        segment = segment_id(file, 1, 3)
        push!(edges, SEA.AdmissionEdge(file, file[1:8], left, right, Float64(left_t),
                                       Float64(right_t), segment, index,
                                       (Float64(fwd), Float64(bwd)),
                                       sha(codeunits(raw_lines[index + 1] * "\n"))))
    end
    hashes = Dict(
        "plan_sha256" => SCI._PLAN_SHA256,
        "t7_review_sha256" => SCI._T7_REVIEW_SHA256,
        "t10_review_sha256" => SCI._T10_REVIEW_FILE_SHA256,
        "candidate_config_sha256" => T10_CANDIDATE_SHA,
        "model_config_sha256" => T10_MODEL_SHA,
        "universe_receipt_sha256" => repeat("a", 64),
        "universe_sha256" => repeat("b", 64),
        "universe_keys_sha256" => key_hash(t10.nodes),
        "edge_receipt_sha256" => sha(read(t10.receipt)),
        "edge_observations_sha256" => sha(read(t10.edge)),
        "node_segments_sha256" => sha(read(t10.node)),
        "edge_source_sha256" => T10_SOURCE_SHA,
        "edge_feature_sha256" => repeat("1", 64),
        "t3_review_sha256" => SCI._T3_REVIEW_SHA256,
        "input_sha256" => repeat("1", 64),
        "source_sha256" => SCI._T11_SOURCE_BUNDLE_SHA256,
        "normalization_sha256" => repeat("d", 64),
    )
    return SEA.AdmissionData(nodes, edges, hashes)
end

function _serializer_fit(t10, data, fit_sha, training_dates, status; reversed=false)
    training_dates = collect(String.(training_dates))
    training_files = sort(unique(node.file for node in data.nodes if node.date in training_dates))
    identities = [(node.file, node.date, node.lobe, node.t_nm) for node in data.nodes]
    probabilities = fill(0.8, length(identities))
    node_order_sha = SEA._node_order_sha256(identities)
    probability_sha = SEA._matrix_hash(reshape(probabilities, :, 1))
    unary = SEA.UnaryFitResult(status, status == "PASS" ? :ok : :insufficient_support,
                               "C1", probabilities, SEA.UnaryViewModel[],
                               training_dates, training_files, identities,
                               node_order_sha, repeat("a", 64), probability_sha)
    training_edges = [edge for edge in data.edges if edge.date in training_dates]
    training_sha = repeat("c", 64)
    residualizer = nothing
    null_fit = nothing
    conditional_fit = nothing
    if status == "PASS"
        coefficients = zeros(15, 2)
        evaluation_order = reversed ? vcat(1, collect(9:15), collect(2:8)) : collect(1:15)
        residualizer = SEA.ResidualizerFit(0.0, coefficients, length(training_edges),
                                           training_sha, evaluation_order)
        null_fit = SEA.NullFit(:ok, :ok, [1.0, 1.0], Matrix{Float64}(I, 2, 2),
                               0.0, [0.0], true, 1, 0.0)
        starts = SEA.ConditionalStartTrace[]
        means = fill(0.0, 4, 2)
        means[:, 1] .= 0.30
        means[:, 2] .= 0.10
        scale = Matrix{Float64}(I, 2, 2)
        for index in 1:5
            push!(starts, SEA.ConditionalStartTrace(index,
                (0.0, 0.25, 0.5, 0.75, 1.0)[index], [0.0], true, true,
                true, :ok, 0, 0.0, means, scale))
        end
        conditional_fit = SEA.ConditionalFit(:ok, :ok, means, scale, 0.0,
                                              starts, 1, true)
    end
    samples = SEA.EdgeSample[]
    transforms = SEA.EdgeTransformState[]
    for edge in data.edges
        predictors = zeros(14)
        transform_status = status == "PASS" ? "eligible" : "unavailable"
        transform_reason = status == "PASS" ? :none : :residualizer_unavailable
        status == "PASS" && push!(samples, SEA.EdgeSample(edge, predictors,
                                      collect(edge.observation),
                                      (0.04, 0.16, 0.16, 0.64), edge.raw_row_sha256))
        push!(transforms, SEA.EdgeTransformState(edge, predictors,
                                                  transform_status, transform_reason))
    end
    support = status == "PASS" ?
        SEA.SupportEvidence(length(training_dates), length(training_files), length(training_edges),
                            fill(Float64(length(training_edges)), 4),
                            fill(isempty(training_dates) ? 0 : 1, 4), 9,
                            length(training_edges) / 9.0, true, :ok) :
        SEA.SupportEvidence(0, 0, 0, fill(0.0, 4), fill(0, 4), 9, 0.0, false,
                            status == "FAIL" ? :numerical_failure : :insufficient_support)
    return SEA.PartitionFitResult(status,
        status == "PASS" ? :ok : status == "FAIL" ? :numerical_failure : :insufficient_support,
        "C1", training_dates, reversed, unary, residualizer, null_fit,
        conditional_fit, samples, transforms, support, training_sha, fit_sha)
end

function _serializer_score_hash(fit_sha, target, date_mean, control_mean,
                                shuffle_quantile, shuffle_pass, scan_scores,
                                gains, shuffle_values)
    lines = [
        "schema=structured-edge-partition-score-v1",
        "fit=$(fit_sha)",
        "target=$(target)",
        "date_mean=$(SEA._format_float(date_mean))",
        "control_mean=$(SEA._format_float(control_mean))",
        "shuffle_quantile=$(SEA._format_float(shuffle_quantile))",
        "shuffle_pass=$(shuffle_pass)",
    ]
    append!(lines, "scan=$(scan):$(SEA._format_float(scan_scores[scan]))" for scan in sort!(collect(keys(scan_scores))))
    append!(lines, "gain=$(identity):$(SEA._format_float(gains[identity]))" for identity in sort!(collect(keys(gains))))
    append!(lines, "shuffle=$(SEA._format_float(value))" for value in shuffle_values)
    return SEA._hash_lines(lines)
end

function _serializer_evaluation(fit, data, target; run_shuffle=true)
    fit.status != "PASS" && return SEA._invalid_evaluation(fit.status, fit.reason, fit, target)
    evaluation = SEA.evaluate_partition(fit, target; run_shuffle=run_shuffle)
    !run_shuffle && return evaluation
    shuffle_values = zeros(500)
    shuffle_quantile = 0.0
    shuffle_pass = evaluation.control_date_mean > shuffle_quantile
    score_hash = _serializer_score_hash(
        fit.fit_sha256, target, evaluation.date_mean, evaluation.control_date_mean,
        shuffle_quantile, shuffle_pass, evaluation.scan_scores, evaluation.gains,
        shuffle_values,
    )
    return SEA.PartitionEvaluation(
        evaluation.status, evaluation.reason, evaluation.model_id,
        evaluation.target_date, evaluation.training_dates, evaluation.fit,
        evaluation.scan_scores, evaluation.date_mean, evaluation.control_scan_scores,
        evaluation.control_date_mean, shuffle_values, evaluation.shuffle_fit_sha256,
        shuffle_quantile, shuffle_pass, evaluation.reversal_pass, evaluation.gains,
        score_hash,
    )
end

function _serializer_model_result(outer_folds, final_gate, full_fit)
    states = String[]
    for fold in outer_folds
        push!(states, fold.inner_gate.status, fold.outer_score.status)
    end
    push!(states, final_gate.status, full_fit.status)
    status = SEA.combine_terminal_states(states)
    reason = status == "PASS" ? :ok : status == "FAIL" ? :numerical_failure :
             status == "BLOCKED" ? :integrity_failure : :conditional_mechanism_not_admitted
    lines = [
        "schema=structured-edge-model-admission-v1",
        "model=C1",
        "status=$(status)",
        "reason=$(String(reason))",
        "final_gate=$(final_gate.gate_sha256)",
        "full_fit=$(full_fit.fit_sha256)",
    ]
    append!(lines, join(("outer", fold.outer_date, fold.inner_gate.gate_sha256,
                         fold.outer_score.score_sha256, fold.fit_sha256), '\t')
            for fold in outer_folds)
    return SEA.ModelAdmissionResult(status, reason, "C1", outer_folds,
                                    final_gate, full_fit, SEA._hash_lines(lines))
end

function make_todo11_serializer(root, t10; reversed=false, terminal_status="PASS")
    output = joinpath(root, "todo11")
    mkpath(output)
    data = _serializer_admission_data(t10)
    dates = sort(t10.dates)
    fit_status = terminal_status
    training_sets = Dict{String,Vector{String}}()
    fit_hashes = Dict{String,String}()
    fit_index = 0
    subsets = Vector{String}[]
    for i in 1:4, j in (i + 1):4
        push!(subsets, [dates[i], dates[j]])
    end
    for heldout in 1:4
        push!(subsets, [dates[i] for i in 1:4 if i != heldout])
    end
    push!(subsets, collect(dates))
    for training_dates in subsets
        key = join(training_dates, ',')
        haskey(training_sets, key) && continue
        fit_index += 1
        training_sets[key] = training_dates
        fit_hashes[key] = repeat(string("0123456789abcdef"[fit_index]), 64)
    end
    fits = Dict{Tuple{String,Bool},SEA.PartitionFitResult}()
    for (key, training_dates) in training_sets
        fits[(key, false)] = _serializer_fit(t10, data, fit_hashes[key], training_dates,
                                             fit_status; reversed=false)
        fits[(key, true)] = _serializer_fit(t10, data, fit_hashes[key], training_dates,
                                            fit_status; reversed=true)
    end
    fit_for(training_dates, reversed_fit=false) =
        fits[(join(sort(training_dates), ','), reversed_fit)]
    outer_folds = SEA.OuterFoldEvidence[]
    for outer_date in dates
        outer_training = [date for date in dates if date != outer_date]
        inner_evaluations = SEA.PartitionEvaluation[]
        for target_date in outer_training
            inner_training = [date for date in outer_training if date != target_date]
            observed = _serializer_evaluation(fit_for(inner_training), data, target_date;
                                               run_shuffle=true)
            reversed_score = _serializer_evaluation(fit_for(inner_training, true), data,
                                                    target_date; run_shuffle=false)
            push!(inner_evaluations, SEA._with_reversal(observed, reversed_score))
        end
        inner_gate = SEA.evaluate_gate(inner_evaluations)
        observed_outer = _serializer_evaluation(fit_for(outer_training), data, outer_date;
                                                run_shuffle=false)
        reversed_outer = _serializer_evaluation(fit_for(outer_training, true), data,
                                                outer_date; run_shuffle=false)
        outer_score = SEA._with_reversal(observed_outer, reversed_outer)
        push!(outer_folds, SEA.OuterFoldEvidence(outer_date, inner_gate, outer_score,
                                                 fit_for(outer_training).fit_sha256))
    end
    final_evals = SEA.PartitionEvaluation[]
    for target_date in dates
        training_dates = [date for date in dates if date != target_date]
        observed = _serializer_evaluation(fit_for(training_dates), data, target_date;
                                           run_shuffle=true)
        reversed_score = _serializer_evaluation(fit_for(training_dates, true), data,
                                                target_date; run_shuffle=false)
        push!(final_evals, SEA._with_reversal(observed, reversed_score))
    end
    final_gate = SEA.evaluate_gate(final_evals)
    model_result = _serializer_model_result(outer_folds, final_gate, fit_for(dates))
    top_status = model_result.status
    top_reason = model_result.reason
    result_lines = [
        "schema=structured-edge-admission-report-v1",
        "status=$(top_status)",
        "reason=$(String(top_reason))",
        "normalization=$(data.hashes["normalization_sha256"])",
        "model=C1:$(model_result.result_sha256)",
    ]
    for key in sort!(collect(keys(data.hashes)))
        key == "normalization_sha256" && continue
        push!(result_lines, "input=$(key):$(data.hashes[key])")
    end
    report = SEA.AdmissionReport(top_status, top_reason, Dict("C1" => model_result),
                                 data.hashes, SEA._hash_lines(result_lines))
    files = SEA.report_files(report)
    for (name, bytes) in files
        write(joinpath(output, name), bytes)
    end
    return (models=joinpath(output, "fitted_edge_models.tsv"),
            unary=joinpath(output, "fitted_unary_nodes.tsv"),
            transforms=joinpath(output, "fitted_edge_transforms.tsv"),
            receipt=joinpath(output, "receipt.toml"), output=output,
            fit=fit_for(dates).fit_sha256,
            partition_fit=fit_for(dates[3:4]).fit_sha256)
end

function make_todo11(root, t10; shared_references=false, reversed=false, terminal_status="PASS")
    shared_references && return make_todo11_serializer(root, t10;
                                                       reversed=reversed,
                                                       terminal_status=terminal_status)
    return make_todo11_legacy(root, t10; shared_references=shared_references,
                               reversed=reversed, terminal_status=terminal_status)
end

function make_fixture(; dates=("20260101", "20260102", "20260103", "20260104"), shared_references=true, reversed=false, terminal_status="PASS")
    root = mktempdir(; prefix="stmfit-structured-chain-")
    copy_authority(root)
    t10 = make_todo10(root; dates=dates)
    effective_shared = shared_references && length(dates) >= 4
    t11 = make_todo11(root, t10; shared_references=effective_shared, reversed=reversed, terminal_status=terminal_status)
    return (root=root, t10=t10, t11=t11)
end

function make_producer_todo10(root; dates=[string("2026010", index) for index in 1:7])
    masks = (52, 11, 38, 25)
    feature = repeat("1", 64)
    forward = repeat("2", 64)
    backward = repeat("3", 64)
    source = SCI._T10_SOURCE_SHA256
    nodes = SEA.AdmissionNode[]
    for date in dates, scan in 1:4
        file = string(date, "_scan", lpad(scan, 2, '0'), ".sxm")
        date_index = parse(Int, date[8:8])
        for lobe in 1:6
            mode = (masks[scan] >> (lobe - 1)) & 1
            signal = mode == 1 ? 2.3 : -2.3
            predictors = ntuple(feature_index ->
                signal * (1.0 + 0.025 * feature_index) +
                0.17 * sin(0.4 * date_index + 0.8 * scan) +
                0.11 * sin(0.31 * lobe + 0.17 * feature_index + 0.13 * scan), 7)
            amplitude = 1.0 + 2.0 * mode + 0.03 * cos(0.5 * lobe + scan)
            push!(nodes, SEA.AdmissionNode(file, date, lobe, 0.5 * (lobe - 1),
                                           amplitude, predictors))
        end
    end
    sort!(nodes; by=node -> (node.file, node.t_nm, node.lobe))
    normalized = SEA.normalize_admission_data(
        SEA.AdmissionData(nodes, SEA.AdmissionEdge[], Dict{String,String}()),
    )
    conditional_means = Dict(
        (0, 0) => [-0.31, 0.12],
        (0, 1) => [0.34, -0.10],
        (1, 0) => [0.34, -0.10],
        (1, 1) => [-0.25, 0.15],
    )
    coefficients = zeros(15, 2)
    coefficients[1, :] .= (0.035, -0.025)
    for row in 2:15
        coefficients[row, 1] = 0.0025 * sin(0.37 * row)
        coefficients[row, 2] = 0.0020 * cos(0.29 * row)
    end
    edge_rows = String[]
    edges = SEA.AdmissionEdge[]
    for file in sort!(unique(node.file for node in nodes))
        chain = filter(node -> node.file == file, nodes)
        segment = segment_id(file, 1, 6)
        scan = parse(Int, match(r"scan(\d+)", file).captures[1])
        for ordinal in 1:5
            left, right = chain[ordinal], chain[ordinal + 1]
            left_index = normalized.node_index[(file, left.lobe)]
            right_index = normalized.node_index[(file, right.lobe)]
            endpoint = vcat(collect(normalized.normalized_predictors[left_index, :]),
                            collect(normalized.normalized_predictors[right_index, :]))
            baseline = vec(vcat(1.0, endpoint)' * coefficients)
            left_mode = (masks[scan] >> (left.lobe - 1)) & 1
            right_mode = (masks[scan] >> (right.lobe - 1)) & 1
            observation = baseline .+ conditional_means[(left_mode, right_mode)] .+
                [0.018 * sin(0.71 * ordinal + 0.19 * scan),
                 0.017 * cos(0.53 * ordinal + 0.23 * scan)]
            fields = [file, string(left.lobe), string(right.lobe), fmt(left.t_nm),
                      fmt(right.t_nm), fmt(right.t_nm - left.t_nm), segment, segment,
                      "eligible", "none", fmt(observation[1]), fmt(observation[2]),
                      SCI._GRID_SHA256, forward, backward, feature, SCI._MODEL_SHA256,
                      source]
            row_bytes = codeunits(join(fields, '\t') * "\n")
            push!(edge_rows, join(fields, '\t'))
            push!(edges, SEA.AdmissionEdge(file, left.date, left.lobe, right.lobe,
                                           left.t_nm, right.t_nm, segment, ordinal,
                                           (observation[1], observation[2]), sha(row_bytes)))
        end
    end
    sort!(edges; by=edge -> (edge.file, edge.left_t_nm, edge.right_t_nm,
                             edge.left_lobe, edge.right_lobe))
    node_rows = [join([node.file, string(node.lobe), fmt(node.t_nm),
                       segment_id(node.file, 1, 6), "connected",
                       node.lobe == 1 ? "start" : "eligible",
                       node.lobe == 6 ? "end" : "eligible", feature,
                       SCI._MODEL_SHA256, source], '\t') for node in nodes]
    edge_bytes = Vector{UInt8}(codeunits(join(SCI._T10_EDGE_HEADER, '\t') * "\n" *
                                           join(edge_rows, "\n") * "\n"))
    node_bytes = Vector{UInt8}(codeunits(join(SCI._T10_NODE_HEADER, '\t') * "\n" *
                                           join(node_rows, "\n") * "\n"))
    input = joinpath(root, "todo10")
    mkpath(input)
    edge_path = joinpath(input, "edge_observations.tsv")
    node_path = joinpath(input, "node_segments.tsv")
    write(edge_path, edge_bytes)
    write(node_path, node_bytes)
    keys = key_hash([(node.file, node.lobe, node.t_nm) for node in nodes])
    receipt_pairs = [
        "schema" => "stmfit-structured-edge-feature-receipt-v1",
        "schema_version" => 1, "status" => "PASS",
        "edge_observations_file" => "edge_observations.tsv",
        "node_segments_file" => "node_segments.tsv",
        "edge_header" => join(SCI._T10_EDGE_HEADER, '\t'),
        "node_header" => join(SCI._T10_NODE_HEADER, '\t'),
        "edge_header_sha256" => sha(codeunits(join(SCI._T10_EDGE_HEADER, '\t') * "\n")),
        "node_header_sha256" => sha(codeunits(join(SCI._T10_NODE_HEADER, '\t') * "\n")),
        "edge_observations_sha256" => sha(edge_bytes),
        "node_segments_sha256" => sha(node_bytes),
        "edge_row_count" => length(edge_rows), "node_row_count" => length(node_rows),
        "feature_sha256" => feature, "candidate_config_sha256" => T10_CANDIDATE_SHA,
        "model_config_sha256" => T10_MODEL_SHA, "universe_receipt_sha256" => repeat("a", 64),
        "keys_sha256" => keys, "grid_sha256" => SCI._GRID_SHA256,
        "forward_patch_sha256" => forward, "backward_patch_sha256" => backward,
        "forward_patch_receipt_sha256" => repeat("4", 64),
        "backward_patch_receipt_sha256" => repeat("5", 64),
        "forward_producer_sha256" => SCI._FORWARD_PRODUCER_SHA256,
        "backward_producer_sha256" => SCI._BACKWARD_PRODUCER_SHA256,
        "forward_command" => "synthetic-forward", "backward_command" => "synthetic-backward",
        "source_sha256" => source, "t2_review_sha256" => SCI._T2_REVIEW_SHA256,
        "t3_review_sha256" => SCI._T3_REVIEW_SHA256,
    ]
    receipt_path = write_text(joinpath(input, "receipt.toml"),
        join([key * " = " * (value isa Integer ? string(value) : repr(String(value)))
              for (key, value) in receipt_pairs], "\n") * "\n")
    hashes = Dict(
        "plan_sha256" => SCI._PLAN_SHA256, "t7_review_sha256" => SCI._T7_REVIEW_SHA256,
        "t10_review_sha256" => SCI._T10_REVIEW_FILE_SHA256,
        "candidate_config_sha256" => T10_CANDIDATE_SHA, "model_config_sha256" => T10_MODEL_SHA,
        "universe_receipt_sha256" => repeat("a", 64), "universe_sha256" => repeat("b", 64),
        "universe_keys_sha256" => keys, "edge_receipt_sha256" => sha(read(receipt_path)),
        "edge_observations_sha256" => sha(edge_bytes), "node_segments_sha256" => sha(node_bytes),
        "edge_source_sha256" => source, "edge_feature_sha256" => feature,
        "t3_review_sha256" => SCI._T3_REVIEW_SHA256, "input_sha256" => feature,
        "source_sha256" => SCI._T11_SOURCE_BUNDLE_SHA256,
    )
    data = SEA.AdmissionData(nodes, edges, hashes)
    return (edge=edge_path, node=node_path, receipt=receipt_path, dates=dates,
            files=sort!(unique(node.file for node in nodes)), nodes=nodes, edges=edges,
            data=data)
end

function make_producer_todo11(root, t10; require_pass=true)
    output = joinpath(root, "todo11")
    mkpath(output)
    report = SEA.evaluate_admission(t10.data)
    require_pass && report.status == "PASS" ||
        (!require_pass || error("producer fixture did not reach PASS: $(report.status)"))
    files = SEA.report_files(report)
    for (name, bytes) in files
        write(joinpath(output, name), bytes)
    end
    return (models=joinpath(output, "fitted_edge_models.tsv"),
            unary=joinpath(output, "fitted_unary_nodes.tsv"),
            transforms=joinpath(output, "fitted_edge_transforms.tsv"),
            receipt=joinpath(output, "receipt.toml"), output=output,
            report=report)
end

function make_producer_fixture(; dates=[string("2026010", index) for index in 1:7],
                                require_pass=true)
    root = mktempdir(; prefix="stmfit-structured-chain-correction6-")
    copy_authority(root)
    t10 = make_producer_todo10(root; dates=dates)
    t11 = make_producer_todo11(root, t10; require_pass=require_pass)
    return (root=root, t10=t10, t11=t11)
end

function _stage_model_result(result, outer_folds; final_gate=result.final_gate)
    states = String[fold.inner_gate.status for fold in outer_folds]
    append!(states, fold.outer_score.status for fold in outer_folds)
    push!(states, final_gate.status, result.full_fit.status)
    status = SEA.combine_terminal_states(states)
    reason = status == "PASS" ? :ok : status == "FAIL" ? :numerical_failure :
             status == "BLOCKED" ? :integrity_failure : :conditional_mechanism_not_admitted
    lines = [
        "schema=structured-edge-model-admission-v1",
        "model=$(result.model_id)",
        "status=$(status)",
        "reason=$(String(reason))",
        "final_gate=$(final_gate.gate_sha256)",
        "full_fit=$(result.full_fit.fit_sha256)",
    ]
    for fold in outer_folds
        push!(lines, join(("outer", fold.outer_date, fold.inner_gate.gate_sha256,
                           fold.outer_score.score_sha256, fold.fit_sha256), '\t'))
    end
    return SEA.ModelAdmissionResult(status, reason, result.model_id, outer_folds,
                                    final_gate, result.full_fit,
                                    SEA._hash_lines(lines))
end

function _stage_report_with_failure(report, stage::Symbol; failed_fit=nothing)
    original = report.models["C1"]
    fold = original.outer_folds[1]
    fit = fold.outer_score.fit
    failed_fit = failed_fit === nothing && stage == :null_only ?
        SEA._invalid_partition("FAIL", :conditional_fit_failed, fit.model_id,
                               copy(fit.training_dates), fit.reversed, fit.unary;
                               residualizer=fit.residualizer, samples=fit.samples,
                               transform_edges=fit.transform_edges, support=fit.support,
                               null_fit=fit.null_fit) : failed_fit
    failed_fit = failed_fit === nothing && stage == :conditional_failure ? begin
        conditional = fit.conditional_fit::SEA.ConditionalFit
        failed_starts = [SEA.ConditionalStartTrace(
            start.start_index, start.alpha, start.objective_trace, false,
            start.monotone, false, :conditional_fit_failed,
            start.rejected_iteration, start.rejected_objective,
            start.final_means, start.final_scale,
        ) for start in conditional.starts]
        failed_conditional = SEA.ConditionalFit(
            :fail, :conditional_fit_failed, conditional.means, conditional.scale,
            last(failed_starts[1].objective_trace), failed_starts, 1, false,
        )
        SEA._invalid_partition("FAIL", :conditional_fit_failed, fit.model_id,
                               copy(fit.training_dates), fit.reversed, fit.unary;
                               residualizer=fit.residualizer, samples=fit.samples,
                               transform_edges=fit.transform_edges, support=fit.support,
                               null_fit=fit.null_fit,
                               conditional_fit=failed_conditional)
    end : failed_fit
    failed_fit === nothing && begin
        error("unknown correction8 stage: $stage")
    end
    failed_evaluation = SEA.evaluate_partition(failed_fit, fold.outer_date;
                                               run_shuffle=false)
    folds = copy(original.outer_folds)
    folds[1] = SEA.OuterFoldEvidence(fold.outer_date, fold.inner_gate,
                                     failed_evaluation, failed_fit.fit_sha256)
    final_partitions = copy(original.final_gate.partitions)
    final_index = findfirst(partition -> partition.target_date == fold.outer_date,
                            final_partitions)
    final_index === nothing && error("stage final gate target is absent")
    final_partitions[final_index] = failed_evaluation
    final_gate = SEA.evaluate_gate(final_partitions)
    models = copy(report.models)
    models["C1"] = _stage_model_result(original, folds; final_gate=final_gate)
    statuses = [result.status for result in values(models)]
    status = SEA.combine_terminal_states(statuses)
    reason = status == "PASS" ? :ok : status == "FAIL" ? :numerical_failure :
             status == "BLOCKED" ? :integrity_failure : :conditional_mechanism_not_admitted
    normalization_hash = report.hashes["normalization_sha256"]
    lines = [
        "schema=structured-edge-admission-report-v1",
        "status=$(status)",
        "reason=$(String(reason))",
        "normalization=$(normalization_hash)",
    ]
    for model_id in sort!(collect(keys(models)))
        push!(lines, "model=$(model_id):$(models[model_id].result_sha256)")
    end
    for key in sort!(collect(keys(report.hashes)))
        key == "normalization_sha256" && continue
        push!(lines, "input=$(key):$(report.hashes[key])")
    end
    return SEA.AdmissionReport(status, reason, models, copy(report.hashes),
                               SEA._hash_lines(lines))
end

function make_producer_stage_fixture(kind::Symbol)
    base = make_producer_fixture()
    report = base.t11.report
    if kind == :zero_support
        original = report.models["C1"].outer_folds[1].outer_score.fit
        unavailable = [SEA.EdgeTransformState(state.edge, state.endpoint_predictors,
                                               "unavailable", :residualizer_unavailable)
                       for state in original.transform_edges]
        failed = SEA._invalid_partition("SKIPPED", :insufficient_dates,
                                        original.model_id, copy(original.training_dates),
                                        original.reversed, original.unary;
                                        transform_edges=unavailable,
                                        support=SEA._empty_support(:insufficient_dates))
        report = _stage_report_with_failure(report, :zero_support; failed_fit=failed)
    elseif kind == :insufficient_support
        original = report.models["C1"].outer_folds[1].outer_score.fit
        kept = Set((sample.edge.file, sample.edge.left_lobe, sample.edge.right_lobe)
                   for date in original.training_dates
                   for file in first(sort!(unique(sample.edge.file for sample in original.samples
                                                if sample.edge.date == date)), 4)
                   for sample in [first(filter(item -> item.edge.date == date &&
                                               item.edge.file == file, original.samples))])
        filtered_samples = [sample for sample in original.samples
                            if (sample.edge.file, sample.edge.left_lobe, sample.edge.right_lobe) in kept]
        filtered_transforms = [state for state in original.transform_edges
                               if (state.edge.file, state.edge.left_lobe, state.edge.right_lobe) in kept]
        unavailable = [state for state in original.transform_edges
                       if (state.edge.file, state.edge.left_lobe, state.edge.right_lobe) ∉ kept]
        support = SEA._support_evidence(
            [sample for sample in filtered_samples if sample.edge.date in original.training_dates],
            original.training_dates,
        )
        residualizer = original.residualizer::SEA.ResidualizerFit
        staged_residualizer = SEA.ResidualizerFit(
            residualizer.ridge, residualizer.coefficients, support.edges,
            residualizer.training_sha256, residualizer.evaluation_order,
        )
        failed = SEA._invalid_partition("SKIPPED", support.reason, original.model_id,
                                        copy(original.training_dates), original.reversed,
                                        original.unary; residualizer=staged_residualizer,
                                        samples=filtered_samples,
                                        transform_edges=vcat(filtered_transforms, [
                                            SEA.EdgeTransformState(state.edge,
                                                state.endpoint_predictors,
                                                "unavailable", :residualizer_unavailable)
                                            for state in unavailable]),
                                        support=support)
        report = _stage_report_with_failure(report, :insufficient_support; failed_fit=failed)
    elseif kind in (:null_only, :conditional_failure)
        report = _stage_report_with_failure(report, kind)
    else
        error("unknown correction8 stage fixture: $kind")
    end
    output = base.t11.output
    # Replace only the serializer-produced report bytes in this private fixture.
    for (name, bytes) in SEA.report_files(report)
        write(joinpath(output, name), bytes)
    end
    return (root=base.root, t10=base.t10, t11=base.t11, report=report, kind=kind)
end

function copy_producer_fixture(source)
    root = mktempdir(; prefix="stmfit-structured-chain-correction7-case-")
    copy_authority(root)
    cp(joinpath(source.root, "todo10"), joinpath(root, "todo10"); force=true)
    cp(joinpath(source.root, "todo11"), joinpath(root, "todo11"); force=true)
    t10 = (edge=joinpath(root, "todo10", "edge_observations.tsv"),
           node=joinpath(root, "todo10", "node_segments.tsv"),
           receipt=joinpath(root, "todo10", "receipt.toml"),
           dates=source.t10.dates, files=source.t10.files,
           nodes=source.t10.nodes, edges=source.t10.edges)
    output = joinpath(root, "todo11")
    t11 = (models=joinpath(output, "fitted_edge_models.tsv"),
           unary=joinpath(output, "fitted_unary_nodes.tsv"),
           transforms=joinpath(output, "fitted_edge_transforms.tsv"),
           receipt=joinpath(output, "receipt.toml"), output=output)
    return (root=root, t10=t10, t11=t11)
end

function refresh_t11_wrappers!(fixture; graph_lists=false)
    output = fixture.t11.output
    receipt_path = fixture.t11.receipt
    receipt = TOML.parse(read(receipt_path, String))
    artifact_names = vcat(SCI._REPORT_FILES, SCI._GRAPH_HANDOFF_FILES)
    bytes = Dict(name => read(joinpath(output, name)) for name in artifact_names)
    for name in artifact_names
        receipt["artifacts"][name] = sha(bytes[name])
        receipt["row_counts"][name] = length(split(chomp(String(copy(bytes[name]))), '\n')) - 1
    end
    receipt["graph_handoff"]["model_row_count"] = receipt["row_counts"]["fitted_edge_models.tsv"]
    receipt["graph_handoff"]["unary_node_row_count"] = receipt["row_counts"]["fitted_unary_nodes.tsv"]
    receipt["graph_handoff"]["edge_transform_row_count"] = receipt["row_counts"]["fitted_edge_transforms.tsv"]
    if graph_lists
        model_lines = split(read(joinpath(output, "fitted_edge_models.tsv"), String), '\n')
        model_rows = [split(line, '\t'; keepempty=true) for line in model_lines[2:end-1]]
        fit_col = findfirst(==("fit_sha256"), SCI._GRAPH_MODEL_HEADER)
        role_col = findfirst(==("fit_role"), SCI._GRAPH_MODEL_HEADER)
        receipt["graph_handoff"]["fit_hashes"] = sort(unique(row[fit_col] for row in model_rows))
        receipt["graph_handoff"]["full_refit_hashes"] = sort(unique(
            row[fit_col] for row in model_rows if row[role_col] == "full_refit"))
    end
    receipt["graph_handoff_replay_sha256"] = graph_replay_hash(
        Dict(name => bytes[name] for name in SCI._GRAPH_HANDOFF_FILES))
    write(receipt_path, sprint(TOML.print, receipt))
    return nothing
end

function tsv_rows(path)
    lines = split(read(path, String), '\n'; keepempty=true)
    return lines[1], [split(line, '\t'; keepempty=true) for line in lines[2:end-1]]
end

function write_tsv_rows(path, header, rows)
    write(path, header * "\n" * join([join(row, '\t') for row in rows], "\n") * "\n")
end

function mutate_tsv_row!(path, row_index, mutator)
    header, rows = tsv_rows(path)
    mutator(rows[row_index])
    write_tsv_rows(path, header, rows)
end

mutate_tsv_row!(mutator::Function, path, row_index) =
    mutate_tsv_row!(path, row_index, mutator)

function mutate_tsv_rows!(path, mutator)
    header, rows = tsv_rows(path)
    mutator(rows)
    write_tsv_rows(path, header, rows)
end

mutate_tsv_rows!(mutator::Function, path) = mutate_tsv_rows!(path, mutator)

function refresh_todo10_and_t11!(fixture)
    refresh_t10_bindings!(fixture)
    refresh_t11_input_bindings!(fixture)
    refresh_t11_wrappers!(fixture)
end

function matrix_result_row(name, family, fixture, mutation, refreshed, layer,
                           expected_status, expected_reason, report)
    return [name, family, fixture, mutation, refreshed, layer,
            string(expected_status), string(expected_reason), string(report.status),
            string(report.reason), string(SCI._SCAN_DECODE_CALL_COUNT[]),
            string(SCI._FACTOR_BUILD_CALL_COUNT[]), string(SCI._DP_CALL_COUNT[])]
end

function rename_model_identity!(fixture, old::String, new::String)
    for name in vcat(SCI._REPORT_FILES, SCI._GRAPH_HANDOFF_FILES)
        path = joinpath(fixture.t11.output, name)
        header, rows = tsv_rows(path)
        model_column = findfirst(==( "model"), split(header, '\t'))
        model_column === nothing && continue
        for row in rows
            row[model_column] == old && (row[model_column] = new)
            name == "fitted_edge_models.tsv" &&
                (row[16] = replace(row[16], "|$(old)|" => "|$(new)|"))
        end
        write_tsv_rows(path, header, rows)
    end
    receipt = TOML.parse(read(fixture.t11.receipt, String))
    receipt["models"][new] = receipt["models"][old]
    delete!(receipt["models"], old)
    write(fixture.t11.receipt, sprint(TOML.print, receipt))
    refresh_t11_wrappers!(fixture; graph_lists=true)
end

function split_fit_alias!(fixture, old_ref::String, split_ref::String)
    model_path = fixture.t11.models
    header, rows = tsv_rows(model_path)
    positions = Dict(name => findfirst(==(name), SCI._GRAPH_MODEL_HEADER)
                     for name in SCI._GRAPH_MODEL_HEADER)
    source_index = only(findall(row -> occursin(split_ref, row[positions["references"]]), rows))
    source = rows[source_index]
    old_fit = source[positions["fit_sha256"]]
    references = split(source[positions["references"]], ';')
    deleteat!(references, findfirst(==(split_ref), references))
    source[positions["reference_count"]] = string(length(references))
    source[positions["references"]] = join(references, ';')
    source[positions["references_sha256"]] = hash_lines([join(split(ref, '|'), '\t') for ref in references])
    new_fit = repeat("9", 64)
    clone = copy(source)
    clone[positions["fit_sha256"]] = new_fit
    clone[positions["fit_role"]] = "partition"
    fields = split(split_ref, '|')
    clone[positions["scope"]] = fields[3]
    clone[positions["outer_date"]] = fields[4]
    clone[positions["target_date"]] = fields[5]
    clone[positions["reference_count"]] = "1"
    clone[positions["references"]] = split_ref
    clone[positions["references_sha256"]] = hash_lines([join(split(split_ref, '|'), '\t')])
    push!(rows, clone)
    sort!(rows; by=row -> row[positions["fit_sha256"]])
    write_tsv_rows(model_path, header, rows)

    for name in ("fitted_unary_nodes.tsv", "fitted_edge_transforms.tsv")
        path = joinpath(fixture.t11.output, name)
        h, values = tsv_rows(path)
        for row in values
            row[2] == old_fit && (row[2] = new_fit)
        end
        write_tsv_rows(path, h, values)
    end
    for (name, fit_column) in (("partitions.tsv", 12), ("ess.tsv", 16), ("starts.tsv", 17))
        path = joinpath(fixture.t11.output, name)
        h, values = tsv_rows(path)
        for row in values
            row[fit_column] == old_fit && (row[fit_column] = new_fit)
        end
        write_tsv_rows(path, h, values)
    end
    # The graph validator intentionally reaches the alias check before the
    # diagnostic rows; retargeting the canonical fit rows is sufficient to
    # keep the mutation outside shallow artifact-hash failure.
    refresh_t11_wrappers!(fixture; graph_lists=true)
end

function infer(fixture)
    return SCI.infer_structured_chains(fixture.root;
        edge_observations=fixture.t10.edge,
        node_segments=fixture.t10.node,
        todo10_receipt=fixture.t10.receipt,
        fitted_edge_models=fixture.t11.models,
        fitted_unary_nodes=fixture.t11.unary,
        fitted_edge_transforms=fixture.t11.transforms,
        todo11_receipt=fixture.t11.receipt)
end

function exhaustive(unaries, factors)
    n = length(unaries)
    logweights = Float64[]
    labels = Vector{Vector{Int}}()
    for code in 0:(2^n - 1)
        sequence = [(code >> (i - 1)) & 1 for i in 1:n]
        value = sum(unaries[i][sequence[i] + 1] for i in 1:n)
        for i in 1:(n - 1)
            value += factors[i][sequence[i] * 2 + sequence[i + 1] + 1]
        end
        push!(logweights, value)
        push!(labels, sequence)
    end
    m = maximum(logweights)
    logz = m + log(sum(exp(value - m) for value in logweights))
    marginal = zeros(n, 2)
    for (weight, sequence) in zip(logweights, labels)
        for i in 1:n
            marginal[i, sequence[i] + 1] += exp(weight - logz)
        end
    end
    best = argmax(logweights)
    return logz, marginal, labels[best], logweights[best]
end

function strict_probability_invariants(report)
    nodes = [node for reference in report.references for block in reference.blocks
             for node in block.nodes]
    @test !isempty(nodes)
    for node in nodes
        p0, p1 = node.marginals
        l0, l1 = node.log_marginals
        @test isfinite(p0) && isfinite(p1)
        @test 0.0 <= p0 <= 1.0
        @test 0.0 <= p1 <= 1.0
        @test p0 + p1 == 1.0
        @test !(p0 > p1) || l0 > l1
        @test !(p1 > p0) || l1 > l0
        l0 == l1 && @test (p0, p1) == (0.5, 0.5)
    end
    return nothing
end

function hard_output_projection(report)
    return (
        report.status,
        report.reason,
        report.log_evidence,
        report.viterbi_score,
        [(
            reference.fit_sha256,
            reference.model,
            reference.fit_role,
            reference.scope,
            reference.outer_date,
            reference.target_date,
            reference.status,
            reference.reason,
            [(
                block.segment_id,
                block.log_evidence,
                block.viterbi_score,
                block.labels,
                [(
                    factor.file,
                    factor.date,
                    factor.left_lobe,
                    factor.right_lobe,
                    factor.segment_id,
                    factor.factor,
                    factor.factor_sha256,
                    factor.raw_row_sha256,
                ) for factor in block.factors],
                [(
                    node.file,
                    node.date,
                    node.lobe,
                    node.t_nm,
                    node.log_marginals,
                    node.label,
                ) for node in block.nodes],
            ) for block in reference.blocks],
        ) for reference in report.references],
    )
end

function exact_bigfloat_binary_softmax(l0::Float64, l1::Float64, precision::Int)
    return setprecision(BigFloat, precision) do
        b0, b1 = BigFloat(l0), BigFloat(l1)
        shift = max(b0, b1)
        weight0 = exp(b0 - shift)
        weight1 = exp(b1 - shift)
        total = weight0 + weight1
        (weight0 / total, weight1 / total)
    end
end

function nonnegative_ulp_distance(left::Float64, right::Float64)
    @assert isfinite(left) && isfinite(right) && left >= 0.0 && right >= 0.0
    return abs(Int128(reinterpret(UInt64, left)) - Int128(reinterpret(UInt64, right)))
end

function assert_bigfloat_softmax_oracle(l0::Float64, l1::Float64)
    exact256 = exact_bigfloat_binary_softmax(l0, l1, 256)
    exact512 = exact_bigfloat_binary_softmax(l0, l1, 512)
    rounded256 = Float64.(exact256)
    rounded512 = Float64.(exact512)
    @test rounded256 == rounded512

    actual = SCI._binary_softmax(l0, l1)
    relative_limit = BigFloat(4) * BigFloat(eps(Float64))
    denominator_floor = BigFloat(floatmin(Float64))
    for index in 1:2
        @test nonnegative_ulp_distance(actual[index], rounded512[index]) <= 2
        relative_error = abs(BigFloat(actual[index]) - exact512[index]) /
                         max(abs(exact512[index]), denominator_floor)
        @test relative_error <= relative_limit
    end
    @test actual[1] + actual[2] == 1.0
    return actual
end

function float_neighborhood(center::Float64, radius::Int=4)
    values = Float64[center]
    lower = center
    upper = center
    for _ in 1:radius
        lower = prevfloat(lower)
        upper = nextfloat(upper)
        push!(values, lower, upper)
    end
    return sort!(unique!(values))
end

@testset "T12 deterministic binary softmax normalization" begin
    # RED witness: direct exponentiation can produce an invalid probability.
    red_pair = (1.7763568394002505e-15, -38.36886073604529)
    old_pair = exp.(red_pair)
    @test old_pair[1] > 1.0
    @test !(all(value -> 0.0 <= value <= 1.0, old_pair) &&
            old_pair[1] + old_pair[2] == 1.0)
    red_actual = assert_bigfloat_softmax_oracle(red_pair...)
    @test red_actual[1] == 1.0
    @test 0.0 < red_actual[2] < 1.0

    @test SCI._binary_softmax(0.0, 0.0) == (0.5, 0.5)
    @test SCI._binary_softmax(0.0, -Inf) == (1.0, 0.0)
    @test SCI._binary_softmax(-Inf, 0.0) == (0.0, 1.0)
    @test SCI._binary_softmax(0.0, -1000.0) == (1.0, 0.0)
    @test SCI._binary_softmax(1000.0, 0.0) == (1.0, 0.0)
    @test SCI._binary_softmax(-1000.0, 0.0) == (0.0, 1.0)

    half_min_subnormal_log = setprecision(BigFloat, 256) do
        Float64(log(BigFloat(nextfloat(0.0)) / BigFloat(2)))
    end
    boundary_centers = (
        0.0,
        log(eps(Float64) / 2),
        log(floatmin(Float64)),
        log(nextfloat(0.0)),
        half_min_subnormal_log,
    )
    for center in boundary_centers, boundary in float_neighborhood(center, 4)
        assert_bigfloat_softmax_oracle(0.0, boundary)
        assert_bigfloat_softmax_oracle(boundary, 0.0)
    end

    one_quarter = SCI._binary_softmax(log(0.25), log(0.75))
    three_quarters = SCI._binary_softmax(log(0.75), log(0.25))
    @test one_quarter == (0.25, 0.75)
    @test three_quarters == (0.75, 0.25)
    @test SCI._binary_softmax(0.25, -0.75) ==
          reverse(SCI._binary_softmax(-0.75, 0.25))

    offset = 1.0e12
    @test SCI._binary_softmax(offset + 0.25, offset) ==
          SCI._binary_softmax(0.25, 0.0)
    huge_offset = 1.0e300
    @test SCI._binary_softmax(huge_offset + 1.0e285, huge_offset) == (1.0, 0.0)
    @test SCI._binary_softmax(floatmax(Float64), -floatmax(Float64)) == (1.0, 0.0)
    @test SCI._binary_softmax(-floatmax(Float64), floatmax(Float64)) == (0.0, 1.0)

    for pair in ((nextfloat(0.0), 0.0), (0.0, nextfloat(0.0)),
                 (prevfloat(0.0), 0.0), (0.0, prevfloat(0.0)))
        first = SCI._binary_softmax(pair...)
        second = SCI._binary_softmax(pair...)
        @test first == second
        @test first[1] + first[2] == 1.0
    end

    for (l0, l1) in ((-Inf, -Inf), (NaN, 0.0), (0.0, NaN),
                     (Inf, 0.0), (0.0, Inf), (-Inf, Inf), (Inf, -Inf))
        error = try
            SCI._binary_softmax(l0, l1)
        catch caught
            caught
        end
        @test error isa SCI._InferenceFailed
        @test error.code == :numerical_failure
    end

    rng = MersenneTwister(0x5354464954)
    for _ in 1:10_000
        l0, l1 = 24.0 * randn(rng), 24.0 * randn(rng)
        actual = assert_bigfloat_softmax_oracle(l0, l1)
        swapped = SCI._binary_softmax(l1, l0)
        @test actual == reverse(swapped)
        @test !(actual[1] > actual[2]) || l0 > l1
        @test !(actual[2] > actual[1]) || l1 > l0

        delta_hi, delta_lo = SCI._two_diff(l0, l1)
        @test isfinite(delta_hi)
        @test BigFloat(delta_hi) + BigFloat(delta_lo) == BigFloat(l0) - BigFloat(l1)
    end
end

function infer_with_paths(fixture; root=fixture.root, edge=fixture.t10.edge,
                          node=fixture.t10.node, t10_receipt=fixture.t10.receipt,
                          models=fixture.t11.models, unary=fixture.t11.unary,
                          transforms=fixture.t11.transforms, receipt=fixture.t11.receipt)
    return SCI.infer_structured_chains(root; edge_observations=edge, node_segments=node,
        todo10_receipt=t10_receipt, fitted_edge_models=models,
        fitted_unary_nodes=unary, fitted_edge_transforms=transforms,
        todo11_receipt=receipt)
end

function rebind_todo10_source!(fixture, source; t11_source=source)
    mutate_tsv_rows!(fixture.t10.node) do rows
        for row in rows
            row[10] = source
        end
    end
    mutate_tsv_rows!(fixture.t10.edge) do rows
        for row in rows
            row[18] = source
        end
    end
    todo10_receipt = TOML.parse(read(fixture.t10.receipt, String))
    todo10_receipt["source_sha256"] = source
    write(fixture.t10.receipt, sprint(TOML.print, todo10_receipt))
    refresh_t10_bindings!(fixture)

    todo11_receipt = TOML.parse(read(fixture.t11.receipt, String))
    todo11_receipt["hashes"]["edge_source_sha256"] = t11_source
    write(fixture.t11.receipt, sprint(TOML.print, todo11_receipt))
    refresh_t11_input_bindings!(fixture)
    return nothing
end

function set_todo10_receipt_source!(fixture, source)
    receipt = TOML.parse(read(fixture.t10.receipt, String))
    receipt["source_sha256"] = source
    write(fixture.t10.receipt, sprint(TOML.print, receipt))
    return nothing
end

function set_todo11_edge_source!(fixture, source)
    receipt = TOML.parse(read(fixture.t11.receipt, String))
    receipt["hashes"]["edge_source_sha256"] = source
    write(fixture.t11.receipt, sprint(TOML.print, receipt))
    return nothing
end

function validation_registry(fixture)
    authority = SCI._authority_snapshots(fixture.root)
    supplied = [SCI._snapshot(fixture.root, path, context) for (path, context) in [
        (fixture.t10.edge, "edge"), (fixture.t10.node, "node"),
        (fixture.t10.receipt, "receipt"), (fixture.t11.models, "models"),
        (fixture.t11.unary, "unary"), (fixture.t11.transforms, "transforms"),
        (fixture.t11.receipt, "admission"),
    ]]
    return Dict(snapshot.relative => snapshot for snapshot in vcat(authority, supplied))
end

function todo10_validation_error(fixture)
    registry = validation_registry(fixture)
    try
        SCI._validate_todo10(fixture.root, fixture.t10.edge, fixture.t10.node,
                             fixture.t10.receipt, registry)
    catch error
        error isa SCI._ChainBlocked || rethrow()
        return error
    end
    return nothing
end

function todo11_validation_error(fixture)
    registry = validation_registry(fixture)
    todo10 = SCI._validate_todo10(fixture.root, fixture.t10.edge, fixture.t10.node,
                                  fixture.t10.receipt, registry)
    try
        SCI._validate_todo11(fixture.root, fixture.t11.models, fixture.t11.unary,
                             fixture.t11.transforms, fixture.t11.receipt, registry,
                             todo10)
    catch error
        error isa SCI._ChainBlocked || rethrow()
        return error
    end
    return nothing
end

function graph_replay_hash(files::Dict{String,Vector{UInt8}})
    io = IOBuffer()
    write(io, codeunits("schema=structured-edge-graph-handoff-replay-v1\n"))
    for name in SCI._GRAPH_HANDOFF_FILES
        bytes = files[name]
        write(io, codeunits("file=$(name)\nlength=$(length(bytes))\n"))
        write(io, bytes)
        write(io, UInt8('\n'))
    end
    return sha(take!(io))
end

function fixture_graph_replay(fixture)
    receipt = TOML.parse(read(fixture.t11.receipt, String))
    bytes = Dict(name => read(joinpath(fixture.t11.output, name))
                 for name in SCI._GRAPH_HANDOFF_FILES)
    return receipt["graph_handoff_replay_sha256"], graph_replay_hash(bytes)
end

function refresh_t11_bindings!(fixture; include_handoff=true)
    output = fixture.t11.output
    receipt_path = fixture.t11.receipt
    receipt = TOML.parse(read(receipt_path, String))
    artifact_names = vcat(SCI._REPORT_FILES, SCI._GRAPH_HANDOFF_FILES)
    bytes = Dict(name => read(joinpath(output, name)) for name in artifact_names)
    for name in artifact_names
        receipt["artifacts"][name] = sha(bytes[name])
        receipt["row_counts"][name] = length(split(chomp(String(copy(bytes[name]))), '\n')) - 1
    end
    receipt["graph_handoff"]["model_row_count"] = receipt["row_counts"]["fitted_edge_models.tsv"]
    receipt["graph_handoff"]["unary_node_row_count"] = receipt["row_counts"]["fitted_unary_nodes.tsv"]
    receipt["graph_handoff"]["edge_transform_row_count"] = receipt["row_counts"]["fitted_edge_transforms.tsv"]
    include_handoff && (receipt["graph_handoff_replay_sha256"] = graph_replay_hash(Dict(name => bytes[name] for name in SCI._GRAPH_HANDOFF_FILES)))
    write(receipt_path, sprint(TOML.print, receipt))
    return nothing
end

function refresh_t11_input_bindings!(fixture)
    receipt = TOML.parse(read(fixture.t11.receipt, String))
    receipt["hashes"]["edge_receipt_sha256"] = sha(read(fixture.t10.receipt))
    receipt["hashes"]["edge_observations_sha256"] = sha(read(fixture.t10.edge))
    receipt["hashes"]["node_segments_sha256"] = sha(read(fixture.t10.node))
    result_lines = ["schema=structured-edge-admission-report-v1",
        "status=$(receipt["status"])", "reason=$(receipt["reason"])",
        "normalization=$(receipt["hashes"]["normalization_sha256"])"]
    for model_id in sort!(collect(keys(receipt["models"])))
        binding = split(String(receipt["models"][model_id]), '|'; keepempty=true)
        push!(result_lines, "model=$(model_id):$(binding[3])")
    end
    for key in sort!(collect(keys(receipt["hashes"])))
        key == "normalization_sha256" && continue
        push!(result_lines, "input=$(key):$(receipt["hashes"][key])")
    end
    receipt["result_sha256"] = hash_lines(result_lines)
    write(fixture.t11.receipt, sprint(TOML.print, receipt))
end

function refresh_t10_bindings!(fixture; tables=true)
    receipt = TOML.parse(read(fixture.t10.receipt, String))
    edge_bytes = read(fixture.t10.edge)
    node_bytes = read(fixture.t10.node)
    tables && begin
        receipt["edge_observations_sha256"] = sha(edge_bytes)
        receipt["node_segments_sha256"] = sha(node_bytes)
        receipt["edge_row_count"] = length(split(chomp(String(copy(edge_bytes))), '\n')) - 1
        receipt["node_row_count"] = length(split(chomp(String(copy(node_bytes))), '\n')) - 1
    end
    write(fixture.t10.receipt, sprint(TOML.print, receipt))
    return nothing
end

function assert_blocked_zero(fixture, report)
    @test report.status == :BLOCKED
    @test SCI._DP_CALL_COUNT[] == 0
    @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
    @test SCI._SCAN_DECODE_CALL_COUNT[] == 0
end

function make_split_fixture()
    fixture = make_fixture(dates=("20260101",), shared_references=false)
    file = fixture.t10.file
    old_segment = segment_id(file, 1, 1)
    new_segment = segment_id(file, 2, 3)
    node_lines = split(read(fixture.t10.node, String), '\n'; keepempty=true)
    node_rows = [split(node_lines[index], '\t'; keepempty=true) for index in 2:4]
    node_rows[1][3] = fmt(1.0); node_rows[1][4] = old_segment; node_rows[1][5] = "isolated"; node_rows[1][7] = "gap_out_of_range"
    node_rows[2][3] = fmt(2.0); node_rows[2][4] = new_segment; node_rows[2][6] = "gap_out_of_range"
    node_rows[3][3] = fmt(2.4); node_rows[3][4] = new_segment
    write(fixture.t10.node, join(vcat([node_lines[1]], [join(row, '\t') for row in node_rows], [""]), '\n'))
    edge_lines = split(read(fixture.t10.edge, String), '\n'; keepempty=true)
    edge_rows = [split(edge_lines[index], '\t'; keepempty=true) for index in 2:3]
    edge_rows[1][5] = fmt(2.0); edge_rows[1][6] = fmt(1.0); edge_rows[1][7] = old_segment; edge_rows[1][9] = "split"; edge_rows[1][10] = "gap_out_of_range"
    edge_rows[1][8] = new_segment; edge_rows[1][11] = "NA"; edge_rows[1][12] = "NA"
    edge_rows[2][4] = fmt(2.0); edge_rows[2][5] = fmt(2.4); edge_rows[2][6] = fmt(2.4 - 2.0)
    edge_rows[2][7] = new_segment; edge_rows[2][8] = new_segment
    write(fixture.t10.edge, join(vcat([edge_lines[1]], [join(row, '\t') for row in edge_rows], [""]), '\n'))
    refresh_t10_bindings!(fixture)
    refresh_t11_input_bindings!(fixture)
    transform_lines = split(read(fixture.t11.transforms, String), '\n'; keepempty=true)
    transform_row = split(transform_lines[3], '\t'; keepempty=true)
    transform_row[11] = fmt(2.0); transform_row[12] = fmt(2.4); transform_row[13] = new_segment
    transform_row[15] = join((file, "2", "3", new_segment), '|')
    transformed_edge_line = join(edge_rows[2], '\t') * "\n"
    transform_row[16] = sha(codeunits(transformed_edge_line))
    write(fixture.t11.transforms, join(vcat([transform_lines[1]], [join(transform_row, '\t')], [""]), '\n'))
    ess_header = join(SCI._ESS_HEADER, '\t')
    ess_rows = [["1", "C1", "full_refit", "NA", "NA", string(category), "1", "1", "1", "1", "1", "1", "1", "true", "ok", fixture.t11.fit] for category in 0:3]
    write(joinpath(fixture.t11.output, "ess.tsv"), ess_header * "\n" * join([join(row, '\t') for row in ess_rows], "\n") * "\n")
    residualizer = hash_lines(vcat([
        "schema=structured-edge-residualizer-v1", "training_sha256=$(repeat("e", 64))",
        "training_row_count=1", "ridge=0", "evaluation_order=$(join(collect(1:15), ','))",
    ], ["$(name)_$(output)=0" for name in SCI._GRAPH_COEFFICIENT_NAMES for output in SCI._GRAPH_OUTPUT_NAMES]))
    rewrite_tsv_row!(fixture.t11.models, 1, row -> (row[19] = residualizer))
    transform_lines = split(read(fixture.t11.transforms, String), '\n'; keepempty=true)
    for index in 2:(length(transform_lines) - 1)
        fields = split(transform_lines[index], '\t'; keepempty=true)
        fields[17] = residualizer
        transform_lines[index] = join(fields, '\t')
    end
    write(fixture.t11.transforms, join(transform_lines, '\n'))
    refresh_t11_bindings!(fixture)
    return fixture
end

function rewrite_tsv_row!(path, row_number, mutator)
    lines = split(read(path, String), '\n'; keepempty=true)
    row = split(lines[row_number + 1], '\t'; keepempty=true)
    mutator(row)
    lines[row_number + 1] = join(row, '\t')
    write(path, join(lines, '\n'))
end

function rebind_unary_probability!(fixture)
    lines = split(read(fixture.t11.unary, String), '\n'; keepempty=true)
    rows = [split(lines[index], '\t'; keepempty=true) for index in 2:(length(lines) - 1)]
    by_fit = Dict{String,Vector{Vector{String}}}()
    for row in rows
        push!(get!(by_fit, row[2], Vector{Vector{String}}()), row)
    end
    probabilities = Dict(fit => matrix_hash([parse(Float64, row[11]) for row in fit_rows])
                         for (fit, fit_rows) in by_fit)
    for row in rows
        row[13] = probabilities[row[2]]
    end
    write(fixture.t11.unary, join(vcat([lines[1]], [join(row, '\t') for row in rows], [""]), '\n'))
    model_lines = split(read(fixture.t11.models, String), '\n'; keepempty=true)
    for index in 2:(length(model_lines) - 1)
        row = split(model_lines[index], '\t'; keepempty=true)
        row[17] = probabilities[row[3]]
        model_lines[index] = join(row, '\t')
    end
    write(fixture.t11.models, join(model_lines, '\n'))
    refresh_t11_bindings!(fixture)
end

@testset "public path-only graph inference" begin
    fixture = make_fixture()
    try
        report = infer(fixture)
        @test report.status == :PASS
        @test report.reason == :graph_inference
        @test length(report.references) == 21
        @test length(report.references[1].blocks) == 1
        @test length(report.references[1].blocks[1].nodes) == 3
        @test length(report.references[1].blocks[1].factors) == 2
        @test [node.lobe for node in report.references[1].blocks[1].nodes] == [1, 2, 3]
        @test [factor.left_lobe for factor in report.references[1].blocks[1].factors] == [1, 2]
        @test all(factor.raw_row_sha256 == repeat("", 0) ||
                  occursin(r"^[0-9a-f]{64}$", factor.raw_row_sha256)
                  for factor in report.references[1].blocks[1].factors)
        @test length(report.references[1].blocks[1].labels) == 3
        @test all(isapprox(sum(node.marginals), 1.0; atol=1e-12, rtol=0)
                  for node in report.references[1].blocks[1].nodes)
        strict_probability_invariants(report)
        report_again = infer(fixture)
        @test report_again.result_sha256 == report.result_sha256
        @test report_again.provenance_sha256 == report.provenance_sha256
        @test hard_output_projection(report_again) == hard_output_projection(report)
        @test report_again.log_evidence == report.log_evidence
        @test report_again.viterbi_score == report.viterbi_score
        authority = SCI._authority_snapshots(fixture.root)
        supplied = [SCI._snapshot(fixture.root, path, context) for (path, context) in [
            (fixture.t10.edge, "edge"), (fixture.t10.node, "node"),
            (fixture.t10.receipt, "receipt"), (fixture.t11.models, "models"),
            (fixture.t11.unary, "unary"), (fixture.t11.transforms, "transforms"),
            (fixture.t11.receipt, "admission"),
        ]]
        registry = Dict(snapshot.relative => snapshot for snapshot in vcat(authority, supplied))
        validated_t10 = SCI._validate_todo10(fixture.root, fixture.t10.edge,
                                              fixture.t10.node, fixture.t10.receipt, registry)
        validated_t11 = SCI._validate_todo11(fixture.root, fixture.t11.models,
                                              fixture.t11.unary, fixture.t11.transforms,
                                              fixture.t11.receipt, registry,
                                              validated_t10)
        model = first(validated_t11.models)
        edge = first(edge for edge in validated_t10.edges if edge.date == report.references[1].target_date)
        transform = validated_t11.transforms[(model.fit_sha256, edge.file, edge.date,
                                              edge.left_lobe, edge.right_lobe)]
        @test SCI._FACTOR_BUILD_CALL_COUNT[] == 40
        @test SCI._SCAN_DECODE_CALL_COUNT[] == 20
        residual = [edge.corr_fwd - transform.pred_fwd,
                    edge.corr_bwd - transform.pred_bwd]
        null_density = SCI._student_log_density(residual, model.null_mean, model.null_scale)
        hand = ntuple(index -> SCI._student_log_density(
            residual, vec(model.conditional_means[index, :]), model.conditional_scale) -
            null_density, 4)
        @test SCI._factor(model, edge, transform) == hand
        @test transform.raw_sha256 == edge.raw_sha256
        factor_result = report.references[1].blocks[1].factors[1]
        @test factor_result.factor == hand
        @test factor_result.factor_sha256 == SCI._factor_sha(edge, model, hand)
        @test factor_result.raw_row_sha256 == edge.raw_sha256
        @test factor_result.fit_sha256 == model.fit_sha256
        @test SCI._DP_CALL_COUNT[] == 48
        println("structured_chain_inference_result_sha256=$(report.result_sha256)")
        println("structured_chain_inference_provenance_sha256=$(report.provenance_sha256)")
        @test length(methods(SCI.infer_structured_chains)) == 1
        @test fieldtype(SCI.ChainInferenceReport, :status) == Symbol
    finally
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "return-time snapshot mutation is blocked" begin
    fixture = make_fixture()
    try
        original = read(fixture.t10.edge)
        SCI._SNAPSHOT_RETURN_HOOK[] = () -> write(fixture.t10.edge, vcat(original, UInt8('#')))
        report = infer(fixture)
        @test report.status == :BLOCKED
        @test report.reason == :input_changed
        @test SCI._DP_CALL_COUNT[] > 0
    finally
        SCI._SNAPSHOT_RETURN_HOOK[] = nothing
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "every fixed authority mutation blocks before inference" begin
    fixture = make_fixture()
    authority = SCI._authority_snapshots(fixture.root)
    try
        for authority_snapshot in authority
            isolated = make_fixture()
            try
                isolated_snapshot = only(filter(snapshot -> snapshot.relative == authority_snapshot.relative,
                                                SCI._authority_snapshots(isolated.root)))
                write(isolated_snapshot.path, vcat(read(isolated_snapshot.path), UInt8('#')))
                assert_blocked_zero(isolated, infer(isolated))
            finally
                rm(isolated.root; recursive=true, force=true)
            end
        end
        expected_paths = Set([
            ".omo/plans/structured-label-free-unit-assignment.md",
            "config/unit_assignment_structured_candidate.toml",
            "config/unit_assignment_structured_model.toml",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase3-t3-universe/correction/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase3-t3-universe/correction/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/correction/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/correction/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/claim-correction/ClaimCorrection.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/canonical-t11-source-bundle.bytes",
            "test/lib/structured_assignment/edge_admission.jl",
            "test/evaluate_structured_edge_admission.jl",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/canonical-t11-source-bundle.bytes",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/canonical-t11-source-bundle.bytes",
        ])
        @test length(authority) == 25
        @test Set(snapshot.relative for snapshot in authority) == expected_paths
        @test length(Set(snapshot.relative for snapshot in authority)) == length(authority)
    finally
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "dynamic per-report graph replay boundary" begin
    pass_fixture = make_fixture()
    skipped_fixture = make_fixture(
        dates=("20260101", "20260102", "20260103", "20260104"),
        terminal_status="SKIPPED",
    )
    try
        pass_report = infer(pass_fixture)
        skipped_report = infer(skipped_fixture)
        pass_receipt_replay, pass_recomputed = fixture_graph_replay(pass_fixture)
        skipped_receipt_replay, skipped_recomputed = fixture_graph_replay(skipped_fixture)
        @test pass_report.status == :PASS
        @test skipped_report.status == :SKIPPED
        strict_probability_invariants(pass_report)
        strict_probability_invariants(skipped_report)
        @test pass_receipt_replay == pass_recomputed
        @test skipped_receipt_replay == skipped_recomputed
        @test pass_receipt_replay != skipped_receipt_replay
        @test pass_receipt_replay !=
              "dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474"
        @test skipped_receipt_replay !=
              "dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474"
    finally
        rm(pass_fixture.root; recursive=true, force=true)
        rm(skipped_fixture.root; recursive=true, force=true)
    end

    literal_fixture = make_fixture()
    try
        receipt = TOML.parse(read(literal_fixture.t11.receipt, String))
        original_replay = receipt["graph_handoff_replay_sha256"]
        @test original_replay !=
              "dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474"
        receipt["graph_handoff_replay_sha256"] =
            "dc15a94bb8bf373a62dcb1fba1598711ece0f99882082de7568eac6aab6aa474"
        write(literal_fixture.t11.receipt, sprint(TOML.print, receipt))
        SCI._DP_CALL_COUNT[] = 0
        blocked = infer(literal_fixture)
        @test blocked.status == :BLOCKED
        @test SCI._SCAN_DECODE_CALL_COUNT[] == 0
        @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
        @test SCI._DP_CALL_COUNT[] == 0
    finally
        rm(literal_fixture.root; recursive=true, force=true)
    end

    stale_fixture = make_fixture()
    try
        write(stale_fixture.t11.models, vcat(read(stale_fixture.t11.models), UInt8('#')))
        SCI._DP_CALL_COUNT[] = 0
        blocked = infer(stale_fixture)
        @test blocked.status == :BLOCKED
        @test SCI._SCAN_DECODE_CALL_COUNT[] == 0
        @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
        @test SCI._DP_CALL_COUNT[] == 0
    finally
        rm(stale_fixture.root; recursive=true, force=true)
    end

    semantic_fixture = make_fixture()
    try
        mutate_tsv_row!(semantic_fixture.t11.unary, 1) do row
            row[9] = "1,0"
        end
        refresh_t11_wrappers!(semantic_fixture)
        SCI._DP_CALL_COUNT[] = 0
        blocked = infer(semantic_fixture)
        @test blocked.status == :BLOCKED
        @test SCI._SCAN_DECODE_CALL_COUNT[] == 0
        @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
        @test SCI._DP_CALL_COUNT[] == 0
    finally
        rm(semantic_fixture.root; recursive=true, force=true)
    end
end

@testset "path, member and snapshot identity adversaries block before DP" begin
    cases = (:duplicate, :escape, :symlink, :nonregular, :missing_lf, :sibling_mutation, :inode_replace)
    for case in cases
        fixture = make_fixture()
        try
            report = if case == :duplicate
                infer_with_paths(fixture; node=fixture.t10.edge)
            elseif case == :escape
                infer_with_paths(fixture; edge=joinpath(fixture.root, "..", "outside.tsv"))
            elseif case == :symlink
                real = fixture.t10.edge * ".real"
                mv(fixture.t10.edge, real)
                symlink(real, fixture.t10.edge)
                infer(fixture)
            elseif case == :nonregular
                original = fixture.t10.edge * ".original"
                mv(fixture.t10.edge, original)
                mkdir(fixture.t10.edge)
                infer(fixture)
            elseif case == :missing_lf
                write(fixture.t10.edge, chomp(read(fixture.t10.edge, String)))
                infer(fixture)
            elseif case == :sibling_mutation
                SCI._POST_VALIDATION_HOOK[] = () -> write(joinpath(fixture.t11.output, "models.tsv"), vcat(read(joinpath(fixture.t11.output, "models.tsv")), UInt8('#')))
                infer(fixture)
            else
                sibling = joinpath(fixture.t11.output, "models.tsv")
                original = read(sibling)
                SCI._POST_VALIDATION_HOOK[] = () -> begin
                    moved = sibling * ".old"
                    mv(sibling, moved)
                    write(sibling, original)
                end
                infer(fixture)
            end
            assert_blocked_zero(fixture, report)
        finally
            SCI._POST_VALIDATION_HOOK[] = nothing
            rm(fixture.root; recursive=true, force=true)
        end
    end
end

@testset "terminal receipt statuses do not fall through to graph DP" begin
    for terminal in ("FAIL", "SKIPPED")
        fixture = make_fixture(terminal_status=terminal)
        try
            report = infer(fixture)
            @test report.status == (terminal == "FAIL" ? :FAIL : :SKIPPED)
            @test SCI._DP_CALL_COUNT[] == 0
            @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
    for terminal in ("FAIL", "SKIPPED")
        fixture = make_fixture(dates=("20260101", "20260102", "20260103", "20260104"), shared_references=true, terminal_status=terminal)
        try
            report = infer(fixture)
            @test report.status == (terminal == "FAIL" ? :FAIL : :SKIPPED)
            @test SCI._DP_CALL_COUNT[] == 0
            @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
            @test SCI._SCAN_DECODE_CALL_COUNT[] == (terminal == "FAIL" ? 0 : 20)
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
end

@testset "serialized references remain distinct while caches reuse fit work" begin
    fixture = make_fixture(dates=("20260101", "20260102", "20260103", "20260104"), shared_references=true)
    try
        report = infer(fixture)
        @test report.status == :PASS
        @test length(report.references) == 21
        roles = [(reference.fit_role, reference.scope) for reference in report.references]
        @test count(==( ("partition", "outer_inner")), roles) == 12
        @test count(==( ("partition", "outer_score")), roles) == 4
        @test count(==( ("partition", "final_lodo")), roles) == 4
        @test count(==( ("full_refit", "full_refit")), roles) == 1
        @test [sum(length(block.nodes) for block in reference.blocks) for reference in report.references[1:20]] == fill(3, 20)
        @test sum(length(block.nodes) for block in report.references[21].blocks) == 12
        strict_probability_invariants(report)
        @test SCI._FACTOR_BUILD_CALL_COUNT[] == 40
        @test SCI._SCAN_DECODE_CALL_COUNT[] == 20
        @test SCI._DP_CALL_COUNT[] == 48
    finally
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "out-of-range gaps split topology and preserve unary blocks" begin
    fixture = make_split_fixture()
    try
        report = infer(fixture)
        @test report.status == :PASS
        @test length(report.references) == 1
        @test length(report.references[1].blocks) == 2
        @test sort([length(block.nodes) for block in report.references[1].blocks]) == [1, 2]
        @test sort([length(block.factors) for block in report.references[1].blocks]) == [0, 1]
        strict_probability_invariants(report)
    finally
        rm(fixture.root; recursive=true, force=true)
    end
    tampered = make_split_fixture()
    try
        rewrite_tsv_row!(tampered.t10.edge, 1, row -> (row[11] = "0.1"))
        refresh_t10_bindings!(tampered)
        refresh_t11_input_bindings!(tampered)
        assert_blocked_zero(tampered, infer(tampered))
    finally
        rm(tampered.root; recursive=true, force=true)
    end
end

@testset "Todo10 and Todo11 adversarial bindings fail closed" begin
    todo10_cases = [
        ("edge_header", fixture -> begin
            lines = split(read(fixture.t10.edge, String), '\n'; keepempty=true)
            fields = split(lines[1], '\t'; keepempty=true); fields[1] = "wrong"; lines[1] = join(fields, '\t')
            write(fixture.t10.edge, join(lines, '\n'))
        end),
        ("edge_duplicate", fixture -> begin
            text = read(fixture.t10.edge, String); lines = split(text, '\n'; keepempty=true)
            insert!(lines, length(lines) - 1, lines[2]); write(fixture.t10.edge, join(lines, '\n')); refresh_t10_bindings!(fixture)
        end),
        ("edge_endpoint_t", fixture -> begin
            rewrite_tsv_row!(fixture.t10.edge, 1, row -> (row[4] = fmt(1.1))); refresh_t10_bindings!(fixture)
        end),
        ("edge_correlation_range", fixture -> begin
            rewrite_tsv_row!(fixture.t10.edge, 1, row -> (row[11] = "1.1")); refresh_t10_bindings!(fixture)
        end),
        ("node_status", fixture -> begin
            rewrite_tsv_row!(fixture.t10.node, 1, row -> (row[5] = "isolated")); refresh_t10_bindings!(fixture)
        end),
        ("receipt_extra_key", fixture -> begin
            receipt = TOML.parse(read(fixture.t10.receipt, String)); receipt["unexpected"] = "x"
            write(fixture.t10.receipt, sprint(TOML.print, receipt))
        end),
        ("receipt_header_binding", fixture -> begin
            receipt = TOML.parse(read(fixture.t10.receipt, String)); receipt["edge_header"] = "wrong"
            write(fixture.t10.receipt, sprint(TOML.print, receipt))
        end),
    ]
    for (name, mutate) in todo10_cases
        fixture = make_fixture()
        try
            mutate(fixture)
            assert_blocked_zero(fixture, infer(fixture))
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
    todo11_cases = [
        ("models_report", fixture -> begin
            rewrite_tsv_row!(joinpath(fixture.t11.output, "models.tsv"), 1, row -> (row[4] = "tampered")); refresh_t11_bindings!(fixture; include_handoff=false)
        end),
        ("partition_report", fixture -> begin
            fixture = fixture
        end),
        ("ess_report", fixture -> begin
            fixture = fixture
        end),
        ("handoff_model_header", fixture -> begin
            lines = split(read(fixture.t11.models, String), '\n'; keepempty=true)
            fields = split(lines[1], '\t'; keepempty=true); fields[1] = "wrong"; lines[1] = join(fields, '\t')
            write(fixture.t11.models, join(lines, '\n')); refresh_t11_bindings!(fixture)
        end),
        ("handoff_unary_row", fixture -> begin
            rewrite_tsv_row!(fixture.t11.unary, 1, row -> (row[10] = "0.9")); refresh_t11_bindings!(fixture)
        end),
        ("receipt_extra_key", fixture -> begin
            receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["unexpected"] = "x"
            write(fixture.t11.receipt, sprint(TOML.print, receipt))
        end),
        ("receipt_hash_key", fixture -> begin
            receipt = TOML.parse(read(fixture.t11.receipt, String)); delete!(receipt["hashes"], "plan_sha256")
            write(fixture.t11.receipt, sprint(TOML.print, receipt))
        end),
        ("receipt_pass_fail_inconsistent", fixture -> begin
            receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["status"] = "FAIL"
            write(fixture.t11.receipt, sprint(TOML.print, receipt))
        end),
        ("receipt_pass_skipped_inconsistent", fixture -> begin
            receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["status"] = "SKIPPED"
            write(fixture.t11.receipt, sprint(TOML.print, receipt))
        end),
    ]
    for (name, mutate) in todo11_cases
        fixture = make_fixture(shared_references=(name == "partition_report" || name == "ess_report"), dates=("20260101", "20260102", "20260103", "20260104"))
        try
            if name == "partition_report"
                rewrite_tsv_row!(joinpath(fixture.t11.output, "partitions.tsv"), 1, row -> (row[12] = repeat("0", 64))); refresh_t11_bindings!(fixture; include_handoff=false)
            elseif name == "ess_report"
                rewrite_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1, row -> (row[16] = repeat("0", 64))); refresh_t11_bindings!(fixture; include_handoff=false)
            else
                mutate(fixture)
            end
            assert_blocked_zero(fixture, infer(fixture))
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
end

@testset "factor report hash is deterministic and sensitive" begin
    baseline_fixture = make_fixture()
    mutated_fixture = make_fixture()
    try
        baseline = infer(baseline_fixture)
        rewrite_tsv_row!(mutated_fixture.t11.transforms, 1, row -> (row[19] = "0.02"))
        refresh_t11_bindings!(mutated_fixture)
        mutated = infer(mutated_fixture)
        @test baseline.status == :PASS
        @test mutated.status == :PASS
        @test baseline.result_sha256 != mutated.result_sha256
        @test length(mutated.references[1].blocks[1].factors) == 2
    finally
        rm(baseline_fixture.root; recursive=true, force=true)
        rm(mutated_fixture.root; recursive=true, force=true)
    end
end

@testset "public unary -Inf and zero-factor semantics" begin
    infeasible = make_fixture()
    zero_factor = make_fixture()
    try
        for row_number in 1:12
            rewrite_tsv_row!(infeasible.t11.unary, row_number,
                            row -> begin row[10] = "0"; row[11] = "0" end)
        end
        rebind_unary_probability!(infeasible)
        failed = infer(infeasible)
        @test failed.status == :FAIL
        @test failed.reason == :unary_failure
        @test SCI._DP_CALL_COUNT[] == 0
        rewrite_tsv_row!(zero_factor.t11.models, 1, row -> begin
            for category in SCI._GRAPH_CATEGORY_ORDER, output in SCI._GRAPH_OUTPUT_NAMES
                row[findfirst(==("conditional_$(category)_mean_$(output)"), SCI._GRAPH_MODEL_HEADER)] = "1"
            end
        end)
        refresh_t11_bindings!(zero_factor)
        report = infer(zero_factor)
        @test report.status == :BLOCKED
        @test SCI._DP_CALL_COUNT[] == 0
        @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
        @test SCI._SCAN_DECODE_CALL_COUNT[] == 0
    finally
        rm(infeasible.root; recursive=true, force=true)
        rm(zero_factor.root; recursive=true, force=true)
    end
end

@testset "public reversed model returns natural-order results" begin
    natural_fixture = make_fixture()
    reversed_fixture = make_fixture(reversed=true)
    try
        natural = infer(natural_fixture)
        reversed = infer(reversed_fixture)
        @test natural.status == :PASS
        @test reversed.status == :PASS
        strict_probability_invariants(natural)
        strict_probability_invariants(reversed)
        natural_nodes = natural.references[1].blocks[1].nodes
        reversed_nodes = reversed.references[1].blocks[1].nodes
        @test [node.lobe for node in reversed_nodes] == [1, 2, 3]
        @test [node.label for node in reversed_nodes] == [node.label for node in natural_nodes]
        @test all(all(isapprox(reversed_nodes[index].marginals[state], natural_nodes[index].marginals[state];
                               atol=1e-12, rtol=1e-12) for state in 1:2)
                  for index in eachindex(natural_nodes))
    finally
        rm(natural_fixture.root; recursive=true, force=true)
        rm(reversed_fixture.root; recursive=true, force=true)
    end
end

@testset "correction3 diagnostic TSV semantics fail before inference" begin
    cases = [
        ("models.tsv", false, row -> (row[8] = repeat("0", 64))),
        ("partitions.tsv", true, row -> (row[11] = repeat("0", 64))),
        ("ess.tsv", false, row -> (row[6] = "4")),
        ("starts.tsv", false, row -> (row[7] = "1")),
        ("scores.tsv", true, row -> (row[7] = "missing.sxm")),
        ("shuffle.tsv", true, row -> (row[6] = "500")),
        ("bootstrap.tsv", false, row -> (row[13] = repeat("b", 64))),
    ]
    for (name, shared, mutate) in cases
        fixture = make_fixture(shared_references=shared,
                               dates=shared ? ("20260101", "20260102", "20260103", "20260104") : ("20260101",))
        try
            rewrite_tsv_row!(joinpath(fixture.t11.output, name), 1, mutate)
            refresh_t11_bindings!(fixture; include_handoff=false)
            assert_blocked_zero(fixture, infer(fixture))
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
end

@testset "correction3 status, transform and held-out topology adversaries" begin
    mixed = make_fixture()
    try
        rewrite_tsv_row!(mixed.t11.models, 1, row -> begin
            row[10] = "FAIL"
            row[11] = "numerical_failure"
        end)
        refresh_t11_bindings!(mixed)
        assert_blocked_zero(mixed, infer(mixed))
    finally
        rm(mixed.root; recursive=true, force=true)
    end

    for receipt_status in ("SKIPPED", "FAIL")
        fixture = make_fixture()
        try
            receipt = TOML.parse(read(fixture.t11.receipt, String))
            receipt["status"] = receipt_status
            receipt["reason"] = receipt_status == "SKIPPED" ? "conditional_mechanism_not_admitted" : "numerical_failure"
            write(fixture.t11.receipt, sprint(TOML.print, receipt))
            assert_blocked_zero(fixture, infer(fixture))
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end

    duplicate = make_fixture(shared_references=true,
                             dates=("20260101", "20260102", "20260103", "20260104"))
    try
        lines = split(read(duplicate.t11.transforms, String), '\n'; keepempty=true)
        divergent = split(lines[2], '\t'; keepempty=true)
        divergent[4] = "divergent_scope"
        write(duplicate.t11.transforms, join(vcat(lines[1:2], join(divergent, '\t'), [""]), '\n'))
        refresh_t11_bindings!(duplicate)
        assert_blocked_zero(duplicate, infer(duplicate))
    finally
        rm(duplicate.root; recursive=true, force=true)
    end

    leakage = make_fixture(shared_references=true,
                           dates=("20260101", "20260102", "20260103", "20260104"))
    try
        rewrite_tsv_row!(leakage.t11.models, 1, row -> (row[8] = "20260101,20260103,20260104"))
        refresh_t11_bindings!(leakage)
        assert_blocked_zero(leakage, infer(leakage))
    finally
        rm(leakage.root; recursive=true, force=true)
    end

    segment = make_fixture()
    try
        rewrite_tsv_row!(segment.t10.node, 2, row -> (row[4] = segment_id(segment.t10.file, 2, 2)))
        refresh_t10_bindings!(segment)
        assert_blocked_zero(segment, infer(segment))
    finally
        rm(segment.root; recursive=true, force=true)
    end
end

@testset "correction3 dependency hash bindings" begin
    for field in ("universe_sha256", "input_sha256", "normalization_sha256")
        fixture = make_fixture()
        try
            receipt = TOML.parse(read(fixture.t11.receipt, String))
            receipt["hashes"][field] = repeat("0", 64)
            write(fixture.t11.receipt, sprint(TOML.print, receipt))
            assert_blocked_zero(fixture, infer(fixture))
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
end

@testset "correction4 real serializer fixtures and staged evidence" begin
    happy = make_fixture()
    sparse = make_fixture(terminal_status="SKIPPED")
    staged_fail = make_fixture(terminal_status="FAIL")
    try
        @test infer(happy).status == :PASS
        @test infer(sparse).status == :SKIPPED
        @test infer(staged_fail).status == :FAIL
        @test length(readlines(joinpath(sparse.t11.output, "starts.tsv"))) == 1
        model_lines = split(read(happy.t11.models, String), '\n'; keepempty=true)
        model_rows = [split(line, '\t'; keepempty=true) for line in model_lines[2:(end - 1)]]
        @test length(model_rows) == 11
        @test count(row[5] == "outer_inner" for row in model_rows) == 6
        @test count(row[5] == "final_lodo" for row in model_rows) == 4
        @test sum(parse(Int, row[14]) for row in model_rows) == 21
        @test any(row[5] == "outer_inner" && row[6] == "20260101" && row[7] == "20260102" && row[8] == "20260103,20260104" for row in model_rows)
        @test any(occursin("outer_inner|20260102|20260101|", row[16]) for row in model_rows)
    finally
        rm(happy.root; recursive=true, force=true)
        rm(sparse.root; recursive=true, force=true)
        rm(staged_fail.root; recursive=true, force=true)
    end

    varying = make_fixture()
    try
        rewrite_tsv_row!(joinpath(varying.t11.output, "ess.tsv"), 1, row -> (row[10] = "2"))
        rewrite_tsv_row!(joinpath(varying.t11.output, "ess.tsv"), 2, row -> (row[11] = "0"))
        refresh_t11_bindings!(varying)
        @test infer(varying).status == :PASS
    finally
        rm(varying.root; recursive=true, force=true)
    end

    semantic_cases = [
        ("ess_invariant_dates", "ess.tsv", row -> (row[7] = "99")),
        ("ess_duplicate_category", "ess.tsv", row -> (row[6] = "1")),
        ("sparse_ess_marked_sufficient", "ess.tsv", row -> (row[14] = "true")),
        ("start_trace_reason", "starts.tsv", row -> (row[10] = "tampered_reason")),
        ("score_gain", "scores.tsv", row -> (row[11] = "2")),
        ("shuffle_value", "shuffle.tsv", row -> (row[7] = "2")),
        ("bootstrap_value", "bootstrap.tsv", row -> (row[6] = "2")),
        ("model_result", "models.tsv", row -> (row[8] = repeat("0", 64))),
    ]
    for (name, artifact, mutate) in semantic_cases
        fixture = make_fixture(terminal_status=(name == "sparse_ess_marked_sufficient" ? "SKIPPED" : "PASS"))
        try
            rewrite_tsv_row!(joinpath(fixture.t11.output, artifact), 1, mutate)
            refresh_t11_bindings!(fixture; include_handoff=false)
            assert_blocked_zero(fixture, infer(fixture))
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
end

@testset "correction6 producer boundary and named publication matrix" begin
    producer = make_producer_fixture()
    try
        report = infer(producer)
        @test report.status == :PASS
        @test report.reason == :graph_inference
        @test length(report.references) == 114
        @test length(producer.t11.report.models) == 2
        for model_id in ("C1", "C2")
            result = producer.t11.report.models[model_id]
            @test result.status == "PASS"
            @test length(result.outer_folds) == 7
            @test all(length(fold.inner_gate.partitions) == 6 for fold in result.outer_folds)
            @test length(result.final_gate.partitions) == 7
            @test length(result.full_fit.training_dates) == 7
        end
        model_rows = [split(line, '\t'; keepempty=true)
                      for line in split(read(producer.t11.models, String), '\n')[2:end-1]]
        @test count(row -> row[4] == "C1", model_rows) == 29
        @test count(row -> row[4] == "C2", model_rows) == 29
        @test all(parse(Int, row[14]) == 1 for row in model_rows if row[4] == "C1" && row[2] == "full_refit")
        @test length(readlines(joinpath(producer.t11.output, "partitions.tsv"))) - 1 == 112
        @test length(readlines(joinpath(producer.t11.output, "ess.tsv"))) - 1 == 456
        @test length(readlines(joinpath(producer.t11.output, "starts.tsv"))) - 1 == 684
        @test length(readlines(joinpath(producer.t11.output, "scores.tsv"))) - 1 == 2800
        @test length(readlines(joinpath(producer.t11.output, "shuffle.tsv"))) - 1 == 56000
        @test length(readlines(joinpath(producer.t11.output, "bootstrap.tsv"))) - 1 == 8000
        fit_hashes = Dict(row[3] => row for row in model_rows)
        for model_id in ("C1", "C2")
            refs = [split(value, '|'; keepempty=true)
                    for row in model_rows if row[4] == model_id
                    for value in split(row[16], ';')]
            @test count(row -> row[3] == "outer_inner", refs) == 42
            @test count(row -> row[3] == "outer_score", refs) == 7
            @test count(row -> row[3] == "final_lodo", refs) == 7
            @test count(row -> row[1] == "full_refit", refs) == 1
            @test length(unique(row[6] for row in refs)) == 29
            @test only(row[6] for row in refs if row[3] == "outer_inner" &&
                       row[4] == "20260101" && row[5] == "20260102") ==
                  only(row[6] for row in refs if row[3] == "outer_inner" &&
                       row[4] == "20260102" && row[5] == "20260101")
            for date in producer.t10.dates
                @test only(row[6] for row in refs if row[3] == "outer_score" && row[4] == date) ==
                      only(row[6] for row in refs if row[3] == "final_lodo" && row[5] == date)
            end
        end
    finally
        rm(producer.root; recursive=true, force=true)
    end

    shape = make_fixture(dates=("20260101", "20260102", "20260103", "20260104"),
                         shared_references=true)
    try
        @test length(readlines(shape.t11.models)) - 1 == 11
        @test length(readlines(joinpath(shape.t11.output, "partitions.tsv"))) - 1 == 20
        @test count(row -> occursin("outer_inner|20260101|20260102|",
                                    row), split(read(shape.t11.models, String), '\n')) == 1
    finally
        rm(shape.root; recursive=true, force=true)
    end

    skipped = make_producer_fixture(dates=["20260101", "20260102", "20260103", "20260104"],
                                    require_pass=false)
    try
        report = infer(skipped)
        @test report.status == :SKIPPED
        @test report.reason == :unary_only
        @test SCI._DP_CALL_COUNT[] == 0
        @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
        @test SCI._SCAN_DECODE_CALL_COUNT[] == 32
        @test length(readlines(joinpath(skipped.t11.output, "partitions.tsv"))) == 1
        @test length(readlines(joinpath(skipped.t11.output, "ess.tsv"))) - 1 == 8
        @test all(startswith(line, "1\t")
                  for line in readlines(joinpath(skipped.t11.output, "bootstrap.tsv"))[2:end])
    finally
        rm(skipped.root; recursive=true, force=true)
    end

    impossible = make_producer_fixture(dates=["20260101", "20260102", "20260103", "20260104"],
                                       require_pass=false)
    try
        rewrite_tsv_row!(joinpath(impossible.t11.output, "models.tsv"), 1,
                         row -> (row[3] = "PASS"; row[4] = "ok"))
        refresh_t11_bindings!(impossible)
        assert_blocked_zero(impossible, infer(impossible))
    finally
        rm(impossible.root; recursive=true, force=true)
    end

    matrix_cases = [
        "fixed_authority", "seven_path_distinctness", "path_escape", "symlink",
        "file_set", "snapshot", "pre_inference_toctou", "todo10_headers_types_order",
        "todo10_raw_hashes_topology_gaps_splits_segments_boundaries_receipt",
        "c1_c2_complete_graph", "shared_fit_reuse", "full_refit_distinct",
        "k_lt_7_skipped", "k_lt_7_impossible_pass", "handoff_coverage",
        "category_order_tied_means_scale", "prediction_status_references",
        "ess_support_kish_ratio_reason", "staged_start_aliases",
        "score_identity_gain_counts_order_raw_scan_date_control_sparse_hash",
        "shuffle_scope_seed_order_run_state_control_quantile_pass_score_alias",
        "bootstrap_values_lower_predicates_status_reason_gate_aggregation",
        "model_tuple_count_final_full_status_reason_result",
        "receipt_status_reason_result_input_universe_artifact",
        "factor_category_density_cache_unary_reduction", "exhaustive_dp_reversal",
        "split_singleton_api_label_free", "density_failure", "dp_failure",
        "unexpected_internal_failure",
    ]
    @test length(matrix_cases) == 30
    @test length(unique(matrix_cases)) == length(matrix_cases)
end

@testset "correction7 executable deep matrix" begin
    canonical = make_producer_fixture()
    matrix_rows = Vector{Vector{String}}()
    cases = NamedTuple[]
    add_case(action, name, family, mutation, refreshed, layer, status, reason; invoke=infer) =
        push!(cases, (; name, family, mutation, refreshed, layer, status, reason,
                       action, invoke))

    _, shuffle_rows = tsv_rows(joinpath(canonical.t11.output, "shuffle.tsv"))
    running_shuffle = filter(row -> row[7] != "NA", shuffle_rows)
    @test length(unique(row[8] for row in running_shuffle)) > 1
    @test length(unique(row[8] for row in running_shuffle)) < length(running_shuffle)

    add_case("authority_unexpected_sibling", "authority_paths",
             "add todo11/unexpected.tsv", "none", "sibling_set", :BLOCKED, :sibling_set) do fixture
        write(joinpath(fixture.t11.output, "unexpected.tsv"), "unexpected\n")
    end
    add_case(f -> nothing, "authority_duplicate_supplied_member", "authority_paths",
             "supply node_segments as edge_observations", "none", "path identity", :BLOCKED, :source_collision;
             invoke=f -> infer_with_paths(f; node=f.t10.edge))
    add_case(f -> nothing, "authority_path_escape", "authority_paths",
             "supply path outside root", "none", "path resolution", :BLOCKED, :path_escape;
             invoke=f -> infer_with_paths(f; edge=joinpath(f.root, "..", "outside.tsv")))
    add_case("authority_symlink_input", "authority_paths",
             "replace edge input with symlink", "none", "path identity", :BLOCKED, :symlink_rejected) do fixture
        real = fixture.t10.edge * ".real"
        mv(fixture.t10.edge, real)
        symlink(real, fixture.t10.edge)
    end
    add_case("authority_directory_input", "authority_paths",
             "replace edge input with directory", "none", "path identity", :BLOCKED, :missing_input) do fixture
        real = fixture.t10.edge * ".original"
        mv(fixture.t10.edge, real)
        mkdir(fixture.t10.edge)
    end
    add_case("authority_missing_final_lf", "authority_paths",
             "remove edge table final LF", "none", "snapshot parser", :BLOCKED, :missing_final_lf) do fixture
        write(fixture.t10.edge, chomp(read(fixture.t10.edge, String)))
    end
    add_case("authority_pre_inference_snapshot_replacement", "authority_paths",
             "replace edge after validation", "none", "pre-inference resnapshot", :BLOCKED, :input_changed) do fixture
        SCI._POST_VALIDATION_HOOK[] = () -> write(fixture.t10.edge, vcat(read(fixture.t10.edge), UInt8('#')))
    end
    add_case("authority_return_time_snapshot_race", "authority_paths",
             "replace edge at return snapshot", "none", "return-time resnapshot", :BLOCKED, :input_changed) do fixture
        SCI._SNAPSHOT_RETURN_HOOK[] = () -> write(fixture.t10.edge, vcat(read(fixture.t10.edge), UInt8('#')))
    end

    add_case("todo10_invalid_canonical_type", "todo10",
             "node lobe 1 -> 1.0", "T10 receipt/input", "Todo10 canonical parser", :BLOCKED, :invalid_key) do fixture
        mutate_tsv_row!(fixture.t10.node, 1) do row
            row[2] = "1.0"
        end
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_row_order_swap", "todo10",
             "swap first two edge rows", "T10 receipt/input", "Todo10 order", :BLOCKED, :topology_mismatch) do fixture
        header, rows = tsv_rows(fixture.t10.edge)
        rows[1], rows[2] = rows[2], rows[1]
        write_tsv_rows(fixture.t10.edge, header, rows)
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_raw_row_hash_rebound", "todo10",
             "change edge observation and refresh outer receipts only", "T10/T11 artifact", "transform raw binding", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(fixture.t10.edge, 1) do row
            row[11] = fmt(parse(Float64, row[11]) + 0.01)
        end
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_duplicate_node", "todo10",
             "duplicate first node row", "T10 receipt", "Todo10 duplicate node", :BLOCKED, :duplicate_key) do fixture
        header, rows = tsv_rows(fixture.t10.node)
        insert!(rows, 2, copy(rows[1]))
        write_tsv_rows(fixture.t10.node, header, rows)
        refresh_t10_bindings!(fixture)
        refresh_t11_input_bindings!(fixture)
    end
    add_case("todo10_missing_edge", "todo10",
             "remove first edge row", "T10 receipt/input", "Todo10 topology", :BLOCKED, :topology_mismatch) do fixture
        header, rows = tsv_rows(fixture.t10.edge)
        deleteat!(rows, 1)
        write_tsv_rows(fixture.t10.edge, header, rows)
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_nonconsecutive_endpoint", "todo10",
             "retarget first edge to nonconsecutive lobe", "T10 receipt/input", "Todo10 topology", :BLOCKED, :topology_mismatch) do fixture
        mutate_tsv_row!(fixture.t10.edge, 1) do row
            row[3] = "3"
            row[5] = fmt(2.0)
            row[6] = fmt(1.0)
        end
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_invalid_gap", "todo10",
             "change eligible gap without endpoint change", "T10 receipt/input", "Todo10 gap", :BLOCKED, :topology_mismatch) do fixture
        mutate_tsv_row!(fixture.t10.edge, 1) do row
            row[6] = fmt(0.1)
        end
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_split_precedence", "todo10",
             "make coordinated out-of-range split with wrong reason", "T10 receipt/input", "Todo10 split precedence", :BLOCKED, :split_reason) do fixture
        node_header, node_rows = tsv_rows(fixture.t10.node)
        edge_header, edge_rows = tsv_rows(fixture.t10.edge)
        file = fixture.t10.files[1]
        file_nodes = filter(row -> row[1] == file, node_rows)
        left_segment = segment_id(file, 1, 1)
        right_segment = segment_id(file, 2, 6)
        shifted = Dict(row[2] => (parse(Float64, row[3]) + (parse(Int, row[2]) >= 2 ? 0.5 : 0.0))
                       for row in file_nodes)
        for row in file_nodes
            row[3] = fmt(shifted[row[2]])
            row[4] = parse(Int, row[2]) == 1 ? left_segment : right_segment
        end
        first_node = only(filter(row -> row[2] == "1", file_nodes))
        second_node = only(filter(row -> row[2] == "2", file_nodes))
        first_node[5] = "isolated"; first_node[7] = "gap_out_of_range"
        second_node[6] = "gap_out_of_range"
        file_edges = filter(row -> row[1] == file, edge_rows)
        for row in file_edges
            left, right = row[2], row[3]
            row[4] = fmt(shifted[left]); row[5] = fmt(shifted[right]); row[6] = fmt(shifted[right] - shifted[left])
            row[7] = right == "2" ? right_segment : right_segment
            row[8] = right_segment
        end
        first_edge = only(filter(row -> row[2] == "1" && row[3] == "2", file_edges))
        first_edge[6] = fmt(1.0)
        first_edge[7] = left_segment
        first_edge[8] = right_segment
        first_edge[9] = "split"
        first_edge[10] = "other"
        write_tsv_rows(fixture.t10.node, node_header, node_rows)
        write_tsv_rows(fixture.t10.edge, edge_header, edge_rows)
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_segment_id", "todo10",
             "change edge segment id", "T10 receipt/input", "Todo10 segment", :BLOCKED, :segment_mismatch) do fixture
        mutate_tsv_row!(fixture.t10.edge, 1) do row
            row[7] = repeat("0", 16)
        end
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_boundary_reason", "todo10",
             "change first node boundary reason", "T10 receipt/input", "Todo10 boundary", :BLOCKED, :boundary_mismatch) do fixture
        mutate_tsv_row!(fixture.t10.node, 1) do row
            row[6] = "eligible"
        end
        refresh_todo10_and_t11!(fixture)
    end
    add_case("todo10_receipt_count", "todo10",
             "decrement edge receipt row count", "none", "Todo10 receipt", :BLOCKED, :row_count) do fixture
        receipt = TOML.parse(read(fixture.t10.receipt, String))
        receipt["edge_row_count"] -= 1
        write(fixture.t10.receipt, sprint(TOML.print, receipt))
    end
    add_case("todo10_receipt_header", "todo10",
             "change edge header receipt binding", "none", "Todo10 receipt", :BLOCKED, :receipt_binding) do fixture
        receipt = TOML.parse(read(fixture.t10.receipt, String))
        receipt["edge_header"] = "wrong"
        write(fixture.t10.receipt, sprint(TOML.print, receipt))
    end
    add_case("todo10_receipt_input_cross_binding", "todo10",
             "change model config receipt binding", "none", "Todo10 receipt", :BLOCKED, :hash_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t10.receipt, String))
        receipt["model_config_sha256"] = repeat("0", 64)
        write(fixture.t10.receipt, sprint(TOML.print, receipt))
    end

    add_case("graph_exact_model_set", "model_graph",
             "coherently rename C2 report/handoff/receipt identities to C3", "artifact/graph/result", "C1/C2 model set", :BLOCKED, :reference_mismatch) do fixture
        rename_model_identity!(fixture, "C2", "C3")
    end
    add_case("graph_missing_outer_inner", "model_graph",
             "remove A/B outer-inner logical reference and rehash graph wrappers", "artifact/graph", "nested graph completeness", :BLOCKED, :row_count) do fixture
        header, rows = tsv_rows(fixture.t11.models)
        ref_col = findfirst(==("references"), SCI._GRAPH_MODEL_HEADER)
        count_col = findfirst(==("reference_count"), SCI._GRAPH_MODEL_HEADER)
        hash_col = findfirst(==("references_sha256"), SCI._GRAPH_MODEL_HEADER)
        index = only(findall(row -> occursin("partition|C1|outer_inner|20260101|20260102|", row[ref_col]), rows))
        refs = split(rows[index][ref_col], ';')
        remove = only(filter(ref -> occursin("partition|C1|outer_inner|20260101|20260102|", ref), refs))
        filter!(ref -> ref != remove, refs)
        rows[index][ref_col] = join(refs, ';')
        rows[index][count_col] = string(length(refs))
        rows[index][hash_col] = hash_lines([join(split(ref, '|'), '\t') for ref in refs])
        rows[index][findfirst(==( "outer_date"), SCI._GRAPH_MODEL_HEADER)] = "20260102"
        rows[index][findfirst(==( "target_date"), SCI._GRAPH_MODEL_HEADER)] = "20260101"
        write_tsv_rows(fixture.t11.models, header, rows)
        for name in ("partitions.tsv", "ess.tsv", "starts.tsv", "scores.tsv", "shuffle.tsv")
            path = joinpath(fixture.t11.output, name)
            h, values = tsv_rows(path)
            filter!(row -> !(row[2] == "C1" && row[3] == "outer_inner" &&
                             row[4] == "20260101" && row[5] == "20260102"), values)
            write_tsv_rows(path, h, values)
        end
        refresh_t11_wrappers!(fixture; graph_lists=true)
    end
    add_case("graph_ab_shared_fit_split", "fit_aliases",
             "clone B/A cross-pair fit and retarget its handoff rows", "artifact/graph", "A/B shared-fit reuse", :BLOCKED, :reference_mismatch) do fixture
        header, rows = tsv_rows(fixture.t11.models)
        ref_col = findfirst(==( "references"), SCI._GRAPH_MODEL_HEADER)
        ref = only(filter(ref -> occursin("partition|C1|outer_inner|20260102|20260101|", ref),
                          vcat([split(row[ref_col], ';') for row in rows]...)))
        split_fit_alias!(fixture, "", String(ref))
    end
    add_case("graph_outer_final_alias_split", "fit_aliases",
             "clone final-L ODO fit and retarget its handoff rows", "artifact/graph", "outer/final fit reuse", :BLOCKED, :reference_mismatch) do fixture
        header, rows = tsv_rows(fixture.t11.models)
        ref_col = findfirst(==("references"), SCI._GRAPH_MODEL_HEADER)
        ref = only(filter(ref -> occursin("partition|C1|final_lodo|NA|20260101|", ref),
                          vcat([split(row[ref_col], ';') for row in rows]...)))
        split_fit_alias!(fixture, "", String(ref))
    end
    add_case("graph_outer_tuple", "model_graph",
             "make one outer-inner reference target equal its outer date", "artifact/graph", "outer-fold tuple", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            index = findfirst(==( "target_date"), SCI._GRAPH_MODEL_HEADER)
            outer = row[findfirst(==( "outer_date"), SCI._GRAPH_MODEL_HEADER)]
            row[index] = outer
        end
        refresh_t11_wrappers!(fixture; graph_lists=true)
    end

    add_case("handoff_missing_unary", "handoff",
             "remove one fitted unary node row", "artifact/graph", "unary coverage", :BLOCKED, :unary_coverage) do fixture
        path = fixture.t11.unary
        header, rows = tsv_rows(path)
        deleteat!(rows, 1)
        write_tsv_rows(path, header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_duplicate_unary", "handoff",
             "duplicate one fitted unary node row", "artifact/graph", "unary duplicate key", :BLOCKED, :duplicate_key) do fixture
        path = fixture.t11.unary
        header, rows = tsv_rows(path)
        insert!(rows, 2, copy(rows[1]))
        write_tsv_rows(path, header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_missing_transform", "handoff",
             "remove one fitted transform edge row", "artifact/graph", "transform coverage", :BLOCKED, :transform_coverage) do fixture
        path = fixture.t11.transforms
        header, rows = tsv_rows(path)
        deleteat!(rows, 1)
        write_tsv_rows(path, header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_duplicate_transform", "handoff",
             "duplicate one fitted transform edge row", "artifact/graph", "transform duplicate key", :BLOCKED, :duplicate_key) do fixture
        path = fixture.t11.transforms
        header, rows = tsv_rows(path)
        insert!(rows, 2, copy(rows[1]))
        write_tsv_rows(path, header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_category_order", "handoff",
             "change graph category order", "artifact/graph", "model schema", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "category_order"), SCI._GRAPH_MODEL_HEADER)] = "00,10,01,11"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_state_order", "handoff",
             "change graph state order", "artifact/graph", "model schema", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "state_order"), SCI._GRAPH_MODEL_HEADER)] = "1,0"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_tied_mixed_mean", "handoff",
             "break exact 01/10 conditional mean tie", "artifact/graph", "density handoff", :BLOCKED, :model_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "conditional_01_mean_corr_fwd"), SCI._GRAPH_MODEL_HEADER)] = "0.123"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_scale_floor", "handoff",
             "set null scale below positive floor", "artifact/graph", "density handoff", :BLOCKED, :model_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "null_scale_ff"), SCI._GRAPH_MODEL_HEADER)] = "0"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_conditional_scale_floor", "handoff",
             "set conditional scale below positive floor", "artifact/graph", "conditional density handoff", :BLOCKED, :model_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "conditional_scale_ff"), SCI._GRAPH_MODEL_HEADER)] = "0"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_invalid_condition", "handoff",
             "change fixed Student-t condition", "artifact/graph", "model condition", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "student_t_nu"), SCI._GRAPH_MODEL_HEADER)] = "7"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_selected_start", "handoff",
             "change selected start index", "artifact/graph", "selected-start binding", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "selected_start"), SCI._GRAPH_MODEL_HEADER)] = "2"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_unary_probability", "handoff",
             "change one unary probability without alias rebind", "artifact/graph", "unary probability binding", :BLOCKED, :unary_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.unary, 1) do row
            row[11] = fmt(parse(Float64, row[11]) + 0.01)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_unary_identity", "handoff",
             "change one unary node identity", "artifact/graph", "unary identity", :BLOCKED, :unary_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.unary, 1) do row
            row[8] = replace(row[8], "|1|" => "|2|")
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_unary_fit_alias", "handoff",
             "change one unary fit digest", "artifact/graph", "unary fit alias", :BLOCKED, :unary_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.unary, 1) do row
            row[12] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_transform_raw_digest", "handoff",
             "change transform raw-row digest", "artifact/graph", "transform raw binding", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.transforms, 1) do row
            row[16] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_transform_status_reason", "handoff",
             "mark finite transform ineligible", "artifact/graph", "transform status/reason", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.transforms, 1) do row
            row[19] = "NA"
            row[20] = "NA"
            row[21] = "unavailable"
            row[22] = "none"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("handoff_transform_prediction_na", "handoff",
             "blank one eligible transform prediction", "artifact/graph", "transform prediction finite/NA", :BLOCKED, :transform_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.transforms, 1) do row
            row[19] = "NA"
        end
        refresh_t11_wrappers!(fixture)
    end

    add_case("ess_missing_category", "ess",
             "remove one ESS row", "artifact/graph", "ESS category coverage", :BLOCKED, :row_count) do fixture
        path = joinpath(fixture.t11.output, "ess.tsv")
        header, rows = tsv_rows(path); deleteat!(rows, 1)
        write_tsv_rows(path, header, rows); refresh_t11_wrappers!(fixture)
    end
    add_case("ess_duplicate_category", "ess",
             "duplicate one ESS category row", "artifact/graph", "ESS duplicate key", :BLOCKED, :duplicate_key) do fixture
        path = joinpath(fixture.t11.output, "ess.tsv")
        header, rows = tsv_rows(path); insert!(rows, 2, copy(rows[1]))
        write_tsv_rows(path, header, rows); refresh_t11_wrappers!(fixture)
    end
    add_case("ess_invariant_dates", "ess",
             "change one ESS date count", "artifact/graph", "ESS invariant", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[7] = "99"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("ess_invariant_scans_edges", "ess",
             "change one ESS scan count", "artifact/graph", "ESS scan/edge invariant", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[8] = "99"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("ess_kish", "ess",
             "change one Kish ESS value", "artifact/graph", "ESS Kish", :BLOCKED, :support_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[10] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("ess_minimum_support", "ess",
             "change minimum supporting scans", "artifact/graph", "ESS minimum support", :BLOCKED, :support_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[11] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("ess_free_parameters", "ess",
             "change free parameter count", "artifact/graph", "ESS free parameters", :BLOCKED, :support_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[12] = "10"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("ess_effective_ratio", "ess",
             "change effective observation ratio", "artifact/graph", "ESS ratio", :BLOCKED, :support_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[13] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("ess_sufficient_flag", "ess",
             "flip sufficient flag", "artifact/graph", "ESS sufficiency", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[14] = "false"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("ess_first_failing_reason", "ess",
             "change first-failing reason", "artifact/graph", "ESS reason", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "ess.tsv"), 1) do row
            row[15] = "insufficient_edges"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_missing", "starts",
             "remove one start row", "artifact/graph", "start coverage", :BLOCKED, :row_count) do fixture
        path = joinpath(fixture.t11.output, "starts.tsv")
        header, rows = tsv_rows(path); deleteat!(rows, 1)
        write_tsv_rows(path, header, rows); refresh_t11_wrappers!(fixture)
    end
    add_case("starts_duplicate", "starts",
             "duplicate one start row", "artifact/graph", "start duplicate key", :BLOCKED, :duplicate_key) do fixture
        path = joinpath(fixture.t11.output, "starts.tsv")
        header, rows = tsv_rows(path); insert!(rows, 2, copy(rows[1]))
        write_tsv_rows(path, header, rows); refresh_t11_wrappers!(fixture)
    end
    add_case("starts_null_alpha", "starts",
             "change null alpha from NA", "artifact/graph", "start schema", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 1) do row
            row[8] = "0"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_null_index", "starts",
             "change null start index", "artifact/graph", "null start metadata", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 1) do row
            row[7] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_conditional_index", "starts",
             "change conditional index outside 1:5", "artifact/graph", "start schema", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 2) do row
            row[7] = "6"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_conditional_alpha", "starts",
             "change conditional start alpha", "artifact/graph", "conditional start alpha", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 2) do row
            row[8] = "0.1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_status_reason", "starts",
             "change PASS start reason", "artifact/graph", "start status/reason", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 1) do row
            row[10] = "tampered"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_alias_trace", "starts",
             "change opaque trace digest on one alias", "artifact/graph", "same-fit start aliases", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 1) do row
            row[16] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_iterations", "starts",
             "change published iteration count", "artifact/graph", "start iteration summary", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 1) do row
            row[11] = "2"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_objective_initial", "starts",
             "change initial objective summary", "artifact/graph", "start objective summary", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 1) do row
            row[12] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_rejected_summary", "starts",
             "change rejected iteration summary", "artifact/graph", "start rejected summary", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 1) do row
            row[14] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_selected_alpha", "handoff",
             "change selected start alpha", "artifact/graph", "selected start binding", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "selected_alpha"), SCI._GRAPH_MODEL_HEADER)] = "0.25"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_selected_objective", "handoff",
             "change selected start objective", "artifact/graph", "selected objective binding", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "selected_objective"), SCI._GRAPH_MODEL_HEADER)] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_selected_convergence", "handoff",
             "change selected convergence flag", "artifact/graph", "selected convergence binding", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.models, 1) do row
            row[findfirst(==( "selected_converged"), SCI._GRAPH_MODEL_HEADER)] = "false"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_conditional_status", "starts",
             "change conditional start status", "artifact/graph", "conditional start status", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "starts.tsv"), 2) do row
            row[9] = "SKIPPED"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("starts_null_only_coverage", "starts",
             "remove five conditional starts and blank conditional handoff", "artifact/graph", "null-only versus five-start coverage", :BLOCKED, :model_mismatch) do fixture
        start_path = joinpath(fixture.t11.output, "starts.tsv")
        start_header, start_rows = tsv_rows(start_path)
        fit = start_rows[1][17]
        filter!(row -> !(row[17] == fit && row[6] == "conditional"), start_rows)
        write_tsv_rows(start_path, start_header, start_rows)
        model_header, model_rows = tsv_rows(fixture.t11.models)
        model_index = only(findall(row -> row[3] == fit, model_rows))
        for name in vcat(["conditional_$(category)_mean_$(output)"
                          for category in SCI._GRAPH_CATEGORY_ORDER for output in SCI._GRAPH_OUTPUT_NAMES],
                         ["conditional_scale_ff", "conditional_scale_fb", "conditional_scale_bb",
                          "selected_start", "selected_alpha", "selected_objective", "selected_converged"])
            model_rows[model_index][findfirst(==(name), SCI._GRAPH_MODEL_HEADER)] = "NA"
        end
        write_tsv_rows(fixture.t11.models, model_header, model_rows)
        refresh_t11_wrappers!(fixture)
    end

    add_case("score_edge_gain", "scores",
             "change published edge gain", "artifact/graph", "edge gain reconstruction", :BLOCKED, :score_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "scores.tsv"), 1) do row
            row[11] = fmt(parse(Float64, row[11]) + 0.1)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_edge_count", "scores",
             "change edge score count", "artifact/graph", "edge count", :BLOCKED, :row_count) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "scores.tsv"), 1) do row
            row[12] = "2"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_edge_identity", "scores",
             "change held-out score edge file identity", "artifact/graph", "score edge identity", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "scores.tsv"), 1) do row
            row[7] = "20260199_scan01.sxm"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_edge_segment", "scores",
             "change held-out score segment", "artifact/graph", "score segment identity", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "scores.tsv"), 1) do row
            row[10] = repeat("0", 16)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_edge_order", "scores",
             "swap two edge score rows", "artifact/graph", "score key order", :BLOCKED, :row_count) do fixture
        header, rows = tsv_rows(joinpath(fixture.t11.output, "scores.tsv"))
        rows[1], rows[2] = rows[2], rows[1]
        write_tsv_rows(joinpath(fixture.t11.output, "scores.tsv"), header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_scan_mean", "scores",
             "change scan mean", "artifact/graph", "scan aggregation", :BLOCKED, :score_mismatch) do fixture
        header, rows = tsv_rows(joinpath(fixture.t11.output, "scores.tsv"))
        index = findfirst(row -> row[6] == "scan", rows)
        rows[index][11] = fmt(parse(Float64, rows[index][11]) + 0.1)
        write_tsv_rows(joinpath(fixture.t11.output, "scores.tsv"), header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_date_count", "scores",
             "change date score count", "artifact/graph", "date aggregation", :BLOCKED, :row_count) do fixture
        header, rows = tsv_rows(joinpath(fixture.t11.output, "scores.tsv"))
        index = findfirst(row -> row[6] == "date", rows)
        rows[index][12] = "1"
        write_tsv_rows(joinpath(fixture.t11.output, "scores.tsv"), header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_date_mean", "scores",
             "change date mean", "artifact/graph", "date aggregation", :BLOCKED, :score_mismatch) do fixture
        header, rows = tsv_rows(joinpath(fixture.t11.output, "scores.tsv"))
        index = findfirst(row -> row[6] == "date", rows)
        rows[index][11] = fmt(parse(Float64, rows[index][11]) + 0.1)
        write_tsv_rows(joinpath(fixture.t11.output, "scores.tsv"), header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_raw_heldout_hash", "scores",
             "change partition raw-heldout digest", "artifact/graph", "raw-heldout score binding", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "partitions.tsv"), 1) do row
            row[14] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_scan_identity", "scores",
             "change scan score file identity", "artifact/graph", "scan score identity", :BLOCKED, :reference_mismatch) do fixture
        header, rows = tsv_rows(joinpath(fixture.t11.output, "scores.tsv"))
        index = findfirst(row -> row[6] == "scan", rows)
        rows[index][7] = "20260199_scan01.sxm"
        write_tsv_rows(joinpath(fixture.t11.output, "scores.tsv"), header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_scan_count", "scores",
             "change scan score edge count", "artifact/graph", "scan score count", :BLOCKED, :row_count) do fixture
        header, rows = tsv_rows(joinpath(fixture.t11.output, "scores.tsv"))
        index = findfirst(row -> row[6] == "scan", rows)
        rows[index][12] = "1"
        write_tsv_rows(joinpath(fixture.t11.output, "scores.tsv"), header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_sparse_zero_coverage", "scores",
             "set one edge score count to zero with zero gain", "artifact/graph", "sparse score formula", :BLOCKED, :invalid_key) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "scores.tsv"), 1) do row
            row[11] = "0"; row[12] = "0"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("score_hash", "scores",
             "change published score digest", "artifact/graph", "score hash", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "scores.tsv"), 1) do row
            row[13] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_control_mean", "shuffle",
             "change one observed permutable control mean", "artifact/graph", "shuffle control", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[9] = fmt(parse(Float64, row[9]) + 0.1)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_scope", "shuffle",
             "change one shuffle scope", "artifact/graph", "shuffle reference", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[3] = "wrong_scope"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_seed_range", "shuffle",
             "set seed to 500", "artifact/graph", "shuffle seed", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[6] = "500"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_seed_order", "shuffle",
             "duplicate seed 1", "artifact/graph", "shuffle seed order", :BLOCKED, :order_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[6] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_run_state", "shuffle",
             "set one running date gain to NA", "artifact/graph", "shuffle run-state", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[7] = "NA"
            row[8] = "not_run"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_nonrun_field_combination", "shuffle",
             "put date gain on a non-running outer-score row", "artifact/graph", "non-run shuffle schema", :BLOCKED, :invalid_hash) do fixture
        header, rows = tsv_rows(joinpath(fixture.t11.output, "shuffle.tsv"))
        index = findfirst(row -> row[3] == "outer_score", rows)
        rows[index][7] = "0"
        write_tsv_rows(joinpath(fixture.t11.output, "shuffle.tsv"), header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_running_value", "shuffle",
             "make running date gain nonfinite", "artifact/graph", "running shuffle value", :BLOCKED, :invalid_number) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[7] = "NaN"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_quantile", "shuffle",
             "change Type-7 upper quantile", "artifact/graph", "shuffle quantile", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[10] = fmt(parse(Float64, row[10]) + 0.1)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_pass", "shuffle",
             "flip shuffle pass", "artifact/graph", "shuffle predicate", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[11] = "false"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_score_binding", "shuffle",
             "change shuffle score binding", "artifact/graph", "shuffle score hash", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[12] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("shuffle_same_fit_seed_alias", "shuffle",
             "change one conditional-fit digest for a shared base-fit/seed", "artifact/graph", "shuffle fit alias", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "shuffle.tsv"), 1) do row
            row[8] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_missing_seed", "bootstrap",
             "remove one bootstrap row", "artifact/graph", "bootstrap seed coverage", :BLOCKED, :order_mismatch) do fixture
        path = joinpath(fixture.t11.output, "bootstrap.tsv")
        header, rows = tsv_rows(path); deleteat!(rows, 1)
        write_tsv_rows(path, header, rows); refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_duplicate_seed", "bootstrap",
             "duplicate one bootstrap row", "artifact/graph", "bootstrap duplicate seed", :BLOCKED, :order_mismatch) do fixture
        path = joinpath(fixture.t11.output, "bootstrap.tsv")
        header, rows = tsv_rows(path); insert!(rows, 2, copy(rows[1]))
        write_tsv_rows(path, header, rows); refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_seed_order", "bootstrap",
             "swap first two bootstrap seeds", "artifact/graph", "bootstrap seed order", :BLOCKED, :order_mismatch) do fixture
        path = joinpath(fixture.t11.output, "bootstrap.tsv")
        header, rows = tsv_rows(path); rows[1][5], rows[2][5] = rows[2][5], rows[1][5]
        write_tsv_rows(path, header, rows); refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_value", "bootstrap",
             "change bootstrap mean", "artifact/graph", "bootstrap value", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[6] = fmt(parse(Float64, row[6]) + 0.1)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_lower_quantile", "bootstrap",
             "change bootstrap lower quantile", "artifact/graph", "bootstrap Type-7 lower", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[7] = fmt(parse(Float64, row[7]) + 0.1)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_predicates", "bootstrap",
             "flip every-date-positive predicate", "artifact/graph", "bootstrap predicates", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[8] = "false"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_gate_hash", "bootstrap",
             "change bootstrap gate hash", "artifact/graph", "bootstrap gate hash", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[13] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_partition_date_aggregation", "bootstrap",
             "change one gate partition/date identity", "artifact/graph", "bootstrap gate aggregation", :BLOCKED, :row_count) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[4] = "20260199"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_every_shuffle", "bootstrap",
             "flip every-shuffle predicate", "artifact/graph", "bootstrap shuffle predicate", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[9] = "false"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_reversal", "bootstrap",
             "flip reversal predicate", "artifact/graph", "bootstrap reversal predicate", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[10] = "false"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_status", "bootstrap",
             "change bootstrap terminal status", "artifact/graph", "bootstrap status", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[11] = "SKIPPED"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("bootstrap_reason", "bootstrap",
             "change bootstrap terminal reason", "artifact/graph", "bootstrap reason", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "bootstrap.tsv"), 1) do row
            row[12] = "tampered"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("model_outer_count", "model_result",
             "change outer fold count", "artifact/graph", "model outer tuple", :BLOCKED, :row_count) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "models.tsv"), 1) do row
            row[5] = "6"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("model_final_gate", "model_result",
             "change final gate binding", "artifact/graph", "model final gate", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "models.tsv"), 1) do row
            row[6] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("model_full_fit", "model_result",
             "change full-fit binding", "artifact/graph", "model full fit", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "models.tsv"), 1) do row
            row[7] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("model_distinct_full_refit", "model_result",
             "swap one full-refit and six-date fit identity coherently", "artifact/graph/result", "distinct-full-refit validator", :BLOCKED, :reference_mismatch) do fixture
        model_header, model_rows = tsv_rows(fixture.t11.models)
        fit_col = findfirst(==( "fit_sha256"), SCI._GRAPH_MODEL_HEADER)
        role_col = findfirst(==( "fit_role"), SCI._GRAPH_MODEL_HEADER)
        model_col = findfirst(==( "model"), SCI._GRAPH_MODEL_HEADER)
        full_index = only(findall(row -> row[model_col] == "C1" && row[role_col] == "full_refit", model_rows))
        full_fit = model_rows[full_index][fit_col]
        outer_fit = first(split(reference, '|')[6]
                          for row in model_rows if row[model_col] == "C1"
                          for reference in split(row[16], ';')
                          if split(reference, '|')[3] == "outer_score")
        for index in eachindex(model_rows)
            model_rows[index][fit_col] = model_rows[index][fit_col] == full_fit ? outer_fit :
                model_rows[index][fit_col] == outer_fit ? full_fit : model_rows[index][fit_col]
            references = model_rows[index][16]
            references = replace(references, full_fit => "__FULL_REFIT_SWAP__")
            references = replace(references, outer_fit => full_fit)
            model_rows[index][16] = replace(references, "__FULL_REFIT_SWAP__" => outer_fit)
        end
        write_tsv_rows(fixture.t11.models, model_header, model_rows)
        for (name, fit_column) in (("fitted_unary_nodes.tsv", 2),
                                   ("fitted_edge_transforms.tsv", 2),
                                   ("partitions.tsv", 12),
                                   ("ess.tsv", 16),
                                   ("starts.tsv", 17))
            path = joinpath(fixture.t11.output, name)
            header, rows = tsv_rows(path)
            for row in rows
                row[fit_column] = row[fit_column] == full_fit ? outer_fit :
                    row[fit_column] == outer_fit ? full_fit : row[fit_column]
            end
            write_tsv_rows(path, header, rows)
        end
        result_header, result_rows = tsv_rows(joinpath(fixture.t11.output, "models.tsv"))
        result_rows[1][7] = outer_fit
        write_tsv_rows(joinpath(fixture.t11.output, "models.tsv"), result_header, result_rows)
        refresh_t11_wrappers!(fixture; graph_lists=true)
        receipt = TOML.parse(read(fixture.t11.receipt, String))
        _, result_values = tsv_rows(joinpath(fixture.t11.output, "models.tsv"))
        result_row = only(filter(row -> row[2] == "C1", result_values))
        _, graph_values = tsv_rows(fixture.t11.models)
        model_row = only(filter(row -> row[model_col] == "C1" && row[role_col] == "full_refit", graph_values))
        references = vcat([split(row[16], ';') for row in graph_values if row[model_col] == "C1"]...)
        outer_dates = sort(unique(split(ref, '|')[4] for ref in references
                                  if split(ref, '|')[3] == "outer_inner"))
        _, bootstrap_values = tsv_rows(joinpath(fixture.t11.output, "bootstrap.tsv"))
        _, partition_values = tsv_rows(joinpath(fixture.t11.output, "partitions.tsv"))
        outer_lines = String[]
        for date in outer_dates
            inner_gate = first(row[13] for row in bootstrap_values
                               if row[2] == "C1" && row[3] == "outer_inner" && row[4] == date)
            outer_reference = only(split(ref, '|') for ref in references
                                   if split(ref, '|')[3] == "outer_score" && split(ref, '|')[4] == date)
            outer_score = only(row for row in partition_values
                               if row[2] == "C1" && row[3] == "outer_score" && row[4] == date)
            push!(outer_lines, join(("outer", date, inner_gate, outer_score[13], outer_reference[6]), '\t'))
        end
        result_hash = hash_lines([
            "schema=structured-edge-model-admission-v1", "model=C1",
            "status=$(result_row[3])", "reason=$(result_row[4])",
            "final_gate=$(result_row[6])", "full_fit=$(result_row[7])",
            outer_lines...,
        ])
        result_row[8] = result_hash
        write_tsv_rows(joinpath(fixture.t11.output, "models.tsv"), result_header, [result_row])
        receipt["models"]["C1"] = "$(result_row[3])|$(result_row[4])|$(result_hash)"
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
        refresh_t11_input_bindings!(fixture)
        refresh_t11_wrappers!(fixture; graph_lists=true)
    end
    add_case("model_status_precedence", "model_result",
             "change overall model status", "artifact/graph", "model status precedence", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "models.tsv"), 1) do row
            row[3] = "SKIPPED"
        end
        receipt = TOML.parse(read(fixture.t11.receipt, String))
        binding = split(String(receipt["models"]["C1"]), '|'; keepempty=true)
        binding[1] = "SKIPPED"
        receipt["models"]["C1"] = join(binding, '|')
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
        refresh_t11_wrappers!(fixture)
    end
    add_case("model_reason_precedence", "model_result",
             "change overall model reason with receipt rebind", "artifact/graph/result", "model reason precedence", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "models.tsv"), 1) do row
            row[4] = "tampered_reason"
        end
        receipt = TOML.parse(read(fixture.t11.receipt, String))
        binding = split(String(receipt["models"]["C1"]), '|'; keepempty=true)
        binding[2] = "tampered_reason"
        receipt["models"]["C1"] = join(binding, '|')
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
        refresh_t11_wrappers!(fixture)
    end
    add_case("model_result_hash", "model_result",
             "change model result hash", "artifact/graph", "model result hash", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(joinpath(fixture.t11.output, "models.tsv"), 1) do row
            row[8] = repeat("0", 64)
        end
        receipt = TOML.parse(read(fixture.t11.receipt, String))
        binding = split(String(receipt["models"]["C1"]), '|'; keepempty=true)
        binding[3] = repeat("0", 64)
        receipt["models"]["C1"] = join(binding, '|')
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
        refresh_t11_wrappers!(fixture)
    end
    add_case("receipt_status", "receipt",
             "change receipt status", "none", "receipt status", :BLOCKED, :status_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["status"] = "SKIPPED"
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_historical_t11_source", "receipt",
             "replace effective source bundle with immutable correction2 bundle",
             "none", "receipt source authority", :BLOCKED, :hash_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String))
        receipt["hashes"]["source_sha256"] = SCI._PRIOR_SCALE_T11_SOURCE_BUNDLE_SHA256
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_reason", "receipt",
             "change receipt reason", "none", "receipt reason", :BLOCKED, :status_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["reason"] = "tampered"
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_result_hash", "receipt",
             "change receipt result hash", "none", "receipt result", :BLOCKED, :hash_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["result_sha256"] = repeat("0", 64)
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_input_authority", "receipt",
             "change universe input authority", "none", "receipt input hash", :BLOCKED, :hash_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["hashes"]["universe_sha256"] = repeat("0", 64)
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_input_hash", "receipt",
             "change Todo 10 input authority", "none", "receipt input authority", :BLOCKED, :hash_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["hashes"]["input_sha256"] = repeat("0", 64)
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_normalization", "receipt",
             "change normalization authority", "none", "receipt normalization authority", :BLOCKED, :hash_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["hashes"]["normalization_sha256"] = repeat("0", 64)
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_artifact_hash", "receipt",
             "change one artifact hash", "none", "receipt artifact hash", :BLOCKED, :hash_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["artifacts"]["scores.tsv"] = repeat("0", 64)
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_artifact_count", "receipt",
             "change one artifact row count", "none", "receipt artifact count", :BLOCKED, :row_count) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["row_counts"]["scores.tsv"] -= 1
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_artifact_map", "receipt",
             "remove one artifact map entry", "none", "receipt artifact map", :BLOCKED, :schema_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); delete!(receipt["artifacts"], "scores.tsv")
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end
    add_case("receipt_model_map", "receipt",
             "change C1 model-map binding", "none", "receipt model map", :BLOCKED, :reference_mismatch) do fixture
        receipt = TOML.parse(read(fixture.t11.receipt, String)); receipt["models"]["C1"] = "PASS|ok|$(repeat("0",64))"
        write(fixture.t11.receipt, sprint(TOML.print, receipt))
    end

    add_case("density_injected_failure", "public_failures",
             "inject density failure at factor construction", "private hook", "density stage", :FAIL, :density_failure) do fixture
        SCI._DENSITY_FAILURE_HOOK[] = () -> nothing
    end
    add_case("dp_injected_failure", "public_failures",
             "inject DP failure at first DP block", "private hook", "DP stage", :FAIL, :dp_failure) do fixture
        SCI._DP_FAILURE_HOOK[] = () -> nothing
    end
    add_case("unexpected_internal_failure", "public_failures",
             "throw ordinary exception from post-validation hook", "private existing hook", "internal inference", :FAIL, :internal_inference) do fixture
        SCI._POST_VALIDATION_HOOK[] = () -> error("correction7 injected unexpected exception")
    end
    add_case("reversal_failure", "public_failures",
             "perturb reversal result after validation", "private existing hook", "reversal check", :FAIL, :reversal_mismatch) do fixture
        SCI._REVERSAL_CHECK_HOOK[] = (natural, reversed) -> (reversed[3][1, 1] += 0.1)
    end

    for specification in cases
        fixture = copy_producer_fixture(canonical)
        SCI._POST_VALIDATION_HOOK[] = nothing
        SCI._SNAPSHOT_RETURN_HOOK[] = nothing
        try
            specification.action(fixture)
            report = specification.invoke(fixture)
            @test report.status == specification.status
            @test report.reason == specification.reason
            push!(matrix_rows, matrix_result_row(specification.name, specification.family,
                "seven-date-producer", specification.mutation,
                specification.refreshed, specification.layer, specification.status,
                specification.reason, report))
        finally
            SCI._POST_VALIDATION_HOOK[] = nothing
            SCI._SNAPSHOT_RETURN_HOOK[] = nothing
            SCI._DENSITY_FAILURE_HOOK[] = nothing
            SCI._DP_FAILURE_HOOK[] = nothing
            rm(fixture.root; recursive=true, force=true)
        end
    end
    skipped_fixture = make_producer_fixture(
        dates=["20260101", "20260102", "20260103", "20260104"],
        require_pass=false,
    )
    try
        skipped_report = infer(skipped_fixture)
        @test skipped_report.status == :SKIPPED
        @test skipped_report.reason == :unary_only
        push!(matrix_rows, matrix_result_row(
            "k_lt_7_producer_skipped", "model_graph", "four-date-producer-skipped",
            "evaluate producer early-SKIPPED report", "producer report_files",
            "unary-only public branch", :SKIPPED, :unary_only, skipped_report))
        impossible = copy_producer_fixture(skipped_fixture)
        try
            mutate_tsv_rows!(joinpath(impossible.t11.output, "models.tsv")) do rows
                for row in rows
                    row[3] = "PASS"; row[4] = "ok"
                end
            end
            receipt = TOML.parse(read(impossible.t11.receipt, String))
            for model_id in ("C1", "C2")
                binding = split(String(receipt["models"][model_id]), '|'; keepempty=true)
                binding[1] = "PASS"; binding[2] = "ok"
                receipt["models"][model_id] = join(binding, '|')
            end
            write(impossible.t11.receipt, sprint(TOML.print, receipt))
            refresh_t11_wrappers!(impossible)
            impossible_report = infer(impossible)
            @test impossible_report.status == :BLOCKED
            @test impossible_report.reason == :status_mismatch
            push!(matrix_rows, matrix_result_row(
                "k_lt_7_impossible_pass", "model_graph", "four-date-producer-skipped",
                "coherently artifact-rebound PASS overall rows", "artifact/graph",
                "early branch status", :BLOCKED, :status_mismatch, impossible_report))
        finally
            rm(impossible.root; recursive=true, force=true)
        end
    finally
        rm(skipped_fixture.root; recursive=true, force=true)
    end
    rm(canonical.root; recursive=true, force=true)
    @test length(matrix_rows) == length(cases) + 2
    evidence_dir = SCALE_REBIND_CORRECTION_EVIDENCE_DIR
    mkpath(evidence_dir)
    matrix_header = ["name", "family", "fixture", "mutation", "refreshed_bindings",
                     "intended_layer", "expected_status", "expected_reason",
                     "actual_status", "actual_reason", "scan_count", "factor_count", "dp_count"]
    write(joinpath(evidence_dir, "correction7-executable-matrix.tsv"),
          join(matrix_header, '\t') * "\n" *
           join([join(row, '\t') for row in matrix_rows], "\n") * "\n")
    write(joinpath(evidence_dir, "correction7-retained-assertions.log"), join([
        "not-counted-as-deep-mutations",
        "public path-only graph inference: PASS replay and factor/cache-key assertions",
        "validated reversal mismatch is an inference failure",
        "private binary recurrence equals exhaustive enumeration",
        "strict provenance blocks before DP",
        "static and forbidden-boundary audits",
        "public unary -Inf and zero-factor semantics",
        "out-of-range gaps split topology and preserve unary blocks",
        "seven-path API rejection of extra matrices/weights and label-free source checks",
        "shuffle conditional-fit rule: different seeds may differ, and conditional-fit digests may also coincide",
        "t11-exclusive-opaque-preimages: producer serializer reports, nested partitions, gate/bootstrap traces, and reversal records",
        "t11-exclusive-opaque-preimages are covered by test_structured_edge_admission.jl and are not counted in the Todo12 deep matrix",
    ], "\n") * "\n")
end

@testset "correction8 semantic replacements and admissions" begin
    correction7_matrix = joinpath(ROOT, ".omo", "evidence",
        "structured-label-free-unit-assignment", "t12", "correction7",
        "executable-matrix.tsv")
    retained_header, retained_values = tsv_rows(correction7_matrix)
    rejected = Set(["score_sparse_zero_coverage", "shuffle_nonrun_field_combination"])
    retained_rows = [row for row in retained_values if row[1] ∉ rejected]
    @test length(retained_rows) == 126
    @test length(unique(row[1] for row in retained_rows)) == length(retained_rows)

    canonical = make_producer_fixture()
    stage_sources = Dict(
        kind => make_producer_stage_fixture(kind)
        for kind in (:zero_support, :insufficient_support, :null_only, :conditional_failure)
    )
    cases = NamedTuple[]
    add_case(name, family, source, mutation, refreshed, layer, status, reason, action;
             invoke=infer) = push!(cases, (; name, family, source, mutation, refreshed,
                                           layer, status, reason, action, invoke))
    add_case(action::Function, name, family, source, mutation, refreshed, layer, status, reason;
             invoke=infer) = add_case(name, family, source, mutation, refreshed, layer,
                                       status, reason, action; invoke=invoke)

    add_case("r1_sparse_terminal_positive", "sparse_scores", :zero_support,
             "serializer-produced sparse terminal publication", "report_files",
             "sparse terminal stage", :SKIPPED, :unary_only, f -> nothing)
    add_case("r1_sparse_numeric_gain_zero_count", "sparse_scores", :zero_support,
             "set sparse date gain numeric while edge count remains zero", "artifact/graph",
             "sparse score semantics", :BLOCKED, :schema_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "scores.tsv")) do rows
            row = first(filter(row -> row[3] == "outer_score" && row[6] == "date", rows))
            row[11] = "0"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r1_sparse_non_na_identity_zero_count", "sparse_scores", :zero_support,
             "set sparse date identity non-NA while edge count remains zero", "artifact/graph",
             "sparse score semantics", :BLOCKED, :row_count) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "scores.tsv")) do rows
            row = first(filter(row -> row[3] == "outer_score" && row[6] == "date", rows))
            row[7] = "20260101_scan01.sxm"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r1_sparse_nonzero_count", "sparse_scores", :zero_support,
             "set sparse date edge count to one while status remains SKIPPED", "artifact/graph",
             "sparse score semantics", :BLOCKED, :row_count) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "scores.tsv")) do rows
            row = first(filter(row -> row[3] == "outer_score" && row[6] == "date", rows))
            row[12] = "1"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r1_sparse_invalid_score_formula", "sparse_scores", :zero_support,
             "rebind partition and score hashes to a wrong sparse invalid-score digest", "artifact/graph/result",
             "sparse score formula", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "partitions.tsv")) do rows
            row = first(filter(row -> row[2] == "C1" && row[3] == "outer_score", rows))
            row[13] = repeat("0", 64)
        end
        mutate_tsv_rows!(joinpath(fixture.t11.output, "scores.tsv")) do rows
            row = first(filter(row -> row[2] == "C1" && row[3] == "outer_score" && row[6] == "date", rows))
            row[13] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r1_residualizer_insufficient_positive", "staged_fits", :insufficient_support,
             "serializer-produced residualizer-present insufficient-support stage", "report_files",
             "insufficient-support stage", :SKIPPED, :unary_only, f -> nothing)
    add_case("r1_insufficient_support_ess_mutation", "staged_fits", :insufficient_support,
             "change ESS support edge count in a residualizer-present stage", "artifact/graph",
             "ESS support semantics", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "ess.tsv")) do rows
            row = first(filter(row -> row[2] == "C1" && row[3] == "outer_score", rows))
            row[9] = "999"
        end
        refresh_t11_wrappers!(fixture)
    end

    add_case("r2_nonrun_control", "shuffle_nonrun", :canonical,
             "change finite control on an actual outer-score non-run row", "artifact/graph",
             "non-run control semantics", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "shuffle.tsv")) do rows
            row = first(filter(row -> row[3] == "outer_score" && row[7] == "NA" && row[8] == "not_run", rows))
            row[9] = fmt(parse(Float64, row[9]) + 0.1)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r2_nonrun_pass", "shuffle_nonrun", :canonical,
             "flip passed on an actual outer-score non-run row", "artifact/graph",
             "non-run predicate semantics", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "shuffle.tsv")) do rows
            row = first(filter(row -> row[3] == "outer_score" && row[7] == "NA" && row[8] == "not_run", rows))
            row[11] = "false"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r2_nonrun_quantile", "shuffle_nonrun", :canonical,
             "publish a numeric upper quantile on an actual non-run row", "artifact/graph",
             "non-run schema semantics", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "shuffle.tsv")) do rows
            row = first(filter(row -> row[3] == "outer_score" && row[7] == "NA" && row[8] == "not_run", rows))
            row[10] = "0"
        end
        refresh_t11_wrappers!(fixture)
    end

    for (name, column_start, values, description) in [
        ("r3_null_subfloor", 27, [fmt(0.5), fmt(0.49999), fmt(0.5)], "null minimum eigenvalue 1e-5"),
        ("r3_conditional_subfloor", 38, [fmt(0.5), fmt(0.49999), fmt(0.5)], "conditional minimum eigenvalue 1e-5"),
        ("r3_null_overcap", 27, [fmt(1.00005), fmt(0.99995), fmt(1.00005)], "null condition about 2e4"),
        ("r3_conditional_overcap", 38, [fmt(1.00005), fmt(0.99995), fmt(1.00005)], "conditional condition about 2e4"),
    ]
        add_case(name, "scale_feasibility", :canonical, description, "artifact/graph",
                 "scale feasibility", :BLOCKED, :model_mismatch) do fixture
            mutate_tsv_row!(fixture.t11.models, 1) do row
                row[column_start:column_start + 2] = values
            end
            refresh_t11_wrappers!(fixture)
        end
    end

    add_case("r4_fitted_model_missing", "fit_coverage", :canonical,
             "remove one fitted model row and refresh graph counts/lists", "artifact/graph",
             "fitted-model coverage", :BLOCKED, :row_count) do fixture
        header, rows = tsv_rows(fixture.t11.models)
        deleteat!(rows, 1)
        write_tsv_rows(fixture.t11.models, header, rows)
        refresh_t11_wrappers!(fixture; graph_lists=true)
    end
    add_case("r4_fitted_model_duplicate", "fit_coverage", :canonical,
             "duplicate one fitted model row and refresh graph counts/lists", "artifact/graph",
             "fitted-model duplicate", :BLOCKED, :duplicate_key) do fixture
        header, rows = tsv_rows(fixture.t11.models)
        insert!(rows, 2, copy(rows[1]))
        write_tsv_rows(fixture.t11.models, header, rows)
        refresh_t11_wrappers!(fixture; graph_lists=true)
    end
    add_case("r4_transform_residualizer", "transform_bindings", :canonical,
             "change transform residualizer to another canonical 64-hex digest", "artifact/graph",
             "transform residualizer binding", :BLOCKED, :hash_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.transforms, 1) do row
            row[17] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r4_transform_fit_reference", "transform_bindings", :canonical,
             "retarget transform fit to an unknown 64-hex digest", "artifact/graph",
             "transform fit-reference binding", :BLOCKED, :reference_mismatch) do fixture
        mutate_tsv_row!(fixture.t11.transforms, 1) do row
            row[2] = repeat("0", 64)
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r4_zero_support_positive", "staged_fits", :zero_support,
             "serializer-produced zero-support terminal stage", "report_files",
             "zero-support stage", :SKIPPED, :unary_only, f -> nothing)
    add_case("r4_zero_support_status_mutation", "staged_fits", :zero_support,
             "change zero-support ESS reason", "artifact/graph",
             "zero-support ESS/status", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "ess.tsv")) do rows
            row = first(filter(row -> row[2] == "C1" && row[3] == "outer_score", rows))
            row[15] = "ok"
        end
        refresh_t11_wrappers!(fixture)
    end
    add_case("r4_null_only_positive", "staged_fits", :null_only,
             "serializer-produced null-only failure stage", "report_files",
             "null-only stage", :FAIL, :todo11_fail, f -> nothing)
    add_case("r4_null_only_coverage", "staged_fits", :null_only,
             "remove null start from a serializer-produced null-only stage", "artifact/graph",
             "null-only start coverage", :BLOCKED, :row_count) do fixture
        path = joinpath(fixture.t11.output, "starts.tsv")
        header, rows = tsv_rows(path)
        deleteat!(rows, findfirst(row -> row[2] == "C1" && row[3] == "outer_score" && row[6] == "null", rows))
        write_tsv_rows(path, header, rows)
        refresh_t11_wrappers!(fixture)
    end
    add_case("r4_conditional_failure_positive", "staged_fits", :conditional_failure,
             "serializer-produced failed-conditional stage", "report_files",
             "conditional-failure stage", :FAIL, :todo11_fail, f -> nothing)
    add_case("r4_conditional_failure_status", "staged_fits", :conditional_failure,
             "change a failed conditional start status", "artifact/graph",
             "conditional-failure start semantics", :BLOCKED, :status_mismatch) do fixture
        mutate_tsv_rows!(joinpath(fixture.t11.output, "starts.tsv")) do rows
            row = first(filter(row -> row[2] == "C1" && row[3] == "outer_score" && row[6] == "conditional", rows))
            row[9] = "PASS"
        end
        refresh_t11_wrappers!(fixture)
    end

    for specification in cases
        source = specification.source == :canonical ? canonical : stage_sources[specification.source]
        fixture = copy_producer_fixture(source)
        SCI._POST_VALIDATION_HOOK[] = nothing
        SCI._SNAPSHOT_RETURN_HOOK[] = nothing
        try
            specification.action(fixture)
            report = specification.invoke(fixture)
            @test report.status == specification.status
            @test report.reason == specification.reason
            push!(retained_rows, matrix_result_row(specification.name, specification.family,
                String(specification.source), specification.mutation,
                specification.refreshed, specification.layer, specification.status,
                specification.reason, report))
        finally
            SCI._POST_VALIDATION_HOOK[] = nothing
            SCI._SNAPSHOT_RETURN_HOOK[] = nothing
            SCI._DENSITY_FAILURE_HOOK[] = nothing
            SCI._DP_FAILURE_HOOK[] = nothing
            rm(fixture.root; recursive=true, force=true)
        end
    end
    @test length(unique(row[1] for row in retained_rows)) == length(retained_rows)
    @test length(retained_rows) == 126 + length(cases)
    @test all(length(row) == 13 for row in retained_rows)

    matrix_header = ["name", "family", "fixture", "mutation", "refreshed_bindings",
                     "intended_layer", "expected_status", "expected_reason",
                     "actual_status", "actual_reason", "scan_count", "factor_count", "dp_count"]
    correction8_matrix_path = joinpath(ROOT, ".omo", "evidence",
        "structured-label-free-unit-assignment", "t12", "correction8",
        "executable-matrix.tsv")
    @test isfile(correction8_matrix_path)
    rm(canonical.root; recursive=true, force=true)
    for source in values(stage_sources)
        rm(source.root; recursive=true, force=true)
    end
end

@testset "correction8 exact scale boundaries" begin
    fixture = make_producer_fixture()
    try
        _, rows = tsv_rows(fixture.t11.models)
        base = first(row for row in rows if row[4] == "C1" && row[2] == "full_refit")
        columns = Dict(name => index for (index, name) in enumerate(SCI._GRAPH_MODEL_HEADER))
        for (start, values) in ((columns["null_scale_ff"], [fmt(1e-4), "0", fmt(1e-4)]),
                                (columns["conditional_scale_ff"], [fmt(1e-4), "0", fmt(1e-4)]),
                                (columns["null_scale_ff"], [fmt(1e-4), "0", "1"]),
                                (columns["conditional_scale_ff"], [fmt(1e-4), "0", "1"]))
            row = String.(copy(base))
            row[start:start + 2] = values
            parsed = SCI._parse_model(row, columns)
            @test parsed.status == "PASS"
        end
    finally
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "correction exact T11/T12 scale eigensolver parity" begin
    counterexample = [
        0.32262332547297057 -0.46740750707094925
        -0.46740750707094925 0.6774766745270292
    ]
    decomposition = eigen(Symmetric(counterexample))
    @test decomposition.values == [9.999999999998899e-5, 0.9999999999999998]
    @test maximum(decomposition.values) / minimum(decomposition.values) == 10000.000000001099
    @test SEA._stored_scale_certificate(counterexample).status == :floor
    disposition = try
        SCI._validate_scale_feasibility(counterexample, "exact non-diagonal floor")
        :ok
    catch error
        error isa SCI._ChainBlocked || rethrow()
        error.code
    end
    @test disposition == :model_mismatch

    parity_cases = [
        ("interior", [0.25 0.0; 0.0 0.5], :ok),
        ("exact diagonal boundary", [1.0e-4 0.0; 0.0 1.0e-4], :ok),
        ("rotated certified boundary",
         [0.25004999999999994 -0.43292609935184084;
          -0.43292609935184084 0.7499500000000001], :ok),
        ("non-diagonal floor failure", counterexample, :floor),
        ("non-diagonal cap failure",
         [0.25007500002499994 -0.43296940066533135;
          -0.43296940066533135 0.7500250000750002], :condition),
    ]
    for (name, matrix, expected) in parity_cases
        authority = SEA._stored_scale_certificate(matrix)
        todo12 = try
            SCI._validate_scale_feasibility(Matrix{Float64}(matrix), name)
            :ok
        catch error
            error isa SCI._ChainBlocked || rethrow()
            error.code == :model_mismatch || rethrow()
            :model_mismatch
        end
        @test authority.status == expected
        @test todo12 == (expected == :ok ? :ok : :model_mismatch)
        @test (authority.status == :ok) == (todo12 == :ok)
    end
end

@testset "correction9 exact frozen scale boundaries" begin
    source_text = read(joinpath(ROOT, "test", "lib", "structured_assignment", "edge_model.jl"), String)
    @test !occursin("SCALE_BOUNDARY", source_text)
    @test !occursin("ROUND_OFF", source_text)

    correction8_matrix = joinpath(ROOT, ".omo", "evidence",
        "structured-label-free-unit-assignment", "t12", "correction8",
        "executable-matrix.tsv")
    matrix_header, retained_rows = tsv_rows(correction8_matrix)
    @test length(retained_rows) == 150
    canonical = make_producer_fixture()
    _, boundary_base_rows = tsv_rows(canonical.t11.models)
    boundary_base = String.(first(row for row in boundary_base_rows
                                   if row[4] == "C1" && row[2] == "full_refit"))
    boundary_columns = Dict(name => index for (index, name) in enumerate(SCI._GRAPH_MODEL_HEADER))
    for (column_start, values) in (
        (boundary_columns["null_scale_ff"], [fmt(1.0e-4), "0", fmt(1.0e-4)]),
        (boundary_columns["conditional_scale_ff"], [fmt(1.0e-4), "0", fmt(1.0e-4)]),
        (boundary_columns["null_scale_ff"], [fmt(1.0e-4), "0", "1"]),
        (boundary_columns["conditional_scale_ff"], [fmt(1.0e-4), "0", "1"]),
    )
        row = copy(boundary_base)
        row[column_start:column_start + 2] = values
        @test SCI._parse_model(row, boundary_columns).status == "PASS"
    end
    cases = [
        ("c9_null_prevfloat_floor", 27,
         [fmt(prevfloat(1.0e-4)), "0", fmt(prevfloat(1.0e-4))],
         "null diagonal minimum prevfloat(1e-4)"),
        ("c9_conditional_prevfloat_floor", 38,
         [fmt(prevfloat(1.0e-4)), "0", fmt(prevfloat(1.0e-4))],
         "conditional diagonal minimum prevfloat(1e-4)"),
        ("c9_null_nextfloat_condition", 27,
         [fmt(1.0e-4), "0", fmt(nextfloat(1.0))],
         "null diagonal maximum nextfloat(1.0)"),
        ("c9_conditional_nextfloat_condition", 38,
         [fmt(1.0e-4), "0", fmt(nextfloat(1.0))],
         "conditional diagonal maximum nextfloat(1.0)"),
        ("c9_non_diagonal_exact_floor_eigensolver_parity", 27,
         ["0.32262332547297057", "-0.46740750707094925", "0.6774766745270292"],
         "non-diagonal exact T11 floor rejection and T12 eigensolver parity"),
    ]
    rows = copy(retained_rows)
    boundary_header = "name\tkind\tvalue_1\tvalue_1_bits\tvalue_2\tvalue_2_bits\tvalue_3\tvalue_3_bits\tdisposition"
    boundary_rows = String[]
    for (name, column_start, values, description) in cases
        fixture = copy_producer_fixture(canonical)
        try
            mutate_tsv_row!(fixture.t11.models, 1) do row
                row[column_start:column_start + 2] = values
            end
            refresh_t11_wrappers!(fixture)
            report = infer(fixture)
            @test report.status == :BLOCKED
            @test report.reason == :model_mismatch
            @test SCI._SCAN_DECODE_CALL_COUNT[] == 0
            @test SCI._FACTOR_BUILD_CALL_COUNT[] == 0
            @test SCI._DP_CALL_COUNT[] == 0
            push!(rows, matrix_result_row(name, "scale_feasibility_adjacent", "seven-date-producer",
                description, "artifact/graph", "exact frozen scale feasibility",
                :BLOCKED, :model_mismatch, report))
            push!(boundary_rows, join((name, column_start == 27 ? "null" : "conditional",
                values[1], bitstring(parse(Float64, values[1])), values[2],
                bitstring(parse(Float64, values[2])), values[3],
                bitstring(parse(Float64, values[3])), "BLOCKED/model_mismatch/0/0/0"), '\t'))
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
    @test length(rows) == 155
    @test length(unique(row[1] for row in rows)) == 155
    evidence_dir = SCALE_REBIND_CORRECTION_EVIDENCE_DIR
    mkpath(evidence_dir)
    write(joinpath(evidence_dir, "correction9-executable-matrix.tsv"), matrix_header * "\n" *
          join([join(row, '\t') for row in rows], "\n") * "\n")
    write(joinpath(evidence_dir, "correction9-boundary-values.tsv"), boundary_header * "\n" *
          join(boundary_rows, "\n") * "\n")
    rm(canonical.root; recursive=true, force=true)
end

@testset "validated reversal mismatch is an inference failure" begin
    fixture = make_fixture()
    try
        SCI._REVERSAL_CHECK_HOOK[] = (natural, reversed) -> (reversed[3][1, 1] += 0.1)
        report = infer(fixture)
        @test report.status == :FAIL
        @test report.reason == :reversal_mismatch
    finally
        SCI._REVERSAL_CHECK_HOOK[] = nothing
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "private binary recurrence equals exhaustive enumeration" begin
    rng = MersenneTwister(12)
    for n in 1:10
        unaries = [(rand(rng), rand(rng)) for _ in 1:n]
        factors = [ntuple(_ -> rand(rng) - 0.5, 4) for _ in 1:(n - 1)]
        dp = SCI._dp_block(unaries, factors)
        exact = exhaustive(unaries, factors)
        @test dp[1] ≈ exact[1] atol=1e-12 rtol=1e-12
        @test dp[3] ≈ exact[2] atol=1e-12 rtol=1e-12
        @test dp[4] == exact[3]
        @test dp[5] ≈ exact[4] atol=1e-12 rtol=1e-12
    end
    tied = SCI._dp_block([(0.0, 0.0), (0.0, 0.0)], [(0.0, 0.0, 0.0, 0.0)])
    @test tied[4] == [0, 0]
    @test SCI._dp_block([(log(0.25), log(0.75))], NTuple{4,Float64}[])[3] ≈ [0.25 0.75] atol=1e-12 rtol=0
    @test SCI._dp_block([(0.0, -Inf)], NTuple{4,Float64}[])[3] == [1.0 0.0]
    @test SCI._dp_block([(-Inf, 0.0)], NTuple{4,Float64}[])[3] == [0.0 1.0]
    @test_throws SCI._InferenceFailed SCI._dp_block([(-Inf, -Inf)], NTuple{4,Float64}[])
    original = SCI._dp_block([(0.1, -0.2), (0.3, 0.4), (-0.1, 0.2)],
                             [(0.2, -0.1, 0.4, 0.0), (0.1, 0.5, -0.2, 0.3)])
    reversed = SCI._dp_block([( -0.1, 0.2), (0.3, 0.4), (0.1, -0.2)],
                             [(0.1, -0.2, 0.5, 0.3), (0.2, 0.4, -0.1, 0.0)])
    @test original[1] ≈ reversed[1] atol=1e-12 rtol=1e-12
    @test reversed[3][end:-1:1, :] ≈ original[3] atol=1e-12 rtol=1e-12
end

@testset "strict provenance blocks before DP" begin
    fixture = make_fixture()
    try
        baseline = infer(fixture)
        baseline_calls = SCI._DP_CALL_COUNT[]
        edge_bytes = read(fixture.t10.edge, String)
        write(fixture.t10.edge, replace(edge_bytes, "0.40000000000000002" => "0.41"; count=1))
        blocked = infer(fixture)
        @test blocked.status == :BLOCKED
        @test SCI._DP_CALL_COUNT[] == 0
        write(fixture.t10.edge, edge_bytes)
        receipt = read(fixture.t11.receipt, String)
        old_replay = match(r"graph_handoff_replay_sha256 = \"([0-9a-f]+)\"", receipt).captures[1]
        write(fixture.t11.receipt, replace(receipt, old_replay => repeat("0", 64); count=1))
        blocked_receipt = infer(fixture)
        @test blocked_receipt.status == :BLOCKED
        @test SCI._DP_CALL_COUNT[] == 0
        write(fixture.t11.receipt, receipt)
        real_edge = fixture.t10.edge * ".real"
        mv(fixture.t10.edge, real_edge)
        symlink(real_edge, fixture.t10.edge)
        blocked_symlink = infer(fixture)
        @test blocked_symlink.status == :BLOCKED
        @test SCI._DP_CALL_COUNT[] == 0
        rm(fixture.t10.edge; force=true)
        mv(real_edge, fixture.t10.edge)
        @test baseline.status == :PASS
        @test SCI._DP_CALL_COUNT[] == 0
    finally
        rm(fixture.root; recursive=true, force=true)
    end
end

@testset "static and forbidden-boundary audits" begin
    chain_source = read(joinpath(@__DIR__, "lib", "structured_assignment", "chain_inference.jl"), String)
    edge_source = read(joinpath(@__DIR__, "lib", "structured_assignment", "edge_model.jl"), String)
    @test occursin("module StructuredChainInference", chain_source)
    @test occursin("include(joinpath(@__DIR__, \"edge_model.jl\"))", chain_source)
    @test !occursin("type_posterior_inference.jl", chain_source * edge_source)
    @test !occursin("include(joinpath(@__DIR__, \"edge_admission.jl\"))", chain_source * edge_source)
    @test !occursin("expected_n", lowercase(chain_source * edge_source))
    @test !occursin("nknnkn", lowercase(chain_source * edge_source))
    @test !occursin("grader", lowercase(chain_source * edge_source))
    @test !occursin("benchmark", lowercase(chain_source * edge_source))
    @test !occursin("pair_prior::", lowercase(chain_source * edge_source))
end

@testset "T12 Todo10/T11 source authority identity correction" begin
    dynamic_source, source_snapshots = SEA.StructuredEdgeFeatures._source_snapshots()
    source_members = Dict(relpath(snapshot.path, ROOT) => snapshot.sha256
                          for snapshot in source_snapshots)
    expected_members = [
        "test/build_label_free_edge_features.jl",
        "test/lib/structured_assignment/edge_features.jl",
        "test/lib/structured_assignment/firewall.jl",
        "test/lib/structured_assignment/universe.jl",
    ]
    @test SCI._T10_SOURCE_SHA256 == dynamic_source
    @test dynamic_source == CANONICAL_T10_SOURCE_SHA
    @test sort(collect(keys(source_members))) == expected_members
    dynamic_lines = ["structured-edge-source-v1"]
    append!(dynamic_lines, "$(member)=$(source_members[member])"
            for member in sort(expected_members))
    @test hash_lines(dynamic_lines) == dynamic_source
    @test source_members["test/lib/structured_assignment/edge_features.jl"] ==
          OBSOLETE_T10_SOURCE_SHA
    @test dynamic_source != source_members["test/lib/structured_assignment/edge_features.jl"]

    canonical = make_fixture()
    try
        rebind_todo10_source!(canonical, CANONICAL_T10_SOURCE_SHA;
                              t11_source=CANONICAL_T10_SOURCE_SHA)
        report = infer(canonical)
        @test report.status == :PASS
        @test report.reason == :graph_inference
        error = todo10_validation_error(canonical)
        @test error === nothing
    finally
        rm(canonical.root; recursive=true, force=true)
    end

    obsolete = make_fixture()
    try
        set_todo10_receipt_source!(obsolete, OBSOLETE_T10_SOURCE_SHA)
        report = infer(obsolete)
        @test report.status == :BLOCKED
        @test report.reason == :hash_mismatch
        error = todo10_validation_error(obsolete)
        @test error !== nothing
        error !== nothing && @test error.message == "Todo 10 source binding differs"
    finally
        rm(obsolete.root; recursive=true, force=true)
    end

    forged_t11 = make_fixture()
    try
        before = TOML.parse(read(forged_t11.t11.receipt, String))
        set_todo11_edge_source!(forged_t11, OBSOLETE_T10_SOURCE_SHA)
        after = TOML.parse(read(forged_t11.t11.receipt, String))
        @test before["result_sha256"] == after["result_sha256"]
        @test before["graph_handoff_replay_sha256"] == after["graph_handoff_replay_sha256"]
        @test before["artifacts"] == after["artifacts"]
        report = infer(forged_t11)
        @test report.status == :BLOCKED
        @test report.reason == :hash_mismatch
        error = todo11_validation_error(forged_t11)
        @test error !== nothing
        error !== nothing && @test error.message == "Todo 10 source binding differs"
    finally
        rm(forged_t11.root; recursive=true, force=true)
    end
end

@testset "T12 distinguishes local admission input from Todo 10 features" begin
    local_input = make_fixture()
    try
        baseline = infer_with_paths(local_input)
        @test baseline.status == :PASS
        baseline_references = repr(baseline.references)
        baseline_log_evidence = baseline.log_evidence
        baseline_viterbi_score = baseline.viterbi_score

        receipt = TOML.parse(read(local_input.t11.receipt, String))
        todo10_feature = receipt["hashes"]["edge_feature_sha256"]
        local_feature = repeat("6", 64)
        @test local_feature != todo10_feature
        receipt["hashes"]["input_sha256"] = local_feature
        write(local_input.t11.receipt, sprint(TOML.print, receipt))
        refresh_t11_input_bindings!(local_input)

        distinct = infer_with_paths(local_input)
        @test distinct.status == :PASS
        @test distinct.reason == baseline.reason
        @test repr(distinct.references) == baseline_references
        @test distinct.log_evidence == baseline_log_evidence
        @test distinct.viterbi_score == baseline_viterbi_score
        @test distinct.result_sha256 != baseline.result_sha256
        @test distinct.provenance_sha256 != baseline.provenance_sha256
    finally
        rm(local_input.root; recursive=true, force=true)
    end

    forged_local = make_fixture()
    try
        receipt = TOML.parse(read(forged_local.t11.receipt, String))
        receipt["hashes"]["input_sha256"] = repeat("0", 64)
        write(forged_local.t11.receipt, sprint(TOML.print, receipt))
        report = infer_with_paths(forged_local)
        assert_blocked_zero(forged_local, report)
    finally
        rm(forged_local.root; recursive=true, force=true)
    end

    forged_edge_feature = make_fixture()
    try
        receipt = TOML.parse(read(forged_edge_feature.t11.receipt, String))
        receipt["hashes"]["edge_feature_sha256"] = repeat("7", 64)
        write(forged_edge_feature.t11.receipt, sprint(TOML.print, receipt))
        refresh_t11_input_bindings!(forged_edge_feature)
        report = infer_with_paths(forged_edge_feature)
        assert_blocked_zero(forged_edge_feature, report)
    finally
        rm(forged_edge_feature.root; recursive=true, force=true)
    end
end
