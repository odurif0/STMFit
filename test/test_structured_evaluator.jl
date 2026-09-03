#!/usr/bin/env julia

module StructuredEvaluatorTests

using Test
using TOML
using SHA
using Random
using Statistics

# Inclusion is deliberately the only evaluator entry here.  Its PROGRAM_FILE
# guard therefore keeps the production CLI dormant during this pure-seam run.
include(joinpath(@__DIR__, "evaluate_structured_unit_assignment.jl"))
const E = StructuredUnitAssignmentEvaluator
Test.TESTSET_PRINT_ENABLE[] = false

# Reuse the authoritative Todo 10 synthetic fixture builder without changing
# the frozen evaluator source.  The prefix ends immediately after
# make_cli_fixture and before the edge-admission test products.
const EDGE_FIXTURE_SUPPORT = let
    support = Module(:T13EdgeFixtureSupport)
    Core.eval(support, :(include(path::AbstractString) = Base.include(@__MODULE__, path)))
    source = joinpath(@__DIR__, "test_structured_edge_admission.jl")
    lines = readlines(source)
    prefix = join(lines[1:1065], "\n")
    Base.include_string(support, prefix, source)
    under_source = join(lines[835:1064], "\n")
    under_source = replace(under_source,
                           "function make_cli_fixture(;" =>
                           "function make_cli_fixture_under(source_root::String;",
                           count=1)
    under_source = replace(under_source,
                           "root = mktempdir(; prefix=\"stmfit-edge-admission-cli-\", cleanup=cleanup)" =>
                           """root = if get(ENV, \"STMFIT_T13_PUBLIC_REVALIDATION\", \"0\") == \"1\"
                               deterministic_root = joinpath(
                                   source_root, \".tmp-t13-public-fixture-runtime-v3\")
                               ispath(deterministic_root) &&
                                   error(\"deterministic public fixture already exists\")
                               mkdir(deterministic_root)
                               deterministic_root
                           else
                               mktempdir(source_root;
                                   prefix=\".tmp-t13-public-fixture-\", cleanup=cleanup)
                           end""",
                           count=1)
    Base.include_string(support, under_source, source)
    support
end

const VIEW_NAMES = (
    "base_local",
    "base_local+bwd_neg_com_t",
    "base_local+bwd_neg_diag45",
    "base_local+split_log_skew",
)
const ARTIFACTS = (
    "bootstrap.tsv", "dates.tsv", "events.tsv", "folds.tsv", "nodes.tsv",
    "receipt.toml", "scans.tsv", "signs.tsv", "summary.tsv",
)
const SCAN_DATES = ["20240101", "20240102", "20240103", "20240104", "20240105"]
const AUTHORITY = "a"^64
const UNIVERSE = "b"^64
const HISTORICAL_EVALUATOR_CONFIG = "config/unit_assignment_structured_evaluator.toml"
const RUNTIME_EVALUATOR_CONFIG = "config/unit_assignment_structured_evaluator_runtime.toml"
const RUNTIME_V2_EVALUATOR_CONFIG = "config/unit_assignment_structured_evaluator_runtime_v2.toml"
const RUNTIME_V3_EVALUATOR_CONFIG = "config/unit_assignment_structured_evaluator_runtime_v3.toml"

function view_set(c2_gap::Real; c1_gap::Real=0.0,
                 c2_missing::Bool=false, names=VIEW_NAMES)
    c1 = (0.0, Float64(c1_gap))
    c2 = c2_missing ? (NaN, NaN) : (0.0, Float64(c2_gap))
    return [E.ViewJointEvidence(name, c1, c2) for name in names]
end

function node_for(file::String, date::String, lobe::Int, c2_gap::Real;
                  c1_gap::Real=0.0, c2_missing::Bool=false,
                  names=VIEW_NAMES, t_nm::Real=lobe)
    return E.NodeEvidence(file, date, lobe, Float64(t_nm),
                          view_set(c2_gap; c1_gap, c2_missing, names))
end

function graph_for(model::String, c1_fit::String, c2_fit::String;
                   common::String="", enabled::Bool=false,
                   logz=Dict{String,Float64}(),
                   marginals=Dict{Tuple{String,Int},Float64}(),
                   t11::Symbol=enabled ? :PASS : :SKIPPED,
                   t12::Symbol=enabled ? :PASS : :SKIPPED,
                   t12_reason::Symbol=enabled ? :ok : :unary_only,
                   t12_result::String=enabled ? "t12-result" : "")
    selected = model == "C1" ? c1_fit : c2_fit
    return E.GraphEvidence(model, selected, "t11-$model", common, t11,
                           "t11-gate-$model", t12, t12_reason, t12_result,
                           Dict{String,Float64}(logz),
                           Dict{Tuple{String,Int},Float64}(marginals))
end

function inner_for(heldout::String, files::Vector{String}; c2_gap::Real=2.0,
                   c1_status::Symbol=:PASS, c1_reason::Symbol=:ok,
                   c2_status::Symbol=:PASS, c2_reason::Symbol=:ok,
                   nodes=nothing)
    nodes === nothing && (nodes = [node_for(file, heldout, index, c2_gap)
                                   for (index, file) in enumerate(files)])
    counts = Dict(file => 0 for file in files)
    return E.InnerFoldEvidence(
        heldout, String[], c1_status, c1_reason, "inner-c1-$heldout",
        c2_status, c2_reason, "inner-c2-$heldout", nodes, counts,
    )
end

function outer_for(index::Int; outer_gap::Real=2.0, inner_gap::Real=2.0,
                   c2_status::Symbol=:PASS, c2_reason::Symbol=:ok,
                   inner_status::Symbol=:PASS,
                   inner_reason::Symbol=:ok,
                   common::String="", eligible_edges::Int=0,
                   outer_nodes=nothing, inner_nodes=nothing,
                   graphs=nothing)
    date = SCAN_DATES[index]
    file = "scan_$date.sxm"
    training = [value for (j, value) in enumerate(SCAN_DATES) if j != index]
    outer_nodes === nothing &&
        (outer_nodes = [node_for(file, date, 1, outer_gap)])
    inner = E.InnerFoldEvidence[]
    for heldout in training
        heldfile = "scan_$heldout.sxm"
        nodes = inner_nodes === nothing ?
            [node_for(heldfile, heldout, 1, inner_gap)] : inner_nodes(heldout)
        push!(inner, E.InnerFoldEvidence(
            heldout,
            [value for value in training if value != heldout],
            inner_status, inner_reason, "inner-c1-$date-$heldout",
            :PASS, :ok, "inner-c2-$date-$heldout", nodes,
            Dict(heldfile => 0),
        ))
    end
    if graphs === nothing
        graphs = Dict(
            "C1" => graph_for("C1", "outer-c1-$date", "outer-c2-$date"; common),
            "C2" => graph_for("C2", "outer-c1-$date", "outer-c2-$date"; common),
        )
    end
    return E.OuterFoldEvidence(
        date, training, inner, :PASS, :ok, "outer-c1-$date",
        c2_status, c2_reason, "outer-c2-$date", common, eligible_edges,
        outer_nodes, E.EdgeNullEvidence[], graphs,
    )
end

function happy_input(; outer_gap::Real=2.0, mutate_index::Int=0,
                     graph_logz::Real=0.0)
    folds = E.OuterFoldEvidence[]
    for index in eachindex(SCAN_DATES)
        gap = index == mutate_index ? outer_gap : 2.0
        push!(folds, outer_for(index; outer_gap=gap,
                               graphs=Dict(
                                   "C1" => graph_for("C1", "outer-c1-$(SCAN_DATES[index])",
                                                     "outer-c2-$(SCAN_DATES[index])";
                                                     logz=Dict("scan_$(SCAN_DATES[index]).sxm" => graph_logz)),
                                   "C2" => graph_for("C2", "outer-c1-$(SCAN_DATES[index])",
                                                     "outer-c2-$(SCAN_DATES[index])";
                                                     logz=Dict("scan_$(SCAN_DATES[index]).sxm" => graph_logz)),
                               )))
    end
    scans = [("scan_$date.sxm", date) for date in SCAN_DATES]
    return E.EvaluationInput(AUTHORITY, UNIVERSE, scans, folds)
end

function replace_outer_node(input::E.EvaluationInput, index::Int, node::E.NodeEvidence)
    old = input.folds[index]
    replacement = E.OuterFoldEvidence(
        old.outer_date, old.training_dates, old.inner_folds,
        old.c1_status, old.c1_reason, old.c1_fit_sha256,
        old.c2_status, old.c2_reason, old.c2_fit_sha256,
        old.common_null_sha256, old.eligible_edge_count, [node], old.edges, old.graphs,
    )
    folds = copy(input.folds)
    folds[index] = replacement
    return E.EvaluationInput(input.authority_sha256, input.universe_sha256,
                             input.scan_dates, folds)
end

function replace_inner_status(input::E.EvaluationInput, fold_index::Int,
                             status::Symbol, reason::Symbol)
    old = input.folds[fold_index]
    first_inner = old.inner_folds[1]
    changed_inner = E.InnerFoldEvidence(
        first_inner.heldout_date, first_inner.training_dates,
        status, reason, first_inner.c1_fit_sha256,
        first_inner.c2_status, first_inner.c2_reason, first_inner.c2_fit_sha256,
        first_inner.nodes, first_inner.eligible_edge_counts,
    )
    inners = copy(old.inner_folds)
    inners[1] = changed_inner
    changed_fold = E.OuterFoldEvidence(
        old.outer_date, old.training_dates, inners,
        old.c1_status, old.c1_reason, old.c1_fit_sha256,
        old.c2_status, old.c2_reason, old.c2_fit_sha256,
        old.common_null_sha256, old.eligible_edge_count,
        old.nodes, old.edges, old.graphs,
    )
    folds = copy(input.folds)
    folds[fold_index] = changed_fold
    return E.EvaluationInput(input.authority_sha256, input.universe_sha256,
                             input.scan_dates, folds)
end

function logsum_pair(a::Real, b::Real)
    x, y = Float64(a), Float64(b)
    m = max(x, y)
    return m + log(exp(x - m) + exp(y - m))
end

function independent_unary(node::E.NodeEvidence, model::Symbol)
    eta0 = Float64[]
    eta1 = Float64[]
    active = 0
    for view in node.views
        c1 = view.c1_log_joint
        c2 = view.c2_log_joint
        all(isnan, c1) && continue
        joints = model == :C1 ? c1 : c2
        all(isfinite, joints) || continue
        a = logsum_pair(joints[1], joints[2])
        p = exp(joints[2] - a)
        p = clamp(p, 1.0e-12, 1.0 - 1.0e-12)
        push!(eta0, a + log1p(-p))
        push!(eta1, a + log(p))
        active += 1
    end
    active == 0 && return (u=0.0, p1=0.5, views=0)
    e0 = mean(eta0)
    e1 = mean(eta1)
    u = logsum_pair(e0, e1)
    return (u=u, p1=exp(e1 - u), views=active)
end

function type7_025(values::Vector{Float64})
    ordered = sort(copy(values))
    h = 1.0 + (length(ordered) - 1) * 0.025
    lo, hi = floor(Int, h), ceil(Int, h)
    lo == hi && return ordered[lo]
    gamma = h - lo
    return (1.0 - gamma) * ordered[lo] + gamma * ordered[hi]
end

function independent_bootstrap(by_date::Dict{String,Vector{Float64}})
    result = Float64[]
    dates = sort(collect(keys(by_date)))
    for seed in 0:499
        rng = MersenneTwister(seed)
        means = Float64[]
        for date in dates
            values = by_date[date]
            sampled = [values[rand(rng, 1:length(values))] for _ in eachindex(values)]
            push!(means, mean(sampled))
        end
        push!(result, mean(means))
    end
    return result
end

function independent_sign_test(means::Vector{Float64})
    k = length(means)
    t_obs = sum(means) / k
    rows = NamedTuple[]
    tail = 0
    for mask in 0:(2^k - 1)
        signs = [((mask >> (index - 1)) & 1) == 1 ? 1 : -1
                 for index in eachindex(means)]
        t = sum(signs[index] * means[index] for index in eachindex(means)) / k
        in_tail = t >= t_obs
        in_tail && (tail += 1)
        push!(rows, (mask=mask, signs=join(signs, ','), t_epsilon=t,
                     in_upper_tail=in_tail))
    end
    return (t_obs=t_obs, rows=rows, tail=tail, denominator=2^k,
            p=tail / (2^k))
end

function independent_binding(provenance::String, tsv::Dict{String,Vector{UInt8}})
    io = IOBuffer()
    for name in sort(collect(keys(tsv)))
        nb = Vector{UInt8}(codeunits(name))
        write(io, UInt64(length(nb))); write(io, nb)
        write(io, UInt64(length(tsv[name]))); write(io, tsv[name])
    end
    tsv_hash = bytes2hex(sha256(take!(io)))
    records = [
        "schema=structured-evaluator-result-binding-v1",
        "provenance_sha256=$provenance",
        "result_tsv_sha256=$tsv_hash",
    ]
    io = IOBuffer()
    for record in records
        bytes = Vector{UInt8}(codeunits(record))
        write(io, UInt64(length(bytes))); write(io, bytes)
    end
    return bytes2hex(sha256(take!(io))), tsv_hash
end

function report_stdout(report::E.EvaluatorReport)
    io = IOBuffer()
    println(io, "status=", report.status)
    println(io, "reason=", report.reason)
    println(io, "result_sha256=", report.result_sha256)
    println(io, "provenance_sha256=", report.provenance_sha256)
    return take!(io)
end

function tab_rows(bytes::Vector{UInt8})
    lines = split(chomp(String(copy(bytes))), '\n')
    return [split(line, '\t') for line in lines[2:end]]
end

function assert_blocked(input::E.EvaluationInput)
    report = E.evaluate(input)
    @test report.status == :BLOCKED
    @test isempty(report.folds)
    @test isempty(report.nodes)
    @test isempty(report.scans)
    @test isempty(report.bootstrap_values)
    @test isempty(report.signs)
    @test report.metrics.eligible_edge_count == 0
    @test report.metrics.sign_denominator == 0
    return report
end

function snapshot_context(input::E.EvaluationInput, manifest::String)
    path = abspath(joinpath(@__DIR__, "evaluate_structured_unit_assignment.jl"))
    bytes = Vector{UInt8}(read(path))
    info = stat(path)
    snap = E._Snapshot(path, bytes, bytes2hex(sha256(bytes)),
                       UInt64(info.device), UInt64(info.inode), UInt64(info.nlink))
    context = E._ProductionContext(pwd(), [snap],
                                   Tuple{String,Vector{String}}[], manifest)
    return context, snap
end

function mixed_status_input(; reverse_folds::Bool=false, fail_only::Bool=false,
                            skipped_only::Bool=false)
    folds = [
        outer_for(1; inner_status=fail_only ? :FAIL : skipped_only ? :PASS : :FAIL,
                  inner_reason=fail_only ? :unary_numerical_failure : :unary_numerical_failure),
        outer_for(2; inner_status=skipped_only ? :PASS : fail_only ? :PASS : :BLOCKED,
                  inner_reason=:authority_or_factor_mismatch),
        outer_for(3),
        outer_for(4; inner_status=skipped_only ? :SKIPPED : :PASS,
                  inner_reason=:scope_unavailable),
        outer_for(5),
    ]
    reverse_folds && reverse!(folds)
    scans = [("scan_$date.sxm", date) for date in SCAN_DATES]
    return E.EvaluationInput(AUTHORITY, UNIVERSE, scans, folds)
end

function local_context(root::String)
    gate = joinpath(root, "gate.json")
    review = joinpath(root, "review.json")
    boulder = joinpath(root, "boulder.json")
    source = joinpath(root, "evaluator.jl")
    input_dir = joinpath(root, "inputs")
    mkpath(input_dir)
    contents = Dict(
        gate => Vector{UInt8}(codeunits("gate")),
        review => Vector{UInt8}(codeunits("review")),
        boulder => Vector{UInt8}(codeunits("boulder")),
        source => Vector{UInt8}(codeunits("source")),
    )
    for (path, bytes) in contents
        write(path, bytes)
    end
    write(joinpath(input_dir, "one.tsv"), "one\n")
    snapshots = E._Snapshot[]
    for path in (gate, review, boulder, source)
        bytes = Vector{UInt8}(read(path))
        info = stat(path)
        push!(snapshots, E._Snapshot(path, bytes, bytes2hex(sha256(bytes)),
                                    UInt64(info.device), UInt64(info.inode), UInt64(info.nlink)))
    end
    context = E._ProductionContext(root, snapshots,
                                   [("inputs", ["one.tsv"])], "local-context")
    return context, (gate=gate, review=review, boulder=boulder, source=source,
                     input_dir=input_dir)
end

function isolated_authority_fixture(root::String;
                                    config_rel::String=HISTORICAL_EVALUATOR_CONFIG)
    source_root = abspath(joinpath(@__DIR__, ".."))
    config = TOML.parse(read(joinpath(source_root, config_rel), String))
    paths = String[config_rel,
                   ".omo/boulder.json",
                   ".omo/plans/structured-label-free-unit-assignment.md",
                   "test/evaluate_structured_unit_assignment.jl",
                   ".omo/evidence/structured-label-free-unit-assignment/t13/evaluator-policy-v1/correction6-GateClosure.json",
                   ".omo/evidence/structured-label-free-unit-assignment/t13/evaluator-policy-v1/correction6-review/AdversarialVerify.json",
                   E._CLOSURE_V1_PUBLICATION_RECEIPT_PATH,
                   E._CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_PATH]
    append!(paths, String(value) for (key, value) in config["authority"]
            if endswith(String(key), "_path"))
    append!(paths, [
        "test/lib/structured_assignment/universe.jl",
        "test/lib/structured_assignment/robust_emissions.jl",
        "test/test_structured_robust_emissions.jl",
        "test/lib/structured_assignment/edge_features.jl",
        "test/build_label_free_edge_features.jl",
        "test/test_structured_edge_features.jl",
        ".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/review/AdversarialVerify.json",
    ])
    if config_rel == RUNTIME_V3_EVALUATOR_CONFIG
        append!(paths, [
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/products/v3/edge_model.jl",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/products/v3/chain_inference.jl",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/products/v3/test_structured_chain_inference.jl",
        ])
    end
    for relative in unique(paths)
        source = joinpath(source_root, relative)
        destination = joinpath(root, relative)
        isfile(source) || error("authority fixture source is absent: $source")
        mkpath(dirname(destination))
        cp(source, destination; force=true)
        chmod(destination, 0o644)
    end
    if config_rel == RUNTIME_V3_EVALUATOR_CONFIG
        root_source = joinpath(source_root, E._RUNTIME_V3_T12_ROOT)
        root_destination = joinpath(root, E._RUNTIME_V3_T12_ROOT)
        ispath(root_destination) && rm(root_destination; recursive=true, force=true)
        cp(root_source, root_destination; force=true)
        for (directory, dirs, files) in walkdir(root_destination)
            chmod(directory, 0o555)
            for child in dirs
                chmod(joinpath(directory, child), 0o555)
            end
            for file in files
                chmod(joinpath(directory, file), 0o444)
            end
        end
        for relative in (
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/root-manifest.sha256",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/DoneClaim.json",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/review/AdversarialVerify.json",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3-publication-receipt.json",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/products/v3/edge_model.jl",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/products/v3/chain_inference.jl",
            ".omo/evidence/structured-label-free-unit-assignment/t12/stable-binary-marginal-normalization-v3/products/v3/test_structured_chain_inference.jl",
        )
            chmod(joinpath(root, relative), 0o444)
        end
    end
    closure_source = joinpath(source_root, E._CLOSURE_V1_ROOT)
    closure_destination = joinpath(root, E._CLOSURE_V1_ROOT)
    cp(closure_source, closure_destination; force=true)
    for (directory, dirs, files) in walkdir(closure_destination)
        chmod(directory, 0o555)
        for child in dirs
            chmod(joinpath(directory, child), 0o555)
        end
        for file in files
            chmod(joinpath(directory, file), 0o444)
        end

    end
    chmod(joinpath(root, E._CLOSURE_V1_PUBLICATION_RECEIPT_PATH), 0o444)
    chmod(joinpath(root, E._CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_PATH), 0o444)
    return (config=joinpath(root, config_rel),
            plan=joinpath(root, ".omo/plans/structured-label-free-unit-assignment.md"),
            boulder=joinpath(root, ".omo/boulder.json"),
            gate=joinpath(root, ".omo/evidence/structured-label-free-unit-assignment/t13/evaluator-policy-v1/correction6-GateClosure.json"),
            prerequisite_review=joinpath(root, ".omo/evidence/structured-label-free-unit-assignment/t13/evaluator-policy-v1/correction6-review/AdversarialVerify.json"),
            claim=joinpath(root, E._CLOSURE_V1_CLAIM_PATH),
            claim_sidecar=joinpath(root, E._CLOSURE_V1_CLAIM_SIDECAR_PATH),
            review=joinpath(root, E._CLOSURE_V1_REVIEW_PATH),
            review_sidecar=joinpath(root, E._CLOSURE_V1_REVIEW_SIDECAR_PATH),
            root_manifest=joinpath(root, E._CLOSURE_V1_ROOT_MANIFEST_PATH),
            publication=joinpath(root, E._CLOSURE_V1_PUBLICATION_RECEIPT_PATH),
            publication_sidecar=joinpath(root, E._CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_PATH),
            noncontrol=joinpath(root, E._CLOSURE_V1_ROOT, "scope-audit.json"),
            closure_root=joinpath(root, E._CLOSURE_V1_ROOT),
            closure_guards=joinpath(root, E._CLOSURE_V1_ROOT, "guards"),
            closure_review=joinpath(root, E._CLOSURE_V1_ROOT, "review"))
end

function rewrite_readonly(path::String, bytes::Vector{UInt8})
    old_mode = UInt(stat(path).mode) & UInt(0o777)
    chmod(path, 0o644)
    write(path, bytes)
    chmod(path, old_mode)
    return nothing
end

function remove_readonly(path::String)
    parent = dirname(path)
    parent_mode = UInt(stat(parent).mode) & UInt(0o777)
    chmod(parent, 0o755)
    rm(path; recursive=true, force=true)
    chmod(parent, parent_mode)
    return nothing
end

function restore_readonly(path::String, bytes::Vector{UInt8})
    parent = dirname(path)
    parent_mode = UInt(stat(parent).mode) & UInt(0o777)
    chmod(parent, 0o755)
    write(path, bytes)
    chmod(parent, parent_mode)
    chmod(path, 0o444)
    return nothing
end

function authority_blocked(f::Function)
    error = try
        f()
        nothing
    catch caught
        caught
    end
    return error isa E.EvaluatorError && error.status == :BLOCKED
end

function make_tree_writable(root::String)
    for (directory, _, files) in walkdir(root)
        chmod(directory, 0o755)
        for file in files
            path = joinpath(directory, file)
            islink(path) || chmod(path, 0o644)
        end
    end
    return nothing
end

function fake_fit_for_edge(; training_date="20240101", heldout_date="20240201",
                           status="unavailable", reason=:residualizer_unavailable,
                           sample=nothing, fit_sha="fit")
    edge = (file="heldout.sxm", date=heldout_date, left_lobe=1,
            right_lobe=2, segment_id="seg")
    state = (status=status, reason=reason, endpoint_predictors=fill(0.0, 14), edge=edge)
    unary = (status=status == "eligible" ? "PASS" : "SKIPPED",
             reason=reason,
             node_identities=[("heldout.sxm", heldout_date, 1, 1.0)])
    samples = sample === nothing ? Any[] : [sample]
    return (training_dates=[training_date], fit_sha256=fit_sha,
            unary=unary, residualizer=nothing, null_fit=nothing,
            transform_edges=[state], samples=samples)
end

function synthetic_selection_context(; full_c2_status="PASS",
                                     partition_c1_status="PASS",
                                     partition_c2_status="PASS")
    dates = ["20240101", "20240102"]
    nodes = [E.T11.AdmissionNode("scan_$date.sxm", date, 1, 1.0, 1.0,
                                 ntuple(_ -> 0.0, 7))
             for date in dates]
    data = E.T11.AdmissionData(nodes, E.T11.AdmissionEdge[], Dict{String,String}())
    normalized = (data=data, normalized_predictors=ones(Float64, length(nodes), 1))
    views(model) = [E.T11.UnaryViewModel(
        name, [1], model == "C1" ? :C1 : :C2, :ok, :none, 2,
        model == "C1" ? fill(5.0, 2, 1) : zeros(2, 1), ones(2, 1))
        for name in E._VIEWS]
    unary(model, date; status="PASS", sha="unary-$model-$date") =
        E.T11.UnaryFitResult(status, :ok, model, [0.5], views(model),
                             [value for value in dates if value != date],
                             ["scan_$date.sxm"],
                             [("scan_$date.sxm", date, 1, 1.0)],
                             "node-$model-$date", sha, "prob-$model-$date")
    fit(model, date; status="PASS", sha="edge-$model-$date") = (
        model_id=model, status=status, reason=:ok,
        training_dates=[value for value in dates if value != date],
        unary=unary(model, date; status, sha="unary-$model-$date"),
        fit_sha256=sha, transform_edges=Any[])
    partitions(model) = [
        (target_date=date,
         training_dates=[value for value in dates if value != date],
         fit=fit(model, date;
                 status=model == "C1" ? partition_c1_status : partition_c2_status))
        for date in dates]
    full_fit(model) = (
        model_id=model, status="PASS", reason=:ok, training_dates=copy(dates),
        unary=unary(model, "full";
                    status=model == "C2" ? full_c2_status : "PASS",
                    sha="unary-$model-full"),
        fit_sha256="edge-$model-full", transform_edges=Any[])
    result(model) = (model_id=model,
                     final_gate=(partitions=partitions(model),),
                     full_fit=full_fit(model))
    metadata = (admission=(models=Dict("C1" => result("C1"),
                                       "C2" => result("C2")),),
                normalized=normalized,
                t12_result_sha256="t12-result",
                t12_provenance_sha256="t12-provenance")
    outer = outer_for(1; c2_status=:SKIPPED,
                      c2_reason=:residualizer_unavailable)
    input = E.EvaluationInput("authority", "universe",
                              [("scan_$date.sxm", date) for date in dates], [outer])
    bindings = (evaluator_config_sha256="evaluator",
                t12_result_sha256="t12-result",
                t12_provenance_sha256="t12-provenance")
    return input, metadata, bindings
end

Test.TESTSET_PRINT_ENABLE[] = false

const RESULTS = @testset "Todo13 evaluator pure seam" begin
    @testset "happy formulas, artifacts, and deterministic replay" begin
        input = happy_input()
        first_report = E.evaluate(input)
        second_report = E.evaluate(input)
        @test first_report.status == :PASS
        @test first_report.reason == :ok
        @test first_report.result_sha256 == second_report.result_sha256
        files1 = E.report_files(first_report)
        files2 = E.report_files(second_report)
        @test files1 == files2
        @test sort(collect(keys(files1))) == sort(collect(ARTIFACTS))
        @test report_stdout(first_report) == report_stdout(second_report)
        @test length(first_report.folds) == 5
        @test length(first_report.nodes) == 5
        @test length(first_report.scans) == 5
        @test length(first_report.dates) == 5
        @test length(first_report.bootstrap_values) == 500
        @test length(first_report.signs) == 2^5
        @test [parse(Int, row[2]) for row in tab_rows(files1["bootstrap.tsv"])] == collect(0:499)

        for name in ARTIFACTS
            @test !isempty(files1[name])
            if !isempty(files1[name])
                @test files1[name][end] == UInt8('\n')
                @test !any(==(UInt8('\r')), files1[name])
            end
        end
        receipt = TOML.parse(String(copy(files1["receipt.toml"])))
        @test receipt["schema"] == "structured_label_free_unit_assignment_evaluator_v1_receipt_v2"
        @test receipt["schema_version"] == 2
        @test receipt["file_count"] == 9
        @test receipt["result_sha256"] == first_report.result_sha256
        @test receipt["provenance_sha256"] == first_report.provenance_sha256
        @test receipt["evaluator_source_sha256"] == first_report.source_sha256
        for key in ("evaluator_config_sha256", "preclosure_plan_sha256",
                    "gateclosure_sha256", "preclosure_boulder_sha256",
                    "review_sha256", "closure_v1_claim_sha256",
                    "closure_v1_review_sha256", "closure_v1_root_manifest_sha256",
                    "closure_v1_publication_receipt_sha256",
                    "closure_v1_publication_receipt_sidecar_sha256")
            @test occursin(r"^[0-9a-f]{64}$", receipt[key])
        end
        @test !haskey(receipt, "plan_sha256")
        @test !haskey(receipt, "boulder_sha256")
        @test receipt["preclosure_plan_sha256"] == E._PLAN_SHA256
        @test receipt["preclosure_boulder_sha256"] == E._BOULDER_SHA256
        @test receipt["live_plan_boulder_runtime_authority"] === false
        @test receipt["closure_v1_claim_sha256"] == E._CLOSURE_V1_CLAIM_SHA256
        @test receipt["closure_v1_review_sha256"] == E._CLOSURE_V1_REVIEW_SHA256
        @test receipt["closure_v1_root_manifest_sha256"] == E._CLOSURE_V1_ROOT_MANIFEST_SHA256
        @test receipt["closure_v1_publication_receipt_sha256"] ==
              E._CLOSURE_V1_PUBLICATION_RECEIPT_SHA256
        @test receipt["closure_v1_publication_receipt_sidecar_sha256"] ==
              E._CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_SHA256
        tsv = Dict(name => files1[name] for name in ARTIFACTS if name != "receipt.toml")
        combined, tsv_hash = independent_binding(first_report.provenance_sha256, tsv)
        @test first_report.result_sha256 == combined
        @test receipt["result_tsv_sha256"] == tsv_hash
        @test first_report.result_sha256 != receipt["result_tsv_sha256"]
        for name in keys(tsv)
            entry = receipt["artifacts"][name]
            @test entry["sha256"] == bytes2hex(sha256(tsv[name]))
            @test entry["bytes"] == length(tsv[name])
        end
        @test [row.outer_date for row in first_report.folds] == SCAN_DATES
        @test [row.date for row in first_report.dates] == SCAN_DATES
        @test [(row.file, row.date) for row in first_report.scans] ==
              [("scan_$date.sxm", date) for date in SCAN_DATES]
        @test [(row.file, row.date) for row in first_report.nodes] ==
              [("scan_$date.sxm", date) for date in SCAN_DATES]
        @test first_report.events[end].ordinal == length(first_report.events)
        @test first_report.events[end].terminal_status == :PASS

        node = input.folds[1].nodes[1]
        want_c1 = independent_unary(node, :C1)
        want_c2 = independent_unary(node, :C2)
        want_d = want_c2.u - want_c1.u
        @test isapprox(first_report.scans[1].d_s, want_d; rtol=1e-13, atol=1e-13)
        @test isapprox(first_report.nodes[1].c1_u, want_c1.u; rtol=1e-13)
        @test isapprox(first_report.nodes[1].selected_u, want_c2.u; rtol=1e-13)
        @test isapprox(first_report.nodes[1].p1, want_c2.p1; rtol=1e-13)
        @test all(isapprox(row.d_k, want_d; rtol=1e-13) for row in first_report.dates)
        by_date = Dict(date => [want_d] for date in SCAN_DATES)
        bootstrap = independent_bootstrap(by_date)
        @test first_report.bootstrap_values == bootstrap
        @test first_report.metrics.bootstrap_lower == type7_025(bootstrap)
        @test first_report.metrics.t_obs == want_d
        @test first_report.metrics.sign_tail_count == 1
        @test first_report.metrics.sign_denominator == 2^5
        @test first_report.metrics.sign_p == 1 / 2^5
        @test all(row.d_k > 0.0 for row in first_report.dates)
        @test first_report.metrics.bootstrap_lower > 0.0
        @test first_report.metrics.sign_p < 0.05
        @test first_report.metrics.coverage == 1.0
        p = want_c2.p1
        want_entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p))
        @test isapprox(first_report.metrics.entropy, want_entropy; rtol=1e-13)
        @test first_report.metrics.agreement_numerator == 6 * 5
        @test first_report.metrics.agreement_denominator == 6 * 5
        @test first_report.metrics.agreement == 1.0
        @test first_report.metrics.unary_parameter_count == 76
        @test first_report.metrics.residualizer_parameter_count == 30
        @test first_report.metrics.null_parameter_count == 5
        @test first_report.metrics.conditional_edge_parameter_count == 9
        @test first_report.metrics.unary_meta_total == 111
        @test first_report.metrics.graph_meta_total == 120
        @test all(row.graph_logz == 0.0 for row in first_report.scans)
        @test all(row.graph_enabled == false for row in first_report.scans)
        @test all(row.in_upper_tail == (row.mask == 31) for row in first_report.signs)
    end

    @testset "T13 unary-selection evidence seam (red contract)" begin
        @test isdefined(E, :UNARY_SELECTION_REASON_CODES)
        @test isdefined(E, :UnarySelectionReference)
        @test isdefined(E, :UnarySelectionFoldEvidence)
        @test isdefined(E, :UnarySelectionScanEvidence)
        @test isdefined(E, :UnarySelectionDateEvidence)
        @test isdefined(E, :UnarySelectionBootstrapEvidence)
        @test isdefined(E, :UnarySelectionDecision)
        @test isdefined(E, :UnarySelectionEvidence)
        @test isdefined(E, :produce_unary_selection_evidence)
        evidence = E._inner_gate_evidence([
            inner_for("20240101", ["one.sxm"]; c2_gap=2.0),
        ])
        @test evidence.decision.model == "C2"
        @test length(evidence.scans) == 1
        @test length(evidence.bootstrap) == 500
    end

    @testset "T13 unary gate evidence is literal and replayable" begin
        expected_reasons = Set([
            :ok, :c2_unavailable_fallback_c1, :nonpositive_date,
            :bootstrap_not_positive, :scope_unavailable,
            :unary_numerical_failure, :denominator_nonpositive,
            :bootstrap_replay_mismatch, :authority_or_factor_mismatch,
            :schema_mismatch, :incomplete_or_stale_evidence,
            :all_views_missing_incident_edge, :selection_reference_mismatch,
        ])
        @test Set(E.UNARY_SELECTION_REASON_CODES) == expected_reasons

        c2 = [inner_for("20240101", ["z.sxm"]; c2_gap=2.0),
              inner_for("20240102", ["a.sxm"]; c2_gap=2.0)]
        reversed_evidence = E._inner_gate_evidence(reverse(c2))
        evidence = E._inner_gate_evidence(c2)
        @test evidence.decision.status == :PASS
        @test evidence.decision.model == "C2"
        @test evidence.decision.reason == :ok
        @test [row.date for row in reversed_evidence.scans] == ["20240101", "20240102"]
        @test [row.date for row in reversed_evidence.dates] == ["20240101", "20240102"]
        @test [row.seed for row in evidence.bootstrap] == collect(0:499)
        expected_by_date = Dict(row.date => [row.delta]
                                for row in evidence.scans)
        @test [row.mean_delta for row in evidence.bootstrap] ==
              independent_bootstrap(expected_by_date)
        second_evidence = E._inner_gate_evidence(c2)
        first_decision = getfield(evidence, :decision)
        second_decision = getfield(second_evidence, :decision)
        @test getproperty(first_decision, :evidence_hash) ==
              getproperty(second_decision, :evidence_hash)

        negative = E._inner_gate_evidence([
            inner_for("20240101", ["negative.sxm"]; c2_gap=-2.0),
        ])
        @test negative.decision.status == :PASS
        @test negative.decision.reason == :nonpositive_date
        @test negative.decision.model == "C1"
        @test negative.decision.selected_reference.model == "C1"

        noisy = E._inner_gate_evidence([
            inner_for("20240101", ["n1.sxm", "n2.sxm"]; c2_gap=-100.0,
                      nodes=[node_for("n1.sxm", "20240101", 1, -100.0),
                             node_for("n2.sxm", "20240101", 2, 3.0)]),
        ])
        @test noisy.decision.status == :PASS
        @test noisy.decision.reason == :bootstrap_not_positive
        @test noisy.decision.model == "C1"

        unavailable = E._inner_gate_evidence([
            inner_for("20240101", ["u.sxm"]; c2_status=:SKIPPED,
                      c2_reason=:residualizer_unavailable),
        ])
        @test unavailable.decision.status == :PASS
        @test unavailable.decision.reason == :c2_unavailable_fallback_c1
        @test unavailable.decision.model == "C1"
        @test isempty(unavailable.bootstrap)

        for (status, want_status, want_reason) in ((:SKIPPED, :SKIPPED, :scope_unavailable),
                                                    (:FAIL, :FAIL, :unary_numerical_failure),
                                                    (:BLOCKED, :BLOCKED,
                                                     :authority_or_factor_mismatch))
            evidence_status = E._inner_gate_evidence([
                inner_for("20240101", ["status.sxm"]; c1_status=status,
                          c1_reason=want_reason),
            ])
            @test evidence_status.decision.status == want_status
            @test evidence_status.decision.reason == want_reason
            @test evidence_status.decision.selected_reference === nothing
        end

        denominator = E._inner_gate_evidence([
            inner_for("20240101", ["zero.sxm"]; nodes=[node_for(
                "zero.sxm", "20240101", 1, 2.0)],),
        ])
        old = denominator.folds[1]
        # The typed seam must retain the selected candidate identities even in
        # the hand fixture; the arithmetic guard itself is checked through the
        # legacy input constructor below.
        @test old.c1_reference.fit_role == "partition"
        bad_counts = E.InnerFoldEvidence(
            "20240101", String[], :PASS, :ok, "c1", :PASS, :ok, "c2",
            [node_for("zero.sxm", "20240101", 1, 2.0)], Dict("zero.sxm" => -1))
        bad_denominator = E._inner_gate_evidence([bad_counts])
        @test bad_denominator.decision.status == :FAIL
        @test bad_denominator.decision.reason == :denominator_nonpositive

        Random.seed!(1234)
        first_draw = rand()
        second_draw = rand()
        Random.seed!(1234)
        @test rand() == first_draw
        E._inner_gate_evidence(c2)
        @test rand() == second_draw

        p0 = E.NodeEvidence("p.sxm", "20240101", 1, 1.0,
                           [E.ViewJointEvidence("base_local", (0.0, 0.0),
                                                (0.0, 0.0))])
        p1 = E.NodeEvidence("p.sxm", "20240101", 1, 1.0,
                           [E.ViewJointEvidence("base_local", (0.0, 0.0),
                                                (5.0, 5.0))])
        d0 = E._inner_gate_evidence([inner_for("20240101", ["p.sxm"];
                                               nodes=[p0])])
        d1 = E._inner_gate_evidence([inner_for("20240101", ["p.sxm"];
                                               nodes=[p1])])
        @test d0.scans[1].delta != d1.scans[1].delta
        @test d0.scans[1].c1_log_evidence == d1.scans[1].c1_log_evidence

        seam_call(features) = E.produce_unary_selection_evidence(pwd();
            evaluator_config=RUNTIME_EVALUATOR_CONFIG,
            features=features,
            candidate_config="config/unit_assignment_structured_candidate.toml",
            model_config="config/unit_assignment_structured_model.toml",
            universe_dir="missing-universe", edge_dir="missing-edge",
            forward_receipt="missing-forward.toml",
            backward_receipt="missing-backward.toml", admission_dir="missing-admission")
        @test_throws E.EvaluatorError seam_call("truth.tsv")
        @test_throws E.EvaluatorError seam_call("expected_count.tsv")
        @test_throws E.EvaluatorError seam_call("/tmp/features.tsv")

        @testset "S2 reference and serialization defects (red)" begin
            outer = outer_for(1)
            outer_reference = E._selection_outer_reference(outer, "C1")
            @test outer_reference.fit_sha256 == "t11-C1"
            @test outer_reference.unary_fit_sha256 == "outer-c1-20240101"
            inner_reference = E._inner_reference(
                inner_for("20240101", ["inner.sxm"]), "C1", "outer", "20240109")
            @test inner_reference.fit_sha256 != inner_reference.unary_fit_sha256
            @test E._selection_effective_outer_model(outer_for(1; c2_status=:SKIPPED),
                                                     "C2") == "C1"
            outer_c1_skipped = E.OuterFoldEvidence(
                outer.outer_date, outer.training_dates, outer.inner_folds,
                :FAIL, :unary_numerical_failure, outer.c1_fit_sha256,
                outer.c2_status, outer.c2_reason, outer.c2_fit_sha256,
                outer.common_null_sha256, outer.eligible_edge_count,
                outer.nodes, outer.edges, outer.graphs)
            @test_throws E.EvaluatorError E._selection_effective_outer_model(
                outer_c1_skipped, "C2")

            scan = E.UnarySelectionScanEvidence(
                "outer", "20240109", "20240101", "scan.sxm", 1, 0, 1,
                0.0, 1.0, 1.0)
            date = E.UnarySelectionDateEvidence("outer", "20240109", "20240101",
                                                1, 1.0, true)
            @test scan.outer_date == "20240109"
            @test scan.date == "20240101"
            @test date.outer_date == "20240109"

            c1a = E.UnarySelectionReference("C1", "partition", "outer_inner",
                                             "20240109", "20240101", "edge",
                                             "unary-a")
            c1b = E.UnarySelectionReference("C1", "partition", "outer_inner",
                                             "20240109", "20240101", "edge",
                                             "unary-b")
            fold_a = E.UnarySelectionFoldEvidence(
                "outer", "20240109", "20240101", String[], :PASS, :ok,
                c1a, :PASS, :ok, c1a)
            fold_b = E.UnarySelectionFoldEvidence(
                "outer", "20240109", "20240101", String[], :PASS, :ok,
                c1b, :PASS, :ok, c1b)
            @test E._selection_evidence_hash([fold_a], [scan], [date],
                                              E.UnarySelectionBootstrapEvidence[]) !=
                  E._selection_evidence_hash([fold_b], [scan], [date],
                                              E.UnarySelectionBootstrapEvidence[])

            fallback_decision = getfield(E._inner_gate_evidence([
                inner_for("20240101", ["fallback.sxm"];
                          c2_status=:SKIPPED,
                          c2_reason=:residualizer_unavailable)]), :decision)
            @test fallback_decision.hash == fallback_decision.decision_sha256
            @test !isempty(fallback_decision.decision_sha256)
        end

        @testset "S2 public-path firewall and root defects (red)" begin
            arguments = [RUNTIME_EVALUATOR_CONFIG,
                         "features.tsv", "config/unit_assignment_structured_candidate.toml",
                         "config/unit_assignment_structured_model.toml", "universe",
                         "edge", "forward.toml", "backward.toml", "admission"]
            no_output_options = String[value for (flag, value) in zip(
                ("--root", "--evaluator-config", "--features", "--candidate-config",
                 "--model-config", "--universe-dir", "--edge-dir", "--forward-receipt",
                 "--backward-receipt", "--admission-dir"),
                (pwd(), arguments...)) for value in (flag, value)]
            no_output_options_map = E._selection_cli_options(no_output_options)
            @test no_output_options_map["--root"] == pwd()
            @test !haskey(no_output_options_map, "--out-dir")
            root_error = try
                E.produce_unary_selection_evidence(joinpath(pwd(), "absent-root");
                    evaluator_config=arguments[1], features=arguments[2],
                    candidate_config=arguments[3], model_config=arguments[4],
                    universe_dir=arguments[5], edge_dir=arguments[6],
                    forward_receipt=arguments[7], backward_receipt=arguments[8],
                    admission_dir=arguments[9])
                nothing
            catch error
                error
            end
            @test root_error isa E.EvaluatorError
            @test root_error.reason == :authority_path_invalid

            encoded_error = try
                seam_call(raw"features\bad")
                nothing
            catch error
                error
            end
            @test encoded_error isa E.EvaluatorError
            @test encoded_error.reason == :cli_error
            for bad_value in ("features%2Fbad", "featuresé.tsv")
                bad_error = try
                    seam_call(bad_value)
                    nothing
                catch error
                    error
                end
                @test bad_error isa E.EvaluatorError
                @test bad_error.reason == :cli_error
            end
            missing_context = E._ProductionContext(
                pwd(), E._Snapshot[], Tuple{String,Vector{String}}[], "manifest")
            missing_error = try
                E._selection_binding_digest(missing_context,
                                            joinpath(pwd(), "missing-input"), "test")
                nothing
            catch error
                error
            end
            @test missing_error isa E.EvaluatorError
            @test missing_error.reason == :incomplete_or_stale_evidence
        end

        @testset "S2 synthetic production-context success and terminal state" begin
            synthetic_input, synthetic_metadata, synthetic_bindings =
                synthetic_selection_context()
            synthetic = E._produce_unary_selection_evidence_from_context(
                synthetic_input, synthetic_metadata, synthetic_bindings)
            @test synthetic.status == :PASS
            @test synthetic.reason == :ok
            @test synthetic.final_decision.status == :PASS
            @test synthetic.final_decision.selected_model == "C2"
            @test synthetic.final_decision.selected_reference.fit_role == "full_refit"
            @test synthetic.final_decision.selected_reference.fit_sha256 == "edge-C2-full"
            @test synthetic.final_decision.selected_reference.unary_fit_sha256 ==
                  "unary-C2-full"
            @test synthetic.outer_decisions[1].selected_model == "C1"
            @test synthetic.outer_decisions[1].selected_reference.fit_sha256 == "t11-C1"
            @test synthetic.outer_decisions[1].selected_reference.unary_fit_sha256 ==
                  "outer-c1-20240101"
            @test synthetic.outer_decisions[1].reason == :c2_unavailable_fallback_c1
            @test all(row.outer_date == "20240101" for row in synthetic.scans
                      if row.decision_scope == "outer")
            @test all(row.outer_date == "NA" for row in synthetic.scans
                      if row.decision_scope == "final")
            @test synthetic.t12_result_sha256 == "t12-result"
            @test synthetic.t12_provenance_sha256 == "t12-provenance"
            old_hash = synthetic.evidence_hash
            changed = E._selection_top_hash(synthetic.outer_decisions,
                                            synthetic.final_decision,
                                            merge(synthetic.bindings,
                                                  (extra="different",)))
            @test old_hash != changed
            selected = synthetic.outer_decisions[1].selected_reference
            changed_selected = E.UnarySelectionReference(
                selected.model, selected.fit_role, selected.scope,
                selected.outer_date, selected.target_date, "changed-edge",
                "changed-unary")
            changed_outer = E._selection_rebind_decision(
                synthetic.outer_decisions[1],
                synthetic.outer_decisions[1].c1_reference,
                synthetic.outer_decisions[1].c2_reference,
                changed_selected; selected_model="C1")
            @test changed_outer.evidence_hash != synthetic.outer_decisions[1].evidence_hash
            @test E._selection_top_hash(
                [changed_outer], synthetic.final_decision,
                synthetic.bindings) != old_hash

            unavailable_input, unavailable_metadata, unavailable_bindings =
                synthetic_selection_context(; full_c2_status="SKIPPED")
            unavailable = E._produce_unary_selection_evidence_from_context(
                unavailable_input, unavailable_metadata, unavailable_bindings)
            @test unavailable.status == :SKIPPED
            @test unavailable.reason == :scope_unavailable
            @test unavailable.final_decision.selected_reference === nothing
            @test unavailable.final_decision.upstream_reason == :ok

            terminal_input, terminal_metadata, terminal_bindings =
                synthetic_selection_context(; partition_c1_status="SKIPPED")
            terminal_synthetic = E._produce_unary_selection_evidence_from_context(
                terminal_input, terminal_metadata, terminal_bindings)
            @test terminal_synthetic.status == :SKIPPED
            @test terminal_synthetic.reason == :scope_unavailable
            @test terminal_synthetic.final_decision.selected_reference === nothing

            @testset "final candidate references bind every terminal decision" begin
                c1 = E.UnarySelectionReference("C1", "partition", "final_lodo",
                                               "NA", "NA", "c1-fit", "c1-unary")
                c2 = E.UnarySelectionReference("C2", "partition", "final_lodo",
                                               "NA", "NA", "c2-fit", "c2-unary")
                bindings = (evaluator_config_sha256="historical",)
                terminal_cases = (
                    (:PASS, :ok, "C1", c1),
                    (:SKIPPED, :scope_unavailable, "", nothing),
                    (:FAIL, :unary_numerical_failure, "", nothing),
                    (:BLOCKED, :authority_or_factor_mismatch, "", nothing),
                )
                for (status, reason, model, selected) in terminal_cases
                    base = E._selection_decision(
                        "final", "NA", status, reason, model, selected, c1, c2,
                        2, 2, 2, 0.25, true, "", "", :upstream)
                    base_top = E._selection_top_hash(
                        E.UnarySelectionDecision[], base, bindings)
                    @test !isempty(base.decision_hash)
                    @test !isempty(base_top)
                    for (candidate_index, reference) in ((1, c1), (2, c2))
                        for field_index in eachindex(fieldnames(E.UnarySelectionReference))
                            values = [getfield(reference, field)
                                      for field in fieldnames(E.UnarySelectionReference)]
                            values[field_index] *= "-changed"
                            changed_reference = E.UnarySelectionReference(values...)
                            changed_c1 = candidate_index == 1 ? changed_reference : c1
                            changed_c2 = candidate_index == 2 ? changed_reference : c2
                            changed_selected = status == :PASS ?
                                (model == "C1" ? changed_c1 : changed_c2) : nothing
                            changed = E._selection_decision(
                                "final", "NA", status, reason, model, changed_selected,
                                changed_c1, changed_c2, 2, 2, 2, 0.25, true, "", "",
                                :upstream)
                            changed_top = E._selection_top_hash(
                                E.UnarySelectionDecision[], changed, bindings)
                            @test changed.decision_hash != base.decision_hash
                            @test changed_top != base_top
                        end
                    end
                end
            end
        end

        @testset "T13 public producer temporary-fixture boundary" begin
            source_root = realpath(joinpath(@__DIR__, ".."))
            if VERSION != v"1.12.6"
                public_error = try
                    E.produce_unary_selection_evidence(source_root;
                        evaluator_config=RUNTIME_V3_EVALUATOR_CONFIG,
                        features=".tmp-t13-public-fixture-missing/features.tsv",
                        candidate_config="config/unit_assignment_structured_candidate.toml",
                        model_config="config/unit_assignment_structured_model.toml",
                        universe_dir=".tmp-t13-public-fixture-missing/universe",
                        edge_dir=".tmp-t13-public-fixture-missing/edges",
                        forward_receipt=".tmp-t13-public-fixture-missing/forward.toml",
                        backward_receipt=".tmp-t13-public-fixture-missing/backward.toml",
                        admission_dir=".tmp-t13-public-fixture-missing/admission")
                    nothing
                catch error
                    error
                end
                @test public_error isa E.EvaluatorError
                @test public_error.reason == :runtime_version_mismatch
            elseif get(ENV, "STMFIT_T13_PUBLIC_REVALIDATION", "0") != "1"
                @test get(ENV, "STMFIT_T13_PUBLIC_REVALIDATION", "0") == "0"
            else
                fixture = EDGE_FIXTURE_SUPPORT.make_cli_fixture_under(source_root;
                                                                       cleanup=false)
                root = fixture.root
                try
                    feature_sha = EDGE_FIXTURE_SUPPORT.file_sha256(fixture.original_features)
                    keys_sha = bytes2hex(sha256(
                        E.T11.StructuredEdgeFeatures.StructuredUniverse._key_identity_bytes(
                            Tuple(E.T11.StructuredEdgeFeatures.StructuredUniverse.LobeKey(
                                node.file, node.lobe) for node in fixture.data.nodes),
                        ),
                    ))
                    commands = E.T11.StructuredEdgeFeatures.expected_producer_commands(
                        source_root;
                        features=fixture.original_features,
                        data_dir=fixture.producer_data,
                        forward_patches=fixture.forward_patches,
                        backward_patches=fixture.backward_patches,
                    )
                    write(fixture.forward_receipt,
                        EDGE_FIXTURE_SUPPORT.synthetic_patch_receipt_text(
                            source_root, :forward, fixture.original_features,
                            fixture.producer_data, fixture.forward_patches, feature_sha,
                            keys_sha, length(fixture.data.nodes), commands.forward))
                    write(fixture.backward_receipt,
                        EDGE_FIXTURE_SUPPORT.synthetic_patch_receipt_text(
                            source_root, :backward, fixture.original_features,
                            fixture.producer_data, fixture.backward_patches, feature_sha,
                            keys_sha, length(fixture.data.nodes), commands.backward))
                    edge_receipt = joinpath(fixture.edges, "receipt.toml")
                    EDGE_FIXTURE_SUPPORT.replace_toml_value!(
                        edge_receipt, "forward_patch_receipt_sha256",
                        EDGE_FIXTURE_SUPPORT.file_sha256(fixture.forward_receipt))
                    EDGE_FIXTURE_SUPPORT.replace_toml_value!(
                        edge_receipt, "backward_patch_receipt_sha256",
                        EDGE_FIXTURE_SUPPORT.file_sha256(fixture.backward_receipt))
                    EDGE_FIXTURE_SUPPORT.replace_toml_value!(
                        edge_receipt, "forward_command", commands.forward)
                    EDGE_FIXTURE_SUPPORT.replace_toml_value!(
                        edge_receipt, "backward_command", commands.backward)
                    admission_dir = joinpath(root, "synthetic", "admission")
                    mkpath(admission_dir)
                    data = E.T11.load_admission_data(source_root;
                        features=relpath(fixture.features, source_root),
                        candidate_config="config/unit_assignment_structured_candidate.toml",
                        model_config="config/unit_assignment_structured_model.toml",
                        universe_dir=relpath(fixture.universe, source_root),
                        edge_dir=relpath(fixture.edges, source_root),
                        forward_receipt=relpath(fixture.forward_receipt, source_root),
                        backward_receipt=relpath(fixture.backward_receipt, source_root))
                    admission = E.T11.evaluate_admission(data)
                    for (name, bytes) in E.T11.report_files(admission)
                        write(joinpath(admission_dir, name), bytes)
                    end
                    public_call() = E.produce_unary_selection_evidence(source_root;
                        evaluator_config=RUNTIME_V3_EVALUATOR_CONFIG,
                        features=relpath(fixture.features, source_root),
                        candidate_config="config/unit_assignment_structured_candidate.toml",
                        model_config="config/unit_assignment_structured_model.toml",
                        universe_dir=relpath(fixture.universe, source_root),
                        edge_dir=relpath(fixture.edges, source_root),
                        forward_receipt=relpath(fixture.forward_receipt, source_root),
                        backward_receipt=relpath(fixture.backward_receipt, source_root),
                        admission_dir=relpath(admission_dir, source_root))
                    first_evidence = try
                        public_call()
                    catch error
                        error
                    end
                    if !(first_evidence isa E.UnarySelectionEvidence)
                        println(stderr, "PUBLIC_PRODUCER_FAILURE: ",
                                sprint(showerror, first_evidence))
                    end
                    @test first_evidence isa E.UnarySelectionEvidence
                    if first_evidence isa E.UnarySelectionEvidence
                        @test first_evidence.status isa Symbol
                        @test first_evidence.status in (:PASS, :SKIPPED, :FAIL)
                        @test !isempty(first_evidence.outer_decisions)
                        @test first_evidence.final_decision isa E.UnarySelectionDecision
                        @test !isempty(first_evidence.t12_result_sha256)
                        @test !isempty(first_evidence.t12_provenance_sha256)
                        second_evidence = public_call()
                        @test second_evidence isa E.UnarySelectionEvidence
                        @test second_evidence.evidence_hash == first_evidence.evidence_hash
                        @test second_evidence.bindings == first_evidence.bindings
                        println("PUBLIC_REPLAY_V3",
                                " status=", first_evidence.status,
                                " evidence_hash=", first_evidence.evidence_hash,
                                " authority_sha256=", first_evidence.authority_sha256,
                                " t12_result_sha256=", first_evidence.t12_result_sha256,
                                " t12_provenance_sha256=", first_evidence.t12_provenance_sha256,
                                " final_decision_sha256=",
                                first_evidence.final_decision.decision_hash)
                    end
                finally
                    make_tree_writable(root)
                    rm(root; recursive=true, force=true)
                end
            end
            @test isempty(filter(name -> startswith(name, ".tmp-t13-public-fixture-"),
                                 readdir(source_root)))
        end
    end

    @testset "pooled date gate and fallback evidence" begin
        date_a = [node_for("a1.sxm", "20240101", 1, -100.0),
                  node_for("a2.sxm", "20240101", 2, -100.0),
                  node_for("a3.sxm", "20240101", 3, 3.0)]
        date_b = [node_for("b1.sxm", "20240102", 1, 2.0)]
        pooled = [
            E.InnerFoldEvidence("20240101", String[], :PASS, :ok, "i1",
                                :PASS, :ok, "j1", date_a,
                                Dict("a1.sxm" => 0, "a2.sxm" => 0, "a3.sxm" => 0)),
            E.InnerFoldEvidence("20240102", String[], :PASS, :ok, "i2",
                                :PASS, :ok, "j2", date_b, Dict("b1.sxm" => 0)),
        ]
        decision = E._inner_gate(pooled)
        @test decision.model == "C2"
        @test decision.reason == :ok
        @test length(decision.hash) == 64
        negative = deepcopy(pooled)
        negative[1] = E.InnerFoldEvidence("20240101", String[], :PASS, :ok, "i1",
                                          :PASS, :ok, "j1",
                                          [node_for("a1.sxm", "20240101", 1, -100.0)],
                                          Dict("a1.sxm" => 0))
        negative_decision = E._inner_gate(negative)
        @test negative_decision.reason == :nonpositive_date
        @test negative_decision.hash != decision.hash
        noisy = [E.InnerFoldEvidence(
            "20240101", String[], :PASS, :ok, "i3", :PASS, :ok, "j3",
            [node_for("n1.sxm", "20240101", 1, -100.0),
             node_for("n2.sxm", "20240101", 2, 3.0)],
            Dict("n1.sxm" => 0, "n2.sxm" => 0),
        )]
        noisy_decision = E._inner_gate(noisy)
        @test noisy_decision.reason == :bootstrap_not_positive
        @test noisy_decision.hash != decision.hash
        unavailable = E.InnerFoldEvidence(
            "20240101", String[], :PASS, :ok, "i4", :SKIPPED,
            :residualizer_unavailable, "j4", [node_for("u.sxm", "20240101", 1, 2.0)],
            Dict("u.sxm" => 0))
        @test E._inner_gate([unavailable]).reason == :c2_unavailable_fallback_c1
        outer = outer_for(1; c2_status=:SKIPPED,
                          c2_reason=:residualizer_unavailable,
                          inner_status=:PASS)
        @test E._score_fold(outer, "C2").reason == :c2_unavailable_fallback_c1
    end

    @testset "terminal precedence and exhaustive signs" begin
        normal = E.evaluate(mixed_status_input())
        reversed = E.evaluate(mixed_status_input(; reverse_folds=true))
        for report in (normal, reversed)
            @test report.status == :BLOCKED
            @test report.reason == :authority_or_factor_mismatch
            @test isempty(report.folds) && isempty(report.nodes) && isempty(report.scans)
            @test isempty(report.bootstrap_values) && isempty(report.signs)
            @test report.metrics.eligible_edge_count == 0
            @test report.metrics.sign_denominator == 0
            @test report.metrics.graph_enabled_fold_count == 0
            @test length(report.events) == 1
            event = only(report.events)
            @test event.ordinal == 1 && event.consumed == false
            @test event.terminal_status == :BLOCKED
            @test event.outer_date == "20240102"
            @test event.reason == :authority_or_factor_mismatch
            @test event.consequence == :authority_or_factor_mismatch
            @test event.formula_count == 0 && event.bootstrap_count == 0
            @test event.sign_count == 0 && event.graph_count == 0
        end
        @test normal.result_sha256 == reversed.result_sha256
        winning_receipt = TOML.parse(String(E._blocker_files(
            E.EvaluatorError(:BLOCKED, normal.reason, "winning"))["receipt.toml"]))
        @test winning_receipt["reason"] == String(normal.reason)

        fail_report = E.evaluate(mixed_status_input(; fail_only=true))
        @test fail_report.status == :FAIL
        @test fail_report.reason == :unary_numerical_failure
        @test length(fail_report.events) == 1
        @test fail_report.events[1].outer_date == "20240101"
        @test fail_report.events[1].terminal_status == :FAIL
        @test fail_report.events[1].reason == :unary_numerical_failure
        @test isempty(fail_report.bootstrap_values)

        skipped_report = E.evaluate(mixed_status_input(; skipped_only=true))
        @test skipped_report.status == :SKIPPED
        @test skipped_report.reason == :scope_unavailable
        @test length(skipped_report.events) == 1
        @test skipped_report.events[1].outer_date == "20240104"
        @test skipped_report.events[1].terminal_status == :SKIPPED
        @test skipped_report.events[1].consequence == :scope_unavailable
        malformed = E.evaluate(E.EvaluationInput("", "", Tuple{String,String}[],
                                                 E.OuterFoldEvidence[]))
        @test malformed.status == :BLOCKED
        @test length(malformed.events) == 1
        @test malformed.events[1].ordinal == 1
        @test malformed.events[1].outer_date == ""
        @test malformed.events[1].consumed == false

        for k in 1:13
            means = [index == 2 ? 0.0 : (isodd(index) ? -0.37 * index : 0.21 * index)
                     for index in 1:k]
            want = independent_sign_test(means)
            got = E._exhaustive_sign_test(means)
            @test got.t_obs == want.t_obs
            @test got.tail == want.tail
            @test got.denominator == 2^k
            @test got.p == want.p
            @test length(got.rows) == length(want.rows)
            for (actual, expected) in zip(got.rows, want.rows)
                @test actual.mask == expected.mask
                @test actual.signs == expected.signs
                @test actual.t_epsilon == expected.t_epsilon
                @test actual.in_upper_tail == expected.in_upper_tail
            end
        end
        for k in (1, 4, 13)
            tied = E._exhaustive_sign_test(zeros(k))
            @test tied.tail == 2^k
            @test tied.denominator == 2^k
            @test tied.p == 1.0
            @test all(row.in_upper_tail for row in tied.rows)
        end
        all4 = E._exhaustive_sign_test(ones(4))
        all5 = E._exhaustive_sign_test(ones(5))
        @test all4.tail == 1 && all4.p == 1 / 16 && !(all4.p < 0.05)
        @test all5.tail == 1 && all5.p == 1 / 32 && all5.p < 0.05
    end

    @testset "shared masks and fixed edge lifecycle" begin
        structural = E.NodeEvidence("x.sxm", "20240101", 1, 1.0,
                                    [E.ViewJointEvidence(name, (NaN, NaN),
                                                         (0.0, 2.0)) for name in VIEW_NAMES])
        c1_missing = E._unary(structural, :C1)
        c2_missing = E._unary(structural, :C2)
        @test c1_missing.views == 0 && c2_missing.views == 0
        @test c1_missing.u == 0.0 && c2_missing.p1 == 0.5
        @test E._hard(0.5) == "?"
        active_missing = E.NodeEvidence("x.sxm", "20240101", 1, 1.0,
                                       [E.ViewJointEvidence("base_local", (0.0, 1.0),
                                                             (NaN, NaN))])
        @test E._unary(active_missing, :C1).views == 1
        c2_state = E._unary(active_missing, :C2)
        @test c2_state.views == 0 && c2_state.c2_unavailable
        partial = E.NodeEvidence("x.sxm", "20240101", 1, 1.0,
                                 [E.ViewJointEvidence("base_local", (NaN, 1.0), (0.0, 1.0))])
        partial_error = try E._unary(partial, :C1); nothing catch err; err end
        @test partial_error isa E.EvaluatorError
        @test partial_error.status == :FAIL
        infinite = E.NodeEvidence("x.sxm", "20240101", 1, 1.0,
                                   [E.ViewJointEvidence("base_local", (Inf, 1.0), (0.0, 1.0))])
        @test (try E._unary(infinite, :C1); false catch err; err isa E.EvaluatorError && err.status == :FAIL end)
        edge = E.EdgeNullEvidence("x.sxm", "20240101", 1, 2, "seg", 7.0)
        missing_nodes = [E.NodeEvidence("x.sxm", "20240101", lobe, Float64(lobe),
                                        [E.ViewJointEvidence(name, (NaN, NaN), (NaN, NaN))
                                         for name in VIEW_NAMES]) for lobe in (1, 2)]
        blocked_fold = E.OuterFoldEvidence(
            "20240101", String[], E.InnerFoldEvidence[], :PASS, :ok, "c1",
            :PASS, :ok, "c2", "null", 1, missing_nodes, [edge],
            Dict("C1" => graph_for("C1", "c1", "c2"; common="null"),
                 "C2" => graph_for("C2", "c1", "c2"; common="null")),
        )
        blocked = try E._score_fold(blocked_fold, "C1"); nothing catch err; err end
        @test blocked isa E.EvaluatorError && blocked.status == :BLOCKED
        unavailable_fit = fake_fit_for_edge()
        @test E._edge_counts(unavailable_fit, "20240201")["heldout.sxm"] == 1
        _, _, available = E._common_null(unavailable_fit, unavailable_fit, "20240201")
        @test !available
        null_missing = outer_for(1; common="", eligible_edges=1)
        skipped = E._score_fold(null_missing, "C1")
        @test skipped.status == :SKIPPED && skipped.reason == :common_null_unavailable
        @test isempty(skipped.nodes) && isempty(skipped.scans)
        unary_only = outer_for(1; common="", eligible_edges=0)
        unary_score = E._score_fold(unary_only, "C1")
        @test unary_score.status == :PASS
        @test unary_score.scans[1].eligible_edge_count == 0
        one_edge = outer_for(1; common="null", eligible_edges=1,
                             outer_nodes=[node_for("scan_20240101.sxm", "20240101", 1, 2.0),
                                          node_for("scan_20240101.sxm", "20240101", 2, 2.0)])
        edge0 = E.OuterFoldEvidence(
            one_edge.outer_date, one_edge.training_dates, one_edge.inner_folds,
            one_edge.c1_status, one_edge.c1_reason, one_edge.c1_fit_sha256,
            one_edge.c2_status, one_edge.c2_reason, one_edge.c2_fit_sha256,
            one_edge.common_null_sha256, 0, one_edge.nodes, E.EdgeNullEvidence[],
            one_edge.graphs,
        )
        edge1 = E.EdgeNullEvidence("scan_20240101.sxm", "20240101", 1, 2, "seg", 7.0)
        with_edge = E.OuterFoldEvidence(
            one_edge.outer_date, one_edge.training_dates, one_edge.inner_folds,
            one_edge.c1_status, one_edge.c1_reason, one_edge.c1_fit_sha256,
            one_edge.c2_status, one_edge.c2_reason, one_edge.c2_fit_sha256,
            one_edge.common_null_sha256, 1, one_edge.nodes, [edge1], one_edge.graphs,
        )
        s0, s1 = E._score_fold(edge0, "C2"), E._score_fold(with_edge, "C2")
        numerator0 = s0.scans[1].meta_loglik - s0.scans[1].c1_loglik
        numerator1 = s1.scans[1].meta_loglik - s1.scans[1].c1_loglik
        @test isapprox(numerator0, numerator1; rtol=1e-13)
        @test isapprox(s0.scans[1].d_s, numerator0 / 2; rtol=1e-13)
        @test isapprox(s1.scans[1].d_s, numerator0 / 3; rtol=1e-13)
        @test s1.scans[1].denominator == 3
        @test s0.scans[1].denominator == 2
        @test s0.scans[1].eligible_edge_count == 0
        @test s1.scans[1].eligible_edge_count == 1
        @test s0.scans[1].common_null_loglik == 0.0
        @test s1.scans[1].common_null_loglik == 7.0

        null_fit = E.T11.NullFit(:ok, :ok, [0.0, 0.0], [1.0 0.0; 0.0 1.0],
                                 0.0, Float64[], true, 1, NaN)
        train_edge = (file="train.sxm", date="20240101", left_lobe=1,
                      right_lobe=2, segment_id="seg")
        train_state = (status="eligible", reason=:none,
                       endpoint_predictors=fill(0.0, 14), edge=train_edge)
        train_unary = (status="PASS", reason=:ok,
                       node_identities=[("train.sxm", "20240101", 1, 1.0)])
        sample1 = (edge=train_edge, endpoint_predictors=fill(0.0, 14),
                   residual=[0.1, 0.2], raw_sha256="raw")
        sample2 = (edge=train_edge, endpoint_predictors=fill(0.0, 14),
                   residual=[0.3, 0.2], raw_sha256="raw")
        c1 = (training_dates=["20240101"], fit_sha256="same-c1", unary=train_unary,
              residualizer=nothing, null_fit=null_fit, transform_edges=[train_state],
              samples=[sample1])
        c1_changed = merge(c1, (samples=[sample2],))
        c2_state = (status="unavailable", reason=:residualizer_unavailable,
                    endpoint_predictors=fill(0.0, 14), edge=train_edge)
        c2 = (training_dates=["20240101"], fit_sha256="same-c2",
              unary=(status="SKIPPED", reason=:residualizer_unavailable,
                     node_identities=Tuple{String,String,Int,Float64}[]),
              residualizer=nothing, null_fit=nothing, transform_edges=[c2_state],
              samples=Any[])
        cache = Dict{String,Any}()
        E._common_null(c1, c2, "20240201"; cache)
        E._common_null(c1_changed, c2, "20240201"; cache)
        @test length(cache) == 2
    end

    @testset "graph, leakage, and topology boundaries" begin
        input = happy_input()
        base = E.evaluate(input)
        extreme_node = node_for("scan_20240101.sxm", "20240101", 1, 20.0)
        changed = E.evaluate(replace_outer_node(input, 1, extreme_node))
        @test [row.inner_hash for row in base.folds] == [row.inner_hash for row in changed.folds]
        @test [row.selected_fit for row in base.folds] == [row.selected_fit for row in changed.folds]
        @test base.scans[1].d_s != changed.scans[1].d_s
        @test base.scans[2:end] == changed.scans[2:end]

        node = node_for("g.sxm", "20240101", 1, 2.0)
        disabled = outer_for(1; outer_nodes=[node])
        disabled_score = E._score_fold(disabled, "C2")
        @test disabled_score.scans[1].graph_enabled == false
        @test disabled_score.nodes[1].p1 == independent_unary(node, :C2).p1
        @test disabled_score.scans[1].graph_logz == 0.0
        enabled_graphs = Dict(
            "C1" => graph_for("C1", "outer-c1-20240101", "outer-c2-20240101";
                               enabled=true, logz=Dict("g.sxm" => 2.5),
                               marginals=Dict(("g.sxm", 1) => 0.25)),
            "C2" => graph_for("C2", "outer-c1-20240101", "outer-c2-20240101";
                               enabled=true, logz=Dict("g.sxm" => 2.5),
                               marginals=Dict(("g.sxm", 1) => 0.25)),
        )
        enabled = outer_for(1; outer_nodes=[node], graphs=enabled_graphs)
        enabled_score = E._score_fold(enabled, "C2")
        @test enabled_score.scans[1].graph_enabled
        @test enabled_score.nodes[1].p1 == 0.25
        @test enabled_score.scans[1].graph_logz == 2.5
        @test enabled_score.scans[1].meta_loglik ==
              enabled_score.scans[1].selected_node_loglik + 2.5
        graph_nodes = [
            node_for("g.sxm", "20240101", 1, 2.0),
            node_for("g.sxm", "20240101", 2, 2.0),
        ]
        chain_node1 = E.T12.ChainNodeResult(
            "g.sxm", "20240101", 1, 1.0, (log(0.75), log(0.25)),
            (0.25, 0.75), 2, "u",
        )
        chain_node2 = E.T12.ChainNodeResult(
            "g.sxm", "20240101", 2, 2.0, (log(0.60), log(0.40)),
            (0.40, 0.60), 2, "u",
        )
        block = E.T12.ChainBlockResult("s1", (chain_node1, chain_node2), (), 4.0, 0.0, (2, 2))
        reference = E.T12.ChainReferenceResult(
            "t11-C1", "C1", "partition", "outer_score", "20240101", "20240101",
            :PASS, :ok, (block,), 4.0, 0.0,
        )
        graph_report = (result_sha256="graph-result",)
        rebuilt = E._graph_from_reference(
            enabled_graphs["C1"], reference, graph_report, graph_nodes,
            E.EdgeNullEvidence[],
        )
        @test rebuilt.scan_logz["g.sxm"] == 4.0
        @test rebuilt.node_p1[("g.sxm", 1)] == 0.75
        @test rebuilt.node_p1[("g.sxm", 2)] == 0.60
        duplicate_reference = E.T12.ChainReferenceResult(
            "t11-C1", "C1", "partition", "outer_score", "20240101", "20240101",
            :PASS, :ok, (block, block), 8.0, 0.0,
        )
        duplicate_error = try
            E._graph_from_reference(enabled_graphs["C1"], duplicate_reference,
                                    graph_report, graph_nodes, E.EdgeNullEvidence[])
            nothing
        catch err
            err
        end
        @test duplicate_error isa E.EvaluatorError && duplicate_error.status == :BLOCKED
        missing_graph = Dict(
            "C1" => graph_for("C1", "outer-c1-20240101", "outer-c2-20240101";
                               enabled=true, logz=Dict("g.sxm" => 2.5)),
            "C2" => enabled_graphs["C2"],
        )
        missing = outer_for(1; outer_nodes=[node], graphs=missing_graph)
        missing_error = try E._score_fold(missing, "C1"); nothing catch err; err end
        @test missing_error isa E.EvaluatorError && missing_error.status == :BLOCKED
        fail_t11 = outer_for(1; outer_nodes=[node], graphs=Dict(
            "C1" => graph_for("C1", "outer-c1-20240101", "outer-c2-20240101"; t11=:FAIL),
            "C2" => graph_for("C2", "outer-c1-20240101", "outer-c2-20240101"),
        ))
        t11_error = try E._score_fold(fail_t11, "C1"); nothing catch err; err end
        @test t11_error isa E.EvaluatorError && t11_error.status == :FAIL
        fail_t12 = outer_for(1; outer_nodes=[node], graphs=Dict(
            "C1" => graph_for("C1", "outer-c1-20240101", "outer-c2-20240101"; t12=:FAIL),
            "C2" => graph_for("C2", "outer-c1-20240101", "outer-c2-20240101"),
        ))
        t12_error = try E._score_fold(fail_t12, "C1"); nothing catch err; err end
        @test t12_error isa E.EvaluatorError && t12_error.status == :FAIL

        no_scan = E.EvaluationInput(input.authority_sha256, input.universe_sha256,
                                    input.scan_dates[1:end-1], input.folds)
        substitute = E.EvaluationInput(input.authority_sha256, input.universe_sha256,
                                       [("other.sxm", SCAN_DATES[1]); input.scan_dates[2:end]],
                                       input.folds)
        duplicate = E.EvaluationInput(input.authority_sha256, input.universe_sha256,
                                      vcat(input.scan_dates, [first(input.scan_dates)]),
                                      input.folds)
        unsorted = E.EvaluationInput(input.authority_sha256, input.universe_sha256,
                                     reverse(input.scan_dates), input.folds)
        for malformed in (no_scan, substitute, duplicate, unsorted)
            assert_blocked(malformed)
        end
        omitted_inner = copy(input.folds)
        old = omitted_inner[1]
        omitted_inner[1] = E.OuterFoldEvidence(
            old.outer_date, old.training_dates, old.inner_folds[1:end-1],
            old.c1_status, old.c1_reason, old.c1_fit_sha256,
            old.c2_status, old.c2_reason, old.c2_fit_sha256,
            old.common_null_sha256, old.eligible_edge_count, old.nodes, old.edges, old.graphs,
        )
        assert_blocked(E.EvaluationInput(input.authority_sha256, input.universe_sha256,
                                         input.scan_dates, omitted_inner))
    end

    @testset "provenance, authority, firewall, and publication" begin
        input = happy_input()
        base = E.evaluate(input)
        view_changed = replace_outer_node(input, 1,
                                          node_for("scan_20240101.sxm", "20240101", 1, 2.0,
                                                   c1_gap=1.0))
        view_report = E.evaluate(view_changed)
        @test view_report.provenance_sha256 != base.provenance_sha256
        @test view_report.result_sha256 != base.result_sha256
        logz_changed = E.evaluate(happy_input(graph_logz=9.0))
        @test logz_changed.provenance_sha256 != base.provenance_sha256
        @test logz_changed.result_sha256 != base.result_sha256
        plain_base = E.report_files(base)
        plain_logz = E.report_files(logz_changed)
        for name in ARTIFACTS
            name == "receipt.toml" && continue
            @test plain_base[name] == plain_logz[name]
        end

        context, snapshot = snapshot_context(input, "context-pass")
        E._PRODUCTION_CONTEXT[input] = context
        bound = E.evaluate(input)
        @test bound.status == :PASS
        @test bound.source_sha256 == snapshot.sha256
        @test bound.provenance_sha256 == "context-pass"
        skipped_base = happy_input()
        skipped_old = skipped_base.folds[1]
        skipped_fold = E.OuterFoldEvidence(
            skipped_old.outer_date, skipped_old.training_dates, skipped_old.inner_folds,
            skipped_old.c1_status, skipped_old.c1_reason, skipped_old.c1_fit_sha256,
            skipped_old.c2_status, skipped_old.c2_reason, skipped_old.c2_fit_sha256,
            skipped_old.common_null_sha256, 1, skipped_old.nodes, skipped_old.edges,
            skipped_old.graphs,
        )
        skipped_folds = copy(skipped_base.folds)
        skipped_folds[1] = skipped_fold
        skipped_input = E.EvaluationInput(
            skipped_base.authority_sha256, skipped_base.universe_sha256,
            skipped_base.scan_dates, skipped_folds,
        )
        skipped_context, skipped_snapshot = snapshot_context(skipped_input, "context-skip")
        E._PRODUCTION_CONTEXT[skipped_input] = skipped_context
        skipped_report = E.evaluate(skipped_input)
        @test skipped_report.status == :SKIPPED
        @test skipped_report.reason == :common_null_unavailable
        @test skipped_report.source_sha256 == skipped_snapshot.sha256
        @test skipped_report.provenance_sha256 == "context-skip"
        fail_input = replace_inner_status(happy_input(), 1, :FAIL, :unary_numerical_failure)
        fail_context, fail_snapshot = snapshot_context(fail_input, "context-fail")
        E._PRODUCTION_CONTEXT[fail_input] = fail_context
        fail_report = E.evaluate(fail_input)
        @test fail_report.status == :FAIL
        @test fail_report.source_sha256 == fail_snapshot.sha256
        @test fail_report.provenance_sha256 == "context-fail"
        delete!(E._PRODUCTION_CONTEXT, input)
        delete!(E._PRODUCTION_CONTEXT, skipped_input)
        delete!(E._PRODUCTION_CONTEXT, fail_input)

        authority = E._authority_snapshots(pwd(), RUNTIME_EVALUATOR_CONFIG)
        paths = [snapshot.path for snapshot in authority.snapshots]
        @test any(endswith(path, "GateClosure.json") for path in paths)
        @test any(endswith(path, "AdversarialVerify.json") for path in paths)
        @test any(endswith(path, "closure-v1/DoneClaim.json") for path in paths)
        @test any(endswith(path, "closure-v1/review/AdversarialVerify.json") for path in paths)
        @test any(endswith(path, "closure-v1/root-manifest.sha256") for path in paths)
        @test any(endswith(path, "closure-v1-publication-receipt.json") for path in paths)
        @test !any(endswith(path, "boulder.json") for path in paths)
        @test !any(endswith(path, "structured-label-free-unit-assignment.md") for path in paths)
        @test any(endswith(path, "evaluate_structured_unit_assignment.jl") for path in paths)
        @test occursin(r"^[0-9a-f]{64}$", authority.authority_sha256)

        # Rebind checks use a complete isolated authority fixture.  The live
        # Plan and Boulder are deliberately changed in the fixture only; the
        # immutable closure authorities must still be the sole postclosure
        # lifecycle proof.
        authority_root = mktempdir()
        try
            fixture = isolated_authority_fixture(authority_root;
                                                 config_rel=RUNTIME_EVALUATOR_CONFIG)
            initial = E._authority_snapshots(authority_root,
                                              RUNTIME_EVALUATOR_CONFIG)
            closure_members = [snapshot for snapshot in initial.snapshots
                               if E._is_closure_path(snapshot.path)]
            @test length(closure_members) == 27
            @test length(initial.closure.snapshots) == 27
            @test Set(snapshot.path for snapshot in closure_members) ==
                  Set(snapshot.path for snapshot in initial.closure.snapshots)
            @test length(initial.inventories) == 3
            @test sort(collect(length(names) for (_, names) in initial.inventories)) == [2, 6, 21]
            @test Set(first.(initial.inventories)) ==
                  Set(relpath(path, authority_root) for path in
                      (fixture.closure_root, fixture.closure_guards, fixture.closure_review))
            @test all(UInt(stat(snapshot.path).mode) & UInt(0o222) == 0
                      for snapshot in closure_members)
            @test all(UInt(stat(joinpath(authority_root, directory)).mode) & UInt(0o222) == 0
                      for (directory, _) in initial.inventories)
            context = E._ProductionContext(authority_root, initial.snapshots,
                                            initial.inventories, "fixture")
            @test E._verify_context(context) === nothing
            write(fixture.plan, "later administrative Plan revision\n")
            write(fixture.boulder, "later administrative Boulder revision\n")
            later = E._authority_snapshots(authority_root,
                                           RUNTIME_EVALUATOR_CONFIG)
            @test later.authority_sha256 == initial.authority_sha256
            @test !any(endswith(path, "boulder.json") for path in
                       (snapshot.path for snapshot in later.snapshots))
            @test !any(endswith(path, "structured-label-free-unit-assignment.md") for path in
                       (snapshot.path for snapshot in later.snapshots))

            # Every immutable semantic control and sidecar is independently
            # tamper-bound.  All changes stay inside this isolated fixture.
            for path in (fixture.gate, fixture.prerequisite_review,
                         fixture.claim, fixture.claim_sidecar, fixture.review,
                         fixture.review_sidecar, fixture.root_manifest,
                         fixture.publication, fixture.publication_sidecar)
                original = Vector{UInt8}(read(path))
                rewrite_readonly(path, vcat(original, UInt8('\n')))
                @test authority_blocked(() ->
                    E._authority_snapshots(authority_root,
                                           RUNTIME_EVALUATOR_CONFIG))
                rewrite_readonly(path, original)
            end

            original_noncontrol = Vector{UInt8}(read(fixture.noncontrol))
            rewrite_readonly(fixture.noncontrol, vcat(original_noncontrol, UInt8('\n')))
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))
            rewrite_readonly(fixture.noncontrol, original_noncontrol)
            rewrite_readonly(fixture.noncontrol, vcat(original_noncontrol, UInt8('x')))
            @test authority_blocked(() -> E._verify_context(context))
            rewrite_readonly(fixture.noncontrol, original_noncontrol)
            @test E._verify_context(context) === nothing

            chmod(fixture.noncontrol, 0o644)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))
            @test authority_blocked(() -> E._verify_context(context))
            chmod(fixture.noncontrol, 0o444)
            chmod(fixture.closure_guards, 0o755)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))
            @test authority_blocked(() -> E._verify_context(context))
            chmod(fixture.closure_guards, 0o555)

            extra_file = joinpath(fixture.closure_root, "extra.txt")
            chmod(fixture.closure_root, 0o755)
            write(extra_file, "extra\n")
            chmod(extra_file, 0o444)
            chmod(fixture.closure_root, 0o555)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))
            remove_readonly(extra_file)

            extra_directory = joinpath(fixture.closure_root, "extra-directory")
            chmod(fixture.closure_root, 0o755)
            mkdir(extra_directory)
            chmod(extra_directory, 0o555)
            chmod(fixture.closure_root, 0o555)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))
            remove_readonly(extra_directory)

            symlink_target = joinpath(authority_root, "symlink-target")
            write(symlink_target, "target\n")
            original_noncontrol = Vector{UInt8}(read(fixture.noncontrol))
            remove_readonly(fixture.noncontrol)
            chmod(fixture.closure_root, 0o755)
            symlink(symlink_target, fixture.noncontrol)
            chmod(fixture.closure_root, 0o555)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))
            remove_readonly(fixture.noncontrol)
            restore_readonly(fixture.noncontrol, original_noncontrol)

            hardlink_target = fixture.claim
            original_noncontrol = Vector{UInt8}(read(fixture.noncontrol))
            remove_readonly(fixture.noncontrol)
            chmod(fixture.closure_root, 0o755)
            hardlink(hardlink_target, fixture.noncontrol)
            chmod(fixture.closure_root, 0o555)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))
            remove_readonly(fixture.noncontrol)
            restore_readonly(fixture.noncontrol, original_noncontrol)

            refreshed = E._authority_snapshots(authority_root,
                                               RUNTIME_EVALUATOR_CONFIG)
            context = E._ProductionContext(authority_root, refreshed.snapshots,
                                            refreshed.inventories, "fixture-refreshed")
            @test E._verify_context(context) === nothing
            inventory_extra = joinpath(fixture.closure_guards, "context-extra.txt")
            chmod(fixture.closure_guards, 0o755)
            write(inventory_extra, "context extra\n")
            chmod(inventory_extra, 0o444)
            chmod(fixture.closure_guards, 0o555)
            @test authority_blocked(() -> E._verify_context(context))
            remove_readonly(inventory_extra)
            @test E._verify_context(context) === nothing

            remove_readonly(fixture.noncontrol)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root,
                                       RUNTIME_EVALUATOR_CONFIG))

            history_root = mktempdir()
            try
                history_fixture = isolated_authority_fixture(history_root)
                seen = Dict{Tuple{UInt64,UInt64},String}()
                # Historical-only preclosure marker reconstruction; runtime
                # authority deliberately does not call _check_plan.
                plan_snapshot = E._check_plan(
                    history_root, "./.omo/plans/structured-label-free-unit-assignment.md", seen)
                @test plan_snapshot.sha256 == bytes2hex(sha256(read(history_fixture.plan)))
                write(history_fixture.plan, "malformed historical Plan\n")
                history_error = try
                    E._check_plan(history_root,
                                  "./.omo/plans/structured-label-free-unit-assignment.md",
                                  Dict{Tuple{UInt64,UInt64},String}())
                    nothing
                catch error
                    error
                end
                @test history_error isa E.EvaluatorError && history_error.status == :BLOCKED
            finally
                isdir(history_root) && make_tree_writable(history_root)
                rm(history_root; recursive=true, force=true)
            end
        finally
            isdir(authority_root) && make_tree_writable(authority_root)
            rm(authority_root; recursive=true, force=true)
        end

        snapshot_root = mktempdir()
        try
            local_context_value, local_paths = local_context(snapshot_root)
            @test E._verify_context(local_context_value) === nothing
            for path in (local_paths.gate, local_paths.review,
                         local_paths.boulder, local_paths.source)
                original = Vector{UInt8}(read(path))
                write(path, vcat(original, UInt8('x')))
                changed_error = try E._verify_context(local_context_value); nothing catch err; err end
                @test changed_error isa E.EvaluatorError
                @test changed_error.status == :BLOCKED
                write(path, original)
                @test E._verify_context(local_context_value) === nothing
            end
            extra = joinpath(local_paths.input_dir, "two.tsv")
            write(extra, "two\n")
            member_error = try E._verify_context(local_context_value); nothing catch err; err end
            @test member_error isa E.EvaluatorError && member_error.status == :BLOCKED
            rm(extra)
            @test E._verify_context(local_context_value) === nothing

            target = joinpath(snapshot_root, "symlink-target")
            write(target, "target")
            victim = local_paths.gate
            original = Vector{UInt8}(read(victim))
            rm(victim)
            symlink(target, victim)
            symlink_error = try E._verify_context(local_context_value); nothing catch err; err end
            @test symlink_error isa E.EvaluatorError && symlink_error.status == :BLOCKED
            rm(victim)
            write(victim, original)
            refreshed, _ = local_context(snapshot_root)
            @test E._verify_context(refreshed) === nothing
        finally
            rm(snapshot_root; recursive=true, force=true)
        end

        make_sensitive(parts) = join(parts, "")
        cli_cases = [
            ["positional"],
            [make_sensitive(["--", "mystery"])],
            ["--help", "--help"],
            [make_sensitive(["--", "syn", "thetic"])],
            [make_sensitive(["--", "t", "ru", "th"])],
        ]
        for args in cli_cases
            error = try E.parse_cli(args); nothing catch err; err end
            @test error isa E.EvaluatorError
            @test error.status == :BLOCKED
            @test error.reason == :cli_error
        end
        diagnostic = join(["forbidden_cli: forbidden ", join(["bench", "mark"], ""),
                           "/", join(["con", "trol"], ""), " CLI value"], "")
        value_cases = [
            ("--root", join(["/tmp/", "BENCH", "MARK", ".tsv"], "")),
            ("--features", join(["/tmp/", "Gr", "Ad", "E", ".tsv"], "")),
            ("--root", join(["/tmp/a/../", "Be", "Nc", "Hm", "Ar", "K", ".tsv"], "")),
            ("--features", join(["/tmp/nested/", "Co", "Nt", "Ro", "L", ".tsv"], "")),
        ]
        for (flag, value) in value_cases
            error = try E.parse_cli([flag, value]); nothing catch err; err end
            @test error isa E.EvaluatorError
            @test error.message == diagnostic
        end

        files = E.report_files(base)
        root = mktempdir()
        try
            @test E._publish_atomic(root, "out", files) == :committed_verified
            @test E._publish_atomic(root, "out", files) == :committed_verified
            @test readdir(root) == ["out"]
            different = copy(files)
            different["summary.tsv"] = vcat(different["summary.tsv"], UInt8('x'))
            collision = try E._publish_atomic(root, "out", different); nothing catch err; err end
            @test collision isa E.EvaluatorError
            @test collision.reason == :publication_collision
            @test readdir(root) == ["out"]
            rename_stage = mktempdir(root; prefix="rename-check-", cleanup=false)
            rename_error = try
                E._rename_noreplace(rename_stage, joinpath(root, "out"))
                nothing
            catch err
                err
            end
            @test rename_error isa E._PublicationError
            @test rename_error.state == :not_committed
            ispath(rename_stage) && rm(rename_stage; recursive=true, force=true)
        finally
            rm(root; recursive=true, force=true)
        end

        pre_root = mktempdir()
        try
            pre_context, pre_paths = local_context(pre_root)
            @test E._publish_atomic(pre_root, "verified", files;
                                    context=pre_context) == :committed_verified
            @test E._publish_atomic(pre_root, "verified", files;
                                    context=pre_context) == :committed_verified
            write(pre_paths.gate, "mutated")
            pre_error = try
                E._publish_atomic(pre_root, "out", files; context=pre_context)
                nothing
            catch err
                err
            end
            @test pre_error isa E._PublicationError
            @test pre_error.state == :not_committed
            @test !ispath(joinpath(pre_root, "out"))
            @test all(!startswith(name, ".stmfit-t13-evaluator-phase1-")
                      for name in readdir(pre_root))
        finally
            rm(pre_root; recursive=true, force=true)
        end
        @test E._publication_error(:ambiguous, :publication_failure, "x").state == :ambiguous
        @test E._publication_error(:not_committed, :publication_failure, "x").state == :not_committed
        @test_throws ArgumentError E._publication_error(:committed_verified,
                                                        :publication_failure, "x")

        @testset "receipt-gated renameat2 fallback" begin
            with_hook(f, hook) = begin
                prior = E._PUBLICATION_HOOK[]
                E._PUBLICATION_HOOK[] = hook
                try
                    f()
                finally
                    E._PUBLICATION_HOOK[] = prior
                end
            end
            force_errno(errno) = (event, args...) ->
                event == :rename_noreplace ? errno : nothing

            fallback_root = mktempdir()
            try
                @test with_hook(force_errno(22)) do
                    E._publish_atomic(fallback_root, "normal", files)
                end == :committed_verified
                @test E._publish_atomic(fallback_root, "normal", files) == :committed_verified
                @test E._same_publication(joinpath(fallback_root, "normal"), files)
                @test readdir(fallback_root) == ["normal"]
                different = copy(files)
                different["summary.tsv"] = vcat(different["summary.tsv"], UInt8('x'))
                collision = try
                    E._publish_atomic(fallback_root, "normal", different)
                    nothing
                catch caught
                    caught
                end
                @test collision isa E.EvaluatorError
                @test collision.reason == :publication_collision
            finally
                rm(fallback_root; recursive=true, force=true)
            end

            for errno in (22, 38, 95)
                root = mktempdir()
                try
                    @test with_hook(force_errno(errno)) do
                        E._publish_atomic(root, "out", files)
                    end == :committed_verified
                    @test E._same_publication(joinpath(root, "out"), files)
                finally
                    rm(root; recursive=true, force=true)
                end
            end

            blocker_files = E._blocker_files(E.EvaluatorError(
                :BLOCKED, :cli_error, "blocked"))
            blocker_root = mktempdir()
            try
                @test with_hook(force_errno(22)) do
                    E._publish_atomic(blocker_root, "blocker", blocker_files)
                end == :committed_verified
                @test E._same_publication(joinpath(blocker_root, "blocker"), blocker_files)
            finally
                rm(blocker_root; recursive=true, force=true)
            end

            for errno in (18, 5, 1, 13, 9, 45)
                root = mktempdir()
                try
                    error = with_hook(force_errno(errno)) do
                        try
                            E._publish_atomic(root, "out", files)
                            nothing
                        catch caught
                            caught
                        end
                    end
                    @test error isa E._PublicationError
                    @test error.state == :not_committed
                    @test !ispath(joinpath(root, "out"))
                    @test all(!startswith(name, ".stmfit-t13-evaluator-phase1-")
                              for name in readdir(root))
                finally
                    rm(root; recursive=true, force=true)
                end
            end

            for failure_event in (:before_fsync_file, :before_fsync_directory)
                root = mktempdir()
                try
                    error = with_hook((event, args...) -> begin
                        event == :rename_noreplace && return 22
                        event == failure_event &&
                            throw(E._publication_error(
                                :not_committed, :injected_fsync_failure,
                                String(failure_event)))
                        nothing
                    end) do
                        try
                            E._publish_atomic(root, "out", files)
                            nothing
                        catch caught
                            caught
                        end
                    end
                    @test error isa E._PublicationError
                    @test error.state == :not_committed
                    @test !ispath(joinpath(root, "out"))
                finally
                    make_tree_writable(root)
                    rm(root; recursive=true, force=true)
                end
            end

            root = mktempdir()
            try
                observed_mode = Ref{UInt}(0)
                result = with_hook((event, args...) -> begin
                    event == :rename_noreplace && return 22
                    if event == :after_reserve_create
                        observed_mode[] = UInt(stat(args[2]).mode) & UInt(0o777)
                    end
                    nothing
                end) do
                    E._publish_atomic(root, "mode", files)
                end
                @test result == :committed_verified
                @test observed_mode[] == UInt(0o700)
            finally
                rm(root; recursive=true, force=true)
            end

            for (target_kind, expected_state) in ((:payload, :not_committed),
                                                   (:receipt, :ambiguous))
                root = mktempdir()
                try
                    destination_bytes = target_kind == :payload ?
                        Vector{UInt8}(codeunits("existing payload\n")) : files["receipt.toml"]
                    error = with_hook((event, args...) -> begin
                        event == :rename_noreplace && return 22
                        if event == :before_hardlink && basename(args[2]) ==
                           (target_kind == :payload ? "bootstrap.tsv" : "receipt.toml")
                            write(args[2], destination_bytes)
                        end
                        nothing
                    end) do
                        try
                            E._publish_atomic(root, String(target_kind), files)
                            nothing
                        catch caught
                            caught
                        end
                    end
                    @test error isa E._PublicationError
                    @test error.state == expected_state
                    target = joinpath(root, String(target_kind),
                                      target_kind == :payload ? "bootstrap.tsv" : "receipt.toml")
                    @test read(target) == destination_bytes
                finally
                    make_tree_writable(root)
                    rm(root; recursive=true, force=true)
                end
            end

            for (mutation, expected_state) in ((:payload, :not_committed),
                                                (:receipt, :not_committed))
                root = mktempdir()
                try
                    error = with_hook((event, args...) -> begin
                        event == :rename_noreplace && return 22
                        if event == :before_hardlink
                            name = basename(args[1])
                            if mutation == :payload && name == "bootstrap.tsv"
                                write(args[1], "mutated payload\n")
                            end
                        elseif mutation == :receipt && event == :before_receipt_link
                            write(joinpath(args[1], "receipt.toml"), "mutated receipt\n")
                        end
                        nothing
                    end) do
                        try
                            E._publish_atomic(root, String(mutation), files)
                            nothing
                        catch caught
                            caught
                        end
                    end
                    @test error isa E._PublicationError
                    @test error.state == expected_state
                    @test !ispath(joinpath(root, String(mutation), "receipt.toml"))
                finally
                    make_tree_writable(root)
                    rm(root; recursive=true, force=true)
                end
            end

            root = mktempdir()
            try
                error = with_hook((event, args...) -> begin
                    event == :rename_noreplace && return 22
                    if event == :after_receipt_link
                        error("receipt visible")
                    end
                    nothing
                end) do
                    try
                        E._publish_atomic(root, "receipt-visible", files)
                        nothing
                    catch caught
                        caught
                    end
                end
                @test error isa E._PublicationError
                @test error.state == :ambiguous
                @test isfile(joinpath(root, "receipt-visible", "receipt.toml"))
            finally
                make_tree_writable(root)
                rm(root; recursive=true, force=true)
            end

            for target in (:parent_after_reserve, :destination_precommit,
                           :stage_precommit, :destination_post_receipt,
                           :parent_post_receipt)
                root = mktempdir()
                try
                    error = with_hook((event, args...) -> begin
                        event == :rename_noreplace && return 22
                        if event == :before_fsync_directory
                            path, state = args[1], args[2]
                            destination = joinpath(root, "out")
                            should_fail = target == :parent_after_reserve &&
                                path == root && state == :not_committed
                            should_fail |= target == :destination_precommit &&
                                path == destination && state == :not_committed
                            should_fail |= target == :stage_precommit &&
                                startswith(path, root * "/.stmfit-t13-evaluator-phase1-") &&
                                state == :not_committed && isdir(destination)
                            should_fail |= target == :destination_post_receipt &&
                                path == destination && state == :ambiguous
                            should_fail |= target == :parent_post_receipt &&
                                path == root && state == :ambiguous
                            should_fail && throw(E._publication_error(
                                state, :injected_fsync_failure, String(target)))
                        end
                        nothing
                    end) do
                        try
                            E._publish_atomic(root, "out", files)
                            nothing
                        catch caught
                            caught
                        end
                    end
                    @test error isa E._PublicationError
                    expected = target in (:destination_post_receipt, :parent_post_receipt) ?
                               :ambiguous : :not_committed
                    @test error.state == expected
                    @test (expected == :ambiguous) ==
                          isfile(joinpath(root, "out", "receipt.toml"))
                finally
                    make_tree_writable(root)
                    rm(root; recursive=true, force=true)
                end
            end

            root = mktempdir()
            try
                preferred_state = Ref{Symbol}(:not_committed)
                error = with_hook((event, args...) -> begin
                    if event == :after_preferred_rename
                        preferred_state[] = args[2]
                        throw(E._publication_error(:ambiguous, :injected_failure,
                                                    "preferred receipt visible"))
                    end
                    nothing
                end) do
                    try
                        E._publish_atomic(root, "preferred-visible", files)
                        nothing
                    catch caught
                        caught
                    end
                end
                @test error isa E._PublicationError
                @test error.state == :ambiguous
                @test preferred_state[] == :ambiguous
                @test isfile(joinpath(root, "preferred-visible", "receipt.toml"))
            finally
                make_tree_writable(root)
                rm(root; recursive=true, force=true)
            end

            context_root = mktempdir()
            publication_root = mktempdir()
            try
                context, context_paths = local_context(context_root)
                for (event_to_mutate, expected_state) in
                    ((:before_receipt_link, :not_committed),
                     (:after_receipt_link, :ambiguous))
                    hook = (event, args...) -> begin
                        event == :rename_noreplace && return 22
                        event == event_to_mutate && write(context_paths.gate, "mutated")
                        nothing
                    end
                    error = with_hook(hook) do
                        try
                            E._publish_atomic(publication_root, String(event_to_mutate), files;
                                               context=context)
                            nothing
                        catch caught
                            caught
                        end
                    end
                    @test error isa E._PublicationError
                    @test error.state == expected_state
                    write(context_paths.gate, "gate")
                    context, context_paths = local_context(context_root)
                end
            finally
                make_tree_writable(context_root)
                rm(context_root; recursive=true, force=true)
                make_tree_writable(publication_root)
                rm(publication_root; recursive=true, force=true)
            end

            for residue in (:empty, :payload)
                root = mktempdir()
                try
                    destination = joinpath(root, "out")
                    mkdir(destination)
                    chmod(destination, 0o700)
                    residue == :payload && write(joinpath(destination, "summary.tsv"), "x")
                    collision = try
                        E._publish_atomic(root, "out", files)
                        nothing
                    catch caught
                        caught
                    end
                    @test collision isa E.EvaluatorError
                    @test collision.reason == :publication_collision
                    @test isdir(destination) && !ispath(joinpath(destination, "receipt.toml"))
                finally
                    make_tree_writable(root)
                    rm(root; recursive=true, force=true)
                end
            end

            failure_events = (:after_reserve, :after_payload_link,
                              :before_receipt_link, :after_receipt_link,
                              :before_final_verify)
            for event_to_fail in failure_events
                root = mktempdir()
                try
                    error = with_hook((event, args...) -> begin
                        event == :rename_noreplace && return 22
                        event == event_to_fail &&
                            throw(E._publication_error(
                                event_to_fail in (:after_receipt_link, :before_final_verify) ?
                                :ambiguous : :not_committed,
                                :injected_failure, String(event_to_fail)))
                        nothing
                    end) do
                        try
                            E._publish_atomic(root, "out", files;
                                               context=nothing)
                            nothing
                        catch caught
                            caught
                        end
                    end
                    @test error isa E._PublicationError
                    expected_state = event_to_fail in (:after_receipt_link, :before_final_verify) ?
                                     :ambiguous : :not_committed
                    @test error.state == expected_state
                    if expected_state == :ambiguous
                        @test isdir(joinpath(root, "out"))
                    end
                finally
                    make_tree_writable(root)
                    rm(root; recursive=true, force=true)
                end
            end

            root = mktempdir()
            try
                race = with_hook((event, args...) -> begin
                    event == :rename_noreplace && return 22
                    if event == :before_reserve
                        destination = args[2]
                        mkdir(destination)
                        chmod(destination, 0o700)
                    end
                    nothing
                end) do
                    try
                        E._publish_atomic(root, "out", files)
                        nothing
                    catch caught
                        caught
                    end
                end
                @test race isa E.EvaluatorError
                @test race.reason == :publication_collision
            finally
                make_tree_writable(root)
                rm(root; recursive=true, force=true)
            end

            root = mktempdir()
            try
                hook = (event, args...) -> begin
                    event == :rename_noreplace && return 18
                    if event == :before_stage_cleanup
                        stage = args[1]
                        ispath(stage) && rm(stage; recursive=true, force=true)
                        mkdir(stage)
                    end
                    nothing
                end
                error = with_hook(hook) do
                    try
                        E._publish_atomic(root, "out", files)
                        nothing
                    catch caught
                        caught
                    end
                end
                @test error isa E._PublicationError
                @test any(startswith(name, ".stmfit-t13-evaluator-phase1-")
                          for name in readdir(root))
            finally
                make_tree_writable(root)
                rm(root; recursive=true, force=true)
            end

            root = mktempdir()
            try
                @test with_hook(force_errno(22)) do
                    E._publish_atomic(root, "out", files)
                end == :committed_verified
                destination = joinpath(root, "out")
                chmod(joinpath(destination, "summary.tsv"), 0o600)
                @test !E._same_publication(destination, files)
                chmod(joinpath(destination, "summary.tsv"), 0o644)
                hardlink(joinpath(destination, "summary.tsv"),
                         joinpath(destination, "hardlinked.tsv"))
                @test !E._same_publication(destination, files)
                rm(joinpath(destination, "hardlinked.tsv"))
                rm(joinpath(destination, "receipt.toml"))
                @test !E._same_publication(destination, files)
                write(joinpath(destination, "receipt.toml"), files["receipt.toml"])
                chmod(joinpath(destination, "receipt.toml"), 0o644)
                symlink(joinpath(destination, "summary.tsv"),
                        joinpath(destination, "linked.tsv"))
                @test !E._same_publication(destination, files)
            finally
                make_tree_writable(root)
                rm(root; recursive=true, force=true)
            end
        end

        blocker_error = E.EvaluatorError(:BLOCKED, :cli_error, "blocked")
        blocker = E._blocker_files(blocker_error)
        @test collect(keys(blocker)) == ["receipt.toml"]
        blocker_receipt = TOML.parse(String(blocker["receipt.toml"]))
        @test blocker_receipt["schema"] ==
              "structured_label_free_unit_assignment_evaluator_v1_blocker_receipt_v2"
        @test blocker_receipt["schema_version"] == 2
        @test !haskey(blocker_receipt, "plan_sha256")
        @test !haskey(blocker_receipt, "boulder_sha256")
        @test blocker_receipt["preclosure_plan_sha256"] == E._PLAN_SHA256
        @test blocker_receipt["preclosure_boulder_sha256"] == E._BOULDER_SHA256
        @test blocker_receipt["evaluator_config_sha256"] ==
              "ec0546096b3c4742cd86d8c3d40788a5894d20fc5a4702b0318581f10f8b0b90"
        @test blocker_receipt["live_plan_boulder_runtime_authority"] === false
        @test blocker_receipt["closure_v1_claim_sha256"] == E._CLOSURE_V1_CLAIM_SHA256
        @test blocker_receipt["closure_v1_review_sha256"] == E._CLOSURE_V1_REVIEW_SHA256
        @test blocker_receipt["closure_v1_root_manifest_sha256"] ==
              E._CLOSURE_V1_ROOT_MANIFEST_SHA256
        @test blocker_receipt["closure_v1_publication_receipt_sha256"] ==
              E._CLOSURE_V1_PUBLICATION_RECEIPT_SHA256
        @test blocker_receipt["closure_v1_publication_receipt_sidecar_sha256"] ==
              E._CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_SHA256
        @test blocker_receipt["status"] == "BLOCKED"
        @test blocker_receipt["work_rows"] == 0
        @test blocker_receipt["bootstrap_rows"] == 0
        @test blocker_receipt["sign_rows"] == 0
        @test !haskey(blocker_receipt, "files")

        function config_context(config_rel::String)
            path = normpath(joinpath(pwd(), config_rel))
            bytes = Vector{UInt8}(read(path))
            info = stat(path)
            snapshot = E._Snapshot(path, bytes, bytes2hex(sha256(bytes)),
                                   UInt64(info.device), UInt64(info.inode),
                                   UInt64(info.nlink))
            return E._ProductionContext(pwd(), [snapshot],
                                        Tuple{String,Vector{String}}[], "config-context")
        end
        runtime_blocker = TOML.parse(String(E._blocker_files(
            blocker_error; context=config_context(RUNTIME_EVALUATOR_CONFIG))["receipt.toml"]))
        @test runtime_blocker["evaluator_config_sha256"] ==
              "4a355c7ec5572976f5f54eaf14f2a01e8c91a9304d789e60b844388ac20939a4"
        historical_context_blocker = TOML.parse(String(E._blocker_files(
            blocker_error; context=config_context(HISTORICAL_EVALUATOR_CONFIG))["receipt.toml"]))
        @test historical_context_blocker["evaluator_config_sha256"] ==
              "ec0546096b3c4742cd86d8c3d40788a5894d20fc5a4702b0318581f10f8b0b90"
        runtime_snapshot = only(config_context(RUNTIME_EVALUATOR_CONFIG).snapshots)
        tampered_snapshot = E._Snapshot(
            runtime_snapshot.path, runtime_snapshot.bytes, "0"^64,
            runtime_snapshot.device, runtime_snapshot.inode, runtime_snapshot.nlink)
        @test_throws E.EvaluatorError E._blocker_files(
            blocker_error; context=E._ProductionContext(
                pwd(), [tampered_snapshot], Tuple{String,Vector{String}}[], "tampered-config"))
        @test_throws E.EvaluatorError E._blocker_files(
            blocker_error; context=E._ProductionContext(
                pwd(), E._Snapshot[], Tuple{String,Vector{String}}[], "missing-config"))

        @test [event.ordinal for event in base.events] == collect(1:length(base.events))
        @test base.events[end].consumed == false
        @test base.events[end].terminal_status == :PASS
        unexpected_input = happy_input()
        bad_context = E._ProductionContext(pwd(), E._Snapshot[],
                                            Tuple{String,Vector{String}}[], "bad")
        E._PRODUCTION_CONTEXT[unexpected_input] = bad_context
        @test_throws ArgumentError E.evaluate(unexpected_input)
        delete!(E._PRODUCTION_CONTEXT, unexpected_input)
    end

    @testset "T13 runtime authority correction rebind" begin
        historical_bytes = Vector{UInt8}(read(joinpath(pwd(), HISTORICAL_EVALUATOR_CONFIG)))
        runtime_bytes = Vector{UInt8}(read(joinpath(pwd(), RUNTIME_EVALUATOR_CONFIG)))
        @test bytes2hex(sha256(historical_bytes)) == E._CONFIG_SHA256
        @test bytes2hex(sha256(runtime_bytes)) == E._RUNTIME_CONFIG_SHA256
        @test E._authority_config_kind(pwd(), HISTORICAL_EVALUATOR_CONFIG) == :historical
        @test E._authority_config_kind(pwd(), RUNTIME_EVALUATOR_CONFIG) == :runtime
        historical_document = TOML.parse(String(historical_bytes))
        runtime_document = TOML.parse(String(runtime_bytes))
        for section in setdiff(Set(String.(keys(historical_document))), Set(["authority"]))
            @test runtime_document[section] == historical_document[section]
        end
        for key in keys(historical_document["authority"])
            key == "t12_edge_model_sha256" && continue
            @test runtime_document["authority"][key] == historical_document["authority"][key]
        end
        @test runtime_document["authority"]["t12_edge_model_sha256"] ==
              E._RUNTIME_T12_EDGE_MODEL_SHA256

        for other in ("config/unit_assignment_structured_evaluator_runtime-copy.toml",
                      "config/not-an-evaluator.toml")
            @test authority_blocked(() -> E._authority_snapshots(pwd(), other))
        end

        authority_root = mktempdir()
        try
            fixture = isolated_authority_fixture(authority_root;
                                                 config_rel=RUNTIME_EVALUATOR_CONFIG)
            authority = E._authority_snapshots(authority_root, RUNTIME_EVALUATOR_CONFIG)
            required = [
                joinpath(authority_root, RUNTIME_EVALUATOR_CONFIG),
                joinpath(authority_root, E._RUNTIME_T12_ROOT_MANIFEST_PATH),
                joinpath(authority_root, E._RUNTIME_T12_CLAIM_PATH),
                joinpath(authority_root, E._RUNTIME_T12_REVIEW_PATH),
                joinpath(authority_root, E._RUNTIME_T12_PUBLICATION_PATH),
                joinpath(authority_root, "test/evaluate_structured_unit_assignment.jl"),
            ]
            snapshot_paths = [snapshot.path for snapshot in authority.snapshots]
            @test authority.config_kind == :runtime
            @test all(count(==(path), snapshot_paths) == 1 for path in required)
            @test all(snapshot.nlink == 1 for snapshot in authority.snapshots)
            @test length(unique((snapshot.device, snapshot.inode)
                               for snapshot in authority.snapshots)) ==
                  length(authority.snapshots)
            @test length(authority.authority_sha256) == 64
            @test !isempty(fixture.config)

            correction_paths = [
                joinpath(authority_root, E._RUNTIME_T12_ROOT_MANIFEST_PATH),
                joinpath(authority_root, E._RUNTIME_T12_CLAIM_PATH),
                joinpath(authority_root, E._RUNTIME_T12_REVIEW_PATH),
                joinpath(authority_root, E._RUNTIME_T12_PUBLICATION_PATH),
            ]
            for path in correction_paths
                original = Vector{UInt8}(read(path))
                rewrite_readonly(path, vcat(original, UInt8('\n')))
                @test authority_blocked(() ->
                    E._authority_snapshots(authority_root, RUNTIME_EVALUATOR_CONFIG))
                rewrite_readonly(path, original)
            end

            missing_path = first(correction_paths)
            missing_bytes = Vector{UInt8}(read(missing_path))
            remove_readonly(missing_path)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_EVALUATOR_CONFIG))
            write(missing_path, missing_bytes)
            chmod(missing_path, 0o644)

            config_path = fixture.config
            original_config = Vector{UInt8}(read(config_path))
            mutated_config = replace(String(copy(original_config)),
                E._RUNTIME_T12_EDGE_MODEL_SHA256 => "0"^64; count=1)
            rewrite_readonly(config_path, Vector{UInt8}(codeunits(mutated_config)))
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_EVALUATOR_CONFIG))
            rewrite_readonly(config_path, original_config)

            initial = E._authority_snapshots(authority_root, RUNTIME_EVALUATOR_CONFIG)
            semantic_error(path::String, bytes::Vector{UInt8}) = begin
                snapshots = copy(initial.snapshots)
                index = findfirst(snapshot -> snapshot.path == path, snapshots)
                index === nothing && error("semantic test snapshot is absent")
                old = snapshots[index]
                snapshots[index] = E._Snapshot(old.path, bytes, old.sha256,
                                               old.device, old.inode, old.nlink)
                try
                    E._validate_runtime_correction_semantics(authority_root, snapshots, :runtime)
                    nothing
                catch error
                    error
                end
            end
            claim_path = joinpath(authority_root, E._RUNTIME_T12_CLAIM_PATH)
            review_path = joinpath(authority_root, E._RUNTIME_T12_REVIEW_PATH)
            publication_path = joinpath(authority_root, E._RUNTIME_T12_PUBLICATION_PATH)
            root_manifest_path = joinpath(authority_root, E._RUNTIME_T12_ROOT_MANIFEST_PATH)
            for (path, needle, replacement) in (
                (claim_path, "\"status\": \"PASS\"", "\"status\": \"FAIL\""),
                (review_path, "\"verdict\": \"PASS\"", "\"verdict\": \"FAIL\""),
                (publication_path, "\"publication_state\": \"committed_verified\"",
                 "\"publication_state\": \"staged\""),
            )
                mutated = replace(String(read(path)), needle => replacement; count=1)
                error = semantic_error(path, Vector{UInt8}(codeunits(mutated)))
                @test error isa E.EvaluatorError
                @test error.status == :BLOCKED
                @test error.reason == :authority_binding_mismatch
            end
            manifest_error = semantic_error(
                root_manifest_path,
                vcat(Vector{UInt8}(read(root_manifest_path)),
                     Vector{UInt8}(codeunits("malformed\n"))))
            @test manifest_error isa E.EvaluatorError
            @test manifest_error.status == :BLOCKED
            @test manifest_error.reason == :authority_bundle_mismatch

            symlink_target = joinpath(authority_root, "correction-symlink-target")
            write(symlink_target, "substitution\n")
            victim = claim_path
            original = Vector{UInt8}(read(victim))
            remove_readonly(victim)
            symlink(symlink_target, victim)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_EVALUATOR_CONFIG))
            remove_readonly(victim)
            write(victim, original)
            chmod(victim, 0o644)

            victim = review_path
            original = Vector{UInt8}(read(victim))
            remove_readonly(victim)
            hardlink(claim_path, victim)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_EVALUATOR_CONFIG))
            remove_readonly(victim)
            write(victim, original)
            chmod(victim, 0o644)
        finally
            make_tree_writable(authority_root)
            rm(authority_root; recursive=true, force=true)
        end

        production_options(config::String) = Dict{String,String}(
            "--root" => pwd(), "--evaluator-config" => config,
            "--features" => "missing/features.tsv",
            "--candidate-config" => "missing/candidate.toml",
            "--model-config" => "missing/model.toml",
            "--universe-dir" => "missing/universe", "--edge-dir" => "missing/edge",
            "--forward-receipt" => "missing/forward.toml",
            "--backward-receipt" => "missing/backward.toml",
            "--admission-dir" => "missing/admission")
        if VERSION == v"1.12.6"
            historical_error = try
                E._load_production_input(production_options(HISTORICAL_EVALUATOR_CONFIG))
                nothing
            catch error
                error
            end
            @test historical_error isa E.EvaluatorError
            @test historical_error.reason == :authority_path_mismatch
        else
            @test VERSION != v"1.12.6"
        end
    end

    @testset "T13 runtime authority correction rebind v2" begin
        runtime_v2_bytes = Vector{UInt8}(read(joinpath(pwd(), RUNTIME_V2_EVALUATOR_CONFIG)))
        @test bytes2hex(sha256(runtime_v2_bytes)) == E._RUNTIME_V2_CONFIG_SHA256
        @test E._authority_config_kind(pwd(), RUNTIME_V2_EVALUATOR_CONFIG) == :runtime_v2
        authority = E._authority_snapshots(pwd(), RUNTIME_V2_EVALUATOR_CONFIG)
        @test authority.config_kind == :runtime_v2
        @test authority.config_sha256 == E._RUNTIME_V2_CONFIG_SHA256
        @test length(authority.authority_sha256) == 64
        contract = E._runtime_contract(:runtime_v2)
        @test contract.edge_model_sha256 ==
              "5e3b1c0371bb0387024714d420d25a96a9d076aa3f0ba6637a10899837a6c8b1"
        @test contract.root_manifest_sha256 ==
              "d079f2cee7e0558a93b1bfda156c979e86deb988beb3ef0cb553759fe4bb9626"

        runtime_v2_path = joinpath(pwd(), RUNTIME_V2_EVALUATOR_CONFIG)
        runtime_v2_bytes_for_context = Vector{UInt8}(read(runtime_v2_path))
        runtime_v2_info = stat(runtime_v2_path)
        runtime_v2_snapshot = E._Snapshot(
            runtime_v2_path, runtime_v2_bytes_for_context,
            bytes2hex(sha256(runtime_v2_bytes_for_context)),
            UInt64(runtime_v2_info.device), UInt64(runtime_v2_info.inode),
            UInt64(runtime_v2_info.nlink))
        runtime_v2_context = E._ProductionContext(
            pwd(), [runtime_v2_snapshot], Tuple{String,Vector{String}}[], "runtime-v2-test")
        blocker = E.EvaluatorError(:BLOCKED, :runtime_v2_test, "runtime-v2")
        blocker_receipt = TOML.parse(String(E._blocker_files(
            blocker; context=runtime_v2_context)["receipt.toml"]))
        @test blocker_receipt["evaluator_config_sha256"] == E._RUNTIME_V2_CONFIG_SHA256

        authority_root = mktempdir()
        try
            isolated_authority_fixture(authority_root; config_rel=RUNTIME_V2_EVALUATOR_CONFIG)
            initial = E._authority_snapshots(authority_root, RUNTIME_V2_EVALUATOR_CONFIG)
            paths = [snapshot.path for snapshot in initial.snapshots]
            for relative in (contract.root_manifest_path, contract.claim_path,
                             contract.review_path, contract.publication_path)
                @test count(==(joinpath(authority_root, relative)), paths) == 1
            end
            for relative in (contract.root_manifest_path, contract.claim_path,
                             contract.review_path, contract.publication_path)
                path = joinpath(authority_root, relative)
                original = Vector{UInt8}(read(path))
                rewrite_readonly(path, vcat(original, UInt8('\n')))
                @test authority_blocked(() ->
                    E._authority_snapshots(authority_root, RUNTIME_V2_EVALUATOR_CONFIG))
                rewrite_readonly(path, original)
            end
        finally
            make_tree_writable(authority_root)
            rm(authority_root; recursive=true, force=true)
        end

        if VERSION == v"1.12.6"
            v1_options = Dict{String,String}(
                "--root" => pwd(), "--evaluator-config" => RUNTIME_EVALUATOR_CONFIG,
                "--features" => "missing/features.tsv",
                "--candidate-config" => "missing/candidate.toml",
                "--model-config" => "missing/model.toml",
                "--universe-dir" => "missing/universe", "--edge-dir" => "missing/edge",
                "--forward-receipt" => "missing/forward.toml",
                "--backward-receipt" => "missing/backward.toml",
                "--admission-dir" => "missing/admission")
            v1_error = try
                E._load_production_input(v1_options)
                nothing
            catch error
                error
            end
            @test v1_error isa E.EvaluatorError
            @test v1_error.reason == :authority_path_mismatch
        end
    end

    @testset "T13 runtime authority correction rebind v3" begin
        v2_document = TOML.parse(read(joinpath(pwd(), RUNTIME_V2_EVALUATOR_CONFIG), String))
        v3_path = joinpath(pwd(), RUNTIME_V3_EVALUATOR_CONFIG)
        v3_bytes = Vector{UInt8}(read(v3_path))
        v3_document = TOML.parse(String(copy(v3_bytes)))
        changed = Set{String}()
        for section in union(Set(String.(keys(v2_document))),
                             Set(String.(keys(v3_document))))
            old_table = get(v2_document, section, nothing)
            new_table = get(v3_document, section, nothing)
            old_table isa AbstractDict && new_table isa AbstractDict || begin
                old_table == new_table || push!(changed, section)
                continue
            end
            for key in union(Set(String.(keys(old_table))),
                             Set(String.(keys(new_table))))
                old_value = get(old_table, key, nothing)
                new_value = get(new_table, key, nothing)
                old_value == new_value || push!(changed, "$section.$key")
            end
        end
        allowed = Set([
            "authority.t12_edge_model_sha256",
            "authority.t12_chain_sha256",
            "authority.t12_test_sha256",
            "authority.t12_authority_root_manifest_path",
            "authority.t12_authority_root_manifest_sha256",
            "authority.t12_authority_claim_path",
            "authority.t12_authority_claim_sha256",
            "authority.t12_authority_review_path",
            "authority.t12_authority_review_sha256",
            "authority.t12_authority_publication_receipt_path",
            "authority.t12_authority_publication_receipt_sha256",
        ])
        @test changed == allowed
        @test bytes2hex(sha256(v3_bytes)) == E._RUNTIME_V3_CONFIG_SHA256
        @test E._authority_config_kind(pwd(), RUNTIME_V3_EVALUATOR_CONFIG) == :runtime_v3
        @test E._RUNTIME_V3_AUTHORITY_SCHEMA == "schema=structured-evaluator-authority-v5"

        authority = E._authority_snapshots(pwd(), RUNTIME_V3_EVALUATOR_CONFIG)
        @test authority.config_kind == :runtime_v3
        @test authority.config_sha256 == E._RUNTIME_V3_CONFIG_SHA256
        @test occursin(r"^[0-9a-f]{64}$", authority.authority_sha256)
        contract = E._runtime_contract(:runtime_v3)
        @test contract.edge_model_sha256 == E._RUNTIME_V3_T12_EDGE_MODEL_SHA256
        @test contract.chain_sha256 == E._RUNTIME_V3_T12_CHAIN_SHA256
        @test contract.corrected_test_sha256 == E._RUNTIME_V3_T12_CORRECTED_TEST_SHA256
        @test contract.root_manifest_sha256 ==
              "45f846e35f2965b50d13f4242add20a7d3e4e25e9ada94e619fa12003ed3b03f"
        @test contract.manifest_entries == 92
        @test contract.edge_model_product_path ==
              "products/v3/edge_model.jl"
        @test contract.chain_product_path == "products/v3/chain_inference.jl"
        @test contract.test_product_path ==
              "products/v3/test_structured_chain_inference.jl"
        @test authority.config_sha256 ==
              E._context_evaluator_config_sha256(E._ProductionContext(
                  pwd(), [authority.snapshots[1]], Tuple{String,Vector{String}}[], "v3"))

        for (config, expected_kind) in ((RUNTIME_EVALUATOR_CONFIG, :runtime),
                                        (RUNTIME_V2_EVALUATOR_CONFIG, :runtime_v2))
            audit_root = mktempdir()
            try
                isolated_authority_fixture(audit_root; config_rel=config)
                audited = E._authority_snapshots(audit_root, config)
                @test audited.config_kind == expected_kind
                @test length(audited.authority_sha256) == 64
            finally
                isdir(audit_root) && make_tree_writable(audit_root)
                rm(audit_root; recursive=true, force=true)
            end
        end

        authority_root = mktempdir()
        try
            fixture = isolated_authority_fixture(authority_root;
                                                 config_rel=RUNTIME_V3_EVALUATOR_CONFIG)
            initial = E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG)
            correction_paths = [
                joinpath(authority_root, contract.root_manifest_path),
                joinpath(authority_root, contract.claim_path),
                joinpath(authority_root, contract.review_path),
                joinpath(authority_root, contract.publication_path),
                joinpath(authority_root, E._RUNTIME_V3_T12_ROOT,
                         contract.edge_model_product_path),
                joinpath(authority_root, E._RUNTIME_V3_T12_ROOT,
                         contract.chain_product_path),
                joinpath(authority_root, E._RUNTIME_V3_T12_ROOT,
                         contract.test_product_path),
            ]
            live_paths = [
                joinpath(authority_root, "test/lib/structured_assignment/edge_model.jl"),
                joinpath(authority_root, "test/lib/structured_assignment/chain_inference.jl"),
                joinpath(authority_root, "test/test_structured_chain_inference.jl"),
            ]
            @test all(count(==(path), [snapshot.path for snapshot in initial.snapshots]) == 1
                      for path in vcat(correction_paths, live_paths))
            bundle_files = String[]
            bundle_directories = String[]
            for (directory, _, files) in walkdir(joinpath(authority_root, E._RUNTIME_V3_T12_ROOT))
                push!(bundle_directories, directory)
                append!(bundle_files, joinpath(directory, file) for file in files)
            end
            bundle_snapshot_paths = [
                snapshot.path for snapshot in initial.snapshots
                if startswith(snapshot.path,
                             joinpath(authority_root, E._RUNTIME_V3_T12_ROOT) * "/")
            ]
            @test length(bundle_files) == E._RUNTIME_V3_T12_ROOT_FILE_COUNT
            @test length(bundle_directories) == E._RUNTIME_V3_T12_ROOT_DIRECTORY_COUNT
            @test length(bundle_snapshot_paths) == E._RUNTIME_V3_T12_ROOT_FILE_COUNT
            @test Set(bundle_snapshot_paths) == Set(bundle_files)
            @test count(==(joinpath(authority_root, E._RUNTIME_V3_T12_PUBLICATION_PATH)),
                        [snapshot.path for snapshot in initial.snapshots]) == 1
            @test all(UInt(stat(path).mode) & UInt(0o777) == UInt(0o444)
                      for path in bundle_files)
            @test all(UInt(stat(path).mode) & UInt(0o777) == UInt(0o555)
                      for path in bundle_directories)
            for path in vcat(correction_paths, live_paths)
                original = Vector{UInt8}(read(path))
                rewrite_readonly(path, vcat(original, UInt8('\n')))
                @test authority_blocked(() ->
                    E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
                rewrite_readonly(path, original)
            end

            live_symlink_target = joinpath(authority_root, "live-symlink-target")
            write(live_symlink_target, "live substitution\n")
            victim = first(live_paths)
            original = Vector{UInt8}(read(victim))
            remove_readonly(victim)
            symlink(live_symlink_target, victim)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            remove_readonly(victim)
            write(victim, original)
            chmod(victim, 0o644)

            victim = live_paths[2]
            original = Vector{UInt8}(read(victim))
            remove_readonly(victim)
            hardlink(live_paths[1], victim)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            remove_readonly(victim)
            write(victim, original)
            chmod(victim, 0o644)

            nonselected = first(filter(path -> !(path in correction_paths), bundle_files))
            original = Vector{UInt8}(read(nonselected))
            rewrite_readonly(nonselected, vcat(original, UInt8('\n')))
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            rewrite_readonly(nonselected, original)

            chmod(nonselected, 0o644)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            chmod(nonselected, 0o444)

            nonselected_directory = dirname(nonselected)
            chmod(nonselected_directory, 0o755)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            chmod(nonselected_directory, 0o555)

            symlink_target = joinpath(authority_root, "bundle-symlink-target")
            write(symlink_target, "substitution\n")
            nonselected_parent = dirname(nonselected)
            remove_readonly(nonselected)
            chmod(nonselected_parent, 0o755)
            symlink(symlink_target, nonselected)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            remove_readonly(nonselected)
            write(nonselected, original)
            chmod(nonselected, 0o444)
            chmod(nonselected_parent, 0o555)

            hardlink_source = first(filter(path -> path != nonselected, bundle_files))
            remove_readonly(nonselected)
            chmod(nonselected_parent, 0o755)
            hardlink(hardlink_source, nonselected)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            remove_readonly(nonselected)
            write(nonselected, original)
            chmod(nonselected, 0o444)
            chmod(nonselected_parent, 0o555)

            extra = joinpath(authority_root, E._RUNTIME_V3_T12_ROOT, "extra-tamper.txt")
            bundle_root = dirname(joinpath(authority_root, E._RUNTIME_V3_T12_ROOT,
                                           "root-manifest.sha256"))
            chmod(bundle_root, 0o755)
            write(extra, "extra\n")
            chmod(extra, 0o444)
            chmod(bundle_root, 0o555)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            remove_readonly(extra)

            missing = nonselected
            missing_bytes = Vector{UInt8}(read(missing))
            remove_readonly(missing)
            @test authority_blocked(() ->
                E._authority_snapshots(authority_root, RUNTIME_V3_EVALUATOR_CONFIG))
            chmod(dirname(missing), 0o755)
            write(missing, missing_bytes)
            chmod(missing, 0o444)
            chmod(dirname(missing), 0o555)

            toctou_authority = E._authority_snapshots(
                authority_root, RUNTIME_V3_EVALUATOR_CONFIG)
            toctou = E._ProductionContext(
                authority_root, toctou_authority.snapshots,
                toctou_authority.inventories, "runtime-v3-toctou")
            rewrite_readonly(nonselected, vcat(original, UInt8('\n')))
            @test authority_blocked(() -> E._verify_context(toctou))
            rewrite_readonly(nonselected, original)
            chmod(nonselected, 0o644)
            @test authority_blocked(() -> E._verify_context(toctou))
            chmod(nonselected, 0o444)
            external = joinpath(authority_root, E._RUNTIME_V3_T12_PUBLICATION_PATH)
            external_bytes = Vector{UInt8}(read(external))
            rewrite_readonly(external, vcat(external_bytes, UInt8('\n')))
            @test authority_blocked(() -> E._verify_context(toctou))
            rewrite_readonly(external, external_bytes)
            chmod(external, 0o644)
            @test authority_blocked(() -> E._verify_context(toctou))
            chmod(external, 0o444)
            chmod(nonselected_directory, 0o755)
            @test authority_blocked(() -> E._verify_context(toctou))
            chmod(nonselected_directory, 0o555)

            toctou_link_target = joinpath(authority_root, "toctou-link-target")
            write(toctou_link_target, "toctou link\n")
            remove_readonly(nonselected)
            chmod(nonselected_parent, 0o755)
            symlink(toctou_link_target, nonselected)
            @test authority_blocked(() -> E._verify_context(toctou))
            remove_readonly(nonselected)
            write(nonselected, original)
            chmod(nonselected, 0o444)
            chmod(nonselected_parent, 0o555)

            remove_readonly(nonselected)
            chmod(nonselected_parent, 0o755)
            hardlink(hardlink_source, nonselected)
            @test authority_blocked(() -> E._verify_context(toctou))
            remove_readonly(nonselected)
            write(nonselected, original)
            chmod(nonselected, 0o444)
            chmod(nonselected_parent, 0o555)

            missing_toctou = Vector{UInt8}(read(nonselected))
            remove_readonly(nonselected)
            @test authority_blocked(() -> E._verify_context(toctou))
            chmod(nonselected_parent, 0o755)
            write(nonselected, missing_toctou)
            chmod(nonselected, 0o444)
            chmod(nonselected_parent, 0o555)

            external_link_target = joinpath(authority_root, "external-link-target")
            write(external_link_target, "external link\n")
            external_parent = dirname(external)
            remove_readonly(external)
            chmod(external_parent, 0o755)
            symlink(external_link_target, external)
            @test authority_blocked(() -> E._verify_context(toctou))
            remove_readonly(external)
            write(external, external_bytes)
            chmod(external, 0o444)

            chmod(bundle_root, 0o755)
            write(extra, "extra\n")
            chmod(extra, 0o444)
            chmod(bundle_root, 0o555)
            @test authority_blocked(() -> E._verify_context(toctou))
            remove_readonly(extra)
            extra_directory = joinpath(bundle_root, "toctou-extra-directory")
            chmod(bundle_root, 0o755)
            mkdir(extra_directory)
            chmod(extra_directory, 0o555)
            chmod(bundle_root, 0o555)
            @test authority_blocked(() -> E._verify_context(toctou))
            chmod(bundle_root, 0o755)
            rm(extra_directory; recursive=true, force=true)
            chmod(bundle_root, 0o555)
        finally
            make_tree_writable(authority_root)
            rm(authority_root; recursive=true, force=true)
        end

        production_options(config::String) = Dict{String,String}(
            "--root" => pwd(), "--evaluator-config" => config,
            "--features" => "missing/features.tsv",
            "--candidate-config" => "missing/candidate.toml",
            "--model-config" => "missing/model.toml",
            "--universe-dir" => "missing/universe", "--edge-dir" => "missing/edge",
            "--forward-receipt" => "missing/forward.toml",
            "--backward-receipt" => "missing/backward.toml",
            "--admission-dir" => "missing/admission")
        for config in (RUNTIME_EVALUATOR_CONFIG, RUNTIME_V2_EVALUATOR_CONFIG)
            old_error = try
                E._load_production_input(production_options(config))
                nothing
            catch error
                error
            end
            @test old_error isa E.EvaluatorError
            @test old_error.reason == :authority_path_mismatch
        end
        v3_error = try
            E._load_production_input(production_options(RUNTIME_V3_EVALUATOR_CONFIG))
            nothing
        catch error
            error
        end
        @test v3_error isa E.EvaluatorError
        @test v3_error.reason == :authority_path_invalid
        @test occursin("candidate config is absent", v3_error.message)
    end

    @testset "T13 runtime-v2 main config preflight and blocker publication" begin
        function main_args(root::String, config::String, output::String)
            return String[
                "--root", root, "--evaluator-config", config,
                "--features", "missing/features.tsv",
                "--candidate-config", "missing/candidate.toml",
                "--model-config", "missing/model.toml",
                "--universe-dir", "missing/universe",
                "--edge-dir", "missing/edge",
                "--forward-receipt", "missing/forward.toml",
                "--backward-receipt", "missing/backward.toml",
                "--admission-dir", "missing/admission",
                "--out-dir", relpath(output, root),
            ]
        end
        function run_main(arguments::Vector{String})
            stdout_path, stdout = mktemp()
            stderr_path, stderr = mktemp()
            status = try
                redirect_stdout(stdout) do
                    redirect_stderr(stderr) do
                        E.main(arguments)
                    end
                end
            finally
                close(stdout)
                close(stderr)
            end
            output = (status, read(stdout_path, String), read(stderr_path, String))
            rm(stdout_path; force=true)
            rm(stderr_path; force=true)
            return output
        end
        function blocker_config_sha(output::String)
            receipt = TOML.parse(read(joinpath(output, "receipt.toml"), String))
            return receipt["evaluator_config_sha256"]
        end

        sandbox = mktempdir(pwd(); prefix=".tmp-t13-main-preflight-", cleanup=false)
        try
            v2_output = joinpath(sandbox, "v2-blocker")
            v2_options = Dict{String,String}(
                "--root" => pwd(), "--evaluator-config" => RUNTIME_V2_EVALUATOR_CONFIG)
            config_only = E._preflight_config_context(v2_options)
            @test length(config_only.snapshots) == 1
            @test isempty(config_only.inventories)
            @test E._context_evaluator_config_sha256(config_only) ==
                  E._RUNTIME_V2_CONFIG_SHA256
            @test config_only.manifest_sha256 ==
                  E._preflight_config_context(v2_options).manifest_sha256

            v3_output = joinpath(sandbox, "v3-blocker")
            v3_options = Dict{String,String}(
                "--root" => pwd(), "--evaluator-config" => RUNTIME_V3_EVALUATOR_CONFIG)
            v3_config_only = E._preflight_config_context(v3_options)
            @test length(v3_config_only.snapshots) == 1
            @test E._context_evaluator_config_sha256(v3_config_only) ==
                  E._RUNTIME_V3_CONFIG_SHA256
            @test v3_config_only.manifest_sha256 ==
                  E._preflight_config_context(v3_options).manifest_sha256
            v3_status, _, v3_stderr = run_main(
                main_args(pwd(), RUNTIME_V3_EVALUATOR_CONFIG, v3_output))
            @test v3_status == 2
            @test occursin("candidate config is absent", v3_stderr)
            @test isfile(joinpath(v3_output, "receipt.toml"))
            @test readdir(v3_output) == ["receipt.toml"]
            @test blocker_config_sha(v3_output) == E._RUNTIME_V3_CONFIG_SHA256

            v2_status, _, v2_stderr = run_main(
                main_args(pwd(), RUNTIME_V2_EVALUATOR_CONFIG, v2_output))
            @test v2_status == 2
            @test occursin("evaluator config path differs", v2_stderr)
            @test isfile(joinpath(v2_output, "receipt.toml"))
            @test readdir(v2_output) == ["receipt.toml"]
            @test blocker_config_sha(v2_output) == E._RUNTIME_V2_CONFIG_SHA256

            v1_output = joinpath(sandbox, "v1-blocker")
            v1_status, _, v1_stderr = run_main(
                main_args(pwd(), RUNTIME_EVALUATOR_CONFIG, v1_output))
            @test v1_status == 2
            @test occursin("evaluator config path differs", v1_stderr)
            @test blocker_config_sha(v1_output) == E._RUNTIME_V1_CONFIG_SHA256

            historical_output = joinpath(sandbox, "historical-blocker")
            historical_status, _, historical_stderr = run_main(
                main_args(pwd(), HISTORICAL_EVALUATOR_CONFIG, historical_output))
            @test historical_status == 2
            @test occursin("evaluator config path differs", historical_stderr)
            @test blocker_config_sha(historical_output) == E._CONFIG_SHA256

            unknown_output = joinpath(sandbox, "unknown-no-blocker")
            unknown_status, _, unknown_stderr = run_main(main_args(
                pwd(), "config/unit_assignment_structured_evaluator_unknown.toml",
                unknown_output))
            @test unknown_status == 2
            @test occursin("accepted authority path", unknown_stderr)
            @test !ispath(unknown_output)

            mutable_root = mktempdir(; prefix="stmfit-t13-main-mutable-", cleanup=false)
            try
                mkpath(joinpath(mutable_root, "config"))
                mutable_config = joinpath(mutable_root, RUNTIME_V2_EVALUATOR_CONFIG)
                cp(joinpath(pwd(), RUNTIME_V2_EVALUATOR_CONFIG), mutable_config; force=true)
                mutable_output = joinpath(mutable_root, "blocker-after-mutation")
                mutable_options = Dict{String,String}(
                    "--root" => mutable_root,
                    "--evaluator-config" => RUNTIME_V2_EVALUATOR_CONFIG)
                @test E._context_evaluator_config_sha256(
                    E._preflight_config_context(mutable_options)) ==
                      E._RUNTIME_V2_CONFIG_SHA256
                original = Vector{UInt8}(read(mutable_config))
                mutated = replace(String(copy(original)),
                    "t12_edge_model_sha256" => "t12_edge_model_sha256_tampered";
                    count=1)
                write(mutable_config, mutated)
                mutation_status, _, mutation_stderr = run_main(
                    main_args(mutable_root, RUNTIME_V2_EVALUATOR_CONFIG, mutable_output))
                @test mutation_status == 2
                @test occursin("bytes differ", mutation_stderr)
                @test !ispath(mutable_output)
            finally
                isdir(mutable_root) && rm(mutable_root; recursive=true, force=true)
            end
        finally
            isdir(sandbox) && rm(sandbox; recursive=true, force=true)
        end
    end

end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    counts = Test.get_test_counts(RESULTS)
    println("assertions=", counts.cumulative_passes,
            " fails=", counts.cumulative_fails,
            " errors=", counts.cumulative_errors,
            " broken=", counts.cumulative_broken)
end

end
