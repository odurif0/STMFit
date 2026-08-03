#!/usr/bin/env julia

using Random
using Statistics
using SHA
using Printf
using TOML
using STMSXMIO

include(joinpath(@__DIR__, "lib", "joint_proxy", "synthetic_validation.jl"))
using .JointProxySyntheticValidation

using Main.JointProxySimulator: SyntheticCase, SimulatorConfig, default_proxy_ensemble, generate_batch, CONTROL_NORMAL,
    CONTROL_NO_MOLECULE, CONTROL_IDENTICAL_MOLDS, CONTROL_CORRUPTED_VIEW, CONTROL_MISSING_BWD, CONTROL_SWAPPED_TYPES
using Main.JointProxyCountCalibration: CountCalibrationCase, calibrate_count_temperature, predict_count
using Main.JointProxyTypeCalibration: TypeCalibrationCase, calibrate_type_temperature, predict_type
using Main.JointProxyTypePosterior: TypePosteriorResult

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_EVIDENCE = joinpath(ROOT, ".omo", "evidence", "joint-proxy-mold-inference-v1")

include(joinpath(@__DIR__, "lib", "joint_proxy", "synthetic_validation_cli_support.jl"))

function run_cli(args::AbstractVector{<:AbstractString}=String[])
    opt = parse_cli(args)
    bundle = _bundle(opt.config)
    reg = bundle.registry
    type_reg = opt.swap_type_mapping ? swap_registry(reg) : reg
    cal_seed = opt.seed + 101
    eval_seed = opt.seed + 202
    type_seed = opt.seed + 303

    cal_n = max(20, cld(opt.cases, 2))
    hold_n = cal_n
    sim_ensemble = default_proxy_ensemble()
    cal_cases = generate_batch(MersenneTwister(cal_seed), _quiet_cfg(), sim_ensemble; n_cases=cal_n, controls=_controls(:normal, cal_n))
    hold_cases = generate_batch(MersenneTwister(cal_seed + 1), _noise_cfg(), sim_ensemble; n_cases=hold_n, controls=_controls(:normal, hold_n))

    count_cases = vcat(_count_cases(cal_cases, :calibration, cal_seed, bundle.registry, bundle.config_hash),
                       _count_cases(hold_cases, :heldout, cal_seed, bundle.registry, bundle.config_hash))
    count_diag_lines = ["stage\tcase_id\tcontrol\ttotal\tfinite\tvalid\tfinite_valid\tbest_n\tbest_gcv"]
    for c in count_cases
        case = first(filter(x -> x.case_id == c.case_id, vcat(cal_cases, hold_cases)))
        push!(count_diag_lines, sprint(io -> _print_case_summary(io, String(c.split), case, c.report)))
    end
    count_model = calibrate_count_temperature(count_cases; seed=cal_seed, config_hash=bundle.config_hash)

    type_cases = vcat(_type_cases(cal_cases, :calibration, type_seed, bundle.registry, bundle.config_hash, bundle.source_hash, bundle.payload_hash),
                      _type_cases(hold_cases, :heldout, type_seed, bundle.registry, bundle.config_hash, bundle.source_hash, bundle.payload_hash))
    type_diag_lines = ["stage\tcase_id\tcontrol\tfamily_reports\ttrue_types"]
    for c in type_cases
        push!(type_diag_lines, join((String(c.split), c.case_id, String(c.control), string(length(c.family_reports)), string(length(c.true_types))), '\t'))
    end
    type_model = calibrate_type_temperature(type_cases; seed=type_seed, config_hash=bundle.config_hash,
        manifest_hash=bundle.payload_hash, source_hash=bundle.source_hash)

    evidence_dir = DEFAULT_EVIDENCE
    mkpath(evidence_dir)
    cal_path = joinpath(evidence_dir, "task-10-synthetic-e2e-calibration.toml")
    write_calibration_toml(cal_path, bundle, count_model, type_model; type_adapter=oracle_geometry_adapter_name())
    if opt.wrong_proxy_hash
        txt = read(cal_path, String)
        txt = replace(txt, bundle.source_hash => repeat("0", 64); count=1)
        write(cal_path, txt)
    end
    _validate_calibration_provenance(cal_path, bundle)

    groups = [(:noiseless, _quiet_cfg(), CONTROL_NORMAL), (:low_noise, _noise_cfg(), CONTROL_NORMAL),
        (:null, _noise_cfg(), CONTROL_NO_MOLECULE), (:identical, _noise_cfg(), CONTROL_IDENTICAL_MOLDS),
        (:corrupted, _noise_cfg(), CONTROL_CORRUPTED_VIEW)]
    sizes = _balanced_sizes(opt.cases, length(groups))
    eval_cases = SyntheticCase[]
    for (i, (label, cfg, _)) in enumerate(groups)
        append!(eval_cases, generate_batch(MersenneTwister(eval_seed + i), cfg, sim_ensemble; n_cases=sizes[i],
            controls=_controls(label == :noiseless || label == :low_noise ? :normal : label, sizes[i])))
    end

    artifact_dir = joinpath(evidence_dir, "task-10-synthetic-e2e-artifacts")
    mkpath(artifact_dir)
    count_rows = NamedTuple[]
    pred_rows = NamedTuple[]
    eval_records = NamedTuple[]
    for case in eval_cases
        report = physical_candidate_report(case, reg)
        count_post = predict_count(report, count_model)
        raw = oracle_geometry_type_result(case, type_reg; swap_type_mapping=false)
        type_post = predict_type(TypeCalibrationCase(case_id=case.case_id, split=:heldout, true_types=copy(case.truth.sequence),
            family_reports=["proxy" => raw], seed=UInt64(type_seed), config_hash=bundle.config_hash,
            manifest_hash=bundle.payload_hash, source_hash=bundle.source_hash, control=_control_flag(case.control)), type_model)
        raw_count_n = selected_candidate(report).n
        raw_type_map = copy(raw.map_sequence)
        push!(eval_records, (; case, count=count_post, type=type_post, raw_count_n, raw_type_map,
            raw_count_correct = !ismissing(raw_count_n) && Int(raw_count_n) == case.truth.N,
            raw_type_correct = mean(raw_type_map .== case.truth.sequence),
            raw_exact_chain = raw_type_map == case.truth.sequence))

        summ = candidate_report_summary(report)
        push!(count_rows, (; file=case.case_id, split="heldout", control=String(case.control),
            candidate_total=string(summ.total), candidate_finite=string(summ.finite), candidate_valid=string(summ.valid),
            candidate_finite_valid=string(summ.finite_valid), candidate_best_n=summ.best_n === missing ? "NA" : string(summ.best_n),
            candidate_best_gcv=@sprintf("%.6f", summ.best_gcv), predicted_n=count_post.predicted_n === missing ? "?" : string(count_post.predicted_n),
            confidence=@sprintf("%.6f", count_post.confidence), abstained=string(count_post.abstained), source_sha256=bundle.source_hash,
            config_sha256=bundle.config_hash, payload_sha256=bundle.payload_hash))
        for row in type_post.rows
            push!(pred_rows, (; file=case.case_id, split="heldout", control=String(case.control), lobe=string(row.lobe),
                predicted=row.predicted_type === missing ? "?" : string(row.predicted_type), confidence=@sprintf("%.6f", row.confidence),
                abstained=string(type_post.abstained[row.lobe]), flags=type_post.flags, source_sha256=bundle.source_hash,
                config_sha256=bundle.config_hash, payload_sha256=bundle.payload_hash))
        end
    end

    _write_tsv(joinpath(artifact_dir, "candidate_n.tsv"), count_rows)
    _write_tsv(joinpath(artifact_dir, "predictions.tsv"), pred_rows)

    replay_cases = SyntheticCase[]
    for (i, (label, cfg, _)) in enumerate(groups)
        append!(replay_cases, generate_batch(MersenneTwister(eval_seed + i), cfg, sim_ensemble; n_cases=sizes[i],
            controls=_controls(label == :noiseless || label == :low_noise ? :normal : label, sizes[i])))
    end
    replay_ok = replay_cases == eval_cases
    if replay_ok
        for (case, original) in zip(replay_cases, eval_records)
            report = physical_candidate_report(case, reg)
            count_post = predict_count(report, count_model)
            raw = oracle_geometry_type_result(case, type_reg; swap_type_mapping=false)
            type_post = predict_type(TypeCalibrationCase(case_id=case.case_id, split=:heldout, true_types=copy(case.truth.sequence),
                family_reports=["proxy" => raw], seed=UInt64(type_seed), config_hash=bundle.config_hash,
                manifest_hash=bundle.payload_hash, source_hash=bundle.source_hash, control=_control_flag(case.control)), type_model)
            replay_ok &= selected_candidate(report).n == original.raw_count_n
            replay_ok &= isequal(count_post.predicted_n, original.count.predicted_n) && count_post.confidence == original.count.confidence
            replay_ok &= raw.map_sequence == original.raw_type_map
            replay_ok &= isequal(type_post.predicted_types, original.type.predicted_types) && type_post.confidences == original.type.confidences
        end
    end

    bounds = cumsum(vcat(0, sizes))
    group_names = ["noiseless", "low_noise", "null", "identical", "corrupted"]
    group_rows = NamedTuple[]
    for (i, label) in enumerate(group_names)
        subset = eval_records[(bounds[i] + 1):bounds[i + 1]]
        metrics = _group_metrics(subset, Symbol(label))
        count_pass = label == "noiseless" ? metrics.count_recovery >= 1.0 - 1e-12 :
            label in ("low_noise", "corrupted") ? metrics.count_recovery >= 0.9 : metrics.abstain_rate >= 0.9
        type_pass = label == "noiseless" ? metrics.type_recovery >= 1.0 - 1e-12 :
            label in ("low_noise", "corrupted") ? metrics.type_recovery >= 0.9 : metrics.abstain_rate >= 0.9
        push!(group_rows, (; metric="$(label)_count_recovery", value=@sprintf("%.6f", metrics.count_recovery), threshold=label == "noiseless" ? "1.0" : "0.9", pass=string(count_pass), details="count"))
        push!(group_rows, (; metric="$(label)_type_recovery", value=@sprintf("%.6f", metrics.type_recovery), threshold=label == "noiseless" ? "1.0" : "0.9", pass=string(type_pass), details="type"))
        push!(group_rows, (; metric="$(label)_exact_chain_rate", value=@sprintf("%.6f", metrics.exact_chain_rate), threshold="diagnostic", pass="n/a", details="type"))
        if label in ("null", "identical")
            push!(group_rows, (; metric="$(label)_abstention", value=@sprintf("%.6f", metrics.abstain_rate), threshold="0.9", pass=string(metrics.abstain_rate >= 0.9), details="count+type"))
        end
    end

    duplicate = first(filter(c -> c.control == CONTROL_NORMAL, eval_cases))
    dup = SyntheticCase(duplicate.case_id * "_dup", STMSXMIO.SXMImage(duplicate.img.filepath, duplicate.img.header, duplicate.img.width, duplicate.img.height,
        duplicate.img.range_nm, duplicate.img.offset_nm, [duplicate.img.channels[1], STMSXMIO.SXMChannel("Z", "arb", "bwd", duplicate.img.channels[1].data)]), duplicate.truth, duplicate.control, duplicate.case_seed)
    dup_report = physical_candidate_report(dup, reg)
    dup_count = predict_count(dup_report, count_model)
    dup_raw = oracle_geometry_type_result(dup, type_reg; swap_type_mapping=false)
    dup_type = predict_type(TypeCalibrationCase(case_id=dup.case_id, split=:heldout, true_types=copy(dup.truth.sequence), family_reports=["proxy" => dup_raw], seed=UInt64(type_seed), config_hash=bundle.config_hash, manifest_hash=bundle.payload_hash, source_hash=bundle.source_hash, control=:ok), type_model)
    original = first(filter(r -> r.case.case_seed == duplicate.case_seed && r.case.control == duplicate.control, eval_records))
    original_count_confidence = original.count.confidence
    original_type_confidence = first(original.type.rows).confidence
    duplicate_rows = NamedTuple[
        (; metric="duplicate_count_confidence", value=@sprintf("%.6f", dup_count.confidence), threshold="<= original", pass=string(dup_count.confidence <= original_count_confidence + 1e-12), details=@sprintf("orig=%.6f", original_count_confidence)),
        (; metric="duplicate_type_confidence", value=@sprintf("%.6f", first(dup_type.rows).confidence), threshold="<= original", pass=string(first(dup_type.rows).confidence <= original_type_confidence + 1e-12), details=@sprintf("orig=%.6f", original_type_confidence)),
        (; metric="deterministic_replay", value=string(replay_ok), threshold="true", pass=string(replay_ok), details="regenerated cases and replayed posteriors")
    ]

    metrics_path = joinpath(evidence_dir, "task-10-synthetic-e2e.tsv")
    write_metrics_tsv(metrics_path, vcat(group_rows, duplicate_rows))
    cp(metrics_path, opt.out; force=true)

    log_path = joinpath(evidence_dir, "task-10-synthetic-e2e.log")
    open(log_path, "w") do io
        println(io, "synthetic joint proxy validation")
        println(io, "config_sha256\t", bundle.config_hash)
        println(io, "payload_sha256\t", bundle.payload_hash)
        println(io, "source_sha256\t", bundle.source_hash)
        println(io, "count_diagnostics")
        foreach(line -> println(io, line), count_diag_lines)
        println(io, "type_diagnostics")
        foreach(line -> println(io, line), type_diag_lines)
        println(io, "eval_cases\t", length(eval_cases))
        println(io, "artifact_dir\t", artifact_dir)
        println(io, "metrics\t", metrics_path)
        println(io, "count_temperature\t", count_model.temperature)
        println(io, "type_temperature\t", type_model.temperature)
        println(io, "type_adapter\t", oracle_geometry_adapter_name())
        println(io, "duplicate_count_confidence\t", dup_count.confidence)
        println(io, "duplicate_type_confidence\t", first(dup_type.rows).confidence)
        println(io, "noiseless_count_recovery\t", first(group_rows).value)
    end
    println(metrics_path)
    rows = vcat(group_rows, duplicate_rows)
    return all(row -> row.pass in ("true", "n/a"), rows) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(run_cli(ARGS))
end
