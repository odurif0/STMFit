Base.@kwdef struct CliOptions2
    config::String
    seed::Int
    cases::Int
    out::String
    wrong_proxy_hash::Bool = false
    swap_type_mapping::Bool = false
end

_fail(msg) = error(String(msg))

function parse_cli(args::AbstractVector{<:AbstractString})::CliOptions2
    cfg = ""; seed = typemin(Int); cases = 0; out = ""; wrong = false; swap = false; i = 1
    while i <= length(args)
        a = String(args[i])
        if a == "--config"; cfg = String(args[i + 1]); i += 2
        elseif a == "--seed"; seed = parse(Int, String(args[i + 1])); i += 2
        elseif a == "--cases"; cases = parse(Int, String(args[i + 1])); i += 2
        elseif a == "--out"; out = String(args[i + 1]); i += 2
        elseif a == "--wrong-proxy-hash"; wrong = true; i += 1
        elseif a == "--swap-type-mapping"; swap = true; i += 1
        else; _fail("unknown argument: $a")
        end
    end
    isempty(cfg) && _fail("--config is required")
    seed != typemin(Int) || _fail("--seed is required")
    cases > 0 || _fail("--cases must be positive")
    isempty(out) && _fail("--out is required")
    return CliOptions2(config=cfg, seed=seed, cases=cases, out=out, wrong_proxy_hash=wrong, swap_type_mapping=swap)
end

const VALIDATION_SOURCE_FILES = [
    normpath(joinpath(@__DIR__, "..", "..", "synthetic_joint_proxy_validation.jl")),
    joinpath(@__DIR__, "synthetic_validation.jl"),
    joinpath(@__DIR__, "synthetic_validation_cli_support.jl"),
    joinpath(@__DIR__, "synthetic_validation_oracle.jl"),
    joinpath(@__DIR__, "synthetic_validation_count.jl"),
    joinpath(@__DIR__, "synthetic_validation_report.jl"),
]

function _validation_source_hash()
    entries = [basename(path) * ":" * bytes2hex(sha256(read(path))) for path in sort(VALIDATION_SOURCE_FILES)]
    return bytes2hex(sha256(codeunits(join(entries, '\n'))))
end

_bundle(config::String) = (; config_hash=bytes2hex(sha256(read(config))),
    source_hash=_validation_source_hash(),
    payload_hash=Main.JointProxyRegistry.load_registry(config).payload_sha256,
    registry=Main.JointProxyRegistry.load_registry(config))

function _validate_calibration_provenance(path::AbstractString, bundle)
    parsed = TOML.parsefile(path)
    provenance = get(parsed, "provenance", Dict{String,Any}())
    get(provenance, "config_sha256", "") == bundle.config_hash || _fail("calibration config hash mismatch")
    get(provenance, "source_sha256", "") == bundle.source_hash || _fail("calibration source hash mismatch")
    get(provenance, "payload_sha256", "") == bundle.payload_hash || _fail("calibration payload hash mismatch")
    return nothing
end

_balanced_sizes(total::Int, groups::Int) = begin
    base, extra = divrem(total, groups)
    [base + (i <= extra ? 1 : 0) for i in 1:groups]
end

_noise_cfg() = SimulatorConfig()
_quiet_cfg() = SimulatorConfig(noise_sigma=0.0, blur_sigma_px=0.0, affine_scale_jitter=0.0,
    drift_strength=0.0, row_offset_sigma=0.0, correlated_noise_frac=0.0)

_controls(kind::Symbol, n::Int) = kind == :normal ? fill(CONTROL_NORMAL, n) :
    kind == :null ? fill(CONTROL_NO_MOLECULE, n) : kind == :identical ? fill(CONTROL_IDENTICAL_MOLDS, n) :
    kind == :corrupted ? fill(CONTROL_CORRUPTED_VIEW, n) : kind == :missing_bwd ? fill(CONTROL_MISSING_BWD, n) :
    fill(CONTROL_SWAPPED_TYPES, n)

_control_flag(control::Symbol) = control in (CONTROL_NO_MOLECULE, CONTROL_IDENTICAL_MOLDS) ? :null :
    control in (CONTROL_CORRUPTED_VIEW, CONTROL_SWAPPED_TYPES) ? :family_conflict : :ok

function _count_cases(cases, split::Symbol, seed::Int, reg, config_hash::String)
    [CountCalibrationCase(case_id=c.case_id, split=split, true_n=c.truth.N,
        report=physical_candidate_report(c, reg), seed=UInt64(seed), config_hash=config_hash) for c in cases]
end

function _type_cases(cases, split::Symbol, seed::Int, reg, config_hash::String, source_hash::String, payload_hash::String; swap=false)
    use_reg = swap ? swap_registry(reg) : reg
    out = TypeCalibrationCase[]
    for c in cases
        raw = oracle_geometry_type_result(c, use_reg; swap_type_mapping=false)
        push!(out, TypeCalibrationCase(case_id=c.case_id, split=split, true_types=copy(c.truth.sequence),
            family_reports=["proxy" => raw], seed=UInt64(seed), config_hash=config_hash,
            manifest_hash=payload_hash, source_hash=source_hash, control=split == :calibration ? :ok : _control_flag(c.control)))
    end
    return out
end

function _status(case::SyntheticCase, count_post, type_post)
    n_pred = count_post.predicted_n
    count_abstained = count_post.abstained || n_pred === missing
    count_correct = !count_abstained && Int(n_pred) == case.truth.N
    type_emitted_correct = 0
    type_abstained = 0
    exact_nonabstained_chain = true
    for i in eachindex(case.truth.sequence)
        pred = type_post.predicted_types[i]
        if ismissing(pred)
            type_abstained += 1
            exact_nonabstained_chain = false
        else
            type_emitted_correct += Int(pred == case.truth.sequence[i])
            exact_nonabstained_chain &= pred == case.truth.sequence[i]
        end
    end
    type_conf = isempty(type_post.confidences) ? 0.0 : maximum(type_post.confidences)
    return (; count_conf = count_post.confidence, count_abstained, count_correct,
        type_emitted_correct, type_total = length(case.truth.sequence), type_abstained,
        exact_nonabstained_chain, type_conf)
end

function _group_metrics(records, label::Symbol)
    statuses = [_status(r.case, r.count, r.type) for r in records]
    raw_count_correct = mean(r.raw_count_correct for r in records)
    raw_type_correct = mean(r.raw_type_correct for r in records)
    abstain_rate = if label in (:no_molecule, :null)
        mean(s.count_abstained && (s.type_abstained == s.type_total) for s in statuses)
    elseif label in (:identical_molds, :identical)
        mean(s.type_abstained == s.type_total for s in statuses)
    else
        mean(s.count_abstained && (s.type_abstained == s.type_total) for s in statuses)
    end
    return (; count_recovery=raw_count_correct, type_recovery=raw_type_correct,
        abstain_rate=abstain_rate,
        exact_chain_rate=mean(r.raw_exact_chain for r in records),
        count_confidence=mean(s.count_conf for s in statuses), type_confidence=mean(s.type_conf for s in statuses))
end

function _write_tsv(path::AbstractString, rows)
    mkpath(dirname(path))
    fields = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for row in rows
            println(io, join((get(row, f, "") for f in fields), '\t'))
        end
    end
    return path
end

function _print_case_summary(io::IO, stage::AbstractString, case::SyntheticCase, report)
    s = candidate_report_summary(report)
    println(io, join((
        stage,
        case.case_id,
        String(case.control),
        string(s.total),
        string(s.finite),
        string(s.valid),
        string(s.finite_valid),
        s.best_n === missing ? "missing" : string(s.best_n),
        @sprintf("%.6f", s.best_gcv),
    ), '\t'))
end
