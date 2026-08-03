using Random
using Printf
using SHA
using Statistics
using TOML

using Main.JointProxyRegistry: load_registry, sha256_file, RegistryConfigError
using Main.JointProxySimulator: SimulatorConfig, generate_batch, default_proxy_ensemble,
    CONTROL_NORMAL, CONTROL_NO_MOLECULE, CONTROL_IDENTICAL_MOLDS,
    CONTROL_SWAPPED_TYPES, CONTROL_MISSING_BWD, CONTROL_CORRUPTED_VIEW,
    ALL_CONTROLS, SyntheticCase, ControlType, generate_case
using Main.JointProxyCountCalibration: CountCalibrationCase, calibrate_count_temperature
using Main.JointProxyTypePosterior: TypePosteriorResult, TypePosteriorGlobalState
using Main.JointProxyTypeCalibration: TypeCalibrationCase, calibrate_type_temperature

struct CalibrationCliError <: Exception
    msg::String
end

Base.showerror(io::IO, e::CalibrationCliError) = print(io, "CalibrationCliError: ", e.msg)

Base.@kwdef struct CliOptions
    config::String
    seed::Int
    cases::Int
    n_min::Int
    n_max::Int
    fast::Bool
    out::String
end

Base.@kwdef struct CalibrationBundle
    config_hash::String
    source_hash::String
    payload_hash::String
    code_revision::Union{Nothing,String}
    registry::Any
    sim_cfg::SimulatorConfig
end

const FORBIDDEN_KEYS = Set(["truth", "benchmark", "manifest", "files", "files_from",
    "files-from", "expected_n", "expected-n", "target_n", "target-n",
    "data_dir", "data-dir", "grading"])

const CONTROL_CYCLE = (CONTROL_NORMAL, CONTROL_NO_MOLECULE, CONTROL_IDENTICAL_MOLDS,
    CONTROL_SWAPPED_TYPES, CONTROL_MISSING_BWD, CONTROL_CORRUPTED_VIEW)

struct SyntheticCountCandidate
    n::Int
    valid::Bool
    joint_gcv::Float64
end

struct SyntheticCountReport
    candidates::Vector{SyntheticCountCandidate}
    view_count::Int
    bwd_missing::Bool
    residual_corr::Float64
end

_q(s::AbstractString) = '"' * replace(String(s), '"' => "\\\"") * '"'
_fmt(x::AbstractString) = _q(x)
_fmt(x::Bool) = x ? "true" : "false"
_fmt(x::Integer) = string(x)
_fmt(x::Real) = isfinite(Float64(x)) ? @sprintf("%.15g", Float64(x)) : "nan"
_fmt(xs::AbstractVector{<:Real}) = "[" * join((_fmt(x) for x in xs), ", ") * "]"
_fmt(xs::AbstractVector{String}) = "[" * join((_q(x) for x in xs), ", ") * "]"

function _next_arg(args::Vector{String}, i::Int, flag::String)
    i < length(args) || throw(CalibrationCliError("$flag requires a value"))
    return args[i + 1]
end

function _scan_forbidden(node, path::String="root")
    node isa AbstractDict || return
    for (k, v) in node
        key = lowercase(String(k))
        key in FORBIDDEN_KEYS && throw(CalibrationCliError("forbidden truth/benchmark field in config: $path.$k"))
        _scan_forbidden(v, "$path.$k")
    end
end

function parse_cli(args::AbstractVector{<:AbstractString})::CliOptions
    config = ""; seed = typemin(Int); cases = 0; n_min = 2; n_max = 12; out = ""; fast = false
    i = 1
    while i <= length(args)
        a = String(args[i])
        if a == "--fast"
            fast = true; i += 1
        elseif a == "--config"
            config = _next_arg(collect(args), i, a); i += 2
        elseif a == "--seed"
            seed = parse(Int, _next_arg(collect(args), i, a)); i += 2
        elseif a == "--cases"
            cases = parse(Int, _next_arg(collect(args), i, a)); i += 2
        elseif a == "--n-min"
            n_min = parse(Int, _next_arg(collect(args), i, a)); i += 2
        elseif a == "--n-max"
            n_max = parse(Int, _next_arg(collect(args), i, a)); i += 2
        elseif a == "--out"
            out = _next_arg(collect(args), i, a); i += 2
        else
            throw(CalibrationCliError("Unknown argument: $a"))
        end
    end
    isempty(config) && throw(CalibrationCliError("--config is required"))
    seed != typemin(Int) || throw(CalibrationCliError("--seed is required"))
    cases > 0 || throw(CalibrationCliError("--cases must be positive"))
    n_min > 0 || throw(CalibrationCliError("--n-min must be positive"))
    n_max >= n_min || throw(CalibrationCliError("--n-max must be >= --n-min"))
    cases >= 4 || throw(CalibrationCliError("--cases must be at least 4 to split calibration and holdout"))
    isempty(out) && throw(CalibrationCliError("--out is required"))
    return CliOptions(config=config, seed=seed, cases=cases, n_min=n_min, n_max=n_max, fast=fast, out=out)
end

function _code_revision()
    for key in ("SOURCE_REVISION", "GIT_COMMIT", "CI_COMMIT_SHA", "GITHUB_SHA")
        haskey(ENV, key) || continue
        val = strip(ENV[key])
        isempty(val) || return val
    end
    return nothing
end

function _source_hash(paths::Vector{String})
    parts = String[]
    for p in paths
        isfile(p) || continue
        push!(parts, basename(p) * ":" * sha256_file(p))
    end
    return bytes2hex(sha256(codeunits(join(parts, "\n"))))
end

function _load_bundle(opt::CliOptions)
    isfile(opt.config) || throw(CalibrationCliError("config file not found: $(opt.config)"))
    config_hash = sha256_file(opt.config)
    cfg = TOML.parsefile(opt.config)
    _scan_forbidden(cfg)
    registry = load_registry(opt.config)
    source_hash = _source_hash(_calibration_source_files())
    payload_hash = registry.payload_sha256
    prov = get(cfg, "provenance", Dict{String,Any}())
    if haskey(prov, "expected_source_sha256")
        expected = String(prov["expected_source_sha256"])
        expected == source_hash || throw(CalibrationCliError("proxy source hash mismatch: expected $expected, actual $source_hash"))
    end
    if haskey(prov, "expected_payload_sha256")
        expected = String(prov["expected_payload_sha256"])
        expected == payload_hash || throw(CalibrationCliError("proxy payload hash mismatch: expected $expected, actual $payload_hash"))
    end
    if haskey(prov, "expected_config_sha256")
        expected = String(prov["expected_config_sha256"])
        expected == config_hash || throw(CalibrationCliError("config hash mismatch: expected $expected, actual $config_hash"))
    end
    return CalibrationBundle(config_hash=config_hash, source_hash=source_hash, payload_hash=payload_hash,
        code_revision=_code_revision(), registry=registry,
        sim_cfg=SimulatorConfig(n_min=opt.n_min, n_max=opt.n_max, seed=opt.seed))
end
