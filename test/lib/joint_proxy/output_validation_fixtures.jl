module JointProxyOutputValidationFixtures

include("output_validation.jl")
include(joinpath(@__DIR__, "..", "script_utils.jl"))
using .ScriptUtils: _ensure_parent

using SHA
using TOML
using Printf

export FIXTURE_FILES, SOURCE_HASH, PAYLOAD_HASH, fixture_rows, write_fixture,
       read_bytes, copy_fixture, write_tsv, config_hash

const FIXTURE_FILES = ["a.sxm", "b.sxm", "c.sxm"]
const SOURCE_HASH = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
const PAYLOAD_HASH = "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210"

config_hash(path::AbstractString) = bytes2hex(sha256(read(path)))
read_bytes(path::AbstractString) = read(path)

function write_tsv(path::AbstractString, fields::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for row in rows
            println(io, join((string(getproperty(row, Symbol(f))) for f in fields), '\t'))
        end
    end
end

_toml_q(s::AbstractString) = '"' * replace(String(s), '"' => "\\\"") * '"'
_toml_fmt(x::AbstractString) = _toml_q(x)
_toml_fmt(x::Bool) = x ? "true" : "false"
_toml_fmt(x::Integer) = string(x)
_toml_fmt(x::Real) = isfinite(Float64(x)) ? @sprintf("%.15g", Float64(x)) : "nan"
_toml_fmt(x::AbstractVector) = "[" * join((_toml_fmt(v) for v in x), ", ") * "]"
_toml_fmt(x::Tuple) = _toml_fmt(collect(x))
_toml_fmt(x) = x === nothing ? "\"\"" : _toml_q(string(x))

function _write_toml(io::IO, name::String, value)
    println(io, "[$name]")
    if value isa AbstractDict
        for (k, v) in pairs(value)
            if v isa AbstractDict || v isa NamedTuple
                _write_toml(io, "$name.$(String(k))", v)
            elseif v isa AbstractVector && !isempty(v) && (first(v) isa AbstractDict || first(v) isa NamedTuple)
                for item in v
                    println(io, "[[$name.$(String(k))]]")
                    for (ik, iv) in pairs(item)
                        println(io, "$(String(ik)) = $(_toml_fmt(iv))")
                    end
                end
            else
                println(io, "$(String(k)) = $(_toml_fmt(v))")
            end
        end
    else
        for (k, v) in pairs(value)
            println(io, "$(String(k)) = $(_toml_fmt(v))")
        end
    end
    println(io)
end

function fixture_rows()
    candidate_n = NamedTuple[
        (; file="a.sxm", n=2, joint_gcv=1.0, delta_rel=0.10, probability=0.30, rank=2, is_map=false,
           view_count=2, bwd_missing=false, residual_corr=0.20, effective_view_factor=0.83, effective_n_views=1.66,
           p_eff=8, nd_joint=12, source_gcv=1.4, joint_nrmse=0.12),
        (; file="a.sxm", n=3, joint_gcv=0.8, delta_rel=0.00, probability=0.70, rank=1, is_map=true,
           view_count=2, bwd_missing=false, residual_corr=0.20, effective_view_factor=0.83, effective_n_views=1.66,
           p_eff=8, nd_joint=12, source_gcv=1.1, joint_nrmse=0.08),
        (; file="b.sxm", n=2, joint_gcv=1.1, delta_rel=0.20, probability=0.60, rank=1, is_map=true,
           view_count=1, bwd_missing=true, residual_corr=0.00, effective_view_factor=1.0, effective_n_views=1.0,
           p_eff=4, nd_joint=8, source_gcv=1.5, joint_nrmse=0.10),
        (; file="b.sxm", n=4, joint_gcv=1.5, delta_rel=0.30, probability=0.40, rank=2, is_map=false,
           view_count=1, bwd_missing=true, residual_corr=0.00, effective_view_factor=1.0, effective_n_views=1.0,
           p_eff=4, nd_joint=8, source_gcv=1.8, joint_nrmse=0.15),
        (; file="c.sxm", n=2, joint_gcv=0.9, delta_rel=0.00, probability=1.0, rank=1, is_map=true,
           view_count=2, bwd_missing=false, residual_corr=0.10, effective_view_factor=0.90, effective_n_views=1.80,
           p_eff=6, nd_joint=10, source_gcv=1.0, joint_nrmse=0.05),
    ]
    candidate_lobes = NamedTuple[
        (; file="a.sxm", n=2, lobe=1, x_nm=0.1, y_nm=0.2, t_nm=0.0, u_nm=0.1, amplitude=0.7, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.2, type_probability_1=0.8, predicted_type=1, confidence=0.8),
        (; file="a.sxm", n=2, lobe=2, x_nm=0.3, y_nm=0.4, t_nm=0.1, u_nm=0.2, amplitude=0.6, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.7, type_probability_1=0.3, predicted_type=0, confidence=0.7),
        (; file="a.sxm", n=3, lobe=1, x_nm=0.1, y_nm=0.2, t_nm=0.0, u_nm=0.1, amplitude=0.7, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.4, type_probability_1=0.6, predicted_type=1, confidence=0.6),
        (; file="a.sxm", n=3, lobe=2, x_nm=0.2, y_nm=0.3, t_nm=0.1, u_nm=0.2, amplitude=0.6, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.5, type_probability_1=0.5, predicted_type="?", confidence=0.5),
        (; file="a.sxm", n=3, lobe=3, x_nm=0.3, y_nm=0.4, t_nm=0.2, u_nm=0.3, amplitude=0.5, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.6, type_probability_1=0.4, predicted_type=0, confidence=0.6),
        (; file="b.sxm", n=2, lobe=1, x_nm=0.1, y_nm=0.2, t_nm=0.0, u_nm=0.1, amplitude=0.5, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.5, type_probability_1=0.5, predicted_type="?", confidence=0.5),
        (; file="b.sxm", n=2, lobe=2, x_nm=0.2, y_nm=0.3, t_nm=0.1, u_nm=0.2, amplitude=0.4, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.3, type_probability_1=0.7, predicted_type="?", confidence=0.7),
        (; file="b.sxm", n=4, lobe=1, x_nm=0.1, y_nm=0.2, t_nm=0.0, u_nm=0.1, amplitude=0.5, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.2, type_probability_1=0.8, predicted_type=1, confidence=0.8),
        (; file="b.sxm", n=4, lobe=2, x_nm=0.2, y_nm=0.3, t_nm=0.1, u_nm=0.2, amplitude=0.4, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.8, type_probability_1=0.2, predicted_type=0, confidence=0.8),
        (; file="b.sxm", n=4, lobe=3, x_nm=0.3, y_nm=0.4, t_nm=0.2, u_nm=0.3, amplitude=0.3, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.6, type_probability_1=0.4, predicted_type=0, confidence=0.6),
        (; file="b.sxm", n=4, lobe=4, x_nm=0.4, y_nm=0.5, t_nm=0.3, u_nm=0.4, amplitude=0.2, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.4, type_probability_1=0.6, predicted_type=1, confidence=0.6),
        (; file="c.sxm", n=2, lobe=1, x_nm=0.1, y_nm=0.2, t_nm=0.0, u_nm=0.1, amplitude=0.7, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.1, type_probability_1=0.9, predicted_type=1, confidence=0.9),
        (; file="c.sxm", n=2, lobe=2, x_nm=0.2, y_nm=0.3, t_nm=0.1, u_nm=0.2, amplitude=0.6, sigma_parallel_nm=0.2,
           sigma_perp_nm=0.1, skew_ratio=1.0, type_probability_0=0.9, type_probability_1=0.1, predicted_type=0, confidence=0.9),
    ]
    predictions = NamedTuple[
        (; file="a.sxm", lobe=1, N_prediction=3, N_probability=0.70, N_confidence=0.70, predicted=1, confidence=0.80, type_probability_0=0.20, type_probability_1=0.80),
        (; file="a.sxm", lobe=2, N_prediction=3, N_probability=0.70, N_confidence=0.70, predicted=0, confidence=0.70, type_probability_0=0.70, type_probability_1=0.30),
        (; file="a.sxm", lobe=3, N_prediction=3, N_probability=0.70, N_confidence=0.70, predicted=1, confidence=0.60, type_probability_0=0.40, type_probability_1=0.60),
        (; file="b.sxm", lobe=1, N_prediction="?", N_probability=0.60, N_confidence=0.40, predicted="?", confidence=0.50, type_probability_0=0.50, type_probability_1=0.50),
        (; file="b.sxm", lobe=2, N_prediction="?", N_probability=0.60, N_confidence=0.40, predicted="?", confidence=0.50, type_probability_0=0.50, type_probability_1=0.50),
        (; file="c.sxm", lobe=1, N_prediction=2, N_probability=1.00, N_confidence=1.00, predicted=1, confidence=0.90, type_probability_0=0.10, type_probability_1=0.90),
        (; file="c.sxm", lobe=2, N_prediction=2, N_probability=1.00, N_confidence=1.00, predicted=0, confidence=0.90, type_probability_0=0.90, type_probability_1=0.10),
    ]
    summaries = NamedTuple[
        (; file="a.sxm", N_prediction=3, N_probability=0.70, N_confidence=0.70, N_abstained=false, candidate_count=2, lobe_count=5, view_count=2, bwd_missing=false),
        (; file="b.sxm", N_prediction="?", N_probability=0.60, N_confidence=0.40, N_abstained=true, candidate_count=2, lobe_count=6, view_count=1, bwd_missing=true),
        (; file="c.sxm", N_prediction=2, N_probability=1.00, N_confidence=1.00, N_abstained=false, candidate_count=1, lobe_count=2, view_count=2, bwd_missing=false),
    ]
    return candidate_n, candidate_lobes, predictions, summaries
end

function write_fixture(dir::AbstractString; config_hash::AbstractString, source_hash::AbstractString, payload_hash::AbstractString,
    config_path::AbstractString, files::Vector{String}=FIXTURE_FILES, chunk::AbstractString="none", data_dir::AbstractString=joinpath(dir, "data"),
    forbidden_column::Union{Nothing,String}=nothing, bad_source_hash::Bool=false)
    mkpath(dir)
    candidate_n, candidate_lobes, predictions, summaries = fixture_rows()
    file_set = Set(files)
    candidate_n = [row for row in candidate_n if row.file in file_set]
    candidate_lobes = [row for row in candidate_lobes if row.file in file_set]
    predictions = [row for row in predictions if row.file in file_set]
    summaries = [row for row in summaries if row.file in file_set]
    write_tsv(joinpath(dir, "candidate_n.tsv"), JointProxyOutputValidation.CANDIDATE_N_FIELDS, candidate_n)
    write_tsv(joinpath(dir, "candidate_lobes.tsv"), JointProxyOutputValidation.CANDIDATE_LOBE_FIELDS, candidate_lobes)
    write_tsv(joinpath(dir, "predictions.tsv"), JointProxyOutputValidation.PREDICTION_FIELDS, predictions)
    write_tsv(joinpath(dir, "chain_summary.tsv"), JointProxyOutputValidation.SUMMARY_FIELDS, summaries)
    prov_source = bad_source_hash ? "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff" : source_hash
    mani = Dict(
        "provenance" => Dict("config_sha256" => config_hash, "source_sha256" => prov_source, "payload_sha256" => payload_hash,
            "calibration_path" => "joint_proxy_calibration.toml", "uncalibrated" => false),
        "inputs" => Dict("config" => config_path, "data_dir" => data_dir, "files" => files, "chunk" => chunk),
        "outputs" => Dict("candidate_n_rows" => length(candidate_n), "candidate_lobes_rows" => length(candidate_lobes),
            "predictions_rows" => length(predictions), "chain_summary_rows" => length(summaries)),
    )
    forbidden_column !== nothing && (mani[forbidden_column] = Dict("expected_N" => 6))
    open(joinpath(dir, "run_manifest.toml"), "w") do io
        for (name, value) in pairs(mani)
            _write_toml(io, String(name), value)
        end
    end
    cal = Dict("provenance" => Dict("config_sha256" => config_hash, "source_sha256" => source_hash, "payload_sha256" => payload_hash),
        "count" => Dict("temperature" => 1.0, "confidence_threshold" => 0.6),
        "type" => Dict("temperature" => 1.0, "confidence_threshold" => 0.5))
    open(joinpath(dir, "joint_proxy_calibration.toml"), "w") do io
        TOML.print(io, cal)
    end
    return dir
end

function copy_fixture(src::AbstractString, dst::AbstractString)
    mkpath(dst)
    for name in ("candidate_n.tsv", "candidate_lobes.tsv", "predictions.tsv", "chain_summary.tsv", "run_manifest.toml", "joint_proxy_calibration.toml")
        cp(joinpath(src, name), joinpath(dst, name); force=true)
    end
    return dst
end

end # module
