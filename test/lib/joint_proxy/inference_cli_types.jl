struct InferenceCliError <: Exception
    msg::String
end

Base.showerror(io::IO, e::InferenceCliError) = print(io, "InferenceCliError: ", e.msg)

Base.@kwdef struct CliOptions
    config::String
    data_dir::String
    files::Vector{String}
    files_from::Union{Nothing,String}
    calibration::Union{Nothing,String}
    outdir::String
    chunk::Union{Nothing,Tuple{Int,Int}}
    uncalibrated::Bool
end

Base.@kwdef struct InferenceCalibration
    count_temperature::Float64
    count_threshold::Float64
    type_temperature::Float64
    type_threshold::Float64
    config_sha256::String
    source_sha256::String
    payload_sha256::String
    path::String
end

Base.@kwdef struct InferenceBundle
    config_hash::String
    source_hash::String
    payload_hash::String
    registry::Any
    calibration::Union{Nothing,InferenceCalibration}
end

Base.@kwdef struct InferenceArtifacts
    candidate_n_rows::Vector{Any}
    candidate_lobes_rows::Vector{Any}
    predictions_rows::Vector{Any}
    chain_summary_rows::Vector{Any}
    manifest::Any
end

const CANDIDATE_N_FIELDS = ["file", "n", "joint_gcv", "delta_rel", "probability", "rank",
    "is_map", "view_count", "bwd_missing", "residual_corr", "effective_view_factor",
    "effective_n_views", "p_eff", "nd_joint", "source_gcv", "joint_nrmse"]

const CANDIDATE_LOBE_FIELDS = ["file", "n", "lobe", "x_nm", "y_nm", "t_nm", "u_nm",
    "amplitude", "sigma_parallel_nm", "sigma_perp_nm", "skew_ratio", "type_probability_0",
    "type_probability_1", "predicted_type", "confidence"]

const PREDICTION_FIELDS = ["file", "lobe", "N_prediction", "N_probability", "N_confidence",
    "predicted", "confidence", "type_probability_0", "type_probability_1"]

const SUMMARY_FIELDS = ["file", "N_prediction", "N_probability", "N_confidence",
    "N_abstained", "candidate_count", "lobe_count", "view_count", "bwd_missing"]

_q(s::AbstractString) = '"' * replace(String(s), '"' => "\\\"") * '"'
_fmt(x::AbstractString) = _q(x)
_fmt(x::Bool) = x ? "true" : "false"
_fmt(x::Integer) = string(x)
_fmt(x::Real) = isfinite(Float64(x)) ? @sprintf("%.15g", Float64(x)) : "nan"
_fmt(xs::AbstractVector) = "[" * join((_fmt(x) for x in xs), ", ") * "]"
_fmt(xs::Tuple) = _fmt(collect(xs))
_fmt(x) = x === nothing ? "\"\"" : _q(string(x))

_tsv_fmt(x::AbstractString) = x
_tsv_fmt(x::Bool) = x ? "true" : "false"
_tsv_fmt(x::Integer) = string(x)
_tsv_fmt(x::Real) = isfinite(Float64(x)) ? @sprintf("%.15g", Float64(x)) : "NaN"
_tsv_fmt(x) = x === nothing ? "" : string(x)

function _field_value(row, field::String)
    sym = Symbol(field)
    if row isa NamedTuple
        hasproperty(row, sym) && return getproperty(row, sym)
        low = Symbol(lowercase(field))
        hasproperty(row, low) && return getproperty(row, low)
    elseif hasproperty(row, sym)
        return getproperty(row, sym)
    elseif hasproperty(row, Symbol(lowercase(field)))
        return getproperty(row, Symbol(lowercase(field)))
    end
    return row[field]
end
