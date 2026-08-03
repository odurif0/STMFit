using SHA
using TOML

export ValidationError, ValidationOptions, ValidationReport,
       CANDIDATE_N_FIELDS, CANDIDATE_LOBE_FIELDS, PREDICTION_FIELDS, SUMMARY_FIELDS,
       parse_cli

const FORBIDDEN_KEYS = Set(["truth", "expected_N", "expected_n", "sequence", "target_N", "control_sequence"])
const FLOAT_TOL = 1e-9

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

struct ValidationError <: Exception
    msg::String
end

Base.showerror(io::IO, e::ValidationError) = print(io, "ValidationError: ", e.msg)

Base.@kwdef struct ValidationOptions
    artifacts::String
    calibration::String
end

Base.@kwdef struct ValidationReport
    files::Int
    candidate_rows::Int
    candidate_lobe_rows::Int
    prediction_rows::Int
    summary_rows::Int
end

_fail(msg::AbstractString) = throw(ValidationError(String(msg)))

function _arg_value(args, i::Int, flag::String)
    i < length(args) || _fail("$flag requires a value")
    return String(args[i + 1])
end

function _parse_cli(args::AbstractVector{<:AbstractString})::ValidationOptions
    artifacts = ""; calibration = ""; i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg == "--artifacts"; artifacts = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--artifacts="); artifacts = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--calibration"; calibration = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--calibration="); calibration = split(arg, "="; limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/validate_joint_proxy_predictions.jl --artifacts DIR --calibration PATH

            Validates candidate_n.tsv, candidate_lobes.tsv, predictions.tsv,
            chain_summary.tsv, and run_manifest.toml.
            """)
            exit(0)
        else
            _fail("Unknown argument: $arg")
        end
    end
    isempty(artifacts) && _fail("--artifacts is required")
    isempty(calibration) && _fail("--calibration is required")
    isdir(artifacts) || _fail("Artifacts directory not found: $artifacts")
    isfile(calibration) || _fail("Calibration file not found: $calibration")
    return ValidationOptions(artifacts=artifacts, calibration=calibration)
end
