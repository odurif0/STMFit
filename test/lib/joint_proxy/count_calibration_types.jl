struct CountCalibrationError <: Exception
    msg::String
end

Base.showerror(io::IO, e::CountCalibrationError) = print(io, "CountCalibrationError: ", e.msg)

Base.@kwdef struct CountCalibrationCase
    case_id::String
    split::Symbol
    true_n::Int
    report::Any
    seed::UInt64
    config_hash::String
end

struct CountPosteriorRow
    n::Int
    joint_gcv::Float64
    delta_rel::Float64
    probability::Float64
    rank::Int
end

struct CountCasePosterior
    case_id::String
    true_n::Int
    rows::Vector{CountPosteriorRow}
    predicted_n::Union{Int,Missing}
    confidence::Float64
    abstained::Bool
    true_rank::Int
    nll::Float64
end

struct CountReliabilityBin
    lo::Float64
    hi::Float64
    count::Int
    correct::Int
    mean_confidence::Float64
    mean_nll::Float64
end

struct CountCalibrationDiagnostics
    temperature_grid::Vector{Float64}
    calibration_nll_grid::Vector{Float64}
    calibration_nll::Float64
    heldout_nll::Float64
    baseline_temp1_nll::Float64
    selected_grid_index::Int
    threshold_source::String
    heldout_cases::Int
    calibration_cases::Int
end

struct CountCalibrationModel
    temperature::Float64
    confidence_threshold::Float64
    calibration_nll::Float64
    heldout_nll::Float64
    baseline_temp1_nll::Float64
    diagnostics::CountCalibrationDiagnostics
    reliability_bins::Vector{CountReliabilityBin}
end

_valid_count_candidate(c) = hasproperty(c, :valid) && hasproperty(c, :joint_gcv) && hasproperty(c, :n) &&
    getproperty(c, :valid) && isfinite(getproperty(c, :joint_gcv)) && isfinite(getproperty(c, :n))

function _ensure_same_seed_and_config(cases::AbstractVector{CountCalibrationCase}, seed::Integer,
                                      config_hash::AbstractString)
    seed_u64 = UInt64(seed)
    for c in cases
        c.seed == seed_u64 || throw(CountCalibrationError(
            "config/seed mismatch: case $(c.case_id) seed=$(c.seed) != requested seed=$seed_u64"))
        c.config_hash == config_hash || throw(CountCalibrationError(
            "config/seed mismatch: case $(c.case_id) config_hash=$(c.config_hash) != requested config_hash=$config_hash"))
    end
end

function _finite_candidates(report)
    hasproperty(report, :candidates) || throw(CountCalibrationError("report missing candidates field"))
    cands = [c for c in getproperty(report, :candidates) if _valid_count_candidate(c)]
    isempty(cands) && throw(CountCalibrationError(
        "all infinite scores: report has no finite valid count candidates"))
    return cands
end
