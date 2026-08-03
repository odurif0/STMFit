struct TypeCalibrationError <: Exception
    msg::String
end

Base.showerror(io::IO, e::TypeCalibrationError) = print(io, "TypeCalibrationError: ", e.msg)

Base.@kwdef struct TypeCalibrationCase
    case_id::String
    split::Symbol
    true_types::Vector{Int}
    family_reports::Vector{Pair{String,TypePosteriorResult}}
    seed::UInt64
    config_hash::String
    manifest_hash::String
    source_hash::String
    control::Symbol = :ok
end

struct TypePosteriorRow
    lobe::Int
    probability_0::Float64
    probability_1::Float64
    predicted_type::Union{Int,Missing}
    confidence::Float64
    rank::Int
end

struct TypeCasePosterior
    case_id::String
    true_types::Vector{Int}
    rows::Vector{TypePosteriorRow}
    predicted_types::Vector{Union{Int,Missing}}
    confidences::Vector{Float64}
    abstained::Vector{Bool}
    nll::Float64
    flags::String
end

struct TypeReliabilityBin
    lo::Float64
    hi::Float64
    count::Int
    correct::Int
    mean_confidence::Float64
    mean_nll::Float64
end

struct TypeCalibrationDiagnostics
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

struct TypeCalibrationModel
    temperature::Float64
    confidence_threshold::Float64
    calibration_nll::Float64
    heldout_nll::Float64
    baseline_temp1_nll::Float64
    diagnostics::TypeCalibrationDiagnostics
    reliability_bins::Vector{TypeReliabilityBin}
    config_hash::String
    manifest_hash::String
    source_hash::String
end

function _valid_type_case(c::TypeCalibrationCase)
    !isempty(c.case_id) && !isempty(c.true_types) && !isempty(c.family_reports) || return false
    all(t -> t in (0, 1), c.true_types) || return false
    for (_, r) in c.family_reports
        hasproperty(r, :lobe_marginals) || return false
        size(r.lobe_marginals, 2) == 2 || return false
        size(r.lobe_marginals, 1) == length(c.true_types) || return false
        all(isfinite, r.lobe_marginals) || return false
    end
    return true
end

function _ensure_same_seed_config_manifest_source(cases::AbstractVector{TypeCalibrationCase}, seed::Integer,
                                                  config_hash::AbstractString,
                                                  manifest_hash::AbstractString,
                                                  source_hash::AbstractString)
    seed_u64 = UInt64(seed)
    for c in cases
        c.seed == seed_u64 || throw(TypeCalibrationError(
            "config/seed mismatch: case $(c.case_id) seed=$(c.seed) != requested seed=$seed_u64"))
        c.config_hash == config_hash || throw(TypeCalibrationError(
            "config/seed mismatch: case $(c.case_id) config_hash=$(c.config_hash) != requested config_hash=$config_hash"))
        c.manifest_hash == manifest_hash || throw(TypeCalibrationError(
            "manifest mismatch: case $(c.case_id) manifest_hash=$(c.manifest_hash) != requested manifest_hash=$manifest_hash"))
        c.source_hash == source_hash || throw(TypeCalibrationError(
            "source mismatch: case $(c.case_id) source_hash=$(c.source_hash) != requested source_hash=$source_hash"))
    end
end

function _validate_calibration_cases(cases::AbstractVector{TypeCalibrationCase})
    isempty(cases) && throw(TypeCalibrationError("zero cases: synthetic calibration requires at least one case"))
    all(_valid_type_case, cases) || throw(TypeCalibrationError("malformed type calibration case"))
    cal = [c for c in cases if c.split == :calibration]
    held = [c for c in cases if c.split == :heldout]
    isempty(cal) && throw(TypeCalibrationError("zero cases: no calibration split provided"))
    isempty(held) && throw(TypeCalibrationError("zero cases: no heldout split provided"))
    any(c.control != :ok for c in cal) && throw(TypeCalibrationError("calibration split may not contain null/disagreement controls"))
    return cal, held
end
