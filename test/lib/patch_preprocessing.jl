module PatchPreprocessing

using TOML
using STMSXMIO: _box_smooth

export PreprocessingSettings, load_patch_preprocessing, load_patch_residual_filter,
       patch_residual

struct PreprocessingSettings
    stride::Int
    flatten::String
    smooth_radius_px::Int
end

"""
Read the three patch-preprocessing fields from a count config. A supplied config
must state all three fields. Omitting the option preserves the legacy extractor
settings exactly: stride=1, flatten="plane+rows", smooth_radius_px=1.
"""
function load_patch_preprocessing(path::Union{Nothing,AbstractString}=nothing)
    path === nothing && return PreprocessingSettings(1, "plane+rows", 1)
    isempty(path) && throw(ArgumentError("--config path must not be empty"))
    isfile(path) || throw(ArgumentError("Config not found: $path"))
    config = TOML.parsefile(path)
    preproc = get(config, "preprocessing", nothing)
    preproc isa AbstractDict || throw(ArgumentError("Config requires [preprocessing]"))
    function integer(key, minimum)
        value = get(preproc, key, nothing)
        value isa Integer && !(value isa Bool) && value >= minimum ||
            throw(ArgumentError("[preprocessing] $key must be an integer >= $minimum"))
        return Int(value)
    end
    stride = integer("stride", 1)
    radius = integer("smooth_radius_px", 0)
    flatten = get(preproc, "flatten", nothing)
    flatten isa AbstractString && lowercase(strip(flatten)) in ("none", "plane", "rows", "plane+rows") ||
        throw(ArgumentError("[preprocessing] flatten must be none, plane, rows, or plane+rows"))
    return PreprocessingSettings(stride, String(flatten), radius)
end

"Load the explicit assignment policy; omitted standalone option keeps legacy extraction."
function load_patch_residual_filter(path::Union{Nothing,AbstractString}=nothing)
    path === nothing && return "smooth_data_only"
    isempty(path) && throw(ArgumentError("Empty assignment config path"))
    cfg = TOML.parsefile(path)
    pre = get(cfg, "preprocessing", nothing)
    pre isa AbstractDict || throw(ArgumentError("Assignment config requires [preprocessing]"))
    mode = get(pre, "patch_residual_filter", nothing)
    mode in ("smooth_data_only", "smooth_residual") ||
        throw(ArgumentError("Explicit patch_residual_filter must be smooth_data_only or smooth_residual"))
    return String(mode)
end

"""
Residual in the preprocessed image's units, before patch interpolation/normalization.
Legacy: S(z) - M. Matched: S(z - M), using the SAME native finite-window smoother.
This does not change flattening, missing-value treatment, geometry or raw patches.
"""
function patch_residual(z::Matrix{Float64}, z_smooth::Matrix{Float64},
                        model::Matrix{Float64}, settings::PreprocessingSettings,
                        mode::AbstractString)
    size(z) == size(z_smooth) == size(model) || throw(DimensionMismatch("Residual image/model grids differ"))
    settings.smooth_radius_px >= 0 || throw(ArgumentError("Negative smoothing radius"))
    mode == "smooth_data_only" && return z_smooth .- model
    mode == "smooth_residual" && return _box_smooth(z .- model, settings.smooth_radius_px)
    throw(ArgumentError("Unknown patch residual filter: $mode"))
end

end # module
