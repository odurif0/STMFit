module PatchPreprocessing

using TOML

export PreprocessingSettings, load_patch_preprocessing

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

end # module
