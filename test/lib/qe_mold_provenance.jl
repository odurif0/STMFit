module QEMoldProvenance

using SHA
using TOML

export write_qe_mold_provenance

# Fixed crossing-policy label recorded for the constant-current observable. It
# names the exact convention implemented by ConstantCurrentCube
# (vacuum -> slab search, first bracketed crossing, linear interpolation,
# ambiguity/absent/nonfinite rejected). Provenance binds the policy string so a
# consumer can detect a convention change without re-reading the engine.
const CC_CROSSING_POLICY = "vacuum_first_bracketed_linear"
const DEFAULT_ISOVALUE_SCAN_INTERVALS = 1024

_sha256_file(path::AbstractString) = open(path, "r") do io
    bytes2hex(SHA.sha256(io))
end

function _artifact_source(path::AbstractString, artifact_sources)
    logical_path = String(path)
    artifact_sources === nothing && return logical_path
    source = haskey(artifact_sources, logical_path) ? artifact_sources[logical_path] :
        haskey(artifact_sources, path) ? artifact_sources[path] : logical_path
    source_path = String(source)
    if source_path != logical_path
        islink(source_path) &&
            error("refusing symlinked staged QE mold artifact: $source_path")
    end
    return source_path
end

function _require_artifact(path::AbstractString, artifact_sources,
        message::AbstractString)
    source = _artifact_source(path, artifact_sources)
    isfile(source) || error("$message: $path")
    return source
end

function _atomic_write(writer::Function, path::AbstractString)
    parent = dirname(abspath(path))
    isdir(parent) || error("provenance output directory does not exist: $parent")
    tmp, io = mktemp(parent)
    committed = false
    try
        writer(io)
        flush(io)
        close(io)
        Base.Filesystem.rename(tmp, path)
        committed = true
    finally
        isopen(io) && close(io)
        !committed && ispath(tmp) && rm(tmp; force=true)
    end
    return path
end

# Validate that every constant-current observable field is present and well
# formed, returning the normalized values so the payload records Float64 /
# Vector{Float64} determinism. Uses error() (ErrorException) for the same
# "named diagnostic, nonzero exit" convention as the rest of this module.
function _entry_value(entry, key::Symbol; label::AbstractString="type frame")
    if entry isa NamedTuple
        haskey(entry, key) || error("$label missing $(String(key))")
        return entry[key]
    elseif entry isa AbstractDict
        if haskey(entry, key)
            return entry[key]
        elseif haskey(entry, String(key))
            return entry[String(key)]
        else
            error("$label missing $(String(key))")
        end
    else
        error("$label entries must be named tuples or dictionaries")
    end
end

function _finite_real(value, name::AbstractString; positive::Bool=false)
    (value isa Real && !(value isa Bool)) ||
        error("$name must be a non-Bool real number")
    out = Float64(value)
    isfinite(out) || error("$name must be finite")
    positive && out <= 0 && error("$name must be positive")
    return out
end

function _positive_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) || error("$name must be a positive integer")
    value > 0 || error("$name must be a positive integer")
    return Int(value)
end

function _type01(value, name::AbstractString)
    (value isa Real && !(value isa Bool)) ||
        error("$name must be a non-Bool numeric type")
    typ = Int(value)
    Float64(value) == Float64(typ) || error("$name must be an integer")
    typ in (0, 1) || error("$name must be 0 or 1")
    return typ
end

function _vec3(value, name::AbstractString)
    value isa AbstractVector || error("$name must be a vector")
    length(value) == 3 || error("$name must have length 3")
    return [_finite_real(v, "$name[$i]") for (i, v) in enumerate(value)]
end

function _normalize_type_frames(type_frames)
    type_frames === nothing &&
        error("constant-current provenance requires type_frames entries for type 0 and type 1")
    frames = Dict{Int,Dict{String,Any}}()
    for entry in type_frames
        typ = _type01(_entry_value(entry, :type), "type frame type")
        haskey(frames, typ) && error("duplicate type frame for type $typ")
        frames[typ] = Dict{String,Any}(
            "type" => typ,
            "origin_nm" => _vec3(_entry_value(entry, :origin_nm),
                                 "type frame $typ origin_nm"),
            "t_axis" => _vec3(_entry_value(entry, :t_axis),
                              "type frame $typ t_axis"),
            "u_axis" => _vec3(_entry_value(entry, :u_axis),
                              "type frame $typ u_axis"),
        )
    end
    Set(keys(frames)) == Set([0, 1]) ||
        error("constant-current provenance requires exactly one frame for type 0 and type 1")
    return Dict{String,Any}[frames[0], frames[1]]
end

function _normalize_type_isovalues(type_isovalues, legacy_isovalue)
    if type_isovalues === nothing
        legacy_isovalue === nothing &&
            error("constant-current provenance requires type_isovalues entries for type 0 and type 1")
        iso = _finite_real(legacy_isovalue, "isovalue")
        return Dict{String,Any}[
            Dict("type" => 0, "isovalue" => iso, "source" => "legacy_common"),
            Dict("type" => 1, "isovalue" => iso, "source" => "legacy_common"),
        ]
    end
    legacy_isovalue === nothing ||
        error("constant-current provenance accepts either type_isovalues or legacy isovalue, not both")
    entries = Dict{Int,Dict{String,Any}}()
    for entry in type_isovalues
        typ = _type01(_entry_value(entry, :type; label="type isovalue"),
            "type isovalue type")
        haskey(entries, typ) && error("duplicate type isovalue for type $typ")
        iso = _finite_real(_entry_value(entry, :isovalue; label="type isovalue"),
            "type isovalue $typ isovalue")
        source = try
            String(_entry_value(entry, :source; label="type isovalue"))
        catch err
            err isa ErrorException || rethrow()
            "unspecified"
        end
        isempty(strip(source)) && error("type isovalue $typ source must not be empty")
        entries[typ] = Dict{String,Any}(
            "type" => typ,
            "isovalue" => iso,
            "source" => source,
        )
    end
    Set(keys(entries)) == Set([0, 1]) ||
        error("constant-current provenance requires exactly one isovalue for type 0 and type 1")
    return Dict{String,Any}[entries[0], entries[1]]
end

function _legacy_type_frames(frame_origin_nm, frame_t_axis, frame_u_axis)
    for (val, name) in ((frame_origin_nm, "frame_origin_nm"),
                        (frame_t_axis, "frame_t_axis"),
                        (frame_u_axis, "frame_u_axis"))
        val === nothing &&
            error("constant-current provenance requires type_frames or complete legacy $name fields")
    end
    origin = _vec3(frame_origin_nm, "frame_origin_nm")
    t_axis = _vec3(frame_t_axis, "frame_t_axis")
    u_axis = _vec3(frame_u_axis, "frame_u_axis")
    return (
        (type=0, origin_nm=origin, t_axis=t_axis, u_axis=u_axis),
        (type=1, origin_nm=origin, t_axis=t_axis, u_axis=u_axis),
    )
end

function _cc_fields(observable, nominal_height_nm, bracket_heights_nm,
                    isovalue, type_isovalues, z_spacing_nm, crossing_policy,
                    invalid_mask, type_frames, frame_origin_nm, frame_t_axis,
                    frame_u_axis, has_templates)
    observable === "constant-current" ||
        error("observable must be \"constant-current\" to bind cc fields; got $(repr(observable))")
    for (val, name) in ((nominal_height_nm, "nominal_height_nm"),
                        (bracket_heights_nm, "bracket_heights_nm"),
                        (z_spacing_nm, "z_spacing_nm"),
                        (crossing_policy, "crossing_policy"),
                        (invalid_mask, "invalid_mask"))
        val === nothing &&
            error("constant-current provenance requires $name; it is the observable contract")
    end
    # Templates are intentionally omitted for the constant-current build, which
    # emits maps+mask only. Connected-template import is a separate step.
    has_templates &&
        error("constant-current provenance binds maps+mask, not templates; drop templates=...")

    nominal = _finite_real(nominal_height_nm, "nominal_height_nm")
    brackets = Float64[_finite_real(h, "bracket_heights_nm[$i]")
                       for (i, h) in enumerate(bracket_heights_nm)]
    typed_isovalues = _normalize_type_isovalues(type_isovalues, isovalue)
    zsp = _finite_real(z_spacing_nm, "z_spacing_nm"; positive=true)
    policy = String(crossing_policy)
    policy == CC_CROSSING_POLICY ||
        error("crossing_policy must be \"$CC_CROSSING_POLICY\" (the fixed engine convention); got $(repr(policy))")
    if type_frames === nothing
        typed = _normalize_type_frames(_legacy_type_frames(
            frame_origin_nm, frame_t_axis, frame_u_axis))
    else
        any(x -> x !== nothing, (frame_origin_nm, frame_t_axis, frame_u_axis)) &&
            error("constant-current provenance accepts either type_frames or legacy common-frame fields, not both")
        typed = _normalize_type_frames(type_frames)
    end
    return nominal, brackets, typed_isovalues, zsp, typed
end

"""
    write_qe_mold_provenance(path; provider, glcn_cube, glcnac_cube, maps,
                             sample_bias_ev, height_nm, half_nm, step_nm,
                             cube_units, templates=nothing, observable=nothing,
                             nominal_height_nm=nothing, bracket_heights_nm=nothing,
                             isovalue=nothing, type_isovalues=nothing,
                             z_spacing_nm=nothing,
                             crossing_policy=nothing, invalid_mask=nothing,
                             type_frames=nothing,
                             frame_origin_nm=nothing, frame_t_axis=nothing,
                             frame_u_axis=nothing,
                             bracket_artifacts=nothing,
                             artifact_sources=nothing) -> path

Write the bound provenance sidecar for an STM/LDOS mold. The schema stays
`stmfit-qe-mold-provenance-v1`.

* Constant-height (the historical default): pass `templates` and omit
  `observable`. The payload is byte-identical to the frozen `stm_dft_v1`
  sidecar (no observable block leaks into it).
* Constant-current (diagnostic): pass `observable="constant-current"` with the
  full observable contract — nominal and bracket heights, per-type
  calibrated/explicit isovalues, cube z spacing, fixed crossing policy, exact
  isovalue scan interval count, the nominal invalid-mask
  path, and the per-type Cu reference frames (origin + t/u axes). `templates` must be
  omitted (the build emits maps+mask only). New callers pass `type_frames` with
  exactly one type-0 and one type-1 frame entry; legacy common-frame callers may
  still pass `frame_origin_nm`, `frame_t_axis`, and `frame_u_axis`, which are
  normalized into identical typed entries. `bracket_artifacts` must bind every
  non-nominal `(height_nm, map_path, invalid_mask_path)` tuple declared by
  `bracket_heights_nm`; an empty vector is required when only the nominal height
  is declared. Transactional callers may pass `artifact_sources`, mapping each
  recorded final artifact path to the same-directory staged file whose bytes
  must be hashed. Recorded paths remain the final paths. Omitting the mapping
  preserves the direct-writer behavior.

Cube hashes, bias, and the map hash are always bound. No nA-to-cube conversion
is claimed for the constant-current observable.
"""
function write_qe_mold_provenance(path::AbstractString;
        provider::AbstractString,
        glcn_cube::AbstractString,
        glcnac_cube::AbstractString,
        maps::AbstractString,
        sample_bias_ev::Real,
        height_nm::Real,
        half_nm::Real,
        step_nm::Real,
        cube_units::AbstractString,
        templates::Union{Nothing,AbstractString}=nothing,
        observable::Union{Nothing,AbstractString}=nothing,
        nominal_height_nm::Union{Nothing,Real}=nothing,
        bracket_heights_nm::Union{Nothing,AbstractVector{<:Real}}=nothing,
        isovalue::Union{Nothing,Real}=nothing,
        type_isovalues::Union{Nothing,AbstractVector}=nothing,
        z_spacing_nm::Union{Nothing,Real}=nothing,
        crossing_policy::Union{Nothing,AbstractString}=nothing,
        invalid_mask::Union{Nothing,AbstractString}=nothing,
        type_frames::Union{Nothing,AbstractVector}=nothing,
        frame_origin_nm::Union{Nothing,AbstractVector{<:Real}}=nothing,
        frame_t_axis::Union{Nothing,AbstractVector{<:Real}}=nothing,
        frame_u_axis::Union{Nothing,AbstractVector{<:Real}}=nothing,
        bracket_artifacts::Union{Nothing,Vector}=nothing,
        isovalue_scan_intervals=DEFAULT_ISOVALUE_SCAN_INTERVALS,
        artifact_sources::Union{Nothing,AbstractDict}=nothing)
    glcn_cube_source = _require_artifact(glcn_cube, artifact_sources,
        "cannot record missing QE mold artifact")
    glcnac_cube_source = _require_artifact(glcnac_cube, artifact_sources,
        "cannot record missing QE mold artifact")
    maps_source = _require_artifact(maps, artifact_sources,
        "cannot record missing QE mold artifact")
    payload = Dict{String,Any}(
        "schema" => "stmfit-qe-mold-provenance-v1",
        "provider" => String(provider),
        "sample_bias_ev" => _finite_real(sample_bias_ev, "sample_bias_ev"),
        "height_nm" => _finite_real(height_nm, "height_nm"),
        "half_nm" => _finite_real(half_nm, "half_nm"; positive=true),
        "step_nm" => _finite_real(step_nm, "step_nm"; positive=true),
        "cube_units" => String(cube_units),
        "glcn_cube_path" => String(glcn_cube),
        "glcnac_cube_path" => String(glcnac_cube),
        "maps_path" => String(maps),
        "glcn_cube_sha256" => _sha256_file(glcn_cube_source),
        "glcnac_cube_sha256" => _sha256_file(glcnac_cube_source),
        "maps_sha256" => _sha256_file(maps_source),
    )
    if templates !== nothing
        templates_source = _require_artifact(templates, artifact_sources,
            "cannot record missing QE mold templates artifact")
        payload["templates_path"] = String(templates)
        payload["templates_sha256"] = _sha256_file(templates_source)
    end

    if observable !== nothing
        nominal, brackets, typed_isovalues, zsp, typed_frames = _cc_fields(
            observable, nominal_height_nm, bracket_heights_nm, isovalue,
            type_isovalues, z_spacing_nm, crossing_policy, invalid_mask,
            type_frames, frame_origin_nm, frame_t_axis, frame_u_axis,
            templates !== nothing)
        invalid_mask_source = _require_artifact(invalid_mask::AbstractString,
            artifact_sources, "cannot record missing constant-current invalid mask")
        payload["observable"] = String(observable)
        payload["nominal_height_nm"] = nominal
        payload["bracket_heights_nm"] = brackets
        payload["type_isovalues"] = typed_isovalues
        payload["z_spacing_nm"] = zsp
        payload["crossing_policy"] = CC_CROSSING_POLICY
        payload["isovalue_scan_intervals"] = _positive_integer(
            isovalue_scan_intervals, "isovalue_scan_intervals")
        payload["type_frames"] = typed_frames
        payload["invalid_mask_path"] = String(invalid_mask)
        payload["invalid_mask_sha256"] = _sha256_file(invalid_mask_source)
        bracket_artifacts === nothing &&
            error("constant-current provenance requires bracket_artifacts to bind every non-nominal height")
        nominal in brackets ||
            error("bracket_heights_nm must include nominal_height_nm")
        length(unique(brackets)) == length(brackets) ||
            error("bracket_heights_nm must not contain duplicates")
        entries = Dict{String,Any}[]
        artifact_heights = Float64[]
        for art in bracket_artifacts
            length(art) == 3 || error("each bracket artifact must contain height, map, and invalid mask")
            h = _finite_real(art[1], "bracket artifact height")
            mp = String(art[2]); ip = String(art[3])
            map_source = _require_artifact(mp, artifact_sources,
                "cannot record missing bracket artifact")
            invalid_source = _require_artifact(ip, artifact_sources,
                "cannot record missing bracket artifact")
            push!(artifact_heights, h)
            push!(entries, Dict{String,Any}(
                "height_nm" => h,
                "map_path" => mp,
                "map_sha256" => _sha256_file(map_source),
                "invalid_mask_path" => ip,
                "invalid_mask_sha256" => _sha256_file(invalid_source),
            ))
        end
        expected_heights = sort(filter(!=(nominal), unique(brackets)))
        sort(artifact_heights) == expected_heights ||
            error("bracket_artifacts heights must exactly match non-nominal bracket_heights_nm")
        payload["bracket_artifacts"] = sort(entries; by=entry -> entry["height_nm"])
    end

    _atomic_write(path) do io
        TOML.print(io, payload; sorted=true)
    end
    return path
end

end
