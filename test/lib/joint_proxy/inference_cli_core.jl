function _read_tsv(path::AbstractString)
    lines = readlines(path)
    data = [l for l in lines if !isempty(strip(l)) && !startswith(strip(l), '#')]
    isempty(data) && return Dict{String,String}[]
    header = split(first(data), '\t'; keepempty=true)
    rows = Dict{String,String}[]
    for line in data[2:end]
        vals = split(line, '\t'; keepempty=true)
        row = Dict{String,String}()
        for (i, h) in enumerate(header)
            row[h] = i <= length(vals) ? vals[i] : ""
        end
        push!(rows, row)
    end
    return rows
end

function _write_tsv(path::AbstractString, fields::Vector{String}, rows)
    _ensure_parent(path)
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for row in rows
            println(io, join((_tsv_fmt(_field_value(row, f)) for f in fields), '\t'))
        end
    end
end

function _write_toml(io::IO, name::String, value)
    println(io, "[$name]")
    if value isa AbstractDict
        for (k, v) in sort(collect(value); by=first)
            if v isa AbstractDict || v isa NamedTuple
                _write_toml(io, "$name.$(String(k))", v)
            elseif v isa AbstractVector && !isempty(v) && (first(v) isa AbstractDict || first(v) isa NamedTuple)
                for item in v
                    println(io, "[[$name.$(String(k))]]")
                    for (ik, iv) in pairs(item)
                        println(io, "$(String(ik)) = $(_fmt(iv))")
                    end
                end
            else
                println(io, "$(String(k)) = $(_fmt(v))")
            end
        end
    else
        for (k, v) in pairs(value)
            println(io, "$(String(k)) = $(_fmt(v))")
        end
    end
    println(io)
end

_source_hash(paths::Vector{String}) = bytes2hex(sha256(codeunits(join([basename(p) * ":" * (isfile(p) ? sha256_file(p) : "missing") for p in paths], '\n'))))

_calibration_source_sha256() = _source_hash(_calibration_source_files())

function _validate_bare_file(name::AbstractString, source::String)
    t = strip(String(name))
    isempty(t) && throw(InferenceCliError("empty file name in $source"))
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*\.sxm$", t) ||
        throw(InferenceCliError("invalid bare SXM name in $source: $t"))
    occursin(r"[/\\\t=,\s]", t) && throw(InferenceCliError("non-bare SXM name in $source: $t"))
    return t
end

function _parse_files_arg(arg::AbstractString)
    files = [_validate_bare_file(tok, "--files") for tok in split(String(arg), ','; keepempty=true)]
    isempty(files) && throw(InferenceCliError("--files is empty"))
    length(unique(files)) == length(files) || throw(InferenceCliError("duplicate file names in --files"))
    return files
end

function _parse_files_from(path::AbstractString)
    isfile(path) || throw(InferenceCliError("files-from list not found: $path"))
    files = String[]
    for (line_no, raw) in enumerate(readlines(path))
        s = strip(raw)
        isempty(s) && continue
        startswith(s, "#") && continue
        tokens = split(s)
        length(tokens) == 1 || throw(InferenceCliError("files-from must contain exactly one bare SXM name per line: $path:$line_no"))
        push!(files, _validate_bare_file(tokens[1], "$path:$line_no"))
    end
    isempty(files) && throw(InferenceCliError("files-from list is empty: $path"))
    length(unique(files)) == length(files) || throw(InferenceCliError("duplicate file names in files-from list"))
    return files
end

function _parse_chunk(token::AbstractString)
    m = match(r"^(\d+)/(\d+)$", strip(String(token)))
    m === nothing && throw(InferenceCliError("--chunk must be I/N"))
    i = parse(Int, m.captures[1]); n = parse(Int, m.captures[2])
    n > 0 && i >= 1 && i <= n || throw(InferenceCliError("invalid chunk $i/$n"))
    return (i, n)
end

function _required_value(args::AbstractVector{<:AbstractString}, i::Int, flag::String)
    i < length(args) || throw(InferenceCliError("$flag requires a value"))
    return String(args[i + 1])
end

function parse_cli(args::AbstractVector{<:AbstractString})::CliOptions
    config = ""; data_dir = ""; files = String[]; files_from = nothing; calibration = nothing
    outdir = ""; chunk = nothing; uncalibrated = false
    i = 1
    while i <= length(args)
        a = String(args[i])
        if a == "--config"
            config = _required_value(args, i, a); i += 2
        elseif a == "--data-dir"
            data_dir = _required_value(args, i, a); i += 2
        elseif a == "--files"
            files_from === nothing || throw(InferenceCliError("use exactly one of --files or --files-from"))
            files = _parse_files_arg(_required_value(args, i, a)); i += 2
        elseif a == "--files-from"
            isempty(files) || throw(InferenceCliError("use exactly one of --files or --files-from"))
            files_from = _required_value(args, i, a); i += 2
        elseif a == "--calibration"
            calibration === nothing || throw(InferenceCliError("use exactly one of --calibration or --uncalibrated"))
            calibration = _required_value(args, i, a); i += 2
        elseif a == "--uncalibrated"
            uncalibrated && throw(InferenceCliError("duplicate --uncalibrated"))
            uncalibrated = true; i += 1
        elseif a == "--outdir"
            outdir = _required_value(args, i, a); i += 2
        elseif a == "--chunk"
            chunk = _parse_chunk(_required_value(args, i, a)); i += 2
        else
            throw(InferenceCliError("unknown argument: $a"))
        end
    end
    isempty(config) && throw(InferenceCliError("--config is required"))
    isempty(data_dir) && throw(InferenceCliError("--data-dir is required"))
    isempty(outdir) && throw(InferenceCliError("--outdir is required"))
    (isempty(files) != (files_from === nothing)) || throw(InferenceCliError("use exactly one of --files or --files-from"))
    uncalibrated || calibration !== nothing || throw(InferenceCliError("--calibration is required unless --uncalibrated"))
    uncalibrated && calibration !== nothing && throw(InferenceCliError("--uncalibrated cannot be combined with --calibration"))
    files = files_from === nothing ? files : _parse_files_from(files_from)
    return CliOptions(config, data_dir, files, files_from, calibration, outdir, chunk, uncalibrated)
end

function _load_calibration(path::AbstractString)
    isfile(path) || throw(InferenceCliError("calibration not found: $path"))
    t = TOML.parsefile(path)
    prov = get(t, "provenance", nothing)
    count = get(t, "count", nothing)
    type = get(t, "type", nothing)
    prov isa AbstractDict || throw(InferenceCliError("calibration missing [provenance]"))
    count isa AbstractDict || throw(InferenceCliError("calibration missing [count]"))
    type isa AbstractDict || throw(InferenceCliError("calibration missing [type]"))
    return InferenceCalibration(Float64(count["temperature"]), Float64(count["confidence_threshold"]),
        Float64(type["temperature"]), Float64(type["confidence_threshold"]),
        String(prov["config_sha256"]), String(prov["source_sha256"]), String(prov["payload_sha256"]), path)
end

function _bundle_for(opts::CliOptions)
    isfile(opts.config) || throw(InferenceCliError("config file not found: $(opts.config)"))
    config_hash = sha256_file(opts.config)
    registry = load_registry(opts.config)
    source_hash = _calibration_source_sha256()
    payload_hash = registry.payload_sha256
    calibration = opts.uncalibrated ? nothing : _load_calibration(opts.calibration::String)
    if calibration !== nothing
        calibration.config_sha256 == config_hash || throw(InferenceCliError("config hash mismatch"))
        calibration.source_sha256 == source_hash || throw(InferenceCliError("source hash mismatch"))
        calibration.payload_sha256 == payload_hash || throw(InferenceCliError("payload hash mismatch"))
    end
    return InferenceBundle(config_hash, source_hash, payload_hash, registry, calibration)
end

_load_bundle_for_test(config::AbstractString) = _bundle_for(CliOptions(config=config, data_dir=".", files=String[], files_from=nothing,
    calibration=nothing, outdir=".", chunk=nothing, uncalibrated=true))

function _chunk_files(files::Vector{String}, chunk::Union{Nothing,Tuple{Int,Int}})
    isempty(files) && throw(InferenceCliError("no files selected"))
    chunk === nothing && return sort(files)
    i, n = chunk
    s = fld((i - 1) * length(files), n) + 1
    e = fld(i * length(files), n)
    return s > e ? String[] : sort(files)[s:e]
end

function _require_files(data_dir::AbstractString, files::Vector{String})
    isdir(data_dir) || throw(InferenceCliError("data-dir not found: $data_dir"))
    resolved = String[]
    for f in files
        path = joinpath(data_dir, f)
        isfile(path) || throw(InferenceCliError("data file not found: $path"))
        push!(resolved, path)
    end
    return resolved
end

function _sort_artifacts(art::InferenceArtifacts)
    cn = sort(art.candidate_n_rows; by=r -> (r.file, r.n))
    cl = sort(art.candidate_lobes_rows; by=r -> (r.file, r.n, r.lobe))
    pr = sort(art.predictions_rows; by=r -> (r.file, r.lobe))
    cs = sort(art.chain_summary_rows; by=r -> r.file)
    return InferenceArtifacts(cn, cl, pr, cs, art.manifest)
end

function _write_outputs(outdir::AbstractString, art::InferenceArtifacts)
    mkpath(outdir)
    _write_tsv(joinpath(outdir, "candidate_n.tsv"), CANDIDATE_N_FIELDS, art.candidate_n_rows)
    _write_tsv(joinpath(outdir, "candidate_lobes.tsv"), CANDIDATE_LOBE_FIELDS, art.candidate_lobes_rows)
    _write_tsv(joinpath(outdir, "predictions.tsv"), PREDICTION_FIELDS, art.predictions_rows)
    _write_tsv(joinpath(outdir, "chain_summary.tsv"), SUMMARY_FIELDS, art.chain_summary_rows)
    open(joinpath(outdir, "run_manifest.toml"), "w") do io
        man = art.manifest
        for (name, value) in sort(collect(pairs(man)); by=first)
            _write_toml(io, String(name), value)
        end
    end
end

function run_cli(args::AbstractVector{<:AbstractString}=String[]; adapter=_real_inference_artifacts)
    opts = parse_cli(args)
    bundle = _bundle_for(opts)
    files = _chunk_files(opts.files, opts.chunk)
    _require_files(opts.data_dir, files)
    art = _sort_artifacts(adapter(opts, bundle, files))
    _write_outputs(opts.outdir, art)
    return opts.outdir
end
