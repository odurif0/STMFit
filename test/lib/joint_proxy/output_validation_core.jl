include(joinpath(@__DIR__, "..", "script_utils.jl"))
using .ScriptUtils: _read_tsv

export validate_artifacts, run_cli

function _parse_bool(text::AbstractString, source::AbstractString)
    t = lowercase(strip(String(text)))
    t in ("true", "false") || _fail("$source has invalid boolean: $text")
    return t == "true"
end

function _parse_int(text::AbstractString, source::AbstractString)
    parsed = tryparse(Int, strip(String(text)))
    parsed === nothing && _fail("$source has invalid integer: $text")
    return parsed
end

function _parse_float(text::AbstractString, source::AbstractString)
    parsed = tryparse(Float64, strip(String(text)))
    parsed === nothing && _fail("$source has invalid numeric value: $text")
    return parsed
end

_parse_predicted(text::AbstractString, source::AbstractString) = begin
    t = strip(String(text)); t in ("0", "1", "?") || _fail("$source has invalid predicted value: $t"); t
end

_basename(text::AbstractString) = basename(strip(String(text)))
_is_prob(x::Real) = isfinite(x) && 0.0 <= x <= 1.0

function _assert_exact_header(header::Vector{String}, expected::Vector{String}, source::AbstractString)
    header == expected || _fail("$source schema mismatch")
end

function _check_forbidden_keys(node, source::AbstractString, path::String="")
    if node isa AbstractDict
        for (k, v) in pairs(node)
            key = String(k)
            key in FORBIDDEN_KEYS && _fail("forbidden expected_N/sequence/truth column in $source: $(isempty(path) ? key : path * "." * key)")
            _check_forbidden_keys(v, source, isempty(path) ? key : path * "." * key)
        end
    elseif node isa AbstractVector
        for (i, v) in enumerate(node)
            _check_forbidden_keys(v, source, "$(path)[$i]")
        end
    end
end

function _load_tsv(path::AbstractString, expected::Vector{String}, source::AbstractString)
    header, rows = _read_tsv(path)
    _assert_exact_header(String.(header), expected, source)
    isempty(rows) && _fail("$source has no rows")
    return rows
end

_manifest_table(path::AbstractString) = isfile(path) ? TOML.parsefile(path) : _fail("run_manifest.toml not found: $path")
_calibration_table(path::AbstractString) = TOML.parsefile(path)

function _require_table(t::AbstractDict, key::AbstractString, source::AbstractString)
    val = get(t, key, nothing)
    val isa AbstractDict || _fail("$source missing [$key]")
    return val
end

function _hash_file(path::AbstractString)
    isfile(path) || _fail("Config file not found: $path")
    return bytes2hex(sha256(read(path)))
end

function _string_hash(value, source::AbstractString, key::AbstractString)
    text = String(value)
    occursin(r"^[0-9a-fA-F]+$", text) || _fail("$source has invalid $key hash")
    return lowercase(text)
end

function _validate_candidate_n(rows, source::AbstractString)
    seen = Set{Tuple{String,Int}}(); by_file = Dict{String,Vector{NamedTuple}}()
    for row in rows
        file = _basename(row["file"]); n = _parse_int(row["n"], "$source file=$file"); n > 0 || _fail("$source file=$file has invalid n")
        key = (file, n); key in seen && _fail("duplicate key in candidate_n.tsv: $(key[1]) n=$(key[2])"); push!(seen, key)
        prob = _parse_float(row["probability"], "$source file=$file n=$n"); rank = _parse_int(row["rank"], "$source file=$file n=$n")
        view_count = _parse_int(row["view_count"], "$source file=$file n=$n"); residual_corr = _parse_float(row["residual_corr"], "$source file=$file n=$n")
        eff_factor = _parse_float(row["effective_view_factor"], "$source file=$file n=$n"); eff_views = _parse_float(row["effective_n_views"], "$source file=$file n=$n")
        p_eff = _parse_float(row["p_eff"], "$source file=$file n=$n"); nd_joint = _parse_int(row["nd_joint"], "$source file=$file n=$n")
        source_gcv = _parse_float(row["source_gcv"], "$source file=$file n=$n"); nrmse = _parse_float(row["joint_nrmse"], "$source file=$file n=$n")
        _is_prob(prob) || _fail("probability drift in candidate_n.tsv for $file")
        rank > 0 && view_count > 0 || _fail("$source file=$file n=$n has invalid metadata")
        -1.0 <= residual_corr <= 1.0 || _fail("$source file=$file n=$n has invalid residual_corr")
        0.0 < eff_factor <= 1.0 && eff_views > 0.0 && p_eff > 0.0 && nd_joint > 0 && source_gcv >= 0.0 && nrmse >= 0.0 ||
            _fail("$source file=$file n=$n has invalid numeric bounds")
        push!(get!(by_file, file, NamedTuple[]), (; n, prob, view_count, bwd_missing=_parse_bool(row["bwd_missing"], "$source file=$file n=$n"), is_map=_parse_bool(row["is_map"], "$source file=$file n=$n")))
    end
    for (file, rows_for_file) in by_file
        isapprox(sum(r.prob for r in rows_for_file), 1.0; atol=FLOAT_TOL, rtol=0.0) || _fail("probability drift in candidate_n.tsv for $file")
        count(r -> r.is_map, rows_for_file) <= 1 || _fail("$source file=$file has multiple MAP candidates")
        v = rows_for_file[1].view_count; m = rows_for_file[1].bwd_missing
        all(r -> r.view_count == v && r.bwd_missing == m, rows_for_file) || _fail("$source file=$file has inconsistent candidate metadata")
    end
    return by_file
end

function _validate_candidate_lobes(rows, source::AbstractString)
    seen = Set{Tuple{String,Int,Int}}(); by_candidate = Dict{Tuple{String,Int},Vector{Int}}()
    for row in rows
        file = _basename(row["file"]); n = _parse_int(row["n"], "$source file=$file"); lobe = _parse_int(row["lobe"], "$source file=$file n=$n")
        n > 0 && lobe > 0 || _fail("$source file=$file has invalid indices")
        key = (file, n, lobe); key in seen && _fail("duplicate key in candidate_lobes.tsv: $(key[1]) n=$(key[2]) lobe=$(key[3])"); push!(seen, key)
        p0 = _parse_float(row["type_probability_0"], "$source file=$file n=$n lobe=$lobe"); p1 = _parse_float(row["type_probability_1"], "$source file=$file n=$n lobe=$lobe")
        _is_prob(p0) && _is_prob(p1) && isapprox(p0 + p1, 1.0; atol=FLOAT_TOL, rtol=0.0) || _fail("type probability drift in candidate_lobes.tsv for $file n=$n lobe=$lobe")
        _parse_predicted(row["predicted_type"], "$source file=$file n=$n lobe=$lobe")
        _is_prob(_parse_float(row["confidence"], "$source file=$file n=$n lobe=$lobe")) || _fail("$source file=$file n=$n lobe=$lobe has invalid confidence")
        push!(get!(by_candidate, (file, n), Int[]), lobe)
    end
    for ((file, n), lobes) in by_candidate
        sort!(lobes); lobes == collect(1:n) || _fail("missing candidate lobe in candidate_lobes.tsv for $file n=$n")
    end
    return by_candidate
end

function _validate_predictions(rows, source::AbstractString, candidate_by_file)
    seen = Set{Tuple{String,Int}}(); by_file = Dict{String,Vector{NamedTuple}}()
    for row in rows
        file = _basename(row["file"]); lobe = _parse_int(row["lobe"], "$source file=$file"); key = (file, lobe)
        key in seen && _fail("duplicate key in predictions.tsv: $(key[1]) lobe=$(key[2])"); push!(seen, key)
        n_raw = strip(row["N_prediction"]); n = n_raw == "?" ? "?" : _parse_int(n_raw, "$source file=$file lobe=$lobe")
        n_prob = _parse_float(row["N_probability"], "$source file=$file lobe=$lobe"); n_conf = _parse_float(row["N_confidence"], "$source file=$file lobe=$lobe")
        pred = _parse_predicted(row["predicted"], "$source file=$file lobe=$lobe"); conf = _parse_float(row["confidence"], "$source file=$file lobe=$lobe")
        p0 = _parse_float(row["type_probability_0"], "$source file=$file lobe=$lobe"); p1 = _parse_float(row["type_probability_1"], "$source file=$file lobe=$lobe")
        _is_prob(n_prob) && _is_prob(n_conf) && _is_prob(conf) || _fail("$source file=$file lobe=$lobe has invalid probability bounds")
        isapprox(p0 + p1, 1.0; atol=FLOAT_TOL, rtol=0.0) || _fail("type probability drift in predictions.tsv for $file")
        n == "?" ? (pred == "?" || _fail("forced final label under uncertain N in predictions.tsv for $file")) : (pred in ("0", "1") || _fail("$source file=$file lobe=$lobe has invalid predicted value"))
        push!(get!(by_file, file, NamedTuple[]), (; lobe, n, n_prob, pred))
    end
    for (file, rows_for_file) in by_file
        nvals = unique(r.n for r in rows_for_file); length(nvals) == 1 || _fail("$source file=$file has inconsistent N_prediction values")
        n = first(nvals); n == "?" ? all(r -> r.pred == "?", rows_for_file) || _fail("forced final label under uncertain N in predictions.tsv for $file") : (
            length(rows_for_file) == Int(n) || _fail("$source file=$file has wrong candidate lobe count");
            match = findfirst(r -> r.n == n, get(candidate_by_file, file, NamedTuple[]));
            match === nothing && _fail("$source file=$file has predictions for missing candidate N");
            isapprox(rows_for_file[1].n_prob, candidate_by_file[file][match].prob; atol=FLOAT_TOL, rtol=0.0) || _fail("probability drift in predictions.tsv for $file")
        )
    end
    return by_file
end

function _validate_summaries(rows, source::AbstractString, candidate_by_file, predictions_by_file, lobes_by_candidate)
    seen = Set{String}()
    for row in rows
        file = _basename(row["file"]); file in seen && _fail("duplicate key in chain_summary.tsv: $file"); push!(seen, file)
        n_raw = strip(row["N_prediction"]); n = n_raw == "?" ? "?" : _parse_int(n_raw, "$source file=$file")
        n_prob = _parse_float(row["N_probability"], "$source file=$file"); n_conf = _parse_float(row["N_confidence"], "$source file=$file")
        n_abstained = _parse_bool(row["N_abstained"], "$source file=$file"); candidate_count = _parse_int(row["candidate_count"], "$source file=$file")
        lobe_count = _parse_int(row["lobe_count"], "$source file=$file"); view_count = _parse_int(row["view_count"], "$source file=$file")
        bwd_missing = _parse_bool(row["bwd_missing"], "$source file=$file")
        _is_prob(n_prob) && _is_prob(n_conf) || _fail("$source file=$file has invalid N probability bounds")
        candidate_rows = get(candidate_by_file, file, NamedTuple[]); prediction_rows = get(predictions_by_file, file, NamedTuple[])
        candidate_count == length(candidate_rows) || _fail("$source file=$file has inconsistent candidate_count")
        expected_lobes = sum((length(v) for ((f, _), v) in lobes_by_candidate if f == file); init=0)
        lobe_count == expected_lobes || _fail("$source file=$file has inconsistent lobe_count")
        if isempty(candidate_rows)
            n_abstained && candidate_count == 0 && lobe_count == 0 && isempty(prediction_rows) ||
                _fail("$source file=$file has invalid empty-candidate abstention")
            continue
        end
        first(candidate_rows).view_count == view_count || _fail("$source file=$file has inconsistent view_count")
        first(candidate_rows).bwd_missing == bwd_missing || _fail("$source file=$file has inconsistent bwd_missing")
        n_abstained == (n == "?") || _fail("$source file=$file has inconsistent abstention semantics")
        n == "?" ? all(r -> r.pred == "?", prediction_rows) || _fail("forced final label under uncertain N in predictions.tsv for $file") :
            length(prediction_rows) == n || _fail("$source file=$file has wrong candidate lobe count")
    end
    return nothing
end

function _validate_manifest(manifest::AbstractDict, calibration::AbstractDict, calibration_path::AbstractString)
    _check_forbidden_keys(manifest, "run_manifest.toml"); _check_forbidden_keys(calibration, "joint_proxy_calibration.toml")
    prov = _require_table(manifest, "provenance", "run_manifest.toml"); inputs = _require_table(manifest, "inputs", "run_manifest.toml"); outputs = _require_table(manifest, "outputs", "run_manifest.toml")
    cprov = _require_table(calibration, "provenance", "joint_proxy_calibration.toml"); _require_table(calibration, "count", "joint_proxy_calibration.toml"); _require_table(calibration, "type", "joint_proxy_calibration.toml")
    cfg = _hash_file(String(inputs["config"]))
    cfg == _string_hash(get(prov, "config_sha256", ""), "run_manifest.toml", "config") == _string_hash(get(cprov, "config_sha256", ""), "joint_proxy_calibration.toml", "config") || _fail("config hash mismatch")
    _string_hash(get(prov, "source_sha256", ""), "run_manifest.toml", "source") == _string_hash(get(cprov, "source_sha256", ""), "joint_proxy_calibration.toml", "source") || _fail("source hash mismatch")
    _string_hash(get(prov, "payload_sha256", ""), "run_manifest.toml", "payload") == _string_hash(get(cprov, "payload_sha256", ""), "joint_proxy_calibration.toml", "payload") || _fail("payload hash mismatch")
    basename(String(get(prov, "calibration_path", ""))) == basename(String(calibration_path)) || _fail("run_manifest.toml calibration_path mismatch")
    files = get(inputs, "files", String[]); files isa AbstractVector || _fail("run_manifest.toml inputs.files has invalid schema")
    manifest_files = sort(String.(files)); manifest_files == sort(unique(manifest_files)) || _fail("run_manifest.toml inputs.files must be unique")
    _parse_int(string(get(outputs, "candidate_n_rows", "")), "run_manifest.toml") > 0 || _fail("run_manifest.toml outputs must be positive")
    _parse_int(string(get(outputs, "candidate_lobes_rows", "")), "run_manifest.toml") > 0 || _fail("run_manifest.toml outputs must be positive")
    _parse_int(string(get(outputs, "predictions_rows", "")), "run_manifest.toml") > 0 || _fail("run_manifest.toml outputs must be positive")
    _parse_int(string(get(outputs, "chain_summary_rows", "")), "run_manifest.toml") > 0 || _fail("run_manifest.toml outputs must be positive")
    return manifest_files
end

function validate_artifacts(opt::ValidationOptions)
    artifact_dir = abspath(opt.artifacts); calibration = _calibration_table(opt.calibration); manifest = _manifest_table(joinpath(artifact_dir, "run_manifest.toml"))
    manifest_files = _validate_manifest(manifest, calibration, opt.calibration)
    candidate_n_rows = _load_tsv(joinpath(artifact_dir, "candidate_n.tsv"), CANDIDATE_N_FIELDS, "candidate_n.tsv")
    candidate_lobe_rows = _load_tsv(joinpath(artifact_dir, "candidate_lobes.tsv"), CANDIDATE_LOBE_FIELDS, "candidate_lobes.tsv")
    prediction_rows = _load_tsv(joinpath(artifact_dir, "predictions.tsv"), PREDICTION_FIELDS, "predictions.tsv")
    summary_rows = _load_tsv(joinpath(artifact_dir, "chain_summary.tsv"), SUMMARY_FIELDS, "chain_summary.tsv")
    candidate_by_file = _validate_candidate_n(candidate_n_rows, "candidate_n.tsv")
    lobes_by_candidate = _validate_candidate_lobes(candidate_lobe_rows, "candidate_lobes.tsv")
    predictions_by_file = _validate_predictions(prediction_rows, "predictions.tsv", candidate_by_file)
    _validate_summaries(summary_rows, "chain_summary.tsv", candidate_by_file, predictions_by_file, lobes_by_candidate)
    summary_files = sort([String(row["file"]) for row in summary_rows])
    summary_files == manifest_files || _fail("run_manifest.toml files do not match artifact files")
    return ValidationReport(files=length(summary_files), candidate_rows=length(candidate_n_rows), candidate_lobe_rows=length(candidate_lobe_rows), prediction_rows=length(prediction_rows), summary_rows=length(summary_rows))
end

function run_cli(args::AbstractVector{<:AbstractString}=String[])
    report = validate_artifacts(_parse_cli(args))
    println("joint proxy validation ok")
    println("  files: ", report.files)
    println("  candidate_n_rows: ", report.candidate_rows)
    println("  candidate_lobe_rows: ", report.candidate_lobe_rows)
    println("  prediction_rows: ", report.prediction_rows)
    println("  summary_rows: ", report.summary_rows)
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cli(ARGS)
end
