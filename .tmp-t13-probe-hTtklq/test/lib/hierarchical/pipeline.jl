# Prediction writer (existing schema + model/version/provenance fields) and
# the top-level pipeline used by the CLI.

const BASE_PREDICTION_COLUMNS = String[
    "file", "lobe", "predicted", "confidence", "amplitude",
    "probability_1", "views_used", "invalid_reason",
]
const EXTRA_PREDICTION_COLUMNS = String["model", "model_version", "provenance_sha256"]
const PREDICTION_COLUMNS = vcat(BASE_PREDICTION_COLUMNS, EXTRA_PREDICTION_COLUMNS)

function _atomic_write(writer::Function, path::AbstractString)
    destination = abspath(path)
    parent = dirname(destination)
    mkpath(parent)
    temp_path, io = mktemp(parent; cleanup=false)
    committed = false
    try
        writer(io)
        flush(io)
        close(io)
        Base.Filesystem.rename(temp_path, destination)
        committed = true
    finally
        isopen(io) && close(io)
        !committed && ispath(temp_path) && rm(temp_path; force=true)
    end
    return nothing
end

function write_predictions(path::AbstractString, records::Vector{LobeRecord},
                           combined::CombinedPrediction, provenance_sha::String)
    n_forced = n_uncertain = 0
    _atomic_write(path) do io
        println(io, join(PREDICTION_COLUMNS, '\t'))
        for (i, rec) in enumerate(records)
            p = combined.probability_1[i]
            reason = combined.invalid_reason[i]
            if isfinite(p)
                pred = p >= 0.5 ? "1" : "0"
                conf = max(p, 1 - p)
                n_forced += 1
                println(io, join([
                    rec.file, rec.lobe, pred,
                    @sprintf("%.8f", conf),
                    @sprintf("%.10g", rec.amplitude),
                    @sprintf("%.8f", p),
                    combined.views_used[i], reason,
                    MODEL_NAME, MODEL_VERSION, provenance_sha,
                ], '\t'))
            else
                n_uncertain += 1
                println(io, join([
                    rec.file, rec.lobe, "?", "0.00000000",
                    @sprintf("%.10g", rec.amplitude), "NA",
                    combined.views_used[i], reason,
                    MODEL_NAME, MODEL_VERSION, provenance_sha,
                ], '\t'))
            end
        end
    end
    return n_forced, n_uncertain
end

# A consumed artifact is bound by a stable semantic role, normalized path
# identity, and the SHA-256 of its bytes.
struct ProvenanceInput
    role::String
    path::String
    sha256::String
end

function ProvenanceInput(role::AbstractString, path::AbstractString)
    normalized_role = strip(String(role))
    isempty(normalized_role) && error("Provenance input role is empty")
    normalized_path = normpath(abspath(String(path)))
    isfile(normalized_path) || error("Provenance input not found: $normalized_path")
    return ProvenanceInput(normalized_role, normalized_path,
                           bytes2hex(SHA.sha256(read(normalized_path))))
end

# Deterministic provenance SHA-256 over every consumed artifact + resolved
# views + fixed model options. Computed once, written on every row.
function _provenance_sha256(inputs::Vector{ProvenanceInput},
                             view_specs::Vector{Pair{String,Vector{String}}},
                             first_seed::Int, n_starts::Int,
                             cov_floor::Float64, max_iter::Int, tol::Float64)
    isempty(inputs) && error("No provenance inputs")
    roles = [input.role for input in inputs]
    length(unique(roles)) == length(roles) ||
        error("Duplicate provenance input role")
    io = IOBuffer()
    for input in sort(inputs; by=input -> input.role)
        println(io, "artifact_role=", input.role)
        println(io, "artifact_path=", input.path)
        println(io, "artifact_sha256=", input.sha256)
    end
    for (name, feats) in view_specs
        println(io, "view=", name, " features=", join(feats, ","))
    end
    println(io, "first_seed=", first_seed)
    println(io, "n_starts=", n_starts)
    println(io, "cov_floor=", cov_floor)
    println(io, "max_iter=", max_iter)
    println(io, "tol=", tol)
    println(io, "model=", MODEL_NAME, " version=", MODEL_VERSION)
    return bytes2hex(SHA.sha256(take!(io)))
end

function run_pipeline(records::Vector{LobeRecord}, inputs::Vector{ProvenanceInput},
                      out_path::AbstractString,
                      view_specs::Vector{Pair{String,Vector{String}}};
                      first_seed::Int=0, n_starts::Int=default_n_starts,
                      cov_floor::Float64=cov_floor_default,
                      max_iter::Int=200, tol::Float64=1e-6)
    isempty(view_specs) && error("No views configured")
    view_fits = ViewFit[]
    for (name, feats) in view_specs
        for fn in feats
            is_forbidden_feature(fn) &&
                error("View $name uses a forbidden feature: $fn")
            any(haskey(r.features, fn) for r in records) ||
                error("View $name uses an unavailable feature: $fn")
        end
        vf = fit_view(records, String[String(f) for f in feats], name;
                      first_seed=first_seed, n_starts=n_starts,
                      cov_floor=cov_floor, max_iter=max_iter, tol=tol)
        @printf("view %-28s rows=%d/%d features=%s one_comp=%s degenerate=%s\n",
                name, length(vf.valid_lobe_idx), length(records),
                join(feats, ","), vf.one_component_evidence, vf.degenerate)
        push!(view_fits, vf)
    end
    combined = combine_views(view_fits, length(records))
    primary_inputs = [input for input in inputs if input.role == "primary_features"]
    length(primary_inputs) == 1 ||
        error("Exactly one primary_features provenance input is required")
    features_path = only(primary_inputs).path
    provenance_sha = _provenance_sha256(inputs, view_specs,
                                        first_seed, n_starts,
                                        cov_floor, max_iter, tol)
    n_forced, n_uncertain = write_predictions(out_path, records, combined, provenance_sha)
    @printf("\nHierarchical equal-prior unit predictions\n")
    @printf("  features:   %s\n", features_path)
    @printf("  out:        %s\n", out_path)
    @printf("  files:      %d\n", length(unique(r.file for r in records)))
    @printf("  rows:       %d\n", length(records))
    @printf("  predicted:  %d\n", n_forced)
    @printf("  uncertain:  %d\n", n_uncertain)
    @printf("  views:      %s\n", join(first.(view_specs), ", "))
    @printf("  model:      %s v%d\n", MODEL_NAME, MODEL_VERSION)
    @printf("  provenance: %s\n", provenance_sha)
    return out_path
end

function run_pipeline(records::Vector{LobeRecord}, features_path::AbstractString,
                      out_path::AbstractString,
                      view_specs::Vector{Pair{String,Vector{String}}}; kwargs...)
    inputs = [ProvenanceInput("primary_features", features_path)]
    return run_pipeline(records, inputs, out_path, view_specs; kwargs...)
end

function run_pipeline(features_path::AbstractString, out_path::AbstractString,
                      view_specs::Vector{Pair{String,Vector{String}}}; kwargs...)
    records = load_records(features_path)
    return run_pipeline(records, features_path, out_path, view_specs; kwargs...)
end
