# Records and deterministic feature loading keyed by (file, lobe).

struct LobeRecord
    file::String
    lobe::Int
    amplitude::Float64
    features::Dict{String,Float64}
end

# Value equality so two LobeRecords parsed from the same TSV compare equal
# regardless of Dict identity / insertion order.
function Base.:(==)(a::LobeRecord, b::LobeRecord)
    return a.file == b.file && a.lobe == b.lobe &&
           a.amplitude == b.amplitude && a.features == b.features
end

function Base.hash(r::LobeRecord, h::UInt)
    return hash((r.file, r.lobe, r.amplitude, r.features), h)
end

function _key(file::AbstractString, lobe::Integer)
    return (basename(strip(file)), Int(lobe))
end

function _parse_lobe(row, source::AbstractString)
    haskey(row, "lobe") || error("$source missing lobe column")
    lobe = tryparse(Int, strip(row["lobe"]))
    lobe === nothing && error("$source has invalid lobe index: $(row["lobe"])")
    lobe >= 1 || error("$source lobe index must be >= 1, got $lobe")
    return lobe
end

function _parse_f(s)
    t = strip(String(s))
    (isempty(t) || t in ("NA", "NaN", "nan", "?")) && return NaN
    return parse(Float64, t)
end

# Build the feature dict from a row, enforcing the firewall. Returns the dict
# and the list of forbidden columns encountered (for an explicit error).
function _record_features(row::Dict{String,String}, path::AbstractString)
    feat = Dict{String,Float64}()
    for (k, v) in row
        key = strip(String(k))
        key in NON_FEATURE_COLUMNS && continue
        is_forbidden_feature(key) &&
            error("Feature TSV $path contains a forbidden column that would " *
                  "leak lobe-order/benchmark/truth information: $key")
        parsed = tryparse(Float64, strip(v))
        parsed === nothing && continue
        feat[key] = parsed
    end
    if !haskey(feat, "integrated") &&
       all(haskey(feat, k) for k in ("amplitude", "sigma_parallel_nm", "sigma_perp_nm"))
        feat["integrated"] = feat["amplitude"] * feat["sigma_parallel_nm"] * feat["sigma_perp_nm"]
    end
    return feat
end

function load_records(path::AbstractString)
    isfile(path) || error("Feature TSV not found: $path")
    header, rows = _read_tsv(path)
    isempty(rows) && error("No rows in $path")
    for col in ("file", "lobe", "amplitude")
        col in header || error("Feature TSV $path missing required column: $col")
    end
    records = LobeRecord[]
    seen = Set{Tuple{String,Int}}()
    for row in rows
        file = basename(strip(row["file"]))
        isempty(file) && error("Feature TSV $path has empty file value")
        lobe = _parse_lobe(row, "Feature TSV $path")
        key = _key(file, lobe)
        key in seen && error("Feature TSV $path has duplicate row for $(key[1]) lobe $(key[2])")
        push!(seen, key)
        amp = _parse_f(row["amplitude"])
        isfinite(amp) || error("Feature TSV $path has non-finite amplitude for $file lobe $lobe: $(row["amplitude"])")
        feats = _record_features(row, path)
        push!(records, LobeRecord(file, lobe, amp, feats))
    end
    sort!(records, by = r -> (r.file, r.lobe))
    return records
end

# Build an n × p feature matrix. `valid[i]` is true iff every feature for
# record i is finite (no NaN). Missing features default to NaN.
function feature_matrix(records::Vector{LobeRecord}, features::Vector{String})
    n = length(records)
    p = length(features)
    X = fill(NaN, n, p)
    valid = trues(n)
    for (i, rec) in enumerate(records), (j, fn) in enumerate(features)
        v = get(rec.features, fn, NaN)
        X[i, j] = v
        isfinite(v) || (valid[i] = false)
    end
    return X, valid
end
