# Print classification stats, exact/tolerant agreement rates, N counts from summary TSV(s).
# Usage:
#   julia --project=. test/summarize.jl results/best_plots/summary_overlap060_hard_combined.tsv
#   julia --project=. test/summarize.jl results/best_plots/summary_overlap060_hard_chunk0*of04.tsv

using DelimitedFiles, Printf, Statistics

summaries = isempty(ARGS) ? ["results/best_plots/summary_harmonized.tsv"] : ARGS
files = String[]
for s in summaries
    if isfile(s)
        push!(files, s)
    elseif endswith(s, ".tsv")
        @warn "$s not found, skipping"
    else
        # Glob-like: treat as prefix and find matching .tsv files in same dir
        dir, pat = dirname(s), basename(s)
        matched = sort([f for f in readdir(dir; join=true) if endswith(f, ".tsv") && startswith(basename(f), pat)])
        isempty(matched) && (@warn "no .tsv files matching '$s'"; continue)
        append!(files, matched)
    end
end
isempty(files) && error("No summary files found")

function parse_set(s)
    s = strip(String(s))
    s == "{}" && return Int[]
    startswith(s, "{") && endswith(s, "}") || return Int[]
    body = s[2:end-1]
    isempty(body) && return Int[]
    return parse.(Int, split(body, ","))
end

function counts(vals)
    d = Dict{Any,Int}()
    for v in vals
        d[v] = get(d, v, 0) + 1
    end
    return sort(collect(d); by=x -> string(x[1]))
end

# Read first file (with header)
global data = readdlm(files[1], '\t', String)
# Append data rows from remaining files
for i in 2:length(files)
    global data
    f = files[i]
    isfile(f) || (println("WARN: $f not found"); continue)
    d2 = readdlm(f, '\t', String)
    data = vcat(data, d2[2:end, :])
end

header = vec(data[1, :])
ncol = length(header)

function idx(name::String)
    i = findfirst(==(name), header)
    i === nothing && error("Missing column: $name")
    return i
end

rows = []
for i in 2:size(data, 1)
    data[i, idx("status")] == "ok" || continue
    push!(rows, data[i, :])
end

if isempty(rows)
    println("No ok rows across $(length(files)) file(s)")
    exit()
end

Nell = [parse(Int, r[idx("N_ell")]) for r in rows]
Ncirc = [parse(Int, r[idx("N_circ")]) for r in rows]
N1d = [parse(Int, r[idx("N_1D")]) for r in rows]
classes = [r[idx("classification")] for r in rows]
common10 = [parse_set(r[idx("common_N_10")]) for r in rows]
commonh = [parse_set(r[idx("common_N_hybrid")]) for r in rows]
support_mismatch_ell = [parse(Float64, r[idx("support_mismatch_ell")]) for r in rows]

dc_e = Ncirc .- Nell
d1_e = N1d .- Nell
d1_c = N1d .- Ncirc

println("=== Enriched summary stats ===")
println("files: $(join(basename.(files), ", "))")
println("n ok = $(length(rows))")
println("\nClassification counts: ", counts(classes))
println("\nExact agreements:")
println("  ell=circ: $(count(==(0), dc_e))/$(length(rows))")
println("  1D=ell:   $(count(==(0), d1_e))/$(length(rows))")
println("  1D=circ:  $(count(==(0), d1_c))/$(length(rows))")
println("  all:      $(count(i -> Nell[i] == Ncirc[i] == N1d[i], eachindex(rows)))/$(length(rows))")
println("\nTolerant agreements:")
println("  common ΔBIC<=10:     $(count(!isempty, common10))/$(length(rows))")
println("  common hybrid:       $(count(!isempty, commonh))/$(length(rows))")
@printf("\nMean |Δ| circ-ell=%.2f, 1D-ell=%.2f, 1D-circ=%.2f\n", mean(abs.(dc_e)), mean(abs.(d1_e)), mean(abs.(d1_c)))
@printf("Mean support mismatch ell=%.1f%%\n", 100mean(support_mismatch_ell))
println("\nN counts:")
println("  ell:  ", counts(Nell))
println("  circ: ", counts(Ncirc))
println("  1D:   ", counts(N1d))
