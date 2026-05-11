#!/usr/bin/env julia

using DelimitedFiles, Printf, Statistics

const SUMMARY = isdefined(Main, :SUMMARY_OVERRIDE) ? Main.SUMMARY_OVERRIDE :
    (length(ARGS) >= 1 ? ARGS[1] : "results/best_plots/summary_harmonized.tsv")

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

isfile(SUMMARY) || error("Summary not found: $SUMMARY")
data = readdlm(SUMMARY, '\t', String)
header = vec(data[1, :])
idx(name) = findfirst(==(name), header) === nothing ? error("Missing column $name") : findfirst(==(name), header)

rows = []
for i in 2:size(data, 1)
    data[i, idx("status")] == "ok" || continue
    push!(rows, data[i, :])
end

if isempty(rows)
    println("No ok rows in $SUMMARY")
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
println("summary: $SUMMARY")
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
