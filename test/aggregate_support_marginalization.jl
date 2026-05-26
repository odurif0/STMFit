#!/usr/bin/env julia

using DelimitedFiles, Printf, Statistics

function _parse_cli(args)
    inputs = String[]
    out_tsv = "results/chitosan_case_audit/support_marginalization_summary.tsv"
    top_rel = 0.05
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--inputs"
            i < length(args) || error("--inputs requires comma-separated TSV paths")
            inputs = split(args[i + 1], ","); i += 2
        elseif startswith(arg, "--inputs=")
            inputs = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--out"
            i < length(args) || error("--out requires a TSV path")
            out_tsv = args[i + 1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--top-rel"
            i < length(args) || error("--top-rel requires a numeric threshold")
            top_rel = parse(Float64, args[i + 1]); i += 2
        elseif startswith(arg, "--top-rel=")
            top_rel = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        else
            push!(inputs, arg); i += 1
        end
    end
    isempty(inputs) && error("No input TSV files. Use --inputs file1.tsv,file2.tsv or positional paths.")
    return inputs, out_tsv, top_rel
end

function _read_tsv(path)
    data, header = readdlm(path, '\t', String, header=true)
    header = vec(String.(header))
    rows = Vector{Dict{String,String}}()
    for i in axes(data, 1)
        push!(rows, Dict(header[j] => data[i, j] for j in eachindex(header)))
    end
    return rows
end

_f(x; default=NaN) = try parse(Float64, x) catch; default end
_i(x; default=0) = try parse(Int, x) catch; default end

function main()
    inputs, out_tsv, top_rel = _parse_cli(ARGS)
    by_file_support = Dict{Tuple{String,String},Vector{Dict{String,String}}}()
    for path in inputs
        for row in _read_tsv(path)
            file = row["file"]
            # Each input file is one support candidate.  Include basename so repeated
            # support lengths from different parameter triples remain separate samples.
            support_id = basename(path)
            push!(get!(by_file_support, (file, support_id), Vector{Dict{String,String}}()), row)
        end
    end

    files = sort(unique(first(k) for k in keys(by_file_support)))
    mkpath(dirname(out_tsv))
    open(out_tsv, "w") do io
        println(io, join(["file", "N", "support_samples", "win_fraction", "top_fraction",
                          "median_regret_rel", "q75_regret_rel", "mean_support_nm",
                          "robust_selected", "ambiguous"], '\t'))
        for file in files
            supports = sort([sid for (f, sid) in keys(by_file_support) if f == file])
            regrets = Dict{Int,Vector{Float64}}()
            wins = Dict{Int,Int}()
            tops = Dict{Int,Int}()
            support_lengths = Float64[]
            for sid in supports
                rows = by_file_support[(file, sid)]
                vals = Tuple{Int,Float64}[]
                for r in rows
                    n = _i(r["N"])
                    s = _f(r["eff_score_N"])
                    isfinite(s) && n > 0 && push!(vals, (n, s))
                end
                isempty(vals) && continue
                best = minimum(s for (_n, s) in vals)
                best_n = sort([n for (n, s) in vals if abs(s - best) <= 1e-15])[1]
                wins[best_n] = get(wins, best_n, 0) + 1
                push!(support_lengths, _f(first(rows)["support_2d_nm"]))
                for (n, s) in vals
                    rel = (s - best) / max(abs(best), eps(Float64))
                    push!(get!(regrets, n, Float64[]), rel)
                    if rel <= top_rel
                        tops[n] = get(tops, n, 0) + 1
                    end
                end
            end
            ns = sort(collect(keys(regrets)))
            isempty(ns) && continue
            med = Dict(n => median(regrets[n]) for n in ns)
            q75 = Dict(n => quantile(regrets[n], 0.75) for n in ns)
            selected = sort(ns; by=n -> (med[n], q75[n]))[1]
            sorted_ns = sort(ns; by=n -> (med[n], q75[n]))
            ambiguous = length(sorted_ns) > 1 && (med[sorted_ns[2]] - med[sorted_ns[1]] <= 0.02 || get(wins, selected, 0) / max(length(supports), 1) < 0.50)
            for n in ns
                println(io, join([file, n, length(supports),
                                  get(wins, n, 0) / max(length(supports), 1),
                                  get(tops, n, 0) / max(length(supports), 1),
                                  med[n], q75[n], mean(filter(isfinite, support_lengths)),
                                  n == selected, ambiguous], '\t'))
            end
        end
    end
    println("Wrote $out_tsv")
end

main()
