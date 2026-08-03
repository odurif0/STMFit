#!/usr/bin/env julia

using Printf
using Statistics

function parse_args(args)
    values = Dict{String,String}()
    i = 1
    while i <= length(args)
        flag = String(args[i])
        flag in ("--synthetic", "--real", "--real-summary", "--out") ||
            error("unknown argument: $flag")
        i < length(args) || error("$flag requires a value")
        values[flag] = String(args[i + 1])
        i += 2
    end
    all(haskey(values, flag) for flag in ("--synthetic", "--real", "--real-summary", "--out")) ||
        error("--synthetic, --real, --real-summary and --out are required")
    return values
end

function read_rows(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("empty table: $path")
    header = split(lines[1], '\t')
    return [Dict(header .=> split(line, '\t')) for line in lines[2:end] if !isempty(strip(line))]
end

parse_float(row, field) = parse(Float64, row[field])
parse_int(row, field) = parse(Int, row[field])

function synthetic_metrics(rows)
    output = NamedTuple[]
    groups = sort(unique((row["injected_family"], row["ablation"], row["noise_sigma"]) for row in rows))
    for (injected, ablation, noise) in groups
        subset = filter(row -> (row["injected_family"], row["ablation"], row["noise_sigma"]) ==
            (injected, ablation, noise), rows)
        recovery = mean(row["correct"] == "true" for row in subset)
        margin = mean(parse_float(row, "true_margin") for row in subset)
        push!(output, (; domain="synthetic", patch_kind="template", injected_family=injected,
            ablation, noise_sigma=noise, metric="recovery", value=recovery,
            n=length(subset), details="synthetic truth only"))
        push!(output, (; domain="synthetic", patch_kind="template", injected_family=injected,
            ablation, noise_sigma=noise, metric="mean_true_margin", value=margin,
            n=length(subset), details="p_true-p_other"))
    end
    return output
end

function real_metrics(rows, summaries)
    output = NamedTuple[]
    groups = sort(unique((row["patch_kind"], row["ablation"]) for row in rows))
    for (patch_kind, ablation) in groups
        subset = filter(row -> (row["patch_kind"], row["ablation"]) == (patch_kind, ablation), rows)
        argmax1 = mean(parse_int(row, "argmax_type") == 1 for row in subset)
        mean_p1 = mean(parse_float(row, "p1") for row in subset)
        confidence = mean(parse_float(row, "confidence") for row in subset)
        for (metric, value, details) in (
            ("argmax1_fraction", argmax1, "label-free raw argmax"),
            ("mean_p1", mean_p1, "uncalibrated family posterior"),
            ("mean_confidence", confidence, "max(p0,p1)"),
        )
            push!(output, (; domain="real", patch_kind, injected_family="none", ablation,
                noise_sigma="none", metric, value, n=length(subset), details))
        end
    end
    for patch_kind in sort(unique(row["patch_kind"] for row in summaries if row["status"] == "ok"))
        subset = filter(row -> row["status"] == "ok" && row["patch_kind"] == patch_kind &&
            row["ablation"] == "combined", summaries)
        isempty(subset) && continue
        masses = [parse_float(row, "stm_dft_v1_mass") for row in subset]
        push!(output, (; domain="real", patch_kind, injected_family="none", ablation="combined",
            noise_sigma="none", metric="mean_stm_dft_v1_mass", value=mean(masses),
            n=length(subset), details="global-state posterior mass"))
    end
    for status in sort(unique(row["status"] for row in summaries))
        subset = filter(row -> row["status"] == status, summaries)
        push!(output, (; domain="real", patch_kind="all", injected_family="none", ablation="all",
            noise_sigma="none", metric="status_$status", value=length(unique(row["file"] for row in subset)),
            n=length(subset), details="unique files"))
    end
    required = ("file", "n", "lobe", "patch_kind", "ablation", "argmax_type")
    if all(row -> all(haskey(row, field) for field in required), rows)
        lookup = Dict((row["file"], row["n"], row["lobe"], row["patch_kind"], row["ablation"]) =>
            parse_int(row, "argmax_type") for row in rows)
        keys3 = unique((row["file"], row["n"], row["lobe"]) for row in rows)
        patch_kinds = ("image", "raw", "residual", "neg_residual")
        for ablation in sort(unique(row["ablation"] for row in rows))
            valid = [key for key in keys3 if all(haskey(lookup,
                (key[1], key[2], key[3], patch_kind, ablation)) for patch_kind in patch_kinds)]
            isempty(valid) && continue
            values(key) = [lookup[(key[1], key[2], key[3], patch_kind, ablation)]
                for patch_kind in patch_kinds]
            metrics = (
                ("all_patch_argmax_agreement", mean(length(unique(values(key))) == 1 for key in valid)),
                ("image_residual_agreement", mean(values(key)[1] == values(key)[3] for key in valid)),
                ("raw_residual_agreement", mean(values(key)[2] == values(key)[3] for key in valid)),
            )
            for (metric, value) in metrics
                push!(output, (; domain="real", patch_kind="cross_representation",
                    injected_family="none", ablation, noise_sigma="none", metric, value,
                    n=length(valid), details="label-free argmax stability"))
            end
        end
    end
    return output
end

function write_metrics(path::AbstractString, rows)
    mkpath(dirname(path))
    fields = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for row in rows
            values = [value isa AbstractFloat ? @sprintf("%.9g", value) : string(value)
                for value in (getproperty(row, field) for field in fields)]
            println(io, join(values, '\t'))
        end
    end
end

function main(args)
    options = parse_args(args)
    synthetic = read_rows(options["--synthetic"])
    real = read_rows(options["--real"])
    summaries = read_rows(options["--real-summary"])
    metrics = vcat(synthetic_metrics(synthetic), real_metrics(real, summaries))
    write_metrics(options["--out"], metrics)
    println("wrote $(options["--out"]) ($(length(metrics)) metrics)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
