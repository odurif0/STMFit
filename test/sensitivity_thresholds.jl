# sensitivity_thresholds.jl — measure how N_selected responds to the
# GCV ambiguity threshold on the chitosan benchmark.
#
# Usage:
#   julia --project=. test/sensitivity_thresholds.jl generate     # write 5 configs
#   julia --project=. test/sensitivity_thresholds.jl submit       # submit 5 batch jobs to HPC
#   julia --project=. test/sensitivity_thresholds.jl local        # run 5 batches locally (slow)
#   julia --project=. test/sensitivity_thresholds.jl compare      # diff N_selected across runs
#
# The "compare" step reads the 5 summary TSVs (one per threshold) and reports
# how many files change N_selected, identifying pivot files. A threshold is
# "robust" if ≤2 files change between adjacent threshold values (e.g. 0.04↔0.06).

using DelimitedFiles, Printf, Statistics

const THRESHOLDS = [0.03, 0.04, 0.05, 0.06, 0.08]
const BASE_CONFIG = joinpath(@__DIR__, "..", "config", "chitosan.toml")
const CONFIG_DIR = joinpath(@__DIR__, "..", "config")
const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "sensitivity")

"Generate 5 derived TOML configs varying gcv_ambiguity_rel_threshold."
function generate_configs()
    mkpath(CONFIG_DIR)
    base = read(BASE_CONFIG, String)
    for thr in THRESHOLDS
        # Replace the threshold line in the [selection] section.
        modified = replace(base,
            r"gcv_ambiguity_rel_threshold\s*=\s*[0-9.]+"m =>
                @sprintf("gcv_ambiguity_rel_threshold = %.2f", thr))
        out = joinpath(CONFIG_DIR, @sprintf("chitosan_sensitivity_%03d.toml", round(Int, thr * 1000)))
        write(out, modified)
        println("wrote $out (threshold=$thr)")
    end
end

"Submit 5 batch jobs to HPC via launch_remote.sh, one config per chunk-group."
function submit_hpc()
    mkpath(RESULTS_DIR)
    for thr in THRESHOLDS
        tag = @sprintf("%03d", round(Int, thr * 1000))
        cfg = "config/chitosan_sensitivity_$tag.toml"
        outdir = "results/sensitivity/thr$tag"
        cmd = `./hpc/launch_remote.sh --no-sync-data --config $cfg --outdir $outdir`
        @info "submitting" thr config=cfg outdir=outdir
        run(cmd)
    end
end

"Run 5 batches locally (slow — use HPC instead if available)."
function run_local()
    mkpath(RESULTS_DIR)
    data_dir = get(ENV, "STMFIT_DATA_DIR", "")
    isempty(data_dir) && error("Set STMFIT_DATA_DIR for local runs")
    for thr in THRESHOLDS
        tag = @sprintf("%03d", round(Int, thr * 1000))
        cfg = "config/chitosan_sensitivity_$tag.toml"
        outdir = "results/sensitivity/thr$tag"
        mkpath(outdir)
        @info "running locally" thr
        run(`julia -t 4 --project=. test/batch_full.jl 48 --config $cfg --data-dir $data_dir --outdir $outdir --exclude-from results/chitosan_exclude.txt`)
    end
end

"Compare N_selected across the 5 threshold runs and write a markdown report."
function compare_runs()
    mkpath(RESULTS_DIR)
    # Collect per-file N_selected for each threshold
    per_thr = Dict{Float64,Dict{String,Int}}()
    for thr in THRESHOLDS
        tag = @sprintf("%03d", round(Int, thr * 1000))
        tsv = joinpath(RESULTS_DIR, "thr$tag", "summary_overlap060_hard.tsv")
        isfile(tsv) || (@warn "missing $tsv — did you run all 5 batches?"; continue)
        d = Dict{String,Int}()
        for row in eachrow(readdlm(tsv, '\t', String, '\n'))
            row[1] == "filepath" && continue   # header
            d[row[1]] = parse(Int, row[9])     # N_selected
        end
        per_thr[thr] = d
    end
    isempty(per_thr) && error("No runs found in $RESULTS_DIR")

    # All files
    files = sort(union((keys(d) for d in values(per_thr))...))

    open(joinpath(RESULTS_DIR, "sensitivity_report.md"), "w") do io
        println(io, "# Sensitivity of N_selected to gcv_ambiguity_rel_threshold\n")
        println(io, "Benchmark: chitosan 240817 (48 files).\n")
        println(io, "Thresholds tested: $(join(THRESHOLDS, ", "))\n")

        # Pivot files: N_selected changes across thresholds
        pivots = String[]
        for f in files
            vals = [get(per_thr[t], f, -1) for t in THRESHOLDS if haskey(per_thr, t)]
            length(unique(vals)) > 1 && push!(pivots, f)
        end
        println(io, "## Summary\n")
        println(io, "- Files where N_selected is **stable** across all thresholds: $(length(files) - length(pivots)) / $(length(files))")
        println(io, "- Files where N_selected **changes** (pivots): $(length(pivots))\n")

        if !isempty(pivots)
            println(io, "## Pivot files\n")
            println(io, "| File | " * join(["thr=$t" for t in THRESHOLDS], " | ") * " |")
            println(io, "|------|" * repeat("------|" , length(THRESHOLDS)))
            for f in pivots
                row = "| $f "
                for t in THRESHOLDS
                    row *= "| $(get(per_thr[t], f, "—")) "
                end
                println(io, row * "|")
            end
        end

        # Adjacent-threshold stability (robustness criterion)
        println(io, "\n## Adjacent-threshold stability\n")
        println(io, "A threshold is robust if ≤2 files change N between adjacent values.\n")
        for i in 1:(length(THRESHOLDS) - 1)
            t1, t2 = THRESHOLDS[i], THRESHOLDS[i + 1]
            haskey(per_thr, t1) && haskey(per_thr, t2) || continue
            d1, d2 = per_thr[t1], per_thr[t2]
            nchanged = count(f -> get(d1, f, -1) != get(d2, f, -1), files)
            verdict = nchanged <= 2 ? "ROBUST" : "SENSITIVE"
            println(io, "- thr $t1 → $t2: $nchanged files change → **$verdict**")
        end
    end
    println("Report written to $(joinpath(RESULTS_DIR, "sensitivity_report.md"))")
    println("Pivot files: $(length(pivots))")
end

# Dispatch
const MODE = length(ARGS) >= 1 ? ARGS[1] : "compare"
if MODE == "generate"
    generate_configs()
elseif MODE == "submit"
    submit_hpc()
elseif MODE == "local"
    run_local()
elseif MODE == "compare"
    compare_runs()
else
    error("Unknown mode '$MODE'. Use: generate | submit | local | compare")
end
