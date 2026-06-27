#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent
include(joinpath(@__DIR__, "lib", "qe_preflight.jl"))
using .QEPreflight: check_dir!

const DEFAULT_OUT = "hpc/qe_molds/qe_input_preflight.tsv"

function _parse_cli(args)
    dirs = String[]
    out = DEFAULT_OUT
    max_total_tasks = 8
    min_mem_mb = 0
    sequential = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--dir"
            push!(dirs, args[i+1]); i += 2
        elseif startswith(arg, "--dir=")
            push!(dirs, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--out"
            out = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--max-total-tasks"
            max_total_tasks = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--max-total-tasks=")
            max_total_tasks = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--min-mem-mb"
            min_mem_mb = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--min-mem-mb=")
            min_mem_mb = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--sequential"
            sequential = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/preflight_qe_mold_inputs.jl [options]

            Options:
              --dir PATH              QE run directory. Repeatable.
              --out PATH              Output report TSV [$(DEFAULT_OUT)]
              --max-total-tasks INT   Max simultaneous Slurm tasks (sum if parallel, max-per-job if --sequential) [8]
              --min-mem-mb INT        Optional minimum #SBATCH --mem per job [0]
              --sequential            Treat dirs as an afterok chain; task budget applies per job

            Full runs use run_qe_mold.sbatch; SCF+PP-only preliminary runs use
            run_scf_pp.sbatch and must not contain a relax step.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(dirs) && (dirs = ["qe/glcn", "qe/glcnac"])
    max_total_tasks > 0 || error("--max-total-tasks must be positive")
    min_mem_mb >= 0 || error("--min-mem-mb must be non-negative")
    return dirs, out, max_total_tasks, min_mem_mb, sequential
end

function main()
    dirs, out, max_total_tasks, min_mem_mb, sequential = _parse_cli(ARGS)
    rows = NamedTuple[]
    total_tasks = 0
    max_tasks_per_job = 0
    for dir in dirs
        tasks = check_dir!(rows, dir; min_mem_mb=min_mem_mb)
        total_tasks += tasks
        max_tasks_per_job = max(max_tasks_per_job, tasks)
    end
    effective_tasks = sequential ? max_tasks_per_job : total_tasks
    effective_tasks <= max_total_tasks || error("Simultaneous Slurm tasks $effective_tasks exceeds max $max_total_tasks")
    push!(rows, (run="all", key="total_tasks", value=string(total_tasks)))
    push!(rows, (run="all", key="max_simultaneous_tasks", value=string(effective_tasks)))
    push!(rows, (run="all", key="max_total_tasks", value=string(max_total_tasks)))
    push!(rows, (run="all", key="min_mem_mb", value=string(min_mem_mb)))
    push!(rows, (run="all", key="sequential", value=string(sequential)))
    push!(rows, (run="all", key="status", value="OK"))
    _ensure_parent(out)
    open(out, "w") do io
        println(io, "run\tkey\tvalue")
        for row in rows
            println(io, join([row.run, row.key, row.value], '\t'))
        end
    end
    println("QE mold input preflight OK")
    println("  dirs:   ", join(dirs, ", "))
    println("  tasks:  ", effective_tasks, " / ", max_total_tasks, sequential ? " simultaneous (sequential)" : " simultaneous")
    println("  report: ", out)
end

main()
