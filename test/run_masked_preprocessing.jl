#!/usr/bin/env julia
# Four saved-geometry cases, one preprocessing diagnostic command each.
module MaskedPreprocessingDriver
using TOML

const FROZEN_FILES = ("260115_016.sxm", "251206_013.sxm", "260215_022.sxm", "251206_034.sxm")
const PATH_OPTIONS = ("--input-dir", "--config", "--acquisition-settings", "--settings", "--outdir")
const HELP = """
Usage (Julia 1.13, Slurm for execution):
  julia --startup-file=no --threads=1 --project=. test/run_masked_preprocessing.jl \
    --input-dir INPUT --config PHYSICAL.toml --acquisition-settings ACQ.toml \
    --settings MASK.toml --outdir NEW_DIR [--dry-run]
All five paths are required. INPUT contains diagnostic_sample.tsv,
selected_summary.tsv, base_geometry.tsv and raw/<file> for exactly the four
original frozen files. The sample's file order is retained, never re-ranked.
Saved N_selected is used only to check geometry coverage, not as an expected N.
Original refined_policy must agree between sample and selected summary.
Four independent one-thread commands run only the masked-preprocessing CLI.
No counting, classifier, feature-extraction or grading stage is invoked.
Per-case outputs are opaque to this driver. Logs and completed stages.tsv rows
are retained, including failures; any failed case makes the overall run fail.
The output path must not exist, including a dangling symlink; its parent must
exist. No overwrite or retry. Execution requires Slurm with at least four CPUs.
--dry-run checks only paths, TOML sections and TSV coverage, then prints commands.
It does not read SXM pixels, invoke child commands, or create any output.
"""

function options(args)
    args == ["--help"] && return nothing
    d = Dict{String,String}(); dry = false; i = 1
    while i <= length(args)
        key = args[i]
        if key == "--dry-run"
            dry && error("Duplicate --dry-run")
            dry = true; i += 1; continue
        end
        key in PATH_OPTIONS || error("Unknown option: $key")
        haskey(d, key) && error("Duplicate option: $key")
        i < length(args) && !startswith(args[i+1], "--") || error("Missing value: $key")
        isempty(strip(args[i+1])) && error("Empty value: $key")
        d[key] = abspath(args[i+1]); i += 2
    end
    all(k -> haskey(d, k), PATH_OPTIONS) || error("All five path arguments are required; use --help")
    return (; d, dry)
end

# Metadata only. Never load an SXM reader or a fitting module in the driver.
function table(path, required)
    lines = readlines(path)
    isempty(lines) && error("Empty table: $path")
    header = String.(strip.(split(first(lines), '\t'; keepempty=true)))
    length(unique(header)) == length(header) || error("Duplicate columns: $path")
    any(isempty, header) && error("Empty column name: $path")
    all(k -> k in header, required) || error("Required columns missing: $path")
    forbidden = ("truth", "sequence", "expected_n", "target_n", "control_sequence",
                 "label", "predicted", "class_count", "class_counts", "composition_prior", "count_prior")
    any(k -> lowercase(k) in forbidden, header) && error("External labels/priors are not diagnostic inputs: $path")
    rows = Dict{String,String}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = String.(strip.(split(line, '\t'; keepempty=true)))
        length(values) == length(header) || error("Malformed table: $path")
        push!(rows, Dict(zip(header, values)))
    end
    return rows
end

function check_output(out)
    (ispath(out) || islink(out)) && error("Output already exists: $out")
    isdir(dirname(out)) || error("Output parent directory must exist: $(dirname(out))")
    return nothing
end

function require_sections(path, names)
    config = TOML.parsefile(path)
    all(k -> get(config, k, nothing) isa AbstractDict, names) || error("Missing TOML sections in $path: $names")
    return nothing
end

function prepare(args; project=dirname(@__DIR__))
    o = options(args); o === nothing && return nothing
    input = o.d["--input-dir"]; out = o.d["--outdir"]
    check_output(out)
    paths = (config=o.d["--config"], acquisition_settings=o.d["--acquisition-settings"],
             settings=o.d["--settings"], sample=joinpath(input, "diagnostic_sample.tsv"),
             selected=joinpath(input, "selected_summary.tsv"), geometry=joinpath(input, "base_geometry.tsv"))
    all(isfile, values(paths)) || error("Missing required input file")
    require_sections(paths.config, ("model", "selection", "preprocessing"))
    require_sections(paths.acquisition_settings, ("acquisition_noise",))
    require_sections(paths.settings, ("masked_robust_preprocessing",))
    sample = table(paths.sample, ["file", "N_selected", "refined_policy"])
    length(sample) == 4 || error("Exactly four frozen diagnostic files are required")
    names = getindex.(sample, "file")
    all(f -> f == basename(f) && endswith(f, ".sxm"), names) || error("Manifest requires SXM basenames")
    length(unique(names)) == 4 || error("Duplicate diagnostic file")
    Set(names) == Set(FROZEN_FILES) || error("Manifest must contain the original frozen four-file set")
    selected = table(paths.selected, ["filepath", "N_selected", "refined_policy"])
    keys = basename.(getindex.(selected, "filepath"))
    all(!isempty, keys) || error("Empty selected-summary filename")
    length(unique(keys)) == length(keys) || error("Duplicate selected-summary basename")
    byfile = Dict(zip(keys, selected))
    geometry = table(paths.geometry, ["file", "N", "lobe"])
    root = abspath(project)
    prefix = [joinpath(Sys.BINDIR, Base.julia_exename()), "--startup-file=no", "--threads=1", "--project=$root"]
    cases = NamedTuple[]
    for row in sample
        file = row["file"]
        raw = joinpath(input, "raw", file)
        isfile(raw) || error("Missing raw file: $file") # Existence only; no pixel read.
        haskey(byfile, file) || error("No saved selection for $file")
        saved = byfile[file]
        n = tryparse(Int, saved["N_selected"])
        n !== nothing && n > 0 || error("Invalid saved N_selected: $file")
        tryparse(Int, row["N_selected"]) == n || error("Sample/selected geometry count mismatch: $file")
        !isempty(saved["refined_policy"]) && row["refined_policy"] == saved["refined_policy"] ||
            error("Original refined_policy missing or inconsistent: $file")
        gs = filter(g -> basename(g["file"]) == file, geometry)
        length(gs) == n && all(g -> tryparse(Int, g["N"]) == n, gs) &&
            Set(tryparse(Int, g["lobe"]) for g in gs) == Set(1:n) ||
            error("Saved geometry/count coverage mismatch: $file")
        stem = splitext(file)[1]
        case_out = joinpath(out, stem)
        argv = [joinpath(root, "test", "diagnose_masked_preprocessing_real.jl"),
                "--file", raw, "--geometry", paths.geometry, "--selected-summary", paths.selected,
                "--config", paths.config, "--acquisition-settings", paths.acquisition_settings,
                "--settings", paths.settings, "--outdir", case_out]
        push!(cases, (file=file, stem=stem, out=case_out, cmd=Cmd(vcat(prefix, argv))))
    end
    return (; o, paths, cases, out)
end

# The executor seam is for synthetic tests. Production always invokes the Cmd
# directly, never through a shell; spaces and quote characters stay in one arg.
function run_command(cmd, io)
    process = run(pipeline(ignorestatus(cmd), stdout=io, stderr=io))
    return process.exitcode
end

function execute(p; executor=run_command)
    v"1.13.0" <= VERSION < v"1.14.0" || error("Julia 1.13 required")
    Threads.nthreads() == 1 || error("Use one Julia thread in the driver")
    check_output(p.out) # Recheck after preparation, including for --dry-run.
    if p.o.dry
        foreach(c -> println(c.cmd), p.cases)
        println("DRY RUN: four frozen files; metadata checked; no subprocesses, outputs or SXM pixels read")
        return nothing
    end
    !isempty(get(ENV, "SLURM_JOB_ID", "")) || error("Run multifile diagnostics only in a Viper Slurm job")
    cpus = tryparse(Int, get(ENV, "SLURM_CPUS_PER_TASK", ""))
    cpus !== nothing && cpus >= 4 || error("At least four requested CPUs are required")
    mkdir(p.out) # Atomic creation also refuses a late collision/symlink.
    logs = joinpath(p.out, "logs"); mkdir(logs)
    stages = joinpath(p.out, "stages.tsv")
    write(stages, "case\tstage\texit_code\telapsed_s\n")
    status_lock = ReentrantLock()
    outcomes = asyncmap(p.cases; ntasks=4) do case
        started = time_ns(); code = 1
        open(joinpath(logs, "$(case.stem)_masked_preprocessing.log"), "w") do io
            println(io, case.cmd); flush(io)
            try
                code = Int(executor(case.cmd, io))
            catch err
                println(io, sprint(showerror, err)); code = 1
            end
        end
        elapsed = max(1, time_ns() - started) / 1e9
        row = (case.file, "masked_preprocessing", code, elapsed)
        # Append/close each completed stage. A later case failure or Slurm timeout
        # must not erase already completed records. No success-dependent filtering.
        lock(status_lock) do
            open(stages, "a") do io
                println(io, join(row, '\t'))
            end
        end
        return row
    end
    all(row -> row[3] == 0, outcomes) || error("Some masked-preprocessing cases failed; see preserved logs/stages.tsv")
    println("Four diagnostic commands completed; this is not scientific validation or a production change")
    return outcomes
end

function main(args=ARGS)
    p = prepare(args)
    return p === nothing ? println(HELP) : execute(p)
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
end
