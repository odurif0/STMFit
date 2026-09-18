#!/usr/bin/env julia
# One bounded four-file diagnostic job. Not a production pipeline or selector.
module LabelFreeExploration
using TOML

const HELP = """
Usage (Julia 1.13, Slurm for execution):
  julia --startup-file=no --threads=1 --project=. test/run_label_free_exploration.jl \
    --input-dir INPUT --config ORIGINAL_MOLECULE.toml \
    --settings config/label_free_exploration.toml --outdir NEW_DIR [--dry-run]
INPUT contains raw/, selected_summary.tsv, base_geometry.tsv, candidate_counts.tsv.
The candidate table has file,candidate_ns,source columns; counts are diagnostics,
not expected counts. Exactly four files run as four independent one-thread tasks.
Each task runs saved-geometry acquisition evidence then native/profile comparisons.
Dry-run checks paths/count metadata/settings and prints commands, reads no SXM pixels
and creates no output. Execution requires Slurm with at least four requested CPUs.
No original inputs, reference outputs or production methods are modified.
"""

function options(args)
    args == ["--help"] && return nothing
    allowed=Set(["--input-dir","--config","--settings","--outdir"])
    d=Dict{String,String}(); dry=false; i=1
    while i<=length(args)
        k=args[i]
        if k=="--dry-run"
            dry && error("Duplicate --dry-run")
            dry=true; i+=1; continue
        end
        k in allowed || error("Unknown option: $k")
        haskey(d,k) && error("Duplicate option: $k")
        i<length(args) && !startswith(args[i+1],"--") || error("Missing value: $k")
        isempty(args[i+1]) && error("Empty value: $k")
        d[k]=abspath(args[i+1]); i+=2
    end
    Set(keys(d))==allowed || error("All four path arguments are required; use --help")
    return (;d,dry)
end

function table(path, required)
    lines=readlines(path); isempty(lines) && error("Empty table: $path")
    h=String.(split(chomp(lines[1]),'\t';keepempty=true))
    length(h)==length(unique(h)) || error("Duplicate columns: $path")
    all(k->k in h,required) || error("Required columns missing: $path")
    rows=Dict{String,String}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        v=String.(split(chomp(line),'\t';keepempty=true))
        length(v)==length(h) || error("Malformed table: $path")
        push!(rows,Dict(zip(h,v)))
    end
    return rows
end

function prepare(args; project=dirname(@__DIR__))
    o=options(args); o===nothing && return nothing
    input=o.d["--input-dir"]; out=o.d["--outdir"]
    (ispath(out)||islink(out)) && error("Output already exists: $out")
    paths=(config=o.d["--config"],settings=o.d["--settings"],
        selected=joinpath(input,"selected_summary.tsv"),
        geometry=joinpath(input,"base_geometry.tsv"),
        candidates=joinpath(input,"candidate_counts.tsv"))
    all(isfile,values(paths)) || error("Missing required input file")
    model=TOML.parsefile(paths.config)
    all(k->haskey(model,k),("model","selection","preprocessing")) || error("Incomplete molecule config")
    settings=TOML.parsefile(paths.settings)
    all(k->haskey(settings,k),("counting_variable_projection","acquisition_noise")) || error("Missing diagnostic settings")
    selected=table(paths.selected,["filepath","N_selected","N_eff","runnerup_N_eff","refined_policy"])
    selected_keys=basename.(getindex.(selected,"filepath"))
    length(selected_keys)==length(unique(selected_keys)) || error("Duplicate selected-summary file")
    byfile=Dict(zip(selected_keys,selected))
    geometry=table(paths.geometry,["file","N","lobe"])
    candidates=table(paths.candidates,["file","candidate_ns","source"])
    length(candidates)==4 || error("This bounded job requires exactly four diagnostic files")
    names=getindex.(candidates,"file")
    length(unique(names))==4 || error("Duplicate diagnostic file")
    cases=NamedTuple[]
    for r in candidates
        file=r["file"]
        file==basename(file) && endswith(lowercase(file),".sxm") || error("Invalid raw filename")
        raw=joinpath(input,"raw",file); isfile(raw) || error("Missing raw file: $file")
        haskey(byfile,file) || error("No saved selection for $file")
        saved=byfile[file]; n=parse(Int,saved["N_selected"])
        expected_candidates=Set([n-1,n,n+1])
        for k in ("N_eff","runnerup_N_eff")
            v=tryparse(Int,saved[k]); v!==nothing && v>0 && push!(expected_candidates,v)
        end
        lo=Int(model["model"]["n_min"]); hi=Int(model["model"]["n_max"])
        filter!(n->lo<=n<=hi,expected_candidates)
        ns=parse.(Int,split(r["candidate_ns"],',';keepempty=true))
        length(unique(ns))==length(ns) && Set(ns)==expected_candidates ||
            error("Candidate list must equal saved selected-neighbors/effective/runner-up union: $file")
        gs=filter(g->basename(g["file"])==file,geometry)
        length(gs)==n && sort(parse.(Int,getindex.(gs,"lobe")))==collect(1:n) &&
            all(g->parse(Int,g["N"])==n,gs) || error("Saved geometry/count mismatch: $file")
        stem=splitext(file)[1]
        common=["--file",raw,"--config",paths.config,"--selected-summary",paths.selected]
        acq=vcat(["--geometry",paths.geometry,"--settings",paths.settings],common,
                 ["--outdir",joinpath(out,stem,"acquisition")])
        vp=vcat(["--diagnostic-config",paths.settings,"--candidate-ns",join(sort(ns),',')],common,
                ["--outdir",joinpath(out,stem,"variable_projection")])
        prefix=[joinpath(Sys.BINDIR,Base.julia_exename()),"--startup-file=no","--threads=1","--project=$(abspath(project))"]
        push!(cases,(file=file,stem=stem,
            acquisition=Cmd(vcat(prefix,[joinpath(project,"test","diagnose_acquisition_noise.jl")],acq)),
            variable_projection=Cmd(vcat(prefix,[joinpath(project,"test","diagnose_counting_variable_projection.jl")],vp))))
    end
    return (;o,paths,cases,out)
end

function execute(p)
    v"1.13.0" <= VERSION < v"1.14.0" || error("Julia 1.13 required")
    Threads.nthreads()==1 || error("Use one Julia thread in the driver")
    if p.o.dry
        for c in p.cases
            println(c.acquisition); println(c.variable_projection)
        end
        println("DRY RUN: four files; metadata consistent; no outputs or pixels read")
        return nothing
    end
    !isempty(get(ENV,"SLURM_JOB_ID","")) || error("Execute multifile diagnostics in a Viper Slurm job")
    parse(Int,get(ENV,"SLURM_CPUS_PER_TASK","0"))>=4 || error("Four CPUs are required")
    mkdir(p.out); mkpath(joinpath(p.out,"logs"))
    outcomes=asyncmap(p.cases;ntasks=4) do c
        mkdir(joinpath(p.out,c.stem))
        rows=Tuple{String,String,Int,Float64}[]
        for stage in (:acquisition,:variable_projection)
            cmd=getproperty(c,stage); started=time_ns(); code=1
            open(joinpath(p.out,"logs","$(c.stem)_$(stage).log"),"w") do io
                println(io,cmd); flush(io)
                try
                    proc=run(pipeline(ignorestatus(cmd),stdout=io,stderr=io))
                    code=proc.exitcode
                catch err
                    println(io,sprint(showerror,err)); code=1
                end
            end
            push!(rows,(c.file,string(stage),code,(time_ns()-started)/1e9))
        end
        rows
    end
    open(joinpath(p.out,"stages.tsv"),"w") do io
        println(io,"file\tstage\texit_code\telapsed_s")
        for rows in outcomes, row in rows; println(io,join(row,'\t')); end
    end
    all(row[3]==0 for rows in outcomes for row in rows) || error("Some diagnostic stages failed; see preserved logs/stages.tsv")
    println("Diagnostic stages completed; numerical/scientific checks remain separate")
end

function main(args=ARGS)
    p=prepare(args)
    p===nothing ? println(HELP) : execute(p)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
