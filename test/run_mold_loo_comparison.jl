#!/usr/bin/env julia
# Two saved-cost arms, identical upstream inputs and classifier, no labels.
module MoldLOOComparison
include(joinpath(@__DIR__,"run_mold_state_comparison.jl"))
include(joinpath(@__DIR__,"lib","mold_leave_one_out.jl"))
const D=MoldStateComparison
const L=MoldLeaveOneOut
using .MoldStateComparison.Pipeline.ReconstructedUnitAssignment, SHA, LinearAlgebra

function write_scores(dir,data)
    diagnostics=Dict{String,String}[]
    for (j,view) in enumerate(("fwd","bwd"))
        path=joinpath(dir,"score_"*view*".tsv"); ispath(path) && error("Score output exists")
        open(path,"w") do io
            println(io,join(L.SCORE_HEADER,'\t'))
            for f in data.files
                records=D.M.records_for(data.base,f)
                decoded=L.decode(records,data.views[j],data.opt)
                L.write_decoded(io,records,decoded)
                for (r,b) in zip(records,decoded), p in 0:1, m in 0:1
                    push!(diagnostics,Dict("file"=>f,"lobe"=>string(r.lobe),"view"=>view,
                        "phase"=>string(p),"mirror"=>string(m),"training_cost"=>string(b.objectives[1+2p+m]),
                        "training_lobes"=>string(b.training_lobes),"selected"=>string((p,m)==(b.phase,b.mirror)),
                        "reason"=>b.reason))
                end
            end
        end
    end
    write_table(joinpath(dir,"state_selection.tsv"),["file","lobe","view","phase","mirror",
        "training_cost","training_lobes","selected","reason"],diagnostics)
end

function execute(o;classifier=D.classify)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    BLAS.set_num_threads(1)
    data=D.load_inputs(o;read_settings=L.settings)
    haskey(o,"--dry-run") && return println("Two cached decoder arms: $(length(data.files)) scans / $(length(data.base)) keys; no inference or output")
    out=abspath(o["--outdir"]); (ispath(out) || islink(out)) && error("Output already exists"); mkpath(out)
    inputs=vcat([joinpath(data.ref,f) for f in (D.INPUTS...,D.REPLAY_FILES...)],data.paths...,[o["--config"],o["--settings"]])
    write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
        [Dict("input"=>abspath(p),"sha256"=>bytes2hex(sha256(read(p)))) for p in inputs])
    cp(o["--config"],joinpath(out,"assignment_config.toml")); cp(o["--settings"],joinpath(out,"decoder_settings.toml"))
    stage="initialization"
    try
        for arm in ("reference","loo")
            stage=arm; dir=joinpath(out,arm); mkpath(joinpath(dir,"logs"))
            for f in D.INPUTS; cp(joinpath(data.ref,f),joinpath(dir,f)); end
            arm=="reference" ? D.write_scores(dir,data,:finite) : write_scores(dir,data)
            if arm=="reference"
                for f in ("score_fwd.tsv","score_bwd.tsv")
                    read(joinpath(dir,f))==read(joinpath(data.ref,f)) || error("Finite control differs: $f")
                end
            end
            classifier(dir,data.cfg,abspath(o["--config"]))
            _,pred=lobe_table(joinpath(dir,"predictions.tsv")); require_same_keys(data.base,pred,"complete cohort")
            if arm=="reference"
                for f in D.REPLAY_FILES
                    read(joinpath(dir,f))==read(joinpath(data.ref,f)) || error("Reference replay differs: $f")
                end
            end
        end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Two complete decoder arms; external grading remains separate")
end
function main(args=ARGS)
    if "--help" in args
        println("run_mold_loo_comparison.jl --reference-dir FINITE_TANGENT_DIR --config TOML --settings TOML --outdir NEW_DIR [--dry-run]")
        return
    end
    execute(D.parse_cli(args))
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && MoldLOOComparison.main()
