#!/usr/bin/env julia
# QE token-order correction only. No STM data, benchmark or production changes.
module BuildQEOrderMolds
include(joinpath(@__DIR__,"lib","qe_cube_molds.jl"))
using .QECubeMolds, TOML, SHA, Printf
const OPTIONS=Set(["--cube0","--cube1","--frame0","--frame1","--config","--settings","--outdir"])

function options(args)
    if args==["--help"]
        println("build_qe_order_molds.jl --cube0 CUBE --cube1 CUBE --frame0 TSV --frame1 TSV --config ASSIGNMENT_TOML --settings ORDER_TOML --outdir NEW_DIR [--dry-run]")
        return nothing
    end
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option: $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in OPTIONS && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden option: $k")
        o[k]=args[i+1]; i+=2
    end
    all(haskey(o,k) for k in OPTIONS) || error("All inputs required")
    for k in setdiff(OPTIONS,Set(["--outdir"])); isfile(o[k]) || error("Missing $k"); end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output already exists")
    QECubeMolds.settings(o["--settings"])
    QECubeMolds.CCMoldNative.load_settings(o["--config"])
    o
end

function execute(o)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    haskey(o,"--dry-run") && return println("QE-order templates: paths/settings checked, no cube parsed or output")
    out=abspath(o["--outdir"]); mkpath(out)
    try
        r=QECubeMolds.CCMoldNative.build_molds(cube0=o["--cube0"],cube1=o["--cube1"],
            frame0=o["--frame0"],frame1=o["--frame1"],config=o["--config"],
            out=joinpath(out,"templates_cc.tsv"),cube_reader=QECubeMolds.read_qe_cube)
        cp(o["--config"],joinpath(out,"assignment_config.toml")); cp(o["--settings"],joinpath(out,"order_settings.toml"))
        provenance=Dict("cube_order"=>"qe_last_axis_fast","julia_version"=>string(VERSION),
            "inputs"=>Dict(k[3:end]=>Dict("path"=>abspath(o[k]),"sha256"=>bytes2hex(sha256(read(o[k]))))
                for k in setdiff(OPTIONS,Set(["--outdir"]))),
            "template_sha256"=>bytes2hex(sha256(read(joinpath(out,"templates_cc.tsv")))))
        open(io->TOML.print(io,provenance),joinpath(out,"provenance.toml"),"w")
        open(joinpath(out,"surface_audit.tsv"),"w") do io
            println(io,"type\tpixel\tisovalue\tmean_height\tvalid_columns\theight\tnormalized")
            for (typ,result) in enumerate((r.type0,r.type1)), i in eachindex(result.heights)
                @printf(io,"%d\t%d\t%.17g\t%.17g\t%d\t%.17g\t%.17g\n",typ-1,i,
                    result.iso,result.mean_height,result.nvalid,result.heights[i],result.normalized[i])
            end
        end
        for (typ,result) in enumerate((r.type0,r.type1))
            println("type",typ-1,": iso=",result.iso," mean_height=",result.mean_height," valid=",result.nvalid,"/",length(result.heights))
        end
        r
    catch err
        open(io->println(io,sprint(showerror,err)),joinpath(out,"failure.txt"),"w")
        rethrow()
    end
end
function main(args=ARGS)
    o=options(args); o===nothing || execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && BuildQEOrderMolds.main()
