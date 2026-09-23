#!/usr/bin/env julia
# The only changed signals are two physical CC scores. No labels or refit.
module TangentMoldScorer
module Connected
include(joinpath(@__DIR__,"score_connected_mold_templates.jl"))
end
include(joinpath(@__DIR__,"lib","tangent_mold_projection.jl"))
include(joinpath(@__DIR__,"lib","reconstructed_unit_assignment.jl"))
include(joinpath(@__DIR__,"lib","patch_preprocessing.jl"))
using .ReconstructedUnitAssignment: lobe_table, require_same_keys, write_table
using STMSXMIO: read_sxm, _coordinate_vectors
using LinearAlgebra, TOML, Printf
const T = TangentMoldProjection
const OPTIONS=Set(["--patches","--templates","--features","--data-dir","--count-config","--config","--settings","--prefix","--out"])

function parse_cli(args)
    if "--help" in args
        println("score_tangent_mold_templates.jl: ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--chunk I/N] [--dry-run]")
        return nothing
    end
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in union(OPTIONS,Set(["--chunk"])) && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden tangent scorer option")
        o[k]=args[i+1]; i+=2
    end
    all(haskey(o,k) for k in OPTIONS) || error("All tangent scorer inputs required")
    for k in setdiff(OPTIONS,Set(["--data-dir","--prefix","--out"]))
        isfile(o[k]) || error("Missing $k")
    end
    isdir(o["--data-dir"]) || error("Missing raw directory")
    o["--prefix"] in ("res","bwd_res") || error("Only the two native residual families are supported")
    for suffix in ("", ".audit.tsv", ".basis.tsv", ".failures.tsv")
        p=o["--out"]*suffix
        (ispath(p)||islink(p)) && error("Output exists: $p")
    end
    chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    length(chunk)==2 && 1<=chunk[1]<=chunk[2] || error("Invalid shard")
    T.settings(o["--settings"])
    return o
end

function execute(o)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    BLAS.set_num_threads(1)
    opt=T.settings(o["--settings"])
    cfg=ReconstructedUnitAssignment.load_config(o["--config"])
    pre=PatchPreprocessing.load_patch_preprocessing(o["--count-config"])
    cfg["preprocessing"]["patch_residual_filter"]=="smooth_residual" || error("Matched residual required")
    cfg["preprocessing"]["pixel_order"]=="u_outer_t_inner" || error("Native patch order required")
    header,base=lobe_table(o["--features"])
    _,patchrows=lobe_table(o["--patches"])
    require_same_keys(base,patchrows,"tangent patch keys")
    for r in values(base); T.validate_geometry(r); end
    templates,pix=Connected._load_templates(o["--templates"])
    Connected._contrast_unary_templates!(templates)
    byfile,patchpix=Connected._load_patches(o["--patches"],o["--prefix"]*"_p")
    coords=collect(-cfg["model"]["mold_half_nm"]:cfg["model"]["mold_step_nm"]:cfg["model"]["mold_half_nm"])
    expected=[@sprintf("p%03d",i) for i in 1:length(coords)^2]
    pix==expected && patchpix==o["--prefix"].*"_".*expected || error("Pixel order/grid mismatch")
    for f in keys(byfile)
        isfile(joinpath(o["--data-dir"],f)) || error("Missing raw $f")
        rs=byfile[f]; n=length(rs)
        [r.lobe for r in rs]==collect(1:n) || error("Noncontiguous lobe keys")
        all(parse(Int,base[(f,r.lobe)]["N"])==n for r in rs) || error("N/key mismatch")
    end
    chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    files=[f for (i,f) in enumerate(sort(collect(keys(byfile)))) if mod1(i,chunk[2])==chunk[1]]
    isempty(files) && error("Empty tangent shard")
    haskey(o,"--dry-run") && return println("Tangent mold scoring: $(length(files)) files; metadata only, no raw pixels or output")
    index=IdDict(v=>k for (k,v) in templates)
    legacy=Connected.Options(o["--patches"],o["--templates"],nothing,o["--out"],o["--prefix"],"ncc","contrast",0.0,1.0)
    audit=Dict{String,String}[]; failures=Dict{String,String}[]
    output=IOBuffer(); println(output,join(Connected.SCORE_HEADER,'\t'))
    mkpath(dirname(o["--out"]))
    # Save the actual pixel-space design for a separate QR/arithmetic check.
    open(o["--out"]*".basis.tsv","w") do bio
        maxcols=opt.adjacent_amplitudes ? 10 : opt.orientation ? 9 : 8
        println(bio,"file\tlobe\tpixel\tbasis_count\t",join(["b$j" for j in 1:maxcols],'\t'))
        for (fi,f) in enumerate(files)
            try
                img=read_sxm(joinpath(o["--data-dir"],f))
                xs,ys=_coordinate_vectors(img;stride=pre.stride)
                scores=Dict{Tuple{Int,Int,Int},Vector{Float64}}()
                for r in byfile[f]
                    neighbors=opt.adjacent_amplitudes ? T.adjacent_rows(base,f,r.lobe,length(byfile[f])) : AbstractDict[]
                    B=T.native_design(xs,ys,base[(f,r.lobe)],coords,pre.smooth_radius_px;orientation=opt.orientation,neighbors)
                    for k in axes(B,1)
                        vals=vcat(B[k,:],fill(NaN,maxcols-size(B,2)))
                        println(bio,join(vcat([f,string(r.lobe),string(k),string(size(B,2))],[@sprintf("%.17g",v) for v in vals]),'\t'))
                    end
                    for parity in 0:1, mirror in 0:1
                        M=hcat(templates[(0,parity,mirror)],templates[(1,parity,mirror)])
                        result=T.projected_score(r.patch,M,B,opt)
                        scores[(r.lobe,parity,mirror)]=result.costs
                        vals=merge((file=f,lobe=r.lobe,parity=parity,mirror=mirror,
                            cost0=result.costs[1],cost1=result.costs[2]),Base.structdiff(result,NamedTuple{(:costs,)}((result.costs,))))
                        push!(audit,Dict(string(k)=>string(v) for (k,v) in pairs(vals)))
                    end
                end
                scorer=(r,t)->begin
                    typ,parity,mirror=index[t]
                    scores[(r.lobe,parity,mirror)][typ+1]
                end
                best=Connected._decode_file(byfile[f],templates,nothing,legacy;scorer,omit_unavailable=opt.omit_unavailable)
                Connected.write_decoded(output,f,byfile[f],best)
                println("[$fi/$(length(files))] $f: $(length(byfile[f])) fixed lobes; tangent CC only"); flush(stdout)
            catch err
                push!(failures,Dict("file"=>f,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')))
                println(stderr,"FAILED $f: ",sprint(showerror,err))
            end
        end
    end
    !isempty(audit) && write_table(o["--out"]*".audit.tsv",sort(collect(keys(first(audit)))),audit)
    if !isempty(failures)
        write_table(o["--out"]*".failures.tsv",["file","reason"],failures)
        error("$(length(failures)) tangent failures; no partial score table")
    end
    write(o["--out"],take!(output))
end
function main(args=ARGS)
    o=parse_cli(args); o===nothing || execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && TangentMoldScorer.main()
