#!/usr/bin/env julia
# Image-only estimation: no molecular count, geometry, template or label is read.
module EstimateImageRegistration2D
include(joinpath(@__DIR__,"estimate_acquisition_shifts.jl"))
include(joinpath(@__DIR__,"lib","image_registration_2d.jl"))
using .EstimateAcquisitionShifts: read_sxm, get_channel, require_direction, PatternConfig, preprocess_channel,
    load_patch_preprocessing, read_table, write_table
using .ImageRegistration2D, Statistics, SHA, TOML
const R=ImageRegistration2D
const OPTIONS=Set(["--raw-file","--base-shifts","--config","--settings","--outdir"])
table(path)=last(read_table(path))
record(x)=Dict(string(k)=>string(v) for (k,v) in pairs(x))
write_records(path,rows)=write_table(path,sort(collect(keys(first(rows)))),rows)

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option")
        if k=="--dry-run"; opts[k]="true"; i+=1; continue; end
        k in union(OPTIONS,Set(["--fold"])) && i<length(args) || error("Missing/forbidden image registration option")
        opts[k]=args[i+1]; i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All image inputs required")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    for k in setdiff(OPTIONS,Set(["--outdir"]))
        isfile(opts[k]) || error("Missing input $k")
    end
    (ispath(opts["--outdir"]) || islink(opts["--outdir"])) && error("Output exists")
    parse(Int,get(opts,"--fold","0")) in 0:2 || error("Invalid fold")
    return opts
end

function inputs(opts)
    s=R.settings(opts["--settings"]); pre=load_patch_preprocessing(opts["--config"])
    basename(opts["--raw-file"])==s.file || error("Wrong scan")
    header,rows=read_table(opts["--base-shifts"])
    header==["file","bwd_sample_dx_px"] || error("Expected image-only integer shift metadata")
    length(unique(getindex.(rows,"file")))==length(rows) || error("Repeated base shift")
    dx=parse(Int,only(r for r in rows if r["file"]==s.file)["bwd_sample_dx_px"])
    return (;s,pre,base=(dx,0))
end

function load_images(path,pre)
    img=read_sxm(path); cf=get_channel(img,"Z";direction="fwd"); cb=get_channel(img,"Z";direction="bwd")
    require_direction(cf,"fwd"); require_direction(cb,"bwd")
    cfg=PatternConfig(filepath=path,channel="Z",direction="fwd",stride=pre.stride,
        flatten=pre.flatten,smooth_radius_px=pre.smooth_radius_px,no_plot=true)
    xs,ys,_,a,_,_,_=preprocess_channel(img,cf,cfg)
    xb,yb,_,b,_,_,_=preprocess_channel(img,cb,cfg)
    xs==xb && ys==yb || error("View grids differ")
    a[.!isfinite.(cf.data[1:pre.stride:end,1:pre.stride:end])].=NaN
    b[.!isfinite.(cb.data[1:pre.stride:end,1:pre.stride:end])].=NaN
    step=(median(diff(xs)),median(diff(ys)))
    all(isapprox.(diff(xs),step[1];rtol=1e-10,atol=1e-12)) &&
        all(isapprox.(diff(ys),step[2];rtol=1e-10,atol=1e-12)) || error("Nonuniform grid")
    return (;a,b,xs,ys,step)
end

function save_result(out,r,im,input,elapsed)
    mkpath(out); s=input.s; base=input.base
    write_records(joinpath(out,"summary.tsv"),[record((file=s.file,fold=r.fold,accepted=r.accepted,status=r.status,
        dx=first(r.peaks).dx,dy=first(r.peaks).dy,base_dx=base[1],base_dy=base[2],
        step_x_nm=im.step[1],step_y_nm=im.step[2],limit_x_px=r.limits[1],limit_y_px=r.limits[2],
        nrows=size(im.a,1),ncols=size(im.a,2),support_pixels=count(r.mask),elapsed_s=elapsed))])
    write_records(joinpath(out,"peaks.tsv"),[merge(record(p),Dict("region"=>string(i-1),"fold"=>string(r.fold))) for (i,p) in enumerate(r.peaks)])
    scores=Dict{String,String}[]
    for region in 0:s.bands,j in eachindex(r.coarse)
        push!(scores,record((;phase="coarse",region,dx=r.coarse[j][1],dy=r.coarse[j][2],correlation=r.scores[region+1,j])))
    end
    for f in r.fine,j in eachindex(f.shifts)
        push!(scores,record((phase="fine",region=f.region,dx=f.shifts[j][1],dy=f.shifts[j][2],correlation=f.values[j])))
    end
    write_records(joinpath(out,"scores.tsv"),scores)
    open(joinpath(out,"support.tsv"),"w") do io
        println(io,"row\tcolumn")
        for x in axes(r.mask,2),y in axes(r.mask,1)
            r.mask[y,x] && println(io,y,'\t',x)
        end
    end
    if r.fold==0
        open(joinpath(out,"images.tsv"),"w") do io
            println(io,"row\tcolumn\tx_nm\ty_nm\tfwd_nm\tbwd_raw_nm")
            for x in axes(im.a,2),y in axes(im.a,1)
                println(io,join((y,x,im.xs[x],im.ys[y],im.a[y,x],im.b[y,x]),'\t'))
            end
        end
    end
end

function execute(opts)
    input=inputs(opts); fold=parse(Int,get(opts,"--fold","0"))
    haskey(opts,"--dry-run") && return println("Image-only dry run: $(input.s.file), base $(input.base), residual ±$(input.s.residual_window_nm) nm, fold $fold; no image scan/output.")
    im=load_images(opts["--raw-file"],input.pre)
    started=time_ns(); r=R.scan(im.a,im.b,im.step,input.base,input.s;fold)
    save_result(opts["--outdir"],r,im,input,(time_ns()-started)/1e9)
    hashes=Dict{String,String}[]
    for key in ("--raw-file","--base-shifts","--config","--settings")
        push!(hashes,Dict("input"=>key,"sha256"=>bytes2hex(sha256(read(opts[key])))))
    end
    write_records(joinpath(opts["--outdir"],"input_hashes.tsv"),hashes)
    cp(opts["--config"],joinpath(opts["--outdir"],"physical.toml")); cp(opts["--settings"],joinpath(opts["--outdir"],"settings.toml"))
    println("fold=$fold, image shift=($(first(r.peaks).dx),$(first(r.peaks).dy)) px; $(r.status)")
end

function merge_chunks(out)
    dirs=[joinpath(out,"fold$f") for f in 0:2]
    all(isdir,dirs) || error("Missing image fold")
    for name in ("input_hashes.tsv","physical.toml","settings.toml")
        all(read(joinpath(dir,name))==read(joinpath(first(dirs),name)) for dir in dirs) || error("Fold inputs differ")
    end
    rows=[only(table(joinpath(d,"summary.tsv"))) for d in dirs]
    parse.(Int,getindex.(rows,"fold"))==[0,1,2] || error("Wrong fold outputs")
    peaks=vcat([table(joinpath(d,"peaks.tsv")) for d in dirs]...)
    s=R.settings(joinpath(first(dirs),"settings.toml"))
    runs=[(fold=parse(Int,r["fold"]),accepted=r["accepted"]=="true",peaks=[(dx=parse(Float64,r["dx"]),dy=parse(Float64,r["dy"]))]) for r in rows]
    eligible=R.eligibility(runs,s)
    for name in ("input_hashes.tsv","physical.toml","settings.toml")
        cp(joinpath(first(dirs),name),joinpath(out,name))
    end
    write_records(joinpath(out,"summary.tsv"),rows); write_records(joinpath(out,"peaks.tsv"),peaks)
    write_records(joinpath(out,"eligibility.tsv"),[record(eligible)])
    println("Image-only registration: $eligible. No molecular fit or recognition grade.")
end

function main(args=ARGS)
    length(args)==2 && first(args)=="--merge" && return merge_chunks(last(args))
    opts=parse_cli(args)
    opts===nothing && return println("estimate_image_registration_2d.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--fold 0|1|2] [--dry-run]; or --merge RUN_DIR")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && EstimateImageRegistration2D.main()
