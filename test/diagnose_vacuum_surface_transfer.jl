#!/usr/bin/env julia
# Raw source-view fitting and target-row holdout; no count or label input.
module DiagnoseVacuumSurfaceTransfer
using STMSXMIO, Statistics, LinearAlgebra, SHA, TOML
include(joinpath(@__DIR__,"lib/vacuum_surface_transfer.jl"))
include(joinpath(@__DIR__,"qe_vacuum_wavefunctions.jl"))
const F=VacuumSurfaceTransfer
const D=QEVacuumWavefunctions
const S=D.S
const V=D.V
sha(path)=S.sha(path)
num(r,k)=parse(Float64,r[k])
integer(r,k)=parse(Int,r[k])

function table(path,rows;header=isempty(rows) ? String[] : string.(keys(first(rows))))
    open(path,"w") do io
        println(io,join(header,'\t'))
        for r in rows
            println(io,join([replace(string(getproperty(r,Symbol(k))),'\n'=>' ','\t'=>' ') for k in header],'\t'))
        end
    end
end

function load_maps(root,c)
    m=c["model"]; run=joinpath(root,m["surface_run"])
    n=1+round(Int,2m["map_half_nm"]/m["map_step_nm"])
    grid=range(-m["map_half_nm"],m["map_half_nm"];length=n)
    heights=Dict{Int,Dict{Int,Matrix{Float64}}}(j=>Dict{Int,Matrix{Float64}}() for j in m["intervals"])
    isos=Dict{Int,Float64}(); refs=NamedTuple[]
    for (type,species) in enumerate(("glcn","glcnac"))
        dir=joinpath(run,species)
        for (kind,rel) in (("surface","continuous/surfaces.tsv"),("xml","data-file-schema.xml"),("frame","frame.tsv"))
            path=joinpath(dir,rel); hash=sha(path)
            hash==m[species*"_"*kind*"_sha256"] || error("Changed $species $kind")
            push!(refs,(;species,kind,sha256=hash))
        end
        rows=V.tsv(joinpath(dir,"continuous/surfaces.tsv"))
        length(rows)==n*n*length(m["intervals"]) || error("Incomplete surface family")
        frame=D.G.C.read_frame(joinpath(dir,"frame.tsv")); axis=frame.t_axis[1:2]; axis/=norm(axis)
        perp=[-axis[2],axis[1]]
        geo=D.xml_geometry(joinpath(dir,"data-file-schema.xml"))
        copper=[a.position_nm[3] for a in geo.geometry.atoms if a.z==29]
        isempty(copper) && error("Missing substrate reference")
        reference=minimum(copper)
        for j in m["intervals"]
            rr=sort(filter(r->integer(r,"interval")==j,rows);by=r->integer(r,"pixel"))
            integer.(rr,"pixel")==collect(1:n*n) || error("Incomplete lateral grid")
            all(r->r["valid"]=="true" && integer(r,"unresolved")==0 && integer(r,"root_count")==1,rr) || error("Uncertified root")
            surface=zeros(n,n)
            iso=num(first(rr),"isovalue")
            haskey(isos,j) && isos[j]!=iso && error("Class-specific isovalue")
            isos[j]=iso
            for (k,r) in enumerate(rr)
                it=1+mod(k-1,n); iu=1+fld(k-1,n)
                xy=frame.origin_nm[1:2]+grid[it]*axis+grid[iu]*perp
                isapprox(num(r,"x_nm"),xy[1];rtol=0,atol=c["preprocessing"]["plane_rank_rtol"]) &&
                    isapprox(num(r,"y_nm"),xy[2];rtol=0,atol=c["preprocessing"]["plane_rank_rtol"]) || error("Wrong surface coordinates")
                num(r,"isovalue")==iso || error("Inconsistent isovalue")
                surface[iu,it]=num(r,"z_nm")-reference
            end
            heights[j][type-1]=surface
        end
    end
    for maps in values(heights); maps[-1]=(maps[0]+maps[1])/2; end
    (;heights,isos,refs)
end

"Header only, without decoding pixels or any benchmark metadata."
function raw_header(path)
    lines=String[]
    open(path) do io
        while !eof(io)
            line=readline(io)
            strip(line)==":SCANIT_END:" && break
            push!(lines,line)
        end
    end
    STMSXMIO._parse_header(join(lines,'\n'))
end

bias(header)=parse(Float64,strip(get(header,"BIAS","NaN")))
matching_bias(header,c)=isapprox(bias(header),c["model"]["sample_bias_v"];rtol=0,atol=c["model"]["bias_atol_v"])

function source_data(img,direction,c)
    p=c["preprocessing"]
    channels=filter(ch->lowercase(ch.name)==lowercase(p["channel"]) && lowercase(ch.direction)==direction,img.channels)
    length(channels)==1 || error("Missing or ambiguous $direction Z channel")
    ch=only(channels)
    obs=STMSXMIO.preprocess_observed_channel(img,ch;stride=p["stride"],flatten=p["fit_flatten"],
        smooth_radius_px=0,plane_rank_rtol=p["plane_rank_rtol"])
    obs.unit=="nm" && obs.status=="observed_only" || error("Unusable physical height channel")
    detection=STMSXMIO.preprocess_observed_channel(img,ch;stride=p["stride"],flatten=p["anchor_flatten"],
        smooth_radius_px=p["anchor_smooth_radius_px"],plane_rank_rtol=p["plane_rank_rtol"])
    aa=F.anchors(obs.xs,obs.ys,obs.raw,detection.z_smooth,c)
    pp=F.patches(obs.xs,obs.ys,obs.z,aa,c)
    (;obs,aa,pp)
end

"Target values can change this scoring vector, never source fit or target calibration."
function targets(source,target,pair,c)
    result=Vector{Vector{Float64}}()
    for patch in source.pp
        vals=fill(NaN,length(patch.inds))
        if pair!==nothing
            lag=pair.reg.applied_dx_px; co=pair.coeff; obs=source.obs
            for (k,i) in enumerate(patch.inds)
                y,x=Tuple(i); xx=x+lag
                pair.masks.test[y] && 1<=xx<=size(target,2) && isfinite(target[y,xx]) || continue
                plane=obs.raw[i]-obs.z[i]
                vals[k]=target[y,xx]-plane-(co[1]+co[2]*obs.xs[x]+co[3]*obs.ys[y])
            end
        end
        push!(result,vals)
    end
    result
end

function one_direction(img,direction,maps,c,out)
    S.newdir(out); source=source_data(img,direction,c)
    target_direction=direction=="fwd" ? "bwd" : "fwd"
    target_channels=filter(ch->lowercase(ch.name)=="z" && ch.direction==target_direction,img.channels)
    table(joinpath(out,"anchors.tsv"),[(;anchor=j,a...,used=any(p->p.anchor==j,source.pp)) for (j,a) in enumerate(source.aa)])
    isempty(source.pp) && return (;direction,status="no_supported_anchor",anchors=length(source.aa),patches=0)
    pair=nothing; reason="ok"; target=fill(NaN,size(source.obs.raw))
    try
        length(target_channels)==1 || error("Missing target direction")
        ch=only(target_channels); scale,unit=STMSXMIO._value_scale(ch.unit)
        unit=="nm" || error("Wrong target units")
        target=ch.data.*scale
        pair=F.calibrate_pair(source.obs.raw,target,source.obs.xs,source.obs.ys,c)
    catch e
        e isa InterruptException && rethrow()
        reason=replace(sprint(showerror,e),'\n'=>' ','\t'=>' ')
    end
    heldout=targets(source,target,pair,c)
    if pair!==nothing
        table(joinpath(out,"registration.tsv"),[(;region=k-1,p...) for (k,p) in enumerate(pair.reg.peaks)])
        table(joinpath(out,"registration_scores.tsv"),[(;region=k-1,lag=lag,correlation=pair.reg.correlations[k,j])
            for k in axes(pair.reg.correlations,1) for (j,lag) in enumerate(pair.reg.lags)])
    end
    open(joinpath(out,"calibration.toml"),"w") do io
        TOML.print(io,Dict("status"=>reason,"identified"=>pair!==nothing && pair.reg.accepted,
            "registration_status"=>pair===nothing ? "unavailable" : pair.reg.status,
            "applied_dx_px"=>pair===nothing ? 0 : pair.reg.applied_dx_px,
            "difference_plane"=>pair===nothing ? fill(NaN,3) : pair.coeff,
            "calibration_pixels"=>pair===nothing ? 0 : pair.calibration_pixels))
    end
    table(joinpath(out,"observations.tsv"),[(;patch=j,anchor=p.anchor,pixel=k,row=i[1],column=i[2],
        dx_nm=p.dx[k],dy_nm=p.dy[k],source_nm=p.values[k],target_nm=heldout[j][k])
        for (j,p) in enumerate(source.pp) for (k,i) in enumerate(p.inds)])
    models=NamedTuple[]; patchrows=NamedTuple[]
    flat=mean(vcat(getproperty.(source.pp,:values)...))
    for interval in c["model"]["intervals"], arm in c["model"]["arms"]
        r=F.fit_dictionary(source.pp,maps.heights[interval],c,arm)
        table(joinpath(out,"costs_$(interval)_$(arm).tsv"),[(;patch=j,state=k,s...,pixels=cost.n,
            residual_mean_nm=cost.mu[k],centered_sse_nm2=cost.v[k],chosen=k==r.fit.choices[j])
            for (j,cost) in enumerate(r.costs) for (k,s) in enumerate(r.ss)])
        train=0.; test=0.; flaterr=0.; copyerr=0.; npix=0
        for (j,p) in enumerate(source.pp)
            chosen=r.ss[r.fit.choices[j]]; pred=r.pred[j]; good=isfinite.(heldout[j])
            trainerr=sum(abs2,p.values-pred); n=count(good)
            err=sum(abs2,heldout[j][good]-pred[good]); flatloss=sum(abs2,heldout[j][good].-flat)
            copyloss=sum(abs2,heldout[j][good]-p.values[good])
            train+=trainerr; test+=err; flaterr+=flatloss; copyerr+=copyloss; npix+=n
            push!(patchrows,(;interval,arm,patch=j,anchor=p.anchor,chosen...,offset_nm=r.fit.offset,
                training_pixels=length(p.values),training_sse_nm2=trainerr,test_pixels=n,
                test_sse_nm2=err,flat_sse_nm2=flatloss,source_copy_sse_nm2=copyloss))
        end
        push!(models,(;interval,isovalue=maps.isos[interval],arm,offset_nm=r.fit.offset,flat_height_nm=flat,
            training_pixels=sum(length(p.values) for p in source.pp),training_sse_nm2=train,
            optimizer_sse_nm2=r.fit.loss,offset_segments=r.fit.segments,test_pixels=npix,
            test_sse_nm2=test,flat_sse_nm2=flaterr,source_copy_sse_nm2=copyerr))
    end
    table(joinpath(out,"models.tsv"),models); table(joinpath(out,"patches.tsv"),patchrows)
    (;direction,status=reason,anchors=length(source.aa),patches=length(source.pp))
end

function process_file(file,maps,c,out)
    S.newdir(out); hash=sha(file); header=raw_header(file)
    open(joinpath(out,"input.toml"),"w") do io
        TOML.print(io,Dict("file"=>basename(file),"raw_sha256"=>hash,"bias_v"=>bias(header)))
    end
    if !matching_bias(header,c)
        table(joinpath(out,"status.tsv"),[(;direction="both",status="unsupported_bias",anchors=0,patches=0)])
        return
    end
    img=read_sxm(file); statuses=NamedTuple[]
    for direction in ("fwd","bwd")
        try
            push!(statuses,one_direction(img,direction,maps,c,joinpath(out,direction)))
        catch e
            e isa InterruptException && rethrow()
            push!(statuses,(;direction,status="failed:"*replace(sprint(showerror,e),'\n'=>' ','\t'=>' '),anchors=0,patches=0))
        end
    end
    sha(file)==hash || error("Raw data changed")
    table(joinpath(out,"status.tsv"),statuses)
end

function main(args=ARGS)
    "--help" in args && return println("diagnose_vacuum_surface_transfer.jl --data-dir DIR --config TOML --root REPO --outdir NEW_DIR [--file ONE_SXM] [--dry-run]")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Duplicate option")
        if k=="--dry-run"; opts[k]="true"; i+=1; continue; end
        k in ("--data-dir","--config","--root","--outdir","--file") && i<length(args) || error("Unexpected/missing option")
        opts[k]=args[i+1]; i+=2
    end
    all(haskey(opts,k) for k in ("--data-dir","--config","--root","--outdir")) || error("Missing input")
    c=F.settings(opts["--config"]); maps=load_maps(opts["--root"],c)
    names=sort(filter(f->endswith(lowercase(f),".sxm") && isfile(joinpath(opts["--data-dir"],f)),readdir(opts["--data-dir"])))
    if haskey(opts,"--file")
        opts["--file"] in names || error("Unknown raw file")
        names=[opts["--file"]]
    end
    isempty(names) && error("No raw images")
    if haskey(opts,"--dry-run")
        eligible=count(f->matching_bias(raw_header(joinpath(opts["--data-dir"],f)),c),names)
        println("Metadata-only: $(length(names)) raw scans, $eligible at the model bias, three fixed isovalues, both directions. No fit/output.")
        return
    end
    out=opts["--outdir"]; S.newdir(out); cp(opts["--config"],joinpath(out,"settings.toml"))
    table(joinpath(out,"surface_hashes.tsv"),maps.refs)
    table(joinpath(out,"files.tsv"),[(;file=f) for f in names])
    BLAS.set_num_threads(1)
    Threads.@threads :static for k in eachindex(names)
        f=names[k]; started=time_ns()
        process_file(joinpath(opts["--data-dir"],f),maps,c,joinpath(out,splitext(f)[1]))
        println("TRANSFER_FILE_COMPLETE ",f," elapsed_s=",(time_ns()-started)/1e9)
        flush(stdout)
    end
    println("SURFACE_TRANSFER_COMPLETE; no external labels, production change or grade")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && DiagnoseVacuumSurfaceTransfer.main()
