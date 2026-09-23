#!/usr/bin/env julia
# No N, molecular geometry, assignment output, template or label is an input.
module DiagnoseDirectionalResponse
include(joinpath(@__DIR__,"estimate_image_registration_2d.jl"))
include(joinpath(@__DIR__,"lib","directional_response.jl"))
using .DirectionalResponse, SHA, TOML, LinearAlgebra
const R=DirectionalResponse
const I=EstimateImageRegistration2D
const OPTIONS=Set(["--data-dir","--base-shifts","--config","--settings","--outdir"])
record(x)=Dict(string(k)=>string(v) for (k,v) in pairs(x))
table(path)=last(I.read_table(path))
write_records(path,rows)=I.write_table(path,sort(collect(keys(first(rows)))),rows)

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option")
        if k=="--dry-run"; opts[k]="true"; i+=1; continue; end
        k in OPTIONS && i<length(args) || error("Missing/forbidden response option")
        opts[k]=args[i+1]; i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All response inputs required")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    for k in ("--base-shifts","--config","--settings")
        isfile(opts[k]) || error("Missing $k")
    end
    isdir(opts["--data-dir"]) || error("Missing image directory")
    (ispath(opts["--outdir"]) || islink(opts["--outdir"])) && error("Refusing response output overwrite")
    return opts
end

function inputs(opts)
    s=R.settings(opts["--settings"]); pre=I.load_patch_preprocessing(opts["--config"])
    pre.stride==1 || error("Pilot requires native pixels")
    head,rows=I.read_table(opts["--base-shifts"])
    head==["file","bwd_sample_dx_px"] || error("Only image shift metadata allowed")
    panel=R.panel(getindex.(rows,"file"),s)
    shifts=Dict(r["file"]=>parse(Int,r["bwd_sample_dx_px"]) for r in rows)
    all(isfile(joinpath(opts["--data-dir"],p.file)) for p in panel) || error("Incomplete acquisition panel")
    return (;s,pre,panel,shifts)
end

function aligned_backward(b,dx)
    ny,nx=size(b); abs(dx)<nx || error("Invalid base displacement")
    out=fill(NaN,ny,nx)
    for x in 1:nx
        source=x+dx
        1<=source<=nx && (out[:,x].=b[:,source])
    end
    return out
end

function save_images(path,im,b)
    open(path,"w") do io
        println(io,"row\tcolumn\tx_nm\ty_nm\tfwd_nm\tbwd_raw_nm\tbwd_registered_nm")
        for x in axes(im.a,2),y in axes(im.a,1)
            println(io,join((y,x,im.xs[x],im.ys[y],im.a[y,x],im.b[y,x],b[y,x]),'\t'))
        end
    end
end

function one_scan(opts,input,p;loader=I.load_images,header_reader=path->I.read_sxm(path).header)
    out=joinpath(opts["--outdir"],splitext(p.file)[1]); mkpath(out)
    path=joinpath(opts["--data-dir"],p.file); base=input.shifts[p.file]
    write_records(joinpath(out,"input_hashes.tsv"),[Dict("input"=>p.file,"sha256"=>bytes2hex(sha256(read(path))))])
    started=time_ns()
    try
        header=header_reader(path)
        times=parse.(Float64,split(header["SCAN_TIME"]))
        length(times)==2 && all(x->isfinite(x) && x>0,times) && times[1]==times[2] || error("Mirrored response needs equal forward/backward line times")
        im=loader(path,input.pre); b=aligned_backward(im.b,base)
        write_records(joinpath(out,"metadata.tsv"),[merge(record(p),record((;base_dx=base,nrows=size(im.a,1),ncols=size(im.a,2),
            step_x_nm=im.step[1],step_y_nm=im.step[2],line_time_fwd_s=times[1],line_time_bwd_s=times[2],
            scan_angle_deg=parse(Float64,header["SCAN_ANGLE"]),scan_direction=header["SCAN_DIR"],bias_v=parse(Float64,header["BIAS"]))))])
        save_images(joinpath(out,"images.tsv"),im,b)
        r=R.scan(im.a,b,input.s;progress=(lag,sign)->begin
            println(p.file," lag=",lag," sign=",sign," scored"); flush(stdout)
        end)
        write_records(joinpath(out,"scores.tsv"),record.(r.records))
        write_records(joinpath(out,"selected.tsv"),record.(r.chosen))
        write_records(joinpath(out,"heldout.tsv"),record.(r.checks))
        write_records(joinpath(out,"summary.tsv"),[merge(record(p),record(r.summary),record((;elapsed_s=(time_ns()-started)/1e9,
            support_pixels=count(r.support.mask),radius_x=r.support.rx,radius_y=r.support.ry,status="completed")))])
        open(joinpath(out,"support.tsv"),"w") do io
            println(io,"row\tcolumn\tfold")
            for x in axes(r.support.mask,2),y in axes(r.support.mask,1)
                r.support.mask[y,x] && println(io,join((y,x,r.support.fold[y]),'\t'))
            end
        end
        println(p.file," complete: ",r.summary); flush(stdout)
    catch err
        reason=replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')
        write_records(joinpath(out,"failure.tsv"),[Dict("file"=>p.file,"reason"=>reason)])
        println(stderr,p.file," FAILED: ",reason); flush(stderr)
    end
end

function summarize(out,panel,s)
    rows=Dict{String,String}[]
    for p in panel
        dir=joinpath(out,splitext(p.file)[1])
        if isfile(joinpath(dir,"failure.tsv"))
            push!(rows,Dict("file"=>p.file,"status"=>"failed","supported"=>"false","reason"=>only(table(joinpath(dir,"failure.tsv")))["reason"]))
        else
            a=only(table(joinpath(dir,"summary.tsv")))
            push!(rows,Dict("file"=>p.file,"status"=>a["status"],"supported"=>a["supported"],"reason"=>"see per-scan heldout and summary"))
        end
    end
    write_records(joinpath(out,"summary.tsv"),rows)
    supporting=count(r->r["supported"]=="true",rows)
    complete=all(r->r["status"]=="completed",rows)
    write_records(joinpath(out,"panel_result.tsv"),[record((;complete,supporting,total=length(rows),
        directional_hypothesis_supported=complete && supporting>=s.min_supporting_scans,
        recognition_grade="not_run_image_only_no_candidate_predictions"))])
    println("Directional response panel: ",supporting,"/",length(rows)," supporting scans; complete=",complete,
        ". No molecular fit, deconvolution output or chemical grade.")
end

function execute(opts;loader=I.load_images,header_reader=path->I.read_sxm(path).header)
    input=inputs(opts)
    if haskey(opts,"--dry-run")
        println("Image-only response panel, no image arrays or output:")
        for p in input.panel; println(p.index,"/",p.cohort_size," ",p.file," base_dx=",input.shifts[p.file]); end
        return
    end
    BLAS.set_num_threads(1)
    out=opts["--outdir"]; mkpath(out)
    cp(opts["--settings"],joinpath(out,"settings.toml")); cp(opts["--config"],joinpath(out,"physical.toml"))
    cp(opts["--base-shifts"],joinpath(out,"base_shifts.tsv"))
    write_records(joinpath(out,"panel.tsv"),record.(input.panel))
    write_records(joinpath(out,"input_hashes.tsv"),[Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))) for k in ("--base-shifts","--config","--settings")])
    Threads.@threads :dynamic for p in input.panel
        one_scan(opts,input,p;loader,header_reader)
    end
    summarize(out,input.panel,input.s)
end

function main(args=ARGS)
    opts=parse_cli(args)
    opts===nothing && return println("diagnose_directional_response.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--dry-run]")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && DiagnoseDirectionalResponse.main()
