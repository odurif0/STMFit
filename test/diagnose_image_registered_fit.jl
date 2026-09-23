#!/usr/bin/env julia
# One fixed native scan: image-calibrated shift versus zero and free translation.
module ImageRegisteredFit
include(joinpath(@__DIR__,"diagnose_paired_shift.jl"))
include(joinpath(@__DIR__,"lib","image_registration_2d.jl"))
using .PairedShiftDiagnostic: table, value, write_records, CASES
using TOML, SHA, DelimitedFiles, LinearAlgebra
const D=PairedShiftDiagnostic
const H=D.H
const R=ImageRegistration2D
const OPTIONS=union(setdiff(D.OPTIONS,Set(["--shift-settings"])),Set(["--fit-settings","--image-dir"]))

function settings(path)
    t=TOML.parsefile(path)
    fields=Dict("model"=>["method","start_seeds","start_unit_box_radius","slsqp_maxeval","max_time_s"],
        "selection"=>["file","count_policy","start_policy","family_policy","image_policy","free_image_agreement_px",
            "free_fold_agreement_px","shift_interior_margin_px","holdout_policy","failure_policy"],
        "preprocessing"=>["input","support","window_and_folds","reference_frame","background_coordinates","noise"])
    Set(keys(t))==Set(keys(fields)) || error("Unknown fit section")
    for (section,keys_) in fields
        Set(keys(t[section]))==Set(keys_) || error("Unknown/missing fit field")
    end
    s=merge(values(t)...)
    for (key,val) in (("method","zero_image_fixed_and_free_translation"),("file","240817_006.sxm"),
        ("count_policy","reuse_native_diagnostic_N"),("start_policy","minimum_training_rss_then_existing_validity"),
        ("family_policy","valid_minimum_full_parameter_gcv"),("image_policy","require_all_image_folds_identified"),
        ("holdout_policy","both_folds_improve_paired_rss"),("failure_policy","explicit_no_retry_no_partial_grade"),
        ("input","verified_saved_registered_fit_pixels"),("support","reuse_all_observed_pixels_without_resampling"),
        ("window_and_folds","reuse_image_registration_settings"),("reference_frame","forward"),
        ("background_coordinates","observed_grid_in_both_views"),("noise","reuse_native_noise"))
        s[key]==val || error("Unexpected $key")
    end
    seeds=s["start_seeds"]
    length(seeds)>=2 && first(seeds)==0 && length(unique(seeds))==length(seeds) &&
        all(v->v isa Integer && !(v isa Bool) && v>=0,seeds) || error("Invalid deterministic starts")
    for key in ("start_unit_box_radius","slsqp_maxeval","max_time_s","free_image_agreement_px","free_fold_agreement_px","shift_interior_margin_px")
        v=s[key]; v isa Real && !(v isa Bool) && isfinite(v) && v>0 || error("Invalid $key")
    end
    s["start_unit_box_radius"]<.5 && s["slsqp_maxeval"] isa Integer || error("Invalid start radius or budget")
    return (;(Symbol(k)=>v for (k,v) in s)...)
end

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); old=String[]; i=1
    while i<=length(args)
        key=args[i]; haskey(opts,key) && error("Repeated option")
        if key=="--dry-run"; opts[key]="true"; push!(old,key); i+=1; continue; end
        key in union(OPTIONS,Set(["--chunk"])) && i<length(args) || error("Missing/forbidden image-fit option")
        opts[key]=args[i+1]
        key in ("--fit-settings","--image-dir","--solver-settings") || append!(old,args[i:i+1]); i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All image-fit inputs required")
    D.D.parse_cli(old); settings(opts["--fit-settings"]); H.S.settings(opts["--solver-settings"])
    for file in ("settings.toml","physical.toml","summary.tsv","eligibility.tsv","input_hashes.tsv","fold0/images.tsv")
        isfile(joinpath(opts["--image-dir"],file)) || error("Missing image calibration $file")
    end
    return opts
end

function inputs(opts)
    input=D.D.inputs(opts); s=settings(opts["--fit-settings"]); dir=opts["--image-dir"]
    image_settings=R.settings(joinpath(dir,"settings.toml")); image_rows=table(joinpath(dir,"summary.tsv"))
    parse.(Int,getindex.(image_rows,"fold"))==[0,1,2] || error("Missing image folds")
    all(r["file"]==input.file && r["accepted"]=="true" for r in image_rows) || error("Image registration is not identified")
    runs=[(fold=parse(Int,r["fold"]),accepted=true,peaks=[(dx=value(r,"dx"),dy=value(r,"dy"))]) for r in image_rows]
    eligibility=R.eligibility(runs,image_settings)
    e=only(table(joinpath(dir,"eligibility.tsv")))
    eligibility.eligible && e["eligible"]=="true" && value(e,"spread_px")==eligibility.spread_px || error("Image folds disagree")
    read(joinpath(dir,"physical.toml"))==read(opts["--physical-config"]) || error("Image/native preprocessing differs")
    old=table(joinpath(opts["--native-dir"],"input_hashes.tsv")); ih=table(joinpath(dir,"input_hashes.tsv"))
    only(r for r in old if r["input"]=="--shifts")["sha256"]==only(r for r in ih if r["input"]=="--base-shifts")["sha256"] || error("Base shift input differs")
    step=H.grid_step(input.pixels); firstrow=first(image_rows)
    all(isapprox.(step,(value(firstrow,"step_x_nm"),value(firstrow,"step_y_nm"));rtol=1e-10)) || error("Image/native grids differ")
    limits=(parse(Int,firstrow["limit_x_px"]),parse(Int,firstrow["limit_y_px"]))
    all(r["base_dx"]==firstrow["base_dx"] && r["base_dy"]=="0" &&
        r["limit_x_px"]==firstrow["limit_x_px"] && r["limit_y_px"]==firstrow["limit_y_px"] for r in image_rows) || error("Image frames differ")
    base=parse(Int,firstrow["base_dx"])
    ny,nx=parse(Int,firstrow["nrows"]),parse(Int,firstrow["ncols"])
    limits==Tuple(min(floor(Int,image_settings.residual_window_nm/step[j]+1e-10),floor(Int,image_settings.max_lag_width_fraction*(j==1 ? nx : ny))) for j in 1:2) || error("Image window differs")
    m=readdlm(joinpath(dir,"fold0","images.tsv"),'\t',Float64;skipstart=1)
    size(m)==(ny*nx,6) && m[:,1]==repeat(1:ny;outer=nx) && m[:,2]==repeat(1:nx;inner=ny) || error("Invalid saved image grid")
    a=reshape(m[:,5],ny,nx); b=reshape(m[:,6],ny,nx)
    rows=parse.(Int,getindex.(input.pixels,"row")); cols=parse.(Int,getindex.(input.pixels,"column"))
    for (j,p) in enumerate(input.pixels)
        y,x=rows[j],cols[j]
        isapprox(a[y,x],value(p,"fwd_nm");atol=1e-12,rtol=1e-10) &&
            isapprox(b[y,x+base],value(p,"bwd_nm");atol=1e-12,rtol=1e-10) || error("Image/native observations differ")
        i=y+(x-1)*ny
        isapprox(m[i,3],value(p,"x_nm");atol=1e-10) && isapprox(m[i,4],value(p,"y_nm");atol=1e-10) || error("Image/native coordinates differ")
    end
    all(parse(Int,r["dx_px"])==base for r in input.fits if r["arm"]=="registered") || error("Native frame shift differs")
    folds=H.folds(rows,cols,image_settings.holdout_block_px,image_settings.holdout_buffer_px)
    solver=merge(H.S.settings(opts["--solver-settings"]),(slsqp_maxeval=s.slsqp_maxeval,max_time_s=s.max_time_s,evaluation_checkpoints=[1,s.slsqp_maxeval]))
    return merge(input,(settings=s,image_settings=image_settings,image_rows=image_rows,step=step,limits=limits,folds=folds,solver=solver))
end

function problem(ctx,mode,model,step,limits,image_shift;indices=eachindex(ctx.data.z))
    mode in ("paired","image_fixed","shifted") || error("Unknown image-fit mode")
    mode=="paired" && return H.problem(ctx,"paired",model,(shift_max_px=1.,),step;indices)
    p=H.problem(ctx,"shifted",model,(shift_max_px=Float64(maximum(limits)),),step;indices)
    p.lo[end-1:end].=-collect(limits); p.hi[end-1:end].=collect(limits)
    mode=="shifted" && return p
    all(-limits[j]<=image_shift[j]<=limits[j] for j in 1:2) || error("Image estimate outside fixed window")
    return merge(p,(q0=p.q0[1:end-2],lo=p.lo[1:end-2],hi=p.hi[1:end-2],
        predict=q->p.predict(vcat(q,collect(image_shift))),allpred=q->p.allpred(vcat(q,collect(image_shift)))))
end

function finish(solved,p,ctx,mode)
    fit=H.finish(solved,p,ctx,"paired")
    # Two image coordinates were estimated from data even though not optimized here.
    np=length(fit.params)+(mode=="image_fixed" ? 2 : 0)
    return merge(fit,(;np,gcv=fit.nd>np ? fit.nd/(fit.nd-np)^2*fit.rss : Inf))
end

function execute(opts)
    BLAS.set_num_threads(1); input=inputs(opts)
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/')); cases=CASES[chunk[1]:chunk[2]:end]
    for c in cases; D.D.context(input,c); end
    haskey(opts,"--dry-run") && return println("Image-calibrated fit dry run: $(input.file), $(length(cases)) families, $(length(input.pixels)) unchanged observed pixels; no fitting/output.")
    out=opts["--outdir"]; mkpath(out); hashes=Dict{String,String}[]
    for key in ("--physical-config","--assignment-config","--model-settings","--settings","--solver-settings","--fit-settings")
        cp(opts[key],joinpath(out,replace(key,"--"=>"")*".toml"))
        push!(hashes,Dict("input"=>key,"sha256"=>bytes2hex(sha256(read(opts[key])))))
    end
    for (tag,key,files) in (("native","--native-dir",["fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv"]),
        ("reference","--reference-dir",["fits.tsv","parameters.tsv","native_hashes.tsv","paired_settings.toml"]),
        ("image","--image-dir",["settings.toml","physical.toml","summary.tsv","eligibility.tsv","input_hashes.tsv","fold0/images.tsv"]))
        for file in files
            push!(hashes,Dict("input"=>"$tag/$file","sha256"=>bytes2hex(sha256(read(joinpath(opts[key],file))))))
        end
    end
    write_records(joinpath(out,"input_hashes.tsv"),hashes)
    cp(joinpath(opts["--image-dir"],"summary.tsv"),joinpath(out,"image_shifts.tsv"))
    cp(joinpath(opts["--image-dir"],"settings.toml"),joinpath(out,"image_settings.toml"))
    fits=Dict{String,String}[]; failures=Dict{String,String}[]
    for case in cases
        ctx=D.D.context(input,case)
        for fold in 0:2
            ix=fold==0 ? eachindex(ctx.data.z) : input.folds[fold].train
            ir=input.image_rows[fold+1]; image_shift=(value(ir,"dx"),value(ir,"dy"))
            zero=problem(ctx,"paired",input.model_options,input.step,input.limits,image_shift;indices=ix)
            scale=max(sum(abs2,zero.predict(zero.q0).-zero.target),input.options.objective_scale_floor)
            for mode in ("paired","image_fixed","shifted")
                p=problem(ctx,mode,input.model_options,input.step,input.limits,image_shift;indices=ix)
                length(p.target)>length(p.q0)+2 || error("Too few training observations")
                for seed in input.settings.start_seeds
                    id="$(case.profile).$(case.family).$mode.$fold.$seed"
                    try
                        println("START $id"); flush(stdout)
                        q0=H.start(p,seed,input.settings.start_unit_box_radius)
                        fit=finish(H.S.solve(merge(p,(;q0)),input.solver,input.options;objective_scale=scale),p,ctx,mode)
                        audit=H.C.stationarity(merge(fit,(initial_rss=scale,)),input.options)
                        delta=mode=="paired" ? (0.,0.) : mode=="image_fixed" ? image_shift : Tuple(fit.params[end-1:end])
                        extra=(shift_x_px=delta[1],shift_y_px=delta[2],image_dx_px=image_shift[1],image_dy_px=image_shift[2],fitted_np=length(fit.params))
                        push!(fits,D.record(out,input,case,ctx,p,fit,audit,mode,fold,seed,scale;extra))
                        println("DONE $id rss=$(fit.rss) valid=$(fit.valid) shift=$delta"); flush(stdout)
                    catch err
                        push!(failures,Dict("id"=>id,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')))
                        println(stderr,"FAILED $id: ",sprint(showerror,err)); flush(stderr)
                    end
                end
            end
        end
    end
    write_records(joinpath(out,"fits.tsv"),fits); write_records(joinpath(out,"failures.tsv"),failures)
    isempty(failures) || error("Image-fit failures; no retry or partial grade")
end

function summaries(fits,s,limits)
    winners=Dict{String,String}[]; pilots=Dict{String,String}[]
    for c in CASES,m in ("paired","image_fixed","shifted"),f in 0:2
        rows=filter(r->r["profile"]==c.profile && r["family"]==c.family && r["mode"]==m && parse(Int,r["fold"])==f,fits)
        length(rows)==length(s.start_seeds) && Set(parse.(Int,getindex.(rows,"seed")))==Set(s.start_seeds) || error("Incomplete start group")
        push!(winners,first(sort(rows;by=r->(value(r,"rss"),parse(Int,r["seed"])))))
    end
    for c in CASES,m in ("image_fixed","shifted")
        row(mode,f)=only(r for r in winners if r["profile"]==c.profile && r["family"]==c.family && r["mode"]==mode && parse(Int,r["fold"])==f)
        controls=[row("paired",f) for f in 0:2]; selected=[row(m,f) for f in 0:2]
        gains=[value(controls[i],"heldout_rss")-value(selected[i],"heldout_rss") for i in 2:3]
        interior=all(abs(value(r,k))<limits[j]-s.shift_interior_margin_px for r in selected for (j,k) in enumerate(("shift_x_px","shift_y_px")))
        spread=maximum(abs(value(a,k)-value(b,k)) for a in selected,b in selected,k in ("shift_x_px","shift_y_px"))
        disagreement=maximum(abs(value(r,k)-value(r,ik)) for r in selected for (k,ik) in (("shift_x_px","image_dx_px"),("shift_y_px","image_dy_px")))
        eligible=first(selected)["valid"]=="true" && value(first(selected),"gcv")<value(first(controls),"gcv") && all(>(0),gains) &&
            interior && spread<=s.free_fold_agreement_px && disagreement<=s.free_image_agreement_px
        push!(pilots,Dict("profile"=>c.profile,"family"=>c.family,"mode"=>m,"selected_id"=>first(selected)["id"],"full_valid"=>first(selected)["valid"],
            "full_gcv_gain"=>string(value(first(controls),"gcv")-value(first(selected),"gcv")),"fold1_rss_gain"=>string(gains[1]),"fold2_rss_gain"=>string(gains[2]),
            "shift_interior"=>string(interior),"shift_spread_px"=>string(spread),"image_disagreement_px"=>string(disagreement),"eligible"=>string(eligible)))
    end
    return (;winners,pilots)
end

function merge_chunks(out)
    dirs=[joinpath(out,"chunk$i") for i in 1:4]; fits=Dict{String,String}[]
    names=["input_hashes.tsv","physical-config.toml","assignment-config.toml","model-settings.toml","settings.toml","solver-settings.toml","fit-settings.toml","image_shifts.tsv","image_settings.toml"]
    for dir in dirs
        isfile(joinpath(dir,"failures.tsv")) && error("Shard failed")
        append!(fits,table(joinpath(dir,"fits.tsv")))
        all(read(joinpath(dir,n))==read(joinpath(first(dirs),n)) for n in names) || error("Shard inputs differ")
    end
    s=settings(joinpath(first(dirs),"fit-settings.toml")); ir=first(table(joinpath(first(dirs),"image_shifts.tsv")))
    limits=(parse(Int,ir["limit_x_px"]),parse(Int,ir["limit_y_px"]))
    length(fits)==length(unique(getindex.(fits,"id")))==36length(s.start_seeds) || error("Incomplete fit schedule")
    result=summaries(fits,s,limits)
    for n in names; cp(joinpath(first(dirs),n),joinpath(out,n)); end
    write_records(joinpath(out,"fits.tsv"),fits); write_records(joinpath(out,"selected_starts.tsv"),result.winners)
    write_records(joinpath(out,"pilot.tsv"),result.pilots)
    for mode in ("image_fixed","shifted")
        eligible=all(any(r["profile"]==profile && r["mode"]==mode && r["eligible"]=="true" for r in result.pilots) for profile in ("gaussian","split"))
        println("All $(length(fits)) endpoints accounted for; $mode full Gaussian/split eligibility=$eligible. No recognition grade.")
    end
end

function main(args=ARGS)
    length(args)==2 && first(args)=="--merge" && return merge_chunks(last(args))
    opts=parse_cli(args)
    opts===nothing && return println("diagnose_image_registered_fit.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--chunk I/N] [--dry-run]; or --merge RUN_DIR")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ImageRegisteredFit.main()
