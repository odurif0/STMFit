#!/usr/bin/env julia
# One bounded raw-image counting diagnostic. No saved counts or benchmark input.
module LocalSigmaCounting
include(joinpath(@__DIR__, "refine_fixed_n_geometry.jl"))
module Pipeline
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
end
using SHA, TOML, LinearAlgebra
const R=FixedNGeometryRefinement
const G=R.G
const VP=R.VP
const E=R.F.Extractor
const P=Pipeline
const ARMS=("global_sigma_max","local_sigma_cap")
const REQUIRED=Set(["--data-dir","--count-config","--settings","--outdir"])
const COLUMNS=["file","repetition","arm","N","family","success","valid","reason",
    "gcv","rss","n_data","p_full","amp_min","amp_range","overlap","pair_overlap",
    "mean_overlap","endpoint_overrun_nm","residual_peak_snr","kappa_max_adj",
    "baseline","tilt_x","tilt_y","elapsed_s","parameters"]
const LOBES=["file","repetition","arm","N","family","lobe","amplitude","x_nm","y_nm",
    "t_nm","u_nm","sigma_parallel_nm","sigma_perp_nm","local_sigma_cap_nm"]
const CONTEXT=["file","raw_sha256","observed_pixels","total_pixels","n_data","noise",
    "origin_x_nm","origin_y_nm","axis_x","axis_y","support_tmin_nm","support_tmax_nm"]
const SELECTED=["file","repetition","arm","status","N_selected","source","gcv"]
const DIAG=["file","repetition","arm","N","family","start","global_status","global_evaluations",
    "global_error","lm_status","lm_converged","lm_iterations","lm_error","rss"]
cell(x)=replace(string(x),'\t'=>' ','\n'=>' ','\r'=>' ')
emit(io,cols,row)=(println(io,join([cell(get(row,k,"NA")) for k in cols],'\t'));flush(io))
rowdict(ps...)=Dict{String,Any}(ps...)
sha(path)=bytes2hex(sha256(read(path)))

function settings(path)
    s=TOML.parsefile(path)
    expected=Dict("model"=>Dict("arms"=>collect(ARMS),"profile"=>"gaussian",
        "native_elliptical_maxiter"=>50,"repetitions"=>2,"n_min"=>2,"n_max"=>14),
        "selection"=>Dict("policy"=>"exhaustive_valid_full_parameter_gcv",
            "tie_break"=>"smaller_N_then_elliptical","failure_policy"=>"retain_all_candidates_no_partial_grade"),
        "preprocessing"=>Dict("policy"=>"unchanged_native_fused","shared_data"=>true))
    s==expected || error("Only the declared counting comparison is supported")
    return s
end

function options(args)
    args==["--help"] && return nothing
    o=Dict{String,String}();i=1
    while i<=length(args)
        k=args[i];haskey(o,k) && error("Repeated option $k")
        if k=="--dry-run";o[k]="true";i+=1;continue;end
        k in union(REQUIRED,Set(["--chunk"])) && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden option $k")
        o[k]=args[i+1];i+=2
    end
    all(haskey(o,k) for k in REQUIRED) || error("All counting inputs required")
    isdir(o["--data-dir"]) && isfile(o["--count-config"]) && isfile(o["--settings"]) || error("Missing input")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    length(chunk)==2 && 1<=chunk[1]<=chunk[2] || error("Invalid chunk")
    settings(o["--settings"])
    o
end

function configs(raw,s,arm)
    arm in ARMS || error("Unknown arm")
    Set(keys(raw))==Set(["model","selection","preprocessing"]) || error("Unexpected count config sections")
    m=raw["model"]
    m["selection_criterion"]==m["cv_method"]=="gcv" || error("Full-parameter GCV required")
    m["chain_tilted_baseline"]===true || error("Frozen tilted baseline required")
    m["selection_policy"]=="support_midpoint_hybrid" || error("Only the frozen original input config is supported; its hybrid policy is NOT used")
    get(m,"peak_profile","gaussian")=="gaussian" && get(m,"shared_sigma_types",0)==0 &&
        get(m,"chain_spacing_model","free")=="free" || error("This comparison requires per-lobe Gaussian widths/free gaps")
    get(raw["preprocessing"],"missing_pixel_policy","median_fill")=="median_fill" || error("Native preprocessing required")
    for k in ("global_maxtime","global_maxiter","max_iter","multistart")
        haskey(m,k) && m[k]>0 || error("Explicit native numerical budget required: $k")
    end
    pcfg,ell,circ=E._configs(m,raw["preprocessing"],"unused")
    ell.sigma_parallel_min_nm==ell.sigma_perp_min_nm && ell.sigma_parallel_max_nm==ell.sigma_perp_max_nm || error("Circular nesting requires identical sigma bounds")
    for cfg in (ell,circ)
        cfg.overlap_constraint=arm;cfg.intelligent_sweep=false
        cfg.n_min=s["model"]["n_min"];cfg.n_max=s["model"]["n_max"]
        isempty(cfg.init_centers_t) && isempty(cfg.init_amplitudes) || error("No 1D initializer")
        G._effective_spacing_min_nm(cfg)
    end
    ell.skip_global=true;ell.multistart=1;ell.max_iter=s["model"]["native_elliptical_maxiter"]
    (;pcfg,ell,circ)
end

function data_bundle(img,pcfg,cfg)
    views=map(("fwd","bwd")) do direction
        found=filter(c->lowercase(c.name)==lowercase(pcfg.roi_channel) && lowercase(c.direction)==direction,img.channels)
        length(found)==1 || error("Exactly one actual $direction channel required; no directional fallback")
        only(found)
    end
    xs,ys,zimg,mask,x,y,z,noise=G._fused_roi_data(img,pcfg)
    axisfull=G._weighted_roi_axis(x,y,z)
    xf,yf,zf,axis,keep,support=G._chain_fit_data(x,y,z,axisfull,cfg)
    # This arm intentionally keeps native imputation; observation counts are
    # descriptive and never eligibility or N-selection inputs.
    f,b=views;stride=max(1,pcfg.stride)
    rawf=f.data[1:stride:end,1:stride:end];rawb=b.data[1:stride:end,1:stride:end]
    observed=count(isfinite.(rawf).&isfinite.(rawb))
    all(isfinite,zf) || error("Nonfinite fit data")
    return (;xs,ys,zimg,x=xf,y=yf,z=zf,zfull=z,noise,axisctx=axis,
        observed,total=length(rawf),support)
end

candidate_ns(d,cfg)=cfg.n_min:min(cfg.n_max,max(1,floor(Int,(d.axisctx.tmax-d.axisctx.tmin)/G._effective_spacing_min_nm(cfg))+1))
eligible(r)=r.success && r.valid && isfinite(r.gcv)
function select_candidate(candidates)
    pool=filter(q->eligible(q.result),candidates)
    isempty(pool) && return nothing
    first(sort(pool;by=q->(q.result.gcv,q.result.n,q.family=="ell" ? 0 : 1)))
end

function fit_family(n,d,cfg;warm_start=nothing,diagnostics=nothing)
    r=G._fit_chain_n(d.xs,d.ys,d.zimg,d.x,d.y,d.z,d.noise,n,d.axisctx,cfg;warm_start,diagnostics)
    VP.finalize_native!(r,d,cfg)
end

function save_fit(ios,key,r,d,cfg,elapsed)
    row=merge(copy(key),rowdict("success"=>r.success,"valid"=>r.valid,"reason"=>r.reason,
        "gcv"=>r.gcv,"rss"=>r.rss,"n_data"=>length(d.z),"p_full"=>G._chain_nparams(r.n,cfg),
        "amp_min"=>r.amp_min,"amp_range"=>r.amp_range,"elapsed_s"=>elapsed,
        "parameters"=>join(r.params,';')))
    for k in ("overlap","endpoint_overrun_nm","residual_peak_snr","kappa_max_adj")
        row[k]=getproperty(r,Symbol(k))
    end
    if r.success
        b,feats,ts,us,sp,sq=G._decode_chain(r.params,r.n,d.axisctx,cfg;amp_min=r.amp_min,amp_range=r.amp_range)
        row["baseline"]=b;row["tilt_x"]=cfg.chain_tilted_baseline ? r.params[2] : 0.
        row["tilt_y"]=cfg.chain_tilted_baseline ? r.params[3] : 0.
        row["pair_overlap"]=G._chain_pair_overlap(feats,sp,sq)
        row["mean_overlap"]=G._chain_overlap(feats,sum(sp)/r.n,sum(sq)/r.n)
        caps=G._chain_local_sigma_caps(r.n,diff(ts),cfg)
        for i in 1:r.n
            emit(ios["lobes"],LOBES,merge(copy(key),rowdict("lobe"=>i,"amplitude"=>feats[i].amplitude,
                "x_nm"=>feats[i].x_nm,"y_nm"=>feats[i].y_nm,"t_nm"=>ts[i],"u_nm"=>us[i],
                "sigma_parallel_nm"=>sp[i],"sigma_perp_nm"=>sq[i],"local_sigma_cap_nm"=>caps[i])))
        end
    end
    emit(ios["candidates"],COLUMNS,row)
end

function execute(o;reader=G.read_sxm,fitter=fit_family)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    Threads.nthreads()==1 || error("One Julia thread per shard required")
    BLAS.set_num_threads(1)
    s=settings(o["--settings"]);raw=TOML.parsefile(o["--count-config"])
    arms=Dict(a=>configs(raw,s,a) for a in ARMS)
    paths=P.raw_index(o["--data-dir"]);chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    files=[f for (i,f) in enumerate(sort(collect(keys(paths)))) if mod1(i,chunk[2])==chunk[1]]
    isempty(files) && error("Empty raw shard")
    haskey(o,"--dry-run") && return println("Metadata only: $(length(files))/$(length(paths)) scans, shard $(join(chunk,'/')); two exhaustive Gaussian arms, two repeats, no saved N or labels")
    out=abspath(o["--outdir"]);(ispath(out)||islink(out)) && error("Output exists")
    mkpath(joinpath(out,"fit_data"))
    cp(o["--settings"],joinpath(out,"settings.toml"));cp(o["--count-config"],joinpath(out,"count.toml"))
    for a in ARMS,f in ("circ","ell")
        cfg=getproperty(arms[a],Symbol(f))
        snapshot=Dict(string(k)=>(getfield(cfg,k) isa Symbol ? string(getfield(cfg,k)) : getfield(cfg,k)) for k in fieldnames(typeof(cfg)))
        open(io->TOML.print(io,snapshot),joinpath(out,"$a.$f.toml"),"w")
    end
    P.write_table(joinpath(out,"inputs.tsv"),["input","sha256"],
        [Dict("input"=>k,"sha256"=>sha(o[k])) for k in ("--count-config","--settings")])
    P.write_table(joinpath(out,"cohort.tsv"),["file","raw_sha256"],
        [Dict("file"=>f,"raw_sha256"=>sha(paths[f])) for f in files])
    schemas=Dict("candidates"=>COLUMNS,"lobes"=>LOBES,"context"=>CONTEXT,"selected"=>SELECTED,"diagnostics"=>DIAG)
    ios=Dict(k=>open(joinpath(out,"$k.tsv"),"w") for k in keys(schemas))
    for (k,cols) in schemas;println(ios[k],join(cols,'\t'));end
    failures=Dict{String,String}[];unavailable=0
    try
        for (i,file) in enumerate(files)
            try
                ac=arms[first(ARMS)];pcfg=deepcopy(ac.pcfg);pcfg.filepath=paths[file]
                img=reader(paths[file]);d=data_bundle(img,pcfg,ac.circ)
                axis=d.axisctx
                emit(ios["context"],CONTEXT,rowdict("file"=>file,"raw_sha256"=>sha(paths[file]),
                    "observed_pixels"=>d.observed,"total_pixels"=>d.total,"n_data"=>length(d.z),"noise"=>d.noise,
                    "origin_x_nm"=>axis.origin[1],"origin_y_nm"=>axis.origin[2],"axis_x"=>axis.axis[1],"axis_y"=>axis.axis[2],
                    "support_tmin_nm"=>axis.tmin,"support_tmax_nm"=>axis.tmax))
                open(joinpath(out,"fit_data",file*".tsv"),"w") do io
                    println(io,"x_nm\ty_nm\tz_nm")
                    for j in eachindex(d.z);println(io,join((d.x[j],d.y[j],d.z[j]),'\t'));end
                end
                for rep in 1:s["model"]["repetitions"],arm in (isodd(rep) ? ARMS : reverse(ARMS))
                    c=arms[arm];candidates=NamedTuple[]
                    for n in candidate_ns(d,c.circ)
                        rc=nothing
                        for family in ("circ","ell")
                            cfg=getproperty(c,Symbol(family));key=rowdict("file"=>file,"repetition"=>rep,"arm"=>arm,"N"=>n,"family"=>family)
                            started=time_ns()
                            callback=q->emit(ios["diagnostics"],DIAG,merge(copy(key),Dict(string(k)=>v for (k,v) in pairs(q))))
                            r=try
                                if family=="ell" && (rc===nothing || !rc.success)
                                    G.ChainModelResult(n=n,success=false,reason="circular initialization unavailable")
                                else
                                    warm=family=="circ" ? nothing : VP.circular_to_elliptical(rc.params,n,cfg)
                                    fitter(n,d,cfg;warm_start=warm,diagnostics=callback)
                                end
                            catch err
                                G.ChainModelResult(n=n,success=false,reason=cell(sprint(showerror,err)))
                            end
                            family=="circ" && (rc=r)
                            save_fit(ios,key,r,d,cfg,(time_ns()-started)/1e9)
                            push!(candidates,(result=r,family=family))
                        end
                    end
                    best=select_candidate(candidates)
                    key=rowdict("file"=>file,"repetition"=>rep,"arm"=>arm,"status"=>best===nothing ? "unavailable" : "ok")
                    if best===nothing
                        unavailable+=1
                    else
                        key["N_selected"]=best.result.n;key["source"]=best.family;key["gcv"]=best.result.gcv
                    end
                    emit(ios["selected"],SELECTED,key)
                    println("[$i/$(length(files))] $file repeat=$rep $arm N=$(get(key,"N_selected","NA")) candidates=$(length(candidates))");flush(stdout)
                end
            catch err
                push!(failures,Dict("file"=>file,"reason"=>cell(sprint(showerror,err))))
            end
        end
    finally
        foreach(close,values(ios))
    end
    isempty(failures) || P.write_table(joinpath(out,"failures.tsv"),["file","reason"],failures)
    isempty(failures) && unavailable==0 || error("Incomplete comparison: $(length(failures)) file failures, $unavailable unavailable selections; outputs retained, no grade")
    println("Complete shard, counting only; no chemical assignment or external grading.")
end

function main(args=ARGS)
    o=options(args)
    o===nothing ? println("diagnose_local_sigma_counting.jl --data-dir RAW --count-config TOML --settings TOML --outdir NEW [--chunk I/N] [--dry-run]\nTwo exhaustive counting arms repeated twice; no labels, saved N, assignment or production change.") : execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && LocalSigmaCounting.main()
