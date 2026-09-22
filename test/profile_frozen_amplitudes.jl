#!/usr/bin/env julia
# Opt-in fixed-geometry experiment. No nonlinear fit, selection or label input.
module FrozenAmplitudeProfile
using GaussianFit2D, TOML, LinearAlgebra, Printf
include(joinpath(@__DIR__, "lib", "counting_variable_projection.jl"))
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment: lobe_table, write_table
const G = GaussianFit2D
const VP = CountingVariableProjection
module Extractor
include(joinpath(@__DIR__, "extract_lobe_features.jl"))
end

const MUTABLE_COLUMNS = Set(["amplitude", "baseline", "tilt_x", "tilt_y", "amp_rel", "gcv"])
const REQUIRED = ["N", "amplitude", "x_nm", "y_nm", "t_nm", "u_nm", "sigma_parallel_nm",
    "sigma_perp_nm", "skew_ratio", "axis_x", "axis_y", "origin_x_nm", "origin_y_nm",
    "baseline", "tilt_x", "tilt_y", "amp_rel", "gcv", "source"]
const PATH_OPTIONS = Set(["--features", "--data-dir", "--config", "--settings", "--out"])
number(row,key) = parse(Float64,row[key])
fmt(x) = @sprintf("%.17g", x)

function parse_cli(args)
    if "--help" in args
        println("""
        Julia 1.13: profile_frozen_amplitudes.jl --features TSV --data-dir DIR
          --config COUNT.toml --settings DIAGNOSTIC.toml --out NEW.tsv
          [--chunk I/N] [--dry-run]
        Solve only the native bounded amplitudes and tilted background of cached
        Gaussian geometry. N, positions, widths, axes, skew and split cache stay
        fixed. Reuse [counting_variable_projection] linear controls; NO outer search.
        All rows must succeed. Failures remain explicit and stop downstream use.
        Multi-file real processing belongs on Viper; dry-run reads metadata only.
        """)
        return nothing
    end
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        key=args[i]
        haskey(opts,key) && error("Repeated option $key")
        if key=="--dry-run"
            opts[key]="true"; i+=1; continue
        end
        key in union(PATH_OPTIONS,Set(["--chunk"])) || error("Unknown or forbidden option $key")
        i<length(args) && !startswith(args[i+1],"--") || error("Missing value $key")
        opts[key]=args[i+1]; i+=2
    end
    all(haskey(opts,k) for k in PATH_OPTIONS) || error("Missing required input")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    isdir(opts["--data-dir"]) || error("Missing raw directory")
    for k in ("--features","--config","--settings")
        isfile(opts[k]) || error("Missing input $k")
    end
    for suffix in ("", ".audit.tsv", ".coefficients.tsv", ".failures.tsv", ".fit_data")
        p=opts["--out"]*suffix
        (ispath(p)||islink(p)) && error("Output exists: $p")
    end
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    length(chunk)==2 && 1<=chunk[1]<=chunk[2] || error("Invalid chunk")
    return opts
end

function validate_chain(rows)
    n=length(rows); n>0 || error("Empty chain")
    all(parse(Int,r["N"])==n for r in rows) || error("Inconsistent N")
    [parse(Int,r["lobe"]) for r in rows]==collect(1:n) || error("Noncontiguous lobes")
    for c in ("axis_x","axis_y","origin_x_nm","origin_y_nm","baseline","tilt_x","tilt_y","gcv","source")
        length(unique(r[c] for r in rows))==1 || error("Inconsistent $c")
    end
    for r in rows, c in setdiff(REQUIRED,["source"])
        isfinite(number(r,c)) || error("Nonfinite $c")
    end
    all(r["source"] in ("circ","ell") for r in rows) || error("Unknown geometry source")
    all(number(r,"skew_ratio")==1 for r in rows) || error("Only cached Gaussian base geometry is profiled")
    all(number(r,c)>0 for r in rows for c in ("amplitude","sigma_parallel_nm","sigma_perp_nm")) || error("Nonpositive geometry/amplitude")
    return n
end

# The cached decimal values are the experiment's fixed geometry. Do not decode
# new centers, renormalize their axis or rerun any geometry optimizer.
function frozen_design(x,y,rows,tilted)
    validate_chain(rows)
    n=length(rows); p=tilted ? 3 : 1
    A=ones(length(x),p+n)
    length(x)==length(y) || throw(DimensionMismatch())
    all(isfinite,x) && all(isfinite,y) || error("Nonfinite fit coordinates")
    if tilted; A[:,2].=x; A[:,3].=y; end
    ax,ay=number(rows[1],"axis_x"),number(rows[1],"axis_y")
    for (j,r) in enumerate(rows)
        cx,cy=number(r,"x_nm"),number(r,"y_nm")
        sp,sq=number(r,"sigma_parallel_nm"),number(r,"sigma_perp_nm")
        A[:,p+j] .= @. exp(-0.5*((((x-cx)*ax+(y-cy)*ay)/sp)^2 + (((x-cx)*(-ay)+(y-cy)*ax)/sq)^2))
    end
    return A
end

# Exact native finite parameter bounds, mapped through the native decoder.
function coefficient_bounds(n,cfg,axis,zimg)
    ampmax=max(maximum(zimg),G.EPS)
    amin=cfg.min_amplitude_fraction*ampmax
    arange=max(ampmax-amin,G.EPS)
    lo,hi=VP.native_raw_bounds(n,cfg)
    lower=VP.linear_coefficients(lo,n,axis,cfg; amp_min=amin,amp_range=arange)
    upper=VP.linear_coefficients(hi,n,axis,cfg; amp_min=amin,amp_range=arange)
    return lower,upper
end

"Half the last written decimal unit, not a fitted feasibility tolerance."
function decimal_radius(s)
    parts=split(lowercase(s),'e'); mant=parts[1]
    exponent=length(parts)==2 ? parse(Int,parts[2]) : 0
    fraction=occursin('.',mant) ? length(split(mant,'.')[2]) : 0
    return 0.5*10.0^(exponent-fraction)
end

function profile_rows(rows,x,y,z,lower,upper,cfg,options)
    n=validate_chain(rows); tilted=cfg.chain_tilted_baseline
    fields=tilted ? ["baseline","tilt_x","tilt_y"] : ["baseline"]
    texts=vcat([rows[1][c] for c in fields],[r["amplitude"] for r in rows])
    initial=parse.(Float64,texts)
    length(initial)==length(lower)==length(upper) || throw(DimensionMismatch())
    # Rounded saved coefficients can be outside a native endpoint by at most
    # half their last printed decimal unit. Do not widen the native bounds.
    for i in eachindex(initial)
        tol=decimal_radius(texts[i])+eps(max(abs(initial[i]),abs(lower[i]),abs(upper[i])))
        lower[i]-tol<=initial[i]<=upper[i]+tol || error("Cached coefficient outside native bounds: $i")
    end
    A=frozen_design(x,y,rows,tilted)
    result=VP.box_lsq(A,z,lower,upper; x0=initial,options)
    result.converged || error("Linear profile did not converge: $(result.status)")
    oldrss=sum(abs2,A*initial-z)
    result.rss<=oldrss+options.mapping_atol+options.mapping_rtol*oldrss || error("Profile increased RSS")
    pfull=G._chain_nparams(n,cfg)
    length(z)>pfull || error("Insufficient fit pixels for full-parameter GCV")
    gcv=length(z)/(length(z)-pfull)^2*result.rss
    out=deepcopy(rows); prefix=length(fields); amps=result.x[prefix+1:end]
    for (i,r) in enumerate(out)
        for (j,c) in enumerate(fields); r[c]=fmt(result.x[j]); end
        r["amplitude"]=fmt(amps[i]); r["amp_rel"]=@sprintf("%.6f",amps[i]/maximum(amps))
        r["gcv"]=fmt(gcv)
        all(r[c]==rows[i][c] for c in setdiff(keys(r),MUTABLE_COLUMNS)) || error("Frozen column changed")
    end
    return (rows=out,result=result,initial=initial,initial_rss=oldrss,pfull=pfull,
        initial_clamp=maximum(abs.(clamp.(initial,lower,upper)-initial)),gcv=gcv)
end

function execute(opts)
    BLAS.set_num_threads(1)
    raw=TOML.parsefile(opts["--config"])
    settings=TOML.parsefile(opts["--settings"])
    options=VP.read_options(settings["counting_variable_projection"])
    pcfg,cfg,_=Extractor._configs(raw["model"],raw["preprocessing"],dirname(opts["--out"]))
    G._chain_peak_profile(cfg)==:gaussian || error("Gaussian base profile required")
    startswith(get(raw["model"],"selection_policy",""),"adaptive_support") && error("Adaptive support cache is outside this experiment")
    header,table=lobe_table(opts["--features"];required=REQUIRED)
    files=sort(unique(first.(collect(keys(table)))))
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    files=[f for (i,f) in enumerate(files) if mod1(i,chunk[2])==chunk[1]]
    isempty(files) && error("Empty chunk")
    chains=Dict(f=>[table[k] for k in sort(collect(keys(table))) if k[1]==f] for f in files)
    for f in files
        validate_chain(chains[f]); isfile(joinpath(opts["--data-dir"],f)) || error("Missing raw $f")
    end
    haskey(opts,"--dry-run") && return println("Frozen base profiles: $(length(files)) files; no pixels read, no fit or output.")
    mkdirpath=opts["--out"]*".fit_data"; mkpath(mkdirpath)
    outputs=Dict{String,String}[]; audit=Dict{String,String}[]; coefficients=Dict{String,String}[]; failures=Dict{String,String}[]
    for (idx,f) in enumerate(files)
        try
            rs=chains[f]; n=length(rs); pcfg.filepath=joinpath(opts["--data-dir"],f)
            img=G.read_sxm(pcfg.filepath)
            xs,ys,zimg,mask,x,y,z,noise=G._fused_roi_data(img,pcfg)
            full=G._weighted_roi_axis(x,y,z)
            xf,yf,zf,axis,keep,support=G._chain_fit_data(x,y,z,full,cfg)
            # Same unmodified source/preprocessing reconstructs native support;
            # reject cache/input mismatch before the linear solve.
            for (key,value,format) in (("axis_x",axis.axis[1],"axis"),("axis_y",axis.axis[2],"axis"),
                ("origin_x_nm",axis.origin[1],"origin"),("origin_y_nm",axis.origin[2],"origin"))
                printed=format=="axis" ? @sprintf("%.8f",value) : @sprintf("%.6f",value)
                printed==rs[1][key] || error("Cached $key differs from current input/support")
            end
            localcfg=deepcopy(cfg); localcfg.chain_circular_sigmas=rs[1]["source"]=="circ"
            lower,upper=coefficient_bounds(n,localcfg,axis,zimg)
            r=profile_rows(rs,xf,yf,zf,lower,upper,localcfg,options)
            open(joinpath(mkdirpath,f*".tsv"),"w") do io
                println(io,"x_nm\ty_nm\tz_nm")
                for j in eachindex(zf); println(io,join(fmt.((xf[j],yf[j],zf[j])),'\t')); end
            end
            append!(outputs,r.rows)
            names=vcat(localcfg.chain_tilted_baseline ? ["baseline","tilt_x","tilt_y"] : ["baseline"], ["amplitude_$j" for j in 1:n])
            for j in eachindex(names)
                push!(coefficients,Dict("file"=>f,"parameter"=>names[j],"initial"=>fmt(r.initial[j]),
                    "value"=>fmt(r.result.x[j]),"lower"=>fmt(lower[j]),"upper"=>fmt(upper[j]),"gradient"=>fmt(r.result.gradient[j])))
            end
            values=(file=f,N=n,n_data=length(zf),p_full=r.pfull,initial_rss=r.initial_rss,rss=r.result.rss,
                gcv=r.gcv,cached_gcv=rs[1]["gcv"],kkt_violation=r.result.kkt_violation,kkt_tolerance=r.result.kkt_tolerance,
                iterations=r.result.iterations,status=r.result.status,initial_clamp=r.initial_clamp,
                active_lower=join(r.result.active_lower,','),active_upper=join(r.result.active_upper,','),
                support_tmin=axis.tmin,support_tmax=axis.tmax,noise=noise,
                amp_max_data=max(maximum(zimg),G.EPS))
            push!(audit,Dict(string(k)=>string(v) for (k,v) in pairs(values)))
            println("[$idx/$(length(files))] $f N=$n RSS $(r.initial_rss) -> $(r.result.rss), KKT $(r.result.status)"); flush(stdout)
        catch err
            push!(failures,Dict("file"=>f,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')))
            println(stderr,"FAILED $f: ",sprint(showerror,err))
        end
    end
    !isempty(audit) && write_table(opts["--out"]*".audit.tsv",sort(collect(keys(first(audit)))),audit)
    !isempty(coefficients) && write_table(opts["--out"]*".coefficients.tsv",["file","parameter","initial","value","lower","upper","gradient"],coefficients)
    if !isempty(failures)
        write_table(opts["--out"]*".failures.tsv",["file","reason"],failures)
        error("$(length(failures)) frozen profiles failed; no incomplete candidate table emitted")
    end
    write_table(opts["--out"],header,outputs)
end

function main(args=ARGS)
    opts=parse_cli(args); opts===nothing || execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && FrozenAmplitudeProfile.main()
