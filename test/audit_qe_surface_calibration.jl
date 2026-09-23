#!/usr/bin/env julia
# Read-only scientific diagnostic: never produces recognition templates.
module AuditQESurfaceCalibration
include(joinpath(@__DIR__, "lib", "qe_cube_molds.jl"))
using .QECubeMolds, LinearAlgebra, Statistics, Printf, SHA, TOML
const C = QECubeMolds.CCMoldNative
const INPUTS = ["--cube0", "--cube1", "--frame0", "--frame1", "--config", "--settings"]

function options(args)
    args == ["--help"] && return nothing
    o = Dict{String,String}()
    i = 1
    while i <= length(args)
        k = args[i]
        haskey(o,k) && error("Repeated option $k")
        if k == "--dry-run"; o[k] = "true"; i += 1; continue; end
        k in vcat(INPUTS,["--outdir"]) && i < length(args) && !startswith(args[i+1],"--") ||
            error("Missing/forbidden option $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o,k) for k in vcat(INPUTS,["--outdir"])) || error("All inputs required")
    all(isfile(o[k]) for k in INPUTS) || error("Missing input")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    C.load_settings(o["--config"]); QECubeMolds.settings(o["--settings"])
    o
end

"Numerical profiles and the full unchanged legacy isovalue response, not a calibration."
function audit(cube, frame, settings)
    z, v = C.surface_grid(cube,frame,settings)
    xy = C._numpy_grid(-settings.half_nm,settings.step_nm,
        round(Int,2settings.half_nm/settings.step_nm)+1)
    normal = C._normal(frame)
    inverse = inv(cube.axes_nm)
    profiles = NamedTuple[]
    for (k,h) in enumerate(z)
        finite = filter(isfinite,v[:,k])
        strict = 0
        for u in xy, t in xy
            q = inverse*(frame.origin_nm+t*frame.t_axis+u*frame.u_axis+h*normal-cube.origin_nm)
            strict += all(0 <= q[a] <= cube.dims[a]-1 for a in 1:3)
        end
        push!(profiles,(height_nm=h, finite=length(finite), strict_inside=strict,
            negative=count(<(0),finite), minimum=isempty(finite) ? NaN : minimum(finite),
            median=isempty(finite) ? NaN : median(finite), maximum=isempty(finite) ? NaN : maximum(finite),
            positive_sum=sum(x for x in finite if x>0; init=0.0),
            negative_abs_sum=-sum(x for x in finite if x<0; init=0.0)))
    end
    # Exactly the original 80-point log grid, including its endpoint convention.
    upper = log10(maximum(filter(isfinite,vec(v)))*settings.isovalue_max_fraction)
    step = (upper-settings.isovalue_min_log10)/(settings.isovalue_count-1)
    response = NamedTuple[]
    for i in 0:settings.isovalue_count-1
        iso = 10.0^(i == settings.isovalue_count-1 ? upper : settings.isovalue_min_log10+i*step)
        hs = Float64[]
        for row in axes(v,1)
            k = findlast(x->x>iso,@view v[row,:])
            k === nothing || push!(hs,z[k])
        end
        push!(response,(isovalue=iso, valid=length(hs),
            mean_height_nm=isempty(hs) ? NaN : mean(hs),
            above_target=count(>(settings.target_height_nm),hs)))
    end
    # Sample the declared target itself, not its nearest z-grid sample.
    points = reduce(vcat,[permutedims(frame.origin_nm+t*frame.t_axis+u*frame.u_axis+
        settings.target_height_nm*normal) for u in xy for t in xy])
    target = C.sample_volume(cube,points)
    finite = filter(isfinite,target)
    summary = Dict("target_height_nm"=>settings.target_height_nm,
        "target_valid"=>length(finite),"target_negative"=>count(<(0),finite),
        "target_minimum"=>(isempty(finite) ? NaN : minimum(finite)),
        "target_median"=>(isempty(finite) ? NaN : median(finite)),
        "target_maximum"=>(isempty(finite) ? NaN : maximum(finite)),
        "legacy_isovalue_lower"=>10.0^settings.isovalue_min_log10,
        "target_positive_sum"=>sum(x for x in finite if x>0;init=0.0),
        "target_negative_abs_sum"=>-sum(x for x in finite if x<0;init=0.0),
        "cell_step_extent_nm"=>collect(cube.dims).*norm.(eachcol(cube.axes_nm)),
        "frame_origin_nm"=>frame.origin_nm,"frame_normal"=>normal,
        "total_pixels"=>length(target),"total_native_negative"=>count(<(0),cube.values),
        "total_native_values"=>length(cube.values),"native_minimum"=>minimum(cube.values))
    (;profiles,response,summary,target)
end

# Diagnostic reproduction of QE 7.4.1 Modules/w0gauss.f90, n=-1. It is NOT
# used to replace any cube or spectrum. In stm.f90 a state below the lower
# bias-window edge has x=(down-energy)/degauss; negative tails are possible.
qe_cold_derivative(x) = exp(-min(200.0,(x-1/sqrt(2.0))^2))*(2-sqrt(2.0)*x)/sqrt(pi)

function table(path,rows)
    open(path,"w") do io
        println(io,join(string.(keys(first(rows))),'\t'))
        for r in rows
            println(io,join([x isa AbstractFloat ? @sprintf("%.17g",x) : string(x) for x in values(r)],'\t'))
        end
    end
end

function execute(o)
    VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 required")
    haskey(o,"--dry-run") && return println("Inputs checked; no cube calculation or output")
    out = abspath(o["--outdir"]); mkpath(out)
    s = C.load_settings(o["--config"])
    summaries = Dict{String,Any}()
    for typ in 0:1
        r = audit(QECubeMolds.read_qe_cube(o["--cube$typ"]),C.read_frame(o["--frame$typ"]),s)
        table(joinpath(out,"type$(typ)_profiles.tsv"),r.profiles)
        table(joinpath(out,"type$(typ)_legacy_response.tsv"),r.response)
        table(joinpath(out,"type$(typ)_target.tsv"),[(pixel=i,value=x) for (i,x) in enumerate(r.target)])
        summaries["type$typ"] = r.summary
        println("type$typ: ",r.summary)
    end
    summaries["scope"] = "diagnostic_only_no_templates_no_new_calibration_no_benchmark"
    summaries["inputs"] = Dict(k[3:end]=>Dict("path"=>abspath(o[k]),
        "sha256"=>bytes2hex(sha256(read(o[k])))) for k in INPUTS)
    open(io->TOML.print(io,summaries),joinpath(out,"summary.toml"),"w")
    summaries
end
function main(args=ARGS)
    o = options(args)
    o === nothing ? println("audit_qe_surface_calibration.jl --cube0 C --cube1 C --frame0 F --frame1 F --config ASSIGNMENT --settings ORDER --outdir NEW_DIR [--dry-run]") : execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && AuditQESurfaceCalibration.main()
