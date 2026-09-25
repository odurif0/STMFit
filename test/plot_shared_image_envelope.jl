#!/usr/bin/env julia
# Saved one-file diagnostic only. Display limits never enter fitting or scores.
using Plots, Statistics, TOML
include(joinpath(@__DIR__,"diagnose_shared_image_envelope.jl"))
const D=DiagnoseSharedImageEnvelope
length(ARGS)==3 || error("plot_shared_image_envelope.jl FOREGROUND_SCAN SAVED_SCAN NEW_PNG")
foreground,saved,out=ARGS; ispath(out) && error("Output exists"); panels=[]
for view in ("fwd","bwd")
    dir=joinpath(saved,view); a=D.load_source(joinpath(foreground,view))
    chosen=TOML.parsefile(joinpath(dir,"selection.toml"))["chosen"]
    ny,nx=size(a.raw); source=D.D.get_array(dir,"source_corrected.f64",Float64,ny,nx)
    envelope=D.D.get_array(dir,chosen*".f64",Float64,ny,nx); residual=source-envelope
    common=Tuple(quantile(source[isfinite.(source)],[.01,.99]))
    scale=quantile(abs.(residual[isfinite.(residual)]),.99)
    for (title,z,limits,colors) in (("source",source,common,:grays),(chosen,envelope,common,:grays),
        ("source - envelope",residual,(-scale,scale),:RdBu))
        p=heatmap(a.meta["xs_nm"],a.meta["ys_nm"],1000z;clims=1000 .* limits,color=colors,
            title="$view / $title",xlabel="x (nm)",ylabel="y (nm)",colorbar_title="pm",aspect_ratio=:equal)
        push!(panels,p)
    end
end
savefig(plot(panels...;layout=(2,3),size=(1500,1000),plot_title="Class-free envelope; display clipping only; no units/chemistry"),out)
