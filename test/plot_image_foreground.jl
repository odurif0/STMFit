#!/usr/bin/env julia
# One saved scan only; display limits never feed support decisions.
using Plots, Statistics
include(joinpath(@__DIR__,"verify_image_foreground.jl"))
const V=VerifyImageForeground
length(ARGS)==3 || error("plot_image_foreground.jl SAVED_SCAN_DIR NEW_PNG TITLE")
dir,out,title=ARGS; ispath(out) && error("Output exists"); panels=[]
for view in ("fwd","bwd")
    path=joinpath(dir,view); data=V.load_direction(path); xs=data.meta["xs_nm"]; ys=data.meta["ys_nm"]
    anchors=V.t(joinpath(path,"anchor_membership.tsv"))
    finite=vcat([z[isfinite.(z)] for z in (data.arrays["initial_smooth.f64"],data.arrays["joint_smooth.f64"])]...)
    limits=Tuple(quantile(finite,[.01,.99]))
    for prefix in ("initial","joint")
        z=data.arrays[prefix*"_smooth.f64"]; mask=data.arrays[prefix*"_core.u8"].==1
        p=heatmap(xs,ys,z;color=:grays,clims=limits,colorbar=false,aspect_ratio=:equal,
            title="$view / $prefix / $(count(mask)) bright pixels",xlabel="x (nm)",ylabel="y (nm)",legend=:topright)
        # Explicit boundary pixels avoid contour colour limits changing the heatmap.
        ny,nx=size(mask)
        edge=[CartesianIndex(y,x) for y in 1:ny for x in 1:nx if mask[y,x] &&
            (y==1 || y==ny || x==1 || x==nx || !all(mask[max(1,y-1):min(ny,y+1),max(1,x-1):min(nx,x+1)]))]
        scatter!(p,[xs[q[2]] for q in edge],[ys[q[1]] for q in edge];color=:cyan,markersize=.5,markerstrokewidth=0,label=false)
        for (retained,color,label) in ((true,:lime,"old anchor inside"),(false,:orange,"old anchor outside"))
            rows=filter(r->V.b(r,prefix*"_available") && V.b(r,prefix*"_core")==retained,anchors)
            scatter!(p,V.n.(rows,"x_nm"),V.n.(rows,"y_nm");color,label,markersize=2,markerstrokewidth=0)
        end
        push!(panels,p)
    end
end
savefig(plot(panels...;layout=(2,2),size=(1150,1000),
    plot_title=title*"; bright support, not molecules; 1-99% display only"),out)
