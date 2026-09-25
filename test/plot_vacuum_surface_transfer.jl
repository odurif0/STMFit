#!/usr/bin/env julia
# One-file visual check of saved diagnostic anchors, never chemical validation.
using STMSXMIO, Plots, Statistics, TOML
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow
length(ARGS)==3 || error("plot_vacuum_surface_transfer.jl SAVED_RUN ONE_RAW_SXM NEW_PNG")
run,raw,out=ARGS
ispath(out) && error("Output exists")
c=TOML.parsefile(joinpath(run,"settings.toml")); img=read_sxm(raw)
dir=joinpath(run,splitext(basename(raw))[1]); panels=[]
for view in ("fwd","bwd")
    ch=only(filter(q->q.name=="Z" && q.direction==view,img.channels))
    obs=STMSXMIO.preprocess_observed_channel(img,ch;stride=1,flatten="plane",smooth_radius_px=0,
        plane_rank_rtol=c["preprocessing"]["plane_rank_rtol"])
    anchors=V.tsv(joinpath(dir,view,"anchors.tsv")); patches=V.tsv(joinpath(dir,view,"patches.tsv"))
    display_limits=Tuple(quantile(obs.z[isfinite.(obs.z)],[.01,.99]))
    for iso in 1:3
        selected=filter(r->r["arm"]=="chemical" && parse(Int,r["interval"])==iso,patches)
        p=heatmap(obs.xs,obs.ys,obs.z;color=:grays,clims=display_limits,colorbar=false,
            aspect_ratio=:equal,title="$view / iso $iso / $(length(selected)) patches",
            xlabel="x (nm)",ylabel="y (nm)",legend=:topright)
        for row in selected
            anchor=anchors[parse(Int,row["anchor"])]; x=parse(Float64,anchor["x"]); y=parse(Float64,anchor["y"])
            type=parse(Int,row["type"]); color=type==0 ? :royalblue : :orangered
            angle=range(0,2pi;length=65); radius=c["selection"]["patch_radius_nm"]
            plot!(p,x.+radius.*cos.(angle),y.+radius.*sin.(angle);color,linewidth=.8,label=false)
        end
        scatter!(p,[NaN],[NaN];color=:royalblue,label="GlcN hypothesis",markersize=3)
        scatter!(p,[NaN],[NaN];color=:orangered,label="GlcNAc hypothesis",markersize=3)
        push!(panels,p)
    end
end
figure=plot(panels...;layout=(2,3),size=(1400,950),
    plot_title="$(basename(raw)): diagnostic extrema, not unit assignments; grayscale 1-99% for display only")
savefig(figure,out)
