#!/usr/bin/env julia
# Descriptive saved-table plot only; no reference selection, trend fit or mask.
using Plots, Statistics
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow
length(ARGS)==3 || error("plot_qe_substrate_reference.jl SAVED_RUN VERIFIED_REPORT NEW_PNG")
run,report,out=ARGS
ispath(out) && error("Output already exists")
n(r,k)=parse(Float64,r[k]); i(r,k)=parse(Int,r[k])
grid=Dict(i(r,"pixel")=>r for r in V.tsv(joinpath(run,"grid.tsv")))
eligible=count(r->parse(Bool,r["eligible"]),values(grid))
rows=filter(r->r["domain"]=="cu_half_gap",V.tsv(joinpath(report,"paired.tsv")))
isempty(rows) && error("No paired Cu-half roots; no plot or substituted domain")
clearance(r)=min(n(grid[i(r,"pixel")],"glcn_margin_nm"),n(grid[i(r,"pixel")],"glcnac_margin_nm"))
y=1000n.(rows,"referenced_contrast_nm"); raw=1000n.(rows,"original_molecular_contrast_nm")
lo,hi=extrema(vcat(y,raw,[0.])); pad=.05*(hi-lo)
panels=[]
for iso in 1:3
    rr=filter(r->i(r,"interval")==iso,rows)
    p=plot(;title="isovalue $iso: $(length(rr)) / $eligible paired sites",xlabel="minimum projected PAW clearance (nm)",
        ylabel="hypothetical referenced contrast (pm)",xlims=(0,1.03maximum(clearance,rows)),
        ylims=(lo-pad,hi+pad),legend=iso==1 ? :bottomright : false,
        left_margin=12Plots.mm,bottom_margin=12Plots.mm,right_margin=5Plots.mm)
    if !isempty(rr)
        scatter!(p,clearance.(rr),1000n.(rr,"referenced_contrast_nm");markersize=2.5,
            markerstrokewidth=0,alpha=.7,label="each valid reference site")
        hline!(p,[1000n(first(rr),"original_molecular_contrast_nm")];color=:black,
            linestyle=:dash,linewidth=1.5,label="original molecular mean contrast")
    end
    hline!(p,[0.];color=:gray,linestyle=:dot,label="zero contrast")
    push!(panels,p)
end
savefig(plot(panels...;layout=(1,3),size=(1650,620),
    plot_title="Reference-site dependence, Cu-half domain only — no adopted calibration; all domains retained in TSVs"),out)
