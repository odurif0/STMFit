#!/usr/bin/env julia
# Descriptive geometry, or every paired Cu-half root. No fitted trend/cutoff.
using Plots
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
include(joinpath(@__DIR__,"lib/molecular_exterior.jl"))
const V=VerifyQESpectralWindow
const E=MolecularExterior
length(ARGS) in (2,3) || error("plot_molecular_exterior.jl GEOMETRY NEW_PNG [REPORT]")
prepared,out=ARGS[1:2]; ispath(out) && error("Output already exists")
n(r,k)=parse(Float64,r[k]); i(r,k)=parse(Int,r[k]); b(r,k)=parse(Bool,r[k])
grid=V.tsv(joinpath(prepared,"membership.tsv")); panels=[]
if length(ARGS)==2
    hulls=V.tsv(joinpath(prepared,"hulls.tsv")); cell=first(hulls)
    for molecule in ("glcn","glcnac")
        p=plot(;title="$molecule: same common support",xlabel="x (nm)",ylabel="y (nm)",
            xlims=(0,n(cell,"cell_x_nm")),ylims=(0,n(cell,"cell_y_nm")),aspect_ratio=:equal,
            legend=:outerbottom,bottom_margin=6Plots.mm,left_margin=8Plots.mm)
        hh=filter(r->r["molecule"]==molecule,hulls); radius=n(first(hh),"radius_nm")
        # Display-only polygon approximation; exact predicates set membership.
        outline=E.hull([(n(r,"x_nm")+radius*cos(t),n(r,"y_nm")+radius*sin(t))
            for r in hh for t in range(0,2pi;length=129)[1:end-1]])
        plot!(p,Shape(Float64.(first.(outline)),Float64.(last.(outline)));
            fillalpha=.12,color=:black,linewidth=1,label="whole molecular envelope")
        for (label,keep,color) in (("outside old support",r->!b(r,"prior_eligible"),:gray75),
            ("removed by envelope",r->b(r,"removed_by_envelope"),:darkorange),
            ("common exterior",r->b(r,"common_exterior"),:royalblue))
            rr=filter(keep,grid)
            scatter!(p,n.(rr,"x_nm"),n.(rr,"y_nm");color,markersize=2,markerstrokewidth=0,label="$label ($(length(rr)))")
        end
        scatter!(p,n.(hh,"x_nm"),n.(hh,"y_nm");color=:black,markersize=3,label="hull vertices")
        push!(panels,p)
    end
    savefig(plot(panels...;layout=(1,2),size=(1300,800),
        plot_title="Geometry-only common exterior; rounded contours for display, exact membership"),out)
else
    allrows=V.tsv(joinpath(ARGS[3],"paired.tsv"))
    rows=filter(r->r["domain"]=="cu_half_gap" && r["group"]=="all_prior",allrows)
    gg=Dict(i(r,"pixel")=>r for r in grid)
    margin(r)=min(n(gg[i(r,"pixel")],"glcn_margin_nm"),n(gg[i(r,"pixel")],"glcnac_margin_nm"))
    ys=vcat(1000n.(rows,"referenced_contrast_nm"),1000n.(rows,"original_molecular_contrast_nm"),[0.])
    lo,hi=extrema(ys); pad=.05*(hi-lo); xmin,xmax=extrema(margin.(rows)); xpad=.04*(xmax-xmin)
    for interval in 1:3
        p=plot(;title="isovalue $interval",xlabel="minimum envelope clearance (nm)",
            ylabel="hypothetical referenced contrast (pm)",ylims=(lo-pad,hi+pad),xlims=(xmin-xpad,xmax+xpad),
            legend=interval==1 ? :bottomright : false,left_margin=12Plots.mm,bottom_margin=12Plots.mm,right_margin=5Plots.mm)
        for (keep,label,color) in ((true,"common exterior",:royalblue),(false,"removed complement",:darkorange))
            rr=filter(r->i(r,"interval")==interval && b(gg[i(r,"pixel")],"common_exterior")==keep,rows)
            scatter!(p,margin.(rr),1000n.(rr,"referenced_contrast_nm");color,markersize=2.5,alpha=.65,markerstrokewidth=0,label)
        end
        hline!(p,[0.];color=:gray,linestyle=:dot,label="zero contrast")
        rr=filter(r->i(r,"interval")==interval,rows)
        hline!(p,[1000n(first(rr),"original_molecular_contrast_nm")];color=:black,linestyle=:dash,label="original molecular mean contrast")
        vline!(p,[0.];color=:gray,linewidth=.7,label="")
        push!(panels,p)
    end
    savefig(plot(panels...;layout=(1,3),size=(1650,620),
        plot_title="All paired Cu-half sites retained — geometric exterior is not a calibrated background"),out)
end
