#!/usr/bin/env julia
# Shows all geometry-tagged PAW-free planes; never selects a matching plane.
using Plots
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow
length(ARGS)==2 || error("plot_qe_xc_components.jl SAVED_RUN NEW_PNG")
run,out=ARGS; ispath(out) && error("Output exists")
number(r,k)=parse(Float64,r[k]); yes(r,k)=parse(Bool,r[k])
panels=[]
components=("xc_total","lda","gga_local","gga_divergence")
labels=("total XC","local-density part","local GGA term","negative FFT divergence")
colors=(:black,:forestgreen,:darkorange,:royalblue)
for molecule in ("glcn","glcnac")
    rows=filter(r->yes(r,"paw_gap") && r["subset"]=="all",V.tsv(joinpath(run,molecule,"analysis/components.tsv")))
    for (stat,title) in (("mean","plane means"),("spatial_sd","lateral spatial SD"))
        p=plot(;title="$molecule — $title",xlabel="accepted-cell z (nm)",ylabel="XC potential (eV)",
            legend=:outerright,left_margin=8Plots.mm,bottom_margin=6Plots.mm)
        for (component,label,color) in zip(components,labels,colors)
            selected=filter(r->r["component"]==component,rows)
            isempty(selected) && error("Missing component")
            plot!(p,number.(selected,"z_nm"),number.(selected,stat).*13.605693122994;
                label,color,linewidth=2)
        end
        push!(panels,p)
    end
    counts=filter(r->yes(r,"paw_gap"),V.tsv(joinpath(run,molecule,"analysis/counts.tsv")))
    z=number.(counts,"z_nm"); ns=number.(counts,"samples")
    p=plot(z,1 .-number.(counts,"gga_active")./ns;label="GGA inactive (source rule)",color=:royalblue,
        title="$molecule — density and source branches",xlabel="accepted-cell z (nm)",
        ylabel="fraction of full lateral plane",legend=:outerright,
        left_margin=8Plots.mm,bottom_margin=6Plots.mm,ylims=(0,1))
    plot!(p,z,number.(counts,"negative_valence")./ns;label="negative valence density",color=:darkorange)
    plot!(p,z,number.(counts,"negative_xc_input")./ns;label="negative XC input density",color=:forestgreen)
    push!(panels,p)
end
savefig(plot(panels[1],panels[4],panels[2],panels[5],panels[3],panels[6];layout=(3,2),size=(1750,1350),
    plot_title="Frozen-state XC attribution — all PAW-free planes, no potential replacement"),out)
