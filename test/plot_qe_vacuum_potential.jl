#!/usr/bin/env julia
# Plots all planes in the predeclared whole-plane PAW-free gap, no fitting.
using Plots
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow
length(ARGS)==2 || error("plot_qe_vacuum_potential.jl SAVED_RUN NEW_PNG")
run,out=ARGS; ispath(out) && error("Output exists")
n(r,k)=parse(Float64,r[k]); b(r,k)=parse(Bool,r[k])
ry_ev=13.605693122994 # QE 7.4.1 energy conversion, not a fitted scale.
panels=[]; geometry=V.tsv(joinpath(run,"geometry.tsv"))
for molecule in ("glcn","glcnac")
    rows=filter(r->b(r,"paw_gap"),V.tsv(joinpath(run,molecule,"analysis/planes.tsv")))
    isempty(rows) && error("No PAW-free planes; no alternate domain")
    g=only(unique(n(r,"midpoint_nm") for r in geometry if r["molecule"]==molecule))
    z=n.(rows,"z_nm"); mu=n.(rows,"mean_total_minus_fermi_ev")
    lo=ry_ev.*(n.(rows,"mean_total_ry").-n.(rows,"min_total_ry"))
    hi=ry_ev.*(n.(rows,"max_total_ry").-n.(rows,"mean_total_ry"))
    p=plot(z,mu;ribbon=(lo,hi),fillalpha=.18,color=:royalblue,label="total: mean and full lateral range",
        title="$molecule — all PAW-free planes",xlabel="accepted-cell z (nm)",ylabel="local potential minus Fermi level (eV)",
        legend=:bottomright,left_margin=10Plots.mm,bottom_margin=8Plots.mm)
    plot!(p,z,mu.+ry_ev.*(n.(rows,"mean_electrostatic_ry").-n.(rows,"mean_total_ry"));
        color=:darkorange,linestyle=:dash,label="electrostatic mean (not total)")
    vline!(p,[g];color=:gray,linestyle=:dot,label="geometric gap midpoint")
    push!(panels,p)
    p=plot(z,ry_ev.*(n.(rows,"max_total_ry").-n.(rows,"min_total_ry"));color=:royalblue,
        label="full lateral range",title="$molecule — lateral variation",xlabel="accepted-cell z (nm)",
        ylabel="total-potential variation (eV)",legend=:topright,
        left_margin=10Plots.mm,bottom_margin=8Plots.mm)
    plot!(p,z,ry_ev.*n.(rows,"std_total_ry");color=:darkorange,label="population spatial SD")
    vline!(p,[g];color=:gray,linestyle=:dot,label="geometric gap midpoint")
    push!(panels,p)
end
savefig(plot(panels[1],panels[3],panels[2],panels[4];layout=(2,2),size=(1350,1150),
    plot_title="Accepted-state vacuum potential — no matching plane or flatness threshold selected"),out)
