module VacuumPotential
using Statistics

function plane_stats(total,electrostatic)
    size(total)==size(electrostatic) && !isempty(total) || error("Different/empty planes")
    all(isfinite,total) && all(isfinite,electrostatic) || error("Nonfinite potentials")
    xc=total.-electrostatic
    (;samples=length(total),mean_total_ry=mean(total),min_total_ry=minimum(total),
        max_total_ry=maximum(total),std_total_ry=std(total;corrected=false),
        mean_electrostatic_ry=mean(electrostatic),min_electrostatic_ry=minimum(electrostatic),
        max_electrostatic_ry=maximum(electrostatic),std_electrostatic_ry=std(electrostatic;corrected=false),
        mean_xc_ry=mean(xc),min_xc_ry=minimum(xc),max_xc_ry=maximum(xc),std_xc_ry=std(xc;corrected=false))
end

"G_parallel=0 local-barrier descriptors, not a propagated state or error bound."
function barrier(stats,energy_ev,ry_ev,bohr_nm)
    lo=stats.min_total_ry-energy_ev/ry_ev; hi=stats.max_total_ry-energy_ev/ry_ev
    avg=stats.mean_total_ry-energy_ev/ry_ev
    allpositive=lo>0
    (;barrier_min_ev=lo*ry_ev,barrier_mean_ev=avg*ry_ev,barrier_max_ev=hi*ry_ev,
        all_sampled_barriers_positive=allpositive,
        lateral_span_over_mean_barrier=allpositive ? (hi-lo)/avg : NaN,
        kappa_min_nm_inv=allpositive ? sqrt(lo)/bohr_nm : NaN,
        kappa_mean_potential_nm_inv=avg>0 ? sqrt(avg)/bohr_nm : NaN,
        kappa_max_nm_inv=hi>0 ? sqrt(hi)/bohr_nm : NaN)
end
end
