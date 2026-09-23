module MoldLeaveOneOut
using TOML, Printf

function settings(path)
    c=TOML.parsefile(path)
    c==Dict("model"=>Dict("state_objective"=>"leave_target_out_unary_minima",
            "chemical_states"=>"independent_per_view","view_weights"=>[1.0,1.0]),
        "selection"=>Dict("missing_cost"=>"omit_both_infinite_keep_unavailable",
            "state_tie"=>"first_phase_then_mirror","transition_penalty"=>0.0,
            "empty_training"=>"unavailable"),
        "preprocessing"=>Dict("cost_source"=>"saved_target_gaussian_affine_audit")) ||
        error("Only fixed leave-target-out mold settings are supported")
    # Same interface as the independent finite control. Views are never pooled;
    # unit weights are retained explicitly, not fitted using chemical outcomes.
    return (weights=Float64.(c["model"]["view_weights"]),transition_penalty=0.0)
end

"Exclude the target before accumulation, not by subtracting from a full sum."
function decode(records, costs, opt)
    opt.transition_penalty==0 || error("Unary-only state selection required")
    n=length(records); n>0 || error("Empty chain")
    getproperty.(records,:lobe)==collect(1:n) || error("Noncontiguous lobes")
    length(unique(getproperty.(records,:file)))==1 || error("Mixed chains")
    states=[(p,m) for p in 0:1 for m in 0:1]
    values=Array{Float64}(undef,n,2,4)
    available=falses(n)
    for (s,(p,m)) in enumerate(states), (i,r) in enumerate(records)
        c=costs[(r.file,r.lobe,mod(r.lobe-1+p,2),m)]
        length(c)==2 && (all(isfinite,c) || all(==(Inf),c)) || error("Invalid unary pair")
        values[i,:,s]=c
        observed=all(isfinite,c)
        if s==1
            available[i]=observed
        else
            observed==available[i] || error("State-dependent missing support")
        end
    end
    output=NamedTuple[]
    for i in 1:n
        train=[j for j in 1:n if j!=i && available[j]]
        objectives=[sum((min(values[j,1,s],values[j,2,s]) for j in train);init=0.0) for s in 1:4]
        if isempty(train)
            push!(output,(phase=-1,mirror=-1,training_cost=0.0,training_lobes=0,
                costs=[Inf,Inf],label=-1,reason="no_remaining_observations",objectives=objectives))
            continue
        end
        s=argmin(objectives); p,m=states[s]; c=values[i,:,s]
        label=available[i] ? Int(c[2]<c[1]) : -1
        push!(output,(phase=p,mirror=m,training_cost=objectives[s],training_lobes=length(train),
            costs=c,label=label,reason=available[i] ? "ok" : "insufficient_native_support",objectives=objectives))
    end
    return output
end

# Per-target states are not falsely called a global chain reconstruction.
const SCORE_HEADER=["file","lobe","predicted","amplitude","physical_label","cost_GlcN","cost_GlcNAc",
    "cost_margin","state_phase","state_mirror","training_cost","training_lobes","state_reason"]
function write_decoded(io, records, decoded)
    length(records)==length(decoded) || error("Decoded row count differs")
    for (r,b) in zip(records,decoded)
        c0,c1=b.costs; physical=b.label<0 ? "?" : b.label==1 ? "GlcNAc" : "GlcN"
        println(io,join([r.file,r.lobe,b.label,@sprintf("%.8g",r.amplitude),physical,
            @sprintf("%.8g",c0),@sprintf("%.8g",c1),@sprintf("%.8g",abs(c0-c1)),
            b.phase,b.mirror,@sprintf("%.17g",b.training_cost),b.training_lobes,b.reason],'\t'))
    end
end
end
