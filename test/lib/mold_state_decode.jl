module MoldStateDecode
module Connected
include(joinpath(@__DIR__, "..", "score_connected_mold_templates.jl"))
end
include(joinpath(@__DIR__, "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment: read_table, lobe_table, require_same_keys
using TOML

function settings(path)
    c = TOML.parsefile(path)
    c == Dict("model" => Dict("state_objective" => "sum_of_per_view_unary_minima",
            "chemical_states" => "independent_per_view", "view_weights" => [1.0, 1.0]),
        "selection" => Dict("missing_cost" => "omit_both_infinite_keep_unavailable",
            "state_tie" => "first_phase_then_mirror", "transition_penalty" => 0.0),
        "preprocessing" => Dict("cost_source" => "saved_target_gaussian_affine_audit")) ||
        error("Only the fixed, label-free mold-state comparison settings are supported")
    return (weights=Float64.(c["model"]["view_weights"]),
            transition_penalty=Float64(c["selection"]["transition_penalty"]))
end

"Load full-precision state costs, never rounded selected margins or predictions."
function load_costs(paths, base)
    costs = Dict{Tuple{String,Int,Int,Int},Vector{Float64}}()
    for path in paths
        _, rows = read_table(path)
        for r in rows
            f = basename(r["file"]); l = parse(Int, r["lobe"])
            p = parse(Int, r["parity"]); m = parse(Int, r["mirror"])
            (f,l) in keys(base) && p in 0:1 && m in 0:1 || error("Unexpected state key")
            k = (f,l,p,m); haskey(costs,k) && error("Duplicate state cost")
            v = parse.(Float64, [r["cost0"],r["cost1"]])
            if r["reason"] == "ok"
                all(isfinite,v) || error("Nonfinite available state")
            elseif r["reason"] == "insufficient_native_support"
                all(==(Inf),v) || error("Missing state costs must remain infinite")
            else
                error("Unsupported saved score failure: $(r["reason"])")
            end
            costs[k] = v
        end
    end
    length(costs) == 4length(base) || error("Incomplete state-cost cohort")
    for (f,l) in keys(base)
        available = [all(isfinite,costs[(f,l,p,m)]) for p in 0:1 for m in 0:1]
        all(==(first(available)),available) || error("State-dependent missing support")
    end
    return costs
end

function records_for(base, file)
    keys_for_file = sort([k for k in keys(base) if k[1] == file])
    n = length(keys_for_file)
    last.(keys_for_file) == collect(1:n) || error("Noncontiguous records")
    all(parse(Int,base[k]["N"]) == n for k in keys_for_file) || error("N/key mismatch")
    return [(file=file,lobe=k[2],amplitude=parse(Float64,base[k]["amplitude"]),patch=Float64[]) for k in keys_for_file]
end

function unary(costs, records, phase, mirror)
    reduce(vcat, permutedims(costs[(r.file,r.lobe,mod(r.lobe-1+phase,2),mirror)]) for r in records)
end

function decode(records, views, opt; mode)
    mode in (:reference,:finite,:shared) || error("Unknown decoder ablation")
    length(views) == 2 && length(opt.weights) == 2 || error("Two views required")
    legacy = Connected.Options("","",nothing,"","","ncc","contrast",opt.transition_penalty,1.0)
    if mode != :shared
        # Preserve the existing eight-state enumeration, arithmetic and tie order.
        templates = Dict((t,p,m)=>[Float64(t),Float64(p),Float64(m)] for t in 0:1 for p in 0:1 for m in 0:1)
        return [Connected._decode_file(records,templates,nothing,legacy;
            scorer=(r,t)->costs[(r.file,r.lobe,Int(t[2]),Int(t[3]))][Int(t[1])+1],
            omit_unavailable=mode==:finite) for costs in views]
    end
    # With unary-only scores, reversed direction duplicates phase according to
    # phase_eff = mod(phase + direction*(N-1), 2). Enumerate four distinct states.
    best = nothing; best_total = Inf
    for phase in 0:1, mirror in 0:1
        candidates = [Connected._decode_state(unary(v,records,phase,mirror),records,0,phase,mirror,
            nothing,legacy;omit_unavailable=true) for v in views]
        total = sum(opt.weights[j]*candidates[j].total for j in 1:2)
        if best === nothing || total < best_total
            best = candidates; best_total = total
        end
    end
    return best
end

function audit_paths(dir, view)
    paths = sort([joinpath(dir,f) for f in readdir(dir) if
        occursin(Regex("^score_"*view*"(?:_chunk[1-9][0-9]*)?\\.tsv\\.audit\\.tsv\$"),f)])
    isempty(paths) && error("Missing $view state-cost audits")
    return paths
end
end
