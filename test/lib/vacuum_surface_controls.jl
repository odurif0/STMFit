module VacuumSurfaceControls
using LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__, "vacuum_surface_transfer.jl"))
const F = VacuumSurfaceTransfer

function settings(path)
    c = TOML.parsefile(path)
    fields = Dict(
        "model" => ["input_run", "input_settings_sha256", "input_files_sha256", "arms",
            "contrast_mean", "height_gain", "plane_basis", "chemical_geometry"],
        "selection" => ["policy", "ties", "rank_failure"],
        "preprocessing" => ["observations", "target_use", "plane_rank_rtol", "replay_rtol", "replay_atol_nm", "replay_atol_nm2"])
    Set(keys(c)) == Set(keys(fields)) || error("Unexpected sections")
    for (section, keys_expected) in fields
        Set(keys(c[section])) == Set(keys_expected) || error("Unknown/missing $section setting")
    end
    m, s, p = c["model"], c["selection"], c["preprocessing"]
    m["arms"] == ["height_contrast", "local_constant", "local_plane", "plane_common", "plane_chemical"] || error("All controls required")
    m["contrast_mean"] == "equal_weight_all_fixed_map_nodes" && m["height_gain"] == 1 || error("Changed contrast convention")
    m["plane_basis"] == "1_dx_nm_dy_nm" && m["chemical_geometry"] == "frozen_source_plane_common_argmin" || error("Changed nuisance/geometry")
    s["policy"] == "all_cases_directions_isovalues_no_grade" &&
        s["ties"] == "first_exact_minimum_no_chemical_confidence" &&
        s["rank_failure"] == "retain_failed_view_no_support_change" || error("Changed selection policy")
    p["observations"] == "saved_native_pixels_no_refiltering" && p["target_use"] == "score_only_never_fit" || error("Changed observation boundary")
    for key in ("plane_rank_rtol", "replay_rtol", "replay_atol_nm", "replay_atol_nm2")
        v = p[key]
        v isa Real && !(v isa Bool) && isfinite(v) && 0 < v < 1 || error("Invalid numerical tolerance")
    end
    for key in ("input_settings_sha256", "input_files_sha256")
        occursin(r"^[0-9a-f]{64}$", m[key]) || error("Invalid input hash")
    end
    c
end

"Projection basis uses source coordinates only; no target or chemical input."
function plane_basis(p, rank_rtol)
    X = hcat(ones(length(p.values)), p.dx, p.dy)
    all(isfinite, X) && all(isfinite, p.values) || error("Nonfinite source observation")
    fact = svd(X; full=false)
    length(fact.S) == 3 && minimum(fact.S) > rank_rtol * maximum(fact.S) || error("Unidentifiable source plane")
    (; X, U=fact.U, S=fact.S, V=fact.V)
end

function remove_plane(residual, basis)
    q = basis.U' * residual
    beta = basis.V * (q ./ basis.S)
    remainder = residual - basis.U * q
    (; beta, loss=sum(abs2, remainder))
end

"Mean contrast is fixed by the map grid, not by source/target pixels or labels."
mean_contrast(maps) = mean((maps[1] - maps[0]) / 2)

function fit_patches(pp, maps, base, c)
    isempty(pp) && error("No source patches")
    ss = F.states(base, "common")
    hs = F.states(base, "chemical")
    delta = mean_contrast(maps)
    bases = [plane_basis(p, c["preprocessing"]["plane_rank_rtol"]) for p in pp]
    common_costs = NamedTuple[]; plane_costs = NamedTuple[]; selected = Int[]
    constants = Float64[]; flatplanes = Vector{Float64}[]
    chemical_costs = NamedTuple[]; chemical_choices = Int[]
    for (p, basis) in zip(pp, bases)
        mu = Float64[]; v = Float64[]; losses = Float64[]; betas = Vector{Float64}[]
        for state in ss
            residual = p.values - F.prediction(p, maps[-1], state, base)
            center = mean(residual)
            push!(mu, center); push!(v, sum(abs2, residual .- center))
            plane = remove_plane(residual, basis)
            push!(losses, plane.loss); push!(betas, plane.beta)
        end
        push!(common_costs, (; n=length(p.values), mu, v))
        push!(plane_costs, (; losses, betas))
        k = argmin(losses); push!(selected, k)
        # Chemistry cannot refit rotation/translation or choose another anchor.
        state = ss[k]; chem_loss = Float64[]; chem_beta = Vector{Float64}[]
        for type in (0, 1)
            residual = p.values - F.prediction(p, maps[type], state, base)
            plane = remove_plane(residual, basis)
            push!(chem_loss, plane.loss); push!(chem_beta, plane.beta)
        end
        push!(chemical_costs, (; losses=chem_loss, betas=chem_beta))
        push!(chemical_choices, argmin(chem_loss))
        push!(constants, mean(p.values)); push!(flatplanes, remove_plane(p.values, basis).beta)
    end
    # Subtracting +/- delta changes residual means, never centered costs.
    height_costs = [(; n=x.n, mu=vcat(x.mu .+ delta, x.mu .- delta), v=vcat(x.v, x.v)) for x in common_costs]
    replay = F.common_offset(common_costs)
    height = F.common_offset(height_costs)
    fits = Dict{String,Vector{NamedTuple}}(arm => NamedTuple[] for arm in c["model"]["arms"])
    for (j, p) in enumerate(pp)
        basis = bases[j]; state = ss[selected[j]]
        k = height.choices[j]; hs0 = hs[k]
        hpred = height.offset .+ (2hs0.type-1)*delta .+ F.prediction(p, maps[-1], hs0, base)
        push!(fits["height_contrast"], (; state=hs0, choice=k, beta=[height.offset, 0., 0.], pred=hpred))
        dummy = (; type=-1, angle=0., tx=0., ty=0.)
        push!(fits["local_constant"], (; state=dummy, choice=0, beta=[constants[j], 0., 0.], pred=fill(constants[j], length(p.values))))
        push!(fits["local_plane"], (; state=dummy, choice=0, beta=flatplanes[j], pred=basis.X*flatplanes[j]))
        beta = plane_costs[j].betas[selected[j]]
        pred = F.prediction(p, maps[-1], state, base) + basis.X*beta
        push!(fits["plane_common"], (; state, choice=selected[j], beta, pred))
        ch = chemical_choices[j]; beta = chemical_costs[j].betas[ch]
        chem_state = merge(state, (; type=ch-1))
        pred = F.prediction(p, maps[ch-1], chem_state, base) + basis.X*beta
        push!(fits["plane_chemical"], (; state=chem_state, choice=ch, beta, pred))
    end
    (; fits, ss, hs, delta, common_costs, plane_costs, chemical_costs, height_costs, height, replay)
end
end
