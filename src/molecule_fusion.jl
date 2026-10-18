# Label-free fusion of per-scan unit calls over repeated scans of one molecule.
#
# Each physical lobe of a molecule track (scans linked by molecule consensus,
# same final count) has a hidden type z. Every scan's call is a noisy binary
# observation: P(call 1 | z=1) = theta1 (detection rate), P(call 1 | z=0) =
# theta0 (false-call rate), P(z=1) = pi. All three are estimated by EM over
# all physical lobes (binomial latent-class model, one repeated annotator).
# Unlike a majority vote, a type detected in a minority of scans is kept when
# false calls are rare. No label, expected sequence or class count is read;
# pi is a fitted mixing weight, and the class-1 component is the one with the
# higher call rate (theta1 > theta0).
module MoleculeFusion

using Statistics

export physical_lobes, fit_latent_class, FusionSettings, load_fusion_settings

struct FusionSettings
    min_scans::Int
    max_iter::Int
    tol::Float64
    init_pi::Float64
    init_theta0::Float64
    init_theta1::Float64
end

function load_fusion_settings(cfg::AbstractDict)
    s = cfg["selection"]; m = cfg["model"]
    get(s, "fusion", "") == "latent_class" || error("Only latent_class molecule fusion is defined")
    out = FusionSettings(Int(s["fusion_min_scans"]), Int(m["fusion_em_max_iter"]), Float64(m["fusion_em_tol"]),
        Float64(m["fusion_init_pi"]), Float64(m["fusion_init_theta0"]), Float64(m["fusion_init_theta1"]))
    out.min_scans >= 2 && out.max_iter > 0 && out.tol > 0 && 0 < out.init_pi < 1 &&
        0 < out.init_theta0 < out.init_theta1 < 1 || error("Invalid molecule-fusion settings")
    return out
end

"""
    physical_lobes(order, tracks, counts, lobes_abs)

`order`: acquisition order; `tracks`: file => track id; `counts`: file => final
lobe count; `lobes_abs`: file => absolute lobe positions (lobe order). For each
track, scans with the track's modal count (ties: earliest in acquisition order)
are ranked along the chain axis of the first such scan. Ranking along a fixed
absolute direction is invariant to drift. Returns (file, lobe) => key.
"""
function physical_lobes(order, tracks, counts, lobes_abs)
    keys_out = Dict{Tuple{String,Int},String}()
    members = Dict{Int,Vector{String}}()
    for f in order; push!(get!(members, tracks[f], String[]), f); end
    for (t, fs) in members
        tally = Dict{Int,Int}(); for f in fs; tally[counts[f]] = get(tally, counts[f], 0) + 1; end
        top = maximum(values(tally))
        nmode = first(counts[f] for f in fs if tally[counts[f]] == top)
        ref = first(f for f in fs if counts[f] == nmode)
        P = lobes_abs[ref]; length(P) >= 2 || continue
        ax = (P[end][1] - P[1][1], P[end][2] - P[1][2]); L = hypot(ax...); L > 0 || continue
        ax = (ax[1] / L, ax[2] / L)
        for f in fs
            counts[f] == nmode || continue
            proj = [p[1] * ax[1] + p[2] * ax[2] for p in lobes_abs[f]]
            rank = invperm(sortperm(proj))
            for (lobe, r) in enumerate(rank)
                keys_out[(f, lobe)] = "track$(t)_$(r)"
            end
        end
    end
    return keys_out
end

"""
    fit_latent_class(k, n, st) -> (posterior, pi, theta0, theta1, iterations)

`k[i]` of `n[i]` scans call physical lobe i class 1.
"""
function fit_latent_class(k::AbstractVector{<:Integer}, n::AbstractVector{<:Integer}, st::FusionSettings)
    length(k) == length(n) && !isempty(k) || error("Latent-class fit needs matching non-empty counts")
    all(0 .<= k .<= n) && all(n .>= 1) || error("Invalid call counts")
    kf = Float64.(k); nf = Float64.(n)
    pi, t0, t1 = st.init_pi, st.init_theta0, st.init_theta1
    w = fill(0.5, length(k)); it = 0
    clampp(x) = clamp(x, 1e-6, 1 - 1e-6)
    for i in 1:st.max_iter
        it = i
        l1 = log(pi) .+ kf .* log(t1) .+ (nf .- kf) .* log(1 - t1)
        l0 = log(1 - pi) .+ kf .* log(t0) .+ (nf .- kf) .* log(1 - t0)
        mx = max.(l0, l1)
        w = exp.(l1 .- mx) ./ (exp.(l1 .- mx) .+ exp.(l0 .- mx))
        npi = clampp(mean(w))
        nt1 = clampp(sum(w .* kf) / max(sum(w .* nf), eps()))
        nt0 = clampp(sum((1 .- w) .* kf) / max(sum((1 .- w) .* nf), eps()))
        delta = max(abs(npi - pi), abs(nt0 - t0), abs(nt1 - t1))
        pi, t0, t1 = npi, nt0, nt1
        delta < st.tol && break
    end
    t1 > t0 || error("Latent-class fit did not separate call rates")
    return (posterior=w, pi=pi, theta0=t0, theta1=t1, iterations=it)
end

end
