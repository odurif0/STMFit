# Label-free reliability of per-unit signals across repeated scans of one molecule.
#
# A run of run_molecule_consensus_chitosan.jl gives, for every scan, its lobes
# (scan_geometry.tsv), the physical lobe each one belongs to (predictions.tsv)
# and the molecule track of the scan (consensus/consensus.tsv). Scans of one
# track show the same molecule, so a physical lobe is measured once per scan.
#
# Signals are residuals: the log of a lobe feature minus its scan median (scan
# gain) minus the mean at the lobe's distance class d to the nearest chain end
# (end effect). The reliability of a residual is its repeatability over the
# scans of a track: a one-way random-effects model over physical lobes gives
# the between-lobe variance (what repeats), the within-lobe variance (what does
# not) and ICC = between / (between + within). Intervals come from a bootstrap
# over tracks, since lobes of one molecule are not independent.
#
# Nothing here reads a label, an expected count or a sequence. The reading of
# the numbers (the 6-mer as a fixed-sequence negative control) lives in the
# journal and in scripts/data/report_unit_reliability.jl.
module UnitReliability

using Statistics, Random, Printf, Dates
using ..ReconstructedUnitAssignment: read_table

export SIGNALS, Lobe, RunUnits, load_run_units, signal_value, distance_profile, track_groups,
       Moments, moments, variance_components, bootstrap_components, ratio_interval,
       pair_correlations, lag_summary, persistent_profiles, profile_autocorrelation,
       autocorrelation_null, count_dispersion, capped_files, centred, normal_cdf

"Signals a lobe offers (log scale): amplitude, integrated amplitude, widths."
const SIGNALS = (:amplitude, :integrated, :width_along, :width_across)

struct Lobe
    k::Int               # lobe index in its scan
    d::Int               # distance to the nearest chain end (1 = terminal lobe)
    physical::String     # physical lobe of the track
    amplitude::Float64
    sigma_par::Float64
    sigma_perp::Float64
end

struct RunUnits
    name::String
    tracks::Dict{String,Vector{String}}   # track => scan files, in acquisition order
    track_of::Dict{String,String}
    lobes::Dict{String,Vector{Lobe}}      # scan file => lobes ordered by k
    n_scan::Dict{String,Int}              # per-scan count of the counting stage
    n_final::Dict{String,Int}             # count after the molecule consensus
end

"Log of one lobe feature (NaN when the feature is not positive)."
function signal_value(l::Lobe, kind::Symbol)
    v = kind == :amplitude ? l.amplitude :
        kind == :integrated ? l.amplitude * l.sigma_par * l.sigma_perp :
        kind == :width_along ? l.sigma_par :
        kind == :width_across ? l.sigma_perp : error("Unknown signal $kind")
    return v > 0 ? log(v) : NaN
end

"""
    load_run_units(run_dir; name, times)

Read a run of the promoted method. `times` (file => any sortable value, for
example a DateTime) orders the scans of a track by acquisition (files without
a time keep the file-name order); without `times` the file names are used,
which order the scans of one session.
"""
function load_run_units(run_dir::AbstractString; name::AbstractString=basename(abspath(run_dir)), times=nothing)
    _, cons = read_table(joinpath(run_dir, "consensus", "consensus.tsv"))
    _, geo = read_table(joinpath(run_dir, "scan_geometry.tsv"))
    _, pred = read_table(joinpath(run_dir, "predictions.tsv"))
    physical = Dict((r["file"], parse(Int, r["lobe"])) => r["physical_lobe"] for r in pred)
    track_of = Dict(r["file"] => r["track"] for r in cons)
    n_scan = Dict(r["file"] => parse(Int, r["N_scan"]) for r in cons)
    n_final = Dict(r["file"] => parse(Int, r["N_final"]) for r in cons)
    rows = Dict{String,Vector{Dict{String,String}}}()
    for r in geo
        push!(get!(rows, r["file"], Dict{String,String}[]), r)
    end
    lobes = Dict{String,Vector{Lobe}}()
    for (f, rs) in rows
        sort!(rs; by=r -> parse(Int, r["lobe"]))
        n = length(rs)
        lobes[f] = [Lobe(parse(Int, r["lobe"]), min(parse(Int, r["lobe"]), n + 1 - parse(Int, r["lobe"])),
                         get(physical, (f, parse(Int, r["lobe"])), ""),
                         parse(Float64, r["amplitude"]), parse(Float64, r["sigma_parallel_nm"]),
                         parse(Float64, r["sigma_perp_nm"])) for r in rs]
    end
    tracks = Dict{String,Vector{String}}()
    for f in keys(track_of)
        haskey(lobes, f) && push!(get!(tracks, track_of[f], String[]), f)
    end
    order(f) = times === nothing ? (typemin(DateTime), f) : (get(times, f, typemin(DateTime)), f)
    for fs in values(tracks)
        sort!(fs; by=order)
    end
    return RunUnits(name, tracks, track_of, lobes, n_scan, n_final)
end

"At most `cap` files, evenly spaced over the list (first and last kept)."
function capped_files(files::AbstractVector, cap::Integer)
    n = length(files)
    n <= cap && return collect(files)
    return [files[round(Int, 1 + (i - 1) * (n - 1) / (cap - 1))] for i in 1:cap]
end

"Scan-centred signal of every lobe of a scan (value minus the scan median)."
function centred(lobes::AbstractVector{Lobe}, kind::Symbol)
    v = [signal_value(l, kind) for l in lobes]
    ok = filter(!isnan, v)
    isempty(ok) && return fill(NaN, length(v))
    return v .- median(ok)
end

"""
    normal_cdf(x)

Standard normal distribution function (Abramowitz and Stegun 7.1.26, absolute
error below 1.5e-7). The reading ceiling of a threshold between two classes
`d` standard deviations apart is `normal_cdf(d / 2)`.
"""
function normal_cdf(x::Real)
    z = abs(x) / sqrt(2)
    t = 1 / (1 + 0.3275911 * z)
    erf = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * exp(-z^2)
    return x >= 0 ? 0.5 * (1 + erf) : 0.5 * (1 - erf)
end

"""
    distance_profile(ru, kind; dcap=4)

Mean scan-centred signal at distance classes d = 1, ..., dcap (dcap stands for
d ≥ dcap): each track gives the mean of all its lobes of a class, and the
profile is the mean over the tracks that have the class (every track counts
once, so a molecule scanned 43 times does not dominate).
"""
function distance_profile(ru::RunUnits, kind::Symbol; dcap::Integer=4)
    sums = zeros(dcap); cnt = zeros(Int, dcap)
    for fs in values(ru.tracks)
        s = zeros(dcap); c = zeros(Int, dcap)
        for f in fs, (l, v) in zip(ru.lobes[f], centred(ru.lobes[f], kind))
            isnan(v) && continue
            j = min(l.d, dcap); s[j] += v; c[j] += 1
        end
        for j in 1:dcap
            c[j] > 0 && (sums[j] += s[j] / c[j]; cnt[j] += 1)
        end
    end
    return [cnt[j] > 0 ? sums[j] / cnt[j] : NaN for j in 1:dcap]
end

"""
    track_groups(ru, kind; min_scans=3, cap=8, dmin=2, profile, min_group=3)

Residuals grouped by physical lobe, for each track of at least `min_scans`
scans (the `cap` scans evenly spaced in acquisition order are used). Only lobes
with d ≥ `dmin` are kept, and only physical lobes seen in at least `min_group`
of the used scans. Returns track => Vector of groups (one residual per scan).
"""
function track_groups(ru::RunUnits, kind::Symbol; min_scans::Integer=3, cap::Integer=8, dmin::Integer=2,
                      profile::AbstractVector=distance_profile(ru, kind), min_group::Integer=3)
    out = Dict{String,Vector{Vector{Float64}}}()
    dcap = length(profile)
    for (t, fs) in ru.tracks
        length(fs) >= min_scans || continue
        groups = Dict{String,Vector{Float64}}()
        for f in capped_files(fs, cap)
            for (l, v) in zip(ru.lobes[f], centred(ru.lobes[f], kind))
                (isnan(v) || l.d < dmin || isempty(l.physical)) && continue
                push!(get!(groups, l.physical, Float64[]), v - profile[min(l.d, dcap)])
            end
        end
        g = [x for x in values(groups) if length(x) >= min_group]
        isempty(g) || (out[t] = g)
    end
    return out
end

"Sufficient statistics of a one-way design: groups, observations, Σ nᵢ², Σx, Σx², Σ (Σ xᵢ)²/nᵢ."
struct Moments
    groups::Int
    n::Int
    sum_n2::Float64
    sum_x::Float64
    sum_x2::Float64
    sum_group::Float64
end
Base.:+(a::Moments, b::Moments) = Moments(a.groups + b.groups, a.n + b.n, a.sum_n2 + b.sum_n2,
                                          a.sum_x + b.sum_x, a.sum_x2 + b.sum_x2, a.sum_group + b.sum_group)
Moments() = Moments(0, 0, 0.0, 0.0, 0.0, 0.0)

function moments(groups::AbstractVector{<:AbstractVector{<:Real}})
    m = Moments()
    for g in groups
        n = length(g); s = sum(g)
        m += Moments(1, n, Float64(n)^2, s, sum(abs2, g), s^2 / n)
    end
    return m
end

"""
    variance_components(m) -> (between, within, icc)

One-way random-effects ANOVA estimates for an unbalanced design: `between` and
`within` are variances (the between one can come out negative from sampling
noise), icc = between / (between + within), NaN when undefined.
"""
function variance_components(m::Moments)
    I, N = m.groups, m.n
    (I >= 2 && N > I) || return (NaN, NaN, NaN)
    ssb = m.sum_group - m.sum_x^2 / N
    ssw = m.sum_x2 - m.sum_group
    msb = ssb / (I - 1); msw = ssw / (N - I)
    n0 = (N - m.sum_n2 / N) / (I - 1)
    between = (msb - msw) / n0
    return (between, msw, between / (between + msw))
end

"""
    bootstrap_components(track_moments; draws=1000, rng)

Variance components of `draws` resamplings of the tracks (with replacement).
Returns (between, within, icc) vectors; undefined draws are dropped.
"""
function bootstrap_components(tm::AbstractVector{Moments}; draws::Integer=1000, rng::AbstractRNG=Random.default_rng())
    B = Float64[]; W = Float64[]; I = Float64[]
    n = length(tm)
    for _ in 1:draws
        m = Moments()
        for _ in 1:n
            m += tm[rand(rng, 1:n)]
        end
        b, w, i = variance_components(m)
        isnan(i) || (push!(B, b); push!(W, w); push!(I, i))
    end
    return B, W, I
end

"Percentile interval (default 95%) of a sample."
interval(x::AbstractVector; level=0.95) = isempty(x) ? (NaN, NaN) :
    (quantile(x, (1 - level) / 2), quantile(x, 1 - (1 - level) / 2))

"""
    ratio_interval(a, b; draws, rng, level)

Point estimate and percentile interval of between(a) / between(b), where `a`
and `b` are vectors of per-track Moments of two independent sets of tracks.
Draws with a non-positive denominator are dropped.
"""
function ratio_interval(a::AbstractVector{Moments}, b::AbstractVector{Moments}; draws::Integer=1000,
                        rng::AbstractRNG=Random.default_rng(), level=0.95)
    Ba, _, _ = bootstrap_components(a; draws, rng)
    Bb, _, _ = bootstrap_components(b; draws, rng)
    k = min(length(Ba), length(Bb))
    r = [Ba[i] / Bb[i] for i in 1:k if Bb[i] > 0]
    pa = variance_components(sum(a; init=Moments()))[1]; pb = variance_components(sum(b; init=Moments()))[1]
    return (pb > 0 ? pa / pb : NaN, interval(r; level)...)
end

"""
    pair_correlations(ru, kind, times; dmin=2, min_shared=4, profile)

For every pair of scans of a track (tracks of at least 3 scans), the Pearson
correlation of their residual profiles over the physical lobes they share, and
the time between the two scans in minutes. `times` maps a file to a DateTime.
Returns a vector of (track, lag_minutes, r).
"""
function pair_correlations(ru::RunUnits, kind::Symbol, times::AbstractDict; dmin::Integer=2, min_shared::Integer=4,
                           profile::AbstractVector=distance_profile(ru, kind))
    out = Tuple{String,Float64,Float64}[]
    dcap = length(profile)
    for (t, fs) in ru.tracks
        length(fs) >= 3 || continue
        vec = Dict{String,Dict{String,Float64}}()
        for f in fs
            d = Dict{String,Float64}()
            for (l, v) in zip(ru.lobes[f], centred(ru.lobes[f], kind))
                (isnan(v) || l.d < dmin || isempty(l.physical)) && continue
                d[l.physical] = v - profile[min(l.d, dcap)]
            end
            vec[f] = d
        end
        for i in 1:length(fs)-1, j in i+1:length(fs)
            a, b = vec[fs[i]], vec[fs[j]]
            shared = [p for p in keys(a) if haskey(b, p)]
            length(shared) >= min_shared || continue
            xa = [a[p] for p in shared]; xb = [b[p] for p in shared]
            (std(xa) > 0 && std(xb) > 0) || continue
            lag = abs(minutes_between(times[fs[i]], times[fs[j]]))
            push!(out, (t, lag, cor(xa, xb)))
        end
    end
    return out
end

"Minutes from timestamp b to timestamp a."
minutes_between(a::DateTime, b::DateTime) = Dates.value(a - b) / 60000

"""
    lag_summary(pairs, bins; draws, rng)

Mean correlation by lag bin (Fisher z within a track, then over tracks) and a
bootstrap interval over tracks. `bins` is a vector of (low, high) in minutes.
Returns a vector of (low, high, mean_r, lo, hi, tracks).
"""
function lag_summary(pairs, bins; draws::Integer=1000, rng::AbstractRNG=Random.default_rng())
    out = Tuple{Float64,Float64,Float64,Float64,Float64,Int}[]
    z(r) = atanh(clamp(r, -0.999, 0.999))
    for (lo, hi) in bins
        per = Dict{String,Vector{Float64}}()
        for (t, lag, r) in pairs
            lo <= lag < hi && push!(get!(per, t, Float64[]), z(r))
        end
        isempty(per) && continue
        v = [mean(x) for x in values(per)]
        boot = [tanh(mean(v[rand(rng, 1:length(v), length(v))])) for _ in 1:draws]
        l, h = interval(boot)
        push!(out, (lo, hi, tanh(mean(v)), l, h, length(v)))
    end
    return out
end

"""
    persistent_profiles(ru, kind; min_scans=3, dmin=2, min_units=4, profile)

For every track of at least `min_scans` scans, the mean residual of each
physical lobe (d ≥ `dmin`, seen in at least `min_scans` scans), ordered along
the chain. Tracks with fewer than `min_units` such lobes are dropped.
"""
function persistent_profiles(ru::RunUnits, kind::Symbol; min_scans::Integer=3, dmin::Integer=2, min_units::Integer=4,
                             profile::AbstractVector=distance_profile(ru, kind))
    out = Dict{String,Vector{Float64}}()
    dcap = length(profile)
    for (t, fs) in ru.tracks
        length(fs) >= min_scans || continue
        acc = Dict{String,Vector{Float64}}()
        for f in fs, (l, v) in zip(ru.lobes[f], centred(ru.lobes[f], kind))
            (isnan(v) || l.d < dmin || isempty(l.physical)) && continue
            push!(get!(acc, l.physical, Float64[]), v - profile[min(l.d, dcap)])
        end
        keys_ = sort([p for (p, x) in acc if length(x) >= min_scans]; by=p -> parse(Int, split(p, '_')[end]))
        length(keys_) >= min_units && (out[t] = [mean(acc[p]) for p in keys_])
    end
    return out
end

"Autocorrelation of a profile at a lag (sum of products over the sum of squares)."
function profile_autocorrelation(p::AbstractVector, lag::Integer)
    m = mean(p); den = sum(abs2, p .- m)
    den > 0 || return 0.0
    return sum((p[i] - m) * (p[i+lag] - m) for i in 1:length(p)-lag) / den
end

"""
    autocorrelation_null(profiles, lag; draws, rng)

Mean autocorrelation over the tracks (observed) and its distribution when each
profile is shuffled along the chain (the null of no order along the chain).
"""
function autocorrelation_null(profiles::AbstractVector{<:AbstractVector}, lag::Integer; draws::Integer=2000,
                              rng::AbstractRNG=Random.default_rng())
    obs = mean(profile_autocorrelation(p, lag) for p in profiles)
    null = [mean(profile_autocorrelation(shuffle(rng, p), lag) for p in profiles) for _ in 1:draws]
    return obs, null
end

"""
    count_dispersion(ru; min_scans=2)

Per track of at least `min_scans` scans, the range (max - min) of the per-scan
counts of the counting stage. Returns a vector of (track, scans, range).
"""
function count_dispersion(ru::RunUnits; min_scans::Integer=2)
    out = Tuple{String,Int,Int}[]
    for (t, fs) in ru.tracks
        length(fs) >= min_scans || continue
        n = [ru.n_scan[f] for f in fs]
        push!(out, (t, length(fs), maximum(n) - minimum(n)))
    end
    return sort(out; by=first)
end

end # module
