# Joint count of the repeated scans of one molecule from the per-scan GCV curves
# of a counting run (campaign J, stage 1; docs/src/journal.md, 2026-10-26).
#
# batch_full.jl writes, for every scan, the score of each candidate N for the
# elliptical and circular models (<scan>/ell_scores.tsv, circ_scores.tsv). The
# effective curve of a scan is the minimum of the valid elliptical and circular
# scores at each N, the Level 1 of the per-scan hierarchy (docs/src/selection.md).
# The joint score of a molecule adds the scans' log ratios to their own minimum,
# so every scan weighs by the shape of its curve, not by its intensity units.
# Label-free: score tables and molecule tracks only.
module JointCount

export read_score_table, read_score_curves, effective_curve, joint_count, vote_count, split_halves

"""
    read_score_table(path)

Valid rows of an `ell_scores.tsv` or `circ_scores.tsv` table of batch_full.jl, as
`N => score` (second column, the configured criterion). The header repeats the
criterion name, so columns are read by position.
"""
function read_score_table(path::AbstractString)
    out = Dict{Int,Float64}()
    lines = readlines(path)
    isempty(lines) && return out
    h = split(lines[1], '\t')
    length(h) >= 8 && h[1] == "N" && h[8] == "valid" || error("Unexpected score table header: $path")
    for l in lines[2:end]
        isempty(strip(l)) && continue
        v = split(l, '\t')
        v[8] == "true" || continue
        n = parse(Int, v[1]); s = parse(Float64, v[2])
        isfinite(s) || continue
        s > 0 || error("Non-positive score for N=$n: $path")
        out[n] = min(get(out, n, Inf), s)
    end
    return out
end

"Effective curve of one scan: at each N, the smaller of the valid elliptical and circular scores."
effective_curve(ell::AbstractDict, circ::AbstractDict) = mergewith(min, Dict{Int,Float64}(ell), Dict{Int,Float64}(circ))

"""
    read_score_curves(root)

Effective curves of every scan under `root` (a directory `<scan>/` holding
`ell_scores.tsv` and, if present, `circ_scores.tsv`, at any depth), keyed by
`<scan>.sxm`.
"""
function read_score_curves(root::AbstractString)
    curves = Dict{String,Dict{Int,Float64}}()
    for (d, _, names) in walkdir(root)
        "ell_scores.tsv" in names || continue
        f = basename(d) * ".sxm"
        haskey(curves, f) && error("Duplicate scan $f under $root")
        ell = read_score_table(joinpath(d, "ell_scores.tsv"))
        circ = "circ_scores.tsv" in names ? read_score_table(joinpath(d, "circ_scores.tsv")) : Dict{Int,Float64}()
        curves[f] = effective_curve(ell, circ)
    end
    return curves
end

"""
    joint_count(curves) -> (N, J)

Joint count of the scans of one molecule. `J(N) = Σ_s [log c_s(N) − min log c_s]`
over the N valid in at least one curve; a curve without that N contributes its
own largest log ratio. `N` minimises `J` (ties: the smaller N). Empty curves do
not contribute; without any curve the result is `(0, empty)`.
"""
function joint_count(curves)
    cs = [c for c in curves if !isempty(c)]
    isempty(cs) && return 0, Dict{Int,Float64}()
    ns = sort(unique(reduce(vcat, [collect(keys(c)) for c in cs])))
    J = Dict(n => 0.0 for n in ns)
    for c in cs
        lo = log(minimum(values(c)))
        worst = log(maximum(values(c))) - lo
        for n in ns
            J[n] += haskey(c, n) ? log(c[n]) - lo : worst
        end
    end
    return first(sort(ns; by=n -> (J[n], n))), J
end

"Strict-majority count of `ns`, else the lower median."
function vote_count(ns::AbstractVector{<:Integer})
    isempty(ns) && throw(ArgumentError("no counts"))
    tally = Dict{Int,Int}()
    for n in ns
        tally[n] = get(tally, n, 0) + 1
    end
    n, k = first(sort(collect(tally); by=p -> (-last(p), first(p))))
    2k > length(ns) && return n
    return sort(ns)[(length(ns) + 1) ÷ 2]
end

"Odd and even members of a time-ordered list of scans."
split_halves(files::AbstractVector) = (files[1:2:end], files[2:2:end])

end
