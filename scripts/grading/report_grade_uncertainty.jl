#!/usr/bin/env julia
# Molecule- and session-level view of a benchmark grade, with bootstrap
# intervals over molecules, and paired comparisons of grades (grading side: it
# reads the grades of grade_consensus_run.jl and the benchmark manifest).
#
# Usage:
#   julia --project=. scripts/grading/report_grade_uncertainty.jl --entry NAME=GRADE_DIR,RUN_DIR[,PROFILE] \
#       [--entry NAME=GRADE_DIR,RUN_DIR[,PROFILE] ...] --manifest TOML --outdir NEW_DIR [--draws 2000] [--seed 20261026]
#
# Per scan (GRADE_DIR/units/grades/PROFILE.tsv, default final): exact N (N_pred = N_truth),
# correct units, emitted errors (false positives + false negatives), exact chain.
# Molecule = track of RUN_DIR/consensus/consensus.tsv; session = first folder of
# the manifest's relative_path. Repeat scans of one molecule are not
# independent, so a molecule counts once and every interval resamples
# molecules (with replacement, the scans of a molecule together; 95% percentile
# intervals):
#   molecule-weighted (primary)   the mean over molecules of each molecule's rate;
#   scan-weighted   the usual totals (rates over scans or units), for continuity;
#   per session     the molecule-weighted rates of its molecules.
# With two or more entries, each entry is compared with the first on their
# shared scans, paired by the first entry's molecules: differences of the rates
# with bootstrap intervals. Outputs: report.md, scans.tsv, molecules.tsv,
# sessions.tsv.
using Printf, Random, Statistics, TOML

function uncertainty_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,Any}("--entry" => NTuple{4,String}[], "--draws" => "2000", "--seed" => "20261026"); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--entry", "--manifest", "--outdir", "--draws", "--seed") && i < length(args) || error("Unknown option or missing value: $k")
        if k == "--entry"
            name, rest = split(args[i+1], "="; limit=2)
            parts = split(rest, ",")
            length(parts) in (2, 3) || error("--entry NAME=GRADE_DIR,RUN_DIR[,PROFILE]")
            push!(o[k], (String(name), String(parts[1]), String(parts[2]), length(parts) == 3 ? String(parts[3]) : "final"))
        else
            o[k] = args[i+1]
        end
        i += 2
    end
    isempty(o["--entry"]) && error("Required: --entry NAME=GRADE_DIR,RUN_DIR")
    for k in ("--manifest", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

"Plain TSV reader (grading files)."
function read_tsv(path)
    lines = filter(l -> !isempty(strip(l)) && !startswith(l, "#"), readlines(path))
    header = split(first(lines), '\t')
    return [Dict(String(h) => String(v) for (h, v) in zip(header, split(l, '\t'; keepempty=true))) for l in lines[2:end]]
end

"Graded scan: molecule, session, and the four per-scan quantities."
struct ScanGrade
    file::String
    molecule::String
    session::String
    exact_n::Bool
    correct::Int
    possible::Int
    errors::Int
    exact_chain::Bool
end

function load_entry(grade_dir, run_dir, sessions; profile="final")
    track = Dict(r["file"] => r["track"] for r in read_tsv(joinpath(run_dir, "consensus", "consensus.tsv")))
    out = ScanGrade[]
    for r in read_tsv(joinpath(grade_dir, "units", "grades", profile * ".tsv"))
        f = basename(r["file"])
        push!(out, ScanGrade(f, get(track, f, "scan:" * f), get(sessions, f, ""), r["N_pred"] == r["N_truth"],
                             parse(Int, r["phys_correct"]), parse(Int, r["N_truth"]),
                             parse(Int, r["phys_fp"]) + parse(Int, r["phys_fn"]), r["phys_seq_match"] == "true"))
    end
    return out
end

"Session of every manifest file: the first folder of its relative_path."
manifest_sessions(path) = Dict(f => (p = get(v, "relative_path", ""); occursin('/', p) ? String(split(p, '/')[1]) : "")
                               for (f, v) in TOML.parsefile(path)["files"])

const METRICS = (:exact_n, :correct_units, :errors_per_scan, :exact_chains)

"Scan-weighted rates of a set of scans."
function rates(s::AbstractVector{ScanGrade})
    n = length(s); n == 0 && return (exact_n=NaN, correct_units=NaN, errors_per_scan=NaN, exact_chains=NaN)
    return (exact_n=count(x -> x.exact_n, s) / n, correct_units=sum(x -> x.correct, s) / sum(x -> x.possible, s),
            errors_per_scan=sum(x -> x.errors, s) / n, exact_chains=count(x -> x.exact_chain, s) / n)
end

"Scans grouped by molecule, in a fixed order."
function by_molecule(s::AbstractVector{ScanGrade})
    g = Dict{String,Vector{ScanGrade}}()
    for x in s
        push!(get!(g, x.molecule, ScanGrade[]), x)
    end
    keys_ = sort(collect(keys(g)))
    return keys_, [g[k] for k in keys_]
end

"Mean over molecules of each molecule's rates."
function molecule_rates(groups)
    r = [rates(g) for g in groups]
    return NamedTuple{METRICS}(Tuple(mean(getproperty(x, m) for x in r) for m in METRICS))
end

"Bootstrap of a statistic of the molecule groups (resampled with replacement); 95% percentile intervals per metric."
function bootstrap(stat, groups; draws, rng)
    isempty(groups) && return NamedTuple{METRICS}(Tuple((NaN, NaN) for _ in METRICS))
    vals = [stat(groups[rand(rng, 1:length(groups), length(groups))]) for _ in 1:draws]
    return NamedTuple{METRICS}(Tuple((quantile([getproperty(v, m) for v in vals], 0.025),
                                     quantile([getproperty(v, m) for v in vals], 0.975)) for m in METRICS))
end

pooled(groups) = rates(reduce(vcat, groups; init=ScanGrade[]))

"Paired difference (b − a) of the rates of two grades on their shared scans, grouped by `a`'s molecules."
function paired_groups(a::AbstractVector{ScanGrade}, b::AbstractVector{ScanGrade})
    bb = Dict(x.file => x for x in b)
    shared = [x for x in a if haskey(bb, x.file)]
    _, ga = by_molecule(shared)
    return [(g, [bb[x.file] for x in g]) for g in ga]
end

function paired_stat(pairs; weight=:scan)
    if weight == :scan
        ra = rates(reduce(vcat, first.(pairs); init=ScanGrade[])); rb = rates(reduce(vcat, last.(pairs); init=ScanGrade[]))
    else
        ra = molecule_rates(first.(pairs)); rb = molecule_rates(last.(pairs))
    end
    return NamedTuple{METRICS}(Tuple(getproperty(rb, m) - getproperty(ra, m) for m in METRICS))
end

pct(x) = isnan(x) ? "NA" : @sprintf("%.1f", 100x)
interval_pct(ci) = "$(pct(ci[1])) to $(pct(ci[2]))"
fmt2(x) = isnan(x) ? "NA" : @sprintf("%.2f", x)
show_metric(m, v) = m == :errors_per_scan ? fmt2(v) : pct(v)
show_ci(m, ci) = m == :errors_per_scan ? "$(fmt2(ci[1])) to $(fmt2(ci[2]))" : interval_pct(ci)
const LABEL = Dict(:exact_n => "exact N (%)", :correct_units => "correct units (%)",
                   :errors_per_scan => "errors per scan", :exact_chains => "exact chains (%)")

function uncertainty_main(args=ARGS)
    o = uncertainty_options(args)
    o === nothing && return println("report_grade_uncertainty.jl --entry NAME=GRADE_DIR,RUN_DIR [--entry ...] --manifest TOML --outdir NEW_DIR")
    draws = parse(Int, o["--draws"]); rng = MersenneTwister(parse(Int, o["--seed"]))
    sessions = manifest_sessions(o["--manifest"])
    entries = [(n, load_entry(g, r, sessions; profile=p)) for (n, g, r, p) in o["--entry"]]
    out = o["--outdir"]; mkpath(out)
    md = IOBuffer()
    println(md, "# Grade by molecule and session, with intervals\n")
    println(md, "Molecule = consensus track. A molecule counts once: the primary rates are the means over molecules of each molecule's rate; the scan-weighted rates are the usual totals, where a molecule scanned 43 times counts 43 times. Intervals: 95%, bootstrap over molecules ($draws draws).\n")
    open(joinpath(out, "scans.tsv"), "w") do io
        println(io, "entry\tfile\tmolecule\tsession\texact_n\tcorrect\tpossible\terrors\texact_chain")
        for (n, s) in entries, x in s
            println(io, join([n, x.file, x.molecule, x.session, x.exact_n, x.correct, x.possible, x.errors, x.exact_chain], '\t'))
        end
    end
    mol_io = open(joinpath(out, "molecules.tsv"), "w")
    println(mol_io, "entry\tmolecule\tsession\tscans\texact_n\tcorrect_units\terrors_per_scan\texact_chains")
    ses_io = open(joinpath(out, "sessions.tsv"), "w")
    println(ses_io, "entry\tsession\tmolecules\tscans\texact_n\texact_n_lo\texact_n_hi\tcorrect_units\tcorrect_units_lo\tcorrect_units_hi\terrors_per_scan\texact_chains")
    for (n, s) in entries
        keys_, groups = by_molecule(s)
        println(md, "## $n: $(length(s)) scans, $(length(groups)) molecules\n")
        sc = pooled(groups); scci = bootstrap(pooled, groups; draws, rng)
        mo = molecule_rates(groups); moci = bootstrap(molecule_rates, groups; draws, rng)
        println(md, "| metric | by molecule (primary) | 95% interval | by scan | 95% interval |\n|---|---:|---|---:|---|")
        for m in METRICS
            println(md, "| $(LABEL[m]) | $(show_metric(m, getproperty(mo, m))) | $(show_ci(m, getproperty(moci, m))) | $(show_metric(m, getproperty(sc, m))) | $(show_ci(m, getproperty(scci, m))) |")
        end
        sizes = sort(length.(groups); rev=true)
        println(md, "\nLargest molecules: $(join(sizes[1:min(5, end)], ", ")) scans (the largest holds $(pct(sizes[1] / length(s)))% of the scans).\n")
        for (k, g) in zip(keys_, groups)
            r = rates(g)
            println(mol_io, join([n, k, g[1].session, length(g), r.exact_n, r.correct_units, r.errors_per_scan, r.exact_chains], '\t'))
        end
        println(md, "| session | molecules | scans | exact N (%) | 95% interval | correct units (%) | 95% interval | exact chains (%) |\n|---|---:|---:|---:|---|---:|---|---:|")
        for ses in sort(unique(x.session for x in s))
            sg = [g for g in groups if g[1].session == ses]
            r = molecule_rates(sg)
            ci = length(sg) >= 2 ? bootstrap(molecule_rates, sg; draws, rng) : NamedTuple{METRICS}(Tuple((NaN, NaN) for _ in METRICS))
            println(md, "| $(isempty(ses) ? "(none)" : ses) | $(length(sg)) | $(sum(length, sg)) | $(pct(r.exact_n)) | $(interval_pct(ci.exact_n)) | $(pct(r.correct_units)) | $(interval_pct(ci.correct_units)) | $(pct(r.exact_chains)) |")
            println(ses_io, join([n, ses, length(sg), sum(length, sg), r.exact_n, ci.exact_n..., r.correct_units, ci.correct_units..., r.errors_per_scan, r.exact_chains], '\t'))
        end
        println(md)
    end
    close(mol_io); close(ses_io)
    if length(entries) >= 2
        a_name, a = entries[1]
        println(md, "## Paired differences against `$a_name` (shared scans, its molecules)\n")
        println(md, "| entry | shared scans | molecules | metric | difference by molecule (primary) | 95% interval | difference by scan | 95% interval |\n|---|---:|---:|---|---:|---|---:|---|")
        for (b_name, b) in entries[2:end]
            pairs = paired_groups(a, b)
            isempty(pairs) && continue
            ds = paired_stat(pairs); dm = paired_stat(pairs; weight=:molecule)
            cis = bootstrap(p -> paired_stat(p), pairs; draws, rng); cim = bootstrap(p -> paired_stat(p; weight=:molecule), pairs; draws, rng)
            nscans = sum(length(first(p)) for p in pairs)
            for m in METRICS
                println(md, "| $b_name | $nscans | $(length(pairs)) | $(LABEL[m]) | $(show_metric(m, getproperty(dm, m))) | $(show_ci(m, getproperty(cim, m))) | $(show_metric(m, getproperty(ds, m))) | $(show_ci(m, getproperty(cis, m))) |")
            end
        end
        println(md, "\nA difference whose interval contains 0 is within the molecule-to-molecule variation of this benchmark.")
    end
    write(joinpath(out, "report.md"), take!(md))
    print(read(joinpath(out, "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && uncertainty_main()
