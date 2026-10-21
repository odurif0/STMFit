#!/usr/bin/env julia
# Campaign A: grade the fully frozen, label-free assignment (grading side:
# reads the benchmark labels; the model itself reads none).
#
# Usage:
#   julia -t N --project=. scripts/grading/grade_frozen_assignment.jl --run RUN_DIR \
#       --config config/unit_assignment_frozen_labelfree.toml --outdir NEW_DIR \
#       [--manifest TOML] [--other-run RUN_DIR] [--ablate cap,weights,fisher]
#
# On a run of the 6-mer benchmark (default manifest chitosan_6mer_benchmark_v2.toml,
# control NKNNKN):
#   1. held-out: each development session (cohorts full145, added_2024) is graded
#      with a model built without it (reference rule on the other sessions,
#      emission rates from all their scans), next to the promoted transductive
#      calls of the run (loso.tsv, loso_files.tsv);
#   2. transfer: the model built on all development sessions is applied to the
#      later sessions (cohort added_2025) (transfer.tsv, transfer_files.tsv;
#      model_development.tsv);
#   3. with --other-run: invariance. The development model is applied to the
#      scans the two runs share; per-scan calls must agree wherever a scan's
#      inputs (assignment features other than emp_fisher, forward and backward
#      patches) agree (invariance.tsv).
# report.md summarises all three, with one vote per molecule, a paired
# bootstrap over molecules (frozen minus transductive), the GlcNAc call rate
# and the errors by chain position (six-lobe scans).
# --ablate (exploratory diagnostics, 2026-10-21) undoes choices of the config:
# cap (no per-track cap on the reference), weights (learned GMM weights),
# fisher (the run's transductive emp_fisher instead of a frozen Fisher).
module Search
include(joinpath(@__DIR__, "search_reference_cohort.jl"))
end
using STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment
using LinearAlgebra, Printf, TOML

const ROOT = dirname(dirname(@__DIR__))

function frozen_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--config", "--outdir", "--manifest", "--other-run", "--ablate") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        haskey(o, args[i]) && error("Repeated option $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--config", "--outdir")) || error("Required: --run --config --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    all(in(("cap", "weights", "fisher")), split(get(o, "--ablate", ""), ","; keepempty=false)) ||
        error("--ablate takes a comma list of cap, weights, fisher")
    return o
end

rows_of(S, sessions) = sort(reduce(vcat, [S.session_rows[s] for s in sessions]))

"Model trained on `sessions` with the config's reference rule."
function session_model(S, R, cfg, st, sessions; freeze_fisher=true)
    files = Set(S.data.file[rows_of(S, sessions)])
    ref = reference_files([r for r in R.consensus if basename(r["filepath"]) in files], st.max_per_track)
    ref_rows = findall(in(ref), R.data.file)
    return train_reference_model(R.data, R.fisher, ref_rows, rows_of(S, sessions), R.physical, cfg; freeze_fisher), length(ref)
end

"GlcNAc call rate and, on six-lobe scans, errors by position (positions 1/6, 2/5, 3/4)."
function call_profile(S, files_calls, truth)
    called = classified = 0
    err = Dict(k => [0, 0] for k in ("ends_1_6", "glcnac_2_5", "middle_3_4"))   # [false GlcNAc, missed GlcNAc]
    n6 = 0
    for (f, calls) in files_calls
        classified += count(>=(0), calls); called += count(==(1), calls)
        length(calls) == 6 || continue
        n6 += 1; t = truth[f]
        for i in 1:6
            calls[i] < 0 && continue
            k = i in (1, 6) ? "ends_1_6" : i in (2, 5) ? "glcnac_2_5" : "middle_3_4"
            calls[i] == 1 && t[i] == 0 && (err[k][1] += 1)
            calls[i] == 0 && t[i] == 1 && (err[k][2] += 1)
        end
    end
    return (rate=called / max(classified, 1), scans6=n6, errors=err)
end

function grade_rows(S, sessions, calls_of, truth)
    grades = Any[]
    for s in sessions
        append!(grades, Search.grade_session(S, s, calls_of(s), truth))
    end
    return grades
end

stats_row(extra, s) = merge(extra, Dict(Search.stat_strings(s)...))

"Per-scan inputs of every lobe of a file, for the invariance check."
function file_inputs(R, f)
    rows = findall(==(f), R.data.file)
    g = copy(R.data.gmm_raw[rows, :]); g[:, end] .= 0      # emp_fisher (last GMM column) is recomputed
    return (n=length(rows), gmm=g, kmeans=[V.X[:, rows] for V in R.data.kmeans], fwd=R.fisher.table.X[rows, :])
end
same_inputs(a, b) = a.n == b.n && isequal(a.gmm, b.gmm) && isequal(a.kmeans, b.kmeans) && isequal(a.fwd, b.fwd)

function frozen_main(args=ARGS)
    o = frozen_options(args)
    o === nothing && return println("grade_frozen_assignment.jl --run RUN_DIR --config TOML --outdir NEW_DIR [--manifest TOML] [--other-run RUN_DIR]")
    BLAS.set_num_threads(1)
    protocol = Dict("data" => Dict("run" => abspath(o["--run"]), "config" => abspath(o["--config"]),
                                   "manifest" => abspath(get(o, "--manifest", joinpath(ROOT, "benchmarks", "chitosan_6mer_benchmark_v2.toml"))),
                                   "control_sequence" => "010010"),
                    "split" => Dict("development_cohorts" => ["full145", "added_2024"], "sealed_cohorts" => ["added_2025"]))
    S = Search.load_setup(protocol)
    cfg = load_config(o["--config"])
    ablate = split(get(o, "--ablate", ""), ","; keepempty=false)
    "cap" in ablate && (cfg["selection"]["reference_max_scans_per_track"] = typemax(Int32))
    "weights" in ablate && (cfg["selection"]["frozen_gmm_weights"] = "learned")
    freeze = !("fisher" in ablate)
    st = load_frozen_settings(cfg)
    R = load_run_assignment(o["--run"], cfg)
    R.data.file == S.data.file && R.data.lobe == S.data.lobe || error("Run tables are not aligned")
    truth = Search.benchmark_truth(S, "010010")
    out = o["--outdir"]; mkpath(out)
    cp(o["--config"], joinpath(out, "config.toml"))
    loso = Dict{String,String}[]; lfiles = Dict{String,String}[]; heldout_calls = Tuple{String,Vector{Int8}}[]
    for s in S.development
        train = filter(!=(s), S.development)
        model, nref = session_model(S, R, cfg, st, train; freeze_fisher=freeze)
        fused, _ = apply_reference_model(model, R.data, R.fisher, S.session_rows[s], R.physical, cfg; freeze_fisher=freeze)
        append!(heldout_calls, [(f, fused[[findfirst(==(r), S.session_rows[s]) for r in S.file_rows[f]]]) for f in S.session_files[s]])
        for (name, calls) in (("frozen", fused), ("transductive_run", S.transductive[S.session_rows[s]]))
            g = Search.grade_session(S, s, calls, truth)
            push!(loso, stats_row(Dict("session" => s, "cohort" => name, "reference_scans" => name == "frozen" ? string(nref) : "",
                                       "theta0" => name == "frozen" ? @sprintf("%.6f", model.theta0) : "",
                                       "theta1" => name == "frozen" ? @sprintf("%.6f", model.theta1) : ""), Search.summarize(g, truth)))
            append!(lfiles, Search.file_rows_out(g, Dict("outer" => s, "cohort" => name)))
        end
        @printf("held-out %s: frozen %s / transductive %s correct\n", s, loso[end-1]["correct"], loso[end]["correct"]); flush(stdout)
    end
    head = vcat(["session", "cohort", "reference_scans", "theta0", "theta1"], Search.STAT_COLUMNS)
    write_table(joinpath(out, "loso.tsv"), head, loso)
    write_table(joinpath(out, "loso_files.tsv"), vcat(["outer", "cohort"], Search.FILE_COLUMNS), lfiles)
    model, nref = session_model(S, R, cfg, st, S.development; freeze_fisher=freeze)
    write_frozen_model(joinpath(out, "model_development.tsv"), model)
    later = rows_of(S, S.sealed)
    fused, _ = apply_reference_model(model, R.data, R.fisher, later, R.physical, cfg; freeze_fisher=freeze)
    pos = Dict(r => k for (k, r) in enumerate(later))
    transfer = Dict{String,String}[]; tfiles = Dict{String,String}[]
    for (name, calls_of) in (("frozen", s -> fused[[pos[r] for r in S.session_rows[s]]]), ("transductive_run", s -> S.transductive[S.session_rows[s]]))
        g = grade_rows(S, S.sealed, calls_of, truth)
        push!(transfer, stats_row(Dict("cohort" => name, "reference_scans" => name == "frozen" ? string(nref) : ""), Search.summarize(g, truth)))
        append!(tfiles, Search.file_rows_out(g, Dict("cohort" => name)))
    end
    write_table(joinpath(out, "transfer.tsv"), vcat(["cohort", "reference_scans"], Search.STAT_COLUMNS), transfer)
    write_table(joinpath(out, "transfer_files.tsv"), vcat(["cohort"], Search.FILE_COLUMNS), tfiles)
    inv = nothing
    if haskey(o, "--other-run")
        R2 = load_run_assignment(o["--other-run"], cfg)
        common = sort(collect(intersect(Set(R.data.file), Set(R2.data.file))))
        f1, p1 = apply_reference_model(model, R.data, R.fisher, collect(eachindex(R.data.file)), R.physical, cfg; freeze_fisher=freeze)
        f2, p2 = apply_reference_model(model, R2.data, R2.fisher, collect(eachindex(R2.data.file)), R2.physical, cfg; freeze_fisher=freeze)
        _, t2 = lobe_table(joinpath(o["--other-run"], "fusion", "predictions_fused.tsv"))
        call(x) = x == "?" ? Int8(-1) : parse(Int8, x)
        rows = Dict{String,String}[]
        for f in common
            i1 = findall(==(f), R.data.file); i2 = findall(==(f), R2.data.file)
            same = same_inputs(file_inputs(R, f), file_inputs(R2, f))
            tr2 = [call(t2[(f, l)]["predicted"]) for l in R2.data.lobe[i2]]
            push!(rows, Dict("file" => f, "same_inputs" => string(same), "N_run" => string(length(i1)), "N_other" => string(length(i2)),
                "frozen_scan_calls_equal" => string(p1.call[i1] == p2.call[i2]), "frozen_fused_equal" => string(f1[i1] == f2[i2]),
                "transductive_equal" => string(length(i1) == length(i2) && S.transductive[i1] == tr2)))
        end
        write_table(joinpath(out, "invariance.tsv"), ["file", "same_inputs", "N_run", "N_other", "frozen_scan_calls_equal",
                                                      "frozen_fused_equal", "transductive_equal"], rows)
        inv = rows
    end
    tot(rows, c, k) = sum(parse(Int, r[k]) for r in rows if r["cohort"] == c)
    open(joinpath(out, "report.md"), "w") do io
        println(io, "# Campaign A: fully frozen, label-free assignment\n")
        println(io, "Config `$(basename(o["--config"]))`; run `$(o["--run"])`." *
                    (isempty(ablate) ? "" : " **Exploratory ablation: $(join(ablate, ", ")) undone.**") * "\n")
        for (title, rows, files) in (("Held-out 2024 sessions (model built without the graded session)", loso, lfiles),
                                     ("Transfer: model built on 2024, applied to 2025", transfer, tfiles))
            println(io, "## $title\n")
            println(io, "| assignment | correct | classified | errors | exact chains | molecules | molecule accuracy | molecules majority exact |")
            println(io, "|---|---|---|---|---|---|---|---|")
            for c in ("frozen", "transductive_run")
                m = Search.molecule_stats([r for r in files if r["cohort"] == c], truth)
                @printf(io, "| %s | %d/%d | %d | %d | %d | %d | %.3f | %d |\n", c, tot(rows, c, "correct"), tot(rows, c, "possible"),
                        tot(rows, c, "classified"), tot(rows, c, "errors"), tot(rows, c, "exact"), m.tracks, m.accuracy, m.majority_exact)
            end
            pm = Search.paired_molecules(files, "transductive_run", "frozen")
            println(io, "\nFrozen minus transductive, per molecule: $(pm.better) better, $(pm.worse) worse, $(pm.same) unchanged; " *
                        "total $(pm.total) units, bootstrap 95% interval [$(pm.lo), $(pm.hi)]; largest $(pm.largest).\n")
        end
        tcalls(files) = [(f, S.transductive[S.file_rows[f]]) for f in files]
        dev_files = [f for s in S.development for f in S.session_files[s]]
        later_files = [f for s in S.sealed for f in S.session_files[s]]
        pos = Dict(r => k for (k, r) in enumerate(later))
        profiles = (("held-out 2024, frozen", call_profile(S, heldout_calls, truth)),
                    ("held-out 2024, transductive", call_profile(S, tcalls(dev_files), truth)),
                    ("2025, frozen", call_profile(S, [(f, fused[[pos[r] for r in S.file_rows[f]]]) for f in later_files], truth)),
                    ("2025, transductive", call_profile(S, tcalls(later_files), truth)))
        println(io, "## Call profile (truth NKNNKN: GlcNAc fraction 1/3 at positions 2 and 5)\n")
        println(io, "| calls | GlcNAc call rate | six-lobe scans | false GlcNAc at ends 1,6 | false GlcNAc at middle 3,4 | missed GlcNAc at 2,5 |")
        println(io, "|---|---|---|---|---|---|")
        for (name, p) in profiles
            e = p.errors
            @printf(io, "| %s | %.3f | %d | %d | %d | %d |\n", name, p.rate, p.scans6, e["ends_1_6"][1], e["middle_3_4"][1],
                    e["glcnac_2_5"][2])
        end
        println(io)
        println(io, "| held-out session | files | frozen | transductive_run | reference scans |")
        println(io, "|---|---|---|---|---|")
        for s in S.development
            r = Dict(x["cohort"] => x for x in loso if x["session"] == s)
            println(io, "| $s | $(r["frozen"]["files"]) | $(r["frozen"]["correct"]) | $(r["transductive_run"]["correct"]) | $(r["frozen"]["reference_scans"]) |")
        end
        if inv !== nothing
            same = [r for r in inv if r["same_inputs"] == "true"]
            println(io, "\n## Invariance (development model applied in both runs)\n")
            println(io, "$(length(inv)) shared scans; $(length(same)) with identical inputs. Among those: per-scan frozen calls " *
                        "identical in $(count(r -> r["frozen_scan_calls_equal"] == "true", same)), fused frozen calls identical in " *
                        "$(count(r -> r["frozen_fused_equal"] == "true", same)), transductive calls identical in " *
                        "$(count(r -> r["transductive_equal"] == "true", same)).")
        end
    end
    print(read(joinpath(out, "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && frozen_main()
