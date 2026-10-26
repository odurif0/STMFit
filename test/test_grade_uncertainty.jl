# Molecule- and session-level grade with bootstrap intervals, on synthetic grades.
#   julia --project=. test/test_grade_uncertainty.jl
using Test, Random
include(joinpath(@__DIR__, "..", "scripts", "grading", "report_grade_uncertainty.jl"))

sg(f, m, s; exact=true, correct=6, errors=0, chain=true) = ScanGrade(f, m, s, exact, correct, 6, errors, chain)

@testset "rates, molecule weighting and bootstrap" begin
    # one big perfect molecule (8 scans) and two poor molecules (1 scan each)
    s = vcat([sg("a$i", "A", "S1") for i in 1:8],
             [sg("b", "B", "S1"; exact=false, correct=3, errors=3, chain=false), sg("c", "C", "S2"; correct=4, errors=2, chain=false)])
    r = rates(s)
    @test r.exact_n == 0.9 && r.exact_chains == 0.8
    @test isapprox(r.correct_units, (8 * 6 + 3 + 4) / 60) && isapprox(r.errors_per_scan, 0.5)
    keys_, groups = by_molecule(s)
    @test keys_ == ["A", "B", "C"]
    m = molecule_rates(groups)
    @test isapprox(m.exact_chains, 1 / 3) && isapprox(m.exact_n, 2 / 3)    # every molecule counts once
    ci = bootstrap(pooled, groups; draws=500, rng=MersenneTwister(1))
    @test ci.exact_chains[1] < r.exact_chains <= ci.exact_chains[2]
    @test 0 <= ci.exact_chains[1] && ci.exact_chains[2] <= 1
end

@testset "paired differences" begin
    a = [sg("x1", "M1", "S"), sg("x2", "M1", "S"), sg("y1", "M2", "S"; correct=4, errors=2, chain=false)]
    b = [sg("x1", "other", "S"; correct=5, errors=1, chain=false), sg("x2", "other", "S"), sg("y1", "M9", "S"; correct=4, errors=2, chain=false),
         sg("z", "M3", "S")]                                                   # z is not shared
    p = paired_groups(a, b)
    @test length(p) == 2 && sum(length(first(x)) for x in p) == 3                # a's molecules, shared scans only
    d = paired_stat(p)
    @test isapprox(d.correct_units, -1 / 18) && isapprox(d.exact_chains, -1 / 3) && d.exact_n == 0
    same = paired_groups(a, a)
    @test all(iszero, values(paired_stat(same))) && all(iszero, values(paired_stat(same; weight=:molecule)))
end

@testset "command line on grade and run folders" begin
    mktempdir() do tmp
        g = mkpath(joinpath(tmp, "grade", "units", "grades")); r = mkpath(joinpath(tmp, "run", "consensus"))
        write(joinpath(g, "final.tsv"), "file\tN_pred\tN_truth\tphys_correct\tphys_fp\tphys_fn\tphys_seq_match\n" *
              "s1.sxm\t6\t6\t6\t0\t0\ttrue\ns2.sxm\t6\t6\t5\t1\t0\tfalse\ns3.sxm\t5\t6\t4\t0\t1\tfalse\n")
        write(joinpath(r, "consensus.tsv"), "file\ttrack\ns1.sxm\t1\ns2.sxm\t1\ns3.sxm\t2\n")
        man = joinpath(tmp, "m.toml")
        write(man, "[files.\"s1.sxm\"]\nrelative_path = \"sesA/s1.sxm\"\n[files.\"s2.sxm\"]\nrelative_path = \"sesA/s2.sxm\"\n[files.\"s3.sxm\"]\nrelative_path = \"sesB/s3.sxm\"\n")
        out = joinpath(tmp, "out")
        redirect_stdout(devnull) do
            uncertainty_main(["--entry", "x=$(joinpath(tmp, "grade")),$(joinpath(tmp, "run"))", "--entry", "y=$(joinpath(tmp, "grade")),$(joinpath(tmp, "run"))",
                              "--manifest", man, "--outdir", out, "--draws", "200"])
        end
        rep = read(joinpath(out, "report.md"), String)
        @test occursin("x: 3 scans, 2 molecules", rep) && occursin("| sesA | 1 | 2 |", rep) && occursin("Paired differences", rep)
        @test length(readlines(joinpath(out, "molecules.tsv"))) == 1 + 4 && length(readlines(joinpath(out, "sessions.tsv"))) == 1 + 4
        @test_throws ErrorException uncertainty_main(["--entry", "x=a,b", "--manifest", man, "--outdir", out])
    end
end
