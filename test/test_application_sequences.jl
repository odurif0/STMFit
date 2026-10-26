# Sequence listing of an application run with per-unit robustness, on a synthetic calls table.
#   julia --project=. test/test_application_sequences.jl
using Test
include(joinpath(@__DIR__, "..", "scripts", "data", "list_application_sequences.jl"))

@testset "sequences and robustness from a calls table" begin
    mktempdir() do tmp
        # scan a: units out of chain order in the table; scan b: one unclassified unit
        rows = ["file\tlobe\tphysical_lobe\td\tmodels_classified\tbright_share\tfull\tother",
                "a.sxm\t2\tt1_2\t2\t10\t0.900\t1\t1",
                "a.sxm\t1\tt1_1\t1\t10\t0.000\t0\t0",
                "a.sxm\t3\tt1_3\t2\t10\t0.500\t0\t1",
                "b.sxm\t1\tt2_1\t1\t10\t0.200\t0\t0",
                "b.sxm\t2\tt2_2\t2\t0\tNA\t\t",
                "c.sxm\t1\tt3_1\t1\t10\t1.000\t1\t1"]
        calls = joinpath(tmp, "calls.tsv"); write(calls, join(rows, "\n") * "\n")
        scans = joinpath(tmp, "scans.tsv"); write(scans, "file\tdecision\nb.sxm\taccept\na.sxm\taccept\nmissing.sxm\taccept\n")
        out = joinpath(tmp, "out")
        redirect_stdout(devnull) do
            sequences_main(["--calls", calls, "--scans", scans, "--outdir", out])
        end
        lines = readlines(joinpath(out, "sequences.tsv"))
        @test lines[1] == "file\tunits\tcalls\trobustness\tbright\trobust_bright\tunstable"
        @test lines[2] == "b.sxm\t2\t0?\t.?\t0\t0\t0"             # scans follow the table order; NA is unknown
        @test lines[3] == "a.sxm\t3\t010\t.#~\t1\t1\t1"            # lobes are put in chain order
        @test length(lines) == 3                                   # a scan absent from the calls table is not listed
        @test occursin("absent from the calls table", read(joinpath(out, "report.md"), String))
        # another model column, all scans, and a stricter threshold
        out2 = joinpath(tmp, "out2")
        redirect_stdout(devnull) do
            sequences_main(["--calls", calls, "--model", "other", "--threshold", "0.95", "--outdir", out2])
        end
        l2 = readlines(joinpath(out2, "sequences.tsv"))
        @test l2[2] == "a.sxm\t3\t011\t.~~\t2\t0\t2" && l2[end] == "c.sxm\t1\t1\t#\t1\t1\t0"
        @test_throws ErrorException sequences_main(["--calls", calls, "--outdir", out2])   # output exists
        @test_throws ErrorException sequences_main(["--calls", calls, "--model", "nope", "--outdir", joinpath(tmp, "out3")])
    end
end
