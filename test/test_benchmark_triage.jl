# Kept/discarded triage applied to three raw benchmark scans.
#   julia --project=. test/test_benchmark_triage.jl
# Skips when the raw scans are absent (STMFIT_DATA_DIR or the default raw link).
using Test
include(joinpath(@__DIR__, "report_benchmark_triage.jl"))

const RAW = get(ENV, "STMFIT_DATA_DIR", joinpath(@__DIR__, "..", "results", "reconstructed_cc_soft_v1", "full146_raw"))
const SCANS = ["240817_002.sxm", "240817_003.sxm", "240817_004.sxm"]

if all(f -> isfile(joinpath(RAW, f)), SCANS)
    @testset "triage" begin
        mktempdir() do tmp
            a = mkpath(joinpath(tmp, "rootA", "20240817_LHe_Cu100")); b = mkpath(joinpath(tmp, "rootB", "copy"))
            for f in SCANS; cp(joinpath(RAW, f), joinpath(a, f)); end
            cp(joinpath(RAW, SCANS[1]), joinpath(b, "same_bytes.sxm"))
            cfg = ["--consensus-config", joinpath(@__DIR__, "..", "config", "molecule_consensus.toml"),
                   "--count-config", joinpath(@__DIR__, "..", "config", "chitosan.toml")]
            redirect_stdout(devnull) do
                main(vcat(["--root", joinpath(tmp, "rootA"), "--outdir", joinpath(tmp, "c1")], cfg))
                main(vcat(["--root", joinpath(tmp, "rootB"), "--outdir", joinpath(tmp, "c2")], cfg))
            end
            cats = [joinpath(tmp, "c1", "catalog.tsv"), joinpath(tmp, "c2", "catalog.tsv")]
            rows = scan_universe(cats)
            @test sort([r["file"] for r in rows]) == SCANS            # same bytes counted once, first catalog wins
            tri = joinpath(tmp, "tri.tsv")
            write(tri, "# test\nfile\tsession\tdecision\tmotif\tsource\tnote\tsha256\n" *
                  join(["$(r["file"])\t$(r["session"])\t$(i == 3 ? "ecarte" : "garde")\t$(i == 3 ? "incertain" : "")\ttest\t\t$(r["sha256"])\n"
                        for (i, r) in enumerate(sort(rows; by=r -> r["file"]))]))
            dest = joinpath(tmp, "dest")
            redirect_stdout(devnull) do
                triage_main(["--catalog", cats[1], "--catalog", cats[2], "--count-config", cfg[4], "--triage", tri, "--dest", dest])
            end
            @test sort(readdir(joinpath(dest, "gardes"))) == SCANS[1:2]
            @test readdir(joinpath(dest, "ecartes")) == SCANS[3:3]
            @test isfile(joinpath(dest, "planches", "gardes", "page_001.png"))
            @test isfile(joinpath(dest, "planches", "ecartes", "incertain", "page_001.png"))
            _, applied = read_table(joinpath(dest, "tri.tsv"))
            @test all(r -> isfile(r["source_path"]), applied)
            short = joinpath(tmp, "short.tsv"); write(short, join(readlines(tri)[1:end-1], "\n") * "\n")
            @test_throws ErrorException triage_main(["--catalog", cats[1], "--count-config", cfg[4],
                                                     "--triage", short, "--dest", joinpath(tmp, "dest2")])
        end
    end
else
    @info "Raw scans absent; triage test skipped" RAW
end
