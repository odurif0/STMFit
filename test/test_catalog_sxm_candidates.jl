# Catalog triage and the incremental central folder, on three raw benchmark scans.
#   julia --project=. test/test_catalog_sxm_candidates.jl
# Skips when the raw scans are absent (STMFIT_DATA_DIR or the default raw link).
using Test
include(joinpath(@__DIR__, "catalog_sxm_candidates.jl"))

const RAW = get(ENV, "STMFIT_DATA_DIR", joinpath(@__DIR__, "..", "results", "reconstructed_cc_soft_v1", "full146_raw"))
const SCANS = ["240817_002.sxm", "240817_003.sxm", "240817_004.sxm"]

@testset "component areas" begin
    m = Bool[1 1 0 0; 0 0 0 1; 1 0 0 1]
    @test component_areas(m) == [2, 2, 1]
end

if all(f -> isfile(joinpath(RAW, f)), SCANS)
    @testset "central folder" begin
        mktempdir() do tmp
            a = mkpath(joinpath(tmp, "rootA", "20240817_LHe_Cu100")); b = mkpath(joinpath(tmp, "rootB", "session"))
            cp(joinpath(RAW, SCANS[1]), joinpath(a, SCANS[1])); cp(joinpath(RAW, SCANS[2]), joinpath(a, SCANS[2]))
            cp(joinpath(RAW, SCANS[1]), joinpath(b, "renamed_copy.sxm"))        # same bytes, other name
            cp(joinpath(RAW, SCANS[3]), joinpath(b, SCANS[3]))
            ex = joinpath(tmp, "exclude.tsv"); write(ex, "file\treason\n$(SCANS[2])\tnot chitosan\n")
            cfg = ["--consensus-config", joinpath(@__DIR__, "..", "config", "molecule_consensus.toml"),
                   "--count-config", joinpath(@__DIR__, "..", "config", "chitosan.toml")]
            central = joinpath(tmp, "central")
            redirect_stdout(devnull) do
                main(vcat(["--roots", joinpath(tmp, "rootA"), "--outdir", joinpath(tmp, "c1"), "--exclude", ex,
                           "--copy-to", central], cfg))
            end
            _, rows = read_table(joinpath(tmp, "c1", "catalog.tsv"))
            st = Dict(r["file"] => r["status"] for r in rows)
            @test st[SCANS[1]] == "usable" && st[SCANS[2]] == "excluded_by_review"
            @test only(r for r in rows if r["file"] == SCANS[1])["session"] == "20240817_LHe_Cu100"
            @test sort(filter(f -> endswith(f, ".sxm"), readdir(central))) == [SCANS[1]]
            redirect_stdout(devnull) do
                main(vcat(["--root", joinpath(tmp, "rootB"), "--outdir", joinpath(tmp, "c2"), "--copy-to", central], cfg))
            end
            _, rows = read_table(joinpath(tmp, "c2", "catalog.tsv"))
            dup = only(r for r in rows if r["file"] == "renamed_copy.sxm")
            @test dup["status"] == "duplicate" && dup["duplicate_of"] == joinpath(central, SCANS[1])
            @test sort(filter(f -> endswith(f, ".sxm"), readdir(central))) == sort([SCANS[1], SCANS[3]])
            _, man = read_table(joinpath(central, "MANIFEST.tsv"))
            @test [m["central_name"] for m in man] == sort([SCANS[1], SCANS[3]])
            @test all(m -> m["sha256"] == bytes2hex(open(sha256, joinpath(central, m["central_name"]))), man)
            mkpath(joinpath(tmp, "bare"))
            @test_throws ErrorException catalog_options(vcat(["--roots", a, "--outdir", joinpath(tmp, "c3"),
                                                              "--copy-to", joinpath(tmp, "bare")], cfg))
        end
    end
else
    @info "Raw scans absent; central-folder test skipped" RAW
end
