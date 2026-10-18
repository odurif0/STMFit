# Data index and working sets built from tables, on three raw benchmark scans.
#   julia --project=. test/test_data_sets.jl
# Skips when the raw scans are absent (STMFIT_DATA_DIR or the default raw link).
using Test, SHA
include(joinpath(@__DIR__, "report_data_sets.jl"))

const RAW = get(ENV, "STMFIT_DATA_DIR", joinpath(@__DIR__, "..", "results", "reconstructed_cc_soft_v1", "full146_raw"))
const SCANS = ["240817_002.sxm", "240817_003.sxm", "240817_004.sxm"]

if all(f -> isfile(joinpath(RAW, f)), SCANS)
    @testset "data sets" begin
        mktempdir() do tmp
            data = joinpath(tmp, "data"); sess = mkpath(joinpath(data, "raw", "NKN2", "20240817_LHe_Cu100"))
            for f in SCANS; cp(joinpath(RAW, f), joinpath(sess, f); follow_symlinks=true); end
            cp(joinpath(RAW, SCANS[1]), joinpath(sess, "240817_002 - Copy.sxm"); follow_symlinks=true)
            sha(f) = bytes2hex(open(sha256, joinpath(sess, f)))
            table(rows) = (p = joinpath(tmp, "t$(rand(UInt32)).tsv");
                write(p, "# test\nfile\tset\tmotif\tsource\tnote\tsha256\traw\n" * join(rows .* "\n")); p)
            good = table(["$(SCANS[1])\tbench\t\ttest\t\t$(sha(SCANS[1]))\tNKN2/20240817_LHe_Cu100/$(SCANS[1])",
                          "$(SCANS[2])\t\tincertain\ttest\twhy\t$(sha(SCANS[2]))\tNKN2/20240817_LHe_Cu100/$(SCANS[2])"])
            redirect_stdout(devnull) do; data_main(["--data", data, "--table", good]); end
            link = joinpath(data, "sets", "bench", SCANS[1])
            @test islink(link) && !isabspath(readlink(link)) && read(link) == read(joinpath(sess, SCANS[1]))
            @test readdir(joinpath(data, "sets")) == ["bench"]             # discarded scans get no folder
            _, idx = read_table(joinpath(data, "INDEX.tsv"))
            st = Dict(r["file"] => r for r in idx)
            @test st[SCANS[1]]["status"] == "bench" && st[SCANS[2]]["status"] == "discarded" && st[SCANS[2]]["motif"] == "incertain"
            @test st[SCANS[3]]["status"] == "not_triaged" && st[SCANS[3]]["substrate"] == "Cu(100)" && st[SCANS[3]]["frame_nm"] == "7.00x7.00"
            @test st["240817_002 - Copy.sxm"]["status"] == "duplicate"
            redirect_stdout(devnull) do; data_main(["--data", data, "--table", good]); end   # rebuild in place
            @test islink(link)
            wrong = table(["$(SCANS[1])\tbench\t\ttest\t\t$(sha(SCANS[1]))\tNKN2/20240817_LHe_Cu100/$(SCANS[3])"])
            @test_throws ErrorException data_main(["--data", data, "--table", wrong])
            nomotif = table(["$(SCANS[2])\t\t\ttest\t\t$(sha(SCANS[2]))\tNKN2/20240817_LHe_Cu100/$(SCANS[2])"])
            @test_throws ErrorException data_main(["--data", data, "--table", nomotif])
            write(joinpath(data, "sets", "bench", "real_file.sxm"), "x")
            @test_throws ErrorException data_main(["--data", data, "--table", good])
        end
    end
else
    @info "Raw scans absent; data-set test skipped" RAW
end
