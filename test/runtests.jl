# Every test of the repository, each file in its own Julia process (the
# scripts the tests load define globals that would clash in one session).
#
#   julia --project=. test/runtests.jl          # all tests, end-to-end run included
#   julia --project=. test/runtests.jl --fast   # without the end-to-end run
#   julia --project=. -e 'using Pkg; Pkg.test()'
#
# A failing file prints its output; the run ends with an error listing them.
using Printf

const ROOT = dirname(@__DIR__)

files = vcat([joinpath(ROOT, "packages", p, "test", "runtests.jl") for p in ("STMFitCore.jl", "STMSXMIO.jl", "GaussianFit2D.jl")],
             sort(filter(f -> startswith(basename(f), "test_"), readdir(@__DIR__; join=true))),
             sort(filter(f -> endswith(f, ".jl"), readdir(joinpath(@__DIR__, "molecule_consensus"); join=true))))
runs = [(f, String[]) for f in files]
"--fast" in ARGS || push!(runs, (joinpath(@__DIR__, "test_reconstructed_pipeline.jl"), ["--e2e"]))

failed = String[]
for (f, args) in runs
    out = IOBuffer()
    t = @elapsed proc = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --project=$ROOT $f $args`); stdout=out, stderr=out))
    name = strip(relpath(f, ROOT) * " " * join(args, " "))
    @printf("%-4s %5.0f s  %s\n", success(proc) ? "PASS" : "FAIL", t, name)
    success(proc) || (push!(failed, name); print(String(take!(out))))
end
isempty(failed) || error("$(length(failed)) test file(s) failed: " * join(failed, ", "))
println("All $(length(runs)) test runs passed.")
