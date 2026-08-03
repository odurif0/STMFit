using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
const FINALIZER = joinpath(@__DIR__, "finalize_qe_mold_workflow.jl")

@testset "finalize QE workflow accepts regex-derived prefixes" begin
    mktempdir() do tmp
        index_dir = joinpath(tmp, "indices")
        mkpath(index_dir)

        run_dirs = String[]
        for prefix in ("glcn_central", "glcnac_central")
            run_dir = joinpath(tmp, prefix)
            mkpath(run_dir)
            push!(run_dirs, run_dir)
            write(joinpath(run_dir, "pp_ldos.in"), "&INPUTPP\n  prefix = '$prefix'\n/\n")
            write(joinpath(run_dir, "$(prefix)_relaxed.xyz"), "1\nfixture\nH 0 0 0\n")
            write(joinpath(run_dir, "$(prefix)_ldos.cube"), "fixture\n")
            write(joinpath(index_dir, "$(prefix)_trimer_indices.tsv"), """
key\tvalue
n_atoms\t1
origin_indices_bare\t1
axis_from_index_bare\t1
axis_to_index_bare\t1
plane_index_bare\t1
""")
        end

        cmd = `$(Base.julia_cmd()) --project=$ROOT $FINALIZER --glcn-dir $(run_dirs[1]) --glcnac-dir $(run_dirs[2]) --index-dir $index_dir --height-nm 0.5 --dry-run`
        result = run(ignorestatus(cmd))
        @test result.exitcode == 0
    end
end
