using Test, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
const CONFIG = joinpath(@__DIR__, "..", "config", "unit_assignment_reconstructed.toml")

@testset "Explicit transverse half-plane descriptor" begin
    coords = collect(-4:4) .* 0.08
    p = zeros(81)
    p[5] = -1.0                # u=-4, t=0
    p[77] = 2.0                # u=+4, t=0
    @test transverse_asymmetry(p, coords; zero_l1=1e-12) == (1.0, "ok")
    @test first(transverse_asymmetry(7p, coords; zero_l1=1e-12)) == 1.0
    @test first(transverse_asymmetry(-p, coords; zero_l1=1e-12)) == -1.0
    reflected_u = vcat([p[(i-1)*9+1:i*9] for i in 9:-1:1]...)
    reflected_t = vcat([reverse(p[(i-1)*9+1:i*9]) for i in 1:9]...)
    @test first(transverse_asymmetry(reflected_u, coords; zero_l1=1e-12)) == -1.0
    @test first(transverse_asymmetry(reflected_t, coords; zero_l1=1e-12)) == 1.0
    t_only = zeros(81); t_only[45] = 1.0
    @test first(transverse_asymmetry(t_only, coords; zero_l1=1e-12)) == 0.0
    @test first(transverse_asymmetry(ones(81), coords; zero_l1=1e-12)) == 0.0
    @test last(transverse_asymmetry(zeros(81), coords; zero_l1=1e-12)) == "zero_patch_mass"
    bad = copy(p); bad[1] = NaN
    @test last(transverse_asymmetry(bad, coords; zero_l1=1e-12)) == "nonfinite_patch"
    @test_throws ErrorException transverse_asymmetry(ones(80), coords; zero_l1=1e-12)
    @test_throws ErrorException transverse_asymmetry(p, coords; zero_l1=-1)
    # This reconstruction is half-plane parity, NOT the u-weighted first moment.
    near = zeros(81); near[50] = 1.0  # u=+1, t=0
    far = zeros(81); far[77] = 1.0    # u=+4, t=0
    @test first(transverse_asymmetry(near, coords; zero_l1=1e-12)) ==
          first(transverse_asymmetry(far, coords; zero_l1=1e-12)) == 1.0
end

@testset "Descriptor table, explicit invalid rows and no benchmark inputs" begin
    mktempdir() do dir
        features = joinpath(dir, "features_full145_name_is_allowed.tsv")
        patches = joinpath(dir, "patches.tsv")
        output = joinpath(dir, "augmented.tsv")
        rows = [Dict("file" => "unknown.sxm", "lobe" => string(i), "amplitude" => "0.1") for i in 1:2]
        write_table(features, ["file", "lobe", "amplitude"], rows)
        columns = ["bwd_res_p" * lpad(string(i), 3, '0') for i in 1:81]
        patchrows = [merge(copy(row), Dict(c => "0.0" for c in columns)) for row in rows]
        patchrows[1]["bwd_res_p077"] = "1.0"
        write_table(patches, vcat(["file", "lobe"], columns), patchrows)
        augment_descriptor(features, patches, output, CONFIG)
        header, actual = lobe_table(output)
        @test "patch_u_asym_reconstructed" in header
        @test actual[("unknown.sxm", 1)]["patch_u_asym_reconstructed"] == "1"
        @test actual[("unknown.sxm", 2)]["patch_u_asym_reconstructed"] == "NA"
        @test actual[("unknown.sxm", 2)]["descriptor_reason"] == "zero_patch_mass"
        @test length(actual) == 2
        @test_throws ErrorException augment_descriptor(features, patches, output, CONFIG)
        for field in ("expected_N", "sequence", "truth", "control_sequence")
            forbidden = joinpath(dir, field * ".tsv")
            write_table(forbidden, ["file", "lobe", field], [Dict("file" => "x.sxm", "lobe" => "1", field => "x")])
            @test_throws ErrorException read_table(forbidden)
        end
        duplicates = joinpath(dir, "duplicate.tsv")
        write_table(duplicates, ["file", "lobe"], [rows[1], rows[1]])
        @test_throws ErrorException lobe_table(duplicates)
        gap = joinpath(dir, "gap.tsv")
        write_table(gap, ["file", "lobe"], [rows[2]])
        @test_throws ErrorException lobe_table(gap)
        @test_throws ErrorException require_same_keys(actual, Dict(("other.sxm", 1) => rows[1]), "fixture")
    end
end

@testset "Soft vote keeps all lobes and unchanged finite-pair arithmetic" begin
    mktempdir() do dir
        features, km, gm, out = [joinpath(dir, name * ".tsv") for name in ("features", "km", "gm", "out")]
        rows = [Dict("file" => "unknown.sxm", "lobe" => string(i)) for i in 1:3]
        write_table(features, ["file", "lobe"], rows)
        kmrows = [merge(rows[i], Dict("predicted" => "1", "probability_1" => ["0.2", "0.8", "0.5"][i])) for i in 1:3]
        gmrows = [merge(rows[i], Dict("predicted" => "1", "probability_1" => ["0.4", "0.2"][i])) for i in 1:2]
        write_table(km, ["file", "lobe", "predicted", "probability_1"], kmrows)
        write_table(gm, ["file", "lobe", "predicted", "probability_1"], gmrows)
        write_soft_vote(features, km, gm, out, CONFIG)
        _, actual = lobe_table(out)
        @test length(actual) == 3
        @test actual[("unknown.sxm", 1)]["predicted"] == "0"
        @test actual[("unknown.sxm", 1)]["confidence"] == "0.40000000"
        @test actual[("unknown.sxm", 2)]["predicted"] == "1"  # exact tie, frozen threshold
        @test actual[("unknown.sxm", 3)]["predicted"] == "?"
        @test actual[("unknown.sxm", 3)]["invalid_reason"] == "missing_gmm"
    end
end
