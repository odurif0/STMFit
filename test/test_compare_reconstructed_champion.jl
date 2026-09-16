using Test
include(joinpath(@__DIR__, "compare_reconstructed_champion.jl"))

@testset "External output comparison, not a production selector" begin
    mktempdir() do dir
        reference = joinpath(dir, "ref.tsv")
        predictions = joinpath(dir, "pred.tsv")
        header = ["file", "lobe", "predicted", "confidence"]
        base = [Dict("file" => "a.sxm", "lobe" => string(i), "predicted" => string(i-1), "confidence" => "0.5") for i in 1:2]
        write_table(reference, header, base)
        changed = deepcopy(base)
        changed[1]["predicted"] = "1"
        changed[2]["predicted"] = "0"
        changed[2]["confidence"] = "0.6"
        write_table(predictions, header, changed)
        result = compare_predictions(predictions, reference, joinpath(dir, "different"))
        @test result["same_assignment"] == "0"
        @test result["identical_key_set"] == "true"
        @test result["identical_assignments"] == "false"
        @test parse(Float64, result["max_confidence_abs_difference"]) ≈ .1
        identical = compare_predictions(reference, reference, joinpath(dir, "same"))
        @test identical["identical_reported_confidences"] == "true"
        subset = joinpath(dir, "subset.tsv")
        write_table(subset, header, base[1:1])
        partial = compare_predictions(subset, reference, joinpath(dir, "partial"))
        @test partial["missing_lobes"] == "1"
        @test partial["identical_assignments"] == "false"
        @test_throws ErrorException compare_predictions(reference, reference, joinpath(dir, "same"))
        manifest = joinpath(dir, "external.toml")
        write(manifest, "[files.\"a.sxm\"]\nexpected_N = 2\n")
        compare_predictions(reference, reference, joinpath(dir, "grade_rows"); benchmark_manifest=manifest)
        @test isfile(joinpath(dir, "grade_rows", "predictions_for_external_grade.tsv"))
    end
end
