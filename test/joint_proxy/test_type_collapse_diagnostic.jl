#!/usr/bin/env julia

using Test
using Statistics

include(joinpath(@__DIR__, "..", "diagnose_joint_proxy_type_collapse.jl"))
include(joinpath(@__DIR__, "..", "merge_joint_proxy_type_ablation.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")

@testset "joint-proxy type-collapse diagnostic" begin
    registry = load_registry(CONFIG)
    mktempdir() do dir
        output = joinpath(dir, "synthetic.tsv")
        options = DiagnosticOptions(config=CONFIG, synthetic_out=output,
            seed=20260721, replicates=2, noise_sigmas=[0.0, 0.10])
        rows_a = synthetic_rows(registry, options)
        rows_b = synthetic_rows(registry, options)

        @test rows_a == rows_b
        @test length(rows_a) == 576
        @test Set(row.ablation for row in rows_a) == Set([
            "geometric", "stm_dft_v1", "combined",
            "geometric_matched", "stm_dft_v1_matched", "combined_matched",
            "geometric_generative", "stm_dft_v1_generative", "combined_generative",
        ])
        @test Set(row.injected_type for row in rows_a) == Set([0, 1])

        dft_self = filter(row -> row.injected_family == "stm_dft_v1" &&
            row.ablation == "stm_dft_v1" && row.noise_sigma == 0.0, rows_a)
        @test count(row -> row.injected_type == 0 && row.correct, dft_self) == 8
        @test count(row -> row.injected_type == 1 && row.correct, dft_self) == 8
        @test all(row -> row.true_margin > 0, dft_self)

        dft_matched = filter(row -> row.injected_family == "stm_dft_v1" &&
            row.ablation == "stm_dft_v1_matched" && row.noise_sigma == 0.0, rows_a)
        @test count(row -> row.injected_type == 0 && row.correct, dft_matched) == 8
        @test count(row -> row.injected_type == 1 && row.correct, dft_matched) == 8
        @test all(row -> row.true_margin > 0, dft_matched)

        dft_generative = filter(row -> row.injected_family == "stm_dft_v1" &&
            row.ablation == "stm_dft_v1_generative" && row.noise_sigma == 0.0, rows_a)
        @test count(row -> row.injected_type == 0 && row.correct, dft_generative) == 8
        @test count(row -> row.injected_type == 1 && row.correct, dft_generative) == 8
        @test all(row -> row.true_margin > 0, dft_generative)

        coords = collect(-registry.grid_half_nm:registry.grid_step_nm:registry.grid_half_nm)
        points = [(; t, u) for u in coords for t in coords]
        patch = [0.3 + 2.0 * exp(-0.5 * ((point.t / 0.20)^2 +
            (point.u / 0.14)^2)) for point in points]
        @test maximum(abs, backbone_residual(patch, points, 0.20, 0.14)) < 1e-10
        @test backbone_residual(-patch, points, 0.20, 0.14) ≈
            -backbone_residual(patch, points, 0.20, 0.14)

        dft_entry = only(entry for entry in registry.entries if entry.source.family == "stm_dft_v1")
        template0 = only(template for template in dft_entry.templates if
            template.type == 0 && template.parity == 0 && template.mirror == 0)
        template1 = only(template for template in dft_entry.templates if
            template.type == 1 && template.parity == 0 && template.mirror == 0)
        gaussian = [exp(-0.5 * ((point.t / 0.20)^2 + (point.u / 0.14)^2)) for point in points]
        joint_patch = 0.3 .+ 2.0 .* gaussian .+ 0.8 .* template0.pixels
        @test joint_generative_score(joint_patch, template0.pixels, points, 0.20, 0.14) >
            joint_generative_score(joint_patch, template1.pixels, points, 0.20, 0.14)
        gaussian_only = 0.3 .+ 2.0 .* gaussian
        @test joint_generative_score(gaussian_only, template0.pixels, points, 0.20, 0.14) ≈ 0.0 atol=1e-10
        inverse_patch = gaussian_only .- 0.8 .* template0.pixels
        @test joint_generative_score(inverse_patch, template0.pixels, points, 0.20, 0.14) ≈ 0.0 atol=1e-10

        combined = last(family_ensembles(registry)).second
        synthetic_candidate = (; lobes=[(; sigma_parallel_nm=0.20, sigma_perp_nm=0.14)],
            patch_tu=points)
        two_view_result = generative_type_posterior(synthetic_candidate,
            [Dict("fwd" => joint_patch, "bwd" => joint_patch)], combined; effective_factor=1.0)
        @test two_view_result !== nothing
        @test two_view_result.lobe_marginals[1, 1] > two_view_result.lobe_marginals[1, 2]
        null_result = generative_type_posterior(synthetic_candidate,
            [Dict("fwd" => gaussian_only, "bwd" => gaussian_only)], combined; effective_factor=1.0)
        @test null_result !== nothing
        @test null_result.lobe_marginals[1, :] ≈ [0.5, 0.5] atol=1e-12
        sparse_patch = fill(NaN, length(gaussian_only))
        sparse_patch[1:2] .= gaussian_only[1:2]
        @test generative_type_posterior(synthetic_candidate,
            [Dict("fwd" => sparse_patch)], combined; effective_factor=1.0) === nothing
        @test_throws ErrorException generative_type_posterior(synthetic_candidate,
            [Dict{String,Vector{Float64}}()], combined; effective_factor=1.0)
        view_sets = generative_patch_sets([
            TypePosteriorLobeEvidence(Dict("fwd" => joint_patch, "bwd" => joint_patch)),
        ])
        @test first.(view_sets) == ["image", "image_bwd", "image_fwd"]
        @test generative_view_factor(last(view_sets[1]), 0.4) == 0.4
        @test generative_view_factor(last(view_sets[2]), 0.4) == 1.0

        affine_bwd = 1.7 .* joint_patch .- 0.25
        affine_metrics = view_transfer_metrics(joint_patch, affine_bwd, registry.grid_n, 1)
        @test affine_metrics.n_pairs == registry.npix
        @test affine_metrics.corr ≈ 1.0 atol=1e-12
        @test affine_metrics.slope ≈ 1.7 atol=1e-12
        @test affine_metrics.intercept ≈ -0.25 atol=1e-12
        @test affine_metrics.affine_nrmse < 1e-12
        @test (affine_metrics.best_dt_steps, affine_metrics.best_du_steps) == (0, 0)

        inverse_metrics = view_transfer_metrics(joint_patch,
            -2.0 .* joint_patch .+ 0.4, registry.grid_n, 0)
        @test inverse_metrics.corr ≈ -1.0 atol=1e-12
        @test inverse_metrics.slope ≈ -2.0 atol=1e-12

        shifted = fill(NaN, registry.npix)
        for iu in 1:registry.grid_n, it in 1:(registry.grid_n - 1)
            source = (iu - 1) * registry.grid_n + it
            target = (iu - 1) * registry.grid_n + it + 1
            shifted[target] = joint_patch[source]
        end
        shift_metrics = view_transfer_metrics(joint_patch, shifted, registry.grid_n, 1)
        @test (shift_metrics.best_dt_steps, shift_metrics.best_du_steps) == (1, 0)
        @test abs(shift_metrics.best_corr) > abs(shift_metrics.corr)
        degenerate = [TypePosteriorLobeEvidence(Dict("fwd" => zeros(registry.npix)))]
        @test diagnostic_posterior(degenerate, combined; effective_factor=1.0) === nothing

        write_tsv(output, rows_a)
        header = first(readlines(output))
        @test !occursin("expected_N", header)
        @test !occursin("benchmark", lowercase(header))
        @test !occursin("sequence", lowercase(header))
    end

    @test_throws ErrorException parse_options(["--config", CONFIG])
    @test_throws ErrorException parse_options([
        "--config", CONFIG, "--real-out", "x.tsv", "--data-dir", "/tmp",
    ])
    @test_throws ErrorException parse_options([
        "--config", CONFIG, "--real-out", "x.tsv", "--real-summary", "s.tsv",
        "--data-dir", "/tmp", "--files-from", "f.txt", "--view-out", "v.tsv",
    ])
    parsed_view = parse_options([
        "--config", CONFIG, "--real-out", "x.tsv", "--real-summary", "s.tsv",
        "--data-dir", "/tmp", "--files-from", "f.txt", "--view-out", "v.tsv",
        "--view-max-shift-steps", "1",
    ])
    @test parsed_view.view_out == "v.tsv"
    @test parsed_view.view_max_shift_steps == 1
    @test diagnostic_file_seed(20260721, "a.sxm") == diagnostic_file_seed(20260721, "a.sxm")
    @test diagnostic_file_seed(20260721, "a.sxm") != diagnostic_file_seed(20260721, "b.sxm")
    deterministic_chain = diagnostic_chain_config(diagnostic_file_seed(20260721, "a.sxm"))
    @test deterministic_chain.skip_global
    @test deterministic_chain.rng_seed == diagnostic_file_seed(20260721, "a.sxm")

    mktempdir() do dir
        duplicate_list = joinpath(dir, "files.txt")
        write(duplicate_list, "duplicate.sxm\nduplicate.sxm\n")
        duplicate_options = DiagnosticOptions(config=CONFIG, data_dir=dir,
            files_from=duplicate_list, real_out=joinpath(dir, "real.tsv"),
            real_summary=joinpath(dir, "summary.tsv"))
        @test_throws ErrorException read_files(duplicate_options)

        empty_table = joinpath(dir, "empty.tsv")
        write_tsv(empty_table, NamedTuple[]; fields=REAL_LOBE_FIELDS)
        @test first(readlines(empty_table)) == join(REAL_LOBE_FIELDS, '\t')
    end

    mktempdir() do dir
        shards = [joinpath(dir, "shard_01"), joinpath(dir, "shard_02")]
        for (i, shard) in enumerate(shards)
            write_tsv(joinpath(shard, LOBE_NAME), [
                (; file="same.sxm", n=i == 1 ? 10 : 2, lobe=1,
                    patch_kind="residual", ablation="combined"),
            ])
            write_tsv(joinpath(shard, SUMMARY_NAME), [
                (; file="file_$i.sxm", patch_kind="residual", ablation="combined"),
            ])
        end
        outdir = joinpath(dir, "merged")
        @test merge_table(shards, LOBE_NAME, ["file", "n", "lobe", "patch_kind", "ablation"],
            joinpath(outdir, LOBE_NAME)) == 2
        @test merge_table(shards, SUMMARY_NAME, ["file", "patch_kind", "ablation"],
            joinpath(outdir, SUMMARY_NAME)) == 2
        @test length(readlines(joinpath(outdir, LOBE_NAME))) == 3
        merged_lobes = readlines(joinpath(outdir, LOBE_NAME))
        @test split(merged_lobes[2], '\t')[2] == "2"
        @test split(merged_lobes[3], '\t')[2] == "10"
        @test_throws ErrorException merge_table([shards[1], shards[1]], LOBE_NAME,
            ["file", "n", "lobe", "patch_kind", "ablation"], joinpath(outdir, "duplicate.tsv"))
    end
end
