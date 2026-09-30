#!/usr/bin/env julia
# Focused regression: julia --project=. test/test_selected_support_context.jl
# Synthetic metadata and coordinate arrays only. No optimizer, SXM, or grading.
using Test
using TOML
using GaussianFit2D

module SelectedSupportExtractor
include(joinpath(@__DIR__, "extract_lobe_features.jl"))
end
const Extractor = SelectedSupportExtractor
const SUPPORT_CONFIG = joinpath(@__DIR__, "..", "config", "chitosan_10_20mer_adaptive_support_rescue.toml")
const ADAPTIVE_POLICY = "adaptive_support_rescue"

# All outcomes emitted by batch_full.jl, including its optional down-only guard.
const SUPPORT_MODES = [(policy=base * suffix, use_rescue=(base == ADAPTIVE_POLICY))
    for base in (ADAPTIVE_POLICY, "adaptive_support_rescue_keep",
                 "adaptive_support_rescue_failed", "adaptive_support_rescue_reject_no_improvement",
                 "adaptive_support_rescue_reject_circ_ell_incoherent", "adaptive_support_rescue_reject_infeasible")
    for suffix in ("", "_robust_guard", "_robust_guard_audit_only")]

function write_summary(dir, header, rows)
    path = joinpath(dir, "selected.tsv")
    open(path, "w") do io
        println(io, join(header, '\t'))
        for row in rows
            println(io, join(row, '\t'))
        end
    end
    return path
end

function mode_contexts(dir)
    rows = [begin
        # A successful guard overwrites selection_source/refined_source. The
        # support choice must survive that change and must not use N_eff.
        source = endswith(mode.policy, "_robust_guard") ? "ell_robust_aicc" :
                 mode.use_rescue ? ADAPTIVE_POLICY : "ell"
        ["/synthetic/scan_$i.sxm", "9", ADAPTIVE_POLICY, source, mode.policy, source, "12", "9"]
    end for (i, mode) in enumerate(SUPPORT_MODES)]
    path = write_summary(dir,
        ["filepath", "N_selected", "selection_policy", "selection_source",
         "refined_policy", "refined_source", "N_eff", "N_refined"], rows)
    return Extractor._read_selected_context(path; selection_policy=ADAPTIVE_POLICY)
end

field_values(config) = Dict(name => deepcopy(getfield(config, name)) for name in fieldnames(typeof(config)))
changed_fields(a, b) = Set(name for name in fieldnames(typeof(a)) if !isequal(getfield(a, name), getfield(b, name)))

@testset "selected support context: Julia 1.13" begin
    @test v"1.13.0" <= VERSION < v"1.14.0"
    @test isempty(Extractor._read_selected_context(""))
    mktempdir() do dir
        @test_throws ErrorException Extractor._read_selected_context(joinpath(dir, "absent.tsv"))
        contexts = mode_contexts(dir)
        @test Set(keys(contexts)) == Set("scan_$i.sxm" for i in eachindex(SUPPORT_MODES))
        for (i, mode) in enumerate(SUPPORT_MODES)
            @test contexts["scan_$i.sxm"] == (n=9, use_rescue=mode.use_rescue, roi=nothing)
        end

        # Legacy nonadaptive summaries still need only these two columns.
        path = write_summary(dir, ["N_selected", "filepath"],
            [["5", "/synthetic/alpha.sxm"], ["11", "beta.sxm"]])
        expected = Dict("alpha.sxm" => (n=5, use_rescue=false, roi=nothing), "beta.sxm" => (n=11, use_rescue=false, roi=nothing))
        @test Extractor._read_selected_context(path) == expected
        @test Extractor._read_selected_context(path; selection_policy="gcv") == expected

        # A known refined outcome is enough even if other policy columns are
        # absent. Do not infer rescue solely from selection_source.
        for mode in SUPPORT_MODES
            path = write_summary(dir, ["filepath", "N_selected", "refined_policy"],
                [["alpha.sxm", "7", mode.policy]])
            @test Extractor._read_selected_context(path)["alpha.sxm"] == (n=7, use_rescue=mode.use_rescue, roi=nothing)
        end
        for policy in ("gcv", "gcv_with_robust_aicc_guard", "support_midpoint_hybrid")
            path = write_summary(dir, ["filepath", "N_selected", "selection_policy", "refined_policy"],
                [["alpha.sxm", "7", policy, "overfit_guard_down_only"]])
            @test Extractor._read_selected_context(path; selection_policy=policy)["alpha.sxm"] == (n=7, use_rescue=false, roi=nothing)
        end
        for mode in filter(m -> !m.use_rescue, SUPPORT_MODES)
            path = write_summary(dir, ["filepath", "N_selected", "refined_policy", "selection_source"],
                [["alpha.sxm", "7", mode.policy, ADAPTIVE_POLICY]])
            @test Extractor._read_selected_context(path)["alpha.sxm"] == (n=7, use_rescue=false, roi=nothing)
        end
    end
end

@testset "ambiguous adaptive metadata fails closed" begin
    mktempdir() do dir
        # The active model alone makes a legacy summary insufficient.
        path = write_summary(dir, ["filepath", "N_selected"], [["alpha.sxm", "9"]])
        @test_throws ErrorException Extractor._read_selected_context(path; selection_policy=ADAPTIVE_POLICY)

        for evidence_column in ("selection_policy", "selection_source", "refined_source")
            path = write_summary(dir, ["filepath", "N_selected", evidence_column],
                [["alpha.sxm", "9", ADAPTIVE_POLICY]])
            @test_throws ErrorException Extractor._read_selected_context(path)
            for unresolved in ("", "NA", "gcv", "unknown_policy", "adaptive_support_rescue_future",
                               "adaptive_support_rescue_reject_unknown", "adaptive_support_rescue_robust_guard_extra",
                               "adaptive_support_rescue_robust_guard_robust_guard", "Adaptive_support_rescue")
                path = write_summary(dir, ["filepath", "N_selected", evidence_column, "refined_policy"],
                    [["alpha.sxm", "9", ADAPTIVE_POLICY, unresolved]])
                @test_throws ErrorException Extractor._read_selected_context(path)
            end
        end
        for unresolved in ("", "NA", "gcv", "unknown_policy", "adaptive_support_rescue_future",
                           "adaptive_support_rescue_reject_unknown", "adaptive_support_rescue_robust_guard_extra")
            path = write_summary(dir, ["filepath", "N_selected", "refined_policy"],
                [["alpha.sxm", "9", unresolved]])
            @test_throws ErrorException Extractor._read_selected_context(path; selection_policy=ADAPTIVE_POLICY)
        end
        for unknown in ("adaptive_support_rescue_future", "adaptive_support_rescue_reject_unknown",
                        "adaptive_support_rescue_keep_robust_guard_robust_guard_audit_only")
            path = write_summary(dir, ["filepath", "N_selected", "refined_policy"],
                [["alpha.sxm", "9", unknown]])
            @test_throws ErrorException Extractor._read_selected_context(path)
        end
    end
end

@testset "malformed selected rows and duplicate basenames fail" begin
    mktempdir() do dir
        for (header, rows) in (
            (["filepath"], [["alpha.sxm"]]),
            (["N_selected"], [["9"]]),
            (["filepath", "N_selected", "N_selected"], [["alpha.sxm", "9", "9"]]),
            (["filepath", "filepath", "N_selected"], [["alpha.sxm", "beta.sxm", "9"]]),
            (["filepath", "N_selected"], [["alpha.sxm"]]),
            (["filepath", "N_selected"], [["alpha.sxm", "9", "extra"]]),
            (["filepath", "N_selected", "refined_policy"], [["alpha.sxm", "9"]]))
            path = write_summary(dir, header, rows)
            @test_throws ErrorException Extractor._read_selected_context(path)
        end
        for n in ("0", "-1", "9.0", "9.5", "NA", "", "nine", "999999999999999999999999999999")
            path = write_summary(dir, ["filepath", "N_selected"], [["alpha.sxm", n]])
            @test_throws ErrorException Extractor._read_selected_context(path)
        end
        for file in ("", "   ", "/")
            path = write_summary(dir, ["filepath", "N_selected"], [[file, "9"]])
            @test_throws ErrorException Extractor._read_selected_context(path)
        end
        for rows in ([["alpha.sxm", "9"], ["alpha.sxm", "9"]],
                     [["alpha.sxm", "9"], ["alpha.sxm", "7"]],
                     [["/synthetic/a/alpha.sxm", "9"], ["/synthetic/b/alpha.sxm", "9"]])
            path = write_summary(dir, ["filepath", "N_selected"], rows)
            @test_throws ErrorException Extractor._read_selected_context(path)
        end
        # A later invalid row must not silently yield a partial context map.
        path = write_summary(dir, ["filepath", "N_selected"], [["alpha.sxm", "9"], ["beta.sxm", "NA"]])
        @test_throws ErrorException Extractor._read_selected_context(path)
    end
end

@testset "metadata reaches both chain configs without leaking" begin
    mktempdir() do dir
        contexts = mode_contexts(dir)
        for profile in ("gaussian", "split")
            cfg = TOML.parsefile(SUPPORT_CONFIG)
            model, preproc = cfg["model"], cfg["preprocessing"]
            model["peak_profile"] = profile
            model_before, preproc_before = deepcopy(model), deepcopy(preproc)
            pcfg_base, ell_base, circ_base = Extractor._configs(model, preproc, dir)
            @test ell_base.n_min == circ_base.n_min == model["n_min"]
            @test ell_base.n_max == circ_base.n_max == model["n_max"]
            @test !ell_base.chain_circular_sigmas && circ_base.chain_circular_sigmas
            for (i, mode) in enumerate(SUPPORT_MODES)
                context = contexts["scan_$i.sxm"]
                pcfg, ell, circ = Extractor._configs(model, preproc, dir; selected_context=context)
                expected_noise = model[mode.use_rescue ? "adaptive_rescue_support_noise_k" : "support_noise_k"]
                expected_padding = model[mode.use_rescue ? "adaptive_rescue_support_padding_nm" : "support_padding_nm"]
                @test isequal(field_values(pcfg), field_values(pcfg_base))
                @test !ell.chain_circular_sigmas && circ.chain_circular_sigmas
                for (chain, base) in ((ell, ell_base), (circ, circ_base))
                    @test chain.n_min == chain.n_max == context.n == 9
                    @test chain.peak_profile == Symbol(profile)
                    @test chain.support_noise_k == expected_noise
                    @test chain.support_padding_nm == expected_padding
                    expected_changed = mode.use_rescue ?
                        Set((:n_min, :n_max, :support_noise_k, :support_padding_nm)) : Set((:n_min, :n_max))
                    @test changed_fields(chain, base) == expected_changed
                end
                @test isequal(model, model_before)
                @test isequal(preproc, preproc_before)
            end
            # A later unselected file uses the original sweep and support; the
            # returned circular/elliptical configs must also be independent.
            _, rescued_ell, rescued_circ = Extractor._configs(model, preproc, dir;
                selected_context=contexts["scan_1.sxm"])
            rescued_ell.support_padding_nm = 99.0
            rescued_ell.n_min = 1
            @test rescued_circ.support_padding_nm == model["adaptive_rescue_support_padding_nm"]
            @test rescued_circ.n_min == 9
            pcfg_after, ell_after, circ_after = Extractor._configs(model, preproc, dir)
            @test isequal(field_values(pcfg_after), field_values(pcfg_base))
            @test isempty(changed_fields(ell_after, ell_base))
            @test isempty(changed_fields(circ_after, circ_base))
            @test isequal(model, model_before)
            @test isequal(preproc, preproc_before)
        end
    end
end

@testset "accepted rescue requires explicit finite physical settings" begin
    mktempdir() do dir
        contexts = mode_contexts(dir)
        rescue, keep = contexts["scan_1.sxm"], contexts["scan_4.sxm"]
        cfg = TOML.parsefile(SUPPORT_CONFIG)
        model, preproc = cfg["model"], cfg["preprocessing"]
        for key in ("adaptive_rescue_support_noise_k", "adaptive_rescue_support_padding_nm")
            invalid_models = Dict{String,Any}[]
            absent = deepcopy(model); delete!(absent, key); push!(invalid_models, absent)
            invalid_values = key == "adaptive_rescue_support_noise_k" ?
                (NaN, Inf, -Inf, 0.0, -0.1, "bad", true) : (NaN, Inf, -Inf, -0.1, "bad", true)
            for value in invalid_values
                bad = deepcopy(model); bad[key] = value; push!(invalid_models, bad)
            end
            for bad in invalid_models
                before, preproc_before = deepcopy(bad), deepcopy(preproc)
                @test_throws ErrorException Extractor._configs(bad, preproc, dir; selected_context=rescue)
                # Inactive rescue settings are not needed for ordinary support.
                for context in (nothing, keep)
                    _, ell, circ = Extractor._configs(bad, preproc, dir; selected_context=context)
                    @test ell.support_noise_k == circ.support_noise_k == model["support_noise_k"]
                    @test ell.support_padding_nm == circ.support_padding_nm == model["support_padding_nm"]
                end
                @test isequal(bad, before)
                @test isequal(preproc, preproc_before)
            end
        end
        # The producer copies the configured values directly: no hidden floor,
        # padding minimum, or replacement with the usual 1.5/0.75 values.
        explicit = deepcopy(model)
        explicit["adaptive_rescue_support_noise_k"] = 0.125
        explicit["adaptive_rescue_support_padding_nm"] = 0.0
        _, ell, circ = Extractor._configs(explicit, preproc, dir; selected_context=rescue)
        @test ell.support_noise_k == circ.support_noise_k == 0.125
        @test ell.support_padding_nm == circ.support_padding_nm == 0.0
    end
end

@testset "fixed N feasibility uses unchanged physical limits, without optimization" begin
    mktempdir() do dir
        contexts = mode_contexts(dir)
        for profile in ("gaussian", "split")
            cfg = TOML.parsefile(SUPPORT_CONFIG)
            model, preproc = cfg["model"], cfg["preprocessing"]
            model["peak_profile"] = profile
            for context in (contexts["scan_1.sxm"], contexts["scan_4.sxm"])
                _, ell, circ = Extractor._configs(model, preproc, dir; selected_context=context)
                for chain in (ell, circ)
                    @test chain.n_min == chain.n_max == 9
                    @test GaussianFit2D._chain_min_span(9, chain) ≈ 4.115845881365984 atol=1e-12
                    for length_nm in (2.86, 1.93)
                        @test !GaussianFit2D._chain_can_fit_support(9, (tmin=0.0, tmax=length_nm), chain)
                    end
                    for length_nm in (5.552874796122559, 5.90756511515228)
                        @test GaussianFit2D._chain_can_fit_support(9, (tmin=0.0, tmax=length_nm), chain)
                    end
                end
            end
        end
    end
end

@testset "synthetic support endpoints replay through the actual fit-data helper" begin
    mktempdir() do dir
        contexts = mode_contexts(dir)
        # A flat central signal gives a contiguous support. Off-axis high values
        # must not enter the tube. Rescue padding alone makes this N=9 feasible.
        t = collect(range(-4.0, 4.0; length=801))
        outside_t = t[1:80:end]
        x = vcat(t, outside_t)
        y = vcat(zeros(length(t)), fill(0.5, length(outside_t)))
        z = vcat(Float64.(abs.(t) .<= 1.5), fill(100.0, length(outside_t)))
        axis = (origin=(0.0, 0.0), axis=[1.0, 0.0], perp=[0.0, 1.0], tmin=-4.0, tmax=4.0)
        inputs_before = deepcopy((x, y, z, axis))
        reference = nothing
        for profile in ("gaussian", "split")
            cfg = TOML.parsefile(SUPPORT_CONFIG)
            model, preproc = cfg["model"], cfg["preprocessing"]
            model["peak_profile"] = profile
            _, base_ell, base_circ = Extractor._configs(model, preproc, dir; selected_context=contexts["scan_4.sxm"])
            _, rescue_ell, rescue_circ = Extractor._configs(model, preproc, dir; selected_context=contexts["scan_1.sxm"])
            for (base, rescue) in ((base_ell, rescue_ell), (base_circ, rescue_circ))
                before = (field_values(base), field_values(rescue))
                # Synthetic off-molecule background is 0 (needed by half_maximum_cap).
                ordinary = GaussianFit2D._chain_fit_data(x, y, z, axis, base; background=0.0)
                expanded = GaussianFit2D._chain_fit_data(x, y, z, axis, rescue; background=0.0)
                for (result, chain) in ((ordinary, base), (expanded, rescue))
                    xf, yf, zf, fitted_axis, keep, meta = result
                    @test meta.support_method == "auto_axis_profile_support" && meta.fallback == "none"
                    @test meta.noise_k == chain.support_noise_k
                    @test meta.padding_nm == chain.support_padding_nm
                    @test meta.final_t_min_nm == fitted_axis.tmin
                    @test meta.final_t_max_nm == fitted_axis.tmax
                    @test meta.final_support_length_nm == fitted_axis.tmax - fitted_axis.tmin
                    @test fitted_axis.origin == axis.origin && fitted_axis.axis == axis.axis && fitted_axis.perp == axis.perp
                    @test keep == ((abs.(y) .<= chain.fit_width_nm) .& (x .>= fitted_axis.tmin) .& (x .<= fitted_axis.tmax))
                    @test (xf, yf, zf) == (x[keep], y[keep], z[keep])
                    @test meta.tube_pixels == length(t) && meta.fit_mask_pixels == count(keep)
                    @test chain.n_min == chain.n_max == 9
                end
                base_axis, rescue_axis = ordinary[4], expanded[4]
                @test rescue_axis.tmin ≈ base_axis.tmin - 0.5 atol=1e-12
                @test rescue_axis.tmax ≈ base_axis.tmax + 0.5 atol=1e-12
                @test expanded[6].threshold_noise < ordinary[6].threshold_noise
                @test !GaussianFit2D._chain_can_fit_support(9, base_axis, base)
                @test GaussianFit2D._chain_can_fit_support(9, rescue_axis, rescue)
                @test count(expanded[5]) > count(ordinary[5])
                @test isequal((field_values(base), field_values(rescue)), before)
                endpoints = (base_axis.tmin, base_axis.tmax, rescue_axis.tmin, rescue_axis.tmax)
                reference === nothing && (reference = endpoints)
                @test endpoints == reference # same reconstruction in both profiles and chain shapes
            end
        end
        @test isequal((x, y, z, axis), inputs_before)
    end
end


@testset "actual CLI rejects unusable selected support before SXM/output" begin
    mktempdir() do dir
        data_dir = mkdir(joinpath(dir, "unreadable_inputs"))
        for file in ("alpha.sxm", "beta.sxm")
            dummy = joinpath(data_dir, file)
            write(dummy, "not an SXM file\n")
            chmod(dummy, 0o000)
        end
        missing_key_config = joinpath(dir, "missing_rescue_key.toml")
        invalid_config = TOML.parsefile(SUPPORT_CONFIG)
        delete!(invalid_config["model"], "adaptive_rescue_support_noise_k")
        open(missing_key_config, "w") do io
            TOML.print(io, invalid_config)
        end
        root = dirname(@__DIR__)
        script = joinpath(@__DIR__, "extract_lobe_features.jl")
        cases = (
            (name="adaptive_n_only", config=SUPPORT_CONFIG, requested="alpha.sxm",
             header=["filepath", "N_selected"], row=["alpha.sxm", "9"],
             error="refined_policy"),
            (name="missing_rescue_key", config=missing_key_config, requested="alpha.sxm",
             header=["filepath", "N_selected", "refined_policy"], row=["alpha.sxm", "9", ADAPTIVE_POLICY],
             error="Accepted adaptive support requires model.adaptive_rescue_support_noise_k"),
            (name="absent_selected_file", config=SUPPORT_CONFIG, requested="beta.sxm",
             header=["filepath", "N_selected", "refined_policy"], row=["alpha.sxm", "9", "adaptive_support_rescue_keep"],
             error="No selected context for requested file: beta.sxm"))
        for case in cases
            path = write_summary(dir, case.header, [case.row])
            output = joinpath(dir, case.name, "features.tsv")
            log = joinpath(dir, case.name * ".log")
            cmd = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$root $script
                   --config $(case.config) --data-dir $data_dir --selected-summary $path
                   --files $(case.requested) --out $output`
            process = open(log, "w") do io
                run(pipeline(ignorestatus(addenv(cmd, "GKSwstype" => "100")), stdout=io, stderr=io))
            end
            diagnostic = read(log, String)
            @test !success(process)
            @test occursin(case.error, diagnostic)
            @test !ispath(output)
            @test !ispath(dirname(output))
        end
    end
end

@testset "registered ROI refit support parsing" begin
    mktempdir() do dir
        cols = ["filepath", "N_selected", "refit_support", "roi_x0_nm", "roi_y0_nm", "roi_x1_nm", "roi_y1_nm"]
        path = write_summary(dir, cols, [["alpha.sxm", "6", "registered_roi", "0.5", "0.0", "3.4", "5.0"]])
        @test Extractor._read_selected_context(path)["alpha.sxm"] == (n=6, use_rescue=false, roi=(0.5, 0.0, 3.4, 5.0))
        path = write_summary(dir, cols, [["alpha.sxm", "6", "scan", "NaN", "NaN", "NaN", "NaN"]])
        @test Extractor._read_selected_context(path)["alpha.sxm"].roi === nothing
        for bad in (["alpha.sxm", "6", "registered_roi", "NaN", "0", "1", "1"], ["alpha.sxm", "6", "registered_roi", "2", "0", "1", "1"],
                    ["alpha.sxm", "6", "elsewhere", "0", "0", "1", "1"])
            path = write_summary(dir, cols, [bad])
            @test_throws ErrorException Extractor._read_selected_context(path)
        end
    end
end
