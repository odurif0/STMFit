# Frozen unit assignment (src/frozen_assignment.jl) and the campaign-B search
# helpers (scripts/grading/search_reference_cohort.jl).
#   julia --project=. test/test_frozen_assignment.jl
using Test, Random, Printf, LinearAlgebra
using STMFit, STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment, STMFit.MoleculeFusion, STMFit.EmpiricalFisherNative
module KM
include(joinpath(@__DIR__, "..", "scripts", "pipeline", "build_labelfree_unit_predictions.jl"))
end
module GMM
include(joinpath(@__DIR__, "..", "scripts", "pipeline", "build_labelfree_gmm_predictions.jl"))
end
module Search
include(joinpath(@__DIR__, "..", "scripts", "grading", "search_reference_cohort.jl"))
end
const ROOT = dirname(@__DIR__)
const CORROBORATED = joinpath(ROOT, "config", "unit_assignment_corroborated_training.toml")
const FROZEN = joinpath(ROOT, "config", "unit_assignment_frozen_reference.toml")
const LABELFREE = joinpath(ROOT, "config", "unit_assignment_frozen_labelfree.toml")
const FEATURES = ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel",
                  "patch_u_asym_reconstructed", "mold_cc_fwd", "mold_cc_bwd", "emp_fisher", "split_log_skew"]
BLAS.set_num_threads(1)

"Synthetic run: 50 six-lobe scans, two lobe types, a few missing values."
function synthetic_run(dir)
    rng = MersenneTwister(11)
    feats = joinpath(dir, "features_predictor.tsv"); patches = joinpath(dir, "patches_bwd17.tsv")
    summary = joinpath(dir, "consensus_summary_final.tsv"); fwd = joinpath(dir, "patches_fwd17.tsv")
    open(feats, "w") do f; open(patches, "w") do p; open(summary, "w") do s; open(fwd, "w") do fw
        println(f, join(vcat(["file", "lobe", "amplitude"], FEATURES), '\t'))
        println(p, join(vcat(["file", "lobe"], [@sprintf("bwd_res_p%03d", i) for i in 1:289]), '\t'))
        println(fw, join(vcat(["file", "lobe"], [@sprintf("res_p%03d", i) for i in 1:289]), '\t'))
        println(s, "filepath\tN_selected\tcount_rule\ttrack")
        for scan in 1:50
            file = @sprintf("s%02d.sxm", scan)
            println(s, file, "\t6\t", scan % 3 == 0 ? "consensus_applied" : "agrees", "\t", cld(scan, 10))
            for lobe in 1:6
                hi = lobe in (2, 5)
                vals = [(hi ? 1.0 : 0.0) + 0.6randn(rng) for _ in FEATURES]
                scan == 7 && lobe == 3 && (vals[5] = NaN)
                println(f, join(vcat([file, string(lobe), @sprintf("%.6f", (hi ? 0.2 : 0.1) + 0.02randn(rng))],
                                     [isfinite(v) ? @sprintf("%.6f", v) : "NaN" for v in vals]), '\t'))
                px = [@sprintf("%.5f", randn(rng) - (hi ? 0.3 : 0.0)) for _ in 1:289]
                scan == 9 && lobe == 1 && (px[5] = "NaN")
                println(p, join(vcat([file, string(lobe)], px), '\t'))
                fx = [@sprintf("%.5f", randn(rng) + (hi ? 0.5 * exp(-((i - 145)^2) / 400) : 0.0)) for i in 1:289]
                scan == 11 && lobe == 2 && (fx[145] = "NaN")
                println(fw, join(vcat([file, string(lobe)], fx), '\t'))
            end
        end
    end; end; end; end
    return feats, patches, summary, fwd
end

@testset "training on the corroborated scans reproduces the production scripts" begin
    mktempdir() do dir
        feats, patches, summary, _ = synthetic_run(dir)
        base4 = join(FEATURES[1:4], ",")
        KM.main(["--features", feats, "--patches", patches, "--out", joinpath(dir, "km.tsv"), "--first-seed", "0",
                 "--seeds", "20", "--interactions", "--training-scans", summary,
                 "--view", "v_base=$base4", "--view", "v_split=$base4,split_log_skew",
                 "--view", "v_comt=$base4,bwd_neg_com_t", "--view", "v_diag45=$base4,bwd_neg_diag45"])
        GMM.main(["--features", feats, "--config", CORROBORATED, "--out", joinpath(dir, "gmm.tsv"), "--first-seed", "0",
                  "--seeds", "10", "--selftrain", "2", "--interactions", "--training-scans", summary,
                  "--view", "v_cc=$base4,patch_u_asym_reconstructed,mold_cc_fwd,mold_cc_bwd,emp_fisher"])
        write_soft_vote(feats, joinpath(dir, "km.tsv"), joinpath(dir, "gmm.tsv"), joinpath(dir, "vote.tsv"), CORROBORATED)
        data = load_assignment_data(feats, patches)
        cfg = load_config(FROZEN)
        keep = read_corroborated_scans(summary)
        model = train_frozen(data, findall(in(keep), data.file), cfg)
        rows = collect(eachindex(data.file))
        pred = predict_frozen(model, data, rows)
        f8(p) = isfinite(p) ? @sprintf("%.8f", p) : "NA"
        for (path, v) in (("km.tsv", pred.p_kmeans), ("gmm.tsv", pred.p_gmm), ("vote.tsv", pred.p_vote))
            _, t = lobe_table(joinpath(dir, path); contiguous=false)
            @test all(f8(v[i]) == t[(data.file[i], data.lobe[i])]["probability_1"] for i in rows)
        end
        _, t = lobe_table(joinpath(dir, "vote.tsv"))
        @test all((pred.call[i] < 0 ? "?" : string(pred.call[i])) == t[(data.file[i], data.lobe[i])]["predicted"] for i in rows)
        @test count(<(0), pred.call) == 1                     # the lobe with a missing GMM feature

        # a frozen model applies unchanged to a subset and survives a write/read round trip
        m2 = STMFit.FrozenAssignment.with_fusion(model, 0.05, 0.8, 0.5)
        write_frozen_model(joinpath(dir, "model.tsv"), m2)
        m3 = read_frozen_model(joinpath(dir, "model.tsv"))
        sub = rows[13:60]
        @test predict_frozen(m3, data, sub).call == pred.call[13:60]
        @test (m3.theta0, m3.theta1, m3.pi) == (0.05, 0.8, 0.5)
        @test_throws ErrorException write_frozen_model(joinpath(dir, "model.tsv"), m2)
    end
end

@testset "fully frozen model: Fisher per parity, reference rule, application" begin
    mktempdir() do dir
        feats, patches, summary, fwd = synthetic_run(dir)
        cfg = load_config(LABELFREE)
        st = load_frozen_settings(cfg)
        @test st.max_per_track == 5 && st.equal_weights && st.pi == 0.5
        data = load_assignment_data(feats, patches)
        fd = load_fisher_data(fwd, cfg, data)
        rows = collect(eachindex(data.file))
        # whole cohort as reference = the production lobe-parity CV, as the files carry it
        cv = cv_scores(load_patches(fwd, "res", load_fisher_config(cfg)), load_fisher_config(cfg))
        fisher = train_fisher(fd, rows)
        expect(x) = isfinite(x.score) ? parse(Float64, @sprintf("%.8g", -parse(Float64, @sprintf("%.6f", x.score)))) : NaN
        @test isequal([fisher_feature(fd, fisher, i) for i in rows], expect.(cv))
        @test isnan(fisher_feature(fd, fisher, findfirst(i -> data.file[i] == "s11.sxm" && data.lobe[i] == 2, rows)))
        # reference rule: corroborated scans, at most 5 per track, first by file name
        _, cons = read_table(summary)
        ref = reference_files(cons, 5)
        @test ref == Set(["s01.sxm", "s02.sxm", "s04.sxm", "s05.sxm", "s07.sxm", "s11.sxm", "s13.sxm", "s14.sxm", "s16.sxm",
                          "s17.sxm", "s22.sxm", "s23.sxm", "s25.sxm", "s26.sxm", "s28.sxm", "s31.sxm", "s32.sxm", "s34.sxm",
                          "s35.sxm", "s37.sxm", "s41.sxm", "s43.sxm", "s44.sxm", "s46.sxm", "s47.sxm"])
        @test all(length([f for f in ref if cld(parse(Int, f[2:3]), 10) == t]) <= 5 for t in 1:5)
        physical = ["track$(cld(parse(Int, f[2:3]), 10))_$(l)" for (f, l) in zip(data.file, data.lobe)]
        ref_rows = findall(in(ref), data.file)
        model = train_reference_model(data, fd, ref_rows, rows, physical, cfg)
        @test length(model.fisher) == 2 && model.pi == 0.5 && model.theta1 > model.theta0
        @test all(g.weights == [0.5, 0.5] for g in model.gmm)
        fused, p = apply_reference_model(model, data, fd, rows, physical, cfg)
        # the model never sees the other scans: applying to a subset gives the same per-scan calls
        sub = findall(f -> f in ("s03.sxm", "s30.sxm"), data.file)
        _, psub = apply_reference_model(model, data, fd, sub, physical, cfg)
        @test psub.call == p.call[sub]
        write_frozen_model(joinpath(dir, "m.tsv"), model)
        m2 = read_frozen_model(joinpath(dir, "m.tsv"))
        @test apply_reference_model(m2, data, fd, rows, physical, cfg)[1] == fused
        @test [x === nothing ? nothing : (x.w_p, x.mid) for x in m2.fisher] == [x === nothing ? nothing : (x.w_p, x.mid) for x in model.fisher]
    end
end

@testset "fixed-parameter fusion equals the latent-class posterior" begin
    calls = Int8[1, 1, 0, 1, 0, 0, -1, 1, 0]
    phys = ["a", "a", "a", "b", "b", "b", "b", "c", ""]
    st = FusionSettings(3, 500, 1e-10, 0.5, 0.1, 0.9)
    fit = fit_emission_rates(calls, phys, st)
    fused = fuse_calls(calls, phys, fit.theta0, fit.theta1, fit.pi, 3)
    @test fused[8:9] == calls[8:9]                            # fewer than 3 calls, or no track: unchanged
    @test fused[7] == fused[4]                                # "?" takes its physical lobe's call
    x3 = ["x", "x", "x"]
    @test fuse_calls(Int8[1, 1, 0], x3, 0.05, 0.8, 0.5, 3) == Int8[1, 1, 1]
    @test fuse_calls(Int8[1, 0, 0], x3, 0.05, 0.8, 0.5, 3) == Int8[0, 0, 0]
    @test fuse_calls(Int8[1, 0, 0], x3, 0.01, 0.8, 0.5, 3) == Int8[1, 1, 1]   # rare false calls: one detection suffices
end

@testset "search helpers" begin
    t = Int8[0, 1, 0, 0, 1, 0]
    @test Search.grade_file(Int8[0, 1, 0, 0, 1, 0], t) == (n=6, correct=6, classified=6, exact=true)
    @test Search.grade_file(Int8[0, 1, 0, 0, -1, 0], t) == (n=6, correct=5, classified=5, exact=false)
    g = Search.grade_file(Int8[1, 0, 0, 1, 0], t)             # N=5: the reverse 0,1,0,0,1 matches the first 5 truth units
    @test (g.correct, g.classified) == (5, 5)
    @test Search.grade_file(Int8[-1, -1], t).correct == 0
    # floating search: removing 3, then 5 reaches the target; nothing else gains
    target = Set([1, 2, 4, 6])
    obj(sets) = [10 - length(symdiff(Set(c), target)) - (4 in c ? 0 : (2 in c ? 0 : 3)) for c in sets]
    C, value, trace = Search.floating_search([1, 2, 3, 4, 5, 6], obj)
    @test Set(C) == target && value == 10
    @test [(x.action, x.unit) for x in trace] == [("start", 0), ("remove", 3), ("remove", 5)]
    C, _, _ = Search.floating_search([1, 2, 3], s -> zeros(Int, length(s)))
    @test C == [1, 2, 3]                                      # no gain: nothing removed
    # re-addition: after removing 1, 2 and 3, putting 1 back is the best step
    v = Dict(Set([2, 3, 4]) => 5, Set([3, 4]) => 6, Set([4]) => 7, Set([1, 4]) => 9, Set([2, 4]) => 2, Set([2, 3]) => 2,
             Set([3]) => 2, Set([1, 3, 4]) => 1, Set([1, 2, 4]) => 1, Set([1, 2, 3]) => 1)
    C, value, trace = Search.floating_search([1, 2, 3, 4], sets -> [get(v, Set(c), 0) for c in sets])
    @test C == [1, 4] && value == 9
    @test [(x.action, x.unit) for x in trace] == [("start", 0), ("remove", 1), ("remove", 2), ("remove", 3), ("add", 1)]
    # local search: flips reach the target from a random start; swaps escape a flip optimum
    target = Set([1, 2, 4, 6])
    f(sets) = [10 - length(symdiff(Set(c), target)) for c in sets]
    C, value, trace = Search.local_search([1, 2, 3, 4, 5, 6], [3, 4], f)
    @test Set(C) == target && value == 10
    @test [x.action for x in trace] == ["start", "add", "add", "remove", "add"]
    w = Dict(Set([1, 2]) => 3, Set([1, 3]) => 5)          # from {1,2}: no flip gains, the swap 2 -> 3 does
    C, value, trace = Search.local_search([1, 2, 3], [1, 2], sets -> [get(w, Set(c), 0) for c in sets])
    @test C == [1, 3] && value == 5 && (trace[end].action, trace[end].unit, trace[end].added) == ("swap", 2, 3)
    C, _, _ = Search.local_search([1, 2, 3], [1, 2], sets -> [get(w, Set(c), 0) for c in sets]; swap_top=0)
    @test C == [1, 2]
    sp = Dict("start_seeds" => [1, 2], "start_keep_probability" => 0.75)
    @test Search.start_set([5, 1, 3], sp, 0) == [1, 3, 5]
    @test Search.start_set(collect(1:40), sp, 1) == Search.start_set(collect(1:40), sp, 1) != Search.start_set(collect(1:40), sp, 2)
    @test_throws ErrorException Search.start_set([1, 2], sp, 3)
end
