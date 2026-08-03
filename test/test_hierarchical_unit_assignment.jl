#!/usr/bin/env julia

# T4 tests for the composition-free hierarchical equal-prior assignment model.
#
# Two test sets:
#   1. baseline_schema        — characterization of the EXISTING label-free
#                               prediction builder (`build_labelfree_unit_predictions.jl`):
#                               pins the prediction TSV schema and the
#                               higher-amplitude -> GlcNAc (1) physical mapping.
#                               Must run GREEN on unchanged code before the new
#                               hierarchical model is exercised.
#   2. hierarchical_model     — failing-first tests for the new
#                               `test/lib/hierarchical_unit_assignment.jl` model
#                               and the `test/build_hierarchical_unit_predictions.jl`
#                               CLI. Composition-free: equal priors, no lobe
#                               order, no transitions, no run lengths, no
#                               class-count/top-k/occupancy controls, no
#                               benchmark/truth columns.
#
# Run: julia --project=. test/test_hierarchical_unit_assignment.jl

using Printf
using Random
using SHA
using Statistics
using Test

const ROOT = dirname(@__DIR__)
const LIB = joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl")
const NEW_CLI = joinpath(@__DIR__, "build_hierarchical_unit_predictions.jl")
const OLD_CLI = joinpath(@__DIR__, "build_labelfree_unit_predictions.jl")

# --------------------------------------------------------------------------
# Fixture helpers (deterministic synthetic feature TSVs).
# --------------------------------------------------------------------------

const BASE_FEATURES = ["amp_prominence", "amp_neighbor_ratio",
                       "integrated_prominence", "amp_rel"]

# Write a synthetic feature TSV. Rows are Dicts; missing columns are skipped.
function _write_feature_tsv(path::AbstractString, rows, extra_cols=String[])
    open(path, "w") do io
        cols = vcat(["file", "N", "lobe", "amplitude", "x_nm", "y_nm",
                     "sigma_parallel_nm", "sigma_perp_nm", "integrated"],
                    BASE_FEATURES, extra_cols)
        println(io, join(cols, '\t'))
        for r in rows
            vals = [string(get(r, c, "")) for c in cols]
            println(io, join(vals, '\t'))
        end
    end
    return path
end

# Build a deterministic two-cluster scan. The HIGH cluster (cluster=1) has the
# larger amplitude and is the one that must map to GlcNAc (1) under the
# physical convention. feature_shift_high/low set the per-feature location.
# A small multiplicative lobe-dependent jitter is applied so per-scan MAD is
# non-zero even when several lobes share a cluster tier.
function _two_cluster_scan(; file::AbstractString,
                             n_low::Int, n_high::Int,
                             amp_low::Real, amp_high::Real,
                             feature_shift_low::Real, feature_shift_high::Real,
                             lobe_start::Int=1,
                             amp_jitter::Real=0.01,
                             feat_jitter::Real=0.05)
    rows = Dict{String,Any}[]
    lobe = lobe_start
    for _ in 1:n_low
        amp = Float64(amp_low) * (1.0 + amp_jitter * sin(lobe))
        fmul = 1.0 + feat_jitter * sin(lobe)
        fv = Float64(feature_shift_low) * fmul
        push!(rows, Dict(
            "file" => file, "N" => n_low + n_high, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", amp),
            "x_nm" => @sprintf("%.6f", 0.1 * lobe),
            "y_nm" => @sprintf("%.6f", 0.05 * lobe),
            "sigma_parallel_nm" => "0.35", "sigma_perp_nm" => "0.40",
            "integrated" => @sprintf("%.10g", amp * 0.35 * 0.40),
            "amp_prominence" => @sprintf("%.6f", fv),
            "amp_neighbor_ratio" => @sprintf("%.6f", fv * 0.5),
            "integrated_prominence" => @sprintf("%.6f", fv * 0.8),
            "amp_rel" => @sprintf("%.6f", fv * 0.6)))
        lobe += 1
    end
    for _ in 1:n_high
        amp = Float64(amp_high) * (1.0 + amp_jitter * sin(lobe))
        fmul = 1.0 + feat_jitter * sin(lobe)
        fv = Float64(feature_shift_high) * fmul
        push!(rows, Dict(
            "file" => file, "N" => n_low + n_high, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", amp),
            "x_nm" => @sprintf("%.6f", 0.1 * lobe),
            "y_nm" => @sprintf("%.6f", 0.05 * lobe),
            "sigma_parallel_nm" => "0.45", "sigma_perp_nm" => "0.50",
            "integrated" => @sprintf("%.10g", amp * 0.45 * 0.50),
            "amp_prominence" => @sprintf("%.6f", fv),
            "amp_neighbor_ratio" => @sprintf("%.6f", fv * 0.5),
            "integrated_prominence" => @sprintf("%.6f", fv * 0.8),
            "amp_rel" => @sprintf("%.6f", fv * 0.6)))
        lobe += 1
    end
    return rows
end

# Scan with affine nuisance applied to features: each feature f -> a*f + b.
# This must NOT change the assignment because per-scan median+MAD profiling
# removes affine nuisance. Includes lobe-dependent multiplicative jitter so
# the within-scan MAD is non-zero.
function _affine_scan(; file::AbstractString,
                       n_low::Int, n_high::Int,
                       amp_low::Real, amp_high::Real,
                       f_low::Real, f_high::Real,
                       scale::Real, shift::Real,
                       lobe_start::Int=1,
                       feat_jitter::Real=0.05)
    rows = Dict{String,Any}[]
    lobe = lobe_start
    for _ in 1:n_low
        amp = Float64(amp_low) * (1.0 + 0.01 * sin(lobe))
        fmul = 1.0 + feat_jitter * sin(lobe)
        s = scale * (Float64(f_low) * fmul) + shift
        push!(rows, Dict(
            "file" => file, "N" => n_low + n_high, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", amp),
            "x_nm" => @sprintf("%.6f", 0.1 * lobe),
            "y_nm" => @sprintf("%.6f", 0.05 * lobe),
            "sigma_parallel_nm" => "0.35", "sigma_perp_nm" => "0.40",
            "integrated" => @sprintf("%.10g", amp * 0.35 * 0.40),
            "amp_prominence" => @sprintf("%.6f", s),
            "amp_neighbor_ratio" => @sprintf("%.6f", s * 0.5),
            "integrated_prominence" => @sprintf("%.6f", s * 0.8),
            "amp_rel" => @sprintf("%.6f", s * 0.6)))
        lobe += 1
    end
    for _ in 1:n_high
        amp = Float64(amp_high) * (1.0 + 0.01 * sin(lobe))
        fmul = 1.0 + feat_jitter * sin(lobe)
        s = scale * (Float64(f_high) * fmul) + shift
        push!(rows, Dict(
            "file" => file, "N" => n_low + n_high, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", amp),
            "x_nm" => @sprintf("%.6f", 0.1 * lobe),
            "y_nm" => @sprintf("%.6f", 0.05 * lobe),
            "sigma_parallel_nm" => "0.45", "sigma_perp_nm" => "0.50",
            "integrated" => @sprintf("%.10g", amp * 0.45 * 0.50),
            "amp_prominence" => @sprintf("%.6f", s),
            "amp_neighbor_ratio" => @sprintf("%.6f", s * 0.5),
            "integrated_prominence" => @sprintf("%.6f", s * 0.8),
            "amp_rel" => @sprintf("%.6f", s * 0.6)))
        lobe += 1
    end
    return rows
end

# Scan where every lobe is drawn from a single feature population (no
# bimodality). The hierarchical 2-component model must detect this and abstain.
function _one_component_scan(; file::AbstractString, n::Int, amp::Real,
                              feature_value::Real, lobe_start::Int=1)
    rows = Dict{String,Any}[]
    for k in 1:n
        lobe = lobe_start + k - 1
        a = Float64(amp) * (1.0 + 0.005 * sin(lobe))
        fv = Float64(feature_value) + 0.01 * cos(lobe)
        push!(rows, Dict(
            "file" => file, "N" => n, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", a),
            "x_nm" => @sprintf("%.6f", 0.1 * lobe),
            "y_nm" => @sprintf("%.6f", 0.05 * lobe),
            "sigma_parallel_nm" => "0.40", "sigma_perp_nm" => "0.45",
            "integrated" => @sprintf("%.10g", a * 0.40 * 0.45),
            "amp_prominence" => @sprintf("%.6f", fv),
            "amp_neighbor_ratio" => @sprintf("%.6f", fv * 0.5),
            "integrated_prominence" => @sprintf("%.6f", fv * 0.8),
            "amp_rel" => @sprintf("%.6f", fv * 0.6)))
    end
    return rows
end

# Scan with degenerate features: the same value for every lobe (zero variance,
# zero MAD). Must be invalidated for that scan/view.
function _degenerate_scan(; file::AbstractString, n::Int, amp::Real,
                           feature_value::Real, lobe_start::Int=1)
    rows = Dict{String,Any}[]
    for k in 1:n
        lobe = lobe_start + k - 1
        push!(rows, Dict(
            "file" => file, "N" => n, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", Float64(amp)),
            "x_nm" => @sprintf("%.6f", 0.1 * lobe),
            "y_nm" => @sprintf("%.6f", 0.05 * lobe),
            "sigma_parallel_nm" => "0.40", "sigma_perp_nm" => "0.45",
            "integrated" => @sprintf("%.10g", Float64(amp) * 0.40 * 0.45),
            "amp_prominence" => @sprintf("%.6f", Float64(feature_value)),
            "amp_neighbor_ratio" => @sprintf("%.6f", Float64(feature_value) * 0.5),
            "integrated_prominence" => @sprintf("%.6f", Float64(feature_value) * 0.8),
            "amp_rel" => @sprintf("%.6f", Float64(feature_value) * 0.6)))
    end
    return rows
end

# Seed-sensitive scan: a uniform amp_prominence ramp plus two large-magnitude
# outliers. The deterministic EM multistart is seed-sensitive on this shape:
# the narrowest quantile start (delta 0.25, multistart #1) recovers the ramp
# body, while later starts seed the means on the outliers and converge to a
# different, lower-likelihood basin. Used by the B2 regression to prove
# --first-seed actually reaches the EM multistart seam.
function _seed_sensitive_scan(; file::AbstractString)
    rows = Dict{String,Any}[]
    ramp = collect(range(-2.0, 2.0; length=40))
    for (i, fv) in enumerate(ramp)
        lobe = i
        amp = 0.05 + 0.001 * i
        push!(rows, Dict(
            "file" => file, "N" => length(ramp) + 2, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", amp),
            "x_nm" => "0.0", "y_nm" => "0.0",
            "sigma_parallel_nm" => "0.40", "sigma_perp_nm" => "0.45",
            "integrated" => @sprintf("%.10g", amp * 0.40 * 0.45),
            "amp_prominence" => @sprintf("%.6f", fv),
            "amp_neighbor_ratio" => @sprintf("%.6f", fv * 0.5),
            "integrated_prominence" => @sprintf("%.6f", fv * 0.8),
            "amp_rel" => @sprintf("%.6f", fv * 0.6)))
    end
    for (k, fv) in enumerate((8.0, 9.0))
        lobe = length(ramp) + k
        amp = 0.20
        push!(rows, Dict(
            "file" => file, "N" => length(ramp) + 2, "lobe" => lobe,
            "amplitude" => @sprintf("%.10g", amp),
            "x_nm" => "0.0", "y_nm" => "0.0",
            "sigma_parallel_nm" => "0.40", "sigma_perp_nm" => "0.45",
            "integrated" => @sprintf("%.10g", amp * 0.40 * 0.45),
            "amp_prominence" => @sprintf("%.6f", fv),
            "amp_neighbor_ratio" => @sprintf("%.6f", fv * 0.5),
            "integrated_prominence" => @sprintf("%.6f", fv * 0.8),
            "amp_rel" => @sprintf("%.6f", fv * 0.6)))
    end
    return rows
end

function _run(cmd::Cmd)
    out = IOBuffer()
    err = IOBuffer()
    p = run(pipeline(cmd; stdout=out, stderr=err), wait=false)
    wait(p)
    return p.exitcode, String(take!(out)), String(take!(err))
end

_julia_cmd(script) = `$(Base.julia_cmd()) --project=$ROOT $script`

function _read_predictions(path::AbstractString)
    lines = readlines(path)
    data = filter(l -> !isempty(strip(l)) && !startswith(strip(l), '#'), lines)
    @test length(data) >= 1
    header = split(data[1], '\t'; keepempty=true)
    rows = Dict{String,String}[]
    for line in data[2:end]
        vals = split(line, '\t'; keepempty=true)
        row = Dict{String,String}()
        for (i, h) in enumerate(header)
            row[h] = i <= length(vals) ? vals[i] : ""
        end
        push!(rows, row)
    end
    return String.(header), rows
end

# ===========================================================================
@testset "baseline_schema" begin
    # Characterize the EXISTING label-free prediction builder so any drift in
    # the frozen prediction schema or the physical mapping convention is
    # caught before the new hierarchical model runs.
    mktempdir() do tmp
        feats = joinpath(tmp, "features.tsv")
        rows = vcat(
            _two_cluster_scan(;
                  file="scan_A.sxm",
                  n_low=3, n_high=3,
                  amp_low=0.05, amp_high=0.10,
                  feature_shift_low=-1.0, feature_shift_high=1.0),
            _two_cluster_scan(;
                  file="scan_B.sxm",
                  n_low=3, n_high=3,
                  amp_low=0.04, amp_high=0.09,
                  feature_shift_low=-1.2, feature_shift_high=0.8))
        _write_feature_tsv(feats, rows)

        out = joinpath(tmp, "predictions.tsv")
        cmd = `$(_julia_cmd(OLD_CLI)) --features $feats --out $out --seeds 5 --first-seed 0`
        code, stdout_str, stderr_str = _run(cmd)
        @test code == 0

        header, pred_rows = _read_predictions(out)

        @testset "prediction TSV schema is frozen" begin
            @test header == ["file", "lobe", "predicted", "confidence",
                             "amplitude", "probability_1", "views_used",
                             "invalid_reason"]
        end

        @testset "higher-amplitude cluster maps to GlcNAc (1)" begin
            for row in pred_rows
                amp = parse(Float64, row["amplitude"])
                if amp > 0.07
                    @test row["predicted"] == "1"
                else
                    @test row["predicted"] == "0"
                end
            end
        end

        @testset "prediction values are in {0,1,?}" begin
            for row in pred_rows
                @test row["predicted"] in ("0", "1", "?")
            end
        end

        @testset "no benchmark/grader columns leak into the schema" begin
            for forbidden in ("sequence", "expected_N", "target_N", "control_sequence",
                              "truth", "nknnkn", "010010", "grade")
                @test !any(occursin(forbidden, lowercase(h)) for h in header)
            end
        end
    end
end

# ===========================================================================
# Hierarchical model tests. These FAIL until test/lib/hierarchical_unit_assignment.jl
# and test/build_hierarchical_unit_predictions.jl exist.
# ===========================================================================

# Holds the loaded module so multiple testsets can share it.
const _HUA = Ref{Any}(nothing)

# Load the library at top level (not inside a function) to avoid Julia's
# world-age warning from referencing a just-defined module. If the file is
# absent (red phase), the include is skipped and hierarchical_model tests
# will fail with a clear diagnostic.
if isfile(LIB)
    include(LIB)
    _HUA[] = HierarchicalUnitAssignment
end

_load_hua() = _HUA[]

@testset "hierarchical_model" begin
    @test isfile(LIB)
    @test isfile(NEW_CLI)
    HUA = _load_hua()
    @test HUA !== nothing
    HUA === nothing && return

    # ----------------------------------------------------------------------
    @testset "module surface and constants" begin
        @test isdefined(HUA, :MODEL_NAME)
        @test getproperty(HUA, :MODEL_NAME) == "hierarchical_equalprior"
        @test isdefined(HUA, :MODEL_VERSION)
        @test getproperty(HUA, :MODEL_VERSION) isa Integer
        @test isdefined(HUA, :EQUAL_PRIOR_WEIGHTS)
        w = getproperty(HUA, :EQUAL_PRIOR_WEIGHTS)
        @test length(w) == 2
        @test all(w .== 0.5)
    end

    # ----------------------------------------------------------------------
    @testset "feature-name firewall denies lobe-order/benchmark leaks" begin
        # Anything that encodes lobe order/position, sequence, expected N,
        # transitions, run lengths, class-count/top-k/occupancy controls, or
        # benchmark/truth columns must be rejected before it can become a view.
        for bad in ["centered_pos", "edge_distance_norm", "lobe_order",
                    "lobe_index", "position_along_chain", "sequence",
                    "control_sequence", "expected_N", "target_N",
                    "transition_prob", "run_length", "nknnkn", "010010",
                    "101101", "truth_label", "top_k", "class_count",
                    "occupancy"]
            @test HUA.is_forbidden_feature(bad)
        end
        # Ordinary physical features must pass.
        for good in BASE_FEATURES
            @test !HUA.is_forbidden_feature(good)
        end
        for good in ["amplitude", "sigma_parallel_nm", "integrated",
                     "bwd_neg_com_t", "split_log_skew"]
            @test !HUA.is_forbidden_feature(good)
        end
    end

    # ----------------------------------------------------------------------
    @testset "load_records is deterministic keyed by (file,lobe)" begin
        mktempdir() do tmp
            feats = joinpath(tmp, "features.tsv")
            rows = vcat(
                _two_cluster_scan(; file="s1.sxm", n_low=2, n_high=2,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0),
                _two_cluster_scan(; file="s2.sxm", n_low=2, n_high=2,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0,
                                  lobe_start=1))
            _write_feature_tsv(feats, rows)

            recs1 = HUA.load_records(feats)
            recs2 = HUA.load_records(feats)
            @test length(recs1) == length(rows)
            @test recs1 == recs2  # deterministic
            @test issorted([(r.file, r.lobe) for r in recs1])
        end
    end

    # ----------------------------------------------------------------------
    @testset "load_records rejects forbidden/benchmark columns" begin
        mktempdir() do tmp
            feats = joinpath(tmp, "features.tsv")
            rows = vcat(
                _two_cluster_scan(; file="s1.sxm", n_low=2, n_high=2,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0),
                _two_cluster_scan(; file="s2.sxm", n_low=2, n_high=2,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0))
            # Append a benchmark/truth column.
            for (i, r) in enumerate(rows)
                r["sequence"] = string(mod(i, 2))
            end
            _write_feature_tsv(feats, rows, ["sequence"])
            @test_throws Exception HUA.load_records(feats)
        end
    end

    # ----------------------------------------------------------------------
    @testset "median+MAD nuisance profiling removes affine per-scan nuisance" begin
        # Two scans with the SAME underlying cluster structure but different
        # affine feature transforms must yield the SAME normalized features
        # under median+MAD profiling.
        mktempdir() do tmp
            base = _two_cluster_scan(; file="s1.sxm", n_low=4, n_high=4,
                                     amp_low=0.05, amp_high=0.10,
                                     feature_shift_low=-1.0, feature_shift_high=1.0)
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, base)
            recs = HUA.load_records(feats)

            # Manually build an affine-shifted copy and check that median+MAD
            # normalization is invariant (up to sign for negative scale).
            X, valid = HUA.feature_matrix(recs, ["amp_prominence"])
            @test all(valid)
            prof = HUA.profile_scan(X)
            # MAD is invariant to x -> a*x + b for a != 0 (|a|*MAD); the
            # normalized value (x - median)/MAD is sign-stable for a>0.
            shift = 7.5
            scale = 3.0
            X2 = scale .* X .+ shift
            prof2 = HUA.profile_scan(X2)
            @test prof2.mad[1] ≈ scale * prof.mad[1] atol = 1e-9
            z1 = (X[:, 1] .- prof.median[1]) ./ prof.mad[1]
            z2 = (X2[:, 1] .- prof2.median[1]) ./ prof2.mad[1]
            @test z1 ≈ z2 atol = 1e-9
        end
    end

    # ----------------------------------------------------------------------
    @testset "zero-MAD scan/view features are invalidated" begin
        # A scan where every lobe has the SAME feature value has MAD = 0.
        # That feature must be invalidated for the entire scan/view under
        # per-scan median+MAD normalization.
        mktempdir() do tmp
            deg = _degenerate_scan(; file="deg.sxm", n=5, amp=0.08,
                                   feature_value=2.0)
            norm = _two_cluster_scan(; file="norm.sxm", n_low=3, n_high=3,
                                     amp_low=0.05, amp_high=0.10,
                                     feature_shift_low=-1.0, feature_shift_high=1.0)
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, vcat(deg, norm))
            recs = HUA.load_records(feats)
            X, _ = HUA.feature_matrix(recs, ["amp_prominence"])
            Z = HUA.normalize_per_scan(recs, X, ["amp_prominence"])
            # All lobes of the degenerate scan must be NaN after normalization.
            deg_idx = findall(r -> r.file == "deg.sxm", recs)
            @test all(isnan, Z[deg_idx, 1])
            # The normal scan is unaffected.
            norm_idx = findall(r -> r.file == "norm.sxm", recs)
            @test all(isfinite, Z[norm_idx, 1])
        end
    end

    # ----------------------------------------------------------------------
    @testset "two-component EM recovers two shared emissions" begin
        # A clearly bimodal 1D feature: lobes at -3 and +3.
        X = vcat(fill(-3.0, 30), fill(3.0, 30))
        fit = HUA.fit_em_two_component(X, 0)
        @test fit.converged
        @test fit.monotone
        means = sort([fit.means[1], fit.means[2]])
        @test means[1] <= -1.0
        @test means[2] >= 1.0
        # Priors remain EXACTLY 0.5/0.5.
        @test fit.weights ≈ [0.5, 0.5]
    end

    # ----------------------------------------------------------------------
    @testset "EM priors stay 0.5/0.5 regardless of generated occupancy" begin
        # Deliberately unequal occupancy: 50 lobes in one cluster, 4 in the
        # other. Equal priors must NOT follow the empirical class proportions.
        X = vcat(fill(-3.0, 50), fill(3.0, 4))
        fit = HUA.fit_em_two_component(X, 0)
        @test fit.converged
        @test fit.weights ≈ [0.5, 0.5]
        # And the means must still be recoverable.
        means = sort([fit.means[1], fit.means[2]])
        @test means[1] <= -1.0
        @test means[2] >= 1.0
    end

    # ----------------------------------------------------------------------
    @testset "EM log-likelihood is monotone non-decreasing" begin
        # Pathological input: tight clusters + a couple of outliers.
        X = vcat(fill(-2.0, 20), fill(2.0, 20), [-10.0, 10.0])
        fit = HUA.fit_em_two_component(X, 0)
        @test fit.monotone
        @test issorted(fit.loglik_trace)
    end

    # ----------------------------------------------------------------------
    @testset "EM covariance floor prevents variance collapse" begin
        # Identical values in each cluster (zero within-cluster variance).
        # The floor must keep variances strictly positive.
        X = vcat(fill(-1.0, 10), fill(1.0, 10))
        fit = HUA.fit_em_two_component(X, 0; cov_floor=1e-4)
        @test all(fit.vars .> 0.0)
        @test all(fit.vars .>= 1e-4 - 1e-12)
    end

    # ----------------------------------------------------------------------
    @testset "EM is deterministic across repeated calls (multi-start)" begin
        X = vcat(fill(-2.5, 25), fill(2.5, 25))
        f1 = HUA.fit_em_two_component(X, 0)
        f2 = HUA.fit_em_two_component(X, 0)
        @test f1.means ≈ f2.means
        @test f1.vars ≈ f2.vars
        @test f1.loglik ≈ f2.loglik atol = 1e-9
    end

    # ----------------------------------------------------------------------
    @testset "one-component companion likelihood is finite and below 2-comp" begin
        X = vcat(fill(-3.0, 30), fill(3.0, 30))
        fit2 = HUA.fit_em_two_component(X, 0)
        fit1 = HUA.fit_one_component(X)
        @test isfinite(fit1.loglik)
        @test isfinite(fit2.loglik)
        @test fit2.loglik >= fit1.loglik - 1e-9
    end

    # ----------------------------------------------------------------------
    @testset "physical mapping: higher-amplitude cluster is GlcNAc (1)" begin
        # Build records where the high-amplitude lobes are obvious, then
        # verify the cluster id that maps to label 1 is the one whose members
        # have higher mean raw amplitude.
        mktempdir() do tmp
            rows = _two_cluster_scan(; file="s1.sxm", n_low=5, n_high=5,
                                     amp_low=0.04, amp_high=0.10,
                                     feature_shift_low=-1.0, feature_shift_high=1.0)
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            recs = HUA.load_records(feats)

            X, _ = HUA.feature_matrix(recs, ["amp_prominence"])
            fit = HUA.fit_em_two_component(X[:, 1], 0)
            resp = HUA.responsibilities(fit, X[:, 1])
            assignments = [argmax(resp[i, :]) for i in 1:length(recs)]
            high_cluster = HUA.physical_high_amplitude_cluster(recs, assignments)
            # The high-amp cluster's mean amplitude must exceed the other's.
            amps = [r.amplitude for r in recs]
            mean_high = mean(amps[assignments .== high_cluster])
            other = 3 - high_cluster  # cluster ids are 1 or 2
            mean_other = mean(amps[assignments .== other])
            @test mean_high > mean_other
        end
    end

    # ----------------------------------------------------------------------
    @testset "affine per-scan nuisance leaves assignments invariant" begin
        mktempdir() do tmp
            # Baseline: two scans with the same cluster structure.
            base_rows = vcat(
                _two_cluster_scan(; file="s1.sxm", n_low=4, n_high=4,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0),
                _two_cluster_scan(; file="s2.sxm", n_low=4, n_high=4,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0))

            # Affine-nuanced version: scan s2 has feature values shifted+scaled.
            affine_rows = vcat(
                _two_cluster_scan(; file="s1.sxm", n_low=4, n_high=4,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0),
                _affine_scan(; file="s2.sxm", n_low=4, n_high=4,
                             amp_low=0.05, amp_high=0.10,
                             f_low=-1.0, f_high=1.0,
                             scale=4.0, shift=10.0))

            f1 = joinpath(tmp, "features_base.tsv")
            f2 = joinpath(tmp, "features_affine.tsv")
            _write_feature_tsv(f1, base_rows)
            _write_feature_tsv(f2, affine_rows)

            r1 = HUA.load_records(f1)
            r2 = HUA.load_records(f2)

            # Per-scan median+MAD normalization must remove the affine nuisance:
            # the normalized feature matrices for s2 must coincide (because
            # (a*x+b - median)/(MAD) is invariant in a>0). Compare the s2 block.
            feats = ["amp_prominence"]
            X1, _ = HUA.feature_matrix(r1, feats)
            X2, _ = HUA.feature_matrix(r2, feats)
            Z1 = HUA.normalize_per_scan(r1, X1, feats)
            Z2 = HUA.normalize_per_scan(r2, X2, feats)
            # Restrict to the s2 lobes.
            s2_idx1 = findall(r -> r.file == "s2.sxm", r1)
            s2_idx2 = findall(r -> r.file == "s2.sxm", r2)
            @test length(s2_idx1) == length(s2_idx2)
            # Floating-point round-off from scale*x+shift propagates through
            # median+MAD; compare element-wise with a tolerance well above
            # Float64 round-off (~1e-7 here) but far below the normalized
            # scale (~1). `≈` on matrices uses Frobenius norm, so use
            # element-wise `all(.≈)` to compare per-entry.
            @test all(isapprox.((Z1[s2_idx1, :] .- Z2[s2_idx2, :]), 0.0; atol=1e-6))
        end
    end

    # ----------------------------------------------------------------------
    @testset "assignments are invariant under deterministic row permutation" begin
        mktempdir() do tmp
            rows = vcat(
                _two_cluster_scan(; file="s1.sxm", n_low=4, n_high=4,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0),
                _two_cluster_scan(; file="s2.sxm", n_low=4, n_high=4,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0))
            f_sorted = joinpath(tmp, "features_sorted.tsv")
            f_shuffled = joinpath(tmp, "features_shuffled.tsv")
            _write_feature_tsv(f_sorted, rows)

            # Deterministic non-trivial permutation (seed-fixed).
            rng = MersenneTwister(12345)
            perm = shuffle(rng, collect(1:length(rows)))
            _write_feature_tsv(f_shuffled, rows[perm])

            out1 = joinpath(tmp, "pred_sorted.tsv")
            out2 = joinpath(tmp, "pred_shuffled.tsv")
            c1 = `$(_julia_cmd(NEW_CLI)) --features $f_sorted --out $out1 --first-seed 0`
            c2 = `$(_julia_cmd(NEW_CLI)) --features $f_shuffled --out $out2 --first-seed 0`
            @test _run(c1)[1] == 0
            @test _run(c2)[1] == 0

            h1, r1 = _read_predictions(out1)
            h2, r2 = _read_predictions(out2)
            @test h1 == h2
            # Join on (file, lobe) and confirm identical predictions.
            d1 = Dict((row["file"], parse(Int, row["lobe"])) => row["predicted"] for row in r1)
            d2 = Dict((row["file"], parse(Int, row["lobe"])) => row["predicted"] for row in r2)
            @test Set(keys(d1)) == Set(keys(d2))
            for k in keys(d1)
                @test d1[k] == d2[k]
            end
        end
    end

    # ----------------------------------------------------------------------
    @testset "CLI preserves prediction schema + model/version/provenance fields" begin
        mktempdir() do tmp
            rows = vcat(
                _two_cluster_scan(; file="s1.sxm", n_low=3, n_high=3,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0),
                _two_cluster_scan(; file="s2.sxm", n_low=3, n_high=3,
                                  amp_low=0.04, amp_high=0.09,
                                  feature_shift_low=-1.1, feature_shift_high=0.9))
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out = joinpath(tmp, "pred.tsv")
            cmd = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out --first-seed 0`
            code, stdout_str, stderr_str = _run(cmd)
            @test code == 0
            header, pred_rows = _read_predictions(out)

            # Existing frozen schema must be preserved as a subset.
            for col in ("file", "lobe", "predicted", "confidence", "amplitude",
                        "probability_1", "views_used", "invalid_reason")
                @test col in header
            end
            # New model/version/provenance fields must be present.
            @test "model" in header
            @test "model_version" in header
            @test "provenance_sha256" in header

            # Higher-amplitude physical convention.
            for row in pred_rows
                amp = parse(Float64, row["amplitude"])
                if row["predicted"] in ("0", "1")
                    if amp > 0.07
                        @test row["predicted"] == "1"
                    else
                        @test row["predicted"] == "0"
                    end
                end
            end

            # Model identity + provenance populated.
            @test all(row["model"] == "hierarchical_equalprior" for row in pred_rows)
            @test all(!isempty(row["provenance_sha256"]) for row in pred_rows)
            @test all(length(row["provenance_sha256"]) == 64 for row in pred_rows)

            # No benchmark/truth content anywhere in the output.
            blob = read(out, String)
            for forbidden in ("sequence", "expected_N", "target_N", "control_sequence",
                              "nknnkn", "010010", "truth")
                @test !occursin(forbidden, lowercase(blob))
            end
        end
    end

    # ----------------------------------------------------------------------
    @testset "CLI --help exits 0 and documents label-free contract" begin
        code, out, err = _run(`$(_julia_cmd(NEW_CLI)) --help`)
        @test code == 0
        joined = lowercase(out * "\n" * err)
        @test occursin("label-free", joined)
        @test occursin("equal prior", joined)
    end

    # ----------------------------------------------------------------------
    @testset "CLI rejects forbidden benchmark/lobe-order flags" begin
        mktempdir() do tmp
            rows = _two_cluster_scan(; file="s1.sxm", n_low=2, n_high=2,
                                     amp_low=0.05, amp_high=0.10,
                                     feature_shift_low=-1.0, feature_shift_high=1.0)
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out = joinpath(tmp, "pred.tsv")
            for flag in ["--truth", "--control-sequence", "--manifest",
                         "--expected-N", "--expected-n", "--target-N",
                         "--transitions", "--transition", "--run-length",
                         "--top-k", "--class-count", "--occupancy",
                         "--lobe-index", "--lobe-order", "--full145",
                         "--control", "--benchmark-manifest"]
                cmd = `$(_julia_cmd(NEW_CLI)) $flag x --features $feats --out $out`
                code, = _run(cmd)
                @test code != 0
            end
        end
    end

    # ----------------------------------------------------------------------
    @testset "CLI abstains on one-component fixture" begin
        mktempdir() do tmp
            # Single feature population across all scans -> 2-comp model should
            # collapse or be flagged as one-component and abstain.
            rows = vcat(
                _one_component_scan(; file="s1.sxm", n=8, amp=0.07, feature_value=0.5),
                _one_component_scan(; file="s2.sxm", n=8, amp=0.07, feature_value=0.5))
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out = joinpath(tmp, "pred.tsv")
            cmd = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out --first-seed 0`
            code, = _run(cmd)
            @test code == 0
            header, pred_rows = _read_predictions(out)
            # The invalid_reason column must mention one-component evidence,
            # and at least one row must abstain with "?".
            reasons = join([row["invalid_reason"] for row in pred_rows], " ")
            @test occursin("one_component", lowercase(reasons))
            @test any(row["predicted"] == "?" for row in pred_rows)
        end
    end

    # ----------------------------------------------------------------------
    @testset "CLI abstains on degenerate-covariance fixture" begin
        mktempdir() do tmp
            # Every feature value identical within each scan: MAD = 0 across
            # the board. Every view must be invalidated; all rows abstain.
            rows = vcat(
                _degenerate_scan(; file="s1.sxm", n=6, amp=0.07, feature_value=1.0),
                _degenerate_scan(; file="s2.sxm", n=6, amp=0.07, feature_value=1.0))
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out = joinpath(tmp, "pred.tsv")
            cmd = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out --first-seed 0`
            code, = _run(cmd)
            @test code == 0
            header, pred_rows = _read_predictions(out)
            @test all(row["predicted"] == "?" for row in pred_rows)
            reasons = join([row["invalid_reason"] for row in pred_rows], " ")
            @test occursin("degenerate", lowercase(reasons))
        end
    end

    # ----------------------------------------------------------------------
    @testset "CLI reports missing-view abstention" begin
        mktempdir() do tmp
            # Build a feature TSV with NaN/NA values on every feature for one
            # lobe: no view can produce a responsibility -> missing-view
            # abstention.
            rows = _two_cluster_scan(; file="s1.sxm", n_low=3, n_high=3,
                                     amp_low=0.05, amp_high=0.10,
                                     feature_shift_low=-1.0, feature_shift_high=1.0)
            # Blank out feature columns for lobe 1.
            for col in BASE_FEATURES
                rows[1][col] = "NA"
            end
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out = joinpath(tmp, "pred.tsv")
            cmd = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out --first-seed 0`
            code, = _run(cmd)
            @test code == 0
            header, pred_rows = _read_predictions(out)
            blank = first(row for row in pred_rows if parse(Int, row["lobe"]) == 1)
            @test blank["predicted"] == "?"
            @test occursin("missing", lowercase(blank["invalid_reason"]))
        end
    end

    # ----------------------------------------------------------------------
    @testset "CLI handles unstable-EM fixture with abstention" begin
        mktempdir() do tmp
            # Two perfectly interleaved identical single-point clusters with
            # extreme outliers - designed to make naive EM unstable. The model
            # must either converge deterministically (multi-start) or abstain
            # with an explicit reason rather than fabricate labels.
            rows = Dict{String,Any}[]
            for lobe in 1:4
                push!(rows, Dict(
                    "file" => "s1.sxm", "N" => 6, "lobe" => lobe,
                    "amplitude" => "0.07", "x_nm" => "0.0", "y_nm" => "0.0",
                    "sigma_parallel_nm" => "0.40", "sigma_perp_nm" => "0.45",
                    "integrated" => "0.0126",
                    "amp_prominence" => "0.0",
                    "amp_neighbor_ratio" => "0.0",
                    "integrated_prominence" => "0.0",
                    "amp_rel" => "0.0"))
            end
            for lobe in 5:6
                push!(rows, Dict(
                    "file" => "s1.sxm", "N" => 6, "lobe" => lobe,
                    "amplitude" => "1e6", "x_nm" => "0.0", "y_nm" => "0.0",
                    "sigma_parallel_nm" => "0.40", "sigma_perp_nm" => "0.45",
                    "integrated" => "1e6",
                    "amp_prominence" => "1e6",
                    "amp_neighbor_ratio" => "1e6",
                    "integrated_prominence" => "1e6",
                    "amp_rel" => "1e6"))
            end
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out = joinpath(tmp, "pred.tsv")
            cmd = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out --first-seed 0`
            code, stdout_str, stderr_str = _run(cmd)
            @test code == 0
            header, pred_rows = _read_predictions(out)
            # Either abstains on every lobe OR produces {0,1,?} values only.
            for row in pred_rows
                @test row["predicted"] in ("0", "1", "?")
            end
        end
    end

    # ----------------------------------------------------------------------
    @testset "deterministic output: repeated runs produce byte-identical TSVs" begin
        mktempdir() do tmp
            rows = vcat(
                _two_cluster_scan(; file="s1.sxm", n_low=3, n_high=3,
                                  amp_low=0.05, amp_high=0.10,
                                  feature_shift_low=-1.0, feature_shift_high=1.0),
                _two_cluster_scan(; file="s2.sxm", n_low=3, n_high=3,
                                  amp_low=0.04, amp_high=0.09,
                                  feature_shift_low=-1.1, feature_shift_high=0.9))
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out1 = joinpath(tmp, "pred1.tsv")
            out2 = joinpath(tmp, "pred2.tsv")
            c1 = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out1 --first-seed 0`
            c2 = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out2 --first-seed 0`
            @test _run(c1)[1] == 0
            @test _run(c2)[1] == 0
            @test read(out1, String) == read(out2, String)
        end
    end

    # ----------------------------------------------------------------------
    # T4-correction B1: bare ordinal / lobe-position feature names must be
    # denied at the parse boundary, not only the compound spellings.
    # ----------------------------------------------------------------------
    @testset "firewall denies bare ordinal/lobe-position names (B1)" begin
        for bad in ["order", "Order", "ORDER",
                    "index", "Index", "INDEX",
                    "rank", "Rank", "RANK",
                    "lobe_pos", "Lobe_Pos", "LOBE_POS",
                    "lobe-pos", "LOBE-POS", "lobe.pos", "lobe pos",
                    "lobepos", "LOBEPOS",
                    "o.r.d.e.r", "i/n/d/e/x", "r-a-n-k",
                    "l.o.b.e/p.o.s", "οrder", "lobe–pos",
                    "lobe／pos", "ｏｒｄｅｒ"]
            @test HUA.is_forbidden_feature(bad)
        end
        # Legitimate physical features must still pass the firewall.
        for good in [BASE_FEATURES;
                     ["amplitude", "sigma_parallel_nm", "sigma_perp_nm",
                      "integrated", "bwd_neg_com_t", "bwd_neg_diag45",
                      "bwd_neg_diag135", "split_log_skew", "x_nm", "y_nm"]]
            @test !HUA.is_forbidden_feature(good)
        end
    end

    @testset "load_records rejects bare ordinal columns (B1)" begin
        mktempdir() do tmp
            for col in ("order", "index", "rank", "lobe_pos")
                feats = joinpath(tmp, "features_$col.tsv")
                rows = _two_cluster_scan(; file="s1.sxm", n_low=2, n_high=2,
                                          amp_low=0.05, amp_high=0.10,
                                          feature_shift_low=-1.0,
                                          feature_shift_high=1.0)
                for (i, r) in enumerate(rows)
                    r[col] = string(i)
                end
                _write_feature_tsv(feats, rows, [col])
                @test_throws Exception HUA.load_records(feats)
            end
        end
    end

    @testset "load_records rejects punctuation/Unicode aliases (B1)" begin
        mktempdir() do tmp
            for (k, col) in enumerate(("o.r.d.e.r", "οrder",
                                       "lobe–pos", "lobe／pos"))
                feats = joinpath(tmp, "variant_$k.tsv")
                rows = _two_cluster_scan(; file="s1.sxm", n_low=2, n_high=2,
                                          amp_low=0.05, amp_high=0.10,
                                          feature_shift_low=-1.0,
                                          feature_shift_high=1.0)
                for (i, r) in enumerate(rows)
                    r[col] = string(i)
                end
                _write_feature_tsv(feats, rows, [col])
                @test_throws Exception HUA.load_records(feats)
            end
        end
    end

    @testset "CLI rejects --view carrying a bare ordinal feature (B1)" begin
        mktempdir() do tmp
            rows = _two_cluster_scan(; file="s1.sxm", n_low=2, n_high=2,
                                      amp_low=0.05, amp_high=0.10,
                                      feature_shift_low=-1.0,
                                      feature_shift_high=1.0)
            for (i, r) in enumerate(rows)
                r["order"] = string(i)
                r["index"] = string(i)
                r["rank"] = string(i)
                r["lobe_pos"] = string(i)
            end
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows, ["order", "index", "rank", "lobe_pos"])
            out = joinpath(tmp, "pred.tsv")
            # `order`, `index`, `rank`, `lobe_pos` are available columns, so
            # only the firewall can reject them: a nonzero exit proves the
            # guard fires at the --view re-check, not an availability error.
            for leak in ("order", "index", "rank", "lobe_pos")
                cmd = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out --view leak=$leak`
                code, = _run(cmd)
                @test code != 0
            end
        end
    end

    @testset "CLI rejects punctuation/Unicode ordinal aliases (B1)" begin
        mktempdir() do tmp
            for (k, leak) in enumerate(("o.r.d.e.r", "οrder",
                                        "lobe–pos", "lobe／pos"))
                rows = _two_cluster_scan(; file="s1.sxm", n_low=2, n_high=2,
                                          amp_low=0.05, amp_high=0.10,
                                          feature_shift_low=-1.0,
                                          feature_shift_high=1.0)
                for (i, r) in enumerate(rows)
                    r[leak] = string(i)
                end
                feats = joinpath(tmp, "variant_$k.tsv")
                out = joinpath(tmp, "variant_$k-pred.tsv")
                _write_feature_tsv(feats, rows, [leak])
                cmd = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out --view leak=$leak`
                code, = _run(cmd)
                @test code != 0
            end
        end
    end

    # ----------------------------------------------------------------------
    # T4-correction B2: --first-seed must reach the EM multistart seam and
    # change the fit (previously parsed + hashed into provenance but silently
    # dropped before fit_em_two_component).
    # ----------------------------------------------------------------------
    @testset "fit_view forwards first_seed to the EM multistart (B2)" begin
        mktempdir() do tmp
            rows = _seed_sensitive_scan(; file="s1.sxm")
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            recs = HUA.load_records(feats)
            vf0 = HUA.fit_view(recs, ["amp_prominence"], "v"; first_seed=0)
            vf3 = HUA.fit_view(recs, ["amp_prominence"], "v"; first_seed=3)
            m0 = sort([vf0.fit2.means[1], vf0.fit2.means[2]])
            m3 = sort([vf3.fit2.means[1], vf3.fit2.means[2]])
            @test !(m0 ≈ m3)
        end
    end

    @testset "CLI --first-seed changes the fit, not only provenance (B2)" begin
        mktempdir() do tmp
            rows = _seed_sensitive_scan(; file="s1.sxm")
            feats = joinpath(tmp, "features.tsv")
            _write_feature_tsv(feats, rows)
            out0  = joinpath(tmp, "pred_s0.tsv")
            out3  = joinpath(tmp, "pred_s3.tsv")
            out0b = joinpath(tmp, "pred_s0b.tsv")
            c0  = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out0  --view v=amp_prominence --first-seed 0`
            c3  = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out3  --view v=amp_prominence --first-seed 3`
            c0b = `$(_julia_cmd(NEW_CLI)) --features $feats --out $out0b --view v=amp_prominence --first-seed 0`
            @test _run(c0)[1] == 0
            @test _run(c3)[1] == 0
            @test _run(c0b)[1] == 0
            # Same seed is byte-deterministic.
            @test read(out0, String) == read(out0b, String)
            h0, r0 = _read_predictions(out0)
            h3, r3 = _read_predictions(out3)
            # Provenance legitimately varies with the seed (faithful binding).
            @test [row["provenance_sha256"] for row in r0] !=
                  [row["provenance_sha256"] for row in r3]
            # The fit must change too: the probability_1 column differs across
            # seeds, not only the provenance column.
            @test [row["probability_1"] for row in r0] !=
                  [row["probability_1"] for row in r3]
        end
    end

    # ----------------------------------------------------------------------
    # T4-correction B3: a monotone-likelihood violation must return a
    # self-consistent last-accepted (means, vars, loglik) record, not a mix of
    # post-violation params with a pre-violation loglik.
    # ----------------------------------------------------------------------
    @testset "EM monotone-violation returns self-consistent params (B3)" begin
        # Forty identical points: within-cluster variance is zero, so the
        # M-step variance is clamped to the cov floor. With means initialized
        # at (-1, 1) and cov_floor=10, the first M-step collapses both means
        # to 0 with floor variance, which DECREASES the total log-likelihood
        # below the initialization value -> _run_em sets monotone=false and
        # breaks. The returned (means, vars) must correspond to trace[end]
        # (the last monotone-accepted record), not the post-violation M-step.
        Z = reshape(zeros(40), :, 1)
        means0 = reshape([-1.0, 1.0], 2, 1)
        vars0  = reshape([1.0, 1.0], 2, 1)
        means, vars, trace, conv, monot =
            HUA._run_em(Z, means0, vars0, HUA.EQUAL_PRIOR_WEIGHTS, 10.0, 200, 1e-6)
        @test !monot
        @test trace[end] ≈ HUA._total_log_likelihood(Z, means, vars,
                                                     HUA.EQUAL_PRIOR_WEIGHTS) atol=1e-9
    end

    @testset "EM fit.loglik matches its means/vars across fixtures (B3)" begin
        # Self-consistency invariant for every public fit: fit.loglik must
        # equal the log-likelihood recomputed from fit.means/vars/weights.
        for X in (vcat(fill(-3.0, 30), fill(3.0, 30)),
                  vcat(fill(-2.0, 20), fill(2.0, 20), [-10.0, 10.0]),
                  vcat(fill(0.0, 40), fill(1.0, 40), [-100.0, 100.0]),
                  collect(range(-2.0, 2.0; length=40)),
                  vcat(fill(0.0, 30), collect(0.5:0.1:2.0)))
            fit = HUA.fit_em_two_component(X, 0)
            @test fit.loglik ≈ HUA.log_likelihood_under(fit, X) atol=1e-9
        end
    end
end
