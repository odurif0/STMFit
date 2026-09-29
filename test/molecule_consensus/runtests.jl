# Synthetic tests for repeated-scan grouping and molecule-level count consensus.
#   julia --project=. test/molecule_consensus/runtests.jl
using Test, Dates, Random, Statistics, TOML
include(joinpath(@__DIR__, "..", "lib", "molecule_consensus.jl"))
using .MoleculeConsensus

const ST = MoleculeConsensus.ConsensusSettings(3.0, 0.04, 0.32, 1.6, 4, 500, 0.5, 6.0, 3, 0.3, "plane+rows", "Z")

frame(name; t=0, off=(10.0, -5.0), ext=(6.0, 6.0), ang=0.0, px=256, session="17.08.2024") =
    ScanFrame(name, DateTime(2024, 8, 17, 12, 0, 0) + Minute(t), session, off, ext, ang, px, px)

# Render absolute-frame Gaussian lobes into a scan's pixel grid (rows = y down).
function render(f::ScanFrame, lobes; sigma=0.25, noise=0.0, rng=MersenneTwister(1))
    z = zeros(f.height, f.width)
    for r in 1:f.height, c in 1:f.width
        x = (c - 1) * f.extent[1] / (f.width - 1); y = (r - 1) * f.extent[2] / (f.height - 1)
        X, Y = absolute_xy(f, x, y)
        z[r, c] = sum(exp(-((X - p[1])^2 + (Y - p[2])^2) / (2sigma^2)) for p in lobes)
    end
    return z .+ noise .* randn(rng, size(z))
end

zigzag(c; n=6, sp=0.6, amp=0.15, ang=0.4) =
    [(c[1] + (i - (n + 1) / 2) * sp * cos(ang) - (isodd(i) ? amp : -amp) * sin(ang),
      c[2] + (i - (n + 1) / 2) * sp * sin(ang) + (isodd(i) ? amp : -amp) * cos(ang)) for i in 1:n]

@testset "absolute frame round trip and rotation consistency" begin
    for ang in (-90.0, -70.0, 0.0, 30.0, 45.0)
        f = frame("a"; ang, off=(3.2, -1.7), ext=(7.0, 7.0))
        for (x, y) in ((0.0, 0.0), (1.3, 5.2), (6.9, 0.4))
            X, Y = absolute_xy(f, x, y); x2, y2 = image_xy(f, X, Y)
            @test isapprox(x2, x; atol=1e-12) && isapprox(y2, y; atol=1e-12)
        end
    end
    # One absolute molecule seen through two rotated frames maps back to itself.
    lobes = zigzag((10.5, -4.6))
    fa = frame("a"; ang=40.0); fb = frame("b"; ang=-60.0, off=(10.3, -4.9))
    for p in lobes
        xa, ya = image_xy(fa, p...); xb, yb = image_xy(fb, p...)
        @test all(isapprox.(absolute_xy(fa, xa, ya), p; atol=1e-12))
        @test all(isapprox.(absolute_xy(fb, xb, yb), p; atol=1e-12))
    end
end

@testset "registration recovers drift and rejects unrelated structure" begin
    lobes = zigzag((10.2, -5.1))
    drift = (0.36, -0.52)
    fa = frame("a"; ang=0.0); fb = frame("b"; t=2, ang=30.0, off=(10.1, -5.0))
    za = render(fa, lobes; noise=0.02)
    zb = render(fb, [(p[1] + drift[1], p[2] + drift[2]) for p in lobes]; noise=0.02, rng=MersenneTwister(2))
    sigma = ST.highpass_sigma_nm / ST.grid_step_nm
    c = (mean(first.(lobes)), mean(last.(lobes)))
    A = highpass(resample_absolute(za, fa, c, ST), sigma)
    B = highpass(resample_absolute(zb, fb, c, ST), sigma)
    reg = register_grids(A, B, ST)
    @test reg.ncc > 0.8
    @test abs(reg.tx - drift[1]) <= ST.grid_step_nm + 1e-9
    @test abs(reg.ty - drift[2]) <= ST.grid_step_nm + 1e-9
    other = [(10.2 + 0.7 * cos(k), -5.1 + 1.1 * sin(2k)) for k in 1:5]
    zc = render(fb, other; noise=0.02, rng=MersenneTwister(3))
    C = highpass(resample_absolute(zc, fb, c, ST), sigma)
    @test register_grids(A, C, ST).ncc < ST.ncc_min
    # Missing pixels stay missing and never become observations.
    zm = copy(za); zm[1:40, :] .= NaN
    M = resample_absolute(zm, fa, c, ST)
    @test any(isnan, M) && all(isfinite, filter(isfinite, M))
end

function synthetic_cohort(counts; linked=trues(length(counts) - 1), drift=(0.1, 0.0), outside=Int[])
    files = ["s$(i).sxm" for i in eachindex(counts)]
    frames = Dict(f => frame(f; t=i) for (i, f) in enumerate(files))
    lobes = Dict(f => zigzag((10.0, -5.0); n=counts[i]) for (i, f) in enumerate(files))
    for i in outside  # molecule footprint partially outside this scan's frame
        frames[files[i]] = frame(files[i]; t=i, off=(12.5, -5.0))
    end
    links = [(a=files[k], b=files[k+1], center_distance_nm=0.1, ncc=linked[k] ? 0.9 : 0.2,
              tx_nm=drift[1], ty_nm=drift[2], linked=linked[k], reason=linked[k] ? "linked" : "low_correlation")
             for k in 1:length(files)-1]
    rows = consensus_counts(files, links, Dict(zip(files, counts)), lobes, frames, ST)
    return Dict(r["file"] => r for r in rows)
end

@testset "strict-majority consensus rules" begin
    r = synthetic_cohort([6, 6, 5, 6])
    @test r["s3.sxm"]["N_final"] == 6 && r["s3.sxm"]["rule"] == "consensus_applied"
    @test r["s3.sxm"]["reference_scan"] in ("s2.sxm", "s4.sxm")
    @test all(r["s$i.sxm"]["rule"] == "agrees" for i in (1, 2, 4))
    @test r["s1.sxm"]["agreement"] == 0.75
    r = synthetic_cohort([6, 7, 7, 6])  # tie: no strict majority, nothing changes
    @test all(v["N_final"] == v["N_scan"] && v["rule"] == "no_strict_majority" for v in values(r))
    r = synthetic_cohort([6, 5])  # two scans are not enough evidence
    @test all(v["rule"] == "track_too_short" && v["N_final"] == v["N_scan"] for v in values(r))
    r = synthetic_cohort([6, 6, 5, 6, 6]; linked=[true, false, true, true])
    @test r["s1.sxm"]["track"] != r["s3.sxm"]["track"]
    @test r["s3.sxm"]["N_final"] == 6 && r["s1.sxm"]["rule"] == "track_too_short"
    r = synthetic_cohort([6, 6, 5, 6]; outside=[3])
    @test r["s3.sxm"]["rule"] == "footprint_outside_frame" && r["s3.sxm"]["N_final"] == 5
    # Majority values other than any benchmark expectation are applied the same way.
    r = synthetic_cohort([9, 9, 8, 9])
    @test r["s3.sxm"]["N_final"] == 9
end

@testset "settings validation" begin
    cfg = Dict("model" => Dict("window_half_nm" => 3.0, "grid_step_nm" => 0.04, "highpass_sigma_nm" => 0.32,
                               "max_shift_nm" => 1.6, "coarse_step_px" => 4),
               "selection" => Dict("ncc_min" => 0.5, "min_overlap_px" => 500, "max_center_distance_nm" => 6.0,
                                   "min_track_scans" => 3, "majority" => "strict", "frame_margin_nm" => 0.3),
               "preprocessing" => Dict("channel" => "Z", "flatten" => "plane+rows"))
    @test load_consensus_settings(cfg).ncc_min == 0.5
    cfg["selection"]["majority"] = "plurality"
    @test_throws ErrorException load_consensus_settings(cfg)
    cfg["selection"]["majority"] = "strict"; cfg["selection"]["min_track_scans"] = 1
    @test_throws ErrorException load_consensus_settings(cfg)
    repo = dirname(dirname(@__DIR__))
    @test load_consensus_settings(TOML.parsefile(joinpath(repo, "config", "molecule_consensus.toml"))).min_track_scans == 3
end

module MCRun
include(joinpath(@__DIR__, "..", "run_molecule_consensus_chitosan.jl"))
end

@testset "consensus counts are kept only when the fixed-N refit succeeds" begin
    mktempdir() do dir
        proposed = joinpath(dir, "consensus_summary.tsv")
        MCRun.MoleculeConsensusRun.write_table(proposed,
            ["filepath", "status", "N_selected", "N_scan", "count_rule"],
            [Dict("filepath" => "a.sxm", "status" => "ok", "N_selected" => "6", "N_scan" => "5", "count_rule" => "consensus_applied"),
             Dict("filepath" => "b.sxm", "status" => "ok", "N_selected" => "6", "N_scan" => "3", "count_rule" => "consensus_applied"),
             Dict("filepath" => "c.sxm", "status" => "ok", "N_selected" => "6", "N_scan" => "6", "count_rule" => "agrees")])
        calls = Vector{Vector{String}}()
        fake(out, name, script, args; threads=1) = begin
            push!(calls, args)
            outpath = args[findfirst(==("--out"), args) + 1]
            open(outpath, "w") do io  # only a.sxm fits at the consensus count
                println(io, "file\tlobe"); for k in 1:6; println(io, "a.sxm\t$k"); end
            end
        end
        final = MCRun.MoleculeConsensusRun.check_consensus_refits(dir, dir, "count.toml", proposed, fake)
        _, rows = MCRun.MoleculeConsensusRun.read_table(final)
        byfile = Dict(r["filepath"] => r for r in rows)
        @test byfile["a.sxm"]["N_selected"] == "6" && byfile["a.sxm"]["count_rule"] == "consensus_applied"
        @test byfile["b.sxm"]["N_selected"] == "3" && byfile["b.sxm"]["count_rule"] == "consensus_fit_failed"
        @test byfile["c.sxm"]["N_selected"] == "6" && byfile["c.sxm"]["count_rule"] == "agrees"
        @test length(calls) == 1 && calls[1][findfirst(==("--files"), calls[1]) + 1] == "a.sxm,b.sxm"
    end
end
