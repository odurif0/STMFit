#!/usr/bin/env julia
# Synthetic/metadata tests only. No raw STM fit and no benchmark labels.
# julia --startup-file=no --threads=1 --project=. test/test_counting_variable_projection.jl
using Test, LinearAlgebra, Random, TOML
include(joinpath(@__DIR__, "diagnose_counting_variable_projection.jl"))
const CLI = CountingVariableProjectionCLI
const VP = CLI.CountingVariableProjection
const G = CLI.G
BLAS.set_num_threads(1)

const BASE_OPTIONS = Dict{String,Any}(
    "outer_maxeval" => 80, "outer_maxtime_s" => 2.0, "outer_xtol_rel" => 1e-7, "outer_ftol_rel" => 1e-10,
    "linear_maxiter" => 1000, "linear_kkt_atol" => 1e-11, "linear_kkt_rtol" => 1e-10,
    "linear_bound_atol" => 1e-12, "linear_svd_rtol" => 1e-12, "mapping_atol" => 1e-11, "mapping_rtol" => 1e-11,
    "native_elliptical_maxiter" => 20)
options(; kwargs...) = VP.read_options(merge(BASE_OPTIONS, Dict(string(k) => v for (k,v) in kwargs)))
const OPTS = options()

function cfg_fixture(; circular=false, tilted=true, spacing="free", shared=0)
    raw = TOML.parsefile(joinpath(@__DIR__, "..", "config", "chitosan.toml"))
    _, cfg, _ = CLI.Extractor._configs(raw["model"], raw["preprocessing"], "unused")
    cfg.chain_circular_sigmas = circular
    cfg.chain_tilted_baseline = tilted
    cfg.chain_spacing_model = spacing
    cfg.shared_sigma_types = shared
    cfg.skip_global = true
    cfg.max_iter = 20
    cfg.multistart = 1
    return cfg
end

# Exhaustive face enumeration provides a separate tiny-problem oracle.
function box_oracle(A, b, lo, hi)
    best = Inf
    for state in Iterators.product(fill((-1, 0, 1), length(lo))...)
        x = zeros(length(lo))
        free = findall(==(0), collect(state))
        for i in eachindex(lo); x[i] = state[i] == -1 ? lo[i] : state[i] == 1 ? hi[i] : 0; end
        if !isempty(free)
            x[free] .= pinv(A[:, free]; rtol=1e-12) * (b-A*x)
        end
        all(lo .- 1e-10 .<= x .<= hi .+ 1e-10) || continue
        best = min(best, sum(abs2, A*x-b))
    end
    return best
end

@testset "explicit diagnostic settings / Julia 1.13" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    @test Threads.nthreads() == 1
    @test VP.read_options(BASE_OPTIONS).outer_maxeval == 80
    @test_throws ErrorException VP.read_options(Dict())
    @test_throws ErrorException VP.read_options(merge(BASE_OPTIONS, Dict("hidden_knob" => 2)))
    for k in keys(BASE_OPTIONS), bad in (0, -1, NaN, Inf, true)
        @test_throws ErrorException VP.read_options(merge(BASE_OPTIONS, Dict(k => bad)))
    end
    @test_throws ErrorException options(outer_maxeval=2.5)
    @test_throws ErrorException options(linear_svd_rtol=1.0)
end

@testset "box LS active bounds, KKT signs and residuals" begin
    A = Matrix{Float64}(I, 4, 4)
    b = [-2.0, 0.4, 3.0, 5.0]
    lo, hi = [0.2, 0.2, 0.2, 1.0], [1.0, 1.0, 1.0, 1.0]
    r = VP.box_lsq(A, b, lo, hi; x0=zeros(4), options=OPTS)
    @test r.converged && r.feasible
    @test r.x ≈ [0.2, 0.4, 1.0, 1.0]
    @test r.active_lower == [1, 4]
    @test r.active_upper == [3, 4]
    @test r.gradient[1] >= 0 && r.gradient[3] <= 0
    @test abs(r.gradient[2]) <= r.kkt_tolerance
    @test r.rss ≈ sum(abs2, A*r.x-b)
    @test r.residual_norm^2 ≈ r.rss
    @test r.kkt_violation <= r.kkt_tolerance
    # Positive finite lower bound matters: nonnegative NNLS gives a different answer.
    nonnegative = VP.box_lsq(A, b, zeros(4), fill(10.0, 4); x0=zeros(4), options=OPTS)
    @test nonnegative.x[1] == 0.0 != r.x[1]
    # Wrong-sign active gradients must release a bound.
    released = VP.box_lsq(A, fill(0.6, 4), zeros(4), ones(4); x0=ones(4), options=OPTS)
    @test released.converged && released.x ≈ fill(0.6, 4)
    limited = VP.box_lsq(A, fill(0.6, 4), zeros(4), ones(4); x0=ones(4), options=options(linear_maxiter=1))
    @test !limited.converged && limited.status == "iteration_limit"
    @test limited.feasible && isfinite(limited.rss)
end

@testset "box LS agrees with all-face oracle, including degeneracy" begin
    rng = MersenneTwister(287)
    for trial in 1:48
        n = 1+mod(trial, 4)
        A = randn(rng, 7, n)
        trial % 3 == 0 && n > 1 && (A[:, end] .= A[:, 1])
        trial % 4 == 0 && (A[:, 1] .= 0)
        b = randn(rng, 7)
        lo = randn(rng, n); hi = lo .+ rand(rng, n) .+ 0.1
        trial % 5 == 0 && (hi[1] = lo[1])
        r = VP.box_lsq(A, b, lo, hi; x0=randn(rng, n), options=OPTS)
        @test r.feasible
        @test r.converged
        @test r.kkt_violation <= r.kkt_tolerance
        @test isapprox(r.rss, box_oracle(A, b, lo, hi); atol=1e-9, rtol=1e-9)
        @test isapprox(r.rss, sum(abs2, A*r.x-b); atol=1e-13)
    end
    for A in (zeros(6, 3), ones(6, 3), [1.0 1.0+1e-14; 2.0 2.0; 3.0 3.0])
        b = collect(1.0:size(A, 1))
        r = VP.box_lsq(A, b, fill(-1.0, size(A, 2)), ones(size(A, 2)); x0=zeros(size(A, 2)), options=OPTS)
        @test r.converged
        @test isapprox(r.rss, box_oracle(A, b, fill(-1.0, size(A, 2)), ones(size(A, 2))); atol=1e-12, rtol=1e-10)
    end
    # SVD truncation does not conceal failed stationarity.
    A = [1.0 0.0; 0.0 1e-4]
    r = VP.box_lsq(A, [0.0, 1.0], [-1.0,-1.0], [1.0,1.0]; x0=zeros(2), options=options(linear_svd_rtol=0.1))
    @test !r.converged && r.status == "stationarity_failure"
    @test_throws ArgumentError VP.box_lsq(zeros(0, 2), Float64[], zeros(2), ones(2); x0=zeros(2), options=OPTS)
    @test_throws ArgumentError VP.box_lsq([NaN;;], [1.0], [0.0], [1.0]; x0=[0.0], options=OPTS)
    @test_throws ArgumentError VP.box_lsq([1.0;;], [1.0], [2.0], [1.0]; x0=[0.0], options=OPTS)
    @test_throws DimensionMismatch VP.box_lsq(ones(2, 2), [1.0], zeros(2), ones(2); x0=zeros(2), options=OPTS)
end

@testset "exact native forward/decode mapping and full parameter count" begin
    rng = MersenneTwister(913)
    axis = (origin=(0.4,-0.3), axis=[0.8,0.6], perp=[-0.6,0.8], tmin=-1.7, tmax=1.8)
    x, y = randn(rng, 71), randn(rng, 71)
    for n in (1, 2, 4), circular in (false,true), tilted in (false,true), spacing in ("free", "uniform", "alternating"), shared in (0,1,2)
        cfg = cfg_fixture(; circular, tilted, spacing, shared)
        lo, hi = VP.native_raw_bounds(n, cfg)
        @test length(lo) == length(hi) == G._chain_nparams(n, cfg)
        @test length(VP.linear_layout(n, cfg).geometry) + n + (tilted ? 3 : 1) == length(lo)
        for p in (lo, hi, lo .+ rand(rng, length(lo)).*(hi-lo))
            kwargs = (amp_min=0.3, amp_range=0.7)
            A = VP.design_matrix(x, y, p, n, axis, cfg; kwargs...)
            coeff = VP.linear_coefficients(p, n, axis, cfg; kwargs...)
            check = VP.mapping_check(A, coeff, x, y, p, n, axis, cfg; kwargs..., options=OPTS)
            @test check.passed
            @test A*coeff ≈ G._chain_model_values(x, y, p, n, axis, cfg; kwargs...)
            recovered = VP.encode_linear(p, coeff, n, axis, cfg; kwargs...)
            @test all(lo .<= recovered .<= hi)
            @test recovered[VP.linear_layout(n,cfg).geometry] == p[VP.linear_layout(n,cfg).geometry]
            @test maximum(abs.(recovered-p)) < 1e-12
        end
        low, high = VP._linear_bounds(zeros(length(lo)), n, axis, cfg; amp_min=0.3, amp_range=0.7)
        @test low[end] ≈ 0.3+0.7*G._rsigmoid(-5.0)
        @test high[end] ≈ 0.3+0.7*G._rsigmoid(5.0)
        @test (low[1],high[1]) == (-5.0,5.0)
        tilted && (@test low[2:3] == [-1.0,-1.0] && high[2:3] == [1.0,1.0])
        @test VP.full_gcv(2.3, 71, n, cfg) ≈ 2.3*71/(71-G._chain_nparams(n,cfg))^2
        @test isinf(VP.full_gcv(2.3, 2, n, cfg))
    end
    cfg = cfg_fixture()
    n=2; p=zeros(G._chain_nparams(n,cfg))
    # The native EPS fallback is exponential, not the finite-fraction branch.
    kw=(amp_min=0.0,amp_range=G.EPS)
    c=VP.linear_coefficients(p,n,axis,cfg;kw...)
    @test c[end] == 1.0
    @test VP.encode_linear(p,c,n,axis,cfg;kw...) ≈ p
    cfg.peak_profile=:split
    @test_throws ErrorException VP.design_matrix(x,y,p,n,axis,cfg;kw...)
end

@testset "raw box matches compiled native fitter, not just copied constants" begin
    # Zero iterations plus skip_global leaves exactly the native clamped warm
    # start. This probes the private bound construction without a real fit or
    # source parsing. It will fail if core changes any private optimizer bound.
    x=collect(range(-2.0,2.0;length=31)); y=zeros(31); z=zeros(31)
    axis=(origin=(0.0,0.0),axis=[1.0,0.0],perp=[0.0,1.0],tmin=-2.0,tmax=2.0)
    cases=((1,false,false,"free",0),(3,true,true,"free",0),
           (4,false,false,"uniform",1),(4,false,true,"alternating",2))
    for (n,circular,tilted,spacing,shared) in cases
        cfg=cfg_fixture(;circular,tilted,spacing,shared)
        cfg.max_iter=0
        lo,hi=VP.native_raw_bounds(n,cfg)
        for side in (-100.0,100.0)
            native=G._fit_chain_n(x,[0.0],ones(1,31),x,y,z,1.0,n,axis,cfg;
                starts=1,warm_start=fill(side,length(lo)))
            @test native.success
            @test native.params == (side<0 ? lo .+ 1e-9 : hi .- 1e-9)
        end
    end
end

function synthetic_data(n, cfg)
    xs = collect(range(-1.5,1.5;length=31)); ys=collect(range(-0.4,0.4;length=9))
    x = [a for b in ys for a in xs]; y = [b for b in ys for a in xs]
    axis = (origin=(0.0,0.0), axis=[1.0,0.0], perp=[0.0,1.0], tmin=-1.5,tmax=1.5)
    p = zeros(G._chain_nparams(n,cfg)); p[1]=0.05
    cfg.chain_tilted_baseline && (p[2:3] .= [0.02,-0.03])
    p[VP.linear_layout(n,cfg).amplitudes] .= range(-1.0,1.0;length=n)
    z = G._chain_model_values(x,y,p,n,axis,cfg;amp_min=0.3,amp_range=0.7)
    z .+= [0.002*sin(0.71*i) for i in eachindex(z)]
    zimg=permutedims(reshape(z,length(xs),length(ys)))
    return (xs=xs,ys=ys,zimg=zimg,x=x,y=y,z=z,zfull=z,noise=0.02,axisctx=axis), p
end

@testset "fixed-N profile and native circular→elliptical refinement" begin
    circ=cfg_fixture(circular=true)
    data, truth=synthetic_data(3,circ)
    native = G._fit_chain_n(data.xs,data.ys,data.zimg,data.x,data.y,data.z,data.noise,3,data.axisctx,circ; starts=1)
    @test native.success
    @test length(native.params) == G._chain_nparams(3,circ)
    VP.finalize_native!(native,data,circ)
    @test native.gcv == VP.full_gcv(native.rss,length(data.z),3,circ)
    kw=(amp_min=native.amp_min,amp_range=native.amp_range,options=OPTS)
    fixed=VP.fixed_geometry_profile(native.params,3,data.x,data.y,data.z,data.axisctx,circ;kw...)
    @test fixed.success && fixed.linear.converged
    @test fixed.rss <= native.rss + 1e-11
    @test fixed.params[VP.linear_layout(3,circ).geometry] == native.params[VP.linear_layout(3,circ).geometry]
    outer=VP.refine_profile(native.params,3,data.x,data.y,data.z,data.axisctx,circ;kw...,initial=fixed)
    @test outer.best.success && outer.best.rss <= fixed.rss
    @test 0 < outer.evaluations <= OPTS.outer_maxeval
    @test outer.failed_evaluations == 0
    @test !isempty(outer.status)
    outer.status in ("MAXEVAL_REACHED","MAXTIME_REACHED") && (@test !outer.converged)
    ell=cfg_fixture(circular=false)
    warm=VP.circular_to_elliptical(native.params,3,ell)
    @test length(warm) == G._chain_nparams(3,ell)
    # Same raw-width duplication as native extraction, actual predictions equal
    # because this fixture has identical parallel/perpendicular width bounds.
    @test G._chain_model_values(data.x,data.y,warm,3,data.axisctx,ell;amp_min=native.amp_min,amp_range=native.amp_range) ≈
          G._chain_model_values(data.x,data.y,native.params,3,data.axisctx,circ;amp_min=native.amp_min,amp_range=native.amp_range)
    @test_throws ErrorException VP.circular_to_elliptical(native.params,3,circ)
    @test_throws ErrorException VP.fixed_geometry_profile(fill(100.0,length(native.params)),3,data.x,data.y,data.z,data.axisctx,circ;kw...)
    impossible=100
    @test_throws ErrorException VP.fixed_geometry_profile(zeros(G._chain_nparams(impossible,circ)),impossible,data.x,data.y,data.z,data.axisctx,circ;kw...)
    # Serialization retains unsuccessful N rows and native convergence unknown.
    ios=(fits=IOBuffer(),params=IOBuffer(),lobes=IOBuffer())
    prep=(file="synthetic.sxm", options=options(outer_maxeval=10,outer_maxtime_s=1.0))
    r, ok = CLI.process_geometry!(ios,prep,data,3,"ell",ell;warm_start=warm)
    @test r !== nothing && r.success
    failed, failedok = CLI.process_geometry!(ios,prep,data,impossible,"circ",circ)
    @test failed !== nothing && !failed.success && !failedok
    lines=split(chomp(String(take!(ios.fits))),'\n')
    @test length(lines) == 6
    rows=[Dict(zip(CLI.FIT_COLUMNS,split(line,'\t';keepempty=true))) for line in lines]
    @test all(row -> length(row)==length(CLI.FIT_COLUMNS), rows)
    @test rows[1]["optimizer_converged"] == "unknown"
    @test rows[1]["optimizer_status"] == "unknown_native_not_exposed"
    @test all(row -> row["status"]=="failed",rows[4:6])
    @test all(row -> row["N"]=="100",rows[4:6])
    @test parse(Float64,rows[2]["rss"]) <= parse(Float64,rows[1]["rss"])+1e-11
    @test all(row -> row["p_full"]==string(G._chain_nparams(3,ell)),rows[1:3])
    @test !isempty(String(take!(ios.params)))
    @test length(split(chomp(String(take!(ios.lobes))),'\n')) == 9
end

@testset "metadata-only CLI, adaptive support, collisions, failure rows" begin
    mktempdir() do dir
        raw=joinpath(dir,"one.sxm"); write(raw,"not an SXM file: dry run must not parse this")
        summary=joinpath(dir,"summary.tsv")
        write(summary,"filepath\tN_selected\trefined_policy\none.sxm\t9\tadaptive_support_rescue_robust_guard\n")
        diag=joinpath(dir,"diagnostic.toml")
        open(io -> TOML.print(io,Dict("counting_variable_projection"=>BASE_OPTIONS)),diag,"w")
        cfg=joinpath(@__DIR__,"..","config","chitosan_10_20mer_adaptive_support_rescue.toml")
        out=joinpath(dir,"new")
        args=["--file",raw,"--config",cfg,"--diagnostic-config",diag,"--selected-summary",summary,"--candidate-ns","3,9,100","--outdir",out]
        before=read(raw)
        p=CLI.prepare(vcat(args,["--dry-run"]))
        @test p.context == (n=9,use_rescue=true)
        @test p.ns == [3,9,100] # neither fixed-six nor auto-range from selected N
        @test p.circ.support_noise_k == p.ell.support_noise_k == p.config["model"]["adaptive_rescue_support_noise_k"]
        @test p.circ.support_padding_nm == p.ell.support_padding_nm == p.config["model"]["adaptive_rescue_support_padding_nm"]
        @test p.ell.skip_global && p.ell.max_iter == OPTS.native_elliptical_maxiter
        @test CLI.main(vcat(args,["--dry-run"])) == 0
        @test !ispath(out) && read(raw)==before
        @test_throws ErrorException CLI.prepare(vcat(args,["--expected-n","6"]))
        @test_throws ErrorException CLI.prepare(vcat(args,["--file",raw]))
        bad=copy(args); bad[findfirst(==("3,9,100"),bad)]="3,3"
        @test_throws ErrorException CLI.prepare(bad)
        bad[findfirst(==("3,3"),bad)]="0"
        @test_throws ErrorException CLI.prepare(bad)
        mkdir(out); @test_throws ErrorException CLI.prepare(args); rm(out)
        write(out,"collision"); @test_throws ErrorException CLI.prepare(args); rm(out)
        symlink(joinpath(dir,"absent"),out); @test_throws ErrorException CLI.prepare(args); rm(out)
        write(summary,"filepath\tN_selected\none.sxm\t9\n")
        @test_throws ErrorException CLI.prepare(args)
        write(summary,"filepath\tN_selected\trefined_policy\none.sxm\t9\tadaptive_support_rescue_keep\n")
        p=CLI.prepare(args)
        @test !p.context.use_rescue
        @test p.circ.support_padding_nm == p.config["model"]["support_padding_nm"]
        # Intentional malformed synthetic raw fixture: no fit happens, every
        # requested candidate/method survives in the output as a failed row.
        @test CLI.run(p) == 1
        fitlines=readlines(joinpath(out,"fits.tsv"))
        @test length(fitlines)==1+length(p.ns)*6
        @test all(line -> occursin("data_preparation_failed",line),fitlines[2:end])
        metadata=TOML.parsefile(joinpath(out,"metadata.toml"))
        @test metadata["status"]=="failed"
        @test metadata["effective_preprocessing"]["initial_sigma_nm"]=="nothing"
        @test metadata["input_sha256"]["file"]==CLI._sha(raw)
        @test metadata["candidate_ns"]==[3,9,100]
        @test_throws ErrorException CLI.prepare(args)
        @test read(raw)==before
    end
end
