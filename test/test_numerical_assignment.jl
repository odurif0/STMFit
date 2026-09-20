using Test, LinearAlgebra, Random, Statistics, TOML, Printf
include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .EmpiricalFisherNative, .ReconstructedUnitAssignment
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
using .GMM.AssignmentCovariance
const ROOT = dirname(@__DIR__)
const BASE_PATH = joinpath(ROOT, "config", "unit_assignment_affine_residual.toml")
const CENTER_PATH = joinpath(ROOT, "config", "unit_assignment_centered_fisher.toml")
const SHRINK_PATH = joinpath(ROOT, "config", "unit_assignment_shrunk_gmm.toml")
const BASE = load_fisher_config(BASE_PATH)
const CENTER = load_fisher_config(CENTER_PATH)
const GRID = fisher_grid(BASE)

@testset "Explicit independent modes, no composition or vote changes" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    @test BASE.score_center == "legacy_centered_mean"
    @test CENTER.score_center == "training_mean"
    control = load_config(BASE_PATH)
    for (path, field, value) in ((CENTER_PATH, "fisher_score_center", "training_mean"),
                               (SHRINK_PATH, "gmm_final_covariance", "ledoit_wolf"))
        candidate = load_config(path)
        @test candidate["model"][field] == value
        candidate["model"][field] = control["model"][field]
        candidate["model"]["name"] = control["model"]["name"]
        @test candidate == control
    end
    for value in (nothing, "guess", true, 1.0)
        cfg = deepcopy(control)
        value === nothing ? delete!(cfg["model"], "fisher_score_center") : (cfg["model"]["fisher_score_center"] = value)
        @test_throws ArgumentError load_fisher_config(cfg)
    end
    for key in ("gmm_final_covariance", "gmm_covariance_ridge"), value in (nothing, "guess", true, -1.0, NaN, Inf)
        cfg = deepcopy(control)
        value === nothing ? delete!(cfg["model"], key) : (cfg["model"][key] = value)
        @test_throws ArgumentError load_covariance_config(cfg)
    end
    @test load_covariance_config(BASE_PATH) == (mode="ridge", ridge=1e-6)
    @test load_covariance_config(SHRINK_PATH) == (mode="ledoit_wolf", ridge=1e-6)
    mktempdir() do dir
        cfg = load_config(SHRINK_PATH); cfg["selection"]["gmm_selftrain"] = 0
        path = joinpath(dir, "invalid.toml")
        open(io -> TOML.print(io, cfg), path, "w")
        @test_throws ErrorException load_config(path)
    end
end

function fisher_fixture()
    rng = MersenneTwister(20260921)
    n = 64
    full = 0.06randn(rng, n, GRID.side^2)
    shape = [exp(-20(u*u+t*t)) + 0.15sin(9u)*cos(11t) for u in GRID.coords for t in GRID.coords]
    background = [0.8 + 0.3u + 0.2t*t for u in GRID.coords for t in GRID.coords]
    for i in 1:n
        full[i,:] .+= (i <= 44 ? -2.0 : 2.0) .* shape .+ (isodd(i) ? 1.5 : 0.5) .* background
    end
    patches = PatchTable([("synthetic.sxm", i) for i in 1:n], full[:,GRID.disk_indices],
        full[:,GRID.center_index], fill("",n), GRID)
    return patches, full
end

@testset "Training mean changes only the score origin, not learned chemistry" begin
    patches, full = fisher_fixture()
    original = copy(patches.X)
    legacy = fit_fisher(patches.X, patches.amplitudes, BASE)
    centered = fit_fisher(patches.X, patches.amplitudes, CENTER)
    @test centered.w_p == legacy.w_p
    @test centered.g0 == legacy.g0 && centered.g1 == legacy.g1
    @test centered.gmm.weights == legacy.gmm.weights
    @test centered.gmm.means == legacy.gmm.means
    @test centered.gmm.assignments == legacy.gmm.assignments
    @test centered.amplitude_means == legacy.amplitude_means
    @test centered.gmm.weights[1] != centered.gmm.weights[2]
    @test centered.mid == vec(mean(patches.X; dims=1))
    @test norm(legacy.mid) < 1e-12
    @test norm(centered.mid) > 0.1
    shift = dot(centered.mid - legacy.mid, legacy.w_p)
    for i in (1, 3, 21, 47, 64)
        x = patches.X[i,:]
        @test score(x, centered) == dot(x - vec(mean(patches.X; dims=1)), legacy.w_p)
        @test score(x, centered) ≈ score(x, legacy) - shift atol=1e-12
        @test maxmirror_score(x, centered, GRID) ≈ maxmirror_score(x, legacy, GRID) - shift atol=1e-12
    end
    @test patches.X == original
    baseline = cv_scores(patches, BASE)
    treatment = cv_scores(patches, CENTER)
    offsets = Float64[]
    for parity in (0,1)
        train = findall(i -> mod(i,2) == parity, 1:64)
        held = findall(i -> mod(i,2) != parity, 1:64)
        model = fit_fisher(patches.X[train,:], patches.amplitudes[train], BASE)
        offset = dot(vec(mean(patches.X[train,:]; dims=1)) - model.mid, model.w_p)
        push!(offsets, offset)
        @test [treatment[i].score for i in held] ≈ [baseline[i].score - offset for i in held] atol=1e-12
    end
    @test !isapprox(offsets[1], offsets[2]; atol=1e-8)
    x = [r.score for r in baseline]; y = [r.score for r in treatment]
    z(v) = (v .- mean(v)) ./ std(v)
    @test z(x .+ 7.0) ≈ z(x) atol=1e-12 # a uniform offset would cancel downstream
    @test norm(z(x)-z(y)) > 1e-3 # opposite-fold offsets can survive
    changed = copy(patches.X); changed[1,:] .+= 0.2
    other = cv_scores(PatchTable(patches.keys,changed,patches.amplitudes,patches.invalid_reasons,GRID), CENTER)
    @test [r.score for r in treatment[3:2:end]] == [r.score for r in other[3:2:end]]
    bad = copy(patches.X); bad[2,1] = NaN
    invalid = cv_scores(PatchTable(patches.keys,bad,patches.amplitudes,patches.invalid_reasons,GRID), CENTER)
    @test [(r.file,r.lobe) for r in invalid] == patches.keys
    @test invalid[2].invalid_reason == "nonfinite_patch_input" && isnan(invalid[2].score)
    renamed = PatchTable([("renamed.sxm",i) for i in 1:64], patches.X,patches.amplitudes,patches.invalid_reasons,GRID)
    @test [r.score for r in cv_scores(renamed,CENTER)] == y
    mktempdir() do dir
        path = joinpath(dir,"patches.tsv")
        open(path,"w") do io
            println(io,join(vcat(["file","lobe"],[@sprintf("res_p%03d",i) for i in 1:289]),'\t'))
            for i in 1:64
                println(io,join(vcat(["synthetic.sxm",string(i)],string.(full[i,:])),'\t'))
            end
        end
        expected, actual = joinpath.(dir,("expected.tsv","actual.tsv"))
        write_scores(expected,treatment)
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_empirical_fisher_native.jl")) --patches $path --prefix res --config $CENTER_PATH --out $actual`)
        @test read(expected) == read(actual)
    end
end

@testset "Ledoit-Wolf algebra, degenerate data and fixed numerical floor" begin
    rng = MersenneTwister(812)
    X = randn(rng,6,90); X[6,:] = X[2,:] .+ 1e-7randn(rng,90)
    C = X .- mean(X; dims=2)
    n,p = size(C,2),size(C,1)
    S = C*C'/n; target = (tr(S)/p)*Matrix{Float64}(I,p,p)
    beta = sum(sum(abs2, x*x' - S) for x in eachcol(C))/n^2
    delta = sum(abs2,S-target)
    expected_lambda = clamp(beta/delta,0,1)
    raw = final_covariance(C; mode="ridge",ridge=1e-6)
    shrunk = final_covariance(C; mode="ledoit_wolf",ridge=1e-6)
    @test raw.covariance == S+1e-6I
    @test raw.shrinkage == 0
    @test shrunk.shrinkage ≈ expected_lambda rtol=1e-13
    @test shrunk.covariance ≈ (1-expected_lambda)*S+expected_lambda*target+1e-6I rtol=1e-13
    @test tr(shrunk.covariance) ≈ tr(raw.covariance)
    @test cond(shrunk.covariance) < cond(raw.covariance)
    @test isposdef(Symmetric(shrunk.covariance))
    @test C == X .- mean(X;dims=2)
    reversed = final_covariance(C[:,end:-1:1];mode="ledoit_wolf",ridge=1e-6)
    @test reversed.shrinkage ≈ shrunk.shrinkage rtol=1e-13
    @test reversed.covariance ≈ shrunk.covariance rtol=1e-13
    rotated = Matrix(qr(randn(rng,p,p)).Q)
    rotation = final_covariance(rotated*C;mode="ledoit_wolf",ridge=1e-6)
    @test rotation.shrinkage ≈ shrunk.shrinkage rtol=1e-13
    @test rotation.covariance ≈ rotated*shrunk.covariance*rotated' rtol=1e-12
    one = final_covariance(zeros(4,1);mode="ledoit_wolf",ridge=1e-6)
    @test one.shrinkage == 0 && one.covariance == 1e-6Matrix{Float64}(I,4,4)
    spherical = hcat(Matrix{Float64}(I,3,3),-Matrix{Float64}(I,3,3))
    @test final_covariance(spherical;mode="ledoit_wolf",ridge=1e-6).shrinkage == 0
    for mode in ("ridge","ledoit_wolf")
        @test_throws ArgumentError final_covariance(zeros(0,2);mode,ridge=1e-6)
        @test_throws ArgumentError final_covariance(zeros(2,0);mode,ridge=1e-6)
        @test_throws ArgumentError final_covariance(fill(NaN,2,2);mode,ridge=1e-6)
        for ridge in (0.,-1.,NaN,Inf,true)
            @test_throws ArgumentError final_covariance(C;mode,ridge)
        end
    end
    @test_throws ArgumentError final_covariance(C;mode="guess",ridge=1e-6)
end

@testset "Shrink only the final hard covariance; retain means, weights and memberships" begin
    rng = MersenneTwister(971)
    data = hcat(randn(rng,4,60) .- 2,randn(rng,4,28) .+ 2)
    fit = GMM._gmm_fit(data,0;ridge=1e-6)
    da,db = [],[]
    a = GMM._mahalanobis_self_train(data,deepcopy(fit)...;iters=2,
        covariance_mode="ridge",ridge=1e-6,diagnostics=da)
    b = GMM._mahalanobis_self_train(data,deepcopy(fit)...;iters=2,
        covariance_mode="ledoit_wolf",ridge=1e-6,diagnostics=db)
    @test a[1] == b[1] && a[3] == b[3]
    @test a[3][1] != a[3][2] && sum(a[3]) ≈ 1
    @test length(da) == length(db) == 2
    for (ra,rb) in zip(da,db)
        @test ra.members == rb.members && ra.mean == rb.mean
        @test ra.sample == rb.sample && ra.weight == rb.weight
        @test 0 < rb.shrinkage <= 1
        @test cond(rb.covariance) <= cond(ra.covariance)
        C = data[:,ra.members] .- ra.mean
        @test ra.covariance == C*C'/length(ra.members)+1e-6I
        @test rb.covariance == final_covariance(C;mode="ledoit_wolf",ridge=1e-6).covariance
    end
    @test_throws ErrorException GMM._mahalanobis_self_train(data,deepcopy(fit)...;
        iters=0,covariance_mode="ledoit_wolf",ridge=1e-6)
end

@testset "GMM CLI selects one configured covariance mode and retains invalid keys" begin
    mktempdir() do dir
        rng = MersenneTwister(157)
        path = joinpath(dir,"features.tsv")
        names = ["f$i" for i in 1:8]
        open(path,"w") do io
            println(io,join(vcat(["file","lobe","amplitude"],names),'\t'))
            for i in 1:80
                amp = 0.05 + 0.1rand(rng)
                values = randn(rng,8); values[1] = amp
                i == 80 && (values[2] = NaN)
                println(io,join(vcat(["synthetic_$(cld(i,8)).sxm",string(mod1(i,8)),string(amp)],string.(values)),'\t'))
            end
        end
        records = GMM._load_records(path)
        for (mode,config) in (("ridge",BASE_PATH),("ledoit_wolf",SHRINK_PATH))
            out = joinpath(dir,mode*".tsv")
            args = ["--config",config,"--features",path,"--out",out,"--seeds","2",
                "--selftrain","2","--interactions","--view","synthetic="*join(names,',')]
            opt = GMM._parse_cli(args)
            @test opt.covariance_mode == mode && opt.covariance_ridge == 1e-6
            diag = []
            probabilities = GMM._view_probability(records,names,opt;diagnostics=diag)
            @test length(probabilities) == 80 && count(isfinite,probabilities) == 79
            @test length(diag) == 2
            @test all(d -> length(d.indices) == 79,diag)
            @test all(c -> size(c.covariance) == (36,36),Iterators.flatten(d.clusters for d in diag))
            @test isequal(probabilities,GMM._view_probability(records,names,opt))
            renamed = [GMM.LobeRecord("renamed_"*r.file,r.lobe,r.amplitude,copy(r.features)) for r in records]
            @test isequal(probabilities,GMM._view_probability(renamed,names,opt))
            expected = joinpath(dir,mode*"_expected.tsv")
            GMM._write_predictions(expected,records,probabilities,Int.(isfinite.(probabilities)))
            @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args`)
            @test read(out) == read(expected)
            _, rows = lobe_table(out)
            @test length(rows) == 80 && count(r -> r["predicted"] == "?",values(rows)) == 1
        end
        for flag in ("--truth","--expected-N","--reference","--manifest","--control-sequence")
            @test_throws ErrorException GMM._parse_cli([flag,"forbidden"])
        end
        @test_throws ErrorException GMM._parse_cli(["--config",SHRINK_PATH,"--features",path,"--selftrain","0"])
    end
end
