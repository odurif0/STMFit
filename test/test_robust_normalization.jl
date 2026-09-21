using Test, Statistics, Random, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
include(joinpath(@__DIR__, "lib", "reconstructed_representation_diagnostics.jl"))
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT, "config", "unit_assignment_patch_support.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_robust_normalization.toml")

function records_for(values; file="synthetic.sxm")
    [GMM.LobeRecord(file, i, 1.0, Dict("f" => Float64(x))) for (i,x) in enumerate(values)]
end
standardize(records; mode="median_iqr", fallback=1.0, kwargs...) =
    GMM._standardized_matrix(records, ["f"]; normalization=mode, scale_fallback=fallback, kwargs...)

@testset "One explicit normalization change, no selection or class prior change" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    a,b = load_config(CONTROL),load_config(CANDIDATE)
    @test load_gmm_normalization(a) == (mode="mean_sample_std", scale_fallback=1.0)
    @test load_gmm_normalization(b) == (mode="median_iqr", scale_fallback=1.0)
    @test b["selection"]["assignment_training_support"] == "all_admissible"
    @test_throws ErrorException ReconstructedRepresentationDiagnostics.load_inputs(Dict(),CANDIDATE)
    b["preprocessing"]["gmm_feature_normalization"] = a["preprocessing"]["gmm_feature_normalization"]
    b["model"]["name"] = a["model"]["name"]
    @test a == b
    mktempdir() do dir
        for field in ("gmm_feature_normalization", "gmm_scale_fallback"),
            value in (nothing, "guess", true, -1.0, 0.0, NaN, Inf)
            cfg = deepcopy(a)
            value === nothing ? delete!(cfg["preprocessing"],field) : (cfg["preprocessing"][field]=value)
            @test_throws ArgumentError load_gmm_normalization(cfg)
            path=joinpath(dir,"bad.toml")
            open(io -> TOML.print(io,cfg),path,"w")
            @test_throws ArgumentError load_config(path)
        end
    end
end

@testset "Exact median and Type-7 IQR, without clipping or a fitted floor" begin
    records = records_for([1,2,3,4,5,1e6])
    z,valid = standardize(records)
    @test all(valid) && size(z)==(6,1)
    @test z[:,1] == ([1,2,3,4,5,1e6] .- 3.5) ./ 2.5
    @test z[end,1] > 1e5 # Robust moments do not bound extreme observations.
    altered=deepcopy(records); altered[end].features["f"]=1e12
    @test first(standardize(altered))[1:5,:] == z[1:5,:]
    for values in ([2.0], fill(3.0,4), [1.,1.,1.,1.,9.], [1.,2.], [1.,2.,3.])
        x=Float64.(values)
        q=quantile(x,[0.25,0.75];alpha=1,beta=1)
        scale=q[2]-q[1]; scale=scale>0 ? scale : 2.0
        actual,valid=standardize(records_for(x);fallback=2.0)
        @test all(valid) && actual[:,1] == (x .- median(x))./scale
    end
    @test first(standardize(records_for([1.,1.,1.,1.,9.]);fallback=2.0))[:,1] == [0,0,0,0,4]
    tiny=[1.,2.,3.,4.].*1e-30
    @test first(standardize(records_for(tiny)))[:,1] ≈ [-1.,-1/3,1/3,1.]
    @test_throws ArgumentError standardize(records;mode="unknown")
    for fallback in (0.,-1.,NaN,Inf,true)
        @test_throws ArgumentError standardize(records;fallback)
    end
end

@testset "Legacy arithmetic, per-feature finite support, masks and scan isolation" begin
    records=vcat(records_for([1,2,3,4,5,1e6]),records_for([2.,NaN,Inf,-Inf];file="second.sxm"))
    for mode in ("mean_sample_std","median_iqr")
        z,valid=standardize(records;mode)
        @test valid == [trues(7);falses(3)]
        @test z[7,1]==0 && all(isnan,z[8:10,1])
        @test [r.lobe for r in records] == [1:6;1:4]
        changed=deepcopy(records); changed[7].features["f"]=9e8
        @test isequal(first(standardize(changed;mode))[1:6,:],z[1:6,:])
        renamed=deepcopy(records)
        for r in renamed; r.file="different_"*r.file; r.lobe+=100; end
        @test isequal(first(standardize(renamed;mode)),z)
        order=[10,3,1,7,4,6,5,2,9,8]
        @test isequal(first(standardize(records[order];mode)),z[order,:])
        mask=trues(10); mask[6]=false
        trained,_=standardize(records;mode,training_mask=mask)
        altered=deepcopy(records); altered[6].features["f"]=-1e12
        @test isequal(first(standardize(altered;mode,training_mask=mask))[mask,:],trained[mask,:])
        _,empty=standardize(records;mode,training_mask=falses(10))
        @test !any(empty)
        @test_throws ArgumentError standardize(records;mode,training_mask=trues(9))
    end
    z,_=standardize(records;mode="mean_sample_std")
    x=[1.,2.,3.,4.,5.,1e6]
    @test z[1:6,1] == (x .- mean(x))./std(x) # Exact legacy operations.
    @test first(standardize(records;training_mask=[trues(5);falses(5)]))[1:6,1] ==
        (x .- 3.0)./2.0
    affine=records_for(7 .+ 3 .* x)
    @test first(standardize(affine)) ≈ first(standardize(records_for(x)))
    for r in records; r.features["g"]=r.lobe==2 ? NaN : 2r.lobe; end
    z,valid=GMM._standardized_matrix(records,["f","g"];normalization="median_iqr",scale_fallback=1.0,interactions=true)
    @test size(z)==(10,3)
    @test isequal(z[:,3],z[:,1].*z[:,2])
    @test !valid[2] && isfinite(z[2,1]) # Missing other features do not discard observed f.
    @test z[1:6,1] == first(standardize(records))[1:6,1]
end

@testset "GMM API and CLI agree; free weights, deterministic output and unavailable rows" begin
    rng=MersenneTwister(719)
    records=GMM.LobeRecord[]
    for file in 1:8, lobe in 1:10
        high=lobe>7
        push!(records,GMM.LobeRecord("synthetic_$file.sxm",lobe,high ? 4.0 : 1.0,
            Dict("f1"=>(high ? 4. : -2.)+0.3randn(rng)+file,
                 "f2"=>(high ? 3. : -1.)+0.4randn(rng)-file)))
    end
    push!(records,GMM.LobeRecord("unavailable.sxm",1,1.0,Dict("f1"=>NaN,"f2"=>NaN)))
    mktempdir() do dir
        features=joinpath(dir,"features.tsv")
        rows=[Dict("file"=>r.file,"lobe"=>string(r.lobe),"amplitude"=>string(r.amplitude),
            "f1"=>string(r.features["f1"]),"f2"=>string(r.features["f2"])) for r in records]
        header=["file","lobe","amplitude","f1","f2"]
        write_table(features,header,rows)
        for (name,cfg) in (("control",CONTROL),("robust",CANDIDATE))
            args=["--features",features,"--config",cfg,"--view","test=f1,f2","--seeds","2","--selftrain","2","--interactions"]
            opt=GMM._parse_cli(args)
            @test opt.normalization == (name=="robust" ? "median_iqr" : "mean_sample_std")
            @test opt.scale_fallback == 1.0
            diag=[]
            p=GMM._view_probability(records,["f1","f2"],opt;diagnostics=diag)
            @test length(p)==81 && all(isfinite,p[1:80]) && isnan(p[81])
            @test isequal(p,GMM._view_probability(records,["f1","f2"],opt))
            @test all(d->d.indices==collect(1:80),diag)
            @test all(d->length(d.clusters)==2 && d.clusters[1].weight != d.clusters[2].weight,diag)
            @test all(d->sum(c.weight for c in d.clusters) ≈ 1,diag)
            expected,out=joinpath.(dir,(name*"_expected.tsv",name*"_cli.tsv"))
            GMM._write_predictions(expected,records,p,Int.(isfinite.(p)))
            @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args --out $out`)
            @test read(expected)==read(out)
            @test_throws ErrorException GMM._parse_cli(vcat(args,["--expected-N","6"]))
            @test_throws ErrorException GMM._parse_cli(vcat(args,["--truth","forbidden.tsv"]))
        end
        forbidden=joinpath(dir,"forbidden.tsv")
        write_table(forbidden,vcat(header,["expected_N"]),rows)
        @test_throws ErrorException lobe_table(forbidden)
    end
end
