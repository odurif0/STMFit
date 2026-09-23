using Test, Random, Printf, TOML
module Patches
include(joinpath(@__DIR__,"test_patch_preprocessing.jl"))
end
module Driver
include(joinpath(@__DIR__,"run_tangent_mold_comparison.jl"))
end
include(joinpath(@__DIR__,"score_tangent_mold_templates.jl"))
const S=TangentMoldScorer
const ROOT=dirname(@__DIR__)
const SETTINGS=joinpath(ROOT,"config/tangent_mold_projection.toml")
const ORIENTATION=joinpath(ROOT,"config/tangent_mold_orientation.toml")
const CONFIG=joinpath(ROOT,"config/unit_assignment_patch_support.toml")
const COUNT=joinpath(ROOT,"config/chitosan.toml")

@testset "Synthetic scorer CLI, shards and complete comparison" begin
    mktempdir() do dir
        raw=joinpath(dir,"raw"); mkdir(raw)
        rng=MersenneTwister(923); geometry=Dict{String,String}[]; patches=Dict{String,String}[]
        for fi in 1:3
            file="synthetic_$fi.sxm"; Patches.synthetic_sxm(joinpath(raw,file))
            for lobe in 1:3
                push!(geometry,Dict("file"=>file,"N"=>"3","lobe"=>string(lobe),"source"=>"ell",
                    "x_nm"=>string(.5+.4lobe),"y_nm"=>"1.2","axis_x"=>"1.0","axis_y"=>"0.0",
                    "sigma_parallel_nm"=>"0.18","sigma_perp_nm"=>"0.14","skew_ratio"=>"1.0"))
                r=Dict("file"=>file,"lobe"=>string(lobe),"amplitude"=>string(.1lobe))
                for i in 1:289; r[@sprintf("res_p%03d",i)]=string(randn(rng)); end
                push!(patches,r)
            end
        end
        features=Driver.write_table(joinpath(dir,"features.tsv"),sort(collect(keys(first(geometry)))),geometry)
        patch=Driver.write_table(joinpath(dir,"patches.tsv"),vcat(["file","lobe","amplitude"],[@sprintf("res_p%03d",i) for i in 1:289]),patches)
        templates=Dict{String,String}[]
        for parity in 0:1, mirror in 0:1
            c=randn(rng,289)
            for typ in 0:1
                r=Dict("name"=>"synthetic","type"=>string(typ),"parity"=>string(parity),"mirror"=>string(mirror))
                for i in 1:289; r[@sprintf("p%03d",i)]=string((2typ-1)*c[i]); end
                push!(templates,r)
            end
        end
        templ=Driver.write_table(joinpath(dir,"templates.tsv"),vcat(["name","type","parity","mirror"],[@sprintf("p%03d",i) for i in 1:289]),templates)
        out=joinpath(dir,"scores.tsv")
        args=["--patches",patch,"--templates",templ,"--features",features,"--data-dir",raw,"--count-config",COUNT,
            "--config",CONFIG,"--settings",SETTINGS,"--prefix","res","--out",out]
        S.main(vcat(args,["--dry-run"]))
        @test !ispath(out) && !ispath(out*".basis.tsv")
        S.main(args)
        _,scores=Driver.lobe_table(out); _,audit=Driver.read_table(out*".audit.tsv")
        @test length(scores)==9 && length(audit)==36
        @test all(r["reason"]=="ok" && parse(Int,r["rank"])==8 for r in audit)
        @test all(0<parse(Float64,r["contrast_fraction"])<=1 for r in audit)
        @test all(isfinite(parse(Float64,r["cost_margin"])) for r in values(scores))
        @test_throws ErrorException S.main(args)
        shard=joinpath(dir,"shard.tsv")
        S.main(vcat(replace.(args,out=>shard),["--chunk","2/2"]))
        _,shards=Driver.lobe_table(shard)
        @test length(shards)==3 && all(k[1]=="synthetic_2.sxm" for k in keys(shards))
        @test all(shards[k]==scores[k] for k in keys(shards))
        oriented=joinpath(dir,"oriented.tsv")
        oa=replace.(replace.(args,out=>oriented),SETTINGS=>ORIENTATION)
        S.main(oa)
        _,oscores=Driver.lobe_table(oriented); _,oaudit=Driver.read_table(oriented*".audit.tsv")
        @test Set(keys(oscores))==Set(keys(scores))
        @test all(r["reason"]=="ok" && parse(Int,r["rank"])==9 for r in oaudit)
        bh,br=Driver.read_table(out*".basis.tsv"); oh,obr=Driver.read_table(oriented*".basis.tsv")
        @test oh==vcat(bh,["b9"])
        @test all(r["basis_count"]=="9" for r in obr)
        @test all(a[k]==b[k] for (a,b) in zip(br,obr) for k in setdiff(bh,["basis_count"]))
        o_shard=joinpath(dir,"oriented_shard.tsv")
        S.main(vcat(replace.(oa,oriented=>o_shard),["--chunk","2/2"]))
        _,oss=Driver.lobe_table(o_shard)
        @test length(oss)==3 && all(oss[k]==oscores[k] for k in keys(oss))
        bad=joinpath(dir,"bad.tsv")
        @test_throws ErrorException S.main(vcat(replace.(args,out=>bad),["--chunk","0/2"]))
        @test !ispath(bad)
        # A numerical failure must retain diagnostics but never emit partial scores.
        broken=deepcopy(geometry); broken[1]["axis_x"]="9.0"
        badfeatures=Driver.write_table(joinpath(dir,"bad_features.tsv"),sort(collect(keys(first(broken)))),broken)
        @test_throws ErrorException S.main(replace.(replace.(args,out=>bad),features=>badfeatures))
        @test !ispath(bad) && !ispath(bad*".basis.tsv")
        comparison=joinpath(dir,"comparison")
        ds=["--data-dir",raw,"--features",features,"--split-features",features,"--templates",templ,
            "--count-config",COUNT,"--config",CONFIG,"--settings",SETTINGS,"--outdir",comparison]
        Driver.tangent_comparison(vcat(ds,["--dry-run"]))
        @test !ispath(comparison)
        calls=[]
        runner=function(out,stage,script,as)
            push!(calls,(stage,script,as)); dest=as[findfirst(==("--outdir"),as)+1]; mkpath(dest)
            Driver.write_table(joinpath(dest,"predictions.tsv"),sort(collect(keys(first(geometry)))),geometry)
        end
        Driver.tangent_comparison(ds;runner)
        @test first.(calls)==["reference","tangent"]
        @test all(c[2]=="run_reconstructed_chitosan.jl" for c in calls)
        @test !("--mold-tangent-settings" in calls[1][3])
        @test calls[2][3][findfirst(==("--mold-tangent-settings"),calls[2][3])+1]==SETTINGS
        for arm in ("reference","tangent")
            @test read(joinpath(comparison,arm,"features.tsv"))==read(features)
            @test read(joinpath(comparison,arm,"features_split.tsv"))==read(features)
        end
        fail=joinpath(dir,"failed")
        empty!(calls)
        oc=joinpath(dir,"orientation_comparison")
        od=vcat(replace.(replace.(ds,comparison=>oc),SETTINGS=>ORIENTATION),["--reference-settings",SETTINGS])
        Driver.tangent_comparison(vcat(od,["--dry-run"]))
        @test !ispath(oc)
        Driver.tangent_comparison(od;runner)
        @test first.(calls)==["reference","tangent"]
        @test calls[1][3][findfirst(==("--mold-tangent-settings"),calls[1][3])+1]==SETTINGS
        @test calls[2][3][findfirst(==("--mold-tangent-settings"),calls[2][3])+1]==ORIENTATION
        _,hs=Driver.read_table(joinpath(oc,"input_hashes.tsv"))
        @test "--reference-settings" in [r["input"] for r in hs]
        @test_throws ErrorException Driver.tangent_comparison(vcat(od,["--reference-settings",SETTINGS]))
        @test_throws ErrorException Driver.tangent_comparison(vcat(od,["--truth","forbidden"]))
        @test_throws ErrorException Driver.tangent_comparison(replace.(ds,comparison=>fail);runner=(a...)->error("intentional"))
        @test isfile(joinpath(fail,"failures.tsv")) && !ispath(joinpath(fail,"tangent"))
        opts=Dict("--config"=>CONFIG,"--count-config"=>COUNT,"--data-dir"=>raw,"--features"=>features,
            "--split-features"=>features,"--mold-tangent-settings"=>SETTINGS,"--outdir"=>joinpath(dir,"no_output"),"--dry-run"=>"true")
        Driver.execute_pipeline(opts)
        for key in ("--features","--split-features")
            o=copy(opts); delete!(o,key); @test_throws ErrorException Driver.execute_pipeline(o)
        end
        for key in ("--acquisition-shifts","--patch-frames","--residual-features-fwd","--residual-features-bwd",
            "--patches-fwd","--patches-bwd","--descriptor-patches")
            @test_throws ErrorException Driver.execute_pipeline(merge(opts,Dict(key=>"invalid.tsv")))
        end
        @test !ispath(opts["--outdir"])
        touch(joinpath(raw,"extra.sxm"))
        @test_throws ErrorException Driver.tangent_comparison(vcat(replace.(ds,comparison=>joinpath(dir,"extra")),["--dry-run"]))
    end
    shell=joinpath(ROOT,"hpc/compare_tangent_molds.sbatch")
    @test success(`bash -n $shell`)
    @test occursin("#SBATCH --time=02:00:00",read(shell,String))
    @test occursin("SLURM_JOB_ID:?",read(shell,String))
    orientation_shell=joinpath(ROOT,"hpc/compare_tangent_orientation.sbatch")
    @test success(`bash -n $orientation_shell`)
    @test occursin("#SBATCH --time=02:00:00",read(orientation_shell,String))
    @test occursin("--reference-settings config/tangent_mold_projection.toml",read(orientation_shell,String))
end
