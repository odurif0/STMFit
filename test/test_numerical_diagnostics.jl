using Test, Random, Statistics
include(joinpath(@__DIR__, "diagnose_numerical_assignment.jl"))
const D = NumericalAssignmentDiagnostics
const E = D.EmpiricalFisherNative
const R = D.ReconstructedUnitAssignment
const G = D.GMM

@testset "Fixed numerical diagnostics on synthetic scans, without truth" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    @test D.main(["--help"]) == 0
    @test_throws ErrorException D.main(["--truth","x","--centered","x","--shrunk","x","--outdir","x"])
    mktempdir() do dir
        rng = MersenneTwister(210921)
        paths = joinpath.(dir,("control","centered","shrunk"))
        foreach(mkpath,paths)
        configs = D.config.(("affine_residual","centered_fisher","shrunk_gmm"))
        options = E.load_fisher_config(configs[1]); grid = E.fisher_grid(options)
        n = 80
        pixels = 0.05randn(rng,n,grid.side^2)
        shape = [exp(-20(u*u+t*t)) for u in grid.coords for t in grid.coords]
        for i in 1:n
            pixels[i,:] .+= (i<=53 ? -2.0 : 2.0) .* shape .+ (isodd(i) ? 0.9 : 0.4)
        end
        pixels[end,1] = NaN
        keys = [("synthetic_$(lpad(cld(i,8),2,'0')).sxm",mod1(i,8)) for i in 1:n]
        patchrows = [Dict{String,Any}("file"=>k[1],"lobe"=>k[2]) for k in keys]
        for i in 1:n, j in 1:grid.side^2
            patchrows[i]["res_p"*lpad(j,3,'0')] = pixels[i,j]
        end
        header = vcat(["file","lobe"],["res_p"*lpad(j,3,'0') for j in 1:grid.side^2])
        patchpath = joinpath(paths[1],"patches_fwd17.tsv")
        R.write_table(patchpath,header,patchrows)
        patches = E.load_patches(patchpath,"res",options)
        scores = [E.cv_scores(patches,E.load_fisher_config(c)) for c in configs[1:2]]
        for m in 1:2; E.write_scores(joinpath(paths[m],"fisher_cv.tsv"),scores[m]); end
        featurevalues = randn(rng,n,7)
        featurevalues[:,1] .= pixels[:,grid.center_index]
        for m in 1:2
            saved = last(R.lobe_table(joinpath(paths[m],"fisher_cv.tsv")))
            rows = Dict{String,Any}[]
            for (i,k) in enumerate(keys)
                row = Dict{String,Any}("file"=>k[1],"lobe"=>k[2],"amplitude"=>0.5+0.1featurevalues[i,1])
                for j in 1:7; row[D.FEATURES[j]]=featurevalues[i,j]; end
                row["emp_fisher"]=saved[k]["score"]
                push!(rows,row)
            end
            R.write_table(joinpath(paths[m],"features_predictor.tsv"),vcat(["file","lobe","amplitude"],D.FEATURES),rows)
        end
        cp(joinpath(paths[1],"features_predictor.tsv"),joinpath(paths[3],"features_predictor.tsv"))
        records = G._load_records(joinpath(paths[1],"features_predictor.tsv"))
        groups = D.scan_partitions(records)
        @test sort(vcat(groups...)) == collect(1:n)
        @test isempty(intersect(Set(r.file for r in records[groups[1]]),Set(r.file for r in records[groups[2]])))
        @test length.(groups) == [40,40]
        @test_throws ErrorException D.scan_partitions(records[1:8])
        for m in (1,3)
            opt = G._parse_cli(["--config",configs[m],"--features",joinpath(paths[m],"features_predictor.tsv"),
                "--view","v_cc="*join(D.FEATURES,','),"--seeds","10","--selftrain","2","--interactions"])
            probs = G._view_probability(records,D.FEATURES,opt)
            G._write_predictions(joinpath(paths[m],"pred_gmm.tsv"),records,probs,Int.(isfinite.(probs)))
        end
        out = joinpath(dir,"diagnostics")
        @test D.main(["--control",paths[1],"--centered",paths[2],"--shrunk",paths[3],"--outdir",out]) == 0
        _,covs = R.read_table(joinpath(out,"covariances.tsv"))
        @test length(covs) == 120
        @test Set(r["partition"] for r in covs) == Set(["full","scan_half_1","scan_half_2"])
        @test all(r -> 0 <= parse(Float64,r["shrinkage"]) <= 1,covs)
        @test all(r -> parse(Float64,r["eigen_min"]) > 0,covs)
        _,z = R.read_table(joinpath(out,"standardized_features.tsv"))
        @test length(z) == n
        _,s = R.read_table(joinpath(out,"scan_withdrawal.tsv"))
        @test length(s) == 2n
        @test Set((r["file"],parse(Int,r["lobe"])) for r in s) == Set(keys)
        for mode in ("ridge","ledoit_wolf")
            @test count(r -> r["mode"]==mode,s) == n
        end
        @test_throws ErrorException D.run_diagnostics(paths...,out)
        tablepath=joinpath(paths[1],"features_predictor.tsv")
        lines=readlines(tablepath)
        lines[1] *= "\texpected_N"
        lines[2:end] .*= "\t6"
        write(tablepath,join(lines,'\n')*"\n")
        @test_throws ErrorException D.run_diagnostics(paths...,joinpath(dir,"forbidden"))
        @test !ispath(joinpath(dir,"forbidden"))
    end
end
