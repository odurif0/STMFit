#!/usr/bin/env julia
# Saved-state application and arithmetic checks only; NEVER fits a learner/image.
module FrozenLearningVerification
using Test, TOML, LinearAlgebra, Printf
include(joinpath(@__DIR__,"diagnose_frozen_learning.jl"))
const D=FrozenLearningDiagnostic
const RU=D.RU

function verify(input,output)
    data=D.inputs(input)
    @testset "Complete frozen-learning saved-output verification" begin
        @test TOML.parsefile(joinpath(output,"assignment.toml"))==data.cfg
        D.settings(joinpath(output,"settings.toml"))
        @test D.Pipeline.selected_counts(joinpath(output,"selected_summary.tsv"))==data.counts
        _,manifest=RU.read_table(joinpath(output,"inputs.tsv"))
        source=only(r["path"] for r in manifest if basename(r["path"])=="selected_summary.tsv") |> dirname
        @test length(manifest)==length(data.hashes)
        @test Dict(normpath(joinpath(input,relpath(r["path"],source)))=>r["sha256"] for r in manifest)==data.hashes
        _,groups=RU.read_table(joinpath(output,"groups.tsv"))
        @test length(groups)==length(data.counts)
        @test Dict(r["file"]=>r["group"] for r in groups)==data.groups
        _,arms=RU.read_table(joinpath(output,"arms.tsv"))
        @test arms==[Dict(string(k)=>v for (k,v) in pairs(a)) for a in D.ARMS]
        fishopts=D.EF.load_fisher_config(data.cfg)
        loaded=Dict(a=>D.load_bank(joinpath(output,"native_"*a,"models.toml")) for a in ("control","observed"))
        patches=Dict(a=>D.EF.load_patches(joinpath(input,a,"patches_fwd17.tsv"),"res",fishopts) for a in ("control","observed"))
        for a in ("control","observed")
            dir=joinpath(output,"native_"*a); bank=loaded[a].bank; fisher=loaded[a].fisher
            valid=isempty.(patches[a].invalid_reasons)
            for p in 0:1
                @test fisher[p].train==findall(i->valid[i] && mod(patches[a].keys[i][2],2)==p,eachindex(valid))
                @test fisher[p].held==findall(i->valid[i] && mod(patches[a].keys[i][2],2)!=p,eachindex(valid))
                @test all(isfinite,fisher[p].w_p) && all(isfinite,fisher[p].mid)
            end
            rec=D.head_records(joinpath(dir,"features_predictor.tsv"),joinpath(input,a,"patches_bwd17.tsv"))
            opt=D.options(joinpath(dir,"features_predictor.tsv"),joinpath(input,a,"patches_bwd17.tsv"),
                joinpath(dir,"unused.tsv"),data.cfgpath,data.cfg)
            for h in (bank.g,values(bank.k)...)
                states=[]
                if h.kind=="gmm"
                    z,valid=D.Gaussian._standardized_matrix(rec.g,h.features;normalization=opt.g.normalization,
                        scale_fallback=opt.g.scale_fallback,interactions=opt.g.interactions,normalization_state=states)
                else
                    z,valid=D.KMeans._standardized_matrix(rec.k,h.features;interactions=opt.k.interactions,normalization_state=states)
                end
                @test isequal(h.scales,states)
                @test [m.seed for m in h.models]==collect(data.cfg["selection"]["first_seed"] .+
                    (0:(h.kind=="gmm" ? data.cfg["selection"]["gmm_seeds"] : data.cfg["selection"]["kmeans_seeds"])-1))
                for m in h.models
                    @test m.high_cluster in (0,1,2)
                    @test m.indices==findall(valid)
                    @test length(m.indices)==length(m.assignments)
                    @test all(c in (1,2) for c in m.assignments)
                    records=h.kind=="gmm" ? rec.g : rec.k
                    naming=[records[i].amplitude for i in m.indices]
                    means=[any(==(c),m.assignments) ? D.mean(naming[m.assignments .== c]) : NaN for c in 1:2]
                    if all(isfinite,means)
                        @test m.high_cluster in (1,2)
                        means[1]==means[2] || @test m.high_cluster==argmax(means)
                    else
                        @test m.high_cluster==0
                    end
                    if h.kind=="gmm"
                        @test all(isfinite,m.means) && all(isfinite,m.weights) && all(>(0),m.weights)
                        for c in m.covariances
                            @test all(isfinite,c) && isapprox(c,c';atol=1e-12,rtol=1e-12)
                            @test isposdef(Symmetric(c)+1e-8I)
                        end
                    else
                        @test all(isfinite,m.centers)
                    end
                end
            end
            for stage in D.REFERENCES
                @test read(joinpath(dir,stage*".tsv"))==read(joinpath(output,a,stage*".tsv"))
            end
        end
        mktempdir() do tmp
            for a in D.ARMS
                dir=joinpath(output,a.name); replay=joinpath(tmp,a.name); mkdir(replay)
                D.EF.write_scores(joinpath(replay,"fisher_cv.tsv"),D.frozen_fisher(patches[a.input],loaded[a.fisher].fisher))
                @test read(joinpath(dir,"fisher_cv.tsv"))==read(joinpath(replay,"fisher_cv.tsv"))
                ft=joinpath(replay,"features_predictor.tsv")
                D.joined(input,a.input,joinpath(replay,"fisher_cv.tsv"),ft,data.cfg)
                @test read(ft)==read(joinpath(dir,"features_predictor.tsv"))
                rec=D.head_records(ft,joinpath(input,a.input,"patches_bwd17.tsv"))
                opt=D.options(ft,joinpath(input,a.input,"patches_bwd17.tsv"),joinpath(replay,"unused.tsv"),data.cfgpath,data.cfg)
                p=D.predict_heads(rec,loaded[a.scales].bank,loaded[a.heads].bank,opt)
                D.Gaussian._write_predictions(joinpath(replay,"pred_gmm.tsv"),rec.g,p.g,Int.(isfinite.(p.g)))
                D.KMeans._write_predictions(joinpath(replay,"pred_kmeans.tsv"),rec.k,p.k,p.kviews)
                for stage in D.REFERENCES
                    D.Pipeline.check_counts(joinpath(dir,stage*".tsv"),data.counts)
                    stage in ("pred_gmm","pred_kmeans") &&
                        @test read(joinpath(dir,stage*".tsv"))==read(joinpath(replay,stage*".tsv"))
                end
                # Final vote is checked independently of the native vote writer.
                g=last(RU.lobe_table(joinpath(dir,"pred_gmm.tsv")))
                k=last(RU.lobe_table(joinpath(dir,"pred_kmeans.tsv")))
                v=last(RU.lobe_table(joinpath(dir,"predictions.tsv")))
                for key in keys(v)
                    @test v[key]["model"]==data.cfg["model"]["name"]
                    if g[key]["predicted"]=="?" || k[key]["predicted"]=="?"
                        @test v[key]["predicted"]=="?" && v[key]["probability_1"]=="NA"
                    else
                        pg=parse(Float64,g[key]["probability_1"]); pk=parse(Float64,k[key]["probability_1"])
                        @test 0<=pg<=1 && 0<=pk<=1
                        meanvote=(pg+pk)/2
                        @test isapprox(parse(Float64,v[key]["probability_1"]),meanvote;atol=5.1e-9,rtol=0)
                        @test isapprox(parse(Float64,v[key]["confidence"]),2abs(meanvote-.5);atol=5.1e-9,rtol=0)
                        @test v[key]["predicted"]==(meanvote>=data.cfg["selection"]["vote_threshold"] ? "1" : "0")
                    end
                end
            end
        end
        for (a,b) in (("control","local_only"),("observed","reverse_local_only"))
            x=last(RU.lobe_table(joinpath(output,a,"predictions.tsv"))); y=last(RU.lobe_table(joinpath(output,b,"predictions.tsv")))
            for k in keys(x)
                data.groups[k[1]]=="fully_observed" || continue
                @test x[k]==y[k]
            end
        end
    end
    nothing
end

function main(args=ARGS)
    args==["--help"] && return println("verify_frozen_learning.jl LOCAL_INPUT RUN_DIR (saved-state application only, no fitting)")
    length(args)==2 || error("Provide local input and one repetition directory")
    verify(args...)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && FrozenLearningVerification.main()
