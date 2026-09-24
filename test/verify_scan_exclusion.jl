#!/usr/bin/env julia
# Saved-state application only. This verifier does not fit any model or image.
module ScanExclusionVerification
using Test, TOML, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_scan_exclusion.jl"))
const S=ScanExclusionDiagnostic
const D=S.D
const RU=D.RU

function verify_bank(bank,train,opt,cfg)
    for h in (bank.g,values(bank.k)...)
        records=h.kind=="gmm" ? train.g : train.k
        states=[]
        if h.kind=="gmm"
            X,valid=D.Gaussian._standardized_matrix(records,h.features;normalization=opt.g.normalization,
                scale_fallback=opt.g.scale_fallback,interactions=opt.g.interactions,normalization_state=states)
        else
            X,valid=D.KMeans._standardized_matrix(records,h.features;interactions=opt.k.interactions,normalization_state=states)
        end
        @test isequal(h.scales,states)
        @test [m.seed for m in h.models]==collect(cfg["selection"]["first_seed"] .+
            (0:(h.kind=="gmm" ? cfg["selection"]["gmm_seeds"] : cfg["selection"]["kmeans_seeds"])-1))
        for m in h.models
            @test m.indices==findall(valid)
            @test length(m.indices)==length(m.assignments)
            @test all(c in (1,2) for c in m.assignments)
            # Independent training-membership application; no EM or kmeans fit.
            memberships=[h.kind=="gmm" ? argmax([D.Gaussian._final_component_score(X[i,:],m.means[:,c],
                m.covariances[c],m.weights[c],opt.g) for c in 1:2]) :
                argmin([sum(abs2,X[i,:]-m.centers[:,c]) for c in 1:2]) for i in m.indices]
            @test memberships==m.assignments
            amps=[records[i].amplitude for i in m.indices]
            means=[any(==(c),m.assignments) ? D.mean(amps[m.assignments .== c]) : NaN for c in 1:2]
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
end

function verify_vote(dir,cfg)
    g=last(RU.lobe_table(joinpath(dir,"pred_gmm.tsv")))
    k=last(RU.lobe_table(joinpath(dir,"pred_kmeans.tsv")))
    v=last(RU.lobe_table(joinpath(dir,"predictions.tsv")))
    RU.require_same_keys(g,k,"vote"); RU.require_same_keys(g,v,"vote")
    for key in keys(v)
        @test v[key]["model"]==cfg["model"]["name"]
        if g[key]["predicted"]=="?" || k[key]["predicted"]=="?"
            @test v[key]["predicted"]=="?" && v[key]["probability_1"]=="NA"
        else
            pg=parse(Float64,g[key]["probability_1"]); pk=parse(Float64,k[key]["probability_1"])
            @test 0<=pg<=1 && 0<=pk<=1
            vote=(pg+pk)/2
            @test isapprox(parse(Float64,v[key]["probability_1"]),vote;atol=5.1e-9,rtol=0)
            @test isapprox(parse(Float64,v[key]["confidence"]),2abs(vote-.5);atol=5.1e-9,rtol=0)
            @test v[key]["predicted"]==(vote>=cfg["selection"]["vote_threshold"] ? "1" : "0")
        end
    end
end

function verify(input,output)
    D.BLAS.set_num_threads(1)
    data=D.inputs(input); S.settings(joinpath(output,"settings.toml"))
    @testset "Complete scan-exclusion saved-output verification" begin
        @test TOML.parsefile(joinpath(output,"assignment.toml"))==data.cfg
        @test D.Pipeline.selected_counts(joinpath(output,"selected_summary.tsv"))==data.counts
        _,manifest=RU.read_table(joinpath(output,"inputs.tsv"))
        source=dirname(only(r["path"] for r in manifest if basename(r["path"])=="selected_summary.tsv"))
        @test length(manifest)==length(data.hashes)
        @test Dict(normpath(joinpath(input,relpath(r["path"],source)))=>r["sha256"] for r in manifest)==data.hashes
        _,groups=RU.read_table(joinpath(output,"groups.tsv"))
        @test length(groups)==length(data.counts)
        @test Dict(r["file"]=>r["group"] for r in groups)==data.groups
        _,ar=RU.read_table(joinpath(output,"arm.tsv")); arm=only(ar)["arm"]
        @test arm in ("control","observed")
        _,folds=RU.read_table(joinpath(output,"folds.tsv"))
        files=sort(collect(keys(data.counts)))
        @test [r["file"] for r in folds]==files
        @test [r["fold"] for r in folds]==[D.@sprintf("fold%04d",i) for i in eachindex(files)]
        opts=D.EF.load_fisher_config(data.cfg)
        patches=D.EF.load_patches(joinpath(input,arm,"patches_fwd17.tsv"),"res",opts)
        jobs=vcat([Dict("fold"=>"native","file"=>"")],folds)
        mktempdir() do tmp
            for fold in jobs
                excluded=fold["file"]; dir=joinpath(output,fold["fold"])
                replay=joinpath(tmp,fold["fold"]); mkdir(replay)
                training=S.training_patches(patches,excluded)
                if !isempty(excluded)
                    @test all(k[1]!=excluded for k in training.keys)
                    @test parse(Int,fold["training_scans"])==length(files)-1==length(unique(first.(training.keys)))
                    @test parse(Int,fold["training_lobes"])==length(training.keys)==sum(values(data.counts))-data.counts[excluded]
                    @test parse(Int,fold["target_lobes"])==data.counts[excluded]
                end
                saved=D.load_bank(joinpath(dir,"models.toml"))
                valid=isempty.(training.invalid_reasons) .& [all(isfinite,x) for x in eachrow(training.X)] .& isfinite.(training.amplitudes)
                for parity in 0:1
                    m=saved.fisher[parity]
                    @test m.train==findall(i->valid[i] && mod(training.keys[i][2],2)==parity,eachindex(valid))
                    @test m.held==findall(i->valid[i] && mod(training.keys[i][2],2)!=parity,eachindex(valid))
                    @test all(isfinite,m.w_p) && all(isfinite,m.mid)
                    @test all(training.keys[i][1]!=excluded for i in m.train)
                end
                D.EF.write_scores(joinpath(replay,"fisher_cv.tsv"),D.frozen_fisher(patches,saved.fisher))
                @test read(joinpath(replay,"fisher_cv.tsv"))==read(joinpath(dir,"fisher_cv.tsv"))
                ft=joinpath(replay,"features_predictor.tsv")
                D.joined(input,arm,joinpath(replay,"fisher_cv.tsv"),ft,data.cfg)
                @test read(ft)==read(joinpath(dir,"features_predictor.tsv"))
                pp=joinpath(input,arm,"patches_bwd17.tsv"); rec=D.head_records(ft,pp)
                train,target=S.split_records(rec,excluded)
                opt=D.options(ft,pp,joinpath(replay,"unused.tsv"),data.cfgpath,data.cfg)
                verify_bank(saved.bank,train,opt,data.cfg)
                if !isempty(excluded)
                    @test all(r.file!=excluded for r in train.g) && all(r.file!=excluded for r in train.k)
                    @test all(r.file==excluded for r in target.g) && all(r.file==excluded for r in target.k)
                end
                header,table=RU.lobe_table(ft)
                RU.write_table(joinpath(replay,"target_features.tsv"),header,[table[D.key(r)] for r in target.g])
                @test read(joinpath(replay,"target_features.tsv"))==read(joinpath(dir,"target_features.tsv"))
                scales=S.local_scaling(target,opt)
                @test S.scaling_rows(scales)==last(RU.read_table(joinpath(dir,"target_scales.tsv")))
                S.write_target(replay,target,scales,saved.bank,opt,data.cfgpath)
                for stage in S.STAGES
                    @test read(joinpath(replay,stage*".tsv"))==read(joinpath(dir,stage*".tsv"))
                    counts=isempty(excluded) ? data.counts : Dict(excluded=>data.counts[excluded])
                    D.Pipeline.check_counts(joinpath(dir,stage*".tsv"),counts)
                end
                verify_vote(dir,data.cfg)
                if isempty(excluded)
                    for stage in D.REFERENCES
                        @test read(joinpath(dir,stage*".tsv"))==read(joinpath(input,arm,stage*".tsv"))
                    end
                end
            end
            S.merge_predictions(tmp,folds,data.counts)
            for stage in S.STAGES
                @test read(joinpath(tmp,stage*".tsv"))==read(joinpath(output,stage*".tsv"))
            end
        end
    end
    nothing
end

function main(args=ARGS)
    args==["--help"] && return println("verify_scan_exclusion.jl LOCAL_INPUT ONE_ARM_RUN (saved-state application only)")
    length(args)==2 || error("Provide local input and one arm/repetition directory")
    verify(args...)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ScanExclusionVerification.main()
