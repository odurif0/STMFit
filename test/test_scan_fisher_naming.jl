# Synthetic only: no benchmark, raw scans or historical predictions.
using Test, Random, Statistics, LinearAlgebra, TOML, SHA, Printf
include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
using .EmpiricalFisherNative
const EF = EmpiricalFisherNative
const RU = EF.ReconstructedUnitAssignment
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
module Attribution
include(joinpath(@__DIR__, "diagnose_fisher_attribution.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT,"config","unit_assignment_patch_support.toml")
const SCAN = joinpath(ROOT,"config","unit_assignment_scan_fisher.toml")
const NAMING = joinpath(ROOT,"config","unit_assignment_relative_naming.toml")
const OPTIONS = load_fisher_config(CONTROL)
const GROUPED = load_fisher_config(SCAN)
const GRID = fisher_grid(OPTIONS)

@testset "Two isolated explicit variants, no combined or label-directed rule" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    base = RU.load_config(CONTROL)
    for (path,field,value) in ((SCAN,"fisher_cv_scheme","scan_hash_twofold"),
                              (NAMING,"gmm_cluster_naming","within_scan_z"))
        candidate = RU.load_config(path)
        @test candidate["selection"][field] == value
        candidate["selection"][field] = base["selection"][field]
        candidate["model"]["name"] = base["model"]["name"]
        @test candidate == base
    end
    for name in readdir(joinpath(ROOT,"config"))
        startswith(name,"unit_assignment_") && endswith(name,".toml") || continue
        c = TOML.parsefile(joinpath(ROOT,"config",name))
        haskey(get(c,"model",Dict()),"gmm_covariance_structure") || continue
        @test RU.load_fisher_cv(c).scheme == (name==basename(SCAN) ? "scan_hash_twofold" : "lobe_parity")
        @test RU.load_gmm_cluster_naming(c) == (name==basename(NAMING) ? "within_scan_z" : "raw_amplitude")
    end
    for field in ("fisher_cv_scheme","fisher_scan_split_seed","gmm_cluster_naming")
        invalid = field == "fisher_scan_split_seed" ? (nothing,-1,true,0.5,Inf,"0") : (nothing,"auto",true,1)
        for value in invalid
            c=deepcopy(base)
            value===nothing ? delete!(c["selection"],field) : (c["selection"][field]=value)
            if field=="gmm_cluster_naming"
                @test_throws ArgumentError RU.load_gmm_cluster_naming(c)
            else
                @test_throws ArgumentError RU.load_fisher_cv(c)
                @test_throws ArgumentError load_fisher_config(c)
            end
        end
    end
    c=RU.load_config(NAMING); c["selection"]["fisher_cv_scheme"]="scan_hash_twofold"
    @test_throws ArgumentError RU.load_gmm_cluster_naming(c)
    for (section,field,value) in (("selection","gmm_resampling","whole_scans"),
            ("selection","gmm_training_weighting","equal_scans"),
            ("model","gmm_covariance_structure","tied"))
        c=RU.load_config(NAMING); c[section][field]=value
        @test_throws ArgumentError RU.load_gmm_cluster_naming(c)
    end
end

function fisher_fixture()
    rng=MersenneTwister(923)
    rowkeys=[(@sprintf("synthetic_%02d.sxm",f),j) for f in 1:12 for j in 1:8]
    full=0.07randn(rng,length(rowkeys),GRID.side^2)
    shape=[exp(-((u-9)^2+(t-9)^2)/18)+.03(t-9) for u in 1:17 for t in 1:17]
    for (i,(file,j)) in enumerate(rowkeys)
        full[i,:] .+= (j<=6 ? -2.0 : 2.0).*shape
    end
    table=PatchTable(rowkeys,full[:,GRID.disk_indices],full[:,GRID.center_index],fill("",length(rowkeys)),GRID)
    table,full
end

@testset "Whole scans have balanced deterministic groups, independent of row contents" begin
    patches,_=fisher_fixture(); rowkeys=patches.keys
    groups=fisher_fold_ids(rowkeys,GROUPED)
    files=sort(unique(first.(rowkeys)))
    # Independent spelling of the documented byte payload and alternating rank rule.
    digests=Dict(f=>bytes2hex(sha256(vcat(collect(codeunits("0")),UInt8[0],collect(codeunits(f))))) for f in files)
    ordered=sort(files;by=f->(digests[f],f))
    expected=Dict(ordered[i]=>isodd(i) ? 0 : 1 for i in eachindex(ordered))
    @test groups == [expected[k[1]] for k in rowkeys]
    @test length(Set(k[1] for (i,k) in enumerate(rowkeys) if groups[i]==0))==6
    @test length(Set(k[1] for (i,k) in enumerate(rowkeys) if groups[i]==1))==6
    @test all(length(unique(groups[findall(k->k[1]==f,rowkeys)]))==1 for f in files)
    @test fisher_fold_ids(reverse(rowkeys),GROUPED)==reverse(groups)
    @test fisher_fold_ids([(k[1],100-k[2]) for k in rowkeys],GROUPED)==groups
    @test fisher_fold_ids([(joinpath("/unrelated",k[1]),k[2]) for k in rowkeys],GROUPED)==groups
    @test fisher_fold_ids(rowkeys,OPTIONS)==[mod(k[2],2) for k in rowkeys]
    @test isempty(fisher_fold_ids(Tuple{String,Int}[],GROUPED))
    # Removing lobes, but not a scan identity, cannot alter another lobe's group.
    subset=collect(1:3:length(rowkeys))
    @test Set(first.(rowkeys[subset]))==Set(files)
    @test fisher_fold_ids(rowkeys[subset],GROUPED)==groups[subset]
    Random.seed!(453); rand(30)
    @test fisher_fold_ids(rowkeys,GROUPED)==groups
end

@testset "All Fisher learning uses only opposite scans; legacy arithmetic unchanged" begin
    table,_=fisher_fixture(); groups=fisher_fold_ids(table.keys,GROUPED)
    diagnostics=[]; rows=cv_scores(table,GROUPED;diagnostics)
    @test length(rows)==length(table.keys) && all(isempty(r.invalid_reason) for r in rows)
    @test length(diagnostics)==2
    for (fold,d) in enumerate(diagnostics)
        group=fold-1
        train=findall(==(group),groups); held=findall(!=(group),groups)
        @test d.train==train && d.held==held && d.status=="fitted"
        @test isempty(intersect(Set(first.(table.keys[train])),Set(first.(table.keys[held]))))
        independent=fit_fisher(table.X[train,:],table.amplitudes[train],GROUPED)
        @test d.model.w_p==independent.w_p && d.model.mid==independent.mid
        @test d.model.gmm.weights==independent.gmm.weights
        @test [rows[i].score for i in held]==[maxmirror_score(table.X[i,:],independent,GRID) for i in held]
        # Every held-out pixel and naming amplitude can change without affecting this model.
        altered=deepcopy(table.X); altered[held,:] .+= 0.13
        amplitudes=copy(table.amplitudes); amplitudes[held] .+= 20
        other=[]; cv_scores(PatchTable(table.keys,altered,amplitudes,table.invalid_reasons,GRID),GROUPED;diagnostics=other)
        @test other[fold].model.w_p==d.model.w_p && other[fold].model.mid==d.model.mid
        @test other[fold].model.amplitude_means==d.model.amplitude_means
    end
    legacy=cv_scores(table,OPTIONS)
    for parity in 0:1
        train=findall(k->mod(k[2],2)==parity,table.keys)
        held=findall(k->mod(k[2],2)!=parity,table.keys)
        independent=fit_fisher(table.X[train,:],table.amplitudes[train],OPTIONS)
        @test [legacy[i].score for i in held]==[maxmirror_score(table.X[i,:],independent,GRID) for i in held]
    end
    order=randperm(MersenneTwister(22),length(rows))
    permuted=PatchTable(table.keys[order],table.X[order,:],table.amplitudes[order],table.invalid_reasons[order],GRID)
    @test [r.score for r in cv_scores(permuted,GROUPED)]==[rows[i].score for i in order]
    invalid=copy(table.invalid_reasons); invalid[1]="missing_pixel"
    bad=cv_scores(PatchTable(table.keys,table.X,table.amplitudes,invalid,GRID),GROUPED)
    @test bad[1].invalid_reason=="missing_pixel" && isnan(bad[1].score)
    @test [r.file=>r.lobe for r in bad]==[k[1]=>k[2] for k in table.keys]
    one=findall(k->k[1]==table.keys[1][1],table.keys)
    absent=cv_scores(PatchTable(table.keys[one],table.X[one,:],table.amplitudes[one],fill("",length(one)),GRID),GROUPED)
    @test length(absent)==length(one)
    @test all(r->startswith(r.invalid_reason,"training_scan_") && isnan(r.score),absent)
    diagnostic=Attribution.FisherAttributionDiagnostic
    own_options=diagnostic.EF.load_fisher_config(SCAN)
    own_table=diagnostic.EF.PatchTable(table.keys,table.X,table.amplitudes,table.invalid_reasons,
        diagnostic.EF.fisher_grid(own_options))
    @test_throws ArgumentError diagnostic.replay_folds(own_table,own_options)
    @test_throws ErrorException diagnostic.RD.load_inputs(Dict(),SCAN)
    @test_throws ErrorException diagnostic.RD.load_inputs(Dict(),NAMING)
end

@testset "Relative naming removes scan offsets/gains, never fixes composition" begin
    amplitudes=[100.,101.,102.,1.,2.,3.]
    records=[GMM.LobeRecord(i<=3 ? "a.sxm" : "b.sxm",mod1(i,3),amplitudes[i],Dict{String,Float64}()) for i in 1:6]
    assigned=[1,1,2,1,2,2]; eligible=trues(6)
    raw=GMM._cluster_amplitude_means(records,1:6,assigned,eligible,nothing)
    relative=GMM._cluster_amplitude_means(records,1:6,assigned,eligible,nothing;mode="within_scan_z",scale_fallback=1.)
    @test raw==Dict(1=>(100+101+1)/3,2=>(102+2+3)/3)
    @test relative[1]≈-2/3 && relative[2]≈2/3
    @test raw[1]>raw[2] && relative[1]<relative[2] # Names can reverse without moving a member.
    transformed=deepcopy(records)
    for i in 1:6; transformed[i].amplitude=(i<=3 ? 2.0 : 3.0)*amplitudes[i]+(i<=3 ? -50.0 : 80.0); end
    renamed=GMM._cluster_amplitude_means(transformed,1:6,assigned,eligible,nothing;mode="within_scan_z",scale_fallback=1.)
    @test all(renamed[k]≈relative[k] for k in keys(relative))
    @test GMM._cluster_amplitude_means(records,1:6,ones(Int,6),eligible,nothing;mode="within_scan_z",scale_fallback=1.)==Dict(1=>0.)
    extra=vcat(records,[GMM.LobeRecord("a.sxm",4,1e9,Dict{String,Float64}())])
    @test GMM._cluster_amplitude_means(extra,1:7,vcat(assigned,2),vcat(eligible,false),nothing;mode="within_scan_z",scale_fallback=1.)==relative
    constant=[GMM.LobeRecord("constant.sxm",i,3.,Dict{String,Float64}()) for i in 1:2]
    @test GMM._cluster_amplitude_means(constant,1:2,[1,2],trues(2),nothing;mode="within_scan_z",scale_fallback=1.)==Dict(1=>0.,2=>0.)
    @test GMM._cluster_amplitude_means(constant,1:1,[1],trues(2),nothing;mode="within_scan_z",scale_fallback=1.)==Dict(1=>0.)
    @test_throws ArgumentError GMM._cluster_amplitude_means(records,1:6,assigned,eligible,nothing;mode="within_scan_z")
    @test_throws ArgumentError GMM._cluster_amplitude_means(records,1:6,assigned,eligible,ones(6);mode="within_scan_z",scale_fallback=1.)
end

@testset "CLI consistency, fixed GMM fitting and label-free inputs" begin
    mktempdir() do dir
        table,full=fisher_fixture()
        patchpath=joinpath(dir,"patches.tsv")
        open(patchpath,"w") do io
            println(io,join(vcat(["file","lobe"],[@sprintf("res_p%03d",i) for i in 1:289]),'\t'))
            for (i,k) in enumerate(table.keys); println(io,join(vcat([k[1],string(k[2])],string.(full[i,:])),'\t')); end
        end
        expected=joinpath(dir,"expected.tsv"); actual=joinpath(dir,"actual.tsv")
        loaded=load_patches(patchpath,"res",GROUPED)
        write_scores(expected,cv_scores(loaded,GROUPED))
        command=`$(Base.julia_cmd()) --project=$ROOT $(joinpath(ROOT,"test","build_empirical_fisher_native.jl")) --patches $patchpath --prefix res --config $SCAN --out $actual`
        @test success(command)
        @test read(actual)==read(expected)
        @test !success(pipeline(ignorestatus(`$command --truth forbidden.tsv`);stdout=devnull,stderr=devnull))
        path=joinpath(dir,"features.tsv"); rng=MersenneTwister(922)
        records=GMM.LobeRecord[]
        open(path,"w") do io
            println(io,"file\tlobe\tamplitude\tsignal\tshape")
            for f in 1:8, j in 1:10
                signal=(j<=7 ? -2.0 : 2.0)+.2randn(rng)
                shape=.1randn(rng)+signal/3
                amp=10.0*f+0.5*signal
                rec=GMM.LobeRecord("synthetic_$f.sxm",j,amp,Dict("signal"=>signal,"shape"=>shape))
                push!(records,rec); println(io,join((rec.file,j,amp,signal,shape),'\t'))
            end
            println(io,"synthetic_8.sxm\t11\t80.0\tNA\t0.0")
        end
        options=[GMM._parse_cli(["--features",path,"--config",c,"--out",joinpath(dir,"gmm_$i.tsv"),
            "--view","synthetic=signal,shape","--seeds","2","--selftrain","2"]) for (i,c) in enumerate((CONTROL,NAMING))]
        loaded_records=GMM._load_records(path); fits=[[],[]]
        probs=[GMM._view_probability(loaded_records,["signal","shape"],o;diagnostics=fits[i]) for (i,o) in enumerate(options)]
        @test fits[1]==fits[2] # Every mean, covariance, membership and free weight, all seeds.
        @test isequal(probs[1][end],NaN) && isequal(probs[2][end],NaN)
        @test all(x->isnan(x)||0<=x<=1,probs[2])
        @test [o.cluster_naming for o in options]==["raw_amplitude","within_scan_z"]
        @test options[1].covariance_structure==options[2].covariance_structure=="full"
        for (i,c) in enumerate((CONTROL,NAMING))
            out=options[i].out_tsv
            command=`$(Base.julia_cmd()) --project=$ROOT $(joinpath(ROOT,"test","build_labelfree_gmm_predictions.jl")) --features $path --config $c --out $out --view synthetic=signal,shape --seeds 2 --selftrain 2`
            @test success(command)
            _,rows=RU.lobe_table(out)
            @test length(rows)==81 && rows[("synthetic_8.sxm",11)]["predicted"]=="?"
            for (j,r) in enumerate(loaded_records)
                isnan(probs[i][j]) && continue
                @test parse(Float64,rows[(r.file,r.lobe)]["probability_1"])≈probs[i][j]
            end
        end
    end
end
