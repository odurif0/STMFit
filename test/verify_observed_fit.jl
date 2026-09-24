#!/usr/bin/env julia
# Saved-output verification/report only. No fit, truth, labels or grading.
module ObservedFitVerification
include(joinpath(@__DIR__,"run_observed_fit_comparison.jl"))
using Test, TOML, SHA, Statistics
const C=ObservedFitComparison
const O=C.ObservedGeometry
const P=C.P
const G=O.G
num(r,k)=parse(Float64,r[k])
int(r,k)=parse(Int,r[k])
table(path)=isfile(path) ? last(P.read_table(path)) : Dict{String,String}[]
close(a,b)=isapprox(a,b;atol=1e-12,rtol=1e-10)

function verify(root; rawdir=nothing)
    counts=P.selected_counts(joinpath(root,"selected_summary.tsv"))
    raw=TOML.parsefile(joinpath(root,"control_count.toml"))
    assignment=TOML.parsefile(joinpath(root,"assignment.toml"))
    settings=O.settings(joinpath(root,"settings.toml"))
    observed=TOML.parsefile(joinpath(root,"observed_count.toml"))
    prefixes=sort(unique([chop(path;tail=length(suffix)) for path in readdir(root;join=true)
        for suffix in (".fits.tsv",".failures.tsv") if startswith(basename(path),"refit") && endswith(path,suffix)]))
    seen=Set{String}(); fitsall=Dict{String,String}[]; failuresall=Dict{String,String}[]
    commonall=Dict{String,String}[]; diagnosticsall=Dict{String,String}[]
    @testset "Observed fits: complete attempts, frozen settings, independent residuals and GCV" begin
        @test observed==C.candidate_config(raw,settings)
        @test !isempty(prefixes)
        hashes=Dict(r["input"]=>r["sha256"] for r in table(joinpath(root,"input_hashes.tsv")))
        for (key,name) in (("--selected-summary","selected_summary.tsv"),("--count-config","control_count.toml"),
                ("--config","assignment.toml"),("--settings","settings.toml"),("--templates","templates.tsv"))
            @test bytes2hex(sha256(read(joinpath(root,name))))==hashes[key]
        end
        rh=table(joinpath(root,"raw_hashes.tsv"))
        @test Set(r["file"] for r in rh)==Set(keys(counts)) && length(rh)==length(counts)
        if rawdir!==nothing
            paths=P.raw_index(rawdir)
            @test Set(keys(paths))==Set(keys(counts))
            @test all(bytes2hex(sha256(read(paths[r["file"]])))==r["sha256"] for r in rh)
        end
        for prefix in prefixes
            fits=table(prefix*".fits.tsv"); failures=table(prefix*".failures.tsv")
            parameters=table(prefix*".parameters.tsv"); families=table(prefix*".families.tsv")
            common=table(prefix*".common.tsv"); diagnostics=table(prefix*".diagnostics.tsv")
            append!(fitsall,fits); append!(failuresall,failures); append!(commonall,common); append!(diagnosticsall,diagnostics)
            files=union(Set(r["file"] for r in fits),Set(r["file"] for r in failures))
            @test isempty(intersect(seen,files)); union!(seen,files)
            @test length(Set((r["file"],r["arm"],r["profile"]) for r in fits))==length(fits)
            @test length(parameters)==sum(int(r,"p_full") for r in fits;init=0)
            @test length(families)==2length(fits)
            for f in files, arm in O.ARMS, profile in O.PROFILES
                successes=[r for r in fits if (r["file"],r["arm"],r["profile"])==(f,arm,profile)]
                reasons=[r for r in failures if r["file"]==f && r["stage"] in (arm*"_"*profile,arm*"_support","load_or_comparison")]
                @test length(successes)==1 || !isempty(reasons)
            end
            tables=Dict((a,p)=>last(P.lobe_table(prefix*".$a.$p.tsv")) for a in O.ARMS for p in O.PROFILES if isfile(prefix*".$a.$p.tsv"))
            @test isempty(failures) ? length(tables)==4 : isempty(tables)
            for t in values(tables)
                @test Set(keys(t))==Set((f,i) for f in files for i in 1:counts[f])
            end
            for a in fits
                file,arm,profile=a["file"],a["arm"],a["profile"]; n=counts[file]
                @test int(a,"N")==n && a["valid"]=="true" && a["reason"]=="ok"
                model=deepcopy(raw["model"]); model["peak_profile"]=profile
                profile=="split" && (model["skew_ratio_max"]=assignment["model"]["split_skew_ratio_max"])
                _,cfg,_=O.F.Extractor._configs(model,raw["preprocessing"],"unused"); cfg.chain_circular_sigmas=a["source"]=="circ"
                pars=sort([r for r in parameters if (r["file"],r["arm"],r["profile"])==(file,arm,profile)];by=r->int(r,"parameter"))
                p=num.(pars,"value"); lo,hi=O.VP.native_raw_bounds(n,cfg)
                @test int.(pars,"parameter")==collect(eachindex(p)) && all(lo.<=p.<=hi)
                @test num.(pars,"lower")==lo && num.(pars,"upper")==hi
                @test length(p)==int(a,"p_full")==G._chain_nparams(n,cfg)
                d=table(joinpath(prefix*".fit_data",file*".$arm.tsv"))
                x=num.(d,"x_nm"); y=num.(d,"y_nm"); z=num.(d,"z_nm")
                @test length(z)==int(a,"n_data")>length(p) && all(isfinite,vcat(x,y,z))
                @test length(Set((int(r,"row"),int(r,"column")) for r in d))==length(d)
                axis=(origin=(num(a,"origin_x_nm"),num(a,"origin_y_nm")),axis=[num(a,"axis_x"),num(a,"axis_y")],
                    perp=[-num(a,"axis_y"),num(a,"axis_x")],tmin=num(a,"support_tmin"),tmax=num(a,"support_tmax"))
                b,feats,ts,us,sp,sq=G._decode_chain(p,n,axis,cfg;amp_min=num(a,"amp_min"),amp_range=num(a,"amp_range"))
                @test close(sum(abs2,axis.axis),1.)
                @test all(cfg.sigma_parallel_min_nm.<=sp.<=cfg.sigma_parallel_max_nm)
                @test all(cfg.sigma_perp_min_nm.<=sq.<=cfg.sigma_perp_max_nm)
                @test all(abs.(us).<=cfg.lateral_max_nm)
                @test all(ts.>=axis.tmin-1e-6) && all(ts.<=axis.tmax+1e-6)
                @test all(diff(ts).>=G._effective_spacing_min_nm(cfg)-1e-12) && all(diff(ts).<=cfg.spacing_max_nm+1e-12)
                pred=b.+p[2].*x.+p[3].*y
                for f in feats
                    dt=(x.-f.x_nm).*axis.axis[1].+(y.-f.y_nm).*axis.axis[2]
                    du=(x.-f.x_nm).*(-axis.axis[2]).+(y.-f.y_nm).*axis.axis[1]
                    widths=profile=="split" ? ifelse.(dt.<0,f.sigma_x_nm/sqrt(f.skew_ratio),f.sigma_x_nm*sqrt(f.skew_ratio)) : fill(f.sigma_x_nm,length(dt))
                    pred .+= f.amplitude .* exp.(-0.5 .* ((dt ./ widths).^2 .+ (du ./ f.sigma_y_nm).^2))
                    @test num(a,"amp_min")<=f.amplitude<=num(a,"amp_min")+num(a,"amp_range")
                    @test 1/cfg.skew_ratio_max<=f.skew_ratio<=cfg.skew_ratio_max
                end
                rss=sum(abs2,pred.-z)
                @test close(rss,num(a,"rss")) && close(length(z)/(length(z)-length(p))^2*rss,num(a,"gcv"))
                peak=maximum(abs.(pred.-z))/max(num(a,"noise"),G.EPS)
                @test close(peak,num(a,"residual_peak_snr")) && peak<=cfg.residual_peak_snr_threshold
                bs=[r for r in families if (r["file"],r["arm"],r["profile"])==(file,arm,profile)]
                @test Set(r["family"] for r in bs)==Set(["circ","ell"])
                @test num(a,"gcv")==minimum(num(r,"gcv") for r in bs if r["valid"]==r["success"]=="true")
                r=G.ChainModelResult(n=n,params=p,amp_min=num(a,"amp_min"),amp_range=num(a,"amp_range"),gcv=num(a,"gcv"))
                serialized=O.R.feature_rows(file,r,(axisctx=axis,),cfg,a["source"])
                haskey(tables,(arm,profile)) && @test all(tables[(arm,profile)][(file,i)]==serialized[i] for i in 1:n)
                reused=arm=="observed" && int(a,"observed_pixels")==int(a,"total_pixels")
                @test (a["reused_control"]=="true")==reused
                if reused
                    c=only(r for r in fits if (r["file"],r["arm"],r["profile"])==(file,"control",profile))
                    @test all(a[k]==c[k] for k in setdiff(keys(a),["arm","reused_control","elapsed_s"]))
                    @test a["elapsed_s"]=="0"
                    @test d==table(joinpath(prefix*".fit_data",file*".control.tsv"))
                    @test all(r["arm"]!="observed" for r in diagnostics if r["file"]==file)
                else
                    @test any((r["file"],r["arm"],r["profile"])==(file,arm,profile) for r in diagnostics)
                end
            end
            for r in common
                d=table(joinpath(prefix*".fit_data",r["file"]*".$(r["profile"]).common.tsv"))
                @test length(d)==int(r,"n_common")>0
                @test all(isfinite,num.(d,"raw_mean_nm"))
                coords=Set((v["x_nm"],v["y_nm"]) for v in d)
                maps=Dict(arm=>Dict((v["x_nm"],v["y_nm"])=>v for v in
                    table(joinpath(prefix*".fit_data",r["file"]*".$arm.tsv"))) for arm in O.ARMS)
                @test coords==intersect(Set(keys(maps["control"])),Set(keys(maps["observed"])))
                for arm in O.ARMS
                    @test all(v[arm*"_target_nm"]==maps[arm][(v["x_nm"],v["y_nm"])]["z_nm"] for v in d)
                    @test close(sum((num(v,arm*"_target_nm")-num(v,arm*"_model_nm"))^2 for v in d),num(r,arm*"_rss"))
                end
            end
        end
        @test seen==Set(keys(counts))
        if !isempty(failuresall)
            @test isfile(joinpath(root,"failures.tsv"))
            @test all(!ispath(joinpath(root,a)) for a in O.ARMS)
            @test !isfile(joinpath(root,"refit.tsv"))
        else
            @test length(fitsall)==4length(counts)
            @test !isfile(joinpath(root,"failures.tsv"))
            for arm in O.ARMS, name in ("features","features_split","features_local","features_descriptor","features_predictor",
                    "patches_fwd17","patches_bwd17","patches_bwd9","training_support","score_fwd","score_bwd","fisher_cv","pred_gmm","pred_kmeans","predictions")
                @test length(P.check_counts(joinpath(root,arm,name*".tsv"),counts))==sum(values(counts))
            end
        end
    end
    return (;counts,fits=fitsall,failures=failuresall,common=commonall,diagnostics=diagnosticsall)
end

function report(r)
    io=IOBuffer(); n=length(r.counts)
    println(io,"# Observed-pixel fixed-N comparison\n\nSaved-output checks passed. No benchmark labels or new fit were used by this verifier.\n")
    println(io,"- Input cohort: $n scans, $(sum(values(r.counts))) predicted lobes.\n- Valid fits: $(length(r.fits))/$(4n).\n- Failure records: $(length(r.failures)).")
    for arm in O.ARMS
        fs=[f for f in keys(r.counts) if all(any(a["file"]==f && a["arm"]==arm && a["profile"]==p for a in r.fits) for p in O.PROFILES)]
        println(io,"- $arm: $(length(fs))/$n scans have both valid profiles.")
    end
    println(io,"- Exact fit reuses: $(count(a->a["reused_control"]=="true",r.fits)).")
    if !isempty(r.failures)
        println(io,"\nIncomplete cohort: no downstream assignment, partial benchmark or promotion.\n\n| File | Stage | Reason |\n|---|---|---|")
        for a in r.failures; println(io,"| $(a["file"]) | $(a["stage"]) | $(replace(a["reason"],'|'=>'/')) |"); end
    else
        println(io,"\nBoth complete assignment tables are available for separate external grading. This is not a reproduction or promotion.")
    end
    println(io,"\n## Diagnostics, not recognition\n")
    for profile in O.PROFILES
        rows=[a for a in r.common if a["profile"]==profile]
        improved=count(a->num(a,"observed_rss")<num(a,"control_rss"),rows)
        identical=count(a->num(a,"observed_rss")==num(a,"control_rss"),rows)
        println(io,"- $profile common-pixel pairs: $(length(rows)); lower candidate RSS $improved, identical $identical, higher $(length(rows)-improved-identical). Unavailable pairs are not counted as improvements.")
    end
    for arm in O.ARMS
        ds=[a for a in r.diagnostics if a["arm"]==arm]
        println(io,"- $arm optimizer records: $(length(ds)); LM converged $(count(a->a["lm_converged"]=="true",ds)); reused candidates have no second optimization record.")
    end
    println(io,"\nPhysical/chemical settings and predicted N are frozen. Native RSS/GCV supports differ; paired RSS includes each arm's background on common observed pixels. No lower-RSS chemical claim. Inherited benchmark-informed calibration is not relabeled strictly label-free.")
    String(take!(io))
end
function main(args=ARGS)
    if isempty(args) || "--help" in args
        return println("verify_observed_fit.jl RUN_DIR NEW_REPORT.md [RAW_DIR]")
    end
    length(args) in (2,3) || error("Expected run, new report and optional raw directory")
    (ispath(args[2]) || islink(args[2])) && error("Report exists")
    r=verify(args[1];rawdir=length(args)==3 ? args[3] : nothing)
    write(args[2],report(r)); println(report(r))
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ObservedFitVerification.main()
