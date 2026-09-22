#!/usr/bin/env julia
# Saved-only checks; no SXM pixels, optimizer, benchmark truth or external grade.
module RegisteredRefitVerification
using Test, TOML, SHA, Statistics, Printf
include(joinpath(@__DIR__,"refit_registered_geometry.jl"))
using .RegisteredRefit: lobe_table, read_table
const RR=RegisteredRefit
const G=RR.G
const ROOT=dirname(@__DIR__)
num(r,k)=parse(Float64,r[k])
int(r,k)=parse(Int,r[k])
close(a,b)=isapprox(a,b;atol=1e-12,rtol=1e-10)

function geometry(root,cache,shiftpath; complete=true)
    _,base=lobe_table(cache)
    shifts=RR.read_shifts(shiftpath,unique(first.(collect(keys(base)))))
    expected=Set(keys(base)); seen=Set{String}()
    raw=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
    assignment=TOML.parsefile(joinpath(ROOT,"config","unit_assignment_patch_support.toml"))
    paths=sort(unique([chop(p;tail=length(suffix)) for p in readdir(root;join=true)
        for suffix in (".fits.tsv",".failures.tsv") if endswith(p,suffix)]))
    combinations=[(a,p) for a in RR.ARMS for p in RR.PROFILES]
    settings_path=joinpath(root,"refit_settings.toml")
    frozen=isfile(settings_path) && RR.settings(settings_path).original_support
    completed_before_failure=Dict("load"=>0,"original_support"=>0,"control_support"=>0,"control_gaussian"=>0,"control_split"=>1,
        "registered_support"=>2,"registered_gaussian"=>2,"registered_split"=>3)
    optional_table(path)=isfile(path) ? last(read_table(path)) : Dict{String,String}[]
    @testset "Observed native refits: parameters, independent RSS/GCV, family choice and zero-shift reuse" begin
        @test !isempty(paths)
        for prefix in paths
            audits=optional_table(prefix*".fits.tsv"); parameters=optional_table(prefix*".parameters.tsv"); bootstrap=optional_table(prefix*".bootstrap.tsv")
            failures=optional_table(prefix*".failures.tsv")
            supports=Dict(r["file"]=>r for r in optional_table(prefix*".support.tsv"))
            frozen || @test isempty(supports)
            complete && @test isempty(failures)
            failed=Dict(r["file"]=>r for r in failures)
            @test length(failed)==length(failures)
            tables=Dict((a,p)=>last(lobe_table(prefix*".$a.$p.tsv")) for a in RR.ARMS for p in RR.PROFILES if isfile(prefix*".$a.$p.tsv"))
            @test isempty(failures) ? length(tables)==4 : isempty(tables)
            files=union(Set(r["file"] for r in audits),Set(keys(failed)))
            @test length(bootstrap)==2length(audits)
            @test Set((r["file"],r["arm"],r["profile"]) for r in parameters)==Set((r["file"],r["arm"],r["profile"]) for r in audits)
            @test length(parameters)==sum(int(a,"p_full") for a in audits;init=0)
            for f in files
                expected_slots=haskey(failed,f) ? combinations[1:completed_before_failure[failed[f]["stage"]]] : combinations
                actual=[(r["arm"],r["profile"]) for r in audits if r["file"]==f]
                @test actual==expected_slots
                if frozen && !isempty(actual)
                    @test haskey(supports,f)
                    support=supports[f]
                    original=last(read_table(joinpath(prefix*".fit_data",f*".original.tsv")))
                    ox=num.(original,"x_nm"); oy=num.(original,"y_nm"); oz=num.(original,"z_native_nm")
                    _,cfg,_=RR.F.Extractor._configs(raw["model"],raw["preprocessing"],"unused")
                    full=G._weighted_roi_axis(ox,oy,oz)
                    _,_,_,axis,keep,meta=G._chain_fit_data(ox,oy,oz,full,cfg)
                    @test keep==[r["in_fit"]=="true" for r in original]
                    @test count(keep)==int(support,"fit_pixels") && length(original)==int(support,"roi_pixels")
                    @test meta.support_method==support["support_method"]
                    for (key,value) in (("axis_x",axis.axis[1]),("axis_y",axis.axis[2]),
                            ("origin_x_nm",axis.origin[1]),("origin_y_nm",axis.origin[2]),
                            ("support_tmin",axis.tmin),("support_tmax",axis.tmax))
                        @test close(num(support,key),value)
                        @test all(a[key]==support[key] for a in audits if a["file"]==f)
                    end
                    RR.check_original_frame((axis=axis,),base[(f,1)])
                    original_pixels=Dict((int(r,"row"),int(r,"column"))=>r for r in original if r["in_fit"]=="true")
                    for arm in RR.ARMS
                        path=joinpath(prefix*".fit_data",f*".$arm.tsv")
                        isfile(path) || continue
                        data=last(read_table(path))
                        @test length(data)<=length(original_pixels)
                        @test all(haskey(original_pixels,(int(r,"row"),int(r,"column"))) for r in data)
                        @test all(r["x_nm"]==original_pixels[(int(r,"row"),int(r,"column"))]["x_nm"] &&
                            r["y_nm"]==original_pixels[(int(r,"row"),int(r,"column"))]["y_nm"] for r in data)
                    end
                end
            end
            @test isempty(intersect(seen,files)); union!(seen,files)
            for table in values(tables)
                @test Set(keys(table))==Set(k for k in expected if first(k) in files)
            end
            fitdata=Dict((f,a)=>last(read_table(joinpath(prefix*".fit_data",f*".$a.tsv"))) for f in files for a in RR.ARMS if isfile(joinpath(prefix*".fit_data",f*".$a.tsv")))
            for a in audits
                file,arm,profile=a["file"],a["arm"],a["profile"]
                n=count(k->first(k)==file,expected)
                @test int(a,"N")==n && a["valid"]=="true" && a["reason"]=="ok"
                @test int(a,"dx_px")==(arm=="control" ? 0 : shifts[file])
                @test (a["reused_zero_shift"]=="true")== (arm=="registered" && shifts[file]==0)
                @test a["optimizer_convergence"]=="unknown_native_not_exposed"
                model=deepcopy(raw["model"]); model["peak_profile"]=profile
                profile=="split" && (model["skew_ratio_max"]=assignment["model"]["split_skew_ratio_max"])
                _,cfg,_=RR.F.Extractor._configs(model,raw["preprocessing"],"unused")
                cfg.chain_circular_sigmas=a["source"]=="circ"
                pars=sort([r for r in parameters if (r["file"],r["arm"],r["profile"])==(file,arm,profile)];by=r->int(r,"parameter"))
                p=num.(pars,"value"); lo,hi=RR.VP.native_raw_bounds(n,cfg)
                @test int.(pars,"parameter")==collect(1:length(p))
                @test all(lo.<=p.<=hi) && all(isfinite,p)
                @test num.(pars,"lower")==lo && num.(pars,"upper")==hi
                pfull=3+n*(cfg.chain_circular_sigmas ? 4 : 5)+(profile=="split" ? n : 0)
                @test length(p)==int(a,"p_full")==pfull
                data=fitdata[(file,arm)]
                x=num.(data,"x_nm"); y=num.(data,"y_nm"); z=num.(data,"z_nm")
                @test length(z)==int(a,"n_data")>pfull
                @test all(isfinite,vcat(x,y,z,num.(data,"fwd_nm"),num.(data,"bwd_nm")))
                @test all(close(num(r,"z_nm"),(num(r,"fwd_nm")+num(r,"bwd_nm"))/2-num(a,"offset")) for r in data)
                @test length(Set((int(r,"row"),int(r,"column")) for r in data))==length(data)
                @test int(a,"observed_pixels")>=int(a,"roi_pixels")>=length(data)
                axis=(origin=(num(a,"origin_x_nm"),num(a,"origin_y_nm")),axis=[num(a,"axis_x"),num(a,"axis_y")],
                    perp=[-num(a,"axis_y"),num(a,"axis_x")],tmin=num(a,"support_tmin"),tmax=num(a,"support_tmax"))
                b,feats,ts,us,sp,sq=G._decode_chain(p,n,axis,cfg;amp_min=num(a,"amp_min"),amp_range=num(a,"amp_range"))
                @test close(sum(abs2,axis.axis),1.)
                @test all(cfg.sigma_parallel_min_nm.<=sp.<=cfg.sigma_parallel_max_nm)
                @test all(cfg.sigma_perp_min_nm.<=sq.<=cfg.sigma_perp_max_nm)
                @test all(abs.(us).<=cfg.lateral_max_nm)
                @test all(ts.>=axis.tmin-1e-6) && all(ts.<=axis.tmax+1e-6)
                @test all(diff(ts).>=G._effective_spacing_min_nm(cfg)-1e-12) && all(diff(ts).<=cfg.spacing_max_nm+1e-12)
                @test close(num(a,"amp_min"),cfg.min_amplitude_fraction*num(a,"amplitude_max_data"))
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
                @test close(rss,num(a,"rss"))
                @test close(length(z)/(length(z)-pfull)^2*rss,num(a,"gcv"))
                overlap=0.
                for i in 1:n-1,j in i+1:n
                    d2=(feats[i].x_nm-feats[j].x_nm)^2+(feats[i].y_nm-feats[j].y_nm)^2
                    overlap=max(overlap,exp(-0.5 * d2/max(mean(sp),mean(sq),G.EPS)^2))
                end
                @test close(overlap,num(a,"overlap")) && overlap<=cfg.max_overlap
                boots=[r for r in bootstrap if (r["file"],r["arm"],r["profile"])==(file,arm,profile)]
                @test Set(r["source"] for r in boots)==Set(["circ","ell"])
                @test num(a,"gcv")==minimum(num(r,"gcv") for r in boots if r["valid"]==r["success"]=="true")
                result=G.ChainModelResult(n=n,params=p,amp_min=num(a,"amp_min"),amp_range=num(a,"amp_range"),gcv=num(a,"gcv"))
                serialized=RR.R.feature_rows(file,result,(axisctx=axis,),cfg,a["source"])
                haskey(tables,(arm,profile)) && @test all(tables[(arm,profile)][(file,i)]==serialized[i] for i in 1:n)
                if arm=="registered" && shifts[file]==0
                    control=only(r for r in audits if (r["file"],r["arm"],r["profile"])==(file,"control",profile))
                    @test all(a[k]==control[k] for k in setdiff(keys(a),["arm","reused_zero_shift","elapsed_s"]))
                    @test int(a,"elapsed_s")==0
                    haskey(tables,(arm,profile)) && @test all(tables[(arm,profile)][(file,i)]==tables[("control",profile)][(file,i)] for i in 1:n)
                    @test fitdata[(file,arm)]==fitdata[(file,"control")]
                end
            end
        end
        @test seen==Set(first.(collect(expected)))
    end
end

function cohort(run,saved,previous)
    _,base=lobe_table(joinpath(run,"cached_features.tsv")); keys0=Set(keys(base))
    files=Set(first.(collect(keys0))); counts=Dict(f=>count(k->first(k)==f,keys0) for f in files)
    tables=("features","features_split","features_local","features_descriptor","features_predictor",
        "patches_fwd17","patches_bwd17","patches_bwd9","training_support","score_fwd","score_bwd","fisher_cv","pred_gmm","pred_kmeans","predictions")
    @testset "Complete cohorts, unchanged registration, exact reference replay, independent votes" begin
        @test !isfile(joinpath(run,"failures.tsv"))
        for suffix in ("",".summary.tsv",".peaks.tsv",".scores.tsv")
            @test read(joinpath(run,"shifts.tsv"*suffix))==read(joinpath(previous,"shifts.tsv"*suffix))
        end
        for arm in ("reference","control","registered")
            dir=joinpath(run,arm)
            @test !isfile(joinpath(dir,"failures.tsv"))
            for table in tables
                path=joinpath(dir,table*".tsv"); _,rows=lobe_table(path)
                @test Set(keys(rows))==keys0
                arm=="reference" && @test read(path)==read(joinpath(saved,table*".tsv"))
            end
            _,summary=read_table(joinpath(dir,"summary.tsv"))
            @test Dict(r["file"]=>int(r,"N_selected") for r in summary)==counts
            @test length(readdir(joinpath(dir,"plots","standalone")))==length(files)
            _,g=lobe_table(joinpath(dir,"pred_gmm.tsv")); _,k=lobe_table(joinpath(dir,"pred_kmeans.tsv")); _,p=lobe_table(joinpath(dir,"predictions.tsv"))
            for key in keys0
                if g[key]["predicted"]=="?" || k[key]["predicted"]=="?"
                    @test p[key]["predicted"]=="?" && p[key]["probability_1"]=="NA"
                else
                    value=(num(g[key],"probability_1")+num(k[key],"probability_1"))/2
                    @test isapprox(num(p[key],"probability_1"),value;atol=5.1e-9,rtol=0)
                    @test p[key]["predicted"]==(value>=.5 ? "1" : "0")
                end
            end
        end
        _,hashes=read_table(joinpath(run,"input_hashes.tsv"))
        settings_path=joinpath(run,"refit_settings.toml")
        frozen=isfile(settings_path) && RR.settings(settings_path).original_support
        paths=Dict("--settings"=>"config/acquisition_registration.toml",
            "--refit-settings"=>(frozen ? "config/registered_refit_original_support.toml" : "config/registered_refit.toml"),
            "--config"=>"config/unit_assignment_patch_support.toml","--count-config"=>"config/chitosan.toml",
            "--features"=>"results/fusion_comparison_20260920/run_v1/symmetric/features.tsv",
            "--split-features"=>"results/fusion_comparison_20260920/run_v1/symmetric/features_split.tsv",
            "--templates"=>"results/reconstructed_cc_soft_v1/full146_v1_inputs/templates_cc.tsv")
        @test Set(r["input"] for r in hashes)==Set(keys(paths))
        for r in hashes; @test r["sha256"]==bytes2hex(sha256(read(joinpath(ROOT,paths[r["input"]])))); end
        isfile(settings_path) && @test read(settings_path)==read(joinpath(ROOT,paths["--refit-settings"]))
    end
    geometry(run,joinpath(run,"cached_features.tsv"),joinpath(run,"shifts.tsv"))
end

"The original geometric support must also match the earlier complete native replay, not merely both new arms."
function original_reference(run, native)
    @test RR.settings(joinpath(run,"refit_settings.toml")).original_support
    prior=Dict{String,Dict{String,String}}(); prior_pixels=Dict{String,String}()
    for path in readdir(native;join=true)
        endswith(path,".fits.tsv") || continue
        prefix=chop(path;tail=length(".fits.tsv"))
        for r in last(read_table(path))
            r["method"]=="native" || continue
            haskey(prior,r["file"]) && error("Repeated native reference")
            prior[r["file"]]=r
            prior_pixels[r["file"]]=joinpath(prefix*".fit_data",r["file"]*".tsv")
        end
    end
    seen=Set{String}()
    @testset "Frozen support exactly replays the earlier native geometric context" begin
        for path in readdir(run;join=true)
            endswith(path,".support.tsv") || continue
            prefix=chop(path;tail=length(".support.tsv"))
            for r in last(read_table(path))
                f=r["file"]; @test !(f in seen); push!(seen,f)
                @test haskey(prior,f)
                for key in ("axis_x","axis_y","origin_x_nm","origin_y_nm","support_tmin","support_tmax")
                    @test r[key]==prior[f][key]
                end
                @test int(r,"fit_pixels")==int(prior[f],"n_data")
                original=filter(r->r["in_fit"]=="true",last(read_table(joinpath(prefix*".fit_data",f*".original.tsv"))))
                old=last(read_table(prior_pixels[f]))
                @test length(original)==length(old)
                @test all(a["x_nm"]==b["x_nm"] && a["y_nm"]==b["y_nm"] && a["z_native_nm"]==b["z_nm"] for (a,b) in zip(original,old))
            end
        end
        @test seen==Set(keys(prior))
    end
end

"Verify recorded failures and every completed fit without turning them into a benchmark."
function recorded_failures(run,saved,previous)
    @testset "Indeterminate comparison: explicit failures, no partial classifier or grade" begin
        _,top=read_table(joinpath(run,"failures.tsv"))
        @test length(top)==1 && top[1]["stage"]=="refit"
        @test !isdir(joinpath(run,"control")) && !isdir(joinpath(run,"registered"))
        failures=reduce(vcat,[last(read_table(p)) for p in readdir(run;join=true) if endswith(p,".failures.tsv")];init=Dict{String,String}[])
        @test !isempty(failures)
        @test length(Set(r["file"] for r in failures))==length(failures)
        _,base=lobe_table(joinpath(run,"cached_features.tsv"))
        @test all(any(first(k)==r["file"] for k in keys(base)) for r in failures)
        for suffix in ("",".summary.tsv",".peaks.tsv",".scores.tsv")
            @test read(joinpath(run,"shifts.tsv"*suffix))==read(joinpath(previous,"shifts.tsv"*suffix))
        end
        for table in ("features","features_split","features_local","features_descriptor","features_predictor",
                "patches_fwd17","patches_bwd17","patches_bwd9","training_support","score_fwd","score_bwd","fisher_cv","pred_gmm","pred_kmeans","predictions")
            @test read(joinpath(run,"reference",table*".tsv"))==read(joinpath(saved,table*".tsv"))
        end
    end
    geometry(run,joinpath(run,"cached_features.tsv"),joinpath(run,"shifts.tsv");complete=false)
    println("Recorded failures checked. The recognition comparison is INDETERMINATE; no candidate grade is produced.")
end

function main(args=ARGS)
    if length(args)==5 && first(args)=="--original-support"
        cohort(args[2:4]...)
        return original_reference(args[2],args[5])
    end
    if length(args)==4 && first(args)=="--recorded-failures"
        return recorded_failures(args[2:end]...)
    end
    length(args) in (1,3) || error("Usage: verify_registered_refit.jl SMOKE_DIR OR [--recorded-failures] RUN_DIR SAVED_SUPPORT_DIR PREVIOUS_REGISTRATION_RUN OR --original-support RUN_DIR SAVED_SUPPORT_DIR PREVIOUS_REGISTRATION_RUN NATIVE_GEOMETRY_RUN")
    length(args)==1 ? geometry(args[1],joinpath(args[1],"features.tsv"),joinpath(args[1],"shifts.tsv")) : cohort(args...)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && RegisteredRefitVerification.main()
