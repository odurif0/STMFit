# Synthetic data only. Also runs the original all-scan exclusion regressions.
include(joinpath(@__DIR__,"test_scan_exclusion.jl"))
const COMPLETE_SETTINGS=joinpath(F.ROOT,"config","scan_exclusion_complete_observed.toml")

@testset "Observation-defined learning, retained partial targets, exact arm invariance" begin
    setting=S.settings(COMPLETE_SETTINGS)
    @test setting["selection"]["training_cohort"]=="fully_observed_fwd_bwd"
    @test !haskey(S.settings(EXCLUSION_SETTINGS)["selection"],"training_cohort")
    mktempdir() do dir
        function replace_fixture_table(path,header,rows)
            # Deliberate mutations only inside this disposable fixture; retain
            # the production writer's refusal to overwrite existing outputs.
            startswith(abspath(path),abspath(dir)*"/") || error("Not a fixture path")
            mktempdir() do staging
                changed=joinpath(staging,"changed.tsv")
                F.RU.write_table(changed,header,rows); cp(changed,path;force=true)
            end
        end
        root=whole_fixture(joinpath(dir,"input");fully_observed=3)
        data=F.inputs(root); files=sort(collect(keys(data.counts)))
        eligible=S.eligible_files(setting,data)
        @test eligible==Set(files[1:3])
        @test S.eligible_files(S.settings(EXCLUSION_SETTINGS),data)==Set(files)
        @test_throws ErrorException S.eligible_files(setting,merge(data,(groups=Dict(f=>"partial" for f in files),)))
        bad=deepcopy(setting); bad["selection"]["training_cohort"]="benchmark_filtered"
        badpath=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),badpath,"w")
        @test_throws ErrorException S.settings(badpath)
        audit=joinpath(root,"refit_chunk1.tsv.fits.tsv")
        header,rows=F.RU.read_table(audit)
        for value in ("101","-1","99")
            altered=deepcopy(rows); altered[1]["observed_pixels"]=value
            replace_fixture_table(audit,header,altered)
            @test_throws ErrorException F.inputs(root)
        end
        altered=deepcopy(rows)
        for r in altered; r["observed_pixels"]="0"; r["total_pixels"]="0"; end
        replace_fixture_table(audit,header,altered)
        @test_throws ErrorException F.inputs(root)
        replace_fixture_table(audit,header,rows)
        native=joinpath(dir,"fixture_native"); D.execute(root,SETTINGS,native)
        for arm in ("control","observed"), stage in D.REFERENCES
            cp(joinpath(native,arm,stage*".tsv"),joinpath(root,arm,stage*".tsv");force=true)
        end
        launcher=joinpath(F.ROOT,"hpc","diagnose_complete_observation.sbatch")
        env=["STMFIT_PROJECT_DIR"=>F.ROOT,"STMFIT_INPUT_DIR"=>root,"STMFIT_OUTDIR"=>joinpath(dir,"hpc"),
            "JULIA_BIN"=>first(Base.julia_cmd().exec)]
        @test success(addenv(`bash $launcher --dry-run`,env...))
        @test !ispath(joinpath(dir,"hpc"))
        @test !success(pipeline(ignorestatus(addenv(`bash $launcher`,env...,"SLURM_JOB_ID"=>""));stdout=devnull,stderr=devnull))
        @test occursin("#SBATCH --no-requeue",read(launcher,String))
        @test occursin("#SBATCH --time=01:00:00",read(launcher,String))
        for arm in ("control","observed")
            out=joinpath(dir,arm); S.execute(root,COMPLETE_SETTINGS,out,arm)
            @test V.verify(root,out)===nothing
            @test length(F.Pipeline.check_counts(joinpath(out,"predictions.tsv"),data.counts))==48
            _,folds=F.RU.read_table(joinpath(out,"folds.tsv"))
            @test [parse(Int,r["training_scans"]) for r in folds]==[2,2,2,3]
            @test [parse(Int,r["training_lobes"]) for r in folds]==[24,24,24,36]
        end
        @test V.verify_complete_pair(root,joinpath(dir,"control"),joinpath(dir,"observed"))===nothing
        data=F.inputs(root); opts=F.EF.load_fisher_config(data.cfg)
        patches=F.EF.load_patches(joinpath(root,"control","patches_fwd17.tsv"),"res",opts)
        @test Set(first.(S.training_patches(patches,files[1];training_files=eligible).keys))==Set(files[2:3])
        @test_throws ErrorException S.training_patches(patches,files[1];training_files=Set(["absent.sxm"]))
        @test_throws ErrorException S.training_patches(patches,files[1];training_files=Set{String}())
        @test_throws ErrorException S.fit_fold(root,"control",joinpath(dir,"bad_native"),"",patches,opts,data;training_files=eligible)
        @test !ispath(joinpath(dir,"bad_native"))
        # A partial scan may be evaluated but cannot teach any shared statistic.
        # A complete held-out scan must not teach either. Mutate these separately.
        x=copy(patches.X); amps=copy(patches.amplitudes)
        descriptor=joinpath(root,"control","features_descriptor.tsv")
        h,rs=F.RU.read_table(descriptor)
        for (step,changed_file) in enumerate((files[4],files[1]))
            for i in eachindex(amps)
                patches.keys[i][1]==changed_file || continue
                x[i,:].+=3.; amps[i]+=17.
            end
            for r in rs
                r["file"]==changed_file || continue
                for f in vcat(["amplitude"],F.BASE,["patch_u_asym_reconstructed"])
                    r[f]=string(7-2parse(Float64,r[f]))
                end
            end
            step==2 && (first(filter(r->r["file"]==files[1],rs))["amp_rel"]="NaN")
            replace_fixture_table(descriptor,h,rs)
            mutated=F.EF.PatchTable(patches.keys,x,amps,patches.invalid_reasons,patches.grid)
            out=joinpath(dir,"mutation$step")
            S.fit_fold(root,"control",out,files[1],mutated,opts,data;training_files=eligible)
            original=joinpath(dir,"control","fold0001")
            @test read(joinpath(out,"models.toml"))==read(joinpath(original,"models.toml"))
            if step==1
                for stage in ("target_features","target_scales",S.STAGES...)
                    @test read(joinpath(out,stage*".tsv"))==read(joinpath(original,stage*".tsv"))
                end
            else
                @test read(joinpath(out,"target_scales.tsv"))!=read(joinpath(original,"target_scales.tsv"))
                predictions=last(F.RU.lobe_table(joinpath(out,"predictions.tsv")))
                @test length(predictions)==data.counts[files[1]]
                @test predictions[(files[1],1)]["predicted"]=="?"
            end
        end
    end
end
