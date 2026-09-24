# Run the frozen-learning regressions and reuse their synthetic fixtures.
include(joinpath(@__DIR__,"test_frozen_learning.jl"))
include(joinpath(@__DIR__,"verify_scan_exclusion.jl"))
const V=ScanExclusionVerification
const S=V.S
const F=S.D
const EXCLUSION_SETTINGS=joinpath(F.ROOT,"config","scan_exclusion.toml")

@testset "Whole-scan exclusion: completeness, no target learning, frozen defaults" begin
    @test S.settings(EXCLUSION_SETTINGS)["selection"]["arms"]==["control","observed"]
    @test S.parse_cli(["--help"])===nothing
    mktempdir() do dir
        root=whole_fixture(joinpath(dir,"input"))
        # Replace the fixture's dummy references by exact native synthetic outputs.
        native=joinpath(dir,"fixture_native"); D.execute(root,SETTINGS,native)
        for arm in ("control","observed"), stage in D.REFERENCES
            cp(joinpath(native,arm,stage*".tsv"),joinpath(root,arm,stage*".tsv");force=true)
        end
        out=joinpath(dir,"control"); args=["--input",root,"--settings",EXCLUSION_SETTINGS,"--outdir",out,"--arm","control"]
        @test S.main(vcat(args,["--dry-run"]))===nothing
        @test !ispath(out)
        for forbidden in ("--truth","--expected-n","--seed","--threshold","--manifest","--training-support")
            @test_throws ErrorException S.parse_cli(vcat(args,[forbidden,"forbidden"]))
        end
        @test_throws ErrorException S.parse_cli(vcat(args,["--arm","observed"]))
        launcher=joinpath(F.ROOT,"hpc","diagnose_scan_exclusion.sbatch")
        env=["STMFIT_PROJECT_DIR"=>F.ROOT,"STMFIT_INPUT_DIR"=>root,"STMFIT_OUTDIR"=>joinpath(dir,"hpc"),
            "JULIA_BIN"=>first(Base.julia_cmd().exec)]
        @test success(addenv(`bash $launcher --dry-run`,env...))
        @test !ispath(joinpath(dir,"hpc"))
        @test !success(pipeline(ignorestatus(addenv(`bash $launcher`,env...,"SLURM_JOB_ID"=>""));stdout=devnull,stderr=devnull))
        @test occursin("#SBATCH --no-requeue",read(launcher,String))
        @test occursin("#SBATCH --time=01:00:00",read(launcher,String))
        for arm in ("control","observed")
            armout=joinpath(dir,arm)
            S.execute(root,EXCLUSION_SETTINGS,armout,arm)
            @test V.verify(root,armout)===nothing
            @test length(F.Pipeline.check_counts(joinpath(armout,"predictions.tsv"),F.inputs(root).counts))==48
        end
        @test_throws ErrorException S.parse_cli(args)
        @test_throws ErrorException S.execute(root,EXCLUSION_SETTINGS,out,"control")
        data=F.inputs(root); opts=F.EF.load_fisher_config(data.cfg)
        patches=F.EF.load_patches(joinpath(root,"control","patches_fwd17.tsv"),"res",opts)
        held=first(sort(collect(keys(data.counts))))
        @test_throws ErrorException S.training_patches(patches,"absent.sxm")
        # Change ALL held-out pixels, center amplitudes, local features and naming
        # amplitudes. Learned parameters and every training statistic must stay exact.
        x=copy(patches.X); amps=copy(patches.amplitudes)
        for i in eachindex(amps)
            patches.keys[i][1]==held || continue
            x[i,:].+=3.; amps[i]+=17.
        end
        altered=F.EF.PatchTable(patches.keys,x,amps,patches.invalid_reasons,patches.grid)
        descriptor=joinpath(root,"control","features_descriptor.tsv")
        header,rows=F.RU.read_table(descriptor)
        for r in rows
            r["file"]==held || continue
            for f in vcat(["amplitude"],F.BASE,["patch_u_asym_reconstructed"])
                r[f]=string(7-2parse(Float64,r[f]))
            end
        end
        first(filter(r->r["file"]==held,rows))["amp_rel"]="NaN"
        changed=joinpath(dir,"changed.tsv"); F.RU.write_table(changed,header,rows)
        cp(changed,descriptor;force=true)
        mutation=joinpath(dir,"mutation")
        S.fit_fold(root,"control",mutation,held,altered,opts,data)
        original=joinpath(out,"fold0001")
        @test read(joinpath(original,"models.toml"))==read(joinpath(mutation,"models.toml"))
        @test read(joinpath(original,"target_scales.tsv"))!=read(joinpath(mutation,"target_scales.tsv"))
        @test read(joinpath(original,"target_features.tsv"))!=read(joinpath(mutation,"target_features.tsv"))
        predictions=last(F.RU.lobe_table(joinpath(mutation,"predictions.tsv")))
        @test length(predictions)==data.counts[held]
        @test predictions[(held,1)]["predicted"]=="?"
    end
end
