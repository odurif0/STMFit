using Test, Random, LinearAlgebra, Printf, SHA
include(joinpath(@__DIR__,"diagnose_fisher_attribution.jl"))
const FA = FisherAttributionDiagnostic
const EF = FA.EF
const RU = FA.RU
const RD = FA.RD
const CONFIG = joinpath(@__DIR__,"..","config","unit_assignment_reconstructed.toml")

"Synthetic forward patches only. No saved cohort table is opened in this suite."
function synthetic_patches(path; n=36, invalid=true, constant=false)
    options=EF.load_fisher_config(CONFIG)
    grid=EF.fisher_grid(options)
    coords=collect(-8:8).*options.step_nm
    t,u=RD.grid_vectors(coords)
    columns=[@sprintf("res_p%03d",i) for i in 1:grid.side^2]
    rng=MersenneTwister(29)
    rows=Dict{String,String}[]
    for i in 1:n
        # Both amplitude groups are represented in both parity folds.
        group=mod(i,4)<2 ? -1.0 : 1.0
        p=constant ? ones(grid.side^2) :
            group .* exp.(-30 .* (t.^2+u.^2)) .+ 0.3 .* sin.(11t .+ i/7) .+
            0.2 .* cos.(13u .+ i/5) .+ 0.06randn(rng,grid.side^2)
        invalid && i==n && (p[1]=NaN)
        row=Dict("file"=>"synthetic.sxm","lobe"=>string(i))
        for (column,value) in zip(columns,p)
            row[column]=isfinite(value) ? @sprintf("%.7g",value) : "NA"
        end
        push!(rows,row)
    end
    RU.write_table(path,vcat(["file","lobe"],columns),rows)
    return EF.load_patches(path,"res",options),options
end

@testset "Retained native fold loop exactly agrees with cv_scores" begin
    mktempdir() do dir
        for (name,n,invalid,constant) in (("normal",36,true,false),("small",8,false,false),("constant",36,false,true))
            patches,options=synthetic_patches(joinpath(dir,name*".tsv");n,invalid,constant)
            projected=EF.load_fisher_config(joinpath(@__DIR__,"..","config","unit_assignment_affine_fisher.toml"))
            @test_throws ArgumentError FA.replay_folds(patches,projected)
            native=EF.cv_scores(patches,options)
            replay=FA.replay_folds(patches,options)
            @test length(replay.scores)==n
            @test [r.invalid_reason for r in replay.scores]==[r.invalid_reason for r in native]
            @test isequal([r.score for r in replay.scores],[r.score for r in native])
            @test [(r.file,r.lobe) for r in replay.scores]==patches.keys
            nativepath,replaypath=joinpath(dir,name*"_native.tsv"),joinpath(dir,name*"_replay.tsv")
            EF.write_scores(nativepath,native);EF.write_scores(replaypath,replay.scores)
            @test read(nativepath)==read(replaypath)
            comparison=FA.compare_saved(nativepath,replaypath)
            @test comparison.all_match
            @test comparison.matching==n
            if name=="normal"
                @test Set(keys(replay.models))==Set([0,1])
                @test count(r->isfinite(r.score),replay.scores)==n-1
                @test replay.scores[end].invalid_reason=="invalid_patch_value:res_p001"
                for info in replay.folds
                    parity=info["training_parity"]=="even" ? 0 : 1
                    @test info["training_rows"]==count(i->isempty(patches.invalid_reasons[i]) && mod(patches.keys[i][2],2)==parity,1:n)
                    @test info["shared_centering_offset"]==-dot(replay.models[parity].mid,replay.models[parity].w_p)
                end
                rows=FA.attribution_rows(patches,replay,comparison)
                @test length(rows)==n
                for (i,row) in enumerate(rows[1:end-1])
                    @test row["training_parity"]==(iseven(i) ? "odd" : "even")
                    @test row["maximum_response"]==native[i].score
                    @test row["even_centered"]+row["absolute_odd_bonus"] ≈ row["maximum_response"] atol=1e-12
                    @test abs(row["analytic_identity_error"])<=1e-12
                end
                @test isnan(rows[end]["maximum_response"])
                @test rows[end]["native_invalid_reason"]==native[end].invalid_reason
            else
                @test isempty(replay.models)
                @test all(r->!isempty(r.invalid_reason),replay.scores)
            end
        end
    end
end

@testset "Comparison exposes score/reason/key differences without tuning" begin
    mktempdir() do dir
        header=["file","lobe","score","invalid_reason"]
        saved=[Dict("file"=>"s.sxm","lobe"=>string(i),"score"=>"1.000000","invalid_reason"=>"") for i in 1:4]
        changed=deepcopy(saved)
        changed[1]["score"]="1.000001"
        changed[2]["invalid_reason"]="synthetic_reason"
        changed[3]["score"]="1.0" # numeric equality alone is insufficient
        changed[4]["lobe"]="5"
        s,r=joinpath(dir,"saved.tsv"),joinpath(dir,"replay.tsv")
        RU.write_table(s,header,saved);RU.write_table(r,header,changed)
        result=FA.compare_saved(s,r)
        @test !result.all_match
        @test result.matching==0
        @test length(result.rows)==5
        @test result.rows[1]["status"]=="score_text_difference"
        @test result.rows[2]["status"]=="reason_difference"
        @test result.rows[3]["status"]=="score_text_difference"
        @test result.rows[4]["status"]=="missing_replay_key"
        @test result.rows[5]["status"]=="missing_saved_key"
        @test result.rows[3]["replayed_minus_saved_serialized"]==0
    end
end

@testset "Synthetic CLI outputs, explicit realization status, input/output boundaries" begin
    mktempdir() do dir
        patches_path=joinpath(dir,"synthetic_patches.tsv")
        patches,options=synthetic_patches(patches_path)
        saved=joinpath(dir,"saved.tsv")
        EF.write_scores(saved,EF.cv_scores(patches,options))
        hashes=Dict(path=>bytes2hex(sha256(read(path))) for path in (patches_path,saved,CONFIG))
        out=joinpath(dir,"replay")
        comparison=FA.run_replay(patches_path,"res",CONFIG,saved,out)
        @test comparison.all_match
        @test all(bytes2hex(sha256(read(path)))==hash for (path,hash) in hashes)
        @test read(saved)==read(joinpath(out,"replayed_fisher_cv.tsv"))
        _,weights=RU.read_table(joinpath(out,"new_fold_weights_mid.tsv"))
        @test length(weights)==2length(patches.grid.disk_indices)
        @test all(haskey(r,"new_centered_mid") for r in weights)
        _,rows=RU.lobe_table(joinpath(out,"replayed_response_attribution.tsv"))
        @test Set(keys(rows))==Set(patches.keys)
        @test all(r["attribution_scope"]=="new_realization_matches_all_saved_exports_not_original_weights" for r in values(rows))
        @test rows[("synthetic.sxm",36)]["maximum_response"]=="NA"
        @test occursin("NOT prove",read(joinpath(out,"report.md"),String))
        @test_throws ErrorException FA.run_replay(patches_path,"res",CONFIG,saved,out)
        empty=joinpath(dir,"empty");mkdir(empty)
        @test_throws ErrorException FA.run_replay(patches_path,"res",CONFIG,saved,empty)
        dangling=joinpath(dir,"dangling");symlink(joinpath(dir,"missing"),dangling)
        @test_throws ErrorException FA.run_replay(patches_path,"res",CONFIG,saved,dangling)
        @test_throws ErrorException FA.run_replay(patches_path,"bwd_res",CONFIG,saved,joinpath(dir,"badprefix"))
        # Fixed replay is NOT silently promoted when saved scores differ.
        changed=readlines(saved)
        fields=split(changed[2],'\t';keepempty=true);fields[3]="123.000000";changed[2]=join(fields,'\t')
        mismatch_saved=joinpath(dir,"mismatch_saved.tsv")
        write(mismatch_saved,join(changed,'\n')*"\n")
        mismatch_out=joinpath(dir,"mismatch")
        result=FA.run_replay(patches_path,"res",CONFIG,mismatch_saved,mismatch_out)
        @test !result.all_match
        _,mismatch=RU.lobe_table(joinpath(mismatch_out,"replayed_response_attribution.tsv"))
        @test all(r["attribution_scope"]=="new_realization_only_saved_exports_differ" for r in values(mismatch))
        @test occursin("Do not attribute",read(joinpath(mismatch_out,"report.md"),String))
        # Real CLI subprocess remains synthetic and uses the same root Julia1.13.
        cliout=joinpath(dir,"cli")
        command=`$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(@__DIR__)) $(joinpath(@__DIR__,"diagnose_fisher_attribution.jl")) --patches $patches_path --prefix res --config $CONFIG --saved-scores $saved --outdir $cliout`
        @test occursin("all_match=true",read(command,String))
        @test_throws ErrorException FA.main(["--patches","x"])
        @test_throws ErrorException FA.main(["--truth","x"])
        @test_throws ErrorException FA.main(["--patches","x","--patches","y"])
    end
end
