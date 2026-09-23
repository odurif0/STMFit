using Test, TOML
include(joinpath(@__DIR__,"run_promoted_count_assignment.jl"))
const P=PromotedCountAssignment
function fixture(dir)
    raw=joinpath(dir,"raw"); mkdir(raw)
    for f in ("a.sxm","b.sxm"); write(joinpath(raw,f),"metadata fixture"); end
    c=joinpath(dir,"control.tsv"); h=joinpath(dir,"hybrid.tsv")
    P.write_table(c,["filepath","N_selected"],[Dict("filepath"=>f,"N_selected"=>string(n)) for (f,n) in (("a.sxm",3),("b.sxm",5))])
    P.write_table(h,["filepath","N_selected","selection_policy"],[Dict("filepath"=>f,"N_selected"=>string(n),
        "selection_policy"=>"support_midpoint_hybrid") for (f,n) in (("a.sxm",4),("b.sxm",5))])
    t=joinpath(dir,"templates.tsv"); write(t,"fixture")
    Dict("--data-dir"=>raw,"--count-config"=>joinpath(P.ROOT,"config/chitosan.toml"),
        "--config"=>joinpath(P.ROOT,"config/unit_assignment_patch_support.toml"),
        "--templates"=>t,"--control-summary"=>c,"--promoted-summary"=>h,"--outdir"=>joinpath(dir,"out"))
end
function fake_runner(calls)
    function runner(out,stage,script,args)
        push!(calls,(stage,script,copy(args)))
        value(k)=args[findfirst(==(k),args)+1]
        counts=P.selected_counts(value("--selected-summary")); dir=value("--outdir"); mkdir(dir)
        rows=[Dict("file"=>f,"lobe"=>string(i),"N"=>string(n)) for (f,n) in sort(collect(counts)) for i in 1:n]
        for name in ("features.tsv","features_split.tsv","predictions.tsv")
            P.write_table(joinpath(dir,name),["file","lobe","N"],rows)
        end
    end
end
@testset "Count-policy handoff is label-free and rebuilds both arms" begin
    @test P.options(["--help"])===nothing
    for k in ("--truth","--expected-N","--manifest","--features","--selection-policy","--settings")
        @test_throws ErrorException P.options([k,"forbidden"])
    end
    mktempdir() do dir
        o=fixture(dir); args=vcat([[k,v] for (k,v) in o]...)
        @test P.options(args)==o
        @test_throws ErrorException P.options(vcat(args,["--config",o["--config"]]))
        P.execute(merge(o,Dict("--dry-run"=>"true")))
        @test !ispath(o["--outdir"])
        calls=[]; P.execute(o;runner=fake_runner(calls))
        @test first.(calls)==["control","hybrid"]
        for (stage,script,a) in calls
            @test script=="run_reconstructed_chitosan.jl"
            @test !any(k in a for k in ("--features","--split-features","--mold-tangent-settings","--manifest"))
            @test a[findfirst(==("--config"),a)+1]==o["--config"]
            @test a[findfirst(==("--templates"),a)+1]==o["--templates"]
        end
        @test P.selected_counts(joinpath(o["--outdir"],"control_counts.tsv"))==Dict("a.sxm"=>3,"b.sxm"=>5)
        @test P.selected_counts(joinpath(o["--outdir"],"hybrid_counts.tsv"))==Dict("a.sxm"=>4,"b.sxm"=>5)
        _,changes=P.read_table(joinpath(o["--outdir"],"count_changes.tsv"))
        @test [r["delta"] for r in changes]==["1","0"]
        @test read(joinpath(o["--outdir"],"original_promoted_summary.tsv"))==read(o["--promoted-summary"])
        @test !isfile(joinpath(o["--outdir"],"failures.tsv"))
        @test_throws ErrorException P.options(args)
        failed=merge(o,Dict("--outdir"=>joinpath(dir,"failed")))
        @test_throws ErrorException P.execute(failed;runner=(args...)->error("fixture failure"))
        _,fs=P.read_table(joinpath(failed["--outdir"],"failures.tsv"))
        @test only(fs)["stage"]=="control"
        @test !ispath(joinpath(failed["--outdir"],"hybrid"))
    end
end
@testset "Reject mixed policy, unavailable counts, labels and cohort mismatch" begin
    for (key,val) in (("selection_policy","gcv"),("refined_policy","adaptive_support_rescue"),
                      ("status","failed"),("N_selected","0"),("expected_N","6"),("filepath","absent.sxm"))
        mktempdir() do dir
            o=fixture(dir); header,rs=P.read_table(o["--promoted-summary"])
            key in header || push!(header,key)
            rs[1][key]=val
            path=joinpath(dir,"bad.tsv"); P.write_table(path,header,rs)
            @test_throws ErrorException P.inputs(merge(o,Dict("--promoted-summary"=>path)))
            @test !ispath(o["--outdir"])
        end
    end
end
@testset "Slurm limits and required exports" begin
    path=joinpath(P.ROOT,"hpc/compare_promoted_counts.sbatch")
    @test success(`bash -n $path`)
    source=read(path,String)
    for word in ("--time=04:00:00","--cpus-per-task=4","--mem=16000MB","SLURM_JOB_ID", "STMFIT_PROMOTED_SUMMARY", "STMFIT_CONTROL_SUMMARY")
        @test occursin(word,source)
    end
end
