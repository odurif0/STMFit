using Test, TOML
include(joinpath(@__DIR__, "run_hybrid_reproduction.jl"))
const H = HybridReproduction

function reproduction_fixture(dir; n=5)
    raw = joinpath(dir, "raw"); mkdir(raw)
    for i in 1:n; write(joinpath(raw, "scan$i.sxm"), "Metadata-only synthetic fixture"); end
    templates = joinpath(dir, "templates.tsv"); write(templates, "fixture")
    return Dict("--data-dir"=>raw, "--count-config"=>joinpath(H.ROOT, "config/chitosan.toml"),
        "--config"=>joinpath(H.ROOT, "config/unit_assignment_patch_support.toml"),
        "--templates"=>templates, "--outdir"=>joinpath(dir, "output"))
end
const REPRO_HEADER = ["filepath", "status", "N_selected", "selection_policy"]
countrow(file, n) = Dict("filepath"=>file, "status"=>"ok", "N_selected"=>string(n),
    "selection_policy"=>"support_midpoint_hybrid")

@testset "Fresh reproduction boundaries" begin
    @test H.options(["--help"]) === nothing
    for k in ("--truth", "--expected-N", "--manifest", "--selected-summary", "--features",
              "--split-features", "--selection-policy", "--promoted-summary", "--reference")
        @test_throws ErrorException H.options([k, "forbidden"])
    end
    mktempdir() do dir
        o = reproduction_fixture(dir)
        args = vcat([[k,v] for (k,v) in o]...)
        @test H.options(args) == o
        @test_throws ErrorException H.options(vcat(args, ["--config", o["--config"]]))
        H.execute(merge(o, Dict("--dry-run"=>"true")))
        @test !ispath(o["--outdir"])
        @test_throws ErrorException H.execute(o; thread_budget=0)
        for (section,key,val) in (("model","selection_policy","gcv"), ("model","selection_criterion","bic"))
            cfg = TOML.parsefile(o["--count-config"]); cfg[section][key] = val
            path = joinpath(dir, key * ".toml"); open(io->TOML.print(io,cfg),path,"w")
            @test_throws ErrorException H.execute(merge(o,Dict("--count-config"=>path,"--dry-run"=>"true")))
        end
        @test !ispath(o["--outdir"])
    end
end

@testset "Complete count merge rejects labels, failures and partial cohorts" begin
    mktempdir() do dir
        a = joinpath(dir,"a.tsv"); b = joinpath(dir,"b.tsv")
        H.write_table(a, REPRO_HEADER, [countrow("b.sxm",3)])
        H.write_table(b, REPRO_HEADER, [countrow("a.sxm",5)])
        output = joinpath(dir,"merged.tsv")
        @test H.merge_counts([a,b], ["a.sxm","b.sxm"], output) == Dict("a.sxm"=>5,"b.sxm"=>3)
        _, rows = H.read_table(output)
        @test [r["filepath"] for r in rows] == ["a.sxm","b.sxm"]
        @test_throws ErrorException H.merge_counts([a,b], ["a.sxm","b.sxm"], output)
        @test_throws ErrorException H.merge_counts([a,a], ["b.sxm"], joinpath(dir,"duplicate.tsv"))
        @test_throws ErrorException H.merge_counts([a], ["a.sxm","b.sxm"], joinpath(dir,"missing.tsv"))
        @test_throws ErrorException H.merge_counts([a,b], ["a.sxm"], joinpath(dir,"extra.tsv"))
        for (i,(key,val)) in enumerate((("status","error"),("N_selected","0"),
                ("selection_policy","gcv"),("expected_N","6"),("refined_policy","adaptive_support_rescue")))
            header = copy(REPRO_HEADER); key in header || push!(header,key)
            row = countrow("a.sxm",4); row[key] = val
            path = joinpath(dir,"bad$i.tsv"); H.write_table(path,header,[row])
            @test_throws ErrorException H.merge_counts([path],["a.sxm"],joinpath(dir,"out$i.tsv"))
            @test !ispath(joinpath(dir,"out$i.tsv"))
        end
    end
end

@testset "Two raw executions, bounded shards, no saved-count substitution" begin
    for n in (1,5)
        mktempdir() do dir
            o = reproduction_fixture(dir;n)
            calls = Channel{Any}(20)
            function runner(out,stage,script,args; threads)
                put!(calls,(stage,script,copy(args),threads))
                value(k) = args[findfirst(==(k),args)+1]
                dest = value("--outdir"); mkdir(dest)
                if script == "batch_full.jl"
                    @test threads == 1
                    @test !ispath(value("--tsv"))
                    @test "--skip-1d" in args
                    @test parse(Int,first(args)) == n
                    i,total = parse.(Int,split(value("--chunk"),'/'))
                    @test total == min(n,4)
                    # Deliberately distinct counts in repeat2: the runner must
                    # retain both, not demand or manufacture identical results.
                    count = basename(out) == "repeat1" ? 3 : 4
                    rows = [countrow("scan$j.sxm",count) for j in 1:n if mod1(j,total)==i]
                    name = total == 1 ? "summary_overlap060_hard.tsv" :
                        "summary_overlap060_hard_chunk$(lpad(i,2,'0'))of$(lpad(total,2,'0')).tsv"
                    H.write_table(joinpath(dest,name),REPRO_HEADER,rows)
                else
                    @test script == "run_reconstructed_chitosan.jl"
                    @test threads == 4
                    @test !any(k in args for k in ("--features","--split-features","--manifest","--truth"))
                    @test startswith(value("--selected-summary"),out)
                    counts = H.selected_counts(value("--selected-summary"))
                    rows = [Dict("file"=>f,"lobe"=>string(i)) for (f,c) in counts for i in 1:c]
                    for table in ("features.tsv","features_split.tsv","predictions.tsv")
                        H.write_table(joinpath(dest,table),["file","lobe"],rows)
                    end
                end
            end
            H.execute(o;runner,thread_budget=4)
            close(calls); events = collect(calls)
            @test length(events) == 2(min(n,4)+1)
            @test [e[2] for e in events][[min(n,4)+1,end]] == fill("run_reconstructed_chitosan.jl",2)
            for i in 1:2
                counts = H.selected_counts(joinpath(o["--outdir"],"repeat$i","counting_summary.tsv"))
                @test length(counts)==n && all(==(i+2),values(counts))
            end
            @test !isfile(joinpath(o["--outdir"],"failures.tsv"))
            @test_throws ErrorException H.execute(o;runner,thread_budget=4)
        end
    end
    mktempdir() do dir
        o = reproduction_fixture(dir)
        runner(args...;kwargs...) = error("Synthetic counting failure")
        @test_throws Exception H.execute(o;runner,thread_budget=4)
        _, rows = H.read_table(joinpath(o["--outdir"],"failures.tsv"))
        @test only(rows)["stage"] == "repeat1_counting"
        @test !ispath(joinpath(o["--outdir"],"repeat2"))
    end
end

@testset "Simple bounded Slurm execution" begin
    path = joinpath(H.ROOT,"hpc/reproduce_hybrid_champion.sbatch")
    @test success(`bash -n $path`)
    source = read(path,String)
    for value in ("--time=04:00:00","--cpus-per-task=4","--mem=16000MB","SLURM_JOB_ID","--dry-run")
        @test occursin(value,source)
    end
    @test !occursin("--selected-summary",source)
end
