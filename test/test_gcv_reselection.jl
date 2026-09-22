using Test
include(joinpath(@__DIR__, "run_gcv_reselection.jl"))
module NativeGCVExtractor
include(joinpath(@__DIR__, "extract_lobe_features.jl"))
end

const RESELECT_JOB = joinpath(ROOT, "hpc", "compare_gcv_reselection.sbatch")

function count_fixture(dir)
    raw = joinpath(dir, "raw"); mkdir(raw)
    for f in ("alpha.sxm", "beta.sxm")
        write(joinpath(raw, f), "Not an SXM: metadata-only synthetic test\n")
    end
    summary = joinpath(dir, "old.tsv")
    write_table(summary, ["filepath", "N_selected"],
        [Dict("filepath"=>joinpath(raw,f), "N_selected"=>string(n)) for (f,n) in (("alpha.sxm",3),("beta.sxm",5))])
    templates = joinpath(dir, "templates.tsv"); write(templates, "synthetic\n")
    return Dict("--data-dir"=>raw, "--count-config"=>joinpath(ROOT,"config","chitosan.toml"),
        "--config"=>joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
        "--selected-summary"=>summary, "--templates"=>templates, "--outdir"=>joinpath(dir,"output"))
end

function synthetic_count_rows(counts)
    [Dict("file"=>f, "lobe"=>string(k), "N"=>string(n), "gcv"=>"0.125", "source"=>"ell")
        for (f,n) in sort(collect(counts)) for k in 1:n]
end
const COUNT_HEADER = ["file", "lobe", "N", "gcv", "source"]

@testset "Original raw GCV rule and frozen settings" begin
    E = NativeGCVExtractor
    result(n,gcv; valid=true, success=true) = (;n,gcv,valid,success,bic=-gcv,aicc=-gcv,cv_nll_mean=-gcv)
    circ = [result(3,.4),result(4,.2),result(5,.1;valid=false),result(2,.01;success=false)]
    ell = [result(3,.3),result(4,.25),result(7,Inf)]
    n,r,source = E._effective_best(E._best_by_n(ell,"gcv"),E._best_by_n(circ,"gcv"),"gcv")
    @test n == 4 && r.gcv == .2 && source == "circ"
    @test E._effective_best(Dict(),Dict(),"gcv") == (0,nothing,"NA")
    n,r,source = E._effective_best(Dict(3=>result(3,.2),4=>result(4,.2)),Dict(3=>result(3,.2)),"gcv")
    @test n == 3 && source == "ell" # Existing deterministic N/source tie rules.
    cfg = TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
    _,free,freec = E._configs(cfg["model"],cfg["preprocessing"],"unused")
    _,fixed,_ = E._configs(cfg["model"],cfg["preprocessing"],"unused";selected_context=(n=9,use_rescue=false))
    @test free.n_min == 2 && free.n_max == 14
    @test free.intelligent_sweep && free.fuse_z_bwd && freec.chain_circular_sigmas
    @test Set(k for k in fieldnames(typeof(free)) if !isequal(getfield(free,k),getfield(fixed,k))) == Set([:n_min,:n_max])
end

@testset "Strict metadata, cohort, and selected-fit checks" begin
    for flag in ("--truth","--expected-N","--reference","--manifest","--features","--split-features","--selection-policy")
        @test_throws ErrorException parse_reselection([flag,"forbidden"])
    end
    mktempdir() do dir
        opts = count_fixture(dir)
        args = vcat([[k,v] for (k,v) in opts]...)
        @test parse_reselection(args) == opts
        @test_throws ErrorException parse_reselection(vcat(args,["--config",opts["--config"]]))
        execute_reselection(merge(opts,Dict("--dry-run"=>"true")))
        @test !ispath(opts["--outdir"])
        path = joinpath(dir,"features.tsv")
        out = joinpath(dir,"new.tsv")
        rows = synthetic_count_rows(Dict("alpha.sxm"=>4,"beta.sxm"=>2))
        write_table(path, COUNT_HEADER, rows)
        expected = Dict("alpha.sxm"=>4,"beta.sxm"=>2)
        @test reselection_summary(path, opts["--data-dir"], keys(expected),out) == expected
        @test_throws ErrorException reselection_summary(path,dir,["alpha.sxm"],joinpath(dir,"cohort_bad.tsv"))
        for (index,(column,value)) in enumerate((("N","8"),("gcv","NaN"),("gcv","-1"),("source","unknown"),("lobe","9")))
            bad = deepcopy(rows); bad[1][column] = value
            badpath = joinpath(dir,"bad_$index.tsv")
            write_table(badpath,COUNT_HEADER,bad)
            @test_throws ErrorException reselection_summary(badpath,dir,keys(expected),joinpath(dir,"bad_$(index)_out.tsv"))
        end
        duplicate = joinpath(dir,"duplicate.tsv")
        write_table(duplicate,COUNT_HEADER,vcat(rows,[rows[1]]))
        @test_throws ErrorException reselection_summary(duplicate,dir,keys(expected),joinpath(dir,"duplicate_out.tsv"))
        write(joinpath(opts["--data-dir"],"extra.sxm"),"synthetic")
        @test_throws ErrorException execute_reselection(merge(opts,Dict("--dry-run"=>"true")))
        @test !ispath(opts["--outdir"])
    end
end

@testset "Two fresh fixed-N arms; sweep contributes counts only" begin
    mktempdir() do dir
        opts = count_fixture(dir)
        calls = []
        function fake_runner(outdir,name,script,args)
            push!(calls,(name,script,copy(args)))
            value(key) = args[findfirst(==(key),args)+1]
            counts = selected_counts(value("--selected-summary"))
            dest = value("--outdir"); mkdir(dest)
            rows = synthetic_count_rows(counts)
            for file in ("features.tsv","features_split.tsv","predictions.tsv")
                write_table(joinpath(dest,file),COUNT_HEADER,rows)
            end
        end
        function fake_exporter(outdir,name,script,args,output; nfiles)
            @test nfiles == 2
            @test name == "gcv_sweep" && script == "extract_lobe_features.jl"
            @test args == ["--data-dir",joinpath(opts["--outdir"],"raw_inputs"),"--config",opts["--count-config"]]
            write_table(output,COUNT_HEADER,synthetic_count_rows(Dict("alpha.sxm"=>4,"beta.sxm"=>2)))
        end
        execute_reselection(opts; exporter=fake_exporter,runner=fake_runner)
        @test [c[1] for c in calls] == ["control","reselected"]
        for (_,script,args) in calls
            @test script == "run_reconstructed_chitosan.jl"
            @test !any(k -> k in args,("--features","--split-features","--patches-fwd","--manifest"))
            @test args[findfirst(==("--config"),args)+1] == opts["--config"]
        end
        @test selected_counts(joinpath(opts["--outdir"],"control_counts.tsv")) == Dict("alpha.sxm"=>3,"beta.sxm"=>5)
        @test selected_counts(joinpath(opts["--outdir"],"reselected_counts.tsv")) == Dict("alpha.sxm"=>4,"beta.sxm"=>2)
        _,changes = read_table(joinpath(opts["--outdir"],"count_changes.tsv"))
        @test [r["delta"] for r in changes] == ["1","-3"]
        @test !isfile(joinpath(opts["--outdir"],"failures.tsv"))
        @test_throws ErrorException execute_reselection(opts;exporter=fake_exporter,runner=fake_runner)
        failopts = merge(opts,Dict("--outdir"=>joinpath(dir,"failed")))
        @test_throws ErrorException execute_reselection(failopts;runner=(args...)->error("synthetic failure"))
        _,failures = read_table(joinpath(failopts["--outdir"],"failures.tsv"))
        @test only(failures)["stage"] == "control"
        @test !ispath(joinpath(failopts["--outdir"],"reselected"))
    end
end

@testset "One bounded Slurm job with metadata-only dry run" begin
    @test success(`bash -n $RESELECT_JOB`)
    mktempdir() do dir
        opts = count_fixture(dir)
        input = joinpath(dir,"inputs"); mkdir(input)
        symlink(opts["--data-dir"],joinpath(input,"full146_raw"))
        cp(opts["--templates"],joinpath(input,"templates_cc.tsv"))
        fake = joinpath(dir,"fake julia")
        write(fake,raw"""#!/usr/bin/env bash
printf '%s\0' "$OPENBLAS_NUM_THREADS" "$JULIA_NUM_THREADS" "$@" > "${RESELECT_CAPTURE:?}"
""")
        chmod(fake,0o755)
        capture = joinpath(dir,"capture")
        env = ["STMFIT_PROJECT_DIR"=>ROOT,"STMFIT_INPUT_DIR"=>input,
            "STMFIT_SELECTED_SUMMARY"=>opts["--selected-summary"],"STMFIT_OUTDIR"=>opts["--outdir"],
            "JULIA_BIN"=>fake,"RESELECT_CAPTURE"=>capture,"SLURM_JOB_ID"=>nothing,"SLURM_CPUS_PER_TASK"=>nothing]
        @test success(addenv(`bash $RESELECT_JOB --dry-run`,env...))
        argv = split(read(capture,String),'\0';keepempty=false)
        @test argv == ["1","4","--startup-file=no","-t","4","--project=.","test/run_gcv_reselection.jl",
            "--data-dir",joinpath(input,"full146_raw"),"--count-config","config/chitosan.toml",
            "--config","config/unit_assignment_patch_support.toml","--selected-summary",opts["--selected-summary"],
            "--templates",joinpath(input,"templates_cc.tsv"),"--outdir",opts["--outdir"],"--dry-run"]
        @test !ispath(opts["--outdir"])
        @test !success(pipeline(addenv(`bash $RESELECT_JOB`,env...);stderr=devnull))
        @test success(addenv(`bash $RESELECT_JOB`,env...,"SLURM_JOB_ID"=>"synthetic","SLURM_CPUS_PER_TASK"=>"4"))
        @test !occursin("--dry-run",read(capture,String))
        @test !success(pipeline(addenv(`bash $RESELECT_JOB`,env...,"SLURM_JOB_ID"=>"synthetic","SLURM_CPUS_PER_TASK"=>"2");stderr=devnull))
        @test !success(pipeline(addenv(`bash $RESELECT_JOB --truth`,env...);stderr=devnull))
        mkdir(opts["--outdir"])
        @test !success(pipeline(addenv(`bash $RESELECT_JOB --dry-run`,env...);stderr=devnull))
    end
end
