using Test, TOML, SHA
include(joinpath(@__DIR__, "run_label_free_exploration.jl"))
const LFE=LabelFreeExploration

function write_table(path, header, rows)
    open(path,"w") do io
        println(io,join(header,'\t'))
        for row in rows; println(io,join(row,'\t')); end
    end
end

@testset "Bounded label-free exploration driver" begin
    @test LFE.options(["--help"])===nothing
    @test_throws ErrorException LFE.options(["--truth","anything"])
    @test_throws ErrorException LFE.options(String[])
    @test_throws ErrorException LFE.options(["--input-dir"])
    mktempdir() do dir
        input=joinpath(dir,"inputs with spaces"); mkpath(joinpath(input,"raw"))
        config=joinpath(dir,"model.toml")
        write(config,"[model]\nn_min=2\nn_max=24\n[selection]\n[preprocessing]\n")
        settings=joinpath(dir,"settings.toml")
        cp(joinpath(@__DIR__,"..","config","label_free_exploration.toml"),settings)
        fishroot=joinpath(input,"fisher"); mkdir(fishroot)
        fishcfg=joinpath(fishroot,"unit_assignment_reconstructed.toml")
        cp(joinpath(@__DIR__,"..","config","unit_assignment_reconstructed.toml"),fishcfg)
        fishinputs=[fishcfg]
        for cohort in ("unknown25","full146")
            mkdir(joinpath(fishroot,cohort))
            for name in ("patches_fwd17.tsv","fisher_cv.tsv")
                path=joinpath(fishroot,cohort,name)
                write(path,"synthetic dry-run input, not fitted patches")
                push!(fishinputs,path)
            end
        end
        chosen=joinpath(input,"candidate_counts.tsv")
        selected=joinpath(input,"selected_summary.tsv")
        geometry=joinpath(input,"base_geometry.tsv")
        names=["case_$i.sxm" for i in 1:4]
        counts=[3,4,5,6]
        sels=[(names[i],counts[i],i==1 ? 5 : counts[i],counts[i]+1,"adaptive_support_rescue_keep") for i in 1:4]
        write_table(selected,["filepath","N_selected","N_eff","runnerup_N_eff","refined_policy"],sels)
        write_table(geometry,["file","N","lobe"],[(names[i],counts[i],j) for i in 1:4 for j in 1:counts[i]])
        candidates=[(names[i],join(i==1 ? [2,3,4,5] : collect(counts[i]-1:counts[i]+1),','),"saved_union") for i in 1:4]
        write_table(chosen,["file","candidate_ns","source"],candidates)
        for file in names; write(joinpath(input,"raw",file),"not an image: dry run must not parse pixels"); end
        out=joinpath(dir,"new outputs")
        args=["--input-dir",input,"--config",config,"--settings",settings,"--outdir",out,"--dry-run"]
        before=Dict(p=>sha256(read(p)) for p in vcat([selected,geometry,chosen,config,settings],fishinputs))
        p=LFE.prepare(args)
        @test length(p.cases)==4
        @test [f.cohort for f in p.fisher]==["unknown25","full146"]
        @test all("res" in f.cmd.exec && "--prefix" in f.cmd.exec for f in p.fisher)
        @test all(any(endswith("patches_fwd17.tsv"),f.cmd.exec) for f in p.fisher)
        @test all(fishcfg in f.cmd.exec for f in p.fisher)
        @test all("--threads=1" in f.cmd.exec for f in p.fisher)
        @test p.o.dry
        @test p.cases[1].file==names[1]
        @test length(p.cases[1].variable_projection.exec)==17
        @test "2,3,4,5" in p.cases[1].variable_projection.exec
        @test "--threads=1" in p.cases[1].acquisition.exec
        @test p.paths.geometry in p.cases[1].acquisition.exec
        @test "--selected-summary" in p.cases[1].variable_projection.exec
        @test all(!occursin("grade",s) for c in p.cases for s in c.variable_projection.exec)
        open(joinpath(dir,"dry.log"),"w") do io
            redirect_stdout(io) do; LFE.execute(p); end
        end
        @test !ispath(out)
        @test all(sha256(read(k))==v for (k,v) in before)
        @test occursin("four files",read(joinpath(dir,"dry.log"),String))
        @test_throws ErrorException LFE.options(vcat(args,["--dry-run"]))
        @test_throws ErrorException LFE.options(vcat(args,["--config",config]))
        nodry=args[1:end-1]
        withenv("SLURM_JOB_ID"=>nothing) do
            @test_throws ErrorException LFE.execute(LFE.prepare(nodry))
        end
        @test !ispath(out)
        withenv("SLURM_JOB_ID"=>"synthetic", "SLURM_CPUS_PER_TASK"=>"2") do
            @test_throws ErrorException LFE.execute(LFE.prepare(nodry))
        end
        @test !ispath(out)
        # Exercise bounded subprocess execution with tiny native Julia stand-ins,
        # not actual fit programs or manufactured scientific results.
        fakeproject=joinpath(dir,"fake project"); mkpath(joinpath(fakeproject,"test"))
        fake_acquisition="j=findfirst(==(\"--outdir\"),ARGS); mkdir(ARGS[j+1]); println(\"synthetic acquisition\"); exit(3)"
        fake_profile="j=findfirst(==(\"--outdir\"),ARGS); mkdir(ARGS[j+1]); println(\"synthetic profile\"); exit(0)"
        write(joinpath(fakeproject,"test","diagnose_acquisition_noise.jl"),fake_acquisition)
        write(joinpath(fakeproject,"test","diagnose_counting_variable_projection.jl"),fake_profile)
        fake_fisher="j=findfirst(==(\"--outdir\"),ARGS); mkdir(ARGS[j+1]); println(\"synthetic Fisher\"); exit(occursin(\"unknown25\",ARGS[j+1]) ? 5 : 0)"
        write(joinpath(fakeproject,"test","diagnose_fisher_attribution.jl"),fake_fisher)
        withenv("SLURM_JOB_ID"=>"synthetic", "SLURM_CPUS_PER_TASK"=>"4") do
            @test_throws ErrorException LFE.execute(LFE.prepare(nodry;project=fakeproject))
        end
        statuses=LFE.table(joinpath(out,"stages.tsv"),["case","stage","exit_code","elapsed_s"])
        @test length(statuses)==10
        @test count(r->r["exit_code"]=="3",statuses)==4
        @test count(r->r["exit_code"]=="0",statuses)==5
        @test count(r->r["exit_code"]=="5",statuses)==1
        @test all(r["stage"]!="fisher_replay" for r in statuses[1:8])
        @test [r["case"] for r in statuses[9:10]]==["unknown25","full146"]
        @test all(isdir(joinpath(out,"fisher",c)) for c in ("unknown25","full146"))
        @test all(isdir(joinpath(out,splitext(n)[1],"acquisition")) for n in names)
        @test all(isdir(joinpath(out,splitext(n)[1],"variable_projection")) for n in names)
        @test all(sha256(read(k))==v for (k,v) in before)
        @test length(readdir(joinpath(out,"logs")))==10
        rm(out;recursive=true)
        missing=fishinputs[2]; mv(missing,missing*".backup")
        @test_throws ErrorException LFE.prepare(args)
        mv(missing*".backup",missing)
        mkdir(out); @test_throws ErrorException LFE.prepare(args); rm(out)
        write(out,"occupied"); @test_throws ErrorException LFE.prepare(args); rm(out)
        symlink(joinpath(dir,"absent"),out); @test_throws ErrorException LFE.prepare(args); rm(out)
        wrong=copy(candidates); wrong[1]=(names[1],"2,3,4","missing_saved_effective")
        write_table(chosen,["file","candidate_ns","source"],wrong)
        @test_throws ErrorException LFE.prepare(args)
        wrong[1]=(names[1],"2,3,4,5,5","duplicate")
        write_table(chosen,["file","candidate_ns","source"],wrong)
        @test_throws ErrorException LFE.prepare(args)
        write_table(chosen,["file","candidate_ns","source"],candidates[1:3])
        @test_throws ErrorException LFE.prepare(args)
        wrong=copy(candidates); wrong[2]=candidates[1]
        write_table(chosen,["file","candidate_ns","source"],wrong)
        @test_throws ErrorException LFE.prepare(args)
        write_table(chosen,["file","candidate_ns","source"],candidates)
        write_table(geometry,["file","N","lobe"],[(names[i],counts[i],1) for i in 1:4])
        @test_throws ErrorException LFE.prepare(args)
    end
    script=read(joinpath(@__DIR__,"..","hpc","label_free_exploration.sbatch"),String)
    @test occursin("--cpus-per-task=4",script)
    @test occursin("--time=02:00:00",script)
    @test occursin("--mem=16000M",script)
    @test occursin("OPENBLAS_NUM_THREADS=1",script)
    @test !occursin("--watch",script)
    @test !occursin("sbatch ",script)
    @test !occursin("--delete",script)
end
