using Test, SHA
include(joinpath(@__DIR__, "run_masked_preprocessing.jl"))
const MPD = MaskedPreprocessingDriver
const ROOT = dirname(@__DIR__)

function write_table(path, header, rows)
    open(path, "w") do io
        println(io, join(header, '\t'))
        for row in rows; println(io, join(row, '\t')); end
    end
end

function fixture(dir)
    input = joinpath(dir, "inputs with spaces and 'quotes'")
    mkpath(joinpath(input, "raw"))
    project = joinpath(dir, "fake project with spaces and 'quotes'")
    mkpath(joinpath(project, "test"))
    config = joinpath(dir, "physical \"quoted\".toml")
    acquisition = joinpath(dir, "acquisition settings.toml")
    settings = joinpath(dir, "masked settings.toml")
    # Synthetic metadata only; no physical estimator or production reader imports.
    write(config, "[model]\n[selection]\n[preprocessing]\n")
    write(acquisition, "[acquisition_noise]\n")
    write(settings, "[masked_robust_preprocessing]\n")
    names = collect(MPD.FROZEN_FILES)
    counts = [3, 4, 5, 6] # Deliberately NOT the original saved N values or truth.
    policies = ["adaptive_support_rescue_keep_robust_guard", "adaptive_support_rescue_keep",
                "adaptive_support_rescue", "adaptive_support_rescue_keep"]
    sample = joinpath(input, "diagnostic_sample.tsv")
    selected = joinpath(input, "selected_summary.tsv")
    geometry = joinpath(input, "base_geometry.tsv")
    sample_header = ["file", "N_selected", "refined_policy", "delta_GCV_rel_eff", "inclusion_reason"]
    sample_rows = [(names[i], counts[i], policies[i], "not used as a score", "frozen") for i in 1:4]
    selected_header = ["filepath", "N_selected", "refined_policy"]
    selected_rows = [(joinpath("original saved paths", names[i]), counts[i], policies[i]) for i in 1:4]
    # The actual input metadata stays at 25 files, not trimmed to the four raws.
    append!(selected_rows, [("unused_$i.sxm", 2, "saved_policy") for i in 1:21])
    geometry_header = ["file", "N", "lobe"]
    geometry_rows = [(names[i], counts[i], j) for i in 1:4 for j in 1:counts[i]]
    append!(geometry_rows, [("unused_$i.sxm", 2, j) for i in 1:21 for j in 1:2])
    write_table(sample, sample_header, sample_rows)
    write_table(selected, selected_header, selected_rows)
    write_table(geometry, geometry_header, geometry_rows)
    raw = [joinpath(input, "raw", name) for name in names]
    foreach(p -> write(p, "SYNTHETIC SENTINEL, NOT AN SXM IMAGE"), raw)
    out = joinpath(dir, "new output with spaces and 'quotes'")
    args = ["--input-dir", input, "--config", config, "--acquisition-settings", acquisition,
            "--settings", settings, "--outdir", out]
    return (; input, project, config, acquisition, settings, names, counts, policies,
            sample, selected, geometry, sample_header, sample_rows, selected_header,
            selected_rows, geometry_header, geometry_rows, raw, out, args)
end

function replaced(f, path, text)
    original = read(path, String)
    try
        write(path, text)
        return f()
    finally
        write(path, original)
    end
end

function changed_rows(f, path, header, rows)
    original = read(path, String)
    try
        write_table(path, header, rows)
        return f()
    finally
        write(path, original)
    end
end

filearg(cmd) = basename(cmd.exec[findfirst(==("--file"), cmd.exec)+1])

@testset "Masked driver options and metadata-only preparation" begin
    @test v"1.13.0" <= VERSION < v"1.14.0"
    @test Threads.nthreads() == 1
    @test MPD.options(["--help"]) === nothing
    for args in (String[], ["--input-dir"], ["--config", "--settings"],
                 ["--truth", "labels.tsv"], ["--expected-N", "6"], ["--sequence", "NKNNKN"],
                 ["--count-prior", "6"], ["--help", "--dry-run"], ["positional"])
        @test_throws ErrorException MPD.options(args)
    end
    mktempdir() do dir
        f = fixture(dir)
        dryargs = vcat(f.args, ["--dry-run"])
        protected = vcat([f.sample, f.selected, f.geometry, f.config, f.acquisition, f.settings], f.raw)
        before = Dict(p => sha256(read(p)) for p in protected)
        plan = MPD.prepare(dryargs; project=f.project)
        @test plan.o.dry
        @test length(plan.cases) == 4
        @test [c.file for c in plan.cases] == f.names
        @test length(MPD.table(f.selected, ["filepath"])) == 25
        @test length(readdir(joinpath(f.input, "raw"))) == 4
        for (i, c) in enumerate(plan.cases)
            expected = [joinpath(Sys.BINDIR, Base.julia_exename()), "--startup-file=no", "--threads=1",
                        "--project=$(f.project)", joinpath(f.project, "test", "diagnose_masked_preprocessing_real.jl"),
                        "--file", f.raw[i], "--geometry", f.geometry, "--selected-summary", f.selected,
                        "--config", f.config, "--acquisition-settings", f.acquisition, "--settings", f.settings,
                        "--outdir", joinpath(f.out, splitext(f.names[i])[1])]
            @test c.cmd.exec == expected
            @test c.out == last(expected)
            @test !any(t -> occursin("diagnose_counting", t) || occursin("diagnose_fisher", t) ||
                           occursin("extract_", t) || occursin("grade", t), c.cmd.exec)
        end
        # A dry run invokes neither the real CLI nor an injected executor. Invalid
        # SXM sentinel bytes and absent fake CLI make accidental science visible.
        calls = Ref(0)
        drylog = joinpath(dir, "metadata dry-run.log")
        withenv("SLURM_JOB_ID" => nothing, "SLURM_CPUS_PER_TASK" => nothing) do
            open(drylog, "w") do io
                redirect_stdout(io) do
                    @test MPD.execute(plan; executor=(cmd, log) -> (calls[] += 1)) === nothing
                end
            end
        end
        @test calls[] == 0
        @test !ispath(f.out)
        text = read(drylog, String)
        @test count(line -> occursin("diagnose_masked_preprocessing_real.jl", line), split(text, '\n')) == 4
        @test occursin("no subprocesses, outputs or SXM pixels read", text)
        @test all(sha256(read(path)) == hash for (path, hash) in before)
        @test_throws ErrorException MPD.options(vcat(dryargs, ["--dry-run"]))
        @test_throws ErrorException MPD.options(vcat(f.args, ["--config", f.config]))
        @test_throws ErrorException MPD.options(vcat(f.args[1:end-1], [""]))
        for key in MPD.PATH_OPTIONS
            index = findfirst(==(key), f.args)
            incomplete = vcat(f.args[1:index-1], f.args[index+2:end])
            @test_throws ErrorException MPD.options(incomplete)
        end
        # Retain frozen manifest order, do not rank by score or inclusion reason.
        changed_rows(f.sample, f.sample_header, reverse(f.sample_rows)) do
            p = MPD.prepare(dryargs; project=f.project)
            @test [c.file for c in p.cases] == reverse(f.names)
        end
        # No real job or CPU allocation is implied by these synthetic guards.
        withenv("SLURM_JOB_ID" => nothing) do
            @test_throws ErrorException MPD.execute(MPD.prepare(f.args; project=f.project))
        end
        for value in (nothing, "", "2", "invalid")
            withenv("SLURM_JOB_ID" => "synthetic-test", "SLURM_CPUS_PER_TASK" => value) do
                @test_throws ErrorException MPD.execute(MPD.prepare(f.args; project=f.project))
            end
        end
        @test !ispath(f.out)
    end
end

@testset "Frozen metadata mismatches and exclusive output" begin
    mktempdir() do dir
        f = fixture(dir)
        check() = MPD.prepare(vcat(f.args, ["--dry-run"]); project=f.project)
        for path in (f.sample, f.selected, f.geometry, f.config, f.acquisition, f.settings, f.raw[1])
            mv(path, path * ".saved")
            try
                @test_throws ErrorException check()
            finally
                mv(path * ".saved", path)
            end
        end
        for (path, text) in ((f.config, "[model]\n[preprocessing]\n"),
                             (f.acquisition, "[wrong_section]\n"), (f.settings, "[wrong_section]\n"))
            replaced(path, text) do; @test_throws ErrorException check(); end
        end
        for rows in (f.sample_rows[1:3], vcat(f.sample_rows, f.sample_rows[1:1]),
                     [f.sample_rows[1], f.sample_rows[1], f.sample_rows[3], f.sample_rows[4]])
            changed_rows(f.sample, f.sample_header, rows) do; @test_throws ErrorException check(); end
        end
        for name in ("subdir/" * f.names[1], "../" * f.names[1], "substitute.sxm")
            rows = copy(f.sample_rows)
            rows[1] = (name, f.counts[1], f.policies[1], "ignored", "frozen")
            changed_rows(f.sample, f.sample_header, rows) do; @test_throws ErrorException check(); end
        end
        for count in (0, 99)
            rows = copy(f.sample_rows); rows[1] = (f.names[1], count, f.policies[1], "ignored", "frozen")
            changed_rows(f.sample, f.sample_header, rows) do; @test_throws ErrorException check(); end
        end
        for policy in ("", "adaptive_support_rescue")
            rows = copy(f.sample_rows); rows[1] = (f.names[1], f.counts[1], policy, "ignored", "frozen")
            changed_rows(f.sample, f.sample_header, rows) do; @test_throws ErrorException check(); end
        end
        for rows in (f.selected_rows[2:end], vcat(f.selected_rows, [("other/" * f.names[1], 3, f.policies[1])]))
            changed_rows(f.selected, f.selected_header, rows) do; @test_throws ErrorException check(); end
        end
        for (n, policy) in ((0, f.policies[1]), (3, ""), (3, "different_saved_policy"))
            rows = copy(f.selected_rows); rows[1] = (f.names[1], n, policy)
            changed_rows(f.selected, f.selected_header, rows) do; @test_throws ErrorException check(); end
        end
        for rows in (f.geometry_rows[2:end], vcat(f.geometry_rows, f.geometry_rows[1:1]))
            changed_rows(f.geometry, f.geometry_header, rows) do; @test_throws ErrorException check(); end
        end
        for row in ((f.names[1], 99, 1), (f.names[1], 3, 2), (f.names[1], 3, 0))
            rows = copy(f.geometry_rows); rows[1] = row
            changed_rows(f.geometry, f.geometry_header, rows) do; @test_throws ErrorException check(); end
        end
        for text in ("", "file\tfile\n", "wrong\nvalue\n", "file\tN\tlobe\nonly_one_field\n",
                     "file\tN\tlobe\ttruth\n", "file\tN\tlobe\tEXPECTED_N\n")
            replaced(f.geometry, text) do; @test_throws ErrorException check(); end
        end
        # Preexisting directory, regular file, dangling symlink and live symlink.
        for kind in (:directory, :file, :dangling, :live)
            if kind == :directory
                mkdir(f.out)
            elseif kind == :file
                write(f.out, "occupied")
            else
                symlink(kind == :live ? f.input : joinpath(dir, "absent"), f.out)
            end
            try
                @test_throws ErrorException check()
            finally
                rm(f.out; recursive=true)
            end
        end
        plan = check()
        symlink(joinpath(dir, "absent"), f.out)
        try
            @test_throws ErrorException MPD.execute(plan)
        finally
            rm(f.out)
        end
        otherargs = copy(f.args); otherargs[end] = joinpath(dir, "missing parent", "new")
        @test_throws ErrorException MPD.prepare(otherargs)
        @test !ispath(f.out)
    end
end

@testset "Every case failure retained and completed rows appended immediately" begin
    mktempdir() do dir
        f = fixture(dir)
        plan = MPD.prepare(f.args; project=f.project)
        seen = String[]
        function fake(cmd, io)
            # This executor does not start any process or read a raw sentinel.
            @test length(readlines(joinpath(f.out, "stages.tsv"))) == length(seen) + 1
            file = filearg(cmd); push!(seen, file)
            println(io, "synthetic executor: ", file)
            file == f.names[1] && return 3
            file == f.names[2] && error("synthetic subprocess spawn failure")
            return 0
        end
        withenv("SLURM_JOB_ID" => "synthetic-test", "SLURM_CPUS_PER_TASK" => "4") do
            @test_throws ErrorException MPD.execute(plan; executor=fake)
        end
        @test seen == f.names
        rows = MPD.table(joinpath(f.out, "stages.tsv"), ["case", "stage", "exit_code", "elapsed_s"])
        @test length(rows) == 4
        @test getindex.(rows, "case") == f.names
        @test getindex.(rows, "exit_code") == ["3", "1", "0", "0"]
        @test all(r -> r["stage"] == "masked_preprocessing", rows)
        @test all(r -> isfinite(parse(Float64, r["elapsed_s"])) && parse(Float64, r["elapsed_s"]) > 0, rows)
        @test length(readdir(joinpath(f.out, "logs"))) == 4
        for c in plan.cases
            text = read(joinpath(f.out, "logs", "$(c.stem)_masked_preprocessing.log"), String)
            @test occursin("synthetic executor: " * c.file, text)
        end
        @test occursin("synthetic subprocess spawn failure",
                       read(joinpath(f.out, "logs", "$(plan.cases[2].stem)_masked_preprocessing.log"), String))
        @test all(!ispath(c.out) for c in plan.cases) # Driver never fabricates CLI outputs.
        @test_throws ErrorException MPD.prepare(f.args; project=f.project)
    end
end

@testset "Four asynchronous one-thread command slots, no more" begin
    mktempdir() do dir
        f = fixture(dir)
        started = Channel{String}(4)
        releases = Dict(name => Channel{Nothing}(1) for name in f.names)
        active = Ref(0); peak = Ref(0); calls = Ref(0)
        function fake(cmd, io)
            active[] += 1; calls[] += 1; peak[] = max(peak[], active[])
            file = filearg(cmd)
            put!(started, file)
            take!(releases[file])
            println(io, "synthetic completed ", file)
            active[] -= 1
            return 0
        end
        withenv("SLURM_JOB_ID" => "synthetic-test", "SLURM_CPUS_PER_TASK" => "8") do
            task = @async MPD.execute(MPD.prepare(f.args; project=f.project); executor=fake)
            names = [take!(started) for _ in 1:4]
            @test Set(names) == Set(f.names)
            @test active[] == 4
            foreach(name -> put!(releases[name], nothing), names)
            outcomes = fetch(task)
            @test length(outcomes) == 4
            @test all(row -> row[3] == 0, outcomes)
        end
        @test peak[] == 4
        @test active[] == 0
        @test calls[] == 4
    end
end

@testset "Native subprocess stand-ins preserve args/logs and nonzero overall status" begin
    mktempdir() do dir
        f = fixture(dir)
        fake = """
        @assert v"1.13.0" <= VERSION < v"1.14.0"
        @assert Threads.nthreads() == 1
        out = ARGS[findfirst(==("--outdir"), ARGS)+1]
        file = basename(ARGS[findfirst(==("--file"), ARGS)+1])
        mkdir(out)
        println("NATIVE SYNTHETIC STDOUT ", file)
        println(stderr, "NATIVE SYNTHETIC STDERR ", file)
        foreach(arg -> println(repr(arg)), ARGS)
        exit(file == $(repr(f.names[1])) ? 7 : 0)
        """
        write(joinpath(f.project, "test", "diagnose_masked_preprocessing_real.jl"), fake)
        # Real native process plumbing, but the only executed CLI is this tiny
        # stand-in in a temporary directory. No actual scientific script is run.
        expression = "include(" * repr(joinpath(@__DIR__, "run_masked_preprocessing.jl")) * "); " *
                     "MaskedPreprocessingDriver.execute(MaskedPreprocessingDriver.prepare(ARGS; project=" * repr(f.project) * "))"
        julia = joinpath(Sys.BINDIR, Base.julia_exename())
        cmd = Cmd(vcat([julia, "--startup-file=no", "--threads=1", "--project=$ROOT", "-e", expression, "--"], f.args))
        log = joinpath(dir, "native stand-in driver.log")
        process = withenv("SLURM_JOB_ID" => "synthetic-test", "SLURM_CPUS_PER_TASK" => "4") do
            open(log, "w") do io
                run(pipeline(ignorestatus(cmd), stdout=io, stderr=io))
            end
        end
        @test process.exitcode != 0
        @test occursin("Some masked-preprocessing cases failed", read(log, String))
        rows = MPD.table(joinpath(f.out, "stages.tsv"), ["case", "stage", "exit_code", "elapsed_s"])
        @test length(rows) == 4
        @test Set(getindex.(rows, "case")) == Set(f.names)
        codes = Dict(row["case"] => row["exit_code"] for row in rows)
        @test codes[f.names[1]] == "7"
        @test all(codes[file] == "0" for file in f.names[2:end])
        for (i, file) in enumerate(f.names)
            stem = splitext(file)[1]
            text = read(joinpath(f.out, "logs", "$(stem)_masked_preprocessing.log"), String)
            @test occursin("NATIVE SYNTHETIC STDOUT " * file, text)
            @test occursin("NATIVE SYNTHETIC STDERR " * file, text)
            @test occursin(repr(f.raw[i]), text)
            @test occursin(repr(f.config), text)
            @test occursin(repr(joinpath(f.out, stem)), text)
            @test isdir(joinpath(f.out, stem))
        end
        @test length(readdir(joinpath(f.out, "logs"))) == 4
    end
end

@testset "SBATCH bounds, shell quoting and cleared-environment handoff (not Slurm export validation)" begin
    script = joinpath(ROOT, "hpc", "masked_preprocessing.sbatch")
    text = read(script, String)
    @test success(`bash -n $script`)
    for directive in ("--job-name=stmfit-masked-background", "--nodes=1", "--ntasks=1",
                      "--cpus-per-task=4", "--mem=16000M", "--time=00:30:00")
        @test occursin("#SBATCH " * directive, text)
    end
    for line in ("GKSwstype=100", "JULIA_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1")
        @test occursin("export " * line, text)
    end
    @test occursin("SBATCH_EXPORT=NONE", text)
    @test occursin("--export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN", text)
    active_lines = join(filter(line -> !startswith(lstrip(line), "#"), split(text, '\n')), '\n')
    @test !occursin(r"\b(sbatch|srun|ssh|rsync|scontrol)\b", active_lines)
    for banned in ("run_label_free_exploration", "diagnose_counting_variable_projection", "diagnose_fisher",
                   "extract_unknown", "--watch", "--array", "--requeue", "11806180", "11812202")
        @test !occursin(banned, text)
    end
    # Driver is stdlib-only and contains exactly one scientific entrypoint path.
    driver = read(joinpath(@__DIR__, "run_masked_preprocessing.jl"), String)
    @test !occursin(r"(?m)^\s*(include\(|using (STMSXMIO|GaussianFit|STMMolecular))", driver)
    @test !occursin("read_sxm(", driver)
    @test count(line -> occursin("diagnose_masked_preprocessing_real.jl", line), split(driver, '\n')) == 1
    mktempdir() do dir
        project = joinpath(dir, "project with spaces and 'quotes'"); mkdir(project)
        input = joinpath(dir, "input with spaces and 'quotes'")
        out = joinpath(dir, "output with spaces and 'quotes'")
        cleared = Dict("PATH" => ENV["PATH"], "HOME" => get(ENV, "HOME", dir), "SBATCH_EXPORT" => "NONE")
        supplied = merge(cleared, Dict("STMFIT_PROJECT_DIR" => project, "STMFIT_INPUT_DIR" => input,
                         "STMFIT_OUTDIR" => out, "JULIA_BIN" => "/bin/echo"))
        function shell_run(env, args, name)
            log = joinpath(dir, name * ".log")
            cmd = Cmd(vcat(["bash", script], args))
            p = open(log, "w") do io
                run(pipeline(ignorestatus(setenv(cmd, env)), stdout=io, stderr=io))
            end
            return p.exitcode, read(log, String)
        end
        code, msg = shell_run(cleared, ["--dry-run"], "empty-env")
        @test code != 0
        @test occursin("STMFIT_PROJECT_DIR: required", msg)
        for key in ("STMFIT_PROJECT_DIR", "STMFIT_INPUT_DIR", "STMFIT_OUTDIR", "JULIA_BIN")
            env = copy(supplied); delete!(env, key)
            code, msg = shell_run(env, ["--dry-run"], "missing_" * key)
            @test code != 0
            @test occursin(key * ":", msg)
        end
        code, msg = shell_run(supplied, ["--dry-run"], "echo")
        @test code == 0
        @test occursin("--threads=1", msg)
        @test occursin("--input-dir " * input, msg)
        @test occursin("--outdir " * out, msg)
        @test endswith(strip(msg), "--dry-run")
        @test !ispath(out)
        for args in (["--help"], ["--retry"], ["--dry-run", "--dry-run"])
            code, msg = shell_run(supplied, args, "invalid_" * string(length(args)) * first(args))
            @test code == 2
            @test occursin("Only optional --dry-run", msg)
        end
        # Unlike echo, this shell stand-in exposes exact argument boundaries and
        # exported thread/plot variables. It does not emulate a compute node.
        stub = joinpath(dir, "argument recorder with 'quotes'")
        write(stub, "#!/bin/bash\nprintf '%s\\n' \"\$@\"\nprintf 'ENV:%s:%s:%s:%s\\n' \"\$JULIA_NUM_THREADS\" \"\$OPENBLAS_NUM_THREADS\" \"\$OMP_NUM_THREADS\" \"\$GKSwstype\"\n")
        chmod(stub, 0o700)
        env = merge(supplied, Dict("JULIA_BIN" => stub))
        expected = ["--startup-file=no", "--threads=1", "--project=$project", "test/run_masked_preprocessing.jl",
                    "--input-dir", input, "--config", "config/chitosan_10_20mer_adaptive_support_rescue.toml",
                    "--acquisition-settings", "config/label_free_exploration.toml", "--settings",
                    "config/masked_robust_preprocessing.toml", "--outdir", out]
        for args in (String[], ["--dry-run"])
            code, msg = shell_run(env, args, "argv_" * string(length(args)))
            @test code == 0
            lines = split(chomp(msg), '\n')
            @test lines[1:end-1] == vcat(expected, args)
            @test lines[end] == "ENV:1:1:1:100"
        end
        @test !ispath(out)
    end
end
