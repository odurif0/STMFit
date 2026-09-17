# julia -t 4 --project=. test/test_resume_reconstructed_unknown25.jl
# Synthetic metadata and fake export only: no SXM read, optimizer or classifier.
using Test, TOML
const TEST_ROOT = dirname(@__DIR__)
const SCRIPT = joinpath(TEST_ROOT, "hpc", "resume_reconstructed_unknown25.sbatch")
module ResumeScript end
const payload = split(read(SCRIPT, String), "--project=. -e '\n"; limit=2)[2]
const entry = "resume_unknown25()\n'\n"
endswith(payload, entry) || error("Unexpected sbatch entrypoint")
# Strip only the final call; resolve includes from the same root as bash -e.
const definitions = replace(chop(payload; tail=length(entry)),
    "include(\"test/" => "include(\"" * TEST_ROOT * "/test/")
Base.include_string(ResumeScript, definitions)
const R = ResumeScript

function fixture(dir)
    input = joinpath(dir, "input"); raw = joinpath(input, "unknown_raw"); mkpath(raw)
    files = vcat(["260215_022.sxm", "260220_083.sxm"], ["synthetic_$i.sxm" for i in 1:23])
    counts = Dict(file => (i == 25 ? 6 : 9) for (i, file) in enumerate(files))
    summary = "filepath\tN_selected\tstatus\trefined_policy\n"
    base = partial = "file\tN\tlobe\tskew_ratio\n"
    for (i, file) in enumerate(files)
        write(joinpath(raw, file), "Not an SXM; metadata must not read this.\n")
        policy = i <= 2 ? "adaptive_support_rescue" : "adaptive_support_rescue_keep"
        summary *= "/synthetic/$file\t$(counts[file])\tok\t$policy\n"
        rows = join("$file\t$(counts[file])\t$lobe\t1.0\n" for lobe in 1:counts[file])
        base *= rows
        i <= 2 && (partial *= rows)
    end
    paths = (summary=joinpath(input, "unknown25_selected_summary.tsv"),
             base=joinpath(input, "unknown25_base_features.tsv"),
             partial=joinpath(input, "unknown25_split_partial.tsv"))
    for (key, text) in ((:summary, summary), (:base, base), (:partial, partial))
        write(getproperty(paths, key), text)
    end
    write(joinpath(input, "templates_cc.tsv"), "Synthetic placeholder; metadata only.\n")
    return (; input, raw, counts, paths..., out=joinpath(dir, "not-created", "out"))
end
snapshot(input) = Dict(joinpath(dir, file) => read(joinpath(dir, file))
                       for (dir, _, files) in walkdir(input) for file in files)
resume(f; dry=true) = withenv("STMFIT_INPUT_DIR" => f.input, "STMFIT_OUTDIR" => f.out,
                            "STMFIT_RESUME_DRY_RUN" => (dry ? "1" : "0")) do
    R.resume_unknown25()
end
function reject(f, path, text, reason)
    original = read(path, String)
    try
        write(path, text)
        err = try resume(f); nothing catch e; e end
        @test err isa ErrorException
        @test occursin(reason, err === nothing ? "no error" : sprint(showerror, err))
        @test !ispath(dirname(f.out))
    finally
        write(path, original)
    end
end
function shell_check(f; out=f.out, args=["--dry-run"])
    cmd = addenv(`bash $SCRIPT $args`, "STMFIT_PROJECT_DIR" => TEST_ROOT,
        "JULIA_BIN" => joinpath(Sys.BINDIR, "julia"), "STMFIT_INPUT_DIR" => f.input,
        "STMFIT_OUTDIR" => out, "SLURM_CPUS_PER_TASK" => "4")
    io = IOBuffer()
    process = run(pipeline(ignorestatus(cmd); stdout=io, stderr=io))
    return success(process), String(take!(io))
end

@testset "Continuation metadata and shell boundaries" begin
    mktempdir() do dir
        f = fixture(dir); before = snapshot(f.input)
        resume(f)
        @test !ispath(dirname(f.out))
        @test snapshot(f.input) == before
        partial = read(f.partial, String); summary = read(f.summary, String)
        for (text, reason) in (
            (replace(partial, "260215_022.sxm\t9\t9\t1.0\n" => ""), "8/9"),
            (replace(partial, "260215_022.sxm\t9\t3\t1.0\n" => ""), "Noncontiguous"),
            (partial * "260215_022.sxm\t9\t10\t1.0\n", "10/9"),
            (partial * "260215_022.sxm\t9\t1\t1.0\n", "Duplicate file/lobe"),
            (partial * "unexpected.sxm\t1\t1\t1.0\n", "both focused files"),
            (join(filter(l -> !startswith(l, "260215_022.sxm"), split(partial, '\n')), '\n'), "both focused files"))
            reject(f, f.partial, text, reason)
        end
        reject(f, f.base, replace(read(f.base, String), "260215_022.sxm\t9\t9\t1.0\n" => ""), "8/9")
        n_only = join((join(split(line, '\t')[1:2], '\t') for line in split(chomp(summary), '\n')), '\n') * "\n"
        for text in (n_only, replace(summary, "\tadaptive_support_rescue\n" => "\t\n"; count=1),
                     replace(summary, "\tadaptive_support_rescue\n" => "\tambiguous\n"; count=1))
            reject(f, f.summary, text, "Missing or unrecognized adaptive-support refined_policy")
        end
        ok, log = shell_check(f)
        @test ok
        @test occursin("25 files / 222 lobes", log)
        @test !ispath(dirname(f.out))
        try
            write(f.summary, n_only)
            ok, log = shell_check(f)
            @test !ok
            @test occursin("supply the original counting summary, not N alone", log)
            @test !ispath(dirname(f.out))
        finally
            write(f.summary, summary)
        end
        # The shell must refuse both an input directory and an input file as
        # output roots, in dry-run and normal mode, before any real work.
        for (path, args) in ((f.input, ["--dry-run"]), (f.base, String[]))
            ok, log = shell_check(f; out=path, args)
            @test !ok
            @test occursin("Output directory already exists", log)
            @test snapshot(f.input) == before
        end
        mkpath(dirname(f.out))
        symlink(joinpath(dir, "absent"), f.out)
        @test_throws ErrorException resume(f)
        rm(f.out); rm(dirname(f.out))
        @test snapshot(f.input) == before
    end
end

# Replace only this process's subprocess boundary. The native export/merge
# helpers still execute; the fake exporter writes only synthetic fixture rows.
@eval ResumeScript begin
    const fake_calls = []
    const fake_lock = ReentrantLock()
    const fake_base = Ref("")
    function run_stage(outdir, name, script, args; threads=min(Threads.nthreads(), 4))
        lock(fake_lock) do
            push!(fake_calls, (script=script, args=copy(args), threads=threads))
        end
        opts = Dict(args[i] => args[i+1] for i in 1:2:length(args))
        if script == "extract_lobe_features.jl"
            excluded = Set(readlines(opts["--exclude-from"]))
            files = sort(setdiff(collect(keys(selected_counts(opts["--selected-summary"]))), excluded))
            chunk, total = parse.(Int, split(opts["--chunk"], '/'))
            wanted = Set(file for (i, file) in enumerate(files) if mod1(i, total) == chunk)
            header, rows = lobe_table(fake_base[])
            write_table(opts["--out"], header, [rows[k] for k in sort(collect(keys(rows))) if first(k) in wanted])
        elseif script != "run_reconstructed_chitosan.jl"
            error("Unexpected scientific subprocess request: $script")
        end
    end
end
@testset "Continuation wiring with fake exports, no fits/classifiers" begin
    mktempdir() do dir
        f = fixture(dir); before = snapshot(f.input); R.fake_base[] = f.base
        count_config = joinpath(TEST_ROOT, "config", "chitosan_10_20mer_adaptive_support_rescue.toml")
        unit_config = joinpath(TEST_ROOT, "config", "unit_assignment_reconstructed.toml")
        configs = Dict(path => read(path, String) for path in (count_config, unit_config))
        resume(f; dry=false)
        complete = joinpath(f.out, "features_split_complete.tsv")
        @test length(R.check_counts(complete, f.counts)) == 222
        @test length(last(R.lobe_table(joinpath(f.out, "features_split_new.tsv")))) == 204
        @test Set(readlines(joinpath(f.out, "cached_split_files.txt"))) == Set(["260215_022.sxm", "260220_083.sxm"])
        expected = TOML.parse(configs[count_config])
        expected["model"]["peak_profile"] = "split"
        expected["model"]["skew_ratio_max"] = TOML.parse(configs[unit_config])["model"]["split_skew_ratio_max"]
        @test TOML.parsefile(joinpath(f.out, "split_config.toml")) == expected
        extracts = filter(c -> c.script == "extract_lobe_features.jl", R.fake_calls)
        @test length(extracts) == 4
        @test all(c -> c.threads == 1, extracts)
        opts = [Dict(c.args[i] => c.args[i+1] for i in 1:2:length(c.args)) for c in extracts]
        @test Set(o["--chunk"] for o in opts) == Set(["1/4", "2/4", "3/4", "4/4"])
        @test all(o -> o["--selected-summary"] == f.summary &&
            o["--exclude-from"] == joinpath(f.out, "cached_split_files.txt") &&
            o["--config"] == joinpath(f.out, "split_config.toml"), opts)
        @test length(R.fake_calls) == 5
        native = only(filter(c -> c.script == "run_reconstructed_chitosan.jl", R.fake_calls))
        @test native.threads == 4
        @test Dict(native.args[i] => native.args[i+1] for i in 1:2:length(native.args)) == Dict(
            "--data-dir" => joinpath(f.out, "raw_inputs"), "--count-config" => count_config,
            "--config" => unit_config, "--selected-summary" => f.summary, "--features" => f.base,
            "--split-features" => complete, "--templates" => joinpath(f.input, "templates_cc.tsv"),
            "--outdir" => joinpath(f.out, "unknown25"))
        @test !ispath(joinpath(f.out, "unknown25"))
        @test snapshot(f.input) == before
        @test all(read(path, String) == content for (path, content) in configs)
    end
end
