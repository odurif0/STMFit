using Test, Random, Statistics, Printf
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
const ASSIGNMENT_CONFIG = joinpath(ROOT, "config", "unit_assignment_reconstructed.toml")
const COUNT_CONFIG = joinpath(ROOT, "config", "chitosan.toml")

@testset "Matched-residual pipeline must regenerate patches" begin
    mktempdir() do dir
        for name in ("matched_residual", "transverse_moment", "affine_residual", "signed_mold", "affine_fisher", "centered_fisher", "shrunk_gmm"),
            key in ("--patches-fwd","--patches-bwd","--descriptor-patches")
            config = joinpath(ROOT,"config","unit_assignment_" * name * ".toml")
            opts = Dict("--config"=>config,"--outdir"=>joinpath(dir,"not_created"),key=>"cached.tsv")
            @test_throws ErrorException execute_pipeline(opts)
            @test !ispath(opts["--outdir"])
        end
    end
end

function pipeline_fixture(dir)
    raw = joinpath(dir, "raw"); mkdir(raw)
    rng = MersenneTwister(20260916)
    geometry = Dict{String,String}[]
    splitrows = Dict{String,String}[]
    fwdrows = Dict{String,String}[]
    bwdrows = Dict{String,String}[]
    descrows = Dict{String,String}[]
    coords17 = collect(-8:8) ./ 8
    coords9 = collect(-4:4) ./ 4
    for fi in 1:5
        file = "synthetic_$(fi).sxm"
        write(joinpath(raw, file), "Fixture only; all extracted artifacts are supplied.\n")
        for lobe in 1:8
            a = 0.05 + 0.05rand(rng)
            row = Dict("file" => file, "lobe" => string(lobe), "N" => "8", "amplitude" => string(a),
                "x_nm" => string(0.5lobe), "y_nm" => string(0.1sin(lobe)),
                "t_nm" => string(0.5lobe), "u_nm" => string(0.1sin(lobe)),
                "sigma_parallel_nm" => string(0.2 + 0.02rand(rng)), "sigma_perp_nm" => "0.10",
                "amp_rel" => string(a / 0.07), "skew_ratio" => "1.0")
            push!(geometry, row)
            push!(splitrows, Dict("file" => file, "lobe" => string(lobe), "skew_ratio" => string(0.8 + 0.5rand(rng))))
            for (prefix, coords, dest) in (("res", coords17, fwdrows), ("bwd_res", coords17, bwdrows), ("bwd_res", coords9, descrows))
                vals = [a * exp(-2(t*t + u*u)) + .02randn(rng) + .01sin(lobe)*u for u in coords for t in coords]
                vals = (vals .- median(vals)) ./ std(vals)
                patch = Dict("file" => file, "lobe" => string(lobe), "amplitude" => string(a))
                for (i, v) in enumerate(vals)
                    patch[@sprintf("%s_p%03d", prefix, i)] = @sprintf("%.7g", v)
                end
                push!(dest, patch)
            end
        end
    end
    paths = Dict{String,String}()
    for (key, name, rows) in (("--features", "features", geometry), ("--split-features", "split", splitrows),
        ("--patches-fwd", "fwd", fwdrows), ("--patches-bwd", "bwd", bwdrows), ("--descriptor-patches", "desc", descrows))
        path = joinpath(dir, name * ".tsv")
        header = vcat(["file", "lobe"], sort(setdiff(collect(keys(first(rows))), ["file", "lobe"])))
        write_table(path, header, rows)
        paths[key] = path
    end
    templates = Dict{String,String}[]
    for typ in 0:1, parity in 0:1, mirror in 0:1
        row = Dict("name" => "fixture_$(typ)_$(parity)_$(mirror)", "type" => string(typ), "parity" => string(parity), "mirror" => string(mirror))
        values = [(typ == 0 ? -1 : 1) * (t + .3u + exp(-3(t*t+u*u))) for u in coords17 for t in coords17]
        for (i, v) in enumerate(values)
            row[@sprintf("p%03d", i)] = string(v)
        end
        push!(templates, row)
    end
    paths["--templates"] = joinpath(dir, "templates.tsv")
    write_table(paths["--templates"], vcat(["name", "type", "parity", "mirror"], [@sprintf("p%03d", i) for i in 1:289]), templates)
    paths["--data-dir"] = raw
    paths["--count-config"] = COUNT_CONFIG
    paths["--config"] = ASSIGNMENT_CONFIG
    paths["--outdir"] = joinpath(dir, "output")
    return paths
end

@testset "Native pipeline boundaries" begin
    @test_throws ErrorException parse_options(["--truth", "forbidden.tsv"])
    @test_throws ErrorException parse_options(["--expected-N", "8"])
    @test_throws ErrorException parse_options(["--reference", "champion.tsv"])
    mktempdir() do dir
        opts = pipeline_fixture(dir)
        parsed = parse_options(vcat([[k, v] for (k, v) in opts]...))
        @test parsed["--data-dir"] == opts["--data-dir"]
        opts["--dry-run"] = "true"
        execute_pipeline(opts)
        @test !ispath(opts["--outdir"])
        counts = Dict("synthetic_$i.sxm" => 8 for i in 1:5)
        @test length(check_counts(opts["--features"], counts)) == 40
        counts["synthetic_1.sxm"] = 7
        @test_throws ErrorException check_counts(opts["--features"], counts)
        bad = joinpath(dir, "empty_summary.tsv")
        write(bad, "filepath\tN_selected\n")
        @test_throws ErrorException selected_counts(bad)
        delete!(opts, "--dry-run")
        rm(joinpath(opts["--data-dir"], "synthetic_5.sxm"))
        @test_throws ErrorException execute_pipeline(opts)
        _, failed = read_table(joinpath(opts["--outdir"], "failures.tsv"))
        @test length(failed) == 5
        @test all(row["stage"] == "input" && row["status"] == "incomplete" for row in failed)
        @test !isfile(joinpath(opts["--outdir"], "predictions.tsv"))
    end
end

if "--e2e" in ARGS
    @testset "Native extracted-input pipeline end to end" begin
        mktempdir() do dir
            opts = pipeline_fixture(dir)
            try
                execute_pipeline(opts)
            catch
                for (root, _, files) in walkdir(opts["--outdir"]), file in files
                    endswith(file, ".log") || continue
                    log = read(joinpath(root, file), String)
                    println("\n--- ", file, " ---\n", last(log, 4000))
                end
                rethrow()
            end
            out = opts["--outdir"]
            _, pred = lobe_table(joinpath(out, "predictions.tsv"))
            @test length(pred) == 40
            @test all(r["predicted"] in ("0", "1", "?") for r in values(pred))
            @test any(r["predicted"] != "?" for r in values(pred))
            @test all(r["model"] == "cc_soft_reconstructed_v1" for r in values(pred))
            @test isfile(joinpath(out, "summary.tsv"))
            @test isfile(joinpath(out, "review_queue.tsv"))
            @test filesize(joinpath(out, "plots", "summary_grid.png")) > 0
            @test length(readdir(joinpath(out, "plots", "standalone"))) == 5
            @test !isfile(joinpath(out, "failures.tsv"))
        end
    end
end

@testset "Bounded feature-export chunks and exact merge" begin
    mktempdir() do dir
        calls = Channel{Tuple{String,Int}}(4)
        function fake_runner(outdir, name, script, args; threads)
            ci = findfirst(==("--chunk"), args)
            chunk = ci === nothing ? "1/1" : args[ci+1]
            put!(calls, (chunk, threads))
            oi = findfirst(==("--out"), args)
            file = "chunk_$(first(split(chunk, '/'))).sxm"
            rows = [Dict("file"=>file, "lobe"=>string(i), "amplitude"=>"1.0") for i in 1:2]
            write_table(args[oi+1], ["file", "lobe", "amplitude"], rows)
        end
        out = joinpath(dir, "merged.tsv")
        export_features(dir, "base", "extract_lobe_features.jl", String[], out;
                        nfiles=5, thread_budget=2, runner=fake_runner)
        _, merged = lobe_table(out)
        @test length(merged) == 4
        @test Set([take!(calls), take!(calls)]) == Set([("1/2", 1), ("2/2", 1)])
        @test_throws ErrorException merge_feature_chunks([out, out], joinpath(dir, "duplicate.tsv"))
        @test !isfile(joinpath(dir, "duplicate.tsv"))
        single = joinpath(dir, "single.tsv")
        export_features(dir, "single", "extract_lobe_features.jl", String[], single;
                        nfiles=1, thread_budget=4, runner=fake_runner)
        @test take!(calls) == ("1/1", 4)
        @test length(last(lobe_table(single))) == 2
        @test_throws ErrorException export_features(dir, "bad", "x", String[], "bad";
                                                    nfiles=0, runner=fake_runner)
    end
end
