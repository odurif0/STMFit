using Test
using SHA
using TOML

include(joinpath(@__DIR__, "AuthorityIdentity.jl"))
using .AuthorityIdentity

function _fixture()
    dir = mktempdir()
    files = Dict{String,Vector{UInt8}}()
    function add(name, bytes)
        p = joinpath(dir, name)
        write(p, bytes)
        files[realpath(p)] = Vector{UInt8}(bytes)
        return realpath(p)
    end
    exe = add("julia", codeunits("fixture executable"))
    sys = add("sys.so", codeunits("fixture sysimage"))
    project = add("Project.toml", codeunits("name = \"fixture\"\n"))
    manifest = add("Manifest.toml", codeunits("manifest_format = \"2.0\"\n"))
    common = DSORecord[]
    versioned = DSORecord[]
    loaded = String[]
    for (name, bytes) in (("libcommon-a.so", codeunits("a")), ("libcommon-b.so", codeunits("b")))
        p = add(name, bytes); push!(loaded, p)
        b = Dict("path" => p, "size" => Int64(length(bytes)), "sha256" => bytes2hex(sha256(bytes)))
        push!(common, DSORecord(name=name, path=p, size=b["size"], sha256=b["sha256"]))
    end
    p = add("libjulia.so", codeunits("julia component")); push!(loaded, p)
    push!(versioned, DSORecord(name=basename(p), path=p, size=length(files[p]),
        sha256=bytes2hex(sha256(files[p])), runtime_label="julia-$(VERSION)",
        runtime_version=string(VERSION)))
    contract = DSOContract(runtime_label="julia-$(VERSION)", runtime_version=string(VERSION),
                           common=common, versioned=versioned)
    env = Dict{String,String}(
        "JULIA_DEPOT_PATH" => dir, "JULIA_LOAD_PATH" => dir,
        "JULIA_PROJECT" => project, "JULIA_HISTORY" => joinpath(dir, "history"),
        "JULIA_STARTUP_FILE" => "default", "JULIA_NUM_THREADS" => "1",
        "TMPDIR" => dir, "JULIA_COMPILED_MODULES" => "1")
    function info(path)
        p = realpath(path); b = files[p]
        Dict{String,Any}("path" => p, "size" => Int64(length(b)),
                         "sha256" => bytes2hex(sha256(b)))
    end
    function probe(; kwargs...)
        vals = Dict{Symbol,Any}(kwargs)
        getv(k, default) = get(vals, k, default)
        return AuthorityIdentity._SystemProbe(
            version=()->VERSION, julia_command=()->exe, sys_bindir=()->dir,
            running_executable=()->exe, canonical_path=p->realpath(p), file_info=info,
            force_numerics=()->nothing, dllist=()->copy(getv(:loaded, loaded)),
            julia_threads=()->getv(:julia_threads, 1), blas_threads=()->getv(:blas_threads, 1),
            sysimage=()->sys, project=()->project, manifest=()->manifest,
            environment=()->copy(getv(:environment, env)), hostname=()->getv(:hostname, "fixture-node"),
            kernel=()->getv(:kernel, "Linux fixture"), cpu=()->getv(:cpu, "fixture-cpu"),
            microcode=()->getv(:microcode, "fixture-microcode"),
            lscpu=()->getv(:lscpu, Vector{UInt8}(codeunits("fixture lscpu"))),
            slurm=()->getv(:slurm, Dict("job_id"=>"job-1", "node"=>"fixture-node",
                "task_id"=>"task-1", "array_job_id"=>"array-1", "array_task_id"=>"3")))
    end
    return (; dir, files, exe, sys, project, manifest, common, versioned, loaded, contract, env, info, probe)
end

@testset "AuthorityIdentity-v2" begin
    @testset "current runtime identity" begin
        @test validate_runtime_identity(runtime_label="julia-$(VERSION)", runtime_version=string(VERSION))
        @test_throws IdentityError validate_runtime_identity(runtime_label="julia-wrong", runtime_version=string(VERSION))
        @test_throws IdentityError validate_runtime_identity(runtime_label="julia-$(VERSION)", runtime_version="0.0.0")
    end

    f = _fixture()
    p = f.probe()
    report = AuthorityIdentity._capture_with_probe(p; dso_contract=f.contract,
        runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
    @test validate_identity(report)
    summary = identity_summary(report)
    @test summary["dso_common_count"] == 2
    @test summary["dso_versioned_count"] == 1
    @test summary["dso_actual_count"] == 3
    @test length(summary["sysimage_sha256"]) == 64

    @testset "strict DSO and process rejection" begin
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe();
            dso_contract=f.contract, runtime_label="julia-wrong", runtime_version=string(VERSION), executable=f.exe)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe();
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.sys)
        altered = DSOContract(runtime_label="julia-$(VERSION)", runtime_version=string(VERSION),
            common=[DSORecord(name=r.name, path=r.path, size=r.size, sha256=repeat("0", 64)) for r in f.common],
            versioned=f.versioned)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe();
            dso_contract=altered, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(loaded=f.loaded[1:2]);
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        extra = copy(f.loaded); push!(extra, f.sys)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(loaded=extra);
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(julia_threads=2);
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(blas_threads=2);
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        badenv = copy(f.env); badenv["JULIA_NUM_THREADS"] = "UNAVAILABLE"
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(environment=badenv);
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
    end

    @testset "system and Slurm rejection" begin
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(lscpu=nothing);
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(microcode="UNAVAILABLE");
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(hostname="other-node");
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        sl = Dict("job_id"=>"job-1", "node"=>"other-node", "task_id"=>"task-1",
                  "array_job_id"=>"array-1", "array_task_id"=>"3")
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(f.probe(slurm=sl);
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe)
        wrong = f.probe()
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(wrong;
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe,
            sysimage_path=f.project)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(wrong;
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe,
            project_path=f.manifest)
        @test_throws IdentityError AuthorityIdentity._capture_with_probe(wrong;
            dso_contract=f.contract, runtime_label="julia-$(VERSION)", runtime_version=string(VERSION), executable=f.exe,
            manifest_path=f.project)
    end

    @testset "strict contract serialization" begin
        json = joinpath(f.dir, "contract.json")
        write_dso_contract(json, f.contract)
        parsed = parse_dso_contract(json)
        @test length(parsed.common) == 2
        @test length(parsed.versioned) == 1
        @test_throws IdentityError parse_dso_contract(joinpath(f.dir, "missing.json"))
    end

    @testset "receipt-last publication" begin
        root = joinpath(f.dir, "identity-root")
        receipt = joinpath(f.dir, "identity.receipt.toml")
        publication = publish_identity(report, root; receipt=receipt)
        @test isfile(joinpath(root, "identity.toml"))
        @test isfile(joinpath(root, "manifest.toml"))
        @test isfile(receipt)
        @test TOML.parsefile(receipt)["schema"] == AuthorityIdentity.RECEIPT_SCHEMA
        filemode = stat(joinpath(root, "identity.toml")).mode & 0o7777
        @test filemode == 0o644
        @test_throws IdentityError publish_identity(report, root; receipt=joinpath(f.dir, "collision.receipt"))
        root2 = joinpath(f.dir, "identity-root-2")
        @test_throws IdentityError publish_identity(report, root2; receipt=receipt)
        @test isdir(root2)
        @test !isfile(joinpath(root2, "receipt.toml"))
        root3 = joinpath(f.dir, "identity-root-3")
        @test_throws IdentityError publish_identity(report, root3;
            receipt=joinpath(f.dir, "fault.receipt"),
            _fault=point -> point == :before_root_rename ? ErrorException("fault") : nothing)
        @test !ispath(root3)
        @test publication["file_count"] == 2
    end
end
