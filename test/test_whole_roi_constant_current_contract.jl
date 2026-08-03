#!/usr/bin/env julia

using Test
using SHA
using TOML
using Printf

const ROOT = dirname(@__DIR__)

include(joinpath(@__DIR__, "lib", "joint_proxy", "whole_roi_real_diagnostic.jl"))
include(joinpath(@__DIR__, "build_constant_current_stm_maps.jl"))

const BUILDER_ROW_WRITE_COUNT = Ref(0)
const BUILDER_ROW_WRITE_INTERRUPT_AFTER = Ref{Union{Nothing,Int}}(nothing)

# A String-specialized test seam leaves the production AbstractString method
# unchanged. Before the set-transaction fix this interrupts immediately after a
# final artifact replacement; afterward the same interruption occurs only while
# writing transaction staging paths.
function _atomic_write_rows(path::String, header::Vector{String}, rows::Vector{String})
    result = invoke(_atomic_write_rows,
        Tuple{AbstractString,Vector{String},Vector{String}}, path, header, rows)
    BUILDER_ROW_WRITE_COUNT[] += 1
    BUILDER_ROW_WRITE_INTERRUPT_AFTER[] == BUILDER_ROW_WRITE_COUNT[] &&
        error("injected builder row-write interruption after $(BUILDER_ROW_WRITE_COUNT[]) outputs")
    return result
end

sha256_file(path::AbstractString) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function toml_string(path::AbstractString)
    return replace(path, "\\" => "\\\\")
end

function write_grid_map_and_mask(map_path::AbstractString, mask_path::AbstractString;
        half_nm::Float64=0.80, step_nm::Float64=0.08)
    coords = collect(-half_nm:step_nm:half_nm)
    open(map_path, "w") do io
        println(io, "type\tt_nm\tu_nm\tvalue")
        for typ in (0, 1), u in coords, t in coords
            value = 0.50 + 0.10typ + 0.01t - 0.02u
            println(io, join((typ, t, u, value), '\t'))
        end
    end
    open(mask_path, "w") do io
        println(io, "type\tt_nm\tu_nm\tstatus")
        for typ in (0, 1), u in coords, t in coords
            println(io, join((typ, t, u, "found"), '\t'))
        end
    end
    return nothing
end

function write_importer_map(path::AbstractString; half_nm::Float64=0.08,
        step_nm::Float64=0.08)
    coords = collect(-half_nm:step_nm:half_nm)
    open(path, "w") do io
        println(io, "type\tt_nm\tu_nm\tvalue")
        for typ in (0, 1), u in coords, t in coords
            value = 1.0 + typ + t + 2u
            println(io, join((typ, t, u, value), '\t'))
        end
    end
    return nothing
end

function write_importer_source_provenance(path::AbstractString, maps_path::AbstractString)
    open(path, "w") do io
        println(io, "schema = \"stmfit-qe-mold-provenance-v1\"")
        println(io, "observable = \"constant-current\"")
        println(io, "maps_sha256 = \"$(sha256_file(maps_path))\"")
    end
    return nothing
end

function run_importer(args::Vector{String}, dir::AbstractString)
    stdout_path = joinpath(dir, "importer.stdout")
    stderr_path = joinpath(dir, "importer.stderr")
    cmd = Cmd(`$(Base.julia_cmd()) --project=$(ROOT) $(joinpath(ROOT, "test", "import_stm_mold_maps.jl")) $(args)`;
        dir=ROOT)
    proc = run(pipeline(ignorestatus(cmd), stdout=stdout_path, stderr=stderr_path))
    return (; code=proc.exitcode, stdout=read(stdout_path, String),
        stderr=read(stderr_path, String), stdout_path, stderr_path)
end

function write_synthetic_cc_cube(path::AbstractString, typ::Int)
    nx, ny, nz = 25, 25, 44
    xstep, ystep, zstep = 0.08, 0.08, 0.05
    phase = typ == 0 ? 0.0 : pi / 3
    open(path, "w") do io
        println(io, "synthetic whole-ROI constant-current cube")
        println(io, "local test fixture only")
        println(io, "1 0.0 0.0 0.0")
        println(io, "$nx $xstep 0.0 0.0")
        println(io, "$ny 0.0 $ystep 0.0")
        println(io, "$nz 0.0 0.0 $zstep")
        println(io, "1 0.0 0.0 0.0 0.0 0.0 0.0")
        for ix in 1:nx, iy in 1:ny, iz in 1:nz
            x = (ix - 1) * xstep
            y = (iy - 1) * ystep
            z = (iz - 1) * zstep
            lateral = 1.2 + 0.15 * cos(2pi * x / 0.64 + phase) +
                0.08 * sin(2pi * y / 0.56 + phase)
            @printf(io, "%.10e\n", lateral * exp(-z / 0.35))
        end
    end
    return nothing
end

const BUILDER_HEIGHTS = (0.40, 0.45, 0.50, 0.55, 0.60)

function builder_output_paths(prefix::AbstractString;
        provenance::AbstractString=string(prefix, ".provenance.toml"))
    maps = [@sprintf("%s_h%03d.tsv", prefix, round(Int, 100height))
        for height in BUILDER_HEIGHTS]
    masks = [replace(path, r"\.tsv$" => ".mask.tsv") for path in maps]
    return (; maps, masks, provenance=String(provenance),
        all=vcat(collect(Iterators.flatten(zip(maps, masks))), [String(provenance)]))
end

function synthetic_builder_args(cube0::AbstractString, cube1::AbstractString,
        prefix::AbstractString; extra_args::Vector{String}=String[])
    args = [
        "--cube0", String(cube0), "--cube1", String(cube1), "--cube-units", "nm",
        "--origin", "0.80,0.80,0.0", "--t-axis", "1,0,0", "--u-axis", "0,1,0",
        "--half-nm", "0.08", "--step-nm", "0.08",
        "--isovalue-scan-intervals", "64", "--out-prefix", String(prefix),
    ]
    append!(args, extra_args)
    return args
end

function run_synthetic_builder(cube0::AbstractString, cube1::AbstractString,
        prefix::AbstractString; extra_args::Vector{String}=String[],
        stdout=devnull, stderr=devnull)
    args = synthetic_builder_args(cube0, cube1, prefix; extra_args)
    cmd = Cmd(`$(Base.julia_cmd()) --project=$(ROOT)
        $(joinpath(ROOT, "test", "build_constant_current_stm_maps.jl")) $(args)`;
        dir=ROOT)
    return run(pipeline(ignorestatus(cmd), stdout=stdout, stderr=stderr))
end

function run_synthetic_builder_inprocess(cube0::AbstractString, cube1::AbstractString,
        prefix::AbstractString; extra_args::Vector{String}=String[],
        failpoint=nothing, interruptpoint=nothing)
    args = synthetic_builder_args(cube0, cube1, prefix; extra_args)
    opt = _parse_cc_cli(args)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            _run_build(opt; failpoint=failpoint, interruptpoint=interruptpoint)
        end
    end
    return nothing
end

function with_builder_row_write_interruption(f::Function, after::Int)
    BUILDER_ROW_WRITE_COUNT[] = 0
    BUILDER_ROW_WRITE_INTERRUPT_AFTER[] = after
    try
        return f()
    finally
        BUILDER_ROW_WRITE_INTERRUPT_AFTER[] = nothing
        BUILDER_ROW_WRITE_COUNT[] = 0
    end
end

const BUILDER_BACKUP_PHASES = [Symbol("after_backup_$index") for index in 1:11]
const BUILDER_INSTALL_PHASES = [Symbol("after_install_$index") for index in 1:10]
const BUILDER_PUBLICATION_PHASES = vcat(
    [:after_marker], BUILDER_BACKUP_PHASES, BUILDER_INSTALL_PHASES,
    [:after_gate, :after_commit])

generation_bytes(paths) = Dict(path => read(path) for path in paths)

function restore_generation(bytes)
    for (path, content) in bytes
        ispath(path) && rm(path; force=true, recursive=true)
        write(path, content)
    end
    return nothing
end

function is_generation(paths, bytes)
    return all(path -> isfile(path) && !islink(path) && read(path) == bytes[path], paths)
end

transaction_residue(dir::AbstractString) =
    filter(name -> occursin("stmfit-txn", name), readdir(dir))

function write_temp_constant_current_config(path::AbstractString, cube0::AbstractString,
        cube1::AbstractString)
    cfg = TOML.parsefile(joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml"))
    cfg["observable"]["glcn_cube_sha256"] = sha256_file(cube0)
    cfg["observable"]["glcnac_cube_sha256"] = sha256_file(cube1)
    open(path, "w") do io
        TOML.print(io, cfg; sorted=true)
    end
    return path
end

function mold_pixel_values(path::AbstractString)
    lines = readlines(path)
    header = split(first(lines), '\t')
    pix_cols = findall(col -> occursin(r"^p\d+$", col), header)
    vals = Float64[]
    for line in lines[2:end]
        fields = split(line, '\t')
        for idx in pix_cols
            parsed = tryparse(Float64, fields[idx])
            parsed === nothing && return nothing
            push!(vals, parsed)
        end
    end
    return length(pix_cols), vals
end

function _write_fixture_type_frames(io; mode::Symbol)
    mode == :none && return nothing
    if mode == :bool_type
        println(io, "\n[[type_frames]]")
        println(io, "type = false")
        println(io, "origin_nm = [1.046771201, 0.8110361601, 0.6923691961]")
        println(io, "t_axis = [0.992443447, -0.02988831232, -0.1190071141]")
        println(io, "u_axis = [0.01607690578, 0.9931928787, -0.1153665414]")
        println(io, "\n[[type_frames]]")
        println(io, "type = 1")
        println(io, "origin_nm = [1.022828209, 0.7711483847, 0.8232730749]")
        println(io, "t_axis = [0.9709379595, 0.02480636745, 0.2380422713]")
        println(io, "u_axis = [0.07523308585, 0.9125603906, -0.4019620831]")
        return nothing
    elseif mode == :bool_frame
        println(io, "\n[[type_frames]]")
        println(io, "type = 0")
        println(io, "origin_nm = [false, 0.8110361601, 0.6923691961]")
        println(io, "t_axis = [0.992443447, -0.02988831232, -0.1190071141]")
        println(io, "u_axis = [0.01607690578, 0.9931928787, -0.1153665414]")
        println(io, "\n[[type_frames]]")
        println(io, "type = 1")
        println(io, "origin_nm = [1.022828209, 0.7711483847, 0.8232730749]")
        println(io, "t_axis = [0.9709379595, 0.02480636745, 0.2380422713]")
        println(io, "u_axis = [0.07523308585, 0.9125603906, -0.4019620831]")
        return nothing
    end
    entries = mode == :one ? (0,) : (0, 1)
    for typ in entries
        println(io, "\n[[type_frames]]")
        println(io, "type = $typ")
        if typ == 0
            println(io, "origin_nm = [1.046771201, 0.8110361601, 0.6923691961]")
            println(io, "t_axis = [0.992443447, -0.02988831232, -0.1190071141]")
            println(io, "u_axis = [0.01607690578, 0.9931928787, -0.1153665414]")
        else
            println(io, "origin_nm = [1.022828209, 0.7711483847, 0.8232730749]")
            println(io, "t_axis = [0.9709379595, 0.02480636745, 0.2380422713]")
            println(io, "u_axis = [0.07523308585, 0.9125603906, -0.4019620831]")
        end
    end
    return nothing
end

function _write_fixture_type_isovalues(io; mode::Symbol)
    mode == :none && return nothing
    if mode == :bool_type
        println(io, "\n[[type_isovalues]]")
        println(io, "type = false")
        println(io, "isovalue = 1.0e-4")
        println(io, "\n[[type_isovalues]]")
        println(io, "type = 1")
        println(io, "isovalue = 1.1e-4")
        return nothing
    elseif mode == :bool_value
        println(io, "\n[[type_isovalues]]")
        println(io, "type = 0")
        println(io, "isovalue = true")
        println(io, "\n[[type_isovalues]]")
        println(io, "type = 1")
        println(io, "isovalue = 1.1e-4")
        return nothing
    end
    entries = mode == :one ? ((0, 1.0e-4),) : ((0, 1.0e-4), (1, 1.1e-4))
    for (typ, iso) in entries
        println(io, "\n[[type_isovalues]]")
        println(io, "type = $typ")
        println(io, "isovalue = $iso")
        println(io, "source = \"calibrated\"")
    end
    return nothing
end

function write_constant_current_fixture(dir::AbstractString;
        type_frames::Symbol=:complete,
        type_isovalues::Symbol=:complete,
        bracket_hashes::Symbol=:complete,
        sidecar_metadata::Symbol=:complete,
        mold_binding::Symbol=:complete)
    map_path = joinpath(dir, "cc_h050.tsv")
    mask_path = joinpath(dir, "cc_h050.mask.tsv")
    provenance_path = joinpath(dir, "cc.provenance.toml")
    mold_path = joinpath(dir, "cc_molds.tsv")
    mold_binding_path = string(mold_path, ".provenance.toml")
    summary_path = joinpath(dir, "summary.tsv")
    file_list_path = joinpath(dir, "file-list.txt")
    data_dir = joinpath(dir, "data")
    mkpath(data_dir)

    map_half_nm = sidecar_metadata == :bad_map_extent ? 0.32 : 0.80
    map_step_nm = sidecar_metadata == :bad_map_step ? 0.16 : 0.08
    write_grid_map_and_mask(map_path, mask_path; half_nm=map_half_nm,
        step_nm=map_step_nm)
    write(mold_path, "name\ttype\tparity\tmirror\tp001\n")
    map_sha = sha256_file(map_path)
    mask_sha = sha256_file(mask_path)
    mold_sha = sha256_file(mold_path)
    height_nm = sidecar_metadata == :bad_height ? 0.55 : 0.50
    half_nm = sidecar_metadata == :bad_grid ? -0.80 : 0.80
    cube_units = sidecar_metadata == :bad_units ? "" : "bohr"
    z_spacing_nm = sidecar_metadata == :bad_z_spacing ? -0.02 : 0.02
    scan_intervals = sidecar_metadata == :bad_scan_intervals ? 0 :
        sidecar_metadata == :float_scan_intervals ? 1024.0 :
        sidecar_metadata == :bool_scan_intervals ? true : 1024
    open(provenance_path, "w") do io
        println(io, "schema = \"stmfit-qe-mold-provenance-v1\"")
        println(io, "provider = \"stm_dft_cc_diag\"")
        println(io, "observable = \"constant-current\"")
        println(io, "sample_bias_ev = -0.300")
        println(io, "height_nm = $height_nm")
        println(io, "half_nm = $half_nm")
        println(io, "step_nm = 0.08")
        println(io, "cube_units = \"$cube_units\"")
        println(io, "z_spacing_nm = $z_spacing_nm")
        sidecar_metadata == :missing_scan_intervals ||
            println(io, "isovalue_scan_intervals = $scan_intervals")
        println(io, "nominal_height_nm = 0.50")
        println(io, "bracket_heights_nm = [0.40, 0.45, 0.50, 0.55, 0.60]")
        println(io, "bracket_policy = \"diagnostics_only\"")
        println(io, "crossing_policy = \"vacuum_first_bracketed_linear\"")
        println(io, "maps_path = \"$(toml_string(map_path))\"")
        println(io, "maps_sha256 = \"$map_sha\"")
        println(io, "invalid_mask_path = \"$(toml_string(mask_path))\"")
        println(io, "invalid_mask_sha256 = \"$mask_sha\"")
        println(io, "glcn_cube_sha256 = \"80cd1d1fde94cf084cc7ea464d2bf065b36b2c015ea0bfaef8cebfee8ff88863\"")
        println(io, "glcnac_cube_sha256 = \"40649ccd9eb6768444b8ff61bf4a639b3940cf926fe3eb42254eb024faf9b5bf\"")
        _write_fixture_type_isovalues(io; mode=type_isovalues)
        _write_fixture_type_frames(io; mode=type_frames)
        for h in (0.40, 0.45, 0.55, 0.60)
            artifact_map_sha = (bracket_hashes == :stale && h == 0.40) ? repeat("a", 64) : map_sha
            println(io, "\n[[bracket_artifacts]]")
            println(io, "height_nm = $h")
            println(io, "map_path = \"$(toml_string(map_path))\"")
            println(io, "map_sha256 = \"$artifact_map_sha\"")
            println(io, "invalid_mask_path = \"$(toml_string(mask_path))\"")
            println(io, "invalid_mask_sha256 = \"$mask_sha\"")
        end
    end
    provenance_sha = sha256_file(provenance_path)
    if mold_binding != :none
        bound_map_sha = mold_binding == :stale_source_map ? repeat("b", 64) : map_sha
        bound_mold_sha = mold_binding == :stale_mold ? repeat("c", 64) : mold_sha
        parity_flip = mold_binding == :bad_parity ? "u" : "t"
        mirror_flip = mold_binding == :bad_mirror ? "t" : "u"
        normalize = mold_binding == :bad_normalize ? "none" : "zscore"
        open(mold_binding_path, "w") do io
            println(io, "schema = \"stmfit-constant-current-mold-binding-v1\"")
            println(io, "source_provenance_path = \"$(toml_string(provenance_path))\"")
            println(io, "source_provenance_sha256 = \"$provenance_sha\"")
            println(io, "source_maps_path = \"$(toml_string(map_path))\"")
            println(io, "source_maps_sha256 = \"$bound_map_sha\"")
            println(io, "molds_path = \"$(toml_string(mold_path))\"")
            println(io, "molds_sha256 = \"$bound_mold_sha\"")
            println(io, "half_nm = 0.80")
            println(io, "step_nm = 0.08")
            mold_binding == :missing_normalize || println(io, "normalize = \"$normalize\"")
            println(io, "parity_flip = \"$parity_flip\"")
            println(io, "mirror_flip = \"$mirror_flip\"")
        end
    end
    write(summary_path, "filepath\tN_selected\n240817_007.sxm\t6\n240817_050.sxm\t6\n")
    write(file_list_path, "240817_007.sxm\n240817_050.sxm\n")
    write(joinpath(data_dir, "240817_007.sxm"), "")
    write(joinpath(data_dir, "240817_050.sxm"), "")
    return (; map_path, mask_path, provenance_path, provenance_sha,
        mold_path, mold_binding_path, summary_path, file_list_path, data_dir)
end

function launcher_env(fx, outdir)
    return Dict(
        "PATH" => get(ENV, "PATH", ""),
        "HOME" => get(ENV, "HOME", ""),
        "STMFIT_SSH_HOST" => "viper.invalid",
        "STMFIT_REMOTE_USER" => "stmfit_t3",
        "STMFIT_LOCAL_DATA" => fx.data_dir,
        "STMFIT_WHOLE_ROI_CONFIG" => "config/joint_proxy_whole_roi_constant_current.toml",
        "STMFIT_WHOLE_ROI_MOLDS" => fx.mold_path,
        "STMFIT_WHOLE_ROI_SUMMARY" => fx.summary_path,
        "STMFIT_WHOLE_ROI_FILE_LIST" => fx.file_list_path,
        "STMFIT_WHOLE_ROI_LOCAL_OUT" => outdir,
        "STMFIT_CONSTANT_CURRENT_PROVENANCE" => fx.provenance_path,
        "STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256" => fx.provenance_sha,
        "STMFIT_CONSTANT_CURRENT_MAPS" => fx.map_path,
        "STMFIT_CONSTANT_CURRENT_INVALID_MASK" => fx.mask_path,
        "STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE" => fx.mold_binding_path,
        "CPUS_PER_TASK" => "4",
        "MEM_PER_CPU" => "4000",
    )
end

function run_launcher_dry_run(env::Dict{String,String}, dir::AbstractString)
    stdout_path = joinpath(dir, "launcher.stdout")
    stderr_path = joinpath(dir, "launcher.stderr")
    cmd = setenv(Cmd(`bash hpc/launch_whole_roi_frozen_contrast.sh --dry-run`;
        dir=ROOT), env)
    proc = run(pipeline(ignorestatus(cmd), stdout=stdout_path, stderr=stderr_path))
    return (; code=proc.exitcode, stdout=read(stdout_path, String),
        stderr=read(stderr_path, String), stdout_path, stderr_path)
end

function run_array_worker(env::Dict{String,String}, dir::AbstractString)
    stdout_path = joinpath(dir, "array.stdout")
    stderr_path = joinpath(dir, "array.stderr")
    cmd = setenv(Cmd(`bash hpc/whole_roi_frozen_contrast_array.sbatch`;
        dir=ROOT), env)
    proc = run(pipeline(ignorestatus(cmd), stdout=stdout_path, stderr=stderr_path))
    return (; code=proc.exitcode, stdout=read(stdout_path, String),
        stderr=read(stderr_path, String), stdout_path, stderr_path)
end

@testset "historical constant-height defaults are still the default surface" begin
    plot = read(joinpath(ROOT, "test", "plot_frozen_contrast_diagnostic.jl"), String)
    launch = read(joinpath(ROOT, "hpc", "launch_whole_roi_frozen_contrast.sh"), String)
    @test occursin("240817_007.sxm", plot)
    @test occursin("240817_050.sxm", plot)
    @test occursin("/tmp/opencode/chitosan_connected_molds_dft_m030_h050_periodic080_half080.tsv", plot)
    @test occursin("DEFAULT_DIAGNOSTIC_CONFIG_PATH", plot)
    @test occursin("240817_007.sxm", launch)
    @test occursin("240817_050.sxm", launch)
end

@testset "constant-current provenance binds map and validity mask before scoring" begin
    mktempdir() do dir
        fx = write_constant_current_fixture(dir)
        env = launcher_env(fx, joinpath(dir, "out"))
        contract = validate_real_observable_contract(
            joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
            env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
        @test contract.mode == "constant-current"
        @test contract.provider == "stm_dft_cc_diag"
        @test contract.maps_sha256 == sha256_file(fx.map_path)
        @test contract.invalid_mask_sha256 == sha256_file(fx.mask_path)
        @test contract.molds_sha256 == sha256_file(fx.mold_path)

        write(fx.mask_path, read(fx.mask_path, String) * "# tamper\n")
        @test_throws ArgumentError validate_real_observable_contract(
            joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
            env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
    end
end

@testset "whole-ROI observable validation binds the actual mold TSV before scoring" begin
    mktempdir() do dir
        fx = write_constant_current_fixture(dir)
        env = launcher_env(fx, joinpath(dir, "out"))
        @test validate_real_observable_contract(
            joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
            env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08).molds_path ==
            abspath(fx.mold_path)

        write(fx.mold_path, read(fx.mold_path, String) * "# stale\n")
        @test_throws ArgumentError validate_real_observable_contract(
            joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
            env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
    end
    mktempdir() do dir
        fx = write_constant_current_fixture(dir; mold_binding=:stale_source_map)
        env = launcher_env(fx, joinpath(dir, "out"))
        @test_throws ArgumentError validate_real_observable_contract(
            joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
            env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
    end
end

@testset "whole-ROI observable validation checks mold-binding importer convention" begin
    for mode in (:bad_parity, :bad_mirror, :bad_normalize, :missing_normalize)
        mktempdir() do dir
            fx = write_constant_current_fixture(dir; mold_binding=mode)
            env = launcher_env(fx, joinpath(dir, "out"))
            @test_throws ArgumentError validate_real_observable_contract(
                joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
                env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
        end
    end
end

@testset "whole-ROI observable validation requires complete typed frame provenance" begin
    mktempdir() do dir
        fx = write_constant_current_fixture(dir; type_frames=:one)
        env = launcher_env(fx, joinpath(dir, "out"))
        @test_throws ArgumentError validate_real_observable_contract(
            joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
            env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
    end
end

@testset "whole-ROI observable validation rejects Bool numeric provenance" begin
    for kwargs in (
            (type_frames=:bool_type,),
            (type_frames=:bool_frame,),
            (type_isovalues=:bool_type,),
            (type_isovalues=:bool_value,))
        mktempdir() do dir
            fx = write_constant_current_fixture(dir; kwargs...)
            env = launcher_env(fx, joinpath(dir, "out"))
            @test_throws ArgumentError validate_real_observable_contract(
                joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
                env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
        end
    end
end

@testset "whole-ROI observable validation checks extraction and bracket metadata" begin
    for kwargs in (
            (type_isovalues=:one,),
            (bracket_hashes=:stale,),
            (sidecar_metadata=:bad_height,),
            (sidecar_metadata=:bad_grid,),
            (sidecar_metadata=:bad_map_extent,),
            (sidecar_metadata=:bad_map_step,),
            (sidecar_metadata=:bad_units,),
            (sidecar_metadata=:bad_z_spacing,),
            (sidecar_metadata=:missing_scan_intervals,),
            (sidecar_metadata=:bad_scan_intervals,),
            (sidecar_metadata=:float_scan_intervals,),
            (sidecar_metadata=:bool_scan_intervals,))
        mktempdir() do dir
            fx = write_constant_current_fixture(dir; kwargs...)
            env = launcher_env(fx, joinpath(dir, "out"))
            @test_throws ArgumentError validate_real_observable_contract(
                joinpath(ROOT, "config", "joint_proxy_whole_roi_constant_current.toml");
                env=env, molds_path=fx.mold_path, half_nm=0.80, step_nm=0.08)
        end
    end
end

@testset "constant-current importer atomically replaces symlinked outputs" begin
    mktempdir() do dir
        maps = joinpath(dir, "maps.tsv")
        source_provenance = joinpath(dir, "maps.provenance.toml")
        molds = joinpath(dir, "molds.tsv")
        bonds = joinpath(dir, "bonds.tsv")
        binding = string(molds, ".provenance.toml")
        victim_molds = joinpath(dir, "victim-molds.tsv")
        victim_bonds = joinpath(dir, "victim-bonds.tsv")
        victim_binding = joinpath(dir, "victim-binding.toml")
        write_importer_map(maps)
        write_importer_source_provenance(source_provenance, maps)
        write(victim_molds, "sentinel molds\n")
        write(victim_bonds, "sentinel bonds\n")
        write(victim_binding, "sentinel binding\n")
        symlink(victim_molds, molds)
        symlink(victim_bonds, bonds)
        symlink(victim_binding, binding)

        result = run_importer([
            "--maps", maps,
            "--out", molds,
            "--bond-out", bonds,
            "--half-nm", "0.08",
            "--step-nm", "0.08",
            "--normalize", "none",
            "--source-provenance", source_provenance,
            "--mold-provenance-out", binding,
        ], dir)

        @test result.code == 0
        @test read(victim_molds, String) == "sentinel molds\n"
        @test read(victim_bonds, String) == "sentinel bonds\n"
        @test read(victim_binding, String) == "sentinel binding\n"
        @test !islink(molds)
        @test !islink(bonds)
        @test !islink(binding)
        @test TOML.parsefile(binding)["molds_sha256"] == sha256_file(molds)
    end
end

@testset "successful constant-current builder bytes and provenance baseline" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_baseline")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        paths = builder_output_paths(prefix)

        first_run = run_synthetic_builder(cube0, cube1, prefix)
        @test first_run.exitcode == 0
        @test all(isfile, paths.all)
        first_bytes = Dict(path => read(path) for path in paths.all)
        provenance = TOML.parsefile(paths.provenance)

        nominal_index = findfirst(==(0.50), BUILDER_HEIGHTS)
        @test provenance["maps_path"] == paths.maps[nominal_index]
        @test provenance["invalid_mask_path"] == paths.masks[nominal_index]
        @test provenance["maps_sha256"] == sha256_file(paths.maps[nominal_index])
        @test provenance["invalid_mask_sha256"] == sha256_file(paths.masks[nominal_index])
        bracket_entries = Dict(entry["height_nm"] => entry
            for entry in provenance["bracket_artifacts"])
        @test Set(keys(bracket_entries)) == Set((0.40, 0.45, 0.55, 0.60))
        for (index, height) in enumerate(BUILDER_HEIGHTS)
            height == 0.50 && continue
            entry = bracket_entries[height]
            @test entry["map_path"] == paths.maps[index]
            @test entry["invalid_mask_path"] == paths.masks[index]
            @test entry["map_sha256"] == sha256_file(paths.maps[index])
            @test entry["invalid_mask_sha256"] == sha256_file(paths.masks[index])
        end

        second_run = run_synthetic_builder(cube0, cube1, prefix)
        @test second_run.exitcode == 0
        @test Dict(path => read(path) for path in paths.all) == first_bytes
    end
end

@testset "staged provenance hashes staged bytes and records final paths" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        final_map = joinpath(dir, "final-map.tsv")
        final_mask = joinpath(dir, "final-map.mask.tsv")
        staged_map = joinpath(dir, "staged-map.tsv")
        staged_mask = joinpath(dir, "staged-map.mask.tsv")
        provenance = joinpath(dir, "staged-provenance.toml")
        write(cube0, "cube zero\n")
        write(cube1, "cube one\n")
        write(staged_map, "staged map bytes\n")
        write(staged_mask, "staged mask bytes\n")
        artifact_sources = Dict(final_map => staged_map, final_mask => staged_mask)

        write_qe_mold_provenance(provenance;
            provider=CC_PROVIDER, glcn_cube=cube0, glcnac_cube=cube1,
            maps=final_map, sample_bias_ev=-0.3, height_nm=0.5,
            half_nm=0.08, step_nm=0.08, cube_units="nm",
            observable="constant-current", nominal_height_nm=0.5,
            bracket_heights_nm=[0.5],
            type_isovalues=[
                (type=0, isovalue=1.0, source="test"),
                (type=1, isovalue=2.0, source="test"),
            ],
            z_spacing_nm=0.05,
            crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=final_mask,
            type_frames=[
                (type=0, origin_nm=[0.0, 0.0, 0.0],
                    t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0]),
                (type=1, origin_nm=[0.0, 0.0, 0.0],
                    t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0]),
            ],
            bracket_artifacts=Tuple{Float64,String,String}[],
            isovalue_scan_intervals=64,
            artifact_sources=artifact_sources)

        payload = TOML.parsefile(provenance)
        @test payload["maps_path"] == final_map
        @test payload["invalid_mask_path"] == final_mask
        @test payload["maps_sha256"] == sha256_file(staged_map)
        @test payload["invalid_mask_sha256"] == sha256_file(staged_mask)
        @test !ispath(final_map)
        @test !ispath(final_mask)
        baseline = read(provenance)

        victim = joinpath(dir, "victim.tsv")
        write(victim, "protected\n")
        rm(staged_map)
        symlink(victim, staged_map)
        @test_throws ErrorException write_qe_mold_provenance(provenance;
            provider=CC_PROVIDER, glcn_cube=cube0, glcnac_cube=cube1,
            maps=final_map, sample_bias_ev=-0.3, height_nm=0.5,
            half_nm=0.08, step_nm=0.08, cube_units="nm",
            observable="constant-current", nominal_height_nm=0.5,
            bracket_heights_nm=[0.5], isovalue=1.0,
            z_spacing_nm=0.05,
            crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=final_mask,
            frame_origin_nm=[0.0, 0.0, 0.0],
            frame_t_axis=[1.0, 0.0, 0.0],
            frame_u_axis=[0.0, 1.0, 0.0],
            bracket_artifacts=Tuple{Float64,String,String}[],
            isovalue_scan_intervals=64,
            artifact_sources=artifact_sources)
        @test read(provenance) == baseline
        @test read(victim, String) == "protected\n"
    end
end

@testset "builder output-set interruption recovery" begin
@testset "builder interruption after a partial map/mask publication preserves the old set" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_partial")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        paths = builder_output_paths(prefix)
        run_synthetic_builder_inprocess(cube0, cube1, prefix)
        old_bytes = Dict(path => read(path) for path in paths.all)

        @test_throws ErrorException with_builder_row_write_interruption(1) do
            run_synthetic_builder_inprocess(cube0, cube1, prefix;
                extra_args=["--step-nm", "0.04"])
        end
        @test all(path -> isfile(path) && read(path) == old_bytes[path], paths.all)
    end
end

@testset "builder interruption after all map/masks but before provenance preserves the old set" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_before_provenance")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        paths = builder_output_paths(prefix)
        run_synthetic_builder_inprocess(cube0, cube1, prefix)
        old_bytes = Dict(path => read(path) for path in paths.all)

        @test_throws ErrorException with_builder_row_write_interruption(10) do
            run_synthetic_builder_inprocess(cube0, cube1, prefix;
                extra_args=["--step-nm", "0.04"])
        end
        @test all(path -> isfile(path) && read(path) == old_bytes[path], paths.all)
    end
end
end


@testset "builder publication phases recover only complete old or new generations" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_phases")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        paths = builder_output_paths(prefix)

        run_synthetic_builder_inprocess(cube0, cube1, prefix)
        old_bytes = generation_bytes(paths.all)
        run_synthetic_builder_inprocess(cube0, cube1, prefix;
            extra_args=["--step-nm", "0.04"])
        new_bytes = generation_bytes(paths.all)
        @test all(path -> old_bytes[path] != new_bytes[path], paths.all)
        restore_generation(old_bytes)

        new_opt = _parse_cc_cli(synthetic_builder_args(cube0, cube1, prefix;
            extra_args=["--step-nm", "0.04"]))
        output_set = _builder_output_set(new_opt)
        gate = paths.provenance
        for phase in BUILDER_PUBLICATION_PHASES
            error_value = try
                run_synthetic_builder_inprocess(cube0, cube1, prefix;
                    extra_args=["--step-nm", "0.04"], interruptpoint=phase)
                nothing
            catch caught
                caught
            end
            @test error_value isa CubeTransactionInterruption
            if error_value isa CubeTransactionInterruption
                @test error_value.point == phase
            end
            @test isfile(_cc_txn_marker(abspath(gate)))

            if phase == :after_marker
                @test is_generation(paths.all, old_bytes)
            elseif phase in BUILDER_BACKUP_PHASES || phase in BUILDER_INSTALL_PHASES
                @test !_cc_txn_present(gate)
            else
                @test is_generation(paths.all, new_bytes)
            end

            _recover_builder_output_set(output_set)
            recovered = phase == :after_commit ? new_bytes : old_bytes
            @test is_generation(paths.all, recovered)
            @test isempty(transaction_residue(dir))
            restore_generation(old_bytes)
        end
    end
end


@testset "interrupted first publication recovers the all-absent generation" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_first_publish")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        paths = builder_output_paths(prefix)
        opt = _parse_cc_cli(synthetic_builder_args(cube0, cube1, prefix))
        output_set = _builder_output_set(opt)

        interruption = try
            run_synthetic_builder_inprocess(cube0, cube1, prefix;
                interruptpoint=:after_install_5)
            nothing
        catch caught
            caught
        end
        @test interruption isa CubeTransactionInterruption
        @test !_cc_txn_present(paths.provenance)
        marker = _cc_txn_marker(abspath(paths.provenance))
        id, _ = _cc_txn_read_marker(marker,
            _cc_txn_signature(output_set.absolute))
        absence_sidecar = _cc_txn_absent(abspath(paths.provenance), id)
        @test isfile(absence_sidecar) && !islink(absence_sidecar)
        absence_bytes = read(absence_sidecar)
        victim = joinpath(dir, "absence-victim")
        write(victim, "protected absence target\n")
        rm(absence_sidecar)
        symlink(victim, absence_sidecar)

        sidecar_error = try
            _recover_builder_output_set(output_set)
            nothing
        catch caught
            caught
        end
        @test sidecar_error isa ErrorException
        @test occursin("symlinked transaction absence sidecar",
            lowercase(sprint(showerror, sidecar_error)))
        @test read(victim, String) == "protected absence target\n"
        rm(absence_sidecar)
        write(absence_sidecar, absence_bytes)
        _recover_builder_output_set(output_set)
        @test !any(_cc_txn_present, paths.all)
        @test isempty(transaction_residue(dir))

        run_synthetic_builder_inprocess(cube0, cube1, prefix)
        @test all(path -> isfile(path) && !islink(path), paths.all)
        @test isempty(transaction_residue(dir))
    end
end


@testset "builder transaction rejects unsafe state and preserves symlink targets" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_symlinks")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        paths = builder_output_paths(prefix)
        run_synthetic_builder_inprocess(cube0, cube1, prefix)
        old_bytes = generation_bytes(paths.all)
        run_synthetic_builder_inprocess(cube0, cube1, prefix;
            extra_args=["--step-nm", "0.04"])
        new_bytes = generation_bytes(paths.all)
        restore_generation(old_bytes)
        new_opt = _parse_cc_cli(synthetic_builder_args(cube0, cube1, prefix;
            extra_args=["--step-nm", "0.04"]))
        output_set = _builder_output_set(new_opt)
        marker = _cc_txn_marker(abspath(paths.provenance))

        write(marker, "malformed marker\n")
        malformed = try
            _run_build(new_opt)
            nothing
        catch caught
            caught
        end
        @test malformed isa ErrorException
        @test occursin("malformed", lowercase(sprint(showerror, malformed)))
        @test is_generation(paths.all, old_bytes)
        rm(marker)

        marker_victim = joinpath(dir, "marker-victim")
        write(marker_victim, "protected marker target\n")
        symlink(marker_victim, marker)
        marker_error = try
            _run_build(new_opt)
            nothing
        catch caught
            caught
        end
        @test marker_error isa ErrorException
        @test occursin("symlinked transaction marker",
            lowercase(sprint(showerror, marker_error)))
        @test is_generation(paths.all, old_bytes)
        @test read(marker_victim, String) == "protected marker target\n"
        rm(marker)

        interruption = try
            run_synthetic_builder_inprocess(cube0, cube1, prefix;
                extra_args=["--step-nm", "0.04"], interruptpoint=:after_marker)
            nothing
        catch caught
            caught
        end
        @test interruption isa CubeTransactionInterruption
        id, phase = _cc_txn_read_marker(marker,
            _cc_txn_signature(output_set.absolute))
        @test phase == "prepared"
        staged_sidecar = _cc_txn_stage(first(output_set.absolute), id)
        staged_bytes = read(staged_sidecar)
        sidecar_victim = joinpath(dir, "sidecar-victim")
        write(sidecar_victim, "protected sidecar target\n")
        rm(staged_sidecar)
        symlink(sidecar_victim, staged_sidecar)
        sidecar_error = try
            _recover_builder_output_set(output_set)
            nothing
        catch caught
            caught
        end
        @test sidecar_error isa ErrorException
        @test occursin("symlinked transaction stage sidecar",
            lowercase(sprint(showerror, sidecar_error)))
        @test is_generation(paths.all, old_bytes)
        @test read(sidecar_victim, String) == "protected sidecar target\n"
        rm(staged_sidecar)
        write(staged_sidecar, staged_bytes)
        _recover_builder_output_set(output_set)
        @test is_generation(paths.all, old_bytes)
        @test isempty(transaction_residue(dir))

        map_victim = joinpath(dir, "map-victim")
        provenance_victim = joinpath(dir, "provenance-victim")
        write(map_victim, "protected map target\n")
        write(provenance_victim, "protected provenance target\n")
        rm(paths.maps[1])
        rm(paths.provenance)
        symlink(map_victim, paths.maps[1])
        symlink(provenance_victim, paths.provenance)
        destination_interruption = try
            run_synthetic_builder_inprocess(cube0, cube1, prefix;
                extra_args=["--step-nm", "0.04"], interruptpoint=:after_install_1)
            nothing
        catch caught
            caught
        end
        @test destination_interruption isa CubeTransactionInterruption
        @test !_cc_txn_present(paths.provenance)
        @test read(map_victim, String) == "protected map target\n"
        @test read(provenance_victim, String) == "protected provenance target\n"
        _recover_builder_output_set(output_set)
        @test islink(paths.maps[1])
        @test islink(paths.provenance)
        @test readlink(paths.maps[1]) == map_victim
        @test readlink(paths.provenance) == provenance_victim
        @test read(map_victim, String) == "protected map target\n"
        @test read(provenance_victim, String) == "protected provenance target\n"
        @test isempty(transaction_residue(dir))

        run_synthetic_builder_inprocess(cube0, cube1, prefix;
            extra_args=["--step-nm", "0.04"])
        @test is_generation(paths.all, new_bytes)
        @test read(map_victim, String) == "protected map target\n"
        @test read(provenance_victim, String) == "protected provenance target\n"
        @test isempty(transaction_residue(dir))
    end
end


@testset "builder failures cannot print success or leave transaction residue" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_failure")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        paths = builder_output_paths(prefix)
        run_synthetic_builder_inprocess(cube0, cube1, prefix)
        old_bytes = generation_bytes(paths.all)
        opt = _parse_cc_cli(synthetic_builder_args(cube0, cube1, prefix;
            extra_args=["--step-nm", "0.04"]))
        output_path = joinpath(dir, "failure.stdout")
        failure = open(output_path, "w") do output
            redirect_stdout(output) do
                try
                    _run_build(opt; failpoint=:after_install_5)
                    nothing
                catch caught
                    caught
                end
            end
        end
        @test failure isa ErrorException
        @test occursin("injected transaction failure", sprint(showerror, failure))
        @test !occursin("Built diagnostic constant-current maps",
            read(output_path, String))
        @test is_generation(paths.all, old_bytes)
        @test isempty(transaction_residue(dir))
    end
end


@testset "builder rejects cross-directory and colliding output sets" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_layout")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        other = mktempdir()
        try
            cross_directory = _parse_cc_cli(synthetic_builder_args(
                cube0, cube1, prefix;
                extra_args=["--provenance", joinpath(other, "provenance.toml")]))
            @test_throws ErrorException _builder_output_set(cross_directory)

            collision = _parse_cc_cli(synthetic_builder_args(cube0, cube1, prefix;
                extra_args=["--provenance", string(prefix, "_h050.tsv")]))
            @test_throws ErrorException _builder_output_set(collision)
            @test !any(_cc_txn_present, builder_output_paths(prefix).all)
        finally
            rm(other; recursive=true, force=true)
        end
    end
end

@testset "generated constant-current map imports as finite whole-ROI 21x21 molds" begin
    mktempdir() do dir
        cube0 = joinpath(dir, "glcn.cube")
        cube1 = joinpath(dir, "glcnac.cube")
        prefix = joinpath(dir, "cc_realgrid")
        write_synthetic_cc_cube(cube0, 0)
        write_synthetic_cc_cube(cube1, 1)
        run(pipeline(`$(Base.julia_cmd()) --project=$(ROOT)
            $(joinpath(ROOT, "test", "build_constant_current_stm_maps.jl"))
            --cube0 $(cube0) --cube1 $(cube1) --cube-units nm
            --origin 0.80,0.80,0.0 --t-axis 1,0,0 --u-axis 0,1,0
            --half-nm 0.80 --step-nm 0.08 --out-prefix $(prefix)`;
            stdout=devnull))

        nominal_map = string(prefix, "_h050.tsv")
        provenance = string(prefix, ".provenance.toml")
        molds = joinpath(dir, "cc_realgrid_h050_connected.tsv")
        binding = string(molds, ".provenance.toml")
        result = run_importer([
            "--maps", nominal_map,
            "--out", molds,
            "--half-nm", "0.80",
            "--step-nm", "0.08",
            "--normalize", "zscore",
            "--source-provenance", provenance,
            "--mold-provenance-out", binding,
        ], dir)
        @test result.code == 0
        pixels = mold_pixel_values(molds)
        @test pixels !== nothing
        npix, vals = pixels
        @test npix == 441
        @test length(vals) == 8 * 441
        @test all(isfinite, vals)

        config = write_temp_constant_current_config(joinpath(dir, "cc_config.toml"),
            cube0, cube1)
        env = Dict(
            "STMFIT_CONSTANT_CURRENT_PROVENANCE" => provenance,
            "STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256" => sha256_file(provenance),
            "STMFIT_CONSTANT_CURRENT_MAPS" => nominal_map,
            "STMFIT_CONSTANT_CURRENT_INVALID_MASK" => string(prefix, "_h050.mask.tsv"),
            "STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE" => binding,
        )
        contract = validate_real_observable_contract(config;
            env=env, molds_path=molds, half_nm=0.80, step_nm=0.08)
        @test contract.molds_sha256 == sha256_file(molds)
        @test contract.maps_sha256 == sha256_file(nominal_map)
    end
end

@testset "forbidden truth and manifest inputs are rejected at runtime boundaries" begin
    @test_throws ArgumentError reject_forbidden_real_runtime_inputs(
        ["--truth", "benchmarks/truth.tsv"]; boundary="CLI")
    @test_throws ArgumentError reject_forbidden_real_runtime_inputs(
        ["--manifest=benchmarks/chitosan_6mer_counting_confirmed.toml"];
        boundary="CLI")
    mktempdir() do dir
        fx = write_constant_current_fixture(dir)
        env = launcher_env(fx, joinpath(dir, "out"))
        env["STMFIT_WHOLE_ROI_TRUTH"] = joinpath(dir, "truth.tsv")
        result = run_launcher_dry_run(env, dir)
        @test result.code != 0
        @test occursin("forbidden", lowercase(result.stderr * result.stdout))
        @test occursin("truth", lowercase(result.stderr * result.stdout))
    end
end

@testset "launcher dry-run propagates explicit no-truth inputs" begin
    mktempdir() do dir
        fx = write_constant_current_fixture(dir)
        outdir = joinpath(dir, "out")
        env = launcher_env(fx, outdir)
        result = run_launcher_dry_run(env, dir)
        @test result.code == 0
        text = result.stdout * result.stderr
        @test occursin("config:  config/joint_proxy_whole_roi_constant_current.toml", text)
        @test occursin("molds:   $(fx.mold_path)", text)
        @test occursin("file-list: $(fx.file_list_path)", text)
        @test occursin("output:  ", text)
        @test occursin(outdir, text)
        @test occursin("array size: 2", text)
        @test occursin("max concurrency: 2", text)
        @test occursin("cpus per task: 4", text)
        @test occursin("julia threads/task: 4", text)
        @test occursin("--array=1-2%2", text)
        @test occursin("STMFIT_WHOLE_ROI_FILE_LIST=", text)
        @test !occursin("scientific success", lowercase(text))
    end
end

@testset "launcher rejects quote and command-substitution path inputs before SSH" begin
    mktempdir() do dir
        fx = write_constant_current_fixture(dir)
        outdir = joinpath(dir, "out")
        env = launcher_env(fx, outdir)
        malicious_mold = joinpath(dir, "cc_molds';touch injected.tsv")
        write(malicious_mold, "x\n")
        env["STMFIT_WHOLE_ROI_MOLDS"] = malicious_mold
        env["STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE"] = string(malicious_mold, ".provenance.toml")
        result = run_launcher_dry_run(env, dir)
        text = result.stdout * result.stderr
        @test result.code != 0
        @test occursin("unsafe", lowercase(text))
        @test !occursin("[dry-run] ssh", text)
        @test !isfile(joinpath(dir, "injected.tsv"))

        env = launcher_env(fx, outdir)
        env["STMFIT_WHOLE_ROI_REMOTE_OUT"] = "/ptmp/stmfit/\$(touch injected)"
        result = run_launcher_dry_run(env, dir)
        text = result.stdout * result.stderr
        @test result.code != 0
        @test occursin("unsafe", lowercase(text))
        @test !occursin("[dry-run] ssh", text)
    end
end

@testset "launcher and array reject unsafe SXM file-list basenames before transport" begin
    mktempdir() do dir
        fx = write_constant_current_fixture(dir)
        outdir = joinpath(dir, "out")
        write(fx.file_list_path, "bad';touch injected.sxm\n")

        env = launcher_env(fx, outdir)
        result = run_launcher_dry_run(env, dir)
        text = result.stdout * result.stderr
        @test result.code != 0
        @test occursin("unsafe", lowercase(text))
        @test occursin("sxm", lowercase(text))
        @test !occursin("[dry-run] ssh", text)
        @test !isfile(joinpath(dir, "injected.sxm"))

        array_env = Dict(
            "PATH" => get(ENV, "PATH", ""),
            "HOME" => get(ENV, "HOME", ""),
            "SLURM_ARRAY_TASK_ID" => "1",
            "STMFIT_REMOTE_PROJECT" => ROOT,
            "STMFIT_DATA_DIR" => fx.data_dir,
            "STMFIT_WHOLE_ROI_MOLDS" => fx.mold_path,
            "STMFIT_WHOLE_ROI_SUMMARY" => fx.summary_path,
            "STMFIT_WHOLE_ROI_CONFIG" => joinpath(ROOT, "config",
                "joint_proxy_whole_roi_constant_current.toml"),
            "STMFIT_WHOLE_ROI_FILE_LIST" => fx.file_list_path,
            "STMFIT_WHOLE_ROI_OUTDIR" => outdir,
            "STMFIT_CONSTANT_CURRENT_PROVENANCE" => fx.provenance_path,
            "STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256" => fx.provenance_sha,
            "STMFIT_CONSTANT_CURRENT_MAPS" => fx.map_path,
            "STMFIT_CONSTANT_CURRENT_INVALID_MASK" => fx.mask_path,
            "STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE" => fx.mold_binding_path,
        )
        result = run_array_worker(array_env, dir)
        text = result.stdout * result.stderr
        @test result.code != 0
        @test occursin("unsafe", lowercase(text))
        @test occursin("sxm", lowercase(text))
        @test !occursin("module", lowercase(text))
        @test !isfile(joinpath(dir, "injected.sxm"))
    end
end
