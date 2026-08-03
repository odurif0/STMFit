#!/usr/bin/env julia

# Focused real fused-fitmask common-registration/frozen-contrast diagnostic.
# Diagnostic only: does not alter fitting, selection, calibration, thresholds,
# mold registry, or production unit assignment.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using Printf
using SHA
using STMSXMIO
using GaussianFit2D

include(joinpath(@__DIR__, "diagnose_joint_proxy_whole_roi_batch.jl"))
include(joinpath(@__DIR__, "lib", "joint_proxy", "whole_roi_real_diagnostic.jl"))
include(joinpath(@__DIR__, "lib", "joint_proxy", "whole_roi_real_plot.jl"))

const REAL_DEFAULT_FILES = ["240817_007.sxm", "240817_050.sxm"]
const REAL_DEFAULT_SUMMARY = joinpath(@__DIR__, "..", "results",
    "best_plots_240817_primary_rerun", "summary_overlap060_hard.tsv")
const REAL_DEFAULT_MOLDS = "/tmp/opencode/chitosan_connected_molds_dft_m030_h050_periodic080_half080.tsv"
const REAL_DEFAULT_CONFIG = Main.WholeRoiDiagnosticConfig.DEFAULT_DIAGNOSTIC_CONFIG_PATH

function _real_usage()
    return """Usage: julia --project=. test/plot_frozen_contrast_diagnostic.jl [options] [files...]
  --out-dir DIR
  --selection-summary PATH
  --molds PATH
  --config PATH
  --file-list PATH
  --observable-provenance PATH
  --observable-provenance-sha256 SHA256
  --observable-map PATH
  --validity-mask PATH
  --observable-mold-provenance PATH
  --half-nm FLOAT
  --step-nm FLOAT
Defaults to 240817_007.sxm and 240817_050.sxm."""
end

function _parse_real_args(args)
    reject_forbidden_real_runtime_inputs(args; boundary="real diagnostic CLI")
    outdir = "/tmp/opencode/frozen_contrast_diagnostic"
    summary, molds, config = REAL_DEFAULT_SUMMARY, REAL_DEFAULT_MOLDS, REAL_DEFAULT_CONFIG
    half_nm, step_nm = 0.80, 0.08
    file_list::Union{Nothing,String} = nothing
    observable_provenance::Union{Nothing,String} = nothing
    observable_provenance_sha256::Union{Nothing,String} = nothing
    observable_map::Union{Nothing,String} = nothing
    validity_mask::Union{Nothing,String} = nothing
    observable_mold_provenance::Union{Nothing,String} = nothing
    files = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        arg in ("-h", "--help") && return (; help=true)
        if arg in ("--out-dir", "--selection-summary", "--molds", "--config",
                   "--file-list", "--observable-provenance",
                   "--observable-provenance-sha256", "--observable-map",
                   "--validity-mask", "--observable-mold-provenance",
                   "--half-nm", "--step-nm")
            i < length(args) || throw(ArgumentError("$arg requires a value"))
            value = args[i+1]
            arg == "--out-dir" && (outdir = value)
            arg == "--selection-summary" && (summary = value)
            arg == "--molds" && (molds = value)
            arg == "--config" && (config = value)
            arg == "--file-list" && (file_list = value)
            arg == "--observable-provenance" && (observable_provenance = value)
            arg == "--observable-provenance-sha256" &&
                (observable_provenance_sha256 = value)
            arg == "--observable-map" && (observable_map = value)
            arg == "--validity-mask" && (validity_mask = value)
            arg == "--observable-mold-provenance" && (observable_mold_provenance = value)
            if arg == "--half-nm"
                parsed = tryparse(Float64, value)
                parsed === nothing && throw(ArgumentError("invalid --half-nm: $value"))
                half_nm = parsed
            elseif arg == "--step-nm"
                parsed = tryparse(Float64, value)
                parsed === nothing && throw(ArgumentError("invalid --step-nm: $value"))
                step_nm = parsed
            end
            i += 2
        elseif startswith(arg, "-")
            throw(ArgumentError("unknown option: $arg"))
        else
            push!(files, String(arg)); i += 1
        end
    end
    half_nm > 0 || throw(ArgumentError("--half-nm must be positive"))
    step_nm > 0 || throw(ArgumentError("--step-nm must be positive"))
    if file_list !== nothing && !isempty(files)
        throw(ArgumentError("--file-list cannot be combined with positional files"))
    end
    selected_files = file_list === nothing ?
        (isempty(files) ? REAL_DEFAULT_FILES : files) : read_real_file_list(file_list)
    return (; help=false, outdir=abspath(outdir), summary=abspath(summary),
        molds=abspath(molds), config=abspath(config), half_nm, step_nm,
        files=selected_files, file_list=file_list === nothing ? nothing : abspath(file_list),
        observable_provenance=observable_provenance === nothing ? nothing :
            abspath(observable_provenance),
        observable_provenance_sha256,
        observable_map=observable_map === nothing ? nothing : abspath(observable_map),
        validity_mask=validity_mask === nothing ? nothing : abspath(validity_mask),
        observable_mold_provenance=observable_mold_provenance === nothing ?
            nothing : abspath(observable_mold_provenance))
end

function _load_fused_fitmask_case(filename, selected_n, data_dir, outdir)
    path = joinpath(data_dir, filename)
    isfile(path) || throw(ArgumentError("SXM file not found: $path"))
    image = STMSXMIO.read_sxm(path)
    pattern = Inference._default_pattern_config(path, outdir)
    chain = Inference._default_chain_config()
    results, _, context = redirect_stdout(devnull) do
        GaussianFit2D.chain_gaussian_sweep(image, pattern, chain)
    end
    views = Inference._build_views(image, pattern,
        Main.JointProxyCandidateViews.build_view_data)
    report = Main.JointProxyCandidateViews.extract_candidate_views(
        results, context, chain, views; patch_half_nm=0.32, patch_step_nm=0.08)
    candidate_index = findfirst(item -> item.n == selected_n, report.candidates)
    candidate_index === nothing && throw(ArgumentError(
        "selected N=$selected_n missing for $filename"))
    candidate = report.candidates[candidate_index]
    fit_index = findfirst(r -> r.n == selected_n && r.success && r.valid, results)
    fit_index === nothing && throw(ArgumentError("valid selected fit missing for $filename"))
    fit = results[fit_index]
    centers = [(Float64(lobe.x_nm), Float64(lobe.y_nm)) for lobe in candidate.lobes]
    margin = 0.90
    ix = findall(x -> minimum(first, centers)-margin <= x <= maximum(first, centers)+margin,
                 context.xs)
    iy = findall(y -> minimum(last, centers)-margin <= y <= maximum(last, centers)+margin,
                 context.ys)
    xs, ys = Float64.(context.xs[ix]), Float64.(context.ys[iy])
    fused = Matrix{Float64}(context.zimg[iy, ix])
    xgrid = repeat(reshape(xs, 1, :), length(ys), 1)
    ygrid = repeat(reshape(ys, :, 1), 1, length(xs))
    model = GaussianFit2D._chain_model_values(vec(xgrid), vec(ygrid),
        fit.params, selected_n, context.axisctx, chain;
        amp_min=fit.amp_min, amp_range=fit.amp_range)
    backbone = reshape(model, length(ys), length(xs))
    t, u = GaussianFit2D._chain_coordinates(xgrid, ygrid, context.axisctx)
    fit_mask = isfinite.(fused) .& isfinite.(backbone) .&
        (abs.(u) .<= context.fit_width_nm) .&
        (t .>= context.axisctx.tmin) .& (t .<= context.axisctx.tmax)
    count(fit_mask) >= 20 || throw(ArgumentError(
        "fused fit mask has fewer than 20 finite pixels for $filename"))
    rows = findall(vec(any(fit_mask; dims=2)))
    cols = findall(vec(any(fit_mask; dims=1)))
    row_range, col_range = first(rows):last(rows), first(cols):last(cols)
    observed = fused[row_range, col_range]
    fit_mask = fit_mask[row_range, col_range]
    observed[.!fit_mask] .= NaN
    theta = atan(context.axisctx.axis[2], context.axisctx.axis[1])
    return WholeRoiObservedCase(observed, backbone[row_range, col_range], centers,
        xs[col_range], ys[row_range], theta)
end

_sequence_text(sequence) = join(sequence)
_float(value) = @sprintf("%.12g", value)
_registration_sha(bundle) = bytes2hex(sha256(bundle.registration_tsv))

struct SimulatedTransactionInterruption <: Exception
    point::Symbol
end

Base.showerror(io::IO, err::SimulatedTransactionInterruption) =
    print(io, "simulated transaction interruption at ", err.point)

_txn_present(path::AbstractString) = ispath(path) || islink(path)
_txn_marker(gate::AbstractString) = string(gate, ".stmfit-txn")
_txn_signature(destinations) = bytes2hex(sha256(join(abspath.(destinations), '\0')))
_txn_stage(path, id) = string(path, ".", id, ".stmfit-txn.new", splitext(path)[2])
_txn_backup(path, id) = string(path, ".", id, ".stmfit-txn.old")
_txn_absent(path, id) = string(path, ".", id, ".stmfit-txn.absent")

function _txn_rename(source::AbstractString, destination::AbstractString)
    rc = ccall(:rename, Cint, (Cstring, Cstring), source, destination)
    rc == 0 || Base.systemerror("rename $source to $destination", true)
    return nothing
end

function _txn_fsync_file(path::AbstractString)
    open(path, "r") do io
        rc = ccall(:fsync, Cint, (Cint,), fd(io))
        rc == 0 || Base.systemerror("fsync $path", true)
    end
    return nothing
end

function _txn_fsync_parent(path::AbstractString)
    directory = dirname(abspath(path))
    directory_fd = ccall(:open, Cint, (Cstring, Cint), directory, 0)
    directory_fd >= 0 || Base.systemerror("open directory $directory", true)
    try
        rc = ccall(:fsync, Cint, (Cint,), directory_fd)
        rc == 0 || Base.systemerror("fsync directory $directory", true)
    finally
        ccall(:close, Cint, (Cint,), directory_fd)
    end
    return nothing
end

function _txn_write_marker(marker, id, phase, signature)
    temp_path, io = mktemp(dirname(marker); cleanup=false)
    try
        println(io, "stmfit-output-transaction-v1")
        println(io, "id\t", id)
        println(io, "phase\t", phase)
        println(io, "signature\t", signature)
        flush(io)
        rc = ccall(:fsync, Cint, (Cint,), fd(io))
        rc == 0 || Base.systemerror("fsync transaction marker", true)
        close(io)
        _txn_rename(temp_path, marker)
        _txn_fsync_parent(marker)
    finally
        isopen(io) && close(io)
        rm(temp_path; force=true)
    end
    return nothing
end

function _txn_read_marker(marker, signature)
    islink(marker) && error("refusing symlinked transaction marker: $marker")
    lines = readlines(marker)
    length(lines) == 4 && lines[1] == "stmfit-output-transaction-v1" ||
        error("malformed output transaction marker: $marker")
    fields = Dict{String,String}()
    for line in lines[2:end]
        parts = split(line, '\t'; limit=2)
        length(parts) == 2 || error("malformed output transaction marker: $marker")
        fields[parts[1]] = parts[2]
    end
    get(fields, "signature", "") == signature ||
        error("output transaction destination set does not match marker: $marker")
    phase = get(fields, "phase", "")
    phase in ("prepared", "committed") ||
        error("invalid output transaction phase in $marker")
    id = get(fields, "id", "")
    occursin(r"^[0-9a-f]{16}$", id) || error("invalid output transaction id in $marker")
    return id, phase
end

function _txn_cleanup(destinations, gate, id; remove_marker=true)
    for destination in destinations
        rm(_txn_stage(destination, id); force=true)
        rm(_txn_backup(destination, id); force=true)
        rm(_txn_absent(destination, id); force=true)
    end
    remove_marker && rm(_txn_marker(gate); force=true)
    _txn_fsync_parent(gate)
    return nothing
end

function _txn_restore_one(destination, id)
    backup = _txn_backup(destination, id)
    absent = _txn_absent(destination, id)
    if _txn_present(backup)
        _txn_present(destination) && rm(destination; force=true)
        _txn_rename(backup, destination)
    elseif isfile(absent)
        _txn_present(destination) && rm(destination; force=true)
        rm(absent; force=true)
    end
    return nothing
end

function _txn_rollback(destinations, gate, id)
    gate_processed = _txn_present(_txn_backup(gate, id)) || isfile(_txn_absent(gate, id))
    gate_processed && _txn_present(gate) && rm(gate; force=true)
    for destination in destinations
        destination == gate || _txn_restore_one(destination, id)
    end
    _txn_restore_one(gate, id)
    _txn_fsync_parent(gate)
    _txn_cleanup(destinations, gate, id)
    return nothing
end

function _txn_recover(destinations, gate)
    marker = _txn_marker(gate)
    _txn_present(marker) || return nothing
    id, phase = _txn_read_marker(marker, _txn_signature(destinations))
    phase == "committed" ? _txn_cleanup(destinations, gate, id) :
        _txn_rollback(destinations, gate, id)
    return nothing
end

function _txn_inject(failpoint, interruptpoint, point)
    interruptpoint == point && throw(SimulatedTransactionInterruption(point))
    failpoint == point && error("injected transaction failure at $point")
    return nothing
end

function _with_staged_outputs(writer::Function, destinations;
        gate=last(destinations), failpoint=nothing, interruptpoint=nothing)
    ordered = unique(abspath.(String.(destinations)))
    isempty(ordered) && throw(ArgumentError("output transaction requires destinations"))
    gate = abspath(String(gate))
    gate in ordered || throw(ArgumentError("transaction gate must be a destination"))
    length(unique(dirname.(ordered))) == 1 ||
        throw(ArgumentError("transaction destinations must share one directory"))
    for destination in ordered
        isdir(destination) && !islink(destination) &&
            throw(ArgumentError("transaction destination is a directory: $destination"))
    end
    _txn_recover(ordered, gate)

    id = bytes2hex(sha256(string(time_ns(), ':', getpid(), ':', objectid(writer))))[1:16]
    staged = Dict(destination => _txn_stage(destination, id) for destination in ordered)
    marker = _txn_marker(gate)
    signature = _txn_signature(ordered)
    for destination in ordered
        for side_path in (staged[destination], _txn_backup(destination, id),
                          _txn_absent(destination, id))
            _txn_present(side_path) && error("transaction side path already exists: $side_path")
        end
    end
    try
        writer(staged)
        all(destination -> isfile(staged[destination]) && !islink(staged[destination]), ordered) ||
            error("one or more staged diagnostic outputs are missing or invalid")
        foreach(_txn_fsync_file, values(staged))
        _txn_write_marker(marker, id, "prepared", signature)
        _txn_inject(failpoint, interruptpoint, :after_marker)

        backup_order = vcat([gate], filter(!=(gate), ordered))
        for (index, destination) in enumerate(backup_order)
            if _txn_present(destination)
                _txn_rename(destination, _txn_backup(destination, id))
            else
                write(_txn_absent(destination, id), "absent\n")
                _txn_fsync_file(_txn_absent(destination, id))
            end
            _txn_fsync_parent(destination)
            _txn_inject(failpoint, interruptpoint, Symbol("after_backup_$index"))
        end

        install_order = vcat(filter(!=(gate), ordered), [gate])
        for (index, destination) in enumerate(install_order)
            _txn_rename(staged[destination], destination)
            _txn_fsync_parent(destination)
            point = destination == gate ? :after_gate : Symbol("after_install_$index")
            _txn_inject(failpoint, interruptpoint, point)
        end
        _txn_write_marker(marker, id, "committed", signature)
        _txn_inject(failpoint, interruptpoint, :after_commit)
        _txn_cleanup(ordered, gate, id)
    catch err
        err isa SimulatedTransactionInterruption && rethrow()
        if _txn_present(marker)
            try
                marker_id, phase = _txn_read_marker(marker, signature)
                marker_id == id || error("output transaction marker changed during commit")
                phase == "committed" ? _txn_cleanup(ordered, gate, id) :
                    _txn_rollback(ordered, gate, id)
            catch rollback_error
                error("output transaction failed and recovery failed: $(sprint(showerror, err)); " *
                      "recovery: $(sprint(showerror, rollback_error))")
            end
        else
            foreach(path -> rm(path; force=true), values(staged))
        end
        rethrow()
    end
    return nothing
end

function _write_summary(io, filename, case, bundle, png)
    reg, score = bundle.nominal.registration, bundle.nominal.scoring
    backbone_only = active_set_nnls(case.observed, [case.backbone])
    common_gain = backbone_only.sse - reg.fit.sse
    sequence = bundle.final_abstain ? "" : _sequence_text(score.best.sequence)
    stability = join(["$(item.label):$(_sequence_text(item.scoring.best.sequence)):" *
        "$(_float(item.scoring.best_runner_margin))" for item in bundle.perturbations], ";")
    negative = only(bundle.controls).scoring
    fields = (filename, length(case.centers), count(reg.mask), reg.direction, reg.phase,
        reg.mirror, join(bundle.nominal.boundary_parameters, ";"),
        _float(reg.shift_t_nm), _float(reg.shift_u_nm),
        _float(reg.rotation_deg), _float(reg.blur_sigma_nm), _float(common_gain),
        _float(reg.fit.sse), _float(score.incremental_contrast_gain),
        _float(score.best_runner_margin), _float(score.complement_margin),
        _float(score.best.contrast_coef), sequence, bundle.final_abstain,
        join(bundle.final_abstention_reasons, ";"), bundle.perturbation_stable,
        bundle.negative_control_beaten, stability,
        "shifted_gain=$(_float(negative.incremental_contrast_gain))", png,
        _registration_sha(bundle))
    println(io, join(fields, '\t'))
end

function _write_controls(io, filename, bundle)
    reg = bundle.nominal.registration
    for item in vcat(bundle.perturbations, bundle.controls)
        score = item.scoring
        println(io, join((filename, item.label, _registration_sha(bundle),
            reg.direction, reg.phase, reg.mirror, _sequence_text(score.best.sequence),
            _float(score.incremental_contrast_gain), _float(score.best_runner_margin),
            _float(score.complement_margin), score.abstain,
            join(score.abstention_reasons, ";")), '\t'))
    end
end

function main(args)
    opts = _parse_real_args(args)
    opts.help && (println(_real_usage()); return 0)
    data_dir = abspath(something(get(ENV, "STMFIT_DATA_DIR", nothing), "."))
    validate_real_diagnostic_paths(opts.summary, opts.molds, opts.config, data_dir)
    validate_real_observable_contract(opts.config;
        provenance_path=opts.observable_provenance,
        provenance_sha256=opts.observable_provenance_sha256,
        maps_path=opts.observable_map,
        invalid_mask_path=opts.validity_mask,
        molds_path=opts.molds,
        mold_provenance_path=opts.observable_mold_provenance,
        half_nm=opts.half_nm,
        step_nm=opts.step_nm)
    for file in opts.files
        isfile(joinpath(data_dir, file)) || throw(ArgumentError(
            "SXM file not found: $(joinpath(data_dir, file))"))
    end
    selected = read_selected_counts(opts.summary)
    all(file -> haskey(selected, file), opts.files) || throw(ArgumentError(
        "one or more files are missing from selection summary"))
    entry = load_converged_entry(opts.molds)
    set = derive_common_contrast(entry)
    config = Main.WholeRoiDiagnosticConfig.load_diagnostic_config(opts.config)
    mkpath(opts.outdir)
    summary_path = joinpath(opts.outdir, "diagnostic_summary.tsv")
    controls_path = joinpath(opts.outdir, "diagnostic_controls.tsv")
    png_paths = [joinpath(opts.outdir, "$(splitext(file)[1])_fused_fitmask.png")
                 for file in opts.files]
    destinations = vcat(png_paths, [controls_path, summary_path])
    _with_staged_outputs(destinations; gate=summary_path) do staged
        open(staged[summary_path], "w") do summary_io
            open(staged[controls_path], "w") do controls_io
                println(summary_io, "file\tn\tpixels\tdirection\tphase\tmirror\tboundary_parameters\tshift_t_nm\tshift_u_nm\trotation_deg\tblur_sigma_nm\tcommon_sse_gain\tcommon_only_sse\tincremental_contrast_gain\tbest_runner_margin\tcomplement_margin\tcontrast_coefficient\tsequence\tabstain\tabstention_reasons\tperturbation_stable\tnegative_control_beaten\tstability_details\tnegative_control_details\tpng_path\tregistration_sha256")
                println(controls_io, "file\tcontrol\tregistration_sha256\tdirection\tphase\tmirror\tbest_sequence\tincremental_contrast_gain\tbest_runner_margin\tcomplement_margin\tabstain\tabstention_reasons")
                for file in opts.files
                    case = _load_fused_fitmask_case(file, selected[file], data_dir, opts.outdir)
                    bundle = run_frozen_real_controls(case, set; config=config,
                        half_nm=opts.half_nm, step_nm=opts.step_nm,
                        control_shift_nm=opts.half_nm)
                    png = joinpath(opts.outdir, "$(splitext(file)[1])_fused_fitmask.png")
                    plot_frozen_real_diagnostic(file, case, set, bundle, staged[png];
                        half_nm=opts.half_nm, step_nm=opts.step_nm)
                    _write_summary(summary_io, file, case, bundle, png)
                    _write_controls(controls_io, file, bundle)
                    flush(summary_io); flush(controls_io)
                    println(file, ": abstain=", bundle.final_abstain,
                        " reasons=", join(bundle.final_abstention_reasons, ";"))
                end
            end
        end
    end
    println("wrote ", summary_path)
    println("wrote ", controls_path)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
