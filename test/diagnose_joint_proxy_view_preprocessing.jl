#!/usr/bin/env julia

include(joinpath(@__DIR__, "diagnose_joint_proxy_type_collapse.jl"))

const VALID_FLATTEN_MODES = ("none", "plane", "rows", "plane+rows")
const PREPROCESSING_TRANSFER_FIELDS = [:file, :flatten_mode, :status, :n, :lobe,
    :n_pairs, :corr, :slope, :intercept, :affine_nrmse, :best_dt_steps,
    :best_du_steps, :best_dt_nm, :best_du_nm, :best_corr, :corr_gain]

Base.@kwdef struct PreprocessingDiagnosticOptions
    config::String
    data_dir::String
    file::String
    out::String
    flatten_modes::Vector{String}
    max_shift_steps::Int
    seed::Int = 20260721
end

function parse_preprocessing_options(args)
    values = Dict{String,String}()
    i = 1
    while i <= length(args)
        flag = String(args[i])
        flag in ("--config", "--data-dir", "--file", "--out", "--flatten-modes",
            "--max-shift-steps", "--seed") || error("unknown argument: $flag")
        i < length(args) || error("$flag requires a value")
        values[flag] = String(args[i + 1])
        i += 2
    end
    for flag in ("--config", "--data-dir", "--file", "--out", "--flatten-modes",
        "--max-shift-steps")
        haskey(values, flag) || error("$flag is required")
    end
    file = values["--file"]
    basename(file) == file && endswith(lowercase(file), ".sxm") ||
        error("--file must be a bare .sxm basename")
    flatten_modes = lowercase.(strip.(split(values["--flatten-modes"], ',')))
    isempty(flatten_modes) && error("--flatten-modes must not be empty")
    length(unique(flatten_modes)) == length(flatten_modes) ||
        error("--flatten-modes contains duplicates")
    all(mode -> mode in VALID_FLATTEN_MODES, flatten_modes) ||
        error("--flatten-modes must use: $(join(VALID_FLATTEN_MODES, ','))")
    max_shift_steps = parse(Int, values["--max-shift-steps"])
    max_shift_steps >= 0 || error("--max-shift-steps must be non-negative")
    return PreprocessingDiagnosticOptions(
        config=values["--config"], data_dir=values["--data-dir"], file=file,
        out=values["--out"], flatten_modes=flatten_modes,
        max_shift_steps=max_shift_steps,
        seed=parse(Int, get(values, "--seed", "20260721")),
    )
end

function pattern_with_flatten(base::GaussianFit2D.PatternConfig, flatten::String)
    flatten in VALID_FLATTEN_MODES || error("invalid flatten mode: $flatten")
    names = fieldnames(typeof(base))
    values = NamedTuple{names}(ntuple(i -> getfield(base, names[i]), length(names)))
    return GaussianFit2D.PatternConfig(; merge(values, (; flatten))...)
end

function preprocessing_transfer_rows(file::String, flatten_mode::String, candidate,
    evidence, registry::ProxyEnsemble, max_shift_steps::Int)
    length(candidate.lobes) == length(evidence) || error("candidate/evidence lobe mismatch")
    rows = NamedTuple[]
    for (lobe, lobe_evidence) in zip(candidate.lobes, evidence)
        patches = lobe_evidence.view_patches
        if !haskey(patches, "fwd") || !haskey(patches, "bwd")
            push!(rows, (; file, flatten_mode, status="bwd_missing", n=candidate.n,
                lobe=lobe.index, n_pairs=0, corr=NaN, slope=NaN, intercept=NaN,
                affine_nrmse=NaN, best_dt_steps=0, best_du_steps=0,
                best_dt_nm=0.0, best_du_nm=0.0, best_corr=NaN, corr_gain=NaN))
            continue
        end
        metrics = view_transfer_metrics(patches["fwd"], patches["bwd"],
            registry.grid_n, max_shift_steps)
        push!(rows, (; file, flatten_mode, status="ok", n=candidate.n,
            lobe=lobe.index, metrics.n_pairs, metrics.corr, metrics.slope,
            metrics.intercept, metrics.affine_nrmse, metrics.best_dt_steps,
            metrics.best_du_steps,
            best_dt_nm=metrics.best_dt_steps * registry.grid_step_nm,
            best_du_nm=metrics.best_du_steps * registry.grid_step_nm,
            metrics.best_corr, metrics.corr_gain))
    end
    return rows
end

function fixed_geometry_case(config::String, data_dir::String, file::String,
    output_dir::String, seed::Int)
    registry = load_registry(config)
    path = joinpath(data_dir, file)
    image = STMSXMIO.read_sxm(path)
    base_pattern = Inference._default_pattern_config(path, output_dir)
    chain = diagnostic_chain_config(diagnostic_file_seed(seed, file))
    results, _best, context = redirect_stdout(devnull) do
        GaussianFit2D.chain_gaussian_sweep(image, base_pattern, chain)
    end
    base_views = Inference._build_views(image, base_pattern,
        Main.JointProxyCandidateViews.build_view_data)
    report = Main.JointProxyCandidateViews.extract_candidate_views(results, context, chain,
        base_views; patch_half_nm=registry.grid_half_nm, patch_step_nm=registry.grid_step_nm)
    isempty(report.candidates) && error("no valid fixed-geometry candidate for $file")
    candidate = argmin(item -> item.joint_gcv, report.candidates)
    return (; registry, image, base_pattern, chain, candidate, context)
end

function run_preprocessing_diagnostic(options::PreprocessingDiagnosticOptions)
    case = fixed_geometry_case(options.config, options.data_dir, options.file,
        dirname(options.out), options.seed)
    registry, image = case.registry, case.image
    base_pattern, candidate = case.base_pattern, case.candidate
    options.max_shift_steps < registry.grid_n ||
        error("--max-shift-steps must be smaller than registry grid_n=$(registry.grid_n)")

    rows = NamedTuple[]
    for flatten_mode in options.flatten_modes
        pattern = pattern_with_flatten(base_pattern, flatten_mode)
        views = Inference._build_views(image, pattern,
            Main.JointProxyCandidateViews.build_view_data)
        evidence = image_lobe_evidence(candidate, views)
        append!(rows, preprocessing_transfer_rows(options.file, flatten_mode, candidate,
            evidence, registry, options.max_shift_steps))
    end
    write_tsv(options.out, rows; fields=PREPROCESSING_TRANSFER_FIELDS)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_preprocessing_diagnostic(parse_preprocessing_options(ARGS))
end
