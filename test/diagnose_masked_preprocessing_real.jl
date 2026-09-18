#!/usr/bin/env julia
# julia --startup-file=no --threads=1 --project=. test/diagnose_masked_preprocessing_real.jl --help
include(joinpath(@__DIR__, "lib", "masked_preprocessing_comparison.jl"))
using .MaskedPreprocessingComparison
using LinearAlgebra

function masked_comparison_help(io=stdout)
    println(io,"""
    Usage (Julia 1.13 only):
      julia --startup-file=no --threads=1 --project=. test/diagnose_masked_preprocessing_real.jl \
        --file RAW.sxm --geometry BASE.tsv --selected-summary SUMMARY.tsv \
        --config PHYSICAL.toml --acquisition-settings ACQ.toml \
        --settings MASK.toml --outdir NEW_DIRECTORY [--dry-run]

    All seven paths are required. Duplicate/unknown flags and existing output
    paths or symlinked output ancestors are rejected. No output overwrite.
    --dry-run reads only config/summary/geometry metadata, checks file identity
    and saved selected/adaptive context, and creates no output. No raw parser,
    pixels, fitting or image diagnostics run in dry-run.

    Real execution belongs on Viper. Four fixed methods: native_reference,
    finite_only_ols, guarded_ols, guarded_huber, each forward/backward. Original
    native ROI, saved geometry and support stay frozen across methods. Native
    plane+rows unsmoothed reference requires nm; raw missingness is not filled
    as observations. Native backward flip occurs once per read, never twice.
    Registration retains ±8 x / ±2 y (85 lags), same support at every lag.
    Own, all-common and native-pair scopes retain coverage losses and failures.
    Grid-best shifts and affine discrepancies are diagnostics, not corrections.
    Optional physical level-aligned RMS uses one shared-background constant per
    view, fixed for all lags. The frozen loader preprocesses auxiliary Current
    arrays, then discards them. No Current evidence is computed or exported by
    this comparison. No saved-fit residual or noise-calibration stage.
    No count fit, chemical prediction, benchmark grade or production change.

    Outputs: view_status.tsv, row_status.tsv, comparison_summary.tsv,
    lag_candidates.tsv (including blocked candidates), local_comparisons.tsv,
    background_summary.tsv (post-fit/descriptive, not noise calibration),
    arrays.jls, metadata.toml. arrays.jls is Julia 1.13 stdlib Serialization of
    plain arrays/scalars/dicts only; z[y,x], raw scaled nm on the analysis grid,
    corrected views/masks/plane/rows, frozen ROI/support/footprint/centers and
    per-scope masks/parameters. Load trusted files only. Native unknowns are NA.
    Inputs/sources and array/TSV hashes are recorded. No winner is selected.
    """)
end

function masked_comparison_options(args)
    paths=("--file","--geometry","--selected-summary","--config","--acquisition-settings","--settings","--outdir")
    opts=Dict{String,String}()
    dry_run=false
    i=1
    while i<=length(args)
        arg=args[i]
        if arg=="--dry-run"
            dry_run && error("Duplicate argument: --dry-run")
            dry_run=true
            i+=1
            continue
        end
        arg in paths || error("Unknown argument: $arg")
        haskey(opts,arg) && error("Duplicate argument: $arg")
        i<length(args) && !startswith(args[i+1],"--") && !isempty(strip(args[i+1])) || error("Missing value for $arg")
        opts[arg]=args[i+1]
        i+=2
    end
    all(k->haskey(opts,k),paths) || error("All seven paths are required; use --help")
    for k in paths[1:6]
        isfile(opts[k]) || error("Input is not a file: $(opts[k])")
    end
    require_new_output(opts["--outdir"])
    return (; paths=opts, dry_run)
end

function masked_comparison_main(args=ARGS)
    require_julia()
    if args in (["--help"],["-h"])
        masked_comparison_help()
        return nothing
    end
    opts=masked_comparison_options(args)
    p=opts.paths
    context=metadata_context(p["--file"],p["--geometry"],p["--config"],p["--selected-summary"],
        p["--acquisition-settings"],p["--settings"])
    if opts.dry_run
        println("Metadata-only dry-run: ",context.file,"; saved N=",context.selected.n,
            "; adaptive rescue=",context.selected.use_rescue,"; 4 methods; 85 lags; no raw read/output")
        return context
    end
    BLAS.set_num_threads(1)
    loaded=load_real_context(p["--file"],p["--geometry"],p["--config"],p["--selected-summary"],context)
    case,rawdata,footprint=loaded.case,loaded.rawdata,loaded.footprint
    views=build_views(case,rawdata,footprint,context.options)
    d=compare_views(case,rawdata,views,footprint,context.s)
    inputs=Dict(k[3:end]=>v for (k,v) in p if k!="--outdir")
    write_comparison(p["--outdir"],case,rawdata,views,footprint,d,context;inputs=inputs)
    println("Diagnostic comparison only: ",p["--outdir"])
    println("View records=",length(d.view_rows),"; all-common pixels=",count(d.all_common),
        "; retained lag rows=",length(d.candidates),"; no accepted shift or method selected")
    return d
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && masked_comparison_main()
