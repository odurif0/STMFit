#!/usr/bin/env julia
# Julia 1.13, root project environment; no Python subprocesses.
# See lib/empirical_fisher_native.jl for the retained legacy conventions and
# the deliberate native initialization/RNG differences from sklearn.
module NativeFisherCLI

include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
using .EmpiricalFisherNative
using .EmpiricalFisherNative.ReconstructedUnitAssignment: load_training_mask

const USAGE = """
Usage: julia --project=. test/build_empirical_fisher_native.jl \\
    --patches PATH --prefix res --config CONFIG --out PATH

Build one cross-validated native Fisher margin TSV from forward patches.
Required output columns: file, lobe, score, invalid_reason. All input keys remain;
invalid patches and missing/degenerate opposite folds receive NA plus a reason.
The output must not exist. No full-data model or prediction TSV is created.

All four options are required. --prefix selects numbered patch columns (res_p001,
etc.). Config supplies PCA, GMM, regularization, seed, layout and patch projection.
The score origin is explicit: legacy centered mean or original training mean.
Config selects lobe_parity or scan_hash_twofold, with an explicit scan split seed.
Scan mode assigns whole scans before invalid-row filtering; all PCA/GMM/Fisher
learning and amplitude naming use the opposite scan group only. Split logs
contain identities/counts, not labels or a benchmark manifest.
With complete_patches training, --training-support PATH is required: only full
forward squares train either parity fold; admissible partial rows are scored.
Affine projection, when selected, applies to training and held-out scoring, not
to the original center amplitudes used to order the unsupervised clusters.
Native Julia and sklearn initialization/RNG differ; byte identity is NOT claimed.
"""

function main(args=ARGS)
    if args == ["--help"] || args == ["-h"]
        print(USAGE)
        return 0
    end
    required = ("--patches", "--prefix", "--config", "--out")
    allowed = (required..., "--training-support")
    options = Dict{String,String}()
    iseven(length(args)) || throw(ArgumentError("expected option/value pairs; use --help"))
    for i in 1:2:length(args)
        name, value = args[i], args[i + 1]
        name in allowed || throw(ArgumentError("unknown option: $name"))
        haskey(options, name) && throw(ArgumentError("duplicate option: $name"))
        isempty(value) && throw(ArgumentError("empty value for $name"))
        options[name] = value
    end
    for name in required
        haskey(options, name) || throw(ArgumentError("missing $name; use --help"))
    end
    ispath(options["--out"]) && throw(ArgumentError("output already exists: $(options["--out"])"))
    config = load_fisher_config(options["--config"])
    patches = load_patches(options["--patches"], options["--prefix"], config)
    mask = load_training_mask(config.training_support, get(options, "--training-support", ""), patches.keys, "fisher")
    groups = fisher_fold_ids(patches.keys, config)
    println("Fisher grouping: ", config.cv_scheme, " split_seed=", config.scan_split_seed)
    if config.cv_scheme == "scan_hash_twofold"
        mapping = Dict(key[1] => groups[i] for (i,key) in enumerate(patches.keys))
        for file in sort(collect(keys(mapping)))
            println("Fisher scan ", file, " group=", mapping[file])
        end
    end
    diagnostics = []
    rows = cv_scores(patches, config; training_mask=mask, diagnostics)
    for d in diagnostics
        train_files = Set(patches.keys[i][1] for i in d.train)
        held_files = Set(patches.keys[i][1] for i in d.held)
        println("Fisher fold ", d.fold, " train_rows=", length(d.train), " held_rows=", length(d.held),
            " train_scans=", length(train_files), " held_scans=", length(held_files),
            " overlap_scans=", length(intersect(train_files,held_files)), " status=", d.status,
            " reason=", isempty(d.reason) ? "none" : d.reason,
            " converged=", d.model === nothing ? "NA" : string(d.model.gmm.converged))
    end
    mask === nothing || println("Fisher complete-square eligibility: $(count(mask))/$(length(mask)) rows before feature validity")
    write_scores(options["--out"], rows)
    invalid = count(row -> !isempty(row.invalid_reason), rows)
    println("Fisher CV: $(options["--out"]) ($(length(rows)) rows, $invalid invalid)")
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(NativeFisherCLI.main())
    catch error
        error isa InterruptException && rethrow()
        print(stderr, "Fisher builder: ")
        showerror(stderr, error)
        println(stderr)
        exit(1)
    end
end
