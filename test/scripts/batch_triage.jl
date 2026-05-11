#!/usr/bin/env julia
# Batch triage: try 2D chain sweep (consensus Z/Current) on every .sxm file,
# recording success/failure in a TSV summary.
# Usage: julia --project=. scripts/batch_triage.jl <directory> [output.tsv] [N_files] [--chunk i/n]

using STMMolecularFit
using GaussianFit2D
using Printf
using Dates

const FWHM_SIGMA = 2.355
const SIGMA_MIN_HARMONIZED_NM = 0.45 / FWHM_SIGMA
const SIGMA_MAX_HARMONIZED_NM = 1.20 / FWHM_SIGMA

function _parse_cli(args)
    length(args) >= 1 || error("Usage: julia batch_triage.jl <directory> [output.tsv] [N_files] [--chunk i/n]")
    data_dir = args[1]
    out_tsv = length(args) >= 2 && !startswith(args[2], "--") ? args[2] : "results/batch_triage_$(Dates.format(now(), "YYYYmmdd_HHMM")).tsv"
    n_files = typemax(Int)
    chunk_idx = 1
    chunk_total = 1
    i = 3
    while i <= length(args)
        arg = args[i]
        if arg == "--chunk"
            i < length(args) || error("--chunk requires i/n")
            chunk_arg = args[i+1]
            i += 2
        elseif startswith(arg, "--chunk=")
            chunk_arg = split(arg, "=", limit=2)[2]
            i += 1
        elseif startswith(arg, "--")
            error("Unknown option: $arg")
        else
            n_files = parse(Int, arg)
            i += 1
            continue
        end
        parts = split(chunk_arg, "/")
        length(parts) == 2 || error("Invalid --chunk '$chunk_arg'; expected i/n")
        chunk_idx = parse(Int, parts[1])
        chunk_total = parse(Int, parts[2])
    end
    1 <= chunk_idx <= chunk_total || error("chunk index must satisfy 1 <= i <= n")
    return data_dir, out_tsv, n_files, chunk_idx, chunk_total
end

const DATA_DIR, OUT_TSV, N_FILES, CHUNK_IDX, CHUNK_TOTAL = _parse_cli(ARGS)
mkpath(dirname(abspath(OUT_TSV)))

function find_sxm(dir::String)
    files = String[]
    for (root, _, names) in walkdir(dir)
        for f in names
            endswith(lowercase(f), ".sxm") && push!(files, joinpath(root, f))
        end
    end
    return sort(files)
end

function has_channel(img, name::String)
    return any(c -> lowercase(c.name) == lowercase(name), img.channels)
end

function safe_fit(filepath::String, pcfg::GaussianFit2D.PatternConfig, ccfg::GaussianFit2D.ChainSweepConfig)
    try
        img = GaussianFit2D.read_sxm(filepath)
        use_consensus = has_channel(img, "Current")
        if use_consensus
            cres = GaussianFit2D.fit_chain_consensus(img, pcfg, ccfg)
            z_best = cres.z.best
            c_best = cres.current !== nothing ? cres.current.best : nothing
            c_n = c_best !== nothing ? c_best.n : missing
            consensus = cres.consensus
            agreement = cres.agreement
        else
            # No Current channel — run Z-only sweep
            results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
            z_best = best
            c_n = missing
            consensus = false
            agreement = "no Current channel"
        end
        return (filepath=basename(filepath), status="ok",
                z_n=z_best.n, z_bic=z_best.bic, z_valid=z_best.valid, z_reason=z_best.reason,
                z_spacing=z_best.mean_spacing_nm, z_spar=z_best.sigma_parallel_nm,
                z_sperp=z_best.sigma_perp_nm, z_chi2=z_best.chi2_reduced,
                c_n=c_n, consensus=consensus, agreement=agreement)
    catch e
        return (filepath=basename(filepath), status="error",
                z_n=missing, z_bic=NaN, z_valid=false, z_reason=string(e),
                z_spacing=NaN, z_spar=NaN, z_sperp=NaN, z_chi2=NaN,
                c_n=missing, consensus=false, agreement="error: $(typeof(e))")
    end
end

# ── Configs ──
pcfg = GaussianFit2D.PatternConfig(
    filepath="", channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1,
    output_dir="results/batch_triage")

ccfg = GaussianFit2D.ChainSweepConfig(
    n_min=2, n_max=14,
    spacing_min_nm=0.35, spacing_max_nm=0.75,
    max_overlap=0.6,
    fit_width_nm=0.15,
    support_threshold_fraction=0.20, support_noise_k=2.5,
    support_padding_nm=0.20,
    sigma_parallel_min_nm=SIGMA_MIN_HARMONIZED_NM,
    sigma_parallel_max_nm=SIGMA_MAX_HARMONIZED_NM,
    sigma_perp_min_nm=SIGMA_MIN_HARMONIZED_NM,
    sigma_perp_max_nm=SIGMA_MAX_HARMONIZED_NM,
    global_maxtime=5.0, global_maxiter=5000, cv_folds=3,
    intelligent_sweep=true, fuse_z_bwd=true)

# ── Find files, skip already-processed ──
sxm_all = find_sxm(DATA_DIR)
sxm_base = sxm_all[1:min(N_FILES, length(sxm_all))]
sxm_files = CHUNK_TOTAL == 1 ? sxm_base : [fp for (i, fp) in enumerate(sxm_base) if mod1(i, CHUNK_TOTAL) == CHUNK_IDX]
CHUNK_TOTAL > 1 && @printf("Chunk %d/%d: %d of %d selected files\n", CHUNK_IDX, CHUNK_TOTAL, length(sxm_files), length(sxm_base))
already = Set{String}()
if isfile(OUT_TSV)
    for line in eachline(OUT_TSV)
        startswith(line, "filepath") && continue
        isempty(strip(line)) && continue
        fname = split(line, '\t')[1]
        isempty(fname) && continue
        push!(already, fname)
    end
end
to_process = [fp for fp in sxm_files if !(basename(fp) in already)]
println("Found $(length(sxm_all)) .sxm files, selected $(length(sxm_files)), $(length(already)) already done, $(length(to_process)) to process")

# ── Process and write progressively ──
println("Writing TSV to $OUT_TSV")
header = ["filepath", "status", "z_n", "z_bic", "z_valid", "z_reason",
          "z_spacing_nm", "z_spar_nm", "z_sperp_nm", "z_chi2",
          "current_n", "consensus", "agreement"]
if !isfile(OUT_TSV) || isempty(read(OUT_TSV))
    open(OUT_TSV, "w") do io
        println(io, join(header, '\t'))
    end
end

function _process_batch(sxm_files, pcfg, ccfg, out_tsv)
    ok = 0; err = 0; valid = 0; cons = 0
    t0 = time()
    for (i, fp) in enumerate(sxm_files)
        @printf("[%3d/%3d] %s ... ", i, length(sxm_files), basename(fp))
        flush(stdout)
        r = safe_fit(fp, pcfg, ccfg)
        open(out_tsv, "a") do io
            line = join([r.filepath, r.status, r.z_n, r.z_bic, r.z_valid, r.z_reason,
                         r.z_spacing, r.z_spar, r.z_sperp, r.z_chi2,
                         r.c_n, r.consensus, r.agreement], '\t')
            println(io, line)
            flush(io)
        end
        if r.status == "ok"
            ok += 1
            println("N=$(r.z_n) valid=$(r.z_valid) consensus=$(r.consensus)")
            r.z_valid && (valid += 1)
            r.consensus && (cons += 1)
        else
            err += 1
            println("FAILED: $(r.z_reason)")
        end
    end
    return (ok=ok, err=err, valid=valid, cons=cons, elapsed=time()-t0)
end

counts = _process_batch(to_process, pcfg, ccfg, OUT_TSV)

# ── Summary ──
println("\n" * "="^60)
println("BATCH TRIAGE SUMMARY: $(length(sxm_files)) files ($(length(already)) already done + $(length(to_process)) new) in $(round(counts.elapsed, digits=1)) s")
println("="^60)
println("  OK (new):    $(counts.ok)")
println("  Errors (new):$(counts.err)")
println("  Valid (new): $(counts.valid)")
println("  Consensus:   $(counts.cons)")
println("\nTSV written to: $OUT_TSV")
