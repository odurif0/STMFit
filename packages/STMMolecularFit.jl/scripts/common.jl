const DEFAULT_DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const FWHM_SIGMA = 2.355
const FWHM_MIN_HARMONIZED_NM = 0.45
const FWHM_MAX_HARMONIZED_NM = 1.20
const SIGMA_MIN_HARMONIZED_NM = FWHM_MIN_HARMONIZED_NM / FWHM_SIGMA
const SIGMA_MAX_HARMONIZED_NM = FWHM_MAX_HARMONIZED_NM / FWHM_SIGMA
const DEFAULT_MAX_OVERLAP = 0.60
const DEFAULT_SPACING_MIN_NM = 0.35
const DEFAULT_SPACING_MAX_NM = 0.75

function list_sxm_files(data_dir::String=DEFAULT_DATA_DIR)
    isdir(data_dir) || error("Data directory not found: $data_dir. Set STMFIT_DATA_DIR or pass a path.")
    files = [f for f in readdir(data_dir) if endswith(lowercase(f), ".sxm")]
    return sort(files)
end

function parse_chunk_arg(args; start_index::Int=1)
    chunk_idx = 1
    chunk_total = 1
    i = start_index
    while i <= length(args)
        arg = args[i]
        if arg == "--chunk"
            i < length(args) || error("--chunk requires i/n")
            chunk_arg = args[i + 1]
            i += 2
        elseif startswith(arg, "--chunk=")
            chunk_arg = split(arg, "=", limit=2)[2]
            i += 1
        else
            i += 1
            continue
        end
        parts = split(chunk_arg, "/")
        length(parts) == 2 || error("Invalid --chunk '$chunk_arg'; expected i/n")
        chunk_idx = parse(Int, parts[1])
        chunk_total = parse(Int, parts[2])
    end
    chunk_total >= 1 || error("chunk total must be >= 1")
    1 <= chunk_idx <= chunk_total || error("chunk index must satisfy 1 <= i <= n")
    return chunk_idx, chunk_total
end

function apply_chunk(items::AbstractVector, chunk_idx::Int, chunk_total::Int)
    chunk_total == 1 && return collect(items)
    return [item for (i, item) in enumerate(items) if mod1(i, chunk_total) == chunk_idx]
end

function make_chain_config(; circular::Bool=false, maxtime::Float64=10.0, maxiter::Int=10000)
    return GaussianFit2D.ChainSweepConfig(
        n_min=2, n_max=14,
        spacing_min_nm=DEFAULT_SPACING_MIN_NM,
        spacing_max_nm=DEFAULT_SPACING_MAX_NM,
        max_overlap=DEFAULT_MAX_OVERLAP,
        fit_width_nm=0.15,
        support_threshold_fraction=0.20,
        support_noise_k=2.5,
        support_padding_nm=0.20,
        global_maxtime=maxtime,
        global_maxiter=maxiter,
        cv_folds=3,
        sigma_parallel_min_nm=SIGMA_MIN_HARMONIZED_NM,
        sigma_parallel_max_nm=SIGMA_MAX_HARMONIZED_NM,
        sigma_perp_min_nm=SIGMA_MIN_HARMONIZED_NM,
        sigma_perp_max_nm=SIGMA_MAX_HARMONIZED_NM,
        intelligent_sweep=true,
        fuse_z_bwd=true,
        chain_circular_sigmas=circular,
    )
end

function best_valid_or_raw(results, raw_best)
    valid = [r for r in results if getproperty(r, :success) && getproperty(r, :valid) && isfinite(getproperty(r, :bic))]
    isempty(valid) && return raw_best
    return sort(valid; by=r -> r.bic)[1]
end
