#!/usr/bin/env julia

using STMMolecularFit, GaussianFit2D, GaussianFit1D
using Printf, Statistics

const DATA_DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"
const BEST_DIR = "results/best_plots"
const OUT_TSV = joinpath(BEST_DIR, "bestN_stats_existing.tsv")
const EXCLUDE = Set(["240817_013.sxm", "240817_014.sxm"])
const FWHM_SIGMA = 2.355
const SIGMA_MIN_HARMONIZED_NM = 0.45 / FWHM_SIGMA
const SIGMA_MAX_HARMONIZED_NM = 1.20 / FWHM_SIGMA

function existing_files()
    files = String[]
    isdir(BEST_DIR) || return files
    for f in readdir(BEST_DIR)
        endswith(f, "_best.png") || continue
        fn = replace(f, "_best.png" => ".sxm")
        fn in EXCLUDE && continue
        push!(files, fn)
    end
    sort!(files)
    return files
end

function bestN_for_file(fn::String)
    fp = joinpath(DATA_DIR, fn)

    scfg = STMMolecularFit.SlideConfig(width_nm=0.70, support_threshold_fraction=0.20,
        support_noise_k=2.5, support_padding_nm=0.20, output_dir=BEST_DIR, no_plot=true)
    fcfg = STMMolecularFit.FitSlideConfig(min_spacing=0.35, max_spacing=0.75,
        output_dir=BEST_DIR, no_plot=true)

    img1 = STMMolecularFit.read_sxm(fp)
    slide = STMMolecularFit.extract_slide(img1, scfg)
    fit_1d = STMMolecularFit.fit_slide(slide, fcfg)
    best1d = GaussianFit1D.best_result(fit_1d.fit_run)

    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=BEST_DIR, no_plot=true)
    ccfg_ell = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14,
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        max_overlap=0.6,
        support_threshold_fraction=0.20, support_noise_k=2.5, support_padding_nm=0.20,
        global_maxtime=10.0, global_maxiter=10000, cv_folds=3,
        sigma_parallel_min_nm=SIGMA_MIN_HARMONIZED_NM,
        sigma_parallel_max_nm=SIGMA_MAX_HARMONIZED_NM,
        sigma_perp_min_nm=SIGMA_MIN_HARMONIZED_NM,
        sigma_perp_max_nm=SIGMA_MAX_HARMONIZED_NM,
        intelligent_sweep=true, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg_ell)
    ccfg_circ.chain_circular_sigmas = true

    img2 = GaussianFit2D.read_sxm(fp)
    res_ell, best_ell_raw, _ctx_ell = GaussianFit2D.chain_gaussian_sweep(img2, pcfg, ccfg_ell)
    res_circ, best_circ_raw, _ctx_circ = GaussianFit2D.chain_gaussian_sweep(img2, pcfg, ccfg_circ)
    valid_ell = [r for r in res_ell if r.success && r.valid && isfinite(r.bic)]
    valid_circ = [r for r in res_circ if r.success && r.valid && isfinite(r.bic)]
    best_ell = isempty(valid_ell) ? best_ell_raw : sort(valid_ell; by=r -> r.bic)[1]
    best_circ = isempty(valid_circ) ? best_circ_raw : sort(valid_circ; by=r -> r.bic)[1]

    return (file=fn, N_ell=best_ell.n, N_circ=best_circ.n, N_1D=best1d.n_peaks,
            d_circ_ell=best_circ.n - best_ell.n,
            d_1d_ell=best1d.n_peaks - best_ell.n,
            d_1d_circ=best1d.n_peaks - best_circ.n,
            bic_ell=best_ell.bic, bic_circ=best_circ.bic, sbic_1d=best1d.student_bic,
            valid_ell=best_ell.valid, valid_circ=best_circ.valid)
end

function counts(vals)
    d = Dict{Int,Int}()
    for v in vals
        d[v] = get(d, v, 0) + 1
    end
    return sort(collect(d); by=x -> x[1])
end

files = existing_files()
println("Files to analyze: $(length(files))")
println("Excluding: $(join(sort(collect(EXCLUDE)), ", "))")

rows = []
open(OUT_TSV, "w") do io
    println(io, join(["file", "N_ell", "N_circ", "N_1D", "circ-ell", "1D-ell", "1D-circ",
                      "BIC_ell", "BIC_circ", "sBIC_1D", "valid_ell", "valid_circ"], '\t'))
    for (i, fn) in enumerate(files)
        @printf("[%2d/%2d] %s ... ", i, length(files), fn)
        flush(stdout)
        try
            r = bestN_for_file(fn)
            push!(rows, r)
            println(io, join([r.file, r.N_ell, r.N_circ, r.N_1D, r.d_circ_ell, r.d_1d_ell, r.d_1d_circ,
                              round(r.bic_ell, digits=1), round(r.bic_circ, digits=1), round(r.sbic_1d, digits=1),
                              r.valid_ell, r.valid_circ], '\t'))
            flush(io)
            println("ell=$(r.N_ell) circ=$(r.N_circ) 1D=$(r.N_1D)")
        catch e
            println(io, join([fn, "ERR", "ERR", "ERR", "ERR", "ERR", "ERR", "ERR", "ERR", "ERR", "ERR", "ERR"], '\t'))
            flush(io)
            println("FAILED: $(sprint(showerror, e))")
        end
    end
end

if isempty(rows)
    println("No rows computed.")
    exit()
end

Nell = [r.N_ell for r in rows]
Ncirc = [r.N_circ for r in rows]
N1d = [r.N_1D for r in rows]
dc_e = [r.d_circ_ell for r in rows]
d1_e = [r.d_1d_ell for r in rows]
d1_c = [r.d_1d_circ for r in rows]

println("\n=== Best-N stats (excluding 013/014) ===")
println("n = $(length(rows))")
println("Exact matches: ell=circ $(count(==(0), dc_e))/$(length(rows)); 1D=ell $(count(==(0), d1_e))/$(length(rows)); 1D=circ $(count(==(0), d1_c))/$(length(rows)); all equal $(count(i -> Nell[i] == Ncirc[i] == N1d[i], eachindex(rows)))/$(length(rows))")
@printf("Mean Δ(circ-ell)=%.2f, mean |Δ|=%.2f\n", mean(dc_e), mean(abs.(dc_e)))
@printf("Mean Δ(1D-ell)=%.2f, mean |Δ|=%.2f\n", mean(d1_e), mean(abs.(d1_e)))
@printf("Mean Δ(1D-circ)=%.2f, mean |Δ|=%.2f\n", mean(d1_c), mean(abs.(d1_c)))
println("\nN_ell counts:  ", counts(Nell))
println("N_circ counts: ", counts(Ncirc))
println("N_1D counts:   ", counts(N1d))
println("\nΔ circ-ell counts: ", counts(dc_e))
println("Δ 1D-ell counts:   ", counts(d1_e))
println("Δ 1D-circ counts:  ", counts(d1_c))
println("\nWrote: $OUT_TSV")
