#!/usr/bin/env julia

using STMMolecularFit
using GaussianFit1D

filepath = length(ARGS) >= 1 ? ARGS[1] : error("Usage: julia extract_and_fit_slide.jl <filepath.sxm> [output_dir]")
outdir = length(ARGS) >= 2 ? ARGS[2] : "results/example_slide_fit"

img = read_sxm(filepath)
slide_cfg = SlideConfig(channel="Z", direction="fwd", width_nm=0.30,
                        support_threshold_fraction=0.08,
                        output_dir=joinpath(outdir, "slide"))
slide = extract_slide(img, slide_cfg)
slide_files = write_slide_outputs(slide, slide_cfg)

fit_cfg = build_config(Dict{String,Any}(
    "filepath" => slide_files.profile,
    "output_dir" => joinpath(outdir, "fit_1d"),
    "min_spacing" => 0.4,
    "max_spacing" => 0.675,
    "fwhm_min" => 0.45,
    "fwhm_max" => 1.2,
    "global_maxtime" => 8.0,
    "no_show" => true,
))

fr = run_fit(fit_cfg; save_cache=true, verbose=true)
best = best_result(fr)
update_model_rankings(fr.all_results, fit_cfg)
export_results(fr.x, fr.y, fr.all_results, fit_cfg)
plot_results(fr.x, fr.y, best, fr.all_results, fit_cfg)

mkpath(outdir)
open(joinpath(outdir, "model_selection.tsv"), "w") do io
    println(io, "n_peaks\tn_params\tbic\tdelta_bic\taic\tr_squared\trss\tselected")
    for r in sort(fr.all_results; by=r -> r.n_peaks)
        println(io, "$(r.n_peaks)\t$(r.n_params)\t$(r.bic)\t$(r.delta_bic)\t$(r.aic)\t$(r.r_squared)\t$(r.rss)\t$(r === best)")
    end
end

println("slide_profile: ", slide_files.profile)
println("support_length_nm: ", slide.support_length_nm)
println("best_n_peaks: ", best.n_peaks)
println("best_bic: ", best.bic)
