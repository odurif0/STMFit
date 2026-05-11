#!/usr/bin/env julia
# Quick 1D vs 2D N comparison for valid files from batch triage.
# Reads existing 2D TSV results, runs 1D extraction + fit, reports comparison.
# Usage: julia --project=. scripts/compare_1d_2d_valid.jl

using STMMolecularFit
using GaussianFit2D
using GaussianFit1D
using DelimitedFiles
using Printf

const DATA_DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"
const TSV = "results/batch_triage_20240817_relaxed.tsv"
const N_TOP = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 15

# Load 2D results, pick top N_TOP by BIC among valid
data = readdlm(TSV)
candidates = Tuple{String,Int,Float64}[]
for i in 2:size(data, 1)
    data[i, 5] isa Bool && data[i, 5] || continue
    push!(candidates, (string(data[i, 1]), Int(data[i, 3]), Float64(data[i, 4])))
end
sort!(candidates, by = x -> x[3])
top = candidates[1:min(N_TOP, length(candidates))]

println("=== 1D vs 2D comparison — top $(length(top)) valid files ===")
println(rpad("File", 24), rpad("N_1D", 6), rpad("sBIC_1D", 10), 
        rpad("N_2D", 6), rpad("sBIC_2D", 10), rpad("ΔN", 6), rpad("spacing", 10))
println(repeat("-", 70))

for (fn, n2d, bic2d) in top
    fp = joinpath(DATA_DIR, fn)
    try
        # 1D slide extraction + fit
        img = STMMolecularFit.read_sxm(fp)
        slide_cfg = STMMolecularFit.SlideConfig(width_nm=0.30, support_threshold_fraction=0.20,
                                                 support_noise_k=2.5, support_padding_nm=0.20,
                                                 output_dir="/tmp/compare_1d2d", no_plot=true)
        slide = STMMolecularFit.extract_slide(img, slide_cfg)
        fit_cfg = STMMolecularFit.FitSlideConfig(min_spacing=0.35, max_spacing=0.75,
                                                  output_dir="/tmp/compare_1d2d")
        fit_1d = STMMolecularFit.fit_slide(slide, fit_cfg)
        best1d = GaussianFit1D.best_result(fit_1d.fit_run)
        n1d = best1d.n_peaks
        sbic1d = best1d.student_bic
        dn = n2d - n1d
        mark = dn == 0 ? "✓" : @sprintf("%+d", dn)
        
        # Also get 2D spacing from TSV
        sp2d = Float64(data[findfirst(row -> string(row[1]) == fn, eachrow(data)), 7])
        
        println(rpad(fn, 24), rpad(n1d, 6), rpad(round(sbic1d, digits=0), 10),
                rpad(n2d, 6), rpad(round(bic2d, digits=0), 10), rpad(mark, 6),
                rpad(round(sp2d, digits=3), 10))
    catch e
        s = sprint(showerror, e)
        println(rpad(fn, 24), "1D FAILED: ", s[1:min(50, length(s))])
    end
end

# Summary stats
matches = 0
for (fn, n2d, bic2d) in top
    n1d_file = joinpath("/tmp/compare_1d2d", fn)
end
