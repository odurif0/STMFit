#!/usr/bin/env julia
using STMMolecularFit, Printf, Statistics

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_006.sxm"
img = STMMolecularFit.read_sxm(FILE)

# Extract raw slide
scfg = STMMolecularFit.SlideConfig(width_nm=0.30, support_noise_k=2.5, support_padding_nm=0.20,
    output_dir="/tmp", no_plot=true)
slide = STMMolecularFit.extract_slide(img, scfg)
@printf("Support: %.2f nm\n", slide.support_length_nm)

# Now manually compute thresholds
dist, y_raw = slide.x, slide.y
ys = STMMolecularFit._smooth1d(y_raw, 5)
baseline = quantile(ys, 0.10)
peak = maximum(ys)
low = ys[ys .<= baseline]
noise = isempty(low) ? STMMolecularFit._mad_std(ys) : STMMolecularFit._mad_std(low)
signal_above = ys .- baseline
signal_pos = signal_above[signal_above .> 0]
otsu_val = STMMolecularFit._otsu_threshold_1d(signal_pos)
threshold_otsu = baseline + otsu_val
threshold_noise = baseline + 2.5 * noise
threshold_frac = baseline + 0.20 * (peak - baseline)

@printf("baseline=%.4f  peak=%.4f  noise=%.4f\n", baseline, peak, noise)
@printf("threshold_otsu=%.4f  threshold_noise=%.4f  threshold_frac(0.20)=%.4f\n",
    threshold_otsu, threshold_noise, threshold_frac)
@printf("Dominant: %s\n", threshold_otsu > threshold_noise ? "Otsu" : "noise_k")
