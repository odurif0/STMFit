#!/usr/bin/env julia
# Signal distribution comparison: 003 vs 034 vs 035.
using GaussianFit2D, Printf, Statistics

const DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"

for fn in ["240817_003.sxm", "240817_034.sxm", "240817_035.sxm"]
    fp = joinpath(DIR, fn)
    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp", no_plot=true)
    img = GaussianFit2D.read_sxm(fp)
    xs, ys, raw, z, z_smooth, unit, noise = GaussianFit2D.preprocess_channel(img,
        GaussianFit2D.get_channel(img, "Z"; direction="fwd"), pcfg)
    
    signal = z_smooth .- minimum(z_smooth)
    maxsig = maximum(signal)
    
    # Quantiles of the signal distribution
    q10 = quantile(vec(signal), 0.10)
    q25 = quantile(vec(signal), 0.25)
    q50 = quantile(vec(signal), 0.50)
    q75 = quantile(vec(signal), 0.75)
    q90 = quantile(vec(signal), 0.90)
    q95 = quantile(vec(signal), 0.95)
    q99 = quantile(vec(signal), 0.99)
    
    threshold_035 = 0.35 * maxsig
    
    # Fraction above threshold
    frac_above = count(signal .>= threshold_035) / length(signal)
    
    println("\n=== $fn ===")
    @printf("  max signal: %.4f nm\n", maxsig)
    @printf("  Signal quantiles: q10=%.4f q25=%.4f q50=%.4f q75=%.4f q90=%.4f q95=%.4f q99=%.4f\n",
        q10, q25, q50, q75, q90, q95, q99)
    @printf("  threshold(0.35×max): %.4f nm → %.1f%% above\n", threshold_035, frac_above*100)
    
    # What fraction would give ~15% ROI?
    for f in [0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
        th = f * maxsig
        frac = count(signal .>= th) / length(signal)
        @printf("    frac=%.2f → threshold=%.4f → %.1f%% above\n", f, th, frac*100)
    end
end
