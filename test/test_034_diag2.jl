#!/usr/bin/env julia
# Proper diagnostics on 034 vs 003 in nm.
using GaussianFit2D, Printf, Statistics

const DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"

for fn in ["240817_003.sxm", "240817_034.sxm"]
    fp = joinpath(DIR, fn)
    img = GaussianFit2D.read_sxm(fp)
    
    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp", no_plot=true)
    
    ch = GaussianFit2D.get_channel(img, "Z"; direction="fwd")
    xs, ys, raw, z, z_smooth, unit, noise = GaussianFit2D.preprocess_channel(img, ch, pcfg)
    
    println("\n=== $fn ===")
    @printf("  Image: %d × %d px, range=(%.2f, %.2f) nm\n", img.width, img.height, img.range_nm...)
    @printf("  After preprocessing (nm): min=%.4f  max=%.4f  std=%.4f  noise=%.5f\n",
        minimum(z_smooth), maximum(z_smooth), std(z_smooth), noise)
    @printf("  Noise/signal: %.4f\n", noise / (maximum(z_smooth) - minimum(z_smooth) + 1e-12))
    
    # Check for double-tip via row autocorrelation in nm
    println("  Row autocorrelation peaks (shift in px → nm):")
    n_rows, n_cols = size(z_smooth)
    shift_votes = Dict{Int,Int}()
    for iy in 1:min(n_rows, 200)
        row = z_smooth[iy, :] .- mean(z_smooth[iy, :])
        s = std(row)
        s < 1e-10 && continue
        best_shift, best_val = 0, -1.0
        for sh in 1:min(50, n_cols÷2)
            c = cor(row[1:end-sh], row[(1+sh):end])
            c > 0.3 && (shift_votes[sh] = get(shift_votes, sh, 0) + 1)
            if c > best_val
                best_val = c
                best_shift = sh
            end
        end
    end
    px = img.range_nm[1] / img.width
    peaks = sort(collect(shift_votes); by=x->-x[2])
    for (sh, cnt) in peaks[1:min(5, end)]
        @printf("    shift=%d px (%.3f nm): %d rows with autocorr > 0.3\n", sh, sh*px, cnt)
    end
    isempty(peaks) && println("    No significant autocorrelation peaks found.")
    
    # fwd/bwd correlation
    ch_bwd = GaussianFit2D.get_channel(img, "Z"; direction="bwd")
    z_bwd = ch_bwd.data .* 1e9
    fwd_v = vec(ch.data .* 1e9)
    bwd_v = vec(z_bwd)
    @printf("  fwd/bwd correlation: %.3f\n", cor(fwd_v, bwd_v))
end
