#!/usr/bin/env julia
# Diagnostics on 034 image: dimensions, noise, double-tip artifacts.
using GaussianFit2D, STMMolecularFit, Printf, Statistics

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_034.sxm"

img = GaussianFit2D.read_sxm(FILE)
@printf("Image: %d × %d px, range=(%.2f, %.2f) nm\n", img.width, img.height, img.range_nm...)

for ch in img.channels
    z = ch.data
    @printf("  %s (%s): min=%.3f max=%.3f std=%.4f\n", ch.name, ch.direction, minimum(z), maximum(z), std(z))
end

# Check for double-tip: look at Z fwd profile along x
z_fwd = GaussianFit2D.get_channel(img, "Z"; direction="fwd").data
z_bwd = GaussianFit2D.get_channel(img, "Z"; direction="bwd").data

# Correlation fwd vs bwd (should be high for good data)
fwd_flat = vec(z_fwd)
bwd_flat = vec(z_bwd)
corr = cor(fwd_flat, bwd_flat)
@printf("\n  fwd/bwd correlation: %.3f\n", corr)

# Look for double-tip signature: cross-correlation of each row with itself shifted
# A double tip would show a peak at a characteristic shift
println("\n  Checking for double-tip artifact (row autocorrelation shift):")
n_rows, n_cols = size(z_fwd)
shifts = 1:20
max_corr_shift = zeros(Int, n_rows)
for iy in 1:min(n_rows, 100)
    row = z_fwd[iy, :] .- mean(z_fwd[iy, :])
    row_std = std(row)
    row_std < 1e-10 && continue
    best_shift = 0
    best_val = -Inf
    for s in shifts
        if s < length(row)
            c = cor(row[1:end-s], row[(1+s):end])
            if c > best_val
                best_val = c
                best_shift = s
            end
        end
    end
    max_corr_shift[iy] = best_shift
end
shift_counts = sort(collect(Dict(s => count(==(s), max_corr_shift) for s in unique(max_corr_shift))); by=first)
for (s, c) in shift_counts
    @printf("    shift=%2d px: %d rows\n", s, c)
end

# Check pixel size and expected double-tip shift
px_nm = img.range_nm[1] / img.width
@printf("\n  Pixel size: %.4f nm/px\n", px_nm)
@printf("  If double-tip at ~0.5 nm shift: would appear at shift=%.0f px\n", 0.5/px_nm)
