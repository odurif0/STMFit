#!/usr/bin/env julia
# Check which threshold dominates in 1D support detection for 006.
using STMMolecularFit, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_006.sxm"

img = STMMolecularFit.read_sxm(FILE)

for (label, has_frac) in [("old (noise only)", false)]
    scfg = STMMolecularFit.SlideConfig(width_nm=0.30,
        support_noise_k=2.5, support_padding_nm=0.20,
        output_dir="/tmp", no_plot=true)
    
    slide = STMMolecularFit.extract_slide(img, scfg)
    
    println("\n=== $label ===")
    println("  support_length: $(slide.support_length_nm) nm")
    println("  n_points: $(length(slide.x))")
    println("  x range: [$(minimum(slide.x)), $(maximum(slide.x))] nm")
end
