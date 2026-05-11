#!/usr/bin/env julia
# Extract a 1D slide/profile from a 2D Nanonis SXM image.
# Uses STMMolecularFit pipeline.
#
# Output: slide_profile.txt (compatible with GaussianFit1D.jl)
#
# Usage:
#   julia --project scripts/extract_slide_1d.jl path/to/image.sxm [options]

using STMMolecularFit
using Printf

function main()
    filepath = length(ARGS) >= 1 ? ARGS[1] : error("Usage: julia extract_slide_1d.jl <filepath>")
    
    cfg = SlideConfig(
        channel="Z", direction="fwd",
        width_nm=0.30,
        support_threshold_fraction=0.20,
        support_noise_k=2.5,
        support_padding_nm=0.20,
        output_dir="results/slide_1d",
        no_plot=false,
    )
    
    img = read_sxm(filepath)
    slide = extract_slide(img, cfg)
    files = write_slide_outputs(slide, cfg)
    
    println("slide_profile: $(files.profile)")
    println("slide_full_profile: $(files.full)")
    println("slide_metadata: $(files.metadata)")
    files.plot !== nothing && println("slide_plot: $(files.plot)")
    @printf("support_length_nm: %.6f\n", slide.support_length_nm)
end

main()
