#!/usr/bin/env julia
using STMMolecularFit, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_006.sxm"
img = STMMolecularFit.read_sxm(FILE)

scfg = STMMolecularFit.SlideConfig(width_nm=0.30, support_noise_k=2.5, support_padding_nm=0.20,
    output_dir="/tmp", no_plot=true)
slide = STMMolecularFit.extract_slide(img, scfg)
@printf("support_length: %.2f nm, n_points: %d\n", slide.support_length_nm, length(slide.x))
@printf("x range: [%.2f, %.2f] nm\n", minimum(slide.x), maximum(slide.x))
