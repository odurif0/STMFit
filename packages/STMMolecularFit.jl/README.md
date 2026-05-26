# STMMolecularFit

Installable Julia package for STM image manipulation.

Current features:

- read Nanonis `.sxm` files;
- list and access channels/directions;
- align backward scans spatially with forward scans;
- STM preprocessing: plane flattening, row flattening, smoothing;
- molecular ROI detection;
- extract a 1D slide/profile from a 2D STM image using `extract_slide`;
- fit the extracted slide with `GaussianFit1D` using `fit_slide`;
- run the complete workflow with `extract_and_fit_slide`.

## Installation from Git/local path

```julia
using Pkg
Pkg.develop(path="/home/durif/Git/GaussianFit1D.jl")
Pkg.develop(path="/home/durif/Git/STMMolecularFit")
```

Then:

```julia
using STMMolecularFit
```

## Basic usage

```julia
using STMMolecularFit

img = read_sxm("path/to/your_image.sxm")
println(channel_names(img))

cfg = SlideConfig(
    channel="Z",
    direction="fwd",
    width_nm=0.30,
    support_noise_k=2.5,
    support_padding_nm=0.25,
    output_dir="results/slide"
)

slide = extract_slide(img, cfg)
files = write_slide_outputs(slide, cfg)

println(files.profile)  # two-column text file: distance_nm, intensity
```

## Complete extraction + fit workflow

```julia
using STMMolecularFit

slide_cfg = SlideConfig(
    channel="Z",
    direction="fwd",
    width_nm=0.30,
    support_noise_k=2.5,
    support_padding_nm=0.25,
    output_dir="results/slide_fit/slide"
)

fit_cfg = FitSlideConfig(
    min_spacing=0.4,
    max_spacing=0.675,
    fwhm_min=0.45,
    fwhm_max=1.2,
    global_maxtime=8.0,
    output_dir="results/slide_fit"
)

result = extract_and_fit_slide(
    "path/to/your_image.sxm";
    slide_config=slide_cfg,
    fit_config=fit_cfg,
)

println(result.fit.best_model.n_peaks)
println(result.fit.best_model.bic)
```

This writes:

```text
slide_profile.txt
slide_full_profile.txt
slide_metadata.tsv
slide_profile.png
model_selection.tsv
best_model.tsv
best_peaks.tsv
fit_1d/*_results.txt
fit_1d/*_n<N>.png
```

## Lower-level fit API

```julia
img = read_sxm("/path/to/file.sxm")
slide = extract_slide(img, SlideConfig())
fit = fit_slide(slide, FitSlideConfig())
write_fit_outputs(fit, FitSlideConfig(output_dir="results/fit"))
```

## Output profile

The exported `slide_profile.txt` is compatible with 1D fitting tools such as `GaussianFit1D`.
