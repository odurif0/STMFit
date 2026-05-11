# Multi-Gaussian Fit GUI

Web GUI for multi-Gaussian fitting of 1D STM/AFM line profiles, built with Stipple.jl + PlotlyJS.

This repository contains **only the GUI layer**. The fitting core is in the separate
[`GaussianFit1D.jl`](https://github.com/od/GaussianFit1D.jl) package.

## Architecture

```text
GaussianFit1D.jl (core package)
  ├── types.jl
  ├── core.jl       ← fitting engine, BIC model selection
  ├── plot.jl       ← static plots (Plots.jl)
  └── cli.jl        ← command-line interface

STMMolecularFitGUI.jl (this repo)
  └── src/STMMolecularFitGUI.jl  ← web GUI (Stipple + PlotlyJS)
```

## Installation

Install the core package first:

```bash
cd /path/to/GaussianFit1D.jl
# Already a working Julia package
```

Then install this GUI:

```bash
cd /path/to/STMMolecularFitGUI.jl
julia --project=. -e 'using Pkg; Pkg.develop(path="/path/to/GaussianFit1D.jl"); Pkg.instantiate()'
```

## Quick start

### Linux / macOS

```bash
bash launch.sh
```

### Windows

Double-click `launch.bat`.

### Manual

```bash
julia --project=. app.jl
```

The GUI opens at `http://localhost:8888`.

## Programmatic use

```julia
using Pkg
Pkg.develop(path="/path/to/GaussianFit1D.jl")
Pkg.develop(path="/path/to/STMMolecularFitGUI.jl")

using STMMolecularFitGUI
STMMolecularFitGUI.run(port=8888)
```
