module Fit2D

# Ordered 2D Gaussian chain fitting for STM molecular patterns.

export fit_2d_chain, Fit2DResult

using ..STMMolecularFit
using GaussianFit2D

Base.@kwdef struct Fit2DResult
    chain_best_n::Int = 0
    chain_best_valid::Bool = false
    chain_best_bic::Float64 = NaN
end

function fit_2d_chain(filepath::AbstractString;
                      ccfg::GaussianFit2D.ChainSweepConfig=GaussianFit2D.ChainSweepConfig())
    img = GaussianFit2D.read_sxm(String(filepath))
    pcfg = GaussianFit2D.PatternConfig(
        filepath=String(filepath), channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1)
    _, best, _ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
    return Fit2DResult(best.n, best.valid, best.bic)
end

fit_2d_chain(img::SXMImage; ccfg=GaussianFit2D.ChainSweepConfig()) =
    fit_2d_chain(img.filepath; ccfg=ccfg)

# Deprecated alias
fit_2d_image(args...; kwargs...) = fit_2d_chain(args...; kwargs...)

end
