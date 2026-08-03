function _gaussian_blur(z::Matrix{Float64}, sigma_px::Float64)
    sigma_px <= 0.0 && return copy(z)
    radius = max(1, ceil(Int, 3.0 * sigma_px))
    ny, nx = size(z)
    kernel = Float64[exp(-0.5 * (i / sigma_px)^2) for i in -radius:radius]
    kernel ./= sum(kernel)
    nk = length(kernel)
    tmp = similar(z)
    @inbounds for iy in 1:ny, ix in 1:nx
        s = 0.0
        for j in 1:nk
            idx = clamp(ix + j - 1 - radius, 1, nx)
            s += kernel[j] * z[iy, idx]
        end
        tmp[iy, ix] = s
    end
    out = similar(z)
    @inbounds for iy in 1:ny, ix in 1:nx
        s = 0.0
        for j in 1:nk
            idx = clamp(iy + j - 1 - radius, 1, ny)
            s += kernel[j] * tmp[idx, ix]
        end
        out[iy, ix] = s
    end
    return out
end

function _add_noise!(z_fwd::Matrix{Float64}, z_bwd::Matrix{Float64},
                     rng::AbstractRNG, noise_sigma::Float64,
                     corr_frac::Float64)
    noise_sigma <= 0.0 && return
    h, w = size(z_fwd)
    if corr_frac > 0.0
        shared = randn(rng, h, w) .* (noise_sigma * sqrt(corr_frac))
        z_fwd .+= shared
        z_bwd .+= shared
    end
    if corr_frac < 1.0
        z_fwd .+= randn(rng, h, w) .* (noise_sigma * sqrt(1.0 - corr_frac))
        z_bwd .+= randn(rng, h, w) .* (noise_sigma * sqrt(1.0 - corr_frac))
    end
end

function _add_row_offsets!(z::Matrix{Float64}, rng::AbstractRNG, sigma::Float64)
    sigma <= 0.0 && return
    h = size(z, 1)
    offsets = randn(rng, h) .* sigma
    @inbounds for iy in 1:h
        z[iy, :] .+= offsets[iy]
    end
end

function _inverse_affine_coords(nuisance::ViewNuisanceTruth, x::Real, y::Real)
    xs = Float64(x) / nuisance.scale_x
    ys = (Float64(y) - nuisance.shear * Float64(x)) / nuisance.scale_y
    return xs, ys
end
