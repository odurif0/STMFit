# Joint shape fit of one molecule on its registered repeat scans (campaign J,
# stage 2; docs/src/journal.md, 2026-10-26).
#
# The scans of a molecule track are mapped onto a reference scan by rigid
# transforms between the physical lobes they share, resampled on one grid of the
# molecule frame (t along the chain, u across), normalised and averaged. The
# stack is fitted with one elliptical Gaussian per unit, each with its own
# position, amplitude, widths and orientation, plus a background plane, by
# Levenberg–Marquardt with an analytic Jacobian over every pixel of the window.
# Label-free: geometry, physical lobes and raw topography only.
module JointShape

using LinearAlgebra, Statistics

export rigid_fit, apply_rigid, invert_rigid, MoleculeFrame, molecule_frame, to_frame, from_frame,
       stack_scans, LobeFit, fit_lobes, lobe_model, SHAPE_BOUNDS

"""
    rigid_fit(src, dst) -> (θ, tx, ty)

Rotation `θ` and translation that map the points `src` onto `dst` in the least
squares sense (two-dimensional Procrustes, no scaling).
"""
function rigid_fit(src::AbstractVector, dst::AbstractVector)
    n = length(src)
    n == length(dst) && n >= 2 || throw(ArgumentError("rigid_fit needs two lists of at least two matched points"))
    cs = (mean(p[1] for p in src), mean(p[2] for p in src)); cd = (mean(p[1] for p in dst), mean(p[2] for p in dst))
    num = 0.0; den = 0.0
    for (s, d) in zip(src, dst)
        ax, ay = s[1] - cs[1], s[2] - cs[2]; bx, by = d[1] - cd[1], d[2] - cd[2]
        num += ax * by - ay * bx; den += ax * bx + ay * by
    end
    θ = atan(num, den); c, s = cos(θ), sin(θ)
    return θ, cd[1] - (c * cs[1] - s * cs[2]), cd[2] - (s * cs[1] + c * cs[2])
end

"Point `p` moved by the rigid transform `T = (θ, tx, ty)`."
apply_rigid(T, p) = (cos(T[1]) * p[1] - sin(T[1]) * p[2] + T[2], sin(T[1]) * p[1] + cos(T[1]) * p[2] + T[3])

"Inverse of a rigid transform."
function invert_rigid(T)
    θ, tx, ty = T; c, s = cos(θ), sin(θ)
    return (-θ, -(c * tx + s * ty), -(-s * tx + c * ty))
end

"Molecule frame: origin, unit axis along the chain and unit perpendicular."
struct MoleculeFrame
    origin::NTuple{2,Float64}
    axis::NTuple{2,Float64}
    perp::NTuple{2,Float64}
end

"""
    molecule_frame(points)

Frame of a chain from its lobe positions in chain order: centroid, principal
axis, oriented from the first to the last point.
"""
function molecule_frame(points::AbstractVector)
    length(points) >= 2 || throw(ArgumentError("a frame needs at least two points"))
    cx, cy = mean(p[1] for p in points), mean(p[2] for p in points)
    sxx = sum((p[1] - cx)^2 for p in points); syy = sum((p[2] - cy)^2 for p in points)
    sxy = sum((p[1] - cx) * (p[2] - cy) for p in points)
    v = eigen(Symmetric([sxx sxy; sxy syy])).vectors[:, 2]
    ax = (v[1], v[2])
    (points[end][1] - points[1][1]) * ax[1] + (points[end][2] - points[1][2]) * ax[2] < 0 && (ax = (-ax[1], -ax[2]))
    return MoleculeFrame((cx, cy), ax, (-ax[2], ax[1]))
end

"Coordinates (t, u) of an absolute point in a molecule frame."
to_frame(F::MoleculeFrame, p) = ((p[1] - F.origin[1]) * F.axis[1] + (p[2] - F.origin[2]) * F.axis[2],
                                 (p[1] - F.origin[1]) * F.perp[1] + (p[2] - F.origin[2]) * F.perp[2])

"Absolute point of the frame coordinates (t, u)."
from_frame(F::MoleculeFrame, t, u) = (F.origin[1] + t * F.axis[1] + u * F.perp[1], F.origin[2] + t * F.axis[2] + u * F.perp[2])

"""
    stack_scans(samplers, ts, us; chain_t, bg_u=0.55, chain_u=0.15, min_share=0.5) -> (image, counts, used)

Mean of the normalised scans on the grid `us × ts` (rows u, columns t). Each
sampler maps (t, u) to a topography value (NaN outside the scan). A scan is
normalised by minus the median of its pixels with |u| > `bg_u` and divided by
the 95th percentile of its chain pixels (|u| < `chain_u`, t within `chain_t`);
scans without at least 20 pixels of each are left out. Pixels seen by fewer than
`min_share` of the used scans are NaN.
"""
function stack_scans(samplers, ts::AbstractVector, us::AbstractVector; chain_t::NTuple{2,Real},
                     bg_u::Real=0.55, chain_u::Real=0.15, min_share::Real=0.5)
    acc = zeros(length(us), length(ts)); cnt = zeros(Int, length(us), length(ts)); used = 0
    for f in samplers
        V = [f(t, u) for u in us, t in ts]
        bg = [V[i, j] for i in eachindex(us), j in eachindex(ts) if abs(us[i]) > bg_u && isfinite(V[i, j])]
        ch = [V[i, j] for i in eachindex(us), j in eachindex(ts)
              if abs(us[i]) < chain_u && chain_t[1] <= ts[j] <= chain_t[2] && isfinite(V[i, j])]
        (length(bg) >= 20 && length(ch) >= 20) || continue
        b = median(bg); s = quantile(ch .- b, 0.95)
        s > 0 || continue
        for k in eachindex(V)
            isfinite(V[k]) && (acc[k] += (V[k] - b) / s; cnt[k] += 1)
        end
        used += 1
    end
    need = max(1, ceil(Int, min_share * used))
    img = [cnt[k] >= need && used > 0 ? acc[k] / cnt[k] : NaN for k in CartesianIndices(acc)]
    return img, cnt, used
end

"Fitted unit: position in the molecule frame (nm), amplitude, widths along and across (nm), orientation (rad)."
struct LobeFit
    t::Float64
    u::Float64
    amplitude::Float64
    sigma_par::Float64
    sigma_perp::Float64
    theta::Float64
end

"Bounds of the joint shape fit: widths (nm), orientation (rad) and displacement from the start (nm)."
const SHAPE_BOUNDS = (sigma=(0.10, 0.60), theta=0.6, shift=0.3)

_sig(x) = 1 / (1 + exp(-x))

function _decode(p, init)
    lo, hi = SHAPE_BOUNDS.sigma
    lobes = LobeFit[]
    for (k, (t0, u0)) in enumerate(init)
        j = 3 + 6(k - 1)
        push!(lobes, LobeFit(t0 + SHAPE_BOUNDS.shift * tanh(p[j+1]), u0 + SHAPE_BOUNDS.shift * tanh(p[j+2]),
                             exp(p[j+3]), lo + (hi - lo) * _sig(p[j+4]), lo + (hi - lo) * _sig(p[j+5]),
                             SHAPE_BOUNDS.theta * tanh(p[j+6])))
    end
    return lobes, (p[1], p[2], p[3])
end

"Model values at the points (T, U) for fitted units and a background plane (b0, bt, bu)."
function lobe_model(T::AbstractVector, U::AbstractVector, lobes::AbstractVector{LobeFit}, bg)
    m = bg[1] .+ bg[2] .* T .+ bg[3] .* U
    for l in lobes
        c, s = cos(l.theta), sin(l.theta)
        @inbounds for i in eachindex(T)
            dt = T[i] - l.t; du = U[i] - l.u
            a = c * dt + s * du; b = -s * dt + c * du
            m[i] += l.amplitude * exp(-0.5 * ((a / l.sigma_par)^2 + (b / l.sigma_perp)^2))
        end
    end
    return m
end

function _jacobian!(J, T, U, p, init)
    lobes, _ = _decode(p, init)
    lo, hi = SHAPE_BOUNDS.sigma
    J[:, 1] .= 1.0; J[:, 2] .= T; J[:, 3] .= U
    for (k, l) in enumerate(lobes)
        j = 3 + 6(k - 1)
        c, s = cos(l.theta), sin(l.theta)
        ip2 = 1 / l.sigma_par^2; in2 = 1 / l.sigma_perp^2
        dtk = SHAPE_BOUNDS.shift * (1 - tanh(p[j+1])^2); duk = SHAPE_BOUNDS.shift * (1 - tanh(p[j+2])^2)
        sp = _sig(p[j+4]); sn = _sig(p[j+5])
        dsp = (hi - lo) * sp * (1 - sp); dsn = (hi - lo) * sn * (1 - sn)
        dth = SHAPE_BOUNDS.theta * (1 - tanh(p[j+6])^2)
        @inbounds for i in eachindex(T)
            dt = T[i] - l.t; du = U[i] - l.u
            a = c * dt + s * du; b = -s * dt + c * du
            g = l.amplitude * exp(-0.5 * (a^2 * ip2 + b^2 * in2))
            J[i, j+1] = g * (a * c * ip2 - b * s * in2) * dtk
            J[i, j+2] = g * (a * s * ip2 + b * c * in2) * duk
            J[i, j+3] = g                                        # d/dlog A
            J[i, j+4] = g * a^2 * ip2 / l.sigma_par * dsp
            J[i, j+5] = g * b^2 * in2 / l.sigma_perp * dsn
            J[i, j+6] = -g * a * b * (ip2 - in2) * dth
        end
    end
    return J
end

"""
    fit_lobes(ts, us, img, init; maxiter=300, tol=1e-10) -> (lobes, background, info)

Fit one elliptical Gaussian per start position `init` (frame coordinates) and a
background plane to the finite pixels of `img` (rows `us`, columns `ts`), by
Levenberg–Marquardt. Widths, orientation and displacement are bounded by
`SHAPE_BOUNDS`. `info` gives the iterations, convergence, and the root mean
square residual relative to the stack's 95th percentile.
"""
function fit_lobes(ts::AbstractVector, us::AbstractVector, img::AbstractMatrix, init::AbstractVector;
                   maxiter::Integer=300, tol::Real=1e-10)
    size(img) == (length(us), length(ts)) || throw(ArgumentError("image size does not match the grid"))
    T = Float64[]; U = Float64[]; Y = Float64[]
    for (i, u) in enumerate(us), (j, t) in enumerate(ts)
        isfinite(img[i, j]) && (push!(T, t); push!(U, u); push!(Y, img[i, j]))
    end
    n = length(init); np = 3 + 6n
    length(Y) > 2np || error("too few pixels for $n units")
    b0 = median(Y[abs.(U) .> 0.55]); lo, hi = SHAPE_BOUNDS.sigma
    p = zeros(np); p[1] = b0
    for (k, (t0, u0)) in enumerate(init)
        j = 3 + 6(k - 1)
        near = argmin((T .- t0) .^ 2 .+ (U .- u0) .^ 2)
        p[j+3] = log(max(Y[near] - b0, 0.05))
        p[j+4] = p[j+5] = log((0.3 - lo) / (hi - 0.3))         # widths start at 0.3 nm
    end
    J = zeros(length(Y), np)
    resid(q) = Y .- lobe_model(T, U, _decode(q, init)...)
    r = resid(p); cost = sum(abs2, r); λ = 1e-3; iters = 0; converged = false
    while iters < maxiter
        iters += 1
        _jacobian!(J, T, U, p, init)
        A = J' * J; g = J' * r
        improved = false
        for _ in 1:20
            δ = try
                (A + λ * Diagonal(diag(A) .+ 1e-12)) \ g
            catch
                λ *= 10; continue
            end
            q = p .+ δ; rq = resid(q); cq = sum(abs2, rq)
            if cq < cost
                rel = (cost - cq) / max(cost, eps())
                p, r, cost = q, rq, cq; λ = max(λ / 3, 1e-9); improved = true
                rel < tol && (converged = true)
                break
            end
            λ *= 4
        end
        (!improved || converged) && (converged = converged || !improved; break)
    end
    lobes, bg = _decode(p, init)
    scale = quantile(Y, 0.95) - b0
    return lobes, bg, (iterations=iters, converged=converged, rms_rel=sqrt(cost / length(Y)) / max(scale, eps()))
end

end
