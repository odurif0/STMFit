# Patch extraction helpers for Todo 2.

# (t,u) offset grid for patches: for u in coords, for t in coords (row-major,
# matching test/extract_lobe_patches_bwd.jl ordering).
function _patch_offset_grid(half_nm::Float64, step_nm::Float64)
    half_nm > 0 || throw(ArgumentError("patch_half_nm must be positive"))
    step_nm > 0 || throw(ArgumentError("patch_step_nm must be positive"))
    coords = collect(-half_nm:step_nm:half_nm)
    isempty(coords) && (coords = [0.0])
    return [(t=t, u=u) for u in coords for t in coords]
end

# Bilinear interpolation at (x, y) on grid (xs, ys) with image z[y, x].
# Returns NaN out of bounds (matches test/extract_lobe_patches_bwd.jl).
function _interp(xs::Vector{Float64}, ys::Vector{Float64},
                 z::Matrix{Float64}, x::Float64, y::Float64)
    ix = searchsortedlast(xs, x)
    iy = searchsortedlast(ys, y)
    (ix < 1 || iy < 1 || ix >= length(xs) || iy >= length(ys)) && return NaN
    x1, x2 = xs[ix], xs[ix + 1]
    y1, y2 = ys[iy], ys[iy + 1]
    dx = x2 - x1
    dy = y2 - y1
    (dx == 0 || dy == 0) && return z[iy, ix]
    tx = (x - x1) / dx
    ty = (y - y1) / dy
    z11 = z[iy, ix];     z21 = z[iy, ix + 1]
    z12 = z[iy + 1, ix]; z22 = z[iy + 1, ix + 1]
    return (1 - tx) * (1 - ty) * z11 + tx * (1 - ty) * z21 +
           (1 - tx) * ty * z12 + tx * ty * z22
end

# Decode the candidate's own geometry and sample aligned residual + raw patches
# at each lobe center, per present view.
function _decode_and_patch(r::ChainModelResult, ccfg::ChainSweepConfig, axisctx,
                           patch_tu, ax, ay, px, py,
                           present_views::Vector{ViewData},
                           resid_imgs::Dict{String,Matrix{Float64}},
                           model_img::Matrix{Float64})
    b0, feats, ts, us, spars, sperps = GaussianFit2D._decode_chain(
        r.params, r.n, axisctx, ccfg; amp_min=r.amp_min, amp_range=r.amp_range)
    lobes = LobePatch[]
    for k in eachindex(feats)
        f = feats[k]
        cx, cy = f.x_nm, f.y_nm
        res_patches = Dict{String,Vector{Float64}}()
        raw_patches = Dict{String,Vector{Float64}}()
        for v in present_views
            rimg = get(resid_imgs, v.label, nothing)
            raw_img = v.z .- model_img
            rp = rimg === nothing ? fill(NaN, length(patch_tu)) :
                 [_interp(v.xs, v.ys, rimg, cx + tu.t * ax + tu.u * px,
                          cy + tu.t * ay + tu.u * py) for tu in patch_tu]
            gp = [_interp(v.xs, v.ys, raw_img, cx + tu.t * ax + tu.u * px,
                          cy + tu.t * ay + tu.u * py) for tu in patch_tu]
            res_patches[v.label] = rp
            raw_patches[v.label] = gp
        end
        push!(lobes, LobePatch(k, cx, cy, ts[k], us[k], f.amplitude,
                               spars[k], sperps[k], f.skew_ratio,
                               res_patches, raw_patches))
    end
    return lobes
end
