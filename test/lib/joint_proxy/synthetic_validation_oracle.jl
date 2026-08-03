calibration_observation_mode() = "synthetic_physical_proxy"

function swap_registry(reg::ProxyEnsemble)
    swapped_entries = ProxyEntry[]
    for entry in reg.entries
        templates = [ProxyTemplate(t.type == 0 ? 1 : 0, t.parity, t.mirror, copy(t.pixels)) for t in entry.templates]
        push!(swapped_entries, ProxyEntry(entry.source, templates))
    end
    return ProxyEnsemble(swapped_entries, reg.grid_half_nm, reg.grid_step_nm, reg.grid_n, reg.npix, reg.payload_sha256, copy(reg.audit))
end

oracle_geometry_adapter_name() = "oracle_geometry_truth_v2"

function _view_data(case::SyntheticCase, direction::String)
    ch = STMSXMIO.get_channel(case.img, "Z"; direction=direction)
    xs = collect(range(0.0, case.img.range_nm[1]; length=case.img.width))
    ys = collect(range(0.0, case.img.range_nm[2]; length=case.img.height))
    return ViewData(direction, xs, ys, ch.data, std(vec(ch.data)), true)
end

function _plane_subtract(z::AbstractMatrix{<:Real}, xs::Vector{Float64}, ys::Vector{Float64})
    X = Float64[]; Y = Float64[]; Z = Float64[]
    for (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        v = Float64(z[iy, ix])
        isfinite(v) || continue
        push!(X, x); push!(Y, y); push!(Z, v)
    end
    isempty(Z) && return Float64.(z)
    A = hcat(ones(length(Z)), X, Y)
    β = A \ Z
    out = similar(Float64.(z))
    for (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        v = Float64(z[iy, ix])
        out[iy, ix] = isfinite(v) ? v - (β[1] + β[2] * x + β[3] * y) : NaN
    end
    return out
end

function _row_median_flatten(z::AbstractMatrix{<:Real})
    out = Matrix{Float64}(undef, size(z, 1), size(z, 2))
    for iy in axes(z, 1)
        row = [Float64(v) for v in z[iy, :] if isfinite(v)]
        row_med = isempty(row) ? 0.0 : median(row)
        for ix in axes(z, 2)
            v = Float64(z[iy, ix])
            out[iy, ix] = isfinite(v) ? v - row_med : NaN
        end
    end
    return out
end

function _preprocess_view(view::ViewData)
    z = _plane_subtract(view.z, view.xs, view.ys)
    z = _row_median_flatten(z)
    return ViewData(view.label, copy(view.xs), copy(view.ys), z, std(vec(z)), view.present)
end

function _axis(view::ViewData)
    w = max.(view.z .- median(view.z), 0.0)
    s = sum(w)
    s > 0 || return (; cx=mean(view.xs), cy=mean(view.ys), ax=1.0, ay=0.0, ux=0.0, uy=1.0)
    xs = repeat(view.xs', length(view.ys), 1)
    ys = repeat(view.ys, 1, length(view.xs))
    cx = sum(xs .* w) / s; cy = sum(ys .* w) / s
    dx = xs .- cx; dy = ys .- cy
    cxx = sum(w .* dx .* dx) / s; cyy = sum(w .* dy .* dy) / s; cxy = sum(w .* dx .* dy) / s
    ang = 0.5 * atan(2 * cxy, cxx - cyy)
    ax, ay = cos(ang), sin(ang)
    return (; cx, cy, ax, ay, ux=-ay, uy=ax)
end

function _sample(view::ViewData, x::Float64, y::Float64)
    xs, ys, z = view.xs, view.ys, view.z
    (x < first(xs) || x > last(xs) || y < first(ys) || y > last(ys)) && return NaN
    ix = clamp(searchsortedlast(xs, x), 1, length(xs) - 1); iy = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
    x1, x2 = xs[ix], xs[ix + 1]; y1, y2 = ys[iy], ys[iy + 1]
    wx = (x - x1) / max(x2 - x1, eps(Float64)); wy = (y - y1) / max(y2 - y1, eps(Float64))
    v11, v21, v12, v22 = z[iy, ix], z[iy, ix + 1], z[iy + 1, ix], z[iy + 1, ix + 1]
    any(!isfinite, (v11, v21, v12, v22)) && return NaN
    return (1 - wx) * (1 - wy) * v11 + wx * (1 - wy) * v21 + (1 - wx) * wy * v12 + wx * wy * v22
end

function _exact_gaussian_backbone(case::SyntheticCase, x::Float64, y::Float64)
    acc = 0.0
    @inbounds for k in eachindex(case.truth.lobe_x_nm)
        dx = x - case.truth.lobe_x_nm[k]
        dy = y - case.truth.lobe_y_nm[k]
        acc += case.truth.amplitudes[k] * exp(-0.5 * ((dx * dx) / (case.truth.sigma_par[k]^2) + (dy * dy) / (case.truth.sigma_perp[k]^2)))
    end
    return acc
end

_view_nuisance(case::SyntheticCase, direction::String) = lowercase(direction) == "fwd" ? case.truth.fwd_nuisance : case.truth.bwd_nuisance

function _molecular_to_image(nuisance, x_mol::Float64, y_mol::Float64)
    x_img = nuisance.scale_x * x_mol
    y_img = nuisance.scale_y * y_mol + nuisance.shear * x_img
    return x_img, y_img
end

function _pixel_scale_nm(view::ViewData)
    dx = length(view.xs) > 1 ? mean(diff(view.xs)) : 1.0
    dy = length(view.ys) > 1 ? mean(diff(view.ys)) : 1.0
    return max(0.5 * (abs(dx) + abs(dy)), eps(Float64))
end

function _blur_broadened_sigma(sigma_nm::Float64, blur_sigma_px::Float64, view::ViewData)
    blur_nm = blur_sigma_px * _pixel_scale_nm(view)
    return sqrt(sigma_nm^2 + blur_nm^2)
end

function _oracle_gaussian_backbone(case::SyntheticCase, view::ViewData, nuisance, x_img::Float64, y_img::Float64)
    x_mol = x_img / nuisance.scale_x
    y_mol = (y_img - nuisance.shear * x_img) / nuisance.scale_y
    acc = 0.0
    @inbounds for k in eachindex(case.truth.lobe_x_nm)
        dx = x_mol - case.truth.lobe_x_nm[k]
        dy = y_mol - case.truth.lobe_y_nm[k]
        σp = _blur_broadened_sigma(case.truth.sigma_par[k], nuisance.blur_sigma_px, view)
        σu = _blur_broadened_sigma(case.truth.sigma_perp[k], nuisance.blur_sigma_px, view)
        amp = case.truth.amplitudes[k] * (case.truth.sigma_par[k] * case.truth.sigma_perp[k]) / (σp * σu)
        acc += amp * exp(-0.5 * ((dx * dx) / (σp * σp) + (dy * dy) / (σu * σu)))
    end
    return acc
end

function _plane_fit_residual(raw_patch::Vector{Float64}, pts)
    finite_idx = Int[]
    X = Float64[]
    Y = Float64[]
    for i in eachindex(raw_patch)
        v = raw_patch[i]
        isfinite(v) || continue
        push!(finite_idx, i)
        push!(X, Float64(getproperty(pts[i], :t)))
        push!(Y, Float64(getproperty(pts[i], :u)))
    end
    length(finite_idx) < 3 && return [isfinite(v) ? v - _finite_mean(raw_patch) : 0.0 for v in raw_patch]
    A = hcat(ones(length(finite_idx)), X, Y)
    β = A \ raw_patch[finite_idx]
    out = Vector{Float64}(undef, length(raw_patch))
    for i in eachindex(raw_patch)
        if isfinite(raw_patch[i])
            t = Float64(getproperty(pts[i], :t))
            u = Float64(getproperty(pts[i], :u))
            out[i] = raw_patch[i] - (β[1] + β[2] * t + β[3] * u)
        else
            out[i] = 0.0
        end
    end
    return out
end

function _exact_backbone_image(case::SyntheticCase, view::ViewData)
    out = Matrix{Float64}(undef, length(view.ys), length(view.xs))
    @inbounds for (iy, y) in enumerate(view.ys), (ix, x) in enumerate(view.xs)
        out[iy, ix] = _exact_gaussian_backbone(case, x, y)
    end
    return out
end

function _sample_matrix(xs::Vector{Float64}, ys::Vector{Float64}, z::Matrix{Float64}, x::Float64, y::Float64)
    (x < first(xs) || x > last(xs) || y < first(ys) || y > last(ys)) && return NaN
    ix = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    iy = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
    x1, x2 = xs[ix], xs[ix + 1]; y1, y2 = ys[iy], ys[iy + 1]
    wx = (x - x1) / max(x2 - x1, eps(Float64)); wy = (y - y1) / max(y2 - y1, eps(Float64))
    v11, v21, v12, v22 = z[iy, ix], z[iy, ix + 1], z[iy + 1, ix], z[iy + 1, ix + 1]
    any(!isfinite, (v11, v21, v12, v22)) && return NaN
    return (1 - wx) * (1 - wy) * v11 + wx * (1 - wy) * v21 + (1 - wx) * wy * v12 + wx * wy * v22
end

function _exact_center_patches(case::SyntheticCase, view_fwd::ViewData, view_bwd::ViewData,
                               reg::ProxyEnsemble; patch_half_nm::Float64=reg.grid_half_nm,
                               patch_step_nm::Float64=reg.grid_step_nm)
    offsets = collect(-patch_half_nm:patch_step_nm:patch_half_nm)
    patch_pts = NamedTuple{(:t,:u),Tuple{Float64,Float64}}[(; t=t, u=u) for u in offsets, t in offsets]
    θ = deg2rad(case.truth.orientation_deg)
    cos_th = cos(θ)
    sin_th = sin(θ)
    lobes = TypePosteriorLobeEvidence[]
    for k in eachindex(case.truth.lobe_x_nm)
        x0 = case.truth.lobe_x_nm[k]
        y0 = case.truth.lobe_y_nm[k]
        fwd_patch = Float64[]
        bwd_patch = Float64[]
        nf = _view_nuisance(case, "fwd")
        nb = _view_nuisance(case, "bwd")
        for u in offsets, t in offsets
            x_mol = x0 + t * cos_th - u * sin_th
            y_mol = y0 + t * sin_th + u * cos_th
            x_img_f, y_img_f = _molecular_to_image(nf, x_mol, y_mol)
            x_img_b, y_img_b = _molecular_to_image(nb, x_mol, y_mol)
            raw_f = _sample(view_fwd, x_img_f, y_img_f)
            raw_b = _sample(view_bwd, x_img_b, y_img_b)
            back_f = _oracle_gaussian_backbone(case, view_fwd, nf, x_img_f, y_img_f)
            back_b = _oracle_gaussian_backbone(case, view_bwd, nb, x_img_b, y_img_b)
            push!(fwd_patch, isfinite(raw_f) ? raw_f - back_f : NaN)
            push!(bwd_patch, isfinite(raw_b) ? raw_b - back_b : NaN)
        end
        push!(lobes, TypePosteriorLobeEvidence(Dict("fwd" => fwd_patch, "bwd" => bwd_patch)))
    end
    return lobes
end

function oracle_geometry_type_result(case::SyntheticCase, reg::ProxyEnsemble; swap_type_mapping::Bool=false)
    use_reg = swap_type_mapping ? swap_registry(reg) : reg
    fwd = _view_data(case, "fwd")
    bwd = _view_data(case, "bwd")
    nf = _view_nuisance(case, "fwd")
    nb = _view_nuisance(case, "bwd")
    lobes = _exact_center_patches(case, fwd, bwd, reg)
    fwd_back = Matrix{Float64}(undef, length(fwd.ys), length(fwd.xs))
    bwd_back = Matrix{Float64}(undef, length(bwd.ys), length(bwd.xs))
    @inbounds for (iy, y_img) in enumerate(fwd.ys), (ix, x_img) in enumerate(fwd.xs)
        fwd_back[iy, ix] = _oracle_gaussian_backbone(case, fwd, nf, x_img, y_img)
    end
    @inbounds for (iy, y_img) in enumerate(bwd.ys), (ix, x_img) in enumerate(bwd.xs)
        bwd_back[iy, ix] = _oracle_gaussian_backbone(case, bwd, nb, x_img, y_img)
    end
    fwd_sub = fwd.z .- fwd_back
    bwd_sub = bwd.z .- bwd_back
    rho = cor(vec(fwd_sub), vec(bwd_sub))
    return _normalize_type_result(infer_type_posterior(lobes, use_reg; rho=rho, effective_factor=effective_view_factor(rho)))
end

