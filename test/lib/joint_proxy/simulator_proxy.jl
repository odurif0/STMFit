function default_proxy_ensemble()
    sites0 = ProxySite[
        ProxySite(0.000, 0.000, 1.00, 0.150, 0.115),
        ProxySite(-0.115, -0.020, 0.42, 0.090, 0.090),
        ProxySite(0.115, 0.020, 0.42, 0.090, 0.090),
        ProxySite(-0.060, 0.120, 0.35, 0.080, 0.080),
        ProxySite(0.075, -0.155, 0.28, 0.090, 0.090),
        ProxySite(-0.090, 0.230, 0.28, 0.075, 0.080),
    ]
    sites1 = ProxySite[
        ProxySite(0.000, 0.000, 1.00, 0.150, 0.115),
        ProxySite(-0.115, -0.020, 0.42, 0.090, 0.090),
        ProxySite(0.115, 0.020, 0.42, 0.090, 0.090),
        ProxySite(-0.060, 0.120, 0.35, 0.080, 0.080),
        ProxySite(0.075, -0.155, 0.28, 0.090, 0.090),
        ProxySite(-0.090, 0.220, 0.22, 0.075, 0.080),
        ProxySite(-0.090, 0.335, 0.48, 0.080, 0.075),
        ProxySite(-0.010, 0.385, 0.42, 0.075, 0.075),
        ProxySite(-0.170, 0.440, 0.36, 0.095, 0.085),
    ]
    return ProxyEnsemble(sites0, sites1)
end

function load_proxy_ensemble(path::String)
    isfile(path) || error("Proxy sites file not found: $path")
    sites0 = ProxySite[]
    sites1 = ProxySite[]
    for line in readlines(path)
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        parts = split(s, '\t')
        length(parts) >= 7 || continue
        tryparse(Int, strip(parts[1])) === nothing && continue
        typ = parse(Int, strip(parts[1]))
        site = ProxySite(
            parse(Float64, strip(parts[3])),
            parse(Float64, strip(parts[4])),
            parse(Float64, strip(parts[5])),
            max(parse(Float64, strip(parts[6])), eps(Float64)),
            max(parse(Float64, strip(parts[7])), eps(Float64)),
        )
        typ == 0 && push!(sites0, site)
        typ == 1 && push!(sites1, site)
        (typ in (0, 1)) || error("Invalid type $typ in proxy sites file: $path")
    end
    isempty(sites0) && error("No GlcN (type=0) proxy sites in $path")
    isempty(sites1) && error("No GlcNAc (type=1) proxy sites in $path")
    return ProxyEnsemble(sites0, sites1)
end

function validate_config(cfg::SimulatorConfig, ensemble::ProxyEnsemble)
    cfg.n_min >= 1 || throw(ArgumentError("n_min must be >= 1, got $(cfg.n_min)"))
    cfg.n_max >= cfg.n_min || throw(ArgumentError("n_max ($(cfg.n_max)) must be >= n_min ($(cfg.n_min))"))
    cfg.width >= 1 || throw(ArgumentError("width must be >= 1, got $(cfg.width)"))
    cfg.height >= 1 || throw(ArgumentError("height must be >= 1, got $(cfg.height)"))
    cfg.range_nm[1] > 0 || throw(ArgumentError("range_nm[1] must be positive, got $(cfg.range_nm[1])"))
    cfg.range_nm[2] > 0 || throw(ArgumentError("range_nm[2] must be positive, got $(cfg.range_nm[2])"))
    cfg.spacing_nm > 0 || throw(ArgumentError("spacing_nm must be positive, got $(cfg.spacing_nm)"))
    cfg.sigma_par_nm > 0 || throw(ArgumentError("sigma_par_nm must be positive, got $(cfg.sigma_par_nm)"))
    cfg.sigma_perp_nm > 0 || throw(ArgumentError("sigma_perp_nm must be positive, got $(cfg.sigma_perp_nm)"))
    cfg.noise_sigma >= 0 || throw(ArgumentError("noise_sigma must be non-negative, got $(cfg.noise_sigma)"))
    cfg.blur_sigma_px >= 0 || throw(ArgumentError("blur_sigma_px must be non-negative, got $(cfg.blur_sigma_px)"))
    cfg.contrast_strength >= 0 || throw(ArgumentError("contrast_strength must be non-negative, got $(cfg.contrast_strength)"))
    cfg.affine_scale_jitter >= 0 || throw(ArgumentError("affine_scale_jitter must be non-negative, got $(cfg.affine_scale_jitter)"))
    cfg.drift_strength >= 0 || throw(ArgumentError("drift_strength must be non-negative, got $(cfg.drift_strength)"))
    cfg.row_offset_sigma >= 0 || throw(ArgumentError("row_offset_sigma must be non-negative, got $(cfg.row_offset_sigma)"))
    0 <= cfg.correlated_noise_frac <= 1 || throw(ArgumentError("correlated_noise_frac must be in [0,1], got $(cfg.correlated_noise_frac)"))
    isempty(ensemble.sites0) && throw(ArgumentError("Empty proxy ensemble: no GlcN (type=0) sites"))
    isempty(ensemble.sites1) && throw(ArgumentError("Empty proxy ensemble: no GlcNAc (type=1) sites"))
    return true
end

_global_state(rng::AbstractRNG) = (
    direction = rand(rng, 0:1),
    phase = rand(rng, 0:1),
    mirror = rand(rng, 0:1),
)

_parity_for_lobe(lobe::Int, n::Int, direction::Int, phase::Int) =
    mod((direction == 0 ? lobe : (n - lobe + 1)) - 1 + phase, 2)

function _transform_site(site::ProxySite, parity::Int, mirror::Int)
    t, u = site.t, site.u
    parity == 1 && (t = -t)
    mirror == 1 && (u = -u)
    return ProxySite(t, u, site.weight, site.sigma_t, site.sigma_u)
end

function _transformed_sites(sites::Vector{ProxySite}, parity::Int, mirror::Int)
    return [_transform_site(site, parity, mirror) for site in sites]
end

function lobe_proxy_sites(case::SyntheticCase, ensemble::ProxyEnsemble, lobe::Int)
    1 <= lobe <= case.truth.N || throw(ArgumentError("lobe index out of range"))
    parity = _parity_for_lobe(lobe, case.truth.N, case.truth.direction, case.truth.phase)
    sites = _sites_for_type(ensemble, case.truth.sequence[lobe], case.control)
    return _transformed_sites(sites, parity, case.truth.mirror)
end

function _sites_for_type(ensemble::ProxyEnsemble, typ::Int, control::ControlType)
    if control == CONTROL_IDENTICAL_MOLDS
        return ensemble.sites0
    elseif control == CONTROL_SWAPPED_TYPES
        return typ == 0 ? ensemble.sites1 : ensemble.sites0
    else
        return typ == 0 ? ensemble.sites0 : ensemble.sites1
    end
end

function _render_view(cfg::SimulatorConfig, xs_nm, ys_nm,
                      N::Int, sequence::Vector{Int},
                      lobe_x::Vector{Float64}, lobe_y::Vector{Float64},
                      amplitudes::Vector{Float64}, sigma_par::Vector{Float64},
                      sigma_perp::Vector{Float64}, theta::Float64,
                      direction::Int, phase::Int, mirror::Int,
                      ensemble::ProxyEnsemble, control::ControlType,
                      contrast_strength::Float64,
                      scale_x::Float64, scale_y::Float64, shear::Float64)
    z = zeros(Float64, cfg.height, cfg.width)
    control == CONTROL_NO_MOLECULE && return z
    cos_th = cos(theta)
    sin_th = sin(theta)
    for k in 1:N
        x0 = lobe_x[k]
        y0 = lobe_y[k]
        A = amplitudes[k]
        sp = sigma_par[k]
        su = sigma_perp[k]
        parity = _parity_for_lobe(k, N, direction, phase)
        sites = _transformed_sites(_sites_for_type(ensemble, sequence[k], control), parity, mirror)
        @inbounds for iy in 1:cfg.height
            y = ys_nm[iy]
            for ix in 1:cfg.width
                x = xs_nm[ix]
                xs = x / scale_x
                ys = (y - shear * x) / scale_y
                dx = xs - x0
                dy = ys - y0
                val = A * exp(-0.5 * (dx * dx / (sp * sp) + dy * dy / (su * su)))
                if contrast_strength > 0.0
                    for s in sites
                        sx = x0 + s.t * cos_th - s.u * sin_th
                        sy = y0 + s.t * sin_th + s.u * cos_th
                        dxs = xs - sx
                        dys = ys - sy
                        val += contrast_strength * s.weight *
                            exp(-0.5 * (dxs * dxs / (s.sigma_t * s.sigma_t) +
                                        dys * dys / (s.sigma_u * s.sigma_u)))
                    end
                end
                z[iy, ix] += val
            end
        end
    end
    return z
end
