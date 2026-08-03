function _generate_case_with_seed(case_seed::UInt64, cfg::SimulatorConfig,
                                  ensemble::ProxyEnsemble;
                                  control::ControlType=CONTROL_NORMAL,
                                  case_id::String="synthetic")
    validate_config(cfg, ensemble)
    control in ALL_CONTROLS || throw(ArgumentError("Unknown control: $control"))
    rng = MersenneTwister(case_seed)

    N = rand(rng, cfg.n_min:cfg.n_max)
    sequence = [rand(rng, 0:1) for _ in 1:N]
    cx = cfg.range_nm[1] / 2.0
    cy = cfg.range_nm[2] / 2.0
    gs = _global_state(rng)
    orientation_deg = (rand(rng) - 0.5) * 20.0
    theta = deg2rad(orientation_deg)
    spacing = cfg.spacing_nm * (1.0 + (rand(rng) - 0.5) * 0.15)
    t_start = -(N - 1) * spacing / 2.0
    lobe_t = [t_start + (k - 1) * spacing + (rand(rng) - 0.5) * 0.06 for k in 1:N]
    cos_th = cos(theta)
    sin_th = sin(theta)
    lobe_x = [cx + t * cos_th for t in lobe_t]
    lobe_y = [cy + t * sin_th for t in lobe_t]
    amplitudes = [cfg.amplitude_nm * (1.0 + (rand(rng) - 0.5) * 0.30) for _ in 1:N]
    sigma_par = [cfg.sigma_par_nm * (1.0 + (rand(rng) - 0.5) * 0.20) for _ in 1:N]
    sigma_perp = [cfg.sigma_perp_nm * (1.0 + (rand(rng) - 0.5) * 0.20) for _ in 1:N]
    xs_nm = collect(range(0.0, cfg.range_nm[1]; length=cfg.width))
    ys_nm = collect(range(0.0, cfg.range_nm[2]; length=cfg.height))

    scale_x_fwd = 1.0 + (rand(rng) - 0.5) * 2.0 * cfg.affine_scale_jitter
    scale_y_fwd = 1.0 + (rand(rng) - 0.5) * 2.0 * cfg.affine_scale_jitter
    shear_fwd = (rand(rng) - 0.5) * 2.0 * cfg.drift_strength
    blur_fwd = cfg.blur_sigma_px * (1.0 + (rand(rng) - 0.5) * 0.50)
    scale_x_bwd = 1.0 + (rand(rng) - 0.5) * 2.0 * cfg.affine_scale_jitter
    scale_y_bwd = 1.0 + (rand(rng) - 0.5) * 2.0 * cfg.affine_scale_jitter
    shear_bwd = (rand(rng) - 0.5) * 2.0 * cfg.drift_strength
    blur_bwd = cfg.blur_sigma_px * (1.0 + (rand(rng) - 0.5) * 0.50)
    tilt_x = (rand(rng) - 0.5) * 0.04
    tilt_y = (rand(rng) - 0.5) * 0.02

    fwd_nuisance = ViewNuisanceTruth(scale_x_fwd, scale_y_fwd, shear_fwd, blur_fwd)
    bwd_nuisance = ViewNuisanceTruth(scale_x_bwd, scale_y_bwd, shear_bwd, blur_bwd)

    z_fwd = _render_view(cfg, xs_nm, ys_nm, N, sequence, lobe_x, lobe_y,
                         amplitudes, sigma_par, sigma_perp, theta,
                         gs.direction, gs.phase, gs.mirror,
                         ensemble, control, cfg.contrast_strength,
                         fwd_nuisance.scale_x, fwd_nuisance.scale_y, fwd_nuisance.shear)
    z_bwd = _render_view(cfg, xs_nm, ys_nm, N, sequence, lobe_x, lobe_y,
                         amplitudes, sigma_par, sigma_perp, theta,
                         gs.direction, gs.phase, gs.mirror,
                         ensemble, control, cfg.contrast_strength,
                         bwd_nuisance.scale_x, bwd_nuisance.scale_y, bwd_nuisance.shear)

    z_fwd = _gaussian_blur(z_fwd, blur_fwd)
    z_bwd = _gaussian_blur(z_bwd, blur_bwd)

    # no_molecule still carries background nuisance when configured: tilt/noise/row offsets.
    if control != CONTROL_NO_MOLECULE || cfg.noise_sigma > 0.0 || cfg.row_offset_sigma > 0.0
        @inbounds for iy in 1:cfg.height, ix in 1:cfg.width
            v = tilt_x * (xs_nm[ix] - cx) + tilt_y * (ys_nm[iy] - cy)
            z_fwd[iy, ix] += v
            z_bwd[iy, ix] += v
        end
    end

    _add_row_offsets!(z_fwd, rng, cfg.row_offset_sigma)
    _add_row_offsets!(z_bwd, rng, cfg.row_offset_sigma)
    _add_noise!(z_fwd, z_bwd, rng, cfg.noise_sigma, cfg.correlated_noise_frac)

    if control == CONTROL_CORRUPTED_VIEW
        y0 = max(1, cfg.height ÷ 3)
        y1 = min(cfg.height, y0 + 5)
        x0 = max(1, cfg.width ÷ 3)
        x1 = min(cfg.width, x0 + 10)
        z_bwd[y0:y1, x0:x1] .= NaN
    end

    channels = SXMChannel[SXMChannel("Z", "arb", "fwd", z_fwd)]
    control != CONTROL_MISSING_BWD && push!(channels, SXMChannel("Z", "arb", "bwd", z_bwd))
    data_info = control == CONTROL_MISSING_BWD ?
        "Channel Z arb fwd" : "Channel Z arb fwd\r\nChannel Z arb bwd"
    header = Dict{String,String}(
        "SCAN_PIXELS" => "$(cfg.width) $(cfg.height)",
        "SCAN_RANGE" => "$(cfg.range_nm[1] * 1e-9) $(cfg.range_nm[2] * 1e-9)",
        "SCAN_OFFSET" => "0 0",
        "DATA_INFO" => data_info,
    )
    img = SXMImage(case_id, header, cfg.width, cfg.height, cfg.range_nm, (0.0, 0.0), channels)
    truth = SyntheticTruth(N, sequence, copy(lobe_x), copy(lobe_y), amplitudes,
                           sigma_par, sigma_perp, orientation_deg,
                           gs.direction, gs.phase, gs.mirror,
                           fwd_nuisance, bwd_nuisance, control)
    return SyntheticCase(case_id, img, truth, control, case_seed)
end

function replay_case(case_seed::UInt64, cfg::SimulatorConfig, ensemble::ProxyEnsemble;
                     control::ControlType=CONTROL_NORMAL,
                     case_id::String="synthetic")
    return _generate_case_with_seed(case_seed, cfg, ensemble; control=control, case_id=case_id)
end

function generate_case(rng::AbstractRNG, cfg::SimulatorConfig, ensemble::ProxyEnsemble;
                       control::ControlType=CONTROL_NORMAL,
                       case_id::String="synthetic")
    case_seed = rand(rng, UInt64)
    return replay_case(case_seed, cfg, ensemble; control=control, case_id=case_id)
end

generate_case(case_seed::UInt64, cfg::SimulatorConfig, ensemble::ProxyEnsemble;
              control::ControlType=CONTROL_NORMAL,
              case_id::String="synthetic") =
    replay_case(case_seed, cfg, ensemble; control=control, case_id=case_id)

function generate_batch(rng::AbstractRNG, cfg::SimulatorConfig, ensemble::ProxyEnsemble;
                        n_cases::Int=32,
                        controls=fill(CONTROL_NORMAL, n_cases))
    length(controls) == n_cases ||
        throw(ArgumentError("controls length ($(length(controls))) must equal n_cases ($n_cases)"))
    cases = SyntheticCase[]
    for i in 1:n_cases
        push!(cases, generate_case(rng, cfg, ensemble; control=controls[i],
                                   case_id=@sprintf("synthetic_%03d", i)))
    end
    return cases
end
