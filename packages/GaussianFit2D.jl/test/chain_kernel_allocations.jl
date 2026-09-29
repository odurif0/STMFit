# Literal pre-barrier arithmetic. Keep separate from the production kernel so
# this regression checks both numerical identity and the allocation failure.
function chain_values_before_barrier(x, y, p, n, ctx, cfg; amp_min=NaN, amp_range=NaN)
    G = GaussianFit2D
    b0, feats, ts, us, _, _ = G._decode_chain(p, n, ctx, cfg; amp_min, amp_range)
    if cfg.chain_tilted_baseline
        bx = p[2]; by = p[3]
        pred = @. b0 + bx*x + by*y
    else
        pred = fill(b0, length(x))
    end
    peak_axes = G._chain_peak_axes(ts, us, ctx, cfg)
    split_profile = G._chain_uses_split_profile(cfg)
    skew_rmax = max(cfg.skew_ratio_max, 1.0 + G.EPS)
    for (k, f) in enumerate(feats)
        ax, ay = peak_axes[k]
        if split_profile
            for i in eachindex(pred)
                dt = (x[i] - f.x_nm) * ax + (y[i] - f.y_nm) * ay
                du = (x[i] - f.x_nm) * (-ay) + (y[i] - f.y_nm) * ax
                pred[i] += f.amplitude * G._chain_split_peak_value(dt, du, f.sigma_x_nm, f.sigma_y_nm, f.skew_ratio, skew_rmax)
            end
        else
            @. pred += f.amplitude * exp(-0.5 * ((((x - f.x_nm)*ax + (y - f.y_nm)*ay)/f.sigma_x_nm)^2 + (((x - f.x_nm)*(-ay) + (y - f.y_nm)*ax)/f.sigma_y_nm)^2))
        end
    end
    return pred
end

@testset "Chain kernel: exact arithmetic and bounded allocations" begin
    G = GaussianFit2D
    ctx = (origin=(0.,0.), axis=[.8,.6], perp=[-.6,.8], tmin=-1.7, tmax=1.7)
    x = collect(range(-1.5,1.5; length=10000)); y = .2sin.(3x)
    for profile in (:gaussian, :split),
        circular in (true, false), tilted in (true, false), bounded in (true, false), n in (1, 6)
        cfg = G.ChainSweepConfig(peak_profile=profile, chain_circular_sigmas=circular,
            chain_tilted_baseline=tilted)
        p = zeros(G._chain_nparams(n, cfg))
        p[1] = .01
        tilted && (p[2:3] .= [.003, -.007])
        first = (tilted ? 3 : 1) + n + G._chain_spacing_param_count(n, cfg)
        p[first+1:first+n] .= [.2,-.3,.1,.4,.25,-.1][1:n]
        profile == :split && (p[end-n+1:end] .= [-.4,.2,.6,-.8,.3,.15][1:n])
        amp_min, amp_range = bounded ? (.03, .08) : (NaN, NaN)
        expected = chain_values_before_barrier(x, y, p, n, ctx, cfg; amp_min, amp_range)
        actual = G._chain_model_values(x, y, p, n, ctx, cfg; amp_min, amp_range)
        @test reinterpret(UInt64, actual) == reinterpret(UInt64, expected)
        if n == 6 && bounded
            # Warmed above. A generous ceiling avoids timing-sensitive tests:
            # the repaired path uses ~0.1 MB, the regressed split path ~19 MB.
            allocated = @allocated G._chain_model_values(x, y, p, n, ctx, cfg; amp_min, amp_range)
            @test allocated < 2_000_000
        end
    end
end
