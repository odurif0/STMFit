#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# synthetic_known_n_validation.jl — Synthetic known-N validation
# for selector comparison.
#
# Generates synthetic Gaussian-chain STM images in memory with known true_N,
# runs either a fast 2D circular sweep or a full circ→ell local-refinement
# sweep, applies core selectors from STMMolecularFit.jl/src/selectors.jl, and
# writes a TSV of selected N and absolute errors.  The known true_N is used only
# for output grading, never to set the fitted candidate range or selector.
#
# Usage:
#   julia --project=. test/synthetic_known_n_validation.jl
#   julia --project=. test/synthetic_known_n_validation.jl --cases 12 --seed 42
#   julia --project=. test/synthetic_known_n_validation.jl --policies gcv,stability_selection
#   julia --project=. test/synthetic_known_n_validation.jl --out results/my_summary.tsv
#   julia --project=. test/synthetic_known_n_validation.jl --noise-scale 0.5 --cases 8
# ──────────────────────────────────────────────────────────────────────────────

using STMMolecularFit, GaussianFit2D, GaussianFit1D
using LinearAlgebra
using DelimitedFiles, Printf, Statistics, Random

# ══════════════════════════════════════════════════════════════════════════════
# Imports from STMMolecularFit (un-exported selectors)
# ══════════════════════════════════════════════════════════════════════════════
import STMMolecularFit: _integrated_robust_aicc_n,
    _support_marginalized_gcv_selection, _stability_selection,
    _laplace_evidence_selection,
    _fwd_bwd_consensus_selection,
    SUPPORT_MARG_REGRET_MARGIN

# ══════════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════════

const SUMMARY_HEADER = ["case_id", "seed", "true_N", "artifact",
                        "policy", "N_eff", "N_selected", "abs_error",
                        "status", "score_or_source", "noise_scale",
                        "mode"]

# Default policies (comma-separated)
const DEFAULT_POLICIES = [
    "gcv",
    "gcv_with_robust_aicc_guard",
    "laplace_evidence_guard",
    "support_marginalized_gcv_guard",
    "stability_selection",
    "fwd_bwd_consensus",
]

# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

function _parse_cli()
    n_cases = 6
    seed = 1234
    outpath = "results/synthetic_known_n/summary.tsv"
    policies_str = join(DEFAULT_POLICIES, ",")
    noise_scale = 1.0
    mode = "circular"

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--cases" && i < length(ARGS)
            n_cases = parse(Int, ARGS[i+1]); i += 2
        elseif arg == "--seed" && i < length(ARGS)
            seed = parse(Int, ARGS[i+1]); i += 2
        elseif arg == "--out" && i < length(ARGS)
            outpath = ARGS[i+1]; i += 2
        elseif arg == "--policies" && i < length(ARGS)
            policies_str = ARGS[i+1]; i += 2
        elseif arg == "--noise-scale" && i < length(ARGS)
            noise_scale = parse(Float64, ARGS[i+1]); i += 2
        elseif arg == "--mode" && i < length(ARGS)
            mode = ARGS[i+1]; i += 2
        else
            i += 1
        end
    end

    policies = split(policies_str, ",")
    return n_cases, seed, outpath, policies, noise_scale, mode
end

# ══════════════════════════════════════════════════════════════════════════════
# Synthetic image generation
# ══════════════════════════════════════════════════════════════════════════════

"""
    _gaussian2d(x, y, x0, y0, A, sx, sy)

Evaluate a 2D Gaussian centred at (x0, y0) with amplitude A and widths
sx (parallel) and sy (perpendicular).  No rotation — the chain is
pre-aligned along the image x-axis.
"""
function _gaussian2d(x, y, x0, y0, A, sx, sy)
    return A * exp(-0.5 * ((x - x0)^2 / sx^2 + (y - y0)^2 / sy^2))
end

"""
    _generate_synthetic_case(case_idx, true_N, rng, noise_scale=1.0;
                             width=128, height=64)

Generate an SXMImage with a Gaussian chain along x.
`noise_scale` multiplies the independent fwd/bwd noise amplitude.
Returns (img, artifact_description).
"""
function _generate_synthetic_case(case_idx, true_N, rng, noise_scale=1.0;
                                  width=128, height=64)
    # Physical extent: ~8 nm × 4 nm
    range_nm = (8.0, 4.0)

    xs_nm = collect(range(0.0, range_nm[1]; length=width))
    ys_nm = collect(range(0.0, range_nm[2]; length=height))

    # Chain centre
    cx, cy = range_nm[1] / 2, range_nm[2] / 2

    # Baseline
    baseline = 0.0
    tilt_x = 0.02   # mild tilt (nm / nm)
    tilt_y = 0.01

    # Chain parameters (nm)
    spacing = 0.62
    sigma_par = 0.20
    sigma_perp = 0.14

    # Generate lobe positions with jitter
    n_eff = true_N
    t_start = -(n_eff - 1) * spacing / 2.0
    lobe_positions = Float64[]
    for k in 1:n_eff
        t = t_start + (k - 1) * spacing + (rand(rng) - 0.5) * 0.08
        push!(lobe_positions, t)
    end

    # Amplitudes with jitter; decrease slightly toward edges
    amplitudes = Float64[]
    for k in 1:n_eff
        base = 0.35 + 0.15 * (1.0 - abs(k - (n_eff + 1) / 2) / n_eff)
        push!(amplitudes, base * (1.0 + (rand(rng) - 0.5) * 0.30))
    end

    # Sigmas with jitter
    spar_vals = Float64[]
    sperp_vals = Float64[]
    for k in 1:n_eff
        push!(spar_vals, sigma_par * (1.0 + (rand(rng) - 0.5) * 0.20))
        push!(sperp_vals, sigma_perp * (1.0 + (rand(rng) - 0.5) * 0.20))
    end

    # Build fwd and bwd z-data matrices
    z_fwd = zeros(height, width)
    z_bwd = zeros(height, width)

    for iy in 1:height, ix in 1:width
        x_nm = xs_nm[ix]
        y_nm = ys_nm[iy]

        # Baseline + tilt
        val_fwd = baseline + tilt_x * (x_nm - cx) + tilt_y * (y_nm - cy)
        val_bwd = val_fwd  # same baseline for both

        # Add each lobe
        for k in 1:n_eff
            x0 = cx + lobe_positions[k]
            y0 = cy
            A = amplitudes[k]
            spar = spar_vals[k]
            sperp = sperp_vals[k]
            val_fwd += _gaussian2d(x_nm, y_nm, x0, y0, A, spar, sperp)
            val_bwd += _gaussian2d(x_nm, y_nm, x0, y0, A, spar, sperp)
        end

        z_fwd[iy, ix] = val_fwd
        z_bwd[iy, ix] = val_bwd
    end

    # Independent Gaussian noise (scaled by noise_scale)
    noise_level = 0.03 * noise_scale
    z_fwd .+= randn(rng, height, width) .* noise_level
    z_bwd .+= randn(rng, height, width) .* noise_level

    # ── Artifacts based on case index ──
    artifact = "none"
    if mod(case_idx, 4) == 1
        # Fwd-only extra blob (false positive lobe)
        artifact = "fwd_extra_blob"
        x_extra = cx + (n_eff ÷ 2) * spacing + rand(rng) * 0.2
        for iy in 1:height, ix in 1:width
            z_fwd[iy, ix] += _gaussian2d(xs_nm[ix], ys_nm[iy],
                                          x_extra, cy,
                                          0.25, sigma_par, sigma_perp)
        end
    elseif mod(case_idx, 4) == 2
        # Stripe artifact (horizontal line in fwd only)
        artifact = "fwd_stripe"
        stripe_y = cy + (rand(rng) - 0.5) * 0.6
        for iy in 1:height, ix in 1:width
            if abs(ys_nm[iy] - stripe_y) < 0.03
                z_fwd[iy, ix] += 0.15 * (1.0 + (rand(rng) - 0.5) * 1.0)
            end
        end
    elseif mod(case_idx, 4) == 3
        # Fwd missing lobe (near-zero amplitude) + mild stripe
        artifact = "fwd_missing_lobe_mild_stripe"
        miss_k = rand(rng, 1:n_eff)
        for iy in 1:height, ix in 1:width
            x_nm = xs_nm[ix]
            y_nm = ys_nm[iy]
            x0 = cx + lobe_positions[miss_k]
            z_fwd[iy, ix] -= _gaussian2d(x_nm, y_nm, x0, cy,
                                          amplitudes[miss_k], spar_vals[miss_k],
                                          sperp_vals[miss_k])
        end
        stripe_y2 = cy + (rand(rng) - 0.5) * 0.4
        for iy in 1:height, ix in 1:width
            if abs(ys_nm[iy] - stripe_y2) < 0.025
                z_fwd[iy, ix] += 0.10
            end
        end
    else
        # No artifact for mod 0
        artifact = "none"
    end

    # Build SXMChannel objects
    # Qualify with GaussianFit2D to avoid ambiguity (same struct in STMMolecularFit)
    ch_fwd = GaussianFit2D.SXMChannel("Z", "arb", "fwd", z_fwd)
    ch_bwd = GaussianFit2D.SXMChannel("Z", "arb", "bwd", z_bwd)

    # Build SXMImage
    case_id = @sprintf("synthetic_case_%03d", case_idx)
    header = Dict{String,String}(
        "SCAN_PIXELS" => "$width $height",
        "SCAN_RANGE" => "$(range_nm[1] * 1e-9) $(range_nm[2] * 1e-9)",
        "SCAN_OFFSET" => "0 0",
        "DATA_INFO" => "Channel Z arb fwd\r\nChannel Z arb bwd",
    )
    img = GaussianFit2D.SXMImage(case_id, header, width, height, range_nm, (0.0, 0.0),
                                  [ch_fwd, ch_bwd])
    return img, artifact
end

# ══════════════════════════════════════════════════════════════════════════════
# Build configs
# ══════════════════════════════════════════════════════════════════════════════

function _build_configs(case_id)
    pcfg = GaussianFit2D.PatternConfig(
        filepath=case_id, channel="Z", direction="fwd",
        roi_channel="Z", stride=1, flatten="none",
        smooth_radius_px=1,
        output_dir="results/synthetic_known_n",
        no_plot=true,
    )

    ccfg = GaussianFit2D.ChainSweepConfig(
        # Fixed synthetic search window.  Do not derive this from true_N: labels
        # are for external grading only, not for fitting or selection.
        n_min=2,
        n_max=10,
        intelligent_sweep=false,
        chain_circular_sigmas=true,
        chain_tilted_baseline=true,
        skip_global=false,
        multistart=1,
        global_maxtime=2.0,
        global_maxiter=3000,
        max_iter=120,
        spacing_min_nm=0.35,
        spacing_max_nm=0.80,
        fit_width_nm=0.35,
        support_noise_k=2.0,
        support_padding_nm=0.4,
        sigma_parallel_min_nm=0.10,
        sigma_parallel_max_nm=0.35,
        sigma_perp_min_nm=0.10,
        sigma_perp_max_nm=0.35,
        max_overlap=0.75,
        kappa_max=10.0,
        selection_criterion="gcv",
        fuse_z_bwd=true,
    )

    return pcfg, ccfg
end

# ══════════════════════════════════════════════════════════════════════════════
# Baseline GCV selection (lowest valid finite-gcv candidate)
# ══════════════════════════════════════════════════════════════════════════════

function _baseline_gcv(results)
    best_n = 0
    best_gcv = Inf
    for r in results
        r.success && r.valid || continue
        isfinite(r.gcv) || continue
        if r.gcv < best_gcv
            best_n = r.n
            best_gcv = r.gcv
        end
    end
    return best_n
end

# ══════════════════════════════════════════════════════════════════════════════
# Per-policy application
# ══════════════════════════════════════════════════════════════════════════════

"""
Apply a single policy and return (N_selected, status, score_or_source).

N_selected is an Int or nothing on failure.
status is "ok" or "error:<message>".
score_or_source is a string description or numeric value.

`ccfg` is the config used for data-tolerant selectors (robust_aicc, etc.).
`results_circ` are always available.  `results_ell` may be the same as
`results_circ` (circular-only mode) or independent elliptical refinements
(circ_ell mode).  Selectors that distinguish ell vs circ receive both sets.
"""
function _apply_policy(policy, img, pcfg, ccfg, results_circ, n_eff, eff_source,
                        results_ell=nothing)
    if results_ell === nothing
        results_ell = results_circ
    end
    # ccfg_ell and ccfg_circ: use the caller's ccfg; in circ_ell mode the
    # caller may have set chain_circular_sigmas=false on ccfg already.
    ccfg_ell = ccfg
    ccfg_circ = deepcopy(ccfg)
    ccfg_circ.chain_circular_sigmas = true

    # Baseline GCV is always available
    if policy == "gcv"
        n_eff > 0 || return (nothing, "error:no_valid_gcv", "NA")
        return (n_eff, "ok", eff_source)
    end

    # gcv_with_robust_aicc_guard: downshift only if robust_n < n_eff
    if policy == "gcv_with_robust_aicc_guard"
        try
            robust_n_raw = _integrated_robust_aicc_n(img, pcfg, ccfg; nu=8.0)
            robust_n, robust_source, robust_score = robust_n_raw
            if robust_n !== nothing && robust_n < n_eff
                return (robust_n, "ok", robust_source)
            else
                return (n_eff, "ok", eff_source)
            end
        catch err
            return (nothing, "error:$err", "NA")
        end
    end

    # laplace_evidence_guard: cap downshift to one lobe
    if policy == "laplace_evidence_guard"
        try
            lap_raw = _laplace_evidence_selection(img, pcfg, ccfg_ell, ccfg_circ,
                                                   results_ell, results_circ)
            lap_n = lap_raw[1]
            if lap_n !== nothing && isfinite(lap_n) && lap_n < n_eff
                # Cap downshift to at most one lobe (n_eff - 1)
                capped = max(lap_n, n_eff - 1)
                return (capped, "ok", @sprintf("laplace_capped_%d", lap_n))
            else
                return (n_eff, "ok", eff_source)
            end
        catch err
            return (nothing, "error:$err", "NA")
        end
    end

    # support_marginalized_gcv_guard: guard similar to batch_full
    if policy == "support_marginalized_gcv_guard"
        try
            sm_raw = _support_marginalized_gcv_selection(img, pcfg, ccfg_ell,
                                                          results_ell, results_circ)
            sm_n = sm_raw[1]
            sm_regret = sm_raw[2]
            if sm_n !== nothing && isfinite(sm_n) && isfinite(sm_regret)
                if sm_regret < SUPPORT_MARG_REGRET_MARGIN
                    # unambiguous: use SM result
                    if sm_n < n_eff
                        return (sm_n, "ok", "support_marginalized_gcv")
                    else
                        return (n_eff, "ok", eff_source)
                    end
                else
                    # ambiguous: guard up to the ambiguity threshold (rel 0.05)
                    # Keep n_eff if the SM regret is too high
                    return (n_eff, "ok", @sprintf("support_marg_ambiguous_%.4f", sm_regret))
                end
            else
                return (n_eff, "ok", eff_source)
            end
        catch err
            return (nothing, "error:$err", "NA")
        end
    end

    # stability_selection
    if policy == "stability_selection"
        try
            stab_raw = _stability_selection(img, pcfg, ccfg_ell,
                                             results_ell, results_circ)
            if stab_raw !== nothing
                stab_n = stab_raw[1]
                if stab_n !== nothing && isfinite(stab_n)
                    return (stab_n, "ok", @sprintf("stability_%.0f%%", stab_raw[2]))
                end
            end
            return (n_eff, "ok", eff_source)
        catch err
            return (nothing, "error:$err", "NA")
        end
    end

    # fwd_bwd_consensus
    if policy == "fwd_bwd_consensus"
        try
            consensus_raw = _fwd_bwd_consensus_selection(img, pcfg, ccfg_ell,
                                                          results_ell, results_circ)
            if consensus_raw !== nothing
                cons_n = consensus_raw[1]
                if cons_n !== nothing && isfinite(cons_n)
                    return (cons_n, "ok", @sprintf("fwd_bwd_consensus_%.2f", consensus_raw[2]))
                end
            end
            return (n_eff, "ok", eff_source)
        catch err
            return (nothing, "error:$err", "NA")
        end
    end

    return (nothing, "error:unknown_policy_$(policy)", "NA")
end

# ══════════════════════════════════════════════════════════════════════════════
# Local elliptical refinement (phase-2)
# ══════════════════════════════════════════════════════════════════════════════

"""
    _circ_ell_refinement(img, pcfg, ccfg, results_circ, ctx)

For each valid N from the circular sweep, refit with elliptical sigmas
(chain_circular_sigmas=false) using the same data and axis context.
Returns the elliptical refinement results.
"""
function _circ_ell_refinement(img, pcfg, ccfg, results_circ, ctx)
    ccfg_ell = deepcopy(ccfg)
    ccfg_ell.chain_circular_sigmas = false
    ccfg_ell.intelligent_sweep = false
    ccfg_ell.multistart = 1  # minimal starts for refinement

    # Reuse the circular sweep's loaded data and axis to avoid re-computation
    override_data = (ctx.xs, ctx.ys, ctx.zimg, ctx.mask,
                     ctx.x, ctx.y, ctx.z, ctx.noise)
    override_axisctx = ctx.axisctx_full

    valid_Ns = sort(unique(r.n for r in results_circ if r.success && r.valid))
    results_ell = GaussianFit2D.ChainModelResult[]

    for n in valid_Ns
        ccfg_n = deepcopy(ccfg_ell)
        ccfg_n.n_min = n
        ccfg_n.n_max = n

        try
            ell_results, _, _ = GaussianFit2D.chain_gaussian_sweep(
                img, pcfg, ccfg_n;
                override_data=override_data,
                override_axisctx=override_axisctx)
            for r in ell_results
                r.success && r.valid && push!(results_ell, r)
            end
        catch err
            @warn "  Elliptical refinement failed for N=$n: $err"
        end
    end

    @printf("  ell refinement: %d / %d valid\n", count(r -> r.valid, results_ell), length(results_ell))
    return results_ell
end

"""
    _best_gcv_n(results)

Find the N with the lowest valid finite GCV in a result set.
Returns (n, gcv) or (0, Inf).
"""
function _best_gcv_n(results)
    best_n = 0
    best_gcv = Inf
    for r in results
        r.success && r.valid && isfinite(r.gcv) || continue
        if r.gcv < best_gcv
            best_n = r.n
            best_gcv = r.gcv
        end
    end
    return best_n, best_gcv
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    n_cases, seed, outpath, policies, noise_scale, mode = _parse_cli()

    @printf("Synthetic known-N validation\n")
    @printf("  mode:        %s\n", mode)
    @printf("  cases:       %d\n", n_cases)
    @printf("  seed:        %d\n", seed)
    @printf("  noise_scale: %g\n", noise_scale)
    @printf("  out:         %s\n", outpath)
    @printf("  policies:    %s\n\n", join(policies, ", "))

    mkpath(dirname(outpath))

    rng = MersenneTwister(seed)
    rows = Vector{Vector{String}}()

    # Keep aggregate counts per policy
    agg_exact = Dict{String,Int}()
    agg_total = Dict{String,Int}()
    agg_abs_err = Dict{String,Float64}()
    for p in policies
        agg_exact[p] = 0
        agg_total[p] = 0
        agg_abs_err[p] = 0.0
    end

    # Cycle true_N 4..8
    true_Ns = [mod1(i, 5) + 3 for i in 1:n_cases]  # cycles 4,5,6,7,8,4,5,...

    for case_idx in 1:n_cases
        true_N = true_Ns[case_idx]
        case_id = @sprintf("synthetic_case_%03d", case_idx)

        @printf("--- Case %d/%d: %s true_N=%d ---\n",
                case_idx, n_cases, case_id, true_N)

        # Generate image (with noise_scale applied)
        img, artifact = _generate_synthetic_case(case_idx, true_N, rng, noise_scale)

        pcfg, ccfg = _build_configs(case_id)

        # ── Circular sweep (always run) ──
        results_circ, best_circ, ctx = try
            GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
        catch err
            @warn "Sweep failed for case $case_id: $err"
            for policy in policies
                row = [case_id, string(seed), string(true_N), artifact,
                       policy, "NA", "NA", "NA",
                       "error:sweep_failed:$err", "NA",
                       @sprintf("%g", noise_scale), mode]
                push!(rows, row)
            end
            continue
        end

        # ── Elliptical refinement (phase-2, only in circ_ell mode) ──
        results_ell = nothing
        if mode == "circ_ell"
            @printf("  circ_ell refinement...\n")
            results_ell = _circ_ell_refinement(img, pcfg, ccfg, results_circ, ctx)
        end

        # ── Determine n_eff ──
        # In circular mode: use best circular GCV N.
        # In circ_ell mode: use the best GCV score across circular and
        # elliptical candidates.  This mirrors the effective min(circ, ell)
        # idea without using the synthetic true_N label.
        n_circ, gcv_circ = _best_gcv_n(results_circ)
        n_eff = n_circ
        eff_source = "gcv_circ"

        if mode == "circ_ell" && results_ell !== nothing && !isempty(results_ell)
            n_ell, gcv_ell = _best_gcv_n(results_ell)
            if n_ell > 0
                if gcv_ell < gcv_circ
                    n_eff = n_ell
                    eff_source = "gcv_ell_min"
                end
            end
        end

        # Fallback if no valid result from either set
        if n_eff <= 0
            if results_circ !== nothing
                n_eff = _baseline_gcv(results_circ)
                n_eff > 0 || (n_eff = 0)
            end
            eff_source = "gcv"
        end

        @printf("  n_eff=%d (%s)\n", n_eff, eff_source)

        # ── Apply each policy ──
        # In circular mode: pass same results for ell and circ (existing behaviour).
        # In circ_ell mode: pass elliptical refinement as results_ell.
        pol_results_ell = (mode == "circ_ell" && results_ell !== nothing) ? results_ell : results_circ
        pol_results_circ = results_circ
        pol_ccfg_ell = (mode == "circ_ell" && results_ell !== nothing) ? deepcopy(ccfg) : ccfg
        if mode == "circ_ell" && results_ell !== nothing
            pol_ccfg_ell.chain_circular_sigmas = false
        end

        for policy in policies
            n_sel, status, score_or_source = _apply_policy(
                policy, img, pcfg, pol_ccfg_ell, pol_results_circ,
                n_eff, eff_source, pol_results_ell)

            n_sel_str = n_sel === nothing ? "NA" : string(n_sel)
            abs_err = (n_sel === nothing || status != "ok") ? "NA" :
                       string(abs(n_sel - true_N))

            row = [case_id, string(seed), string(true_N), artifact,
                   policy, string(n_eff), n_sel_str, abs_err,
                   status, score_or_source, @sprintf("%g", noise_scale),
                   mode]
            push!(rows, row)

            # Aggregate
            p_key = policy
            if haskey(agg_total, p_key)
                agg_total[p_key] += 1
                if status == "ok" && n_sel !== nothing
                    agg_abs_err[p_key] += abs(n_sel - true_N)
                    if n_sel == true_N
                        agg_exact[p_key] += 1
                    end
                end
            end
        end
    end

    # ── Write TSV ──
    open(outpath, "w") do io
        println(io, join(SUMMARY_HEADER, "\t"))
        for row in rows
            println(io, join(row, "\t"))
        end
    end
    @printf("\nWrote %d rows to %s\n", length(rows), outpath)

    # ── Aggregate summary ──
    @printf("\n%s\n", "="^75)
    @printf("Aggregate summary (seed=%d)\n", seed)
    @printf("%s\n", "="^75)
    @printf("%-40s %7s %7s %12s\n", "Policy", "Exact", "Total", "Mean|Δ|")
    @printf("%s\n", "-"^75)
    for policy in policies
        t = agg_total[policy]
        e = agg_exact[policy]
        mae = t > 0 ? agg_abs_err[policy] / t : NaN
        @printf("%-40s %7d %7d %12.4f\n", policy, e, t, mae)
    end
    @printf("%s\n", "="^75)

    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
