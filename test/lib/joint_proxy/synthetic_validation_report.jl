_finite_valid_count_candidate(c) = hasproperty(c, :valid) && hasproperty(c, :joint_gcv) && hasproperty(c, :n) &&
    getproperty(c, :valid) && isfinite(getproperty(c, :joint_gcv)) && isfinite(getproperty(c, :n))

_peaks(ts, prof) = begin
    isempty(prof) && return Float64[]
    thr = maximum(prof) * 0.35
    idx = [i for i in 2:length(prof)-1 if prof[i] >= thr && prof[i] >= prof[i - 1] && prof[i] > prof[i + 1]]
    isempty(idx) ? Float64[] : sort(ts[sort(idx; by=i -> prof[i], rev=true)])
end

_finite_mean(vals; default::Float64=0.0) = begin
    xs = [Float64(v) for v in vals if isfinite(v)]
    isempty(xs) ? default : mean(xs)
end

_finite_max(vals; default::Float64=0.0) = begin
    xs = [Float64(v) for v in vals if isfinite(v)]
    isempty(xs) ? default : maximum(xs)
end

function _backbone_residual(raw_patch::Vector{Float64}, pts, sigma_parallel_nm::Float64, sigma_perp_nm::Float64)
    n = length(raw_patch)
    n == length(pts) || return [isfinite(v) ? v - _finite_mean(raw_patch) : 0.0 for v in raw_patch]
    σt2 = sigma_parallel_nm^2
    σu2 = sigma_perp_nm^2
    σt2 > 0 && σu2 > 0 && isfinite(σt2) && isfinite(σu2) || return [isfinite(v) ? v - _finite_mean(raw_patch) : 0.0 for v in raw_patch]
    basis = Vector{Float64}(undef, n)
    finite_idx = Int[]
    for i in 1:n
        t = Float64(getproperty(pts[i], :t))
        u = Float64(getproperty(pts[i], :u))
        basis[i] = exp(-0.5 * ((t * t) / σt2 + (u * u) / σu2))
        if isfinite(raw_patch[i]) && isfinite(basis[i])
            push!(finite_idx, i)
        end
    end
    if length(finite_idx) < 2
        meanv = _finite_mean(raw_patch)
        return [isfinite(v) ? v - meanv : 0.0 for v in raw_patch]
    end
    X = hcat(ones(length(finite_idx)), basis[finite_idx])
    y = raw_patch[finite_idx]
    β = X \ y
    res = Vector{Float64}(undef, n)
    for i in 1:n
        res[i] = isfinite(raw_patch[i]) && isfinite(basis[i]) ? raw_patch[i] - (β[1] + β[2] * basis[i]) : 0.0
    end
    return res
end

function _patch(view::ViewData, ax, t0::Float64, u0::Float64, half_nm::Float64, step_nm::Float64)
    coords = collect(-half_nm:step_nm:half_nm)
    vals = Float64[]; pts = NamedTuple{(:t,:u),Tuple{Float64,Float64}}[]
    for u in coords, t in coords
        push!(pts, (; t=t, u=u))
        x = ax.cx + (t0 + t) * ax.ax + (u0 + u) * ax.ux
        y = ax.cy + (t0 + t) * ax.ay + (u0 + u) * ax.uy
        push!(vals, _sample(view, x, y))
    end
    return vals, pts
end

function physical_candidate_report(case::SyntheticCase, reg::ProxyEnsemble; swap_type_mapping::Bool=false, n_min::Int=2, n_max::Int=12)
    cfg = Main.JointProxySimulator.SimulatorConfig()
    fwd = _preprocess_view(_view_data(case, "fwd"))
    has_bwd = any(ch -> lowercase(ch.direction) == "bwd", case.img.channels)
    bwd = has_bwd ? _preprocess_view(_view_data(case, "bwd")) : nothing
    bwd_missing = bwd === nothing || !all(isfinite, bwd.z)
    view_count = bwd_missing ? 1 : 2
    ax = _axis(fwd)
    ts, prof, signal = _profile(fwd, ax)
    span = isempty(ts) ? (-1.0, 1.0) : (first(ts), last(ts))
    rho = bwd_missing ? NaN : cor(vec(fwd.z), vec(bwd.z))
    eff = effective_view_factor(rho)
    # A missing/corrupted backward view falls back to the valid forward view;
    # only insufficient molecular signal makes the count geometry ambiguous.
    ambiguous = signal < 0.5
    flat = ambiguous ? 1.0 : 0.05 + 0.05 * exp(-signal)
    step = ambiguous ? 0.0 : max(0.05, (0.20 + 0.80 * eff) * (1.0 - 0.8 * max(clamp(rho, 0.0, 0.95), 0.0)))
    views = ViewRecalibration[ViewRecalibration("fwd", 1.0, 0.0, 0.0, Float64[], count(isfinite, fwd.z), true, "ok")]
    bwd_missing || push!(views, ViewRecalibration("bwd", 1.0, 0.0, 0.0, Float64[], count(isfinite, bwd.z), true, "ok"))
    candidates = CandidateNView[]
    for n in n_min:n_max
        gcv, chosen_centers, is_ambiguous = _candidate_geometry_scan(ts, prof, n)
        chosen_centers = isempty(chosen_centers) ? collect(range(span[1], span[2]; length=n)) : chosen_centers
        lobes = LobePatch[]
        for (i, t0) in enumerate(chosen_centers)
            raw_fwd, pts = _patch(fwd, ax, t0, 0.0, reg.grid_half_nm, reg.grid_step_nm)
            raw_bwd = bwd_missing ? copy(raw_fwd) : _patch(bwd::ViewData, ax, t0, 0.0, reg.grid_half_nm, reg.grid_step_nm)[1]
            res_fwd = _backbone_residual(raw_fwd, pts, cfg.sigma_par_nm, cfg.sigma_perp_nm)
            res_bwd = _backbone_residual(raw_bwd, pts, cfg.sigma_par_nm, cfg.sigma_perp_nm)
            push!(lobes, LobePatch(i, ax.cx + t0 * ax.ax, ax.cy + t0 * ax.ay, t0, 0.0,
                _finite_max(raw_fwd), 0.20, 0.14, 1.0,
                Dict("fwd" => res_fwd, "bwd" => res_bwd), Dict("fwd" => raw_fwd, "bwd" => raw_bwd)))
        end
        patch_pts = NamedTuple{(:t,:u),Tuple{Float64,Float64}}[(; t=t, u=0.0) for t in chosen_centers]
        candidate_ambiguous = ambiguous || is_ambiguous
        final_gcv = candidate_ambiguous ? 1.0 : gcv
        push!(candidates, CandidateNView(n, signal, true, true, candidate_ambiguous ? "ambiguous" : "ok", view_count, bwd_missing, views,
            final_gcv, final_gcv, rho, eff, view_count * eff, 8 * n, 12 * n, lobes, patch_pts))
    end
    return CandidateViewReport(candidates, SkippedCandidate[], view_count, bwd_missing, ["ok"])
end

function candidate_report_summary(report::CandidateViewReport)
    cands = collect(report.candidates)
    total = length(cands)
    finite = count(c -> isfinite(c.joint_gcv) && isfinite(c.n), cands)
    valid = count(c -> getproperty(c, :valid), cands)
    finite_valid = count(c -> _finite_valid_count_candidate(c), cands)
    best = isempty(cands) ? nothing : first(sort(cands; by=c -> (c.joint_gcv, c.n)))
    return (; total, finite, valid, finite_valid,
        best_n = best === nothing ? missing : best.n,
        best_gcv = best === nothing ? NaN : best.joint_gcv,
        ambiguous = total > 1 && !isempty(cands) && all(c -> isapprox(c.joint_gcv, first(cands).joint_gcv; atol=1e-12, rtol=0.0), cands))
end

selected_candidate(report::CandidateViewReport) = first(sort(report.candidates; by=c -> (c.joint_gcv, c.n)))

_normalize_type_result(raw) = Main.JointProxyTypePosterior.TypePosteriorResult(
    raw.lobe_marginals, raw.map_sequence, raw.log_evidence, raw.global_state_posterior, raw.proxy_family_sensitivity)

_abstained_type_result(n::Int) = Main.JointProxyTypePosterior.TypePosteriorResult(
    fill(0.5, n, 2), fill(0, n), 0.0, Main.JointProxyTypePosterior.TypePosteriorGlobalState[], Pair{String,Float64}[])

function type_result_for_candidate(cand, reg::ProxyEnsemble; swap_type_mapping::Bool=false, fallback_n::Int=0)
    use_reg = swap_type_mapping ? swap_registry(reg) : reg
    lobes = [TypePosteriorLobeEvidence(l.residual_patches) for l in cand.lobes]
    try
        return _normalize_type_result(infer_type_posterior(lobes, use_reg; rho=cand.residual_corr, effective_factor=cand.effective_view_factor))
    catch err
        err isa ArgumentError || rethrow()
        return _abstained_type_result(max(length(lobes), fallback_n))
    end
end

function selected_type_result(case::SyntheticCase, report::CandidateViewReport, reg::ProxyEnsemble; swap_type_mapping::Bool=false)
    return type_result_for_candidate(selected_candidate(report), reg; swap_type_mapping=swap_type_mapping, fallback_n=case.truth.N)
end

function write_calibration_toml(path::AbstractString, bundle_hashes::NamedTuple, count_model, type_model; type_adapter::String=oracle_geometry_adapter_name())
    open(path, "w") do io
        println(io, "[provenance]")
        println(io, "config_sha256 = \"$(bundle_hashes.config_hash)\"")
        println(io, "source_sha256 = \"$(bundle_hashes.source_hash)\"")
        println(io, "payload_sha256 = \"$(bundle_hashes.payload_hash)\"")
        println(io, "\n[count]")
        println(io, "temperature = $(count_model.temperature)")
        println(io, "confidence_threshold = $(count_model.confidence_threshold)")
        println(io, "\n[type]")
        println(io, "temperature = $(type_model.temperature)")
        println(io, "confidence_threshold = $(type_model.confidence_threshold)")
        println(io, "\n[type.provenance]")
        println(io, "adapter = \"$type_adapter\"")
    end
    return path
end

function write_metrics_tsv(path::AbstractString, rows)
    mkpath(dirname(path))
    fields = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for row in rows
            println(io, join((get(row, f, "") for f in fields), '\t'))
        end
    end
    return path
end

