module MaskedPreprocessingFixtures

using Statistics, LinearAlgebra, TOML
import GaussianFit2D
using STMSXMIO: SXMImage, SXMChannel

export synthetic_cases, native_reference, study_options, support_metrics, run_signal_study

"""Load published explicit diagnostic settings; no fallback or hidden defaults."""
study_options() = TOML.parsefile(joinpath(@__DIR__, "..", "..", "config",
                                        "masked_robust_preprocessing.toml"))

_compact(t) = abs(t) < 1 ? cospi(t / 2)^2 : 0.0

"""
Five predetermined 73×97 images; no parameter sweep, random seed, raw file,
count selection, molecular fit, class label, or automatic exclusion inference.

The compact signal contains four bumps and a signed localized transverse ramp.
Its rectangular exclusion strictly contains its support. The Gaussian case has
noncompact tails and an intentionally narrow exclusion. Truth is used only to
construct these synthetic inputs and evaluate errors, never by either method.
"""
function synthetic_cases()
    xs = collect(range(0.0, 8.0; length=97))
    ys = collect(range(0.0, 6.0; length=73))
    plane = [0.7 + 0.055x - 0.042y for y in ys, x in xs]
    row_offsets = [0.025sinpi(2y / 6) + 0.012cospi(6y / 6) for y in ys]
    background = plane .+ row_offsets
    centers = (1.8, 2.6, 3.4, 4.2)
    even = [sum(0.22 * _compact((x-4)/0.55) * _compact((y-c)/0.55)
                for c in centers) for y in ys, x in xs]
    gradient_template = [((x-4)/1.2) * _compact((x-4)/1.2) *
                         _compact((y-3)/2.0) for y in ys, x in xs]
    gradient_amplitude_nm = 0.08
    foreground = BitMatrix([2.6 <= x <= 5.4 && 0.8 <= y <= 5.2
                            for y in ys, x in xs])
    # A fixed corner sets ONLY the arbitrary additive level in evaluation.
    # The method-specific offset uses this SAME support for both methods.
    anchor = BitMatrix([x <= 0.75 && y <= 0.5 for y in ys, x in xs])
    structured = BitMatrix([!((y < 1.5 && x >= 3.0) ||
                              (y > 4.5 && x <= 3.0) ||
                              (2.0 <= y <= 4.0 && 6.0 <= x <= 6.5) ||
                              (2.9 <= y <= 3.3 && 3.6 <= x <= 4.0))
                            for y in ys, x in xs])
    descriptions = (
        ("compact_complete", "Compact signal, full observation, covering exclusion"),
        ("compact_structured_missing", "Same compact signal, asymmetric blocks/stripe/hole missing"),
        ("compact_background_outliers", "Same compact signal, fixed positive background outlier patch"),
        ("gaussian_leaky_exclusion", "Gaussian tails and a deliberately too-narrow exclusion"),
        ("compact_background_absent_rows", "Structured missingness plus foreground-only rows and one empty row"),
    )
    cases = NamedTuple[]
    for (name, description) in descriptions
        local_even, local_gradient = copy(even), copy(gradient_template)
        supplied_background = .!foreground
        observed = trues(size(background))
        contamination = zeros(size(background))
        if name in ("compact_structured_missing", "compact_background_absent_rows")
            observed .= structured
        elseif name == "compact_background_outliers"
            for iy in eachindex(ys), ix in eachindex(xs)
                if 1.6 <= ys[iy] <= 4.3 && 6.0 <= xs[ix] <= 6.6
                    contamination[iy, ix] = 1.5
                end
            end
        elseif name == "gaussian_leaky_exclusion"
            local_even = [sum(0.22exp(-0.5*((x-4)/0.75)^2 - 0.5*((y-c)/0.40)^2)
                              for c in centers) for y in ys, x in xs]
            local_gradient = [((x-4)/1.1) * exp(-0.5*((x-4)/1.1)^2 -
                                               0.5*((y-3)/1.5)^2)
                              for y in ys, x in xs]
            supplied_background = BitMatrix([!(3.6 <= x <= 4.4 && 1.1 <= y <= 4.9)
                                              for y in ys, x in xs])
        end
        absent_rows = Int[]
        if name == "compact_background_absent_rows"
            absent_rows = findall(y -> 2.75 <= y <= 3.25, ys)
            for iy in absent_rows
                observed[iy, supplied_background[iy, :]] .= false
            end
            observed[68, :] .= false
        end
        signal = local_even .+ gradient_amplitude_nm .* local_gradient
        raw = background .+ contamination .+ signal
        null_raw = background .+ contamination
        # NaN and signed Inf encode genuinely unobserved values. Finite sentinels
        # are NOT hidden data: the estimator has no separate observation mask.
        for (k, idx) in enumerate(findall(.!observed))
            missing_value = k % 3 == 0 ? Inf : k % 3 == 1 ? NaN : -Inf
            raw[idx] = null_raw[idx] = missing_value
        end
        push!(cases, (; name, description, xs=copy(xs), ys=copy(ys),
                     plane=copy(plane), row_offsets=copy(row_offsets),
                     background=copy(background), signal, even=local_even,
                     gradient_template=local_gradient, gradient_amplitude_nm,
                     contamination, raw, null_raw, observed, supplied_background,
                     foreground=copy(foreground), background_evaluation=.!foreground,
                     anchor=copy(anchor), absent_rows))
    end
    return cases
end

"""Call the actual native entrypoint with native shared types; never read an SXM."""
function native_reference(case, z_nm=case.raw)
    ny, nx = size(z_nm)
    # Exercise native metre-to-nanometre conversion. SXMChannel stores [y,x].
    ch = SXMChannel("Z", "m", "fwd", Matrix{Float64}(z_nm .* 1e-9))
    img = SXMImage("synthetic://$(case.name)", Dict{String,String}(), nx, ny,
                   (last(case.xs), last(case.ys)), (0.0, 0.0), [ch])
    cfg = GaussianFit2D.PatternConfig(stride=1, flatten="plane+rows",
                                     smooth_radius_px=0)
    xs, ys, filled, flattened, smoothed, unit, _ =
        GaussianFit2D.preprocess_channel(img, ch, cfg)
    # Native median-filled values influence its estimator, but are not observed
    # evidence. Restore raw validity before any common-support metric.
    observed = BitMatrix(isfinite.(z_nm))
    corrected = copy(flattened)
    corrected[.!observed] .= NaN
    return (; xs, ys, filled, flattened, smoothed, unit, observed,
             corrected, valid=observed .& isfinite.(corrected),
             status=:NATIVE, reason=:median_imputation_global_row_level)
end

_safe_mean(v) = isempty(v) ? NaN : mean(v)
_safe_rmse(v) = isempty(v) ? NaN : sqrt(mean(abs2, v))

"""
Evaluate one method on an explicitly supplied common observed support. No metric
silently chooses its own easier pixels. Levels are aligned with one scalar per
image on the fixed common anchor; no ramp, row term, gain, or smoothing is removed.

`background_*` measures error in the estimated nuisance surface, after scalar
alignment. `background_residual_*` also contains injected contamination.
`foreground_*` measures total recovered-signal error. Paired response compares
signal-present minus signal-absent runs with identical background/missingness.
The transverse metric projects paired error on the centered known local ramp;
it is a signal-preservation diagnostic, not a chemical score. `metric_status`
and `metric_reason` describe evaluation availability, separately from estimator
status. No support or no fixed anchor blocks the metrics; there is no alternate
anchor. Empty background/foreground subsets and unavailable transverse contrast
are explicit. A calculable metric on less than all observed support is PARTIAL,
not a claim of full-frame recovery; raw missing counts remain separate.
"""
function support_metrics(case, corrected, null_corrected, support)
    support = BitMatrix(support)
    all((.!support) .| (case.observed .& isfinite.(corrected) .&
                       isfinite.(null_corrected))) || error("Metric support is not jointly observed/available")
    a = support .& case.anchor
    bg = support .& case.background_evaluation
    fg = support .& case.foreground
    # Observable-data gauge only: injected truth is NOT used to set the level.
    # A Gaussian tail reaching this corner biases the absolute reference level;
    # the same fixed anchor/support is used by every compared method.
    level = _safe_mean(corrected[a])
    response = corrected .- null_corrected
    response_level = _safe_mean(response[a])
    foreground_error = (corrected .- level .- case.signal)[fg]
    # inferred nuisance = observed_raw - aligned corrected signal
    nuisance_error = (case.raw .- corrected .+ level .- case.background)[bg]
    residual_error = (corrected .- level .- case.signal)[bg]
    paired_error = (response .- response_level .- case.signal)[fg]
    template = case.gradient_template[fg]
    centered_template = isempty(template) ? template : template .- mean(template)
    denom = sum(abs2, centered_template)
    transverse_error = denom > 0 ? dot(paired_error, centered_template) / denom : NaN
    observed_fg = count(case.observed .& case.foreground)
    n_support, n_anchor, n_bg, n_fg = count(support), count(a), count(bg), count(fg)
    missing_observed = count(case.observed .& .!support)
    metric_status, metric_reason = if n_support == 0
        (:UNAVAILABLE, :no_support)
    elseif n_anchor == 0
        (:UNAVAILABLE, :anchor_unavailable)
    elseif n_bg == 0 && n_fg == 0
        (:UNAVAILABLE, :no_evaluation_support)
    elseif n_bg == 0
        (:PARTIAL, :background_unavailable)
    elseif n_fg == 0
        (:PARTIAL, :foreground_unavailable)
    elseif denom == 0
        (:PARTIAL, :transverse_contrast_unavailable)
    elseif missing_observed > 0
        (:PARTIAL, :incomplete_observed_support)
    else
        (:OK, :ok)
    end
    return (; metric_status, metric_reason,
             support_pixels=n_support, anchor_pixels=n_anchor,
             background_pixels=n_bg, foreground_pixels=n_fg,
             observed_pixels=count(case.observed), observed_foreground_pixels=observed_fg,
             raw_missing_pixels=count(.!case.observed),
             unavailable_observed_pixels=count(case.observed .& .!support),
             unavailable_observed_foreground_pixels=count(case.observed .& case.foreground .& .!support),
             coverage=count(support)/count(case.observed),
             foreground_coverage=observed_fg == 0 ? NaN : count(fg)/observed_fg,
             level_offset_nm=level, paired_level_offset_nm=response_level,
             background_bias_nm=_safe_mean(nuisance_error),
             background_rmse_nm=_safe_rmse(nuisance_error),
             background_residual_bias_nm=_safe_mean(residual_error),
             background_residual_rmse_nm=_safe_rmse(residual_error),
             foreground_bias_nm=_safe_mean(foreground_error),
             foreground_rmse_nm=_safe_rmse(foreground_error),
             paired_foreground_bias_nm=_safe_mean(paired_error),
             paired_foreground_rmse_nm=_safe_rmse(paired_error),
             injected_transverse_amplitude_nm=case.gradient_amplitude_nm,
             recovered_transverse_amplitude_nm=case.gradient_amplitude_nm + transverse_error)
end


"""
    run_signal_study(estimator, opts; cases=synthetic_cases())

Run four fixed variants, each with/without injected signal. `estimator` is the
prototype's `masked_background` function; the library never supplies a fake
engine. Return `(;cases,outcomes,metrics)`. Each outcome contains case (name),
method, signal_run, null_run. Metrics use own support, the all-method intersection,
and native-versus-each-prototype intersections. A failed method is retained with
empty availability and NaN errors; it cannot force a successful-looking subset.
No across-case mean or automatic winner is calculated.
"""
function run_signal_study(estimator, opts; cases=synthetic_cases())
    outcomes, metrics = NamedTuple[], NamedTuple[]
    for case in cases
        runs = NamedTuple[]
        for method in ("native_reference", "finite_only_ols", "guarded_ols", "guarded_huber")
            if method == "native_reference"
                signal_run = native_reference(case)
                null_run = native_reference(case, case.null_raw)
            else
                mask = method == "finite_only_ols" ? trues(size(case.raw)) : case.supplied_background
                loss = method == "guarded_huber" ? :huber : :ols
                signal_run = estimator(case.xs, case.ys, case.raw, mask, opts; loss)
                null_run = estimator(case.xs, case.ys, case.null_raw, mask, opts; loss)
            end
            push!(runs, (; case=case.name, method, signal_run, null_run))
        end
        append!(outcomes, runs)
        supports = [case.observed .& r.signal_run.valid .& r.null_run.valid for r in runs]
        common = reduce((a,b) -> a .& b, supports)
        function record(r, comparison, support)
            m = support_metrics(case, r.signal_run.corrected, r.null_run.corrected, support)
            push!(metrics, (; case=case.name, comparison, method=r.method,
                signal_status=r.signal_run.status, signal_reason=r.signal_run.reason,
                null_status=r.null_run.status, null_reason=r.null_run.reason,
                signal_available_pixels=count(r.signal_run.valid),
                null_available_pixels=count(r.null_run.valid), m...))
        end
        for (r, s) in zip(runs, supports)
            record(r, "own", s)
            record(r, "all_common", common)
        end
        for k in 2:length(runs)
            comparison = "native_vs_" * runs[k].method
            shared = supports[1] .& supports[k]
            record(runs[1], comparison, shared)
            record(runs[k], comparison, shared)
        end
    end
    return (; cases, outcomes, metrics)
end

end # module
