# Synthetic only. No SXM read, molecular fit, count/chemical selection, or grade.
# julia --startup-file=no --threads=1 --project=. test/test_masked_preprocessing_signal.jl
# Optional: --outdir <new-directory> saves complete evidence before assertions.
using Test, Statistics, LinearAlgebra, TOML
include(joinpath(@__DIR__, "lib", "masked_robust_preprocessing.jl"))
include(joinpath(@__DIR__, "lib", "masked_preprocessing_fixtures.jl"))
using .MaskedRobustPreprocessing, .MaskedPreprocessingFixtures
BLAS.set_num_threads(1)

function write_named_rows(path, rows)
    isempty(rows) && error("Cannot write a schema-free empty table")
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), '\t'))
        for row in rows
            keys(row) == keys(first(rows)) || error("Inconsistent evidence schema")
            println(io, join(string.(values(row)), '\t'))
        end
    end
end

function save_evidence(study, path, settings)
    (ispath(path) || islink(path)) && error("Refusing to overwrite existing synthetic evidence: $path")
    mkpath(path)
    open(joinpath(path, "settings.toml"), "w") do io
        TOML.print(io, settings)
    end
    write_named_rows(joinpath(path, "metrics.tsv"), study.metrics)
    statuses, rows = NamedTuple[], NamedTuple[]
    for outcome in study.outcomes, (kind, r) in (("signal", outcome.signal_run), ("null", outcome.null_run))
        coefficients = get(r, :coefficients, (intercept_nm=NaN, slope_x=NaN, slope_y=NaN))
        push!(statuses, (; case=outcome.case, method=outcome.method, input=kind,
            status=r.status, reason=r.reason, observed=count(r.observed), valid=count(r.valid),
            converged=get(r, :converged, missing), iterations=get(r, :iterations, missing),
            robust_scale_nm=get(r, :robust_scale, NaN),
            initial_objective_nm2=get(r, :initial_objective_nm2, NaN),
            objective_nm2=get(r, :objective_nm2, NaN),
            stationarity_inf_nm=get(r, :stationarity_inf_nm, NaN),
            design_rank=get(r, :design_rank, missing),
            background_pixels=hasproperty(r, :background_mask) ? count(r.background_mask) : missing,
            coefficient_change=get(r, :coefficient_change, NaN),
            weight_change=get(r, :weight_change, NaN),
            scale_floor_active=get(r, :scale_floor_active, missing), coefficients...))
        if hasproperty(r, :row_status)
            for iy in eachindex(r.row_status)
                push!(rows, (; case=outcome.case, method=outcome.method, input=kind, row=iy,
                    observed=count(r.observed[iy, :]), available=count(r.valid[iy, :]),
                    background_pixels=r.row_background_counts[iy], status=r.row_status[iy],
                    row_offset_nm=r.row_offsets[iy]))
            end
        end
    end
    write_named_rows(joinpath(path, "statuses.tsv"), statuses)
    write_named_rows(joinpath(path, "row_statuses.tsv"), rows)
    open(joinpath(path, "pixels.tsv"), "w") do io
        methods = ("native_reference", "finite_only_ols", "guarded_ols", "guarded_huber")
        base = ("case", "row", "col", "x_nm", "y_nm", "observed", "supplied_background",
                "foreground_evaluation", "level_anchor", "raw_nm", "null_raw_nm",
                "true_background_nm", "true_signal_nm", "true_even_nm", "transverse_template",
                "injected_transverse_amplitude_nm", "contamination_nm")
        extra = ["$(m)_$(kind)_nm" for m in methods for kind in ("signal", "null")]
        println(io, join(vcat(collect(base), extra), '\t'))
        for c in study.cases
            runs = [only(filter(o -> o.case == c.name && o.method == m, study.outcomes)) for m in methods]
            for ix in eachindex(c.xs), iy in eachindex(c.ys)
                fixed = (c.name, iy, ix, c.xs[ix], c.ys[iy], c.observed[iy,ix],
                    c.supplied_background[iy,ix], c.foreground[iy,ix], c.anchor[iy,ix],
                    c.raw[iy,ix], c.null_raw[iy,ix], c.background[iy,ix],
                    c.signal[iy,ix], c.even[iy,ix], c.gradient_template[iy,ix],
                    c.gradient_amplitude_nm, c.contamination[iy,ix])
                outputs = [r.corrected[iy,ix] for o in runs for r in (o.signal_run, o.null_run)]
                println(io, join(vcat(collect(fixed), outputs), '\t'))
            end
        end
    end
    open(joinpath(path, "notes.txt"), "w") do io
        println(io, "Julia $VERSION; threads=$(Threads.nthreads()); BLAS=$(BLAS.get_num_threads()).")
        println(io, "Synthetic evidence only; fixed five cases and four variants, no winner or across-case mean.")
        println(io, "Native: actual GaussianFit2D.preprocess_channel, native SXM types, stride=1, plane+rows, smooth=0.")
        println(io, "Original NaN/Inf validity restored after native median filling. No filled sample is observed evidence.")
        println(io, "Observable level: mean of corrected samples at x<=0.75,y<=0.5 on the same common support.")
        println(io, "Only one scalar is removed per image. Truth is NOT used for this gauge alignment.")
        println(io, "The Gaussian tail can bias that fixed reference level; local ramp contrast is offset invariant.")
        println(io, "background_* evaluates estimated nuisance versus truth; background_residual_* also includes contamination.")
        println(io, "foreground_* is absolute signal error after level alignment; paired_* is injection response, not recovery.")
        println(io, "Guarded masks are trusted input, not inferred. Leaky exclusions remain contaminated.")
        println(io, "Unsupported rows remain NaN. Coverage/own/native-pair common metrics prevent gains by silently dropping rows.")
        println(io, "metric_status/reason are separate from estimator status; no support or no fixed anchor blocks metrics without fallback.")
        println(io, "Partial observed-support metrics are NOT full raw coverage. The background-absent fixture loses206 observed foreground pixels.")
        println(io, "On the structured-missing fixture guarded_huber has WORSE absolute background/foreground RMSE than guarded_ols.")
        println(io, "Sequential plane then rows can bias x slope with unequal per-row x support; y tilt/row terms are not unique.")
        println(io, "Fixed initial OLS MAD Huber scale is not calibrated noise. No controls were tuned after seeing outcomes.")
        println(io, "Evidence is saved before tests, so presence of files does not imply assertions passed. Read the test log.")
        println(io, "settings.toml retains the actually-used config/masked_robust_preprocessing.toml settings snapshot.")
        TOML.print(io, settings)
    end
end

if ARGS == ["--help"]
    println("Synthetic signal preservation tests; optional --outdir <new-directory>.")
    exit(0)
end
outdir = isempty(ARGS) ? nothing : length(ARGS) == 2 && first(ARGS) == "--outdir" ? last(ARGS) :
    error("Usage: test_masked_preprocessing_signal.jl [--outdir <new-directory>]")
outdir === nothing || !(ispath(outdir) || islink(outdir)) ||
    error("Refusing existing synthetic evidence path: $outdir")
settings = study_options()
opts = read_options(settings)
cases = synthetic_cases()
before = deepcopy(cases)
study = run_signal_study(masked_background, opts; cases)
outdir === nothing || save_evidence(study, outdir, settings)
outcome(name, method) = only(filter(o -> o.case == name && o.method == method, study.outcomes))
metric(name, method; comparison="all_common") =
    only(filter(m -> m.case == name && m.method == method && m.comparison == comparison, study.metrics))
const TOL = 1e-12 # Arithmetic equality only, not a scientific noise/quality cutoff.

@testset "Fixed synthetic family and actual native reference" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    @test Threads.nthreads() == 1
    @test BLAS.get_num_threads() == 1
    @test settings["masked_robust_preprocessing"] == Dict(
        "huber_delta"=>1.345, "maxiter"=>50, "coefficient_rtol"=>1e-8,
        "weight_atol"=>1e-8, "scale_floor_nm"=>1e-9, "rank_rtol"=>1e-12,
        "min_background_pixels"=>200, "min_background_fraction"=>0.05,
        "min_row_background_pixels"=>20)
    mktempdir() do dir
        @test_throws ErrorException save_evidence(study, dir, settings)
        dangling = joinpath(dir, "dangling")
        symlink(joinpath(dir, "absent"), dangling)
        @test islink(dangling) && !ispath(dangling)
        @test_throws ErrorException save_evidence(study, dangling, settings)
        @test !ispath(joinpath(dir, "absent"))
    end
    @test length(cases) == 5
    @test length(study.outcomes) == 20
    @test length(study.metrics) == 70
    @test isequal(cases, before)
    @test isequal(cases, synthetic_cases())
    for c in cases
        @test size(c.raw) == (73, 97)
        @test c.observed == isfinite.(c.raw) == isfinite.(c.null_raw)
        @test c.anchor .& c.foreground == falses(size(c.raw))
        @test all(c.observed[c.anchor])
        @test minimum(c.gradient_template) < 0 < maximum(c.gradient_template)
        @test c.signal[c.observed] ≈ (c.raw .- c.null_raw)[c.observed] atol=TOL
        if c.name != "gaussian_leaky_exclusion"
            @test all(iszero, c.signal[c.supplied_background])
        else
            @test any(!iszero, c.signal[c.supplied_background])
        end
        if any(.!c.observed)
            @test any(isnan, c.raw)
            @test any(==(Inf), c.raw)
            @test any(==(-Inf), c.raw)
        end
        native = outcome(c.name, "native_reference").signal_run
        @test native.xs == c.xs && native.ys == c.ys
        @test native.unit == "nm"
        @test native.filled[c.observed] ≈ c.raw[c.observed] atol=TOL
        @test all(isfinite, native.filled)
        @test native.smoothed == native.flattened
        @test native.valid == c.observed
        @test all(isnan, native.corrected[.!c.observed])
        @test all(isapprox(median(native.flattened[iy,:]), median(native.flattened[1,:]); atol=TOL)
                  for iy in axes(c.raw,1))
        if any(.!c.observed)
            @test all(isapprox(v, median(c.raw[c.observed]); atol=TOL) for v in native.filled[.!c.observed])
        else
            # Actual native and finite-only OLS are the same transform up to a
            # constant on complete images; a zero-level change is not a gain.
            finite = outcome(c.name, "finite_only_ols").signal_run
            native_level = mean(native.corrected[c.anchor])
            finite_level = mean(finite.corrected[c.anchor])
            @test native.corrected .- native_level ≈ finite.corrected .- finite_level atol=TOL
        end
        for method in ("finite_only_ols", "guarded_ols", "guarded_huber")
            o = outcome(c.name, method)
            expected_mask = method == "finite_only_ols" ? c.observed : c.observed .& c.supplied_background
            for r in (o.signal_run, o.null_run)
                @test r.observed == c.observed
                @test r.background_mask == expected_mask
                @test all((.!r.valid) .| c.observed)
                @test all(isnan, r.corrected[.!r.valid])
                @test all(isnan, r.background[.!r.valid])
                @test r.status in (:OK, :PARTIAL, :UNAVAILABLE, :NONCONVERGED)
                @test r.iterations <= opts.maxiter
                if !r.converged
                    @test !any(r.valid)
                    @test all(isnan, r.corrected)
                end
            end
        end
    end
end

@testset "Observable gauge and identical comparison denominators" begin
    c = first(cases)
    o = outcome(c.name, "native_reference")
    support = c.observed
    m = support_metrics(c, o.signal_run.corrected, o.null_run.corrected, support)
    shifted = support_metrics(c, o.signal_run.corrected .+ 9.0, o.null_run.corrected .+ 3.0, support)
    altered_truth = support_metrics(merge(c, (; signal=c.signal .+ 10.0)),
                                   o.signal_run.corrected, o.null_run.corrected, support)
    @test shifted.level_offset_nm ≈ m.level_offset_nm + 9 atol=TOL
    @test shifted.paired_level_offset_nm ≈ m.paired_level_offset_nm + 6 atol=TOL
    @test altered_truth.level_offset_nm == m.level_offset_nm
    @test altered_truth.paired_level_offset_nm == m.paired_level_offset_nm
    @test altered_truth.foreground_bias_nm ≈ m.foreground_bias_nm - 10 atol=TOL
    for key in (:background_rmse_nm, :foreground_rmse_nm, :paired_foreground_rmse_nm,
                :recovered_transverse_amplitude_nm)
        @test getproperty(shifted,key) ≈ getproperty(m,key) atol=TOL
    end
    @test_throws ErrorException support_metrics(c, fill(NaN,size(c.raw)), o.null_run.corrected, support)
    empty = support_metrics(c, fill(NaN,size(c.raw)), fill(NaN,size(c.raw)), falses(size(c.raw)))
    @test empty.coverage == 0.0
    @test empty.metric_status == :UNAVAILABLE && empty.metric_reason == :no_support
    @test empty.anchor_pixels == empty.background_pixels == empty.foreground_pixels == 0
    @test empty.unavailable_observed_pixels == count(c.observed)
    @test isnan(empty.foreground_rmse_nm) && isnan(empty.background_rmse_nm)
    @test isnan(empty.level_offset_nm) && isnan(empty.recovered_transverse_amplitude_nm)
    no_anchor = support_metrics(c, o.signal_run.corrected, o.null_run.corrected,
                                support .& .!c.anchor)
    @test no_anchor.metric_status == :UNAVAILABLE
    @test no_anchor.metric_reason == :anchor_unavailable
    @test no_anchor.anchor_pixels == 0
    @test no_anchor.support_pixels == count(support) - count(c.anchor)
    @test no_anchor.background_pixels > 0 && no_anchor.foreground_pixels > 0
    @test isnan(no_anchor.level_offset_nm) && isnan(no_anchor.paired_level_offset_nm)
    @test isnan(no_anchor.background_rmse_nm) && isnan(no_anchor.foreground_rmse_nm)
    @test isnan(no_anchor.paired_foreground_rmse_nm) && isnan(no_anchor.recovered_transverse_amplitude_nm)
    no_foreground = support_metrics(c, o.signal_run.corrected, o.null_run.corrected,
                                    support .& .!c.foreground)
    @test no_foreground.metric_status == :PARTIAL
    @test no_foreground.metric_reason == :foreground_unavailable
    @test no_foreground.foreground_pixels == 0 && no_foreground.background_pixels > 0
    @test isfinite(no_foreground.background_rmse_nm) && isnan(no_foreground.foreground_rmse_nm)
    no_background = support_metrics(c, o.signal_run.corrected, o.null_run.corrected,
                                    support .& c.foreground)
    @test no_background.metric_status == :UNAVAILABLE && no_background.metric_reason == :anchor_unavailable
    @test no_background.background_pixels == 0 && no_background.foreground_pixels > 0
    @test isnan(no_background.background_rmse_nm)
    for c in cases, comparison in unique(m.comparison for m in study.metrics if m.case == c.name)
        records = filter(m -> m.case == c.name && m.comparison == comparison, study.metrics)
        for m in records
            own = outcome(c.name,m.method)
            chosen = comparison == "own" ? [own] :
                [outcome(c.name,k.method) for k in records]
            common = copy(c.observed)
            for r in chosen; common .&= r.signal_run.valid .& r.null_run.valid; end
            @test m.support_pixels == count(common)
            @test m.foreground_pixels == count(common .& c.foreground)
            @test m.background_pixels == count(common .& c.background_evaluation)
            @test m.raw_missing_pixels == count(.!c.observed)
            @test m.unavailable_observed_pixels + m.support_pixels == count(c.observed)
            @test m.metric_status == (m.unavailable_observed_pixels > 0 ? :PARTIAL : :OK)
            @test m.metric_reason == (m.unavailable_observed_pixels > 0 ? :incomplete_observed_support : :ok)
            @test m.unavailable_observed_foreground_pixels + m.foreground_pixels == count(c.observed .& c.foreground)
            if count(common .& c.anchor) > 0
                @test m.level_offset_nm == mean(own.signal_run.corrected[common .& c.anchor])
                err = (own.signal_run.corrected .- m.level_offset_nm .- c.signal)[common .& c.foreground]
                @test m.foreground_rmse_nm == sqrt(sum(abs2,err)/length(err))
            end
        end
    end
end

@testset "Preservation is conditional, not background recovery" begin
    for c in cases
        c.name == "gaussian_leaky_exclusion" && continue
        for method in ("guarded_ols", "guarded_huber")
            o = outcome(c.name, method)
            # Excluded injection cannot change the background samples/estimator.
            @test o.signal_run.coefficients == o.null_run.coefficients
            @test o.signal_run.valid == o.null_run.valid
            if any(o.signal_run.valid)
                m = metric(c.name, method; comparison="own")
                @test m.paired_foreground_rmse_nm < TOL
                @test m.recovered_transverse_amplitude_nm ≈ c.gradient_amplitude_nm atol=TOL
            end
        end
    end
    for method in ("guarded_ols", "guarded_huber")
        c = first(cases)
        o = outcome(c.name, method)
        @test o.signal_run.status == :OK
        @test metric(c.name,method).coverage == 1.0
        @test metric(c.name,method).foreground_rmse_nm < TOL
        @test metric(c.name,method).background_rmse_nm < TOL
    end
    @test metric("compact_complete", "native_reference").paired_foreground_rmse_nm > TOL
    @test metric("compact_structured_missing", "guarded_ols").background_rmse_nm > TOL
    @test metric("compact_structured_missing", "guarded_ols").foreground_rmse_nm > TOL
    # Leakage is a retained adverse result; robust plane fitting is not a
    # guarantee that localized chemical-looking gradients or tails survive.
    for method in ("finite_only_ols", "guarded_ols", "guarded_huber")
        m = metric("gaussian_leaky_exclusion", method; comparison="own")
        if m.support_pixels > 0
            @test m.paired_foreground_rmse_nm > TOL
            @test abs(m.recovered_transverse_amplitude_nm - 0.08) > TOL
        end
    end
end

@testset "Missing observations and background-absent rows stay visible" begin
    c = only(filter(c -> c.name == "compact_background_absent_rows", cases))
    lost_foreground = count(c.observed[c.absent_rows,:])
    @test lost_foreground > 0
    @test !any(c.observed[68,:])
    @test all(c.foreground[c.absent_rows,:][c.observed[c.absent_rows,:]])
    @test metric(c.name,"native_reference";comparison="own").coverage == 1.0
    @test metric(c.name,"finite_only_ols";comparison="own").coverage == 1.0
    for method in ("guarded_ols", "guarded_huber")
        r = outcome(c.name,method).signal_run
        @test r.status == :PARTIAL
        @test r.reason == :insufficient_row_background
        @test all(==(:insufficient_background), r.row_status[c.absent_rows])
        @test all(isnan, r.corrected[c.absent_rows,:])
        @test all(isnan, r.row_offsets[c.absent_rows])
        m = metric(c.name,method;comparison="own")
        @test m.unavailable_observed_pixels == lost_foreground
        @test m.unavailable_observed_foreground_pixels == lost_foreground
        @test m.foreground_coverage < 1.0
        common_native = metric(c.name,"native_reference";comparison="native_vs_"*method)
        @test common_native.support_pixels == m.support_pixels
        @test common_native.foreground_pixels == m.foreground_pixels
        @test common_native.coverage < metric(c.name,"native_reference";comparison="own").coverage
    end
end

println("Synthetic per-case common-support metrics; no across-case mean:")
for m in study.metrics
    m.comparison == "all_common" || continue
    println(m.case, '\t', m.method, '\t', m.signal_status, '\t',
        "metrics=", m.metric_status, '/', m.metric_reason, '\t',
        "n=", m.support_pixels, '/', m.observed_pixels, '\t',
        "fg=", m.foreground_pixels, '/', m.observed_foreground_pixels, '\t',
        "bg_rmse_nm=", m.background_rmse_nm, '\t', "fg_rmse_nm=", m.foreground_rmse_nm, '\t',
        "paired_fg_rmse_nm=", m.paired_foreground_rmse_nm, '\t',
        "transverse_nm=", m.recovered_transverse_amplitude_nm)
end
outdir === nothing || println("Evidence saved to ", outdir)
