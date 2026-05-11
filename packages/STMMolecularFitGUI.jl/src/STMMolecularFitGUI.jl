"""
Multi-Gaussian Fitting GUI (Stipple.jl + PlotlyJS).

Usage:
    julia --project=. app.jl

Or programmatically:
    using STMMolecularFitGUI
    STMMolecularFitGUI.run()
"""

module STMMolecularFitGUI

using Stipple, Stipple.ReactiveTools
using PlotlyJS
using JSON3
using Printf
using GaussianFit1D

# Import Html constructors
const H = Stipple.Html

PEAK_COLORS = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
               "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"]

_peak_color(i::Int) = PEAK_COLORS[mod1(i, length(PEAK_COLORS))]

function _rgba(hex, alpha=0.18)
    h = strip(hex, '#')
    r = parse(Int, h[1:2], base=16)
    g = parse(Int, h[3:4], base=16)
    b = parse(Int, h[5:6], base=16)
    "rgba($r,$g,$b,$alpha)"
end

# ===========================================================================
# Fit plot (data + total fit + peaks, no residuals)
# ===========================================================================
function build_fit_plot_json(x, y, result, cfg)
    n_peaks = result.n_peaks
    popt = result.popt
    x_unit = cfg.x_unit
    asymmetric = cfg.asymmetric_edges && n_peaks >= 2
    FINE = 250

    xx = range(minimum(x), maximum(x), length=FINE)
    y0_val = popt[1]
    y_total = GaussianFit1D.predict_fit(xx, result, cfg)

    n_data = length(x)
    stride = max(1, n_data ÷ 200)
    if n_data > 200
        x_plot = x[1:stride:end]; y_plot = y[1:stride:end]
    else
        x_plot = x; y_plot = y
    end

    traces = Any[]
    # Data points
    push!(traces, Dict(:x=>x_plot, :y=>y_plot, :mode=>"markers",
        :type=>"scatter", :name=>"Data",
        :marker=>Dict(:color=>"gray", :size=>3, :opacity=>0.6),
        :hovertemplate=>"x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>"))
    # Total fit
    push!(traces, Dict(:x=>xx, :y=>y_total, :mode=>"lines",
        :type=>"scatter", :name=>"Total fit",
        :line=>Dict(:color=>"black", :width=>2),
        :hovertemplate=>"x: %{x:.3f}<br>fit: %{y:.3f}<extra></extra>"))
    # Baseline
    push!(traces, Dict(:x=>[x[1], x[end]], :y=>[y0_val, y0_val],
        :type=>"scatter", :mode=>"lines",
        :name=>"y0=$(round(y0_val, digits=3))",
        :line=>Dict(:color=>"gray", :width=>1, :dash=>"dot"),
        :hoverinfo=>"skip"))
    # Individual peaks
    centers = GaussianFit1D._params_to_centers(popt, n_peaks)
    for i in 0:(n_peaks-1)
        A = GaussianFit1D._get_amplitude(popt, i)
        sigma = GaussianFit1D._get_sigma(popt, i)
        mu = centers[i+1]
        fwhm = GaussianFit1D.FWHM_TO_SIGMA * sigma
        color = _peak_color(i+1)
        max_sig = sigma
        if asymmetric && (i == 0 || i == n_peaks - 1)
            sigma_outer = popt[end - (i == 0 ? 1 : 0)]
            max_sig = max(sigma, sigma_outer)
        end
        span = 4.0 * max_sig
        x_min_pk = max(minimum(x), mu - span)
        x_max_pk = min(maximum(x), mu + span)
        n_peak_pts = max(20, Int(ceil(FINE * (x_max_pk - x_min_pk) / (maximum(x) - minimum(x)))))
        xp = range(x_min_pk, x_max_pk, length=n_peak_pts)
        if asymmetric && (i == 0 || i == n_peaks - 1)
            popt_idx = length(popt) - 2 + (i == 0 ? 0 : 1)
            sigma_outer = popt[popt_idx+1]
            fwhm_outer = GaussianFit1D.FWHM_TO_SIGMA * sigma_outer
            z = xp .- mu
            s = i == 0 ? [zv < 0 ? sigma_outer : sigma for zv in z] : [zv < 0 ? sigma : sigma_outer for zv in z]
            g = y0_val .+ A .* exp.(-0.5 .* (z ./ s) .^ 2)
            name = "Peak $(i+1) $(round(mu, digits=2)) $x_unit"
            hover = "Peak $(i+1) (asym)<br>Center: $(round(mu, digits=3)) $x_unit<br>FWHM_in: $(round(fwhm, digits=3))<br>FWHM_out: $(round(fwhm_outer, digits=3))<extra></extra>"
        else
            g = y0_val .+ A .* exp.(-0.5 .* ((xp .- mu) ./ sigma) .^ 2)
            name = "Peak $(i+1) $(round(mu, digits=2)) $x_unit"
            hover = "Peak $(i+1)<br>Center: $(round(mu, digits=3)) $x_unit<br>Amp: $(round(A, digits=3))<br>FWHM: $(round(fwhm, digits=3)) $x_unit<extra></extra>"
        end
        fill_trace = Dict(:x=>xp, :y=>g, :mode=>"lines",
            :type=>"scatter", :showlegend=>false,
            :fill=>"tozeroy", :fillcolor=>_rgba(color, 0.15),
            :line=>Dict(:width=>0), :hoverinfo=>"skip")
        push!(traces, fill_trace)
        push!(traces, Dict(:x=>xp, :y=>g, :mode=>"lines",
            :type=>"scatter", :name=>name,
            :line=>Dict(:color=>color, :width=>1.2, :dash=>"dash"),
            :hovertemplate=>hover))
    end

    x_range = maximum(x) - minimum(x)
    x_pad = 0.02 * x_range
    r2_val = round(result.r_squared, digits=4)
    title_text = "Multi-Gaussian Fit — $n_peaks peaks, BIC=$(round(result.bic, digits=1)), R²=$r2_val"
    layout_dict = Dict(
        :xaxis => Dict(:title => "Distance ($x_unit)", :range => [x[1]-x_pad, x[end]+x_pad],
            :zeroline => false, :gridcolor => "#e0e0e0"),
        :yaxis => Dict(:title => "Intensity", :zeroline => false, :gridcolor => "#e0e0e0"),
        :height => 500, :autosize => true, :hovermode => "x unified",
        :legend => Dict(:font => Dict(:size => 9), :y => 1.0),
        :margin => Dict(:l => 60, :r => 20, :t => 50, :b => 40),
        :plot_bgcolor => "rgba(0,0,0,0)", :paper_bgcolor => "rgba(0,0,0,0)",
        :annotations => [Dict(:text => title_text, :xref => "paper", :yref => "paper",
            :x => 0.5, :y => 1.06, :showarrow => false, :font => Dict(:size => 14))])
    json_str = JSON3.write(_stringify_keys(Dict("data"=>traces, "layout"=>layout_dict,
        "config"=>Dict("displayModeBar"=>false))))
    return json_str
end

# ===========================================================================
# Residuals plot (standalone)
# ===========================================================================
function build_residuals_plot_json(x, y, result, cfg)
    y_fit = GaussianFit1D.predict_fit(x, result, cfg)
    x_unit = cfg.x_unit
    n_data = length(x)
    stride = max(1, n_data ÷ 200)
    if n_data > 200
        x_plot = x[1:stride:end]; y_plot = y[1:stride:end]
    else
        x_plot = x; y_plot = y
    end
    residuals = y .- y_fit
    res_plot = n_data > 200 ? residuals[1:stride:end] : residuals
    x_range = maximum(x) - minimum(x)
    x_pad = 0.02 * x_range
    traces = [
        Dict(:x=>x_plot, :y=>res_plot, :type=>"scatter", :mode=>"markers",
            :name=>"Residuals", :marker=>Dict(:color=>"gray", :size=>3, :opacity=>0.6),
            :hovertemplate=>"x: %{x:.3f}<br>res: %{y:.4f}<extra></extra>"),
        Dict(:x=>[x[1], x[end]], :y=>[0,0], :type=>"scatter", :mode=>"lines",
            :line=>Dict(:color=>"black", :width=>0.5),
            :showlegend=>false, :hoverinfo=>"skip"),
    ]
    layout_dict = Dict{Symbol,Any}(
        :xaxis => Dict(:title => "Distance ($x_unit)", :range => [x[1]-x_pad, x[end]+x_pad]),
        :yaxis => Dict(:title => "Residuals"),
        :height => 200, :autosize => true, :hovermode => "x unified",
        :margin => Dict(:l => 60, :r => 20, :t => 30, :b => 40))
    json_str = JSON3.write(_stringify_keys(Dict("data"=>traces, "layout"=>layout_dict,
        "config"=>Dict("displayModeBar"=>false))))
    println("    [res_plot] JSON: $(round(length(json_str)/1024, digits=1)) KB")
    return json_str
end

# ===========================================================================
# BIC comparison plot (standalone bar chart)
# ===========================================================================
function build_bic_plot_json(all_results, result, cfg)
    threshold = cfg.bic_competition_threshold

    n_list = [r.n_peaks for r in all_results]
    bic_list = [r.bic for r in all_results]
    bi = argmin(bic_list)

    bar_colors = String[]
    for (j, bic) in enumerate(bic_list)
        delta = bic - bic_list[bi]
        if j == bi
            push!(bar_colors, "rgb(220, 50, 50)")
        elseif delta <= threshold
            push!(bar_colors, "rgb(255, 190, 50)")
        else
            push!(bar_colors, "rgb(190, 195, 210)")
        end
    end

    y_min = minimum(bic_list)
    y_max = maximum(bic_list)
    y_pad = max(0.25 * (y_max - y_min), 50)

    fig = PlotlyJS.plot(PlotlyJS.bar(;
        x=n_list, y=bic_list, marker_color=bar_colors,
        text=[Printf.@sprintf("%.0f", b) for b in bic_list],
        textposition="outside", textfont=attr(size=9),
        hovertemplate="n=%{x}<br>BIC=%{y:.1f}<extra></extra>",
        name="BIC", showlegend=false), PlotlyJS.Layout(;
        xaxis_title="Number of Gaussians",
        yaxis_title="BIC",
        yaxis_range=[y_min - y_pad, y_max + y_pad],
        height=250,
        margin=Dict(:l => 60, :r => 20, :t => 40, :b => 40),
        bargap=0.3))

    current_n = result.n_peaks
    if current_n != n_list[bi]
        PlotlyJS.add_trace!(fig, PlotlyJS.scatter(;
            x=[current_n], y=[result.bic], mode="markers",
            name="Current: $current_n peaks",
            marker=attr(color="black", size=10, symbol="diamond"),
            hovertemplate="Current: n=$current_n<extra></extra>"))
    end

    json_str = _plotly_to_json(fig)
    println("    [bic_plot] JSON: $(round(length(json_str) / 1024, digits=1)) KB")
    return json_str
end

# Legacy wrapper
function build_plot_json(x, y, result, all_results, cfg)
    return build_fit_plot_json(x, y, result, cfg)
end

_stringify_keys(d::Dict) = Dict{String,Any}(string(k) => _stringify_keys(v) for (k, v) in d)
_stringify_keys(v::AbstractVector) = [_stringify_keys(x) for x in v]
_stringify_keys(x) = x

function _plotly_to_json(fig)
    p = hasproperty(fig, :plot) ? fig.plot : fig
    traces = Dict{String,Any}[]
    for t in p.data
        d = Dict{String,Any}("type" => t.kind)
        for (k, v) in t.fields
            d[string(k)] = _to_json_val(v)
        end
        push!(traces, d)
    end
    layout = Dict{String,Any}()
    for (k, v) in p.layout.fields
        k == :template && continue
        layout[string(k)] = _to_json_val(v)
    end
    return JSON3.write(Dict("data" => traces, "layout" => layout))
end

_to_json_val(v) = v
_to_json_val(d::Dict) = Dict(string(k) => _to_json_val(v) for (k, v) in d)
_to_json_val(v::AbstractRange) = collect(v)
_to_json_val(v::AbstractVector) = [_to_json_val(x) for x in v]

# ===========================================================================
# Table builder
# ===========================================================================

function build_table_rows(all_results, cfg)
    best_bic_val = minimum(r.bic for r in all_results)
    threshold = cfg.bic_competition_threshold
    rows = Dict{String,Any}[]
    for r in all_results
        delta = r.bic - best_bic_val
        push!(rows, Dict{String,Any}(
            "n_peaks" => r.n_peaks,
            "BIC" => round(r.bic, digits=1),
            "dBIC" => round(delta, digits=1),
            "R2" => round(r.r_squared, digits=5),
            "RSS" => round(r.rss, digits=6),
            "n_params" => r.n_params,
            "best" => delta == 0.0,
            "competitive" => delta > 0.0 && delta <= threshold,
            "success" => r.success,
        ))
    end
    return rows
end

# ===========================================================================
# Reactive Model
# ===========================================================================

const GUI_DEFAULTS = (
    filepath = "",
    output_dir = "",
    min_spacing = GaussianFit1D.DEFAULT_CONFIG.min_spacing,
    max_spacing = GaussianFit1D.DEFAULT_CONFIG.max_spacing,
    fwhm_min = GaussianFit1D.DEFAULT_CONFIG.fwhm_min,
    fwhm_max = GaussianFit1D.DEFAULT_CONFIG.fwhm_max,
    offset_to_zero = GaussianFit1D.DEFAULT_CONFIG.offset_to_zero,
    amplitude_min_fraction = GaussianFit1D.DEFAULT_CONFIG.amplitude_min_fraction,
    asymmetric_edges = GaussianFit1D.DEFAULT_CONFIG.asymmetric_edges,
    edge_sigma_min = GaussianFit1D.DEFAULT_CONFIG.edge_sigma_min,
    edge_sigma_max = GaussianFit1D.DEFAULT_CONFIG.edge_sigma_max,
    global_maxiter = GaussianFit1D.DEFAULT_CONFIG.global_maxiter,
    global_maxtime = GaussianFit1D.DEFAULT_CONFIG.global_maxtime,
    global_tol = GaussianFit1D.DEFAULT_CONFIG.global_tol,
    nlopt_algorithm = string(GaussianFit1D.DEFAULT_CONFIG.nlopt_algorithm),
    curve_fit_maxfev = GaussianFit1D.DEFAULT_CONFIG.curve_fit_maxfev,
    early_stop_patience = GaussianFit1D.DEFAULT_CONFIG.early_stop_patience,
    early_stop_dbic = GaussianFit1D.DEFAULT_CONFIG.early_stop_dbic,
    bic_threshold = GaussianFit1D.DEFAULT_CONFIG.bic_competition_threshold,
)

@app FitModel begin
    # --- File ---
    @in filepath = GUI_DEFAULTS.filepath
    @in output_dir = GUI_DEFAULTS.output_dir
    @in file_content = ""

    # --- Parameters ---
    @in min_spacing = GUI_DEFAULTS.min_spacing
    @in max_spacing = GUI_DEFAULTS.max_spacing
    @in fwhm_min = GUI_DEFAULTS.fwhm_min
    @in fwhm_max = GUI_DEFAULTS.fwhm_max
    @in offset_to_zero = GUI_DEFAULTS.offset_to_zero
    @in amplitude_min_fraction = GUI_DEFAULTS.amplitude_min_fraction
    @in asymmetric_edges = GUI_DEFAULTS.asymmetric_edges
    @in edge_sigma_min = GUI_DEFAULTS.edge_sigma_min
    @in edge_sigma_max = GUI_DEFAULTS.edge_sigma_max
    @in global_maxiter = GUI_DEFAULTS.global_maxiter
    @in global_maxtime = GUI_DEFAULTS.global_maxtime
    @in global_tol = GUI_DEFAULTS.global_tol
    @in nlopt_algorithm = GUI_DEFAULTS.nlopt_algorithm
    @in curve_fit_maxfev = GUI_DEFAULTS.curve_fit_maxfev
    @in early_stop_patience = GUI_DEFAULTS.early_stop_patience
    @in early_stop_dbic = GUI_DEFAULTS.early_stop_dbic
    @in bic_threshold = GUI_DEFAULTS.bic_threshold

    # --- State ---
    @out running = false
    @out fit_ready = false
    @out error_msg = ""
    @out progress_text = ""
    @out progress_percent = 0
    @out y_offset = "\u2014"

    # --- Button triggers ---
    @in runFitTrigger = false
    @in refreshTrigger = 0

    # --- Results ---
    @out plot_json = "{}"
    @out bic_plot_json = "{}"
    @out residuals_plot_json = "{}"
    @out table_json = "[]"
    @out model_labels_json = "[]"
    @out warnings_html = ""
    @in model_idx = 0
    @out best_idx = 0

    # --- Tabs ---
    @in active_tab = "plots"
    @out log_text = ""

    # --- File browser ---
    @out datalist_html = ""

    # --- File content from browser ---
    @onchange file_content begin
        if !isempty(strip(file_content))
            if isempty(_session_dir)
                _session_dir = mktempdir()
                push!(MGF_SESSION_DIRS, _session_dir)
            end
            for old_name in filter(f -> occursin("_mgf_", f) && endswith(f, ".txt"), readdir(_session_dir))
                rm(joinpath(_session_dir, old_name); force=true)
            end
            safe_name = "_mgf_$(replace(split(tempname(), '/')[end], '-' => '_')).txt"
            data_path = joinpath(_session_dir, safe_name)
            write(data_path, file_content)
            _load_path = data_path
            println("[GUI] Session upload saved to: $data_path")
        end
    end

    # --- Internal ---
    @private _load_path = ""
    @private _session_dir = ""
    @private _live_results = GaussianFit1D.FitResult[]
    @private _live_n = 0
    @private x_data = Float64[]
    @private y_data = Float64[]
    @private all_results = GaussianFit1D.FitResult[]
    @private cfg_store = GaussianFit1D.FitConfig()
    @private _table_rows = Dict{String,Any}[]

    # --- Handlers ---

    @onchange filepath begin
        fit_ready = false
        if isempty(file_content); _load_path = String(strip(filepath)); end
        error_msg = ""
        plot_json = "{}"
        bic_plot_json = "{}"
        residuals_plot_json = "{}"
        table_json = "[]"
        warnings_html = ""
        progress_percent = 0
        progress_text = ""

        fp = String(strip(filepath))
        if length(fp) > 2
            dirpath = dirname(fp)
            if isdir(dirpath)
                try
                    entries = readdir(dirpath; join=true, sort=true)
                    candidates = filter(f -> endswith(f, ".txt") || endswith(f, ".dat") ||
                                         endswith(f, ".csv") || endswith(f, ".tsv"), entries)
                    if !isempty(candidates)
                        opts = ["<option value=\"$f\">" for f in candidates]
                        datalist_html = "<datalist id=\"mgf-filelist\">" * join(opts) * "</datalist>"
                    else
                        datalist_html = ""
                    end
                catch
                    datalist_html = ""
                end
            else
                datalist_html = ""
            end
        else
            datalist_html = ""
        end
    end

    @onbutton runFitTrigger begin
        fp = isempty(_load_path) ? String(strip(filepath)) : _load_path
        if isempty(fp)
            error_msg = "Please enter a data file path."
            return
        end
        println("=== [GUI] FIT TRIGGERED ===")
        running = true
        progress_percent = 0
        error_msg = ""
        progress_text = "Loading data..."
        runFitTrigger = false
        @async begin
            try
                cfg = GaussianFit1D.build_config(Dict{String,Any}(
                    "filepath" => fp,
                    "output_dir" => isempty(output_dir) ? nothing : output_dir,
                    "min_spacing" => min_spacing,
                    "max_spacing" => max_spacing,
                    "fwhm_min" => fwhm_min,
                    "fwhm_max" => fwhm_max,
                    "offset_to_zero" => offset_to_zero,
                    "amplitude_min_fraction" => amplitude_min_fraction,
                    "asymmetric_edges" => asymmetric_edges,
                    "edge_sigma_min" => edge_sigma_min,
                    "edge_sigma_max" => edge_sigma_max,
                    "global_maxiter" => global_maxiter,
                    "global_maxtime" => global_maxtime,
                    "global_tol" => global_tol,
                    "nlopt_algorithm" => nlopt_algorithm,
                    "curve_fit_maxfev" => curve_fit_maxfev,
                    "early_stop_patience" => early_stop_patience,
                    "early_stop_dbic" => early_stop_dbic,
                    "bic_competition_threshold" => bic_threshold,
                    "use_student_bic" => false,
                    "no_show" => true,
                ))

                x_raw, y_raw = GaussianFit1D.load_data(fp)
                cfg_store = cfg

                x_data = copy(x_raw)
                y_data = copy(y_raw)

                y_min = minimum(y_raw)
                if cfg.offset_to_zero && abs(y_min) > 1e-10
                    y_offset = Printf.@sprintf("%.4f", abs(y_min))
                    y_data = y_data .- y_min
                else
                    y_offset = "0.0000"
                end

                progress_text = "Fitting models..."
                progress_percent = 5
                log_text = ""
                table_json = "[]"
                model_labels_json = "[]"

                _live_results = GaussianFit1D.FitResult[]
                _live_n = 0
                progress_cb = (direction, n, current, total, result) -> begin
                    push!(_live_results, result)
                    _live_n = n
                    pct = 5 + round(Int, 90 * current / max(total, 1))
                    progress_percent = pct
                    msg = "Fitting n=$n ($direction) [$current/$total]"
                    progress_text = msg
                    log_text = log_text * (isempty(log_text) ? "" : "\n") * msg
                    if !fit_ready && !isempty(_live_results)
                        fit_ready = true
                    end
                    if fit_ready && !isempty(_live_results)
                        live_snapshot = copy(_live_results)
                        fr = GaussianFit1D.FitRunResult(x_data, y_data, live_snapshot, cfg)
                        best_r = GaussianFit1D.best_result(fr)
                        if best_r !== nothing
                            plot_json = build_fit_plot_json(x_data, y_data, best_r, cfg)
                            residuals_plot_json = build_residuals_plot_json(x_data, y_data, best_r, cfg)
                            bic_plot_json = build_bic_plot_json(live_snapshot, best_r, cfg)
                            rows = build_table_rows(live_snapshot, cfg_store)
                            _table_rows = rows
                            table_json = JSON3.write(rows)
                            labels = [r === best_r ? "n=$(r.n_peaks) \u2605" : "n=$(r.n_peaks)" for r in live_snapshot]
                            model_labels_json = JSON3.write(labels)
                        end
                    end
                end

                fit_result = GaussianFit1D.run_fit(x_raw, y_raw, cfg; save_cache=false, progress_callback=progress_cb)
                all_r = fit_result.all_results
                all_results = all_r

                if isempty(all_r)
                    error_msg = "No models could be fitted. Check your data and parameter settings."
                    running = false
                    return
                end

                progress_text = "Rendering results..."

                b = GaussianFit1D.best_result(fit_result)
                bi = findfirst(r -> r === b, all_r)
                best_idx = isnothing(bi) ? 1 : bi
                println("  [GUI] best_result: n=$(b.n_peaks), idx=$best_idx")

                rows = build_table_rows(all_r, cfg_store)
                _table_rows = rows
                table_json = JSON3.write(rows)
                println("  [GUI] table_json done ($(length(rows)) rows)")

                labels = [r === b ? "n=$(r.n_peaks) \u2605" : "n=$(r.n_peaks)" for r in all_r]
                model_labels_json = JSON3.write(labels)
                model_idx = best_idx - 1

                warns = b.warnings
                if isempty(warns)
                    warnings_html = "<em class='text-muted'>No warnings.</em>"
                else
                    items = ["<li>$w</li>" for w in warns]
                    warnings_html = "<ul class='mb-0'>" * join(items) * "</ul>"
                end

                println("  [GUI] building fit + residuals + BIC plots...")
                fit_ready = true
                plot_json = build_fit_plot_json(x_data, y_data, b, cfg)
                residuals_plot_json = build_residuals_plot_json(x_data, y_data, b, cfg)
                bic_plot_json = build_bic_plot_json(all_r, b, cfg)
                println("  [GUI] fit=$(round(length(plot_json)/1024, digits=1)) KB, res=$(round(length(residuals_plot_json)/1024, digits=1)) KB, bic=$(round(length(bic_plot_json)/1024, digits=1)) KB")
                progress_percent = 100
                progress_text = "Done!"
                log_summary = "\n\n=== FIT COMPLETE ===\nBest: n=$(b.n_peaks), BIC=$(round(b.bic,digits=1)), R²=$(round(b.r_squared,digits=5))\nTotal models fitted: $(length(all_r))"
                log_text = log_text * log_summary
                running = false

                try GaussianFit1D.export_results(fit_result.x, fit_result.y, all_r, cfg) catch e @warn "export_results failed" exception=e end

            catch e
                @warn "Fit failed" exception=(e, catch_backtrace())
                error_msg = if e isa ArgumentError || e isa SystemError
                    sprint(showerror, e)
                else
                    try string(e.msg) catch; sprint(showerror, e) end
                end
                running = false
                fit_ready = false
            end
        end
    end

    @onchange refreshTrigger begin
        idx = model_idx + 1
        if idx >= 1 && idx <= length(all_results)
            result = all_results[idx]
            warns = result.warnings
            if isempty(warns)
                warnings_html = "<em class='text-muted'>No warnings.</em>"
            else
                items = ["<li>$w</li>" for w in warns]
                warnings_html = "<ul class='mb-0'>" * join(items) * "</ul>"
            end
            plot_json = build_fit_plot_json(x_data, y_data, result, cfg_store)
            residuals_plot_json = build_residuals_plot_json(x_data, y_data, result, cfg_store)
            bic_plot_json = build_bic_plot_json(all_results, result, cfg_store)
        end
    end
end

const MGF_SESSION_DIRS = Set{String}()
atexit(() -> for d in MGF_SESSION_DIRS; rm(d; recursive=true, force=true); end)

# ===========================================================================
# UI Helpers
# ===========================================================================

function param_div(label, fieldname; kind="number", step=nothing, min_val=nothing)
    if kind == "checkbox"
        return H.div(class="mgf-param", [
            H.label(label, class="mgf-param-label"),
            H.label(class="mgf-checkbox", [
                """<input type="checkbox" v-model="$fieldname"/>""",
            ]),
        ])
    end
    step_attr = step !== nothing ? "step=\"$step\"" : ""
    min_attr = min_val !== nothing ? "min=\"$min_val\"" : ""
    H.div(class="mgf-param", [
        H.label(label, class="mgf-param-label"),
        """<input type="$kind" v-model="$fieldname" $step_attr $min_attr class="mgf-input"/>""",
    ])
end

# ===========================================================================
# UI Definition
# ===========================================================================

function ui()
    [
        H.h2("Multi-Gaussian Fit", style="margin-top:8px;margin-bottom:2px;font-size:1.4rem"),

        H.div(style="display:flex;gap:0;margin-bottom:4px;align-items:flex-end", [
            H.div(style="flex:5;position:relative", [
                H.small("Data file path", class="text-muted", style="display:block;font-size:.75rem;margin-bottom:2px"),
                """<div class="mgf-dropzone"
                  @dragover.prevent="\$event.target.classList.add('drag-over')"
                  @dragleave.prevent="\$event.target.classList.remove('drag-over')"
                  @drop.prevent="mgf_handleDrop(\$event); \$event.target.classList.remove('drag-over')">
                  <div style="display:flex;gap:0">
                  <input id="mgf-filepath-input" v-model="filepath" @input="__fpRaw = \$event.target.value" list="mgf-filelist" placeholder="/path/to/data.txt" style="flex:100;padding:2px 6px;font-size:.82rem;height:28px;border:1px solid #ced4da;border-radius:4px 0 0 4px;box-sizing:border-box"/>
                  <button type="button" id="mgf-browse-btn" style="flex:0 0 32px;padding:0;height:28px;font-size:.82rem;background:#f0f2f5;border:1px solid #ced4da;border-left:0;border-radius:0 4px 4px 0;cursor:pointer" title="Browse for text files">&#128194;</button>
                  </div>
                  <span v-if="false" class="mgf-drop-hint">Drop file here</span>
                </div>""",
                """<div v-html="datalist_html"></div>""",
            ]),
            H.div(style="flex:3;padding-left:8px", [
                H.small("Output dir (optional)", class="text-muted", style="display:block;font-size:.75rem;margin-bottom:2px"),
                """<input v-model="output_dir" placeholder="auto (~/.multigaussianfit)" style="width:100%;padding:2px 6px;font-size:.82rem;height:28px;border:1px solid #ced4da;border-radius:4px;box-sizing:border-box"/>""",
            ]),
            H.div(style="flex:0 0 auto;padding-left:8px", [
                H.small("\u00a0", class="text-muted", style="display:block;font-size:.75rem;margin-bottom:2px"),
                """<button type="button" class="mgf-btn-run" style="width:100%;min-width:90px" :disabled="running" @click="runFitTrigger = true"><span class="mgf-btn-icon">&#9654;</span> Run fit</button>""",
            ]),
        ]),

        H.div(v__if="running", style="margin:6px 0", [
            H.div(class="mgf-progress-bar", [
                """<div class="mgf-progress-fill" :style="{width: progress_percent+'%'}"></div>""",
            ]),
            H.div(style="text-align:center;font-size:.75rem;color:#6c757d;margin-top:4px", [
                "{{progress_text}}",
            ]),
        ]),

        H.div(v__if="error_msg != ''", class="mgf-error", [
            H.span("{{error_msg}}"),
        ]),

        H.div(class="mgf-params", [
            H.span("Boundary limits", class="mgf-group-label"),
            param_div("Min space (nm)", "min_spacing"; step=0.01),
            param_div("Max space (nm)", "max_spacing"; step=0.01),
            param_div("FWHM min (nm)", "fwhm_min"; step=0.01, min_val=0.01),
            param_div("FWHM max (nm)", "fwhm_max"; step=0.01, min_val=0.01),
             H.span(class="mgf-param-cluster", [
                param_div("Offset to 0", "offset_to_zero"; kind="checkbox"),
                H.div(class="mgf-param-offset", [
                    H.small("Offset", class="mgf-param-label"),
                    """<input v-model="y_offset" disabled class="mgf-input" style="width:60px;text-align:center;background:#f0f2f5"/>""",
                ]),
            ]),
            param_div("Amp. min factor", "amplitude_min_fraction"; step=0.05, min_val=0.0),
            H.span(class="mgf-param-cluster", [
                param_div("Asymm.", "asymmetric_edges"; kind="checkbox"),
                param_div("Edge \u03c3 min", "edge_sigma_min"; step=0.1, min_val=0.1),
                param_div("Edge \u03c3 max", "edge_sigma_max"; step=0.1, min_val=0.1),
            ]),
            H.span("", class="mgf-sep"),
            H.span("Convergence", class="mgf-group-label"),
            param_div("NLopt max iter", "global_maxiter"; step=1000, min_val=100),
            param_div("NLopt max time", "global_maxtime"; step=5, min_val=1),
            param_div("NLopt tol", "global_tol"; step=0.0001, min_val=1e-8),
            param_div("LsqFit evals", "curve_fit_maxfev"; step=100, min_val=10),
            param_div("Patience", "early_stop_patience"; step=1, min_val=1),
            param_div("BIC stop", "early_stop_dbic"; step=10, min_val=0),
            param_div("BIC threshold", "bic_threshold"; step=1.0, min_val=0.0),
        ]),

        H.div(v__if="fit_ready", style="margin-top:8px", [
            H.div(class="mgf-tab-bar", [
                """<button class="mgf-tab-btn" :class="{'mgf-tab-active':active_tab=='plots'}" @click="active_tab='plots'">Plots</button>
                <button class="mgf-tab-btn" :class="{'mgf-tab-active':active_tab=='comparison'}" @click="active_tab='comparison'">Comparison</button>
                <button class="mgf-tab-btn" :class="{'mgf-tab-active':active_tab=='warnings'}" @click="active_tab='warnings'">Warnings</button>
                <button class="mgf-tab-btn" :class="{'mgf-tab-active':active_tab=='log'}" @click="active_tab='log'">Log</button>""",
            ]),

            H.div(v__show="active_tab=='plots'", [
                H.div(class="mgf-selector", [
                    H.div(style="text-align:center;font-weight:600;margin-bottom:6px;font-size:.82rem", "Select model"),
                    H.div(style="display:flex;justify-content:center;flex-wrap:wrap;gap:4px 8px", [
                        H.label(v__for="(label,idx) in model_labels_parsed", class="mgf-radio-label", [
                            """<input type="radio" v-model="model_idx" :value="idx" @change="refreshTrigger += 1" />""",
                            H.span("{{ label }}", class="mgf-radio-text"),
                        ]),
                    ]),
                ]),
                H.div(id="plotly-fit-container", style="width:100%;height:500px;border:1px solid #dee2e6;border-radius:4px;overflow:hidden;margin-bottom:8px"),
                H.div(id="plotly-residuals-container", style="width:100%;height:200px;border:1px solid #dee2e6;border-radius:4px;overflow:hidden;margin-bottom:8px"),
                H.div(id="plotly-bic-container", style="width:100%;height:250px;border:1px solid #dee2e6;border-radius:4px;overflow:hidden"),
            ]),

            H.div(v__show="active_tab=='comparison'", [
                H.div(class="mgf-card", [
                    H.h5("Model comparison"),
                    H.table(class="mgf-table", [
                        H.thead(H.tr([
                            H.th("n"), H.th("BIC"), H.th("dBIC"),
                            H.th("R\u00b2"), H.th("RSS"), H.th("params"), H.th(""),
                        ])),
                        H.tbody([
                            H.tr(v__for="row in table_rows_parsed", [
                                H.td("{{ row.n_peaks }}"),
                                H.td("{{ row.BIC }}"),
                                H.td("{{ row.dBIC }}"),
                                H.td("{{ row.R2 }}"),
                                H.td("{{ row.RSS }}"),
                                H.td("{{ row.n_params }}"),
                                H.td("{{ row.best ? '\\u2605' : (row.competitive ? '\\u25cf' : '') }}"),
                            ]),
                        ]),
                    ]),
                ]),
            ]),

            H.div(v__show="active_tab=='warnings'", [
                H.div(style="max-height:500px;overflow:auto;font-size:.82rem;padding:8px;border:1px solid #dee2e6;border-radius:4px;background:#fff",
                    v__html="warnings_html"),
            ]),

            H.div(v__show="active_tab=='log'", [
                H.div(style="max-height:500px;overflow:auto;font-family:monospace;font-size:.78rem;padding:10px;border:1px solid #dee2e6;border-radius:4px;background:#1e1e1e;color:#d4d4d4;white-space:pre-wrap",
                    "{{log_text}}"),
            ]),
        ]),

        H.div(v__if="!fit_ready", class="mgf-placeholder", [
            H.p("No fit results yet. Enter a data file path and click Run fit.", class="text-muted", style="padding-top:250px"),
        ]),
    ]
end

# ===========================================================================
# Vue hooks
# ===========================================================================

Stipple.js_created(::Type{<:FitModel}) = raw"""
this.model_labels_parsed = [];
this.table_rows_parsed = [];
this.__fpRaw = '';
"""

Stipple.js_mounted(::Type{<:FitModel}) = raw"""
var mgf_fitTimer = null;
var mgf_resTimer = null;
var mgf_bicTimer = null;

var mgfBrowseBtn = document.getElementById('mgf-browse-btn');
if (mgfBrowseBtn) {
    var __vue = this;
    mgfBrowseBtn.addEventListener('click', function() {
        var fileInput = document.createElement('input');
        fileInput.type = 'file';
        fileInput.accept = '.txt,.dat,.csv,.tsv';
        fileInput.onchange = function(e) {
            var file = e.target.files[0];
            if (file) {
                var reader = new FileReader();
                reader.onload = function(ev) {
                    __vue.filepath = file.name;
                    __vue.file_content = ev.target.result;
                };
                reader.readAsText(file);
            }
        };
        fileInput.click();
    });
}
this.mgf_handleDrop = function(e) {
    var file = e.dataTransfer.files[0];
    if (file) {
        var self = this;
        var reader = new FileReader();
        reader.onload = function(ev) {
            self.filepath = file.name;
            self.file_content = ev.target.result;
        };
        reader.readAsText(file);
    }
};

this.$watch('plot_json', function(newVal) {
    clearTimeout(mgf_fitTimer);
    var __vue = this;
    mgf_fitTimer = setTimeout(function() {
        if (newVal && newVal !== '{}') {
            try {
                var fig = JSON.parse(newVal);
                fig.config = fig.config || {};
                fig.config.responsive = false;
                var el = document.getElementById('plotly-fit-container');
                if (el) {
                    Plotly.react(el, fig.data, fig.layout, fig.config);
                } else console.warn('plotly-fit-container not found');
            } catch(e) { console.error('Plotly fit error:', e); }
        }
    }, 100);
}, { immediate: true });

this.$watch('residuals_plot_json', function(newVal) {
    clearTimeout(mgf_resTimer);
    mgf_resTimer = setTimeout(function() {
        if (newVal && newVal !== '{}') {
            try {
                var fig = JSON.parse(newVal);
                fig.config = fig.config || {};
                fig.config.responsive = false;
                var el = document.getElementById('plotly-residuals-container');
                if (el) Plotly.react(el, fig.data, fig.layout, fig.config);
                else console.warn('plotly-residuals-container not found');
            } catch(e) { console.error('Plotly residuals error:', e); }
        }
    }, 100);
}, { immediate: true });

this.$watch('bic_plot_json', function(newVal) {
    clearTimeout(mgf_bicTimer);
    mgf_bicTimer = setTimeout(function() {
        if (newVal && newVal !== '{}') {
            try {
                var fig = JSON.parse(newVal);
                fig.config = fig.config || {};
                fig.config.responsive = false;
                var el = document.getElementById('plotly-bic-container');
                if (el) Plotly.react(el, fig.data, fig.layout, fig.config);
                else console.warn('plotly-bic-container not found');
            } catch(e) { console.error('Plotly BIC error:', e); }
        }
    }, 100);
}, { immediate: true });

this.$watch('model_labels_json', function(newVal) {
    try { this.model_labels_parsed = JSON.parse(newVal || '[]'); }
    catch(e) { this.model_labels_parsed = []; }
}, { immediate: true });

this.$watch('table_json', function(newVal) {
    try { this.table_rows_parsed = JSON.parse(newVal || '[]'); }
    catch(e) { this.table_rows_parsed = []; }
}, { immediate: true });
"""

# ===========================================================================
# Layout and route
# ===========================================================================

const MGF_LAYOUT = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <% Stipple.sesstoken() %>
    <title>Multi-Gaussian Fit</title>
    <style>[v-cloak] { display: none }</style>
    <style>
      @keyframes mgf-progress-anim { 0% { background-position: 0 0; } 100% { background-position: 40px 0; } }
      .mgf-progress-bar { width:100%; height:8px; background:#e9ecef; border-radius:4px; overflow:hidden; }
      .mgf-progress-fill { height:100%; border-radius:4px; background:linear-gradient(90deg,#0d6efd,#6610f2); background-size:40px 40px; animation:mgf-progress-anim 0.8s linear infinite; transition:width .3s ease; }
      .mgf-btn-run { display:inline-flex; align-items:center; gap:6px; padding:4px 16px; font-size:.85rem; font-weight:600; color:#fff; background:#0d6efd; border:none; border-radius:6px; cursor:pointer; height:28px; box-shadow:0 2px 4px rgba(13,110,253,.3); transition:all .15s ease; }
      .mgf-btn-run:hover:not(:disabled) { background:#0b5ed7; box-shadow:0 4px 8px rgba(13,110,253,.4); transform:translateY(-1px); }
      .mgf-btn-run:active:not(:disabled) { transform:translateY(0); box-shadow:0 1px 2px rgba(13,110,253,.3); }
      .mgf-btn-run:disabled { opacity:.65; cursor:not-allowed; filter:grayscale(30%); }
      .mgf-btn-icon { font-size:.95rem; line-height:1; }
      .mgf-dropzone { position:relative; }
      .mgf-dropzone.drag-over::after { content:''; position:absolute; inset:0; border:2px dashed #0d6efd; border-radius:4px; background:rgba(13,110,253,.05); pointer-events:none; z-index:10; }
      .mgf-drop-hint { position:absolute; right:8px; top:50%; transform:translateY(-50%); font-size:.75rem; color:#0d6efd; font-weight:600; pointer-events:none; z-index:11; }
      .mgf-error { padding:6px 10px; font-size:.82rem; color:#842029; background:#f8d7da; border:1px solid #f5c2c7; border-radius:4px; margin-bottom:8px; }
      .mgf-params { display:flex; flex-wrap:nowrap; overflow-x:auto; align-items:flex-end; gap:0 4px; padding:6px 0; margin-bottom:10px; }
      .mgf-param { flex:0 0 auto; width:78px; overflow:hidden; }
      .mgf-param-label { display:block; font-size:.7rem; line-height:1.1; margin-bottom:1px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; color:#6c757d; }
      .mgf-param-offset { flex:0 0 auto; text-align:center; }
      .mgf-param-cluster { display:flex; gap:1px; }
      .mgf-param-cluster .mgf-param { width:auto; max-width:74px; }
      .mgf-param-cluster .mgf-param-label { font-size:.60rem; }
      .mgf-group-label { flex:0 0 auto; align-self:center; font-size:.68rem; line-height:1; font-weight:700; color:#6c757d; margin:0 8px 0 2px; white-space:nowrap; }
      .mgf-sep { flex:0 0 auto; border-left:1px solid #ccc; margin:0 6px; align-self:stretch; }
      .mgf-checkbox { display:flex; align-items:center; height:24px; margin-top:1px; cursor:pointer; }
      .mgf-checkbox input { width:16px; height:16px; cursor:pointer; }
      .mgf-input { padding:0 2px; font-size:.75rem; height:24px; width:100%; border:1px solid #ced4da; border-radius:4px; box-sizing:border-box; }
      .mgf-selector { text-align:center; margin-bottom:8px; padding:8px; background:#f0f2f5; border-radius:4px; }
      .mgf-radio-label { display:inline-flex; align-items:center; margin-right:14px; font-size:0.9rem; cursor:pointer; white-space:nowrap; }
      .mgf-radio-label input[type="radio"] { margin-right:4px; cursor:pointer; accent-color:#0d6efd; }
      .mgf-tab-bar { display:flex; border-bottom:2px solid #c0c7ce; margin-bottom:8px; gap:2px; }
      .mgf-tab-btn { padding:8px 18px; font-size:.85rem; background:transparent; border:1px solid transparent; border-bottom:none; border-radius:6px6px00; cursor:pointer; color:#495057; font-weight:600; transition:all .15s; }
      .mgf-tab-btn:hover { color:#0d6efd; background:#e7f1ff; }
      .mgf-tab-active { color:#0d6efd !important; background:#fff !important; border-color:#c0c7ce !important; margin-bottom:-2px; }
      .mgf-radio-text { user-select:none; }
      .mgf-table { width:100%; border-collapse:collapse; font-size:.82rem; text-align:center; }
      .mgf-table th { background:#f4f6f9; font-weight:600; padding:5px 8px; border-bottom:2px solid #dee2e6; font-size:.76rem; }
      .mgf-table td { padding:4px 8px; border-bottom:1px solid #eee; }
      .mgf-card { flex:1; background:#fff; border:1px solid #dee2e6; border-radius:4px; padding:10px; }
      .mgf-card h5 { font-size:.88rem; font-weight:600; margin:0 0 6px; }
      .mgf-placeholder { text-align:center; padding-top:200px; height:700px; }
      .mgf-placeholder p { color:#adb5bd; }
    </style>
    <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
    <% join(Stipple.Layout.theme(; core_theme = true), "\\n    ") %>
  </head>
  <body>
    <% Stipple.page(model, partial = true, v__cloak = true,
          [Stipple.Genie.Renderer.Html.@yield],
          Stipple.@if(:isready),
          class = "container-fluid",
          style = "max-width:1440px; margin:0 auto; padding:10px 24px 40px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif; color:#212529",
          core_theme = true,
          include_themes = false,
          include_deps = true
        )
    %>
  </body>
</html>
"""

Stipple.Pages.Page("/"; view=ui, model=FitModel, layout=MGF_LAYOUT)

function run(; port::Int=8888, host::String="0.0.0.0")
    println("\n  Multi-Gaussian Fit GUI starting at http://localhost:$port")
    println("  Press Ctrl+C to stop.\n")
    Stipple.up(port, host)
    wait()
end

function run(port::Int)
    run(; port=port)
end

export FitModel, run, build_plot_json

end # module
