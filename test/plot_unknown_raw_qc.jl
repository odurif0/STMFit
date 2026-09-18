#!/usr/bin/env julia
# Saved-geometry visual QC only. No fitting, N selection, classifier or labels.
# julia --startup-file=no --threads=1 --project=. test/plot_unknown_raw_qc.jl --help
module UnknownRawQC

using STMSXMIO
using GaussianFit2D: PatternConfig, preprocess_channel
using Statistics, TOML, Printf, SHA
using Plots

include(joinpath(@__DIR__, "lib", "patch_preprocessing.jl"))
using .PatchPreprocessing: load_patch_preprocessing

const ORIENTATION = "native reader orientation: bwd x-reversed; no extra registration"
const RESIDUAL_REASON = "Omitted: rounded decoded TSV geometry is not the exact original model/cache. " *
    "The fitting image used native fused ROI data and a ROI-dependent zero; separate patch views " *
    "use independently flattened/smoothed channels. No exact residual is claimed."
const DISPLAY_QUANTILES = (0.005, 0.995) # display only; never a scientific threshold

function read_columns(path, columns)
    lines = readlines(path)
    isempty(lines) && error("Empty TSV: $path")
    header = String.(split(chomp(lines[1]), '\t'; keepempty=true))
    length(unique(header)) == length(header) || error("Duplicate TSV column: $path")
    all(c -> c in header, columns) || error("Required columns missing in $path: $columns")
    indexes = [findfirst(==(c), header) for c in columns]
    rows = Dict{String,String}[]
    for (i, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        values = split(chomp(line), '\t'; keepempty=true)
        length(values) == length(header) || error("Malformed TSV row $(i+1): $path")
        push!(rows, Dict(c => String(values[j]) for (c,j) in zip(columns,indexes)))
    end
    isempty(rows) && error("No TSV rows: $path")
    return rows
end

function keyed(rows, column)
    out = Dict{String,Dict{String,String}}()
    for row in rows
        file = basename(row[column])
        endswith(lowercase(file), ".sxm") || error("Not an SXM filename: $file")
        haskey(out, file) && error("Duplicate file: $file")
        out[file] = row
    end
    return out
end

function load_inputs(features, selected, summary, review)
    count_cols = ["filepath", "N_selected", "refined_policy", "support_2D_ell_nm",
                  "support_2D_circ_nm", "artifact_fwd_bwd_corr", "artifact_fwd_bwd_nrmse"]
    selected_rows = read_columns(selected, count_cols)
    counts = keyed(selected_rows, "filepath")
    files = [basename(r["filepath"]) for r in selected_rows] # original membership/order
    # Read ONLY count metadata, not saved chemical assignments or class counts.
    summaries = keyed(read_columns(summary, ["file", "N_selected"]), "file")
    reviews = keyed(read_columns(review, ["file", "N_predicted_lobes", "review_status", "review_reasons",
                                          "mean_confidence", "uncertain_fraction"]), "file")
    Set(keys(summaries)) == Set(files) || error("Application summary coverage differs")
    Set(keys(reviews)) == Set(files) || error("Review queue coverage differs")
    geometry_cols = ["file", "N", "lobe", "x_nm", "y_nm", "t_nm", "u_nm", "axis_x", "axis_y",
                     "origin_x_nm", "origin_y_nm", "source"]
    byfile = Dict{String,Vector{Dict{String,String}}}()
    for row in read_columns(features, geometry_cols)
        push!(get!(byfile, basename(row["file"]), Dict{String,String}[]), row)
    end
    Set(keys(byfile)) == Set(files) || error("Geometry coverage differs from selected summary")
    for file in files
        n = parse(Int, counts[file]["N_selected"])
        n > 0 || error("Nonpositive saved N: $file")
        parse(Int, summaries[file]["N_selected"]) == n || error("Application N differs: $file")
        parse(Int, reviews[file]["N_predicted_lobes"]) == n || error("Review N differs: $file")
        rs = sort!(byfile[file]; by=r -> parse(Int,r["lobe"]))
        [parse(Int,r["lobe"]) for r in rs] == collect(1:n) || error("Incomplete/duplicate lobe keys: $file")
        all(r -> parse(Int,r["N"]) == n, rs) || error("Geometry N differs: $file")
        for r in rs, c in geometry_cols[4:11]
            isfinite(parse(Float64,r[c])) || error("Nonfinite geometry: $file $c")
        end
    end
    return (; files, counts, reviews, byfile)
end

# Deliberately do not use get_channel's missing-direction fallback.
function exact_channel(img, name, direction)
    channels = filter(c -> lowercase(c.name) == lowercase(name) &&
                           lowercase(c.direction) == lowercase(direction), img.channels)
    length(channels) <= 1 || error("Ambiguous channel $name $direction")
    return isempty(channels) ? nothing : only(channels)
end

function prepare_view(img, direction, settings)
    ch = exact_channel(img, "Z", direction)
    ch === nothing && return nothing
    scale, unit = STMSXMIO._value_scale(ch.unit)
    # Raw is the actual native-reader matrix: no flatten, smoothing or NaN fill.
    # Raw stays full resolution even if native preprocessing config uses stride.
    raw_x, raw_y = STMSXMIO._coordinate_vectors(img)
    raw = ch.data .* scale
    any(isfinite, raw) || error("No finite Z $direction pixels")
    cfg = PatternConfig(filepath=img.filepath, channel="Z", direction=direction,
                        stride=settings.stride, flatten=settings.flatten,
                        smooth_radius_px=settings.smooth_radius_px, no_plot=true)
    xs, ys, _native_imputed_raw, flat, smooth, scaled_unit, noise = preprocess_channel(img, ch, cfg)
    unit == scaled_unit || error("Unit conversion mismatch")
    return (; direction, raw_x, raw_y, raw, xs, ys, flat, smooth, unit,
            raw_nonfinite=count(!isfinite,raw), noise)
end

function saved_centers(rows)
    return [(lobe=parse(Int,r["lobe"]), x=parse(Float64,r["x_nm"]),
             y=parse(Float64,r["y_nm"])) for r in rows]
end

function geometry_check(rows, img)
    discrepancies = Float64[]
    for r in rows
        t,u,ax,ay,ox,oy = [parse(Float64,r[c]) for c in
            ("t_nm","u_nm","axis_x","axis_y","origin_x_nm","origin_y_nm")]
        # Same inverse transform as _decode_chain and both patch extractors.
        x,y = ox + t*ax - u*ay, oy + t*ay + u*ax
        push!(discrepancies, hypot(x-parse(Float64,r["x_nm"]), y-parse(Float64,r["y_nm"])))
    end
    centers = saved_centers(rows)
    outside = count(c -> !(0 <= c.x <= img.range_nm[1] && 0 <= c.y <= img.range_nm[2]), centers)
    return (; max_projection_error_nm=maximum(discrepancies), centers_outside=outside)
end

function centers_without_finite_raw_pixel(view, centers)
    view === nothing && return length(centers)
    return count(centers) do c
        inframe = first(view.raw_x) <= c.x <= last(view.raw_x) && first(view.raw_y) <= c.y <= last(view.raw_y)
        inframe || return true
        ix = argmin(abs.(view.raw_x .- c.x)); iy = argmin(abs.(view.raw_y .- c.y))
        return !isfinite(view.raw[iy,ix])
    end
end

function display_limits(z)
    good = filter(isfinite, vec(z))
    isempty(good) && error("No finite display pixels")
    lo,hi = quantile(good, collect(DISPLAY_QUANTILES))
    if lo == hi
        pad = max(abs(lo)*1e-6, 1e-12)
        lo -= pad; hi += pad
    end
    return (lo,hi)
end

function missing_panel(title, reason)
    p = plot(; title, axis=false, grid=false, legend=false, xlims=(0,1), ylims=(0,1))
    annotate!(p, 0.5,0.5, text("UNAVAILABLE\n" * reason, 11, :red, :center))
    return p
end

function image_panel(view, mode, range_nm; centers=nothing)
    direction = view === nothing ? "missing" : view.direction
    view === nothing && return missing_panel("Z channel", "Requested direction not recorded; no fallback")
    xs,ys,z = mode == :raw ? (view.raw_x,view.raw_y,view.raw) : (view.xs,view.ys,view.smooth)
    title = mode == :raw ? "RAW Z $(direction) (native reader)" : "PREPROCESSED Z $(direction)"
    if view.raw_nonfinite > 0
        percent = @sprintf("%.2f",100view.raw_nonfinite/length(view.raw))
        title *= mode == :raw ? "\nMissing raw pixels: $(view.raw_nonfinite) ($percent%); NOT imputed" :
                               "\nWARNING: native imputation of $percent% missing raw pixels"
    end
    p = heatmap(xs,ys,z; title, xlabel="local x (nm)", ylabel="local y (nm)",
                color=:grays, clims=display_limits(z), colorbar_title="Z ($(view.unit))",
                aspect_ratio=:equal, xlims=(0,range_nm[1]), ylims=(0,range_nm[2]),
                widen=false, legend=false, colorbar=:right, grid=false, yflip=false,
                titlefontsize=12, guidefontsize=10, tickfontsize=9,
                margin=5Plots.mm, left_margin=12Plots.mm, right_margin=16Plots.mm)
    if centers !== nothing
        scatter!(p,[c.x for c in centers],[c.y for c in centers]; marker=:cross,
                 markersize=5, markerstrokewidth=1.5, color=:cyan, label=false)
        # Numbers are model component indices, not chemical assignments.
        dx,dy = 0.018range_nm[1],0.018range_nm[2]
        for c in centers
            tx = c.x > 0.92range_nm[1] ? c.x-dx : c.x+dx
            ty = c.y > 0.92range_nm[2] ? c.y-dy : c.y+dy
            annotate!(p, tx, ty, text(string(c.lobe),10,:cyan,:center))
        end
    end
    return p
end

function render_figure(file, views, img, settings, rows, selected_n; overlay=false)
    centers = overlay ? saved_centers(rows) : nothing
    panels = [image_panel(v,mode,img.range_nm; centers) for mode in (:raw,:smooth) for v in views]
    title = overlay ? "$file | saved N=$selected_n (algorithmic, NOT truth) | numbers = lobe indices" :
                      "$file | BLIND VIEW: no saved count, centers or chemical labels"
    prep = "Native preprocessing: $(settings.flatten), box radius=$(settings.smooth_radius_px) px, stride=$(settings.stride)."
    foot = "$ORIENTATION.\nLocal scan frame; no offset/angle transform; equal nm axes. Per-panel 0.5%-99.5% color limits only."
    top = plot(;axis=false,grid=false,legend=false,xlims=(0,1),ylims=(0,1))
    annotate!(top,0.5,0.65,text(title,13,:black,:center))
    annotate!(top,0.5,0.2,text(prep,10,:black,:center))
    bottom = plot(;axis=false,grid=false,legend=false,xlims=(0,1),ylims=(0,1))
    annotate!(bottom,0.5,0.55,text(foot,10,:black,:center))
    # Constant physical aspect in every panel, including the rectangular scan.
    body = plot(panels...;layout=(2,2))
    return plot(top,body,bottom;layout=grid(3,1;heights=[0.07,0.87,0.06]),size=(1800,1600))
end

escape_tsv(x) = replace(string(x), "\\"=>"\\\\", "\t"=>"\\t", "\r"=>"\\r", "\n"=>"\\n")
function write_tsv(path, columns, rows)
    open(path,"w") do io
        println(io,join(columns,'\t'))
        for row in rows
            println(io,join([escape_tsv(get(row,c,"")) for c in columns],'\t'))
        end
    end
end

function header_inventory(img)
    h = img.header
    scan = Dict{String,Any}("file"=>basename(img.filepath),"status"=>"ok",
        "bias_header"=>get(h,"BIAS",""),"bias_V"=>something(tryparse(Float64,get(h,"BIAS","")),NaN),
        "width_px"=>img.width,"height_px"=>img.height,
        "range_x_nm"=>img.range_nm[1],"range_y_nm"=>img.range_nm[2],
        "native_dx_nm"=>img.range_nm[1]/(img.width-1),"native_dy_nm"=>img.range_nm[2]/(img.height-1),
        "offset_x_nm"=>img.offset_nm[1],"offset_y_nm"=>img.offset_nm[2],
        "scan_angle_header"=>get(h,"SCAN_ANGLE",""),"scan_direction"=>get(h,"SCAN_DIR",""),
        "record_date"=>get(h,"REC_DATE",""),"record_time"=>get(h,"REC_TIME",""),
        "record_temp_header_unverified"=>get(h,"REC_TEMP",""),
        "scan_time_header"=>get(h,"SCAN_TIME",""),"acquisition_time_header"=>get(h,"ACQ_TIME",""),
        "z_controller_header"=>get(h,"Z-CONTROLLER",""),
        "lockin_status"=>get(h,"Lock-in>Lock-in status",""),
        "lockin_modulated_signal"=>get(h,"Lock-in>Modulated signal",""),
        "lockin_demodulated_signal"=>get(h,"Lock-in>Demodulated signal",""),
        "lockin_amplitude_header"=>get(h,"Lock-in>Amplitude",""),
        "lockin_frequency_Hz_header"=>get(h,"Lock-in>Frequency (Hz)",""),
        "scan_file_header"=>get(h,"SCAN_FILE",""),
        "data_info_header"=>get(h,"DATA_INFO",""), "reader_orientation"=>ORIENTATION)
    channels = Dict{String,Any}[]
    for c in img.channels
        vals = filter(isfinite, vec(c.data))
        push!(channels,Dict("file"=>basename(img.filepath),"channel_name"=>c.name,
            "recorded_unit"=>c.unit,"recorded_direction_expanded"=>c.direction,
            "reader_transform"=>lowercase(c.direction)=="bwd" ? "x-reversed from acquisition bytes" : "no direction reversal",
            "pixels"=>length(c.data),"finite_pixels"=>length(vals),"nonfinite_pixels"=>count(!isfinite,c.data),
            "min_recorded_unit"=>isempty(vals) ? NaN : minimum(vals),
            "median_recorded_unit"=>isempty(vals) ? NaN : median(vals),
            "max_recorded_unit"=>isempty(vals) ? NaN : maximum(vals),
            "sample_sd_recorded_unit"=>length(vals)<2 ? NaN : std(vals;corrected=true)))
    end
    headers = [Dict("file"=>basename(img.filepath),"key"=>k,"value"=>h[k]) for k in sort(collect(keys(h)))]
    return (;scan,channels,headers)
end

const SCAN_COLUMNS = ["file","status","bias_header","bias_V","width_px","height_px","range_x_nm","range_y_nm",
    "native_dx_nm","native_dy_nm","offset_x_nm","offset_y_nm","scan_angle_header","scan_direction",
    "record_date","record_time","record_temp_header_unverified","scan_time_header","acquisition_time_header",
    "z_controller_header","lockin_status","lockin_modulated_signal","lockin_demodulated_signal",
    "lockin_amplitude_header","lockin_frequency_Hz_header","scan_file_header","data_info_header","reader_orientation"]
const CHANNEL_COLUMNS = ["file","channel_name","recorded_unit","recorded_direction_expanded","reader_transform",
                         "pixels","finite_pixels","nonfinite_pixels","min_recorded_unit","median_recorded_unit",
                         "max_recorded_unit","sample_sd_recorded_unit"]
const QC_COLUMNS = ["file","status","error","saved_N_algorithmic","support_policy_saved","support_ell_nm_saved",
    "support_circ_nm_saved","geometry_source","max_projection_error_nm","centers_outside_frame",
    "raw_fwd_nonfinite","raw_bwd_nonfinite","raw_fwd_finite_fraction","raw_bwd_finite_fraction",
    "centers_without_finite_nearest_raw_pixel_fwd","centers_without_finite_nearest_raw_pixel_bwd","data_warning",
    "saved_artifact_fwd_bwd_corr","saved_artifact_fwd_bwd_nrmse",
    "review_status_saved","review_reasons_saved","uncalibrated_mean_margin_saved","uncertain_fraction_saved",
    "blind_png","overlay_png","residual_status","registration"]

function parse_cli(args)
    keys = ["features","selected-summary","summary","review-queue","data-dir","config","out-dir"]
    opt = Dict{String,String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--help","-h")
            println("""
            Saved-geometry raw STM QC. No fitting or classification.
            Use Julia 1.13 --startup-file=no --threads=1 --project=. and GKSwstype=100.
            Required: --features TSV --selected-summary TSV --summary TSV --review-queue TSV
                      --data-dir DIR --config TOML --out-dir NEW_DIR
            Optional: --file NAME.sxm (one selected scan for preview; default ALL selected files)
            Outputs: separate blind and indexed-overlay PNGs, native SXM header/channel/scan
                     inventory, display/QC metadata, source hashes, README and browser index.
            A failed scan remains in the inventory and gets an explicit error figure; exit != 0.
            """)
            return nothing
        end
        startswith(arg,"--") || error("Unknown argument $arg")
        key = arg[3:end]
        key in vcat(keys,["file"]) || error("Unknown option $arg")
        haskey(opt,key) && error("Duplicate $arg")
        i < length(args) && !startswith(args[i+1],"--") || error("Missing value for $arg")
        opt[key] = args[i+1]; i += 2
    end
    all(k -> haskey(opt,k),keys) || error("Missing required options; use --help")
    for key in ["features","selected-summary","summary","review-queue","config"]
        isfile(opt[key]) || error("Missing input $(opt[key])")
    end
    isdir(opt["data-dir"]) || error("Missing raw data directory")
    out = abspath(opt["out-dir"])
    (ispath(out) || islink(out)) && error("Output already exists: $out")
    return opt
end

filehash(path) = bytes2hex(open(sha256,path))
html_escape(x) = replace(string(x),"&"=>"&amp;","<"=>"&lt;",">"=>"&gt;","\""=>"&quot;")
function write_notes(out, scan_rows, channel_rows, qc_rows, settings)
    successful = filter(r -> get(r,"status","")=="ok",scan_rows)
    biases = sort(unique(string(r["bias_V"]) for r in successful))
    channel_types = sort(unique("$(r["channel_name"]) [$(r["recorded_unit"])] $(r["recorded_direction_expanded"])" for r in channel_rows))
    lockin = sort(unique(string(get(r,"lockin_status","")) for r in successful))
    open(joinpath(out,"README.md"),"w") do io
        println(io,"""
        # Raw STM visual QC (saved geometry only)

        Open `index.html`. Inspect `*_blind.png` before `*_overlay.png` to avoid count/label anchoring.
        All $(length(qc_rows)) requested files remain in `qc_index.tsv`, including any failures.
        No chemical assignments or class colors are used. Saved N is algorithmic, not truth.
        Current QC review flags and uncalibrated mean vote margins are copied as metadata only.
        No fitting, count selection, support choice, classifier, thresholds or parameters changed.

        ## Coordinates and preprocessing

        - Data read by STMSXMIO; this is **native reader orientation**, not unaltered file bytes.
          Backward matrices have x reversed. No additional shift, rotation, reflection or drift correction.
          The native reader does not reverse y for SCAN_DIR. Pixel rows map to increasing native y.
        - Coordinates match the saved geometry and patch extractors: matrix indexing is [y,x],
          x spans 0..SCAN_RANGE_x and y spans 0..SCAN_RANGE_y, with endpoints included.
          Native spacing is range/(pixels-1), not range/pixels. This is the software's local scan
          frame, **not a proven absolute laboratory frame**. SCAN_OFFSET/SCAN_ANGLE are inventoried
          but not applied. Equal nm axes; full field retained. Non-square scans are not stretched.
        - Raw panels use full native matrices with only unit conversion (m to nm), no flattening,
          smoothing or NaN imputation. They are not directly displayed acquisition bytes.
        - Preprocessed panels call GaussianFit2D.preprocess_channel, which uses shared STMSXMIO
          helpers: flatten=$(settings.flatten), smooth_radius_px=$(settings.smooth_radius_px),
          stride=$(settings.stride). Any NaN/Inf replacement is native median imputation, recorded
          by direction in qc_index.tsv; it is not applied to the raw panels. Affected figures display
          explicit missing-pixel/imputation warnings. Render status `ok` means the files rendered,
          not that the scan is complete or scientifically valid. `data_warning` is separate.
        - Per-panel color limits are fixed display quantiles 0.5% and 99.5% of finite pixels;
          saturation changes colors only, not data. Exact color limits/counts are in display_limits.tsv.
          Fwd/bwd color limits are independent: compare colorbars before comparing heights.
        - Centers are the saved x_nm/y_nm values, not newly detected peaks or registration landmarks.
          Projection roundoff and out-of-frame centers are recorded. The nearest original raw pixel
          is also checked for finiteness at each center; a finite nearest pixel does NOT guarantee
          complete surrounding patches, trustworthy lobe shape, or valid coverage. No error bars are available;
          decimal precision is not measurement accuracy. The overlay does not establish one chemical
          unit per lobe, a correct N, or chemical identity. A finite QC margin is not calibrated confidence.
        - Saved support policy/lengths are metadata, not a newly inferred ROI. No support is reselected.
          A drift/tip change between directions can cause misregistration even after the native x flip.
          Fwd/bwd images of the same scan are repeated views, not independent molecules.

        ## Residual

        $RESIDUAL_REASON
        In particular, do not subtract a fused-fit baseline from the separately preprocessed raw panels
        and call it the original fit residual. Only centers are overlaid; no model surface is reconstructed.

        ## Native SXM inventory and measurement limits

        Successfully parsed scans: $(length(successful))/$(length(qc_rows)).
        Recorded BIAS values (V): $(join(biases,", ")).
        Recorded channel/direction types: $(join(channel_types,"; ")).
        Lock-in status header values: $(join(lockin,", ")).

        `scans.tsv` records bias, dimensions/ranges, offsets, angle, scan direction, date/time,
        line/acquisition timing, Z feedback/setpoint and lock-in headers. `channels.tsv` lists only
        actually stored DATA_INFO channels, expanded by direction, with finite/nonfinite pixel counts
        and finite-pixel min/median/max/sample standard deviation (n-1 denominator), in recorded units.
        These are descriptive raw-channel checks, not independent-pixel uncertainties or chemical
        features. No imputation is used for these statistics. `headers.tsv` preserves every
        native parsed header field, including calibration/offset rows in DATA_INFO. Multiline fields
        are TSV-escaped. A configured lock-in field is NOT evidence of a recorded lock-in channel.
        REC_TEMP is unverified header metadata, not a validated sample-temperature measurement.

        Recorded Z and Current, when present, can support cross-view artifact/feedback diagnostics.
        Current under active Z feedback is coupled to topography/setpoint, not independent chemistry.
        This inventory alone does not establish multi-bias images of the same molecule, a bias sweep,
        I(V), dI/dV, or spatial spectroscopy. Filenames, offsets or similar images do not establish
        molecular identity or independence. No other raw directory or spectroscopy sidecar was searched.

        Experimental suggestions only (NOT measured here): independently registered same-molecule
        multi-bias imaging or controlled I(V)/dI/dV with appropriate references could test electronic
        contrast. This report recommends no operating-voltage change or chemical interpretation.
        Tip state, drift, feedback coupling, stability and damage must be assessed experimentally;
        repeated images of one molecule must stay grouped, not be counted as independent validation.
        """)
    end
    open(joinpath(out,"index.html"),"w") do io
        println(io,"<!doctype html><meta charset='utf-8'><title>Raw STM QC</title>")
        println(io,"<style>body{font:17px sans-serif;max-width:1100px;margin:2em auto}td,th{padding:8px;text-align:left}img{max-width:100%}</style>")
        println(io,"<h1>Saved-geometry raw STM QC</h1><p>Open the blind view first. Saved N is algorithmic, not truth. No chemical labels.</p><p><a href='README.md'>Conventions and limits</a> | <a href='scans.tsv'>Scan inventory</a> | <a href='channels.tsv'>Recorded channels</a></p>")
        println(io,"<table><tr><th>File</th><th>Blind view</th><th>Numbered centers</th><th>Rendering status</th><th>Raw-data warning</th></tr>")
        for r in qc_rows
            println(io,"<tr><td>$(html_escape(r["file"]))</td><td><a href='$(html_escape(r["blind_png"]))'>No count/model</a></td><td><a href='$(html_escape(r["overlay_png"]))'>Saved geometry</a></td><td>$(html_escape(r["status"]))</td><td>$(html_escape(get(r,"data_warning","")))</td></tr>")
        end
        println(io,"</table>")
    end
end

function main(args=ARGS)
    VERSION.major == 1 && VERSION.minor == 13 || error("Use Julia 1.13")
    opt = parse_cli(args)
    opt === nothing && return
    input = load_inputs(opt["features"],opt["selected-summary"],opt["summary"],opt["review-queue"])
    settings = load_patch_preprocessing(opt["config"])
    files = input.files
    if haskey(opt,"file")
        opt["file"] in files || error("Preview file is not in original selected summary")
        files = [opt["file"]]
    end
    # Snapshot input hashes before reads; verify after rendering. Missing raw stays explicit.
    source_paths = [opt[k] for k in ("features","selected-summary","summary","review-queue","config")]
    append!(source_paths,[joinpath(opt["data-dir"],f) for f in files])
    hashes = Dict(p => isfile(p) ? filehash(p) : "MISSING" for p in source_paths)
    out = opt["out-dir"]; mkpath(out)
    scan_rows = Dict{String,Any}[]; channel_rows = Dict{String,Any}[]; header_rows = Dict{String,Any}[]
    qc_rows = Dict{String,Any}[]; display_rows = Dict{String,Any}[]
    for (i,file) in enumerate(files)
        rows = input.byfile[file]; countrow = input.counts[file]; review = input.reviews[file]
        stem = splitext(file)[1]
        q = Dict{String,Any}("file"=>file,"status"=>"ok","error"=>"",
            "saved_N_algorithmic"=>countrow["N_selected"],"support_policy_saved"=>countrow["refined_policy"],
            "support_ell_nm_saved"=>countrow["support_2D_ell_nm"],"support_circ_nm_saved"=>countrow["support_2D_circ_nm"],
            "geometry_source"=>join(unique(r["source"] for r in rows),","),
            "saved_artifact_fwd_bwd_corr"=>countrow["artifact_fwd_bwd_corr"],
            "saved_artifact_fwd_bwd_nrmse"=>countrow["artifact_fwd_bwd_nrmse"],
            "review_status_saved"=>review["review_status"],"review_reasons_saved"=>review["review_reasons"],
            "uncalibrated_mean_margin_saved"=>review["mean_confidence"],"uncertain_fraction_saved"=>review["uncertain_fraction"],
            "blind_png"=>stem*"_blind.png","overlay_png"=>stem*"_overlay.png",
            "residual_status"=>RESIDUAL_REASON,"registration"=>ORIENTATION,"data_warning"=>"")
        try
            img = read_sxm(joinpath(opt["data-dir"],file))
            inv = header_inventory(img)
            push!(scan_rows,inv.scan); append!(channel_rows,inv.channels); append!(header_rows,inv.headers)
            check = geometry_check(rows,img)
            q["max_projection_error_nm"] = check.max_projection_error_nm
            q["centers_outside_frame"] = check.centers_outside
            views = [prepare_view(img,d,settings) for d in ("fwd","bwd")]
            for (d,v) in zip(("fwd","bwd"),views)
                if v === nothing
                    q["status"] = "missing_view"; q["error"] *= "Z $d not recorded. "
                    continue
                end
                q["raw_$(d)_nonfinite"] = v.raw_nonfinite
                q["raw_$(d)_finite_fraction"] = 1-v.raw_nonfinite/length(v.raw)
                q["centers_without_finite_nearest_raw_pixel_$(d)"] = centers_without_finite_raw_pixel(v,saved_centers(rows))
                if v.raw_nonfinite > 0
                    q["data_warning"] *= "Z $d has nonfinite raw pixels; native preprocessing imputes them. "
                end
                for (mode,z) in (("raw",v.raw),("preprocessed",v.smooth))
                    lo,hi = display_limits(z)
                    push!(display_rows,Dict("file"=>file,"direction"=>d,"view"=>mode,"unit"=>v.unit,
                        "color_low"=>lo,"color_high"=>hi,"finite_pixels"=>count(isfinite,z),
                        "pixels_below_color_low"=>count(x->isfinite(x)&&x<lo,z),
                        "pixels_above_color_high"=>count(x->isfinite(x)&&x>hi,z)))
                end
            end
            for overlay in (false,true)
                fig = render_figure(file,views,img,settings,rows,countrow["N_selected"];overlay)
                savefig(fig,joinpath(out,q[overlay ? "overlay_png" : "blind_png"]))
                closeall()
            end
        catch err
            q["status"] = "error"; q["error"] = sprint(showerror,err)
            any(r -> r["file"] == file,scan_rows) || push!(scan_rows,Dict("file"=>file,"status"=>"error"))
            for key in ("blind_png","overlay_png")
                fig = missing_panel(file, "Read/render failed. See qc_index.tsv; scan NOT omitted.")
                plot!(fig;size=(1200,800)); savefig(fig,joinpath(out,q[key])); closeall()
            end
            @warn "Raw QC error; retained explicit failure artifacts" file exception=(err,catch_backtrace())
        end
        push!(qc_rows,q)
        println("[$i/$(length(files))] $file: $(q["status"]) (saved geometry only)")
        flush(stdout)
    end
    write_tsv(joinpath(out,"qc_index.tsv"),QC_COLUMNS,qc_rows)
    write_tsv(joinpath(out,"scans.tsv"),SCAN_COLUMNS,scan_rows)
    write_tsv(joinpath(out,"channels.tsv"),CHANNEL_COLUMNS,channel_rows)
    write_tsv(joinpath(out,"headers.tsv"),["file","key","value"],header_rows)
    write_tsv(joinpath(out,"display_limits.tsv"),["file","direction","view","unit","color_low","color_high",
        "finite_pixels","pixels_below_color_low","pixels_above_color_high"],display_rows)
    write_tsv(joinpath(out,"source_hashes.tsv"),["path","sha256"],[Dict("path"=>abspath(p),"sha256"=>hashes[p]) for p in source_paths])
    write_notes(out,scan_rows,channel_rows,qc_rows,settings)
    for p in source_paths
        now = isfile(p) ? filehash(p) : "MISSING"
        now == hashes[p] || error("Input changed during QC: $p")
    end
    all(r -> r["status"] == "ok",qc_rows) || error("QC includes failures/missing views; see $out/qc_index.tsv")
    println("Complete: $(length(files)) files; no fit, count, support or assignment change. Output: $out")
    return (;qc_rows,scan_rows,channel_rows)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    UnknownRawQC.main()
end
