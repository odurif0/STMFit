module ResidualFeatures
# Lightweight subtraction-only interface. Main features and sampling stay intact.
export read_residual_models, subtraction_model
const MUTABLE=Set(["amplitude","baseline","tilt_x","tilt_y","amp_rel","gcv"])
key(r)=(basename(r["file"]),parse(Int,r["lobe"]))
number(r,k)=begin
    value=parse(Float64,r[k]); isfinite(value) || error("Nonfinite residual model $k"); value
end

function read_profile(path,base)
    lines=readlines(path); isempty(lines) && error("Empty residual model")
    header=String.(split(first(lines),'\t'))
    length(unique(header))==length(header) || error("Duplicate residual columns")
    Set(header)==union(Set(keys(first(values(base)))),Set(["profile_view"])) || error("Unexpected residual model columns")
    table=Dict{Tuple{String,Int},Dict{String,String}}(); source=""
    for line in lines[2:end]
        fields=split(line,'\t';keepempty=true)
        length(fields)==length(header) || error("Malformed residual model row")
        row=Dict{String,String}(zip(header,fields)); k=key(row)
        haskey(base,k) && !haskey(table,k) || error("Unexpected/duplicate residual model key")
        row["profile_view"] in ("fwd","bwd") || error("Unknown profile source view")
        isempty(source) && (source=row["profile_view"])
        source==row["profile_view"] || error("Mixed source views")
        for c in setdiff(keys(base[k]),MUTABLE)
            row[c]==base[k][c] || error("Residual model changed frozen field $c for $k")
        end
        for c in MUTABLE; number(row,c); end
        number(row,"amplitude")>0 && number(row,"amp_rel")>0 && number(row,"gcv")>=0 || error("Invalid residual model coefficients")
        table[k]=row
    end
    Set(keys(table))==Set(keys(base)) || error("Incomplete residual model coverage")
    byfile=Dict{String,Vector{Dict{String,String}}}()
    for k in sort(collect(keys(table))); push!(get!(byfile,first(k),Dict{String,String}[]),table[k]); end
    for rs in values(byfile)
        for c in ("baseline","tilt_x","tilt_y","gcv")
            length(unique(r[c] for r in rs))==1 || error("Inconsistent chain coefficient $c")
        end
    end
    return byfile,source
end

function read_residual_models(fwd,bwd,rows;acquisition_shifts=nothing,patch_frames=nothing)
    fwd===nothing && bwd===nothing && return nothing
    fwd!==nothing && bwd!==nothing || error("Both residual model tables are required")
    acquisition_shifts===nothing && patch_frames===nothing || error("Residual model experiment excludes shifts and local frames")
    base=Dict{Tuple{String,Int},Dict{String,String}}()
    for row in rows
        any(startswith(c,"model_axis") || c in ("model_orientation","profile_view") for c in keys(row)) && error("Original global-axis main features required")
        k=key(row); haskey(base,k) && error("Duplicate main feature key"); base[k]=row
        number(row,"skew_ratio")==1 || error("Gaussian base required")
        for c in ("amplitude","sigma_parallel_nm","sigma_perp_nm")
            number(row,c)>0 || error("Nonpositive residual geometry")
        end
        for c in ("x_nm","y_nm","axis_x","axis_y","baseline","tilt_x","tilt_y")
            number(row,c)
        end
    end
    isempty(base) && error("Empty main features")
    fm,fs=read_profile(fwd,base); bm,bs=read_profile(bwd,base)
    fs!=bs || error("Need the two distinct source views, not a duplicate model")
    return (fwd=fm,bwd=bm)
end

function subtraction_model(rows,xs,ys,eval_peak)
    ax,ay=number(rows[1],"axis_x"),number(rows[1],"axis_y")
    b0,bx,by=(number(rows[1],c) for c in ("baseline","tilt_x","tilt_y"))
    model=zeros(length(ys),length(xs))
    for iy in eachindex(ys),ix in eachindex(xs)
        model[iy,ix]=b0+bx*xs[ix]+by*ys[iy]
    end
    for r in rows
        amp,cx,cy,sp,sq,sk=(number(r,c) for c in ("amplitude","x_nm","y_nm","sigma_parallel_nm","sigma_perp_nm","skew_ratio"))
        for iy in eachindex(ys),ix in eachindex(xs)
            model[iy,ix]+=eval_peak(xs[ix],ys[iy],cx,cy,ax,ay,amp,sp,sq,sk)
        end
    end
    return model
end
end
