module ImageForeground
using Statistics, TOML
using STMSXMIO

function settings(path)
    c=TOML.parsefile(path)
    fields=Dict(
        "model"=>["transfer_run","controls_run","transfer_settings_sha256","files_sha256","segmentation","background","arms","guard_padding_nm","guard_shape"],
        "selection"=>["policy","connectivity","component_order","noise_multiplier","numeric_height_floor_nm","min_background_pixels","min_background_fraction","min_row_background_pixels","anchor_membership"],
        "preprocessing"=>["channel","stride","initial_flatten","smooth_radius_px","rank_rtol","missing_policy","centering","binary_arrays"])
    Set(keys(c))==Set(keys(fields)) || error("Unexpected sections")
    for (s,ks) in fields; Set(keys(c[s]))==Set(ks) || error("Unknown/missing $s setting"); end
    m,s,p=c["model"],c["selection"],c["preprocessing"]
    m["segmentation"]=="exact_otsu_and_source_hf_floor" &&
        m["background"]=="one_masked_joint_x_slope_and_row_offsets" &&
        m["arms"]==["initial","joint_background"] && m["guard_shape"]=="axis_aligned_box_ceil_pixels" || error("Changed model")
    s["policy"]=="retain_all_cases_and_components_no_count_or_chemical_prior" &&
        s["connectivity"]==8 && s["component_order"]=="first_row_column" &&
        s["anchor_membership"]=="center_pixel_only_no_patch_refit" || error("Changed support policy")
    p["channel"]=="Z" && p["stride"]==1 && p["initial_flatten"]=="plane+rows" &&
        p["missing_policy"]=="observed_only_unavailable_rows_remain_missing" &&
        p["centering"]=="source_finite_median" &&
        p["binary_arrays"]=="little_endian_column_major_float64_or_uint8" || error("Changed observation convention")
    for (d,ks) in ((m,["guard_padding_nm"]),(s,["noise_multiplier","numeric_height_floor_nm"]),(p,["rank_rtol"]))
        all(k->d[k] isa Real && !(d[k] isa Bool) && isfinite(d[k]) && d[k]>0,ks) || error("Invalid positive setting")
    end
    0<s["min_background_fraction"]<1 && 0<p["rank_rtol"]<1 || error("Invalid fraction/tolerance")
    for k in ("min_background_pixels","min_row_background_pixels")
        s[k] isa Int && s[k]>=3 || error("Invalid background count")
    end
    p["smooth_radius_px"] isa Int && p["smooth_radius_px"]>=0 || error("Invalid smoothing radius")
    all(occursin(r"^[0-9a-f]{64}$",m[k]) for k in ("transfer_settings_sha256","files_sha256")) || error("Invalid input hash")
    c
end

"Source horizontal first-difference MAD; descriptive scale, not calibrated noise."
function hf_scale(raw)
    differences=Float64[]
    for y in axes(raw,1),x in 2:size(raw,2)
        isfinite(raw[y,x]) && isfinite(raw[y,x-1]) && push!(differences,raw[y,x]-raw[y,x-1])
    end
    isempty(differences) && return NaN
    1.4826*median(abs.(differences.-median(differences)))/sqrt(2)
end

"Exact finite-observation two-level variance maximum, with no histogram bins."
function otsu(values)
    v=sort(filter(isfinite,vec(values)))
    isempty(v) && return (;threshold=NaN,score=NaN,split=0,n=0)
    first(v)<last(v) || return (;threshold=first(v),score=0.,split=0,n=length(v))
    origin=first(v)+(last(v)-first(v))/2
    w=v.-origin; total=sum(w); prefix=0.; best=-Inf; at=0; n=length(v)
    for j in 1:n-1
        prefix+=w[j]
        v[j]<v[j+1] || continue
        gap=prefix/j-(total-prefix)/(n-j)
        score=(j/n)*((n-j)/n)*gap^2
        if score>best; best=score; at=j; end
    end
    isfinite(best) && at>0 || error("Otsu numerical failure")
    (;threshold=v[at]+(v[at+1]-v[at])/2,score=best,split=at,n)
end

function segment(z,noise,c)
    observed=isfinite.(z); values=z[observed]
    isempty(values) && return (;mask=falses(size(z)),threshold=NaN,otsu_threshold=NaN,otsu_score=NaN,otsu_split=0,
        median_nm=NaN,floor_threshold_nm=NaN,status="unavailable")
    isfinite(noise) && noise>=0 || error("No source horizontal-difference scale")
    s=c["selection"]; o=otsu(values); center=median(values)
    floor_threshold=center+max(s["noise_multiplier"]*noise,s["numeric_height_floor_nm"])
    threshold=max(o.threshold,floor_threshold)
    mask=BitMatrix(observed .& (z .> threshold))
    (;mask,threshold,otsu_threshold=o.threshold,otsu_score=o.score,otsu_split=o.split,
        median_nm=center,floor_threshold_nm=floor_threshold,status=any(mask) ? "bright_support" : "no_bright_support")
end

"All eight-connected components, including small, edge and missing-adjacent ones."
function components(mask,z,xs,ys)
    labels=zeros(Int,size(mask)); rows=NamedTuple[]; ny,nx=size(mask)
    for y in 1:ny,x in 1:nx
        mask[y,x] && labels[y,x]==0 || continue
        id=length(rows)+1; todo=[CartesianIndex(y,x)]; labels[y,x]=id; k=1; edge=false; holes=false
        while k<=length(todo)
            iy,ix=Tuple(todo[k]); k+=1
            edge|=iy==1 || iy==ny || ix==1 || ix==nx
            for yy in max(1,iy-1):min(ny,iy+1),xx in max(1,ix-1):min(nx,ix+1)
                isfinite(z[yy,xx]) || (holes=true)
                if mask[yy,xx] && labels[yy,xx]==0
                    labels[yy,xx]=id; push!(todo,CartesianIndex(yy,xx))
                end
            end
        end
        push!(rows,(;component=id,pixels=length(todo),first_row=y,first_column=x,
            centroid_x_nm=mean(xs[i[2]] for i in todo),centroid_y_nm=mean(ys[i[1]] for i in todo),
            mean_height_nm=mean(z[todo]),max_height_nm=maximum(z[todo]),touches_edge=edge,touches_missing=holes))
    end
    (;labels,rows)
end

"Box guard, rounded outward in each physical image coordinate; never a component filter."
function box_guard(mask,rx,ry)
    ny,nx=size(mask); rows=falses(ny,nx); out=falses(ny,nx)
    for y in 1:ny
        prefix=vcat(0,cumsum(Int.(mask[y,:])))
        for x in 1:nx; rows[y,x]=prefix[min(nx,x+rx)+1]-prefix[max(1,x-rx)]>0; end
    end
    for x in 1:nx
        prefix=vcat(0,cumsum(Int.(rows[:,x])))
        for y in 1:ny; out[y,x]=prefix[min(ny,y+ry)+1]-prefix[max(1,y-ry)]>0; end
    end
    out
end

"One exact joint OLS solve: x slope plus free row levels on supplied source background."
function joint_background(raw_centered,xs,guard,c)
    s=c["selection"]; ny,nx=size(raw_centered)
    size(guard)==(ny,nx) && length(xs)==nx && nx>=2 && all(isfinite,xs) && all(diff(xs).>0) || error("Invalid background coordinates")
    xcenter=first(xs)+(last(xs)-first(xs))/2; xscale=last(xs)-first(xs); xn=(xs.-xcenter)./xscale
    bg=BitMatrix(isfinite.(raw_centered).&.!guard)
    counts=[count(@view bg[y,:]) for y in 1:ny]
    usable=counts.>=s["min_row_background_pixels"]
    used=copy(bg); used[.!usable,:].=false
    n=count(used); xmean=fill(NaN,ny); zmean=fill(NaN,ny); level=fill(NaN,ny)
    xx=0.; xz=0.
    for y in findall(usable)
        inds=findall(@view used[y,:]); x=xn[inds]; z=raw_centered[y,inds]
        xmean[y]=mean(x); zmean[y]=mean(z)
        dx=x.-xmean[y]; dz=z.-zmean[y]
        xx+=sum(abs2,dx); xz+=sum(dx.*dz)
    end
    status = n<max(s["min_background_pixels"],s["min_background_fraction"]*length(bg)) ? "insufficient_background" :
        !(isfinite(xx) && isfinite(xz) && xx>c["preprocessing"]["rank_rtol"]*n) ? "rank_failure" : "ok"
    slope=status=="ok" ? xz/xx : NaN
    corrected=fill(NaN,size(raw_centered))
    if status=="ok"
        for y in findall(usable)
            level[y]=zmean[y]-slope*xmean[y]
            for x in 1:nx
                isfinite(raw_centered[y,x]) && (corrected[y,x]=raw_centered[y,x]-level[y]-slope*xn[x])
            end
        end
    end
    rows=[(;row=y,background_pixels=counts[y],usable=status=="ok" && usable[y],
        xmean=xmean[y],zmean_centered_nm=zmean[y],row_level_centered_nm=level[y]) for y in 1:ny]
    (;corrected,used,rows,slope,xcenter,xscale,within_xx=xx,within_xz=xz,status,background_pixels=n)
end

"Two fixed passes, no target, class, count, mask iteration or post-result retuning."
function analyse(raw,initial_smooth,xs,ys,c)
    size(raw)==size(initial_smooth)==(length(ys),length(xs)) || error("Different image grids")
    length(xs)>=2 && length(ys)>=2 && all(isfinite,xs) && all(isfinite,ys) &&
        all(diff(xs).>0) && all(diff(ys).>0) || error("Invalid coordinates")
    observed=isfinite.(raw); any(observed) || error("No raw observations")
    reference=median(raw[observed]); centered=raw.-reference
    noise=hf_scale(raw); initial=segment(initial_smooth,noise,c)
    pad=c["model"]["guard_padding_nm"]
    rx=ceil(Int,pad/(xs[2]-xs[1])); ry=ceil(Int,pad/(ys[2]-ys[1]))
    guard=box_guard(initial.mask,rx,ry)
    background=joint_background(centered,xs,guard,c)
    final_smooth=STMSXMIO._box_smooth(background.corrected,c["preprocessing"]["smooth_radius_px"])
    final=segment(final_smooth,noise,c)
    initial_components=components(initial.mask,initial_smooth,xs,ys)
    final_components=components(final.mask,final_smooth,xs,ys)
    (;reference,centered,noise,initial,guard,background,final_smooth,final,initial_components,final_components,rx,ry)
end
end
