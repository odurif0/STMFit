module TangentMoldProjection

using LinearAlgebra, TOML
export settings, validate_geometry, native_design, projector, project, projected_score

function settings(path)
    c = TOML.parsefile(path)
    Set(keys(c)) == Set(["model", "selection", "preprocessing"]) || error("Unexpected tangent settings section")
    Set(keys(c["model"])) == Set(["basis", "rank_rtol", "zero_norm_rtol"]) || error("Unexpected tangent model key")
    sk = Set(keys(c["selection"]))
    sk in (Set(["minimum_observed_pixels", "minimum_observed_fraction"]),
           Set(["minimum_observed_pixels", "minimum_observed_fraction", "missing_cost"])) || error("Unexpected tangent selection key")
    c["preprocessing"] == Dict("sampling" => "native_box_bilinear") || error("Native tangent sampling required")
    m, s = c["model"], c["selection"]
    m["basis"] in ("target_gaussian_affine", "target_gaussian_affine_orientation",
                   "target_gaussian_affine_adjacent_amplitudes") || error("Unsupported tangent basis")
    # Absence preserves the historical configurations byte-for-byte. New paired
    # experiments explicitly declare the previously tested missing-cost repair.
    missing_cost = get(s, "missing_cost", "legacy_infinite")
    missing_cost in ("legacy_infinite", "omit_both_infinite") || error("Unsupported missing-cost policy")
    for k in ("rank_rtol", "zero_norm_rtol")
        x = m[k]
        x isa Real && !(x isa Bool) && isfinite(x) && 0 < x < 1 || error("Invalid $k")
    end
    s["minimum_observed_pixels"] === 5 && s["minimum_observed_fraction"] == 0.5 ||
        error("The native connected-mold observation guard is frozen")
    return (rank_rtol=Float64(m["rank_rtol"]), zero_norm_rtol=Float64(m["zero_norm_rtol"]),
            minimum_pixels=5, minimum_fraction=0.5,
            orientation=m["basis"] == "target_gaussian_affine_orientation",
            adjacent_amplitudes=m["basis"] == "target_gaussian_affine_adjacent_amplitudes",
            omit_unavailable=missing_cost == "omit_both_infinite")
end

number(r, k) = parse(Float64, r[k])
function validate_geometry(r)
    all(isfinite(number(r,k)) for k in ("x_nm","y_nm","axis_x","axis_y","sigma_parallel_nm","sigma_perp_nm")) ||
        error("Nonfinite tangent geometry")
    all(number(r,k)>0 for k in ("sigma_parallel_nm","sigma_perp_nm")) || error("Invalid tangent width")
    abs(number(r,"axis_x")^2 + number(r,"axis_y")^2 - 1) <= 3e-8 || error("Nonunit cached axis")
    parse(Float64,get(r,"skew_ratio", "1.0")) == 1.0 || error("Only symmetric Gaussian base geometry is supported")
    get(r,"source", "") in ("circ", "ell") || error("Unknown base Gaussian family")
    r["source"] == "circ" && number(r,"sigma_parallel_nm") != number(r,"sigma_perp_nm") && error("Circular widths differ")
    any(startswith(k,"model_axis") || k in ("model_orientation","profile_view") for k in keys(r)) &&
        error("Only unchanged global-axis geometry is supported")
    return r
end

# Rescaled derivatives span the same tangent space as exact derivatives with
# respect to amplitude, center t/u and log widths. Circular fits have one width.
# The amplitude feature itself remains unchanged in the downstream classifier.
has_orientation(r) = r["source"] == "ell" && number(r,"sigma_parallel_nm") != number(r,"sigma_perp_nm")

function basis_at(x,y,r; orientation::Bool=false)
    dx=x-number(r,"x_nm"); dy=y-number(r,"y_nm")
    t=dx*number(r,"axis_x")+dy*number(r,"axis_y")
    u=-dx*number(r,"axis_y")+dy*number(r,"axis_x")
    a=t/number(r,"sigma_parallel_nm"); b=u/number(r,"sigma_perp_nm")
    g=exp(-0.5*(a*a+b*b))
    affine=[1.0,t,u,g,a*g,b*g]
    base = r["source"] == "circ" ? vcat(affine,(a*a+b*b)*g) : vcat(affine,a*a*g,b*b*g)
    if orientation && has_orientation(r)
        sp=number(r,"sigma_parallel_nm"); su=number(r,"sigma_perp_nm")
        # dG/dtheta at a fixed native pixel and fixed sampling frame. Factoring
        # sp^2-su^2 avoids cancellation for nearly circular ellipses. The exact
        # zero derivative is omitted only at equality, with no anisotropy cutoff.
        push!(base,a*b*g*((sp-su)*(sp+su)/(sp*su)))
    end
    return base
end

"Native clipped box smoothing followed by the extractor's bilinear sampling."
function native_design(xs,ys,r,coords,radius::Int; orientation::Bool=false, neighbors=AbstractDict[])
    validate_geometry(r)
    foreach(validate_geometry, neighbors)
    length(neighbors) <= 2 || error("At most two adjacent amplitude columns")
    orientation && !isempty(neighbors) && error("Orientation/neighbor combination is outside this experiment")
    radius >= 0 || error("Negative smoothing radius")
    length(xs)>=2 && length(ys)>=2 && all(isfinite,xs) && all(isfinite,ys) &&
        all(>(0),diff(xs)) && all(>(0),diff(ys)) || error("Invalid native image grid")
    !isempty(coords) && all(isfinite,coords) || error("Invalid patch coordinates")
    p=(r["source"]=="circ" ? 7 : 8) + Int(orientation && has_orientation(r)) + length(neighbors)
    B=fill(NaN,length(coords)^2,p)
    cache=Dict{Tuple{Int,Int},Vector{Float64}}()
    function smoothed(ix,iy)
        get!(cache,(ix,iy)) do
            v=zeros(p); n=0
            for j in max(1,iy-radius):min(length(ys),iy+radius), i in max(1,ix-radius):min(length(xs),ix+radius)
                b = basis_at(xs[i],ys[j],r;orientation)
                # d(A_neighbor*G_neighbor)/dA_neighbor = G_neighbor. Evaluate
                # on the SAME native pixels, then sample the TARGET's frame.
                for q in neighbors
                    dx=xs[i]-number(q,"x_nm"); dy=ys[j]-number(q,"y_nm")
                    t=dx*number(q,"axis_x")+dy*number(q,"axis_y")
                    u=-dx*number(q,"axis_y")+dy*number(q,"axis_x")
                    push!(b,exp(-0.5*((t/number(q,"sigma_parallel_nm"))^2+(u/number(q,"sigma_perp_nm"))^2)))
                end
                v .+= b; n+=1
            end
            v/n
        end
    end
    k=0
    for u in coords, t in coords
        k+=1
        x=number(r,"x_nm")+t*number(r,"axis_x")-u*number(r,"axis_y")
        y=number(r,"y_nm")+t*number(r,"axis_y")+u*number(r,"axis_x")
        ix=searchsortedlast(xs,x); iy=searchsortedlast(ys,y)
        (1<=ix<length(xs) && 1<=iy<length(ys)) || continue
        tx=(x-xs[ix])/(xs[ix+1]-xs[ix]); ty=(y-ys[iy])/(ys[iy+1]-ys[iy])
        B[k,:]=(1-tx)*(1-ty)*smoothed(ix,iy)+tx*(1-ty)*smoothed(ix+1,iy)+
            (1-tx)*ty*smoothed(ix,iy+1)+tx*ty*smoothed(ix+1,iy+1)
    end
    return B
end

"Topological previous/next lobes, never chemical classes or nearest-by-signal selection."
function adjacent_rows(base, file, lobe, n)
    1 <= lobe <= n || error("Invalid target lobe")
    return [base[(file,j)] for j in (lobe-1,lobe+1) if 1 <= j <= n]
end

"Unit-normalize columns before the fixed numerical SVD rank decision."
function projector(B, opt)
    all(isfinite,B) && size(B,1)>size(B,2) || error("Invalid tangent support")
    norms=[norm(B[:,j]) for j in axes(B,2)]
    all(>(0),norms) || error("Empty tangent column")
    D=B ./ reshape(norms,1,:)
    f=svd(D;full=false)
    rank=count(>(opt.rank_rtol*first(f.S)),f.S)
    Q=f.U[:,1:rank]
    return (;Q,rank,smin=last(f.S),smax=first(f.S))
end
project(v,p)=v-p.Q*(p.Q'*v)

"One common observed support and projector for the observation and both molds."
function projected_score(patch, templates, B, opt)
    n=length(patch)
    size(B,1)==n && size(templates,1)==n && size(templates,2)==2 || error("Tangent dimensions differ")
    all(isfinite,templates) || error("Nonfinite physical template")
    good=findall(isfinite,patch)
    if length(good)<max(opt.minimum_pixels,ceil(Int,opt.minimum_fraction*n))
        return (costs=[Inf,Inf],reason="insufficient_native_support",n=length(good),rank=0,
                patch_fraction=NaN,contrast_fraction=NaN,smin=NaN,smax=NaN)
    end
    p=projector(B[good,:],opt)
    y=patch[good]; M=templates[good,:]
    yp=project(y,p); Mp=project(M,p)
    c=M[:,2]-M[:,1]; cp=Mp[:,2]-Mp[:,1]
    fraction(v,w)=norm(v)>0 ? sum(abs2,w)/sum(abs2,v) : NaN
    norm(c)>0 && norm(cp)>opt.zero_norm_rtol*norm(c) || error("Numerically annihilated physical contrast")
    all(norm(Mp[:,j])>opt.zero_norm_rtol*norm(M[:,j]) for j in 1:2) || error("Numerically annihilated physical mold")
    diag=(n=length(good),rank=p.rank,patch_fraction=fraction(y,yp),contrast_fraction=fraction(c,cp),smin=p.smin,smax=p.smax)
    norm(y)>0 && norm(yp)>opt.zero_norm_rtol*norm(y) ||
        return merge((costs=[Inf,Inf],reason="zero_projected_patch"),diag)
    costs=[-dot(yp,Mp[:,j])/(norm(yp)*norm(Mp[:,j])) for j in 1:2]
    return merge((;costs,reason="ok"),diag)
end
end
