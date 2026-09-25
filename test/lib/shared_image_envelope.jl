module SharedImageEnvelope
using LinearAlgebra, SparseArrays, TOML

function settings(path)
    c=TOML.parsefile(path)
    fields=Dict("model"=>["foreground_run","transfer_run","controls_run","foreground_settings_sha256",
        "transfer_settings_sha256","controls_settings_sha256","files_sha256","basis","minimum_knot_intervals_nm",
        "interval_count","background","probe"],
        "selection"=>["policy","degrees_of_freedom","foreground_use","groups","probe_use"],
        "preprocessing"=>["observations","target","missing","unsupported_basis","rank_pivot_rtol","normal_equation_rtol","binary_arrays"])
    Set(keys(c))==Set(keys(fields)) || error("Unexpected sections")
    for (s,k) in fields; Set(keys(c[s]))==Set(k) || error("Unknown/missing $s setting"); end
    m,s,p=c["model"],c["selection"],c["preprocessing"]
    m["basis"]=="whole_image_open_uniform_tensor_cubic_bspline" &&
        m["interval_count"]=="max_one_floor_image_span_over_minimum_interval" &&
        m["background"]=="reuse_source_masked_x_slope_and_row_levels" &&
        m["probe"]=="saved_common_geometry_glcnac_minus_glcn_zero_outside_patch" || error("Changed model")
    widths=m["minimum_knot_intervals_nm"]
    length(widths)==3 && all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>0,widths) &&
        all(diff(widths).<0) || error("Invalid fixed knot candidates")
    s["policy"]=="minimum_source_conditional_gcv_background_then_coarse_to_fine" &&
        s["degrees_of_freedom"]=="trace_B_plus_P_minus_PB_conditional_on_frozen_mask" &&
        s["foreground_use"]=="report_only_all_finite_source_pixels_fit" &&
        s["groups"]==["all","foreground","background"] &&
        s["probe_use"]=="post_fit_only_no_selection_from_retention" || error("Changed selection")
    p["observations"]=="saved_unsmoothed_source_native_pixels" &&
        p["target"]=="saved_calibration_then_source_background_subtraction_score_only" &&
        p["missing"]=="keep_missing_no_prediction_score_or_probe_imputation" &&
        p["unsupported_basis"]=="drop_exact_zero_columns_no_change_to_observed_function_space" &&
        p["binary_arrays"]=="little_endian_column_major_float64" || error("Changed observation boundary")
    all(occursin(r"^[0-9a-f]{64}$",m[k]) for k in ("foreground_settings_sha256","transfer_settings_sha256","controls_settings_sha256","files_sha256")) || error("Invalid hash")
    all(p[k] isa Real && !(p[k] isa Bool) && 0<p[k]<1 for k in ("rank_pivot_rtol","normal_equation_rtol")) || error("Invalid tolerance")
    c
end

"Open clamped cubic B-splines: iterative Cox--de Boor, no extrapolation."
function basis(axis,width)
    length(axis)>=2 && all(isfinite,axis) && all(diff(axis).>0) && width>0 || error("Invalid axis")
    lo,hi=first(axis),last(axis); segments=max(1,floor(Int,(hi-lo)/width))
    breaks=collect(range(lo,hi;length=segments+1))
    knots=vcat(fill(lo,3),breaks,fill(hi,3)); nb=length(knots)-4
    a=zeros(length(axis),nb)
    for (r,x) in enumerate(axis)
        if x==hi; a[r,end]=1.; continue; end
        v=[knots[j]<=x<knots[j+1] ? 1. : 0. for j in 1:length(knots)-1]
        for degree in 1:3
            w=zeros(length(v)-1)
            for j in eachindex(w)
                left=knots[j+degree]-knots[j]; right=knots[j+degree+1]-knots[j+1]
                left>0 && (w[j]+=(x-knots[j])/left*v[j])
                right>0 && (w[j]+=(knots[j+degree+1]-x)/right*v[j+1])
            end
            v=w
        end
        a[r,:]=v
    end
    (;matrix=sparse(a),knots,segments,step=(hi-lo)/segments)
end

"Background operator on the *same* finite, usable-row source support."
function background_operator(valid,guard,xn)
    size(valid)==size(guard) && length(xn)==size(valid,2) || error("Different background grids")
    inds=findall(valid); rows=sort(unique(q[1] for q in inds)); rowid=Dict(y=>j for (j,y) in enumerate(rows))
    N=length(inds); R=length(rows); R>0 || error("No usable row")
    ii=repeat(collect(1:N),2); jj=vcat([rowid[q[1]] for q in inds],fill(R+1,N))
    Z=sparse(ii,jj,vcat(ones(N),[xn[q[2]] for q in inds]),N,R+1)
    bg=BitVector([!guard[q] for q in inds]); Zb=Z[bg,:]
    counts=vec(Array(sum(Zb[:,1:R];dims=1))); all(counts.>0) || error("Empty background row")
    sx=Vector(transpose(Zb[:,1:R])*Zb[:,end]); means=sx./counts
    xx=sum((xn[inds[j][2]]-means[rowid[inds[j][1]]])^2 for j in findall(bg))
    isfinite(xx) && xx>0 || error("Background slope unidentifiable")
    (;Z,Zb,bg,counts,means,xx,inds,rows,rank=R+1)
end

"Structured solve of the row-level/x-slope background normal matrix."
function bg_solve(b,rhs::AbstractMatrix)
    R=length(b.counts); size(rhs,1)==R+1 || error("Wrong background RHS")
    slope=(rhs[end,:]-vec(transpose(b.means)*rhs[1:R,:]))/b.xx
    levels=rhs[1:R,:]./b.counts .- b.means*transpose(slope)
    vcat(levels,transpose(slope))
end
bg_solve(b,rhs::AbstractVector)=vec(bg_solve(b,reshape(rhs,:,1)))
background(b,v)=b.Z*bg_solve(b,transpose(b.Zb)*v[b.bg])

function system(X,c)
    g=transpose(X)*X; scales=sqrt.(diag(g))
    all(isfinite,scales) && all(scales.>0) || error("Unsupported spline coefficient")
    d=spdiagm(0=>1 ./scales); factor=cholesky(Symmetric(d*g*d))
    pivot=minimum(abs2,diag(sparse(factor.L)))
    pivot>c["preprocessing"]["rank_pivot_rtol"] || error("Spline pivot rank failure")
    (;g,scales,factor,pivot)
end
solve(s,rhs::AbstractVector)=(s.factor\(rhs./s.scales))./s.scales
solve(s,rhs::AbstractMatrix)=(s.factor\(rhs./s.scales))./s.scales
gcv(rss,n,df)=isfinite(df) && 0<=df<n ? (rss/n)/(1-df/n)^2 : Inf

function fit_candidate(z,valid,xs,ys,b,width,c)
    by=basis(ys,width); bx=basis(xs,width)
    full=kron(bx.matrix,by.matrix)[vec(valid),:]
    # An identically zero column has no parameter on observed pixels. Keep its
    # index for provenance, never estimate or predict in its missing region.
    active=findall(vec(sum(abs2,full;dims=1)).>0)
    X=full[:,active]; N,P=size(X)
    P<N || error("Insufficient spline data")
    sys=system(X,c); rhs=transpose(X)*z; beta=solve(sys,rhs); pred=X*beta
    normal=norm(transpose(X)*(z-pred),Inf)/max(norm(rhs,Inf),eps(Float64))
    normal<=c["preprocessing"]["normal_equation_rtol"] || error("Spline normal equation failure")
    TX=bg_solve(b,Matrix(transpose(b.Zb)*X[b.bg,:]))
    GinvXZ=solve(sys,Matrix(transpose(X)*b.Z))
    trace_pb=sum(TX.*transpose(GinvXZ)); df=b.rank+P-trace_pb
    rss=sum(abs2,z-pred)
    (;X,sys,beta,pred,width,active,segments_x=bx.segments,segments_y=by.segments,
        step_x=bx.step,step_y=by.step,knots_x=bx.knots,knots_y=by.knots,
        parameters=P,trace_pb,df,rss,gcv=gcv(rss,N,df),normal)
end

"No target, masks used for scores, chemical map or anchor enters this fit."
function fit_all(z,valid,xs,ys,b,c)
    N=length(z); N==count(valid) || error("Changed finite population")
    null=(;pred=zeros(N),beta=Float64[],active=Int[],width=0.,parameters=0,trace_pb=0.,df=Float64(b.rank),
        rss=sum(abs2,z),gcv=gcv(sum(abs2,z),N,b.rank),normal=0.,segments_x=0,segments_y=0,step_x=NaN,step_y=NaN)
    fits=Any[null]; statuses=["ok"]
    for width in c["model"]["minimum_knot_intervals_nm"]
        try
            result=fit_candidate(z,valid,xs,ys,b,width,c)
            isfinite(result.gcv) || error("Nonfinite conditional GCV")
            push!(fits,result); push!(statuses,"ok")
        catch e
            e isa InterruptException && rethrow()
            push!(fits,nothing); push!(statuses,"failed:"*replace(sprint(showerror,e),'\n'=>' ','\t'=>' '))
        end
    end
    chosen=argmin([f===nothing ? Inf : f.gcv for f in fits])
    isfinite(fits[chosen].gcv) || error("No valid model")
    (;fits,statuses,chosen)
end

"Local sensitivity only: frozen mask, background rule and selected spline resolution."
function probe(delta,b,fit)
    bb=background(b,delta); corrected=delta-bb
    beta=fit.parameters==0 ? Float64[] : solve(fit.sys,transpose(fit.X)*corrected)
    pred=fit.parameters==0 ? zeros(length(delta)) : fit.X*beta
    (;beta,input_sse=sum(abs2,delta),background_sse=sum(abs2,bb),
        corrected_sse=sum(abs2,corrected),envelope_sse=sum(abs2,pred),residual_sse=sum(abs2,corrected-pred))
end
end
