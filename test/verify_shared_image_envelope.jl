#!/usr/bin/env julia
# Saved-coefficient/operator verification, not a new experiment or model selection.
module VerifySharedImageEnvelope
using Test, LinearAlgebra, SparseArrays, Statistics, TOML
include(joinpath(@__DIR__,"diagnose_shared_image_envelope.jl"))
const D=DiagnoseSharedImageEnvelope
const R=D.R
const t=D.t
const n=D.n
const i=D.i
close(a,b)=isequal(a,b) || isapprox(a,b;rtol=1e-8,atol=1e-12)

function recursive_basis(x,t,j,d)
    d==0 && return t[j]<=x<t[j+1] ? 1. : 0.
    a=t[j+d]-t[j]; b=t[j+d+1]-t[j+1]
    (a>0 ? (x-t[j])/a*recursive_basis(x,t,j,d-1) : 0.)+
        (b>0 ? (t[j+d+1]-x)/b*recursive_basis(x,t,j+1,d-1) : 0.)
end
function basis_reference(axis,width)
    segments=max(1,floor(Int,(last(axis)-first(axis))/width))
    knots=vcat(fill(first(axis),3),collect(range(first(axis),last(axis);length=segments+1)),fill(last(axis),3))
    nb=segments+3
    B=[x==last(axis) ? Float64(j==nb) : recursive_basis(x,knots,j,3) for x in axis,j in 1:nb]
    @test all(isapprox.(sum(B;dims=2),1.;atol=1e-13,rtol=0))
    sparse(B)
end

"Four weights, independent of the inference interpolator."
function map_value(map,state,dx,dy,base)
    co,si=cosd(state.angle),sind(state.angle)
    x=co*(dx-state.tx)+si*(dy-state.ty); y=-si*(dx-state.tx)+co*(dy-state.ty)
    half,step=base["model"]["map_half_nm"],base["model"]["map_step_nm"]
    u=(x+half)/step+1; v=(y+half)/step+1; nx=size(map,2); ny=size(map,1)
    1<=u<=nx && 1<=v<=ny || error("Probe outside fixed map")
    ix=min(floor(Int,u),nx-1); iy=min(floor(Int,v),ny-1); a=u-ix; b=v-iy
    (1-a)*(1-b)*map[iy,ix]+a*(1-b)*map[iy,ix+1]+(1-a)*b*map[iy+1,ix]+a*b*map[iy+1,ix+1]
end

function verify_direction(dir,paths,maps,base,c)
    a=D.load_source(paths["foreground"]); other=D.load_source(paths["target"])
    cal=TOML.parsefile(joinpath(dir,"calibration.toml")); oldcal=TOML.parsefile(joinpath(paths["foreground"],"pair_status.toml"))
    @test cal==oldcal
    for h in t(joinpath(dir,"input_hashes.tsv")); @test R.sha(joinpath(paths[h["role"]],h["file"]))==h["sha256"]; end
    meta=a.meta; ny,nx=size(a.raw); good=a.valid; N=count(good)
    arrays=Dict{String,Matrix{Float64}}()
    for row in t(joinpath(dir,"arrays.tsv"))
        path=joinpath(dir,row["file"]); @test filesize(path)==i(row,"bytes")==8ny*nx
        @test R.sha(path)==row["sha256"]
        arrays[row["file"]]=D.D.get_array(dir,row["file"],Float64,ny,nx)
    end
    @test all(isequal.(arrays["source_corrected.f64"][good],a.z))
    @test all(isnan,arrays["source_corrected.f64"][.!good])
    target=arrays["target_corrected.f64"]; masks=R.F.row_masks(ny,base); plane=cal["difference_plane"]; lag=cal["applied_dx_px"]
    @test !any(masks.cal .& masks.test)
    for y in 1:ny,x in 1:nx
        xx=x+lag
        if good[y,x] && masks.test[y] && 1<=xx<=nx && isfinite(other.raw[y,xx])
            bg=n(a.rows[y],"row_level_centered_nm")+meta["slope_normalized_nm"]*a.xn[x]
            want=other.raw[y,xx]+other.meta["source_reference_nm"]-meta["source_reference_nm"]-
                plane[1]-plane[2]*meta["xs_nm"][x]-plane[3]*meta["ys_nm"][y]-bg
            @test close(target[y,x],want)
        else; @test isnan(target[y,x]); end
    end
    models=t(joinpath(dir,"models.tsv")); coeff=t(joinpath(dir,"coefficients.tsv")); selection=TOML.parsefile(joinpath(dir,"selection.toml"))
    @test length(models)==4 && count(r->r["chosen"]=="true",models)==1
    @test selection["chosen"]==only(r for r in models if r["chosen"]=="true")["arm"]
    selected=argmin(n.(models,"gcv")); @test models[selected]["chosen"]=="true"
    # Dense background solve is independent of the production structured row elimination.
    Z=a.b.Z; Zbg=Z[a.b.bg,:]; bgfactor=cholesky(Symmetric(Matrix(transpose(Zbg)*Zbg)))
    @test size(Z,2)==a.b.rank
    @test isapprox(Z*(bgfactor\(transpose(Zbg)*a.raw[good][a.b.bg])),a.bg;rtol=1e-9,atol=1e-12)
    projections=Dict{String,Any}()
    for (k,row) in enumerate(models)
        arm=row["arm"]; @test i(row,"pixels")==N && i(row,"background_rank")==a.b.rank
        row["status"]=="ok" || continue
        pred=arrays[arm*".f64"]; @test all(isnan,pred[.!good])
        P=i(row,"parameters"); rss=sum(abs2,a.z-pred[good])
        @test close(rss,n(row,"source_sse_nm2"))
        if P==0
            @test arm=="background_only" && all(iszero,pred[good])
            @test n(row,"degrees_of_freedom")==a.b.rank && n(row,"trace_pb")==0
            projections[arm]=nothing
        else
            width=c["model"]["minimum_knot_intervals_nm"][k-1]; @test width==n(row,"width_nm")
            Xfull=kron(basis_reference(meta["xs_nm"],width),basis_reference(meta["ys_nm"],width))[vec(good),:]
            active=findall(vec(sum(abs2,Xfull;dims=1)).>0); X=Xfull[:,active]
            @test size(X,2)==P
            beta=n.(filter(r->r["arm"]==arm,coeff),"value_nm"); @test length(beta)==P
            @test i.(filter(r->r["arm"]==arm,coeff),"coefficient")==collect(1:P)
            @test i.(filter(r->r["arm"]==arm,coeff),"basis_column")==active
            replay=X*beta
            @test isapprox(replay,pred[good];rtol=1e-9,atol=1e-12)
            normal=norm(transpose(X)*(a.z-replay),Inf)/max(norm(transpose(X)*a.z,Inf),eps(Float64))
            @test normal<=c["preprocessing"]["normal_equation_rtol"]
            G=transpose(X)*X; factor=cholesky(Symmetric(G))
            TX=bgfactor\Matrix(transpose(Zbg)*X[a.b.bg,:])
            GXZ=factor\Matrix(transpose(X)*Z); overlap=sum(transpose(TX).*GXZ)
            df=size(Z,2)+P-overlap
            @test close(overlap,n(row,"trace_pb")) && close(df,n(row,"degrees_of_freedom"))
            projections[arm]=(;X,factor)
        end
        @test close((rss/N)/(1-n(row,"degrees_of_freedom")/N)^2,n(row,"gcv"))
    end
    rows=t(joinpath(dir,"scores.tsv")); @test length(rows)==3(count(r->r["status"]=="ok",models)+1)
    for row in rows
        arm=row["arm"]; pred=arm=="source_copy" ? arrays["source_corrected.f64"] : arrays[arm*".f64"]
        group=row["group"]; mask=good .& isfinite.(target)
        group=="foreground" && (mask .&=a.mask)
        group=="background" && (mask .&=.!a.mask)
        @test i(row,"test_pixels")==count(mask)
        @test close(n(row,"test_sse_nm2"),sum(abs2,target[mask]-pred[mask]))
        @test parse(Bool,row["identified"])==cal["identified"]
        @test parse(Bool,row["chosen"])==(arm==selection["chosen"])
    end
    probes=t(joinpath(dir,"probes.tsv"))
    if isfile(joinpath(paths["controls"],"patches.tsv"))
        old=filter(r->r["arm"]=="plane_common",t(joinpath(paths["controls"],"patches.tsv")))
        membership=Dict(i(r,"anchor")=>r for r in t(joinpath(paths["foreground"],"anchor_membership.tsv")))
        expected=filter(r->membership[i(r,"anchor")]["joint_core"]=="true",old)
        @test Set((i(r,"interval"),i(r,"patch")) for r in expected)==Set((i(r,"interval"),i(r,"patch")) for r in probes)
        old=Dict((i(r,"interval"),i(r,"patch"))=>r for r in old)
        obs=t(joinpath(paths["transfer"],"observations.tsv")); groups=Dict{Int,Vector{eltype(obs)}}()
        for r in obs; push!(get!(groups,i(r,"patch"),eltype(obs)[]),r); end
        pc=t(joinpath(dir,"probe_coefficients.tsv")); index=zeros(Int,ny,nx); index[good]=collect(1:N)
        coeffindex=Dict{Tuple{Int,Int},Vector{eltype(pc)}}()
        for r in pc; push!(get!(coeffindex,(i(r,"interval"),i(r,"patch")),eltype(pc)[]),r); end
        projection=projections[selection["chosen"]]
        for row in probes
            interval,patch=i(row,"interval"),i(row,"patch"); rr=groups[patch]
            inds=[index[i(r,"row"),i(r,"column")] for r in rr]
            @test length(inds)==i(row,"pixels")
            if any(iszero,inds); @test row["status"]=="missing_source_support"; continue; end
            @test row["status"]=="ok"
            o=old[(interval,patch)]; state=(;angle=n(o,"angle"),tx=n(o,"tx"),ty=n(o,"ty"))
            h=maps.heights[interval]; delta=zeros(N)
            delta[inds]=[map_value(h[1],state,n(r,"dx_nm"),n(r,"dy_nm"),base)-
                map_value(h[0],state,n(r,"dx_nm"),n(r,"dy_nm"),base) for r in rr]
            bg=Z*(bgfactor\(transpose(Zbg)*delta[a.b.bg])); corrected=delta-bg
            pred=zeros(N)
            if projection!==nothing
                beta_rows=coeffindex[(interval,patch)]; beta=n.(beta_rows,"value_nm")
                @test i.(beta_rows,"coefficient")==collect(1:size(projection.X,2))
                pred=projection.X*beta
                normal=transpose(projection.X)*(corrected-pred)
                @test norm(normal,Inf)<=1e-9max(norm(transpose(projection.X)*corrected,Inf),eps(Float64))
            end
            for (key,v) in (("input_sse_nm2",sum(abs2,delta)),("background_sse_nm2",sum(abs2,bg)),
                ("corrected_sse_nm2",sum(abs2,corrected)),("envelope_sse_nm2",sum(abs2,pred)),
                ("residual_sse_nm2",sum(abs2,corrected-pred)))
                @test close(v,n(row,key))
            end
            @test close(sum(abs2,corrected-pred)/sum(abs2,delta),n(row,"remaining_fraction"))
        end
    else
        @test isempty(probes) && selection["probe_status"]=="upstream_no_supported_physical_patches"
    end
    @test selection["probes"]==length(probes)
    (;models,rows,probes,identified=cal["identified"])
end

function main(args=ARGS)
    length(args)==3 || error("verify_shared_image_envelope.jl REPO SAVED_RUN NEW_REPORT_DIR")
    root,run,out=args; c=D.F.settings(joinpath(run,"settings.toml")); R.S.newdir(out); BLAS.set_num_threads(1)
    paths=Dict(k=>joinpath(root,c["model"][k*"_run"]) for k in ("foreground","transfer","controls"))
    base=R.F.settings(joinpath(paths["transfer"],"settings.toml")); maps=R.load_maps(root,base)
    scores=NamedTuple[]; models=NamedTuple[]; probes=NamedTuple[]; statuses=NamedTuple[]
    @testset "Shared image envelope saved evidence" begin
        for (role,path) in paths
            @test R.sha(joinpath(path,"settings.toml"))==c["model"][role*"_settings_sha256"]
            @test R.sha(joinpath(path,"files.tsv"))==c["model"]["files_sha256"]
        end
        files=getindex.(t(joinpath(run,"files.tsv")),"file")
        @test length(unique(files))==length(files)
        for file in files
            stem=splitext(file)[1]; ss=t(joinpath(run,stem,"status.tsv"))
            @test Set(getindex.(ss,"direction"))==Set(["fwd","bwd"]) && length(ss)==2
            for status in ss
                view=status["direction"]; push!(statuses,(;file,view,status=status["status"]))
                status["status"]=="ok" || continue
                p=Dict(k=>joinpath(v,stem,view) for (k,v) in paths)
                p["target"]=joinpath(paths["foreground"],stem,view=="fwd" ? "bwd" : "fwd")
                r=verify_direction(joinpath(run,stem,view),p,maps,base,c)
                for row in r.rows
                    q=(;file,view,arm=row["arm"],group=row["group"],identified=r.identified,
                        test_pixels=i(row,"test_pixels"),test_sse_nm2=n(row,"test_sse_nm2"))
                    push!(scores,q)
                    row["chosen"]=="true" && push!(scores,merge(q,(;arm="gcv_selected")))
                end
                append!(models,[(;file,view,identified=r.identified,arm=row["arm"],status=row["status"],
                    chosen=row["chosen"]=="true",width_nm=n(row,"width_nm"),gcv=n(row,"gcv"),df=n(row,"degrees_of_freedom")) for row in r.models])
                append!(probes,[(;file,view,identified=r.identified,interval=i(row,"interval"),patch=i(row,"patch"),status=row["status"],
                    input_sse_nm2=n(row,"input_sse_nm2"),residual_sse_nm2=n(row,"residual_sse_nm2"),remaining_fraction=n(row,"remaining_fraction")) for row in r.probes])
                @test i(status,"probes")==length(r.probes)
            end
            println("VERIFIED ",file); flush(stdout)
            GC.gc()
        end
    end
    R.table(joinpath(out,"statuses.tsv"),statuses); R.table(joinpath(out,"scores.tsv"),scores)
    R.table(joinpath(out,"models.tsv"),models); R.table(joinpath(out,"probes.tsv"),probes)
    summary=NamedTuple[]; comparisons=NamedTuple[]; retention=NamedTuple[]
    for subset in ("all","identified")
        subsetrows=filter(r->subset=="all" || r.identified,scores)
        for group in c["selection"]["groups"],arm in ("background_only","spline_1","spline_2","spline_3","gcv_selected","source_copy")
            rr=filter(r->r.group==group && r.arm==arm,subsetrows); pixels=sum((r.test_pixels for r in rr);init=0)
            sse=sum((r.test_sse_nm2 for r in rr);init=0.)
            push!(summary,(;subset,group,arm,files=length(unique(r.file for r in rr)),views=length(rr),pixels,
                sse_nm2=sse,rms_pm=pixels>0 ? 1000sqrt(sse/pixels) : NaN))
        end
        for group in c["selection"]["groups"],reference in ("background_only","source_copy")
            picked=Dict((r.file,r.view)=>r for r in subsetrows if r.group==group && r.arm=="gcv_selected")
            refs=filter(r->r.group==group && r.arm==reference,subsetrows)
            deltas=[picked[(r.file,r.view)].test_sse_nm2-r.test_sse_nm2 for r in refs]
            @test all(picked[(r.file,r.view)].test_pixels==r.test_pixels for r in refs)
            push!(comparisons,(;subset,group,reference,views=length(refs),better=count(<(0),deltas),worse=count(>(0),deltas),
                tied=count(iszero,deltas),sse_difference_nm2=sum(deltas),
                mse_gain=-sum(deltas)/sum(r.test_sse_nm2 for r in refs)))
        end
        for interval in 1:3
            rr=filter(r->(subset=="all" || r.identified) && r.interval==interval && r.status=="ok",probes)
            energy=sum((r.input_sse_nm2 for r in rr);init=0.); left=sum((r.residual_sse_nm2 for r in rr);init=0.)
            fractions=[r.remaining_fraction for r in rr if isfinite(r.remaining_fraction)]
            push!(retention,(;subset,interval,probes=length(rr),input_sse_nm2=energy,residual_sse_nm2=left,
                remaining_fraction=energy>0 ? left/energy : NaN,min_fraction=isempty(fractions) ? NaN : minimum(fractions),
                median_fraction=isempty(fractions) ? NaN : median(fractions),max_fraction=isempty(fractions) ? NaN : maximum(fractions)))
        end
    end
    R.table(joinpath(out,"summary.tsv"),summary); R.table(joinpath(out,"comparisons.tsv"),comparisons)
    R.table(joinpath(out,"retention.tsv"),retention)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifySharedImageEnvelope.main()
