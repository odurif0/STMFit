#!/usr/bin/env julia
# Re-read saved pixels and independently check masks and the reported image peaks.
module VerifyImageRegistration2D
using Test, DelimitedFiles, Statistics, SHA
include(joinpath(@__DIR__,"estimate_image_registration_2d.jl"))
const D=EstimateImageRegistration2D
val(r,k)=parse(Float64,r[k])

function read_images(out)
    m=readdlm(joinpath(out,"fold0","images.tsv"),'\t',Float64;skipstart=1)
    ny,nx=Int(maximum(m[:,1])),Int(maximum(m[:,2]))
    size(m,1)==ny*nx || error("Incomplete saved image grid")
    a=fill(NaN,ny,nx); b=copy(a); seen=falses(ny,nx)
    for i in axes(m,1)
        y,x=Int(m[i,1]),Int(m[i,2]); seen[y,x] && error("Duplicate image pixel")
        seen[y,x]=true; a[y,x]=m[i,5]; b[y,x]=m[i,6]
    end
    return (;a,b,grid=m)
end

function independent_score(a,b,mask,base,dx,dy,region,bands)
    aa=0.; bb=0.; ab=0.; ny=size(a,1)
    for y in axes(a,1)
        region==0 || min(bands,1+(y-1)*bands÷ny)==region || continue
        cols=findall(@view mask[y,:]); isempty(cols) && continue
        av=a[y,cols]; bv=Float64[]
        for x in cols
            yy,xx=y+base[2]+dy,x+base[1]+dx
            yl,xl=floor(Int,yy),floor(Int,xx); yh,xh=ceil(Int,yy),ceil(Int,xx)
            fy,fx=yy-yl,xx-xl
            # At integral positions floor==ceil; all four references remain observed.
            push!(bv,(1-fy)*((1-fx)*b[yl,xl]+fx*b[yl,xh])+fy*((1-fx)*b[yh,xl]+fx*b[yh,xh]))
        end
        av=av.-mean(av); bv=bv.-mean(bv)
        aa+=sum(abs2,av); bb+=sum(abs2,bv); ab+=sum(av.*bv)
    end
    return aa>0 && bb>0 ? clamp(ab/sqrt(aa*bb),-1.,1.) : NaN
end

function verify(out)
    s=D.R.settings(joinpath(out,"settings.toml")); im=read_images(out); a,b=im.a,im.b
    summaries=D.table(joinpath(out,"summary.tsv")); allpeaks=D.table(joinpath(out,"peaks.tsv"))
    @testset "Image-only input identity and complete full/fold outputs" begin
        @test parse.(Int,getindex.(summaries,"fold"))==[0,1,2]
        @test length(allpeaks)==3(s.bands+1)
        for fold in 0:2
            dir=joinpath(out,"fold$fold")
            for name in ("input_hashes.tsv","physical.toml","settings.toml")
                @test read(joinpath(dir,name))==read(joinpath(out,name))
            end
            @test only(D.table(joinpath(dir,"summary.tsv")))==summaries[fold+1]
        end
        for (key,name) in (("--config","physical.toml"),("--settings","settings.toml"))
            row=only(r for r in D.table(joinpath(out,"input_hashes.tsv")) if r["input"]==key)
            @test row["sha256"]==bytes2hex(sha256(read(joinpath(out,name))))
        end
    end
    @testset "Fixed observed support, no heldout source reuse, and direct peak arithmetic" begin
        for fold in 0:2
            dir=joinpath(out,"fold$fold"); r=summaries[fold+1]
            base=(parse(Int,r["base_dx"]),parse(Int,r["base_dy"]))
            step=(val(r,"step_x_nm"),val(r,"step_y_nm")); ny,nx=size(a)
            limits=(parse(Int,r["limit_x_px"]),parse(Int,r["limit_y_px"]))
            @test limits==Tuple(min(floor(Int,s.residual_window_nm/step[j]),floor(Int,s.max_lag_width_fraction*size(a,3-j))) for j in 1:2)
            @test parse(Int,r["nrows"])==ny && parse(Int,r["ncols"])==nx
            @test all(isapprox(im.grid[i,3]-im.grid[1,3],(im.grid[i,2]-1)*step[1];atol=1e-10) &&
                isapprox(im.grid[i,4]-im.grid[1,4],(im.grid[i,1]-1)*step[2];atol=1e-10) for i in axes(im.grid,1))
            allowed(y,x)=1<=y<=ny && 1<=x<=nx && (fold==0 ||
                (mod(div(y-1,s.holdout_block_px)+div(x-1,s.holdout_block_px),2)!=fold-1 &&
                 s.holdout_buffer_px<=mod(y-1,s.holdout_block_px)<s.holdout_block_px-s.holdout_buffer_px &&
                 s.holdout_buffer_px<=mod(x-1,s.holdout_block_px)<s.holdout_block_px-s.holdout_buffer_px))
            expected=falses(size(a))
            for x in 1:nx,y in 1:ny
                isfinite(a[y,x]) && allowed(y,x) || continue
                ok=true
                for sy in -limits[2]:limits[2],sx in -limits[1]:limits[1]
                    yy,xx=y+base[2]+sy,x+base[1]+sx
                    if !(1<=yy<=ny && 1<=xx<=nx && isfinite(b[yy,xx]) && allowed(y+sy,x+sx))
                        ok=false; break
                    end
                end
                expected[y,x]=ok
            end
            for y in 1:ny
                count(@view expected[y,:])>=s.min_row_pixels || (expected[y,:].=false)
            end
            mask=falses(size(a)); pairs=readdlm(joinpath(dir,"support.tsv"),'\t',Int;skipstart=1)
            for i in axes(pairs,1)
                y,x=pairs[i,1],pairs[i,2]; @test !mask[y,x]; mask[y,x]=true
            end
            @test mask==expected && count(mask)==parse(Int,r["support_pixels"])
            peaks=D.table(joinpath(dir,"peaks.tsv")); scores=D.table(joinpath(dir,"scores.tsv"))
            @test peaks==filter(p->parse(Int,p["fold"])==fold,allpeaks)
            for region in 0:s.bands
                p=only(p for p in peaks if parse(Int,p["region"])==region)
                coarse=filter(c->c["phase"]=="coarse" && parse(Int,c["region"])==region,scores)
                fine=filter(c->c["phase"]=="fine" && parse(Int,c["region"])==region,scores)
                @test [(val(c,"dx"),val(c,"dy")) for c in coarse]==[(Float64(x),Float64(y)) for y in -limits[2]:limits[2] for x in -limits[1]:limits[1]]
                order(c)=(-val(c,"correlation"),val(c,"dx")^2+val(c,"dy")^2,val(c,"dx"),val(c,"dy"))
                good=filter(c->isfinite(val(c,"correlation")),coarse)
                center=isempty(good) ? (0.,0.) : let c=first(sort(good;by=order)); (val(c,"dx"),val(c,"dy")) end
                ranges=[max(-limits[j],center[j]-s.refinement_radius_px):s.refinement_step_px:min(limits[j],center[j]+s.refinement_radius_px) for j in 1:2]
                @test [(val(c,"dx"),val(c,"dy")) for c in fine]==[(x,y) for y in ranges[2] for x in ranges[1]]
                goodfine=filter(c->isfinite(val(c,"correlation")),fine)
                if !isempty(goodfine)
                    winner=first(sort(goodfine;by=order))
                    @test (val(winner,"dx"),val(winner,"dy"))==(val(p,"dx"),val(p,"dy"))
                    @test val(winner,"correlation")==val(p,"correlation")
                end
                zero=only(c for c in coarse if val(c,"dx")==val(c,"dy")==0)
                for (dx,dy,reported) in ((val(p,"dx"),val(p,"dy"),val(p,"correlation")),(0.,0.,val(zero,"correlation")),
                    (val(p,"distant_dx"),val(p,"distant_dy"),val(p,"correlation")-val(p,"distant_peak_gap")))
                    isfinite(dx) && isfinite(dy) || continue
                    score=independent_score(a,b,mask,base,dx,dy,region,s.bands)
                    @test isnan(score) ? isnan(reported) : isapprox(score,reported;atol=3e-12,rtol=3e-11)
                end
                far=filter(c->isfinite(val(c,"correlation")) && max(abs(val(c,"dx")-val(p,"dx"))*step[1],abs(val(c,"dy")-val(p,"dy"))*step[2])>s.peak_neighborhood_nm,coarse)
                gap=isempty(far) ? NaN : val(p,"correlation")-maximum(val.(far,"correlation"))
                @test isequal(gap,val(p,"distant_peak_gap"))
                relevant=[y for y in 1:ny if region==0 || min(s.bands,1+(y-1)*s.bands÷ny)==region]
                np=sum(count(@view mask[y,:]) for y in relevant); nr=count(y->any(@view mask[y,:]),relevant)
                @test np==parse(Int,p["pixels"]) && nr==parse(Int,p["rows"])
                boundary=abs(val(p,"dx"))>=limits[1] || abs(val(p,"dy"))>=limits[2]
                usable=np>=s.min_band_pixels && nr>=s.min_band_rows && val(p,"correlation")>=s.min_correlation && gap>=s.min_distant_peak_gap && !boundary
                @test (p["usable"]=="true")==usable && (p["boundary"]=="true")==boundary
            end
            globalpeak=only(p for p in peaks if p["region"]=="0")
            accepted=all(p["usable"]=="true" && max(abs(val(p,"dx")-val(globalpeak,"dx"))*step[1],abs(val(p,"dy")-val(globalpeak,"dy"))*step[2])<=s.max_band_disagreement_nm for p in peaks) &&
                all(limits[j]*step[j]>s.peak_neighborhood_nm+step[j] for j in 1:2)
            @test (r["accepted"]=="true")==accepted
            @test val(r,"dx")==val(globalpeak,"dx") && val(r,"dy")==val(globalpeak,"dy")
        end
        spread=maximum(max(abs(val(a,"dx")-val(b,"dx")),abs(val(a,"dy")-val(b,"dy"))) for a in summaries,b in summaries)
        e=only(D.table(joinpath(out,"eligibility.tsv")))
        @test val(e,"spread_px")==spread
        @test (e["eligible"]=="true")== (all(r["accepted"]=="true" for r in summaries) && spread<=s.max_fold_disagreement_px)
    end
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==1 || error("Usage: verify_image_registration_2d.jl RUN_DIR")
    VerifyImageRegistration2D.verify(only(ARGS))
end
