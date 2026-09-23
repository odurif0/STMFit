#!/usr/bin/env julia
# Saved-only checks, with independent convolution, centering and residual sums.
module VerifyDirectionalResponse
using Test, DelimitedFiles, Statistics, LinearAlgebra, SHA
include(joinpath(@__DIR__,"diagnose_directional_response.jl"))
const D=DiagnoseDirectionalResponse
number(r,k)=parse(Float64,r[k])

function images(dir)
    a=readdlm(joinpath(dir,"images.tsv"),'\t',Float64;skipstart=1)
    ny,nx=Int(maximum(a[:,1])),Int(maximum(a[:,2]))
    @test size(a)==(ny*nx,7)
    @test a[:,1]==repeat(collect(1:ny),nx)
    @test a[:,2]==repeat(collect(1:nx);inner=ny)
    return reshape(a[:,5],ny,nx),reshape(a[:,6],ny,nx),reshape(a[:,7],ny,nx)
end

"Independent matrix-vector convolution; no production kernel/transfer calls."
function response(z,mu,sign,epsilon)
    mu==0 && return copy(z)
    a=mu/(1+mu); last=0
    while (1+a)*a^(last+1)>epsilon; last+=1; end
    offs=sign.*([-1;collect(0:last)].-Int(2mu))
    weights=[-a;[(1-a*a)*a^k for k in 0:last]]
    out=fill(NaN,size(z))
    for x in max(1,1-minimum(offs)):min(size(z,2),size(z,2)-maximum(offs))
        out[:,x]=z[:,x.+offs]*weights
    end
    return out
end

function centered_pairs(a,b,pa,pb,mask,fold,dx,dy)
    interpolate(z,y,x)=begin
        yl,yh=floor(Int,y),ceil(Int,y); xl,xh=floor(Int,x),ceil(Int,x)
        fy,fx=y-yl,x-xl
        (1-fy)*((1-fx)*z[yl,xl]+fx*z[yl,xh])+fy*((1-fx)*z[yh,xl]+fx*z[yh,xh])
    end
    out=[(Float64[],Float64[],Float64[],Float64[]) for _ in 1:2]
    for y in axes(a,1)
        cols=findall(view(mask,y,:)); isempty(cols) && continue
        vectors=(collect(a[y,cols]),collect(b[y,cols]),
            [interpolate(pa,y+dy,x+dx) for x in cols],[interpolate(pb,y-dy,x-dx) for x in cols])
        for k in 1:4
            append!(out[fold[y]][k],vectors[k].-mean(vectors[k]))
        end
    end
    return out
end

function independent_scores(v,fold,s)
    tr=fold==0 ? ntuple(k->vcat(v[1][k],v[2][k]),4) : v[fold]
    te=fold==0 ? tr : v[3-fold]
    gf=clamp(dot(tr[1],tr[3])/sum(abs2,tr[3]),s.gain_min,s.gain_max)
    gb=clamp(dot(tr[2],tr[4])/sum(abs2,tr[4]),s.gain_min,s.gain_max)
    tf=sum(abs2,tr[1].-gf.*tr[3]); tb=sum(abs2,tr[2].-gb.*tr[4])
    ef=sum(abs2,te[1].-gf.*te[3]); eb=sum(abs2,te[2].-gb.*te[4])
    return (;gain_fwd=gf,gain_bwd=gb,train_fwd=tf,train_bwd=tb,test_fwd=ef,test_bwd=eb,train_rss=tf+tb,test_rss=ef+eb)
end

"Separable observed-run erosion, independent of production's integral image."
function expected_mask(a,b,rx,ry,s)
    ny,nx=size(a); h=falses(ny,nx); mask=copy(h); fold=zeros(Int,ny)
    for y in 1:ny
        start=1
        for x in 1:nx+1
            if x>nx || !isfinite(a[y,x]) || !isfinite(b[y,x])
                start+rx<=x-1-rx && (h[y,start+rx:x-1-rx].=true)
                start=x+1
            end
        end
    end
    for x in 1:nx
        start=1
        for y in 1:ny+1
            if y>ny || !h[y,x]
                start+ry<=y-1-ry && (mask[start+ry:y-1-ry,x].=true)
                start=y+1
            end
        end
    end
    for y in 1:ny
        offset=mod(y-1,s.holdout_band_rows)
        if !(s.holdout_buffer_rows<=offset<s.holdout_band_rows-s.holdout_buffer_rows) || count(view(mask,y,:))<s.min_pixels_per_row
            mask[y,:].=false
        else
            fold[y]=1+mod(fld(y-1,s.holdout_band_rows),2)
        end
    end
    return mask,fold
end

function verify(out)
    BLAS.set_num_threads(1)
    s=D.R.settings(joinpath(out,"settings.toml")); panel=D.table(joinpath(out,"panel.tsv"))
    base=D.table(joinpath(out,"base_shifts.tsv")); summaries=D.table(joinpath(out,"summary.tsv"))
    @testset "Response inputs and panel identity" begin
        expected=D.record.(D.R.panel(getindex.(base,"file"),s))
        @test panel==expected
        @test getindex.(panel,"file")==getindex.(summaries,"file")
        hashes=D.table(joinpath(out,"input_hashes.tsv"))
        for (key,path) in (("--base-shifts","base_shifts.tsv"),("--config","physical.toml"),("--settings","settings.toml"))
            @test only(r for r in hashes if r["input"]==key)["sha256"]==bytes2hex(sha256(read(joinpath(out,path))))
        end
    end
    supported=0; complete=true
    for p in panel
        file=p["file"]; dir=joinpath(out,splitext(file)[1])
        @testset "Saved directional result $file" begin
            if isfile(joinpath(dir,"failure.tsv"))
                @test only(D.table(joinpath(dir,"failure.tsv")))["file"]==file
                @test only(r for r in summaries if r["file"]==file)["status"]=="failed"
                complete=false
                continue
            end
            a,raw,b=images(dir); meta=only(D.table(joinpath(dir,"metadata.tsv")))
            dx=parse(Int,only(r for r in base if r["file"]==file)["bwd_sample_dx_px"])
            @test meta["base_dx"]==string(dx)
            @test number(meta,"line_time_fwd_s")==number(meta,"line_time_bwd_s")>0
            for x in axes(b,2)
                source=x+dx
                @test 1<=source<=size(b,2) ? isequal(b[:,x],raw[:,source]) : all(isnan,b[:,x])
            end
            summary=only(D.table(joinpath(dir,"summary.tsv")))
            rx,ry=parse(Int,summary["radius_x"]),parse(Int,summary["radius_y"])
            extent=0
            for mu in s.mean_lag_px
                mu==0 && continue
                coefficient=mu/(1+mu); last=0
                while (1+coefficient)*coefficient^(last+1)>s.kernel_tail_l1; last+=1; end
                extent=max(extent,Int(2mu)+1,abs(last-Int(2mu)))
            end
            @test rx==extent+ceil(Int,s.residual_shift_max_px)+1
            @test ry==ceil(Int,s.residual_shift_max_px)+1
            expected,fold=expected_mask(a,b,rx,ry,s)
            pixels=readdlm(joinpath(dir,"support.tsv"),'\t',Int;skipstart=1); mask=falses(size(a))
            @test length(Set(Tuple(pixels[i,:]) for i in axes(pixels,1)))==size(pixels,1)
            for i in axes(pixels,1)
                y,x,f=pixels[i,:]; mask[y,x]=true
                @test f==fold[y] && all(1+mod(fld(yy-1,s.holdout_band_rows),2)==f for yy in y-ry:y+ry)
            end
            @test mask==expected && count(mask)==parse(Int,summary["support_pixels"])
            for f in 1:2
                rows=findall(==(f),fold)
                @test length(rows)>=s.min_rows && sum(mask[rows,:])>=s.min_pixels
            end
            scores=D.table(joinpath(dir,"scores.tsv")); chosen=D.table(joinpath(dir,"selected.tsv"))
            grid=collect(-s.residual_shift_max_px:s.residual_shift_step_px:s.residual_shift_max_px)
            keys=Set((number(r,"lag"),parse(Int,r["sign"]),number(r,"dx"),number(r,"dy"),parse(Int,r["fold"])) for r in scores)
            expected_keys=Set((mu,sign,x,y,fold) for mu in s.mean_lag_px for sign in (mu==0 ? (1,) : (1,-1)) for x in grid for y in grid for fold in 0:2)
            @test keys==expected_keys && length(keys)==length(scores)
            @test length(chosen)==9
            for r in chosen
                fold_id=parse(Int,r["fold"]); mode=r["mode"]; lag=number(r,"lag"); sign=parse(Int,r["sign"])
                pool=filter(x->x["fold"]==r["fold"] && (number(x,"lag")==0 || (mode=="physical" && x["sign"]=="1") || (mode=="reversed" && x["sign"]=="-1")),scores)
                mode=="translation" && filter!(x->number(x,"lag")==0,pool)
                sort!(pool;by=x->(number(x,"train_rss"),number(x,"lag"),abs(number(x,"dx"))+abs(number(x,"dy")),number(x,"dx"),number(x,"dy")))
                @test all(r[k]==v for (k,v) in first(pool))
                pa=response(b,lag,-sign,s.kernel_tail_l1); pb=response(a,lag,sign,s.kernel_tail_l1)
                v=centered_pairs(a,b,pa,pb,mask,fold,number(r,"dx"),number(r,"dy"))
                for (key,value) in pairs(independent_scores(v,fold_id,s))
                    @test isapprox(value,number(r,string(key));atol=1e-10,rtol=1e-9)
                end
            end
            checks=D.table(joinpath(dir,"heldout.tsv"))
            predictive=true; directional=true
            for f in 1:2
                c=only(r for r in chosen if r["fold"]==string(f) && r["mode"]=="translation")
                p=only(r for r in chosen if r["fold"]==string(f) && r["mode"]=="physical")
                n=only(r for r in chosen if r["fold"]==string(f) && r["mode"]=="reversed")
                gain(k)=(number(c,k)-number(p,k))/number(c,k)
                reverse=(number(n,"test_rss")-number(p,"test_rss"))/number(n,"test_rss")
                check=only(r for r in checks if r["fold"]==string(f))
                @test isapprox(gain("test_fwd"),number(check,"gain_fwd");atol=1e-12)
                @test isapprox(gain("test_bwd"),number(check,"gain_bwd");atol=1e-12)
                @test isapprox(reverse,number(check,"reverse_advantage");atol=1e-12)
                predictive &= min(gain("test_fwd"),gain("test_bwd"))>=s.min_heldout_gain_fraction
                directional &= reverse>=s.min_reverse_advantage_fraction
            end
            physical=filter(r->r["mode"]=="physical",chosen)
            interior=all(r->0<number(r,"lag")<last(s.mean_lag_px) && max(abs(number(r,"dx")),abs(number(r,"dy")))<s.residual_shift_max_px,physical)
            gain_interior=all(r->s.gain_min<number(r,"gain_fwd")<s.gain_max && s.gain_min<number(r,"gain_bwd")<s.gain_max,physical)
            spread(k)=maximum(number.(physical,k))-minimum(number.(physical,k))
            lag_spread=spread("lag"); shift_spread=max(spread("dx"),spread("dy"))
            ok=interior && gain_interior && predictive && directional && lag_spread<=s.lag_agreement_px && shift_spread<=s.shift_agreement_px
            for (key,value) in pairs((;interior,gain_interior,predictive,directional,lag_spread,shift_spread,supported=ok))
                @test summary[string(key)]==string(value)
            end
            @test only(r for r in summaries if r["file"]==file)["supported"]==string(ok)
            supported+=ok
        end
    end
    @testset "Complete panel conclusion, no recognition claim" begin
        r=only(D.table(joinpath(out,"panel_result.tsv")))
        @test r["complete"]==string(complete)
        @test r["supporting"]==string(supported)
        @test r["total"]==string(length(panel))
        @test r["directional_hypothesis_supported"]==string(complete && supported>=s.min_supporting_scans)
        @test r["recognition_grade"]=="not_run_image_only_no_candidate_predictions"
    end
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifyDirectionalResponse.verify(only(ARGS))
