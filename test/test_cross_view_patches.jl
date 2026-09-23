#!/usr/bin/env julia
# Includes legacy patch regressions, then opt-in subtraction-only SXM checks.
include(joinpath(@__DIR__,"test_patch_preprocessing.jl"))
const RF=ForwardPatches.ResidualFeatures

function write_rows(path,rows)
    header=sort(collect(keys(first(rows))))
    open(path,"w") do io
        println(io,join(header,'\t'))
        for r in rows; println(io,join([r[c] for c in header],'\t')); end
    end
    return path
end

function independent_patch(xs,ys,image,row)
    vals=Float64[]
    for u in -0.32:0.08:0.32,t in -0.32:0.08:0.32
        x=parse(Float64,row["x_nm"])+t; y=parse(Float64,row["y_nm"])+u
        ix=searchsortedlast(xs,x); iy=searchsortedlast(ys,y)
        tx=(x-xs[ix])/(xs[ix+1]-xs[ix]); ty=(y-ys[iy])/(ys[iy+1]-ys[iy])
        push!(vals,sum(image[iy+j,ix+i]*(i==0 ? 1-tx : tx)*(j==0 ? 1-ty : ty) for j in 0:1,i in 0:1))
    end
    return (vals.-median(vals))./std(vals)
end

@testset "Subtraction-only models preserve main features, raw patches and grids" begin
    mktempdir() do dir
        file="fixture.sxm"; synthetic_sxm(joinpath(dir,file))
        rows=[Dict("file"=>file,"lobe"=>string(i),"N"=>"2","source"=>"ell",
            "amplitude"=>"0.6","baseline"=>"0.1","tilt_x"=>"0.02","tilt_y"=>"-0.01","amp_rel"=>"1.0","gcv"=>"0.01",
            "x_nm"=>string(x),"y_nm"=>"1.2","t_nm"=>string(x),"u_nm"=>"1.2",
            "axis_x"=>"1.0","axis_y"=>"0.0","sigma_parallel_nm"=>"0.2","sigma_perp_nm"=>"0.13","skew_ratio"=>"1.0") for (i,x) in enumerate((1.0,1.5))]
        features=write_rows(joinpath(dir,"features.tsv"),rows)
        mf=[merge(r,Dict("profile_view"=>"fwd")) for r in rows]
        mb=[merge(r,Dict("profile_view"=>"bwd")) for r in rows]
        fwd=write_rows(joinpath(dir,"fwd.tsv"),mf); bwd=write_rows(joinpath(dir,"bwd.tsv"),mb)
        models=RF.read_residual_models(fwd,bwd,rows)
        @test models.fwd[file]==mf && models.bwd[file]==mb
        img=ForwardPatches.read_sxm(joinpath(dir,file))
        pcfg=ForwardPatches.PatternConfig(filepath=joinpath(dir,file),channel="Z",direction="fwd",stride=1,flatten="plane+rows",smooth_radius_px=1,no_plot=true)
        xs,ys,_,zf,_,_,_=ForwardPatches.preprocess_channel(img,ForwardPatches.get_channel(img,"Z";direction="fwd"),pcfg)
        _,_,_,zb,_,_,_=ForwardPatches.preprocess_channel(img,ForwardPatches.get_channel(img,"Z";direction="bwd"),pcfg)
        explicit(rs)=[parse(Float64,rs[1]["baseline"])+parse(Float64,rs[1]["tilt_x"])*x+parse(Float64,rs[1]["tilt_y"])*y+
            sum(parse(Float64,r["amplitude"])*exp(-0.5*((x-parse(Float64,r["x_nm"]))/0.2)^2-0.5*((y-1.2)/0.13)^2) for r in rs) for y in ys,x in xs]
        smooth(z)=[mean(z[max(1,y-1):min(end,y+1),max(1,x-1):min(end,x+1)]) for y in axes(z,1),x in axes(z,2)]
        common=["--features",features,"--data-dir",dir,"--config",COUNT_CONFIG,"--half-nm","0.32","--step-nm","0.08"]
        for (cfg,mode) in ((ASSIGNMENT_CONFIG,"legacy"),(MATCHED_CONFIG,"matched")), (name,mod) in (("fwd",ForwardPatches),("bwd",BackwardPatches))
            args=vcat(common,["--assignment-config",cfg])
            original=joinpath(dir,"$(mode)_$(name)_original.tsv")
            identity=joinpath(dir,"$(mode)_$(name)_identity.tsv")
            extract_quietly(mod,vcat(args,["--out",original]))
            extract_quietly(mod,vcat(args,["--residual-features-fwd",fwd,"--residual-features-bwd",bwd,"--out",identity]))
            @test read(original)==read(identity)
        end
        # Coefficients differ between views, but sampling and metadata still use rows.
        for (i,r) in enumerate(mf); r["amplitude"]=string(0.3+0.4i); r["tilt_x"]="0.06"; end
        for (i,r) in enumerate(mb); r["amplitude"]=string(1.0-0.3i); r["baseline"]="0.02"; end
        write_rows(fwd,mf); write_rows(bwd,mb)
        evaluated=RF.subtraction_model(mf,xs,ys,ForwardPatches._eval_peak)
        # Elementwise rounding bound: the oracle sums lobes before adding the
        # plane, whereas production accumulates them in place. Not a matrix norm.
        @test maximum(abs,evaluated-explicit(mf))<=4eps(maximum(abs,evaluated))
        for (arm,fp,bp,fr,br) in (("same",fwd,bwd,mf,mb),("cross",bwd,fwd,mb,mf)),
            (cfg,mode) in ((ASSIGNMENT_CONFIG,"legacy"),(MATCHED_CONFIG,"matched"))
            rf=mode=="matched" ? smooth(zf-explicit(fr)) : smooth(zf)-explicit(fr)
            rb=mode=="matched" ? smooth(zb-explicit(br)) : smooth(zb)-explicit(br)
            for (name,mod,channels) in (("fwd",ForwardPatches,[("raw",smooth(zf)),("res",rf)]),
                ("bwd",BackwardPatches,[("bwd_raw",smooth(zb)),("bwd_res",rb),("diff_raw",smooth(zf)-smooth(zb)),("diff_res",rf-rb)]))
                out=joinpath(dir,"$(mode)_$(arm)_$(name).tsv")
                extract_quietly(mod,vcat(common,["--assignment-config",cfg,"--residual-features-fwd",fp,"--residual-features-bwd",bp,"--out",out]))
                _,actual=mod._read_tsv(out)
                _,original=mod._read_tsv(joinpath(dir,"$(mode)_$(name)_original.tsv"))
                @test length(actual)==2
                for i in 1:2
                    for col in ("file","lobe","t_nm","u_nm","amplitude"); @test actual[i][col]==rows[i][col]; end
                    for col in keys(actual[i])
                        occursin("raw_p",col) && @test actual[i][col]==original[i][col]
                    end
                    for (prefix,im) in channels
                        expected=independent_patch(xs,ys,im,rows[i])
                        written=[parse(Float64,actual[i][prefix*"_p"*lpad(string(j),3,'0')]) for j in 1:81]
                        @test written≈expected atol=5e-6 rtol=5e-6
                    end
                end
            end
        end
        sentinel=joinpath(dir,"untouched.tsv"); write(sentinel,"untouched")
        for mod in (ForwardPatches,BackwardPatches)
            @test_throws ErrorException mod.main(vcat(common,["--residual-features-fwd",fwd,"--out",sentinel]))
            @test read(sentinel,String)=="untouched"
            @test_throws ErrorException mod._parse_cli(vcat(common,["--residual-features-fwd",fwd,"--residual-features-fwd",fwd]))
        end
        @test_throws ErrorException RF.read_residual_models(fwd,fwd,rows)
        @test_throws ErrorException RF.read_residual_models(fwd,bwd,rows;patch_frames="anything")
        @test_throws ErrorException RF.read_residual_models(fwd,bwd,rows;acquisition_shifts="anything")
        @test_throws ErrorException RF.read_residual_models(fwd,bwd,[merge(r,Dict("model_axis_x"=>"1")) for r in rows])
        for (field,value) in (("N","3"),("x_nm","1.000001"),("amplitude","NaN"),("profile_view","both"),("tilt_x","Inf"),("amplitude","-1"))
            bad=deepcopy(mf); bad[1][field]=value; write_rows(fwd,bad)
            @test_throws ErrorException RF.read_residual_models(fwd,bwd,rows)
        end
        write_rows(fwd,mf[1:1]); @test_throws ErrorException RF.read_residual_models(fwd,bwd,rows)
        write_rows(fwd,vcat(mf,mf[1:1])); @test_throws ErrorException RF.read_residual_models(fwd,bwd,rows)
    end
end
