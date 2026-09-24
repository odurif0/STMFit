using Test, TOML, LinearAlgebra, Statistics, SHA
include(joinpath(@__DIR__,"measure_calibration.jl"))
include(joinpath(@__DIR__,"summarize_calibration_measurements.jl"))
const M=CalibrationMeasurements
const R=CalibrationMeasurementReport
const C=M.settings(M.SETTINGS)

function synthetic_image(path="synthetic.sxm"; backward=true, unit="nm", bad=false, flat=false)
    xs=collect(range(0,6;length=201)); ys=collect(range(0,2;length=101))
    z=[flat ? 0.0 : sum(exp(-((x-t)/.15)^2/2-((y-1)/.10)^2/2) for t in (1.,3.,5.)) for y in ys,x in xs]
    bad && (z[1,1]=NaN)
    channels=[M.SXMChannel("Z",unit,"fwd",z)]
    backward && push!(channels,M.SXMChannel("Z",unit,"bwd",copy(z)))
    M.SXMImage(path,Dict{String,String}(),length(xs),length(ys),(6.,2.),(0.,0.),channels)
end

@testset "Serialized measurement verification and paired report" begin
    mktempdir() do dir
        raw=joinpath(dir,"raw"); mkdir(raw)
        root=joinpath(dir,"results"); mkdir(root)
        for i in 1:3; write(joinpath(raw,"image$i.sxm"),"Synthetic report fixture $i"); end
        for shard in 1:2
            o=M.options(["--data-dir",raw,"--outdir",joinpath(root,"chunk$shard"),"--chunk","$shard/2"])
            M.execute(o;reader=path->synthetic_image(path))
        end
        r=R.inspect(root,raw,M.SETTINGS)
        @test r.scans==3 && length(r.rows)==12 && length(r.pairs)==6
        @test all(p->p["width_median_nm_relative_difference"]==0,r.pairs)
        @test all(p->p["axis_difference_deg"]==0,r.pairs)
        @test all(p->p["peak_count_equal"],r.pairs)
        R.main([root,raw,M.SETTINGS])
        @test occursin("not a count-accuracy grade",read(joinpath(root,"comparison.md"),String))
        @test length(R.read_table(joinpath(root,"paired_measurements.tsv")))==6
        @test_throws ErrorException R.main([root,raw,M.SETTINGS])
        # Output verification fails on a raw-byte or settings change.
        write(joinpath(raw,"image1.sxm"),"changed")
        @test_throws ErrorException R.inspect(root,raw,M.SETTINGS)
        write(joinpath(raw,"image1.sxm"),"Synthetic report fixture 1")
        open(joinpath(root,"chunk1/measurement_settings.toml"),"a") do io
            println(io,"# changed")
        end
        @test_throws ErrorException R.inspect(root,raw,M.SETTINGS)
    end
end

@testset "Bounded ordinary Slurm execution" begin
    path=joinpath(M.ROOT,"hpc/audit_calibration_measurements.sbatch")
    @test success(`bash -n $path`)
    source=read(path,String)
    for token in ("--time=01:00:00","--cpus-per-task=4","--mem=16000MB",
            "SLURM_JOB_ID","STMFIT_DATA_DIR","--dry-run","--chunk","wait")
        @test occursin(token,source)
    end
    for token in ("batch_full.jl","--expected-N","grade_unit_assignment","run_unknown_unit_assignment")
        @test !occursin(token,source)
    end
end
profile(t,z;px=minimum(diff(t)),n=ones(Int,length(t)))=(;t,z,n,px,binw=minimum(diff(t)))

@testset "Peak measurements, not assumed physical lobe widths" begin
    @test M.maxima(Float64[])==Int[]
    @test M.maxima([1.,1.,1.])==Int[]
    @test M.maxima([0.,1.,1.,0.])==[2]
    @test M.maxima([0.,1.,NaN,1.,0.])==Int[]
    @test M.maxima([1.,0.,1.])==Int[]
    t=collect(-2.:.002:2.); sigma=.2; y=exp.(-t.^2/(2sigma^2))
    r=M.peak_measurement(t,y,only(M.maxima(y)))
    @test isapprox(r.width_nm,2sqrt(2log(2))*sigma;atol=2e-5)
    @test isapprox(r.prominence_nm,1.;atol=1e-12)
    @test r.left_crossing_nm<r.t_nm<r.right_crossing_nm
    @test r.left_base_nm<=r.left_crossing_nm<r.right_crossing_nm<=r.right_base_nm
    for scale in (0.25,3.), offset in (-7.,4.)
        q=M.peak_measurement(t,scale.*y .+ offset,only(M.maxima(y)))
        @test q.width_nm≈r.width_nm atol=1e-12
        @test q.prominence_nm≈scale*r.prominence_nm atol=1e-12
    end
    # Independent published topographic-prominence example (one-based k=6).
    q=M.peak_measurement(collect(0.:8.),[0.,1.,0.,3.,1.,3.,0.,4.,0.],6)
    @test q.prominence_nm==3
    @test (q.left_base_nm,q.right_base_nm)==(2.,6.)
    @test q.width_nm==1.25
    # Overlap changes apparent half-prominence width: no universal sigma conversion.
    mixture=exp.(-((t.-.3)./.2).^2/2).+exp.(-((t.+.3)./.2).^2/2)
    @test length(M.maxima(mixture))==2
    for k in M.maxima(mixture)
        q=M.peak_measurement(t,mixture,k)
        @test abs(q.width_nm-2sqrt(2log(2))*.2)>.02
    end
    ob=M.observations(profile(t,y),C,0.)
    @test ob.n_widths==1 && ob.n_spacings==0 && isnan(ob.spacing_median_nm)
    empty=M.observations(profile(collect(0.:4.),zeros(5)),C,0.)
    @test empty.status=="unavailable" && isnan(empty.width_low_nm)
    # Peaks on opposite sides of a gap do not establish an observed spacing.
    gaps=M.observations(profile(collect(0.:8.),[0.,1.,0.,NaN,NaN,NaN,0.,1.,0.];px=.1),C,0.)
    @test gaps.n_widths==2 && gaps.n_spacings==0
    rejected=M.observations(profile(t,y),C,1.)
    @test rejected.n_widths==0 && all(p->p.reason=="low_prominence",rejected.peaks)
end

@testset "Weighted PCA, full profiles and explicit missing bins" begin
    xs=collect(range(-3,3;length=121)); ys=collect(range(-2,2;length=81))
    z=[exp(-((x*cos(.4)+y*sin(.4))/.9)^2/2-((-x*sin(.4)+y*cos(.4))/.15)^2/2) for y in ys,x in xs]
    f=M.axis_frame(xs,ys,z,C)
    @test abs(dot(f.axis,[cos(.4),sin(.4)]))>.9999
    inds=findall(f.mask); points=[[xs[i[2]],ys[i[1]]] for i in inds]
    weights=z[inds].-minimum(z[inds]); centre=sum(points.*weights)/sum(weights)
    cov=sum(w*((p-centre)*transpose(p-centre)) for (w,p) in zip(weights,points))/sum(weights)
    @test f.origin≈centre atol=1e-12
    @test abs(dot(f.axis,eigen(Symmetric(cov)).vectors[:,2]))≈1 atol=1e-12
    @test M.axis_frame(xs,ys,zeros(length(ys),length(xs)),C)===nothing
    transformed=M.axis_frame(xs,ys,3z.+2,C)
    @test abs(dot(transformed.axis,f.axis))≈1 atol=1e-12
    @test transformed.origin≈f.origin atol=1e-12
    p=M.strip_profile(xs,ys,z,f,C)
    old=M.strip_profile(xs,ys,z,f,C;legacy=true)
    @test length(p.t)>=length(old.t)
    @test sum(p.n)>=sum(old.n)
    # This elongated Gaussian's whole narrow strip can lie in its brightest
    # 30%; explicitly mask its tails to exercise the clipping distinction.
    clipped=merge(f,(;mask=f.mask .& (z .> maximum(z)/2)))
    clipped_profile=M.strip_profile(xs,ys,z,clipped,C;legacy=true)
    @test sum(p.n)>sum(clipped_profile.n)
    @test length(p.t)>length(clipped_profile.t)
    @test all(isfinite,p.z[p.n.>0])
    @test all(isnan,p.z[p.n.==0])
    # Artificial sampling gap is retained instead of being zero-filled.
    gap=copy(z); gap[:,55:68].=NaN
    pg=M.strip_profile(xs,ys,gap,f,C)
    @test count(==(0),pg.n)>0
    @test all(isnan,pg.z[pg.n.==0])
    @test issorted(pg.t) && all(>(0),diff(pg.t))
end

@testset "Legacy estimates and placeholders stay separate" begin
    p=profile(collect(0.:.02:2.),[exp(-((x-1)/.08)^2/2) for x in 0.:.02:2.];px=.01)
    ob=M.observations(p,C,0.;legacy=true)
    threshold=median(p.z)+2std(p.z)
    positions=[k for k in 2:length(p.z)-1 if p.z[k]>threshold && p.z[k]>=p.z[k-1] && p.z[k]>=p.z[k+1]]
    widths=Float64[]
    for k in positions
        win=max(1,k-3):min(length(p.z),k+3); amp=p.z[k]-minimum(p.z[win])
        half=minimum(p.z[win])+amp/2; above=win[p.z[win].>=half]
        width=p.t[min(last(above)+1,length(p.t))]-p.t[first(above)]
        amp>0 && .02<width<5 && push!(widths,width)
    end
    @test ob.n_peaks==length(positions)
    @test ob.width_low_nm==quantile(widths,.25)
    @test ob.width_high_nm==quantile(widths,.95)
    r=M.measure_image(synthetic_image(;flat=true),C)
    for row in r.rows
        @test row["status"]=="unavailable"
        @test isnan(row["width_low_nm"]) && isnan(row["spacing_median_nm"])
        @test !row["production_calibration_emitted"]
        if row["method"]=="legacy_bootstrap"
            @test row["legacy_width_fallback_would_apply"]
            @test row["legacy_spacing_fallback_would_apply"]
            @test row["legacy_reported_width_low_nm"]==.3
            @test row["legacy_reported_width_high_nm"]==1.
            @test row["legacy_reported_spacing_nm"]==.5
        else
            @test isnan(row["legacy_reported_width_low_nm"])
        end
    end
end

@testset "Shared I/O, two actual directions, and no imputed measurements" begin
    r=M.measure_image(synthetic_image(),C)
    @test length(r.rows)==4
    f=only(filter(x->x["method"]=="full_profile_prominence" && x["direction"]=="fwd",r.rows))
    b=only(filter(x->x["method"]=="full_profile_prominence" && x["direction"]=="bwd",r.rows))
    @test f["status"]=="apparent_only" && f["width_samples"]==3
    @test f["spacing_median_nm"]≈2 atol=.1
    for k in M.SUMMARY_HEADER
        k=="direction" && continue
        @test isequal(f[k],b[k])
    end
    cfg=M.GaussianFit2D.PatternConfig(filepath="synthetic.sxm",channel="Z",direction="fwd",
        stride=1,flatten="plane+rows",smooth_radius_px=1,no_plot=true)
    img=synthetic_image(); _,_,_,z,zs,_,disp=M.GaussianFit2D.preprocess_channel(img,first(img.channels),cfg)
    @test f["pipeline_dispersion_nm"]==disp
    @test f["hf_mad_nm"]==M.mad_scale(vec(z.-zs))
    @test f["pixel_x_nm"]==6/200 && f["pixel_y_nm"]==2/100
    partial=M.measure_image(synthetic_image(;bad=true),C)
    for row in partial.rows
        @test row["raw_nonfinite"]==1
        @test row["observed_pixels"]==201*101-1
        @test row["smoothed_observed_pixels"]==201*101-4
        @test !row["production_calibration_emitted"]
        @test row["method"]=="legacy_bootstrap" ? startswith(row["status"],"legacy_imputed_") : row["status"]=="apparent_only"
    end
    for (img,reason) in ((synthetic_image(;backward=false),"missing_or_ambiguous_direction"),
            (synthetic_image(;unit="V"),"invalid_units_or_preprocessed_grid"))
        r=M.measure_image(img,C)
        @test length(r.rows)==4
        selected=filter(x->x["direction"]=="bwd",r.rows)
        @test all(x->x["status"]=="unavailable" && x["reason"]==reason,selected)
        @test all(x->!x["production_calibration_emitted"],r.rows)
    end
end

@testset "Observed-only shared preprocessing" begin
    kwargs=(;stride=1,flatten="plane+rows",smooth_radius_px=1,plane_rank_rtol=1e-12)
    img=synthetic_image(); ch=first(img.channels)
    a=M.STMSXMIO.preprocess_observed_channel(img,ch;kwargs...)
    cfg=M.GaussianFit2D.PatternConfig(filepath=img.filepath,channel="Z",direction="fwd",
        stride=1,flatten="plane+rows",smooth_radius_px=1,no_plot=true)
    xs,ys,raw,z,zs,unit,_=M.GaussianFit2D.preprocess_channel(img,ch,cfg)
    @test a.xs==xs && a.ys==ys && a.raw==raw
    @test a.z==z && a.z_smooth==zs && a.unit==unit
    @test all(a.observed) && all(a.smoothed_observed)
    partial=synthetic_image(;bad=true)
    a=M.STMSXMIO.preprocess_observed_channel(partial,first(partial.channels);kwargs...)
    @test count(a.observed)==length(a.raw)-1
    @test count(a.smoothed_observed)==length(a.raw)-4
    @test isnan(a.z[1,1]) && all(isnan,a.z_smooth[1:2,1:2])
    @test all(isfinite,a.z[a.observed])
    @test all(isfinite,a.z_smooth[a.smoothed_observed])
    @test !any(a.smoothed_observed .& .!a.observed)
    # Entire unacquired rows remain missing through plane, row and box operations.
    ch.data[1:12,:].=NaN
    a=M.STMSXMIO.preprocess_observed_channel(img,ch;kwargs...)
    @test all(isnan,a.z[1:12,:]) && all(isnan,a.z_smooth[1:13,:])
    @test any(isfinite,a.z_smooth[14:end,:])
    ch.data.=NaN
    @test M.STMSXMIO.preprocess_observed_channel(img,ch;kwargs...).status=="no_observed_pixels"
    ch.data[20,:].=1
    @test M.STMSXMIO.preprocess_observed_channel(img,ch;kwargs...).status=="background_plane_unidentifiable"
    ch.data.=NaN
    for i in 1:min(size(ch.data)...); ch.data[i,i]=1; end
    @test M.STMSXMIO.preprocess_observed_channel(img,ch;kwargs...).status=="background_plane_unidentifiable"
    @test_throws ErrorException M.STMSXMIO.preprocess_observed_channel(img,ch;kwargs...,stride=0)
end

@testset "CLI, frozen settings, shards and diagnostic-only output" begin
    @test M.options(["--help"])===nothing
    for key in ("--n-lobe","--expected-N","--truth","--benchmark","--manifest","--selected-summary","--templates")
        @test_throws ErrorException M.options([key,"forbidden"])
    end
    mktempdir() do dir
        raw=joinpath(dir,"raw"); mkdir(raw)
        for i in 1:5; write(joinpath(raw,"image$i.sxm"),"Synthetic metadata-only fixture $i"); end
        args=["--data-dir",raw,"--outdir",joinpath(dir,"audit"),"--chunk","2/4"]
        o=M.options(vcat(args,["--dry-run"]))
        M.execute(o;reader=_ -> error("Dry-run must not read pixels"))
        @test !ispath(o["--outdir"])
        @test_throws ErrorException M.options(vcat(args,["--chunk","1/4"]))
        @test_throws ErrorException M.options(["--data-dir",raw,"--outdir",joinpath(dir,"bad"),"--chunk","1/6"])
        bad=deepcopy(C); bad["selection"]["expected_N"]=6
        path=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),path,"w")
        @test_throws ErrorException M.settings(path)
        delete!(o,"--dry-run")
        rows=M.execute(o;reader=path->synthetic_image(path))
        @test length(rows)==4 && all(r->r["file"]=="image2.sxm",rows)
        @test Set(readdir(o["--outdir"]))==Set(["measurements.tsv","profiles.tsv","peaks.tsv","raw_hashes.tsv","measurement_settings.toml"])
        @test read(joinpath(o["--outdir"],"measurement_settings.toml"))==read(M.SETTINGS)
        @test !occursin("N_selected",read(joinpath(o["--outdir"],"measurements.tsv"),String))
        @test !any(occursin("sigma_parallel",f) for f in readdir(o["--outdir"]))
        @test_throws ErrorException M.execute(o;reader=_ -> error("Must reject existing output"))
        hash=bytes2hex(sha256(read(joinpath(raw,"image2.sxm"))))
        @test occursin(hash,read(joinpath(o["--outdir"],"raw_hashes.tsv"),String))
        @test_throws ErrorException M.options(["--file",joinpath(raw,"image1.sxm"),"--data-dir",raw])
        @test isfile(joinpath(raw,"image1.sxm"))
    end
end
