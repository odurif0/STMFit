using Test, TOML, Random
include(joinpath(@__DIR__,"verify_image_registration_2d.jl"))
const V=VerifyImageRegistration2D
const D=V.D
const R=D.R
const ROOT=dirname(@__DIR__)
const CONFIG=joinpath(ROOT,"config","image_registration_2d.toml")
const S=R.settings(CONFIG)

@testset "Strict label-free image settings and fixed masks" begin
    @test S.residual_window_nm==.16 && S.refinement_step_px==.125
    mktempdir() do dir
        for (section,key,value) in (("selection","expected_N",6),("selection","sequence","NKNNKN"),
            ("model","refinement_step_px",.3),("preprocessing","holdout_buffer_px",32),
            ("preprocessing","residual_window_nm",.04))
            t=TOML.parsefile(CONFIG); t[section][key]=value; path=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,t),path,"w")
            @test_throws ErrorException R.settings(path)
        end
    end
    @test_throws ErrorException D.parse_cli(["--labels","forbidden"])
    rng=Xoshiro(131); a=randn(rng,17,23); b=randn(rng,17,23)
    a[6,13]=NaN; b[8,14]=NaN; base=(-3,1); limits=(2,3)
    for fold in 0:2
        allowed=R.training_mask(size(a),fold,16,1)
        mask=R.common_support(a,b,base,limits,2;allowed)
        brute=falses(size(a))
        for y in axes(a,1),x in axes(a,2)
            isfinite(a[y,x]) && allowed[y,x] || continue
            brute[y,x]=all(1<=y+base[2]+dy<=size(b,1) && 1<=x+base[1]+dx<=size(b,2) &&
                1<=y+dy<=size(b,1) && 1<=x+dx<=size(b,2) && allowed[y+dy,x+dx] &&
                isfinite(b[y+base[2]+dy,x+base[1]+dx]) for dy in -limits[2]:limits[2],dx in -limits[1]:limits[1])
        end
        for y in axes(a,1); count(@view brute[y,:])>=2 || (brute[y,:].=false); end
        @test mask==brute
    end
    @test R.sample([1. NaN; NaN NaN],1.,1.)==1
    @test R.sample([1. 3.; 5. 7.],1.25,1.5)==3
end

# A textured, nonperiodic analytic image; no molecular count or template.
texture(x,y)=sin(.17x+.063y)+.6cos(.31x-.12y)+.4sin(.071x+.29y)+
    .3exp(-((x-65)^2+(y-40)^2)/350)+.00008x*y
function fixture(dx,dy;shape=(96,128),base=(-4,0))
    a=[texture(x,y)+.01y for y in 1:shape[1],x in 1:shape[2]]
    b=[1.3texture(x-base[1]-dx,y-base[2]-dy)-.02y+3 for y in 1:shape[1],x in 1:shape[2]]
    return a,b,base
end

@testset "Subpixel sign, row offsets, independent score and conditional folds" begin
    a,b,base=fixture(.375,-.625)
    s=merge(S,(residual_window_nm=.06,peak_neighborhood_nm=.02,max_band_disagreement_nm=.02,min_band_pixels=64,min_band_rows=4,))
    runs=[R.scan(a,b,(.01,.01),base,s;fold=f) for f in 0:2]
    for r in runs
        p=first(r.peaks)
        @test abs(p.dx-.375)<=.125 && abs(p.dy+.625)<=.125
        @test p.correlation>.999
        @test all(q->abs(q.dx-.375)<=.25 && abs(q.dy+.625)<=.25,r.peaks)
        for region in 0:s.bands
            p=r.peaks[region+1]
            @test V.independent_score(a,b,r.mask,base,p.dx,p.dy,region,s.bands)≈p.correlation atol=1e-12
        end
    end
    # No heldout source pixel is read even when backward interpolation moves it.
    fold=1; allowed=R.training_mask(size(a),fold,s.holdout_block_px,s.holdout_buffer_px)
    aa=copy(a); bb=copy(b); aa[.!allowed].+=50
    for x in axes(b,2),y in axes(b,1)
        yy,xx=y-base[2],x-base[1]
        if !(1<=yy<=size(b,1) && 1<=xx<=size(b,2) && allowed[yy,xx]); bb[y,x]-=100; end
    end
    changed=R.scan(aa,bb,(.01,.01),base,s;fold)
    @test isequal(changed.scores,runs[2].scores) && isequal(changed.peaks,runs[2].peaks)
    @test R.eligibility(runs,s).spread_px<=.25
    @test_throws ErrorException R.eligibility(runs[1:2],s)
    # Preserve saved output format and independently re-read all three folds.
    mktempdir() do dir
        ny,nx=size(a); im=(;a,b,step=(.01,.01),xs=.01 .* (0:nx-1),ys=.01 .* (0:ny-1))
        settings=TOML.parsefile(CONFIG)
        for section in values(settings),key in keys(section); section[key]=getproperty(s,Symbol(key)); end
        path=joinpath(dir,"synthetic.toml"); open(io->TOML.print(io,settings),path,"w")
        out=joinpath(dir,"run"); mkpath(out)
        for r in runs
            dest=joinpath(out,"fold$(r.fold)"); D.save_result(dest,r,im,(;s,base),0.)
            cp(path,joinpath(dest,"settings.toml")); cp(joinpath(ROOT,"config","chitosan.toml"),joinpath(dest,"physical.toml"))
            D.write_records(joinpath(dest,"input_hashes.tsv"),[Dict("input"=>key,"sha256"=>bytes2hex(D.sha256(read(joinpath(dest,name)))))
                for (key,name) in (("--config","physical.toml"),("--settings","settings.toml"))])
        end
        D.merge_chunks(out); V.verify(out)
        @test_throws Exception D.merge_chunks(out)
    end
end

@testset "No signal, periodic ambiguity, and boundary are explicit failures" begin
    s=merge(S,(residual_window_nm=.06,peak_neighborhood_nm=.02,min_band_pixels=32,min_band_rows=2,min_row_pixels=4,bands=2))
    a,b,base=fixture(8.,0.;shape=(48,64),base=(0,0))
    r=R.scan(a,b,(.01,.01),base,s)
    @test !r.accepted && occursin("boundary_peak",r.status)
    constant=R.scan(ones(48,64),ones(48,64),(.01,.01),(0,0),s)
    @test !constant.accepted && occursin("constant_signal",constant.status)
    periodic=[sin(pi*x/2) for y in 1:48,x in 1:64]
    repeat_=R.scan(periodic,periodic,(.01,.01),(0,0),s)
    @test !repeat_.accepted && occursin("ambiguous_peak",repeat_.status)
    empty=R.scan(fill(NaN,48,64),ones(48,64),(.01,.01),(0,0),s)
    @test !empty.accepted && occursin("insufficient_support",empty.status)
    @test success(`bash -n $(joinpath(ROOT,"hpc","estimate_image_registration_2d.sbatch"))`)
end
