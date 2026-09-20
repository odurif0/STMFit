using Test, TOML, Statistics, LinearAlgebra, Random
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
const ROOT = dirname(@__DIR__)
const BASE_CONFIG = joinpath(ROOT,"config","unit_assignment_matched_residual.toml")
const VARIANTS = (("transverse_moment","transverse_first_moment"),
                  ("affine_residual","affine_residual_half_plane_asymmetry"))
const COORDS = collect(-4:4) .* .08
const U = repeat(COORDS; inner=9)
const T = repeat(COORDS; outer=9)
const ZERO_L1 = 1e-12
descriptor(p,kind) = transverse_descriptor(p,COORDS,kind;zero_l1=ZERO_L1)
project(p) = ReconstructedUnitAssignment._affine_residual_patch(p,COORDS)

@testset "One explicit descriptor change, no numerical tuning" begin
    base = load_config(BASE_CONFIG)
    for (name,kind) in VARIANTS
        config = load_config(joinpath(ROOT,"config","unit_assignment_"*name*".toml"))
        @test config["model"]["descriptor"] == kind
        @test config["model"]["name"] == "cc_soft_"*name*"_v1"
        config["model"]["descriptor"] = base["model"]["descriptor"]
        config["model"]["name"] = base["model"]["name"]
        @test config == base
    end
    mktempdir() do dir
        bad = deepcopy(base); bad["model"]["descriptor"] = "guess"
        path = joinpath(dir,"bad.toml")
        open(io->TOML.print(io,bad),path,"w")
        @test_throws ErrorException load_config(path)
    end
end

@testset "Position-weighted first moment and unchanged legacy formula" begin
    near = zeros(81); near[50] = 1.0 # u=+0.08, t=0
    far = zeros(81); far[77] = 1.0   # u=+0.32, t=0
    @test descriptor(near,"transverse_first_moment") == (.25,"ok")
    @test descriptor(far,"transverse_first_moment") == (1.0,"ok")
    @test descriptor(near,"transverse_half_plane_asymmetry") == (1.0,"ok")
    @test descriptor(far,"transverse_half_plane_asymmetry") == (1.0,"ok")
    rng=MersenneTwister(20260920)
    for _ in 1:20
        p=randn(rng,81)
        @test descriptor(p,"transverse_half_plane_asymmetry") == transverse_asymmetry(p,COORDS;zero_l1=ZERO_L1)
        expected=sum((U/.32).*p)/sum(abs,p)
        @test first(descriptor(p,"transverse_first_moment")) ≈ expected atol=2e-16
    end
end

@testset "Affine projection: QR reference, removed signal and invariances" begin
    design=hcat(ones(81),T,U)
    shoulder=exp.(-0.5 .* (((U .- 0.12) ./ 0.055).^2 .+ ((T .- 0.04) ./ 0.07).^2))
    original=copy(shoulder)
    # Independently solve the full least-squares problem, not the separable formula.
    expected=shoulder-design*(design\shoulder)
    residual=project(shoulder)
    @test residual ≈ expected atol=5e-16
    @test project(residual) ≈ residual atol=5e-16
    @test maximum(abs,design'*residual) < 3e-14
    @test 0 < norm(residual) < norm(shoulder) # affine signal is removed, not preserved
    @test norm(shoulder)^2 ≈ norm(residual)^2+norm(shoulder-residual)^2 atol=3e-14
    expected_descriptor=sum(sign.(U).*expected)/sum(abs,expected)
    @test first(descriptor(shoulder,"affine_residual_half_plane_asymmetry")) ≈ expected_descriptor atol=2e-15
    @test abs(expected_descriptor) > .01 # retained localized odd structure
    for coefficients in ((.3,1.2,-.7),(-5.,-2.,9.),(0.,0.,3.))
        plane=design*collect(coefficients)
        @test maximum(abs,project(plane)) < ZERO_L1/81
        @test last(descriptor(plane,"affine_residual_half_plane_asymmetry")) == "zero_affine_residual_mass"
        @test maximum(abs,project(shoulder+plane)-residual) < 5e-15 # per-pixel, not aggregate L2 error
        @test first(descriptor(shoulder+plane,"affine_residual_half_plane_asymmetry")) ≈ expected_descriptor atol=2e-14
    end
    @test first(descriptor(shoulder+3U,"transverse_half_plane_asymmetry")) !=
          first(descriptor(shoulder,"transverse_half_plane_asymmetry"))
    for (_,kind) in VARIANTS
        value=first(descriptor(shoulder,kind))
        @test abs(value) <= 1
        @test first(descriptor(7shoulder,kind)) ≈ value atol=2e-15
        @test first(descriptor(-shoulder,kind)) ≈ -value atol=2e-15
        # Julia reshape is t-first, u-second for the actual serialization.
        @test first(descriptor(vec(reverse(reshape(shoulder,9,9);dims=2)),kind)) ≈ -value atol=2e-15
        @test first(descriptor(vec(reverse(reshape(shoulder,9,9);dims=1)),kind)) ≈ value atol=2e-15
        @test last(descriptor(zeros(81),kind)) == "zero_patch_mass"
        @test last(descriptor(fill(1e-16,81),kind)) == "zero_patch_mass"
        for bad in (NaN,Inf,-Inf)
            p=copy(shoulder); p[7]=bad
            @test last(descriptor(p,kind)) == "nonfinite_patch"
        end
        @test_throws ErrorException transverse_descriptor(ones(80),COORDS,kind;zero_l1=ZERO_L1)
        @test_throws ErrorException transverse_descriptor(shoulder,reverse(COORDS),kind;zero_l1=ZERO_L1)
        @test_throws ErrorException transverse_descriptor(shoulder,COORDS,kind;zero_l1=-1)
        @test_throws ErrorException transverse_descriptor(ones(9),zeros(3),kind;zero_l1=ZERO_L1)
    end
    @test shoulder == original
    @test_throws ErrorException descriptor(shoulder,"unregistered")
end

@testset "All-key tables, missing reasons and actual descriptor CLI" begin
    mktempdir() do dir
        features=joinpath(dir,"features.tsv"); patches=joinpath(dir,"patches.tsv")
        header=["bwd_res_p"*lpad(string(i),3,'0') for i in 1:81]
        shape=exp.(-0.5 .* (((U .- 0.12) ./ 0.055).^2 .+ ((T .- 0.04) ./ 0.07).^2))
        inputs=[shape, 1.0 .+ U .- 2T, zeros(81), fill(NaN,81)]
        base=[Dict("file"=>"synthetic.sxm","lobe"=>string(i),"amplitude"=>"0.1") for i in 1:4]
        write_table(features,["file","lobe","amplitude"],base)
        patchrows=[merge(base[i],Dict(c=>string(v) for (c,v) in zip(header,inputs[i]))) for i in 1:4]
        write_table(patches,vcat(["file","lobe"],header),patchrows)
        for (name,kind) in VARIANTS
            config=joinpath(ROOT,"config","unit_assignment_"*name*".toml")
            out=joinpath(dir,name*".tsv"); cliout=joinpath(dir,name*"_cli.tsv")
            augment_descriptor(features,patches,out,config)
            _,rows=lobe_table(out)
            @test Set(keys(rows)) == Set(("synthetic.sxm",i) for i in 1:4)
            for i in 1:4
                value,reason=descriptor(inputs[i],kind)
                @test rows[("synthetic.sxm",i)]["descriptor_reason"] == reason
                stored=rows[("synthetic.sxm",i)]["patch_u_asym_reconstructed"]
                @test isfinite(value) ? parse(Float64,stored) == value : stored == "NA"
            end
            cmd=`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_reconstructed_descriptor.jl")) --features $features --patches $patches --config $config --out $cliout`
            @test success(cmd)
            @test read(out) == read(cliout)
            @test_throws ErrorException augment_descriptor(features,patches,out,config)
        end
    end
end
