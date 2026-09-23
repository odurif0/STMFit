using Test, LinearAlgebra
include(joinpath(@__DIR__,"audit_qe_surface_calibration.jl"))
const A = AuditQESurfaceCalibration
const C = A.C

@testset "Diagnostic surface audit: analytic field and unchanged legacy grid" begin
    dims=(5,5,8)
    # First-axis-fast in-memory values, after QE input adaptation.
    vals=vec([0.3+0.1i+0.02j-0.3k for i in 0:4,j in 0:4,k in 0:7])
    c=C.Cube(zeros(3),Matrix(Diagonal([0.1,0.1,0.1])),dims,vals)
    f=C.Frame([0.2,0.2,0.0],[1.,0.,0.],[0.,1.,0.])
    s=C.MoldSettings(.1,.1,.2,0.,1.,.1,-3.,.8,8,"first_below_target")
    r=A.audit(c,f,s)
    @test r.summary["target_valid"]==9
    @test r.summary["target_negative"]==6
    @test r.summary["target_median"]≈-.06
    @test r.summary["native_minimum"]==minimum(vals)
    @test r.summary["total_native_negative"]==count(<(0),vals)
    @test r.profiles[1].strict_inside==9
    @test r.profiles[end].strict_inside==0
    @test r.profiles[end].finite==0
    @test isnan(r.profiles[end].median)
    @test length(r.response)==8
    z,v=C.surface_grid(c,f,s)
    original=C.iso_for_mean_height_legacy(v,z,s)
    chosen=first(x for x in r.response if x.valid>0 && x.mean_height_nm<s.target_height_nm)
    @test chosen.isovalue==original.iso
    @test chosen.mean_height_nm==original.mean_height
    @test chosen.valid==original.nvalid
    @test r.response[1].isovalue==10.0^s.isovalue_min_log10
    @test r.summary["frame_normal"]==[0.,0.,1.]
    outside=C.MoldSettings(.1,.1,1.2,0.,1.,.1,-3.,.8,8,"first_below_target")
    absent=A.audit(c,f,outside)
    @test absent.summary["target_valid"]==0
    @test isnan(absent.summary["target_median"])
    mktempdir() do dir
        A.table(joinpath(dir,"profile.tsv"),r.profiles)
        @test length(readlines(joinpath(dir,"profile.tsv")))==length(z)+1
    end
end
@testset "QE cold-smearing derivative can be negative" begin
    @test A.qe_cold_derivative(0.)≈2exp(-.5)/sqrt(pi)
    @test abs(A.qe_cold_derivative(sqrt(2.)))<1e-15
    @test A.qe_cold_derivative(2.)<0
    @test A.qe_cold_derivative(3.)<0
    @test A.qe_cold_derivative(-2.)>0
end
@testset "No grading inputs or template-generation interface" begin
    @test A.options(["--help"])===nothing
    for k in ("--truth","--expected-N","--isovalue","--height","--templates","--manifest")
        @test_throws ErrorException A.options([k,"forbidden"])
    end
    @test_throws ErrorException A.options(String[])
end
