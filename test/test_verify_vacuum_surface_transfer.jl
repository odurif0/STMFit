# Existing synthetic workflow also supplies its immutable toy image constructor.
include(joinpath(@__DIR__,"test_vacuum_surface_transfer.jl"))
include(joinpath(@__DIR__,"verify_vacuum_surface_transfer.jl"))
const V=VerifyVacuumSurfaceTransfer

@testset "Independent saved-prediction checker" begin
    grid=collect(-.32:.04:.32)
    a=[1+.03exp(-(t/.12)^2-(u/.12)^2) for u in grid,t in grid]
    b=[1.025+.02exp(-(t/.08)^2-((u+.05)/.1)^2) for u in grid,t in grid]
    maps=(heights=Dict(j=>Dict(-1=>(a+b)/2,0=>a,1=>b) for j in 1:3),isos=Dict(j=>Float64(j) for j in 1:3))
    mktempdir() do dir
        out=joinpath(dir,"saved"); R.one_direction(toy_image(),"fwd",maps,C,out)
        verified=V.verify_direction(out,maps,C)
        @test verified.patches>=2
        @test length(verified.models)==6
    end
end
