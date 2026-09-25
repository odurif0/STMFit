using Test, LinearAlgebra
include(joinpath(@__DIR__,"qe_paw_plane_clearance.jl"))
const P=QEPAWPlaneClearance

@testset "Header-only PAW radius and periodic plane geometry" begin
    @test P.main(["--help"])===nothing
    mktempdir() do d
        cube=joinpath(d,"one.cube")
        header="fixture\nno density payload needed\n1 0 0 0\n10 1 0 0\n10 0 1 0\n10 0 0 1\n6 6 9.8 0 0\n"
        write(cube,header)
        geo=P.cube_atoms(cube)
        @test geo.cell_nm≈fill(10P.G.C.BOHR_NM,3)
        @test only(geo.atoms).z==6
        @test only(geo.atoms).position_nm≈[9.8,0,0]*P.G.C.BOHR_NM
        write(cube,replace(header,"10 1 0 0"=>"10 1 0.1 0"))
        @test_throws ErrorException P.cube_atoms(cube)
        write(cube,replace(header,"1 0 0 0\n10"=>"-1 0 0 0\n10"))
        @test_throws ErrorException P.cube_atoms(cube)
        pp=joinpath(d,"C.UPF")
        upf="""
        <PP_HEADER element=" C" is_paw="true" mesh_size="4"/>
        <PP_R>0.01 1 2 3</PP_R>
        <PP_AUGMENTATION cutoff_r="-1" cutoff_r_index="2"></PP_AUGMENTATION>
        <PP_BETA.1 cutoff_radius_index="3"></PP_BETA.1>
        """
        write(pp,upf)
        radius=P.paw_radius(pp)
        @test radius.z==6 && radius.outer_index==3
        @test radius.radius_nm==2P.G.C.BOHR_NM
        write(pp,replace(upf,"cutoff_r_index=\"2\""=>"cutoff_r_index=\"4\""))
        @test P.paw_radius(pp).radius_nm==3P.G.C.BOHR_NM
        write(pp,replace(upf,"cutoff_r_index=\"2\""=>"cutoff_r_index=\"8\""))
        @test_throws ErrorException P.paw_radius(pp)
        write(pp,replace(upf,"is_paw=\"true\""=>"is_paw=\"false\""))
        @test_throws ErrorException P.paw_radius(pp)
    end
    geometry=(cell_nm=[10.,10.,10.],atoms=[(z=6,position_nm=[9.8,0,0])])
    radii=Dict(6=>.1)
    @test P.clearance([.1 0 0],geometry,radii).minimum_clearance_nm≈.2
    @test P.clearance([19.9 0 0],geometry,radii).minimum_clearance_nm≈0 atol=1e-14
    @test P.clearance([9.8 0 0],geometry,radii).inside_or_boundary==1
    # Independent explicit neighboring-cell enumeration validates minimum-image geometry.
    points=[.1 .3 .5; 9.8 9.9 .1; 4. 4. 4.]
    exact=[minimum(norm(point-atom.position_nm-[i,j,k].*geometry.cell_nm)-radii[atom.z]
            for atom in geometry.atoms for i in -1:1 for j in -1:1 for k in -1:1)
        for point in eachrow(points)]
    result=P.clearance(points,geometry,radii)
    @test result.points==3
    @test result.minimum_clearance_nm≈minimum(exact)
    @test result.inside_or_boundary==count(<=(0),exact)
end
