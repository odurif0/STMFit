using Test, LinearAlgebra, TOML, Printf
include(joinpath(@__DIR__,"qe_vacuum_crossings.jl"))
const Q=QEVacuumCrossings
const X=Q.X
const ROOT=normpath(joinpath(@__DIR__,".."))

@testset "Complete open crossing intervals, including singular levels" begin
    z=collect(0.:3.)
    @test X.intervals(z,[4.,3.,2.,1.])==[(1.,4.)]
    @test X.crossing(z,[4.,3.,2.,1.],3.).z==1.
    @test X.crossing(z,[4.,3.,2.,1.],4.).status==:boundary
    @test isempty(X.intervals(z,[1.,2.,3.,4.]))
    @test X.crossing(z,[1.,2.,3.,4.],2.5).status==:ascending
    @test X.intervals(z,[4.,3.,3.,2.])==[(2.,3.),(3.,4.)]
    @test X.crossing(z,[4.,3.,3.,2.],3.).status==:plateau
    @test X.intervals(z,[4.,1.,2.,0.])==[(0.,1.),(2.,4.)]
    @test X.crossing(z,[4.,1.,2.,0.],1.5).status==:multiple
    @test X.crossing(z,[4.,2.,4.,5.],2.).status==:tangent
    @test X.crossing(z,[4.,3.,2.,1.],0.).status==:absent
    @test isempty(X.intervals(z,fill(2.,4)))
    @test X.intersect_intervals([(0.,1.),(2.,4.)],[(1.,3.)])==[(2.,3.)]
    @test isempty(X.intersect_intervals([(0.,1.)],[(1.,2.)]))
    @test_throws ErrorException X.intervals([0.,0.],[1.,2.])
    @test_throws ErrorException X.intervals([0.,1.],[NaN,2.])
    # Exhaust all 243 five-knot profiles over three density levels. Test
    # singular levels as well as interiors against explicit root counting.
    for tuple in Iterators.product(ntuple(_->0.:2.,5)...)
        v=collect(tuple); zs=collect(0.:4.); ranges=X.intervals(zs,v)
        for iso in 0.:.25:2.
            @test any(lo<iso<hi for (lo,hi) in ranges)==(X.crossing(zs,v,iso).status==:found)
        end
        for (lo,hi) in ranges
            iso=(lo+hi)/2; r=X.crossing(zs,v,iso)
            @test r.status==:found
            k=clamp(floor(Int,r.z)+1,1,4)
            @test v[k]+(r.z-zs[k])*(v[k+1]-v[k])≈iso
        end
    end
end

@testset "End-to-end synthetic periodic field, empty full gap and valid half gap" begin
    mktempdir() do run
        cp(joinpath(ROOT,"config/qe_vacuum_crossings.toml"),joinpath(run,"settings.toml"))
        mkdir(joinpath(run,"pseudo"))
        for element in ("C","Cu")
            write(joinpath(run,"pseudo",element*".UPF"),"""
            <PP_HEADER element="$element" is_paw="true" mesh_size="3"/>
            <PP_R>0.01 0.02 0.1</PP_R>
            <PP_AUGMENTATION cutoff_r_index="3"></PP_AUGMENTATION>
            <PP_BETA.1 cutoff_radius_index="3"></PP_BETA.1>
            """)
        end
        function record(io,items...)
            b=IOBuffer(); foreach(x->write(b,x),items); data=take!(b)
            write(io,Int32(length(data))); write(io,data); write(io,Int32(length(data)))
        end
        wave=joinpath(run,"wfc.dat"); cell=20Matrix{Float64}(I,3,3)
        open(wave,"w") do io
            record(io,Int32(1),zeros(3),Int32(1),Int32(-1),1.)
            record(io,Int32.([2,2,1,3])); record(io,2pi*inv(transpose(cell)))
            record(io,Int32.([0 0; 0 0; 0 1]))
            for _ in 1:3; record(io,ComplexF64[1.,.3]); end
        end
        for m in Q.MOLECULES
            dir=joinpath(run,m); mkdir(dir)
            xml=joinpath(dir,"data-file-schema.xml")
            write(xml,"""
            <qes Units="Hartree atomic units"><creator VERSION="7.4.1"/>
            <input><species name="C"><pseudo_file>C.UPF</pseudo_file></species>
            <species name="Cu"><pseudo_file>Cu.UPF</pseudo_file></species></input>
            <output><convergence_info><wf_collected>true</wf_collected><scf_conv>
            <convergence_achieved>true</convergence_achieved><scf_error>1e-10</scf_error>
            </scf_conv></convergence_info><atomic_structure nat="2" alat="20">
            <atomic_positions><atom name="Cu" index="1">0 0 0</atom>
            <atom name="C" index="2">10 10 3</atom></atomic_positions>
            <cell><a1>20 0 0</a1><a2>0 20 0</a2><a3>0 0 20</a3></cell></atomic_structure>
            <basis_set><gamma_only>true</gamma_only><ecutwfc>25</ecutwfc>
            <fft_grid nr1="32" nr2="32" nr3="32"></fft_grid></basis_set>
            <band_structure><nks>1</nks><nbnd>3</nbnd><lsda>false</lsda>
            <noncolin>false</noncolin><spinorbit>false</spinorbit>
            <ks_energies><k_point weight="2">0 0 0</k_point>
            <eigenvalues size="3">-0.02 0 0.02</eigenvalues></ks_energies>
            <smearing degauss="0.0001">mv</smearing><fermi_energy>0.005</fermi_energy>
            </band_structure><total_energy><etot>-1</etot></total_energy></output></qes>
            """)
            frame=joinpath(dir,"frame.tsv")
            write(frame,"key\tvalue\norigin_nm\t0.5,0.5,0.2\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
            cube=joinpath(dir,"reference_smooth.cube")
            open(cube,"w") do io
                println(io,"analytic field\nnative order\n0 0 0 0")
                println(io,"32 0.625000 0 0\n32 0 0.625000 0\n32 0 0 0.625000")
                for ix in 0:31,iy in 0:31,iz in 0:31
                    @printf(io,"%.4e\n",2*(1+.6cos(2pi*iz/32))^2/8000)
                end
            end
            meta=Dict("config_sha256"=>Q.S.sha(joinpath(run,"settings.toml")),
                "xml_sha256"=>Q.S.sha(xml),"smooth_sha256"=>Q.S.sha(cube),
                "frame_sha256"=>Q.S.sha(frame),"wfc_sha256"=>Q.S.sha(wave),
                "pseudo_sha256"=>Dict(e*".UPF"=>Q.S.sha(joinpath(run,"pseudo",e*".UPF")) for e in ("C","Cu")))
            open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
        end
        @test Q.analyze(run,wave,wave)===nothing
        out=joinpath(run,"analysis"); result=TOML.parsefile(joinpath(out,"summary.toml"))
        @test all(values(result["checks"]))
        @test !result["direct_between_knots_uniqueness_proven"]
        @test result["surface_rows"]==2*289
        rows=Q.D.V.tsv(joinpath(out,"common_intervals.tsv"))
        @test all(r->parse(Int,r["interval"])==0,filter(r->r["domain"]=="full_paw_gap",rows))
        @test all(r->parse(Int,r["interval"])==1,filter(r->r["domain"]=="molecular_half_paw_gap",rows))
        surfaces=Q.D.V.tsv(joinpath(out,"surfaces.tsv"))
        for row in surfaces
            z=parse(Float64,row["z_nm"]); direct=parse(Float64,row["direct_at_surface"])
            @test direct≈2*(1+.6cos(2pi*z/(20Q.D.G.C.BOHR_NM)))^2/8000 rtol=1e-13
        end
        @test_throws ErrorException Q.analyze(run,wave,wave)
    end
end

@testset "Geometry domains and vertical bilinear profiles" begin
    geo=(;cell_nm=[2.,2.,3.],atoms=[(;z=29,position_nm=[0.,0.,0.]),(;z=6,position_nm=[1.,1.,1.])])
    r=Dict(29=>.1,6=>.2); d=X.domains(geo,r,30)
    @test d.lower≈1.2 && d.upper≈2.9 && d.midpoint≈2.05
    @test all(d.lower .< d.knots .< d.upper)
    @test all(d.knots[d.half].<d.midpoint)
    @test d.knots[last(d.half)+1]>=d.midpoint
    @test_throws ErrorException X.domains(geo,Dict(29=>2.,6=>2.),30)
    frame=(;origin_nm=[1.,1.,.7],t_axis=[1.,.1,.3])
    xy=X.lateral_points(frame,.32,.04)
    @test size(xy)==(289,3) && all(xy[:,3].==.7)
    @test xy[145,:]≈frame.origin_nm
    @test norm(xy[2,:]-xy[1,:])≈.04
    @test dot(xy[2,:]-xy[1,:],xy[18,:]-xy[1,:])≈0 atol=1e-16
    cube=(;dims=[4,5,6],grid=[2x+3y+5z for z in 0:5,y in 0:4,x in 0:3])
    cell=[4.,5.,6.]; pts=[.2 .4 0.; 1.3 2.1 0.]; zs=[1.,2.,3.]
    vals=X.profiles(cube,cell,pts,zs)
    @test vals≈[2pts[p,1]+3pts[p,2]+5h for h in zs,p in 1:2]
    @test_throws ErrorException X.profiles(cube,cell,[3.2 0 0],zs)
    @test_throws ErrorException X.profiles(cube,cell,pts,[6.])
    @test_throws ErrorException X.profiles(cube,cell,pts,[1.1,2.])
end

@testset "Frozen scope and all-outcome serialization" begin
    config=joinpath(ROOT,"config/qe_vacuum_crossings.toml"); s=Q.settings(config)
    @test s["selection"]["domains"]==["full_paw_gap","molecular_half_paw_gap"]
    @test !s["preprocessing"]["clip_negative_values"]
    @test Q.main(["--help"])===nothing
    mktempdir() do dir
        changed=deepcopy(s); changed["preprocessing"]["half_nm"]=.64
        file=joinpath(dir,"bad.toml"); open(io->TOML.print(io,changed),file,"w")
        @test_throws ErrorException Q.settings(file)
        file=joinpath(dir,"empty.tsv"); Q.table(file,(:a,:b),NamedTuple[])
        @test read(file,String)=="a\tb\n"
        Q.table(file,(:a,:b),[(;a=NaN,b=0)])
        @test read(file,String)=="a\tb\nNaN\t0\n"
    end
    job=read(joinpath(ROOT,"hpc/qe_vacuum_crossings.sbatch"),String)
    @test occursin("#SBATCH --cpus-per-task=4",job)
    @test occursin("#SBATCH --time=01:00:00",job)
    @test occursin("#SBATCH --no-requeue",job)
    @test !occursin("pw.x",job) && !occursin("pp.x",job) && !occursin("sbatch ",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_vacuum_crossings.sbatch"))`)
    @test true
end
