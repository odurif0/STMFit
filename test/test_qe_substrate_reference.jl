using Test, LinearAlgebra, TOML
include(joinpath(@__DIR__,"verify_qe_substrate_reference.jl"))
const R=VerifyQESubstrateReference.R
const G=R.G
const ROOT=normpath(joinpath(@__DIR__,".."))

@testset "Periodic projection and complete segment geometry" begin
    geo=(;cell_nm=[2.,2.,4.],atoms=[(;z=29,position_nm=[0.,0.,0.]),
        (;z=29,position_nm=[1.,1.,.5]),(;z=6,position_nm=[1.98,.1,1.])])
    radii=Dict(29=>.1,6=>.2)
    @test G.projection_margin([.02,.1],geo,radii)≈-.16
    @test G.projection_margin([2.02,.1],geo,radii)≈-.16
    @test G.projection_margin([1.98,.3],geo,radii)≈0 atol=1e-15
    grid=G.grid([geo,geo],[radii,radii],[32,24])
    @test length(grid)==768 && [p.pixel for p in grid]==collect(1:768)
    @test allunique((p.fx,p.fy) for p in grid)
    @test all(p->p.eligible==(p.glcn_margin_nm>0 && p.glcnac_margin_nm>0),grid)
    for p in grid
        explicit=minimum(norm([p.x_nm,p.y_nm]-a.position_nm[1:2]-[ix,iy].*geo.cell_nm[1:2])-radii[a.z]
            for a in geo.atoms if a.z!=29 for ix in -1:1 for iy in -1:1)
        @test p.glcn_margin_nm≈explicit atol=1e-14
    end
    shifted=merge(geo,(;atoms=[geo.atoms[1:2]...,(;z=6,position_nm=[.8,.8,1.])]))
    pair=G.grid([geo,shifted],[radii,radii],[32,24])
    @test count(p->p.eligible,pair)<count(p->p.eligible,grid)
    @test any(p->p.glcn_margin_nm>0 && p.glcnac_margin_nm<0 && !p.eligible,pair)
    bad=merge(geo,(;cell_nm=[2.1,2.,4.]))
    @test_throws ErrorException G.grid([geo,bad],[radii,radii],[32,24])
    @test_throws ErrorException G.grid([geo,geo],[radii,radii],[0,24])
    old=R.C.Q.X.domains(geo,radii,100)
    ds=G.domains(geo,radii,100,old)
    @test [d.domain for d in ds]==["full_cu_gap","cu_half_gap","saved_molecular_half"]
    @test ds[1].lower_paw_nm≈.6 && ds[1].upper_paw_nm≈3.9
    @test .6<ds[1].lo<ds[2].hi<ds[2].midpoint_nm<ds[1].hi<3.9
    @test ds[3].lo==first(old.knots) && ds[3].hi==old.knots[last(old.half)]
    # Endpoints are clear, but the middle of this column intersects a sphere.
    @test R.D.P.clearance([1.98 .1 .7; 1.98 .1 1.3],geo,radii).inside_or_boundary==0
    @test G.segment_margin([1.98,.1],.7,1.3,geo,radii)≈-.2
    for p in filter(p->p.eligible,grid), d in ds
        v=G.segment_margin([p.x_nm,p.y_nm],d.lo,d.hi,geo,radii)
        @test v>0
        for z in (d.lo,(d.lo+d.hi)/2,d.hi)
            @test R.D.P.clearance([p.x_nm p.y_nm z],geo,radii).minimum_clearance_nm>=v-1e-14
        end
    end
end

function record(io,items...)
    b=IOBuffer(); foreach(x->write(b,x),items); data=take!(b)
    write(io,Int32(length(data))); write(io,data); write(io,Int32(length(data)))
end

@testset "Analytic WFC witness keeps ascending/full-gap failures" begin
    mktempdir() do run
        cp(joinpath(ROOT,"config/qe_substrate_reference.toml"),joinpath(run,"settings.toml"))
        mkdir(joinpath(run,"pseudo"))
        for e in ("C","Cu")
            write(joinpath(run,"pseudo",e*".UPF"),"""
            <PP_HEADER element="$e" is_paw="true" mesh_size="3"/>
            <PP_R>0.01 0.02 0.1</PP_R>
            <PP_AUGMENTATION cutoff_r_index="3"></PP_AUGMENTATION>
            <PP_BETA.1 cutoff_radius_index="3"></PP_BETA.1>
            """)
        end
        for m in R.MOLECULES
            dir=joinpath(run,m); mkdir(dir)
            write(joinpath(dir,"data-file-schema.xml"),"""
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
            write(joinpath(dir,"frame.tsv"),"key\tvalue\norigin_nm\t0.5,0.5,0.2\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
            write(joinpath(dir,"molecular_surfaces.tsv"),"fixture\nnot used by inference\n")
        end
        # An archived but inactive alternate Cu potential must not overwrite
        # the actual XML-referenced Cu radius in the independent reader.
        write(joinpath(run,"pseudo","Cu_other.UPF"),replace(read(joinpath(run,"pseudo","Cu.UPF"),String),
            "0.01 0.02 0.1"=>"0.01 0.02 0.3"))
        independent=VerifyQESubstrateReference.independent_geometry(run,"glcn")
        @test independent.radii["Cu"]==.1R.D.G.C.BOHR_NM
        @test length(independent.radii)==2
        wfc=joinpath(run,"wfc.dat"); cell=20Matrix{Float64}(I,3,3)
        open(wfc,"w") do io
            record(io,Int32(1),zeros(3),Int32(1),Int32(-1),1.)
            record(io,Int32.([2,2,1,3])); record(io,2pi*inv(transpose(cell)))
            record(io,Int32.([0 0;0 0;0 1]))
            for _ in 1:3; record(io,ComplexF64[1.,.3]); end
        end
        s=R.settings(joinpath(run,"settings.toml")); g=R.geometry(run,s)
        R.table(joinpath(run,"grid.tsv"),g.grid)
        R.table(joinpath(run,"domains.tsv"),[merge((;molecule=m),d) for (k,m) in enumerate(R.MOLECULES) for d in g.domains[k]])
        R.table(joinpath(run,"segment_clearance.tsv"),g.margins)
        isos=[.0001,.00015,.00025]
        R.table(joinpath(run,"common_intervals.tsv"),[(;source="joint",domain="molecular_half_paw_gap",interval=j,isovalue=iso)
            for (j,iso) in enumerate(isos)])
        dir=joinpath(run,"glcn")
        meta=Dict("config_sha256"=>R.S.sha(joinpath(run,"settings.toml")),"wfc_sha256"=>R.S.sha(wfc),
            "xml_sha256"=>R.S.sha(joinpath(dir,"data-file-schema.xml")),"frame_sha256"=>R.S.sha(joinpath(dir,"frame.tsv")),
            "molecular_surface_sha256"=>R.S.sha(joinpath(dir,"molecular_surfaces.tsv")),
            "pseudo_sha256"=>Dict(e*".UPF"=>R.S.sha(joinpath(run,"pseudo",e*".UPF")) for e in ("C","Cu")),
            "geometry_sha256"=>Dict(name=>R.S.sha(joinpath(run,name)) for name in
                ("grid.tsv","domains.tsv","segment_clearance.tsv","common_intervals.tsv")))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
        @test R.analyze(run,"glcn",wfc;witness=true)===nothing
        roots=R.V.tsv(joinpath(dir,"witness/roots.tsv"))
        @test length(roots)==9 && all(r->R.i(r,"pixel")==first(filter(p->p.eligible,g.grid)).pixel,roots)
        for r in roots
            if r["domain"]=="full_cu_gap"
                @test !parse(Bool,r["valid"]) && R.i(r,"root_count")==2
                @test R.i(r,"descending")==R.i(r,"ascending")==1
                @test isnan(R.n(r,"z_nm"))
            elseif r["domain"]=="cu_half_gap" && R.n(r,"isovalue")==.00025
                # This analytic root lies exactly at a subdivision boundary.
                # The unchanged strict endpoint-sign test must retain ambiguity.
                @test !parse(Bool,r["valid"]) && R.i(r,"unresolved")==2
                @test R.i(r,"root_count")==0 && isnan(R.n(r,"z_nm"))
                ls=filter(l->l["domain"]==r["domain"] && R.i(l,"interval")==R.i(r,"interval") && l["status"]=="unresolved",
                    R.V.tsv(joinpath(dir,"witness/leaves.tsv")))
                exact=20R.D.G.C.BOHR_NM/4
                @test minimum(R.n(l,"lo") for l in ls)<=exact<=maximum(R.n(l,"hi") for l in ls)
            else
                @test parse(Bool,r["valid"]) && R.i(r,"root_count")==1
                exact=acos((sqrt(R.n(r,"isovalue")*4000)-1)/.6)/(2pi)*20R.D.G.C.BOHR_NM
                @test R.n(r,"z_lower_nm")<=exact<=R.n(r,"z_upper_nm")
            end
            @test R.i(r,"unresolved")== (r["domain"]=="cu_half_gap" && R.n(r,"isovalue")==.00025 ? 2 : 0)
        end
        summary=TOML.parsefile(joinpath(dir,"witness/summary.toml"))
        @test all(values(summary["checks"])) && summary["valid_roots"]==5
        @test summary["unresolved_leaves"]==2
        @test_throws ErrorException R.analyze(run,"glcn",wfc;witness=true)
        changed=deepcopy(s); changed["preprocessing"]["lateral_grid"]=[16,12]
        open(io->TOML.print(io,changed),joinpath(run,"bad.toml"),"w")
        @test_throws ErrorException R.settings(joinpath(run,"bad.toml"))
    end
end

@testset "Bounded CLI" begin
    @test R.main(["--help"])===nothing
    job=read(joinpath(ROOT,"hpc/qe_substrate_reference.sbatch"),String)
    @test occursin("#SBATCH --cpus-per-task=4",job) && occursin("#SBATCH --time=01:00:00",job)
    @test occursin("#SBATCH --no-requeue",job) && !occursin("sbatch ",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_substrate_reference.sbatch"))`)
    @test true
end
