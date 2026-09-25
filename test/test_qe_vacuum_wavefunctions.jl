using Test, LinearAlgebra, TOML, Printf
include(joinpath(@__DIR__,"qe_vacuum_wavefunctions.jl"))
const D=QEVacuumWavefunctions
const W=D.W
const ROOT=normpath(joinpath(@__DIR__,".."))
BLAS.set_num_threads(1)

function write_record(io,items...)
    b=IOBuffer()
    for item in items; write(b,item); end
    payload=take!(b); n=Int32(length(payload))
    write(io,n); write(io,payload); write(io,n)
end
function fixture_wfc(file,miller,coeff,cell;gamma=Int32(-1),scale=1.0)
    open(file,"w") do io
        write_record(io,Int32(1),zeros(3),Int32(1),gamma,scale)
        write_record(io,Int32.([size(miller,2),size(miller,2),1,size(coeff,2)]))
        write_record(io,2pi*inv(transpose(cell)))
        write_record(io,Int32.(miller))
        for col in eachcol(coeff); write_record(io,collect(col)); end
    end
end

@testset "Native Gamma WFC reader and independent analytic Fourier sums" begin
    mktempdir() do dir
        file=joinpath(dir,"wfc.dat")
        miller=[0 1 0; 0 0 -1; 0 0 0]
        coeff=ComplexF64[.3 .2 -.1; .1+.05im .2 .3-.07im; -.04im .1im .1+.03im]
        cell=diagm([10.,12.,14.]); fixture_wfc(file,miller,coeff,cell)
        w=W.read_wfc(file,[1,3]); weights=[2.,2.]
        @test w.selected==[1,3] && w.nbnd==3 && w.ng==3
        @test w.miller==miller && w.coeff==coeff[:,[1,3]] && w.zero==1
        @test w.cell≈cell rtol=1e-15
        @test W.parseval(w,weights)≈sum(2*(abs2(coeff[1,b])+2sum(abs2,coeff[2:end,b])) for b in (1,3))
        fractions=[0. .25 .37 .81 .2; 0. .5 .14 .63 .7; 0. .1 .9 .21 .4]
        explicit=[sum(weights[b]*abs2(w.coeff[1,b]+
            sum(w.coeff[g,b]*cis(2pi*dot(miller[:,g],f))+
                conj(w.coeff[g,b])*cis(-2pi*dot(miller[:,g],f)) for g in 2:3))
            for b in 1:2)/det(cell) for f in eachcol(fractions)]
        dv=W.density(w,fractions,weights;block_points=2)
        @test dv≈explicit rtol=1e-14 atol=1e-18
        @test dv==W.density(w,fractions,weights;block_points=2,parallel=false)
        @test W.density(w,fractions.+[1,2,-3],weights;block_points=3)≈explicit rtol=1e-13
        grid=hcat([[i,j,k]/8 for i in 0:7 for j in 0:7 for k in 0:7]...)
        @test sum(W.density(w,grid,weights;block_points=32))/512*det(cell)≈W.parseval(w,weights) rtol=1e-14
        @test_throws ErrorException W.read_wfc(file,[3,1])
        @test_throws ErrorException W.read_wfc(file,[1,1])
        @test_throws ErrorException W.read_wfc(file,[4])
        @test_throws ErrorException W.density(w,fractions,[-1.,2.];block_points=2)
        @test_throws ErrorException W.density(w,fractions,weights;block_points=0)
        fixture_wfc(file,miller,coeff,cell;gamma=Int32(1))
        @test W.read_wfc(file,[2]).coeff==coeff[:,2:2]
        fixture_wfc(file,miller,coeff,cell;gamma=Int32(0))
        @test_throws ErrorException W.read_wfc(file,[1])
        fixture_wfc(file,miller,coeff,cell;scale=2.)
        @test_throws ErrorException W.read_wfc(file,[1])
        fixture_wfc(file,[0 1 -1; 0 0 0; 0 0 0],coeff,cell)
        @test_throws ErrorException W.read_wfc(file,[1])
        bad=copy(coeff); bad[1,1]+=.1im; fixture_wfc(file,miller,bad,cell)
        @test_throws ErrorException W.read_wfc(file,[1])
        fixture_wfc(file,miller,coeff,cell)
        open(file,"a") do io; write(io,Int32(1)); end
        @test_throws ErrorException W.read_wfc(file,[1])
        @test_throws ErrorException W.record(IOBuffer(reinterpret(UInt8,Int32[43])),44)
    end
end

@testset "Frozen scope and precise-grid support" begin
    config=joinpath(ROOT,"config/qe_vacuum_wavefunctions.toml")
    s=D.settings(config)
    @test s["model"]["sample_bias_ry"]==-0.0220495933
    @test !s["preprocessing"]["normalize_components"] && !s["preprocessing"]["clip_negative_values"]
    @test D.main(["--help"])===nothing
    mktempdir() do dir
        file=joinpath(dir,"bad.toml"); changed=deepcopy(s)
        changed["preprocessing"]["half_nm"]=.64
        open(io->TOML.print(io,changed),file,"w")
        @test_throws ErrorException D.settings(file)
    end
    cell=diagm([10.,10.,10.]); dims=[10,10,10]
    points=transpose([2.3 3.4; 4.1 5.2; 6.7 7.1])*D.G.C.BOHR_NM
    nv=D.native_vertices(cell,dims,points)
    @test length(nv.indices)==15 # The two cells share vertex (3,5,7).
    @test nv.fractions≈[.23 .34; .41 .52; .67 .71]
    @test all(i->all(0 .<= collect(i) .< dims),nv.indices)
    @test_throws ErrorException D.native_vertices(cell,dims,[10. 0 0]*D.G.C.BOHR_NM)
    @test_throws ErrorException D.native_vertices(cell,dims,[-.01 0 0]*D.G.C.BOHR_NM)
    # The source records b_i WITH 2*pi/alat, independent of our parser.
    source=read(joinpath(ROOT,"qe/gamma_source_20260925/PW/src/pw_restart_new.f90"),String)
    source_units=occursin("tpiba*bg(:,1)",replace(source," "=>""))
    @test source_units
    job=read(joinpath(ROOT,"hpc/qe_vacuum_wavefunctions.sbatch"),String)
    @test occursin("#SBATCH --cpus-per-task=4",job)
    @test occursin("#SBATCH --time=01:00:00",job)
    @test occursin("#SBATCH --no-requeue",job)
    @test !occursin("pw.x",job) && !occursin("pp.x",job) && !occursin("sbatch ",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_vacuum_wavefunctions.sbatch"))`)
    @test true
end

@testset "End-to-end independent WFC/cube fixture with unselected flat decay" begin
    mktempdir() do dir
        case_dir=joinpath(dir,"toy"); mkdir(case_dir); mkdir(joinpath(dir,"pseudo"))
        config=joinpath(dir,"settings.toml")
        cp(joinpath(ROOT,"config/qe_vacuum_wavefunctions.toml"),config)
        xml=joinpath(case_dir,"data-file-schema.xml")
        write(xml,"""
        <qes Units="Hartree atomic units"><creator VERSION="7.4.1"/>
        <input><species name="C"><pseudo_file>C.UPF</pseudo_file></species></input>
        <output>
        <convergence_info><wf_collected>true</wf_collected><scf_conv>
        <convergence_achieved>true</convergence_achieved><scf_error>1e-10</scf_error>
        </scf_conv></convergence_info>
        <atomic_structure nat="1" alat="20"><atomic_positions>
        <atom name="C" index="1">0.01 0.01 0.01</atom></atomic_positions>
        <cell><a1>20 0 0</a1><a2>0 20 0</a2><a3>0 0 20</a3></cell></atomic_structure>
        <basis_set><gamma_only>true</gamma_only><ecutwfc>25</ecutwfc>
        <fft_grid nr1="32" nr2="32" nr3="32"></fft_grid></basis_set>
        <band_structure><nks>1</nks><nbnd>3</nbnd><lsda>false</lsda>
        <noncolin>false</noncolin><spinorbit>false</spinorbit>
        <ks_energies><k_point weight="2">0 0 0</k_point>
        <eigenvalues size="3">-0.02 0 0.02</eigenvalues></ks_energies>
        <smearing degauss="0.0001">mv</smearing><fermi_energy>0.005</fermi_energy>
        </band_structure><total_energy><etot>-1</etot></total_energy>
        </output></qes>
        """)
        pseudo=joinpath(dir,"pseudo/C.UPF")
        write(pseudo,"""
        <PP_HEADER element="C" is_paw="true" mesh_size="3"/>
        <PP_R>0.01 0.02 0.1</PP_R>
        <PP_AUGMENTATION cutoff_r_index="3"></PP_AUGMENTATION>
        <PP_BETA.1 cutoff_radius_index="3"></PP_BETA.1>
        """)
        write(joinpath(case_dir,"frame.tsv"),"key\tvalue\norigin_nm\t0.5,0.5,0.2\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
        cell=20Matrix{Float64}(I,3,3)
        miller=[0 1 0; 0 0 1; 0 0 0]
        coeff=ComplexF64[.2 .45 .1; .01 .12+.05im .02; .03 -.02+.03im .04]
        wave=joinpath(dir,"wfc.dat"); fixture_wfc(wave,miller,coeff,cell)
        cube=joinpath(case_dir,"reference_smooth.cube")
        # Independent explicit +/-G series; not the production density function.
        truth(f)=2abs2(coeff[1,2]+sum(coeff[g,2]*cis(2pi*dot(miller[:,g],f))+
            conj(coeff[g,2])*cis(-2pi*dot(miller[:,g],f)) for g in 2:3))/det(cell)
        open(cube,"w") do io
            println(io,"analytic fixture\nnative last-index-fast\n0 0 0 0")
            println(io,"32 0.625000 0 0\n32 0 0.625000 0\n32 0 0 0.625000")
            for ix in 0:31,iy in 0:31,iz in 0:31
                @printf(io,"%.4e\n",truth([ix,iy,iz]/32))
            end
        end
        meta=Dict("config_sha256"=>D.S.sha(config),"xml_sha256"=>D.S.sha(xml),
            "smooth_sha256"=>D.S.sha(cube),"wfc_sha256"=>D.S.sha(wave),
            "pseudo_sha256"=>Dict("C.UPF"=>D.S.sha(pseudo)))
        open(io->TOML.print(io,meta),joinpath(case_dir,"metadata.toml"),"w")
        @test D.compare(case_dir,config,wave)===nothing
        summary=TOML.parsefile(joinpath(case_dir,"comparison/summary.toml"))
        @test all(values(summary["checks"]))
        @test summary["selected_bands"]==[2] && summary["selected_weights"]==[2.]
        @test summary["native_violations"]==0
        @test summary["plane_points"]==3*289
        planes=D.V.tsv(joinpath(case_dir,"comparison/planes.tsv"))
        @test length(planes)==3*289
        for row in planes
            i=parse(Int,row["pixel"])-1
            point=[.5+(i%17-8)*.04,.5+(i÷17-8)*.04,.2+parse(Float64,row["height_nm"])]/D.G.C.BOHR_NM
            @test parse(Float64,row["direct"])≈truth(cell\point) rtol=1e-13 atol=1e-18
        end
        decay=D.V.tsv(joinpath(case_dir,"comparison/decay.tsv"))
        @test length(decay)==2*289
        @test all(row->abs(parse(Float64,row["log_decay_per_nm"]))<1e-12,decay)
        @test all(row->parse(Float64,row["interpolation_relative_l2"])>0,
            D.V.tsv(joinpath(case_dir,"comparison/plane_summary.tsv")))
        @test_throws ErrorException D.compare(case_dir,config,wave) # No overwrite.
    end
end
