using Test, LinearAlgebra, TOML
include(joinpath(@__DIR__,"qe_vacuum_continuous.jl"))
const C=QEVacuumContinuous
const ROOT=normpath(joinpath(@__DIR__,".."))

function record(io,items...)
    buffer=IOBuffer(); foreach(x->write(buffer,x),items); bytes=take!(buffer)
    write(io,Int32(length(bytes))); write(io,bytes); write(io,Int32(length(bytes)))
end

@testset "Complete synthetic saved-profile to continuous-surface workflow" begin
    mktempdir() do run
        config=joinpath(run,"settings.toml"); cp(joinpath(ROOT,"config/qe_vacuum_continuous.toml"),config)
        mkdir(joinpath(run,"pseudo")); dir=joinpath(run,"glcn"); mkdir(dir)
        for element in ("C","Cu")
            write(joinpath(run,"pseudo",element*".UPF"),"""
            <PP_HEADER element="$element" is_paw="true" mesh_size="3"/>
            <PP_R>0.01 0.02 0.1</PP_R>
            <PP_AUGMENTATION cutoff_r_index="3"></PP_AUGMENTATION>
            <PP_BETA.1 cutoff_radius_index="3"></PP_BETA.1>
            """)
        end
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
        wave=joinpath(run,"wfc.dat"); cell=20Matrix{Float64}(I,3,3)
        open(wave,"w") do io
            record(io,Int32(1),zeros(3),Int32(1),Int32(-1),1.)
            record(io,Int32.([2,2,1,3])); record(io,2pi*inv(transpose(cell)))
            record(io,Int32.([0 0; 0 0; 0 1]))
            for _ in 1:3; record(io,ComplexF64[1.,.3]); end
        end
        s=C.settings(config); g=C.Q.geometry(run,"glcn",s); z=g.domain.knots
        truth(h)=2*(1+.6cos(2pi*h/(20C.D.G.C.BOHR_NM)))^2/8000
        profiles=[(;pixel=p,z_nm=h,x_nm=g.xy[p,1],y_nm=g.xy[p,2],direct=truth(h))
            for p in axes(g.xy,1) for h in z]
        C.Q.table(joinpath(run,"glcn_profiles.tsv"),keys(first(profiles)),profiles)
        isos=[.0001,.00015,.00025]
        common=[(;source="joint",domain="molecular_half_paw_gap",interval=j,isovalue=iso)
            for (j,iso) in enumerate(isos)]
        C.Q.table(joinpath(run,"common_intervals.tsv"),keys(first(common)),common)
        saved=[(;molecule="glcn",pixel=p,interval=j,
            z_nm=C.Q.X.crossing(z[g.domain.half],truth.(z[g.domain.half]),iso).z)
            for p in axes(g.xy,1) for (j,iso) in enumerate(isos)]
        C.Q.table(joinpath(run,"surfaces.tsv"),keys(first(saved)),saved)
        meta=Dict("config_sha256"=>C.S.sha(config),"xml_sha256"=>C.S.sha(xml),
            "frame_sha256"=>C.S.sha(frame),"wfc_sha256"=>C.S.sha(wave),
            "references"=>Dict(n=>C.S.sha(joinpath(run,n)) for n in ("glcn_profiles.tsv","surfaces.tsv","common_intervals.tsv")),
            "pseudo_sha256"=>Dict(e*".UPF"=>C.S.sha(joinpath(run,"pseudo",e*".UPF")) for e in ("C","Cu")))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
        @test C.analyze(run,"glcn",wave)===nothing
        out=joinpath(dir,"continuous"); result=TOML.parsefile(joinpath(out,"summary.toml"))
        @test all(values(result["checks"])) && result["certified_complete_maps"]
        @test result["valid_roots"]==result["targets"]==867
        @test result["unresolved_leaves"]==0
        @test result["max_abs_relative_density_error"]<1e-5
        for r in C.D.V.tsv(joinpath(out,"surfaces.tsv"))
            iso=C.number(r,"isovalue"); zl=C.number(r,"z_lower_nm"); zh=C.number(r,"z_upper_nm")
            exact=acos((sqrt(iso*4000)-1)/.6)/(2pi)*20C.D.G.C.BOHR_NM
            @test zl<=exact<=zh
            @test zh-zl<=s["selection"]["root_width_nm"]
            @test C.number(r,"direct")==C.number(r,"repeat")
        end
        @test_throws ErrorException C.analyze(run,"glcn",wave)
        changed=deepcopy(s); changed["selection"]["taylor_order"]=8
        bad=joinpath(run,"bad.toml"); open(io->TOML.print(io,changed),bad,"w")
        @test_throws ErrorException C.settings(bad)
    end
end

@testset "CLI and bounded batch" begin
    @test C.main(["--help"])===nothing
    job=read(joinpath(ROOT,"hpc/qe_vacuum_continuous.sbatch"),String)
    @test occursin("#SBATCH --cpus-per-task=4",job) && occursin("#SBATCH --time=01:00:00",job)
    @test occursin("#SBATCH --no-requeue",job) && !occursin("sbatch ",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_vacuum_continuous.sbatch"))`)
    @test true
end
