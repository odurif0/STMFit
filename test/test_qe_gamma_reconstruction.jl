using Test, TOML, Statistics, Printf
include(joinpath(@__DIR__,"qe_gamma_reconstruction.jl"))
include(joinpath(@__DIR__,"lib/qe_preflight.jl"))
const G=QEGammaReconstruction
const ROOT=dirname(@__DIR__)
const CONFIG=joinpath(ROOT,"config/qe_gamma_reconstruction.toml")

@testset "Frozen Gamma reconstruction and inherited spatial grid" begin
    s=G.settings(CONFIG)
    old=TOML.parsefile(joinpath(ROOT,"config/unit_assignment_patch_support.toml"))["model"]
    @test s["preprocessing"]["half_nm"]==old["mold_half_nm"]==.32
    @test s["preprocessing"]["step_nm"]==old["mold_step_nm"]==.04
    @test s["model"]["scf_acceptance_ry"]==5e-5
    @test !s["preprocessing"]["clip_negative_values"]
    @test G.main(["--help"])===nothing
    @test_throws ErrorException G.main(["--truth","NKNNKN"])
    mktempdir() do d
        for (section,key,value) in (("model","qe_commit","latest"),("model","sample_bias_ry",-.03),
            ("preprocessing","half_nm",.64),("preprocessing","clip_negative_values",true),
            ("selection","expected_N",6),("selection","require_patched_repeat_sha256",false),
            ("selection","cube_significant_digits",4),("selection","numeric_roundoff_atol",1e-6))
            bad=deepcopy(s); bad[section][key]=value
            file=joinpath(d,"bad.toml"); open(io->TOML.print(io,bad),file,"w")
            @test_throws ErrorException G.settings(file)
        end
    end
end

include(joinpath(@__DIR__,"verify_qe_gamma_reconstruction.jl"))
@testset "Synthetic complete PP comparison and independent verification" begin
    mktempdir() do d
        case=joinpath(d,"synthetic"); mkdir(case)
        s=G.settings(CONFIG)
        xml="""
        <qes:espresso Units="Hartree atomic units"><creator VERSION="7.4.1"/>
        <output><basis_set><gamma_only>true</gamma_only></basis_set>
        <convergence_info><scf_conv><convergence_achieved>true</convergence_achieved>
        <scf_error>1e-6</scf_error></scf_conv><wf_collected>true</wf_collected></convergence_info>
        <total_energy><etot>-10</etot></total_energy><band_structure>
        <nks>1</nks><nbnd>9</nbnd><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit>
        <fermi_energy>0</fermi_energy><smearing degauss="0.01">mv</smearing>
        <ks_energies><k_point weight="2">0 0 0</k_point>
        <eigenvalues size="9">-.09 -.04 -.03 -.02 -.01 -.005 .001 .02 .04</eigenvalues>
        </ks_energies></band_structure></output></qes:espresso>
        """
        xmlfile=joinpath(case,"data-file-schema.xml"); write(xmlfile,xml)
        sp=G.S.spectrum(xmlfile,s)
        open(joinpath(case,"metadata.toml"),"w") do io
            TOML.print(io,Dict("config_sha256"=>G.S.sha(CONFIG),"xml_sha256"=>G.S.sha(xmlfile),"ildos_integral_expected"=>4.))
        end
        write(joinpath(case,"frame.tsv"),"key\tvalue\norigin_nm\t0.5,0.5,0.1\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
        step=.5/G.C.BOHR_NM; average=4/(27step^3)
        full=[average*(1+.1cos(2pi*k/3)+.2sin(2pi*i/3)) for i in 0:2 for j in 0:2 for k in 0:2]
        half=(full .+ mean(full))/2
        for name in ("site_control","stock_control","patched_control","stock_ildos","patched_ildos","repeat_ildos","reference")
            density=name=="stock_ildos" ? half : occursin("ildos",name) ? full : -full
            open(joinpath(case,name*".cube"),"w") do io
                println(io,"synthetic cube\nQE last axis fast\n0 0 0 0")
                for k in 1:3
                    axis=zeros(3); axis[k]=step
                    @printf(io,"3 %.17g %.17g %.17g\n",axis...)
                end
                for v in density; @printf(io,"%13.4E\n",v); end
            end
            name=="reference" || write(joinpath(case,name*".in"),G.S.pp_input("synthetic","checkpoint",name,occursin("ildos",name) ? 10 : 5,sp,s))
        end
        @test G.compare(case,CONFIG)===nothing
        @test all(values(TOML.parsefile(joinpath(case,"comparison/summary.toml"))["checks"]))
        VerifyQEGamma.verify_case(case,CONFIG)
        @test_throws ErrorException G.compare(case,CONFIG) # Never overwrite evidence.
    end
end

@testset "Cube-precision error envelopes do not fit tolerances to outcomes" begin
    s=G.settings(CONFIG)
    @test G.quantum(0.,5)==0
    @test G.quantum(.26387e-3,5)≈1e-8
    @test G.quantum(-.19854e-1,5)≈1e-6
    @test G.quantum(1.,5)≈1e-4
    a=[.1,-.025,2e-9,0.]
    @test G.quantized_difference(a,a,s).violations==0
    @test G.quantized_difference(a,a.*1.01,s).violations==3
    @test_throws ErrorException G.quantized_difference(a,[1.],s)
    full=[2+.4cos(2pi*j/17)+.2sin(4pi*j/17) for j in 0:16]
    half=(full .+ mean(full))/2
    @test G.half_pair_relation(half,full,s).violations==0
    @test G.half_pair_relation(round.(half;sigdigits=5),round.(full;sigdigits=5),s).violations==0
    @test G.half_pair_relation(full,full,s).violations>0
    @test G.half_pair_relation(half.+.01,full,s).violations==17
    @test_throws ErrorException G.half_pair_relation([1.],full,s)
end

function fixture(d,scriptname,commands)
    mkdir(joinpath(d,"pseudo")); write(joinpath(d,"pseudo/H.UPF"),"preflight fixture")
    scf="""
    &CONTROL
    prefix = 'h'
    pseudo_dir = './pseudo'
    /
    &SYSTEM
    nat = 1
    ntyp = 1
    /
    ATOMIC_SPECIES
    H 1 H.UPF
    ATOMIC_POSITIONS angstrom
    H 0 0 0
    K_POINTS gamma
    """
    write(joinpath(d,"pw_scf.in"),scf)
    write(joinpath(d,"pp_ldos.in"),"&INPUTPP\n plot_num = 5\n sample_bias = -0.0220495933\n/\n")
    script="""
    #!/bin/bash
    #SBATCH --ntasks-per-node=8
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=96000MB
    QE_NTASKS=8
    srun -n "\$QE_NTASKS" pp.x -in pp_ldos.in
    """*commands
    write(joinpath(d,scriptname),script)
    scf,script
end
@testset "Existing QE preflight explicitly supports PP-only and rejects SCF" begin
    mktempdir() do d
        scf,script=fixture(d,"run_pp_only.sbatch","")
        rows=NamedTuple[]
        @test QEPreflight.check_dir!(rows,d;min_mem_mb=96000)==8
        @test only(r.value for r in rows if r.key=="mode")=="pp_only"
        for forbidden in ("pw.x -in pw_scf.in","pw.x -in pw_relax.in",
                "extract_qe_relaxed_xyz.jl","update_qe_positions_from_xyz.jl")
            write(joinpath(d,"run_pp_only.sbatch"),script*forbidden)
            @test_throws ErrorException QEPreflight.check_dir!(NamedTuple[],d)
        end
        write(joinpath(d,"run_pp_only.sbatch"),script)
        @test_throws ErrorException QEPreflight.check_dir!(NamedTuple[],d;min_mem_mb=100000)
        write(joinpath(d,"run_scf_pp.sbatch"),script*"pw.x -in pw_scf.in\n")
        @test QEPreflight.check_dir!(NamedTuple[],d)==8
        write(joinpath(d,"pw_relax.in"),scf)
        write(joinpath(d,"run_qe_mold.sbatch"),script*"pw.x -in pw_relax.in\nextract_qe_relaxed_xyz.jl\nupdate_qe_positions_from_xyz.jl\n")
        @test QEPreflight.check_dir!(NamedTuple[],d)==8
    end
    job=read(joinpath(ROOT,"hpc/qe_gamma_reconstruction.sbatch"),String)
    @test !occursin(r"\bpw\.x\b",job)
    @test occursin("#SBATCH --time=01:00:00",job) && occursin("#SBATCH --no-requeue",job)
    @test occursin("GIT_CONFIG_VALUE_0=never",job)
    @test occursin("cp -a",job) && !occursin("rm -",job)
    @test occursin("stock_control stock_ildos",job)
    @test occursin("patched_control patched_ildos repeat_ildos",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_gamma_reconstruction.sbatch")) $(joinpath(ROOT,"hpc/launch_qe_molds_remote.sh")) $(joinpath(ROOT,"hpc/submit_qe_molds.sh"))`)
    @test true
end
