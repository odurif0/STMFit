using Test, TOML, Statistics, Printf
include(joinpath(@__DIR__,"qe_paw_decomposition.jl"))
include(joinpath(@__DIR__,"verify_qe_paw_decomposition.jl"))
const D=QEPAWDecomposition
const ROOT=dirname(@__DIR__)
const CONFIG=joinpath(ROOT,"config/qe_paw_decomposition.toml")

@testset "Frozen component policy and rounding-only additivity" begin
    s=D.settings(CONFIG)
    @test !s["preprocessing"]["normalize_components"]
    @test D.main(["--help"])===nothing
    @test_throws ErrorException D.main(["--truth","NKNNKN"])
    mktempdir() do dir
        for (section,key,value) in (("model","sample_bias_ry",-.03),("model","scf_acceptance_ry",1e-4),
                ("preprocessing","half_nm",.64),("preprocessing","normalize_components",true),
                ("preprocessing","clip_negative_values",true),("selection","expected_N",6),
                ("selection","cube_significant_digits",4),("selection","numeric_roundoff_atol",1e-6),
                ("selection","require_corrected_total_sha256",false))
            bad=deepcopy(s); bad[section][key]=value; file=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,bad),file,"w")
            @test_throws ErrorException D.settings(file)
        end
    end
    smooth=[2+.4sin(2pi*k/17) for k in 0:16]
    aug=[3cos(2pi*k/17) for k in 0:16]; total=smooth+aug
    @test any(<(0),total)
    @test D.additivity(total,smooth,aug,s).violations==0
    @test D.additivity(round.(total;sigdigits=5),round.(smooth;sigdigits=5),round.(aug;sigdigits=5),s).violations==0
    @test D.additivity(total,smooth,zero.(aug),s).violations>0
    @test D.additivity(total.+.01,smooth,aug,s).violations==17
    @test_throws ErrorException D.additivity([1.],smooth,aug,s)
    @test_throws ErrorException D.additivity([NaN],[0.],[0.],s)
end

@testset "Signed component fixture and independent complete check" begin
    mktempdir() do dir
        case_dir=joinpath(dir,"synthetic"); mkdir(case_dir)
        s=D.settings(CONFIG)
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
        xmlfile=joinpath(case_dir,"data-file-schema.xml"); write(xmlfile,xml)
        sp=D.S.spectrum(xmlfile,s)
        write(joinpath(case_dir,"frame.tsv"),"key\tvalue\norigin_nm\t0.5,0.5,0.1\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
        step=.5/D.G.C.BOHR_NM; average=4/(27step^3)
        smooth=[average*(1+.1cos(2pi*k/3)+.2sin(2pi*i/3)) for i in 0:2 for j in 0:2 for k in 0:2]
        aug=[2.4average*cos(2pi*k/3) for i in 0:2 for j in 0:2 for k in 0:2]
        full=smooth+aug
        for name in (D.NAMES...,"reference","reference_total")
            values=occursin("smooth",name) ? smooth : occursin("augmentation",name) ? aug : full
            (endswith(name,"control") || name=="reference") && (values=-full)
            open(joinpath(case_dir,name*".cube"),"w") do io
                println(io,"synthetic cube\nQE last axis fast\n0 0 0 0")
                for k in 1:3
                    axis=zeros(3); axis[k]=step
                    @printf(io,"3 %.17g %.17g %.17g\n",axis...)
                end
                for value in values; @printf(io,"%13.4E\n",value); end
            end
            name in D.NAMES && write(joinpath(case_dir,name*".in"),D.S.pp_input("synthetic","checkpoint",name,endswith(name,"control") ? 5 : 10,sp,s))
        end
        meta=Dict("config_sha256"=>D.S.sha(CONFIG),"xml_sha256"=>D.S.sha(xmlfile),"ildos_integral_expected"=>4.,
            "reference_total_sha256"=>D.S.sha(joinpath(case_dir,"reference_total.cube")),
            "reference_stm_sha256"=>D.S.sha(joinpath(case_dir,"reference.cube")))
        open(io->TOML.print(io,meta),joinpath(case_dir,"metadata.toml"),"w")
        @test D.compare(case_dir,CONFIG)===nothing
        @test all(values(TOML.parsefile(joinpath(case_dir,"comparison/summary.toml"))["checks"]))
        VerifyQEPAW.verify_case(case_dir,CONFIG)
        attribution=VerifyQEPAW.V.tsv(joinpath(case_dir,"comparison/attribution.tsv"))
        @test any(parse(Int,r["total_negative"])>0 for r in attribution)
        @test all(parse(Int,r["smooth_negative"])==0 for r in attribution)
        @test all(r["total_negative"]==r["negative_total_with_negative_augmentation"] for r in attribution)
        @test_throws ErrorException D.compare(case_dir,CONFIG)
    end
end

@testset "Minimal reversible source variants and PP-only one-job boundary" begin
    job=read(joinpath(ROOT,"hpc/qe_paw_decomposition.sbatch"),String)
    @test !occursin(r"\bpw\.x\b",job)
    @test occursin("#SBATCH --time=01:00:00",job) && occursin("#SBATCH --no-requeue",job)
    @test occursin("module load",job) && occursin("mkl/2024.0",job)
    @test occursin("GIT_CONFIG_VALUE_0=never",job)
    @test occursin("cmp full_ildos.cube reference_total.cube",job)
    @test !occursin("rm -",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_paw_decomposition.sbatch"))`)
    @test true
    # The actual source is an ignored, pinned download; use a tiny owned copy.
    source=joinpath(ROOT,"qe/gamma_source_20260925/PP/src/local_dos.f90")
    @test D.S.sha(source)=="003d5f92756420115d9393dd849e3e9270c50317abb3564cd2d432ab901da4c9"
    mktempdir() do dir
        mkpath(joinpath(dir,"PP/src")); file=joinpath(dir,"PP/src/local_dos.f90"); cp(source,file)
        run(Cmd(`git apply $(joinpath(ROOT,"hpc/patches/qe-7.4.1-gamma-density.patch"))`;dir=dir))
        full=D.S.sha(file)
        @test full=="8c803e855ebf9bbf5b15a8b4218815f0bfca73d701925d6cc556d593def96263"
        for component in ("smooth","augmentation")
            patch=joinpath(ROOT,"hpc/patches/qe-7.4.1-ildos-$component.patch")
            run(Cmd(`git apply --check $patch`;dir=dir)); run(Cmd(`git apply $patch`;dir=dir))
            @test D.S.sha(file)!=full
            code=read(file,String)
            @test occursin("IF (gamma_only) psic(dfftp%nlm(:))",code)
            @test occursin(component=="smooth" ? "IF (iflag /= 3) CALL addusdens" : "IF (iflag == 3) rho%of_g",code)
            run(Cmd(`git apply -R $patch`;dir=dir))
            @test D.S.sha(file)==full
        end
    end
end
