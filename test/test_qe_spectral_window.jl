using Test, TOML, SHA, LinearAlgebra, Printf
include(joinpath(@__DIR__,"qe_spectral_window.jl"))
const S=QESpectralWindow
const ROOT=dirname(@__DIR__)
const CONFIG=joinpath(ROOT,"config/qe_spectral_window.toml")
const C=S.C

function xml_fixture(;values=[-.09,-.04,-.03,-.02,-.01,-.005,.001,.02,.04])
    """
    <qes:espresso Units="Hartree atomic units">
    <creator NAME="PWSCF" VERSION="7.4.1">PWSCF</creator>
    <output>
    <basis_set><gamma_only>true</gamma_only></basis_set>
    <convergence_info><scf_conv><convergence_achieved>true</convergence_achieved>
    <scf_error>1.0e-6</scf_error></scf_conv><wf_collected>true</wf_collected></convergence_info>
    <total_energy><etot>-10.0</etot></total_energy>
    <band_structure><nks>1</nks><nbnd>$(length(values))</nbnd>
    <lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit>
    <fermi_energy>0.0</fermi_energy><smearing degauss="0.01">mv</smearing>
    <ks_energies><k_point weight="2.0">0.0 0.0 0.0</k_point>
    <eigenvalues size="$(length(values))">$(join(values," "))</eigenvalues></ks_energies>
    </band_structure></output></qes:espresso>
    """
end

include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
@testset "Independent native-order verification reader" begin
    V=VerifyQESpectralWindow
    mktempdir() do d
        file=joinpath(d,"affine.cube")
        open(file,"w") do io
            println(io,"affine field\nlast axis fast\n0 1 2 3\n3 2 0 0\n4 0 3 0\n5 0 0 4")
            for i in 0:2, j in 0:3, k in 0:4
                println(io,100i+10j+k)
            end
        end
        c=V.cube(file)
        @test c.dims==[3,4,5]
        @test c.origin==[1.,2.,3.]
        @test c.grid[5,4,3]==234
        @test V.sample(c,[1.,2.,3.]+c.axes*[.4,1.2,2.3])≈54.3
        @test_throws AssertionError V.sample(c,[0.,0.,0.])
        @test V.element("<a><b>x</b></a>","b")=="x"
        @test_throws AssertionError V.element("<b>x</b><b>y</b>","b")
    end
end

@testset "Frozen spectral policy and explicit rejection" begin
    s=S.settings(CONFIG)
    @test s["model"]["scf_acceptance_ry"]==5e-5
    @test s["model"]["sample_bias_ry"]*S.RY_EV≈-.3 atol=1e-9
    @test !s["preprocessing"]["clip_negative_values"]
    mktempdir() do d
        for (section,key,value) in (("model","candidate_plot_num",3),("model","scf_acceptance_ry",1e-4),
            ("model","sample_bias_ry",-.04),("selection","expected_N",6),
            ("selection","require_archived_control_sha256",false),
            ("preprocessing","clip_negative_values",true),("preprocessing","half_nm",.5))
            bad=deepcopy(s); bad[section][key]=value
            p=joinpath(d,"bad.toml"); open(io->TOML.print(io,bad),p,"w")
            @test_throws ErrorException S.settings(p)
        end
    end
    @test S.main(["--help"])===nothing
    for opt in ("--truth","--expected-N","--composition","--height","--threshold","--manifest")
        @test_throws ErrorException S.main(["prepare-spectrum",opt,"forbidden"])
    end
    @test_throws ErrorException S.main(["unknown"])
    @test_throws ErrorException S.main(["compare","--control"])
    @test_throws ErrorException S.main(["compare","--control","a","--control","b"])
end

@testset "Hartree conversion, spectral weights and accepted checkpoint" begin
    s=S.settings(CONFIG)
    mktempdir() do d
        xml=joinpath(d,"input.xml"); write(xml,xml_fixture())
        original=S.sha(xml); sp=S.spectrum(xml,s)
        @test sp.accepted
        @test sp.scf_error_ry==2e-6
        @test sp.energy_ry==-20.0
        @test sp.fermi_ev==0.0
        @test sp.emin_ev≈s["model"]["sample_bias_ry"]*S.RY_EV
        @test sp.emax_ev==0.0
        @test sp.smearing_ev≈.02*S.RY_EV
        @test sp.first_band==2 && sp.last_band==8
        @test sp.rows[1].legacy_weight==sp.rows[9].legacy_weight==0.0
        @test [r.band for r in sp.rows if r.legacy_weight<0]==[2,3]
        @test [r.band for r in sp.rows if r.ildos_weight>0]==[5,6]
        @test all(r.ildos_weight>=0 for r in sp.rows)
        @test sum(r.ildos_weight for r in sp.rows)==4.0
        @test all(r.legacy_weight==2.0 for r in sp.rows[5:6])
        @test sp.rows[3].energy_ev≈-.03*S.HARTREE_EV
        @test S.cold_derivative(0.0)≈2exp(-.5)/sqrt(pi)
        @test abs(S.cold_derivative(sqrt(2.0)))<1e-15
        @test S.cold_derivative(2.0)<0 && S.cold_derivative(-2.0)>0
        for (old,new) in (("Hartree atomic units","Ry"),("VERSION=\"7.4.1\"","VERSION=\"7.5\""),
            ("<nks>1</nks>","<nks>2</nks>"),("<lsda>false","<lsda>true"),
            ("<gamma_only>true","<gamma_only>false"),("<nbnd>9","<nbnd>8"),
            ("<fermi_energy>0.0","<fermi_energy>NaN"),("degauss=\"0.01\"","degauss=\"0.0\""),
            ("<scf_error>1.0e-6","<scf_error>3.0e-5"),("<wf_collected>true","<wf_collected>false"),
            ("<convergence_achieved>true","<convergence_achieved>false"),("weight=\"2.0\"","weight=\"1.0\""),
            ("0.0 0.0 0.0","0.1 0.0 0.0"),("<fermi_energy>0.0</fermi_energy>",
                "<fermi_energy>0.0</fermi_energy><fermi_energy>0.0</fermi_energy>"))
            bad=joinpath(d,"bad.xml"); write(bad,replace(xml_fixture(),old=>new))
            @test_throws ErrorException S.spectrum(bad,s)
        end
        bad=joinpath(d,"boundary.xml"); write(bad,xml_fixture(values=[-.1,-.04,-.01,0.,.1]))
        @test_throws ErrorException S.spectrum(bad,s)
        bad=joinpath(d,"uncollected.xml"); write(bad,replace(xml_fixture(),"<wf_collected>true"=>"<wf_collected>false"))
        @test !S.spectrum(bad,s;require_accepted=false).accepted
        @test S.sha(xml)==original

        out=joinpath(d,"spectral")
        S.prepare_spectrum(xml,CONFIG,"glcn_central",joinpath(d,"checkpoint"),out)
        meta=TOML.parsefile(joinpath(out,"metadata.toml"))
        @test meta["legacy_negative_bands"]==2
        @test meta["ildos_bands"]==2 && meta["ildos_integral_expected"]==4.0
        @test meta["xml_sha256"]==original
        @test read(joinpath(out,"data-file-schema.xml"))==read(xml)
        @test read(joinpath(out,"settings.toml"))==read(CONFIG)
        control=read(joinpath(out,"control.in"),String); ildos=read(joinpath(out,"ildos.in"),String)
        @test occursin("plot_num = 5",control) && occursin("sample_bias",control) && !occursin("emin",control)
        @test occursin("plot_num = 10",ildos) && occursin("emin",ildos) && occursin("emax",ildos)
        @test !occursin("sample_bias",ildos) && !occursin("degauss",ildos)
        @test length(readlines(joinpath(out,"bands.tsv")))==10
        @test_throws ErrorException S.prepare_spectrum(xml,CONFIG,"glcn_central",d,out)
        symlink(out,joinpath(d,"link")); @test_throws ErrorException S.newdir(joinpath(d,"link"))
        for (prefix,checkpoint,stem,num) in (("a'","x","a",5),("a","x'","a",5),("a","x","a\nb",5),("a","x","a",3))
            @test_throws ErrorException S.pp_input(prefix,checkpoint,stem,num,sp,s)
        end
    end
end

function cube_fixture(path,value;dims=(5,5,5),step=.4)
    open(path,"w") do io
        println(io,"synthetic cube\nQE last-axis-fast")
        println(io,"0 0 0 0")
        for k in 1:3
            a=zeros(3); a[k]=step/C.BOHR_NM
            @printf(io,"%d %.17g %.17g %.17g\n",dims[k],a...)
        end
        for i in 1:prod(dims); @printf(io,"%.17g\n",value); end
    end
end

@testset "Independent cube integration, strict plane bounds and complete comparison" begin
    mktempdir() do d
        config=S.settings(CONFIG)
        v=.25; dims=(5,5,5); step=.4
        raw=joinpath(d,"control.cube"); cand=joinpath(d,"ildos.cube"); ref=joinpath(d,"ref.cube")
        cube_fixture(raw,-v); cube_fixture(cand,v); cp(raw,ref)
        c=S.QECubeMolds.read_qe_cube(cand)
        @test S.cube_stats(c).integral≈prod(dims)*v*(step/C.BOHR_NM)^3
        @test S.cube_stats(c).negative==0
        old=S.QECubeMolds.read_qe_cube(raw)
        @test S.cube_stats(old).negative==prod(dims)
        @test S.cube_stats(old).negative_abs_integral≈S.cube_stats(c).integral
        frame=C.Frame([.8,.8,.2],[1.,0.,0.],[0.,1.,0.])
        @test S.plane(c,frame,.5,config)≈fill(v,289)
        outside=C.Frame([0.,0.,0.],[1.,0.,0.],[0.,1.,0.])
        @test_throws ErrorException S.plane(c,outside,.5,config)
        f=joinpath(d,"frame.tsv")
        write(f,"key\tvalue\norigin_nm\t0.8,0.8,0.2\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
        @test C.read_frame(f).origin_nm==frame.origin_nm
        m=joinpath(d,"meta.toml")
        meta=Dict("config_sha256"=>S.sha(CONFIG),"ildos_integral_expected"=>S.cube_stats(c).integral)
        open(io->TOML.print(io,meta),m,"w")
        out=joinpath(d,"comparison")
        S.compare(raw,cand,ref,f,m,CONFIG,out)
        @test length(readlines(joinpath(out,"planes.tsv")))==1+3*289
        @test length(readlines(joinpath(out,"plane_summary.tsv")))==7
        @test length(readlines(joinpath(out,"cube_summary.tsv")))==3
        @test TOML.parsefile(joinpath(out,"comparison.toml"))["inputs"]["candidate"]["sha256"]==S.sha(cand)
        @test_throws ErrorException S.compare(raw,cand,ref,f,m,CONFIG,out)
        @test_throws ErrorException S.compare(raw,cand,cand,f,m,CONFIG,joinpath(d,"bad_control"))
        @test !ispath(joinpath(d,"bad_control"))
        meta["ildos_integral_expected"]*=2
        open(io->TOML.print(io,meta),m,"w")
        @test_throws ErrorException S.compare(raw,cand,ref,f,m,CONFIG,joinpath(d,"bad_integral"))
        @test !ispath(joinpath(d,"bad_integral"))
    end
end

@testset "QE launcher keeps custom Julia and explicit Slurm export" begin
    launch=read(joinpath(ROOT,"hpc/launch_qe_molds_remote.sh"),String)
    submit=read(joinpath(ROOT,"hpc/submit_qe_molds.sh"),String)
    job=read(joinpath(ROOT,"hpc/qe_spectral_window.sbatch"),String)
    @test occursin("JULIA_BIN=",submit) && occursin("--export=ALL,JULIA_BIN=",submit)
    @test occursin("'\$JULIA_BIN' --project",launch)
    @test occursin("#SBATCH --time=01:00:00",job) && occursin("#SBATCH --no-requeue",job)
    @test occursin("cp -a",job) && !occursin("rm -",job) && !occursin("pw_relax",job)
    @test occursin("cmp control.cube",job) && occursin("cmp glcnac_central_ldos_accept_plain.cube",job)
    run(`bash -n $(joinpath(ROOT,"hpc/launch_qe_molds_remote.sh")) $(joinpath(ROOT,"hpc/submit_qe_molds.sh")) $(joinpath(ROOT,"hpc/qe_spectral_window.sbatch"))`)
    @test true
end
