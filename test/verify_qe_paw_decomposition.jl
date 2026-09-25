#!/usr/bin/env julia
# Independent native-order component check, not the experiment's reader/helpers.
module VerifyQEPAW
using Test, TOML, Statistics, LinearAlgebra, Printf
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow

function verify_case(case_dir,config)
    s=TOML.parsefile(config); meta=TOML.parsefile(joinpath(case_dir,"metadata.toml"))
    @test !s["preprocessing"]["clip_negative_values"] && !s["preprocessing"]["normalize_components"]
    @test meta["config_sha256"]==V.digest(config)
    xml=read(joinpath(case_dir,"data-file-schema.xml"),String)
    output=V.element(xml,"output"); bands=V.element(output,"band_structure")
    convergence=V.element(output,"convergence_info"); scf=V.element(convergence,"scf_conv")
    @test strip(V.element(convergence,"wf_collected"))=="true"
    @test strip(V.element(scf,"convergence_achieved"))=="true"
    @test 0<=2V.value(V.element(scf,"scf_error"))<=s["model"]["scf_acceptance_ry"]
    @test strip(V.element(V.element(output,"basis_set"),"gamma_only"))=="true"
    @test strip(V.element(bands,"lsda"))=="false" && strip(V.element(bands,"noncolin"))=="false"
    eigen=V.value.(split(match(r"<eigenvalues[^>]*>(.*?)</eigenvalues>"s,bands).captures[1]))
    ef=V.value(V.element(bands,"fermi_energy")); lo=ef+s["model"]["sample_bias_ry"]/2
    expected=2count(e->lo<=e<=ef,eigen)
    @test expected==meta["ildos_integral_expected"]
    @test meta["xml_sha256"]==V.digest(joinpath(case_dir,"data-file-schema.xml"))
    result=TOML.parsefile(joinpath(case_dir,"comparison/summary.toml"))
    @test all(values(result["checks"]))
    @test meta["reference_total_sha256"]==V.digest(joinpath(case_dir,"reference_total.cube"))==V.digest(joinpath(case_dir,"full_ildos.cube"))
    @test meta["reference_stm_sha256"]==V.digest(joinpath(case_dir,"reference.cube"))
    for component in ("smooth","augmentation")
        @test V.digest(joinpath(case_dir,component*"_ildos.cube"))==V.digest(joinpath(case_dir,"repeat_"*component*".cube"))
    end
    rows=V.tsv(joinpath(case_dir,"comparison/cube_summary.tsv"))
    pixels=V.tsv(joinpath(case_dir,"comparison/planes.tsv"))
    f=Dict(r["key"]=>V.value.(split(r["value"],',')) for r in V.tsv(joinpath(case_dir,"frame.tsv")))
    normal=cross(f["t_axis"],f["u_axis"]); normal/=norm(normal)
    grid=collect(-8:8)*s["preprocessing"]["step_nm"]
    @test first(grid)==-s["preprocessing"]["half_nm"]
    @test length(rows)==8 && length(pixels)==8*3*289
    hartree_ev=4.3597447222071e-18/1.602176634e-19
    for row in rows
        name=row["observable"]; file=joinpath(case_dir,name*".cube"); c=V.cube(file)
        integral=sum(c.raw)*abs(det(c.axes))
        @test V.digest(file)==result["cube_sha256"][name]
        @test parse(Int,row["values"])==length(c.raw)
        @test parse(Int,row["negative"])==count(<(0),c.raw)
        @test V.value(row["minimum"])==minimum(c.raw)
        @test V.value(row["maximum"])==maximum(c.raw)
        @test V.value(row["integral"])≈integral rtol=1e-12
        if endswith(name,"control")
            @test V.digest(file)==meta["reference_stm_sha256"]
        else
            input=read(joinpath(case_dir,name*".in"),String)
            @test occursin("plot_num = 10",input)
            for (key,target) in (("emin",lo*hartree_ev),("emax",ef*hartree_ev))
                @test V.value(match(Regex(key*"\\s*=\\s*([^\\s]+)"),input).captures[1])≈target rtol=1e-14
            end
        end
        name=="full_ildos" && @test abs(integral-expected)<=s["selection"]["ildos_integral_rtol"]*expected
        for height in s["preprocessing"]["diagnostic_heights_nm"]
            plane=filter(r->r["observable"]==name && V.value(r["height_nm"])==height,pixels)
            @test length(plane)==289
            for (i,(u,t)) in enumerate((u,t) for u in grid for t in grid)
                point=(f["origin_nm"]+t*f["t_axis"]+u*f["u_axis"]+height*normal)/0.05291772109
                @test parse(Int,plane[i]["pixel"])==i
                @test V.value(plane[i]["value"])≈V.sample(c,point) rtol=1e-10 atol=1e-15
            end
        end
        @printf("%s %s integral=%.12g negatives=%d min=%.12g\n",basename(case_dir),name,integral,count(<(0),c.raw),minimum(c.raw))
    end
    cubes=[V.cube(joinpath(case_dir,n*".cube")) for n in ("full_ildos","smooth_ildos","augmentation_ildos")]
    total,smooth,aug=(c.raw for c in cubes)
    @test all(c.dims==cubes[1].dims && c.origin==cubes[1].origin && c.axes==cubes[1].axes for c in cubes)
    digits=s["selection"]["cube_significant_digits"]; atol=s["selection"]["numeric_roundoff_atol"]
    rounding(x)=iszero(x) ? 0. : .5*exp10(floor(log10(abs(x)))+1-digits)
    violations=count(eachindex(total,smooth,aug)) do i
        abs(total[i]-(smooth[i]+aug[i]))>rounding(total[i])+rounding(smooth[i])+rounding(aug[i])+atol
    end
    @test violations==0
    for row in V.tsv(joinpath(case_dir,"comparison/attribution.tsv"))
        height=V.value(row["height_nm"])
        pv=[V.value.([r["value"] for r in pixels if r["observable"]==name && V.value(r["height_nm"])==height])
            for name in ("full_ildos","smooth_ildos","augmentation_ildos")]
        t,p,a=pv
        @test parse(Int,row["total_negative"])==count(<(0),t)
        @test parse(Int,row["smooth_negative"])==count(<(0),p)
        @test parse(Int,row["augmentation_negative"])==count(<(0),a)
        @test parse(Int,row["negative_total_with_nonnegative_smooth"])==count(i->t[i]<0 && p[i]>=0,eachindex(t))
        @test parse(Int,row["negative_total_with_negative_augmentation"])==count(i->t[i]<0 && a[i]<0,eachindex(t))
    end
    println(basename(case_dir),": native component additivity independently verified")
end

function main(args=ARGS)
    @assert VERSION.major==1 && VERSION.minor==13
    args==["--help"] && return println("verify_qe_paw_decomposition.jl RUN_DIR")
    length(args)==1 || error("Expected saved component run")
    run_dir=only(args)
    @testset "Independent PAW component verification" begin
        for m in ("glcn","glcnac")
            verify_case(joinpath(run_dir,m),joinpath(run_dir,"settings.toml"))
        end
    end
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifyQEPAW.main()
