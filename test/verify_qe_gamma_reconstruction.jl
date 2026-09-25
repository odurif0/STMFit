#!/usr/bin/env julia
# Independent saved-output check: native-order reader, no Gamma experiment helpers.
module VerifyQEGamma
using Test, TOML, LinearAlgebra, Statistics, Printf
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow
function verify_case(case,config)
    s=TOML.parsefile(config); meta=TOML.parsefile(joinpath(case,"metadata.toml"))
    @test meta["config_sha256"]==V.digest(config)
    xml=read(joinpath(case,"data-file-schema.xml"),String)
    output=V.element(xml,"output"); bands=V.element(output,"band_structure")
    convergence=V.element(output,"convergence_info")
    scf=V.element(convergence,"scf_conv")
    @test strip(V.element(convergence,"wf_collected"))=="true"
    @test strip(V.element(scf,"convergence_achieved"))=="true"
    @test 0<=2V.value(V.element(scf,"scf_error"))<=s["model"]["scf_acceptance_ry"]
    eigen=V.value.(split(match(r"<eigenvalues[^>]*>(.*?)</eigenvalues>"s,bands).captures[1]))
    ef=V.value(V.element(bands,"fermi_energy"))
    lo=ef+s["model"]["sample_bias_ry"]/2
    expected=2count(e->lo<=e<=ef,eigen)
    @test expected==meta["ildos_integral_expected"]
    @test meta["xml_sha256"]==V.digest(joinpath(case,"data-file-schema.xml"))
    hartree_ev=4.3597447222071e-18/1.602176634e-19
    for name in ("stock_ildos","patched_ildos","repeat_ildos")
        input=read(joinpath(case,name*".in"),String)
        @test occursin("plot_num = 10",input)
        for (key,target) in (("emin",lo*hartree_ev),("emax",ef*hartree_ev))
            @test V.value(match(Regex(key*"\\s*=\\s*([^\\s]+)"),input).captures[1])≈target rtol=1e-14
        end
    end
    result=TOML.parsefile(joinpath(case,"comparison/summary.toml"))
    @test all(values(result["checks"]))
    @test V.digest(joinpath(case,"site_control.cube"))==V.digest(joinpath(case,"reference.cube"))
    @test V.digest(joinpath(case,"stock_control.cube"))==V.digest(joinpath(case,"patched_control.cube"))
    @test V.digest(joinpath(case,"patched_ildos.cube"))==V.digest(joinpath(case,"repeat_ildos.cube"))
    rows=V.tsv(joinpath(case,"comparison/cube_summary.tsv"))
    pixels=V.tsv(joinpath(case,"comparison/planes.tsv"))
    f=Dict(r["key"]=>V.value.(split(r["value"],',')) for r in V.tsv(joinpath(case,"frame.tsv")))
    normal=cross(f["t_axis"],f["u_axis"]); normal/=norm(normal)
    grid=collect(-8:8)*s["preprocessing"]["step_nm"]
    @test first(grid)==-s["preprocessing"]["half_nm"]
    @test length(rows)==6 && length(pixels)==6*3*289
    for row in rows
        name=row["observable"]; file=joinpath(case,name*".cube"); c=V.cube(file)
        integral=sum(c.raw)*abs(det(c.axes))
        @test result["cube_sha256"][name]==V.digest(file)
        @test parse(Int,row["values"])==length(c.raw)
        @test parse(Int,row["negative"])==count(<(0),c.raw)
        @test V.value(row["minimum"])==minimum(c.raw)
        @test V.value(row["maximum"])==maximum(c.raw)
        @test V.value(row["integral"])≈integral rtol=1e-12
        if occursin("ildos",name)
            @test abs(integral-expected)<=s["selection"]["ildos_integral_rtol"]*expected
        end
        for h in s["preprocessing"]["diagnostic_heights_nm"]
            plane=filter(r->r["observable"]==name && V.value(r["height_nm"])==h,pixels)
            @test length(plane)==289
            for (i,(u,t)) in enumerate((u,t) for u in grid for t in grid)
                point=(f["origin_nm"]+t*f["t_axis"]+u*f["u_axis"]+h*normal)/0.05291772109
                @test parse(Int,plane[i]["pixel"])==i
                @test V.value(plane[i]["value"])≈V.sample(c,point) rtol=1e-10 atol=1e-15
            end
        end
        @printf("%s %s integral=%.12g negative=%d min=%.9g max=%.9g\n",
            basename(case),name,integral,count(<(0),c.raw),minimum(c.raw),maximum(c.raw))
    end
    stock=V.cube(joinpath(case,"stock_ildos.cube")).raw
    fixed=V.cube(joinpath(case,"patched_ildos.cube")).raw
    digits=s["selection"]["cube_significant_digits"]
    rounding(x)=iszero(x) ? 0.0 : .5*exp10(floor(log10(abs(x)))+1-digits)
    average=sum(fixed)/length(fixed)
    average_error=sum(rounding(x) for x in fixed)/length(fixed)
    errors=0
    for i in eachindex(stock,fixed)
        allowed=2rounding(stock[i])+rounding(fixed[i])+average_error+2s["selection"]["numeric_roundoff_atol"]
        errors+=abs(2stock[i]-fixed[i]-average)>allowed
    end
    @test errors==0
    reference=V.cube(joinpath(case,"reference.cube")).raw
    rebuilt=V.cube(joinpath(case,"stock_control.cube")).raw
    control_errors=count(eachindex(reference,rebuilt)) do i
        abs(reference[i]-rebuilt[i])>rounding(reference[i])+rounding(rebuilt[i])+s["selection"]["numeric_roundoff_atol"]
    end
    @test control_errors==0
    println(basename(case),": full-volume conjugate-pair relation independently verified")
end
function main(args=ARGS)
    @assert VERSION.major==1 && VERSION.minor==13
    args==["--help"] && return println("verify_qe_gamma_reconstruction.jl RUN_DIR")
    length(args)==1 || error("Expected saved run directory")
    run=only(args)
    @testset "Independent Gamma output verification" begin
        for m in ("glcn","glcnac")
            verify_case(joinpath(run,m),joinpath(run,"settings.toml"))
        end
    end
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifyQEGamma.main()
