using Test, LinearAlgebra, Printf, TOML, SHA
module Legacy
include(joinpath(@__DIR__,"test_cc_mold_native.jl"))
end
include(joinpath(@__DIR__,"build_qe_order_molds.jl"))
const D=BuildQEOrderMolds
const Q=D.QECubeMolds
const C=Q.CCMoldNative
const ROOT=dirname(@__DIR__)
const POLICY=joinpath(ROOT,"config/qe_cube_order.toml")
const CONFIG=joinpath(ROOT,"config/unit_assignment_patch_support.toml")

function last_fast_fixture(path,n,origin,axes,density)
    open(path,"w") do io
        println(io,"Synthetic QE-order cube\nthird coordinate written fastest")
        @printf(io,"0 %.17g %.17g %.17g\n",(origin./C.BOHR_NM)...)
        for a in 1:3
            @printf(io,"%d %.17g %.17g %.17g\n",n[a],(axes[:,a]./C.BOHR_NM)...)
        end
        for i in 0:n[1]-1,j in 0:n[2]-1,k in 0:n[3]-1
            @printf(io,"%.17g\n",density(i,j,k))
        end
    end
    path
end

@testset "QE order, noncubic affine field and unchanged legacy parser" begin
    mktempdir() do dir
        n=(4,3,5); origin=[-.1,.2,-.3]
        A=[.13 .01 .02; .02 .17 .01; 0 .03 .19]
        field(i,j,k)=1+i+10j+100k
        path=last_fast_fixture(joinpath(dir,"qe.cube"),n,origin,A,field)
        digest=sha256(read(path)); raw=C.read_cube(path); qe=Q.read_qe_cube(path)
        @test raw.values[1:5]==[1,101,201,301,401]
        @test qe.values[1:5]==[1,2,3,4,11]
        @test raw.dims==qe.dims==n && raw.origin_nm==qe.origin_nm && raw.axes_nm==qe.axes_nm
        @test sha256(read(path))==digest
        coords=[Float64[i,j,k] for i in 0:3,j in 0:2,k in 0:4]
        for c in coords
            points=reshape(origin+A*c,1,3)
            @test only(C.sample_volume(qe,points))≈field(c...) atol=1e-10
        end
        for c in ([.25,.5,.75],[1.1,1.2,2.3],[2.5,.7,3.2])
            points=reshape(origin+A*c,1,3)
            @test only(C.sample_volume(qe,points))≈field(c...) atol=1e-10
            @test !isapprox(only(C.sample_volume(raw,points)),field(c...);atol=1e-3)
        end
        # Preserve the historical partial-boundary weighting, not a second fix.
        p=reshape(qe.origin_nm+qe.axes_nm*[-.5,0.,0.],1,3)
        @test only(C.sample_volume(qe,p))≈.5 atol=1e-10
        p=reshape(qe.origin_nm+qe.axes_nm*[-1.1,0.,0.],1,3)
        @test isnan(only(C.sample_volume(qe,p)))
        @test C.read_cube(path).values==raw.values
        @test Q.settings(POLICY)["model"]["cube_order"]=="qe_last_axis_fast"
        for (table,key,value) in (("model","cube_order","legacy"),
            ("selection","surface_policy","nearest_target"),("preprocessing","cube_units","angstrom"))
            s=TOML.parsefile(POLICY); s[table][key]=value
            bad=joinpath(dir,"bad.toml"); open(io->TOML.print(io,s),bad,"w")
            @test_throws ArgumentError Q.settings(bad)
        end
    end
end

@testset "Whole-surface equivalence, provenance and isolated builder CLI" begin
    mktempdir() do dir
        c0,c1,f0,f1=Legacy.synthetic_pair(dir)
        newpaths=String[]
        for path in (c0,c1)
            raw=C.read_cube(path); n1,n2,n3=raw.dims
            lookup=(i,j,k)->raw.values[1+i+n1*j+n1*n2*k]
            push!(newpaths,last_fast_fixture(path*".qe",raw.dims,raw.origin_nm,raw.axes_nm,lookup))
        end
        # Legacy fixture in first-fast order and QE fixture in last-fast order
        # describe the same analytic field; every downstream operation is literal.
        expected=joinpath(dir,"expected.tsv")
        old=C.build_molds(cube0=c0,cube1=c1,frame0=f0,frame1=f1,config=CONFIG,out=expected)
        out=joinpath(dir,"new")
        args=["--cube0",newpaths[1],"--cube1",newpaths[2],"--frame0",f0,"--frame1",f1,
            "--config",CONFIG,"--settings",POLICY,"--outdir",out]
        D.main(vcat(args,["--dry-run"]))
        @test !ispath(out)
        now=D.main(args)
        @test read(joinpath(out,"templates_cc.tsv"))==read(expected)
        for (a,b) in zip((old.type0,old.type1),(now.type0,now.type1))
            for k in (:iso,:mean_height,:nvalid,:heights,:normalized,:variants)
                @test isequal(getproperty(a,k),getproperty(b,k))
            end
        end
        p=TOML.parsefile(joinpath(out,"provenance.toml"))
        @test p["cube_order"]=="qe_last_axis_fast"
        @test p["template_sha256"]==bytes2hex(sha256(read(expected)))
        @test length(p["inputs"])==6
        for input in values(p["inputs"]); @test input["sha256"]==bytes2hex(sha256(read(input["path"]))); end
        @test length(readlines(joinpath(out,"surface_audit.tsv")))==579
        @test_throws ErrorException D.main(args)
        link=joinpath(dir,"link"); symlink(out,link)
        @test_throws ErrorException D.main(replace.(args,out=>link))
        for extra in (["--truth","labels.tsv"],["--cube0",c0],["--settings"])
            @test_throws ErrorException D.main(vcat(replace.(args,out=>joinpath(dir,"unused")),extra))
        end
        badframe=joinpath(dir,"badframe.tsv"); write(badframe,"")
        failed=joinpath(dir,"failed")
        @test_throws ArgumentError D.main(replace.(replace.(args,out=>failed),f1=>badframe))
        @test isfile(joinpath(failed,"failure.txt")) && !isfile(joinpath(failed,"templates_cc.tsv"))
        cli=joinpath(ROOT,"test/build_qe_order_molds.jl")
        cliout=joinpath(dir,"cli")
        cliargs=replace.(args,out=>cliout)
        run(pipeline(`$(Base.julia_cmd()) --project=$ROOT $cli $cliargs`;stdout=devnull))
        @test read(joinpath(cliout,"templates_cc.tsv"))==read(expected)
    end
end
