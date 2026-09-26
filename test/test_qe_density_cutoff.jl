#!/usr/bin/env julia
using Test, TOML, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"qe_density_cutoff.jl"))
const C=QEDensityCutoff
const R=C.R; const S=C.S; const N=C.N
const ROOT=normpath(joinpath(@__DIR__,".."))
const CONFIG=joinpath(ROOT,"config/qe_density_cutoff.toml")

@testset "Cutoff scope and all-plane descriptive statistics" begin
    s=C.settings(CONFIG)
    mktempdir() do dir
        for (section,key,value) in (("model","candidate_ecutrho_ry",1080.),("model","ecutwfc_ry",60.),
                ("model","scf_acceptance_ry",1e-4),("model","startingpot","file"),
                ("model","scf_max_seconds",14000),("selection","adopt_reference",true),
                ("preprocessing","clip_density",true),("preprocessing","cross_grid_interpolation",true))
            bad=deepcopy(s); bad[section][key]=value
            file=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),file,"w")
            @test_throws ErrorException C.settings(file)
        end
    end
    rho=reshape([-2.,-1.,1.,2.,-4.,0.,2.,6.,-1.,-1.,-1.,-1.,1.,1.,1.,1.],2,2,4)
    total=2rho; electro=rho
    g=(;geo=(;geometry=(;cell_nm=[2.,2.,4.]),dims=[2,2,4]),lower=.5,upper=3.5,midpoint=2.)
    rows=C.plane_rows(rho,total,electro,g,"glcn",720.)
    @test length(rows)==4 && [r.k for r in rows]==collect(0:3)
    @test [r.z_nm for r in rows]==[0.,1.,2.,3.]
    @test [r.paw_gap for r in rows]==[false,true,true,true]
    @test [r.lower_half for r in rows]==[false,true,false,false]
    @test [r.negative_density_fraction for r in rows]==[.5,.25,1.,0.]
    @test rows[1].mean_negative_density==.75 && rows[1].mean_density==0
    @test rows[1].std_density==sqrt(2.5) && rows[1].std_xc_ry==sqrt(2.5)
    @test rows[3].std_density==0 && rows[3].min_density==rows[3].max_density==-1
    @test_throws ErrorException C.plane_rows(rho[:,:,1:3],total,electro,g,"glcn",720.)
    bad=copy(rho); bad[1]=NaN
    @test_throws ErrorException C.plane_rows(bad,total,electro,g,"glcn",720.)
    @test occursin("plot_num=0",C.pp_input("glcn","density"))
    @test occursin("plot_num=11",C.pp_input("glcnac","electrostatic"))
    @test_throws ErrorException C.pp_input("glcn","ldos")
end

function metadata_tests(root)
    s=C.settings(CONFIG); old=joinpath(root,s["model"]["baseline_run"])
    @testset "Accepted metadata and fixed-geometry SCF inputs" begin
        for (j,m) in enumerate(R.MOLECULES)
            xml=joinpath(old,m,"data-file-schema.xml"); a=C.state(xml,s,360.)
            input=C.input_text(xml,m,s)
            @test S.sha(xml)==R.D.G.XML_SHA[j]
            @test occursin("ecutrho=720.0",input) && occursin("ecutwfc=50.0",input)
            @test occursin("conv_thr="*(j==1 ? "1.0e-7" : "5.0e-5"),input) && occursin("startingpot='atomic'",input)
            @test occursin("startingwfc='file'",input) && occursin("restart_mode='from_scratch'",input)
            @test !occursin("nr1=",input) && !occursin("calculation='relax'",input)
            @test occursin("mixing_ndim="*(j==1 ? "8" : "20"),input)
            # Independent line/card parser checks every coordinate, not a fixture snapshot.
            lines=split(input,'\n'); ia=findfirst(==("ATOMIC_POSITIONS bohr"),lines)
            ic=findfirst(==("CELL_PARAMETERS bohr"),lines)
            for k in 1:3
                @test parse.(Float64,split(lines[ic+k]))==a.geo.cell[:,k]
            end
            atoms=a.geo.geometry.atoms
            @test findfirst(==("K_POINTS gamma"),lines)-ia-1==length(atoms)
            for (k,atom) in enumerate(atoms)
                cols=split(lines[ia+k])
                @test R.D.P.Z[cols[1]]==atom.z
                @test parse.(Float64,cols[2:4])*R.D.G.C.BOHR_NM==atom.position_nm
            end
            mktempdir() do dir
                open(io->TOML.print(io,s),joinpath(dir,"settings.toml"),"w")
                case=joinpath(dir,m); mkpath(joinpath(case,"work",m*"_central.save"))
                cp(realpath(xml),joinpath(case,"accepted.xml"))
                candidate=replace(read(xml,String),"<ecutrho>1.800000000000000E+002</ecutrho>"=>"<ecutrho>3.600000000000000E+002</ecutrho>",
                    "nr1=\"240\" nr2=\"180\" nr3=\"250\""=>"nr1=\"336\" nr2=\"256\" nr3=\"360\"")
                file=joinpath(case,"work",m*"_central.save/data-file-schema.xml")
                write(file,candidate)
                write(joinpath(case,"pw_scf.out"),"Starting wfcs from file\nInitial potential from superposition of free atoms\nconvergence has been achieved\nJOB DONE.\n")
                @test C.check_scf(dir,m)===nothing
                @test read(joinpath(case,"data-file-schema.xml"),String)==candidate
                write(file,replace(candidate,"<convergence_achieved>true</convergence_achieved>"=>"<convergence_achieved>false</convergence_achieved>"))
                @test_throws ErrorException C.check_scf(dir,m)
                write(file,replace(candidate,"<functional>PBE</functional>"=>"<functional>PBE0</functional>"))
                @test_throws ErrorException C.check_scf(dir,m)
            end
        end
    end
end

function verify_saved(root,run,m)
    s=C.settings(joinpath(run,"settings.toml")); dir=joinpath(run,m)
    rows=R.V.tsv(joinpath(dir,"analysis/planes.tsv"))
    grids=R.V.tsv(joinpath(dir,"analysis/grids.tsv"))
    report=TOML.parsefile(joinpath(dir,"analysis/summary.toml"))
    old=joinpath(root,s["model"]["baseline_run"])
    p=R.settings(joinpath(root,"config/qe_vacuum_potential.toml"))["preprocessing"]
    num(r,k)=parse(Float64,r[k]); int(r,k)=parse(Int,r[k])
    @testset "Independent cube summaries for both density cutoffs" begin
        @test report["exact_total_repeat"] && !report["cross_grid_interpolation"] && !report["potential_adopted"]
        @test report["config_sha256"]==S.sha(joinpath(run,"settings.toml"))
        @test report["xml_sha256"]==S.sha(joinpath(dir,"data-file-schema.xml"))
        for ext in ("dat","cube")
            @test S.sha(joinpath(dir,"total.$ext"))==S.sha(joinpath(dir,"repeat_total.$ext"))
        end
        for (name,cutoff,source,geometry_run) in (("baseline",360.,joinpath(old,m,"stock"),old),("candidate",720.,dir,run))
            cubes=Dict(stem=>R.V.cube(joinpath(source,stem*".cube")) for stem in C.STEMS)
            g=R.case_geometry(geometry_run,m); selected=filter(r->num(r,"ecutrho_ry")==cutoff,rows)
            gg=only(filter(r->r["source"]==name,grids)); dims=cubes["density"].dims
            @test all(c.dims==dims for c in values(cubes)) && dims==g.geo.dims
            @test length(selected)==dims[3] && int.(selected,"k")==collect(0:dims[3]-1)
            @test [int(gg,k) for k in ("nx","ny","nz")]==dims
            for stem in C.STEMS, ext in ("dat","cube")
                @test S.sha(joinpath(source,stem*"."*ext))==report["export_sha256"][name*"/"*stem*"."*ext]
            end
            for (k,row) in enumerate(selected)
                @test row["molecule"]==m && int(row,"samples")==dims[1]*dims[2]
                z=(k-1)*g.geo.geometry.cell_nm[3]/dims[3]
                @test num(row,"z_nm")==z
                @test parse(Bool,row["paw_gap"])==(g.lower<z<g.upper)
                @test parse(Bool,row["lower_half"])==(g.lower<z<g.midpoint)
                arrays=Dict(stem=>vec(collect(view(cubes[stem].grid,k,:,:))) for stem in C.STEMS)
                bounds=Dict(stem=>maximum(N.halfquantum(v,p["cube_significant_digits"])+
                    10N.halfquantum(v,p["native_significant_digits"])+p["numeric_roundoff_atol_ry"] for v in arrays[stem]) for stem in C.STEMS)
                arrays["xc"]=arrays["total"].-arrays["electrostatic"]
                bounds["xc"]=bounds["total"]+bounds["electrostatic"]
                for stem in (C.STEMS...,"xc")
                    v=arrays[stem]; mu=sum(v)/length(v); sigma=sqrt(sum((x-mu)^2 for x in v)/length(v))
                    suffix=stem=="density" ? "_density" : "_"*stem*"_ry"
                    for (key,value) in (("min",minimum(v)),("mean",mu),("max",maximum(v)),("std",sigma))
                        @test abs(num(row,key*suffix)-value)<=bounds[stem]
                    end
                end
                rho=arrays["density"]; b=bounds["density"]
                @test count(x->x < -b,rho)<=int(row,"negative_density")<=count(x->x < b,rho)
                @test num(row,"negative_density_fraction")==int(row,"negative_density")/length(rho)
                @test abs(num(row,"mean_negative_density")-sum(x->max(-x,0.),rho)/length(rho))<=b
            end
            rho=cubes["density"].grid; voxel=abs(det(g.geo.cell))/length(rho)
            bound=sum(v->N.halfquantum(v,p["cube_significant_digits"])+
                10N.halfquantum(v,p["native_significant_digits"])+p["numeric_roundoff_atol_ry"],rho)*voxel
            @test abs(num(gg,"charge_electrons")-sum(rho)*voxel)<=bound
            @test abs(num(gg,"negative_charge_electrons")-sum(x->max(-x,0.),rho)*voxel)<=bound
            empty!(cubes); GC.gc()
        end
    end
end

if isempty(ARGS)
    metadata_tests(ROOT) # XML-only local checks; no saved volume is loaded.
elseif length(ARGS)==3
    verify_saved(ARGS...)
else
    error("Usage: test_qe_density_cutoff.jl [ROOT SAVED_RUN MOLECULE]")
end
