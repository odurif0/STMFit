#!/usr/bin/env julia
using Test, Printf, LinearAlgebra, Statistics
include(joinpath(@__DIR__,"qe_vacuum_potential.jl"))
const R=QEVacuumPotential
const N=R.N
const P=R.P

function native_fixture(file,values;dims=(2,3,4),nzp=5,number=1,ibrav=0)
    open(file,"w") do io
        println(io,"synthetic native potential")
        println(io,join((size(values,1),size(values,2),nzp,dims...,1,1),' '))
        println(io,"$ibrav 10.00000000 0 0 0 0 0")
        println(io,"1 0 0\n0 1.5 0\n0 0 2")
        println(io,"100 7.2 50 $number\n1 H 1.00\n1 0.25 0.5 0.75 1")
        for (j,v) in enumerate(values)
            @printf(io,"%17.9E",v)
            j%5==0 && println(io)
        end
        println(io)
    end
end

function verify_saved(run)
    s=R.settings(joinpath(run,"settings.toml")); p=s["preprocessing"]
    geometry=R.V.tsv(joinpath(run,"geometry.tsv"))
    num(r,k)=parse(Float64,r[k]); integer(r,k)=parse(Int,r[k]); bool(r,k)=parse(Bool,r[k])
    @testset "Independent cube-plane summaries and selected-band barriers" begin
        for molecule in R.MOLECULES
            dir=joinpath(run,molecule); out=joinpath(dir,"analysis")
            meta=R.TOML.parsefile(joinpath(dir,"metadata.toml")); result=R.TOML.parsefile(joinpath(out,"summary.toml"))
            @test result["config_sha256"]==meta["config_sha256"]==R.S.sha(joinpath(run,"settings.toml"))
            @test result["xml_sha256"]==meta["xml_sha256"]==R.S.sha(joinpath(dir,"data-file-schema.xml"))
            @test result["exact_total_repeat"] && !result["matching_plane_adopted"] && !result["wavefunctions_propagated"]
            for (file,sha) in result["export_sha256"]; @test R.S.sha(joinpath(dir,file))==sha; end
            for ext in ("dat","cube"); @test R.S.sha(joinpath(dir,"total.$ext"))==R.S.sha(joinpath(dir,"repeat_total.$ext")); end
            cubes=Dict(name=>R.V.cube(joinpath(dir,name*".cube")) for name in ("total","electrostatic"))
            rows=R.V.tsv(joinpath(out,"planes.tsv")); bs=R.V.tsv(joinpath(out,"barriers.tsv"))
            sp=R.S.spectrum(joinpath(dir,"data-file-schema.xml"),s); bands=filter(b->b.ildos_weight>0,sp.rows)
            nz=cubes["total"].dims[3]; dims=cubes["total"].dims
            @test cubes["electrostatic"].dims==dims
            @test length(rows)==nz==result["planes"] && integer.(rows,"k")==collect(0:nz-1)
            @test result["selected_bands"]==length(bands) && length(bs)==nz*length(bands)
            @test Set((integer(r,"k"),integer(r,"band")) for r in bs)==Set((k,b.band) for k in 0:nz-1 for b in bands)
            gg=filter(r->r["molecule"]==molecule,geometry)
            @test length(gg)==nz && integer.(gg,"k")==collect(0:nz-1)
            @test count(r->bool(r,"paw_gap"),rows)==result["paw_free_planes"]
            for (j,row) in enumerate(rows)
                @test row["molecule"]==molecule && integer(row,"samples")==dims[1]*dims[2]
                @test all(row[k]==gg[j][k] for k in ("z_nm","paw_gap","lower_half"))
                arrays=Dict{String,Any}(name=>view(cubes[name].grid,j,:,:) for name in ("total","electrostatic"))
                arrays["xc"]=arrays["total"].-arrays["electrostatic"]
                bounds=Dict(name=>maximum(N.halfquantum(v,p["cube_significant_digits"])+
                    10N.halfquantum(v,p["native_significant_digits"])+p["numeric_roundoff_atol_ry"] for v in arrays[name])
                    for name in ("total","electrostatic"))
                bounds["xc"]=bounds["total"]+bounds["electrostatic"]
                for name in ("total","electrostatic","xc")
                    v=arrays[name]; mu=sum(v)/length(v); sigma=sqrt(sum((x-mu)^2 for x in v)/length(v))
                    # Population SD is 1-Lipschitz in the normalized L2 norm;
                    # the max printing envelope also bounds mean/min/max errors.
                    for (key,value) in (("mean",mu),("min",minimum(v)),("max",maximum(v)),("std",sigma))
                        @test abs(num(row,key*"_"*name*"_ry")-value)<=bounds[name]
                    end
                end
                @test num(row,"mean_total_minus_fermi_ev")==num(row,"mean_total_ry")*R.S.RY_EV-sp.fermi_ev
            end
            for r in bs
                row=rows[integer(r,"k")+1]; band=only(filter(b->b.band==integer(r,"band"),bands))
                @test r["molecule"]==molecule && all(r[k]==row[k] for k in ("z_nm","paw_gap","lower_half"))
                @test num(r,"energy_ev")==band.energy_ev
                lo=num(row,"min_total_ry")-band.energy_ev/R.S.RY_EV
                hi=num(row,"max_total_ry")-band.energy_ev/R.S.RY_EV
                avg=num(row,"mean_total_ry")-band.energy_ev/R.S.RY_EV
                @test num(r,"barrier_min_ev")==lo*R.S.RY_EV && num(r,"barrier_mean_ev")==avg*R.S.RY_EV
                @test num(r,"barrier_max_ev")==hi*R.S.RY_EV
                @test bool(r,"all_sampled_barriers_positive")== (lo>0)
                for (key,value) in (("kappa_min_nm_inv",lo>0 ? sqrt(lo)/R.D.G.C.BOHR_NM : NaN),
                    ("kappa_mean_potential_nm_inv",avg>0 ? sqrt(avg)/R.D.G.C.BOHR_NM : NaN),
                    ("kappa_max_nm_inv",hi>0 ? sqrt(hi)/R.D.G.C.BOHR_NM : NaN),
                    ("lateral_span_over_mean_barrier",lo>0 ? (hi-lo)/avg : NaN))
                    @test isequal(num(r,key),value)
                end
            end
        end
    end
end

@testset "Native QE order, padding, units and descriptive barriers" begin
    s=R.settings(joinpath(@__DIR__,"../config/qe_vacuum_potential.toml")); p=s["preprocessing"]
    @test R.S.HARTREE_EV==2R.S.RY_EV
    # Deliberately nonsquare: swapping x/z cannot pass. Padding is not data.
    storage=[Float64(i+10j+100k)/100 for i in 1:3,j in 1:4,k in 1:4]
    mktempdir() do dir
        file=joinpath(dir,"native.dat"); native_fixture(file,storage)
        a=N.readplot(file)
        @test a.dims==(2,3,4) && a.padded_dims==(3,4,5)
        @test a.storage==storage && a.grid==storage[1:2,1:3,:]
        @test a.plot_num==1 && a.cell_bohr==diagm([10.,15.,20.])
        @test a.atoms[1].position_bohr==[2.5,5.,7.5]
        geo=(;cell=a.cell_bohr,dims=collect(a.dims),geometry=(;atoms=[(;z=1,position_nm=[2.5,5.,7.5]*R.D.G.C.BOHR_NM)]))
        g=(;geo); c=(;dims=collect(a.dims),origin=zeros(3),axes=diagm([5.,5.,5.]),grid=permutedims(a.grid,(3,2,1)))
        check=R.validate_grid(a,c,g,p)
        @test check.voxels==24 && check.max_abs_difference_ry==0
        bad=copy(c.grid); bad[1,1,1]+=.1
        @test_throws ErrorException R.validate_grid(a,merge(c,(;grid=bad)),g,p)
        @test_throws ErrorException R.validate_grid(a,merge(c,(;axes=2c.axes)),g,p)
        @test_throws ErrorException R.validate_grid(a,merge(c,(;origin=ones(3))),g,p)
        native_fixture(joinpath(dir,"bad_lattice.dat"),storage;ibrav=1)
        @test_throws ErrorException N.readplot(joinpath(dir,"bad_lattice.dat"))
        # A nonfinite physical or padding payload is never silently accepted.
        badstorage=copy(storage); badstorage[3,4,4]=NaN
        native_fixture(joinpath(dir,"bad_value.dat"),badstorage)
        @test_throws ErrorException N.readplot(joinpath(dir,"bad_value.dat"))
        write(joinpath(dir,"truncated.dat"),read(file,String)[1:end-18])
        @test_throws Exception N.readplot(joinpath(dir,"truncated.dat"))
        write(joinpath(dir,"extra.dat"),read(file,String)*" 0.0\n")
        @test_throws ErrorException N.readplot(joinpath(dir,"extra.dat"))
    end
    @test N.halfquantum(0.,10)==0
    @test N.halfquantum(-1.2345,5)==.00005
    t=fill(2.,2,3); e=fill(1.5,2,3); r=P.plane_stats(t,e)
    @test r.samples==6 && r.mean_total_ry==r.min_total_ry==r.max_total_ry==2
    @test r.std_total_ry==0 && r.mean_xc_ry==.5 && r.std_xc_ry==0
    # With E=1 Ry and V=2 Ry, kappa=1/Bohr; no Hartree/Ry factor of two.
    b=P.barrier(r,R.S.RY_EV,R.S.RY_EV,R.D.G.C.BOHR_NM)
    @test b.barrier_min_ev==b.barrier_mean_ev==b.barrier_max_ev==R.S.RY_EV
    @test b.kappa_min_nm_inv==b.kappa_mean_potential_nm_inv==b.kappa_max_nm_inv==1/R.D.G.C.BOHR_NM
    @test b.all_sampled_barriers_positive && b.lateral_span_over_mean_barrier==0
    shifted=P.barrier(P.plane_stats(t.+3,e.+3),4R.S.RY_EV,R.S.RY_EV,R.D.G.C.BOHR_NM)
    @test isequal(b,shifted)
    for energy in (2R.S.RY_EV,3R.S.RY_EV)
        bb=P.barrier(r,energy,R.S.RY_EV,R.D.G.C.BOHR_NM)
        @test !bb.all_sampled_barriers_positive && isnan(bb.kappa_min_nm_inv)
        @test isnan(bb.kappa_mean_potential_nm_inv) && isnan(bb.lateral_span_over_mean_barrier)
    end
    t=[1. 2.;3. 4.]; r=P.plane_stats(t,t.-.25)
    @test r.mean_total_ry==2.5 && r.std_total_ry==sqrt(1.25)
    @test r.min_total_ry==1 && r.max_total_ry==4 && r.mean_xc_ry==.25
    b=P.barrier(r,0.,R.S.RY_EV,R.D.G.C.BOHR_NM)
    @test b.kappa_max_nm_inv==2b.kappa_min_nm_inv && b.lateral_span_over_mean_barrier==3/2.5
    @test_throws ErrorException P.plane_stats(fill(NaN,2,2),t)
    @test_throws ErrorException P.plane_stats(t,ones(2,3))
    @test occursin("plot_num=1",R.pp_input("glcn","total",1))
    @test occursin("plot_num=11",R.pp_input("glcnac","electrostatic",11))
    @test_throws ErrorException R.pp_input("other","total",1)
    @test_throws ErrorException R.pp_input("glcn","total",10)
end
isempty(ARGS) || (length(ARGS)==1 ? verify_saved(only(ARGS)) : error("Usage: test_qe_vacuum_potential.jl [SAVED_RUN]"))
