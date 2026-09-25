#!/usr/bin/env julia
# Saved-table check with an independent global crossing-count sweep.
using Test, TOML, LinearAlgebra, Statistics
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow
num(row,key)=parse(Float64,row[key])

function verify(run)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    out=joinpath(run,"analysis"); summary=TOML.parsefile(joinpath(out,"summary.toml"))
    config=TOML.parsefile(joinpath(run,"settings.toml")); p=config["preprocessing"]
    pixels=(1+round(Int,2p["half_nm"]/p["step_nm"]))^2
    geometry=Dict(r["molecule"]=>r for r in V.tsv(joinpath(run,"geometry_before_density.tsv")))
    common=V.tsv(joinpath(out,"common_intervals.tsv")); profiles=Dict{String,Any}()
    @testset "Independent saved vacuum-crossing evidence" begin
        @test all(values(summary["checks"]))
        @test !summary["direct_between_knots_uniqueness_proven"]
        @test V.digest(joinpath(run,"settings.toml"))==summary["config_sha256"]
        for m in ("glcn","glcnac")
            rows=V.tsv(joinpath(out,m*"_profiles.tsv")); g=geometry[m]
            nz=parse(Int,g["full_knots"]); z=num.(rows[1:nz],"z_nm")
            @test length(rows)==pixels*nz==summary["molecules"][m]["profile_points"]
            @test all(>(0),diff(z))
            @test first(z)==num(g,"first_knot_nm") && last(z)==num(g,"last_knot_nm")
            @test num(g,"lower_paw_nm")<first(z)<last(z)<num(g,"upper_paw_nm")
            for (i,r) in enumerate(rows)
                @test parse(Int,r["pixel"])==1+(i-1)÷nz
                @test num(r,"z_nm")==z[1+(i-1)%nz]
            end
            cube=reshape(num.(rows,"cube"),nz,pixels)
            direct=reshape(num.(rows,"direct"),nz,pixels)
            @test all(isfinite,direct) && all(>(0),direct) && all(>(0),cube)
            @test norm(cube-direct)/norm(direct)≈summary["molecules"][m]["interpolation_relative_l2"] rtol=1e-13
            profiles[m]=(;z,cube,direct,rows)
        end
        for domain in config["selection"]["domains"], source in ("cube","direct")
            # Count all segment crossings simultaneously over the UNION of
            # density knots. This does not call the per-column interval code.
            events=Tuple{Float64,Int,Int,Int}[]; singular=Float64[]; column=0
            for m in ("glcn","glcnac")
                pr=profiles[m]; n=domain=="full_paw_gap" ? length(pr.z) : parse(Int,geometry[m]["half_knots"])
                vals=getproperty(pr,Symbol(source))
                for pixel in 1:pixels
                    column+=1; v=vals[1:n,pixel]
                    append!(singular,[first(v),last(v)])
                    for i in 2:n-1
                        ((v[i]-v[i-1])*(v[i+1]-v[i])<=0) && push!(singular,v[i])
                    end
                    for (a,b) in zip(v[1:end-1],v[2:end])
                        a==b && (push!(singular,a); continue)
                        lo,hi=minmax(a,b); d=a>b ? 1 : 0; u=1-d
                        push!(events,(lo,column,d,u)); push!(events,(hi,column,-d,-u))
                    end
                end
            end
            sort!(events;by=first); counts=zeros(Int,column,2); bad=column; i=1
            claimed=[(num(r,"lo"),num(r,"hi")) for r in common if r["domain"]==domain &&
                r["source"]==source && parse(Int,r["interval"])>0]
            while i<=length(events)
                lo=events[i][1]
                while i<=length(events) && events[i][1]==lo
                    _,c,d,u=events[i]; was=counts[c,1]==1 && counts[c,2]==0
                    counts[c,1]+=d; counts[c,2]+=u
                    now=counts[c,1]==1 && counts[c,2]==0
                    bad+=Int(was)-Int(now); i+=1
                end
                i>length(events) && break
                hi=events[i][1]; iso=lo+(hi-lo)/2
                # Adjacent floating-point knots have no representable interior.
                lo<iso<hi || continue
                @test (bad==0)==any(a<iso<b for (a,b) in claimed)
            end
            @test all(==(0),counts)
            @test !any(any(a<x<b for (a,b) in claimed) for x in singular)
            println(domain," ",source,": complete global sweep agrees; intervals=",claimed)
        end
        for domain in config["selection"]["domains"]
            ranges(source)=[(num(r,"lo"),num(r,"hi")) for r in common if r["domain"]==domain &&
                r["source"]==source && parse(Int,r["interval"])>0]
            expected=sort([(max(a,c),min(b,d)) for (a,b) in ranges("cube") for (c,d) in ranges("direct")
                if max(a,c)<min(b,d)])
            @test ranges("joint")==expected
        end
        surfaces=V.tsv(joinpath(out,"surfaces.tsv"))
        @test length(surfaces)==summary["surface_rows"]
        @test length(surfaces)==2pixels*count(r->r["source"]=="joint" && parse(Int,r["interval"])>0,common)
        for r in surfaces
            m=r["molecule"]; pr=profiles[m]; pixel=parse(Int,r["pixel"]); iso=num(r,"isovalue")
            n=parse(Int,geometry[m][r["domain"]=="full_paw_gap" ? "full_knots" : "half_knots"])
            claim=only(filter(c->c["source"]=="joint" && c["domain"]==r["domain"] && c["interval"]==r["interval"],common))
            @test num(claim,"lo")<iso<num(claim,"hi") && iso==num(claim,"isovalue")
            @test iso≈exp((log(num(claim,"lo"))+log(num(claim,"hi")))/2) rtol=1e-14
            position=pr.rows[1+(pixel-1)*length(pr.z)]
            @test num(r,"x_nm")==num(position,"x_nm") && num(r,"y_nm")==num(position,"y_nm")
            for (field,height) in ((pr.direct,"z_nm"),(pr.cube,"cube_z_nm"))
                v=field[1:n,pixel]
                roots=findall(i->min(v[i],v[i+1])<iso<max(v[i],v[i+1]),1:n-1)
                @test length(roots)==1
                k=only(roots); @test v[k]>iso>v[k+1]
                z=pr.z[k]+(iso-v[k])/(v[k+1]-v[k])*(pr.z[k+1]-pr.z[k])
                @test num(r,height)≈z rtol=1e-14
            end
            @test num(r,"relative_density_error")≈(num(r,"direct_at_surface")-iso)/iso rtol=1e-14
            @test num(r,"ring_relative_vertical_nm")≈num(r,"z_nm")-num(geometry[m],"ring_z_nm") rtol=1e-14
        end
        for interval in sort(unique(parse(Int,r["interval"]) for r in surfaces)), m in ("glcn","glcnac")
            rows=filter(r->r["molecule"]==m && parse(Int,r["interval"])==interval,surfaces)
            @test sort(parse.(Int,getindex.(rows,"pixel")))==collect(1:pixels)
            errors=abs.(num.(rows,"relative_density_error")); z=num.(rows,"z_nm")
            delta=abs.(z-num.(rows,"cube_z_nm"))
            println((;interval,molecule=m,isovalue=num(first(rows),"isovalue"),
                min_z=minimum(z),mean_z=mean(z),max_z=maximum(z),
                median_abs_density_error=median(errors),max_abs_density_error=maximum(errors),
                max_cube_height_difference_nm=maximum(delta)))
        end
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    ARGS==["--help"] ? println("verify_qe_vacuum_crossings.jl RUN_DIR (read-only saved tables)") :
        (length(ARGS)==1 ? verify(only(ARGS)) : error("Expected saved run directory"))
end
