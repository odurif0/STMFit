#!/usr/bin/env julia
# Read-only partition/root/output consistency; not a third WFC reconstruction.
module VerifyQEVacuumContinuous
using Test, TOML, Statistics
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const V=VerifyQESpectralWindow
n(r,k)=parse(Float64,r[k])
i(r,k)=parse(Int,r[k])
b(r,k)=parse(Bool,r[k])

function verify(run)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    settings=TOML.parsefile(joinpath(run,"settings.toml")); q=settings["selection"]
    common=filter(r->r["source"]=="joint" && i(r,"interval")>0,V.tsv(joinpath(run,"common_intervals.tsv")))
    ids=sort(i.(common,"interval")); pixels=(1+round(Int,2settings["preprocessing"]["half_nm"]/settings["preprocessing"]["step_nm"]))^2
    @testset "Saved continuous root partitions and refined surfaces" begin
        for m in ("glcn","glcnac")
            out=joinpath(run,m,"continuous"); summary=TOML.parsefile(joinpath(out,"summary.toml"))
            @test all(values(summary["checks"]))
            @test summary["config_sha256"]==V.digest(joinpath(run,"settings.toml"))
            @test summary["interval_arithmetic_version"]==q["interval_arithmetic_version"]
            meta=TOML.parsefile(joinpath(run,m,"metadata.toml"))
            @test summary["wfc_sha256"]==meta["wfc_sha256"]
            for (name,sha) in meta["references"]; @test V.digest(joinpath(run,name))==sha; end
            roots=V.tsv(joinpath(out,"roots.tsv")); leaves=V.tsv(joinpath(out,"leaves.tsv"))
            profiles=V.tsv(joinpath(out,"profiles.tsv")); surfaces=V.tsv(joinpath(out,"surfaces.tsv"))
            @test length(roots)==pixels*length(ids)==summary["targets"]
            @test Set((i(r,"pixel"),i(r,"interval")) for r in roots)==Set((p,j) for p in 1:pixels for j in ids)
            partitions=Dict{Tuple{Int,Int},Vector{typeof(first(leaves))}}()
            for leaf in leaves
                push!(get!(partitions,(i(leaf,"pixel"),i(leaf,"interval")),typeof(first(leaves))[]),leaf)
            end
            @test Set(keys(partitions))==Set((p,j) for p in 1:pixels for j in ids)
            for root in roots
                key=(i(root,"pixel"),i(root,"interval")); rows=partitions[key]
                @test n(first(rows),"lo")==summary["domain_lower_nm"]
                @test n(last(rows),"hi")==summary["domain_upper_nm"]
                @test all(n(rows[k],"hi")==n(rows[k+1],"lo") for k in 1:length(rows)-1)
                iso=n(root,"isovalue")
                @test iso==n(only(filter(r->i(r,"interval")==key[2],common)),"isovalue")
                for leaf in rows
                    @test n(leaf,"lo")<n(leaf,"hi") && n(leaf,"isovalue")==iso
                    status=leaf["status"]
                    if status=="excluded"
                        @test iso<n(leaf,"rho_lo") || iso>n(leaf,"rho_hi")
                    elseif status in ("unique","excluded_monotone","precision_limit")
                        @test (i(leaf,"direction")==-1 && n(leaf,"derivative_hi")<0) ||
                            (i(leaf,"direction")==1 && n(leaf,"derivative_lo")>0)
                        if status!="excluded_monotone"
                            @test n(leaf,"lo")<=n(leaf,"root_lo")<=n(leaf,"root_hi")<=n(leaf,"hi")
                            status=="unique" && @test n(leaf,"root_hi")-n(leaf,"root_lo")<=q["root_width_nm"]
                        end
                    else
                        @test status in ("unresolved","node_limit")
                    end
                end
                unique_rows=filter(r->r["status"]=="unique",rows)
                unresolved=count(r->r["status"] in ("unresolved","node_limit","precision_limit"),rows)
                @test length(unique_rows)==i(root,"root_count") && unresolved==i(root,"unresolved")
                valid=unresolved==0 && length(unique_rows)==1 && i(only(unique_rows),"direction")==-1
                @test valid==b(root,"valid")
                if valid
                    @test n(root,"z_lower_nm")==n(only(unique_rows),"root_lo")
                    @test n(root,"z_upper_nm")==n(only(unique_rows),"root_hi")
                    @test n(root,"z_lower_nm")<=n(root,"z_nm")<=n(root,"z_upper_nm")
                else
                    @test isnan(n(root,"z_nm"))
                end
            end
            @test count(r->b(r,"valid"),roots)==summary["valid_roots"]==length(surfaces)
            lookup=Dict((i(r,"pixel"),i(r,"interval"))=>r for r in roots if b(r,"valid"))
            @test Set((i(r,"pixel"),i(r,"interval")) for r in surfaces)==Set(keys(lookup))
            @test sum(i(r,"unresolved") for r in roots)==summary["unresolved_leaves"]
            @test summary["certified_complete_maps"]==all(r->b(r,"valid"),roots)
            previous=V.tsv(joinpath(run,m*"_profiles.tsv"))
            @test length(profiles)==length(previous)
            for (row,ref) in zip(profiles,previous)
                @test i(row,"pixel")==i(ref,"pixel") && n(row,"z_nm")==n(ref,"z_nm") && n(row,"direct")==n(ref,"direct")
                distance=max(0.,n(row,"lower")-n(row,"direct"),n(row,"direct")-n(row,"upper"))
                @test n(row,"lower")<=n(row,"upper")
                @test distance==n(row,"distance")<=q["numeric_roundoff_atol"]
            end
            for row in surfaces
                root=lookup[(i(row,"pixel"),i(row,"interval"))]
                for key in keys(root); @test row[key]==root[key]; end
                @test n(row,"direct")==n(row,"repeat")
                @test n(row,"distance")==max(0.,n(row,"density_lower")-n(row,"direct"),n(row,"direct")-n(row,"density_upper"))
                @test n(row,"distance")<=q["numeric_roundoff_atol"]
                @test n(row,"relative_density_error")== (n(row,"direct")-n(row,"isovalue"))/n(row,"isovalue")
            end
            println(m,": all targets=",length(roots),", valid=",length(surfaces),
                ", unresolved leaves=",summary["unresolved_leaves"],
                ", profiles outside strict enclosure=",count(r->n(r,"distance")>0,profiles))
            for id in ids
                rows=filter(r->i(r,"interval")==id,surfaces)
                isempty(rows) && continue
                errors=abs.(n.(rows,"relative_density_error")); shifts=abs.(n.(rows,"z_nm")-n.(rows,"old_z_nm"))
                println((;molecule=m,interval=id,pixels=length(rows),max_relative_error=maximum(errors),
                    median_relative_error=median(errors),max_height_shift_nm=maximum(shifts),
                    max_root_width_nm=maximum(n(r,"z_upper_nm")-n(r,"z_lower_nm") for r in rows),
                    outside_strict_enclosure=count(r->n(r,"distance")>0,rows)))
            end
        end
    end
end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    ARGS==["--help"] ? println("verify_qe_vacuum_continuous.jl RUN_DIR (saved tables only)") :
        (length(ARGS)==1 ? VerifyQEVacuumContinuous.verify(only(ARGS)) : error("Expected saved run directory"))
end
