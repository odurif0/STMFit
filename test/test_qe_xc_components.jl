#!/usr/bin/env julia
using Test, Statistics, TOML
include(joinpath(@__DIR__, "qe_xc_components.jl"))
const C=QEXCComponents
const P=C.settings(joinpath(@__DIR__,"../config/qe_xc_components.toml"))["preprocessing"]

@testset "XC decomposition, source branches and empty subsets" begin
    a=reshape([1.,2.,3.,4.],2,2); b=fill(.5,2,2); c=fill(.2,2,2); d=fill(.1,2,2)
    check=C.additivity(a,b,c,d,a.-b.-c.-d,P)
    @test check.voxels==4 && check.max_error_to_bound<1
    @test_throws ErrorException C.additivity(a,b,c,d,a,P)
    @test_throws ErrorException C.additivity(a,b,c,d,fill(NaN,2,2),P)
    @test_throws ErrorException C.additivity(a,b,c,d,zeros(3,2),P)
    @test C.description([1.,2.,3.]) == (;samples=3,minimum=1.,mean=2.,maximum=3.,spatial_sd=sqrt(2/3))
    empty=C.description(Float64[])
    @test empty.samples==0 && all(isnan,(empty.minimum,empty.mean,empty.maximum,empty.spatial_sd))
    rho=[1e-3,-1e-3,1e-8,2e-4,0.]; grad2=[1e-6,1e-6,1e-6,1e-12,0.]
    active=[1.,1.,0.,0.,0.]; local_gga=[1.,2.,0.,0.,0.]
    @test C.source_masks(rho,grad2,active,local_gga,P)==(;definite=5,ambiguous=0)
    @test_throws ErrorException C.source_masks(rho,grad2,zeros(5),local_gga,P)
    @test_throws ErrorException C.source_masks(rho,grad2,[1.,0.,0.,0.,0.],zeros(5),P)
    @test_throws ErrorException C.source_masks(rho,fill(-1.,5),active,local_gga,P)
    @test_throws ErrorException C.source_masks(rho,grad2,fill(.5,5),local_gga,P)
    @test C.source_masks([1e-6],[1.],[0.],[0.],P)==(;definite=0,ambiguous=1)
    @test C.source_masks([1.],[1e-10],[0.],[0.],P)==(;definite=0,ambiguous=1)
    @test occursin("../work/glcn",C.extraction_input("glcn"))
    @test occursin("../../work/glcn",C.extraction_input("glcn",true))
    @test occursin("nfile=1",C.cube_input("lda")) && !occursin("plot_num",C.cube_input("lda"))
    @test_throws ErrorException C.extraction_input("unknown")
    @test_throws ErrorException C.cube_input("other")
end

function verify_saved(root,run; load_settings=C.settings,
        reference=(root,run,s,m)->joinpath(root,s["model"]["source_run"],m),
        control_key="exact_total_reference_and_repeats")
    s=load_settings(joinpath(run,"settings.toml")); p=s["preprocessing"]
    geometry=C.V.tsv(joinpath(run,"geometry.tsv"))
    number(r,k)=parse(Float64,r[k]); boolean(r,k)=parse(Bool,r[k])
    @testset "Independent component cube summaries and branches" begin
        for molecule in C.R.MOLECULES
            dir=joinpath(run,molecule); old=reference(root,run,s,molecule)
            saved=TOML.parsefile(joinpath(dir,"analysis/summary.toml"))
            @test saved[control_key] && !saved["density_modified"] && !saved["potential_adopted"]
            @test saved["config_sha256"]==C.S.sha(joinpath(run,"settings.toml"))
            for (file,sha) in saved["component_sha256"]; @test C.S.sha(joinpath(dir,file))==sha; end
            for field in C.FIELDS
                @test C.S.sha(joinpath(dir,field*".dat"))==C.S.sha(joinpath(dir,"repeat",field*".dat"))
            end
            cubes=Dict(name=>C.V.cube(joinpath(dir,name*".cube")) for name in C.FIELDS)
            total=C.V.cube(joinpath(dir,"xc_total.cube")); electro=C.V.cube(joinpath(old,"electrostatic.cube"))
            @test all(c.dims==total.dims && c.axes==total.axes && c.origin==total.origin for c in values(cubes))
            counts=C.V.tsv(joinpath(dir,"analysis/counts.tsv")); rows=C.V.tsv(joinpath(dir,"analysis/components.tsv"))
            nz=total.dims[3]
            @test length(counts)==nz==saved["planes"] && length(rows)==nz*3*8
            @test saved["native_cube_voxels"]==prod(total.dims)*length(C.FIELDS)
            lookup=Dict((parse(Int,r["k"]),r["subset"],r["component"])=>r for r in rows)
            @test length(lookup)==length(rows)
            gg=filter(r->r["molecule"]==molecule,geometry)
            for k in 1:nz
                fields=Dict{String,Any}(name=>view(cubes[name].grid,k,:,:) for name in C.FIELDS)
                fields["xc_total"]=view(total.grid,k,:,:).-view(electro.grid,k,:,:)
                mask=fields["gga_active"]
                @test all(x->x in (0.,1.),mask)
                @test parse(Int,counts[k]["samples"])==length(mask)
                @test parse(Int,counts[k]["gga_active"])==sum(mask)
                @test all(counts[k][q]==gg[k][q] for q in ("z_nm","paw_gap","lower_half"))
                bounds=Dict(name=>maximum(C.N.halfquantum(v,5)+10C.N.halfquantum(v,10)+1e-13 for v in fields[name]) for name in C.FIELDS)
                bounds["xc_total"]=sum(maximum(C.N.halfquantum(v,5)+10C.N.halfquantum(v,10)+1e-13 for v in view(x.grid,k,:,:)) for x in (total,electro))
                for subset in ("all","gga_active","gga_inactive")
                    selected=subset=="all" ? trues(size(mask)) : mask.==(subset=="gga_active" ? 1. : 0.)
                    for name in (C.FIELDS...,"xc_total")
                        r=lookup[(k-1,subset,name)]; samples=fields[name][selected]; n=length(samples)
                        @test r["molecule"]==molecule && all(r[q]==gg[k][q] for q in ("z_nm","paw_gap","lower_half"))
                        @test parse(Int,r["samples"])==n
                        if n==0
                            @test all(isnan(number(r,q)) for q in ("minimum","maximum","mean","spatial_sd"))
                        else
                            mu=sum(samples)/n; sd=sqrt(sum((v-mu)^2 for v in samples)/n)
                            for (q,value) in (("minimum",minimum(samples)),("maximum",maximum(samples)),("mean",mu),("spatial_sd",sd))
                                @test abs(number(r,q)-value)<=bounds[name]
                            end
                        end
                    end
                end
                # Independent rounded-cube additivity, not a re-call of the native checker.
                reconstructed=fields["lda"].+fields["gga_local"].+fields["gga_divergence"]
                bound=sum(bounds[q] for q in ("lda","gga_local","gga_divergence","xc_total"))
                @test maximum(abs.(reconstructed.-fields["xc_total"]))<=bound
                @test all(iszero,fields["gga_local"][mask.==0])
            end
        end
    end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    isempty(ARGS) || (length(ARGS)==2 ? verify_saved(ARGS...) : error("test_qe_xc_components.jl [ROOT SAVED_RUN]"))
end
