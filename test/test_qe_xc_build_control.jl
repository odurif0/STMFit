#!/usr/bin/env julia
using Test, TOML, Statistics
include(joinpath(@__DIR__,"test_qe_xc_components.jl"))
include(joinpath(@__DIR__,"qe_xc_build_control.jl"))
const B=QEXCBuildControl
const CONFIG=joinpath(@__DIR__,"../config/qe_xc_build_control.toml")

# Independent reference reduction. A naive serial generator sum can lose
# hundreds of ulps on a 43,200-point plane even though the production pairwise
# reduction is accurate. Keep the 16-eps comparison, not a data-fitted tolerance.
function compensated_squares(xs)
    total=0.; correction=0.
    for x in xs
        term=x*x; next=total+term
        correction+=abs(total)>=abs(term) ? (total-next)+term : (term-next)+total
        total=next
    end
    total+correction
end

@testset "Same-build controls and descriptive cross-build differences" begin
    s=B.settings(CONFIG); p=s["preprocessing"]
    @test_throws ErrorException C.settings(CONFIG)
    @test_throws ErrorException B.settings(joinpath(@__DIR__,"../config/qe_xc_components.toml"))
    @test s["preprocessing"]==P && !s["selection"]["adopt_reference"]
    @test occursin("plot_num=0",B.control_input("glcn","density"))
    @test occursin("plot_num=11",B.control_input("glcnac","electrostatic"))
    @test_throws ErrorException B.control_input("unknown","total")
    @test_throws ErrorException B.control_input("glcn","xc")
    r=B.differences([1.,2.,3.,0.],[1.,1.,4.,0.],p)
    @test r.samples==4 && r.different==2 && r.maximum_absolute==1
    @test r.mean_signed==0 && r.mean_absolute==.5 && r.rms==sqrt(.5)
    @test r.relative_l2==sqrt(2/18) && r.within_combined_print_bound==2
    @test B.differences(zeros(2),zeros(2),p).relative_l2==0
    @test isinf(B.differences(ones(2),zeros(2),p).relative_l2)
    @test_throws ErrorException B.differences([NaN],[1.],p)
    @test_throws ErrorException B.differences([1.],[1.,2.],p)
    @test_throws ErrorException B.differences(Float64[],Float64[],p)
    tiny=B.differences([1.0 + 2eps()],[1.],p)
    @test tiny.different==1 && tiny.within_combined_print_bound==1 # not byte identity
    @test compensated_squares([1.,2.,3.])==14
    @test compensated_squares(zeros(43_200))==0
    reference=fill(.3,240,180)
    perturbed=reference.+reshape([1e-9*cos(i/7) for i in 1:length(reference)],size(reference))
    @test isapprox(compensated_squares(reference),Float64(BigFloat(.3)^2*length(reference));rtol=eps())
    fullplane=B.differences(perturbed,reference,p)
    @test isapprox(fullplane.relative_l2,sqrt(compensated_squares(perturbed.-reference)/compensated_squares(reference));rtol=16eps())
    # Match readplot's range-indexed Cartesian view, not just a dense matrix
    # or a colon-indexed contiguous view. The unmaterialized reduction fails
    # this label-free fixture by about 2.78e-13 relatively, with zero padding too.
    for padding in (0,1,2)
        storage=fill(.3,240+padding,180,1)
        shifted=storage.+reshape([1e-9*cos(i/7) for i in eachindex(storage)],size(storage))
        native=view(storage,1:240,1:180,:); candidate=view(shifted,1:240,1:180,:)
        a=view(candidate,:,:,1); b=view(native,:,:,1)
        @test IndexStyle(typeof(b))==IndexCartesian()
        observed=B.differences(a,b,p).relative_l2
        expected=sqrt(compensated_squares(a.-b)/compensated_squares(b))
        @test isapprox(observed,expected;rtol=16eps())
        @test observed==B.differences(copy(a),copy(b),p).relative_l2
    end
    mktempdir() do dir
        cp(CONFIG,joinpath(dir,"settings.toml"))
        for m in B.R.MOLECULES
            for (a,b) in B.identity_pairs(), file in (a,b)
                path=joinpath(dir,m,file); mkpath(dirname(path)); write(path,"identical fixture\n")
            end
            for file in B.CONTROL_FILES
                path=joinpath(dir,m,file); mkpath(dirname(path))
                isfile(path) || write(path,"control fixture\n")
            end
            mkpath(joinpath(dir,m,"controls"))
            open(io->TOML.print(io,Dict("same_build_controls_exact"=>true,
                "config_sha256"=>B.S.sha(CONFIG),"control_exports_sha256"=>Dict(file=>B.S.sha(joinpath(dir,m,file)) for file in B.CONTROL_FILES))),
                joinpath(dir,m,"controls/summary.toml"),"w")
            records=[(;a,b,sha_a=B.S.sha(joinpath(dir,m,a)),sha_b=B.S.sha(joinpath(dir,m,b)),exact=true) for (a,b) in B.identity_pairs()]
            B.S.table(joinpath(dir,m,"controls/identities.tsv"),records)
        end
        @test isnothing(B.require_controls(dir))
        write(joinpath(dir,"glcn/rho_valence.dat"),"different fixture\n")
        @test_throws ErrorException B.require_controls(dir)
        write(joinpath(dir,"glcn/stock/density.dat"),"different fixture\n")
        @test_throws ErrorException B.require_controls(dir) # jointly changed pair is still rejected
    end
end

function verify_controls(root,run)
    s=B.settings(joinpath(run,"settings.toml")); p=s["preprocessing"]
    num(r,k)=parse(Float64,r[k]); integer(r,k)=parse(Int,r[k]); boolean(r,k)=parse(Bool,r[k])
    @testset "Independent saved paired controls and full native differences" begin
        for m in B.R.MOLECULES
            dir=joinpath(run,m); old=joinpath(root,s["model"]["source_run"],m)
            saved=TOML.parsefile(joinpath(dir,"controls/summary.toml"))
            @test saved["same_build_controls_exact"] && !saved["cross_build_acceptance_decision"] && !saved["potential_adopted"]
            @test Set(keys(saved["control_exports_sha256"]))==Set(B.CONTROL_FILES)
            for (file,sha) in saved["control_exports_sha256"]; @test B.S.sha(joinpath(dir,file))==sha; end
            identities=B.V.tsv(joinpath(dir,"controls/identities.tsv"))
            @test Set((r["a"],r["b"]) for r in identities)==Set(B.identity_pairs())
            for r in identities
                @test boolean(r,"exact") && r["sha_a"]==r["sha_b"]==B.S.sha(joinpath(dir,r["a"]))==B.S.sha(joinpath(dir,r["b"]))
            end
            rows=B.V.tsv(joinpath(dir,"controls/differences.tsv")); totals=B.V.tsv(joinpath(dir,"controls/totals.tsv"))
            geometry=filter(r->r["molecule"]==m,B.V.tsv(joinpath(run,"geometry.tsv")))
            pairs=(("same_build_density",joinpath(dir,"stock/density.dat"),joinpath(dir,"rho_valence.dat")),
                ("same_build_total",joinpath(dir,"stock/total.dat"),joinpath(dir,"xc_total.dat")),
                ("cross_build_density",joinpath(dir,"stock/density.dat"),joinpath(dir,"site/density.dat")),
                ("cross_build_total",joinpath(dir,"stock/total.dat"),joinpath(old,"total.dat")),
                ("cross_build_electrostatic",joinpath(dir,"stock/electrostatic.dat"),joinpath(old,"electrostatic.dat")))
            @test length(totals)==5 && length(rows)==5length(geometry)
            for (name,afile,bfile) in pairs
                a=B.N.readplot(afile); b=B.N.readplot(bfile)
                rr=filter(r->r["comparison"]==name,rows)
                t=only(filter(r->r["comparison"]==name,totals))
                @test t["sha_a"]==B.S.sha(afile) && t["sha_b"]==B.S.sha(bfile)
                @test integer.(rr,"k")==collect(0:a.dims[3]-1)
                for (k,r) in enumerate(rr)
                    @test all(r[q]==geometry[k][q] for q in ("z_nm","paw_gap","lower_half"))
                    delta=vec(a.grid[:,:,k]).-vec(b.grid[:,:,k]); ref=vec(b.grid[:,:,k]); n=length(delta)
                    @test integer(r,"samples")==n && integer(r,"different")==sum(x!=0 for x in delta)
                    @test num(r,"maximum_absolute")==maximum(abs,delta)
                    @test isapprox(num(r,"mean_signed"),sum(delta)/n;rtol=16eps(),atol=eps()*sum(abs,delta)/n)
                    @test isapprox(num(r,"mean_absolute"),sum(abs,delta)/n;rtol=16eps())
                    @test isapprox(num(r,"rms"),sqrt(compensated_squares(delta)/n);rtol=16eps())
                    normref=compensated_squares(ref); normdiff=compensated_squares(delta)
                    expected=normref>0 ? sqrt(normdiff/normref) : iszero(normdiff) ? 0. : Inf
                    @test isapprox(num(r,"relative_l2"),expected;rtol=16eps())
                    count=0; ratio=0.
                    for (x,y) in zip(a.grid[:,:,k],b.grid[:,:,k])
                        bound=B.N.halfquantum(x,10)+B.N.halfquantum(y,10); d=abs(x-y)
                        count+=d<=bound; ratio=max(ratio,bound>0 ? d/bound : d==0 ? 0. : Inf)
                    end
                    @test integer(r,"within_combined_print_bound")==count
                    @test num(r,"max_error_to_combined_print_bound")==ratio
                end
                @test integer(t,"different")==sum(integer(r,"different") for r in rr)
                @test num(t,"maximum_absolute")==maximum(num(r,"maximum_absolute") for r in rr)
                @test integer(t,"within_combined_print_bound")==sum(integer(r,"within_combined_print_bound") for r in rr)
            end
            @test saved["comparison_voxels"]==sum(integer(t,"samples") for t in totals)
        end
    end
end

if !isempty(ARGS)
    length(ARGS)==2 || error("test_qe_xc_build_control.jl [ROOT SAVED_RUN]")
    verify_controls(ARGS...)
    verify_saved(ARGS...;load_settings=B.settings,reference=(root,run,s,m)->joinpath(run,m,"stock"),
        control_key="exact_same_build_controls_and_repeats")
end
