using Test, LinearAlgebra, Random, Statistics, TOML, Printf
include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .EmpiricalFisherNative, .ReconstructedUnitAssignment
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONFIGS = Dict(name => joinpath(ROOT,"config","unit_assignment_"*suffix*".toml") for
    (name,suffix) in (("control","centered_fisher"),("support","patch_support"),("volume","gaussian_score")))
const BASE = load_fisher_config(CONFIGS["control"])
const SUPPORT = load_fisher_config(CONFIGS["support"])
const GRID = fisher_grid(BASE)
const COORDS = collect(-4:4) .* .08
const U, T = repeat(COORDS;inner=9), repeat(COORDS;outer=9)
const KIND = "affine_residual_half_plane_asymmetry"
desc(values; mode="complete_disk_symmetric", coords=COORDS) =
    transverse_descriptor(values,coords,KIND;zero_l1=1e-12,patch_support=mode)

@testset "Independent explicit hypotheses, no fitted thresholds" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    base = load_config(CONFIGS["control"])
    for (name,section,field,value) in (("support","preprocessing","assignment_patch_support","complete_disk_symmetric"),
                                     ("volume","model","gmm_final_score","gaussian_density"))
        cfg = load_config(CONFIGS[name])
        @test cfg[section][field] == value
        cfg[section][field] = base[section][field]
        cfg["model"]["name"] = base["model"]["name"]
        @test cfg == base
    end
    mktempdir() do dir
        for (section,field) in (("preprocessing","assignment_patch_support"),("model","gmm_final_score")),
            value in (nothing,"guess",true,0.5)
            cfg = deepcopy(base)
            value === nothing ? delete!(cfg[section],field) : (cfg[section][field] = value)
            path = joinpath(dir,"invalid.toml")
            open(io->TOML.print(io,cfg),path,"w")
            @test_throws ErrorException load_config(path)
            if field == "assignment_patch_support"
                @test_throws ArgumentError load_fisher_config(cfg)
            else
                @test_throws ArgumentError GMM._load_final_score(cfg)
            end
        end
        for (name,section,field,value) in (("support","model","descriptor","transverse_first_moment"),
                                         ("volume","selection","gmm_selftrain",0))
            cfg=load_config(CONFIGS[name]); cfg[section][field]=value
            path=joinpath(dir,"invalid.toml"); open(io->TOML.print(io,cfg),path,"w")
            @test_throws ErrorException load_config(path)
        end
    end
end

@testset "Complete patches retain identical arithmetic" begin
    rng=MersenneTwister(819)
    for _ in 1:20
        p=randn(rng,81)
        @test desc(p) == desc(p;mode="full_square")
    end
    @test isequal(desc(zeros(81)),desc(zeros(81);mode="full_square"))
    @test isequal(desc(ones(81)),desc(ones(81);mode="full_square"))
end

@testset "Observed symmetric support: no imputation, independent QR and invariances" begin
    shape=exp.(-.5 .* (((U .- .12) ./ .055).^2 .+ ((T .- .04) ./ .07).^2))
    p=copy(shape); p[1]=NaN; original=copy(p)
    keep,reason=ReconstructedUnitAssignment._symmetric_patch_support(p,COORDS)
    @test reason == "ok_masked_symmetric"
    @test findall(.!keep) == [1,9,73,81]
    @test count(keep)==77
    @test keep == vec(reverse(reshape(keep,9,9);dims=1))
    @test keep == vec(reverse(reshape(keep,9,9);dims=2))
    disk=[i for (i,(u,t)) in enumerate((u,t) for u in -4:4 for t in -4:4) if u*u+t*t<=16]
    @test length(disk)==49 && all(keep[disk])
    design=hcat(ones(count(keep)),U[keep],T[keep])
    @test rank(design)==3
    residual=p[keep]-design*(design\p[keep])
    expected=sum(sign.(U[keep]).*residual)/sum(abs,residual)
    @test first(desc(p)) ≈ expected atol=2e-15
    @test last(desc(p)) == "ok_masked_symmetric"
    @test last(desc(p;mode="full_square")) == "nonfinite_patch"
    @test isequal(original,p)
    # Values on excluded reflection partners must have no effect at all.
    changed=copy(p); changed[[9,73,81]].=[1e9,-1e6,47.]
    @test desc(changed)==desc(p)
    for bad in (Inf,-Inf,NaN)
        q=copy(p); q[1]=bad
        @test desc(q)==desc(p)
    end
    for coeff in ((.3,1.2,-.7),(-5.,-2.,9.))
        plane=coeff[1] .+ coeff[2]*U .+ coeff[3]*T
        @test first(desc(p+plane)) ≈ expected atol=3e-14
        plane[1]=NaN
        @test last(desc(plane)) == "zero_affine_residual_mass"
    end
    @test first(desc(7p)) ≈ expected atol=2e-15
    @test first(desc(-p)) ≈ -expected atol=2e-15
    @test first(desc(vec(reverse(reshape(p,9,9);dims=2)))) ≈ -expected atol=2e-15
    @test first(desc(vec(reverse(reshape(p,9,9);dims=1)))) ≈ expected atol=2e-15
    # Every one of the 49 disk pixels is mandatory, even an edge pixel.
    for i in disk
        q=copy(shape); q[i]=NaN
        @test last(desc(q)) == "incomplete_descriptor_disk"
        @test isnan(first(desc(q)))
    end
    diskonly=copy(shape); diskonly[setdiff(1:81,disk)].=NaN
    @test isfinite(first(desc(diskonly)))
    @test count(first(ReconstructedUnitAssignment._symmetric_patch_support(diskonly,COORDS)))==49
    @test last(desc(fill(NaN,81))) == "incomplete_descriptor_disk"
    @test_throws ErrorException desc(p;coords=COORDS .+ .01)
    @test_throws ErrorException desc(p;mode="guess")
    @test_throws ErrorException transverse_descriptor(p,COORDS,"transverse_first_moment";
        zero_l1=1e-12,patch_support="complete_disk_symmetric")
end

@testset "Fisher uses its complete scoring disk, keeps invalid keys and schema checks" begin
    rng=MersenneTwister(820)
    full=.06randn(rng,64,289)
    shape=[exp(-20(u*u+t*t))+.15sin(9u)*cos(11t) for u in GRID.coords for t in GRID.coords]
    for i in 1:64
        full[i,:] .+= (i<=44 ? -2. : 2.) .* shape .+ .4
    end
    outside=setdiff(1:289,GRID.disk_indices)
    @test length(GRID.disk_indices)==197
    mktempdir() do dir
        header=vcat(["file","lobe"],[@sprintf("res_p%03d",i) for i in 1:289])
        rows=[merge(Dict("file"=>"synthetic.sxm","lobe"=>string(i)),
            Dict(header[j+2]=>string(full[i,j]) for j in 1:289)) for i in 1:64]
        path=joinpath(dir,"complete.tsv"); write_table(path,header,rows)
        a,b=load_patches(path,"res",BASE),load_patches(path,"res",SUPPORT)
        @test a.keys==b.keys && a.X==b.X && a.amplitudes==b.amplitudes
        expected=cv_scores(a,BASE)
        out1,out2=joinpath.(dir,("legacy.tsv","complete_support.tsv"))
        write_scores(out1,expected); write_scores(out2,cv_scores(b,SUPPORT))
        @test read(out1)==read(out2)
        # Remove ALL outside-disk values; no disk vector, fit or score may change.
        for row in rows, j in outside
            row[header[j+2]] = isodd(j) ? "NA" : "NaN"
        end
        partial=joinpath(dir,"partial.tsv"); write_table(partial,header,rows)
        c=load_patches(partial,"res",SUPPORT)
        @test c.X==a.X && c.amplitudes==a.amplitudes && c.keys==a.keys
        @test all(isempty,c.invalid_reasons)
        @test all(!isempty,load_patches(partial,"res",BASE).invalid_reasons)
        out3=joinpath(dir,"partial_scores.tsv"); write_scores(out3,cv_scores(c,SUPPORT))
        @test read(out3)==read(out1)
        cli=joinpath(dir,"cli.tsv")
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_empirical_fisher_native.jl")) --patches $partial --prefix res --config $(CONFIGS["support"]) --out $cli`)
        @test read(cli)==read(out1)
        # Center and all disk positions remain mandatory; never drop keys.
        rows[1][header[GRID.center_index+2]]="NA"
        rows[2][header[first(GRID.disk_indices)+2]]="Inf"
        invalid=joinpath(dir,"invalid.tsv"); write_table(invalid,header,rows)
        d=load_patches(invalid,"res",SUPPORT)
        @test d.keys==a.keys
        @test !isempty(d.invalid_reasons[1]) && !isempty(d.invalid_reasons[2])
        @test all(isnan,d.X[1:2,:])
        @test d.X[3:end,:]==a.X[3:end,:]
        incomplete=joinpath(dir,"missing_column.tsv")
        write_table(incomplete,filter(!=(header[first(outside)+2]),header),rows)
        @test all(startswith(r,"missing_patch_column:") for r in load_patches(incomplete,"res",SUPPORT).invalid_reasons)
    end
end

@testset "Descriptor CLI: all rows retained, complete/disk-missing/peripheral-missing" begin
    mktempdir() do dir
        shape=exp.(-.5 .* (((U .- .12) ./ .055).^2 .+ ((T .- .04) ./ .07).^2))
        partial=copy(shape); partial[1]=NaN
        bad=copy(shape); bad[41]=NaN
        values=[shape,partial,bad]
        base=[Dict("file"=>"synthetic.sxm","lobe"=>string(i),"amplitude"=>".1") for i in 1:3]
        cols=[@sprintf("bwd_res_p%03d",i) for i in 1:81]
        features,patches,out,cli=joinpath.(dir,("features.tsv","patches.tsv","out.tsv","cli.tsv"))
        write_table(features,["file","lobe","amplitude"],base)
        rows=[merge(base[i],Dict(c=>string(v) for (c,v) in zip(cols,values[i]))) for i in 1:3]
        write_table(patches,vcat(["file","lobe"],cols),rows)
        augment_descriptor(features,patches,out,CONFIGS["support"])
        _,actual=lobe_table(out)
        @test length(actual)==3
        for i in 1:3
            value,reason=desc(values[i]); row=actual[("synthetic.sxm",i)]
            @test row["descriptor_reason"]==reason
            @test isfinite(value) ? parse(Float64,row["patch_u_asym_reconstructed"])==value : row["patch_u_asym_reconstructed"]=="NA"
        end
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_reconstructed_descriptor.jl")) --features $features --patches $patches --config $(CONFIGS["support"]) --out $cli`)
        @test read(out)==read(cli)
        forbidden=joinpath(dir,"forbidden.tsv"); write(forbidden,"file\tlobe\texpected_N\nsynthetic.sxm\t1\t6\n")
        @test_throws ErrorException lobe_table(forbidden)
    end
end

@testset "Final Gaussian score algebra, unchanged fit and independent CLI" begin
    mktempdir() do dir
        rng=MersenneTwister(821)
        path=joinpath(dir,"features.tsv"); names=["f$i" for i in 1:8]
        open(path,"w") do io
            println(io,join(vcat(["file","lobe","amplitude"],names),'\t'))
            for i in 1:80
                amp=.05+.1rand(rng); values=randn(rng,8); values[1]=amp
                i==80 && (values[2]=NaN)
                println(io,join(vcat(["synthetic_$(cld(i,8)).sxm",string(mod1(i,8)),string(amp)],string.(values)),'\t'))
            end
        end
        records=GMM._load_records(path)
        opts=Dict(); diagnostics=Dict()
        for name in ("control","volume")
            args=["--config",CONFIGS[name],"--features",path,"--out",joinpath(dir,name*".tsv"),
                "--seeds","2","--selftrain","2","--interactions","--view","synthetic="*join(names,',')]
            opt=GMM._parse_cli(args); opts[name]=opt
            @test opt.final_score==(name=="control" ? "mahalanobis" : "gaussian_density")
            @test opt.covariance_mode=="ridge" && opt.covariance_ridge==1e-6
            diag=[]; p=GMM._view_probability(records,names,opt;diagnostics=diag); diagnostics[name]=diag
            @test length(p)==80 && count(isfinite,p)==79
            @test isequal(p,GMM._view_probability(records,names,opt))
            expected=joinpath(dir,name*"_expected.tsv")
            GMM._write_predictions(expected,records,p,Int.(isfinite.(p)))
            @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args`)
            @test read(expected)==read(opt.out_tsv)
            _,rows=lobe_table(opt.out_tsv)
            @test length(rows)==80 && rows[("synthetic_10.sxm",8)]["predicted"]=="?"
        end
        @test diagnostics["control"] == diagnostics["volume"] # fitted parameters/members/weights exactly equal
        @test_throws ErrorException GMM._parse_cli(["--config",CONFIGS["volume"],"--features",path,"--selftrain","0"])
        for _ in 1:12
            p=4; A=randn(rng,p,p); C=A*A'+1e-6I; x=randn(rng,p); mu=randn(rng,p); w=.37
            @test isposdef(Symmetric(C))
            guard=C+1e-8I; d=x-mu
            expected=log(w)-.5*(p*log(2π)+logdet(guard)+dot(d,guard\d))
            a=GMM._final_component_score(x,mu,C,w,opts["control"])
            b=GMM._final_component_score(x,mu,C,w,opts["volume"])
            # LU and Cholesky can disagree by O(eps*cond) for nearly singular
            # matrices. This is a test error bound, not a covariance parameter.
            @test b ≈ expected rtol=max(5e-11,8eps(Float64)*cond(guard)) atol=1e-10
            well=C+.1I; wellguard=well+1e-8I
            well_expected=log(w)-.5*(p*log(2π)+logdet(wellguard)+dot(d,wellguard\d))
            @test GMM._final_component_score(x,mu,well,w,opts["volume"]) ≈ well_expected rtol=5e-13 atol=1e-11
            @test b-a ≈ -.5*(p*log(2π)+logdet(guard)) rtol=5e-11 atol=1e-10
            L=cholesky(Symmetric(C)+1e-8I).L; z=L\(x.-mu)
            @test a == log(w)-.5dot(z,z) # bit-identical legacy expression
            # Equal volumes leave the relative score unchanged.
            other_mu=randn(rng,p)
            a2=GMM._final_component_score(x,other_mu,C,1-w,opts["control"])
            b2=GMM._final_component_score(x,other_mu,C,1-w,opts["volume"])
            @test a-a2 ≈ b-b2 rtol=5e-11 atol=1e-10
        end
        # At the common center: the old rule prefers the heavier broad cluster;
        # normalized Gaussian density prefers the compact cluster. No labels.
        cs=[Matrix{Float64}(I,2,2),100Matrix{Float64}(I,2,2)]
        weights=[.4,.6]
        scores(name)=[GMM._final_component_score(zeros(2),zeros(2),cs[c],weights[c],opts[name]) for c in 1:2]
        @test argmax(scores("control"))==2
        @test argmax(scores("volume"))==1
    end
end
