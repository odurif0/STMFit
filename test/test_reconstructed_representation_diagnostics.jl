using Test, LinearAlgebra, Statistics, Printf, SHA, Random
include(joinpath(@__DIR__, "diagnose_reconstructed_representation.jl"))
using .ReconstructedRepresentationCLI.ReconstructedRepresentationDiagnostics
const RD = ReconstructedRepresentationCLI.ReconstructedRepresentationDiagnostics
const RU = RD.ReconstructedUnitAssignment
const EF = RD.EmpiricalFisherNative
const CONFIG = joinpath(@__DIR__, "..", "config", "unit_assignment_reconstructed.toml")
const ATOL = 1e-12
const ZERO = 1e-12
const C9 = collect(-4:4) .* 0.08
const T9, U9 = RD.grid_vectors(C9)

@testset "Mass identity, exact production signs/support and transformations" begin
    samples = [U9, -U9, T9, ones(81), exp.(-12 .* (T9.^2 + U9.^2)),
               U9 .+ 0.3T9 .+ 0.07, sin.(23U9) .* cos.(7T9)]
    for p in samples
        a = patch_audit(p, C9; zero_l1=ZERO, identity_atol=ATOL)
        d, reason = RU.transverse_asymmetry(p, C9; zero_l1=ZERO)
        @test a["d_raw"] == d
        @test a["reason"] == reason == "ok"
        @test a["mass_aligned"] + a["mass_opposed"] + a["mass_central"] ≈ a["mass_l1"]
        @test a["fraction_aligned"] - a["fraction_opposed"] ≈ d atol=ATOL
        @test 2a["fraction_opposed"] + a["fraction_central"] ≈ 1-d atol=ATOL
        @test all(abs(a[k]) <= ATOL for k in ("mass_identity_error", "d_identity_error", "one_minus_d_identity_error", "numerator_attribution_error"))
        for scale in (0.25, 8.0)
            @test patch_audit(scale*p, C9; zero_l1=ZERO, identity_atol=ATOL)["d_raw"] ≈ d atol=ATOL
        end
        @test patch_audit(-p, C9; zero_l1=ZERO, identity_atol=ATOL)["d_raw"] ≈ -d atol=ATOL
        flip_t = vcat([reverse(p[(j-1)*9+1:j*9]) for j in 1:9]...)
        flip_u = vcat([p[(j-1)*9+1:j*9] for j in 9:-1:1]...)
        @test patch_audit(flip_t, C9; zero_l1=ZERO, identity_atol=ATOL)["d_raw"] ≈ d atol=ATOL
        @test patch_audit(flip_u, C9; zero_l1=ZERO, identity_atol=ATOL)["d_raw"] ≈ -d atol=ATOL
    end
    @test patch_audit(U9, C9; zero_l1=ZERO, identity_atol=ATOL)["d_raw"] ≈ 1
    @test patch_audit(-U9, C9; zero_l1=ZERO, identity_atol=ATOL)["d_raw"] ≈ -1
    @test patch_audit(T9, C9; zero_l1=ZERO, identity_atol=ATOL)["d_raw"] ≈ 0 atol=ATOL
    central = zeros(81); central[41] = 5
    a = patch_audit(central, C9; zero_l1=ZERO, identity_atol=ATOL)
    @test a["mass_central"] == 5
    @test a["d_raw"] == 0
    @test a["fraction_central"] == 1
    for p in (zeros(81), fill(1e-15,81))
        a = patch_audit(p, C9; zero_l1=ZERO, identity_atol=ATOL)
        @test a["reason"] == "zero_patch_mass"
        @test isnan(a["d_raw"])
    end
    for value in (NaN, Inf, -Inf)
        p = copy(U9); p[1] = value
        a = patch_audit(p, C9; zero_l1=ZERO, identity_atol=ATOL)
        @test a["reason"] == "nonfinite_patch"
        @test a["finite_pixels"] == 80
        @test isnan(a["gradient_u_per_nm"])
    end
    @test_throws ErrorException patch_audit(ones(80), C9; zero_l1=ZERO, identity_atol=ATOL)
    @test_throws ErrorException patch_audit(U9, C9; zero_l1=-1, identity_atol=ATOL)
    @test_throws ErrorException patch_audit(U9, C9; zero_l1=ZERO, identity_atol=NaN)
    @test_throws ErrorException patch_audit(U9, C9 .+ 0.01; zero_l1=ZERO, identity_atol=ATOL)
end

@testset "Fixed affine projection and diagnostic denominator distinction" begin
    p = 0.4 .+ 0.7T9 .+ 2.3U9
    a = affine_components(p, C9)
    @test a.c ≈ 0.4
    @test a.bt ≈ 0.7
    @test a.bu ≈ 2.3
    @test norm(a.residual) < ATOL
    projection = hcat(ones(81), T9, U9) \ p
    @test [a.c, a.bt, a.bu] ≈ projection
    d = patch_audit(p, C9; zero_l1=ZERO, identity_atol=ATOL)
    @test d["numerator_constant"] ≈ 0 atol=ATOL
    @test d["numerator_t"] ≈ 0 atol=ATOL
    @test d["numerator_u"] ≈ 2.3sum(abs,U9)
    @test d["affine_removed_reason"] == "zero_patch_mass"
    @test isnan(d["d_affine_removed_diagnostic"])
    @test d["gradient_t_energy_fraction"] + d["gradient_u_energy_fraction"] ≈ 1
    # Centered nonlinear remainder prevents a vanishing denominator.
    p2 = p .+ 0.4 .* cos.(21U9) .* sin.(11T9)
    a2 = affine_components(p2, C9)
    @test abs(dot(a2.residual,T9)) < ATOL
    @test abs(dot(a2.residual,U9)) < ATOL
    @test abs(sum(a2.residual)) < ATOL
    d2 = patch_audit(p2,C9;zero_l1=ZERO,identity_atol=ATOL)
    @test sum(d2[k] for k in ("gradient_t_energy_fraction","gradient_u_energy_fraction","affine_residual_energy_fraction")) ≈ 1
    @test d2["d_gradient_removed_diagnostic"] == first(RU.transverse_asymmetry(p2-a2.along_u,C9;zero_l1=ZERO))
    @test d2["u_over_raw_mass"] + d2["affine_residual_over_raw_mass"] ≈ d2["d_raw"] atol=ATOL
    constant = patch_audit(ones(81),C9;zero_l1=ZERO,identity_atol=ATOL)
    @test isnan(constant["gradient_u_energy_fraction"])
    @test constant["gradient_removed_reason"] == "ok"
    @test constant["affine_removed_reason"] == "zero_patch_mass"
end

@testset "Exact native max-reflection identity including common offset" begin
    grid = EF.fisher_grid(EF.load_fisher_config(CONFIG))
    n = length(grid.disk_indices)
    rng = MersenneTwister(4)
    dummy = EF.GMMFit([0.5,0.5],zeros(2,1),[ones(1,1),ones(1,1)],[1,2],true,0,0.0)
    for _ in 1:20
        x,w,mid = randn(rng,n),randn(rng,n),randn(rng,n)
        model = EF.FisherModel(w,mid,zeros(n),ones(n),(0.,1.),dummy,nothing)
        d = response_decomposition(x,w,mid,grid)
        @test d.direct_max ≈ EF.maxmirror_score(x,model,grid) atol=ATOL
        @test d.analytic_max ≈ d.direct_max atol=ATOL
        @test d.shared_centering_offset == -dot(mid,w)
        @test d.even_centered ≈ (d.original+d.mirrored)/2 atol=ATOL
        @test d.odd ≈ (d.original-d.mirrored)/2 atol=ATOL
        @test d.reflection_bonus >= 0
    end
    # Pure even/odd weights on the symmetric sub-support. Shared offset remains.
    x = randn(rng,n); rx = EF.flip_u_disk(x,grid)
    @test !any(iszero,grid.mirror_indices) # saved 17x17 support is fully reflection-closed
    padded = EF.FisherGrid(3, [-1.0,0.0,1.0], [1,2,4,5,6], [0,2,5,4,3], 5)
    px, pw, pm = randn(rng,5), randn(rng,5), randn(rng,5)
    padded_identity = response_decomposition(px,pw,pm,padded)
    @test padded_identity.analytic_max ≈ padded_identity.direct_max atol=ATOL
    @test EF.flip_u_disk(EF.flip_u_disk(px,padded),padded) != px
    symmetric = copy(x)
    for i in eachindex(symmetric)
        grid.mirror_indices[i] == 0 && (symmetric[i] = 0)
    end
    even = (symmetric + EF.flip_u_disk(symmetric,grid))/2
    odd = (symmetric - EF.flip_u_disk(symmetric,grid))/2
    @test EF.flip_u_disk(even,grid) ≈ even
    @test EF.flip_u_disk(odd,grid) ≈ -odd
    for w in (even, odd)
        mid = ones(n)
        d=response_decomposition(x,w,mid,grid)
        @test d.analytic_max ≈ d.direct_max atol=ATOL
    end
    pure=response_decomposition(x,odd,ones(n),grid)
    @test pure.even_uncentered ≈ 0 atol=ATOL
    @test pure.analytic_max ≈ pure.shared_centering_offset+abs(pure.odd) atol=ATOL
    @test_throws ErrorException response_decomposition(x[2:end],odd,ones(n),grid)
    @test_throws ErrorException response_decomposition(fill(NaN,n),odd,ones(n),grid)
    negative = response_decomposition(-ones(n),ones(n),zeros(n),grid)
    @test negative.analytic_max == negative.direct_max == -n
end

@testset "Shared file/cohort gradients are descriptive fixed-basis means" begin
    keys = [("a.sxm",1),("a.sxm",2),("b.sxm",1),("b.sxm",2)]
    patches = Dict(keys[1]=>U9+T9/2, keys[2]=>3U9+T9/2,
                   keys[3]=>-U9+T9/2, keys[4]=>fill(NaN,81))
    diagnostics = Dict(k=>patch_audit(patches[k],C9;zero_l1=ZERO,identity_atol=ATOL) for k in keys)
    RD.add_common!(diagnostics,patches,C9,keys;zero_l1=ZERO)
    @test diagnostics[keys[1]]["common_file_gradient_u_per_nm"] ≈ 2
    @test diagnostics[keys[2]]["common_file_gradient_u_per_nm"] ≈ 2
    @test diagnostics[keys[3]]["common_file_gradient_u_per_nm"] ≈ -1
    @test diagnostics[keys[1]]["common_cohort_gradient_u_per_nm"] ≈ 1
    @test diagnostics[keys[4]]["reason"] == "nonfinite_patch"
    @test isnan(diagnostics[keys[4]]["d_common_file_gradient_removed_diagnostic"])
    for key in keys[1:3]
        r = diagnostics[key]
        expected = first(RU.transverse_asymmetry(patches[key]-r["common_file_gradient_u_per_nm"]*U9,C9;zero_l1=ZERO))
        @test r["d_common_file_gradient_removed_diagnostic"] == expected
    end
end

# Load only the unchanged arithmetic functions for exact source-conformance
# tests. Never include scripts' unconditional main(), imports or fit functions.
module NativePreprocessingFixture
using Statistics
struct LobeRecord
    file::String
    features::Dict{String,Float64}
end
end
for script in ("build_labelfree_gmm_predictions.jl", "build_labelfree_unit_predictions.jl")
    source = read(joinpath(@__DIR__,script),String)
    start = findfirst("function _standardized_matrix",source).start
    finish = findnext("\nfunction ",source,start+1).start
    Base.include_string(NativePreprocessingFixture,source[start:prevind(source,finish)],script)
    @testset "Exact preprocessing from $script" begin
        raw = [10. 0 1; 12. 0 4; 14. NaN 7; 5. 3. NaN; NaN 5. NaN; NaN NaN NaN]
        files = ["a","a","a","b","b","c"]
        features = ["f1","f2","f3"]
        records = [NativePreprocessingFixture.LobeRecord(files[i],Dict(features[j]=>raw[i,j] for j in 1:3)) for i in 1:6]
        for interactions in (false,true)
            native, valid = NativePreprocessingFixture._standardized_matrix(records,features;interactions)
            result=standardize_columns(raw,files;interactions)
            @test isequal(native,result.expanded)
            @test valid == result.valid
            @test result.z[1:3,1] == [-1.,0.,1.]
            @test result.z[1:2,2] == [0.,0.]
            @test result.z[4,1] == 0
            shifted=copy(raw); shifted[1:3,1] .+= 100; shifted[4:5,1] .-= 50
            @test isequal(standardize_columns(shifted,files;interactions).expanded,native)
        end
    end
end
source = read(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl"),String)
start = findfirst("function _negative_moment",source).start
finish = findnext("\nfunction ",source,start+1).start
Base.include_string(NativePreprocessingFixture,source[start:prevind(source,finish)],"native_negative_moment")
@testset "Existing negative moments: no new patch features" begin
    c = collect(range(-1.,1.;length=9))
    for p in (U9, T9, ones(81), zeros(81), [NaN;U9[2:end]])
        native = NativePreprocessingFixture._negative_moment(p,c)
        actual = negative_moments(p)
        @test isequal(actual.com_t,native.com_t)
        @test isequal(actual.diag45,native.diag45)
    end
end

function fixture(dir)
    paths=Dict(name=>joinpath(dir,name*".tsv") for name in ("patches-bwd9","patches-bwd17","patches-fwd17","descriptor","predictor","fisher"))
    ks=[(f,i) for f in ("a.sxm","b.sxm") for i in 1:4]
    descriptor=Dict{String,String}[]
    for (name,prefix,side,step) in (("patches-bwd9","bwd_res",9,.08),("patches-bwd17","bwd_res",17,.04),("patches-fwd17","res",17,.04))
        cols=[@sprintf("%s_p%03d",prefix,i) for i in 1:side^2]
        coords=collect(-(side÷2):side÷2).*step
        t,u=RD.grid_vectors(coords)
        rows=Dict{String,String}[]
        for (file,lobe) in ks
            p=(lobe==1 ? copy(u) : lobe==2 ? 0.3t+2u : lobe==3 ? zeros(length(u)) : copy(t))
            file=="b.sxm" && lobe==4 && (p[1]=NaN)
            row=Dict("file"=>file,"lobe"=>string(lobe))
            for (col,value) in zip(cols,p)
                row[col]=isfinite(value) ? @sprintf("%.7g",value) : "NA"
            end
            push!(rows,row)
            if side==9
                saved=[RD.number(row[col]) for col in cols]
                d,reason=RU.transverse_asymmetry(saved,coords;zero_l1=ZERO)
                push!(descriptor,Dict("file"=>file,"lobe"=>string(lobe),"patch_u_asym_reconstructed"=>RD.text(d),"descriptor_reason"=>reason))
            end
        end
        RU.write_table(paths[name],vcat(["file","lobe"],cols),rows)
    end
    RU.write_table(paths["descriptor"],["file","lobe","patch_u_asym_reconstructed","descriptor_reason"],descriptor)
    features=Dict{String,String}[]; fisher=Dict{String,String}[]
    for (i,key) in enumerate(ks)
        row=copy(descriptor[i])
        for (j,feature) in enumerate(vcat(RD.BASE,["mold_cc_fwd","mold_cc_bwd","split_log_skew"]))
            row[feature]=string(sin(i+j))
        end
        score=i==8 ? "NA" : @sprintf("%.6f",10i/3)
        row["emp_fisher"]=i==8 ? "NaN" : @sprintf("%.8g",-parse(Float64,score))
        push!(features,row)
        push!(fisher,Dict("file"=>key[1],"lobe"=>string(key[2]),"score"=>score,"invalid_reason"=>i==8 ? "invalid_patch_value:res_p001" : ""))
    end
    RU.write_table(paths["predictor"],vcat(["file","lobe"],RD.GMM_FEATURES,["split_log_skew","descriptor_reason"]),features)
    RU.write_table(paths["fisher"],["file","lobe","score","invalid_reason"],fisher)
    return paths
end

@testset "Saved-table end-to-end, support, CLI, collision and source immutability" begin
    mktempdir() do dir
        paths=fixture(dir)
        hashes=Dict(k=>bytes2hex(sha256(read(v))) for (k,v) in paths)
        out=joinpath(dir,"audit")
        result=run_audit(paths;production_config=CONFIG,outdir=out,identity_atol=ATOL,descriptor_atol=ATOL)
        @test result == (keys=8,files=2,descriptor_matches=8)
        _,rows=RU.read_table(joinpath(out,"patch_diagnostics.tsv"))
        @test length(rows)==24
        @test length(Set((r["file"],r["lobe"]) for r in rows))==8
        @test count(r->r["reason"]=="nonfinite_patch",rows)==3
        @test count(r->r["reason"]=="zero_patch_mass",rows)==6
        @test all(r["d_raw"]=="NA" for r in rows if r["reason"]!="ok")
        _,scores=RU.read_table(joinpath(out,"saved_fisher_audit.tsv"))
        @test scores[end]["saved_score"]=="NA"
        @test scores[end]["saved_invalid_reason"]=="invalid_patch_value:res_p001"
        @test all(r["decomposition_status"]=="unavailable_saved_weights_and_mid_not_exported" for r in scores)
        @test occursin("not identifiable",read(joinpath(out,"report.md"),String))
        @test all(bytes2hex(sha256(read(v)))==hashes[k] for (k,v) in paths)
        @test_throws ErrorException run_audit(paths;production_config=CONFIG,outdir=out,identity_atol=ATOL,descriptor_atol=ATOL)
        empty=joinpath(dir,"empty");mkdir(empty)
        @test_throws ErrorException run_audit(paths;production_config=CONFIG,outdir=empty,identity_atol=ATOL,descriptor_atol=ATOL)
        dangling=joinpath(dir,"dangling");symlink(joinpath(dir,"missing"),dangling)
        @test_throws ErrorException run_audit(paths;production_config=CONFIG,outdir=dangling,identity_atol=ATOL,descriptor_atol=ATOL)
        # Missing keys and duplicates are fatal, not silently filtered.
        original=read(paths["fisher"],String)
        lines=readlines(paths["fisher"])
        write(paths["fisher"],join(lines[1:end-1],"\n")*"\n")
        @test_throws ErrorException RD.load_inputs(paths,CONFIG)
        write(paths["fisher"],original*lines[2]*"\n")
        @test_throws ErrorException RD.load_inputs(paths,CONFIG)
        write(paths["fisher"],original)
        # A changed saved descriptor is rejected without any new output directory.
        original_descriptor=read(paths["descriptor"],String)
        changed_lines = split(chomp(original_descriptor), '\n')
        changed_fields = split(changed_lines[2], '\t'; keepempty=true)
        changed_fields[3] = "0.5"
        changed_lines[2] = join(changed_fields, '\t')
        changed_descriptor = join(changed_lines, '\n') * "\n"
        @test changed_descriptor != original_descriptor
        write(paths["descriptor"],changed_descriptor)
        badout=joinpath(dir,"bad")
        @test_throws ErrorException run_audit(paths;production_config=CONFIG,outdir=badout,identity_atol=ATOL,descriptor_atol=ATOL)
        @test !ispath(badout)
        write(paths["descriptor"],original_descriptor)
        # Pixel NA is an invalid retained row; a missing full pixel column is a
        # schema error, matching the production descriptor's boundary.
        original_patch=read(paths["patches-bwd9"],String)
        write(paths["patches-bwd9"],replace(original_patch,"bwd_res_p001"=>"absent_p001";count=1))
        @test_throws ErrorException RD.load_inputs(paths,CONFIG)
        write(paths["patches-bwd9"],original_patch)
        args=reduce(vcat,[["--"*k,v] for (k,v) in sort(collect(paths))])
        append!(args,["--production-config",CONFIG,"--identity-atol","1e-12","--descriptor-atol","1e-12","--outdir",joinpath(dir,"cli")])
        command=`$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(@__DIR__)) $(joinpath(@__DIR__,"diagnose_reconstructed_representation.jl")) $args`
        @test occursin("8 keys",read(command,String))
        @test_throws ErrorException ReconstructedRepresentationCLI.main(["--outdir","x"])
        @test_throws ErrorException ReconstructedRepresentationCLI.main(["--truth","x"])
        @test_throws ErrorException ReconstructedRepresentationCLI.main(["--fisher","x","--fisher","y"])
    end
end
