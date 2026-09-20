using Test, TOML
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))

const CONTROL = joinpath(ROOT, "config", "unit_assignment_affine_residual.toml")
const SIGNAL_CONFIGS = ("signed_mold", "affine_fisher")

@testset "Two independent settings, no composition or voting change" begin
    base = load_config(CONTROL)
    @test base["model"]["mold_margin_mode"] == "absolute_cost_margin"
    @test base["preprocessing"]["fisher_patch_projection"] == "none"
    @test base["model"]["fisher_projection_zero_l1"] == 1e-12
    for name in SIGNAL_CONFIGS
        cfg = load_config(joinpath(ROOT,"config","unit_assignment_"*name*".toml"))
        @test cfg["model"]["name"] == "cc_soft_"*name*"_v1"
        cfg["model"]["name"] = base["model"]["name"]
        if name == "signed_mold"
            @test cfg["model"]["mold_margin_mode"] == "signed_cost_difference"
            cfg["model"]["mold_margin_mode"] = base["model"]["mold_margin_mode"]
        else
            @test cfg["preprocessing"]["fisher_patch_projection"] == "affine_disk"
            cfg["preprocessing"]["fisher_patch_projection"] = base["preprocessing"]["fisher_patch_projection"]
        end
        @test cfg == base
    end
    mktempdir() do dir
        for (section,key) in (("model","mold_margin_mode"),("preprocessing","fisher_patch_projection")),
            value in (nothing,"guess",true)
            cfg=deepcopy(base)
            value === nothing ? delete!(cfg[section],key) : (cfg[section][key]=value)
            path=joinpath(dir,section*key*string(value)*".toml")
            open(io->TOML.print(io,cfg),path,"w")
            @test_throws ErrorException load_config(path)
        end
    end
end

@testset "Signed mold contrast uses physical costs, not decoded labels" begin
    positive=Dict("cost_GlcN"=>"0.4","cost_GlcNAc"=>"-0.2","cost_margin"=>"0.60000001",
                  "predicted"=>"0","physical_label"=>"GlcN")
    @test mold_margin(positive,"absolute_cost_margin") == 0.60000001 # preserve saved rounding
    @test mold_margin(positive,"signed_cost_difference") ≈ 0.6
    negative=merge(positive,Dict("cost_GlcN"=>"-0.2","cost_GlcNAc"=>"0.4"))
    @test mold_margin(negative,"signed_cost_difference") ≈ -0.6
    @test mold_margin(negative,"absolute_cost_margin") == mold_margin(positive,"absolute_cost_margin")
    for label in ("0","1","?","arbitrary")
        changed=merge(positive,Dict("predicted"=>label,"physical_label"=>label))
        @test mold_margin(changed,"signed_cost_difference") == mold_margin(positive,"signed_cost_difference")
    end
    @test mold_margin(Dict("cost_GlcN"=>"0.2","cost_GlcNAc"=>"0.2"),"signed_cost_difference") == 0
    for bad in ("NA","NaN","Inf","-Inf","broken"), key in ("cost_GlcN","cost_GlcNAc")
        @test isnan(mold_margin(merge(positive,Dict(key=>bad)),"signed_cost_difference"))
    end
    @test_throws ErrorException mold_margin(positive,"guess")
    @test_throws KeyError mold_margin(Dict("cost_margin"=>"1"),"signed_cost_difference")
end

@testset "Whole-key feature join and benchmark input rejection" begin
    mktempdir() do dir
        ids=[Dict("file"=>"synthetic.sxm","lobe"=>string(i)) for i in 1:3]
        function table(name,columns,extra)
            path=joinpath(dir,name*".tsv")
            write_table(path,vcat(["file","lobe"],columns),[merge(ids[i],extra[i]) for i in 1:3])
            path
        end
        base=table("base",["patch_u_asym_reconstructed"],[Dict("patch_u_asym_reconstructed"=>x) for x in ("0.4","-0.2","NA")])
        split=table("split",["skew_ratio"],[Dict("skew_ratio"=>"1.0") for _ in 1:3])
        fisher=table("fisher",["score"],[Dict("score"=>x) for x in ("1.2","-0.7","NA")])
        cols=["cost_GlcN","cost_GlcNAc","cost_margin","predicted","physical_label"]
        data=[Dict(zip(cols,values)) for values in (("0.4","-0.2","0.60000001","1","GlcNAc"),
              ("-0.2","0.4","0.60000001","0","GlcN"),("NA","0.2","NA","?","?"))]
        fwd=table("fwd",cols,data)
        bwd=table("bwd",cols,[merge(r,Dict("cost_GlcN"=>r["cost_GlcNAc"],"cost_GlcNAc"=>r["cost_GlcN"])) for r in data])
        outputs=Dict{String,Any}()
        for mode in ("absolute_cost_margin","signed_cost_difference")
            path=joinpath(dir,mode*".tsv")
            join_predictor_features(base,split,fwd,bwd,fisher,path;margin_mode=mode)
            header,rows=lobe_table(path)
            outputs[mode]=(header,rows,path)
            @test Set(keys(rows)) == Set(("synthetic.sxm",i) for i in 1:3)
            @test rows[("synthetic.sxm",3)]["mold_cc_fwd"] == "NaN"
            @test rows[("synthetic.sxm",3)]["mold_cc_bwd"] == "NaN"
            @test rows[("synthetic.sxm",3)]["patch_u_asym_reconstructed"] == "NaN"
            @test rows[("synthetic.sxm",3)]["emp_fisher"] == "NaN"
            # Decoded template labels are irrelevant to both feature modes.
            other=table(mode*"_poison",cols,[merge(r,Dict("predicted"=>"opposite","physical_label"=>"unused")) for r in data])
            repeat=joinpath(dir,mode*"_repeat.tsv")
            join_predictor_features(base,split,other,bwd,fisher,repeat;margin_mode=mode)
            @test read(path) == read(repeat)
        end
        ha,a,_=outputs["absolute_cost_margin"]; hb,b,_=outputs["signed_cost_difference"]
        @test ha == hb
        stable=setdiff(ha,["mold_cc_fwd","mold_cc_bwd"])
        @test all(a[k][c]==b[k][c] for k in keys(a) for c in stable)
        @test [b[("synthetic.sxm",i)]["mold_cc_fwd"] for i in 1:2] == ["0.6","-0.6"]
        @test [b[("synthetic.sxm",i)]["mold_cc_bwd"] for i in 1:2] == ["-0.6","0.6"]
        for field in ("truth","expected_N","control_sequence","benchmark_class_count"),
            mode in ("absolute_cost_margin","signed_cost_difference")
            forbidden=table(field*mode,[field],[Dict(field=>"forbidden") for _ in 1:3])
            out=joinpath(dir,field*mode*"_out.tsv")
            @test_throws ErrorException join_predictor_features(forbidden,split,fwd,bwd,fisher,out;margin_mode=mode)
            @test !ispath(out)
        end
        for flag in ("--truth","--expected-N","--reference","--manifest","--control-sequence")
            @test_throws ErrorException parse_options([flag,"forbidden"])
        end
    end
end
