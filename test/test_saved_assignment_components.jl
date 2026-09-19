using Test, TOML, SHA
include("export_saved_assignment_components.jl")
const SC = SavedAssignmentComponents

function fixture(dir)
    h="file\tlobe\tpredicted\tconfidence\tprobability_1\tinvalid_reason\tmodel\n"
    rows = (
        ["a.sxm\t1\t1\t0\t0.5\tok\tbase", "a.sxm\t2\t1\t0\t0.5\tok\tbase", "a.sxm\t3\t?\t0\tNA\tunavailable_gmm\tbase", "b.sxm\t1\t?\t0\tNA\tunavailable_kmeans\tbase", "b.sxm\t2\t1\t0.7\t0.85\tok\tbase", "c.sxm\t1\t0\t0.5\t0.25\tok\tbase"],
        ["a.sxm\t1\t0\t1\t0\tok\tkm", "a.sxm\t2\t1\t0.5\t0.75\tok\tkm", "a.sxm\t3\t0\t1\t0\tok\tkm", "b.sxm\t1\t?\t0\tNA\tmissing_view\tkm", "b.sxm\t2\t1\t0.6\t0.8\tok\tkm", "c.sxm\t1\t0\t0\t0.50000000\tok\tkm"],
        ["a.sxm\t1\t1\t1\t1\tok\tgmm", "a.sxm\t2\t0\t0.5\t0.25\tok\tgmm", "a.sxm\t3\t?\t0\tNA\tmissing_view\tgmm", "b.sxm\t1\t1\t1\t1\tok\tgmm", "b.sxm\t2\t1\t0.8\t0.9\tok\tgmm", "c.sxm\t1\t0\t1\t0\tok\tgmm"])
    paths=[joinpath(dir,n*".tsv") for n in ("control","km","gmm")]
    for (path,rs) in zip(paths,rows); write(path,h*join(rs,"\n")*"\n"); end
    return paths
end

@testset "Saved component endpoints preserve decisions and common validity" begin
    mktempdir() do dir
        paths=fixture(dir); before=SC.sha.(paths); out=joinpath(dir,"out")
        metadata=SC.export_components(paths...,out)
        @test SC.sha.(paths)==before
        @test read(joinpath(out,"control.tsv"))==read(paths[1])
        @test Set(readdir(out))==Set(["control.tsv","metadata.toml",(n*".tsv" for n in SC.ENDPOINTS)...])
        @test metadata==TOML.parsefile(joinpath(out,"metadata.toml"))
        control=last(SC.lobe_table(paths[1])); km=last(SC.lobe_table(paths[2])); gm=last(SC.lobe_table(paths[3]))
        for (name,source) in zip(SC.ENDPOINTS,(km,gm))
            dest=last(SC.lobe_table(joinpath(out,name*".tsv")))
            @test Set(keys(dest))==Set(keys(control))
            @test metadata["counts"][name]==Dict("rows"=>6,"available"=>4,"unavailable"=>2)
            @test SC.sha(joinpath(out,name*".tsv"))==metadata["outputs"][name*".tsv"]
            for key in keys(control)
                chosen=control[key]["predicted"]=="?" ? control[key] : source[key]
                for field in SC.REQUIRED; @test dest[key][field]==chosen[field]; end
                @test dest[key]["model"]==name
            end
        end
        copied=last(SC.lobe_table(joinpath(out,SC.ENDPOINTS[1]*".tsv")))
        @test copied[("c.sxm",1)]["predicted"]=="0"
        @test copied[("c.sxm",1)]["probability_1"]=="0.50000000"
        @test_throws ErrorException SC.export_components(paths...,out)
        @test SC.sha.(paths)==before
        link=joinpath(dir,"dangling"); symlink(joinpath(dir,"absent"),link)
        @test_throws ErrorException SC.export_components(paths...,link)
        @test islink(link)
        @test_throws ErrorException SC.export_components(paths...,joinpath(dir,"missing_parent","out"))
    end
end

@testset "Reject malformed input and altered availability without writing outputs" begin
    changes = (
        (2, s->s*split(s,'\n')[2]*"\n"), # duplicate key
        (2, s->replace(s,"a.sxm\t2\t"=>"a.sxm\t8\t")), # noncontiguous
        (3, s->replace(s,"c.sxm\t1\t"=>"d.sxm\t1\t")), # wrong keys
        (3, s->replace(s,"a.sxm\t1\t1\t1\t1"=>"a.sxm\t1\t2\t1\t1")),
        (2, s->replace(s,"a.sxm\t1\t0\t1\t0"=>"a.sxm\t1\t0\t1\t1.2")),
        (2, s->replace(s,"a.sxm\t1\t0\t1\t0"=>"a.sxm\t1\t0\tNaN\t0")),
        (1, s->replace(s,"a.sxm\t1\t1\t0\t0.5\tok"=>"a.sxm\t1\t?\t0\tNA\tunavailable")),
        (1, s->replace(s,"a.sxm\t3\t?\t0\tNA"=>"a.sxm\t3\t0\t0\t0")),
        (1, s->replace(s,"a.sxm\t3\t?\t0\tNA"=>"a.sxm\t3\t?\t0\t0.5")),
        (2, s->replace(s,"model"=>"ground_truth")),
        (2, s->replace(s,"model"=>"expected_N")),
        (3, s->first(split(s,'\n'))*"\n"))
    for (index,transform) in changes
        mktempdir() do dir
            paths=fixture(dir); write(paths[index],transform(read(paths[index],String)))
            out=joinpath(dir,"out")
            @test_throws ErrorException SC.export_components(paths...,out)
            @test !ispath(out)
        end
    end
end

@testset "CLI help and arity are metadata-only" begin
    @test SC.main(["--help"])===nothing
    @test_throws ErrorException SC.main(String[])
end
