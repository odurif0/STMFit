using Test
include(joinpath(@__DIR__,"report_observed_assignment_effects.jl"))
const E=ObservedAssignmentEffects
function vote_fixture()
    g=Dict(("a",1)=>Dict("predicted"=>"1","probability_1"=>"0.8"),
           ("b",1)=>Dict("predicted"=>"?","probability_1"=>"NA"))
    k=Dict(("a",1)=>Dict("predicted"=>"0","probability_1"=>"0.4"),
           ("b",1)=>Dict("predicted"=>"0","probability_1"=>"0.2"))
    p=Dict(("a",1)=>Dict("predicted"=>"1","probability_1"=>"0.6","confidence"=>"0.2","model"=>"fixture"),
           ("b",1)=>Dict("predicted"=>"?","probability_1"=>"NA","confidence"=>"NA","model"=>"fixture"))
    g,k,p
end
@testset "Saved-only vote arithmetic and exact input groups" begin
    g,k,p=vote_fixture()
    @test E.verify_vote(g,k,p,.5,"fixture")==2
    for (key,value) in (("model","wrong"),("predicted","0"),("probability_1","0.7"),("confidence","0.1"))
        bad=deepcopy(p); bad[("a",1)][key]=value
        @test_throws ErrorException E.verify_vote(g,k,bad,.5,"fixture")
    end
    bad=deepcopy(p); bad[("b",1)]["predicted"]="0"
    @test_throws ErrorException E.verify_vote(g,k,bad,.5,"fixture")
    bad=deepcopy(g); bad[("a",1)]["probability_1"]="1.2"
    @test_throws ErrorException E.verify_vote(bad,k,p,.5,"fixture")
    @test_throws ErrorException E.verify_vote(Dict(),k,p,.5,"fixture")
    q=deepcopy(p); q[("a",1)]["confidence"]="0.3"; q[("b",1)]["predicted"]="1"
    @test E.differences(p,q,Set(["a"]))==(lobes=1,changed_lobes=1,changed_scans=1)
    @test E.differences(p,q,Set(["a"]);prediction=true)==(lobes=1,changed_lobes=0,changed_scans=0)
    @test E.differences(p,q,Set(["a","b"]);prediction=true)==(lobes=2,changed_lobes=1,changed_scans=1)
    @test E.differences(p,q,Set{String}())==(lobes=0,changed_lobes=0,changed_scans=0)
    @test_throws ErrorException E.differences(p,Dict(),Set(["a"]))
end
