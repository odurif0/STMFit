using Test, LinearAlgebra
include(joinpath(@__DIR__,"lib/vacuum_fourier_intervals.jl"))
const F=VacuumFourierIntervals
const IA=F.IA
const I=F.I

@testset "Enclosed Gamma regrouping and analytic density" begin
    miller=[0 1 0 2; 0 0 1 -1; 0 1 -2 0]
    coeff=ComplexF64[.3 .2; .1+.05im .2; -.04im .1im; .03-.02im .04+.01im]
    w=(;miller,coeff,zero=1,volume=1000.)
    fx,fy=.213,.371; L=2.; weights=[2.,2.]
    c=F.column(w,[I(fx),I(fy)],I(L),weights,12)
    @test all(F.valid,c.a) && all(F.valid,c.b) && all(F.valid,c.omega)
    @test size(c.a)==(3,2)
    # High-precision independent full +/-G evaluation checks cancellation,
    # negative Gz, z-independent lateral terms and the G=0 multiplicity.
    function truth(z;derivative=false)
        setprecision(256) do
            total=BigFloat(0)
            for band in 1:2
                psi=BigFloat(real(coeff[1,band]))
                dpsi=BigFloat(0)
                for g in 2:4
                    phase=2BigFloat(pi)*(miller[1,g]*BigFloat(fx)+miller[2,g]*BigFloat(fy)+miller[3,g]*BigFloat(z)/BigFloat(L))
                    psi+=2real(Complex{BigFloat}(coeff[g,band])*cis(phase))
                    dpsi+=2real(Complex{BigFloat}(coeff[g,band])*cis(phase)*im*2BigFloat(pi)*miller[3,g]/BigFloat(L))
                end
                total+=BigFloat(weights[band])*(derivative ? 2psi*dpsi : psi^2)/BigFloat(w.volume)
            end
            total
        end
    end
    for z in range(.1,1.7;length=101)
        box=F.density(c,I(z)); exact=truth(z)
        @test BigFloat(IA.inf(box))<=exact<=BigFloat(IA.sup(box))
    end
    for (lo,hi) in ((.3,.31),(.8,.85),(.2,1.1))
        e=F.enclosure(c,lo,hi)
        @test F.valid(e.rho) && F.valid(e.drho)
        for z in range(lo,hi;length=21)
            exact=truth(z)
            @test BigFloat(IA.inf(e.rho))<=exact<=BigFloat(IA.sup(e.rho))
            slope=truth(z;derivative=true)
            @test BigFloat(IA.inf(e.drho))<=slope<=BigFloat(IA.sup(e.drho))
        end
    end
    @test_throws ErrorException F.column(w,[I(fx),I(fy)],I(L),[2.,-1.],12)
end

@testset "Continuous root exclusion, uniqueness and explicit ambiguity" begin
    # rho=(1+0.6 cos(2*pi*z))^2: positive, but globally nonmonotonic.
    a=reshape(I.([1.,.6]),2,1); b=fill(I(0),2,1)
    c=F.make_series(a,b,[I(0),I(2)*I(pi)],[I(1)],I(1),12)
    opts=(;root_width=1e-7,max_depth=24,max_nodes=20000,newton_steps=32)
    r=F.crossings(c,.05,.49,.7;opts...)
    @test r.valid_root && r.unresolved==0 && length(r.roots)==1
    expected=acos((sqrt(.7)-1)/.6)/(2pi)
    @test only(r.roots).root_lo<=expected<=only(r.roots).root_hi
    @test only(r.roots).root_hi-only(r.roots).root_lo<=opts.root_width
    @test first(r.leaves).lo==.05 && last(r.leaves).hi==.49
    @test all(r.leaves[i].hi==r.leaves[i+1].lo for i in 1:length(r.leaves)-1)
    allroots=F.crossings(c,.05,.95,.7;opts...)
    @test !allroots.valid_root && allroots.unresolved==0 && length(allroots.roots)==2
    @test sort([x.direction for x in allroots.roots])==[-1,1]
    absent=F.crossings(c,.05,.49,3.;opts...)
    @test !absent.valid_root && isempty(absent.roots) && absent.unresolved==0
    # The exact minimum is deliberately uncertified, never called a crossing.
    tangent_c=F.make_series(reshape(I.([1.,.5]),2,1),b,c.omega,c.weights,c.volume,12)
    tangent=F.crossings(tangent_c,.45,.55,.25;opts...)
    @test !tangent.valid_root
    limited=F.crossings(c,.05,.95,.7;root_width=1e-7,max_depth=0,max_nodes=1,newton_steps=1)
    @test limited.unresolved>0 && !limited.valid_root
    @test IA.isequal_interval(F.density(c,I(.31)),F.density(c,I(.31)))
    @test_throws ErrorException F.crossings(c,1.,0.,.7;opts...)
end
