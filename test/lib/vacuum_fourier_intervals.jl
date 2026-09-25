module VacuumFourierIntervals
using IntervalArithmetic
const IA=IntervalArithmetic
const FI=Interval{Float64}
I(x)=interval(Float64,x)
I(a,b)=interval(Float64,a,b)
valid(x)=isguaranteed(x) && decoration(x)==com && isfinite(inf(x)) && isfinite(sup(x))
signcertain(x)=inf(x)>0 ? 1 : sup(x)<0 ? -1 : 0

"Group the FULL real Gamma series into interval-enclosed vertical harmonics."
function column(w,xy_fraction,L_nm,weights,order)
    length(xy_fraction)==2 && all(valid,xy_fraction) || error("Invalid lateral coordinate enclosures")
    valid(L_nm) && inf(L_nm)>0 || error("Invalid cell height")
    2<=order<=12 || error("Unsupported Taylor order")
    size(w.coeff,2)==length(weights) && all(>(0),weights) || error("Invalid state weights")
    kmax=maximum(abs,w.miller[3,:]); nb=length(weights)
    a=fill(I(0),kmax+1,nb); b=copy(a)
    for band in 1:nb; a[1,band]=-I(real(w.coeff[w.zero,band])); end
    twopi=I(2)*I(pi)
    for g in axes(w.miller,2)
        m=w.miller[3,g]; k=abs(m)+1
        phase=twopi*(I(w.miller[1,g])*xy_fraction[1]+I(w.miller[2,g])*xy_fraction[2])
        sn,cs=sin(phase),cos(phase)
        for band in 1:nb
            cr=I(real(w.coeff[g,band])); ci=I(imag(w.coeff[g,band]))
            a[k,band]+=I(2)*(cr*cs-ci*sn)
            m==0 || (b[k,band]-=I(2sign(m))*(cr*sn+ci*cs))
        end
    end
    omega=[twopi*I(k)/L_nm for k in 0:kmax]
    make_series(a,b,omega,I.(weights),I(w.volume),order)
end

function make_series(a,b,omega,weights,volume,order)
    size(a)==size(b)==(length(omega),length(weights)) || error("Incompatible series")
    all(valid,a) && all(valid,b) && all(valid,omega) && all(valid,weights) && valid(volume) || error("Unguaranteed series")
    inf(volume)>0 && all(x->inf(x)>0,weights) || error("Invalid positive scale")
    all(x->inf(x)>=0,omega) && isequal_interval(first(omega),I(0)) || error("Invalid harmonic frequencies")
    # Global derivative magnitudes bound the Taylor remainder, not the local
    # value itself. Keep all harmonics, including very small coefficients.
    remainder=fill(I(0),2,length(weights))
    for band in eachindex(weights), k in eachindex(omega), d in 1:2
        remainder[d,band]+=(I(mag(a[k,band]))+I(mag(b[k,band])))*omega[k]^(order+d)
    end
    (;a,b,omega,weights,volume,order,remainder)
end

"Direct enclosed scalar density; no FFT, interpolation or polynomial truncation."
function density(c,z)
    rho=I(0)
    sn=sin.(c.omega.*z); cs=cos.(c.omega.*z)
    for band in eachindex(c.weights)
        psi=I(0)
        for k in eachindex(c.omega)
            psi+=c.a[k,band]*cs[k]+c.b[k,band]*sn[k]
        end
        rho+=c.weights[band]*psi^2
    end
    rho/c.volume
end

"Enclose density and its first derivative everywhere in [lo,hi]."
function enclosure(c,lo,hi)
    box=I(lo,hi); center=I(mid(box)); delta=box-center; h=I(mag(delta))
    order=c.order; rho=drho=I(0)
    sn=sin.(c.omega.*center); cs=cos.(c.omega.*center)
    for band in eachindex(c.weights)
        jet=fill(I(0),order+2)
        for k in eachindex(c.omega)
            t0=c.a[k,band]*cs[k]+c.b[k,band]*sn[k]
            t1=-c.a[k,band]*sn[k]+c.b[k,band]*cs[k]
            rotations=(t0,t1,-t0,-t1); power=I(1)
            for j in 0:order+1
                jet[j+1]+=rotations[1+j%4]*power
                power*=c.omega[k]
            end
        end
        psi=jet[order+1]/I(factorial(order))
        dpsi=jet[order+2]/I(factorial(order))
        for j in order-1:-1:0
            psi=psi*delta+jet[j+1]/I(factorial(j))
            dpsi=dpsi*delta+jet[j+2]/I(factorial(j))
        end
        rem0=c.remainder[1,band]*h^(order+1)/I(factorial(order+1))
        rem1=c.remainder[2,band]*h^(order+1)/I(factorial(order+1))
        psi+=I(-sup(rem0),sup(rem0)); dpsi+=I(-sup(rem1),sup(rem1))
        rho+=c.weights[band]*psi^2
        drho+=I(2)*c.weights[band]*psi*dpsi
    end
    rho/=c.volume; drho/=c.volume
    valid(rho) && valid(drho) || error("Invalid Taylor enclosure")
    (;rho,drho)
end

function refine(c,lo,hi,iso,width,maxiter)
    box=I(lo,hi)
    for iteration in 1:maxiter
        diam(box)<=width && return (;box,converged=true,iterations=iteration-1)
        e=enclosure(c,inf(box),sup(box)); signcertain(e.drho)!=0 || error("Lost monotonic proof")
        center=I(mid(box)); newton=center-(density(c,center)-I(iso))/e.drho
        next=intersect_interval(box,newton)
        isempty_interval(next) && error("Existing root excluded by interval Newton")
        box=next
    end
    (;box,converged=diam(box)<=width,iterations=maxiter)
end

"Complete interval subdivision; inconclusive leaves are preserved, never skipped."
function crossings(c,lo,hi,iso;root_width,max_depth,max_nodes,newton_steps)
    all(isfinite,(lo,hi,iso,root_width)) && lo<hi && root_width>0 || error("Invalid root domain")
    max_depth>=0 && max_nodes>0 && newton_steps>0 || error("Invalid numerical limits")
    stack=[(lo,hi,0)]; leaves=NamedTuple[]; nodes=0
    while !isempty(stack)
        a,b,depth=pop!(stack); nodes+=1
        status=:unresolved; root_lo=root_hi=NaN; direction=0
        rho_lo=rho_hi=derivative_lo=derivative_hi=NaN
        if nodes>max_nodes
            status=:node_limit
        else
            e=enclosure(c,a,b); rho_lo,rho_hi=bounds(e.rho); derivative_lo,derivative_hi=bounds(e.drho)
            if iso<rho_lo || iso>rho_hi
                status=:excluded
            else
                direction=signcertain(e.drho)
                if direction!=0
                    sa=signcertain(density(c,I(a))-I(iso)); sb=signcertain(density(c,I(b))-I(iso))
                    if sa!=0 && sb!=0
                        if sa==sb
                            status=:excluded_monotone
                        else
                            r=refine(c,a,b,iso,root_width,newton_steps)
                            root_lo,root_hi=bounds(r.box)
                            status=r.converged ? :unique : :precision_limit
                        end
                    end
                end
            end
            if status==:unresolved && depth<max_depth
                middle=mid(I(a,b))
                if a<middle<b
                    push!(stack,(middle,b,depth+1)); push!(stack,(a,middle,depth+1))
                    continue
                end
            end
        end
        push!(leaves,(;lo=a,hi=b,depth,status,direction,root_lo,root_hi,
            rho_lo,rho_hi,derivative_lo,derivative_hi))
    end
    sort!(leaves;by=x->x.lo)
    roots=filter(x->x.status==:unique,leaves)
    unresolved=count(x->x.status in (:unresolved,:node_limit,:precision_limit),leaves)
    valid_root=unresolved==0 && length(roots)==1 && only(roots).direction==-1
    (;leaves,roots,unresolved,nodes,valid_root)
end

end
