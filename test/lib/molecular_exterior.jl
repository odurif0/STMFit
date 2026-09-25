module MolecularExterior
const Q=Rational{BigInt}
point(p)=(Q(p[1]),Q(p[2]))
sub(a,b)=(a[1]-b[1],a[2]-b[2])
dot2(a,b)=a[1]*b[1]+a[2]*b[2]
orient(a,b,c)=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])

"Andrew chains, exact for the supplied finite Float64 coordinates."
function hull(points)
    all(p->length(p)==2 && all(isfinite,p),points) || error("Invalid planar points")
    ps=sort!(unique(point.(points))); length(ps)>=3 || error("Insufficient centers")
    function chain(xs)
        result=eltype(ps)[]
        for p in xs
            while length(result)>=2 && orient(result[end-1],result[end],p)<=0; pop!(result); end
            push!(result,p)
        end
        result
    end
    lo=chain(ps); hi=chain(reverse(ps)); result=vcat(lo[1:end-1],hi[1:end-1])
    length(result)>=3 || error("Degenerate molecular footprint")
    result
end

function distance2(p,poly)
    all(orient(poly[j],poly[mod1(j+1,length(poly))],p)>=0 for j in eachindex(poly)) && return Q(0)
    minimum(eachindex(poly)) do j
        a=poly[j]; b=poly[mod1(j+1,length(poly))]; edge=sub(b,a); v=sub(p,a)
        t=clamp(dot2(v,edge)/dot2(edge,edge),Q(0),Q(1))
        r=(v[1]-t*edge[1],v[2]-t*edge[2]); dot2(r,r)
    end
end

"Envelope of the complete trimer; Cu radii NEVER set its molecular dilation."
function footprint(geometry,radii)
    atoms=filter(a->a.z!=29,geometry.atoms); isempty(atoms) && error("No molecular atoms")
    coords=[a.position_nm[1:2] for a in atoms]; poly=hull(coords)
    radius=Q(maximum(radii[a.z] for a in atoms)); radius>0 || error("Invalid molecular radius")
    cell=point(geometry.cell_nm[1:2]); all(>(0),cell) || error("Invalid lateral cell")
    # This experiment does not infer how to unwrap a molecule crossing the cell.
    all(minimum(p[k] for p in poly)-radius>0 && maximum(p[k] for p in poly)+radius<cell[k] for k in 1:2) ||
        error("Offset molecular bbox crosses a cell boundary; no automatic reimaging")
    (;poly,radius,cell,atoms)
end

function clearance(xy,f)
    p=point(xy); p=(mod(p[1],f.cell[1]),mod(p[2],f.cell[2]))
    d2=minimum(distance2((p[1]+ix*f.cell[1],p[2]+iy*f.cell[2]),f.poly) for ix in -1:1 for iy in -1:1)
    comparison=sign(d2-f.radius^2)
    # Rounded distance is report-only. Exact squared comparisons define support.
    # Inside the undilated hull this margin saturates at -radius; it is not
    # the full signed penetration depth into the rounded envelope.
    (;distance2=d2,exterior=comparison>0,boundary=comparison==0,
        margin_nm=sqrt(Float64(d2))-Float64(f.radius))
end
end
