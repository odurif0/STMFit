module SubstrateReferenceGeometry
using LinearAlgebra

"Distance to periodically repeated molecular PAW disks in the lateral plane."
function projection_margin(xy,geometry,radii)
    molecules=filter(a->a.z!=29,geometry.atoms)
    isempty(molecules) && error("No molecular atoms")
    minimum(molecules) do atom
        d=xy-atom.position_nm[1:2]
        d-=geometry.cell_nm[1:2].*round.(d./geometry.cell_nm[1:2])
        norm(d)-radii[atom.z]
    end
end

"Complete cell-center grid; eligibility depends on geometry of BOTH states."
function grid(geometries,radii,dims)
    length(geometries)==length(radii)==2 || error("Two paired geometries required")
    length(dims)==2 && all(>(0),dims) || error("Invalid lateral grid")
    all(g->g.cell_nm[1:2]==first(geometries).cell_nm[1:2],geometries) || error("Different lateral cells")
    rows=NamedTuple[]
    for iy in 0:dims[2]-1, ix in 0:dims[1]-1
        fx,fy=(ix+.5)/dims[1],(iy+.5)/dims[2]
        xy=[fx,fy].*first(geometries).cell_nm[1:2]
        margins=[projection_margin(xy,g,r) for (g,r) in zip(geometries,radii)]
        push!(rows,(;pixel=length(rows)+1,ix,iy,fx,fy,x_nm=xy[1],y_nm=xy[2],
            glcn_margin_nm=margins[1],glcnac_margin_nm=margins[2],eligible=all(>(0),margins)))
    end
    rows
end

"Conservative Cu bounding planes, plus the unchanged old molecular-half domain."
function domains(geometry,radii,nz,old)
    cu=filter(a->a.z==29,geometry.atoms)
    isempty(cu) && error("No substrate atoms")
    all(a->0<=a.position_nm[3]<geometry.cell_nm[3],geometry.atoms) || error("Noncanonical z")
    lower=maximum(a.position_nm[3]+radii[a.z] for a in cu)
    upper=minimum(a.position_nm[3]+geometry.cell_nm[3]-radii[a.z] for a in cu)
    0<lower<upper<geometry.cell_nm[3] || error("No Cu PAW gap")
    midpoint=(lower+upper)/2
    knots=[k*geometry.cell_nm[3]/nz for k in 0:nz-1 if lower<k*geometry.cell_nm[3]/nz<upper]
    half=filter(<(midpoint),knots)
    min(length(knots),length(half))>=2 || error("Insufficient native knots")
    [(;domain="full_cu_gap",lo=first(knots),hi=last(knots),lower_paw_nm=lower,upper_paw_nm=upper,midpoint_nm=midpoint),
     (;domain="cu_half_gap",lo=first(half),hi=last(half),lower_paw_nm=lower,upper_paw_nm=upper,midpoint_nm=midpoint),
     (;domain="saved_molecular_half",lo=first(old.knots),hi=old.knots[last(old.half)],
        lower_paw_nm=old.lower,upper_paw_nm=old.upper,midpoint_nm=old.midpoint)]
end

"Conservative margin for EVERY point on a vertical segment, not just endpoints."
function segment_margin(xy,lo,hi,geometry,radii)
    lo<hi && hi-lo<geometry.cell_nm[3] || error("Invalid segment")
    minimum(geometry.atoms) do atom
        d=xy-atom.position_nm[1:2]
        d-=geometry.cell_nm[1:2].*round.(d./geometry.cell_nm[1:2])
        dz=minimum(max(lo-(atom.position_nm[3]+k*geometry.cell_nm[3]),
            atom.position_nm[3]+k*geometry.cell_nm[3]-hi,0.) for k in -1:1)
        hypot(norm(d),dz)-radii[atom.z]
    end
end
end
