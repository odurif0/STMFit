module VacuumCrossings
using LinearAlgebra

"Strict, unique descending crossing of a piecewise-linear profile."
function crossing(z, v, iso)
    length(z)==length(v)>=2 || error("Profile needs at least two paired samples")
    all(isfinite,z) && all(>(0),diff(z)) && all(isfinite,v) && isfinite(iso) || error("Invalid profile")
    roots=Float64[]; directions=Int[]; plateau=false
    for i in 1:length(v)-1
        a,b=v[i],v[i+1]
        a==b==iso && (plateau=true)
        if min(a,b)<iso<max(a,b)
            push!(roots,z[i]+(iso-a)/(b-a)*(z[i+1]-z[i]))
            push!(directions,a>b ? -1 : 1)
        end
    end
    for i in eachindex(v)
        v[i]==iso || continue
        push!(roots,z[i])
        push!(directions,1<i<length(v) && v[i-1]>iso>v[i+1] ? -1 :
            (1<i<length(v) && v[i-1]<iso<v[i+1] ? 1 : 0))
    end
    status=plateau ? :plateau : isempty(roots) ? :absent : length(roots)>1 ? :multiple :
        (v[1]==iso || v[end]==iso) ? :boundary : directions[1]==0 ? :tangent :
        directions[1]==1 ? :ascending : :found
    (;status,z=status==:found ? only(roots) : NaN,roots=length(roots),
        descending=count(==(-1),directions),ascending=count(==(1),directions))
end

"All maximal OPEN isovalue intervals giving exactly one descending crossing."
function intervals(z,v)
    crossing(z,v,0.) # validate, without selecting an isovalue
    events=Dict{Float64,Tuple{Int,Int}}(Float64(x)=>(0,0) for x in v)
    for (a,b) in zip(v[1:end-1],v[2:end])
        a==b && continue
        lo,hi=minmax(a,b); delta=a>b ? (1,0) : (0,1)
        events[lo]=events[lo].+delta; events[hi]=events[hi].-delta
    end
    knots=sort(collect(keys(events))); d=a=0; result=Tuple{Float64,Float64}[]
    for j in 1:length(knots)-1
        lo,hi=knots[j:j+1]; dd,da=events[lo]; d+=dd; a+=da
        d==1 && a==0 || continue
        if !isempty(result) && result[end][2]==lo && crossing(z,v,lo).status==:found
            result[end]=(result[end][1],hi)
        else
            push!(result,(lo,hi))
        end
    end
    result
end

function intersect_intervals(a,b)
    result=Tuple{Float64,Float64}[]; i=j=1
    while i<=length(a) && j<=length(b)
        lo=max(a[i][1],b[j][1]); hi=min(a[i][2],b[j][2])
        lo<hi && push!(result,(lo,hi))
        if a[i][2]<b[j][2]; i+=1; else; j+=1; end
    end
    result
end

"Horizontal projected chain frame; vertical motion is substrate +z."
function lateral_points(frame,half,step)
    half>0 && step>0 && isinteger(2half/step) || error("Invalid lateral grid")
    t=[frame.t_axis[1],frame.t_axis[2],0.]; norm(t)>0 || error("Vertical tangent")
    t/=norm(t); u=[-t[2],t[1],0.]
    grid=range(-half,half;length=1+Int(2half/step))
    reduce(vcat,[permutedims(frame.origin_nm+x*t+y*u) for y in grid for x in grid])
end

"Geometry only: bounding PAW planes and native z knots strictly inside the gap."
function domains(geometry,radii,nz)
    cell=geometry.cell_nm
    all(a->0<=a.position_nm[3]<cell[3],geometry.atoms) || error("Atoms need canonical z cell")
    lower=maximum(a.position_nm[3]+radii[a.z] for a in geometry.atoms)
    upper=minimum(a.position_nm[3]+cell[3]-radii[a.z] for a in geometry.atoms)
    0<lower<upper<cell[3] || error("No single above-molecule PAW-free gap")
    midpoint=(lower+upper)/2
    dz=cell[3]/nz
    knots=[k*dz for k in 0:nz-1 if lower<k*dz<upper]
    length(knots)>=2 || error("Insufficient full-gap knots")
    half=findall(<(midpoint),knots)
    length(half)>=2 || error("Insufficient molecular-half knots")
    (;lower,upper,midpoint,knots,half)
end

"Bilinear lateral interpolation at exact native vertical knots, no z rounding."
function profiles(cube,cell_nm,xy,z)
    nx,ny,nz=cube.dims
    out=zeros(length(z),size(xy,1))
    qz=z./cell_nm[3].*nz
    kz=round.(Int,qz)
    all(abs.(qz-kz).<=4eps.(max.(1.,abs.(qz)))) || error("Expected exact native z knots")
    all(0 .<= kz .< nz) || error("Vertical support outside cube")
    for (p,point) in enumerate(eachrow(xy))
        q=point[1:2]./cell_nm[1:2].*[nx,ny]; lo=floor.(Int,q); f=q-lo
        all(0 .<= lo) && all(lo.+1 .< [nx,ny]) || error("Lateral support outside cube")
        for dx in 0:1,dy in 0:1
            weight=(dx==0 ? 1-f[1] : f[1])*(dy==0 ? 1-f[2] : f[2])
            for (j,k) in enumerate(kz)
                out[j,p]+=weight*cube.grid[k+1,lo[2]+dy+1,lo[1]+dx+1]
            end
        end
    end
    out
end

end
