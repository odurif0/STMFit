module QENativePlot
using LinearAlgebra
number(x)=parse(Float64,replace(x,'D'=>'E','d'=>'e'))

"QE 7.4.1 Modules/plot_io.f90; only explicit orthorhombic cells are supported."
function readplot(file)
    open(file) do io
        title=readline(io); dims=parse.(Int,split(readline(io)))
        length(dims)==8 || error("Bad native dimensions")
        nxp,nyp,nzp,nx,ny,nz,nat,ntyp=dims
        all(>(0),dims) && nxp>=nx && nyp>=ny && nzp>=nz || error("Invalid native dimensions")
        line=split(readline(io)); length(line)==7 || error("Bad lattice header")
        parse(Int,line[1])==0 || error("Only explicit cells supported")
        celldm=number.(line[2:7]); celldm[1]>0 && all(iszero,celldm[2:end]) || error("Unexpected celldm")
        at=hcat([number.(split(readline(io))) for _ in 1:3]...)
        size(at)==(3,3) && all(isfinite,at) && at==diagm(diag(at)) && all(>(0),diag(at)) || error("Nonorthorhombic cell")
        basis=split(readline(io)); length(basis)==4 || error("Bad cutoff record")
        plot_num=parse(Int,basis[4]); cutoffs=number.(basis[1:3])
        species=NamedTuple[]
        for j in 1:ntyp
            row=split(readline(io)); length(row)==3 && parse(Int,row[1])==j || error("Species record")
            push!(species,(;element=String(row[2]),valence=number(row[3])))
        end
        atoms=NamedTuple[]
        for j in 1:nat
            row=split(readline(io)); length(row)==5 && parse(Int,row[1])==j || error("Atom record")
            tau=number.(row[2:4]); type=parse(Int,row[5])
            all(isfinite,tau) && 1<=type<=ntyp || error("Invalid atom")
            push!(atoms,(;tau,position_bohr=celldm[1]*tau,element=species[type].element))
        end
        # QE writes nxp*nyp*nz (not nxp*nyp*nzp); x is the fastest index.
        values=Vector{Float64}(undef,nxp*nyp*nz); k=0
        for line in eachline(io), token in split(line)
            k+=1; k<=length(values) || error("Extra native values")
            values[k]=number(token)
        end
        k==length(values) && all(isfinite,values) || error("Incomplete/nonfinite native payload")
        storage=reshape(values,nxp,nyp,nz)
        (;title,dims=(nx,ny,nz),padded_dims=(nxp,nyp,nzp),cell_bohr=celldm[1]*at,
            alat=celldm[1],at,plot_num,cutoffs,species,atoms,storage,grid=view(storage,1:nx,1:ny,:))
    end
end

"Half the decimal printing quantum; a numerical export bound, not physical error."
halfquantum(x,digits)=iszero(x) ? 0. : .5*10.0^(floor(log10(abs(x)))-digits+1)
end
