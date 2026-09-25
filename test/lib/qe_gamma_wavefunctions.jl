module QEGammaWavefunctions
using LinearAlgebra

# QE 7.4.1 Modules/io_base.f90: little-endian, four-byte Fortran record
# markers. Deliberately reject other formats, spinors and non-Gamma states.
function record(io, expected)
    Base.ENDIAN_BOM == 0x04030201 || error("Only little-endian hosts supported")
    n = read(io,Int32)
    n == expected || error("Unexpected Fortran record size: $n != $expected")
    payload = read(io,n)
    length(payload)==n && read(io,Int32)==n || error("Truncated/mismatched record")
    IOBuffer(payload)
end

function header(io)
    r=record(io,44)
    ik=read(r,Int32); k=read!(r,zeros(3)); spin=read(r,Int32)
    gamma=read(r,Int32); scale=read(r,Float64)
    ik==1 && all(iszero,k) && spin==1 && gamma in (-1,1) && scale==1.0 ||
        error("Not an unscaled scalar Gamma PW checkpoint")
    r=record(io,16); ngw,ng,npol,nbnd=Int.(read!(r,zeros(Int32,4)))
    ngw==ng && ng>0 && npol==1 && nbnd>0 || error("Unsupported WFC dimensions")
    reciprocal=reshape(read!(record(io,72),zeros(9)),3,3)
    all(isfinite,reciprocal) && reciprocal==diagm(diag(reciprocal)) &&
        all(>(0),diag(reciprocal)) || error("Expected the audited orthorhombic cell")
    # write_collected_wfc stores b_i in inverse Bohr, INCLUDING 2*pi.
    cell=2pi*inv(transpose(reciprocal))
    (;ng,nbnd,reciprocal,cell,volume=det(cell))
end

function read_wfc(file, selected)
    issorted(selected) && allunique(selected) && !isempty(selected) || error("Invalid band list")
    open(file) do io
        h=header(io)
        all(b->1<=b<=h.nbnd,selected) || error("Band outside checkpoint")
        miller=reshape(read!(record(io,12h.ng),zeros(Int32,3h.ng)),3,h.ng)
        triples=[Tuple(miller[:,i]) for i in 1:h.ng]
        allunique(triples) || error("Duplicate reciprocal vector")
        zero=only(findall(==((0,0,0)),triples))
        lookup=Set(triples)
        all(t->t==(0,0,0) || !((-t[1],-t[2],-t[3]) in lookup),triples) ||
            error("Both members of a Gamma conjugate pair were stored")
        coeff=zeros(ComplexF64,h.ng,length(selected))
        for band in 1:h.nbnd
            column=findfirst(==(band),selected)
            if column===nothing
                read(io,Int32)==16h.ng || error("Invalid skipped band record")
                skip(io,16h.ng)
                read(io,Int32)==16h.ng || error("Invalid skipped band trailer")
            else
                read!(record(io,16h.ng),view(coeff,:,column))
            end
        end
        eof(io) || error("Unexpected trailing WFC records")
        all(isfinite,coeff) && all(isreal,coeff[zero,:]) || error("Invalid real Gamma state")
        merge(h,(;miller,coeff,zero,selected=collect(selected)))
    end
end

"Direct reciprocal-space norm; no density grid or augmentation normalization."
function parseval(w,weights)
    length(weights)==size(w.coeff,2) || error("Weight count mismatch")
    sum(weights[b]*(2sum(abs2,view(w.coeff,:,b))-abs2(w.coeff[w.zero,b]))
        for b in eachindex(weights))
end

"Evaluate the real Gamma Fourier series at fractional coordinates, then square."
function density(w, fractions, weights; block_points, parallel=true)
    size(fractions,1)==3 && all(isfinite,fractions) || error("Invalid coordinates")
    length(weights)==size(w.coeff,2) && all(>(0),weights) || error("Invalid weights")
    block_points>0 || error("Invalid block size")
    BLAS.get_num_threads()==1 || error("Use one BLAS thread per independent point block")
    n=size(fractions,2); result=zeros(n)
    mill=Float64.(w.miller)
    real_coeff=2real.(w.coeff); imag_coeff=-2imag.(w.coeff)
    function block(first_point)
        indices=first_point:min(n,first_point+block_points-1)
        phases=transpose(mill)*view(fractions,:,indices)
        cs=similar(phases); sn=similar(phases)
        @inbounds for i in eachindex(phases)
            sn[i],cs[i]=sincos(2pi*phases[i])
        end
        psi=transpose(real_coeff)*cs+transpose(imag_coeff)*sn
        psi .-= real.(w.coeff[w.zero,:])
        result[indices]=vec(sum(weights.*abs2.(psi);dims=1))/w.volume
        nothing
    end
    starts=1:block_points:n
    if parallel
        Threads.@threads for j in eachindex(starts); block(starts[j]); end
    else
        for i in starts; block(i); end
    end
    result
end
end
