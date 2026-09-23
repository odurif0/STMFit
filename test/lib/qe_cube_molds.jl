module QECubeMolds

include(joinpath(@__DIR__,"cc_mold_native.jl"))
using .CCMoldNative, TOML

function settings(path)
    expected=Dict("model"=>Dict("cube_order"=>"qe_last_axis_fast"),
        "selection"=>Dict("surface_policy"=>"unchanged_assignment_config"),
        "preprocessing"=>Dict("cube_units"=>"bohr"))
    parsed=TOML.parsefile(path)
    parsed==expected || throw(ArgumentError("Only the declared QE-order input correction is supported"))
    parsed
end

"""Adapt QE's third-index-fast token stream to the unchanged first-index-fast
sampler's in-memory layout. Geometry and numerical interpolation are unchanged.
This is an input permutation, not a rotation/transposition of physical axes.
"""
function read_qe_cube(path::AbstractString)
    c=CCMoldNative.read_cube(path)
    n1,n2,n3=c.dims
    ordered=vec(permutedims(reshape(c.values,n3,n2,n1),(3,2,1)))
    CCMoldNative.Cube(c.origin_nm,c.axes_nm,c.dims,ordered)
end

end
