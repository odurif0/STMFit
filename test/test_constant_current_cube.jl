#!/usr/bin/env julia

# Failing-first tests for the constant-current cube extraction library.
#
# Covers analytic planar/exponential cube recovery, vacuum-side first crossing,
# linear interpolation, absent/ambiguous/nonfinite/low-z cases, masks and
# diagnostics, and the mean-height isovalue policy. QE plot_num=5 values are a
# discrete |psi|^2 sum, never amperes; the isovalue is fixed by the 0.50 nm
# mean-height policy above the Cu reference, not by a benchmark grade.

using Test
using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "lib", "constant_current_cube.jl"))
using .ConstantCurrentCube: CubeGrid, FrameRef, CrossingResult, SurfaceResult,
                            first_vacuum_crossing, constant_current_surface,
                            isovalue_for_mean_height,
                            DEFAULT_ISOVALUE_SCAN_INTERVALS

# --------------------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------------------

"`LDOS = slope * z + intercept` for every lateral column (z in nm, LDOS finite)."
function _planar_cube(; nx::Int=4, ny::Int=4, nz::Int=20, step::Real=0.1,
                     slope::Real=-2.0, intercept::Real=4.0)
    origin = [0.0, 0.0, 0.0]
    axes = Float64[step 0.0 0.0; 0.0 step 0.0; 0.0 0.0 step]
    n = (nx, ny, nz)
    values = Float64[]
    for ix in 0:(nx - 1), iy in 0:(ny - 1), iz in 0:(nz - 1)
        push!(values, slope * (iz * step) + intercept)
    end
    return CubeGrid(origin, axes, n, values)
end

"`LDOS = amplitude * exp(-z/decay)`; physically decreasing with height."
function _exponential_cube(; nx::Int=4, ny::Int=4, nz::Int=40, step::Real=0.05,
                           amplitude::Real=10.0, decay::Real=0.20)
    origin = [0.0, 0.0, 0.0]
    axes = Float64[step 0.0 0.0; 0.0 step 0.0; 0.0 0.0 step]
    n = (nx, ny, nz)
    values = Float64[]
    for ix in 0:(nx - 1), iy in 0:(ny - 1), iz in 0:(nz - 1)
        push!(values, amplitude * exp(-(iz * step) / decay))
    end
    return CubeGrid(origin, axes, n, values)
end

"Cube with a lateral gaussian bump so neighbouring columns differ."
function _bump_cube(; nx::Int=24, ny::Int=24, nz::Int=24, step::Real=0.08,
                   height_amp::Real=5.0, height_decay::Real=0.20,
                   lateral_amp::Real=1.0, lateral_sigma::Real=0.20)
    origin = [0.0, 0.0, 0.0]
    axes = Float64[step 0.0 0.0; 0.0 step 0.0; 0.0 0.0 step]
    n = (nx, ny, nz)
    cx = (nx - 1) * step / 2
    cy = (ny - 1) * step / 2
    values = Float64[]
    for ix in 0:(nx - 1), iy in 0:(ny - 1), iz in 0:(nz - 1)
        x = ix * step
        y = iy * step
        z = iz * step
        lateral = lateral_amp * exp(-((x - cx)^2 + (y - cy)^2) / (2 * lateral_sigma^2))
        push!(values, (height_amp + lateral) * exp(-z / height_decay))
    end
    return CubeGrid(origin, axes, n, values)
end

function _frame_at_origin()
    return FrameRef([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
end

"Build a cube from one equally sampled z profile per lateral x column."
function _column_cube(columns::Vector{Vector{Float64}}; step::Float64=1.0)
    nz = length(first(columns))
    all(length(column) == nz for column in columns) || error("column lengths must match")
    axes = Float64[step 0.0 0.0; 0.0 step 0.0; 0.0 0.0 step]
    return CubeGrid([0.0, 0.0, 0.0], axes, (length(columns), 1, nz), vcat(columns...))
end

function _error_message(f::Function)
    try
        f()
    catch err
        return sprint(showerror, err)
    end
    return ""
end

# --------------------------------------------------------------------------------------
# CubeGrid validation
# --------------------------------------------------------------------------------------

@testset "CubeGrid validation" begin
    @test_nowarn CubeGrid([0.0, 0.0, 0.0], diagm([0.1, 0.1, 0.1]), (3, 3, 3), zeros(27))
    @test_throws ArgumentError CubeGrid([0.0, 0.0], diagm([0.1, 0.1, 0.1]), (3, 3, 3), zeros(27))
    @test_throws ArgumentError CubeGrid([0.0, 0.0, 0.0], zeros(2, 3), (3, 3, 3), zeros(27))
    @test_throws ArgumentError CubeGrid([0.0, 0.0, 0.0], diagm([0.1, 0.1, 0.1]), (3, 3), zeros(27))
    @test_throws ArgumentError CubeGrid([0.0, 0.0, 0.0], diagm([0.1, 0.1, 0.1]), (3, 3, 3), zeros(26))
    @test_throws ArgumentError CubeGrid([0.0, 0.0, 0.0], diagm([0.1, 0.1, 0.1]), (0, 3, 3), zeros(27))
    @test_throws ArgumentError CubeGrid([NaN, 0.0, 0.0], diagm([0.1, 0.1, 0.1]), (3, 3, 3), zeros(27))
    @test_throws ArgumentError CubeGrid([0.0, 0.0, 0.0], diagm([0.1, 0.1, 0.1]), (3, 3, 3), [zeros(26); NaN])
    @test_throws ArgumentError CubeGrid([0.0, 0.0, 0.0], zeros(3, 3), (3, 3, 3), zeros(27))
end

@testset "FrameRef validation" begin
    @test_nowarn FrameRef([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    @test_throws ArgumentError FrameRef([0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    @test_throws ArgumentError FrameRef([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    @test_throws ArgumentError FrameRef([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
    @test_throws ArgumentError FrameRef([NaN, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
end

# --------------------------------------------------------------------------------------
# first_vacuum_crossing
# --------------------------------------------------------------------------------------

@testset "first_vacuum_crossing: planar column exact" begin
    # LDOS = -2z + 4 on a 0.1 nm grid: crossing of iso=1.0 lands exactly at z=1.5.
    z = collect(0.0:0.1:1.9)
    v = [-2.0 * zi + 4.0 for zi in z]
    res = first_vacuum_crossing(z, v, 1.0)
    @test res.status === :found
    @test res.z ≈ 1.5 atol = 1e-12
    @test res.n_brackets == 1
    @test res.bracket == 17   # iz=17 is z=1.6, bracketing pair (1.6, 1.5) -> upper index 17
end

@testset "first_vacuum_crossing: planar column interpolated" begin
    # Same profile, iso=0.5 -> analytic z = (4-0.5)/2 = 1.75, half-grid offset.
    z = collect(0.0:0.1:1.9)
    v = [-2.0 * zi + 4.0 for zi in z]
    res = first_vacuum_crossing(z, v, 0.5)
    @test res.status === :found
    @test res.z ≈ 1.75 atol = 1e-12
    @test res.n_brackets == 1
end

@testset "first_vacuum_crossing: exponential within one grid step" begin
    step = 0.05
    decay = 0.20
    amp = 10.0
    z = collect(0.0:step:1.95)
    v = [amp * exp(-zi / decay) for zi in z]
    iso = 1.0
    analytic_z = -decay * log(iso / amp)
    res = first_vacuum_crossing(z, v, iso)
    @test res.status === :found
    @test res.z ≈ analytic_z atol = step   # within one interpolated-grid tolerance
    @test res.n_brackets == 1
end

@testset "first_vacuum_crossing: absent (all below isovalue)" begin
    z = collect(0.0:0.1:1.0)
    v = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
    res = first_vacuum_crossing(z, v, 5.0)
    @test res.status === :absent
    @test isnan(res.z)
    @test res.n_brackets == 0
    @test res.bracket == 0
end

@testset "first_vacuum_crossing: absent (all above isovalue)" begin
    z = collect(0.0:0.1:1.0)
    v = 10.0 .* ones(length(z))
    res = first_vacuum_crossing(z, v, 1.0)
    @test res.status === :absent
    @test isnan(res.z)
end

@testset "first_vacuum_crossing: ambiguous multiple crossings" begin
    # Oscillating column: [10, 0.5, 2, 0.3, 8] at z=[0, 0.1, 0.2, 0.3, 0.4], iso=1.
    z = [0.0, 0.1, 0.2, 0.3, 0.4]
    v = [10.0, 0.5, 2.0, 0.3, 8.0]
    res = first_vacuum_crossing(z, v, 1.0)
    @test res.status === :ambiguous
    @test res.n_brackets == 4
    @test res.bracket == 5          # first bracket closest to vacuum is the pair (z=0.4, z=0.3)
    @test isnan(res.z)              # ambiguity is rejected, never zero-filled
end

@testset "first_vacuum_crossing: nonfinite column" begin
    z = collect(0.0:0.1:0.5)
    v = [4.0, 3.0, NaN, 1.0, 0.0, -1.0]
    res = first_vacuum_crossing(z, v, 1.5)
    @test res.status === :nonfinite
    @test isnan(res.z)
end

@testset "first_vacuum_crossing: insufficient z resolution" begin
    z = [0.0]
    v = [1.0]
    res = first_vacuum_crossing(z, v, 0.5)
    @test res.status === :insufficient_z
    @test isnan(res.z)
end

@testset "first_vacuum_crossing: non-finite isovalue rejected" begin
    z = collect(0.0:0.1:1.0)
    v = copy(z)
    @test_throws ArgumentError first_vacuum_crossing(z, v, NaN)
    @test_throws ArgumentError first_vacuum_crossing(z, v, Inf)
end

@testset "first_vacuum_crossing: mismatched column lengths" begin
    @test_throws ArgumentError first_vacuum_crossing([0.0, 0.1], [1.0], 0.5)
end

@testset "first_vacuum_crossing: non-ascending z rejected" begin
    @test_throws ArgumentError first_vacuum_crossing([0.1, 0.0], [1.0, 2.0], 0.5)
end

@testset "first_vacuum_crossing: exact touching counts as a bracket" begin
    # v[3] sits exactly on iso; the function must still report a unique crossing
    # when the bracket resolves unambiguously via interpolation.
    z = [0.0, 0.1, 0.2, 0.3, 0.4]
    v = [3.0, 2.0, 1.0, 0.0, -1.0]
    res = first_vacuum_crossing(z, v, 1.0)
    @test res.status === :found
    @test res.z ≈ 0.2 atol = 1e-12
    @test res.n_brackets == 1
end

# --------------------------------------------------------------------------------------
# constant_current_surface
# --------------------------------------------------------------------------------------

@testset "constant_current_surface: planar cube uniform crossing" begin
    cube = _planar_cube()
    res = constant_current_surface(cube, 1.0)
    @test all(res.valid)
    @test all(res.height .≈ 1.5)
    @test all(res.status .=== :found)
    @test res.counts[:found] == prod(cube.n[1:2])
    @test get(res.counts, :absent, 0) == 0
    @test get(res.counts, :ambiguous, 0) == 0
end

@testset "constant_current_surface: exponential cube within tolerance" begin
    step = 0.05
    decay = 0.20
    amp = 10.0
    cube = _exponential_cube(; step=step, decay=decay, amplitude=amp)
    iso = 1.0
    analytic_z = -decay * log(iso / amp)
    res = constant_current_surface(cube, iso)
    @test all(res.valid)
    @test all(isapprox.(res.height, analytic_z; atol=step))
end

@testset "constant_current_surface: diagnostics reflect mixed columns" begin
    # Build a small cube whose left half is absent (vacuum) and right half crosses.
    nx, ny, nz = 4, 1, 20
    step = 0.1
    origin = [0.0, 0.0, 0.0]
    axes = Float64[step 0.0 0.0; 0.0 step 0.0; 0.0 0.0 step]
    values = Float64[]
    for ix in 0:(nx - 1), iy in 0:(ny - 1), iz in 0:(nz - 1)
        if ix >= 2
            push!(values, -2.0 * (iz * step) + 4.0)   # crossing exists for iso in (0, 4)
        else
            push!(values, 10.0)                        # always above iso=1 -> absent
        end
    end
    cube = CubeGrid(origin, axes, (nx, ny, nz), values)
    res = constant_current_surface(cube, 1.0)
    @test res.counts[:found] == 2
    @test res.counts[:absent] == 2
    @test all(res.valid[3:4, 1])
    @test all(.!res.valid[1:2, 1])
    @test all(res.height[3:4, 1] .≈ 1.5)
    @test all(isnan.(res.height[1:2, 1]))
end

@testset "constant_current_surface: non-finite cube values rejected" begin
    nx, ny, nz = 2, 2, 3
    step = 0.1
    origin = [0.0, 0.0, 0.0]
    axes = Float64[step 0.0 0.0; 0.0 step 0.0; 0.0 0.0 step]
    values = Float64[]
    for ix in 0:(nx - 1), iy in 0:(ny - 1), iz in 0:(nz - 1)
        push!(values, ix == 0 ? NaN : -2.0 * (iz * step) + 4.0)
    end
    # CubeGrid constructor rejects non-finite values; build via direct field access
    # is disallowed, so a non-finite cube simply cannot be constructed.
    @test_throws ArgumentError CubeGrid(origin, axes, (nx, ny, nz), values)
end

@testset "constant_current_surface: z-axis must point toward vacuum" begin
    # axes[3,3] < 0: z decreases with iz, which violates the vacuum-at-high-index contract.
    bad_axes = Float64[0.1 0.0 0.0; 0.0 0.1 0.0; 0.0 0.0 -0.1]
    @test_throws ArgumentError constant_current_surface(
        CubeGrid([0.0, 0.0, 0.0], bad_axes, (2, 2, 3), zeros(12)), 1.0)
end

# --------------------------------------------------------------------------------------
# isovalue_for_mean_height
# --------------------------------------------------------------------------------------

@testset "isovalue_for_mean_height: planar cube exact inversion" begin
    cube = _planar_cube(; slope=-2.0, intercept=4.0)
    frame = _frame_at_origin()
    # Target 0.50 nm: analytic iso = 4 - 2 * 0.50 = 3.0.
    iso = isovalue_for_mean_height(cube, frame, 0.50)
    @test iso ≈ 3.0 atol = 1e-6
    # Verify round-trip: surface at this isovalue has mean height 0.50 nm.
    surf = constant_current_surface(cube, iso)
    valid_heights = surf.height[surf.valid]
    @test !isempty(valid_heights)
    @test mean(valid_heights) ≈ 0.50 atol = 1e-6
end

@testset "isovalue_for_mean_height: scan policy is explicit and validated" begin
    cube = _planar_cube(; slope=-2.0, intercept=4.0)
    frame = _frame_at_origin()
    implicit = isovalue_for_mean_height(cube, frame, 0.50)
    explicit = isovalue_for_mean_height(cube, frame, 0.50;
        scan_intervals=DEFAULT_ISOVALUE_SCAN_INTERVALS)
    @test DEFAULT_ISOVALUE_SCAN_INTERVALS == 1024
    @test implicit == explicit
    for invalid in (0, -1, 1.5, true)
        @test_throws ArgumentError isovalue_for_mean_height(cube, frame, 0.50;
            scan_intervals=invalid)
    end
end

@testset "isovalue_for_mean_height: iteration and tolerance parameters are validated" begin
    cube = _planar_cube(; slope=-2.0, intercept=4.0)
    frame = _frame_at_origin()
    for invalid in (0, -1, 1.5, true)
        @test_throws ArgumentError isovalue_for_mean_height(cube, frame, 0.50;
            max_iter=invalid)
    end
    for invalid in (0.0, -1.0, NaN, Inf, -Inf, true)
        @test_throws ArgumentError isovalue_for_mean_height(cube, frame, 0.50;
            height_tol_factor=invalid)
    end
end

@testset "isovalue_for_mean_height: exponential cube recovery" begin
    step = 0.02
    decay = 0.20
    amp = 10.0
    cube = _exponential_cube(; step=step, decay=decay, amplitude=amp, nz=100)
    frame = _frame_at_origin()
    target = 0.50
    iso = isovalue_for_mean_height(cube, frame, target)
    # Analytic: z_cross = -decay * ln(iso/amp). Set z_cross = target -> iso = amp * exp(-target/decay).
    analytic_iso = amp * exp(-target / decay)
    @test iso ≈ analytic_iso rtol = 1e-3
    surf = constant_current_surface(cube, iso)
    valid_heights = surf.height[surf.valid]
    @test mean(valid_heights) ≈ target atol = step
end

@testset "isovalue_for_mean_height: bump cube mean height policy" begin
    cube = _bump_cube()
    frame = _frame_at_origin()
    target = 0.50
    iso = isovalue_for_mean_height(cube, frame, target)
    @test isfinite(iso)
    @test iso > 0
    surf = constant_current_surface(cube, iso)
    valid_heights = surf.height[surf.valid]
    @test !isempty(valid_heights)
    # The calibrated isovalue must put the mean valid height at the policy target.
    @test mean(valid_heights) ≈ target atol = 5e-3
end

@testset "isovalue_for_mean_height: explicit isovalue range option" begin
    cube = _planar_cube()
    frame = _frame_at_origin()
    iso = isovalue_for_mean_height(cube, frame, 0.50; iso_lo=0.5, iso_hi=3.5)
    @test iso ≈ 3.0 atol = 1e-6
end

@testset "isovalue_for_mean_height: unique interior branch survives invalid endpoints" begin
    cube = _column_cube([[3.0, 2.0, 1.0]])
    frame = _frame_at_origin()
    iso1 = isovalue_for_mean_height(cube, frame, 1.1; iso_lo=0.0, iso_hi=4.0)
    iso2 = isovalue_for_mean_height(cube, frame, 1.1; iso_lo=0.0, iso_hi=4.0)
    @test iso1 ≈ 1.9 atol = 1e-6
    @test iso2 == iso1  # fixed 1025-point scan and bisection policy is deterministic
end

@testset "isovalue_for_mean_height: multiple target branches are explicitly ambiguous" begin
    z = collect(0.0:1.0:10.0)
    broad_decreasing = 10.0 .- z
    narrow_increasing = 2.0 .+ 0.2 .* z
    cube = _column_cube([broad_decreasing, narrow_increasing])
    msg = _error_message() do
        isovalue_for_mean_height(cube, _frame_at_origin(), 5.0;
                                 iso_lo=0.1, iso_hi=9.9)
    end
    @test occursin("ambiguous", lowercase(msg))
    @test occursin("multiple", lowercase(msg)) || occursin("nonmonotone", lowercase(msg))
end

@testset "isovalue_for_mean_height: support-changing target bracket is ambiguous" begin
    z = collect(0.0:1.0:10.0)
    broad_decreasing = 10.0 .- z
    entering_increasing = 2.0 .+ 0.1 .* z
    cube = _column_cube([broad_decreasing, entering_increasing])
    msg = _error_message() do
        isovalue_for_mean_height(cube, _frame_at_origin(), 7.0;
                                 iso_lo=1.5, iso_hi=2.5)
    end
    @test occursin("ambiguous discontinuous", lowercase(msg))
    @test occursin("support changes", lowercase(msg))
end

@testset "isovalue_for_mean_height: exact scan endpoints are continuity-ambiguous" begin
    frame = _frame_at_origin()

    # At the lower endpoint iso=2.0, only the first column is valid and its
    # crossing is exactly z=1.0. At the one-sided neighbour iso=2.5, both
    # columns are valid, so the endpoint root lacks two-sided fixed support.
    lower_cube = _column_cube([[4.0, 2.0, 0.0], [4.0, 3.0, 2.4]])
    lower_msg = _error_message() do
        isovalue_for_mean_height(lower_cube, frame, 1.0;
                                 iso_lo=2.0, iso_hi=2.5, scan_intervals=1)
    end
    @test occursin("endpoint", lowercase(lower_msg))
    @test occursin("continuity", lowercase(lower_msg))

    # At the upper endpoint iso=3.0, only the first column remains valid and
    # its crossing is exactly z=1.0. The lower one-sided neighbour has both
    # columns valid, which is the opposite interval endpoint orientation.
    upper_cube = _column_cube([[4.0, 3.0, 2.0], [2.8, 2.4, 2.0]])
    upper_msg = _error_message() do
        isovalue_for_mean_height(upper_cube, frame, 1.0;
                                 iso_lo=2.5, iso_hi=3.0, scan_intervals=1)
    end
    @test occursin("endpoint", lowercase(upper_msg))
    @test occursin("continuity", lowercase(upper_msg))
end

@testset "isovalue_for_mean_height: isolated exact interior response is continuity-ambiguous" begin
    frame = _frame_at_origin()

    # The exact response at iso=2.0 has no finite response at either immediate
    # neighbour: iso=1.0 has two crossings and iso=3.0 has none.
    isolated_cube = _column_cube([[0.0, 2.0, 0.0]])
    isolated_msg = _error_message() do
        isovalue_for_mean_height(isolated_cube, frame, 1.0;
                                 iso_lo=1.0, iso_hi=3.0, scan_intervals=2)
    end
    @test occursin("continuity", lowercase(isolated_msg))
    @test occursin("isolated", lowercase(isolated_msg))

    # The exact response at iso=2.0 has a finite response only on its upper
    # side; the immediate lower neighbour has no crossing.
    one_sided_cube = _column_cube([[1.5, 2.0, 4.0]])
    one_sided_msg = _error_message() do
        isovalue_for_mean_height(one_sided_cube, frame, 1.0;
                                 iso_lo=1.0, iso_hi=3.0, scan_intervals=2)
    end
    @test occursin("continuity", lowercase(one_sided_msg))
    @test occursin("isolated", lowercase(one_sided_msg))
end

@testset "isovalue_for_mean_height: no valid branch remains explicit" begin
    cube = _column_cube([[3.0, 2.0, 1.0]])
    msg = _error_message() do
        isovalue_for_mean_height(cube, _frame_at_origin(), 1.0;
                                 iso_lo=4.0, iso_hi=5.0)
    end
    @test occursin("no valid crossings", lowercase(msg))
end

@testset "isovalue_for_mean_height: unreachable target remains explicit" begin
    cube = _column_cube([[3.0, 2.0, 1.0]])
    msg = _error_message() do
        isovalue_for_mean_height(cube, _frame_at_origin(), 3.0;
                                 iso_lo=1.1, iso_hi=2.9)
    end
    @test occursin("outside the reachable", lowercase(msg))
end

@testset "isovalue_for_mean_height: unattainable target errors" begin
    # Target beyond cube z-extent: no isovalue can reach it.
    cube = _planar_cube()
    frame = _frame_at_origin()
    @test_throws ErrorException isovalue_for_mean_height(cube, frame, 100.0)
end

# --------------------------------------------------------------------------------------
# critic2 vacuum-side convention sanity
# --------------------------------------------------------------------------------------

@testset "vacuum-side first crossing convention" begin
    # Two-bracket column: confirm the bracket closest to vacuum (high index) is
    # reported as the first crossing, matching critic2's vacuum-side convention.
    z = collect(0.0:0.1:0.5)
    v = [5.0, 0.2, 3.0, 0.2, 3.0, 0.2]
    res = first_vacuum_crossing(z, v, 1.0)
    @test res.status === :ambiguous
    @test res.bracket == length(z)          # pair (z[6], z[5]) closest to vacuum
    @test res.n_brackets >= 2
end
