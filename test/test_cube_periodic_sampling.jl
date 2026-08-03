#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "cube_to_stm_maps.jl"))

function periodic_test_cube(; n=(7, 6, 5), skew=0.25)
    axes = [0.10 skew * 0.08 0.0;
            0.00 0.08        0.0;
            0.00 0.00        0.12]
    values = Float64[]
    for ix in 0:(n[1] - 1), iy in 0:(n[2] - 1), iz in 0:(n[3] - 1)
        push!(values, sinpi(2ix / n[1]) + 0.4cospi(2iy / n[2]) + 0.1iz)
    end
    return CubeGrid([0.3, -0.2, 0.5], axes, n, values)
end

@testset "diagnostic periodic cube sampling" begin
    cube = periodic_test_cube()
    q = [2.35, 3.20, 1.40]
    point = cube.origin + cube.axes * q

    @test _sample_cube_periodic(cube, point, (true, true, false)) ≈
          _sample_cube(cube, point) atol=1e-14

    period_x = cube.n[1] .* cube.axes[:, 1]
    period_y = cube.n[2] .* cube.axes[:, 2]
    reference = _sample_cube_periodic(cube, point, (true, true, false))
    @test _sample_cube_periodic(cube, point + period_x, (true, true, false)) ≈ reference atol=1e-13
    @test _sample_cube_periodic(cube, point - 2period_y, (true, true, false)) ≈ reference atol=1e-13

    epsilon = 1e-8
    left = cube.origin + cube.axes * [-epsilon, q[2], q[3]]
    right = cube.origin + cube.axes * [cube.n[1] - epsilon, q[2], q[3]]
    @test _sample_cube_periodic(cube, left, (true, true, false)) ≈
          _sample_cube_periodic(cube, right, (true, true, false)) atol=1e-13

    outside_z = cube.origin + cube.axes * [q[1], q[2], cube.n[3] + 0.1]
    @test isnan(_sample_cube_periodic(cube, outside_z, (true, true, false)))
    @test_throws ArgumentError _sample_cube_periodic(cube, point, (true, true, true))

    @test _parse_periodic_axes("none") == (false, false, false)
    @test _parse_periodic_axes("x") == (true, false, false)
    @test _parse_periodic_axes("xy") == (true, true, false)
    @test _parse_periodic_axes("y,x") == (true, true, false)
    @test_throws ArgumentError _parse_periodic_axes("z")
    @test_throws ArgumentError _parse_periodic_axes("xz")
    @test_throws ArgumentError _parse_periodic_axes("banana")
end

@testset "observable and constant-current CLI parsing" begin
    # --observable is opt-in and must accept only the documented values.
    # Use a real temp file because _parse_cli validates cube existence.
    tmp = tempname()
    write(tmp, "dummy cube\nv\n0 0 0 0\n1 0.1 0 0\n1 0 0.1 0\n1 0 0 0.1\n0.0\n")
    common = ["--cube", "0:$tmp", "--origin", "0,0,0",
              "--t-axis", "1,0,0", "--u-axis", "0,1,0", "--height-nm", "0.5"]
    @test _parse_cli(copy(common)).observable == "constant-height"
    cc = copy(common)
    append!(cc, ["--observable", "constant-current"])
    @test _parse_cli(cc).observable == "constant-current"
    bad = copy(common)
    append!(bad, ["--observable", "bogus"])
    @test_throws Exception _parse_cli(bad)
    iso = copy(common)
    append!(iso, ["--observable", "constant-current", "--isovalue", "1.5"])
    opt_iso = _parse_cli(iso)
    @test opt_iso.isovalue == 1.5
    @test opt_iso.mean_height_nm == 0.50
    mh = copy(common)
    append!(mh, ["--observable", "constant-current", "--mean-height-nm", "0.45"])
    opt_mh = _parse_cli(mh)
    @test opt_mh.isovalue === nothing
    @test opt_mh.mean_height_nm == 0.45
    @test opt_mh.isovalue_scan_intervals == 1024
    scan = copy(common)
    append!(scan, ["--observable", "constant-current",
                   "--isovalue-scan-intervals", "64"])
    @test _parse_cli(scan).isovalue_scan_intervals == 64
    for invalid in ("0", "-1", "1.5")
        bad_scan = copy(common)
        append!(bad_scan, ["--observable", "constant-current",
                           "--isovalue-scan-intervals", invalid])
        @test_throws Exception _parse_cli(bad_scan)
    end
end

@testset "constant-current paired output transaction is symlink-safe and recoverable" begin
    mktempdir() do dir
        map = joinpath(dir, "map.tsv")
        mask = joinpath(dir, "map.mask.tsv")
        victim = joinpath(dir, "victim.tsv")
        write(map, "old map\n")
        write(victim, "protected victim\n")
        symlink(victim, mask)

        @test_throws ErrorException _with_output_transaction(
                [mask, map]; gate=map, failpoint=:after_install_1) do staged
            write(staged[map], "new map\n")
            write(staged[mask], "new mask\n")
        end
        @test read(map, String) == "old map\n"
        @test islink(mask)
        @test readlink(mask) == victim
        @test read(victim, String) == "protected victim\n"

        @test_throws CubeTransactionInterruption _with_output_transaction(
                [mask, map]; gate=map, interruptpoint=:after_gate) do staged
            write(staged[map], "interrupted map\n")
            write(staged[mask], "interrupted mask\n")
        end
        recovered = Ref(false)
        _with_output_transaction([mask, map]; gate=map) do staged
            recovered[] = read(map, String) == "old map\n" && islink(mask)
            write(staged[map], "final map\n")
            write(staged[mask], "final mask\n")
        end
        @test recovered[]
        @test read(map, String) == "final map\n"
        @test read(mask, String) == "final mask\n"
        @test read(victim, String) == "protected victim\n"
        @test isempty(filter(name -> occursin("stmfit-txn", name), readdir(dir)))
    end
end

@testset "eleven-output provenance-gated transaction covers every phase" begin
    mktempdir() do dir
        artifacts = [joinpath(dir, "artifact-$(lpad(index, 2, '0')).tsv")
            for index in 1:10]
        provenance = joinpath(dir, "artifacts.provenance.toml")
        destinations = vcat(artifacts, [provenance])
        old_bytes = Dict(path => "old $(basename(path))\n" for path in destinations)
        new_bytes = Dict(path => "new $(basename(path))\n" for path in destinations)
        restore!(bytes) = foreach(path -> write(path, bytes[path]), destinations)
        complete(bytes) = all(path -> isfile(path) && !islink(path) &&
            read(path, String) == bytes[path], destinations)
        residue() = filter(name -> occursin("stmfit-txn", name), readdir(dir))
        backup_phases = [Symbol("after_backup_$index") for index in 1:11]
        install_phases = [Symbol("after_install_$index") for index in 1:10]
        phases = vcat([:after_marker], backup_phases, install_phases,
            [:after_gate, :after_commit])
        restore!(old_bytes)

        for point in phases
            interruption = try
                _with_output_transaction(destinations; gate=provenance,
                        interruptpoint=point) do staged
                    for path in destinations
                        write(staged[abspath(path)], new_bytes[path])
                    end
                end
                nothing
            catch caught
                caught
            end
            @test interruption isa CubeTransactionInterruption
            if interruption isa CubeTransactionInterruption
                @test interruption.point == point
            end
            @test isfile(_cc_txn_marker(abspath(provenance)))
            if point == :after_marker
                @test complete(old_bytes)
            elseif point in backup_phases || point in install_phases
                @test !_cc_txn_present(provenance)
            else
                @test complete(new_bytes)
            end

            _cc_txn_recover(abspath.(destinations), abspath(provenance))
            @test complete(point == :after_commit ? new_bytes : old_bytes)
            @test isempty(residue())
            restore!(old_bytes)
        end

        writer_called = Ref(false)
        marker = _cc_txn_marker(abspath(provenance))
        write(marker, "malformed\n")
        @test_throws ErrorException _with_output_transaction(
                destinations; gate=provenance) do _
            writer_called[] = true
        end
        @test !writer_called[]
        @test complete(old_bytes)
        rm(marker)

        marker_victim = joinpath(dir, "marker-victim")
        write(marker_victim, "protected\n")
        symlink(marker_victim, marker)
        @test_throws ErrorException _with_output_transaction(
                destinations; gate=provenance) do _
            writer_called[] = true
        end
        @test !writer_called[]
        @test complete(old_bytes)
        @test read(marker_victim, String) == "protected\n"
        rm(marker)
        @test isempty(residue())
    end
end
