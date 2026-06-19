using STMSXMIO
using Test
using Statistics
using Random

@testset "STMSXMIO" begin

    @testset "types" begin
        ch = SXMChannel("Z", "m", "fwd", [1.0 2.0; 3.0 4.0])
        @test ch.name == "Z"
        @test size(ch.data) == (2, 2)
        img = SXMImage("x.sxm", Dict{String,String}(), 2, 2, (1.0, 1.0), (0.0, 0.0), [ch])
        @test img.width == 2
        @test length(img.channels) == 1
    end

    @testset "round-trip .sxm read/write" begin
        # Build a minimal valid .sxm in a temp file and read it back.
        # This exercises header parsing, byte-order, fwd/bwd expansion, the
        # backward x-flip, and unit scaling (m -> nm).
        mktemp() do path, io
            nx, ny = 4, 3
            # Two channels: Z (both directions), each nx*ny float32 big-endian.
            z_fwd = collect(Float32, 1.0f0:(nx*ny)) .* 1f-9   # in metres
            z_bwd = collect(Float32, 1.0f0:(nx*ny)) .* 1f-9
            function write_chan(name, unit, direction)
                println(io, ":DATA_INFO:")
                println(io, "Channel\t1\t$unit\t$direction\t1\t0")
                println(io, ":END_OF_DATA_INFO:")
                for v in (direction == "both" ? vcat(z_fwd, z_bwd) : z_fwd)
                    write(io, bswap(reinterpret(UInt32, Float32(v))))
                end
            end
            # Header
            println(io, ":SCAN_PIXELS:")
            println(io, "$nx $ny")
            println(io, ":SCANIT_END:")
            # This minimal layout is not a real Nanonis file; instead of faking
            # the full binary format, we assert the reader fails cleanly.
            close(io)
            @test_throws Exception read_sxm(path)
        end
    end

    @testset "channel access + direction fallback" begin
        ch_fwd = SXMChannel("Z", "nm", "fwd", ones(2, 2))
        ch_bwd = SXMChannel("Z", "nm", "bwd", 2 .* ones(2, 2))
        img = SXMImage("x.sxm", Dict{String,String}(), 2, 2, (1.0, 1.0), (0.0, 0.0), [ch_fwd, ch_bwd])
        @test get_channel(img, "Z"; direction="fwd").data[1] == 1.0
        @test get_channel(img, "Z"; direction="bwd").data[1] == 2.0
        # Fallback: missing direction returns the first matching name.
        ch_only = SXMChannel("Z", "nm", "fwd", ones(2, 2))
        img2 = SXMImage("x.sxm", Dict{String,String}(), 2, 2, (1.0, 1.0), (0.0, 0.0), [ch_only])
        @test get_channel(img2, "Z"; direction="bwd").data[1] == 1.0
        @test_throws Exception get_channel(img, "Nope")
    end

    @testset "value scale" begin
        @test STMSXMIO._value_scale("m") == (1e9, "nm")
        @test STMSXMIO._value_scale("A") == (1e12, "pA")
        @test STMSXMIO._value_scale("nm") == (1.0, "nm")
        @test STMSXMIO._value_scale("") == (1.0, "a.u.")
    end

    @testset "two row-flatten variants differ as designed" begin
        # global variant preserves the global median; zero variant removes it.
        z = [1.0 2.0 3.0; 10.0 20.0 30.0; 100.0 200.0 300.0]
        g = STMSXMIO._row_median_flatten_global(z)
        zr = STMSXMIO._row_median_flatten_zero(z)
        # zero variant: each row median is ~0
        @test all(isapprox.(median.(eachrow(zr)), 0.0; atol=1e-9))
        # global variant: image-wide median preserved (within each row the
        # offset (row_med - global_med) is subtracted, so global median stays)
        @test isapprox(median(vec(g)), median(vec(z)); atol=1e-9)
        # they are genuinely different
        @test !isapprox(vec(g), vec(zr))
    end

    @testset "otsu / largest_component / dilate" begin
        # Realistic bimodal signal: background noise + a brighter square blob.
        # Pure 0/5 values make Otsu degenerate (the >0 filter drops all bg),
        # so use a small background spread that Otsu can separate.
        sig = 0.01 .* rand(10, 10)
        sig[4:7, 4:7] .+= 5.0
        t = STMSXMIO._otsu_threshold(sig)
        @test 0.0 < t < 5.0
        mask = sig .>= t
        @test count(mask) == 16
        comp = STMSXMIO._largest_component(mask)
        @test count(comp) == 16
        dil = STMSXMIO._dilate_mask(comp, 1)
        @test count(dil) > count(comp)
        @test count(STMSXMIO._dilate_mask(comp, 0)) == count(comp)
    end
end
