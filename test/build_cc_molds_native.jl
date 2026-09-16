#!/usr/bin/env julia
# Native reconstruction of cc_mold_builder.py --legacy. No Python subprocesses.
include(joinpath(@__DIR__, "lib", "cc_mold_native.jl"))
using .CCMoldNative
using Printf

function cc_mold_main(args=ARGS)
    if args == ["--help"] || args == ["-h"]
        println("""
        Usage: julia --project=. test/build_cc_molds_native.jl \\
          --cube0 PATH --cube1 PATH --frame0 PATH --frame1 PATH \\
          --config config/unit_assignment_reconstructed.toml --out PATH

        All six options are required. Reads scalar LDOS cubes (bohr) and frame
        TSVs (origin_nm, t_axis, u_axis). Reads [model] mold_* settings only.
        Writes eight scorer-compatible rows (GlcN/GlcNAc × parity × mirror).
        Existing output is never replaced. No benchmark labels are read.

        This preserves the Python --legacy first-axis-fast cube indexing,
        highest-occupied-height rule, and first mean height below target.
        It is not the strict first-vacuum-crossing diagnostic builder.
        """)
        return 0
    end
    options = Dict{String,String}()
    required = ("--cube0", "--cube1", "--frame0", "--frame1", "--config", "--out")
    i = 1
    while i <= length(args)
        pair = split(args[i], '='; limit=2)
        key = pair[1]
        key in required || throw(ArgumentError("unknown option: $key"))
        haskey(options, key) && throw(ArgumentError("duplicate option: $key"))
        if length(pair) == 2
            value = pair[2]
            i += 1
        else
            i < length(args) && !startswith(args[i + 1], "--") ||
                throw(ArgumentError("missing value for $key"))
            value = args[i + 1]
            i += 2
        end
        isempty(value) && throw(ArgumentError("empty value for $key"))
        options[key] = value
    end
    missing = filter(k -> !haskey(options, k), required)
    isempty(missing) || throw(ArgumentError("missing required options: $(join(missing, ", "))"))
    result = build_molds(cube0=options["--cube0"], cube1=options["--cube1"],
                         frame0=options["--frame0"], frame1=options["--frame1"],
                         config=options["--config"], out=options["--out"])
    for (label, item) in (("GlcN", result.type0), ("GlcNAc", result.type1))
        @printf("%s: iso=%.5g mean_h=%.4f valid=%d/%d\n",
                label, item.iso, item.mean_height, item.nvalid, length(item.heights))
    end
    println("templates: ", options["--out"])
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(cc_mold_main())
    catch err
        err isa InterruptException && rethrow()
        println(stderr, "cc mold builder: ", sprint(showerror, err))
        exit(1)
    end
end
