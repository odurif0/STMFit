#!/usr/bin/env julia

# Validate connected mold inputs before running scientific decoding.
# Checks the aligned site/proxy TSV, unary template TSV, bond template TSV,
# and patch/template pixel compatibility.

using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _read_tsv

const DEFAULT_ATOMS = "templates/chitosan_geometric_sites.tsv"
const DEFAULT_TEMPLATES = "templates/chitosan_connected_molds.tsv"
const DEFAULT_BONDS = "templates/chitosan_connected_bond_molds.tsv"
const DEFAULT_PATCHES = "results/unit_separability/lobe_patches_selectedN_primary.tsv"
const DEFAULT_REPORT = "results/unit_assignment/connected_mold_validation.txt"

struct Options
    atoms::String
    templates::String
    bonds::String
    patches::String
    prefix::String
    report::String
    require_all::Bool
end

function _parse_cli(args)
    atoms = DEFAULT_ATOMS
    templates = DEFAULT_TEMPLATES
    bonds = DEFAULT_BONDS
    patches = DEFAULT_PATCHES
    prefix = "res_p"
    report = DEFAULT_REPORT
    require_all = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--atoms"; atoms = args[i+1]; i += 2
        elseif startswith(arg, "--atoms="); atoms = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--templates"; templates = args[i+1]; i += 2
        elseif startswith(arg, "--templates="); templates = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--bond-templates"; bonds = args[i+1]; i += 2
        elseif startswith(arg, "--bond-templates="); bonds = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--patches"; patches = args[i+1]; i += 2
        elseif startswith(arg, "--patches="); patches = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--prefix"; prefix = args[i+1]; i += 2
        elseif startswith(arg, "--prefix="); prefix = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--report"; report = args[i+1]; i += 2
        elseif startswith(arg, "--report="); report = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--require-all"; require_all = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/validate_connected_molds.jl [options]

            Options:
              --atoms PATH           Geometric/proxy site TSV [$(DEFAULT_ATOMS)]
              --templates PATH       Unary connected mold TSV [$(DEFAULT_TEMPLATES)]
              --bond-templates PATH  Sliding bond mold TSV [$(DEFAULT_BONDS)]
              --patches PATH         Patch TSV for pixel-count compatibility [$(DEFAULT_PATCHES)]
              --prefix STR           Patch prefix raw_p or res_p [res_p]
              --report PATH          Validation report [$(DEFAULT_REPORT)]
              --require-all          Treat missing optional upstream files as errors

            The validator never reads truth labels and never checks composition.
            It only checks file format, required type/parity/mirror combinations,
            and template/patch dimensions.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    return Options(atoms, templates, bonds, patches, prefix, report, require_all)
end

function _has_columns(header, required)
    missing = [c for c in required if !(c in header)]
    return isempty(missing), missing
end

function _numeric_values(rows, cols)
    vals = Float64[]
    for row in rows, col in cols
        v = tryparse(Float64, get(row, col, ""))
        v === nothing || push!(vals, v)
    end
    return vals
end

function _push_file_check!(lines, errors, label, path; required::Bool=true)
    if isempty(strip(path)) && !required
        push!(lines, "[SKIP] $label: not supplied")
        return false
    end
    if isfile(path)
        push!(lines, "[OK] $label: $path")
        return true
    end
    msg = "[MISSING] $label: $path"
    push!(lines, msg)
    required && push!(errors, msg)
    return false
end

function _validate_atoms!(lines, errors, opt)
    isfile(opt.atoms) || return
    header, rows = _read_tsv(opt.atoms)
    required = ["type", "atom", "t_nm", "u_nm", "weight", "sigma_t_nm", "sigma_u_nm"]
    ok, missing = _has_columns(header, required)
    ok || (push!(errors, "Atom TSV missing columns: $(join(missing, ", "))"); return)
    counts = Dict(0 => 0, 1 => 0)
    for row in rows
        typ = tryparse(Int, row["type"])
        typ in (0, 1) && (counts[typ] += 1)
    end
    for typ in (0, 1)
        counts[typ] > 0 || push!(errors, "Atom TSV has no rows for type=$typ")
    end
    vals = _numeric_values(rows, ["weight", "sigma_t_nm", "sigma_u_nm"])
    any(v -> !isfinite(v) || v <= 0, vals) && push!(errors, "Atom TSV contains non-positive/non-finite weights or sigmas")
    push!(lines, "  atom counts: type0=$(counts[0]) type1=$(counts[1])")
end

function _validate_unary_templates!(lines, errors, opt)
    isfile(opt.templates) || return nothing
    header, rows = _read_tsv(opt.templates)
    required = ["name", "type", "parity", "mirror"]
    ok, missing = _has_columns(header, required)
    ok || (push!(errors, "Unary template TSV missing columns: $(join(missing, ", "))"); return nothing)
    pix = [c for c in header if occursin(r"^p\d+$", c)]
    isempty(pix) && push!(errors, "Unary template TSV has no pNNN pixel columns")
    keys = Set{Tuple{Int,Int,Int}}()
    for row in rows
        typ = tryparse(Int, row["type"])
        parity = tryparse(Int, row["parity"])
        mirror = tryparse(Int, row["mirror"])
        typ in (0, 1) && parity in (0, 1) && mirror in (0, 1) && push!(keys, (typ, parity, mirror))
    end
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        (typ, parity, mirror) in keys || push!(errors, "Unary template missing type=$typ parity=$parity mirror=$mirror")
    end
    push!(lines, "  unary templates: rows=$(length(rows)) pixels=$(length(pix))")
    return pix
end

function _validate_bond_templates!(lines, errors, opt)
    isfile(opt.bonds) || return nothing
    header, rows = _read_tsv(opt.bonds)
    required = ["name", "left_type", "right_type", "parity", "mirror"]
    ok, missing = _has_columns(header, required)
    ok || (push!(errors, "Bond template TSV missing columns: $(join(missing, ", "))"); return nothing)
    lpix = [c for c in header if occursin(r"^l_p\d+$", c)]
    rpix = [c for c in header if occursin(r"^r_p\d+$", c)]
    isempty(lpix) && push!(errors, "Bond template TSV has no l_pNNN pixel columns")
    length(lpix) == length(rpix) || push!(errors, "Bond template left/right pixel counts differ")
    keys = Set{Tuple{Int,Int,Int,Int}}()
    for row in rows
        lt = tryparse(Int, row["left_type"])
        rt = tryparse(Int, row["right_type"])
        parity = tryparse(Int, row["parity"])
        mirror = tryparse(Int, row["mirror"])
        lt in (0, 1) && rt in (0, 1) && parity in (0, 1) && mirror in (0, 1) && push!(keys, (lt, rt, parity, mirror))
    end
    for lt in (0, 1), rt in (0, 1), parity in (0, 1), mirror in (0, 1)
        (lt, rt, parity, mirror) in keys || push!(errors, "Bond template missing left=$lt right=$rt parity=$parity mirror=$mirror")
    end
    push!(lines, "  bond templates: rows=$(length(rows)) left_pixels=$(length(lpix)) right_pixels=$(length(rpix))")
    return lpix
end

function _validate_patch_compat!(lines, errors, opt, unary_pix, bond_pix)
    isfile(opt.patches) || return
    header, rows = _read_tsv(opt.patches)
    pix = [c for c in header if startswith(c, opt.prefix)]
    isempty(pix) && (push!(errors, "Patch TSV has no columns with prefix $(opt.prefix)"); return)
    unary_pix !== nothing && length(unary_pix) != length(pix) && push!(errors, "Unary pixel count $(length(unary_pix)) != patch pixel count $(length(pix))")
    bond_pix !== nothing && length(bond_pix) != length(pix) && push!(errors, "Bond pixel count $(length(bond_pix)) != patch pixel count $(length(pix))")
    files = Set(get(row, "file", "") for row in rows)
    push!(lines, "  patches: rows=$(length(rows)) files=$(length(files)) $(opt.prefix)_pixels=$(length(pix))")
end

function main()
    opt = _parse_cli(ARGS)
    lines = String[]
    errors = String[]

    push!(lines, "Connected Mold Validation")
    push!(lines, "=========================")
    push!(lines, "")

    _push_file_check!(lines, errors, "geometric/proxy site TSV", opt.atoms; required=opt.require_all)
    _push_file_check!(lines, errors, "unary template TSV", opt.templates; required=true)
    _push_file_check!(lines, errors, "bond template TSV", opt.bonds; required=false)
    _push_file_check!(lines, errors, "patch TSV", opt.patches; required=false)
    push!(lines, "")

    _validate_atoms!(lines, errors, opt)
    unary_pix = _validate_unary_templates!(lines, errors, opt)
    bond_pix = _validate_bond_templates!(lines, errors, opt)
    _validate_patch_compat!(lines, errors, opt, unary_pix, bond_pix)

    push!(lines, "")
    if isempty(errors)
        push!(lines, "RESULT: OK")
    else
        push!(lines, "RESULT: FAIL ($(length(errors)) issue(s))")
        append!(lines, ["  - $e" for e in errors])
    end

    mkpath(dirname(opt.report))
    open(opt.report, "w") do io
        for line in lines
            println(io, line)
        end
    end
    println(join(lines, '\n'))
    println("\nWrote: ", opt.report)
    isempty(errors) || exit(1)
end

main()
