#!/usr/bin/env julia

# Checker for config/unit_assignment_candidate.toml.
#
# Enforces the T0 contract from `.omo/plans/improve-unit-assignment-benchmark.md`:
#   1. FIREWALL — forbidden benchmark motif/encoding, truth/count column names,
#      and benchmark truth/grade paths may appear only inside [grader_only].
#   2. STRUCTURE — every T0 decision is present with the exact frozen value.
#   3. STATE — grade_status in {locked, frozen_once, graded}; provenance.status
#      and candidate.frozen_hash are consistent with grade_status.
#   4. HASH BINDING — when grade_status = "frozen_once" or "graded",
#      candidate.frozen_hash must equal SHA-256 of the exact non-grader source
#      bytes, excluding only [grader_only] plus the frozen_hash and grade_status
#      assignment text needed for self-reference and the frozen_once -> graded
#      lifecycle transition. Other comments and line endings remain bound.
#
# CLI:
#   julia --project=. test/check_unit_assignment_candidate_manifest.jl \
#       --config config/unit_assignment_candidate.toml [--expect-locked | --expect-frozen-once]
#
# Exit 0 on PASS with a parsed-status summary on stdout; exit nonzero on FAIL
# with diagnostics on stderr. PASS is decided from parsed values and exit code,
# not from log wording.

using SHA, TOML

const FORBIDDEN_TOKENS = ["NKNNKN", "010010", "sequence", "expected_N", "target_N"]
const FORBIDDEN_PATH_FRAGMENTS = [
    "benchmarks/",
    "results/benchmark_grades",
    "report_unit_assignment_benchmark",
    "grade_unit_assignment",
    "_unit_sequences.tsv",
]

struct Options
    config::String
    expect_locked::Bool
    expect_frozen_once::Bool
end

function _fail(msg)
    println(stderr, "Unit-assignment candidate manifest: FAIL")
    println(stderr, msg)
    return 1
end

function _parse_cli(args)
    config = ""
    expect_locked = false
    expect_frozen_once = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i < length(args) || error("--config requires a value")
            config = args[i+1]; i += 2
        elseif startswith(arg, "--config=")
            config = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--expect-locked"
            expect_locked = true; i += 1
        elseif arg == "--expect-frozen-once"
            expect_frozen_once = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/check_unit_assignment_candidate_manifest.jl --config PATH [--expect-locked | --expect-frozen-once]

            Validates the T0 firewall, structure, state, and hash-binding contract
            of config/unit_assignment_candidate.toml. PASS is decided from parsed
            values and exit code, not from log wording. Mutations after a frozen
            hash are rejected via the exact-source SHA-256 comparison.
            """)
            return nothing
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(config) && error("--config is required")
    return Options(config, expect_locked, expect_frozen_once)
end

# Walk a dotted path through nested Dicts; returns (value, ok).
function _get(parsed, path::AbstractString)
    parts = split(path, '.')
    cur = parsed
    for p in parts
        cur isa AbstractDict && haskey(cur, p) || return nothing, false
        cur = cur[p]
    end
    return cur, true
end

# Split a single-line TOML source line at its first comment marker outside a
# basic or literal string. Multiline strings are rejected separately, making
# table and assignment boundaries unambiguous without maintaining a TOML lexer.
function _split_comment(line::AbstractString)
    code = IOBuffer()
    comment = IOBuffer()
    in_basic = false
    in_literal = false
    escaped = false
    found_comment = false
    for c in line
        if found_comment
            write(comment, c)
        elseif in_basic
            write(code, c)
            if escaped
                escaped = false
            elseif c == '\\'
                escaped = true
            elseif c == '"'
                in_basic = false
            end
        elseif in_literal
            write(code, c)
            c == '\'' && (in_literal = false)
        elseif c == '#'
            found_comment = true
            write(comment, c)
        else
            write(code, c)
            c == '"' && (in_basic = true)
            c == '\'' && (in_literal = true)
        end
    end
    return String(take!(code)), String(take!(comment))
end

function _line_ending(line::String)
    endswith(line, "\r\n") && return "\r\n"
    endswith(line, "\n") && return "\n"
    endswith(line, "\r") && return "\r"
    return ""
end

function _assignment_key(code::String)
    m = match(r"^\s*(?:([A-Za-z0-9_-]+)|\"([^\"]+)\"|'([^']+)')\s*=", code)
    m === nothing && return nothing
    return something(m.captures[1], m.captures[2], m.captures[3])
end

# A deliberately conservative source contract. TOML multiline strings and
# non-simple table headers are rejected so text that resembles a table or
# assignment can never alter the firewall/hash boundary.
function _source_contract(text::String)
    section = ""
    grader_tables = 0
    frozen_hash_assignments = 0
    violations = String[]
    hash_buf = IOBuffer()

    for (lineno, raw) in enumerate(eachline(IOBuffer(text); keep=true))
        ending = _line_ending(raw)
        body = isempty(ending) ? raw : chop(raw; tail=length(ending))
        code, comment = _split_comment(body)

        if occursin("\"\"\"", code) || occursin("'''", code)
            push!(violations, "line $lineno: TOML multiline strings are not permitted in this manifest")
        end

        stripped = strip(code)
        header = match(r"^\[\s*([A-Za-z0-9_.-]+)\s*\]$", stripped)
        if startswith(stripped, '[') && endswith(stripped, ']') && header === nothing
            push!(violations, "line $lineno: unsupported table-header syntax; use simple [name] headers")
            section = "__invalid_header__"
        elseif header !== nothing
            section = header.captures[1]
            if section == "grader_only"
                grader_tables += 1
            elseif startswith(section, "grader_only.")
                push!(violations, "line $lineno: nested [grader_only.*] tables are forbidden")
            end
        end

        in_grader = section == "grader_only"
        if !in_grader
            lowered = lowercase(string(code, comment))
            for tok in FORBIDDEN_TOKENS
                occursin(lowercase(tok), lowered) && push!(violations, "line $lineno: forbidden token '$tok' outside [grader_only]: $(strip(body))")
            end
            for frag in FORBIDDEN_PATH_FRAGMENTS
                occursin(lowercase(frag), lowered) && push!(violations, "line $lineno: forbidden path fragment '$frag' outside [grader_only]: $(strip(body))")
            end
        end

        key = _assignment_key(code)
        if section == "candidate" && key == "frozen_hash"
            frozen_hash_assignments += 1
        end

        in_grader && continue
        if section == "candidate" && key in ("frozen_hash", "grade_status")
            # Exclude only lifecycle assignment text. Retain an inline comment
            # and the exact line ending so all other source bytes stay bound.
            write(hash_buf, comment)
            write(hash_buf, ending)
        else
            write(hash_buf, raw)
        end
    end

    grader_tables == 1 || push!(violations, "manifest must contain exactly one real [grader_only] table; found $grader_tables")
    return violations, frozen_hash_assignments, bytes2hex(sha256(take!(hash_buf)))
end

function _semantic_string_violations(value::AbstractString, path::String, kind::String)
    lowered = lowercase(value)
    violations = String[]
    for tok in FORBIDDEN_TOKENS
        occursin(lowercase(tok), lowered) &&
            push!(violations, "parsed $kind at $path contains forbidden token '$tok'")
    end
    for frag in FORBIDDEN_PATH_FRAGMENTS
        occursin(lowercase(frag), lowered) &&
            push!(violations, "parsed $kind at $path contains forbidden path fragment '$frag'")
    end
    return violations
end

# TOML.parse has already decoded \uXXXX and \UXXXXXXXX escapes. Walk every key
# and string value outside the sole top-level [grader_only] table so semantic
# tokens cannot evade the exact-source firewall by changing their spelling.
function _semantic_firewall(parsed::AbstractDict)
    violations = String[]

    function walk(value, path::String; root::Bool=false)
        if value isa AbstractDict
            for (key, child) in pairs(value)
                key_text = string(key)
                root && key_text == "grader_only" && continue
                key_path = isempty(path) ? key_text : "$path.$key_text"
                append!(violations, _semantic_string_violations(key_text, key_path, "key"))
                walk(child, key_path)
            end
        elseif value isa AbstractVector
            for (index, child) in pairs(value)
                walk(child, "$path[$index]")
            end
        elseif value isa AbstractString
            append!(violations, _semantic_string_violations(value, path, "string value"))
        end
    end

    walk(parsed, ""; root=true)
    return violations
end

const EXACT_CHECKS = [
    ("candidate.name", "hierarchical_equalprior_constant_current"),
    ("candidate.plan", ".omo/plans/improve-unit-assignment-benchmark.md"),
    ("candidate.version", 1),
    ("denominator.files", 145),
    ("denominator.control_positions_per_file", 6),
    ("denominator.control_positions_total", 870),
    ("features.policy", "equal_priors"),
    ("views.policy", "equal_weights"),
    ("bootstrap.count", 500),
    ("bootstrap.seeds_range_start", 0),
    ("bootstrap.seeds_range_stop", 499),
    ("bootstrap.resample_unit", "whole_scan_within_training_dates"),
    ("date_parser.rule", "leading_yyyymmdd_token"),
    ("date_parser.ambiguous", "fail"),
    ("date_parser.missing", "fail"),
    ("constant_current.mean_height_policy_nm", 0.50),
    ("constant_current.sensitivity_bracket_nm", [0.40, 0.45, 0.50, 0.55, 0.60]),
    ("constant_current.bracket_policy", "diagnostics_only"),
    ("real_gates.registration_interior_to_bounds", true),
    ("real_gates.beat_shifted_and_rotated_controls", true),
    ("real_gates.forward_backward_agreement", true),
    ("real_gates.fixed_perturbation_agreement", true),
    ("real_gates.abstain_on_invalid_isosurface", true),
    ("real_gates.abstain_on_registration_boundary", true),
    ("real_gates.do_not_widen_bounds_or_select_bracket_from_outputs", true),
    ("promotion.denominator_files", 145),
    ("promotion.denominator_control_positions", 870),
    ("promotion.honest_correct_min", 677),
    ("promotion.physical_accuracy_classified_min_fraction", 0.789),
    ("promotion.exact_chains_min", 18),
]

const CONTAIN_CHECKS = [
    ("features.hierarchical_emission_policy", "two-component"),
    ("features.hierarchical_emission_policy", "0.5"),
    ("features.hierarchical_emission_policy", "no occupancy regularizer"),
    ("one_vs_two_gate.pass_condition", "positive"),
    ("one_vs_two_gate.pass_condition", "95%"),
]

const REQUIRED_KEYS = String[
    "candidate", "denominator", "firewall", "provenance", "features", "views",
    "bootstrap", "date_parser", "constant_current", "one_vs_two_gate",
    "real_gates", "ranking", "promotion", "grader_only",
]
const MANDATORY_PROVENANCE_PATHS = String[
    "provenance.feature_tsves",
    "provenance.constant_current_cubes",
    "provenance.generated_maps",
    "provenance.molds",
    "provenance.configs",
]

_is_sha256_hex(s) = typeof(s) == String && length(s) == 64 && all(c -> c in "0123456789abcdef", s)

# Exact source for the frozen hash. Source-contract validation has already made
# section boundaries unambiguous. Only grader-only bytes and lifecycle
# assignment text are excluded; comments and line endings remain bound.
function _exact_source_hash(text::String)
    _, _, digest = _source_contract(text)
    return digest
end

function _require_ranking(parsed)
    order, ok = _get(parsed, "ranking.order")
    ok && order isa AbstractVector && !isempty(order) ||
        return ["ranking.order must be a non-empty list"]
    return String[]
end

function _check_view_weights(parsed)
    # Equal weights: views.policy must be equal_weights AND views.weight is a
    # finite number in (0,1]. The actual per-view weight is pinned at T0; T4 may
    # raise the count only by adding views that share this same weight.
    errs = String[]
    pol, ok = _get(parsed, "views.policy")
    (ok && pol == "equal_weights") || push!(errs, "views.policy must be \"equal_weights\", got $(ok ? pol : "<missing>")")
    w, ok2 = _get(parsed, "views.weight")
    (ok2 && w isa Number && isfinite(w) && w > 0.0 && w <= 1.0) || push!(errs, "views.weight must be a finite number in (0,1], got $(ok2 ? w : "<missing>")")
    return errs
end

function _check_locked(parsed, status, frozen_hash_assignments::Int)
    errs = String[]
    prov, ok = _get(parsed, "provenance.status")
    (ok && prov == "pending") || push!(errs, "locked candidate must have provenance.status = \"pending\", got $(ok ? prov : "<missing>")")
    fh, has_fh = _get(parsed, "candidate.frozen_hash")
    !has_fh || push!(errs, "locked candidate must not yet declare candidate.frozen_hash")
    frozen_hash_assignments == 0 || push!(errs, "locked candidate must contain no real candidate.frozen_hash assignment")
    return errs
end

function _check_frozen_provenance(parsed, status::String)
    errs = String[]
    for path in MANDATORY_PROVENANCE_PATHS
        value, ok = _get(parsed, path)
        if !ok
            push!(errs, "$status candidate is missing $path; every mandatory provenance field must be a lowercase 64-hex SHA-256")
        elseif !_is_sha256_hex(value)
            push!(errs, "$path must be a lowercase 64-hex SHA-256 for a $status candidate, got $(repr(value))")
        end
    end
    return errs
end

function _check_frozen_binding(parsed, text::String, status::String, frozen_hash_assignments::Int)
    errs = String[]
    prov, ok = _get(parsed, "provenance.status")
    (ok && prov == "frozen") || push!(errs, "$status candidate must have provenance.status = \"frozen\", got $(ok ? prov : "<missing>")")
    append!(errs, _check_frozen_provenance(parsed, status))
    fh, has_fh = _get(parsed, "candidate.frozen_hash")
    has_fh || (push!(errs, "$status candidate must declare candidate.frozen_hash"); return errs)
    frozen_hash_assignments == 1 || push!(errs, "$status candidate must contain exactly one real candidate.frozen_hash assignment; found $frozen_hash_assignments")
    _is_sha256_hex(fh) || push!(errs, "candidate.frozen_hash must be a 64-character lowercase-hex SHA-256, got: $fh")
    recomputed = _exact_source_hash(text)
    fh == recomputed || push!(errs, "Frozen hash mismatch: stored=$fh recomputed=$recomputed. Mutation after the freeze is rejected.")
    return errs
end

function _run(opt::Options)
    isfile(opt.config) || return _fail("Config not found: $(opt.config)")
    text = read(opt.config, String)
    parsed = try
        TOML.parse(text)
    catch e
        return _fail("TOML parse error in $(opt.config): $(sprint(showerror, e))")
    end

    # 1. Unambiguous source boundary and firewall
    violations, frozen_hash_assignments, _ = _source_contract(text)
    append!(violations, _semantic_firewall(parsed))
    isempty(violations) || return _fail(join(vcat(["Source/firewall violations:"], violations), "\n  "))

    # 2. Required top-level sections
    missing_sections = [k for k in REQUIRED_KEYS if !haskey(parsed, k)]
    isempty(missing_sections) || return _fail("Missing required sections: $(join(missing_sections, ", "))")

    # 3. Required exact and substring field checks
    errs = String[]
    for (path, expected) in EXACT_CHECKS
        v, ok = _get(parsed, path)
        ok && v == expected || push!(errs, "$path must equal $(repr(expected)), got $(ok ? repr(v) : "<missing>")")
    end
    for (path, needle) in CONTAIN_CHECKS
        v, ok = _get(parsed, path)
        ok && v isa AbstractString && occursin(needle, v) || push!(errs, "$path must contain $(repr(needle)), got $(ok ? repr(v) : "<missing>")")
    end
    append!(errs, _require_ranking(parsed))
    append!(errs, _check_view_weights(parsed))

    # 4. State consistency
    status, ok = _get(parsed, "candidate.grade_status")
    ok && status in ("locked", "frozen_once", "graded") || push!(errs, "candidate.grade_status must be one of locked|frozen_once|graded, got $(ok ? repr(status) : "<missing>")")
    if ok && status == "locked"
        append!(errs, _check_locked(parsed, status, frozen_hash_assignments))
    elseif ok && status in ("frozen_once", "graded")
        append!(errs, _check_frozen_binding(parsed, text, status, frozen_hash_assignments))
    end

    # 5. Expectation flags
    if opt.expect_locked
        (ok && status == "locked") || push!(errs, "--expect-locked requires grade_status = \"locked\", got $(ok ? repr(status) : "<missing>")")
    end
    if opt.expect_frozen_once
        (ok && status == "frozen_once") || push!(errs, "--expect-frozen-once requires grade_status = \"frozen_once\", got $(ok ? repr(status) : "<missing>")")
    end

    isempty(errs) || return _fail(join(vcat(["Contract violations:"], errs), "\n  "))

    # PASS
    prov, _ = _get(parsed, "provenance.status")
    println("Unit-assignment candidate manifest: OK")
    println("  config:        ", opt.config)
    println("  grade_status:  ", status)
    println("  provenance:    ", prov)
    println("  firewall:      pass (forbidden tokens confined to [grader_only])")
    if status in ("frozen_once", "graded")
        fh, _ = _get(parsed, "candidate.frozen_hash")
        println("  frozen_hash:   ", fh)
    end
    println("  status:        ok")
    return 0
end

function main(args=ARGS)
    opt = try
        _parse_cli(args)
    catch e
        println(stderr, "Argument error: ", sprint(showerror, e))
        return 2
    end
    opt === nothing && return 0
    return _run(opt)
end

exit(main())
