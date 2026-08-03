#!/usr/bin/env julia

const DEFAULT_DOCS = ["docs/src/chitosan_runbook.md", "docs/src/unit_assignment.md"]
const START_MARK = "<!-- UNKNOWN-CHITOSAN-WORKFLOW:START -->"
const END_MARK = "<!-- UNKNOWN-CHITOSAN-WORKFLOW:END -->"
const FORBIDDEN = ["fill the ground truth", "--truth", "--full145", "NKNNKN", "expected_N"]
const TERMINAL_START_MARK = "<!-- T9-TERMINAL-STATUS:START -->"
const TERMINAL_END_MARK = "<!-- T9-TERMINAL-STATUS:END -->"
const TERMINAL_REQUIRED = [
    "terminal `BLOCKED`",
    "13 held-out dates",
    "500-seed",
    "`0.55621530226471905`",
    "`NO_ELIGIBLE_CHALLENGER`",
    "`SKIPPED_NO_ELIGIBLE_CHALLENGER`",
    "grader invocation count was zero",
    "one-shot grade budget remains unused",
    "`grade_status = \"locked\"`",
    "`provenance.status = \"pending\"`",
    "no frozen hash",
    "`9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`",
    "`78.9%` classified physical accuracy",
    "`17/145` exact chains",
    "fixed-denominator honest view",
    "`hierarchical_equalprior`",
    "diagnostic, not frozen or promoted",
    "URI-path",
    "partial-view QC",
    "unstable/nonmonotone fits abstain",
    "records-based pipeline",
    "every consumed primary",
    "stable role, normalized path, and byte SHA-256",
    "resolved views and model options",
    "`--isovalue-scan-intervals` policy (default 1024)",
    "`isovalue_scan_intervals` in provenance",
    "at that declared resolution",
    "one continuous fixed-support root branch",
    "does not unblock T3",
    "rejects multiline TOML strings",
    "exactly one real top-level `[grader_only]` table",
    "case-insensitively",
    "every other non-grader byte",
    "`candidate.frozen_hash`",
    "`candidate.grade_status`",
    "`frozen_once` and `graded`",
    "lifecycle transition does not rehash",
    "not stateless historical proof",
]
const TERMINAL_FORBIDDEN = [
    "candidate was promoted",
    "challenger was promoted",
    "benchmark validated",
    "grade budget used",
    "grade_unit_assignment.jl",
]

function _parse_cli(args)
    docs = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--doc"
            i < length(args) || error("--doc requires a value")
            push!(docs, args[i+1]); i += 2
        elseif startswith(arg, "--doc=")
            push!(docs, split(arg, "="; limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/check_unknown_workflow_docs.jl [--doc PATH]

            Checks UNKNOWN-CHITOSAN-WORKFLOW marked sections for benchmark-only
            strings. With no --doc, checks the maintained runbook and unit docs.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    return isempty(docs) ? DEFAULT_DOCS : docs
end

function _unknown_section(path::String)
    text = read(path, String)
    start = findfirst(START_MARK, text)
    stop = findfirst(END_MARK, text)
    start === nothing && error("Missing start marker in $path")
    stop === nothing && error("Missing end marker in $path")
    first(stop) > last(start) || error("End marker precedes start marker in $path")
    return text[last(start)+1:first(stop)-1]
end

function _marked_section(path::String, start_mark::String, end_mark::String)
    text = read(path, String)
    start = findfirst(start_mark, text)
    stop = findfirst(end_mark, text)
    start === nothing && error("Missing terminal-status start marker in $path")
    stop === nothing && error("Missing terminal-status end marker in $path")
    first(stop) > last(start) || error("Terminal-status end marker precedes start marker in $path")
    return text[last(start)+1:first(stop)-1]
end

function _check_doc(path::String)
    section = _unknown_section(path)
    hits = String[]
    lowered = lowercase(section)
    for token in FORBIDDEN
        haystack = token == "fill the ground truth" ? lowered : section
        needle = token == "fill the ground truth" ? token : token
        occursin(needle, haystack) && push!(hits, token)
    end
    isempty(hits) || error("$path unknown workflow section contains forbidden benchmark-only strings: $(join(hits, ", "))")
    terminal = _marked_section(path, TERMINAL_START_MARK, TERMINAL_END_MARK)
    terminal_flat = replace(terminal, r"\s+" => " ")
    missing = filter(token -> !occursin(token, terminal_flat), TERMINAL_REQUIRED)
    isempty(missing) || error("$path terminal-status section is missing required statements: $(join(missing, ", "))")
    terminal_hits = filter(token -> occursin(token, lowercase(terminal_flat)), TERMINAL_FORBIDDEN)
    isempty(terminal_hits) || error("$path terminal-status section contains a promotion or grading instruction: $(join(terminal_hits, ", "))")
    println("ok\t", path)
    return nothing
end

function main(args=ARGS)
    docs = _parse_cli(args)
    for path in docs
        isfile(path) || error("Doc not found: $path")
        _check_doc(path)
    end
end

main()
