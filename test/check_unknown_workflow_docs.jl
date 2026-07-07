#!/usr/bin/env julia

const DEFAULT_DOCS = ["docs/src/chitosan_runbook.md", "docs/src/unit_assignment.md"]
const START_MARK = "<!-- UNKNOWN-CHITOSAN-WORKFLOW:START -->"
const END_MARK = "<!-- UNKNOWN-CHITOSAN-WORKFLOW:END -->"
const FORBIDDEN = ["fill the ground truth", "--truth", "--full145", "NKNNKN", "expected_N"]

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
