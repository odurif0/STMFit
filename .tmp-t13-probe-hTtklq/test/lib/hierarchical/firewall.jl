# Firewall for the hierarchical equal-prior unit-assignment model.
#
# Benchmark/truth/lobe-order/sequence/transition/run-length/class-count/
# top-k/occupancy fields may never enter features or options. The patterns
# below are matched case-insensitively against the lowercased feature name
# and are curated so ordinary physical features (`amp_prominence`,
# `sigma_parallel_nm`, `bwd_neg_com_t`, `split_log_skew`) never trigger a
# false positive.

const FORBIDDEN_FEATURE_SUBSTRINGS = String[
    "sequence", "transition", "run_length", "run-length",
    "expected_n", "expected-n", "target_n", "target-n",
    "truth", "control_sequence", "control-sequence",
    "nknnkn", "010010", "101101",
    "centered_pos", "centered-pos", "edge_distance", "edge-distance",
    "lobe_order", "lobe-order", "lobe_index", "lobe-index",
    "lobe_position", "lobe-position", "position_along", "position-along",
    "_position", "-position",
    "class_count", "class-count", "occupancy", "top_k", "top-k", "topk",
    "order", "index", "rank", "lobe_pos", "lobepos",
]

const FORBIDDEN_FEATURE_COMPACT = String[
    replace(pat, r"[^a-z0-9]+" => "") for pat in FORBIDDEN_FEATURE_SUBSTRINGS
]

# Feature columns that are NEVER features even though they appear in the TSV.
# `N` is the chain count (expected N) and must never enter features.
const NON_FEATURE_COLUMNS = Set(["file", "N", "lobe", "source"])

const FORBIDDEN_FLAGS = Set([
    "--truth", "--control-sequence", "--manifest", "--benchmark-manifest",
    "--expected-n", "--expected-n", "--target-n", "--target-n",
    "--full145", "--control",
    "--transitions", "--transition", "--run-length", "--run-lengths",
    "--top-k", "--topk", "--class-count", "--class-counts",
    "--occupancy", "--occupancy-regularizer",
    "--lobe-index", "--lobe-order", "--lobe-position",
    "--sequence",
])

function is_forbidden_feature(name::AbstractString)
    n = lowercase(strip(String(name)))
    isempty(n) && return true
    all(isascii, n) || return true
    compact = replace(n, r"[^a-z0-9]+" => "")
    isempty(compact) && return true
    for (pat, compact_pat) in zip(FORBIDDEN_FEATURE_SUBSTRINGS,
                                  FORBIDDEN_FEATURE_COMPACT)
        (occursin(pat, n) || occursin(compact_pat, compact)) && return true
    end
    return false
end

function is_forbidden_flag(arg::AbstractString)
    name = split(String(arg), '='; limit=2)[1]
    return name in FORBIDDEN_FLAGS
end
