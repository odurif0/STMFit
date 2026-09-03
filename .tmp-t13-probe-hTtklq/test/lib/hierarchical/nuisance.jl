# Per-scan/per-feature median+MAD nuisance profiling without class assignments.
# Zero-MAD scan/view features are invalidated (set to NaN).

struct ScanProfile
    median::Vector{Float64}
    mad::Vector{Float64}
    zero_mad::Vector{Bool}
end

function _median_finite(v::AbstractVector)
    good = filter(isfinite, Float64.(v))
    isempty(good) && return NaN
    return median(good)
end

function _mad_finite(v::AbstractVector, med::Real)
    good = filter(isfinite, Float64.(v))
    isempty(good) && return NaN
    return median(abs.(good .- med))
end

# Global median+MAD profile (per feature column). Used by tests and as a
# building block; the pipeline uses normalize_per_scan.
function profile_scan(X::AbstractMatrix)
    n, p = size(X)
    meds = Vector{Float64}(undef, p)
    mads = Vector{Float64}(undef, p)
    zero = Vector{Bool}(undef, p)
    for j in 1:p
        m = _median_finite(view(X, :, j))
        d = _mad_finite(view(X, :, j), m)
        medv = isfinite(m) ? m : 0.0
        madv = isfinite(d) ? d : 0.0
        meds[j] = medv
        mads[j] = madv
        zero[j] = !isfinite(m) || !isfinite(d) || d <= 0.0
    end
    return ScanProfile(meds, mads, zero)
end

# Per-scan median+MAD normalization. For each (file, feature), compute the
# median and MAD over that scan's finite values; zero-MAD features in a scan
# are invalidated for every lobe of that scan (set to NaN). Returns a matrix
# of the same shape as X.
function normalize_per_scan(records::Vector{LobeRecord}, X::AbstractMatrix,
                            features::Vector{String})
    n, p = size(X)
    @assert length(records) == n
    Z = similar(X)
    files = sort(unique(r.file for r in records))
    for file in files
        idxs = findall(r -> r.file == file, records)
        for j in 1:p
            col = view(X, idxs, j)
            med = _median_finite(col)
            mad = _mad_finite(col, med)
            if !isfinite(med) || !isfinite(mad) || mad <= 0.0
                # Zero-MAD scan/view feature: invalidate for every lobe in scan.
                for i in idxs
                    Z[i, j] = NaN
                end
            else
                for i in idxs
                    v = X[i, j]
                    Z[i, j] = isfinite(v) ? (v - med) / mad : NaN
                end
            end
        end
    end
    return Z
end
