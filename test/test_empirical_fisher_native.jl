# Focused synthetic tests, Julia 1.13:
#   julia --project=. test/test_empirical_fisher_native.jl
# The numerical reference uses the documented system `python3` with NumPy only.
# It does not fit sklearn, read benchmark labels, or compare saved predictions.
# Production library/CLI never call Python. All fixtures are temporary/synthetic.
using Test
using LinearAlgebra
using Statistics
using Random
using Printf
using TOML
using DelimitedFiles

include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
using .EmpiricalFisherNative
const EF = EmpiricalFisherNative
const ROOT = dirname(@__DIR__)
const CONFIG_PATH = joinpath(ROOT, "config", "unit_assignment_reconstructed.toml")
const CONFIG = TOML.parsefile(CONFIG_PATH)
const OPTIONS = load_fisher_config(CONFIG)
const GRID = fisher_grid(OPTIONS)

function fixture_table(; n=64)
    rng = MersenneTwister(27)
    full = 0.08 .* randn(rng, n, GRID.side^2)
    pattern = [exp(-((row - 9)^2 + (col - 9)^2) / 18) + 0.03 * (col - 9)
               for row in 1:GRID.side for col in 1:GRID.side]
    # Deliberately unbalanced components; both occur in each parity fold.
    for i in 1:n
        full[i, :] .+= (i <= 44 ? -3.0 : 3.0) .* pattern
    end
    table = PatchTable([("synthetic.sxm", i) for i in 1:n], full[:, GRID.disk_indices],
                       full[:, GRID.center_index], fill("", n), GRID)
    return table, full
end

function write_patch_fixture(path, full; keys=[("synthetic.sxm", i) for i in axes(full, 1)],
                             replacements=Dict{Tuple{Int,Int},String}(), missing_column=0,
                             prefix="res", reversed=false)
    pixel_indices = [i for i in axes(full, 2) if i != missing_column]
    open(path, "w") do io
        println(io, join(vcat(["file", "lobe"], [@sprintf("%s_p%03d", prefix, i) for i in pixel_indices]), '\t'))
        order = reversed ? reverse(collect(axes(full, 1))) : axes(full, 1)
        for i in order
            pixels = [get(replacements, (i, j), string(full[i, j])) for j in pixel_indices]
            println(io, join(vcat([keys[i][1], string(keys[i][2])], pixels), '\t'))
        end
    end
end

function reason_of(f)
    try
        f()
        return "NO_ERROR"
    catch error
        error isa EF.FisherFitError || rethrow()
        return error.reason
    end
end

@testset "Explicit Fisher config and Julia runtime" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    @test OPTIONS.pca_components == 10
    @test OPTIONS.noise_regularization == 0.01
    @test OPTIONS.gmm_regularization == 0.001
    @test GRID.side == 17
    @test length(GRID.disk_indices) == 197
    @test GRID.center_index == 145
    for name in ("fisher_noise_regularization", "fisher_gmm_regularization",
                 "fisher_gmm_tolerance", "mold_half_nm", "mold_step_nm")
        for bad in (0.0, -1.0, NaN, Inf, "0.1", true)
            config = deepcopy(CONFIG)
            config["model"][name] = bad
            @test_throws ArgumentError load_fisher_config(config)
        end
    end
    for name in ("fisher_pca_components", "fisher_gmm_maxiter")
        for bad in (0, -1, 2.5, true)
            config = deepcopy(CONFIG)
            config["model"][name] = bad
            @test_throws ArgumentError load_fisher_config(config)
        end
    end
    for bad in (-1, 0.5, true)
        config = deepcopy(CONFIG)
        config["model"]["fisher_seed"] = bad
        @test_throws ArgumentError load_fisher_config(config)
    end
    config = deepcopy(CONFIG)
    config["preprocessing"]["fisher_layout"] = "physical_u"
    @test_throws ArgumentError load_fisher_config(config)
    config = deepcopy(CONFIG)
    config["model"]["mold_step_nm"] = 0.03
    @test_throws ArgumentError load_fisher_config(config)
    config = deepcopy(CONFIG)
    delete!(config["model"], "fisher_noise_regularization")
    @test_throws ArgumentError load_fisher_config(config)
    config = deepcopy(CONFIG)
    config["model"]["fisher_pca_components"] = 1000
    @test_throws ArgumentError load_fisher_config(config)
end

@testset "Forward disk and legacy Python flattening" begin
    # This is exactly the extractor order: physical u outer, t inner.
    flat = Float64[100u + t for u in 1:GRID.side for t in 1:GRID.side]
    disk = flat[GRID.disk_indices]
    flipped = flip_u_disk(disk, GRID)
    @test flip_u_disk(flipped, GRID) == disk
    @test flipped != disk
    @test GRID.mirror_indices[findfirst(==(GRID.center_index), GRID.disk_indices)] ==
          findfirst(==(GRID.center_index), GRID.disk_indices)
    # An independently constructed full row-major reflection keeps u and flips t.
    reflected_full = Float64[100u + (GRID.side + 1 - t) for u in 1:GRID.side for t in 1:GRID.side]
    @test flipped == reflected_full[GRID.disk_indices]
    @test_throws DimensionMismatch flip_u_disk(disk[1:end-1], GRID)

    python = raw"""
import numpy as np
coords = np.arange(-0.32, 0.321, 0.04)
t, u = np.meshgrid(coords, coords, indexing='ij')
disk = t*t + u*u <= 0.32**2
indices = np.arange(17*17).reshape(17, 17)[disk]
flat = np.array([100*u+t for u in range(1,18) for t in range(1,18)], dtype=float)
x = flat.reshape(17,17)[disk]
g = np.zeros(17*17)
g[indices] = x
flipped = g.reshape(17,17)[:,::-1].ravel()[indices]
for row in (coords, indices+1, x, flipped):
    print('\t'.join(format(float(v), '.17g') for v in row))
"""
    reference = split(chomp(read(`python3 -c $python`, String)), '\n')
    values = [parse.(Float64, split(line, '\t')) for line in reference]
    @test GRID.coords == values[1]
    @test GRID.disk_indices == Int.(values[2])
    @test disk == values[3]
    @test flipped == values[4]
    mktempdir() do dir
        path = joinpath(dir, "forward.tsv")
        write_patch_fixture(path, reshape(flat, 1, :); keys=[("/synthetic/path/scan.sxm", 1)])
        loaded = load_patches(path, "res", OPTIONS)
        @test loaded.keys == [("scan.sxm", 1)]
        @test loaded.X[1, :] == disk
        @test loaded.amplitudes == [909.0]
        @test loaded.invalid_reasons == [""]
    end
end

@testset "Full-covariance EM formulas and free weights" begin
    Z = [0.0 0.0; 2.0 1.0; 4.0 -1.0; 6.0 2.0]
    R = [1.0 0.0; 0.9 0.1; 0.5 0.5; 0.4 0.6]
    reg = OPTIONS.gmm_regularization
    weights, means, covariance = EF._m_step(Z, R, reg)
    mass = [2.8, 1.2]
    @test weights ≈ mass ./ 4
    @test weights != [0.5, 0.5]
    for c in 1:2
        mu = sum(R[i, c] .* Z[i, :] for i in 1:4) / mass[c]
        expected = sum(R[i, c] .* ((Z[i, :] - mu) * transpose(Z[i, :] - mu)) for i in 1:4) / mass[c] + reg * I
        @test means[c, :] ≈ mu
        @test covariance[c] ≈ expected
        @test covariance[c][1, 2] != 0
    end
    posterior, average_loglikelihood = EF._e_step(Z, weights, means, covariance)
    densities = [weights[c] * exp(-dot(Z[i, :] - means[c, :], covariance[c] \ (Z[i, :] - means[c, :])) / 2) /
                 sqrt((2pi)^2 * det(covariance[c])) for i in 1:4, c in 1:2]
    total = sum(densities; dims=2)
    @test posterior ≈ densities ./ total
    @test average_loglikelihood ≈ mean(log.(total))
    @test vec(sum(posterior; dims=2)) ≈ ones(4)
end

@testset "Opt-in transverse mirror: same fit, physical u symmetry" begin
    config = TOML.parsefile(joinpath(ROOT, "config", "unit_assignment_transverse_fisher.toml"))
    @test config["model"]["name"] == "cc_soft_transverse_fisher_v1"
    @test config["preprocessing"]["fisher_layout"] == "physical_u_outer_t_inner"
    options = load_fisher_config(config)
    grid = fisher_grid(options)
    config["model"]["name"] = CONFIG["model"]["name"]
    config["preprocessing"]["fisher_layout"] = CONFIG["preprocessing"]["fisher_layout"]
    @test config == CONFIG # No second scientific setting differs.
    @test grid.coords == GRID.coords
    @test grid.disk_indices == GRID.disk_indices
    @test grid.center_index == GRID.center_index
    @test all(>(0), grid.mirror_indices)
    @test grid.mirror_indices[grid.mirror_indices] == collect(eachindex(grid.disk_indices))
    @test grid.mirror_indices != GRID.mirror_indices
    full = Float64[100u + t for u in 1:grid.side for t in 1:grid.side]
    reflected = Float64[100(grid.side + 1 - u) + t for u in 1:grid.side for t in 1:grid.side]
    disk = full[grid.disk_indices]
    @test flip_u_disk(disk, grid) == reflected[grid.disk_indices]
    @test flip_u_disk(flip_u_disk(disk, grid), grid) == disk
    @test sum(abs2, flip_u_disk(disk, grid)) == sum(abs2, disk)

    table, patches = fixture_table()
    legacy_model = fit_fisher(table.X, table.amplitudes, OPTIONS)
    model = fit_fisher(table.X, table.amplitudes, options)
    @test model.w_p == legacy_model.w_p
    @test model.mid == legacy_model.mid
    @test model.gmm.assignments == legacy_model.gmm.assignments
    @test model.gmm.weights == legacy_model.gmm.weights
    for i in (1, 2, 47, 64)
        x = table.X[i, :]
        @test maxmirror_score(x, model, grid) == maxmirror_score(flip_u_disk(x, grid), model, grid)
    end
    # Independent asymmetric probe: legacy t symmetry is not physical u symmetry.
    weight = zeros(length(disk))
    slot = findfirst(i -> grid.mirror_indices[i] != i && GRID.mirror_indices[i] != i,
                     eachindex(disk))
    weight[slot] = 1
    probe = FisherModel(weight, zeros(length(disk)), model.g0, model.g1,
                        model.amplitude_means, model.gmm, nothing)
    x = zeros(length(disk)); x[slot] = 2
    @test maxmirror_score(x, probe, grid) == 2
    @test maxmirror_score(flip_u_disk(x, grid), probe, grid) == 2
    @test maxmirror_score(flip_u_disk(x, grid), probe, GRID) == 0

    native = PatchTable(table.keys, table.X, table.amplitudes, table.invalid_reasons, grid)
    scores = cv_scores(native, options)
    for parity in (0, 1)
        train = findall(k -> mod(k[2], 2) == parity, table.keys)
        held = findall(k -> mod(k[2], 2) != parity, table.keys)
        fold = fit_fisher(table.X[train, :], table.amplitudes[train], OPTIONS)
        @test [scores[i].score for i in held] ==
              [max(score(table.X[i, :], fold), score(flip_u_disk(table.X[i, :], grid), fold)) for i in held]
    end
    bad = copy(table.X); bad[2, 1] = NaN
    invalid = cv_scores(PatchTable(table.keys, bad, table.amplitudes, table.invalid_reasons, grid), options)
    @test length(invalid) == length(scores)
    @test invalid[2].invalid_reason == "nonfinite_patch_input"
    @test isnan(invalid[2].score)
end

@testset "PCA/sample-noise/Fisher algebra versus fixed-partition NumPy" begin
    rng = MersenneTwister(203)
    X = randn(rng, 40, 18)
    X[1:13, :] .-= 1.3
    X[14:end, :] .+= 2.1
    X .+= transpose(collect(1.0:18.0)) # nonzero original mean, intentionally
    amps = vcat(fill(0.7, 13), fill(2.2, 27))
    partition = vcat(fill(2, 13), fill(1, 27))
    Xc = X .- mean(X; dims=1)
    V = svd(Xc; full=false).V[:, 1:OPTIONS.pca_components]
    Z = Xc * V
    linear = EF._fisher_linear(Xc, Z, V, amps, partition, OPTIONS.noise_regularization)
    @test collect(linear.amplitude_means) ≈ [0.7, 2.2]
    @test linear.g0 ≈ vec(mean(Xc[1:13, :]; dims=1))
    @test linear.g1 ≈ vec(mean(Xc[14:end, :]; dims=1))
    @test norm(linear.mid) < 1e-13
    # Full latent sample covariance, not pooled residual covariance or ML /n.
    noise = cov(Z; dims=1, corrected=true) + OPTIONS.noise_regularization * I
    analytic = V * (noise \ (transpose(V) * (linear.g1 - linear.g0)))
    @test linear.w_p ≈ analytic
    swapped = EF._fisher_linear(Xc, Z, V, amps, 3 .- partition, OPTIONS.noise_regularization)
    @test swapped.w_p == linear.w_p
    @test swapped.amplitude_means == linear.amplitude_means
    @test reason_of(() -> EF._fisher_linear(Xc, Z, V, ones(40), partition, OPTIONS.noise_regularization)) == "amplitude_tie"
    @test reason_of(() -> EF._fisher_linear(Xc, Z, V, amps, ones(Int, 40), OPTIONS.noise_regularization)) == "empty_hard_component"
    python = raw"""
import numpy as np, sys
X = np.loadtxt(sys.argv[1], delimiter='\t')
a = np.loadtxt(sys.argv[2], dtype=int)
amps = np.loadtxt(sys.argv[3])
k, reg = int(sys.argv[4]), float(sys.argv[5])
Xc = X-X.mean(0)
_, _, Vt = np.linalg.svd(Xc, full_matrices=False)
V = Vt[:k]
Z = Xc @ V.T
mu = np.array([Xc[a==c].mean(0) for c in (1,2)])
ma = [amps[a==c].mean() for c in (1,2)]
order = np.argsort(ma)
g0, g1 = mu[order[0]], mu[order[1]]
mid = Xc.mean(0)
cov = np.cov(Z.T) + reg*np.eye(k)
z0, z1 = (g0-mid) @ V.T, (g1-mid) @ V.T
w = V.T @ np.linalg.solve(cov,z1-z0)
for row in (w, mid, g0, g1, (X-mid) @ w):
    print('\t'.join(format(float(v), '.17g') for v in row))
"""
    mktempdir() do dir
        xfile, afile, ampfile = [joinpath(dir, name) for name in ("x.tsv", "a.tsv", "amps.tsv")]
        writedlm(xfile, X, '\t')
        writedlm(afile, partition)
        writedlm(ampfile, amps)
        text = read(`python3 -c $python $xfile $afile $ampfile $(OPTIONS.pca_components) $(OPTIONS.noise_regularization)`, String)
        numerical = [parse.(Float64, split(line, '\t')) for line in split(chomp(text), '\n')]
        @test linear.w_p ≈ numerical[1] rtol=1e-11 atol=1e-12
        @test linear.mid ≈ numerical[2] atol=1e-13
        @test linear.g0 ≈ numerical[3] rtol=1e-12 atol=1e-12
        @test linear.g1 ≈ numerical[4] rtol=1e-12 atol=1e-12
        @test (X .- transpose(linear.mid)) * linear.w_p ≈ numerical[5] rtol=1e-11 atol=1e-11
    end
end

@testset "Deterministic fit, amplitude convention, raw midpoint and reflection" begin
    table, full = fixture_table()
    one = fit_fisher(table.X, table.amplitudes, OPTIONS)
    two = fit_fisher(table.X, table.amplitudes, OPTIONS)
    @test one.w_p == two.w_p
    @test one.mid == two.mid
    @test one.gmm.assignments == two.gmm.assignments
    @test one.gmm.weights == two.gmm.weights
    @test one.gmm.converged
    @test one.amplitude_means[1] < one.amplitude_means[2]
    @test one.amplitude_means[1] < 0 < one.amplitude_means[2]
    @test sort(one.gmm.weights) ≈ [20 / 64, 44 / 64] atol=1e-8
    @test length(one.w_p) == 197
    @test length(one.gmm.covariances) == 2
    @test size(one.gmm.covariances[1]) == (10, 10)
    @test norm(one.mid) < 1e-12
    @test norm(vec(mean(table.X; dims=1))) > 1
    @test score(table.X[1, :], one) < 0
    @test score(table.X[end, :], one) > 0
    for i in (1, 2, 47, 64)
        x = table.X[i, :]
        @test score(x, one) == dot(x - one.mid, one.w_p)
        @test !isapprox(score(x, one), dot(x - vec(mean(table.X; dims=1)), one.w_p))
        expected = max(score(x, one), score(flip_u_disk(x, GRID), one))
        @test maxmirror_score(x, one, GRID) == expected
        @test maxmirror_score(flip_u_disk(x, GRID), one, GRID) == expected
    end
    # maxmirror keeps a negative maximum; it is NOT max(abs(original), abs(mirror)).
    slot = findfirst(i -> GRID.mirror_indices[i] != i, eachindex(GRID.disk_indices))
    negative = zeros(length(GRID.disk_indices))
    negative[slot] = -1.0
    negative[GRID.mirror_indices[slot]] = -2.0
    weight = zeros(length(negative))
    weight[slot] = 1.0
    simple = FisherModel(weight, zeros(length(weight)), one.g0, one.g1, one.amplitude_means, one.gmm, nothing)
    @test score(negative, simple) == -1.0
    @test score(flip_u_disk(negative, GRID), simple) == -2.0
    @test maxmirror_score(negative, simple, GRID) == -1.0

    limited_config = deepcopy(CONFIG)
    limited_config["model"]["fisher_gmm_maxiter"] = 1
    limited_options = load_fisher_config(limited_config)
    limited_model = fit_fisher(table.X, table.amplitudes, limited_options)
    @test !limited_model.gmm.converged
    @test limited_model.gmm.iterations == 1
    limited_rows = @test_logs (:warn, r"^Fisher GMM reached") (:warn, r"^Fisher GMM reached") begin
        cv_scores(table, limited_options)
    end
    @test all(row -> isfinite(row.score) && isempty(row.invalid_reason), limited_rows)
    @test reason_of(() -> fit_fisher(table.X[1:9, :], table.amplitudes[1:9], OPTIONS)) == "too_few_training_rows"
    @test reason_of(() -> fit_fisher(ones(20, 197), ones(20), OPTIONS)) == "no_patch_variation"
    @test reason_of(() -> fit_fisher(table.X, ones(64), OPTIONS)) == "amplitude_tie"
    bad = copy(table.X)
    bad[1, 1] = NaN
    @test reason_of(() -> fit_fisher(bad, table.amplitudes, OPTIONS)) == "nonfinite_training_input"
    @test_throws DimensionMismatch fit_fisher(table.X, [1.0], OPTIONS)
end

@testset "Even/odd CV independence and missing folds retain all keys" begin
    table, _ = fixture_table()
    rows = cv_scores(table, OPTIONS)
    @test length(rows) == length(table.keys)
    @test all(row -> isempty(row.invalid_reason) && isfinite(row.score), rows)
    odd = findall(isodd, last.(table.keys))
    even = findall(iseven, last.(table.keys))
    odd_model = fit_fisher(table.X[odd, :], table.amplitudes[odd], OPTIONS)
    even_model = fit_fisher(table.X[even, :], table.amplitudes[even], OPTIONS)
    @test [rows[i].score for i in even] == [maxmirror_score(table.X[i, :], odd_model, GRID) for i in even]
    @test [rows[i].score for i in odd] == [maxmirror_score(table.X[i, :], even_model, GRID) for i in odd]
    # Mutating one even patch must not change other even scores: neither its PCA
    # nor its amplitude mapping may enter the model trained on odd rows.
    changed_X, changed_amp = copy(table.X), copy(table.amplitudes)
    changed_X[2, :] .+= 0.75 .* collect(range(-1.0, 1.0; length=size(changed_X, 2)))
    changed_amp[2] += 20
    changed = PatchTable(table.keys, changed_X, changed_amp, table.invalid_reasons, GRID)
    changed_rows = cv_scores(changed, OPTIONS)
    @test [rows[i].score for i in even if i != 2] == [changed_rows[i].score for i in even if i != 2]
    @test rows[2].score != changed_rows[2].score
    # Training amplitudes can reverse mapping; held-out amplitudes cannot.
    reversed_amp = copy(table.amplitudes)
    reversed_amp[odd] .*= -1
    flipped = fit_fisher(table.X[odd, :], reversed_amp[odd], OPTIONS)
    @test flipped.w_p ≈ -odd_model.w_p
    @test flipped.gmm.assignments == odd_model.gmm.assignments
    @test flipped.gmm.weights == odd_model.gmm.weights

    odd_only = PatchTable(table.keys[odd], table.X[odd, :], table.amplitudes[odd], fill("", length(odd)), GRID)
    absent = cv_scores(odd_only, OPTIONS)
    @test length(absent) == length(odd)
    @test all(row -> isnan(row.score) && row.invalid_reason == "training_even:too_few_training_rows", absent)
    too_small = PatchTable(table.keys[1:12], table.X[1:12, :], table.amplitudes[1:12], fill("", 12), GRID)
    too_small_rows = cv_scores(too_small, OPTIONS)
    @test length(too_small_rows) == 12
    @test all(row -> isnan(row.score) && endswith(row.invalid_reason, ":too_few_training_rows"), too_small_rows)
    constant = PatchTable(table.keys, zeros(size(table.X)), zeros(64), fill("", 64), GRID)
    degenerate = cv_scores(constant, OPTIONS)
    @test length(degenerate) == 64
    @test all(row -> isnan(row.score) && endswith(row.invalid_reason, ":no_patch_variation"), degenerate)
    tied = PatchTable(table.keys, table.X, ones(64), fill("", 64), GRID)
    @test all(row -> endswith(row.invalid_reason, ":amplitude_tie"), cv_scores(tied, OPTIONS))
    # Invalid held-out rows never contaminate fitting or disappear from output.
    bad = copy(table.X)
    bad[2, 1] = NaN
    invalid = PatchTable(table.keys, bad, table.amplitudes, fill("", 64), GRID)
    badrows = cv_scores(invalid, OPTIONS)
    @test badrows[2].invalid_reason == "nonfinite_patch_input"
    @test isnan(badrows[2].score)
    @test [badrows[i].score for i in even if i != 2] == [rows[i].score for i in even if i != 2]
end

@testset "Patch TSV validation, all-key coverage and six-decimal serialization" begin
    table, full = fixture_table()
    mktempdir() do dir
        input = joinpath(dir, "patches.tsv")
        # An invalid full-patch corner outside the disk still invalidates its row.
        replacements = Dict((1, 1) => "NA", (2, 2) => "NaN", (3, 3) => "Inf", (4, 4) => "bad")
        write_patch_fixture(input, full; replacements=replacements, reversed=true)
        loaded = load_patches(input, "res", OPTIONS)
        @test loaded.keys == table.keys
        @test startswith(loaded.invalid_reasons[1], "invalid_patch_value")
        @test startswith(loaded.invalid_reasons[2], "nonfinite_patch_value")
        @test startswith(loaded.invalid_reasons[3], "nonfinite_patch_value")
        @test startswith(loaded.invalid_reasons[4], "invalid_patch_value")
        @test loaded.X[5:end, :] == table.X[5:end, :]
        result = cv_scores(loaded, OPTIONS)
        @test length(result) == 64
        @test all(row -> isnan(row.score) && !isempty(row.invalid_reason), result[1:4])
        @test all(row -> isfinite(row.score) && isempty(row.invalid_reason), result[5:end])
        output = joinpath(dir, "scores.tsv")
        write_scores(output, result)
        lines = readlines(output)
        @test lines[1] == "file\tlobe\tscore\tinvalid_reason"
        @test length(lines) == 65
        @test all(occursin("\tNA\t", line) for line in lines[2:5])
        @test all(occursin(r"\t-?\d+\.\d{6}\t$", line) for line in lines[6:end])
        @test_throws ArgumentError write_scores(output, result)
        @test readlines(output) == lines

        missing = joinpath(dir, "missing-column.tsv")
        write_patch_fixture(missing, full; missing_column=200)
        missing_rows = cv_scores(load_patches(missing, "res", OPTIONS), OPTIONS)
        @test length(missing_rows) == 64
        @test all(row -> row.invalid_reason == "missing_patch_column:res_p200", missing_rows)
        short = joinpath(dir, "short-row.tsv")
        write_patch_fixture(short, full[1:1, :])
        header = first(readlines(short))
        write(short, header * "\nshort.sxm\t7\t1.0\n")
        short_rows = cv_scores(load_patches(short, "res", OPTIONS), OPTIONS)
        @test length(short_rows) == 1
        @test (short_rows[1].file, short_rows[1].lobe) == ("short.sxm", 7)
        @test short_rows[1].invalid_reason == "invalid_patch_value:res_p002"
        @test isnan(short_rows[1].score)
        empty = joinpath(dir, "empty.tsv")
        write_patch_fixture(empty, zeros(0, 289))
        @test isempty(cv_scores(load_patches(empty, "res", OPTIONS), OPTIONS))
        schema = joinpath(dir, "schema.tsv")
        write(schema, "file\tlobe\n")
        @test isempty(load_patches(schema, "res", OPTIONS).keys)
        write(schema, "")
        @test_throws ArgumentError load_patches(schema, "res", OPTIONS)
        write(schema, "file\nscan.sxm\n")
        @test_throws ArgumentError load_patches(schema, "res", OPTIONS)
        write(schema, "file\tlobe\tlobe\n")
        @test_throws ArgumentError load_patches(schema, "res", OPTIONS)
        for bad_lobe in ("0", "-1", "NaN", "")
            write(schema, "file\tlobe\nscan.sxm\t$bad_lobe\n")
            @test_throws ArgumentError load_patches(schema, "res", OPTIONS)
        end
        write(schema, "file\tlobe\n\t1\n")
        @test_throws ArgumentError load_patches(schema, "res", OPTIONS)
        write(schema, "file\tlobe\n/a/scan.sxm\t1\n/b/scan.sxm\t1\n")
        @test_throws ArgumentError load_patches(schema, "res", OPTIONS)
        write(schema, "file\tlobe\nscan.sxm\t1\textra\n")
        @test_throws ArgumentError load_patches(schema, "res", OPTIONS)
        @test_throws ArgumentError load_patches(input, "", OPTIONS)

        roundfile = joinpath(dir, "rounded.tsv")
        rounding = [FisherScore("s.sxm", 1, 1.23456789, ""), FisherScore("s.sxm", 2, -0.0000001, ""),
                    FisherScore("s.sxm", 3, NaN, ""), FisherScore("s.sxm", 4, 8.0, "invalid_fixture")]
        write_scores(roundfile, rounding)
        @test readlines(roundfile)[2:end] == ["s.sxm\t1\t1.234568\t", "s.sxm\t2\t-0.000000\t",
                                            "s.sxm\t3\tNA\tnonfinite_fisher_score", "s.sxm\t4\tNA\tinvalid_fixture"]
        # Tie-to-even and negative zero agree with Python's historical .6f output.
        python = "print(format(1.23456789, '.6f')); print(format(-0.0000001, '.6f')); print(format(0.0078125, '.6f'))"
        rounded = split(chomp(read(`python3 -c $python`, String)), '\n')
        @test rounded == [@sprintf("%.6f", x) for x in (1.23456789, -0.0000001, 0.0078125)]
    end
end

@testset "Native Julia CLI output, reproducibility and error messages" begin
    mktempdir() do dir
        table, full = fixture_table()
        input = joinpath(dir, "patches.tsv")
        write_patch_fixture(input, full; reversed=true)
        cli = joinpath(@__DIR__, "build_empirical_fisher_native.jl")
        julia = Base.julia_cmd()
        output = joinpath(dir, "first.tsv")
        second = joinpath(dir, "second.tsv")
        help = read(`$julia --project=$ROOT $cli --help`, String)
        @test occursin("byte identity is NOT claimed", help)
        @test occursin("--patches PATH --prefix res --config CONFIG --out PATH", help)
        message = read(`$julia --project=$ROOT $cli --patches $input --prefix res --config $CONFIG_PATH --out $output`, String)
        @test occursin("64 rows, 0 invalid", message)
        read(`$julia --project=$ROOT $cli --patches $input --prefix res --config $CONFIG_PATH --out $second`, String)
        @test read(output) == read(second)
        reference = joinpath(dir, "inprocess.tsv")
        write_scores(reference, cv_scores(table, OPTIONS))
        @test read(output) == read(reference)
        log = IOBuffer()
        process = run(pipeline(ignorestatus(`$julia --project=$ROOT $cli --patches $input --prefix res --config $CONFIG_PATH --out $output`), stdout=log, stderr=log))
        @test process.exitcode == 1
        @test occursin("output already exists", String(take!(log)))
        @test read(output) == read(second)
        config = deepcopy(CONFIG)
        config["model"]["mold_step_nm"] = 0.0
        badconfig = joinpath(dir, "bad.toml")
        open(io -> TOML.print(io, config), badconfig, "w")
        badoutput = joinpath(dir, "invalid.tsv")
        process = run(pipeline(ignorestatus(`$julia --project=$ROOT $cli --patches $input --prefix res --config $badconfig --out $badoutput`), stdout=log, stderr=log))
        @test process.exitcode == 1
        @test occursin("mold_step_nm must be positive and finite", String(take!(log)))
        @test !ispath(badoutput)
        process = run(pipeline(ignorestatus(`$julia --project=$ROOT $cli --bogus x`), stdout=log, stderr=log))
        @test process.exitcode == 1
        @test occursin("unknown option", String(take!(log)))
    end
end
