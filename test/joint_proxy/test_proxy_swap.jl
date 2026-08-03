#!/usr/bin/env julia

using Test
using TOML
using SHA
using Statistics
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "proxy_registry.jl"))
include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "type_posterior.jl"))

using .JointProxyRegistry:
    ProxySource, ProxyTemplate, ProxyEntry, ProxyEnsemble,
    RegistryConfigError, load_registry, load_stm_template_tsv,
    canonical_payload_hash, sha256_file
using .JointProxyTypePosterior:
    TypePosteriorLobeEvidence, infer_type_posterior

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")
const GEOMETRIC_SITES = joinpath(ROOT, "templates", "chitosan_geometric_sites.tsv")

_height_from_path(path::AbstractString) = parse(Float64, match(r"h(\d{3})", basename(String(path))).captures[1]) / 100

function _config_prelim_specs()
    cfg = TOML.parsefile(CONFIG)
    sources = cfg["stm_prelim"]["sources"]
    return [(Float64(src["height_nm"]), String(src["path"]), Float64(src["weight"])) for src in sources]
end

function _local_prelim_specs()
    template_dir = joinpath(ROOT, "templates")
    paths = sort(filter(p -> occursin(r"^chitosan_connected_molds_stm_prelim_h\d{3}_half032\.tsv$", basename(p)),
                        readdir(template_dir; join=true));
                 by=_height_from_path)
    return [( _height_from_path(path), path, float(i) ) for (i, path) in enumerate(paths)]
end

function _write_registry_config(path::AbstractString, prelim_specs::Vector{<:Tuple};
                               sites_tsv::AbstractString=GEOMETRIC_SITES)
    cfg = Dict(
        "grid" => Dict("half_nm" => 0.32, "step_nm" => 0.08),
        "geometric" => Dict("sites_tsv" => sites_tsv, "parity_flip" => "t",
                             "mirror_flip" => "u", "normalize" => "zscore"),
        "stm_prelim" => Dict(
            "family_weight" => 0.5,
            "sources" => [Dict("height_nm" => Float64(h), "path" => String(p), "weight" => Float64(w))
                           for (h, p, w) in prelim_specs],
        ),
        "provenance" => Dict("include_bond_templates" => false),
    )
    open(path, "w") do io
        TOML.print(io, cfg)
    end
    return path
end

function _template_mix(a::Vector{Float64}, b::Vector{Float64}, α::Float64)
    v = a .+ α .* b
    n = norm(v)
    n > 0 || throw(ArgumentError("template mix produced a zero vector"))
    return v ./ n
end

function _make_lobes_from_templates(templates::Vector{ProxyTemplate}; n::Int=4)
    lobes = TypePosteriorLobeEvidence[]
    for i in 1:n
        typ = isodd(i) ? 0 : 1
        parity = mod(i - 1, 2)
        mirror = 0
        idx = 1 + typ + 2 * parity + 4 * mirror
        ref = templates[idx].pixels
        alt = templates[1 + (1 - typ) + 2 * parity + 4 * mirror].pixels
        push!(lobes, TypePosteriorLobeEvidence(Dict(
            "fwd" => _template_mix(ref, alt, 0.20 + 0.01 * i),
            "bwd" => _template_mix(ref, alt, 0.10 + 0.01 * i),
        )))
    end
    return lobes
end

function _write_synthetic_template_tsv(path::AbstractString, templates::Vector{ProxyTemplate})
    open(path, "w") do io
        cols = vcat(["name", "type", "parity", "mirror"], [@sprintf("p%03d", i) for i in 1:length(first(templates).pixels)])
        println(io, join(cols, '\t'))
        for t in templates
            name = @sprintf("pair_t%d_p%d_m%d", t.type, t.parity, t.mirror)
            vals = [@sprintf("%.8g", x) for x in t.pixels]
            println(io, join(vcat([name, string(t.type), string(t.parity), string(t.mirror)], vals), '\t'))
        end
    end
    return path
end

function _validate_explicit_pair_provenance(prov)
    for key in (:provider_family, :status, :audited, :source_sha256, :payload_sha256)
        hasproperty(prov, key) || throw(RegistryConfigError("explicit-pair provenance missing $(String(key))"))
    end
    prov.provider_family == "explicit_pair" || throw(RegistryConfigError("explicit-pair provenance family must be explicit_pair"))
    prov.status == "converged" || throw(RegistryConfigError("explicit-pair provenance status must be converged"))
    prov.audited === true || throw(RegistryConfigError("explicit-pair provenance must be audited"))
    for key in (:source_sha256, :payload_sha256)
        val = String(getproperty(prov, key))
        occursin(r"^[0-9a-f]{64}$", val) || throw(RegistryConfigError("explicit-pair provenance $(String(key)) must be sha256"))
    end
    return prov
end

function _make_explicit_pair_fixture(dir::AbstractString, reference::ProxyEntry)
    pair_tsv = joinpath(dir, "synthetic_explicit_pair_converged.tsv")
    _write_synthetic_template_tsv(pair_tsv, reference.templates)
    templates = load_stm_template_tsv(pair_tsv; npix=81)
    src_sha = sha256_file(pair_tsv)
    entry = ProxyEntry(ProxySource("explicit_pair", "synthetic_converged", pair_tsv, src_sha, 0.50, 1.0, true), templates)
    payload = canonical_payload_hash([entry], 0.32, 0.08, 81)
    prov = (; provider_family="explicit_pair", status="converged", audited=true,
            source_sha256=src_sha, payload_sha256=payload)
    _validate_explicit_pair_provenance(prov)
    audit = [
        "[OK] explicit_pair source: $pair_tsv",
        "     sha256=$src_sha",
        "[OK] explicit_pair provenance: converged/audited",
        "     payload_sha256=$payload",
    ]
    return ProxyEnsemble([entry], 0.32, 0.08, 9, 81, payload, audit), prov
end

function _present_stm_source_heights(ens::ProxyEnsemble)
    return [e.source.height_nm for e in ens.entries if e.source.family == "stm_prelim"]
end

function _present_stm_source_weights(ens::ProxyEnsemble)
    return [e.source.weight for e in ens.entries if e.source.family == "stm_prelim"]
end

@testset "proxy_swap: Todo 11" begin
    @testset "geometric-only fresh fixture works with all optional local prelims absent" begin
        mktempdir() do dir
            missing_specs = [(h, joinpath(dir, basename(p)), w) for (h, p, w) in _config_prelim_specs()]
            cfg = _write_registry_config(joinpath(dir, "geo_only.toml"), missing_specs)
            ens = load_registry(cfg)
            @test length(ens.entries) == 1
            @test ens.entries[1].source.family == "geometric"
            @test isapprox(ens.entries[1].source.weight, 1.0; atol=1e-12)
            @test any(occursin("missing optional stm_prelim source", lowercase(a)) for a in ens.audit)
            @test any(occursin("bond/pair templates excluded", lowercase(a)) for a in ens.audit)

            lobes = _make_lobes_from_templates(ens.entries[1].templates)
            got = infer_type_posterior(lobes, ens)
            @test all(isapprox.(sum(got.lobe_marginals, dims=2), ones(length(lobes), 1); atol=1e-12))
            @test isapprox(sum(last.(got.proxy_family_sensitivity)), 1.0; atol=1e-12)
            @test got.proxy_family_sensitivity[1].first == "geometric"
            @test isapprox(got.proxy_family_sensitivity[1].second, 1.0; atol=1e-12)
        end
    end

    @testset "local preliminary STM heights are included, missing gap is audited, and weights renormalize in place" begin
        local_specs = _local_prelim_specs()
        if isempty(local_specs)
            @test true
        else
            mktempdir() do dir
                gap_path = joinpath(dir, "synthetic_missing_h047.tsv")
                augmented = Tuple{Float64,String,Float64}[]
                for (h, p, w) in local_specs
                    push!(augmented, (h, p, w))
                    if isapprox(h, 0.45; atol=1e-8)
                        push!(augmented, (0.47, gap_path, 99.0))
                    end
                end
                cfg = _write_registry_config(joinpath(dir, "with_locals.toml"), augmented)
                ens = load_registry(cfg)

                @test ens.entries[1].source.family == "geometric"
                @test _present_stm_source_heights(ens) == [h for (h, p, w) in local_specs]
                @test isapprox(sum(e.source.weight for e in ens.entries), 1.0; atol=1e-12)
                @test isapprox(ens.entries[1].source.weight, 0.5; atol=1e-12)
                @test any(occursin("missing optional stm_prelim source", lowercase(a)) for a in ens.audit)

                present_weights = [w for (h, p, w) in augmented if isfile(p)]
                expected = [0.5; 0.5 .* (present_weights ./ sum(present_weights))]
                got_weights = [e.source.weight for e in ens.entries]
                @test isapprox.(got_weights, expected; atol=1e-12) |> all

                lobes = _make_lobes_from_templates(ens.entries[1].templates)
                got = infer_type_posterior(lobes, ens)
                @test all(isapprox.(sum(got.lobe_marginals, dims=2), ones(length(lobes), 1); atol=1e-12))
                @test isapprox(sum(last.(got.proxy_family_sensitivity)), 1.0; atol=1e-12)
                @test Set(first.(got.proxy_family_sensitivity)) == Set(["geometric", "stm_prelim"])
            end
        end
    end

    @testset "temporary explicit-pair fixture accepts the same posterior engine and rejects pair-likelihood / malformed provenance" begin
        mktempdir() do dir
            geo_cfg = _write_registry_config(joinpath(dir, "geo_for_pair.toml"), [(h, joinpath(dir, basename(p)), w) for (h, p, w) in _config_prelim_specs()])
            geo = load_registry(geo_cfg)
            reference = geo.entries[1]
            lobes = _make_lobes_from_templates(reference.templates)

            pair_ens, prov = _make_explicit_pair_fixture(dir, reference)
            @test prov.provider_family == "explicit_pair"
            @test pair_ens.entries[1].source.family == "explicit_pair"
            @test any(occursin("explicit_pair source", lowercase(a)) for a in pair_ens.audit)

            got = infer_type_posterior(lobes, pair_ens)
            @test all(isapprox.(sum(got.lobe_marginals, dims=2), ones(length(lobes), 1); atol=1e-12))
            @test isapprox(sum(last.(got.proxy_family_sensitivity)), 1.0; atol=1e-12)
            @test got.proxy_family_sensitivity[1].first == "explicit_pair"
            @test isapprox(got.proxy_family_sensitivity[1].second, 1.0; atol=1e-12)

            @test_throws ArgumentError infer_type_posterior(lobes, pair_ens; pair_bond_evidence=[1.0])
            @test_throws RegistryConfigError _validate_explicit_pair_provenance((; provider_family="explicit_pair", status="draft", audited=false,
                source_sha256="bad", payload_sha256=prov.payload_sha256))
        end
    end
end

function _any_nonpass(ts)
    for f in (:n_fail, :n_failed, :n_errors, :n_error, :n_non_pass)
        if hasproperty(ts, f)
            v = getproperty(ts, f)
            v isa Number && v > 0 && return true
        end
    end
    if hasproperty(ts, :results)
        for c in ts.results
            c isa Test.DefaultTestSet && _any_nonpass(c) && return true
        end
    end
    return false
end

const _ts = Test.get_testset_depth() > 0 ? Test.get_testset() : nothing
exit(_ts === nothing || !_any_nonpass(_ts) ? 0 : 1)
