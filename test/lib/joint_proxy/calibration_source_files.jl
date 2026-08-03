function _calibration_source_files(helper_root::AbstractString=@__DIR__)
    names = [
        "calibration_source_files.jl",
        "calibration_cli.jl", "calibration_cli_types.jl", "calibration_cli_core.jl", "calibration_cli_adapter.jl",
        "proxy_registry.jl", "joint_proxy_types.jl", "joint_proxy_geometry.jl", "joint_proxy_registry_loading.jl",
        "simulator.jl", "simulator_types.jl", "simulator_generation.jl", "simulator_nuisance.jl", "simulator_proxy.jl",
        "candidate_views.jl", "candidate_views_types.jl", "candidate_views_patches.jl", "candidate_views_scoring.jl",
        "count_calibration.jl", "count_calibration_types.jl", "count_calibration_scoring.jl", "count_calibration_fitting.jl",
        "type_posterior.jl", "type_posterior_types.jl", "type_posterior_scoring.jl", "type_posterior_inference.jl",
        "type_calibration.jl", "type_calibration_types.jl", "type_calibration_scoring.jl", "type_calibration_fitting.jl",
    ]
    paths = [joinpath(helper_root, name) for name in names]
    push!(paths, normpath(joinpath(helper_root, "..", "..", "calibrate_joint_proxy_molds.jl")))
    return paths
end
