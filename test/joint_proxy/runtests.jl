#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TESTS = [
    "test/joint_proxy/test_proxy_registry.jl",
    "test/joint_proxy/test_candidate_views.jl",
    "test/joint_proxy/test_simulator.jl",
    "test/joint_proxy/test_type_posterior.jl",
    "test/joint_proxy/test_count_calibration.jl",
    "test/joint_proxy/test_type_calibration.jl",
    "test/joint_proxy/test_calibration_cli.jl",
    "test/joint_proxy/test_inference_cli.jl",
    "test/joint_proxy/test_output_validation.jl",
    "test/joint_proxy/test_proxy_swap.jl",
    "test/joint_proxy/test_type_collapse_diagnostic.jl",
    "test/joint_proxy/test_view_preprocessing_diagnostic.jl",
    "test/joint_proxy/test_view_row_drift_diagnostic.jl",
    "test/joint_proxy/test_view_registration_diagnostic.jl",
    "test/joint_proxy/test_mold_support_audit.jl",
    "test/joint_proxy/test_whole_roi_mold_model.jl",
    "test/joint_proxy/test_whole_roi_diagnostic_config.jl",
    "test/joint_proxy/test_whole_roi_common_contrast.jl",
    "test/joint_proxy/test_whole_roi_common_registration.jl",
    "test/joint_proxy/test_whole_roi_frozen_contrast.jl",
    "test/joint_proxy/test_whole_roi_end_to_end.jl",
    "test/joint_proxy/test_frozen_contrast_diagnostic.jl",
    "test/joint_proxy/test_frame_compatibility_audit.jl",
    "test/joint_proxy/test_frame_convergence_audit.jl",
    "test/joint_proxy/test_type_ablation_summary.jl",
]

for script in TESTS
    println("==> ", script)
    run(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(ROOT, script))`)
end
