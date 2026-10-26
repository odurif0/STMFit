"""
    STMFit

Library shared by the command-line programs in `scripts/`: molecule
consensus, latent-class fusion, unit assignment (transductive and frozen), patch preprocessing,
connected-component molds, empirical Fisher views, benchmark-free reliability of
per-unit signals, joint counts and joint shape fits of repeated scans, and
script helpers. The fit
engines live in the packages under `packages/`. Scripts load a part with, for
example, `using STMFit.MoleculeConsensus`.
"""
module STMFit

include("script_utils.jl")                    # ScriptUtils
include("reconstructed_unit_assignment.jl")   # ReconstructedUnitAssignment
include("unit_assignment_report.jl")          # UnitAssignmentReport
include("empirical_fisher_native.jl")         # EmpiricalFisherNative
include("cc_mold_native.jl")                  # CCMoldNative
include("patch_preprocessing.jl")             # PatchPreprocessing
include("molecule_consensus.jl")              # MoleculeConsensus
include("molecule_fusion.jl")                 # MoleculeFusion
include("frozen_assignment.jl")               # FrozenAssignment
include("unit_reliability.jl")                # UnitReliability
include("joint_count.jl")                     # JointCount
include("joint_shape.jl")                     # JointShape

"Root folder of the STMFit repository (config/, benchmarks/, dft/, data/)."
root() = pkgdir(@__MODULE__)

end # module
