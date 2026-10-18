"""
    STMFit

Library shared by the command-line programs in `scripts/`: molecule
consensus, latent-class fusion, unit assignment, patch preprocessing,
connected-component molds, empirical Fisher views and script helpers. The fit
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

"Root folder of the STMFit repository (config/, benchmarks/, dft/, data/)."
root() = pkgdir(@__MODULE__)

end # module
