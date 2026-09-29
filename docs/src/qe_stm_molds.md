# DFT molds for GlcN/GlcNAc

The assignment compares every lobe patch with simulated STM molds of a central
GlcN or GlcNAc unit in an adsorbed chain. STM probes the local density of states
(LDOS) in the bias window, not a molecular contour. The molds therefore come
from Quantum ESPRESSO LDOS calculations on Cu(100). The scientific details and
their limits are in the [DFT calculation note](dft_calculation_note.md).

## What is used in production

| Item | Location |
|---|---|
| LDOS cubes (accepted, relaxed trimers) | `qe/glcn_restart5/glcn_central_ldos.cube`, `qe/glcnac/glcnac_central_ldos.cube` (untracked; local and Viper) |
| Central-unit frames | `templates/qe_frames/glcn_central_frame.tsv`, `glcnac_central_frame.tsv` |
| Molds | `templates/chitosan_cc_molds_native_v1.tsv` |
| Provenance (all SHA-256) | `templates/chitosan_cc_molds_native_v1.provenance.toml` |
| QE structures and input templates | `hpc/qe_molds/` |

The two systems are GlcN–GlcN–GlcN and GlcN–GlcNAc–GlcN on a three-layer
8×6 Cu(100) slab, with the same substrate registry and central-ring frame, so
the only chemical variable is the acetylation of the central unit.

## Cube → molds

```bash
julia --project=. test/build_cc_molds_native.jl \
    --cube0 qe/glcn_restart5/glcn_central_ldos.cube --frame0 templates/qe_frames/glcn_central_frame.tsv \
    --cube1 qe/glcnac/glcnac_central_ldos.cube --frame1 templates/qe_frames/glcnac_central_frame.tsv \
    --config config/unit_assignment_corroborated_training.toml --out NEW.tsv
```

For each type, the builder samples the cube in the central-unit frame (origin at
the central ring, `t` along the β-(1→4) backbone, `u` transverse towards the C2
substituent). It scans 80 isovalues and keeps the first whose constant-current
surface has a mean height below the 0.50 nm target. Each surface point is the
highest occupied z sample. The result is a 17×17 height map (half-width
0.32 nm, step 0.04 nm), normalized, for each pyranose parity and surface
mirror. It is byte-identical to the Python reference `test/lib/cc_mold_builder.py`
(`test/test_cc_mold_native.jl --real-cubes`). Mean heights: GlcN 0.496 nm
(72/289 supported samples), GlcNAc 0.381 nm (235/289).

`test/score_connected_mold_templates.jl` then decodes each chain against these
molds ([Unit assignment](unit_assignment.md)).

**Known conventions: these are empirical templates, not a validated physical
mapping.**

1. The builder keeps the historical reader's first-index-fast cube addressing,
   so that it replays the champion lineage. QE 7.4.1 writes the last index
   fastest. A correctly addressed variant was graded once (676/870 correct,
   27/145 exact chains, against 679/33 for the legacy molds) and was not adopted.
2. The `plot_num = 5` quantity is a smearing-weighted, potentially signed sum.
   The GlcNAc cube has 131,984 negative values out of 10.8 million.
3. The legacy isovalue bracket includes the rise toward the periodic substrate
   image.

The molds work as data-derived templates on the benchmark. Their coordinates
and heights are not a physical STM simulation
([DFT note](dft_calculation_note.md#Qualification-of-the-STM-observable)).

## Rules

- The sample bias (−0.300 V) comes from the STM experiment. The 0.50 nm height
  and the isovalue rule were declared before any grade.
- The β-(1→4) linkage constrains orientation and parity only. It never imposes
  a GlcNAc count or an alternating sequence.
- Unit truth (`010010`) stays outside mold construction and scoring.
- New cubes (other contexts, larger slabs, converged k-points) make a new,
  versioned mold file with its own provenance. They never overwrite
  `chitosan_cc_molds_native_v1`.

## Producing new cubes

The QE preparation, preflight, relaxation hand-off and frame-extraction scripts
of the 2026 campaign were retired with it. Recover them from the tag
`archive/pre-cleanup-20260929` (`test/prepare_qe_mold_inputs.jl`,
`test/preflight_qe_mold_inputs.jl`, `test/extract_qe_relaxed_xyz.jl`,
`test/extract_qe_mold_frame.jl`, `hpc/submit_qe_molds.sh`,
`hpc/launch_qe_molds_remote.sh`). A new DFT campaign needs a human-approved
scope and time box. Run QE on compute nodes only, keep the two systems
sequential under the one-node QOS limit, and never run it on login nodes.
