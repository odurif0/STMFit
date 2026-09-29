# QE inputs of the GlcN/GlcNAc mold calculations

Record of the exact structures and Quantum ESPRESSO input templates behind the
accepted LDOS cubes used for the production molds. The calculation, its
parameters and its limits are described in `docs/src/dft_calculation_note.md`;
the cube → mold step is in `docs/src/qe_stm_molds.md`.

| File | Content |
|---|---|
| `glcn_central_trimer.xyz`, `glcnac_central_trimer.xyz` | Initial (unrelaxed) GlcN–GlcN–GlcN and GlcN–GlcNAc–GlcN trimers |
| `*_trimer_atoms.tsv`, `*_trimer_indices.tsv` | Atom tables and central-unit frame indices |
| `*_trimer_slab_pilot.xyz`, `*_slab_pilot_meta.tsv` | Trimer on the three-layer 8×6 Cu(100) pilot slab (the model actually computed) |
| `*_trimer_slab.xyz`, `*_slab_meta.tsv` | Earlier 8×8×4 slab (abandoned: memory) |
| `*.xyz.template` | Placeholder layout of the expected XYZ files |
| `pw_relax.in.template`, `pw_scf.in.template`, `pp_ldos.in.template` | QE input templates (relaxation, SCF, `pp.x` STM LDOS) |

Validation of the initial structures: GlcN `C18 H35 N3 O13` (69 atoms, no
acetyl), GlcNAc `C20 H37 N3 O14` (74 atoms, one central acetyl, no acetylated
neighbour).

The preparation, preflight, relaxation hand-off and frame-extraction scripts
of the completed campaign are archived under the tag
`archive/pre-cleanup-20260929`. Relaxed geometries, QE outputs and cubes are
large and untracked (`qe/`). QE runs on compute nodes only, sequentially under
the one-node QOS limit.
