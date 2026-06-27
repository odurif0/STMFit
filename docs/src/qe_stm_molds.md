# Quantum ESPRESSO STM Molds

This page describes how to produce STM/LDOS molds for GlcN/GlcNAc unit
assignment using Quantum ESPRESSO on an HPC system such as MPCDF Raven/Viper.

## Objective

The desired mold is not a molecular contour. STM is sensitive to the local density
of states, the bias window, the tip, the adsorbed geometry, and the Cu(100)
substrate. The useful input for STMFit is therefore a simulated STM/LDOS map in
the same local frame as the experimental lobe patches:

```text
type    t_nm    u_nm    value
0       -0.48   -0.48   ...
1       -0.48   -0.48   ...
```

`type=0` is GlcN and `type=1` is GlcNAc. The `t` axis is along the beta-(1->4)
backbone, `u` is transverse, and the origin is the central fitted lobe.

## Recommended First Calculation

Start with two linked trimers on Cu(100):

```text
X-GlcN-X
X-GlcNAc-X
```

Use the same neighbor choice `X` in both calculations. This tests whether the
central unit has a detectable STM signature without yet paying for all context
combinations. If the imported DFT-STM mold improves patch-level scoring, extend
to all eight left/central/right contexts.

## QE Workflow

For the 240817 chitosan scans, the SXM headers contain a uniform bias of
`-0.300 V` over the folder. The first prepared QE input directories therefore use
`emin=-0.3 eV`, `emax=0.0 eV` for occupied-state LDOS export. These run
directories live under `qe/` and are gitignored because QE outputs can be large.

Before submission, validate the prepared inputs with:

```bash
julia --project=. test/preflight_qe_mold_inputs.jl \
    --dir qe/glcn \
    --dir qe/glcnac \
    --out hpc/qe_molds/qe_input_preflight.tsv \
    --max-total-tasks 8 \
    --sequential
```

`--sequential` treats the two run directories as an `afterok` chain and enforces
the task budget as the **maximum simultaneous** count (8), not the sum (16), so
two 8-task jobs fit within a single-node QOS group limit.

Then submit on the HPC system with:

```bash
bash hpc/submit_qe_molds.sh --watch --sequential
```

From the local workstation, after configuring `hpc/remote.env`, use:

```bash
bash hpc/launch_qe_molds_remote.sh --dry-run
bash hpc/launch_qe_molds_remote.sh --watch
```

This syncs code and the prepared QE input files, runs remote preflight, then
submits the two QE jobs. The remote launcher uses `SSH_CONNECT_TIMEOUT=180` by
default; increase it in `hpc/remote.env` if the gateway or password/OTP flow is
slow.

After the QE jobs finish, convert the relaxed geometries and cubes to STMFit mold
templates with an explicit sampling height:

```bash
julia --project=. test/finalize_qe_mold_workflow.jl \
    --height-nm HEIGHT \
    --glcn-dir qe/glcn \
    --glcnac-dir qe/glcnac
```

The script extracts the two typed frames, samples the two LDOS cubes, writes
`templates/chitosan_stm_maps.tsv`, and imports connected unary/bond templates.

Scripts and artifacts involved (the operational runbook with job history and
Raven/Viper module details lives in `hpc/qe_molds/README.md`):

| Artifact / script | Role |
|---|---|
| `test/build_initial_chitosan_trimer_xyz.jl` | Deterministic initial X-GlcN-X / X-GlcNAc-X trimer XYZs (starting point for relaxation, not optimized). |
| `test/build_qe_slab_trimer_xyz.jl` | Add a reproducible Cu(100) slab below an oriented trimer XYZ. |
| `test/prepare_qe_mold_inputs.jl` | Serialize a vetted slab+trimer XYZ into `pw_relax.in`, `pw_scf.in`, `pp_ldos.in`, and `run_qe_mold.sbatch`. |
| `test/preflight_qe_mold_inputs.jl` | Validate prepared run dirs (nat/ntyp, pseudos, Slurm tasks/memory, `--sequential` budget). |
| `hpc/submit_qe_molds.sh` | Preflight + submit the prepared run dirs on the cluster. |
| `hpc/launch_qe_molds_remote.sh` | Local sync + remote preflight/submission wrapper. |
| `test/finalize_qe_mold_workflow.jl` | Post-QE: relaxed XYZ + cube → STMFit mold templates (canonical end-to-end path). |
| `hpc/qe_mods/*_trimer.xyz`, `*_slab*.xyz` | Generated geometries (pilot: `8×6×3` slab). |

Before running real QE jobs, the local helper chain can be checked without QE:

```bash
julia --project=. test/smoke_qe_mold_workflow.jl
```

This creates a temporary synthetic trimer/slab, prepares QE inputs via
`--cell-metadata`, simulates the relax-to-SCF handoff, extracts a frame, samples
two synthetic cube files with typed `--frame` arguments, and imports connected STM
molds. It does not validate chemistry or run Quantum ESPRESSO.

The initial trimer geometries can be regenerated with:

```bash
julia --project=. test/build_initial_chitosan_trimer_xyz.jl \
    --out-dir hpc/qe_molds
```

Those files are deterministic starting geometries for QE relaxation. They are not
DFT-optimized structures. Both use GlcN neighbors so the first comparison is
`GlcN-GlcN-GlcN` versus `GlcN-GlcNAc-GlcN`. Companion `*_indices.tsv` files give
the central-ring frame indices.

Validate the generated structures with:

```bash
julia --project=. test/validate_chitosan_trimer_structures.jl \
    --dir hpc/qe_molds \
    --out hpc/qe_molds/structure_validation.tsv
```

The current report confirms that `glcn` has 0 acetyl units and `glcnac` has 1
central acetyl unit, with GlcN neighbors in both structures. After the default
`8×8×4` slab, the validated frame command is:

```text
--origin-indices 268,269,270,271,272,273 --axis-from 271 --axis-to 268 --plane-index 277
```

Once a chemically vetted XYZ exists, generate the concrete QE input set with:

If the trimer XYZ is separate from the Cu(100) slab, first assemble a simple slab
model around it:

```bash
julia --project=. test/build_qe_slab_trimer_xyz.jl \
    --molecule hpc/qe_molds/glcn_central_trimer.xyz \
    --out hpc/qe_molds/glcn_central_trimer_slab.xyz \
    --metadata hpc/qe_molds/glcn_central_trimer_slab_meta.tsv \
    --nx 8 \
    --ny 8 \
    --layers 4 \
    --center-indices 12,13,14,15,16,17 \
    --height-above-top 2.6 \
    --vacuum 18.0
```

`--center-indices 12,13,14,15,16,17` centers the central pyranose ring over the
same Cu(100) registry in both GlcN and GlcNAc calculations. The slab helper does
not rotate the beta-(1->4) trimer; it only places an already oriented molecule
above a reproducible Cu(100) slab.

Then generate QE inputs:

```bash
julia --project=. test/prepare_qe_mold_inputs.jl \
    --xyz hpc/qe_molds/glcn_central_trimer_slab.xyz \
    --cell-metadata hpc/qe_molds/glcn_central_trimer_slab_meta.tsv \
    --out-dir qe/glcn \
    --prefix glcn_central \
    --fix-below-z ZCUT \
    --emin-ev EMIN \
    --emax-ev EMAX
```

Repeat with `type=1` / `glcnac_central` geometry. The helper writes
`pw_relax.in`, `pw_scf.in`, `pp_ldos.in`, and `run_qe_mold.sbatch`.
`--cell-metadata` consumes the `cell_a/b/c` rows written by
`build_qe_slab_trimer_xyz.jl`; explicit `--cell-a/--cell-b/--cell-c` overrides
remain available.

The conceptual QE-level sequence (the generated `run_qe_mold.sbatch` performs
this automatically with explicit MPI task counts) is:

```bash
module load intel/2024.0 impi/2021.11 qe/7.4.1
export OMP_NUM_THREADS=1
QE_NTASKS=${SLURM_NTASKS:-8}
srun -n "$QE_NTASKS" --cpu-bind=cores pw.x -in pw_relax.in > glcn_central_relax.out
julia --project=/path/to/STMFit /path/to/STMFit/test/extract_qe_relaxed_xyz.jl \
    --qe-out glcn_central_relax.out \
    --out glcn_central_relaxed.xyz \
    --metadata glcn_central_relaxed_meta.tsv
julia --project=/path/to/STMFit /path/to/STMFit/test/update_qe_positions_from_xyz.jl \
    --input pw_scf.in \
    --xyz glcn_central_relaxed.xyz \
    --out pw_scf_relaxed.in
srun -n "$QE_NTASKS" --cpu-bind=cores pw.x -in pw_scf_relaxed.in > glcn_central_scf.out
srun -n "$QE_NTASKS" --cpu-bind=cores pp.x -in pp_ldos.in > glcn_central_pp.out
```

Repeat for the GlcNAc-central trimer. The generated `run_qe_mold.sbatch` sets
`#SBATCH --ntasks-per-node=8 --cpus-per-task=1`, derives `QE_NTASKS` from
`$SLURM_NTASKS`, and runs the relax → XYZ extract → SCF-update → SCF → PP chain.
After launch, verify the QE header reports `Number of MPI processes: 8`. Set
`STMFIT_ROOT=/path/to/STMFit` if the default relative path is not correct for
your QE run directory. If the QE output reports `ATOMIC_POSITIONS (alat)` or
`CELL_PARAMETERS (alat)`, pass `--alat-angstrom` to
`extract_qe_relaxed_xyz.jl`; angstrom and bohr units are converted automatically.

## Frame Extraction

After relaxation, choose central-unit atom indices in the relaxed XYZ to define
the local frame used for map extraction:

```bash
julia --project=. test/extract_qe_mold_frame.jl \
    --xyz qe/glcn/glcn_central_relaxed.xyz \
    --origin-indices I,J \
    --axis-from I \
    --axis-to J \
    --plane-index K \
    --height-nm HEIGHT \
    --out qe/glcn/frame.tsv
```

The script writes `origin_nm`, `t_axis`, `u_axis`, `normal_axis`, and
`height_nm`. Pass the resulting file directly to `cube_to_stm_maps.jl` with a
type prefix:

```text
--frame 0:qe/glcn/frame.tsv
```

The legacy explicit form `--origin OX,OY,OZ --t-axis TX,TY,TZ --u-axis UX,UY,UZ
--height-nm HEIGHT` is still accepted. The atom indices are 1-based XYZ atom line
indices. `origin-indices` should identify the central lobe/ring center, the axis
atoms should follow the beta-(1->4) backbone, and `plane-index` should lie on the
positive C2-substituent side.

## Cube To STMFit Map TSV

After `pp.x` writes cube files, extract the central-unit patch plane with:

```bash
julia --project=. test/cube_to_stm_maps.jl \
    --cube 0:qe/glcn/glcn_central_ldos.cube \
    --frame 0:qe/glcn/frame.tsv \
    --cube 1:qe/glcnac/glcnac_central_ldos.cube \
    --frame 1:qe/glcnac/frame.tsv \
    --cube-units bohr \
    --out templates/chitosan_stm_maps.tsv
```

The origins are in nm in the cube coordinate system, and axes are unit vectors in
that same coordinate system. `height_nm` samples a constant-height plane along
`normal = t_axis × u_axis`. If `pp.x` already exports a physical STM-like plane
rather than a volumetric LDOS cube, set the frame origin and height so the
sampled plane matches that map.

## Import And Score

Convert the map TSV into connected templates:

```bash
julia --project=. test/import_stm_mold_maps.jl \
    --maps templates/chitosan_stm_maps.tsv \
    --out templates/chitosan_connected_molds_stm.tsv \
    --bond-out templates/chitosan_connected_bond_molds_stm.tsv \
    --half-nm 0.48 \
    --step-nm 0.08
```

Then decode the experimental patches:

```bash
julia --project=. test/score_connected_mold_templates.jl \
    --patches results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --templates templates/chitosan_connected_molds_stm.tsv \
    --bond-templates templates/chitosan_connected_bond_molds_stm.tsv \
    --prefix raw_p \
    --template-mode contrast \
    --out results/unit_assignment/stm_mold_predictions.tsv
```

Grade only after the mold has been frozen:

```bash
julia --project=. test/grade_unit_assignment.jl \
    --predictions results/unit_assignment/stm_mold_predictions.tsv \
    --truth benchmarks/chitosan_240817_unit_sequences.tsv
```

## Rules

- Choose the bias window from the STM experiment, not from benchmark performance.
- Choose the sampling height from the simulated/experimental STM setup, not from
  sequence accuracy.
- The beta-(1->4) linkage constrains orientation and parity only; it must not
  impose a GlcNAc count or an alternating sequence.
- Keep `010010` and any unit truth outside mold construction and scoring.
