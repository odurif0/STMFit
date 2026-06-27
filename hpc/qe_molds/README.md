# Quantum ESPRESSO DFT-STM molds for chitosan

This folder contains a minimal Quantum ESPRESSO workflow for producing STM/LDOS
molds usable by STMFit unit assignment. The target output is not a molecular
contour; it is a simulated STM/LDOS map exported to:

```text
templates/chitosan_stm_maps.tsv
```

## Scientific target

Start with two beta-(1->4)-linked trimers:

```text
GlcN   central unit: X-GlcN-X
GlcNAc central unit: X-GlcNAc-X
```

Use the same neighbor choice `X` in both first calculations. If this improves the
patch-level score, extend to all eight left/central/right contexts.

The DFT maps should represent the central unit in the adsorbed chain context on
Cu(100), at the experimental STM bias.

## Minimal QE sequence

1. Build Cu(100) slab + chitosan trimer geometry as XYZ.
2. Generate QE inputs with `test/prepare_qe_mold_inputs.jl`.
3. Relax with `pw.x`.
4. Run a converged SCF calculation.
5. Export an integrated LDOS/STM-like cube with `pp.x`.
6. Convert the cube to STMFit's map TSV with `test/cube_to_stm_maps.jl`.
7. Import the map TSV into connected templates with `test/import_stm_mold_maps.jl`.
8. Score with `test/score_connected_mold_templates.jl`.

`glcn_central_trimer_slab.xyz.template` and
`glcnac_central_trimer_slab.xyz.template` document the two XYZ structures to
provide. They are placeholders, not scientific geometries.

Initial, unoptimized trimer geometries can be regenerated from the repository
root with:

```bash
julia --project=. test/build_initial_chitosan_trimer_xyz.jl \
    --out-dir hpc/qe_molds
```

This writes `glcn_central_trimer.xyz` and `glcnac_central_trimer.xyz`, plus atom
tables and frame-index TSVs. The generated structures use GlcN neighbors in both
calculations and must be relaxed before production STM/LDOS maps are extracted.

Validate the generated structures with:

```bash
julia --project=. test/validate_chitosan_trimer_structures.jl \
    --dir hpc/qe_molds \
    --out hpc/qe_molds/structure_validation.tsv
```

Current validation summary:

| structure | formula | central acetyl units | neighbor acetyl units | bare atoms | slab atoms |
|---|---:|---:|---:|---:|---:|
| `glcn` | `C18 H35 N3 O13` | 0 | 0 | 69 | 325 |
| `glcnac` | `C20 H37 N3 O14` | 1 | 0 | 74 | 330 |

The validated default-frame indices after the `8×8×4` slab are:

```text
--origin-indices 268,269,270,271,272,273 --axis-from 271 --axis-to 268 --plane-index 277
```

## Visualization

Static previews can be regenerated with:

```bash
python3 test/render_xyz_preview.py hpc/qe_molds/glcn_central_trimer.xyz \
    --out hpc/qe_molds/previews/glcn_central_trimer.png \
    --title "GlcN central trimer: GlcN-GlcN-GlcN" \
    --molecule-only
python3 test/render_xyz_preview.py hpc/qe_molds/glcnac_central_trimer.xyz \
    --out hpc/qe_molds/previews/glcnac_central_trimer.png \
    --title "GlcNAc central trimer: GlcN-GlcNAc-GlcN" \
    --molecule-only
```

For interactive inspection, Avogadro is available locally:

```bash
avogadro2 hpc/qe_molds/glcn_central_trimer.xyz
avogadro2 hpc/qe_molds/glcnac_central_trimer_slab.xyz
```

## MPCDF module sketch

For the 240817 benchmark scans, the SXM headers report a uniform bias of
`-0.300 V` and setpoint `2.0 pA`. The prepared local run directories use
`emin=-0.3 eV`, `emax=0.0 eV` for the first occupied-state LDOS export:

```text
qe/glcn
qe/glcnac
```

These `qe/` run directories are intentionally gitignored because real QE runs
produce large scratch and cube files.

Before submitting, run the input preflight:

```bash
julia --project=. test/preflight_qe_mold_inputs.jl \
    --dir qe/glcn \
    --dir qe/glcnac \
    --out hpc/qe_molds/qe_input_preflight.tsv \
    --max-total-tasks 8 \
    --sequential \
    --min-mem-mb 96000
```

On the HPC system, submit both QE jobs with:

```bash
bash hpc/submit_qe_molds.sh --watch --sequential
```

Use `--dry-run` first if you only want validation and the exact `sbatch`
commands. By default `--sequential` is on, so the two QE run directories
(`qe/glcn`, `qe/glcnac`) are submitted as an `afterok` chain and the task budget
is enforced as the *maximum simultaneous* task count, not the sum; this keeps a
single 8-task job within the Raven `n0001` group limit. It loads the configured
Julia module for its preflight step if `julia` is not already on `PATH`, and the
generated compute-job scripts load the same Julia module for the relax-to-SCF
handoff. Generated Slurm scripts set `--ntasks-per-node=8 --cpus-per-task=1` and
launch QE with `srun -n "$QE_NTASKS" --cpu-bind=cores`; memory is
`ntasks × mem-per-task-mb` (`96000 MB` for the validated `8 × 12000 MB` pilot),
because `--mem=0` is rejected on shared Raven nodes.

On Raven, the tested QE module stack is:

```bash
module load intel/2024.0 impi/2021.11 qe/7.4.1
```

Override this with `QE_COMPILER_MODULE`, `QE_MPI_MODULE`, and `QE_MODULE` in
`hpc/remote.env` if the cluster module tree changes. The prepared runs include
local PSLibrary/KJPAW pseudos in `qe/<run>/pseudo/*.UPF`; the remote launcher
syncs these files explicitly, and `test/preflight_qe_mold_inputs.jl` verifies
that every `ATOMIC_SPECIES` pseudo file exists before submission.

The original full pilot (`8×8×4` slab, Cu `spn` PAW, `2×2×1` k-points,
`ecutwfc=80`, `ecutrho=640`) exceeded the Raven QOS memory limit: QE estimated
`~159 GB` dynamic RAM per MPI process. The active Raven pilot therefore uses an
`8×6×3` slab with `12 Å` vacuum, the lighter Cu `dn` PAW pseudo
(`Cu.pbe-dn-kjpaw_psl.1.0.0.UPF`, `z_valence=11`), `K_POINTS gamma`
(`--kpoints 1,1,1`), `ecutwfc=50`, `ecutrho=360`, `8` MPI tasks, and
`96000MB` per job.
This is a resource-constrained first DFT-STM mold calculation; record any later
convergence upgrades separately from benchmark grading.

The Raven `n0001` QOS can also enforce a one-node group limit, and pending
dependent jobs can still count against it. Submit GlcN alone first with
`--dir qe/glcn`; submit GlcNAc only after GlcN succeeds. Use the sequential or
parallel multi-dir modes only when your QOS allows the extra pending allocation.
The 2026-06-25 GlcN pilot job `28288215` showed why the clean fallback is needed:
after an in-place reduction from 4 to 2 MPI tasks, it started but failed with
Slurm `OUT_OF_MEMORY` at a 16 GB node limit; QE estimated `34.30 GB` total
dynamic RAM. The clean 2-task/48 GB replacement `28292263` was useful as a
memory probe, but was too slow: after ~7 h it had completed only the initial SCF
(`bfgs steps = 0`) and was still in the next SCF. Regenerate the run directories
as 8-task, 96 GB, 24 h jobs before resubmission:

```bash
julia --project=. test/prepare_qe_mold_inputs.jl \
  --xyz hpc/qe_molds/glcn_central_trimer_slab_pilot.xyz \
  --cell-metadata hpc/qe_molds/glcn_central_trimer_slab_pilot_meta.tsv \
  --out-dir qe/glcn --prefix glcn_central \
  --fix-below-z 1.807501 --emin-ev -0.3 --emax-ev 0.0 \
  --ecutwfc 50 --ecutrho 360 --kpoints 1,1,1 \
  --ntasks 8 --mem-per-task-mb 12000 --walltime 24:00:00
```

Regenerate `qe/glcnac` with the analogous `glcnac_central` pilot XYZ/metadata and
submit it only after GlcN completes successfully. Generated sbatch files now set
`#SBATCH --cpus-per-task=1`, derive `QE_NTASKS=${SLURM_NTASKS:-8}`, and call
`srun -n "$QE_NTASKS" --cpu-bind=cores ...`; verify the QE header says
`Number of MPI processes: 8` after launch.
Set `QE_MIN_MEM_MB=96000` in `hpc/remote.env` or pass `--min-mem-mb 96000` to
make remote preflight reject stale 48 GB sbatch files before submission.
The first optimized GlcN-only resubmission is job `28303162`
(`8` CPUs, `96000M`, `24:00:00`).

### Restarting the relax from the best-so-far geometry (preserves BFGS progress)

The generated sbatch wipes `qe_tmp` and restarts from the original geometry, so a
naive resubmission after a walltime timeout loses all BFGS steps. To preserve the
relaxation progress, restart from the **latest geometry in the relax output**:

```bash
# 1. Fetch the (timed-out or running) relax output
rsync -avz -e "ssh ..." \
  raven:/u/oldu/code/STMFit/qe/glcn/glcn_central_relax.out qe/glcn/

# 2. Extract the last ATOMIC_POSITIONS block as XYZ
julia --project=. test/extract_qe_relaxed_xyz.jl \
    --qe-out qe/glcn/glcn_central_relax.out \
    --out qe/glcn/glcn_central_best.xyz \
    --metadata qe/glcn/glcn_central_best_meta.tsv

# 3. Regenerate a fresh relax dir starting from that geometry
julia --project=. test/prepare_qe_mold_inputs.jl \
    --xyz qe/glcn/glcn_central_best.xyz \
    --cell-metadata hpc/qe_molds/glcn_central_trimer_slab_pilot_meta.tsv \
    --out-dir qe/glcn_restart --prefix glcn_central \
    --fix-below-z 1.807501 --emin-ev -0.3 --emax-ev 0.0 \
    --ecutwfc 50 --ecutrho 360 --kpoints 1,1,1 \
    --ntasks 8 --mem-per-task-mb 12000 --walltime 24:00:00
cp qe/glcn/pseudo/*.UPF qe/glcn_restart/pseudo/

# 4. Submit (single dir, sequential irrelevant)
bash hpc/launch_qe_molds_remote.sh --dir qe/glcn_restart --max-total-tasks 8 --min-mem-mb 96000
```

This is a fresh relax (not a QE `restart_mode`), but starting from the
partially-converged geometry — far fewer ionic steps than from the hand-built
structure.

### Preliminary mold from an unconverged geometry (pipeline de-risking)

While a relax is still running (or restarting), a preliminary LDOS cube can be
produced from the best-so-far geometry to validate the downstream
finalize/score/grade pipeline end-to-end. This is an **approximate** mold, not a
production map:

```bash
# Build a SCF+PP-only dir from the best geometry (no relax step)
julia --project=. test/update_qe_positions_from_xyz.jl \
    --input qe/glcn/pw_scf.in --xyz qe/glcn/glcn_central_best.xyz \
    --out qe/glcn_prelim/pw_scf.in
cp qe/glcn/pp_ldos.in qe/glcn_prelim/ && cp qe/glcn/pseudo/*.UPF qe/glcn_prelim/pseudo/
# run_scf_pp.sbatch = run_qe_mold.sbatch with the relax/extract/update steps removed
bash hpc/launch_qe_molds_remote.sh --dir qe/glcn_prelim --max-total-tasks 8 --min-mem-mb 48000
```

The preliminary cube lives in `qe/glcn_prelim/glcn_central_ldos.cube`. Use it only
to exercise `cube_to_stm_maps.jl` → `import_stm_mold_maps.jl` → scoring; do not
treat its GlcN mold as final. The first such preliminary run is job `28354566`.

From the local workstation, after configuring `hpc/remote.env`, you can sync and
submit in one step:

```bash
bash hpc/launch_qe_molds_remote.sh --dry-run
bash hpc/launch_qe_molds_remote.sh --watch
```

The remote launcher syncs the repository while excluding local `qe/` outputs,
then separately syncs only the prepared QE input files and `pseudo/*.UPF` files
for `qe/glcn` and `qe/glcnac`.
It uses `SSH_CONNECT_TIMEOUT=180` by default; increase this in `hpc/remote.env`
if the MPCDF gateway or password/OTP flow needs more time.

After both QE jobs finish and the cubes are present, convert them to STMFit molds
with an explicitly chosen sampling height:

```bash
julia --project=. test/finalize_qe_mold_workflow.jl \
    --height-nm HEIGHT \
    --glcn-dir qe/glcn \
    --glcnac-dir qe/glcnac
```

This extracts `frame.tsv` files from the relaxed geometries, samples the LDOS
cubes with typed frames, writes `templates/chitosan_stm_maps.tsv`, and imports
connected unary/bond mold templates.
`finalize_qe_mold_workflow.jl` now derives the relaxed-slab frame indices from
`hpc/qe_molds/*_trimer_indices.tsv` plus the atom-count offset of the relaxed
XYZ. This avoids hardcoded `8×8×4` indices and works for the current `8×6×3`
pilot slabs.

Check exact module names on Raven/Viper:

```bash
find-module qe
module load intel/2024.0 impi/2021.11 qe/7.4.1
```

Then generate inputs from a vetted XYZ structure:

If you have only an oriented trimer XYZ, first add a Cu(100) slab:

```bash
julia --project=. test/build_qe_slab_trimer_xyz.jl \
    --molecule hpc/qe_molds/glcn_central_trimer.xyz \
    --out hpc/qe_molds/glcn_central_trimer_slab.xyz \
    --metadata hpc/qe_molds/glcn_central_trimer_slab_meta.tsv \
    --nx 8 --ny 8 --layers 4 \
    --center-indices 12,13,14,15,16,17 \
    --height-above-top 2.6 --vacuum 18.0
```

The central-ring indices keep GlcN and GlcNAc on the same Cu(100) registry. The
generated `*_indices.tsv` files also provide the post-slab frame command for a
default `8×8×4` slab.

Then write QE inputs:

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

Repeat for `glcnac_central_trimer_slab.xyz`. Adapt the generated
`run_qe_mold.sbatch` if the MPCDF module name differs.

The generated Slurm script extracts relaxed coordinates from `pw.x` output and
updates `pw_scf.in` automatically before the production SCF. If the run directory
is not under the STMFit checkout, set:

```bash
export STMFIT_ROOT=/path/to/STMFit
```

After relaxation, extract the local frame for cube sampling:

```bash
julia --project=. test/extract_qe_mold_frame.jl \
    --xyz qe/glcn/glcn_central_relaxed.xyz \
    --origin-indices I,J \
    --axis-from I --axis-to J --plane-index K \
    --height-nm HEIGHT \
    --out qe/glcn/frame.tsv
```

## Important constraints

- Do not choose the DFT height, bias window, or geometry because it improves the
  benchmark sequence. Choose them from the experimental bias/setpoint and DFT
  convergence/physics.
- The beta-(1->4) linkage constrains orientation and parity only. It must not
  impose an alternating sequence or a fixed GlcNAc count.
- Freeze the imported mold before running `grade_unit_assignment.jl`.

## Output conversion example

After `pp.x` produced two cube files, one for GlcN and one for GlcNAc:

```bash
julia --project=. test/cube_to_stm_maps.jl \
    --cube 0:qe/glcn/glcn_central_ldos.cube \
    --frame 0:qe/glcn/frame.tsv \
    --cube 1:qe/glcnac/glcnac_central_ldos.cube \
    --frame 1:qe/glcnac/frame.tsv \
    --cube-units bohr \
    --out templates/chitosan_stm_maps.tsv
```

The frame files contain `origin_nm`, `t_axis`, `u_axis`, and `height_nm` in the
cube coordinate system. `t` must point along the molecular backbone. `u` should
point toward the C2 substituent side for the reference orientation. `height_nm`
samples a constant-height plane above the central unit; set it from the intended
Tersoff-Hamann/STM map.

Then import and score:

```bash
julia --project=. test/import_stm_mold_maps.jl \
    --maps templates/chitosan_stm_maps.tsv \
    --out templates/chitosan_connected_molds_stm.tsv \
    --bond-out templates/chitosan_connected_bond_molds_stm.tsv \
    --half-nm 0.48 --step-nm 0.08

julia --project=. test/score_connected_mold_templates.jl \
    --patches results/unit_separability/lobe_patches_selectedN_primary_half048.tsv \
    --templates templates/chitosan_connected_molds_stm.tsv \
    --bond-templates templates/chitosan_connected_bond_molds_stm.tsv \
    --prefix raw_p \
    --template-mode contrast \
    --out results/unit_assignment/stm_mold_predictions.tsv
```
