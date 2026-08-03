# Quantum ESPRESSO STM Molds

This page describes how to produce STM/LDOS molds for GlcN/GlcNAc unit
assignment using Quantum ESPRESSO on an HPC system such as MPCDF Raven/Viper.

For the scientist-facing description of the atomic systems, approximations,
parameters, quantities calculated, completed jobs, and interpretation limits,
see [DFT-STM Calculation Note for GlcN and GlcNAc](dft_calculation_note.md).

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
`-0.300 V` over the folder. QE `pp.x` with `plot_num=5` uses `sample_bias` in
Ry, not `emin`/`emax`; the prepared production inputs therefore use
`sample_bias=-0.0220495933 Ry`, which QE reports as `-0.3000 eV`. These run
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
molds. The constant-current smoke leg uses a dedicated strictly monotone
`value = 3z` cube so its target has one continuous fixed-support root branch;
the constant-height leg retains its original synthetic cube. This fixture
replacement preserves all smoke assertions. It does not validate chemistry or
run Quantum ESPRESSO.

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
    --sample-bias-ev -0.3
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

## Diagnostic constant-current transformation

`test/build_constant_current_stm_maps.jl` provides a separate, diagnostic-only
transformation of the accepted LDOS cubes. It finds the first bracketed crossing
when sampling from vacuum toward the slab and writes both the resulting height
map and an explicit validity mask. Missing, ambiguous, non-finite, or
insufficiently sampled crossings remain invalid; they are never replaced by
zero or extrapolated heights.

```bash
julia --project=. test/build_constant_current_stm_maps.jl \
    --cube0 qe/glcn/glcn_central_ldos.cube \
    --frame0 qe/glcn/frame.tsv \
    --cube1 qe/glcnac/glcnac_central_ldos.cube \
    --frame1 qe/glcnac/frame.tsv \
    --cube-units bohr \
    --half-nm 0.80 \
    --step-nm 0.08 \
    --isovalue-scan-intervals 1024 \
    --out-prefix /tmp/opencode/chitosan_cc_diag
```

The nominal `0.50 nm` mean-height policy and fixed
`0.40, 0.45, 0.50, 0.55, 0.60 nm` bracket are declared before grading. Each
height produces a typed map and mask; the sidecar binds the observable, both
cube hashes, sample bias, per-type Cu reference frames, nominal and bracket
heights, typed nominal isovalues, cube-normal spacing, crossing policy, nominal
map/mask hashes, and the remaining bracket artifact paths and hashes. The
`type_frames` and `type_isovalues` provenance tables each contain exactly one
entry for type `0` and one for type `1`, so relaxed GlcN and GlcNAc cubes may
keep their distinct local frames and separately calibrated nominal isovalues
instead of copying, averaging, or selecting one value for both geometries.
Historical common-frame or single-explicit-isovalue synthetic calls are still
normalized into two typed entries. The constant-current map and validity mask
are published as a recoverable gate-last transaction: same-directory staged
files and a durable `prepared` marker precede destination changes, the map gate
is backed up first and installed last, and a durable `committed` marker records
the completed generation. On rerun, prepared recovery restores the old complete
set; committed recovery retains the new complete set and removes transaction
residue. Destination symlinks are replaced as directory entries rather than
followed. The same bounded protocol publishes the frozen diagnostic summary,
controls, and PNG set, with its summary TSV as the gate. It covers the tested
writer failures and explicit interruption points, but does not claim recovery
from filesystem corruption or hostile concurrent mutation of sidecars.
Provenance creation fails unless every declared non-nominal bracket
height has a bound map and invalid-mask artifact. The two cubes must also
resolve to the same normal-axis sampling spacing. Periodic cube wrapping is
rejected in this provenance-bound builder because that transform is not part of
the sidecar contract; it remains available only in the lower-level T1 diagnostic
sampler. The real whole-ROI diagnostic imports the nominal map over
`--half-nm 0.80 --step-nm 0.08`, producing the 21×21 (`441` pixel) mold support
expected by the whole-ROI consumer; the builder invocation for that chain must
use the same explicit support. Whole-ROI real diagnostic validation rejects a constant-current
sidecar unless both typed frame and typed isovalue entries are present, the
height/grid/units/z-spacing fields and `isovalue_scan_intervals` are well formed,
the nominal map extents and
step match those fields, the bracket artifact hashes match their paths, and the
connected mold TSV consumed by scoring has its own binding sidecar tying it back
to the nominal map hash and importer convention. The matching registration
config is
`config/joint_proxy_whole_roi_constant_current.toml`; its registration ranges
and numerical-tie tolerances are copied unchanged from the frozen whole-ROI
diagnostic config.

For each typed cube, mean-height calibration uses the explicit validated
`--isovalue-scan-intervals` policy, default `1024`. The sidecar records that value
as `isovalue_scan_intervals`, and whole-ROI validation requires it to match. The
default therefore evaluates 1025 equally spaced isovalues including both
endpoints of the declared interval. It partitions the valid response into
continuous fixed-column-support segments and accepts the target only when exactly
one segment contains exactly one root, which is then refined by bisection.
The public calibration controls are validated before any target conversion,
frame work, scan construction, or bisection: `max_iter` must be a positive
non-`Bool` integer and `height_tol_factor` must be a positive finite non-`Bool`
real. Valid custom values feed the existing bisection tolerance and iteration
limit; the defaults and root policy are unchanged.
An exact root at either endpoint of the declared scan range is rejected because
two-sided continuity and a crossing within the declared interval cannot be
established. An exact interior root requires both immediate adjacent scan
responses to be finite and their valid-column support to remain unchanged;
an isolated exact sample, including one with only one finite neighbour, is
rejected as continuity-ambiguous. Multiple target roots and any target crossing
at a support discontinuity are also explicit ambiguity errors; the code does
not choose the first, lowest, highest, or best-covered branch. A target outside
the reachable range or an interval with no valid crossings also fails
explicitly. Root uniqueness is assessed at the declared scan resolution, not
claimed independently of resolution.

This stronger rejection contract does not make the accepted real cube gate pass.
The accepted GlcNAc response remains multibranch, so T3 remains terminal
`BLOCKED` without a separately predeclared and validated physical branch policy.
No support threshold or branch preference was retuned from that failure.

This transformation does **not** calibrate current. QE `plot_num=5` is a
discrete occupied-state `|psi_n(r)|^2` sum over the bias window, not amperes, so
no nA-to-cube conversion is claimed. The temporary provider name is
`stm_dft_cc_diag`; it is not registered as `stm_dft_v1`, does not modify
`config/joint_proxy_molds.toml`, and cannot enter fitting, `N_selected`,
calibration, thresholds, or production abstention. Its synthetic gate validates
the transformation and scoring mechanics only, not GlcN/GlcNAc transfer on real
STM scans.

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

For the diagnostic constant-current chain, import the nominal map with its
constant-current sidecar so the connected mold TSV gets a child binding sidecar:

```bash
julia --project=. test/import_stm_mold_maps.jl \
    --maps /tmp/opencode/chitosan_cc_diag_h050.tsv \
    --out /tmp/opencode/chitosan_cc_diag_h050_connected.tsv \
    --half-nm 0.80 \
    --step-nm 0.08 \
    --source-provenance /tmp/opencode/chitosan_cc_diag.provenance.toml
```

This writes
`/tmp/opencode/chitosan_cc_diag_h050_connected.tsv.provenance.toml` by default.
Whole-ROI constant-current validation receives the actual `--molds` path and
checks this child sidecar before loading the TSV for scoring. A stale mold,
stale source-map hash, malformed numeric provenance, or missing bracket artifact
is rejected before any scan is scored.

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

- Choose the sample bias from the STM experiment, not from benchmark performance.
- Choose the sampling height from the simulated/experimental STM setup, not from
  sequence accuracy.
- The beta-(1->4) linkage constrains orientation and parity only; it must not
  impose a GlcNAc count or an alternating sequence.
- Keep `010010` and any unit truth outside mold construction and scoring.

## Joint posterior provider contract

`config/joint_proxy_molds.toml` exposes STM molds through a versioned registry
rather than embedding them in inference code. Each provider supplies unary
templates for both physical unit identities across parity and mirror states on
the configured 9×9 patch grid. The active `stm_dft_v1` source is required and
hash-pinned; geometric and DFT each carry weight `0.5`. Preliminary heights stay
listed but disabled and the tracked geometric provider remains reproducible.

The consumer writes five auditable artifacts: candidate-N posterior rows,
candidate-lobe geometry/type rows, selected-candidate `0/1/?` predictions, chain
summaries, and a run manifest containing config/source/payload SHA-256 hashes.
Neither provider nor consumer accepts an expected N, composition count, control
sequence, or grader input.

The production promotion performed the following steps:

1. extract both maps with the same physical height, field of view, orientation,
   normalization, parity, mirror, and grid convention;
2. register the generated TSVs and provenance sidecar as the non-preliminary
   `stm_dft_v1` source family without
   changing hard-label rules or production `N_selected`;
3. regenerate the synthetic calibration so provenance binds the new source and
   payload hashes;
4. run paired synthetic non-regression/QC first, retaining `?` calls; and
5. use any known sequence only afterward for external reporting, never to choose
   height, family weight, flip, threshold, or abstention.

The two-file 2026-07-11 smoke used the geometric/preliminary contract and
validated file/schema flow only. Both real chains abstained on hard count, so it
does not validate chemical identity or counting accuracy.
