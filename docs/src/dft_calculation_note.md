# DFT-STM Calculation Note for GlcN and GlcNAc

## Purpose and scientific scope

These Quantum ESPRESSO calculations are intended to generate local
Tersoff--Hamann-like occupied-state LDOS molds for distinguishing a central
glucosamine (GlcN) unit from a central N-acetylglucosamine (GlcNAc) unit in an
adsorbed chitosan chain on Cu(100). They are not adsorption-energy calculations
and do not constitute a complete simulation of the experimental tip or tunneling
current.

The comparison uses the same chain context for both central units:

```text
GlcN case:   GlcN--GlcN--GlcN
GlcNAc case: GlcN--GlcNAc--GlcN
```

Thus, the chemical variable is the acetylation state of the central unit. No
experimental unit sequence, expected composition, or benchmark label enters the
DFT setup.

## Atomic models

| quantity | central GlcN | central GlcNAc |
|---|---:|---:|
| isolated trimer formula | `C18 H35 N3 O13` | `C20 H37 N3 O14` |
| trimer atoms | 69 | 74 |
| Cu atoms in active pilot slab | 144 | 144 |
| total atoms in QE calculation | 213 | 218 |
| central acetylated units | 0 | 1 |
| acetylated neighbor units | 0 | 0 |

The active resource-constrained model uses a three-layer `8 x 6` Cu(100) slab.
The lower 96 Cu atoms are frozen during ionic relaxation and the upper Cu layer
plus the molecular trimer are relaxed. Both molecules use the same lateral
Cu(100) registry and central-ring frame.

The orthorhombic cells are:

| system | `a` (A) | `b` (A) | `c` (A) |
|---|---:|---:|---:|
| GlcN | 20.449528 | 15.337146 | 21.486251 |
| GlcNAc | 20.449528 | 15.337146 | 21.811770 |

The slab top is at approximately `z = 3.615 A`. The pilot geometry was built
with about 12 A of vacuum. The different `c` values reflect the different
molecular extents, not a different in-plane substrate model.

## Electronic-structure parameters

The two systems use the same electronic-structure settings:

- Quantum ESPRESSO `pw.x` / `pp.x`, Raven module `qe/7.4.1`;
- PBE exchange-correlation functional;
- Grimme D3 dispersion correction;
- scalar, non-spin-polarized calculation (`nspin` is not set);
- PAW pseudopotentials from the PSLibrary/KJPAW family;
- Cu `pbe-dn-kjpaw` potential with 11 valence electrons;
- wavefunction cutoff `ecutwfc = 50 Ry`;
- charge-density cutoff `ecutrho = 360 Ry`;
- Gamma-point-only Brillouin-zone sampling;
- Marzari--Vanderbilt smearing, `degauss = 0.02 Ry`;
- 8 MPI ranks, one OpenMP thread per rank;
- 96 GB memory and a 24-hour walltime allocation on Raven.

This is the reduced pilot model adopted after the earlier `8 x 8 x 4`,
`80/640 Ry`, `2 x 2 x 1` setup exceeded the available memory. The pilot is
appropriate for producing first local LDOS molds, but substrate-size,
k-point, and cutoff convergence remain future scientific checks.

## Geometry relaxation and final SCF

Ionic positions are optimized with BFGS. The current resource-aware
"relaxed-enough" thresholds are:

```text
forc_conv_thr = 6.0e-3 Ry/Bohr
etot_conv_thr = 1.0e-3 Ry
```

The relaxation SCF threshold is `1.0e-6 Ry`; the final fixed-geometry SCF target
is stricter, `1.0e-7 Ry`. These thresholds were chosen from observed force
plateaus and available walltime, not from agreement with unit-assignment labels.
The same ionic policy is applied to GlcN and GlcNAc.

### GlcN result

The GlcN geometry accumulated BFGS progress over several 24-hour restarts.
Restart5 began from the final restart4 geometry and already satisfied the
relaxed-enough criterion, requiring no additional BFGS step. Raven job
`28784933` completed successfully:

- relaxation-stage SCF converged in 34 iterations;
- final fixed-geometry SCF converged in 37 iterations;
- `pw.x` and `pp.x` ended normally;
- `glcn_central_relaxed.xyz` was written;
- `glcn_central_ldos.cube` was produced (approximately 136 MB).

This is the first completed production GlcN LDOS calculation in the workflow.

### GlcNAc result and current retry

Raven job `28790944` completed the ionic relaxation and wrote
`glcnac_central_relaxed.xyz`. Its final SCF, however, did not reach the
`1.0e-7 Ry` target within the default 100 electronic iterations. Consequently,
QE wrote configuration-only restart metadata rather than collected final
wavefunctions, and `pp.x` could not produce an LDOS cube.

The relaxation is not being repeated. SCF+PP-only retry `28811340` used the
relaxed geometry and saved charge density with:

```text
electron_maxstep = 300
mixing_mode = local-TF
mixing_beta = 0.1
mixing_ndim = 16
startingpot = file
startingwfc = atomic+random
```

This retry also stopped without electronic convergence after 300 iterations.
Its best and final estimated SCF accuracies were `3.461e-5 Ry` and
`4.056e-5 Ry`, respectively, so the convergence guard correctly withheld
`pp.x`. Memory use was about 30.7 GB, well below the 96 GB request.

Numerical retry2 `28851882` keeps the relaxed geometry, Hamiltonian, cutoffs,
Gamma sampling, smearing, and strict `1.0e-7 Ry` target unchanged. It tests the
plain mixing path that converged the GlcN control, with a longer Pulay history:

```text
electron_maxstep = 300
mixing_mode = plain
mixing_beta = 0.3
mixing_ndim = 20
startingpot = file
startingwfc = atomic+random
```

Retry2 also stopped after 300 iterations, with minimum/final estimated accuracy
`3.239e-5 / 3.897e-5 Ry`. The explicit convergence guard again prevented
`pp.x`, and no cube was written. Because two distinct mixing paths reproduce the
same approximate `3e-5--4e-5 Ry` floor, no further mixing-only retry is currently
accepted as a production calculation.

The documented common acceptance criterion is `5e-5 Ry`, below `0.7 meV` for
the complete cell; the electronic smearing remains `0.02 Ry` and the STM
sample bias is `-0.300 V`. Account-wide Raven limits blocked several acceptance
submissions before start, so the exact 3.30 GB retry2 density, distributed
wavefunctions, and mixing checkpoint were transferred to Viper. With the same
QE 7.4.1 stack and eight-MPI decomposition, plain job `10640236` and TF job
`10640237` each converged in one iteration at `4.027e-5 Ry`; both had energy
`-31763.44437132 Ry`, Fermi energy `0.3038 eV`, and successful `pp.x` output.
Comparison job `10640238` passed: the complete cubes are byte-identical, their
0.50 nm normalized maps have zero pointwise difference and correlation 1, and
their energy difference is zero. The plain cube is the canonical GlcNAc result;
the TF cube is retained as independent numerical-path validation. GlcN already
satisfies the same acceptance criterion with its stricter `6e-8 Ry` result.

## LDOS quantity exported for STM molds

`pp.x` uses `plot_num = 5` for the Tersoff-Hamann STM quantity. For this plot
mode QE 7.4.1 ignores `emin`, `emax`, and `degauss_ldos`; inspection of
`PP/src/stm.f90` confirms that it integrates states from the Fermi level to
`E_F + sample_bias`. The production input is therefore

```text
sample_bias = -0.0220495933 Ry
```

QE prints `Sample bias = -0.3000 eV`, matching the bias recorded in the 240817
STM scans. PP-only Viper jobs `10640445` and `10640446` regenerated and compared
the corrected cubes. The GlcNAc plain/TF cubes are byte-identical (SHA-256
`40649ccd9eb6768444b8ff61bf4a639b3940cf926fe3eb42254eb024faf9b5bf`); the
corrected GlcN cube has SHA-256
`80cd1d1fde94cf084cc7ea464d2bf065b36b2c015ea0bfaef8cebfee8ff88863`.
The output is a three-dimensional Gaussian cube. STMFit later
samples constant-height planes above the central ring; the nominal physical
height is `0.50 nm`, with `0.40--0.60 nm` retained as a sensitivity bracket.

The resulting maps are local LDOS proxies. They do not include an explicit tip
orbital, tip relaxation, finite-temperature transport, solvent, or a calibrated
constant-current feedback loop.

## Diagnostic constant-current observable

The accepted GlcN and GlcNAc cubes can also be transformed into a
constant-current-like **diagnostic** with
`test/build_constant_current_stm_maps.jl`. For each lateral column, the
transform searches from vacuum toward the slab and linearly interpolates the
first bracketed isovalue crossing. The isovalue is fixed by the same `0.50 nm`
mean-height-above-Cu policy used to anchor the existing mold workflow. The
predeclared `0.40:0.05:0.60 nm` heights remain sensitivity diagnostics; none may
be selected from benchmark performance.

Every map is paired with a validity mask, and the provenance sidecar binds the
two accepted cube hashes, `-0.300 eV` bias, distinct per-type Cu reference
frames, nominal and bracket heights, typed per-type nominal isovalues, z
spacing, crossing convention, and map/mask hashes.
The transformation leaves the common `5e-5 Ry` electronic acceptance criterion
unchanged. It also does not turn the QE quantity into an absolute tunneling
current: `plot_num=5` remains a discrete sum of `|psi_n(r)|^2` over the bias
window, with no justified conversion to nA.

The transformed maps use the temporary provider identity `stm_dft_cc_diag` and
the separate config `config/joint_proxy_whole_roi_constant_current.toml`. They
are not promoted into `config/joint_proxy_molds.toml`; the active
constant-height `stm_dft_v1` artifacts, registry payload, and pinned hashes stay
unchanged. Exact synthetic cases and fixed noise, drift, blur, and frame
perturbations pass at all five heights, while the common-only control abstains.
That result is a software/physics-fixture gate, not evidence that chemical unit
identity transfers to experimental constant-current scans.

## Preliminary versus production molds

Preliminary GlcN/GlcNAc mold TSVs at heights `0.40`, `0.45`, `0.50`, `0.55`,
and `0.60 nm` remain available for diagnostics but are disabled in
`config/joint_proxy_molds.toml`.
They were generated from unconverged or unrelaxed geometries to validate the
conversion and inference pipeline. They must not be described as converged DFT
identity references.

The production replacement gate is:

1. obtain successful final SCF and LDOS cubes for both relaxed systems;
2. sample each cube in its own relaxed-geometry frame, with the same
   normalization, grid, and physical height policy;
3. import them as a new versioned mold provider;
4. regenerate synthetic calibration so provenance hashes bind the new maps;
5. perform unlabeled sensitivity and visual QC before any external grading.

Both relaxed systems have crossed the cube-generation and numerical-stability
gate. The accepted cubes now feed the versioned 9×9 `stm_dft_v1` provider at
`0.50 nm`. Its tracked map, unary-template TSV, and provenance sidecar are bound
to pinned cube/map/template hashes by the registry loader. The five preliminary
height sources were not relabeled: they remain listed with `enabled=false`.
Synthetic paired A/B checks found no recovery or abstention regression relative
to geometric-only across three fixed seeds; this is a software/provenance gate,
not evidence of real chemical classification accuracy.

## Reproducibility map

| item | location |
|---|---|
| validated starting structures | `hpc/qe_molds/structure_validation.tsv` |
| pilot cell metadata | `hpc/qe_molds/*_slab_pilot_meta.tsv` |
| prepared QE inputs and outputs | `qe/glcn_restart5/`, `qe/glcnac/` |
| submitted job ledger | `hpc/qe_molds/qe_jobs.tsv` |
| operational history | `hpc/qe_molds/README.md` |
| mold conversion workflow | `docs/src/qe_stm_molds.md` |
| inference provider configuration | `config/joint_proxy_molds.toml` |
| production synthetic calibration | `config/joint_proxy_calibration_dft_m030_h050_v1.toml` |

Large QE outputs and cube files under `qe/` are intentionally not tracked by
Git. The scientific state should therefore be reconstructed from this note,
the tracked inputs/job ledger, and the archived remote/local QE outputs.
