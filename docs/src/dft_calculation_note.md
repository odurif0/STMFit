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
mode QE 7.4.1 ignores `emin`, `emax`, and `degauss_ldos`; `PP/src/stm.f90` uses
the window between the Fermi level and `E_F + sample_bias`, with smearing-weighted
states outside that window (see the September 24 qualification below).
The production input is therefore

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

**September 26 saved-export XC revalidation.** Separate analysis-only job
**11989810** completes in **21m13s** with **22,642** independent control/norm
and **70,062** component-cube checks, retaining all **26** byte identities and
the original **16-eps** norm bound. Source archive hashes remain unchanged;
new reports are fetched and checksum-verified. No QE or SCF is rerun. Native
cross-build maximum total-potential differences are **8.3e-9/8.7e-9 Ry**,
reported descriptively rather than accepted against a fitted tolerance.

The decomposition at the accepted density is qualified. Local-density XC has
the largest individual spatial SD on **93/122 GlcN and 88/115 GlcNAc**
PAW-free planes and every **61/57** lower-half plane; divergence is largest
on the remaining full-gap planes. SD rankings are not explained-variance
fractions. Native negative valence fractions reach **50.07/49.78%** within a
plane, so small global negative-charge integrals do not establish local-tail
accuracy. This is not density/cutoff convergence, a GGA-only explanation,
a justified potential substitution or a new mold/recognition result.

**September 26 paired-build XC follow-up (historically failed).** Job
**11988093** runs for **25m44s**, passing all **26** same-build byte identities
and internal native/cube/additivity checks for both accepted states. It then
ends **FAILED 1:0**: independent checking passes **22,289** assertions and
fails **313**, all on the descriptive relative-L2 norm reduction. Independent
component-cube checking is not reached. The Cartesian-view summation issue
is reproduced synthetically and corrected locally without changing the
**16-eps** criterion, any exported field or any physical parameter. The full
saved exports subsequently pass the separate revalidation above; its original
component tables and plots remain provisional. All outputs and both
executables are fetched; post-failure checks confirm the pinned original
density/wavefunction and copied checkpoint hashes are unchanged. No second
job, SCF, threshold change, mold replacement or champion promotion follows.

**September 26 XC component attempt (not qualified).** Job **11982255**
compiles the isolated PP helper and produces the first GlcN native components,
then stops at the required byte comparison with the site density control.
The first density difference is **1e-11 electrons/Bohr^3**; native/cube
total controls also differ between builds. These first witnesses are not
whole-volume accuracy bounds. No repeat, GlcNAc export or component analysis
is reached. All partial outputs and the executable are fetched/checksum-
verified; all original and copied checkpoint hashes remain unchanged. The
guard is not loosened and no retry occurs in this cycle. A paired stock and
instrumented executable from the same build is needed before attributing
the potential variation. The accepted states, SCF gate and molds are unchanged.

**September 26 local effective potential.** PP-only job **11981276** completes
in **4m07s** on complete copies of the accepted states. Total local potential
(`plot_num=1`) and electrostatic control (`plot_num=11`) are exported; the
total repeats byte for byte in native and cube formats. All **43.2 million**
native/cube comparisons lie within declared printing bounds, **63,784**
independent saved-output checks pass, and original/copy state hashes remain
unchanged. All new exports and reports are fetched and checksum-identical.

Across the **122/115** complete PAW-free native planes, GlcN/GlcNAc lateral
total-potential spans range from **1.2921--13.9877 / 1.2939--12.3596 eV**.
Minimum spatial SD is **0.2038/0.2021 eV**. Even where every sampled barrier
is positive, lateral span/mean barrier is at least **0.3283/0.3072**.
XC SD exceeds electrostatic SD on **235/237** planes. This identifies a
component to investigate, not a demonstrated density-tail artifact or
cutoff-convergence result. No flatness criterion, matching plane, wavefunction
continuation or potential replacement is adopted. In particular the
electrostatic control is not a substitute Hamiltonian. Both the common
**5e-5 Ry** electronic acceptance and production molds remain unchanged.

**September 25 off-projection reference diagnostic.** Job **11980577**
completes in **19m02s**, with **582,693** independent saved-output checks.
Geometry alone retains **474/768** common cell-grid sites outside both
molecular PAW projections. All full segments clear every periodically repeated
sphere (minimum **0.0008697363 nm**). All **8,532** outcomes remain; there are
no unresolved leaves, searches/direct queries repeat exactly, and all **19,300**
direct queries lie inside the grouped-series enclosures. Both first-site
witnesses are byte-identical between local and Viper.

Full Cu-gap, Cu-half-gap and old molecular-half outcomes are separate. Every
full-gap case has multiple roots, including ascent toward the periodic slab;
two GlcNAc third-isovalue sites also have multiple descending roots. Cu-half
paired valid supports are **424/436/453**. Their inter-state background
differences have spatial SD **64.51/66.74/72.66 pm**, larger than the original
mean molecular contrast near **21 pm**. Hypothetical reference choice can
reverse its sign. This is not a supported scalar Cu reference, current
calibration or chemical assignment. Atomic PAW projection clearance does not
exclude internal molecular spaces or electronic tails. A whole-footprint/
common-exterior follow-up retains 325 candidates but does not reduce the
reference ambiguity: Cu-half paired spatial SD remains 65.42--73.94 pm, with
both referenced contrast signs. No distance cutoff is selected from the
plotted contrasts. The accepted states,
**5e-5 Ry** criterion and production sources are unchanged.

**September 25 shared-image envelope.** Job **11980183** fits a chemistry-free
source-only B-spline envelope after the frozen masked image background, with
resolution selected by the complete conditional GCV trace. It processes all
146 scans / 292 views; **132,131** independent grouped checks pass. On 141
identified-registration views, foreground cross-view MSE improves **24.537%**
against source copy. This is an image-prediction result, not a DFT or chemical
recognition validation. On identified physical-probe support, only
**52.08--52.49%** of the injected GlcNAc-minus-GlcN energy survives at the three
fixed isovalues. Probe geometry/support is inherited, and the mask/resolution
is not reselected under perturbation. The finest envelope grid is chosen in
every view, without demonstrating an optimal physical scale. A consistent
image/DFT background reference and joint envelope/contrast model remain open;
no off-molecule reference, current calibration, accepted state, SCF criterion
or production mold is changed.

**September 25 height/envelope controls.** Job **11979008** keeps the fixed
surface-transfer cohort, patches, calibration and all three isovalues. A common
shape with just the fixed mean class-height contrast reproduces
**99.52--99.57%** of the original chemical hypotheses on identified-registration
views. With the same source-only affine nuisance and chemical geometry frozen
at the common minimum, the physical shapes give **5.995--6.044 pm** RMS,
worse than the chemistry-free local plane (**5.447 pm**). Independent saved
checks pass **2,584,957/2,584,957**. The near-single-class hypotheses after plane
profiling are not chemical composition or unit counts. These conditional,
reused-cohort controls expose the current observation/support limitation;
they do not prove that chemistry is absent or change the numerically verified
maps, accepted DFT states, SCF criterion, current conversion or production molds.
Foreground and a shared molecular envelope are needed before a new chemical
comparison, without selecting pixels or settings from benchmark labels.

**September 25 raw-image transfer.** Job **11978725** tests all six continuous
maps on source-only diagnostic patches, retaining their relative physical
heights with unit gain and one common offset per source scan. Separate target
row blocks calibrate registration and measure prediction. All 146 raw scans
remain accounted for: 141 scored, one without sufficient native patch support
and four incompatible biases. On 138 identified-registration views, the chemical
dictionary lowers MSE by **54.48--55.06%** relative to the common surface across
the three fixed isovalues. Its **44.80--45.36 pm** RMS is nevertheless far worse
than source copy (**6.63 pm**). The first eligible scan's overlay associates
class hypotheses with molecular brightness versus background and includes many
nonmolecular anchors. This visual observation is not a chemical annotation;
the gain cannot establish chemical transfer. An explicit background, shared
molecular envelope and height-only null are the next relevant controls. No map,
current conversion, calibration, production provider or champion is promoted.

**September 25 continuous representative surfaces.** After the sampled-profile
job **11977798** found three common open isovalue intervals in the
geometry-defined molecular half of the vacuum gap, job **11978192** tests all
three fixed representatives using the untruncated finite-basis Fourier signal.
All **1,734** roots are unique and descending with no unresolved region; six
17x17 maps are complete. Search and independent 3D query repetitions are exact.
The maximum relative density residual is **1.77867e-7** (**0.00001779%**),
compared with **2.089%** for the earlier linearly interpolated roots. The largest
height correction is **0.668041 pm**. All **68,493** saved native-profile values
and **1,734** new surface values lie inside their strict interval enclosures.

This resolves the numerical crossing question at three representative
isovalues, not throughout the earlier intervals or in the isolated-slab limit.
The full periodic gap still has no common sampled-profile interval. The
molecular-half boundary remains a geometry convention, and finite-basis,
periodic-image, tip and experimental-current uncertainties remain. No new QE
calculation, SCF criterion, state normalization, production mold or recognition
claim follows. The next question is a shared label-free experimental observation
model, not additional root precision.

**September 25 independent wavefunction reconstruction.** Job **11977502**
reads the accepted collected coefficients and evaluates their Gamma Fourier
series directly in Julia, independently of QE's density FFT and `local_dos`.
Every one of the **13,872** native vertices required by the six declared
planes agrees with the saved smooth density within its frozen printing bounds.
Serial/parallel repetitions are identical. Minimum PAW clearances for the
native vertices are **0.161923 / 0.103670 nm**, and all **1,734** direct plane
queries are strictly outside all spheres and positive. This is pointwise
verification in the audited domain, not a new whole-volume reconstruction.

The reciprocal-space smooth norms are **18.1758333630 / 20.0011288529**;
exact-cell cube integrals differ by only **2.45e-8 / 5.83e-8**. Earlier cube
integrals used the rounded six-decimal axes. Both conventions are retained;
no component is rescaled to 24/26. Relative-L2 interpolation errors at
0.40/0.50/0.60 nm are GlcN **0.260/0.311/0.521%**, GlcNAc
**0.402/0.495/0.526%**. These are plane-level norms, not per-pixel maxima.
At 0.60 nm the largest pointwise relative interpolation errors are
**9.06% / 4.70%**; the direct density, not interpolation, supplies the following
decay diagnostic.

The mean density decreases strongly for both molecules. At individual lateral
points GlcN decreases over both intervals, but GlcNAc rises at **1/289** points
over 0.40–0.50 nm and **30/289** over 0.50–0.60 nm, up to **11.51-fold** in the
latter interval. A correct positive density therefore does not by itself
establish a unique molecular-side constant-current surface. The cause of the
rising tails and a physically supported branch remain open; no column is
masked, made monotonic or removed by a new height cutoff. All **293,483,676
bytes** of this run are fetched and checksum-identical. Production molds,
SCF acceptance, current calibration and recognition remain unchanged.

**September 25 PAW component attribution.** Job **11976981** reproduces the
Gamma-corrected full ILDOS and all legacy controls byte for byte. Separate
smooth and augmentation exports from unchanged projector/band sums each
repeat exactly. Their sum reproduces every voxel of both full cubes within
the predeclared printing precision; independent verification passes
**27,992/27,992** checks. Neither component is separately normalized.

The smooth cubes have no negative native samples. Every negative total sample
on the six predeclared outside-sphere planes has positive smooth density and
negative augmentation. This establishes component attribution, not numerical
convergence of a vacuum observable or a recognition improvement. No clipping,
new mold, height selection, calibration or production change follows.

The PAW transformation equals the identity outside its augmentation regions
([Blöchl 1994, section II A](https://ims.ut.ee/images/3/32/Phys-Rev-B_1994_50_17953_Blochl_DFT-PAW.pdf)).
This motivates a next, separately bounded test of a vacuum-only wavefunction
observable with an independently checked reconstruction and geometric support.
It does not authorize dropping augmentation inside atomic spheres or forcing
the smooth whole-cell integral to equal the full state weight. Finite-G
leakage remains a plausible mechanism, not a completed cutoff-convergence test;
vacuum decay, constant-current extraction and experimental transferability
still need to be established for the candidate.

**September 25 Gamma reconstruction result.** Job **11976185** completed the
paired stock/patched QE 7.4.1 test on the same accepted collected states.
Filling the missing conjugate reciprocal-space entries removes the additive
floor: `stock = (patched + mean(patched))/2` holds for all **21.6 million**
voxels within the previously fixed cube-printing error bounds. All legacy STM
controls are byte-identical to their archives, corrected ILDOS repeats are
exact, and corrected integrals are **24.0001613 / 26.0000836** for selected
weights **24 / 26**. Independent verification passes **20,968** checks; all
outputs and both binaries are fetched and hash-verified. This establishes the
isolated reconstruction defect, not a new production STM observable.

The corrected sharp-window density still has negative samples on the unchanged
0.40/0.50/0.60 nm ring-relative planes: GlcN **0/1/45** and GlcNAc
**45/111/137** out of 289 at each height. At 0.50 nm the GlcNAc negative plane
sum is **14.78%** of its positive sum. Geometry-only checks using the active
UPF outer augmentation indices and periodic cell show all sampled points are
outside the declared PAW spheres, with minimum clearance **0.112 nm**.
The outer radius follows the pinned
[QE UPF reader](https://github.com/QEF/q-e/blob/qe-7.4.1/upflib/read_upf_new.f90);
the [UPF format](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
defines radial lengths in Bohr.
The subsequent component diagnostic above separates smooth and augmentation
contributions on these unchanged states and verifies their sum against the
corrected total. Finite-reciprocal-grid leakage is not yet established as the
mechanism. No clipping, height choice, SCF criterion change, mold import or
recognition grade follows.

**Earlier September 25 spectral follow-up (incomplete).** The missing accepted GlcNAc
collected state was recovered on copies with the same plain SCF input and
**5e-5 Ry** criterion. Both archived STM cubes reproduced byte for byte; the
recovered 2.46 GB GlcNAc checkpoint is now local/hash-verified. Job **11975572**
then computed GlcN plot_num=10 ILDOS in the sharp experimental −0.300 V window,
with **12** selected states, total weight **24**. Its cube integrates to
**24.0001607365**, without negative samples, but these checks do not establish
spatial correctness.

Its nearly constant vacuum floor is **2.6387e-4**, almost exactly half its
cell-average density **5.27747e-4**. Official
[`local_dos.f90`](https://github.com/QEF/q-e/blob/qe-7.4.1/PP/src/local_dos.f90)
fills only `nl` before its final inverse transform, unlike the normal Gamma
[`density helper`](https://github.com/QEF/q-e/blob/qe-7.4.1/FFTXlib/src/fft_helper_subroutines.f90),
which also fills conjugate `nlm`. Missing these partners would give
`(rho + mean(rho))/2`, conserving the integral but corrupting spatial contrast.
That hypothesis is now validated by the paired rebuilt-binary result above.
No algebraic offset correction of rounded cubes or new production observable
is adopted.

The separate paired build/PP job **11975972** reproduced both site controls
exactly but failed during stock compilation because the MKL build paths were
missing. Neither a rebuilt stock nor a patched binary was produced, so it
does not itself test the conjugate-fill hypothesis. Checkpoints were unchanged
and outputs were fetched. Its diagnosed MKL build correction enabled the
separately bounded continuation **11976185**, without altering the failed run
or the frozen scientific settings.

The earlier spectral job **11975572** stopped on the enlarged diagnostic grid: preparation accidentally used
half-width/step **0.64/0.08 nm**, twice the existing **0.32/0.04 nm** grid. A
header-only preflight now detects its out-of-volume samples before preparation.
GlcNAc ILDOS was not run; there was no retry, clipping, domain shrinkage or
recognition grade. This error does not invalidate the narrower historical plane
statistics below. See the dated journal entry for complete partial evidence.

**September 24 physical qualification.** The accepted GlcNAc cube contains
**131,984 negative native values / 10,800,000** (minimum **−0.0037643**); GlcN
contains none. On the existing 17×17 ring-frame plane at 0.50 nm, **35/289**
GlcNAc samples are negative. This precedes normalization and is not created by
trilinear interpolation. The [QE 7.4.1 STM source](https://github.com/QEF/q-e/blob/qe-7.4.1/PP/src/stm.f90)
weights out-of-window states with `w0gauss`, using the SCF smearing type/width.
Its [cold-smearing derivative](https://github.com/QEF/q-e/blob/qe-7.4.1/Modules/w0gauss.f90)
is negative for x>√2. This observable is a potentially **signed weighted sum**,
not a guaranteed nonnegative sum restricted to the bias window. This identifies
a mechanism consistent with the cube; no wavefunction-level reproduction has
yet attributed every negative value to it.

The legacy surface search also includes the rise toward the periodic substrate
image, then out-of-domain samples. At the 0.50 nm target the GlcN maximum is
**2.48346e−7**, below the legacy isovalue lower bound **3.16228e−6**; retaining
that bracket and interpolating cannot recover the molecular-side target.
Ring-relative and Cu-relative heights are not interchangeable. No clipped or
absolute-value density, arbitrary branch cutoff, new template, SCF criterion,
registry replacement or recognition result follows this audit. The **5e−5 Ry**
acceptance and plain/TF repeatability remain true, but neither guarantees a
nonnegative STM observable or calibrated isosurface. Evidence:
`results/qe_surface_audit_20260924/` and the dated journal entry.

**Separate reconstructed-assignment reader warning (September 24).** Its
historical `cc_mold_native.jl` preserves a first-index-fast token interpretation
to replay the lost champion lineage. The QE 7.4.1 writer emits the third index
fastest. An opt-in correction in `qe_cube_molds.jl` now adapts that input order;
`config/qe_cube_order.toml` leaves all other legacy surface settings unchanged.
This does not alter the accepted cubes, SCF criterion, registry, or the separate
constant-height/constant-current helpers below, which already use last-index-
fast addresses. Numerical Python/Julia parity of the historical reader was not
validation of its physical coordinates. Benchmark comparison is separate from
this file-format correctness decision; see the September 24 journal entry.
The completed fixed-calibration comparison gives **676/870 correct, 27/145
exact**, versus exact legacy **679/33**, both coverage **852**. This does not
promote the changed assignment pipeline or validate the legacy coordinates.
The corrected reader remains opt-in; no SCF, height or support retuning follows.

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
current: `plot_num=5` remains a smearing-weighted (potentially signed) sum of
`|psi_n(r)|^2` around the bias window, with no justified conversion to nA.

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
