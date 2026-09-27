# Calibration: measurements, assumptions and limits

## Fixed-density diagonalization is not current calibration

The bounded follow-up in `config/qe_diagonal_precision.toml` compares
**1e-10 / 1e-12 Ry** NSCF eigenvalue tolerances at the same accepted
**60/720-Ry** density for each molecule. Both use full occupied/empty-state
accuracy and start from independent identical copies. No SCF acceptance,
geometry, PAW data, physical query, spectral bias or normalization changes.
The two NSCF results, not just their differences from the earlier SCF,
provide the fixed-density tolerance comparison. Eigenvalue convergence is
not a pointwise orbital-error bound. An unchanged or changed vacuum tail
alone cannot establish physical accuracy, calibrated current or chemical
recognition. No new mold or champion is adopted by this diagnostic.

## Wavefunction-basis sensitivity is not a new calibration

After the 360/720-Ry density-cutoff comparison, the user approves a separate
**50 → 60 Ry wavefunction-cutoff** test at fixed **720 Ry** density cutoff
on both accepted geometries. Two sequential jobs, each at most **8 CPU /
96 GB / 4 h**, compare the same bias and three physical planes. The baseline
must replay exactly; candidate SCFs keep the same physical settings, original
per-state targets and common **5e-5-Ry** acceptance. They start from atomic
potential and atomic+random orbitals in the enlarged basis, with no
cross-cutoff restart or mixing history. This different numerical starting
path and finite SCF tolerance limit causal attribution of a two-point change.

The previously observed boundary maxima and persistent nonmonotone columns
are retained, not removed by cropping, normalization or height choice. The
new test measures sensitivity of the smooth bias-window densities; it does
not by itself establish convergence, a monotone physical tail, a calibrated
current-height relation or improved recognition. Frozen molds and champion
remain unchanged. Both jobs **11993413/11994503** complete successfully,
with exact baseline replay, all four parallel/serial repeats exact, and
**59,312** independent saved-output assertions per molecule passing on
Viper and locally. Original-state and fetched checkpoint hashes pass.
Relative-L2 changes across the three heights are **0.486/1.142/3.527%**
for GlcN and **0.921/1.716/3.606%** for GlcNAc, while individual 0.60-nm
ratios span **0.596–1.235 / 0.613–1.427**. GlcNAc still has **1/289** rising
columns across 0.40–0.50 nm and **30 → 28/289** across 0.50–0.60 nm.
Every sampled density maximum remains on the diagnostic window boundary.
This completed sensitivity test does not qualify a new calibration or
recognition improvement; see the full paired results in the
[DFT calculation note](dft_calculation_note.md).

## Potential qualification precedes vacuum-tail propagation

The separately authorized **360-to-720-Ry** density-cutoff test holds PBE,
`ecutwfc=50 Ry`, the two accepted geometries and other physical settings fixed.
It measures signed density and total/electrostatic/XC plane statistics, not a
new STM observable. Native grids differ; every plane is retained, with the
same geometric gap boundaries and no interpolation or favorable-height choice.
Both authorized jobs complete successfully and their full local copies are
checksum-verified. Smaller minimum total/XC variations and negative-charge
integrals establish sensitivity, not a converged tail, current calibration or
permission to replace the molds. Negative samples persist and GlcNAc's
sampled maximum total-potential variation increases; the change is not a
uniform improvement across the native domain. See the paired results in the
[DFT calculation note](dft_calculation_note.md). The common **5e-5-Ry** SCF
acceptance and the original per-state electronic targets remain unchanged.
Exactly two sequential eight-CPU jobs are used, each within four hours; no
additional SCF, calibration or production change follows.

The completed saved-orbital follow-up confirms that this cutoff sensitivity
also changes the bias-window smooth density outside PAW spheres. Job
**11992917** evaluates the same **0.40/0.50/0.60-nm** planes directly, with
no new SCF, changed bias, matched band-count prior, normalization or selected
height.
All columns remain visible. The two 360-Ry references reproduce exactly,
all four parallel/serial controls pass and **129,555** independent saved-table
checks pass. Relative-L2 changes are **1.736/6.252/12.145%** for GlcN and
**1.308/0.712/3.689%** for GlcNAc across the three heights. GlcNAc still has
**1/289 and 30/289** rising columns across the two adjacent-height intervals;
the upper interval's largest ratio increases **11.514 → 12.491**. A small
global norm change therefore does not qualify local tail accuracy, and
raising the density cutoff does not remove this sampled branch ambiguity.
Nor does the comparison prove that every rise is numerical. Implementation
checks do not establish cutoff convergence, calibrated experimental current,
a unique constant-current branch or a recognition gain. The existing frozen
molds remain unchanged; no height or normalization is selected from this test.

A geometrically empty region is not necessarily a region of constant local
effective potential. The completed PP-only diagnostic exports the accepted-state
total potential and its electrostatic control, retaining every native lateral
pixel and height. It tags the geometry-only all-atom PAW-free gap without
selecting a matching height or flatness threshold. This supplies necessary
physical information for a possible tail model, not a new height/current
calibration or proof that the present reference spread is numerical noise.

In Ry/Bohr units the zero-lateral-momentum constant-potential equation gives
`kappa^2=V_Ry-E_Ry`; XML energies must first be converted from Hartree.
Lateral potential variation couples Fourier modes, vertical variation can
invalidate a constant barrier, and a suitable potential alone does not prove
the boundary wavefunctions are accurate. Report these limitations and the
sampled variation before attempting a continuation. The material-specific
height and flatness criterion in [Park and Yoon's supplementary method](https://www.nature.com/articles/s41598-022-10870-0)
are not calibration values for these chitosan cells. No production change is
made by this diagnostic.

Job **11981276** exports and checks both accepted states without changing
their density or wavefunctions. Across **122/115** PAW-free planes, the full
lateral total-potential range is never below **1.2921/1.2939 eV** and spatial
SD is never below **0.2038/0.2021 eV**. Among band/plane pairs positive at
every sampled lateral point, the range/mean-barrier ratio is at least
**0.3283/0.3072**. These are descriptive extrema, not a new physical
acceptance threshold or a bound between grid nodes. No constant-potential
matching region is qualified or selected by this result.

XC spatial SD is larger than the electrostatic SD on **235/237** planes.
The PP logs also report negative-density warnings, but these global messages
do not establish where those values occur or explain the local XC variation.
Do not silently substitute the smoother electrostatic potential for the
effective potential, clip density, or impose a borrowed vacuum level. The
density/XC treatment and numerical-cutoff dependence need separate examination;
the accepted electronic convergence criterion is not changed by this finding.

The attribution diagnostic keeps the accepted density and PBE functional
fixed and separates local-density XC, local GGA and negative FFT divergence.
It also exports the actual valence-plus-core input and original GGA activity
mask. The internal driver tests absolute density and **squared** gradient;
this does not authorize taking the absolute value of a density used as an STM
observable. An inactive local branch may still receive a nonzero FFT-divergence
term from elsewhere. Attribution and reproducible sums cannot establish a
cutoff-converged tail, an artifact's size or a new physical barrier. In
particular no component may be dropped merely because that makes the potential
flatter. No new calibration or recognition result belongs to this diagnostic.

Its first attempt **11982255** fails the required cross-build density identity
control before repeats or GlcNAc processing. Small first-difference witnesses
do not establish a whole-volume error bound or validate the physical
decomposition. The next control must separate the build effect from the
instrumentation effect; the guard is not relaxed and no component is adopted.

The separately bounded paired-build follow-up holds configure options and
compiled libraries fixed, changing only the PP extraction hook/helper.
Stock density/total controls and component repeats must agree byte for byte
inside that build. Cross-build discrepancies are reported at every native
voxel and plane without choosing an acceptance tolerance from their size.
Only after both states pass their same-build checks can component sums be
compared with that build's total minus electrostatic control. Even a successful
attribution would be a numerical diagnostic, not density/cutoff convergence,
a repaired vacuum Hamiltonian or a recognition improvement.

That follow-up, **11988093**, passes the within-build identities and internal
native/cube/additivity checks, but fails **313** independent descriptive
relative-L2 checks. A Cartesian-view summation defect is reproduced without
physical data and corrected locally, without relaxing the **16-eps** check.
The independent component-cube suite was not reached; that run and its derived
tables/figures remain provisional. Separate saved-export recheck **11989810**
passes all **92,704** independent control/norm/component checks with unchanged
source hashes and precision bounds. No new QE computation is used.

Attribution is now qualified at the fixed accepted density. Local-density XC
has the largest individual component SD on all **61/61 and 57/57** lower-gap
planes (minima **0.14727/0.14941 eV**). Across the full gaps it is largest on
**93/122 and 88/115** planes, with divergence largest on the remainder. These
are SD rankings, not variance fractions; cancellations prevent interpreting
their ratios as explained variance. Native valence-density negative fractions
reach **50.07/49.78%** within a plane despite small whole-cell negative-charge
integrals. Neither a GGA-only cause nor accurate local tails follow from those
integrals. This does not establish cutoff convergence or justify a matching
plane, density clipping, calibration or potential replacement.

## Atomic z origin is not an observed STM background

Subtracting the lowest (or highest) Cu atomic coordinate from a calculated
constant-density tip height does not supply the apparent substrate height at
that density. The experimental masked background and this atomic origin are
different conventions. The bounded off-projection diagnostic therefore queries
the same accepted wavefunctions and three saved isovalues at a fixed lateral
grid, outside the union of both molecular PAW projections, with independently
audited PAW-free vertical segments. It retains full-gap, Cu-half and old
molecular-half outcomes without choosing a favorable branch/domain.

The motivation follows the tip-position LDOS dependence of the
[Tersoff--Hamann model](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.31.805);
the proposed background comparison is an inference to test, not a measured
current conversion. Molecular tails and neighboring periodic adsorbates may
still dominate off-projection points. A scalar background is only defensible
if its spatial and state dependence are understood; no spread cutoff, class
offset or chemical confidence is inferred here. Even a supported reference
cannot cure the shape mismatch left after free local affine profiling. No
production map, current, count, assignment or champion is changed.

The completed diagnostic **11980577** passes **582,693** independent checks.
On paired Cu-half valid sites, background-difference spatial SD is
**64.51--72.66 pm**, larger than the original mean molecular contrast
(**20.63--21.01 pm**). Hypothetically referencing each map to the same lateral
site can produce either contrast sign. The first site gives **+47.78 to +55.49
pm**, while averages across their paired valid supports give
**-20.94 to -16.75 pm**. Neither is adopted as chemical calibration. These
different counterfactual references are not experimental uncertainties or
evidence for composition. The old molecular-half domain has different and
smaller supports, and all full-gap cases have multiple roots.

The geometric screen excludes atomic PAW disks, not the entire molecular
footprint or its electronic tails. This motivated a saved-output follow-up
with a declared whole-footprint/common-exterior geometry rule, keeping all
original outcomes. The observed distance dependence does not authorize picking
a favorable clearance threshold or a desired chemical-contrast sign.

The saved-output follow-up defines the whole footprint as the convex hull of
all non-Cu atomic centers plus a disk of the largest active molecular PAW
radius. Strict periodic exterior common to both states is fixed before
conditional heights are aggregated; the removed complement and every
domain/isovalue remain. This deliberately covers intramolecular spaces, but
PAW radius is not an electronic decay scale. Exterior membership or lower
spatial spread alone cannot establish a clean-Cu plateau, transferability,
chemical calibration or physical uncertainty.

This follow-up is now complete: **325** common exterior candidates yield
**307/312/317** paired valid Cu-half roots. Background-difference spatial SD
remains **65.42/68.17/73.94 pm**, and hypothetical referenced contrasts still
span both signs. The geometric exclusion alone does not stabilize the
reference. No extra clearance, scalar class offset, domain or isovalue is
selected. Testing a potential-qualified vacuum continuation is a distinct
future hypothesis, not a demonstrated cause or an adopted replacement.

## A shared envelope must not silently remove chemical contrast

The September 25 shared-image diagnostic fits one source-only spline surface
to the complete usable image after the frozen masked background. Its GCV
resolution is selected before target scoring or DFT perturbation probes.
The trace includes both linear stages and their overlap, conditional on the
fixed mask; source-derived masks and correlated residuals still prevent an
unconditional risk claim. No effective-sample-size placeholder is reinterpreted.
The trace formula follows the general linear-smoother GCV expression in
[Golub, Heath and Wahba (1979)](https://pages.stat.wisc.edu/~wahba/stat860/pdf1/golub.heath.wahba.pdf);
it does not establish the chosen knot scales as physically correct.

Probes add the saved physical GlcNAc-minus-GlcN difference on one old native
patch and zero elsewhere. Both background and envelope responses are computed
with the selected resolution and mask frozen. Their residual energy can reveal
signal absorption, but is not chemical classification accuracy. The masked
background is oblique, so ratios are not clipped to [0,1]. Probe boundaries,
unvalidated anchors and the absence of re-selection under perturbation limit
any physical interpretation. A good predictive envelope can still destroy
chemical information. Conversely, retaining a perturbation does not show that
real-image chemistry is identifiable. No new current, height, count, composition
or confidence calibration is introduced.

The completed job **11980183** passes **132,131** independent grouped checks.
On 141 identified-registration views, the class-free envelope lowers matched
bright-support prediction MSE by **24.537%** against source copy, improving
120 views and worsening 21. Yet it removes approximately **48%** of the
injected contrast energy; some probes lose over **90%**. Every source-GCV
choice is at the finest declared grid. Better image prediction therefore
does not identify a chemical-preserving resolution or validate subtracting
the envelope before assignment. A jointly constrained observation model is
still needed. In particular, the DFT lowest-Cu-atom reference and the fitted
experimental STM background are different conventions; whether a supported
off-molecule reference at the same isovalue resolves part of this mismatch
is an open physical question, not an adopted calibration or a promised fix.

## Chemical choice versus one fixed physical shape

A gain from a two-shape dictionary over a plain local plane does not by itself
establish useful chemical selection. The fixed-shape replay compares both
saved candidate predictions with the saved source-selected prediction, at the
same geometry, gain, nuisance coefficients and pixel support. Neither fixed
candidate is an estimated composition or an assignment to the real molecule.
Masks are post-scoring annotations only; source/target pixels, geometry,
isovalue and type are not chosen again. Exact prior selected/control error
replay and independent interpolation check the comparison. Both fixed arms,
both masks, all isovalues, all failed input statuses and the unfiltered and
rejected groups remain. No new physical/chemical calibration is inferred from
this reused-cohort diagnostic, and no confidence cutoff is introduced.

The completed replay (job **11979724**, **3,434,090** checks passing) shows that
fixed GlcN alone explains **93.78%, 99.18%, 100%** of the earlier retained-support
gain over a plain plane. The selector adds only **0.07417%, 0.01299%, 0%** in
MSE relative to that fixed shape, improving **2, 2, 0** of 138 identified views.
At the third isovalue, all four different GlcNAc choices lack heldout target
pixels; equal target error is not validation of those choices. The full
identified comparison slightly favors fixed GlcN at every isovalue. This
supports a generic shape contribution, not robust chemical identification or
an all-GlcN composition. Shared molecular-envelope/background modelling remains
an open observation-model question. The shared-image experiment above tests
prediction and contrast absorption, not chemical calibration.

## Bright support is not molecular identity

The two-pass image-foreground diagnostic estimates a source-only bright mask,
guards it and fits one background x slope jointly with supported-row offsets.
Within-row centering avoids the sequential plane/row x-slope ambiguity under
unequal coverage, conditional on the estimated exclusion mask. It does not
prove that the exclusion contains the whole molecular signal. Centered heights
and coefficients are saved, without treating a tiny numerical floor as noise.
The exact Otsu split plus raw horizontal-difference floor is an explicit
engineering detector, not physical count or chemistry calibration. All
components, weak/empty outcomes and unsupported pixels remain reported;
there is no inferred chemical composition or largest-component prior.
Full-row signal cannot be separated from unrestricted row offsets by this
model, and bright substrate/contamination may survive its mask. Cross-view
overlap is descriptive; common observed coverage and the unfiltered predictions
must remain visible before interpreting any conditional fit improvement.

The completed source-only comparison (job **11979519**, all 146 scans) improves
pooled mask agreement **0.94707 -> 0.95909** on 141 identified-registration
views at equal common coverage. It loses only 45 additional smoothed pixels,
without losing any old patch center. This is a stability measurement, not
foreground ground truth. On the 138 previously physically scored identified
views, the unchanged chemical-plus-plane predictions on retained anchors give
**5.884--5.924 pm** RMS versus **5.930 pm** for a plain local plane: a weak
**0.22--1.55%** pooled-MSE gain, present in only about half the views. The same
model chooses the GlcN shape for **99.34--99.74%** of these anchors. The completed
fixed-shape control above explains almost all of the gain without useful
broad chemical selection. No expected composition is used to reject
or repair these frequencies. All cases, masks, rejected anchors and isovalues
remain reported; the unfiltered result still favors the plain local plane.
Neither counting calibration nor the frozen application/champion changes.

## Distinguishing chemical shape from local height trends

The surface-control diagnostic holds the raw-derived observation/calibration
boundary fixed. Its height-only contrast uses the same common map with two
constant shifts, fixed by the physical maps' mean difference, and a shared
source-scan offset. It does not calibrate shifts from experimental composition.
Separate local constant and plane controls measure how much the source alone
predicts without any chemical template. Both plane-template arms use the same
source-only nuisance basis; the common arm fixes geometry before the two
chemical shapes are compared. Removing each candidate's affine component is
an explicit diagnostic loss of height/slope information, not a new physical
normalization for production. Evaluation pixels never determine a source
parameter. Better cross-view prediction would still not prove foreground,
correct count, chemical labels, independent noise or a calibrated current.
These controls follow inspection of the same cohort's earlier predictions;
the reused heldout rows are a diagnostic comparison, not fresh validation.

**Completed result:** job **11979008** preserves all 146 cases without a new
failure. On identified-registration views, height-only contrast retains
**99.52--99.57%** of the old chemical state choices, despite explaining only
about **40%** of their predictive gain over the original common surface.
Plane-plus-chemical shapes improve on plane-plus-common, but remain
**21.13--23.12% worse in MSE** than a local plane with no chemistry
(**5.995--6.044 pm** versus **5.447 pm** RMS). After plane profiling,
**99.48--99.72%** of all diagnostic patches choose the same GlcN hypothesis;
this is not a composition estimate. Neither the new comparison nor its
**2,584,957** passing arithmetic checks justify unit assignment or physical
current calibration. Foreground and a shared molecular envelope remain to be
established without a desired count; removing affine height here does not prove
chemical information absent from better-localized observations.

## Experimental transfer without inherited count calibration

The bounded surface-transfer diagnostic uses only the raw source image to
locate bright extrema and fit each of the three fixed surface families. These
disjoint local patches do not estimate `N_selected`, splice height maps into an
additive chain model, or use the old benchmark-informed width/support settings.
The source plane is removed using observed pixels; row flattening and smoothing
locate anchors only. Fits retain unsmoothed physical nm values. Both DFT species
share unit gain and one source-scan height offset, preserving their predicted
relative height and corrugation; neither class is individually standardized.
Lateral evaluation is bilinear in the saved 17x17 height grids, not a fresh
continuous-WFC root calculation at every experimental coordinate.

Only alternating buffered target-row blocks enter translation and height-plane
calibration. Complementary blocks evaluate the frozen source prediction. Target
evaluation pixels cannot affect source anchors, discrete chemical/geometry
states, the shared offset or target corrections. Both directions, all isovalues,
the chemistry-free mean surface, a constant and source-copy predictions remain
reported. An unresolved registration is flagged, not silently declared physical.

This tests a local height-shape hypothesis with a finite exploratory geometry
grid. It does not establish tip/current conversion, context transfer from the
DFT trimer, physically accurate calibration, independent noise or chemical
recognition. Cross-view improvement must not be read as a benchmark grade, and
image-derived anchors must not be relabeled as validated molecular counts.

**Completed result:** job **11978725** scores 141 scans, retaining one
unsupported-resolution case and four bias mismatches. On the 138 views with
identified registration, all three chemical dictionaries improve on the common
surface (**44.80--45.36 pm** versus **66.83--67.23 pm** RMS), but source copy
remains far more accurate (**6.63 pm**). A fixed first-scan overlay places many
anchors on background/borders and associates the chemical hypothesis with
molecular brightness. This is evidence of a possible background/height confound,
not chemical truth. The next comparison needs an explicit background and a
class-independent molecular envelope, including a height-only contrast null.
Neither the predictive gain nor the **7,281,567** passing saved-output arithmetic
checks validate chemical assignment, justify a count, or calibrate a setpoint.
No production calibration is changed.

## Continuous roots do not choose an experimental setpoint

The follow-up at the three saved midpoint isovalues replaces vertical linear
interpolation by the finite Fourier series of the same accepted states. It
keeps the earlier molecular-half domain and every lateral pixel. Interval
enclosures include rounding in lateral phase evaluation, coefficient regrouping,
vertical harmonics and density/derivative evaluation. All harmonics contribute;
Taylor polynomials bound ranges with an explicit remainder, not a modified
physical density. The implementation uses
[IntervalArithmetic.jl](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/usage/).

An unresolved, multiple, ascending or boundary crossing cannot yield a valid
surface pixel. Refinement cannot turn a failed map into a complete map by
discarding pixels. The 1e-7 nm root-width limit is a numerical target, not an
experimental height uncertainty. A continuous-root check validates the
finite-basis series on this fixed domain, not the adequacy of periodic images,
the tip, spectral convergence or inherited experimental calibration. No one
of the three isovalues is chosen using a recognition grade.

**Completed result:** job **11978192** certifies all **1,734** representative
roots as unique and descending, with no unresolved region and exact
serial/parallel repetition. Six 17x17 maps remain complete. The independent
3D signal's maximum relative departure from the isovalue falls from **2.089%**
to **0.00001779%**; the largest change from the earlier linear direct-profile
height is **0.668041 pm**. Every root enclosure is at most **1e-7 nm** wide.
These are numerical errors in the fixed model, not physical height accuracy
or experimental recognition. Only the three representative isovalues have
been certified, not the whole open intervals containing them.

The next scientific question is experimental transfer with a shared
label-free observation model. Separately normalizing each chemical surface
would remove part of their predicted contrast. Adding constant-current height
patches also requires justification: heights are roots of a density equation,
not additive density contributions. The existing height-template assembly
cannot silently supply that justification. No mold or calibration is promoted
by this numerical result.

## Surface families without a selected experimental current

The geometry-defined crossing diagnostic changes the sampling direction from
the earlier tilted ring-normal planes to the substrate normal. This is an
explicit new diagnostic coordinate convention, not a silent production change.
In the [Tersoff-Hamann model](https://doi.org/10.1103/PhysRevB.31.805), tunnelling
current is related to sample density at the tip. A common isovalue is therefore
a useful model observable, but its conversion to the experimental setpoint
still requires tip and spectral assumptions not supplied by this calculation.

The vacuum search bounds come from the accepted atomic geometry and PAW radii,
before reading densities: all space between the two periodic PAW envelopes,
and its molecular-side half. The half-gap convention is not evidence that the
periodic image's wavefunctions are negligible. Comparing both domains retains
that limitation explicitly; neither can be adjusted after a failed result.

All 578 columns enter each common-isovalue intersection. A missing or ambiguous
column is not removed to obtain a surface. The complete piecewise-linear
interval family and empty results are preserved, with direct Fourier checks
of the knot profiles and representative surfaces. No experimental height,
per-class normalization, current-to-density factor, threshold or benchmark
label is used. Uniqueness between continuous Fourier samples and transfer to
experimental chain images remain distinct questions, even if surfaces exist.

**Completed result:** job **11977798** finds three joint open intervals in the
molecular half-gap and none in the full gap. The intervals include every pixel
and overlap the independently reconstructed and interpolated-profile results.
Six representative maps are exported, with exact repeat and no missing pixel.
Their largest direct-density departure from the requested isovalue is **2.089%**;
the maximum cube/direct-profile height difference is **1.427 pm**. Every half-gap
profile nevertheless contains an increasing segment, and four GlcNAc columns
have disjoint admissible intervals. This result does not impose monotonicity.
The geometric half-gap, finite z sampling, finite-basis periodic slab and unknown
tip/current conversion remain limitations. No one of the three intervals is
chosen as a chemical model or experimental setpoint.

## Direct wavefunctions check the vacuum observable, not the current setpoint

The bounded September 25 WFC diagnostic independently reconstructs the smooth
sharp-window density from the accepted reciprocal-space coefficients, without
calling QE's `local_dos` or a density FFT. Its amplitude follows the stored
wavefunction normalization, spin/k-point weights and cell volume; neither a
whole-cell renormalization nor a fit to experimental contrast is permitted.
All queries and their native interpolation vertices must lie outside the
periodically repeated PAW spheres. This never removes augmentation from an
inside-sphere observable.

Full-precision XML geometry, not the six-decimal printed cube coordinates,
defines native-grid point locations. Their numerical agreement is checked
before scientific comparison; both interpolation conventions are retained in
the diagnostic. All three previously declared ring-relative heights are
reported, including zero or rising density if present. Agreement with native
QE samples and plausible vacuum decay would not establish cutoff convergence,
a tip model, an absolute current conversion or a benchmark improvement. No
constant-current isovalue or experimental height is calibrated in this cycle.

**Completed result:** job **11977502** verifies all 13,872 required native
vertices, the Parseval integral and exact serial/parallel repetition. All 1,734
direct plane samples are positive. GlcN decreases at every lateral point over
both sampled intervals, but GlcNAc increases at **1/289** points over 0.40–0.50 nm
and **30/289** over 0.50–0.60 nm, up to **11.51-fold**. These columns are retained
unchanged. This verifies the finite-basis reconstruction in the audited domain,
not a monotonic vacuum tail, unique isovalue crossing, absolute current or
strict label-free experimental calibration. No mold or inference change follows.

## PAW component attribution precedes any new observable calibration

The September 25 component test compares the unchanged corrected total with
separately calculated smooth and augmentation contributions. It is diagnostic
only: neither component is a new mold, and there is no per-component
renormalization, clipping, spectral/height/cutoff sweep or label access.
The full-density reference and source-state hashes are fixed before the job.
Job **11976981** verifies additive reconstruction throughout both cubes,
exact reference/repeat reproduction and **27,992** independent assertions.
The smooth cubes have no negative samples; every negative total-plane sample
has negative augmentation and positive smooth density. This positive smooth
component does not alone validate a calibrated STM observable or justify
removing augmentation inside its atomic support. A vacuum-only candidate
still needs independently checked reconstruction, geometric support, decay
and constant-current diagnostics before any experimental comparison.

## Gamma reconstruction is a software/observable check, not recalibration

The separate September 25 paired PP-only comparison tests the predicted
`stock = (patched + mean(patched))/2` relation from missing Gamma Fourier
partners. It rebuilds the same QE source with and without one conjugate-fill
instruction and uses unchanged accepted wavefunctions. Controls cover the site
binary, rebuilt legacy STM, exact patched repeats and selected-state integrals.
The whole-volume mean in the relation is fixed by the candidate cube, not fitted
to vacuum pixels or STM labels. Numerical bounds come from QE's printed
precision plus the explicit roundoff allowance.

The actual historical 17×17 lateral grid is **0.32/0.04 nm** half-width/step,
validated against the source config and both volume domains. All three declared
heights remain diagnostics; no height, current, isovalue or chemical threshold
is calibrated here. Job **11976185** passes the frozen software checks, but
does not establish a nonnegative vacuum observable, a physically converged STM
model, label-free historical calibration or improved experimental recognition.
At 0.50 nm the corrected GlcN/GlcNAc planes still contain **1/111** negative
samples out of 289 each, outside all declared PAW augmentation spheres. No
height or clipping rule is selected from this result; the component attribution
above supplies a diagnosis, not a new physical calibration.

## Explicit spectral-window diagnostic

The September 25 QE diagnostic compares the archived plot_num=5 STM proxy with
plot_num=10 ILDOS from the same electronic states. The candidate uses a sharp
zero-temperature window at the existing −0.300 V experimental bias, with no
fitted broadening. It does not claim that the Gamma-only slab or zero-temperature
approximation is spectrally converged, or that an ILDOS value is a current in pA.
The code checks the factor-two Hartree/Ry conversion and uses the full-precision
XML Fermi energy, not its rounded printed value.

Unlike legacy STM's cold-smearing derivative tails, the selected ILDOS state
weights are nonnegative. However, QE's two paths also differ in PAW/US
augmentation and Fourier interpolation. Positive weights alone do not certify
nonnegative grid samples; the diagnostic preserves and counts every negative
value, checks the integrated density against the selected state weights, and
exports all three existing local height planes without selecting a height.
The integral tolerance is numerical, not a chemical-confidence threshold.

Missing accepted GlcNAc collected wavefunctions are restored only in an isolated
copy of the surviving pre-acceptance snapshot, with the unchanged **5e-5 Ry**
criterion. Byte-identical archived STM controls are required before comparing
ILDOS. No numerical tolerance or SCF criterion is relaxed after a failure.
There is no new isovalue calibration, template, recognition result or production
change. See [configuration](config.md) and the dated journal entry.

**Observed outcome:** both legacy cubes reproduce exactly and the accepted
GlcNAc collected checkpoint is recovered. The candidate comparison is incomplete:
preparation mistakenly doubled the existing grid's half-width and spacing,
causing an out-of-domain refusal after GlcN ILDOS. A new header-only preflight
catches this before computation. The smaller historical assignment domain has
not thereby been shown invalid. No GlcNAc ILDOS was run and no retry followed.
GlcN ILDOS integrates to **24.0001607365** versus expected **24**, but its vacuum
floor matches half its mean density. Missing conjugate-G completion in QE's
Gamma post-processing was a source-supported hypothesis at that stage; the
later paired-binary result above validates it. A conserved integral does not
establish spatial correctness; no rounded-cube offset subtraction or template
promotion follows.

## Numerical background conditioning is not recalibration

The September 25 opt-in `background_conditioning.toml` comparison changes only
coordinates used by the numerical optimizer. A centered, QR-normalized basis
represents the same tilted background, and exact transformed inequalities retain
its original feasible set. Molecular amplitudes, positions, widths, overlap,
saved count/family, pixels, support, noise and full parameter count do not change
as model definitions. Their fitted endpoint values may change. The target's
population standard deviation supplies a numerical coordinate scale, not a
physical amplitude bound or independently measured noise.

Four cases are chosen by predetermined ranks of saved numerical repeatability,
not external counting or chemical errors. Same-start repetitions and a separate
forward/gradient check measure numerical stability and stationarity. A smaller
RSS or a better-conditioned Jacobian does not establish better counting or
recognition. Existing inherited calibration provenance remains benchmark-informed;
this diagnostic cannot certify strict label-free calibration or an optimum.
No full-cohort benchmark or production change is included in this experiment.

Completed result: the two SLSQP bases reach the same RSS, with no endpoint
passing the unchanged stationarity audit. Centering improves that gradient on
one scan, worsens it on three and costs 21.39% more optimizer time; it is not
retained. SLSQP's 10.42–29.25% RSS improvement over time-capped LM is a numerical
result on four saved-count cases, not calibration or recognition validation.
All 24 fits are natively valid and both SLSQP repetitions match exactly. Source
**a706593**, job **11974388**, **0:0 in 14m43s**; 1,063 output checks pass.
No parameter, production or champion change follows.

## Width-conditioned counting is not recalibration

The September 24 counting diagnostic keeps all historical numerical physical
bounds. Its hypothesis is narrower than independent calibration: the existing
gap floor protects against **maximum allowed** width even when the fitted
lobes are narrower. With `c=sqrt(-2log(max_overlap))`, the alternative floor is
`max(spacing_min, c*sigma_min)` and each lobe's parallel/perpendicular sigma is
bounded above by its nearest adjacent axial gap divided by `c`, as well as the
original maximum. Shared width types use the tightest member cap.

For Gaussians this guarantees that the radial envelope evaluated at another
center is at most the unchanged overlap limit. It is conservative about
anisotropy and lateral displacement, not an integrated overlap probability or
a measured STM response. The local arm uses actual pair widths for this
diagnostic; legacy mean-width reporting stays unchanged in the default arm.
Both arms use common pixels/support and exhaustive full-parameter GCV counting.
The bounded two-repeat test does not run chemical assignment or promote a
champion. Historical benchmark-informed width/support provenance remains an
independent limitation; no apparent profile median or expected N calibrates
this coupling. Settings and failure rules are in [configuration](config.md).

The completed September 25 comparison is negative: both repeats give **42/145
exact N**, versus exhaustive control **101/145**; summed count error increases
**58→165**. The envelope checks pass, but the candidate mostly increases
overcounting under the fixed solver budgets. No LM fit declares convergence.
Reject it as a counting replacement; no post-grade bounds or calibration are
introduced. `results/local_sigma_counting_20260924/report.md`.

### What crystallography can and cannot supply

A September 24 literature check finds an independent structural reference,
not ready-to-use STM bounds. Hydrated chitosan has a twofold helical structure
with an axial crystal repeat of 10.34(4) Å
([Okuyama et al., 1997](https://doi.org/10.1021/ma970509n)); the anhydrous
structure has a 10.311(7) Å axial repeat
([Naito et al., 2016](https://doi.org/10.1002/bip.22818)). Dividing the twofold
repeat by two gives an **inferred crystal-axis advance of about 0.516–0.517 nm
per residue**, not a measured spacing between STM maxima on Cu(100).
Chitosan salts also exhibit distinct helical conformations
([Lertworasirikul et al., 2003](https://doi.org/10.1016/S0008-6215(03)00145-9)).

The transfer from atomic structure to apparent STM shape needs separate
justification. For example, the surface-specific sucrose study combines
subunit-resolved Cu(100) imaging and molecular modeling
([Abb et al., 2019](https://doi.org/10.1002/anie.201901340)); it does not
calibrate chitosan lobe widths. None of these papers supplies this pipeline's
Gaussian sigma range, support padding or uncertainty. No crystal period,
ad hoc tolerance around it, or known chain length is inserted into the frozen
comparison. An independently justified surface/contrast calibration remains
open.

## Two-view measurement audit (2026-09-24)

`test/measure_calibration.jl` now reports observed apparent widths/spacings,
missing measurements and the former bootstrap's fallback usage. **It no longer
emits a ready-to-use production TOML.** A local-window or half-prominence width
is not automatically an isolated monomer FWHM; overlap, smoothing and baseline
affect it. The existing `config/chitosan_auto.toml` remains historical evidence,
not a certified measurement-derived calibration. The full146 measurement audit
is complete; no new recognition score or production parameter is claimed.

**Result:** old widths require fallback on **245/292 views (83.9%)**, spacings
on **261/292 (89.4%)**. The observed-only diagnostic yields apparent widths on
292 views and spacings on 245. This greater availability does not establish
better agreement: on the same 18 available width pairs, median relative
forward/backward difference is **14.97% legacy versus 35.26% new**; on the same
nine spacing pairs it is **28.57% versus 31.72%**. Relative difference is
`2|fwd-bwd|/(fwd+bwd)`. On all 146 new width pairs the median is 17.33%, and
on 117 new spacing pairs it is 12.50%; these larger subsets are not a direct
comparison to the sparse legacy pairs. Forty-seven views yield a single
apparent peak, not isolated monomer evidence. **No automatic physical
calibration is justified by this experiment.**

Source **76934a6**, job **11962552**, **0:0 in 55s**; all 146 scans, 584 rows and
25 outputs/logs are local/hash-verified. Measurements and comparison read no
benchmark labels, expected N or saved geometry. Historical choice provenance
is not erased by this diagnostic. Complete results and limits:
`results/calibration_measurements_20260924/report.md` and the dated journal.

```bash
julia --project=. test/measure_calibration.jl scan.sxm \
    --config config/calibration_measurements.toml --outdir results/my_measurement
```

The directory must be new. Without `--outdir`, a single-file run uses
`results/calibration_measurements/<scan-stem>/`. Directory input requires
`--data-dir RAW --outdir NEW`; `--chunk I/N` supports the four-shard HPC audit.
`--dry-run` reads metadata/config only. Expected counts, truth, benchmark
manifests, old features and predictions are not accepted inputs; the formerly
unused `--n-lobe` argument is rejected.

Outputs are `measurements.tsv` (both actual directions and both methods),
`profiles.tsv` (bin values and pixel counts), `peaks.tsv` (accepted and rejected
measurements), raw hashes and an exact copy of measurement settings. Missing
observations remain `NA`; all scan/direction rows remain present. No substitution
of the forward channel for a missing backward channel is allowed. The new
method fits the background to observed pixels only, preserves the missing-pixel
mask, and accepts a smoothed pixel only when its entire in-image box footprint
is observed. No interpolation fills holes. Empty rows stay missing; an
unidentifiable background plane or wholly unobserved view stays unavailable.
This diagnostic helper belongs to shared `STMSXMIO`; fitter preprocessing is
unchanged. `observed_pixels` and `smoothed_observed_pixels` describe this
observational coverage for both comparison rows, not the legacy imputed support.

The legacy method preserves the old bright-pixel-only, local-window arithmetic
for audit, including the old preprocessing's median filling of missing pixels.
Rows affected by that filling are explicitly marked `legacy_imputed_*`; they
are not measurements based solely on observed pixels. Its `legacy_reported_*`
columns explicitly show what the old writer
would have returned, including defaults; the observed columns never contain
those fallback substitutions. The new method uses weighted covariance once,
all observed pixels in the axis strip, missing empty bins, plateau-centred peaks,
and interpolated widths at half topographic prominence. This follows the
definitions of [peak prominence](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_prominences.html)
and [relative peak width](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_widths.html),
without depending on SciPy. Spacing is between accepted apparent peaks in the
same finite profile segment, not an established chemical repeat length.

`pipeline_dispersion_nm` preserves the original fit preprocessor's smoothed-image
dispersion as a reference, including for the new method; it is not an
observed-mask estimate. `hf_mad_nm` measures 1.4826·MAD of flattened minus
smoothed pixels in each method, using only complete observed footprints in the
new method.
The previous tool mislabeled the former as the latter. Neither is certified as
pure background noise; the new contrast filter is descriptive, not a 3-sigma
false-positive guarantee. Pixel steps follow the shared I/O coordinate grid
(`range/(pixels−1)` before stride), not the old printed `range/pixels` estimate.

All choices are explicit in `config/calibration_measurements.toml`; see
[configuration](config.md#Diagnostic-calibration-measurements). Descriptive
thresholds preceded raw measurements; missing-data handling was amended after
the first-file preflight and then frozen before cohort execution. No setting
is optimized to recover known counts. Apparent
availability and forward/backward agreement still do not establish physical
Gaussian bounds, fit-tube width, support padding or calibrated uncertainty.

After fetching the four shards, run the label-free comparison and serialized
measurement checks (no image fitting):

```bash
julia --project=. test/summarize_calibration_measurements.jl \
    results/my_measurement_batch /path/to/raw config/calibration_measurements.toml
```

This writes `comparison.md` and `paired_measurements.tsv`, retaining unavailable
pairs instead of substituting defaults or selecting only agreeing scans.

**Historical provenance audit (2026-09-24).** Label-free inference does not
establish label-free calibration. The archived support-padding and fit-width
choices explicitly used known `N_ell=6` gains. The hybrid rule, its 0.30 upshift
threshold and the gap>=2 extension were also compared using full146 count
grades. See [the historical record](journal_archive.md) and
[the current audit](journal.md#2026-09-24-—-Fresh-hybrid-reproduction-and-historical-label-use-audit).
These settings remain frozen for reproducibility, not certified as
benchmark-independent calibration. A strict-label-free successor needs a
separately justified calibration/selection procedure without those labels;
resetting arbitrary defaults or reproducing old scores would not provide one.
No numerical parameter is changed by this audit.

Two independent raw-to-prediction repetitions now give **679/870 correct,
29/145 exact chains, 848/870 coverage**, with **123/145 exact N**. Both final
prediction files are byte-identical, despite slight geometry variation; the
saved-count 694/43/855 result is not recovered. This verifies current-run
repeatability, not strict-label-free calibration. No parameter sweep or promotion
follows. `results/hybrid_reproduction_20260924/report.md`.

The opt-in promoted-count assignment comparison uses two saved
count vectors and rebuilds both geometries at fixed N. It changes neither GCV,
the counting thresholds nor the assignment calibration. It is not a fresh
count-policy benchmark; timed-optimizer variability is reported separately.
The completed comparison gives **694/870 correct, 43/145 exact, coverage 855**,
versus fresh control **676/34/852** and history **677/36/854**: a saved-count
development result, not a raw reproduction or calibration validation. The repaired
split kernel preserves model arithmetic, but time-limited fits are not assumed
identical. No new physical parameter, threshold or composition prior is used.
Full evidence: `results/promoted_counts_20260924/report.md`.

The September 24 read-only surface audit identifies a periodic-substrate branch,
an isovalue lower bound above the entire GlcN target plane, and negative values
already present in the accepted GlcNAc cube (35/289 target samples). A new root
solver or zero-filling does not establish a physical calibration. QE's inherited
cold-smearing weights can be signed; see [the calculation note](dft_calculation_note.md).
No new calibration or template is promoted; existing recognition grades stand.

The opt-in QE cube-order correction changes physical addressing, not calibration:
the writer emits the third coordinate fastest, unlike the historical reader.
Accepted cubes/frames reproduce old templates exactly under the old reader.
`config/qe_cube_order.toml` selects only the corrected input layout. Keep the
0.50 nm target, existing first-below-target rule, grid and normalization; no
height, support or isovalue search is selected by benchmark labels. This does
not establish an explicit tip, calibrated current or physical adequacy of the
remaining legacy surface policy. See the dated journal comparison.
Completed: **676/870 correct, 27/145 exact**, versus legacy **679/33**, coverage
**852** in both. Correct-order mean heights **0.255375/0.292328 nm** and support
**40/289, 58/289** expose a discontinuous remaining calibration before any
grade. No height/branch/support adjustment follows; retain the correct reader,
not a recognition promotion. `results/cube_order_20260924/report.md`.

The September 24 leave-target-out comparison changes only which lobes select
phase/mirror, not the physical state costs or calibration. It uses other lobes
in the same view without weighting, temperature, class-count prior or target
fallback. Exact ties retain first phase/mirror; no remaining evidence gives an
unavailable CC score. Geometry and overlapping patches remain shared, so this
is not independent cross-validation or noise calibration. Settings are explicit
in `config/mold_leave_one_out.toml`; no threshold is learned from the benchmark.
Completed: **670/870 correct, 32/145 exact**, versus tangent control **679/33**,
both coverage **852**. Reject this state-selection ablation without changing
physical calibration or the vote rule. `results/mold_loo_20260924/report.md`.

The adjacent-amplitude experiment uses only the previous/next fitted lobes,
without chemical labels, neighbor fitting or a signal-dependent selection rule.
Their amplitude derivatives are sampled on the same native grid and target
patch as the original tangent columns. `config/tangent_mold_neighbors.toml`
differs from `config/tangent_mold_finite.toml` only by this basis extension; both
explicitly retain the tested missing-cost omission. All numerical tolerances,
physical widths, N and classifier settings remain frozen. Projected contrast
energy diagnoses numerical loss, not retained experimental chemical information.
This does not repeat the failed wider-patch experiment or change calibration.
Completed: **671/870 correct, 29/145 exact**, versus target-only **679/33**, both
coverage **852**. Reject the adjacent-amplitude basis without a new cutoff or
strength adjustment. Retaining over 93.97% of template contrast energy did not
guarantee retained recognition. `results/tangent_neighbors_20260923/report.md`.

The September 23 mold-state comparison changes no physical calibration. It uses
saved tangent costs, omits both-infinite observations from geometric-state
selection without filling their scores, then separately tests a common
parity/mirror using fixed equal view weights. Chemical states remain independent
between views. These NCC costs are not likelihoods; sharing a state is a testable
hypothesis, not proof of identical acquisition contrast. No weight, missing-data
threshold, geometry, composition or vote rule may be tuned against labels.
Settings are explicit in `config/mold_state_comparison.toml`. Completed result:
omission preserves all predictions at **679/870 correct, 33/145 exact**, while
sharing gives **675/32**, both coverage **852**. Retain the opt-in correctness
fix, reject sharing without weight/threshold adjustment. Production defaults
remain unchanged; the complete comparison is available for human review at
`results/mold_states_20260923/report.md`.

The subsequent orientation test adds only the analytic ellipse-angle derivative
at fixed pixels/frame (`config/tangent_mold_orientation.toml`). Exactly isotropic
rows have no angle derivative; every nonzero derivative is normalized with the
same rank tolerance, without an ellipticity threshold. Near-isotropic behavior
is therefore a nuisance-span convention, not calibrated angular uncertainty.
No physical bound, fitted angle, patch support or classifier setting changes.
It gives **678/870 correct, 30/145 exact**, versus tangent **679/33** at the same
**852** coverage. Reject this orientation ablation without adjusting calibration
or introducing a near-isotropy cutoff. Neither result beats historical **677/36**
and coverage **854** on all criteria.

The September 23 tangent-mold experiment is not a new physical calibration.
`config/tangent_mold_projection.toml` fixes the affine/target-Gaussian derivative
span before grading, respecting circular versus elliptical widths. The two
`1e-12` tolerances concern numerical rank/annihilation, not chemical confidence.
Native smoothing/interpolation, finite-pixel requirements, physical templates
and classifier calibration remain unchanged. Projecting out shape errors can
also remove chemistry; retained contrast energy is diagnostic only. No rank,
support, threshold or tangent column may be chosen using benchmark labels.
The fixed comparison gives **679/870 correct, 33/145 exact**, against support
**676/34**, both with coverage **852/870**. It exceeds historical **677/36** only
per unit, not on exact chains or coverage (854); no calibration/default changes
follow this result. Full evidence is in the dated journal entry.

## Parameter classification

The former bootstrap's "measured / principled / free" table overstated what a
single STM image established. The useful distinction is now:

| Quantity | What the current evidence establishes |
|---|---|
| Pixel coordinates | Shared SXM header/grid convention, not independent instrument calibration |
| Profile width and spacing | Apparent, preprocessing- and peak-selection-dependent measurements, or explicitly unavailable |
| High-frequency MAD / image dispersion | Different descriptive scales, not automatically background-noise sigma |
| Gaussian `sigma_parallel_*` | FWHM/2.355 only for an appropriate isolated Gaussian; not justified by arbitrary local or overlapping peak spans |
| Fit width, support padding, spacing bounds, overlap | Physical/model choices requiring separate justification; ±30%, three repeats or a fixed padding multiplier are assumptions, not measurements |
| GCV, model form, optimizer budget | Explicit selection/numerical choices, not physically calibrated from the benchmark |

Historical auto-calibration experiments remain in the journal. Their old output
files are not rewritten and their externally graded results are not new
measurement validation. No automatic replacement config is generated here.

### Structured diagnostic policy

Todo 9's structured follow-up diagnostics are not calibrated from benchmark
labels. Their estimator semantics are frozen in the
`[selection.diagnostics.policy]` block of
`config/unit_assignment_structured_model.toml`; missing or differing keys
fail closed as `BLOCKED`. The policy fixes the fixed-ν=8 Student-t residual and
one-component fallback, date-centered unbalanced ICC(1,1), equal-feature
summed-squares pooling, the configured view contrast, explicit Hyndman–Fan
Type 7 quantiles, inclusive finite-sample tails, strict zero/equality handling,
and deterministic Holm ties.

Channel dropout is a separate follow-up statistic: it removes exactly
`bwd_neg_com_t` and `bwd_neg_diag45`, freshly refits `base_local` using only
the same inner-training partition, and rescored the untouched held-out rows.
It must never reuse or alias the original view contrast. The proposed
featurewise-standardized residual tests and hierarchical reliability model are
deferred to a separately preregistered v3 study and do not alter v2.

## Effective sample size — why GCV is the canonical criterion

The reported STM image/residual correlations are strong (ρ ≈ 0.9–0.95 at
lag 1; reported range 17–100 px). A narrow fit tube (~10 px across) may contain
too little independent spatial information to estimate a correlation model
reliably. Correlation measured over the molecular image also includes signal
structure: it is not automatically a background-noise covariance estimate.

The `n÷9` heuristic is a fixed placeholder, not a measured number of independent
observations. Replacing it with another ad hoc effective sample size can greatly
change absolute BIC/AICc values without improving the noise model. Their iid
absolute values must therefore not be interpreted as calibrated model evidence;
they remain secondary diagnostics/guards, not a replacement for canonical GCV.

**GCV** (`RSS·n/(n−p)²`) avoids inserting an arbitrary `n_eff` into the
per-candidate score and remains the canonical practical criterion. For a linear
smoother, generalized cross-validation approximates leave-one-out error by
replacing individual leverage corrections with their average; it is not the
exact leave-one-out identity in general. This code also uses a parameter-count
approximation in a nonlinear, constrained fit. Spatial correlation, active
bounds and model mismatch can therefore affect its predictive interpretation.
GCV does not by itself calibrate count uncertainty or establish chemical truth.
BIC/AICc remain secondary diagnostics/guards; the `n ÷ 9` placeholder is unchanged.
Numerically profiling amplitudes/background must not reduce the model parameter
count used in GCV.

## Calibrating a new molecule

1. Inspect raw scans and both acquisition directions with the diagnostic above.
2. Separate unavailable values, apparent measurements and independent physical
   evidence; do not take a median of silently substituted defaults.
3. Establish a model linking observed shape to the intended physical parameters,
   including overlap, tip/preprocessing effects and measurement uncertainty.
4. Record independently justified bounds in a separate explicit TOML, without
   benchmark labels. The diagnostic does not perform this step automatically.
5. Freeze that candidate before fitting and external grading; visual QC and
   repeated old benchmark scores do not establish unknown-chain chemistry.

Noise scale and correlation depend on acquisition, tip/feedback state and
preprocessing. They should not be assumed identical across sessions merely
because the instrument and flatten/smooth settings are unchanged. Estimates
from an image containing the molecule can also depend on its structure.
Recheck these quantities from suitable background or independent acquisition
evidence when transferring calibration, as well as molecule-specific FWHM,
spacing and `n_max`. This caution does not change the current production config.

## Structured evaluator-v1 is not physical calibration (correction3 pending review)

The evaluator config is a policy/static/synthetic correction3 prerequisite, not a
physical calibration and not an application or benchmark result. It becomes
authoritative only after parent acceptance, a fresh independent Oracle PASS, and
reviewer-owned `GateClosure`. Its
canonical path/hash/runtime/member checks must pass before any formula,
bootstrap, or graph work; Todo 13 remains blocked until that gate.

The exact graph reference is one partition/outer-score T12 block per held-out
date, with per-scan `logZ` summed only from that scan's selected blocks. T12
marginals, not Viterbi labels or report-wide sums, provide enabled probabilities;
disabled output uses selected unary `q`. Status rows distinguish `BLOCKED`,
`FAIL`, consumed `SKIPPED`, and `PASS`/`SKIPPED` final gates. Type 7 uses
`h=1+(n-1)p`, with 500-value lower quantiles exactly
`0.525*x_(13)+0.475*x_(14)`; sign masks use sorted dates, inclusive `>=`,
`count/2^K`, and no `+1`.

Evidence is limited to static policy checks and hand-computed/synthetic fixtures.
It creates no 10–20mer application claim and no benchmark claim. External labels
or grader outputs, if ever used, are permitted only after the independent policy
gate and outside this label-free prerequisite.

Correction2 distinguishes consumed diagnostic `SKIPPED` rows from terminal
statuses, resolves exactly one selected unary fit and matching T11/T12 partition
reference per held-out date, and computes entropy and view agreement over their
pooled frozen node/pair populations. Descriptor-based authority snapshots reject
path substitution, symlinks, hardlinks, identity collisions, and changed bytes
before any scientific work.

Correction3 integration binds worker-produced static and synthetic evidence; it
does not execute the absent Todo13 evaluator, use real data, use labels or
composition priors, or make a benchmark/10–20mer claim. The bound mutation
projection is 680 policy/semantic rows plus 163 authority rows (843 total),
not an integrator reimplementation. The validator uses the exact policy
bootstrap contract (500 Mersenne-Twister seeds, Type 7 lower quantile and 32
exhaustive sign masks), descriptor snapshots, final revalidation, and the
unchanged GCV/`n_eff` authority boundary. `n_eff` remains an authority-bound
policy statement rather than a newly calibrated physical quantity.

Correction4 is only publication-provenance remediation. It regenerates the
correction3 static/synthetic bytes twice and publishes them once with
descriptor-relative, exclusive no-replace creation. It introduces no physical
calibration, threshold, model, T8/T11/T12, GCV, `n_eff`, label, benchmark,
application, or Todo change; the live configuration remains a candidate and
Todo13 stays absent and blocked pending parent acceptance, a fresh Oracle, and
GateClosure.

Correction5 is a narrow administrative terminalization caveat. The six
correction4 canonical bytes remain valid by reference and are not republished;
the correction4 claim is non-authoritative because its predecessor hash and
finalizer/Boulder closure were defective. Correction5 records the failure and
cleans only the captured residue. It adds no physical calibration or science,
and changes no configuration, threshold, GCV, `n_eff`, T8/T11/T12, label,
benchmark, application, or Todo behavior. Todo13 remains blocked pending parent
acceptance, a fresh Oracle PASS, and reviewer-owned GateClosure.

Correction6 is an administrative publication boundary, not a calibration
change. Correction5 cleanup remains valid, but Correction5 is non-authoritative
because its evidence closure omitted the cleanup receipt and its terminal replay
omitted all six required per-path canonical bindings. Correction6 references the
existing Correction4 canonical bytes without republishing them. The parent owns
the external checkpoint and performs atomic no-replace directory publication;
publication is established only by the parent's external receipt after the
staged 19-file root and both manifest namespaces have been verified.

Correction6 never authorizes Todo13 itself. Todo13 remains blocked pending
parent acceptance of the external receipt, a fresh independent Oracle PASS, and
reviewer-owned `GateClosure`. No policy, configuration, calibration, threshold,
GCV, `n_eff`, T8, T11, T12, label, benchmark, application, or Todo behavior
changes.

## Reconstruction versus calibration (2026-09-16)

The user-approved Julia reconstruction has its own explicit numerical settings
in `config/unit_assignment_reconstructed.toml`; see the config and unit-assignment
references. The backward residual half-plane descriptor is a new definition,
not a recovered calibration or an inferred original formula. Its 9×9 grid and
normalization are fixed before comparison. The mold/Fisher settings port the
existing method; no benchmark labels, expected counts, or frozen predictions
are used to calibrate these quantities. The counting calibration is unchanged.

The September 20 opt-in `unit_assignment_transverse_fisher.toml` tests physical-u
instead of physical-t reflection when scoring the empirical Fisher mold. Its
`fisher_layout` is a coordinate convention, not a calibrated length, noise
estimate or fitted threshold. Patch sampling, normalization, disk support,
Fisher training, seeds and voting remain fixed. Neither reflection symmetry nor
a synthetic test alone establishes improved chemical recognition; the complete
external benchmark comparison is required. Original numerical settings are unchanged.
That fixed comparison regresses to 669/870 correct and 27/145 exact, versus
671/870 and 28/145 for the control at identical coverage; the transverse
variant is not retained. No recalibration or threshold adjustment follows.

The separate September 20 matched-residual comparison tests
`patch_residual_filter = "smooth_residual"`: `S(data - model)` rather than
`S(data) - model`. This makes an exact model image yield zero residual and
retains an injected residual as `S(signal)`, not as an unsmoothed signal.
It changes no smoothing radius, physical calibration, geometry, selected count,
normalization, classifier setting or missing-value treatment. It removes the
deterministic mismatch between smoothed data and unsmoothed model, not an
arbitrary background or transverse gradient. Synthetic consistency alone does
not establish better recognition; the full external comparison is required.
Both older assignment configs now state `smooth_data_only` explicitly, changing
their file hashes but not their effective settings.
The completed paired grade gives **672/849 (79.2%), 29 exact**, versus
671/849 and 28 for the control. This small development-benchmark gain supports
keeping the opt-in candidate, not a claim of calibrated chemistry or independent
generalization; historical 677/854 and 36 exact remains better. No threshold or
other calibration adjustment follows.

The subsequent two-descriptor comparison changes neither calibration nor
residual extraction. Its first moment uses the known transverse coordinate,
normalized by the existing patch half-width. Its affine-residual asymmetry uses
the orthogonal least-squares projection away from `[1,t,u]` on the full 9×9 grid.
This removes the affine component of molecular signal too: only the orthogonal
remainder is retained. It is not an estimated physical background, an acquisition
correction or a guarantee of preserving chemistry. Synthetic plane invariance
and agreement with an independent QR projection establish the algebra only.
All geometry, physical scales, missingness rules, classifier settings and voting
stay fixed; only the one descriptor supplied to GMM changes. Both variants are
declared before a single full-cohort comparison and external grade.
The completed comparison gives **672/849, 28 exact** for the first moment and
**673/849 (79.3%), 33 exact** for affine residuals, versus **672/849, 29 exact**
for the replayed control. Retain the affine config as a working candidate only;
history still has four more correct positions and three more exact chains.
The sixteen changed decisions all end at zero vote margin, and the benchmark
has been reused extensively. Neither this small gain nor the projection algebra
establishes a physical background correction or calibrated chemical confidence.

The independent signed-CC and affine-Fisher comparisons introduce no new
molecular calibration. Signed CC uses the existing physical template cost
ordering, not an experimental class or an oracle sign flip. Affine Fisher removes
the geometric `[1,t,u]` subspace from its 197-pixel disk at training and scoring;
the original unprojected center amplitude still names the learned groups.
Its explicit `fisher_projection_zero_l1 = 1e-12` numerical floor marks a vanished
projected signal unavailable, without missing-pixel imputation. This also removes
affine molecular signal. Templates, physical scales, mixture composition freedom
and count selection are unchanged. The completed September 21 comparison is
negative: signed CC gives **668/849 (78.7%), 32 exact** and affine Fisher
**672/849 (79.2%), 32 exact**, against the exactly replayed **673/849, 33 exact**
control, all at 849/870 coverage. Neither replaces the affine-descriptor working
candidate. Some incorrect new decisions have voting margin one; this is not
evidence of calibrated chemical confidence. No post-grade sign flip, threshold,
seed or hybrid is fitted to the benchmark.


The September 21 numerical comparison also introduces no molecular calibration.
`fisher_score_center` explicitly chooses the historical near-zero centered mean
or the original opposite-fold training mean; it does not impose equal cluster
populations. `gmm_final_covariance` chooses historical `ridge` or analytical
`ledoit_wolf` shrinkage on the final hard covariance only. The explicit
`gmm_covariance_ridge = 1e-6` is the unchanged numerical floor, not a measured
noise variance. No direction, sign, composition, threshold or seed is calibrated
from truth. Whole-scan withdrawal diagnoses sensitivity to training data, not
held-out recognition; Fisher features are not retrained within those halves.
Neither covariance conditioning nor stability calibrates chemical confidence.
The completed comparison improves **673 → 675 correct / 870**, with unchanged
**33/145 exact** and 849/870 coverage, using training-mean centering. This is the
opt-in working candidate at that stage, still below historical 677/36. Shrinkage gives
**665 correct / 34 exact** and is not retained as the replacement: its median
full-cohort covariance condition number falls from about 23,138 to 191, but
recognition loses eight correct positions. Whole-scan withdrawal changes 370/893
GMM decisions for ridge versus 350/893 for shrinkage: substantial sensitivity
remains despite modest improvement of that diagnostic. Better conditioning and
stability are not a substitute for the external recognition grade. No combined
variant, threshold, seed or parameter is tuned afterward.

### Patch support and covariance-volume score (2026-09-21)

`assignment_patch_support = "complete_disk_symmetric"` changes availability, not
the finite-pixel median/sample-standard-deviation patch normalization. The Fisher
disk must be completely observed. The affine backward descriptor additionally
uses reflection-closed observed support around a complete central disk; discarded
corner partners are not estimated or filled. Complete patches retain identical
arithmetic. Newly usable rows can change learned Fisher/GMM features throughout
the cohort, so this is a full-cohort experiment, not a local repair of only the
previously missing predictions. Missing disk pixels still require abstention.

`gmm_final_score = "gaussian_density"` adds the covariance-volume normalization
only at final assignment. It does not refit means, covariances or mixture weights,
nor change the two distance-based hard self-training iterations. The raw-amplitude
cluster naming convention and hard seed-vote aggregation remain, but their
resulting assignments may change. Dependent lobes, per-file standardization and
learned clusters mean that neither normalized density nor vote margin establishes
chemical confidence. No benchmark-derived sign, class ratio or threshold enters
either candidate. These are two independent, uncombined alternatives to 675/33;
no normalization or scan-weighting experiment is included.

The completed comparison gives **676/870 correct and 34/145 exact** for support,
versus 675/33: three recovered correct predictions are offset by two regressions
among previously available rows. Coverage rises **849 → 852/870**, while
classified accuracy falls **79.5% → 79.3%**. The gain on the fixed denominator,
exact chains and coverage supports retaining it as the latest opt-in candidate,
not as a calibrated-confidence improvement or a new champion. Historical 677/36
still leads. The Gaussian-volume score gives **666/870 and 32/145**, with all
19 changed final decisions at zero vote margin, and is rejected. No threshold,
class count, per-file choice or parameter is adjusted after these results.

### Separate complete training from partial prediction (2026-09-21)

The sole new `[selection] assignment_training_support = "complete_patches"`
mode restricts Fisher/GMM learning by observed patch completeness, never by
predicted type or external correctness. Fisher requires its full forward square;
GMM requires all three patch families it consumes. The support candidate's
disk/symmetry rule still governs whether partial patches can be scored. This
does not change extraction, patch normalization, physical calibration, selected
N, DFT molds, vote thresholds or k-means.

The GMM's per-file feature-normalization moments, fits, weights and raw-amplitude
group naming are estimated from complete rows only, then frozen. Moments still
use finite values separately per feature; there is no robust estimator, global
fallback or scan reweighting. A file with no finite complete training value for
a feature cannot provide that view. Opposite-lobe-parity Fisher scoring remains;
this is not a new held-out-scan validation design. Complete training rows can
still be noisy and fewer rows can destabilize learning, so improvement is a
hypothesis, not a guarantee. The full 146-file output and all regressions must
be checked before the separate 145-file external grade is considered.

The completed test is **negative: 674/870 correct, 33/145 exact**, versus
676/34, with coverage falling **852 → 850/870**. All four fitted lobes of
`240818_019` have partial squares, so that file has no complete-row moments;
the prescribed abstention loses its two previously usable predictions. Elsewhere,
one gain is offset by one loss, including loss of exactness for `240818_020`.
This is evidence against this complete-training variant as a replacement, not
permission to add a benchmark-directed fallback. Keep the 676/34 support
candidate; normalization rules, counts and the frozen application stay unchanged.

### Robust GMM feature normalization (2026-09-21)

The newly authorized comparison changes only GMM per-file feature moments from
mean/sample-std to median/IQR, before pairwise products. This is not physical
pixel calibration and does not affect Fisher or k-means. The support candidate
at 676/34 is the control; complete-case training is not reused. All finite
admissible values remain eligible, including usable partial patches.

`[preprocessing] gmm_feature_normalization` selects the explicit formula, and
`gmm_scale_fallback = 1.0` records the degenerate-scale guard in every native
config. Q25/Q75 use Type 7 interpolation, with no normal-consistency rescaling.
An empty feature support stays unavailable; a zero/nonfinite IQR uses one.
Positive small IQRs are not floored. Median/IQR resists changes in extreme
observations when estimating moments, but does not bound their transformed
values; short scans or interactions can therefore become less stable, not more.
The completed comparison is negative: **667/870 correct, 667/852 (78.3%),
10/145 exact**, versus support control 676/34 at unchanged coverage. It loses
24 exact chains and gains none. No IQR fallback occurs in 1,168 scan/feature
pairs, but the expanded maximum grows from 5.44 to 335.47. This observed tail
amplification does not, by itself, establish the cause of the recognition loss.
Retain the mean/std support candidate; no clipping or alternative quantiles are
introduced after the grade.

All choices precede grading. No expected count, sequence, composition, benchmark
label, class balancing, scan weighting or benchmark-tuned fallback enters this
experiment. Counting, DFT and the frozen unknown-chain application are untouched.

### Equal-scan GMM weighting (2026-09-21)

The next authorized comparison leaves physical calibration and normalization
unchanged. It tests whether giving each scan equal total training influence
helps the 676/34 support candidate. `[selection] gmm_training_weighting` selects
`equal_lobes` (existing behavior) or `equal_scans` (one opt-in candidate).
Weights are proportional to inverse usable training-row count, normalized to
mean one to retain the objective's numerical scale. They affect GMM learning
and physical amplitude naming consistently, not the separate Fisher/k-means
heads or final vote. Within-scan finite feature moments stay unchanged.

This is not class balancing, an expected-N prior, or an effective-sample-size
estimate. Both chemical populations remain free. A scan with few usable lobes
can gain influence even if noisy; equal influence is a hypothesis, not a
guarantee of correct chemistry. The fixed comparison is one control
and one variant on Viper with a 30-minute ceiling. No count refit, bootstrap,
DFT change, unknown-chain rerun or post-grade adjustment is included.

The completed result is mixed: **675/870 correct, 675/852 (79.2%), 36/145 exact**,
versus the replayed support control **676/852, 34 exact**. It gains two exact
chains but loses one correct position; all seven changed decisions end at zero
vote margin. Retain the 676/34 support candidate as the primary reference, with
historical 677/36 still unexceeded. This result does not establish that short
or noisy scans caused the loss, nor justify a new physical calibration or
post-grade weight adjustment. The reused benchmark is development evidence.

### Continuous GMM seed aggregation (2026-09-21)

One authorized comparison retains the GMM's normalized per-seed scores instead
of replacing each by an argmax decision before averaging seeds. The explicit
`selection.gmm_seed_aggregation` modes are `hard_vote` (legacy) and
`mean_membership` (one opt-in candidate). Only aggregation changes, not the
learned groups, their amplitude-based physical names or the final threshold.
There is no temperature, label-fitted probability calibration or mixture-size
constraint. Scores may already be saturated, and smooth votes need not improve
chemical correctness. Their confidence remains an uncalibrated vote margin.

Physical calibration, selected N, missing-data handling, eight-decimal output
precision and the two-head vote stay fixed. Scope: support control plus one
candidate, one Viper allocation capped at 30 minutes; no scan weighting,
resampling, new naming rule, count refit, DFT change, unknown-chain rerun or
post-grade adjustment.

The completed result is **negative: 671/870 correct, 671/852 (78.8%), 26/145
exact**, versus the exactly replayed 676/34 support control at the same coverage.
All 24 changed decisions leave old exact ties (GMM=1, k-means=0) and become 0
as the continuous GMM score drops below 1. The GMM still emits class 1 for these
24 rows; the loss arises in the unchanged final mean-of-two-heads decision.
One exact chain is gained and nine lost. Reject this variant, retain the 676/34
support config and leave calibration, precision and thresholds unchanged.
Historical 677/36 remains the target. This is evidence about this specific
uncalibrated fusion, not a general rejection of continuous scores or independent
validation on unknown chains.

### Whole-scan GMM bootstrap (2026-09-21)

One fixed comparison resamples **whole scans**, not individual lobes or chemical
classes. Twenty replicates each draw S usable scans with replacement, using
seeds 0–19; all admissible rows of a drawn scan travel together. The ten existing
initialization seeds are reused within each replicate. Row multiplicity enters
learning and the unchanged raw-amplitude naming rule, not an imposed class prior.
Per-scan scaling and physical calibration remain unchanged. Unsampled scans
still receive predictions, so this is not held-out validation.

The explicit settings are `selection.gmm_resampling="whole_scans"`,
`gmm_bootstrap_replicates=20` and `gmm_bootstrap_seed=0` in the separately named
opt-in config. Earlier configs retain `none`, `1`, `0`. The candidate averages
hard seed votes within each replicate and then valid replicate means equally;
it does not select a best seed or tune confidence. Scores remain uncalibrated
frequencies. Final threshold, precision, N, unavailable-row handling, k-means,
Fisher, CC and all upstream observations remain fixed.

Scope is one support control and one candidate, one Viper job limited to one
hour, four requested CPUs and 16 GB. No retry, combined variant, new naming
rule, count refit, DFT change, unknown-chain rerun or post-grade adjustment.
Job **11925188** completes **0:0 in 7m17s**. All 200 seed fits are accepted;
bagging gives **665/870 correct, 665/852 (78.1%), 10/145 exact**, versus the
exactly replayed **676/852, 34 exact** support control. Coverage is unchanged;
24 exact chains are lost and none gained. Of 62 final 1→0 changes, 57 leave
old GMM=1 / k-means=0 ties under the unchanged fusion. Reject the variant;
retain the **676/34 support candidate**, still below historical 677/36. No
post-grade change to fusion, seeds, replicate count, naming, precision or
threshold follows. This negative result concerns this specified ensemble and
fusion, not all bagging methods; the reused grade is not independent validation.

### Tied GMM covariance throughout learning (2026-09-22)

The user authorizes one comparison of the 676/34 support control with a GMM
whose two components share covariance throughout learning. This is not the
previously rejected final-only Ledoit-Wolf estimate or final Gaussian-volume
score. The same features, per-scan scaling, interactions, ten seeds, two hard
updates, physical group-naming rule and final vote stay fixed. N, Fisher, CC,
support, acquisition preprocessing and unknown25 are not changed.

`model.gmm_covariance_structure="tied"` pools within-component scatter using
the learned responsibilities in every EM update and hard indicators during
initialization/self-training; divide by the number of usable rows and add the
unchanged ridge. Means and mixture masses stay free, with no class-count prior.
The covariance has 666 independent entries in the 36-dimensional view instead
of 1,332 across two separate matrices. This reduces model freedom, not a claim
of calibrated noise or chemical certainty. Better numerical conditioning alone
would not establish better recognition.

One Viper job is capped at 30 minutes, four requested CPUs and 16 GB, without
automatic retry. Exact control replay and all-key/N/upstream/vote checks precede
external grading. Labels, expected N, benchmark membership and historical
predictions are not learning or inference inputs. No naming, grouped-Fisher,
fusion, threshold, resampling or post-grade search accompanies this test.
The completed result is **negative: 666/870 correct, 666/852 (78.2%), 6/145
exact**, versus the exactly replayed 676/34 support control at the same
coverage. No exact chain is gained and 28 are lost. Reducing covariance freedom
does not improve recognition in this fixed representation/fusion. Learned
smaller-component masses are about 13.2–14.4%, not an imposed class proportion;
final class-1 emissions decrease from 223 to 156 across all 900 rows. These
post-hoc observations do not authorize a composition prior, new naming rule or
threshold adjustment. Reject this variant and retain **676/34**; historical
677/36 remains the target. This reused grade is not independent validation.

### Whole-scan Fisher and relative GMM naming (2026-09-22)

Two independent tests retain the 676/34 support control and all physical/count
calibration. The first changes only Fisher's cross-fitting unit: entire scans
are kept in one of two deterministic SHA-256-ranked groups, with an explicit
split seed zero. Every fitted transformation, center and naming operation
uses only the opposite group. Invalid rows do not choose the grouping, and
neither chemical labels nor a required number of lobes is used. This tests
sensitivity to shared scan structure; it is not proof that the parity feature
was chemically wrong or that scans are independent molecules.

The second changes only the global names of the already fitted GMM components.
`gmm_cluster_naming="within_scan_z"` centers/scales raw amplitudes within scans
using eligible, feature-valid training rows; the higher global group mean names
class 1. The existing explicit `gmm_scale_fallback` applies to constant or
singleton training scans. Positive scan gains/offsets then cancel for nonzero
scales, but no chemical accuracy or invariance under arbitrary tip changes is
claimed. No class quota, composition prior, per-scan forced group or DFT anchor
is introduced; the underlying higher-amplitude physical convention remains an
assumption. GMM means, covariances, memberships and free masses stay unchanged.

The two candidates are not combined. Labels belong only to the external grade
after fixed outputs; one Viper job is capped at 30 minutes with no retry.
The completed comparison gives **668/870 correct, 668/852 (78.4%), 34/145 exact**
for whole-scan Fisher: eight fewer correct positions than the replayed support
control, with four exact chains gained and four lost. Both disjoint 73-scan
fits converge. All 26 final flips enter existing zero-margin vote ties; no
threshold or split-seed adjustment follows. Relative naming gives **676/870,
676/852 (79.3%), 34/145**, with unchanged fitted parameters, component names,
scores and decisions across all ten seeds. These naming conventions agree on
this cohort, not necessarily under other acquisition conditions. Retain
**676/34**; historical 677/36 remains the target. This negative recognition
result does not establish that parity cross-fitting is independent validation,
and the reused benchmark remains development evidence.

### Fixed alternative mixture learning (2026-09-22)

Two independently tested profiles retain the support control's 8 features,
28 pairwise products, per-scan scaling, amplitude naming, seeds, selected N,
Fisher and final vote. Rank four (`gmm_factor_rank`) and Student df five
(`gmm_student_df`) are fixed before real inference, with no benchmark-driven
choice. The iteration, tolerance, minimum-mass and numerical-guard parameters
are explicit in the [configuration reference](config.md); Gaussian learning
is unchanged. Neither alternative estimates a chemical class proportion.

The factor-analyzer model initializes each component spectrally, using the
mean of the bottom `p-rank` covariance eigenvalues for isotropic residual
noise. EM jointly updates the component center and loadings using augmented
latent moments; separate diagonal noise is floored at the existing ridge.
This follows [Ghahramani and Hinton's MFA equations](https://www.cs.toronto.edu/~hinton/absps/tr-96-1.pdf),
including their component-specific-noise option. The two pure-Mahalanobis
hard steps each retain one such conditional parameter update.

For Student learning, latent precisions are `(df+p)/(df+distance²)`.
Centers use responsibility-times-precision weights, while each scale scatter
is divided by the **responsibility mass**, not the precision-weighted mass;
the existing ridge is then added. See
[Peel and McLachlan](https://people.smp.uq.edu.au/GeoffMcLachlan/pm_sc00.pdf).
The same update is retained during hard reassignment. The stored scale is
not multiplied by `df/(df-2)`. Final scoring deliberately remains the
control's log-mass minus half squared distance, not a Student density.
Thus neither scores nor seed votes are calibrated chemical probabilities.

No preprocessing, DFT calibration, count selection, abstention threshold or
unknown-chain claim changes. The repeatedly reused benchmark can compare
these frozen choices, but cannot provide independent validation.

The completed comparison is negative: factors **620/870 correct, 23/145
exact**, Student **631/870, 25/145**, versus exactly replayed **676/34** support,
at **852/870 coverage** throughout. All ten Student fits converge in 65–86
updates; all factor fits reach 200 updates without satisfying tolerance.
Their smallest noise variance is about 7.59e-5, above the 1e-6 floor, so no
floor saturation or collapsed seed is observed. The capped MFA result does
not characterize a fully converged optimum. False positives increase more
than false negatives decrease; almost every changed final decision enters a
zero-margin vote tie. Retain support without post-grade cap/rank/df,
threshold, naming or fusion adjustment. Source **e2205cb**, job **11935072**;
full losses: `results/factor_student_mixtures_20260922/report.md`.

### Covariance regularization during learning and Student density (2026-09-22)

`unit_assignment_em_shrinkage.toml` applies spherical covariance shrinkage at
initialization, each Gaussian EM M-step and both hard updates. For centered
vectors `y_i`, responsibility weights `a_i = r_i / sum(r)`, define
`S = sum(a_i y_i y_i')`, `T = tr(S)/p I`, and
`b = sum(a_i^2 ||y_i y_i' - S||_F^2)`. The fixed coefficient is
`alpha = clamp(b / ||S-T||_F^2, 0, 1)`, or zero if `S == T`. Use
`(1-alpha) S + alpha T + 1e-6 I`. Hard memberships reduce to the previously
tested equal-weight formula. Means and free mixture masses use unchanged
updates; neither group is forced to contain a specified number of lobes.

This is a **declared fixed-weight plug-in extension**, not a reproduced
implementation or optimality claim for
[Halbe, Bortman and Aladjem's regularized GMM](https://cris.bgu.ac.il/en/publications/regularized-mixture-density-estimation-with-an-analytical-setting/).
The latter motivates regularization during learning. The extension is checked
against direct weighted outer products and the hard-membership limit.
Learned responsibilities, dependent lobes and estimated centering invalidate
any automatic interpretation as independent measurement-error variance. No
`n_eff`, physical noise estimate or chemically calibrated confidence is inferred.
Shrinkage need not increase unpenalized likelihood monotonically; decreases
and bounded-fit convergence are reported without filtering seeds by grade.

`unit_assignment_student_density.toml` retains fixed df `nu=5` and the
[Student latent-precision updates](https://people.smp.uq.edu.au/GeoffMcLachlan/pm_sc00.pdf),
but uses the same log mass-density in the E-step, both hard assignments and
final scoring: `log(pi_c) - log(det(scale_c + guard I))/2 -
(nu+p)/2 * log1p(distance_c^2/nu)`. Only the common fixed-dimension/fixed-df
normalizing constant is omitted. Scale is not converted to covariance.
Hard assignment now maximizes this score, rather than minimizing distance
alone; the parameter update still uses responsibility-times-precision means
and responsibility-mass scale denominators. This is a classification-style
hard-update extension of the soft Student mixture, not pure soft EM throughout.

Both arms preserve raw-amplitude naming, ten hard seed votes, k-means fusion,
thresholds, support and selected N. Internal model memberships and final vote
fractions are not calibrated chemical probabilities. Benchmark labels are
used only in external grading.

Completed outcome: **637/870 correct, 24/145 exact** for all-update shrinkage,
**627/870, 25/145** for coherent Student, versus exactly replayed **676/34**
support, all at **852/870 coverage**. All ten fits per arm meet tolerance:
25 or 43 updates for shrinkage, 65–86 for Student. Shrinkage coefficients
range about 0.02993–0.06040; seeds 1 and 3 each show eleven ordinary-likelihood
decreases, with maximum drop 1.14377. None is excluded or retried. Student
latent precisions range about 0.003045–5.447092. Numerical convergence does
not establish good recognition or a global optimum. Reject these two fixed
versions, retain support and leave the historical 677/36 target open. No
post-grade tuning or unknown-chain claim follows. Source **bb6943b**, job
**11936229**, full report `results/em_shrinkage_student_density_20260922/report.md`.

### Fresh raw-GCV counts with unchanged calibration (2026-09-22)

Repeating the original direct-extractor GCV sweep after symmetric filtering
gives **661/870 correct, 32/145 exact**, versus fresh fixed-N control
**675/34** and saved support **676/34**. Exact counts fall **106 → 101/145**;
coverage falls **852 → 849/870**. No physical bound, GCV/guard, `n_eff`,
threshold, optimizer limit or DFT source changes. This is not the promoted
batch hybrid counting lineage or an exhaustive sweep.

The fresh control changes three decisions despite identical N/settings.
Same-N base/split centers vary by up to 0.00845/0.07637 nm between fresh arms;
their maximum relative GCV differences are 0.115/1.049%. The existing
time-limited global search is a plausible, unproven source of variability,
not evidence for a new calibrated noise model. Whole-cohort learning also
couples predictions across scans. Reject this version without retuning;
the result is end-to-end development evidence, not an isolated causal effect
of N or independent validation. See `results/gcv_reselection_20260922/report.md`.

## Opt-in label-free exploration (2026-09-18)

The September 23 cross-view residual experiment profiles only amplitudes and
tilted backgrounds independently on the two acquisitions. Original fused
geometry, support, q05 offset and amplitude bounds remain shared. Profile tables
change subtraction only: raw patches, windows, main amplitudes and split features
are fixed. Same-view and opposite-view subtraction use the exact same profile
pair. KKT convergence and native validity of their fused mean are required;
individual-view errors are diagnostics, not a new calibrated noise or chemical
validity criterion. This is conditional cross subtraction, **not independent
cross-validation**: both views already informed geometry/support/bounds and
shared chemical contrast can be fitted away in both. No new physical threshold
or class prior; settings are in `config/cross_view_residual.toml`. The completed
test gives **675/870 correct, 30/145 exact** in same-view and **673/30** in
cross-view, versus saved **676/34**, all at **852/870 coverage**. All 292 solves
converge and 146 fused means pass validity, but neither variant improves
recognition. Retain support without calibration changes; journal and
`results/cross_view_residual_20260923/report.md` retain all failures and limitations.

The September 23 local-Gaussian follow-up changes the **fitted base model**, not
the patch frame. A quadratic tangent is derived at each evaluation from current
fitted centers, with no independently optimized angles. Model axes are exported
for both residual subtractions; the global sampling frame and split cache remain
unchanged. Shared full-precision starts and equal 300-iteration LM limits compare
global/local ellipses; only physically valid minimum-GCV models are eligible,
with the original full parameter counts. Existing spacing/width/amplitude boxes,
support, residual guard and orientation-independent overlap/κ proxies are retained,
not reinterpreted as exact anisotropic overlaps. Circular Gaussians keep the
exact global arithmetic. No degree/angle calibration from benchmark labels is
allowed. Settings: `config/local_gaussian_orientation.toml`; see the journal.
The completed test gives **664/870 correct, 30/145 exact**, versus matched
global **674/31** and saved **676/34**, with unchanged **852/870 coverage**.
All 584 fits are natively valid but iteration-limited, not converged optima.
Lower GCV on 87/146 scans is not improved chemical recognition. Reject this
version without changing bounds, degree, angles, iteration limits or thresholds;
retain saved support. See `results/local_gaussian_20260923/report.md`.

The September 23 `local_patch_orientation.toml` experiment derives a quadratic
centerline tangent from frozen Gaussian centers and rotates only patch sampling.
It does **not** calibrate or refit molecular orientation: the Gaussian residual
model retains its original global axis, N and all saved coefficients. Frame
rotation preserves the rounded global-axis norm and applies identically to both
scan directions. Existing finite-patch support can change with window coverage;
it must be reported, not repaired through labels or adjusted abstention. The
global-axis profile remains the default: the complete comparison gives
**670/870 correct, 30/145 exact**, versus global replay **676/34**, with unchanged
**852/870 coverage**. This rejects patch-frame rotation alone, not a consistently
refitted locally oriented molecular model. No post-grade angle/degree adjustment
or calibration change follows; `results/local_orientation_20260923/report.md`.

The September 23 `directional_response.toml` pilot tests relative scan-response
phase on five acquisition-order scans. Equal forward/backward line times are
required for the same-time-constant, mirrored first-order hypothesis. The
normalized relative transfer has unit spectral magnitude before the common
bilinear translation; its mean lag is removed to separate shape from displacement.
It cannot calibrate absolute/common blur, identify a full controller, or recover
latent molecular amplitudes. Its **0–4 pixel** lag range is exploratory.

Training and checking use buffered scan-line bands and shared observed sources.
Native preprocessing/base shifts still use the full image, and within-row
means are conditioned on in the heldout score: this is conditional predictive
evidence, not independent validation or a noise covariance estimate. A reverse
response direction is a matched negative control. No molecular geometry, N,
unit assignment, physical calibration or champion parameter changes follow
automatically, even from a positive image-only result.

The completed panel supports **no correction (0/5)**. Four scans select the
zero-response control; the fifth's tiny full-image gain is not stable across
folds. Reverse-sign fits saturate the 4-pixel limit and are not a calibration.
No physical parameter, direction convention or threshold is changed. See
`results/directional_response_20260923/report.md`.

The September 23 `image_registration_2d.toml` follow-up estimates displacement
from observed images without a molecular shape, count or template. Its fixed
**±0.16 nm** residual window is an exploratory hypothesis, not a newly measured
vertical calibration. All candidate interpolation sources share an observed
mask; each training fold also excludes every possible heldout backward source.
Global/band/fold agreement and peak separation are engineering checks, not
confidence intervals. Native full-image preprocessing and saved molecular
initialization prevent claims of independent validation. Image agreement and
fit improvement are not recognition improvement; no production calibration,
noise, residual threshold or class prior changes in this pilot.

The measured full/fold peaks are interior but ambiguous and disagree by up
to **2 px**. The previous model-preferred (-1,-1) px decreases direct image
correlation on all three supports. This does not identify a physical correction
or distinguish model compensation from local view differences. No calibration
is accepted and the conditional real fit is not run; report:
`results/image_registration_2d_20260923/report.md`.

The conditional `image_registered_fit.toml` arm holds the image-estimated
translation fixed, alongside zero/free controls at unchanged physical bounds.
Its parameter-count GCV includes both calibration coordinates; this is a
plug-in diagnostic, not an unbiased effective-degrees-of-freedom derivation.
It is synthetically tested but **not run on real data**, because the direct
image estimates fail the declared identifiability/consistency checks.

The September 23 `paired_shift.toml` pilot retains native N, geometry bounds,
amplitude scale, observed pixels and noise. Four deterministic starts test the
local minimum before a nested backward molecular translation, limited to one
pixel in either direction. The forward view defines the molecular frame; view
planes remain in observed coordinates to avoid trading a shifted plane for an
intercept. At zero translation the previous paired model is reproduced exactly.
Validity uses the predicted **observed-view mean**, not an unshifted latent
molecule against shifted data; the original mean/view residual limit is intact.
GCV counts all acquisition parameters. Block-heldout scoring is conditional on
the already saved N, support and native initialization, not independent
validation; correlated noise and shared initialization prohibit that claim.
No parameter, noise, threshold or class composition is calibrated to labels.

Completed result: paired RSS and all reserved-block scores improve, but all
selected shifts hit **−1 pixel in both axes** and full paired validity still
fails. This is evidence of a direction worth diagnosing, not an identified
registration correction or grounds for expanding bounds after observation.
Native calibration and recognition reference **676/34** remain unchanged;
`results/paired_shift_20260923/report.md`.

The September 23 `paired_solver.toml` experiment compares fresh native LM with
SLSQP on the same fixed one-scan mean/paired objectives and native starts.
An anchored affine unit-box parameterization and constant initial-RSS scaling
change numerical coordinates, not physical boxes or the minimizer set. All
molecular/acquisition coefficients remain jointly fitted; no plane/amplitude
profiling or physical calibration is added. Numerical budgets and tolerances
are explicit in [configuration](config.md). Native validity and the previous
independent stationarity audit are retained, separately from solver stop codes.
No recognition or chemical accuracy claim follows from solver convergence.

SLSQP reduces the projected gradients by **24,000–440,000-fold**, without
passing the frozen strict stationarity audit or adding a fully valid family.
Gaussian elliptical maxima become **3.559 fused / 3.443 paired mean**, with
paired views **4.291/3.542**, still failing 3.5. The result supports numerical
progress, not relaxed physical calibration or recognition improvement. Report:
`results/paired_solver_20260923/report.md`.

The September 23 `paired_convergence.toml` follow-up changes only the local
iteration ceiling (300 versus 10,000), retaining all saved starts, objective
pixels, physical bounds and native stopping tolerances on the single blocked
scan. Each fit has a 600-second internal LM limit; the one Viper job is capped
at one hour. An independent two-step finite-difference audit measures
box-constrained stationarity after scaling parameters to the unit box and RSS
by its initial value. Its explicit numerical tolerances are documented in
[configuration](config.md), never used to relax validity, choose N or infer
chemistry. Neither a native small-step stop nor a stationary point proves a
global optimum. No physical calibration or recognition claim follows merely
from reduced residuals; this experiment does not grade the benchmark.

All eight long fits hit 10,000 iterations and remain nonstationary. Although
paired Gaussian mean validity now passes (**3.491 < 3.5**), forward/backward
maxima **4.349/3.575** still fail. No complete comparison or calibration change
follows; see `results/paired_convergence_20260923/report.md`.

The September 22 paired-acquisition feasibility experiment preserves the native
Gaussian/split mean shape and all molecular bounds at saved N. It adds
equal/opposite view gain and background-plane terms, fixing mean gain to one
to remove the trivial gain/amplitude ambiguity. `paired_acquisition.toml` bounds
gains to `[0.5,1.5]`, the differential intercept to `[-5,5] nm` and differential
tilts to `[-1,1]`; these are exploratory numerical boxes, not a new physical
calibration. The matched fused and paired continuations share their start and
300-iteration budget. Native mean validity is never relaxed; paired residuals
must additionally pass the unchanged guard in both acquisitions. Stacked GCV
is a heuristic with four added coefficients, not an independence/noise claim.
The first real check is confined to the previously failed scan, without labels.

That check does not clear native validity: Gaussian elliptical mean maximum/
noise is **3.651 fused / 3.518 paired**, above 3.5, with capped nonconverged
continuations. Paired split passes mean validity but fails the forward-view
guard. This is not proof of a wrong N or inadequate physical widths; numerical
and model limitations remain unresolved. No calibration is changed and no
full-cohort recognition grade follows. Report:
`results/paired_acquisition_20260922/report.md`.

The newly authorized `registered_refit_original_support.toml` follow-up freezes
the original unregistered native ROI, axis/origin, tube and axial bounds for
both refit arms. It changes no physical constant and does not expand support
to accommodate saved N. The legacy preprocessing's imputation influence remains
in this original support; only actually observed samples enter either new
objective. The offset is estimated on the observed original ROI. Frame/cache
consistency is checked, with a saved-only replay against prior native geometry
before grading. No new N, noise estimate or independence claim. The three
previous span failures disappear, but `240817_006.sxm` fails the native Gaussian
residual guard after registration; 145/146 scans complete all fits. No partial
grade is produced, and no guard is relaxed. All original supports replay
exactly; retain support **676/34**, with recognition still inconclusive. See
`results/registered_original_support_20260922/report.md`.

The authorized September 22 registered-refit follow-up keeps physical bounds,
registration gates, selected N and assignment settings unchanged. It refits
both Gaussian and split profiles with native optimizers, selecting their
family by GCV at the saved N; it does not reselect N. The observed unsmoothed
mean supplies the fit, with native finite-only ROI statistics and no fabricated
fit pixels. Nearest-observed sampling is initialization-only. Masks cannot undo
earlier native flattening's imputation influence. The zero-shift arm follows
the identical rules and its fits are reused exactly when registration is zero.
ROIs/axes/support may differ for nonzero shifts; lower RSS across these arms is
not itself an improvement in a common objective. Native noise estimates and
the `n_eff` placeholder are unchanged; no independence or convergence claim is
made. Settings: `registered_refit.toml`. The result is **inconclusive**: two
zero-shift controls and one registered scan cannot fit their saved N within
the new observed support under the unchanged conservative parametrization.
143/146 scans complete all fits. This does not prove incorrect counts or worse
recognition; no candidate grade is available. No bound, support parameter or N
is adjusted afterward. Retain support **676/34**; full diagnostic:
`results/registered_refit_20260922/report.md`.

The September 22 acquisition-translation experiment changes neither molecular
calibration nor cached N/geometry. `acquisition_registration.toml` declares all
search and acceptance settings: up to 2 nm / one quarter width in x, four
contiguous bands, signed correlation ≥0.60, distant-peak gap ≥0.01 and 0.08 nm
peak-neighborhood/band-agreement distances. At least 256 common pixels / eight
rows per band and 16 pixels per row are required. These are experimental
identifiability gates, not a noise estimate or chemical validation. Forward is
the fixed reference; only backward observations move before subtracting the
saved Gaussian model. Native flattening stays unchanged; restored raw masks
cannot undo earlier imputation influence on that background. A separate
zero-shift mask control isolates this effect. No N, optimizer, GCV, `n_eff`,
class threshold, composition prior or missing-head fallback changes.
The completed comparison is **negative: 671/870 correct, 26/145 exact**, versus
both controls **676/34**, coverage **850 versus 852/870**. Stronger local
image agreement does not validate this assignment change. Keep saved support,
without post-grade calibration changes; report:
`results/acquisition_registration_20260922/report.md`.

The subsequent September 22 geometric-profile comparison keeps N, physical
calibration, preprocessing and assignment settings fixed. Joint and profiled
LN_BOBYQA arms share one fresh native initialization chosen by valid GCV at the
saved N. Both use the native κ-penalized RSS and native validity checks, unlike
the old raw-RSS-only prototype. Equal evaluation/time ceilings are not equal
computational cost and do not prove convergence; actual stops and costs are
saved. Split skew remains cached. No new noise covariance, effective sample
size, composition prior or class threshold is introduced. The completed external
comparison is **negative: 632/870 correct, 20/145 exact**, versus joint and saved
support **676/34**, at unchanged **852/870 coverage**. All 146 profiled pixel
objectives improve but recognition worsens; 103 searches are evaluation-capped.
Do not replace the saved support geometry or retune calibration from this result.
Settings remain in `geometry_profile_comparison.toml`; full report:
`results/geometry_profile_20260922/report.md`.

The authorized September 22 fixed-geometry amplitude comparison changes no
physical calibration. It solves only the Gaussian amplitudes and tilted plane
inside their native finite bounds, using the existing linear-profile solver.
The literal saved decimal geometry and N stay fixed; the split shape is reused.
GCV keeps every original model parameter and is diagnostic only. A lower RSS
or satisfied KKT condition does not establish better chemical recognition.
The completed comparison demonstrates that distinction: RSS decreases on all
146 scans (median **0.681%**), but recognition falls to **675/870 correct,
32/145 exact**, versus exactly replayed support **676/34**, at unchanged
**852/870 coverage**. N and every geometric field remain identical. Reject
this variant without retuning any bound, calibration, threshold or vote rule;
historical 677/36 remains unexceeded. See the dated [journal](journal.md) and
`results/frozen_amplitude_profile_20260922/report.md` for full diagnostics.

`config/label_free_exploration.toml` is diagnostic configuration, not a new
molecular calibration. Its empty `[model]`, `[selection]` and `[preprocessing]`
sections do not replace the explicit original molecule config. Native support,
amplitude/geometry bounds, and selection settings remain fixed. Unknown adaptive
inputs still require their original selected summary.

`[counting_variable_projection]` supplies bounded optimizer budgets, linear
solver/KKT/SVD tolerances, model-mapping tolerances and the explicit native
elliptical iteration budget. They are numerical controls, not fitted physical
quantities or externally graded hyperparameters. `[representation]` supplies
arithmetic/serialization comparison tolerances only. `[acquisition_noise]`
declares a small lag window, background exclusion, sample sufficiency checks and
correlation/block diagnostics. These are exploratory validity rules, not
calibrated uncertainty guarantees or replacement production thresholds.

An empirical supervised classifier score is not an information-theoretic upper
bound on the signal. Conversely, stable label-free clusters, lower residuals,
view agreement or reproducible votes do not prove chemical identity. Parameters
must not be selected by repeatedly reading external benchmark grades.


The implemented auxiliary `Current` diagnostic reports direct coupled channel
scales/correlations only; the optional nuisance projection was not implemented.
It does not calibrate chemistry, modify fit weights or supply an independent
likelihood. A fixed Fisher replay exports new fold weights under unchanged
reconstruction settings and checks saved scores; it is neither chemical
calibration nor recovery of unique historical coefficients. Raw-QC plates also
retain and flag nonfinite samples that native preprocessing imputes.


## Masked-background synthetic prototype (2026-09-18)

The user approved starting the first alternative after the completed label-free
exploration. The initial deliverable is an opt-in prototype and synthetic
signal-preservation evidence, not a new production calibration or a real-data
result. Controls are explicit in `config/masked_robust_preprocessing.toml`.

Native median imputation precedes plane/row flattening. Restoring missing masks
afterward cannot remove the imputation's influence on valid pixels. The prototype
instead estimates the background from finite supplied-background pixels only.
It does not remove a transverse component from each molecular patch by default.
A correct foreground exclusion can protect an injected localized signal, but
this conditional property does not prove that a real exclusion is correct.
Gaussian tails, unrecognized molecules or instrumental structure can remain in
the chosen background.

Sequential plane fitting followed by row medians also has a specific limitation:
unequal x support across rows can make row offsets bias the fitted x slope.
Separate y tilt and arbitrary row offsets are not uniquely identifiable.
Signal-preservation checks must therefore include absolute background error,
not only the difference between signal-present and signal-absent images.

The predeclared synthetic comparison separates four operations: native
imputation/reference flattening; finite-only OLS without foreground exclusion;
finite guarded OLS; and finite guarded Huber. No control is selected using
benchmark outcomes. Coverage loss and unsupported/nonconverged cases remain
visible alongside same-support errors. Huber's fixed initial-residual MAD and
its numerical floor are not calibrated STM noise, and improved agreement does
not establish chemical identity or correct counts.

A possible later real-data comparison would retain the original four diagnostic
cases and saved geometry/support. Such geometry already depends on the original
preprocessing; freezing it makes comparisons controlled, not independent. No
new real-data execution or cluster submission is part of this synthetic-first
step. Registration search boundaries, GCV and production defaults stay unchanged.


Synthetic verification is complete: 441 engine assertions, 1,116 signal-study
assertions and 12,641 independent saved-table checks pass. No Huber control was
retuned. The tests demonstrate conditional preservation with a correct supplied
mask, but also sequential x-slope bias, outlier-sensitive initial scale, signal
attenuation with a leaky exclusion and explicit loss of unsupported foreground
coverage. Huber is not uniformly better than guarded OLS. There is no new STM
calibration, real-data result, calibrated uncertainty or production default.


The authorized four-scan real-data pilot keeps the synthetic prototype's
controls unchanged. It measures background changes, observed coverage and
forward/backward agreement, not molecular accuracy. Comparisons use fixed common
observed supports over the original lag grid; dropping an unavailable method or
silently comparing different pixel populations is not allowed. A failed method
can leave the all-method comparison unavailable while native-versus-other-method
comparisons remain separately reportable.

The level convention is one median per view on the same observable guarded
background support, held fixed over all lags. It is not truth alignment or gain
calibration. Better post-fit background statistics or image agreement would not
show that molecular contrast was preserved on real images. The exclusion and
fit support still depend on the old native preprocessing and saved geometry.
No new count, chemical fit, noise covariance, independent-sample estimate,
benchmark grade or production-default change is part of this pilot.

## Observed pixels from background correction through fitting (2026-09-24)

The full-profile measurement audit found nonfinite samples in 72 of 292 actual
acquisition views. Restoring their masks after native median-filled background
subtraction does not undo the fill's influence on acquired pixels. The authorized
next comparison therefore carries the raw observation mask through background
correction, smoothing, fixed-N fitting and residual patches. Its explicit
settings are in `config/observed_fit.toml`; production defaults are unchanged.

This is not a physical recalibration. Sigma, spacing, overlap, native optimizer
budgets, GCV selection, the `n_eff = max(10, ndata ÷ 9)` diagnostic placeholder,
and all chemical thresholds remain unchanged. N is taken from the first of the
two identical latest raw-count reproductions, not from expected chain length.
The input physical configuration still has benchmark-informed historical
provenance: label-free inference in this experiment does not establish a
strictly label-free calibration lineage.

The experiment tests consistent observation handling, including the background
used for residual patches; it is not an isolated optimizer comparison. On fully
acquired images preprocessing and fit inputs must be identical, and the fit is
reused exactly. On partial images background, ROI, geometry and descriptors may
change together. Shared-pixel residuals and optimizer termination are retained
as diagnostics, without selecting a method or retuning controls from benchmark
labels. A missing valid fit is a negative completeness result, not permission to
relax physical guards or grade a favorable subset.

Both complete fixed-N executions give **680/870 correct, 35/145 exact chains,
848/870 coverage**. This does not validate a physical recalibration: all bounds
are inherited, common-pixel RSS is lower on 17 partial scans and higher on 19
for each profile, and no LM fit declares convergence. The net exact-chain gain
occurs on fully observed scans with identical local inputs between arms; the
downstream Fisher/classifier dependence on the rest of the cohort must be
separated from local geometry quality. All 863 candidate decisions repeat,
but sixteen confidences change by 0.4. No champion or application promotion.

The subsequent frozen-learning diagnostic changes no physical calibration.
Its captured per-scan centers/scales and amplitude-based cluster names are
statistical transformations, not independently validated physical bounds or
chemical labels. Holding them fixed separates a local-input intervention from
cohort re-estimation in a declared order; it does not certify strict label-free
calibration, independent validation or a transferable optimum.

Its completed ordered comparison gives **678/28 → 671/25 → 665/17 → 674/21 →
680/35**, coverage **848/870** throughout. New local inputs alone lose seven
correct positions and three exact chains under the fixed control bank; changing
classifiers after Fisher/scaling gains fourteen exact chains. The reverse
frozen-input control is **680/36**. These order-dependent results do not validate
new geometry or any physical recalibration. Exact same-input repeatability of
captured predictive parameters is established for two executions, not calibration
independence, generalization or stability to changed measurements.

The subsequent whole-scan exclusion diagnostic removes each target from
Fisher/PCA, mixture, k-means and amplitude-based class naming. Its own local
normalization remains input-dependent preprocessing. This tests transfer of
the current statistical learners to an excluded scan, not an independent
physical recalibration or a newly collected holdout. Inherited sigma, spacing,
overlap, support and count-selection choices are unchanged; successful replay
or an improved external grade cannot certify their label-free provenance.

That exclusion test is now complete and negative: control **678/28 → 672/19**,
observed **680/35 → 673/29** (correct /870, exact /145), coverage **848/870**.
The two repetitions agree byte for byte, including captured predictive states.
No exact chain is gained. This demonstrates sensitivity to target participation
in unlabelled cohort learning, not evidence of label leakage in the learner or
a physical recalibration. The historical calibration limitation remains separate.

The subsequent diagnostic restricts all shared learning to fully observed
acquisitions, still excluding the target. Eligibility requires both raw views
finite at every acquisition-grid pixel and comes from the saved observation
audit, not known sequence, count or composition. It is not a guarantee of
artifact-free images, complete extracted boundary patches or independent
chemical calibration. All target scans keep their own unlabelled local
normalization. This differs from the earlier negative complete-patch rule,
which used eligible-row normalization and left k-means unchanged. The new
`selection.training_cohort` policy does not modify any physical/selection
parameter or repair inherited benchmark-informed calibration provenance.

Its completed, repeated comparison is negative: **670/870 correct, 16/145
exact** on control inputs and **670/15** on observed inputs, versus native
**678/28**, **680/35**, at the same **848/870** coverage. No native exact chain
is gained. Acquisition completeness alone does not improve this learning
cohort, nor validate physical calibration. Counts remain fixed: their 20
missing benchmark positions cap abstention-free coverage at **850/870**, below
historical **854/870**. This external diagnostic motivates independently
justified calibration/counting, not selecting parameters from expected N.
