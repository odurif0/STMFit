# DFT calculation note: GlcN and GlcNAc molds

## Purpose and scope

These Quantum ESPRESSO calculations produce local Tersoff–Hamann-like
occupied-state LDOS maps to distinguish a central glucosamine (GlcN) from a
central N-acetylglucosamine (GlcNAc) in an adsorbed chitosan chain on Cu(100).
They are not adsorption-energy calculations and do not simulate the tip or the
tunnelling current.

```text
GlcN case:   GlcN–GlcN–GlcN
GlcNAc case: GlcN–GlcNAc–GlcN
```

The chemical variable is the acetylation of the central unit. No experimental
sequence, composition or benchmark label enters the DFT setup.

## Atomic models

| Quantity | Central GlcN | Central GlcNAc |
|---|---:|---:|
| Trimer formula | `C18 H35 N3 O13` | `C20 H37 N3 O14` |
| Trimer atoms | 69 | 74 |
| Cu atoms (three-layer 8×6 Cu(100) slab) | 144 | 144 |
| Total atoms | 213 | 218 |
| Cell `a × b × c` (Å) | 20.4495 × 15.3371 × 21.4863 | 20.4495 × 15.3371 × 21.8118 |

The lower 96 Cu atoms are frozen; the top Cu layer and the trimer relax. Both
systems share the lateral Cu(100) registry and the central-ring frame. The slab
top is at z ≈ 3.615 Å, with about 12 Å of vacuum. Structures and input
templates: `hpc/qe_molds/`.

## Electronic structure

Both systems use QE 7.4.1 `pw.x`/`pp.x`, PBE with Grimme D3, scalar
non-spin-polarized PAW (PSLibrary/KJPAW; Cu `pbe-dn-kjpaw`, 11 valence
electrons), `ecutwfc = 50 Ry`, `ecutrho = 360 Ry`, Γ-only sampling, and
Marzari–Vanderbilt smearing `degauss = 0.02 Ry`, on 8 MPI ranks with 96 GB.
This reduced pilot model replaced an `8×8×4`, 80/640 Ry, 2×2×1 setup that
exceeded memory. Substrate-size, k-point and cutoff convergence remain open.

## Relaxation and SCF acceptance

BFGS with `forc_conv_thr = 6e-3 Ry/Bohr` and `etot_conv_thr = 1e-3 Ry`
("relaxed enough", chosen from force plateaus and walltime, not from labels).
Relaxation SCF threshold 1e-6 Ry; final SCF target 1e-7 Ry.

- **GlcN** (Raven 28784933, restart 5): relaxation and final SCF converged
  (37 iterations, 6e-8 Ry).
- **GlcNAc** (Raven 28790944): relaxation completed, but the final SCF stalled
  at 3–4e-5 Ry under two different mixing schemes (300 iterations each).

The **common acceptance criterion is 5e-5 Ry** (below 0.7 meV for the cell).
Continuing the exact GlcNAc retry state on Viper, plain (10640236) and
local-TF (10640237) mixing each converged in one iteration at 4.027e-5 Ry, with
identical energies (−31763.44437132 Ry). Their cubes are byte-identical. The
plain cube is canonical; the TF cube is independent numerical-path validation.

## LDOS quantity

`pp.x` with `plot_num = 5` (Tersoff–Hamann STM mode). QE uses the window
between E_F and E_F + bias, and ignores `emin`, `emax` and `degauss_ldos`.
`sample_bias = −0.0220495933 Ry` (−0.300 eV), the bias of the 240817 scans.
Accepted cubes (PP jobs 10640445/10640446):

| Cube | SHA-256 |
|---|---|
| `qe/glcn_restart5/glcn_central_ldos.cube` | `80cd1d1fde94cf084cc7ea464d2bf065b36b2c015ea0bfaef8cebfee8ff88863` |
| `qe/glcnac/glcnac_central_ldos.cube` | `40649ccd9eb6768444b8ff61bf4a639b3940cf926fe3eb42254eb024faf9b5bf` |

The maps are local LDOS proxies: no explicit tip orbital, tip relaxation,
finite-temperature transport or calibrated constant-current feedback.

## Qualification of the STM observable

A September 24–29 audit asked whether the accepted cubes support a physical
constant-current observable. Its findings (details: archived journal, tag
`archive/pre-cleanup-20260929`):

- **Signed quantity.** QE weights out-of-window states with the smearing
  function, whose cold-smearing derivative is negative. The GlcNAc cube has
  131,984/10,800,000 negative values (minimum −0.0038); GlcN has none. On the
  0.50 nm ring-frame plane, 35/289 GlcNAc samples are negative.
- **Γ reconstruction defect.** Stock QE omits conjugate reciprocal-space
  entries at Γ and adds a constant floor. A patched build removes it exactly
  (job 11976185). Negative samples remain and are attributed to PAW
  augmentation. Direct wavefunction reconstruction (job 11977502) passes its
  checks and shows non-monotonic GlcNAc vacuum tails.
- **Isovalue bracket.** The legacy search includes the rise toward the periodic
  substrate image. At 0.50 nm the GlcN maximum is below the legacy lower bound.
- **Cutoff and basis sensitivity.** Paired 360/720 Ry density cutoffs
  (11990241/11991069) and 50/60 Ry wavefunction cutoffs change vacuum-tail
  quantities without removing rising columns or negative fractions. These are
  sensitivities at two points, not convergence.
- **Fixed-density precision.** Davidson (11998561/11999206) and CG
  (12018453/12023139) NSCF comparisons at 1e-10/1e-12 Ry did not complete:
  unconverged eigenvalues, then time limits. There is no precision
  qualification.
- **Image controls.** On raw scans, a chemistry-free local plane fits patches
  as well as the chemical molds (RMS 5.45 pm against 6.0 pm), so height and
  envelope confound chemistry at the patch level.

**Conclusion.** The accepted SCF states and the 5e-5 Ry criterion stand. No
physically calibrated constant-current observable was established, and no
mold was replaced.

## Production molds

`templates/chitosan_cc_molds_native_v1.tsv` is built from the two accepted
cubes and the tracked frames (`templates/qe_frames/`) by
`test/build_cc_molds_native.jl` ([DFT molds](qe_stm_molds.md)). It keeps the
historical first-index-fast cube addressing and the legacy isovalue rule, which
replay the champion lineage. The molds are therefore empirical templates:
their benchmark value is measured by external grading, not by the physics
above. Changing the criterion, the addressing or the cubes makes a new
versioned mold file and a new graded experiment.

## Reproducibility map

| Item | Location |
|---|---|
| Structures, slab metadata, QE input templates | `hpc/qe_molds/` |
| QE outputs and cubes (large, untracked) | `qe/` locally, `/ptmp/oldu/...` and `/u/oldu/...` on the cluster |
| Frames, molds and provenance | `templates/` |
| Retired preparation, audit and diagnostic scripts | tag `archive/pre-cleanup-20260929` |
