# STMFit

STMFit reconstructs molecular chains from raw STM images without labels. It
counts the units of each chain, locates them, and assigns each unit to one of
two physical classes (0 = dimmer, 1 = brighter) with an explicit posterior
uncertainty. The chemical identity of the classes is not established. The
reference system is chitosan on Cu(100).

## Current state (2026-10-25)

The promoted method is per-scan GCV counting, then repeated-scan molecule
consensus, corroborated per-scan assignment and latent-class fusion. Since
2026-10-25 the count uses the capped axial support and the later fixed-N
refits the legacy one (decoupled supports). A fresh raw run (Viper job
12380719, source `70a9baa`) was graded externally on the 145-scan NKNNKN 6-mer
benchmark:

| Profile | Exact N /145 | Correct /870 | Exact chains /145 | Classified /870 | Errors |
|---|---:|---:|---:|---:|---:|
| Previous record | 129 | 694 | 43 | 855 | 161 |
| Legacy support (job 12042058) | 137 | 771 | 88 | 865 | 94 |
| **Promoted method** (decoupled supports) | **137** | **776** | **85** | **865** | **89** |

By molecule (41 molecules, each counted once; one holds 43 of the 145 scans):
exact N 87.8% [76.8, 96.3], correct units 79.5% [74.6, 84.3], exact chains
26.8% [14.6, 41.5], 1.13 errors per scan (95% intervals, bootstrap over
molecules). On the 492-scan full set (117 molecules): correct units 72.8%
[69.1, 76.7], exact chains 22.3% [15.4, 29.9]. Decisions use the molecule-level
figures (user, 2026-10-26).

The decoupled supports were adopted by the user's decision after the grade:
the declared rule required no regression on any figure, and exact chains went
from 88 to 85 while errors went from 94 to 89. They were adopted for the long
chains, whose truncated supports they remove; their counts are closer to 17
experimenter length readings (mean difference 1.0 unit, 1.5 before), a
comparison only, since 10–20mer lengths are not known. Identical runs give identical outputs.

What this does and does not establish:

- Inference is label-free end to end. Benchmark membership and the `NKNNKN`
  control are read only by the grading scripts.
- The counting parameters in `config/chitosan.toml` (the base of
  the promoted `config/chitosan_unified_hybrid.toml`) were historically chosen
  with benchmark grades. They are frozen, but they are not independently
  calibrated ([Calibration](calibration.md)).
- NKNNKN (N = GlcNAc, K = GlcN) places its two GlcN next to the chain ends.
  The brighter class falls on them, which contradicts the former naming rule
  "brighter = GlcNAc", and edge adjacency and chemistry are confounded. A
  control molecule with another sequence is needed for chemical validation.
- The unknown 10–20mer application runs. Its chain lengths are not known
  (the experimenters' length readings are a comparison, not a truth), and its
  unit calls do not match their acetylation readings beyond chance. Nothing
  there is validated yet ([Runbook](chitosan_runbook.md)).
- Benchmark-free checks on the 10–20mer (2026-10-25, [journal](journal.md)):
  the per-unit amplitude repeats across scans, but no more than in the
  fixed-sequence 6-mer, and reads K against N with a ceiling of about 55%. The
  6-mer cannot validate an N/K reading (it is a position code). The 6-mer
  subset the assignment is learned on moves one call in five, and 18 of 65
  repeated molecules differ by more than 2 units in count.

## Where to go

| Page | Content |
|---|---|
| [Pipeline](pipeline.md) | Data flow, packages, stages and outputs |
| [Model selection](selection.md) | Count selection per scan, molecule consensus |
| [Unit assignment](unit_assignment.md) | Features, molds, classifiers, vote, fusion, uncertainty |
| [Calibration](calibration.md) | What each parameter is based on, and its limits |
| [Configuration](config.md) | Every config file and key |
| [Runbook](chitosan_runbook.md) | Running, verifying and grading; 10–20mer application |
| [DFT molds](qe_stm_molds.md), [DFT note](dft_calculation_note.md) | GlcN/GlcNAc LDOS molds from Quantum ESPRESSO |
| [Mathematical background](math.md) | Model, GCV, κ penalty, diagnostics |
| [Research journal](journal.md) | Decisions, rejected approaches, open questions |
| [API](api.md) | Package reference |
