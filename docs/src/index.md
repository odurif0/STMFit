# STMFit

STMFit reconstructs molecular chains from raw STM images without labels. It
counts the units of each chain, locates them, and assigns each unit GlcN (0) or
GlcNAc (1) with an explicit posterior uncertainty. The reference system is
chitosan on Cu(100).

## Current state (2026-09-29)

The promoted method is per-scan GCV counting, then repeated-scan molecule
consensus, corroborated per-scan assignment and latent-class fusion. A fresh
raw run (Viper job 12025539, source `3c414a4`) was graded externally on the
145-scan NKNNKN 6-mer benchmark:

| Profile | Exact N /145 | Correct /870 | Exact chains /145 | Classified /870 | Errors |
|---|---:|---:|---:|---:|---:|
| Previous record | 129 | 694 | 43 | 855 | 161 |
| **Promoted method** | **137** | **772** | **88** | **865** | **93** |

Four runs gave 772, 772, 772 and 769 correct (88 exact chains and 137/145
exact N each); the spread comes from the time-limited global search.

What this does and does not establish:

- Inference is label-free end to end. Benchmark membership and the `NKNNKN`
  control are read only by the grading scripts.
- The counting parameters in `config/chitosan.toml` were historically chosen
  with benchmark grades. They are frozen, but they are not independently
  calibrated ([Calibration](calibration.md)).
- NKNNKN places GlcNAc next to the chain ends, so edge adjacency and chemistry
  are confounded. A control molecule with another sequence is needed for
  chemical validation.
- The unknown 10–20mer application runs, but its counts and calls are not
  validated ([Runbook](chitosan_runbook.md)).

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
