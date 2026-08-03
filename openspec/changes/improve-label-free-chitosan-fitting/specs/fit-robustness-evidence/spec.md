## ADDED Requirements

### Requirement: Candidate comparison is predeclared
Before running comparative fits, the evaluation SHALL freeze the candidate configurations, scan list, perturbation grid, label-free metrics, tie rules, failure policy, and promotion criteria. The frozen declaration SHALL exclude benchmark truth and grade data.

#### Scenario: Comparison starts from a complete declaration
- **WHEN** a candidate comparison is launched
- **THEN** validation confirms that every candidate, scan, perturbation, metric, tie rule, and failure rule is declared and hash-bound before the first fit runs

#### Scenario: Candidate is added after results are observed
- **WHEN** an undeclared candidate or perturbation is introduced after comparison output exists
- **THEN** it is treated as a new versioned experiment and cannot be merged into the frozen comparison

### Requirement: Evidence covers representative scans and perturbations
The comparison SHALL include representative externally unlabeled short-chain and 10–20mer scans and SHALL measure stability under predeclared image perturbations relevant to STM acquisition and preprocessing. Labels SHALL NOT be used to choose scans, perturbation strengths, or winners.

#### Scenario: Perturbed scan remains within the declared operating range
- **WHEN** a scan is evaluated under each predeclared noise, drift, blur, crop, or flattening perturbation
- **THEN** the evidence reports fit success, selected count stability, GCV change, residual diagnostics, physical guards, and abstention for every candidate

### Requirement: Promotion uses label-free aggregate evidence
A challenger SHALL be eligible for promotion only when it passes all hard scientific and integrity gates and improves the predeclared aggregate label-free objective without unacceptable regressions in stability, abstention, physical validity, or runtime. A benchmark grade SHALL NOT determine eligibility.

#### Scenario: Challenger improves external benchmark accuracy only
- **WHEN** a challenger has no qualifying label-free improvement but obtains a better post-hoc benchmark grade
- **THEN** it remains ineligible for promotion

#### Scenario: Challenger improves GCV but destabilizes counts
- **WHEN** aggregate GCV improves while perturbation or held-out-scan stability crosses a predeclared regression limit
- **THEN** the challenger is rejected or retained as diagnostic according to the frozen failure policy

### Requirement: Comparison output is reproducible and auditable
The comparison SHALL emit deterministic machine-readable results, command logs, environment and Julia version information, input/output hashes, per-scan failure reasons, and a terminal PASS, FAIL, or BLOCKED verdict. Interrupted or partial runs SHALL NOT be represented as complete evidence.

#### Scenario: One shard fails or times out
- **WHEN** a Viper comparison shard exits nonzero, times out, or is interrupted
- **THEN** the merged evidence is marked incomplete or BLOCKED and identifies the missing shard without silently dropping its scans

#### Scenario: Identical frozen comparison is rerun
- **WHEN** the same code, configuration, inputs, seeds, and environment are used
- **THEN** normalized machine-readable outputs and the terminal verdict are reproducible

### Requirement: Heavyweight evaluation respects project execution boundaries
Focused single-file checks SHALL run locally, while multi-scan, exhaustive, or heavyweight comparisons SHALL run on Viper with a dry-run archived before submission. Concurrent CPU use SHALL respect the configured account limit, and compute SHALL NOT run on login nodes.

#### Scenario: Viper comparison is prepared
- **WHEN** a heavyweight comparison is ready for submission
- **THEN** an archived dry-run identifies exact inputs, outputs, array shape, CPU total, memory, walltime, and commands before any Slurm job is submitted

### Requirement: Scientific history is updated
Every experiment, failed approach, bug fix, parameter decision, and promotion verdict SHALL add a dated entry to `docs/src/journal.md`. Parameter additions or renames SHALL also update configuration and calibration documentation.

#### Scenario: Candidate experiment reaches a terminal verdict
- **WHEN** a comparison ends in PASS, FAIL, or BLOCKED
- **THEN** the journal records what ran, why, evidence paths, observed limitations, and the next decision without overstating benchmark or application validation
