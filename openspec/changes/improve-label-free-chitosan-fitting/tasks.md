## 1. Freeze the Comparison Contract

- [ ] 1.1 Read `docs/src/journal.md` and inventory the current production preprocessing, 2D fit, selector, abstention, and HPC entrypoints in a dated baseline journal entry.
- [ ] 1.2 Add a versioned comparison-manifest template that declares candidate configs, scan list, perturbations, seeds, metrics, tie/failure rules, runtime budget, and output paths without label fields.
- [ ] 1.3 Add failing contract tests for truth/expected-`N`/sequence/grade inputs, hidden selection defaults, missing hashes, tampered artifacts, and post-result manifest mutation.
- [ ] 1.4 Implement manifest parsing, forbidden-input rejection, completeness checks, and SHA-256 binding until the contract tests pass.

## 2. Establish the Immutable Baseline

- [ ] 2.1 Add a thin candidate adapter that invokes shared `STMSXMIO` preprocessing, `GaussianFit2D` fitting, and `STMMolecularFit` selection without duplicating their logic.
- [ ] 2.2 Define the current production TOML as candidate `baseline` and emit a normalized per-scan record containing GCV, `N_selected`, physical guards, abstention, residual diagnostics, elapsed time, and failure reason.
- [ ] 2.3 Add focused regression tests proving the adapter preserves current 2D selection and that enabling `--no-skip-1d` adds diagnostics without changing `N_selected`.
- [ ] 2.4 Add tests proving no physically valid candidate yields a machine-readable abstention instead of a forced count.

## 3. Diagnose Label-Free Failure Modes

- [ ] 3.1 Implement residual-audit output for spatial structure, support/edge truncation, overlap rejection, non-finite optimization, multi-start disagreement, and preprocessing sensitivity.
- [ ] 3.2 Run focused local audits on the predeclared clean short-chain scan and one long-chain scan, archive exact commands, and record observed failure classes in the journal.
- [ ] 3.3 Freeze separate challenger TOMLs for only the observed mechanisms, keeping preprocessing, model, and selection parameters in their required config sections.
- [ ] 3.4 Add mechanism-specific tests showing each challenger changes only its declared behavior and leaves GCV as the base selector.

## 4. Build Robustness Evidence

- [ ] 4.1 Implement deterministic noise, drift, blur, crop, and flattening perturbations from manifest-declared strengths and seeds without reading labels.
- [ ] 4.2 Implement paired baseline/challenger records for every scan and perturbation, including failed and abstaining runs in the declared denominator.
- [ ] 4.3 Implement whole-scan bootstrap confidence intervals for paired relative GCV and baseline-relative non-inferiority checks for stability, physical failures, abstention, and residual structure.
- [ ] 4.4 Add deterministic synthetic and focused real-file tests for perturbation replay, paired denominators, bootstrap reproducibility, and unconditional integrity failures.

## 5. Produce Auditable Comparison Bundles

- [ ] 5.1 Implement atomic shard output with manifest/config/input hashes, environment and Julia versions, exact commands, normalized per-scan records, and explicit timeout/interruption receipts.
- [ ] 5.2 Implement a merge validator that requires every declared shard and scan, rejects duplicate or incompatible records, and emits terminal `PASS`, `FAIL`, or `BLOCKED` JSON.
- [ ] 5.3 Add adversarial tests for missing shards, silently dropped scans, duplicate scans, hash mismatch, non-finite metrics, timeout, interruption, and deterministic rerun.
- [ ] 5.4 Verify protected production configs and selectors remain unchanged throughout comparison tests.

## 6. Add the Viper Execution Lane

- [ ] 6.1 Add a Viper array launcher and worker that shard by scan, use at most four Julia threads per task, keep total concurrent CPUs at or below eight, and never run fits on the login node.
- [ ] 6.2 Add `--dry-run` output covering exact manifest, configs, scan list, output directory, array shape, CPU total, memory, walltime, sync, submission, merge, and fetch commands.
- [ ] 6.3 Archive and independently review a no-submit dry-run for the frozen comparison before enabling submission.
- [ ] 6.4 Execute the frozen Viper comparison, fetch all shards, run the merge validator, and archive cleanup evidence for jobs and temporary remote state.

## 7. Decide and Document

- [ ] 7.1 Independently review the firewall, candidate isolation, denominator completeness, perturbation evidence, statistical decision, adversarial failures, and cleanup receipt.
- [ ] 7.2 On `FAIL` or `BLOCKED`, leave production unchanged and journal the terminal mechanism and next hypothesis; on confirmed `PASS`, freeze all candidate and evidence hashes.
- [ ] 7.3 For a confirmed winner only, switch the explicit production config reference and verify rollback by restoring the prior reference in a focused test.
- [ ] 7.4 Run external benchmark grading only after the winner is frozen, then update every affected README/runbook/method headline while keeping 6mer benchmark and 10–20mer application claims separate.
- [ ] 7.5 Update `docs/src/journal.md`, `docs/src/config.md`, and `docs/src/calibration.md` for all experiments and any added or renamed parameters.
