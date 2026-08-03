## Context

The production path combines preprocessing from `STMSXMIO`, the 2D engine in `GaussianFit2D`, selection in `STMMolecularFit`, and orchestration in `test/batch_full.jl`. Its current benchmark headline is useful external evidence, but it cannot be used to choose fitting behavior. The repository also contains many standalone research scripts, so a challenger can appear promising without sharing one candidate contract, one failure policy, or one reproducible comparison bundle.

The scientific constraints are strict: GCV drives base-model selection; spatially correlated residuals make absolute iid information criteria unreliable; `N_selected` cannot use expected counts or sequences; 1D remains diagnostic-only; parameters live in TOML; and heavyweight comparisons run on Viper rather than login nodes or the local workstation. The objective is therefore not to maximize a benchmark score directly. It is to select the strongest frozen candidate under predeclared image-only evidence, then measure its benchmark effect separately.

## Goals / Non-Goals

**Goals:**

- Give every baseline and challenger the same typed inputs, outputs, failure semantics, and artifact-integrity checks.
- Compare candidates with paired, scan-level evidence led by GCV and guarded by physical validity, perturbation stability, residual structure, abstention, and runtime.
- Make partial, failed, timed-out, and abstaining runs visible rather than silently excluding them.
- Produce a terminal, reviewable decision that can promote a candidate by changing an explicit production config reference and can be rolled back by restoring that reference.
- Preserve a hard separation between candidate development and external benchmark grading.

**Non-Goals:**

- Replacing GCV with BIC/AICc or retuning `n_eff = n ÷ 9`.
- Using unit-assignment, expected-`N`, sequence, class-count, or benchmark-grade information to fit or rank candidates.
- Making the 1D diagnostic influence `N_selected`.
- Promoting constant-current unit-assignment diagnostics into lobe-count selection.
- Treating visually plausible 10–20mer output as benchmark validation.
- Repairing or promoting `STMMolecularFitGUI` as a production entrypoint.

## Decisions

### 1. Use a frozen comparison manifest

A versioned TOML comparison manifest will name the baseline and challenger config files, immutable scan-list artifact, perturbation matrix, deterministic seeds, runtime budget, output schema, metrics, tie policy, and terminal decision rule. A preflight command will reject forbidden label fields and hash every declared input before execution.

This is preferred over ad hoc command lines because the manifest makes the scientific decision surface reviewable before results exist. A general hyperparameter optimizer was rejected: optimizing against a benchmark leaks labels, while optimizing against one unlabeled scan risks overfitting acquisition artifacts.

### 2. Keep the current production path as an immutable baseline adapter

The comparison runner will call the existing preprocessing, 2D fit, and selector entrypoints through a thin candidate adapter. The current production configuration is candidate `baseline`; challengers use separate TOML files and output directories. No comparison run mutates production configs, selectors, registries, or cached baseline output.

This is preferred over copying fit logic into a research runner, which would make a favorable result impossible to attribute to the actual production path.

### 3. Build challengers from observed label-free failure modes

Before adding model complexity, a residual-audit stage will classify failures visible without truth: structured residual autocorrelation, edge/support truncation, overlap infeasibility, non-finite optimization, multi-start instability, and preprocessing sensitivity. Challenger families will be limited to changes that address a recorded failure class, such as preprocessing/background treatment, initialization/multi-start policy, or physical support constraints. Each family gets its own config and one stated mechanism.

Brute-force combinations are rejected because they multiply researcher degrees of freedom and make post-hoc selection likely. High-`N` disappearance must first trigger inspection of `max_overlap`, not selector retuning.

### 4. Use paired GCV improvement as the primary aggregate objective

For each scan and perturbation, the runner will emit the baseline and challenger GCV for the selected physically valid 2D base model. The comparison will use the paired relative difference

`delta_gcv = (gcv_challenger - gcv_baseline) / max(abs(gcv_baseline), eps(Float64))`,

where lower is better. Promotion requires the scan-bootstrap 95% upper confidence bound of the mean paired `delta_gcv` to be below zero. Resampling is by whole scan, never by pixel, because residual pixels are spatially correlated.

Absolute BIC/AICc ranking was rejected because its iid assumptions are unsuitable here. An unpaired aggregate was rejected because scan difficulty would dominate model differences.

### 5. Treat stability and physical behavior as non-inferiority gates

The runner will report, per candidate and scan:

- disagreement of `N_selected` across predeclared perturbations and deterministic reruns;
- physical-constraint rejection and no-valid-candidate rates;
- abstention and non-finite-output rates;
- residual autocorrelation/structure diagnostics;
- elapsed time and peak-memory evidence when available.

A challenger passes a guard only when its paired regression relative to baseline stays within the baseline's scan-bootstrap 95% uncertainty envelope. Physical invariant violations, hidden defaults, truth-boundary violations, missing outputs, or silently dropped scans are unconditional failures. Runtime has a separately frozen operational ceiling; it cannot be traded against benchmark accuracy.

Using uncertainty-derived non-inferiority bounds is preferred over choosing numerical tolerances after seeing candidate results.

### 6. Evaluate in dependency-ordered gates

1. Static firewall and manifest validation.
2. Focused local synthetic and single-file contract tests.
3. Deterministic rerun and perturbation smoke test on a predeclared small unlabeled panel.
4. Viper comparison over the complete frozen scan list, sharded by file with at most four Julia threads per task and at most eight CPUs concurrently.
5. Merge validation that requires every declared scan and shard.
6. Independent review of code, evidence, failure classes, and cleanup.
7. Optional external grade of the frozen winner in a separate process.

This staged design fails cheaply before expensive compute and prevents an incomplete HPC run from becoming promotion evidence.

### 7. Emit one normalized evidence bundle

The bundle will contain the manifest and hashes, environment metadata, exact commands, per-scan machine-readable fit/selection records, perturbation records, shard receipts, merge validation, summary statistics, logs, plots used only for diagnosis, and a terminal `PASS`, `FAIL`, or `BLOCKED` receipt. Normalized TSV/JSON content will exclude timestamps and machine-specific paths from reproducibility hashes.

Plots alone are rejected as evidence because they cannot prove denominator completeness, deterministic behavior, or leakage boundaries.

### 8. Promote by explicit config switch with rollback

No code default changes during comparison. After a confirmed PASS, promotion changes the explicit production config reference and updates all affected benchmark/application claims only after the separate external grade. Rollback restores the prior config reference and archived hash; evidence remains immutable.

## Risks / Trade-offs

- **[A GCV improvement can preserve the same count while changing shape quality]** → retain residual-structure, perturbation, and physical guards and archive visual diagnostics without using them as labels.
- **[Uncertainty-derived gates can be inconclusive on a small panel]** → return `BLOCKED` and expand only through a new predeclared scan-list version, never by selectively adding favorable scans.
- **[Candidate-family restriction may miss a useful architecture]** → allow new families only after a journaled residual mechanism and a new manifest version.
- **[Perturbations may not match all real acquisition artifacts]** → document the physical rationale and keep raw, unperturbed held-out scans as a separate gate.
- **[HPC interruption can create partial evidence]** → require complete shard manifests, explicit timeout/interruption receipts, and atomic merge output.
- **[Runtime and storage grow with candidate count]** → gate candidates on focused checks first, shard on Viper, and retain only declared evidence artifacts.
- **[External benchmark results may tempt post-hoc tuning]** → freeze once, grade separately, and require a new versioned change for any subsequent hypothesis.

## Migration Plan

1. Add manifest and output schemas plus failing firewall/integrity contract tests.
2. Implement the baseline adapter and reproduce current focused outputs byte-for-byte or within predeclared numerical tolerances.
3. Add residual-audit output and one challenger family at a time, each tied to a journaled failure mechanism.
4. Run local focused verification, archive a Viper dry-run, then execute and merge the frozen comparison.
5. Obtain independent evidence review. On `FAIL` or `BLOCKED`, leave production unchanged.
6. On confirmed `PASS`, freeze candidate hashes, switch the explicit production config reference, run the separate external grade, and update every affected headline and method document.
7. Roll back by restoring the prior production config reference and hash if post-promotion application checks reveal a regression.

## Open Questions

- Which first residual failure class is most prevalent under the frozen baseline audit: support/overlap, preprocessing/background, or optimizer instability?
- What operational runtime ceiling should be frozen for the Viper comparison based on the current baseline dry-run and account quota?
- Which long-chain scans form the first application panel without making claims of benchmark validation?
