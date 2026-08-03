## ADDED Requirements

### Requirement: Fitting remains independent of external labels
The fitting and selection pipeline SHALL NOT read or derive decisions from expected lobe counts, unit sequences, benchmark truth, benchmark grades, class counts, or benchmark-specific filenames. External labels MAY be consumed only by a separate reporting command after the candidate configuration and artifacts are frozen.

#### Scenario: Forbidden truth input is supplied
- **WHEN** a fitting, calibration, candidate-ranking, thresholding, or abstention command receives a truth, expected-`N`, sequence, grade, or benchmark-manifest input
- **THEN** the command exits nonzero before producing fit or selection output and names the forbidden input boundary

#### Scenario: Frozen candidate is externally graded
- **WHEN** a candidate has been frozen without label access and a separate grader is invoked
- **THEN** the grader may report benchmark metrics without modifying the candidate configuration, selected model, thresholds, or abstention policy

### Requirement: Candidate behavior is configuration-bound
Every candidate SHALL declare preprocessing, physical model constraints, selection parameters, and abstention rules in a versioned TOML file using the established `[preprocessing]`, `[model]`, and `[selection]` sections. The evidence bundle SHALL bind the configuration and executable input artifacts by SHA-256.

#### Scenario: Hidden candidate default affects selection
- **WHEN** a candidate relies on a selection-affecting parameter that is neither present in its TOML nor inherited from an explicitly named frozen baseline
- **THEN** candidate validation fails before fitting

#### Scenario: Bound candidate input changes
- **WHEN** a configuration or input artifact differs from the SHA-256 recorded in the evidence bundle
- **THEN** comparison or promotion fails with an artifact-integrity error

### Requirement: GCV drives base model selection
The pipeline SHALL use GCV as the primary base-model selection criterion. BIC, AICc, residual diagnostics, and effective-sample-size calculations MAY act only as diagnostics or predeclared guards and SHALL NOT reinterpret or retune the existing `n_eff = n ÷ 9` placeholder.

#### Scenario: Candidate models are all physically valid
- **WHEN** multiple candidate fits satisfy the configured physical constraints
- **THEN** the base selection is determined by the declared GCV rule before any diagnostic guard is applied

#### Scenario: Information criterion disagrees with GCV
- **WHEN** BIC or AICc prefers a different base model from GCV
- **THEN** the evidence records the disagreement and applies only the predeclared guard behavior without replacing GCV as the base selector

### Requirement: Physical infeasibility produces abstention
The pipeline SHALL reject fits that violate configured spacing, overlap, support, or finite-output constraints. If no candidate remains physically valid, the pipeline SHALL abstain with machine-readable reasons instead of forcing an `N_selected` value.

#### Scenario: Every count candidate violates max overlap
- **WHEN** all fitted count candidates exceed the configured physical `max_overlap` constraint
- **THEN** the result contains no forced lobe count and records a physical-infeasibility abstention

### Requirement: One-dimensional fitting remains diagnostic-only
The 1D fitting path SHALL remain skipped by default and SHALL NOT feed `N_selected`. Enabling the 1D path SHALL add only cross-check diagnostics.

#### Scenario: One-dimensional diagnostic is enabled
- **WHEN** the user supplies `--no-skip-1d`
- **THEN** the 2D-derived `N_selected` is unchanged and the output adds a separately identified 1D cross-check
