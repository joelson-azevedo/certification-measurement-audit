# Analysis implementation guide

The `src/` package contains a deliberately selected set of validated numerical modules. All core module bytes are unchanged from the blind audit. [File provenance](../docs/FILE_PROVENANCE.csv) records exact original and staged hashes. The new `verify_release.py` provides portable local orchestration for a declared verification subset; it does not replace the complete original analysis pipeline.

| Purpose | Modules |
| --- | --- |
| Integrity and scoring | `data_integrity`, `scoring` |
| Classical statistics and alternatives | `ctt`, `distractors` |
| Binary structure | `tetrachoric`, `dimensionality`, `parallel_analysis` |
| Signed item-response estimation and diagnostics | `irt_2pl`, `irt_estimation`, `irt_diagnostics` |
| Residual dependence | `local_dependence`, `q3_reference_simulation` |
| Conditional regional screening | `dif`, `dif_matching`, `dif_sensitivity` |
| Common-scale research | `linking_transformations`, `linking_stocking_lord`, `linking_haebara`, `linking_anchor_diagnostics`, `linking_uncertainty` |
| Score correspondence | `score_distribution`, `irt_true_score_equating`, `irt_observed_score_equating`, `equating_decision_impact`, `equating_uncertainty` |
| Hypothetical decision behavior | `classification_probabilities`, `classification_posterior_predictive`, `classification_uncertainty`, `decision_consistency` |

The adjacent `configs/` directory preserves the exact recorded protocols from classical analysis through decision consistency. It includes seeds, item specifications, missingness rules, optimizer safeguards, reference grids, multiplicity and uncertainty conventions. Preserve item ordering and explicit unavailable statuses when using these APIs. The [methods appendix](../METHODOLOGY.md) explains the definitions and limitations.

Packaged tests use small controlled fixtures and independent numerical identities, rather than the undisclosed generation process for the delivered examination. One cache test's source-root path is adapted to this `analysis/` layout; its statistical assertions and all numerical modules are unchanged. `test_release.py` is new release verification. Internal workspace preservation tests, the complete stage runners, presentation-source tests, obsolete outputs and candidate-level predictions are not copied.

The source is a research implementation with explicit failure handling and conditional inference, not a deployment-ready scoring service. In particular, cache APIs expecting a source root should receive this `analysis/` directory. The [full replay plan](../docs/FULL_REPLAY_PLAN.md) explains why complete stage replay is unavailable in the curated candidate. No network access or operational scoring changes are performed by the quick verifier.

