# Controlled full-replay plan and missing components

**Tier 2 status: unavailable in this curated distribution.** This is a reference plan for a complete approved audit workspace, not a claim that the commands below can run here. No hidden simulator or generating truth is a prerequisite. The missing components are audit orchestration/provenance outputs and full scientific uncertainty archives.

## Process order and dependencies

| Original command in a complete workspace | Inputs and outputs | Principal cost/uncertainty |
| --- | --- | --- |
| `python -m src.run_stage1` | Delivered brief/manifest/CSV files to independently scored administrations and administered-item long data. | Exhaustive integrity/history/scoring checks; no latent estimation. |
| `python -m src.run_stage2` | Released processed data to classical statistics, distractors and candidate alpha bootstrap. | Within-form row bootstrap; retain undefined values and original population. [MP02] |
| `python -m src.run_stage3` | First-attempt matrices to raw/conditioned associations, parallel references and factor diagnostics. | Repeated tetrachoric inversion, matched binary references and factor optimization. [MP03] |
| `python -m src.run_stage4` | Item-response matrices to full signed calibrations, sensitivities, observed information and EAP. | Multi-start marginal likelihood, derivative/information and quadrature checks. [MP04] |
| `python -m src.run_stage5` | Saved calibrations/missing masks to residual reference draws and pair evidence. | Fixed-calibration response simulation and repeated posterior/residual calculations. [MP05] |
| `python -m src.run_stage6` | Saved calibration and first responses to matching, conditional logistic models, support and sensitivities. | Target-specific rest matching, nested regressions and bounded purification. [MP06] |
| `python -m src.run_stage7` | Saved item models/covariances and anchor evidence to common-scale links and paired uncertainty. | Candidate bootstrap recalibration and joint-covariance propagation. [MP07] |
| `python -m src.run_stage8` | Exact paired calibration/link caches to score distributions, correspondence and uncertainty. | Converged marginal integration for each retained parameter/link pair. [MP08] |
| `python -m src.run_stage9` | Saved source models, correspondence and paired caches to hypothetical joint decisions and posterior predictions. | Higher-resolution integration and paired uncertainty propagation. [MP09] |
| `python -m src.run_stage10` | Final reports and structured evidence to issue, action, gap and claim registers. | Source verification and synthesis; no new psychometric estimation. [E09] |

Each stage first validates earlier artifacts and gate status. In the original workflow, full project tests and separate stage verification scripts close the inventories after execution. Later stages depend on those checks, not just the existence of similarly named CSV files. This is why copying only runner modules would not create a supported full replay.

## Exact protocols and seeds

The staged [configuration directory](../analysis/README.md) preserves the recorded configurations without altering methods or parameterization. Alpha bootstrap uses 2,000 resamples and form seeds 20261021/20261022. Parallel analysis uses 100 valid references per form/specification with the documented seed formula. Signed calibration retains the initial 31/41-node diagnostics and final 81/101-node resolution checks. [MP02] [MP03] [MP04]

Residual reference simulation uses 500 valid main draws per form, with seeds 20261125/20261126 and separate pilot streams. Primary linking uncertainty uses 200 paired candidate bootstraps, seed 20261127; covariance sensitivity uses 1,000 paired draws, seed 20261227. Downstream equating and decision uncertainty reuse those exact pairs; they do not add independent samples. The controlled numerical Monte Carlo fixture has its separate seed and is not the certification-data simulator. [MP05] [MP07] [MP08] [MP09]

Full replay requires the original stage executors, stage-specific README/verification machinery, historical inventory files and checks, released processed matrices, full saved model/posterior/covariance artifacts, and attempt-matched uncertainty arrays/cache manifests. These are intentionally not all bundled. A portable replacement must either rebuild and validate the provenance chain safely in a separate workspace or define and independently test a new replay contract. Frozen audit sources must not be edited to bypass integrity checks.

## Environment and runtime

Core numerical requirements are pinned in `environment/requirements-quick.txt`. Full plotting and executive-PDF regeneration additionally require the inspected plotting/PDF toolchain and licensed fonts; that full environment is not distributed. No full-pipeline duration is claimed. The expensive steps are repeated tetrachoric references, multi-start calibration and paired bootstrap recalibration. Quick verification reports its actual measured runtime, which must not be extrapolated to full replay or other machines.

Providing omitted audit components, establishing rights to redistribute them, and running a complete replay in a clean workspace are open release-work items. Tier 1 success verifies a useful declared subset; it does not justify calling the curated package fully reproducible.

[E09]: ../docs/claims.md#e09
[MP02]: ../docs/claims.md#mp02
[MP03]: ../docs/claims.md#mp03
[MP04]: ../docs/claims.md#mp04
[MP05]: ../docs/claims.md#mp05
[MP06]: ../docs/claims.md#mp06
[MP07]: ../docs/claims.md#mp07
[MP08]: ../docs/claims.md#mp08
[MP09]: ../docs/claims.md#mp09

The author approved this explicitly limited portfolio scope during final release preparation. Full Tier 2 replay remains unavailable; it is not an unresolved human decision or an implemented development plan. Separate exact-candidate publication authorization remains required.
