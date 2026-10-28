# Selected retained results

These are explicit selections from the blind audit, copied without changing original cells, ordering, field names or missing values. They are expected evidence for verification, not fresh Stage 13 estimates. [Provenance](../docs/FILE_PROVENANCE.csv) gives original/staged hashes. [Claims](../docs/claims.md) identify actual rows/fields and interpretation restrictions.

| Directory | Selection and reason |
| --- | --- |
| `stage1/` | Delivered counts, scoring summary, form membership, observed outcomes and original structured data dictionary. |
| `stage2/` | All primary form-item classical statistics, form reliability and distractor results. |
| `stage3/` | Population/coding sensitivities, raw/conditioned matrix diagnostics and parallel-analysis values. |
| `stage4/` | Saved parameter, uncertainty and fit tables including lower-resolution checks and final reporting models; explicit model IDs prevent accidental substitution. |
| `stage5/` | Fixed-calibration residual maxima, thresholds and simulation completion evidence. |
| `stage6/` | Primary within-form conditional regional results. |
| `stage7/` | Anchor specifications, exact scale-linking constants and conditional link uncertainty. |
| `stage8/` | Declared provisional score correspondence and common-reference marginal distributions. |
| `stage9/` | Same-form/cross-form and near-reference conditional agreement, parameter uncertainty and plotted conditional curves. |
| `stage10/` | Integrated issue, action, gap and claim registers plus their source-defined summary. |

`source_excerpts.json` contains selected literal original methodological passages, keyed by evidence ID. The registry records original source hashes and line spans. It allows the reader to inspect explanations without inaccessible machine-local paths or publishing the entire internal working archive.

Candidate-level proficiency/prediction files, counterfactual person decisions, complete covariance/posterior arrays and full uncertainty caches are intentionally omitted. A weak item's unavailable difficulty remains unavailable. Scoring facts, diagnostics, conditional model estimates and proposed actions are different evidence categories. Retained results do not constitute an operational conversion table or demonstrate historical certification error. [E07] [MP08] [MP09] [E09]

Quick verification regenerates aggregate scoring, classical point statistics and selected fixed-estimate model calculations in `verification/`, then compares them against these retained results. It does not overwrite this directory or claim to reproduce every original estimation.

[E07]: ../docs/claims.md#e07
[E09]: ../docs/claims.md#e09
[MP08]: ../docs/claims.md#mp08
[MP09]: ../docs/claims.md#mp09
