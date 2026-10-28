# Can You Defend Your Pass/Fail Decisions?

**An Independent Measurement Audit of a Professional Certification Examination**

**Joelson Azevêdo · Independent Researcher**  
Measurement Science · Psychometrics · Statistical Validation · AI Evaluation

**Independent blind analysis of a synthetic professional certification dataset.** No real organization, examination, candidate or certification decision is represented. This independent portfolio case examines what evidence can support consequential pass/fail decisions. [E01]

[Executive memo](reports/EXECUTIVE_ONE_PAGE_SUMMARY.pdf) · [Executive brief](reports/EXECUTIVE_MEASUREMENT_BRIEF.pdf) · [Technical case study](CASE_STUDY.md) · [Run Tier 1](REPRODUCIBILITY.md)

## The measurement problem

Correct implementation of a passing rule is necessary, but it is only part of the evidence needed to defend certification. This audit reconstructed scoring, investigated item behavior, evaluated provisional form comparisons and examined how modeled decisions behave near a consequential boundary. Its central conclusion is: **scoring was reproducible within the supplied records, while substantive decision validity and operational form comparability remained unestablished.** [E10] [E08]

The fictional examination uses two forms with common and form-specific multiple-choice items. The included delivery contains synthetic administration histories/responses, the answer key, form/domain metadata and the original disclosure/manifest. [Data dictionary](DATA_DICTIONARY.md) and [data provenance](data/proposed/README.md) explain their roles. Statistical demonstration in synthetic records cannot independently validate certification conclusions for a real population. [E01]

## Principal findings

- **The recorded scoring process was reproducible.** Across 2,400 administrations, there were 0 raw-score and 0 pass/fail discrepancies against the supplied key and rule. This establishes implementation consistency within the delivery, not content or standard validity. [Q01] [Q05] [Q06]
- **Different items require different expert inquiries.** ANC_14 showed negative functioning, ANC_18 showed conditional regional differences, and B_11 had very weak discrimination with unstable model difficulty. None establishes a wrong key, confirmed bias or automatic invalidity. [E02] [E06] [E07]
- **The form-comparison question remains open.** First-attempt pass rates were A 73.40% and B 61.94%, from different, unadjusted groups. Under a provisional common-item link and the same assumed population, A=42 corresponds to B=38.70, with a 95% bootstrap interval 38.15-39.13. **This is not an authorized conversion or recommended new passing cut.** [Q09] [Q10] [Q12] [Q13] [Q14] [Q15]
- **An overall average conceals boundary uncertainty.** Under original rules, modeled cross-form agreement is 86.59% overall (95% bootstrap interval 85.66%-87.64%) and 46.30% at the fixed modeled proficiency where A's expected score equals 42. The latter is neither the observed raw-score group nor a validated mastery threshold. These are hypothetical independent administrations, not observed retest results or accuracy estimates. [Q16] [Q17] [Q18] [Q19] [Q20]


## What the author developed

Under the author's research direction, the project developed independent scoring reconstruction, classical item and distractor analysis, binary dimensionality assessment, signed-slope item-response estimation, residual diagnostics, conditional regional screening, provisional common-item linking, discrete score correspondence and hypothetical decision-consistency calculations. Numerical fixtures, source-cell checks and staged documentation support implementation review; they do not establish substantive certification validity. [MP02] [MP03] [MP04] [MP05] [MP06] [MP07] [MP08] [MP09] [E08]

## Inspect the work

| Reader need | Destination |
| --- | --- |
| Executive interpretation | [One-page memo](reports/EXECUTIVE_ONE_PAGE_SUMMARY.pdf), [five-page brief](reports/EXECUTIVE_MEASUREMENT_BRIEF.pdf) and [report attribution](reports/README.md) |
| Technical reasoning | [Case study](CASE_STUDY.md), [methods](METHODOLOGY.md) and [limitations](LIMITATIONS.md) |
| Implementation and evidence | [Code guide](analysis/README.md), [selected results](results/README.md), [figures/captions](figures/README.md) and [claim sources](docs/claims.md) |
| Local verification | [Reproduction commands](REPRODUCIBILITY.md), [requirements](environment/requirements-quick.txt), [data dictionary](DATA_DICTIONARY.md) and [provenance](docs/PROVENANCE.md) |

## Reproducibility scope

**Tier 1 — Verified reproducibility:** checks included data/manifests, reconstructs scoring and first attempts, recomputes classical point statistics and selected model calculations using frozen parameters, resolves claim cells and runs the packaged numerical tests. It includes local metadata/license checks and disables Python socket operations during verification. Run from this repository's root with the pinned environment; see [exact commands](REPRODUCIBILITY.md).

**Tier 2 — Full analytical replay unavailable:** complete stage executors, processed/model evidence, historical gates and matched uncertainty caches are not all distributed. The package does not refit the entire audit or regenerate its bootstrap intervals. [The full replay gap](docs/FULL_REPLAY_PLAN.md) is disclosed and accepted for this portfolio release. The project is not a fully reproducible archive of every analytical stage and supplies no independent empirical replication.

## Authorship and AI-assisted research

This independent measurement science project was conducted under the direction of Joelson Azevêdo, with AI assistance in software development, statistical implementation, automated testing, technical documentation, and editorial production.

Research objectives, methodological decisions, interpretation of findings, and publication decisions were reviewed and directed by the author.

The project includes numerical verification, automated testing, and reproducibility checks. These procedures do not constitute independent empirical validation of the certification decisions.

The project is an independent synthetic portfolio case, with no institutional affiliation beyond the stated Independent Researcher designation, no sponsorship/commissioning claim and no real-client engagement. The approved author attribution accompanies the unchanged executive PDFs without altering their original generic metadata label. Cite using [CITATION.cff](CITATION.cff).

## Licensing and preparation status

Original project code is **MIT**. Original documentation, reports, figures and synthetic data are **CC BY 4.0**, with third-party exceptions. Read [LICENSE_POLICY.md](LICENSE_POLICY.md), [MIT text](LICENSES/MIT.txt), [CC notice](LICENSES/CC-BY-4.0-NOTICE.md) and the [per-file scope](docs/LICENSE_SCOPE.csv).

This candidate is prepared locally for the author-approved account/repository and intended public visibility. The destination is not asserted to exist. **PUBLICATION STATUS: NOT PUBLISHED — AWAITING EXPLICIT HUMAN AUTHORIZATION.** No operational certification change is authorized.

[E01]: docs/claims.md#e01
[E02]: docs/claims.md#e02
[E06]: docs/claims.md#e06
[E07]: docs/claims.md#e07
[E08]: docs/claims.md#e08
[E10]: docs/claims.md#e10
[MP02]: docs/claims.md#mp02
[MP03]: docs/claims.md#mp03
[MP04]: docs/claims.md#mp04
[MP05]: docs/claims.md#mp05
[MP06]: docs/claims.md#mp06
[MP07]: docs/claims.md#mp07
[MP08]: docs/claims.md#mp08
[MP09]: docs/claims.md#mp09
[Q01]: docs/claims.md#q01
[Q05]: docs/claims.md#q05
[Q06]: docs/claims.md#q06
[Q09]: docs/claims.md#q09
[Q10]: docs/claims.md#q10
[Q12]: docs/claims.md#q12
[Q13]: docs/claims.md#q13
[Q14]: docs/claims.md#q14
[Q15]: docs/claims.md#q15
[Q16]: docs/claims.md#q16
[Q17]: docs/claims.md#q17
[Q18]: docs/claims.md#q18
[Q19]: docs/claims.md#q19
[Q20]: docs/claims.md#q20
