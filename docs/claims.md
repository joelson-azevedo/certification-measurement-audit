# Evidence guide

Every cited claim resolves to retained data, configuration or a literal source excerpt. Original evidence IDs are retained where available.

<a id="q03"></a>
## Q03

Scored/administered items in each form; independently corroborated by both membership rows.

**Qualification:** Metadata-defined form membership, not proof of content validity.

- [results/stage1/consistency_by_form.csv](../results/stage1/consistency_by_form.csv) - CSV data row 1, field `long_rows_per_sitting_min`.

<a id="q04"></a>
## Q04

Common bank-item count; form observations remain separate.

**Qualification:** Common designation is not independent anchor invariance evidence.

- [results/stage10/stage10_summary.json](../results/stage10/stage10_summary.json) - JSON field `shared_bank_items`.

<a id="q01"></a>
## Q01

All supplied administrations covered by scoring reconstruction.

**Qualification:** Internal consistency with supplied records and key does not authenticate external records or validate the standard.

- [results/stage1/operational_pass_rates.csv](../results/stage1/operational_pass_rates.csv) - CSV data row 1, field `administrations`.

<a id="q02"></a>
## Q02

First observed attempts, distinct candidates across forms.

**Qualification:** Observed administrative fact; separate original samples, no adjustment for assignment or cohort composition.

- [results/stage1/operational_pass_rates.csv](../results/stage1/operational_pass_rates.csv) - CSV data row 4, field `administrations`.

<a id="q11"></a>
## Q11

Unchanged original raw cut used for both forms; source cut_second is also verified.

**Qualification:** A recorded operational rule, not a validated mastery criterion.

- [results/stage9/cross_form_decision_consistency.csv](../results/stage9/cross_form_decision_consistency.csv) - CSV data row 1, field `cut_first`.

<a id="q05"></a>
## Q05

Raw-score reconstruction discrepancies over all supplied administrations.

**Qualification:** Correct supplied-key implementation, not confirmation that the key or standard is substantively valid.

- [results/stage1/audit_summary.json](../results/stage1/audit_summary.json) - JSON field `raw_score_mismatches`.

<a id="q06"></a>
## Q06

Recorded pass/fail reconstruction discrepancies.

**Qualification:** No inference about actual competence or historical decision correctness.

- [results/stage1/audit_summary.json](../results/stage1/audit_summary.json) - JSON field `pass_fail_mismatches`.

<a id="q07"></a>
## Q07

Form A integrated conditional North-minus-South credit probability difference, displayed in percentage points.

**Qualification:** Estimated matching proxy, supported within-form reference; not an overall pass-rate difference, cross-form drift or confirmed bias.

- [results/stage6/dif_primary_results.csv](../results/stage6/dif_primary_results.csv) - CSV data row 18, field `signed_reference_mean_difference`.

<a id="q08"></a>
## Q08

Form B integrated conditional North-minus-South credit probability difference, displayed in percentage points.

**Qualification:** Same qualification as Q07; different form-specific matching reference.

- [results/stage6/dif_primary_results.csv](../results/stage6/dif_primary_results.csv) - CSV data row 78, field `signed_reference_mean_difference`.

<a id="q09"></a>
## Q09

Observed first-attempt Form A pass rate.

**Qualification:** Observed administrative fact; separate original samples, no adjustment for assignment or cohort composition.

- [results/stage1/operational_pass_rates.csv](../results/stage1/operational_pass_rates.csv) - CSV data row 5, field `pass_rate`.

<a id="q10"></a>
## Q10

Observed first-attempt Form B pass rate.

**Qualification:** Observed administrative fact; separate original samples, no adjustment for assignment or cohort composition.

- [results/stage1/operational_pass_rates.csv](../results/stage1/operational_pass_rates.csv) - CSV data row 6, field `pass_rate`.

<a id="q12"></a>
## Q12

Reference A score in provisional research correspondence.

**Qualification:** Reference in a model-based mapping; not an authorization to convert scores.

- [results/stage8/provisional_cut_correspondence.csv](../results/stage8/provisional_cut_correspondence.csv) - CSV data row 5, field `reference_A_score`.

<a id="q13"></a>
## Q13

Continuous interpolated observed-score B correspondence to reference A score.

**Qualification:** Research-only: not an authorized conversion or recommended new passing cut. Bootstrap interval must appear adjacent; fixed anchor/model/reference assumptions.

- [results/stage8/provisional_cut_correspondence.csv](../results/stage8/provisional_cut_correspondence.csv) - CSV data row 5, field `observed_inverse_B_continuous`.

<a id="q14"></a>
## Q14

Lower endpoint of source 95% pointwise bootstrap interval for correspondence.

**Qualification:** Parameter/link uncertainty conditional on fixed selection/model/reference; not candidate competence or a recommended cutoff.

- [results/stage8/provisional_cut_correspondence.csv](../results/stage8/provisional_cut_correspondence.csv) - CSV data row 5, field `observed_inverse_bootstrap_low_95`.

<a id="q15"></a>
## Q15

Upper endpoint of same source 95% bootstrap interval.

**Qualification:** Same qualification as Q14.

- [results/stage8/provisional_cut_correspondence.csv](../results/stage8/provisional_cut_correspondence.csv) - CSV data row 5, field `observed_inverse_bootstrap_high_95`.

<a id="q16"></a>
## Q16

Cross-form A/B hypothetical agreement averaged over the common Form A standard-normal reference, original cuts.

**Qualification:** Independent responses conditional on fixed modeled proficiency, no learning/exposure; not observed retest agreement or accuracy.

- [results/stage9/cross_form_decision_consistency.csv](../results/stage9/cross_form_decision_consistency.csv) - CSV data row 1, field `agreement`.

<a id="q17"></a>
## Q17

Lower 95% bootstrap endpoint for global cross-form agreement.

**Qualification:** Parameter/link uncertainty, not repeated-exam outcome variation; model/design uncertainty remains.

- [results/stage9/classification_parameter_uncertainty.csv](../results/stage9/classification_parameter_uncertainty.csv) - CSV data row 96, field `low_95`.

<a id="q18"></a>
## Q18

Upper 95% bootstrap endpoint for global cross-form agreement.

**Qualification:** Same qualification as Q17.

- [results/stage9/classification_parameter_uncertainty.csv](../results/stage9/classification_parameter_uncertainty.csv) - CSV data row 96, field `high_95`.

<a id="q20"></a>
## Q20

Descriptive expected-score reference for the conditional estimate.

**Qualification:** Not an observed-score group or independently validated competence cut.

- [results/stage9/near_cut_classification_uncertainty.csv](../results/stage9/near_cut_classification_uncertainty.csv) - CSV data row 17, field `A_expected_score_location`.

<a id="q19"></a>
## Q19

Conditional A/B agreement at fixed proficiency where Form A's model-expected score equals 42.

**Qualification:** Descriptive proficiency location, not people observed scoring 42 or a validated mastery threshold; same hypothetical pairing as global, different conditioning.

- [results/stage9/near_cut_classification_uncertainty.csv](../results/stage9/near_cut_classification_uncertainty.csv) - CSV data row 17, field `agreement`.

<a id="vq01"></a>
## VQ01

Form A correlation of ANC_14 credit with the other-item total.

**Qualification:** Sample-specific negative alignment; does not establish a wrong key.

- [results/stage2/ctt_item_statistics.csv](../results/stage2/ctt_item_statistics.csv) - CSV data row 14, field `corrected_item_total_r`.

<a id="vq02"></a>
## VQ02

Form B counterpart.

**Qualification:** Same qualification as VQ01.

- [results/stage2/ctt_item_statistics.csv](../results/stage2/ctt_item_statistics.csv) - CSV data row 74, field `corrected_item_total_r`.

<a id="vq03"></a>
## VQ03

B_11 correlation of credit with the other-item total.

**Qualification:** Near-zero sample alignment does not establish an invalid item; model difficulty remains unavailable.

- [results/stage2/ctt_item_statistics.csv](../results/stage2/ctt_item_statistics.csv) - CSV data row 91, field `corrected_item_total_r`.

<a id="d01"></a>
## D01

Retake administrations; not independent repeat validation.

**Qualification:** Retake administrations; not independent repeat validation.

- [results/stage1/audit_summary.json](../results/stage1/audit_summary.json) - JSON field `retakes`.

<a id="g01"></a>
## G01

Integrated issues, not confirmed defects.

**Qualification:** Integrated issues, not confirmed defects.

- [results/stage10/stage10_summary.json](../results/stage10/stage10_summary.json) - JSON field `distinct_integrated_issues`.

<a id="g02"></a>
## G02

Proposed actions; no implementation authority.

**Qualification:** Proposed actions; no implementation authority.

- [results/stage10/stage10_summary.json](../results/stage10/stage10_summary.json) - JSON field `proposed_actions`.

<a id="g03"></a>
## G03

Evidence gaps; external availability remains unknown.

**Qualification:** Evidence gaps; external availability remains unknown.

- [results/stage10/stage10_summary.json](../results/stage10/stage10_summary.json) - JSON field `external_evidence_gaps`.

<a id="d02"></a>
## D02

First observed Form A group size

**Qualification:** Different observed candidates; assignment not independently established.

- [results/stage1/operational_pass_rates.csv](../results/stage1/operational_pass_rates.csv) - CSV data row 5, field `administrations`.

<a id="d03"></a>
## D03

First observed Form B group size

**Qualification:** Different observed candidates; assignment not independently established.

- [results/stage1/operational_pass_rates.csv](../results/stage1/operational_pass_rates.csv) - CSV data row 6, field `administrations`.

<a id="d05"></a>
## D05

Primary complete-case dimensionality population A

**Qualification:** Selection and missingness assumptions remain unresolved.

- [results/stage3/dimensionality_sensitivity.csv](../results/stage3/dimensionality_sensitivity.csv) - CSV data row 1, field `n_candidates`.

<a id="d06"></a>
## D06

Primary complete-case dimensionality population B

**Qualification:** Selection and missingness assumptions remain unresolved.

- [results/stage3/dimensionality_sensitivity.csv](../results/stage3/dimensionality_sensitivity.csv) - CSV data row 4, field `n_candidates`.

<a id="r01"></a>
## R01

First-attempt within-form alpha A

**Qualification:** Score consistency; not dimensionality, classification accuracy or validity.

- [results/stage2/test_statistics_by_form.csv](../results/stage2/test_statistics_by_form.csv) - CSV data row 1, field `cronbach_alpha`.

<a id="r02"></a>
## R02

First-attempt within-form alpha B

**Qualification:** Score consistency; not dimensionality, classification accuracy or validity.

- [results/stage2/test_statistics_by_form.csv](../results/stage2/test_statistics_by_form.csv) - CSV data row 2, field `cronbach_alpha`.

<a id="l01"></a>
## L01

Stocking-Lord B-to-A scale_A

**Qualification:** Provisional S18 link, fixed common-item selection; no operational interchangeability.

- [results/stage7/linking_constants_by_method.csv](../results/stage7/linking_constants_by_method.csv) - CSV data row 5, field `scale_A`.

<a id="l02"></a>
## L02

Stocking-Lord B-to-A shift_B

**Qualification:** Provisional S18 link, fixed common-item selection; no operational interchangeability.

- [results/stage7/linking_constants_by_method.csv](../results/stage7/linking_constants_by_method.csv) - CSV data row 5, field `shift_B`.

<a id="l03"></a>
## L03

Haebara-directional-A B-to-A scale_A

**Qualification:** Provisional S18 link, fixed common-item selection; no operational interchangeability.

- [results/stage7/linking_constants_by_method.csv](../results/stage7/linking_constants_by_method.csv) - CSV data row 6, field `scale_A`.

<a id="l04"></a>
## L04

Haebara-directional-A B-to-A shift_B

**Qualification:** Provisional S18 link, fixed common-item selection; no operational interchangeability.

- [results/stage7/linking_constants_by_method.csv](../results/stage7/linking_constants_by_method.csv) - CSV data row 6, field `shift_B`.

<a id="c01"></a>
## C01

Population-integrated same-form modeled agreement A

**Qualification:** Hypothetical independent administrations, fixed proficiency, original cut; not observed retest agreement.

- [results/stage9/same_form_decision_consistency.csv](../results/stage9/same_form_decision_consistency.csv) - CSV data row 1, field `agreement`.

<a id="c02"></a>
## C02

Population-integrated same-form modeled agreement B

**Qualification:** Hypothetical independent administrations, fixed proficiency, original cut; not observed retest agreement.

- [results/stage9/same_form_decision_consistency.csv](../results/stage9/same_form_decision_consistency.csv) - CSV data row 5, field `agreement`.

<a id="c03"></a>
## C03

Same-form agreement bootstrap interval A low_95

**Qualification:** Pointwise parameter/link uncertainty conditional on fixed model, population, independent administrations and original rules.

- [results/stage9/classification_parameter_uncertainty.csv](../results/stage9/classification_parameter_uncertainty.csv) - CSV data row 1, field `low_95`.

<a id="c04"></a>
## C04

Same-form agreement bootstrap interval A high_95

**Qualification:** Pointwise parameter/link uncertainty conditional on fixed model, population, independent administrations and original rules.

- [results/stage9/classification_parameter_uncertainty.csv](../results/stage9/classification_parameter_uncertainty.csv) - CSV data row 1, field `high_95`.

<a id="c05"></a>
## C05

Same-form agreement bootstrap interval B low_95

**Qualification:** Pointwise parameter/link uncertainty conditional on fixed model, population, independent administrations and original rules.

- [results/stage9/classification_parameter_uncertainty.csv](../results/stage9/classification_parameter_uncertainty.csv) - CSV data row 77, field `low_95`.

<a id="c06"></a>
## C06

Same-form agreement bootstrap interval B high_95

**Qualification:** Pointwise parameter/link uncertainty conditional on fixed model, population, independent administrations and original rules.

- [results/stage9/classification_parameter_uncertainty.csv](../results/stage9/classification_parameter_uncertainty.csv) - CSV data row 77, field `high_95`.

<a id="e02"></a>
## E02

Negative item-rest correlations, one-factor loadings and signed 2PL slopes converge in both forms. No wrong key or cause is established.

**Qualification:** substantive defect, causal mechanism, fairness or operational validity not established

- [results/stage10/claim_traceability.csv](../results/stage10/claim_traceability.csv) - CSV data row 1, field `claim`.

<a id="e06"></a>
## E06

Primary North-South conditional contrasts are positive and screened in each form. Within-form DIF is not cross-form drift or confirmed bias.

**Qualification:** substantive defect, causal mechanism, fairness or operational validity not established

- [results/stage10/claim_traceability.csv](../results/stage10/claim_traceability.csv) - CSV data row 2, field `claim`.

<a id="e07"></a>
## E07

Very weak item-rest alignment and near-flat signed 2PL functioning make the difficulty transformation unstable. This does not establish an incorrect item.

**Qualification:** substantive defect, causal mechanism, fairness or operational validity not established

- [results/stage10/claim_traceability.csv](../results/stage10/claim_traceability.csv) - CSV data row 3, field `claim`.

<a id="e08"></a>
## E08

No independently established competent/not-competent criterion or validated mastery threshold is supplied. Model classifications are not ground truth.

**Qualification:** substantive defect, causal mechanism, fairness or operational validity not established

- [results/stage10/claim_traceability.csv](../results/stage10/claim_traceability.csv) - CSV data row 15, field `claim`.

<a id="e10"></a>
## E10

Delivered scores and classifications reproduce the supplied key and original threshold; external record authenticity and standard validity remain separate.

**Qualification:** certification validity

- [results/stage10/claim_traceability.csv](../results/stage10/claim_traceability.csv) - CSV data row 21, field `claim`.

<a id="mp02"></a>
## MP02

Declared Stage 2 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE2_CONFIG.json](../analysis/configs/STAGE2_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP02`.

<a id="mp03"></a>
## MP03

Declared Stage 3 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE3_CONFIG.json](../analysis/configs/STAGE3_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP03`.

<a id="mp04"></a>
## MP04

Declared Stage 4 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE4_CONFIG.json](../analysis/configs/STAGE4_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP04`.

<a id="mp05"></a>
## MP05

Declared Stage 5 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE5_CONFIG.json](../analysis/configs/STAGE5_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP05`.

<a id="mp06"></a>
## MP06

Declared Stage 6 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE6_CONFIG.json](../analysis/configs/STAGE6_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP06`.

<a id="mp07"></a>
## MP07

Declared Stage 7 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE7_CONFIG.json](../analysis/configs/STAGE7_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP07`.

<a id="mp08"></a>
## MP08

Declared Stage 8 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE8_CONFIG.json](../analysis/configs/STAGE8_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP08`.

<a id="mp09"></a>
## MP09

Declared Stage 9 analytical protocol, definitions and numerical choices

**Qualification:** Protocol choices and numerical checks are conditional; not universal psychometric standards or external validation.

- [analysis/configs/STAGE9_CONFIG.json](../analysis/configs/STAGE9_CONFIG.json) - JSON field `$`.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `MP09`.

<a id="e01"></a>
## E01

This case study uses fully synthetic data and a fictional professional certification program. It was created solely to demonstrate measurement and psychometric methodology. No real organization, examination, candidate or certification decision is represented.

**Qualification:** The delivered fictional dataset supports investigation and internal verification, not real-client outcomes or independent external replication.

- [data/proposed/ANALYST_BRIEF.md](../data/proposed/ANALYST_BRIEF.md) - JSON field `lines:11-13`.

<a id="e03"></a>
## E03

**GO WITH OBSERVATIONS.** A — evidence supports a dominant common dimension, with numerical and item-functioning observations. Separate-form response evidence supports a strong leading dimension; this is not proof of strict unidimensionality or satisfaction of all IRT assumptions. The prespecified sensitivities evaluate whether the conclusion depends on operational omission scoring or the previously flagged common item ANC_14.

**Qualification:** The delivered fictional dataset supports investigation and internal verification, not real-client outcomes or independent external replication.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `E03`.

<a id="e04"></a>
## E04

The reference is conditional on fixed fitted item parameters, the N(0,1) latent distribution and the observed missing mask. Statistical screening prioritizes content/response-process investigation; it cannot identify the cause of association. Neither form exceeds a family-wise threshold in this run. The A maximum is well inside its reference; the B maximum is close to its 95% threshold. Pointwise-only counts are exploratory and compatible with chance exceedances among many pairs; they do not establish systematic excess dependence or confirm any individual pair. Existing negative functioning of ANC_14 and weak discrimination of B_11 remain separate unresolved findings. No response, key, item, raw score or 42/60 rule changed.

**Qualification:** The delivered fictional dataset supports investigation and internal verification, not real-client outcomes or independent external replication.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `E04`.

<a id="e09"></a>
## E09

Statistical strength, substantive strength, severity if confirmed and investigative priority are separate ordinal judgments. No numerical composite risk score is used. These categories express evidence and potential consequence, not probabilities that a defect exists. All current causes and substantive defect statuses remain unresolved.

**Qualification:** The delivered fictional dataset supports investigation and internal verification, not real-client outcomes or independent external replication.

- [results/source_excerpts.json](../results/source_excerpts.json) - JSON field `E09`.
