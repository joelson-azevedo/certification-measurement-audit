# Can You Defend Your Pass/Fail Decisions?

**Joelson Azevêdo · Independent Researcher**

## An Independent Measurement Audit of a Professional Certification Examination

**Independent blind analysis of a synthetic professional certification dataset.** This fictional case represents no real client, examination or certified professional. It demonstrates an audit workflow and its evidentiary limits. The package is a local review candidate, not a published research release. [E01]

## The problem: a correct score can still support an insufficiently defended decision

A certification examination makes a consequential claim when it turns a response record into a pass or fail. The computational rule can be implemented correctly even when important questions about that claim remain unanswered. Does the key represent defensible content? Does the passing standard reflect the required competence? Would comparable performance on another form receive a comparable interpretation? How much instability should be expected near the boundary? Those questions require different evidence; reproducing a total does not answer them collectively. [E10] [E08]

This audit therefore treated scoring integrity as an entry condition rather than a verdict on the whole program. It sought to determine what could be defended from the supplied records, what statistical patterns warranted investigation, and what additional evidence would be needed before changing policy. The central result was deliberately bounded: **scoring was reproducible within the supplied records, while substantive decision validity and operational form comparability remained unestablished.** That conclusion preserves the useful findings without converting uncertainty into a declaration of program failure. [E10] [E08]

The practical concern is not confined to average score reliability. A negatively functioning item can contribute to a score whose overall internal consistency is high. Different groups can have different observed pass rates without the difference being caused solely by their forms. A large population-level agreement estimate can coexist with uncertainty at a consequential boundary. The audit's value lies in separating these questions and linking each proposed action to the evidence it actually requires. [E02] [Q09] [Q10] [Q16] [Q19]

## The delivered examination and the units of analysis

The fictional examination uses Forms A and B, each with 60 scored items and 20 designated common items. The original rule classifies a raw score of 42 or more out of 60 as PASS. The delivery contains 2,400 administrations for 2,160 first observed candidates and 240 retake administrations. Primary item and form investigations use the first observed sitting per candidate: 1,109 on A and 1,051 on B. Those form groups contain distinct candidates. All supplied sittings remain part of the scoring-integrity investigation. [Q03] [Q04] [Q11] [Q01] [Q02] [D01] [D02] [D03]

An administration is a sitting, not a person. An item response belongs to a metadata-defined form, not simply to a column whose name resembles that form. An empty response on an administered item is an omission; an empty response on an item exclusive to the other form is structural absence. Operational scoring gives omissions zero credit and retains the full denominator. Structural absences do not enter that denominator. Preserving this distinction prevents a wide response file from silently becoming an incorrectly defined test. [E01] [Q03]

First observed attempts reduce dependence from repeated candidates in the main analyses. They do not establish that a person had never taken an examination outside the delivery window. The available retakes are selected events, potentially affected by learning, exposure, selection and changing proficiency. They were not treated as an independent validation study of stable repeated decisions. The same discipline applies to form assignment: observed background similarity cannot establish randomized or equivalent-group assignment when its mechanism has not been documented. [MP03] [MP07] [MP09]

## A blind investigation with an auditable chain of evidence

The analyst worked only with the delivered blind-audit package. The supplied brief explicitly describes synthetic records and a fictional program. It provides operational definitions and a checksum manifest, while the investigation is separated from private generation code, hidden parameters and any intended anomaly list. That separation matters because an analyst should discover patterns from observable evidence rather than confirm a known answer. It does not turn synthetic observations into real empirical validation. [E01]

Checksums were verified before substantive processing. Original inputs were preserved, reconstructed scores were stored separately, and later stages recorded their configurations, source hashes, numerical diagnostics and preservation checks. Stages used explicit progression gates: a usable analytical result could proceed with observations while unresolved scientific questions remained open. A successful gate described readiness for the stated research purpose, not certification validity. The curated package retains a manifest, selected original tables and code provenance so readers can examine the chain without navigating the complete internal archive. [E10] [E09]

This distinction also governs the word independent. Scoring was reconstructed separately from the recorded totals. Numerical implementations were checked against alternative identities, controlled fixtures and derivative or integration checks. Those are valuable internal controls. Multiple analyses of the same delivered responses, however, do not create a new external cohort. Candidate bootstrap samples and reference simulations quantify particular uncertainties under stated assumptions; they are not independently observed professionals. [MP04] [MP05] [MP07] [MP09]

## Establishing the computational foundation

The first substantive question was whether the supplied outcomes followed the supplied key and passing rule. The implementation matched the official key to each form's metadata-defined membership, reconstructed credit, rebuilt totals, and independently reapplied the inclusive boundary. Across 2,400 administrations, it found 0 raw-score mismatches and 0 pass/fail mismatches. These are concrete administrative findings. They do not authenticate the delivery against an external authority or establish that the key and standard are substantively correct. [Q01] [Q05] [Q06] [E10]

Preservation was as important as agreement. The audit did not overwrite an original response, key, score or decision. An apparent mismatch would have required diagnosis, with original and reconstructed values retained. Likewise, undefined quantities in later analyses remained unavailable rather than being filled with zero. This rule is particularly important when an apparently simple mathematical transformation is unstable: a finite software output can still be an unsupported measurement interpretation. [E10] [MP04]

Classical analysis then asked how each item aligned with performance on the rest of its form. The corrected item-total correlation excluded the focal response from the total, avoiding the most direct part-whole inflation. Difficulty proportions and distractor choices used explicit administered-item denominators, including operational omission zeros. Reliability used a consistent sample-variance convention so KR-20 and alpha agreed. The resulting alpha values were 0.9015 for A and 0.9101 for B. These describe score consistency in the observed form populations; they do not demonstrate dimensionality, decision accuracy or validity of the passing standard. [MP02] [R01] [R02]

Distractor analysis added a response-process question to the summary statistics. Selection patterns were related to rest-test performance using ties-preserving groups and declared opportunity/count screens. A rare alternative is not automatically defective, and a frequently selected alternative is not automatically misleading. The substantive interpretation depends on the intended skill, wording, alternatives and key rationale, none of which can be recovered from response letters alone. This makes distractors useful evidence for targeted expert review rather than a mechanical deletion rule. [MP02]

## Investigating the measurement model without hiding inconvenient behavior

Internal consistency was not used to justify a latent model by itself. Dimensionality was investigated from binary responses with tetrachoric associations, raw eigenvalue patterns, a matched binary parallel-analysis reference and exploratory common-factor solutions. The primary complete-case samples were 1,015 in A and 957 in B. Complete-case selection, operational-zero sensitivity and same-candidate item exclusion were kept distinct because they answer different missingness and item-set questions. [D05] [D06] [MP03]

The evidence favored a dominant common dimension in both forms, with qualifications. Form A's raw tetrachoric matrix required numerical conditioning for factor analysis. Raw matrices remained the basis for scree and parallel analysis; a separately conditioned copy supported exploratory factor extraction. Additional factors could reduce residual summaries while producing inadmissible communalities. Those solutions were not accepted merely because an optimizer converged or because the metadata listed content domains. A broad dimension was a reasonable provisional hypothesis, not proof of strict unidimensionality. [E03] [MP03]

![Raw tetrachoric parallel analysis for Form A](figures/dimensionality_parallel_analysis_A.png)

*The retained figure compares observed raw tetrachoric roots with exact-marginal binary permutation references for the primary complete-case A sample. The secondary panel shows smaller roots. This is a conditional dimensionality diagnostic, not an estimate of the true number of professional competencies. Figure data and the finite-reference protocol are retained.* [MP03]

The response-model calibration used a signed logistic 2PL, estimating intercepts and discrimination directly. It did not force every item to discriminate positively, reverse a negatively functioning response, or select a model by removing an inconvenient item. Each form was separately identified against a standard-normal proficiency distribution. Consequently, apparently similar proficiency coordinates from the two calibrations were not automatically a common metric. Numerical convergence, identification and substantive interpretation were treated as separate conditions. [MP04]

Quadrature, multiple starts, gradients, boundaries and observed-information behavior were checked explicitly. Initial lower-resolution checks were retained when their EAP differences exceeded the declared tolerance; reporting resolution was increased and independently checked without relaxing that criterion. Parameter uncertainty remained conditional on the response model, latent distribution, omission treatment and independence assumptions. The calibrated model was a tool for examining measurement behavior, not a device for turning missing validity evidence into a positive conclusion. [MP04]

## Three item inquiries, rather than one generic anomaly label

ANC_14 showed negative alignment with the rest of the examination in both forms. Its corrected correlations were -0.3628 in A and -0.4189 in B, accompanied by negative one-factor loadings and signed response-model slopes. The methods use related information from the same delivery, but their consistent direction makes the item a strong investigative priority. The appropriate next step is a documented review of wording, alternatives, official key rationale, content alignment and delivery. The statistics alone do not prove an incorrect key or authorize reverse scoring. [VQ01] [VQ02] [E02]

ANC_18 presented a different question. After conditioning on estimated proficiency from other responses, modeled credit probability was higher for North than South by 11.51 percentage points in A and 11.11 in B, averaged over each form's supported comparison range. Logistic matching excluded the focal response and the previously negative item, retained a primary omnibus multiplicity screen, and withheld unsupported ordinary inference when numerical conditions failed. The result warrants regional interpretation, job-relevance and response-process review. It is not an overall pass-rate difference, demonstrated bias or direct evidence of cross-form drift. [Q07] [Q08] [MP06] [E06]

B_11 had very weak item-rest alignment, with correlation 0.0264, and unstable response-model difficulty. A nearly flat response curve makes the difficulty quotient fragile: small discrimination changes can produce large or misleading location changes. The audit retained that location as unavailable when its support conditions failed. Review should examine the intended competence and retention rationale, rather than treating an unstable quotient as an extreme but precise difficulty or declaring the item automatically invalid. These three cases illustrate why investigative priority should reflect the type of evidence and the question that remains unresolved. [VQ03] [E07] [MP04]

## Testing residual assumptions and conditional regional differences

A dominant dimension does not establish local independence. Residual analysis compared observed responses with probabilities evaluated at the fitted EAP and examined pairwise Q3 within each form. Mean adjustment addressed a known finite-test feature of this residual construction, while a fitted-design reference reproduced the same missing mask, EAP calculation and centering procedure. The primary family-wise threshold came from maxima over the within-form pair family, not from a generic correlation cutoff. [MP05]

No pair exceeded the prespecified primary family-wise threshold in the retained run. That result is narrower than proof of local independence. The reference held item calibration fixed, omitted item-parameter uncertainty and used a finite simulation count. Form B's global maximum was near its screening boundary, with Monte Carlo uncertainty relevant to the tail interpretation. Exploratory pointwise elevations remained review priorities rather than confirmed shared-stimulus mechanisms or a basis for an invented precision correction. [E04] [MP05]

![Fixed-calibration residual reference distributions by form](figures/residual_reference_by_form.png)

*Each panel compares its own observed maximum adjusted Q3 with its own fixed-calibration reference. The missing-response pattern and residual procedure are reproduced. No family-wise exceedance in this run is not proof of response independence; the B result remains close to its screening boundary. The reference maxima and thresholds are packaged separately.* [MP05] [E04]

Regional screening likewise retained its conditioning assumptions. Estimated rest proficiency is a proxy, not an independently observed matching variable. Its parameters were fitted to the same response sample; focal exclusion does not remove that earlier reuse. Common-support restrictions prevented released contrasts in sparse or unsupported regions, and alternate matching and purification were labeled sensitivities. A stable exclusion set cannot establish invariance, just as a small p-value cannot supply the item wording needed to assess a bias hypothesis. [MP06] [E06]

## Comparing forms on a provisional common scale

Observed first-attempt passing rates were 73.40% for A and 61.94% for B. The groups were distinct and unadjusted. Without independent assignment and population evidence, the difference cannot be attributed solely to form difficulty, interpreted as demonstrated unfairness, or used to select a more attractive passing rule. The audit instead asked a conditional research question: what follows if the fitted models can be placed on a defensible provisional common scale? [Q09] [Q10] [MP07]

Common-item selection was evidence-driven and fixed before interpreting the linking results. The full designated pool, a set excluding ANC_14, and a primary set additionally excluding ANC_18 were retained as distinct specifications. The regional exclusion was precautionary; within-form DIF does not itself prove cross-form drift. Excluding an anchor from a linking objective did not delete that item from either full examination model. Broad domain representation remained visible, while detailed content representation and higher-proficiency coverage remained unresolved. [MP07] [E06]

Stocking-Lord and directional Haebara addressed related but distinct discrepancies: aggregate anchor characteristic curves versus item-level probability curves. Both transformed B onto the A metric with positive scale. For the primary set, Stocking-Lord gave scale 1.0398 and shift -0.0160; directional Haebara gave 1.0421 and -0.0179. Their broad compatibility supports a provisional numerical transformation, not independent anchor invariance or operational interchangeability. Leave-one-anchor-out, grid and anchor-set sensitivities reduced some forms of circularity without providing an independent sample. [L01] [L02] [L03] [L04] [MP07]

Scale linking was followed by explicitly defined score correspondence. True-score correspondence pairs expected scores at the same modeled proficiency; observed-score correspondence pairs percentiles of marginal score distributions under one common reference population. The latter used exact Bernoulli score recursion and a declared discrete mid-CDF interpolation convention. It did not replace the marginal distribution with a normal score approximation or pretend that fractional values were attainable raw scores. These methodological distinctions matter even when particular numerical answers happen to be close. [MP08]

Under the primary provisional link and a common standard-normal population on the Form A metric, an A reference score of 42 corresponds to approximately **B=38.70**, with a **95% bootstrap interval 38.15-39.13**. **This is not an authorized conversion or recommended new passing cut.** The interval propagates paired candidate-bootstrap calibration and linking uncertainty with common-item selection fixed. It does not cover uncertain content, assignment or model structure, and it is not an interval for candidate competence. Research-only counterfactual changes therefore cannot be counted as historical certification errors. [Q12] [Q13] [Q14] [Q15] [MP08]

## Agreement is not accuracy, especially near a boundary

Decision consistency was calculated for hypothetical independent administrations at unchanged proficiency. Same-form agreement and cross-form agreement were kept separate. Under the original rules and common A reference, modeled same-form agreement was 89.76% for A (95% bootstrap interval 89.32%-90.27%) and 88.13% for B (87.40%-88.85%). These are conditional model results, not observed agreements from the selected retakes. Memory, learning, exposure or response dependence could alter real repeated performance. [C01] [C02] [C03] [C04] [C05] [C06] [MP09]

Modeled cross-form agreement under original rules was **86.59%**, with a **95% bootstrap interval 85.66%-87.64%**, integrated across the same assumed population. At the fixed modeled proficiency where A's expected score equals 42, agreement was **46.30%**. This is neither the observed group scoring that raw total nor a validated mastery location. The population average includes many proficiency locations far from the transition region; it should not be used as a summary of stability for every individual close to the boundary. [Q16] [Q17] [Q18] [Q19] [Q20]

![Conditional A/B agreement under original and research-only rules](figures/conditional_cross_form_agreement.png)

*The original-rule curve is distinct from the retained research-only boundary sensitivities. The vertical reference is the original modeled A expected-score location, not an observed-score group or competence threshold. The curve shows hypothetical shared-proficiency calculations. Alternative boundaries are not recommendations and were not chosen to maximize agreement.* [MP09] [Q20]

No independently validated competent/not-competent criterion or substantive mastery threshold was supplied. Classification accuracy was therefore not identifiable. A model can predict consistency with its own response assumptions without establishing whether a certification is substantively correct. Neither a narrow bootstrap interval nor a favorable agreement statistic supplies the missing criterion. The audit withheld claims about wrongly accepted or rejected candidates and retained external standard-setting and validation work as prerequisites. [E08]

## Turning findings into a defensible governance sequence

Evidence synthesis produced 20 integrated issues, 27 proposed actions and 14 evidence gaps. Those are organizational outputs, not counts of confirmed defects. Statistical strength, substantive support, possible severity and investigative priority were kept separate. The register avoided a fabricated numerical risk score that could make a strong statistical signal look like a known probability of harm. Missing evidence about the standard received high priority because of its decision relevance, not because the existing standard had been proved wrong. [G01] [G02] [G03] [E09]

The recommended sequence begins with authority and evidence for the passing standard, followed by independent content and regional review, a defensible form-comparison and anchor-validation design, and prospective evidence about actual decision behavior. Actions require different inputs and cannot be completed merely by generating another favorable model output. Future item, scoring or equating proposals need separate authorization after their substantive dependencies are resolved. Nothing in the audit changes the delivered key, responses, original rule or candidate decisions. [E09] [E08]

## Limits, reproduction and the value of the audit

The most consequential limitations are substantive and design-related. Actual item wording, official rationales, standard-setting evidence, delivery/assignment records and independent validation were not supplied. Statistical limitations include finite samples, omitted-response assumptions, provisional latent structure, possible dependence, estimated DIF matching, fixed anchor selection and finite bootstrap tails. Computational stability addresses a different problem: it supports confidence that a declared calculation was implemented consistently, while leaving the adequacy of its assumptions open. [MP04] [MP05] [MP06] [MP07] [MP08] [MP09] [E08]

This curated candidate supports a practical inspection path. Readers can start with the approved executive PDFs, follow claims to retained result cells, inspect the unchanged numerical modules and run a quick verification that reconstructs scoring and selected calculations. The numerical tests use controlled fixtures independent of this case's undisclosed generation process. Full end-to-end replay is explicitly unavailable here because the curated distribution omits the internal stage orchestration and full uncertainty archive. The missing components, original stage order and exact recorded configurations are documented rather than concealed behind a claim of complete reproducibility.

The author has confirmed project-specific synthetic generation under their direction and approved authorship and category licenses. Synthetic status is also supported by the original brief, observed schemas and identifier patterns; those technical observations are distinct from the author’s rights attestation. No private generating truth is included. Separate authorization is still required to publish the exact prepared candidate, and no research correspondence is authorized for operational use. [E01]

The defensible outcome is a bounded evidence program. The audit establishes a reproducible administrative foundation, identifies distinct measurement investigations, quantifies selected conditional uncertainties and explains which stronger conclusions remain unsupported. Its demonstrated capability is the connection between statistical implementation, critical model assessment and decision restraint. Independent measurement auditing is useful precisely when it distinguishes a result that can be computed from a decision that can be justified. [E10] [E08] [E09]

See [methods](METHODOLOGY.md), [limitations](LIMITATIONS.md), [reproduction](REPRODUCIBILITY.md), [source code](analysis/README.md), [evidence index](docs/claims.md) and [review status](docs/RIGHTS_STATUS.md).

[C01]: docs/claims.md#c01
[C02]: docs/claims.md#c02
[C03]: docs/claims.md#c03
[C04]: docs/claims.md#c04
[C05]: docs/claims.md#c05
[C06]: docs/claims.md#c06
[D01]: docs/claims.md#d01
[D02]: docs/claims.md#d02
[D03]: docs/claims.md#d03
[D05]: docs/claims.md#d05
[D06]: docs/claims.md#d06
[E01]: docs/claims.md#e01
[E02]: docs/claims.md#e02
[E03]: docs/claims.md#e03
[E04]: docs/claims.md#e04
[E06]: docs/claims.md#e06
[E07]: docs/claims.md#e07
[E08]: docs/claims.md#e08
[E09]: docs/claims.md#e09
[E10]: docs/claims.md#e10
[G01]: docs/claims.md#g01
[G02]: docs/claims.md#g02
[G03]: docs/claims.md#g03
[L01]: docs/claims.md#l01
[L02]: docs/claims.md#l02
[L03]: docs/claims.md#l03
[L04]: docs/claims.md#l04
[MP02]: docs/claims.md#mp02
[MP03]: docs/claims.md#mp03
[MP04]: docs/claims.md#mp04
[MP05]: docs/claims.md#mp05
[MP06]: docs/claims.md#mp06
[MP07]: docs/claims.md#mp07
[MP08]: docs/claims.md#mp08
[MP09]: docs/claims.md#mp09
[Q01]: docs/claims.md#q01
[Q02]: docs/claims.md#q02
[Q03]: docs/claims.md#q03
[Q04]: docs/claims.md#q04
[Q05]: docs/claims.md#q05
[Q06]: docs/claims.md#q06
[Q07]: docs/claims.md#q07
[Q08]: docs/claims.md#q08
[Q09]: docs/claims.md#q09
[Q10]: docs/claims.md#q10
[Q11]: docs/claims.md#q11
[Q12]: docs/claims.md#q12
[Q13]: docs/claims.md#q13
[Q14]: docs/claims.md#q14
[Q15]: docs/claims.md#q15
[Q16]: docs/claims.md#q16
[Q17]: docs/claims.md#q17
[Q18]: docs/claims.md#q18
[Q19]: docs/claims.md#q19
[Q20]: docs/claims.md#q20
[R01]: docs/claims.md#r01
[R02]: docs/claims.md#r02
[VQ01]: docs/claims.md#vq01
[VQ02]: docs/claims.md#vq02
[VQ03]: docs/claims.md#vq03
