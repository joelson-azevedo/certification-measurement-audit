# Methodology and implementation appendix

**Joelson Azevêdo · Independent Researcher**

**Independent blind analysis of a synthetic professional certification dataset.** These methods describe the delivered fictional case and its approved analyses. Implementations are research code with numerical safeguards, not an operational certification service. The [code guide](analysis/README.md) separates unchanged numerical modules from release-verification orchestration. [E01]

## Population, scoring and missingness

Primary item/model investigations use one first observed administration per candidate, separately within each form. Chronological histories are reconstructed independently of file order and the supplied attempt sequence. This reduces repeated-person dependence, but does not rule out pre-window attempts. Retakes are retained for administrative reconstruction and explicitly appropriate descriptive context; their selection and potential proficiency change preclude treating them as independent external retest validation. [MP03] [MP07] [MP09]

Form membership is determined by metadata: common items plus items exclusive to the delivered form. Correct administered responses receive credit; operational omissions receive zero, with denominator 60 and original passing boundary 42. Structural nonmembers do not become omission zeros. Original outcomes are retained alongside independently reconstructed values. Percentage comparison uses decimal arithmetic with an explicitly declared export-precision tolerance, rather than changing the raw-score boundary through percentage rounding. [Q03] [Q11] [E10]

The primary latent-response analyses treat administered omissions as missing response likelihood terms, not observed incorrect answers. This assumes the missingness mechanism can be ignored for the declared model inference; the delivery does not establish MCAR, MAR or an omission cause. Complete-case dimensionality, all-first-attempt operational-zero sensitivities and same-candidate missing-aware IRT sensitivities are different designs. Similar results cannot establish that those missingness mechanisms are equivalent. [MP03] [MP04]

## Item statistics, reliability and distractors

Difficulty is the correct proportion among administered candidates, including operational omission zeros. Wilson intervals use z=1.959963984540054 and describe itemwise uncertainty under independent first candidates. Corrected discrimination is corr(item, total minus item); constant item/rest variance and insufficient observations remain undefined. Alpha and KR-20 use consistent sample-variance conventions: the pq sum is multiplied by n/(n-1) when paired with the sample total variance. Deletion alpha recalculates the reduced rest-score variance. Values are not clipped into a conventional range. [MP02]

Alpha uncertainty uses 2,000 within-form candidate-row resamples, seeds 20261021 for A and 20261022 for B, and empirical percentile endpoints. Items and the observed form population remain fixed in interpretation. Secondary all-administration summaries do not attach intervals that treat repeat candidates as independent. Neither alpha nor KR-20 establishes dimensionality, local independence, classification consistency at the boundary or substantive validity. [MP02]

Difficulty, discrimination and option screens are declared exploratory conventions, not universal standards. Distractors use 5 ties-preserving performance groups, minimum band and option counts, and an explicit lower-band error-opportunity condition. Alternative selection is related to rest performance without including the focal response in the grouping total. A rare option does not automatically fail; available opportunity, conditional selection and content intent matter. The exact thresholds are retained in the configuration and source. [MP02]

## Binary dimensionality and numerical conditioning

The primary dimensionality analysis uses complete cases on the full form. The ANC_14 exclusion sensitivity retains those same people; the operational-zero sensitivity instead retains all first candidates, changing coding and selection together. Tetrachoric associations assume a latent bivariate-normal threshold representation. Marginal thresholds and the joint probability determine each correlation. Constant margins are unidentifiable rather than rescued by continuity correction. [MP03]

When any joint cell is zero, 0.5 is added to all four cells and thresholds are recalculated. Positive cells below 5 are flagged, not automatically shifted. The deterministic Plackett-integral calculation uses transformed 32-node Gauss-Legendre quadrature with doubled-node verification, safeguarded inversion and a documented adaptive fallback. Each pair retains its original counts, correction, convergence and numerical checks. A finite matrix condition number is not proof of positive definiteness. [MP03]

Raw tetrachoric matrices remain the inputs for scree and binary parallel analysis. A separate EFA copy may be conditioned using Higham positive-semidefinite/unit-diagonal projections with Dykstra correction, an eigenvalue floor of 1e-6 and diagonal rescaling. Maximum absolute adjustment and relative Frobenius adjustment must each remain at most 0.05. These are numerical safeguards, not psychometric standards. Conditioning changes the estimated matrix and can affect secondary structure even when its overall norm is small. [MP03]

Parallel analysis independently permutes binary columns, preserving exact observed marginals and sample size, and applies the same raw tetrachoric estimator. It uses 100 valid references and the order-matched 95th percentile. The seed rule is 20261023 plus 1000 times form index plus 100 times specification index. Leading consecutive exceedances are retained until the first nonexceedance; later isolated exceedances do not silently add dimensions. The reference is not a continuous-normal or mixed-estimator substitute. [MP03]

Common-factor extraction uses ULS/MINRES uniqueness optimization rather than PCA loadings, with deterministic starts and oblique Promax for multifactor solutions. Optimizer convergence and admissibility are separate. Heywood cases, lower-bound hits, communalities and factor correlations remain visible. Eigenvalue trace fractions are PCA-style descriptions, not common-factor variance explained. Broad content domains do not establish distinct latent traits or defensible domain subscores. [MP03]

## Signed-slope 2PL calibration and diagnostic fit

The response model is P(X=1|theta)=logistic(d+a theta), with logistic scale D=1 and genuinely signed discrimination a. Each form independently fixes theta to a standard normal distribution. Positive aggregate slope fixes global reflection; it does not reverse any individual response. The direct intercept/slope representation preserves weak, zero and negative behavior. Difficulty b=-d/a is withheld when |a| is below 0.1 or the supported slope interval includes zero. A negative-slope crossing location is not ordinary increasing-item difficulty. [MP04]

Marginal maximum likelihood uses stable log probabilities, log-sum-exp integration, analytic gradients and deterministic multiple starts. Initial 31/41-node checks were retained; final reporting used 81 nodes checked at 101 after the unchanged EAP tolerance required escalation. Slope/intercept safeguards are respectively +/-6 and +/-10. Supported convergence requires the projected gradient per candidate to be no more than 1e-6 in addition to optimizer success. Boundaries and failed fits remain explicit. [MP04]

Observed marginal information is computed with the Louis identity, including posterior score covariance. Symmetry, eigenvalues, condition number and directional derivative checks precede a supported covariance calculation. Unsupported information or boundary behavior withholds uncertainty rather than triggering an unquestioned matrix inverse. Conditional item intervals omit model misspecification, unresolved content, selection and later linking uncertainty. Posterior proficiency SD and inverse-square-root information SEM are different quantities. [MP04]

Descriptive item checks use focal-response-excluded posterior grouping and integrated response predictions; original calibration still reuses the sample. Fixed-parameter model-replicated envelopes use 200 response/latent samples, seeds beginning at 20261024, and the observed missing mask. They include response/latent sampling, not item-parameter uncertainty, and are not calibrated item-fit p-values or external validation. Better fit after adding flexibility does not establish a valid passing interpretation. [MP04]

## Residual local-dependence screening

Residuals are observed credit minus item probability at the saved EAP. Raw Q3 correlates residual columns on pair-specific observed subsets. Adjusted Q3 subtracts the unweighted mean of unique finite off-diagonal pairs within the same form/specification. The mean being zero after adjustment is a construction, not a fit result. Focal response reuse and finite-test shrinkage remain relevant. [MP05]

The primary fixed-calibration reference targets 500 valid main simulations after a separate pilot, with attempt ceilings and retained failures. Main seeds are 20261125 for A and 20261126 for B; attempt-specific SeedSequence streams reproduce each draw. New standard-normal proficiencies and conditionally independent responses are generated at the saved item parameters, the original missing mask is restored, EAP is recalculated at 81 nodes, and the complete residual procedure is repeated. This is not a refitting bootstrap. [MP05]

The positive maximum across 1,770 pairs within each form defines its own family-wise 95% and 99% reference thresholds using the fixed higher order-statistic convention. Exceedance is strict; no combined cross-form error-control claim is made. Plus-one empirical tail areas, finite-reference resolution and Monte Carlo uncertainty are retained. Exploratory pairwise elevations and leave-two-out joint-cell diagnostics are not confirmed dependence causes. No unsupported LD-X2, dependence-adjusted SEM or testlet remedy is invented. [MP05]

## Regional DIF: matching, support and multiplicity

South is the reference and North the focal category. Forms are analyzed separately. Primary rest EAP uses fixed signed calibrations, available other responses and an 81-node posterior, excluding both the focal item and ANC_14. The matching variable is centered over matching-eligible candidates in the form, including omitted-target candidates; target-valid rows are then used consistently across nested regressions. Focal exclusion does not remove original calibration reuse or establish an invariant matching set. [MP06]

Nested logits start with matching proficiency, add a group term, and then an interaction. The primary discovery test is the omnibus comparison; additive and interaction increments support interpretation. Ordinary unpenalized ML requires convergence, full rank, supported information and no complete/quasi-complete separation. Unsupported fits are unavailable; penalized estimates are not inserted into ordinary likelihood-ratio tests. Deterministic scaling/backtracking and explicit failure states are part of the implementation. [MP06]

BH q=0.05 applies to 60 primary omnibus tests per form; a 120-test cross-form family is a separate sensitivity. It is not pooled calibration or an additional unadjusted discovery route. Reference support intersects group 5th-95th percentile ranges, uses 10 bins and requires at least 20 target responses per group in a bin. Integrated effects use an equal-group-weighted empirical matching distribution restricted to supported bins. Each form/item therefore has its own supported reference. [MP06]

Wald/delta-method intervals condition on estimated matching and fixed item parameters; they are not simultaneous or matching-uncertainty-adjusted. With interaction, the group coefficient and additive odds ratio do not describe a constant effect everywhere. Bounded cumulative purification, prorated rest-score matching, operational-zero matching and uniform-only Mantel-Haenszel are sensitivities. Stable purification is not proof of invariance. DIF remains a conditional statistical signal requiring substantive item and response-process review. [MP06]

## Anchor assessment and scale linking

The fixed anchor specifications are S20, S19 excluding ANC_14, and primary S18 additionally excluding ANC_18. The latter precaution reflects regional evidence, not a cross-form drift test. Broad content coverage, curve discrepancies, signed source behavior, supported proficiency overlap and leave-one-anchor-out behavior are examined together. Easier retained anchors and limited upper-range information restrict interpretation; broad-domain representation is not detailed-content validation. Full-form response models keep every original item. [MP07]

The direction is theta_A=scale_A theta_B+shift_B, with positive scale. Therefore a_B_on_A=a_B/scale_A and d_B_on_A=d_B-a_B shift_B/scale_A. These expressions remain defined for signed or near-zero slopes. The difficulty transformation is reported only where the original difficulty is supported. Information and posterior dispersion transformations have their own scale factors. Probability preservation and inverse round trips are checked directly. [MP07]

Stocking-Lord minimizes weighted anchor test-characteristic-curve differences; directional Haebara minimizes item probability differences on the A grid. The primary grid has 601 nodes on [-3,3] with normalized standard-normal trapezoidal weights; fine and broad grids are fixed sensitivities. Multiple starts, log-scale positivity, derivative checks and interior solutions govern numerical support. Agreement between methods does not establish invariance or justify a universal preferred method. [MP07]

Primary uncertainty uses 200 valid paired, within-form candidate bootstraps with full signed-model recalibration and fixed S18 links, seed 20261127. Secondary uncertainty uses 1,000 paired full-covariance normal parameter draws, seed 20261227, preserving within-form intercept/slope and cross-item covariance. The forms use independent streams while each pair and both linking methods remain coupled. Failure ceilings and every attempt are retained. Selection remains fixed, so these intervals do not include anchor-selection uncertainty or substantive/design uncertainty. [MP07]

## True-score and observed-score correspondence

True-score correspondence solves expected-score functions at the same linked proficiency. Observed-score correspondence instead integrates exact conditional Bernoulli distributions over the same declared population for both forms. The primary reference is standard normal on the A metric; transformed B-identification and equal-mixture populations are sensitivities, not empirically verified population distributions. Lord-Wingersky recursion supports signed response probabilities without imposing positive item slopes. Conditional and marginal probability, moment and total-variance identities are checked. [MP08]

The discrete convention uses integer mid-CDF ranks with linear interpolation between knots; B-to-A maps the B rank through the inverse A mid-CDF. Flat/nonunique intervals and out-of-range tails remain unavailable rather than clipped. Fractional values are continuous research correspondences, not observed integer scores. Separate floor, nearest-half-up and ceiling conventions illustrate quantization dependence without recommending an operational boundary. [MP08]

Distribution integration checks 81, 161, 321 and if needed 641 nodes against prespecified probability and mapping tolerances; this is separate from the saved calibration quadrature. Paired uncertainty reuses the saved calibration/link replicates without relinking or additional random draws. Bootstrap intervals remain pointwise, conditional on the fixed anchor set and model/reference assumptions. Counterfactual aggregate classifications hold the original impact sample fixed; they neither identify historical errors nor measure real repeated-test behavior. [MP08]

## Consistency, accuracy and validity

At fixed proficiency, same-form hypothetical agreement is p squared plus (1-p) squared. Cross-form agreement is p_A p_B+(1-p_A)(1-p_B). Joint cells are integrated over the common population; multiplying marginal pass rates would discard shared proficiency. Original-rule and research-only scenarios remain separate. The fixed expected-score reference is not the observed raw-score group, whose predictions would require a response-pattern posterior. [MP09]

The calculations assume independent responses conditional on fixed proficiency, no practice/memory/exposure change, adequate response models and provisional linking. Exact score tails and quadrature convergence support the computation, not those assumptions. Bootstrap intervals describe parameter/link estimation uncertainty, distinct from posterior proficiency uncertainty, future Bernoulli variation or candidate-to-candidate response-pattern differences. Repeating the same items independently is an especially strong hypothetical assumption. [MP09]

Accuracy requires an independently defensible competence criterion or mastery threshold. None was supplied. Accordingly, no substantively correct-certification percentage or wrongly accepted/rejected count is estimated. Validity additionally concerns the evidence supporting the intended interpretation and use. Consistency can contribute to that argument, but cannot replace content, standard-setting, assignment and external decision evidence. [E08]

## Evidence integration and responsible use

The synthesis distinguishes administrative facts, statistical diagnostics, model-dependent estimates, substantive interpretations requiring external evidence, and proposed investigations. Priority, potential severity and certainty of a defect remain separate. Same-dataset convergence, internal numerical checking, candidate bootstrap uncertainty and external empirical replication are explicitly different forms of evidence. [E09]

This appendix points to [packaged claim sources](docs/claims.md), [exact configurations](analysis/README.md), [selected results](results/README.md), [reproduction tiers](REPRODUCIBILITY.md) and [limits](LIMITATIONS.md). Full original stage reports remain controlled; selected literal methodological passages are retained in `results/source_excerpts.json` so no inaccessible local path is needed to inspect the explanations used here.

[E01]: docs/claims.md#e01
[E08]: docs/claims.md#e08
[E09]: docs/claims.md#e09
[E10]: docs/claims.md#e10
[MP02]: docs/claims.md#mp02
[MP03]: docs/claims.md#mp03
[MP04]: docs/claims.md#mp04
[MP05]: docs/claims.md#mp05
[MP06]: docs/claims.md#mp06
[MP07]: docs/claims.md#mp07
[MP08]: docs/claims.md#mp08
[MP09]: docs/claims.md#mp09
[Q03]: docs/claims.md#q03
[Q11]: docs/claims.md#q11
