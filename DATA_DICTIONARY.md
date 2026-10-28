# Proposed synthetic release data dictionary

**Synthetic project data under CC BY 4.0.** The original brief states that all records and the examination are synthetic and fictional. These proposed files are retained without data modification so the delivered manifest and analytical populations can be checked. The author has confirmed project-specific generation under their direction and approved the data license; actual external publication remains separately unauthorized. [E01]

## Files and units

| Proposed file | Unit and contents |
| --- | --- |
| `data/proposed/input/exam_administrations.csv` | One examination administration per row, including repeated candidate histories; response letters, supplied outcomes and background categories. |
| `data/proposed/input/item_metadata.csv` | One bank item per row; form membership, broad domain, keyed alternative and flags. |
| `data/proposed/input/answer_key.csv` | One bank item and keyed alternative per row, supplied separately for reconstruction. |
| `data/proposed/ANALYST_BRIEF.md` | Original fictional operational definitions and synthetic disclosure. |
| `data/proposed/input_manifest.sha256` | Original hashes for the brief and supplied CSVs. It does not authenticate itself or external provenance. |

The dataset has 2,400 administrations, 2,160 first observed candidates and 240 retake administrations. Each form contains 60 administered/scored items; 20 bank items are common. The unchanged passing rule is raw score at least 42 out of 60. These are delivered-data facts, not a validated competence threshold. [Q01] [Q02] [D01] [Q03] [Q04] [Q11]

## Administration fields

| Field | Meaning and handling |
| --- | --- |
| `administration_id` | Synthetic `ADM_` sitting code; unique within delivery. Necessary for response/history joins and reproducible ordering. |
| `candidate_id` | Synthetic `CAN_` person code; repeats link legitimate attempts. Necessary to reconstruct first observed attempts; it is not a real-person identifier. |
| `administration_date` | Fictional ISO calendar date, without time of day. Used to order observed histories. |
| `form` | A or B; drives metadata membership. |
| `attempt_number` | Supplied integer sequence, independently compared with chronological order. |
| `region` | Fictional North/South category. Not inferred ethnicity, legal protection or a causal explanation. |
| `experience_band` | Fictional grouped experience descriptor; broad labels do not supply a complete background-variable codebook. |
| `raw_score` | Supplied number-correct total. It is checked, not overwritten. |
| `percent_score` | Supplied percentage text. Original precision is preserved. |
| `pass_fail` | Supplied PASS/FAIL; independently reconstructed from the original boundary. |
| Item-ID response columns | A, B, C, D or empty string. Membership determines whether an empty field is an administered omission or structural absence. Read with `keep_default_na=False`. |

## Metadata and key fields

`item_id` identifies a bank item. Metadata `form` is AB for common membership or A/B for a form-exclusive item. `domain` is a broad supplied content category, not evidence of latent dimensionality. `correct_option` is the supplied keyed alternative; agreement between metadata and separate key is checked. `anchor_item` and `scored_item` are exported True/False flags, whose designations do not establish invariance or substantive suitability. [E01] [MP03] [MP07]

Operational omissions receive zero credit; structural nonmembers do not enter a form matrix. Primary latent-response matrices treat administered omissions as missing likelihood terms under explicit assumptions. First observed attempts are derived from dates and IDs, not read as an unquestioned analytical flag. No timeout cause, detailed item wording, full group codebook or external competence criterion is included. [MP03] [MP04] [E08]

## Results and provenance

The original field-by-field dictionary is retained as [structured dictionary evidence](results/stage1/data_dictionary.csv), including observed levels, units and provenance. [Selected results](results/README.md) retain source column names and unavailable values. In particular, an empty model difficulty is not zero; a source-defined screen is not confirmed invalidity; linked and predictive values are research quantities. Candidate-level posterior or counterfactual exports are deliberately omitted.

Synthetic identifier syntax, limited categories, fictional dates and the original disclosure support the proposed release assessment. They do not independently prove data origin or authorize redistribution. The author’s rights/data attestation and applied scope are documented in [rights status](docs/RIGHTS_STATUS.md) and [license policy](LICENSE_POLICY.md). Technical inspection does not independently verify the private generating process.

[D01]: docs/claims.md#d01
[E01]: docs/claims.md#e01
[E08]: docs/claims.md#e08
[MP03]: docs/claims.md#mp03
[MP04]: docs/claims.md#mp04
[MP07]: docs/claims.md#mp07
[Q01]: docs/claims.md#q01
[Q02]: docs/claims.md#q02
[Q03]: docs/claims.md#q03
[Q04]: docs/claims.md#q04
[Q11]: docs/claims.md#q11
