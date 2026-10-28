# Professional Technical Inspector Certification

## Data delivery brief

The fictional certification body requests an independent review of the examination
records supplied in this package. The Professional Technical Inspector Certification
is a technical professional credential. Please use the supplied operational records,
item metadata and answer key to establish what the evidence supports and what
additional information the certification body would need to provide.

This case study uses fully synthetic data and a fictional professional certification
program. It was created solely to demonstrate measurement and psychometric methodology.
No real organization, examination, candidate or certification decision is represented.

## Examination operations

- Two operational examination forms are used: A and B.
- Each form contains 60 scored multiple-choice items, with four alternatives per item
  (A, B, C and D) and one keyed correct alternative.
- The forms share 20 common items. Each form also contains 40 exclusive items, giving
  a bank of 100 unique items.
- The five competence domains are Inspection Fundamentals, Safety & Procedures,
  Technical Interpretation, Equipment & Methods and Professional Practice.
- The operational passing rule is at least 42 correct responses out of 60: a raw
  score of 42 or higher is PASS; a raw score below 42 is FAIL.
- An omitted administered response receives zero credit. The denominator remains 60.
- The records cover 18 months, from January 1, 2024 through June 30, 2025. These dates
  are fictional.
- Candidates may take the examination more than once. The same candidate_id links
  their attempts; each administration_id identifies a separate sitting.

## Supplied records

All three CSV files use UTF-8, a comma delimiter, a header row, ISO dates (YYYY-MM-DD)
and a decimal point. Empty response fields mean NA. No spreadsheet software or
connection to another system is required to use these files.

**input/exam_administrations.csv** has one row per examination administration.
Its administrative fields are administration_id, candidate_id, administration_date,
form, attempt_number, region and experience_band. The supplied outcome fields are
raw_score, percent_score and pass_fail. percent_score is on the 0–100 scale. The
remaining 100 columns contain item responses identified by item_id; the values are
A, B, C, D or an empty field. A candidate's attempt_number records their attempt
sequence. Region and experience band are supplied background characteristics.

**input/item_metadata.csv** contains item_id, form, domain, correct_option,
anchor_item and scored_item. The item form label AB identifies a common item;
A or B identifies an item exclusive to that form. The Boolean fields use True/False.
All 60 operational items on each form are scored.

**input/answer_key.csv** contains item_id and correct_option. The key is supplied
separately so that recorded scores can be checked against the response records.

To determine the items administered in a sitting, select metadata records whose
form is either AB or the sitting's form. This identifies 60 administered items.
An empty response for one of those items is an omission and receives zero credit.
An empty response for an item belonging exclusively to the other form is structural
NA: the item was not administered. This rule distinguishes the two types of missing
response without an additional response file. Individual reasons for an omission
are not included in this delivery.

**input_manifest.sha256** lists SHA-256 checksums for the three CSV inputs and this
brief. After extraction, checksums can be verified from the package root using a
compatible SHA-256 utility, for example `sha256sum -c input_manifest.sha256`.

## Questions from the certification body

1. Are the exam data structurally sound?
2. Are the scores reproducible?
3. Are the exam forms being administered consistently?
4. Are there anomalies that require investigation?
5. What evidence is needed to assess the quality of certification decisions?

Please distinguish observations supported by the supplied records from questions
that require additional operational documentation or evidence from the certification
body. The package is self-contained and may be transferred to a separate workspace.
