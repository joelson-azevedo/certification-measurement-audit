"""Validate public delivery only. Never edit inputs or silently repair records."""

from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path

import pandas as pd

ADMIN_FIELDS = [
    "administration_id", "candidate_id", "administration_date", "form",
    "attempt_number", "region", "experience_band", "raw_score", "percent_score", "pass_fail",
]
META_FIELDS = ["item_id", "form", "domain", "correct_option", "anchor_item", "scored_item"]
KEY_FIELDS = ["item_id", "correct_option"]
DOMAINS = [
    "Inspection Fundamentals", "Safety & Procedures", "Technical Interpretation",
    "Equipment & Methods", "Professional Practice",
]
OPTIONS = {"A", "B", "C", "D"}
SUPPLIED = ["ANALYST_BRIEF.md", "input/answer_key.csv",
            "input/exam_administrations.csv", "input/item_metadata.csv"]
EXCEPTION_FIELDS = [
    "exception_id", "severity", "category", "affected_record", "affected_variable",
    "observed_value", "expected_rule", "description", "proposed_action", "resolution_status",
]


class Findings:
    def __init__(self):
        self.rows = []

    def add(self, severity, category, record, variable, observed, rule, description,
            action="Request clarification or corrected delivery from the certification body.",
            status="OPEN"):
        self.rows.append(dict(zip(EXCEPTION_FIELDS, [
            f"DQ-{len(self.rows) + 1:04d}", severity, category, str(record), variable,
            str(observed), rule, description, action, status,
        ])))

    @property
    def critical(self):
        return any(r["severity"] == "CRITICAL" for r in self.rows)

    def frame(self):
        return pd.DataFrame(self.rows, columns=EXCEPTION_FIELDS)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_manifest(root, findings):
    """Validate manifest coverage and hashes before reading substantive records."""
    root = Path(root).resolve()
    rows, seen = [], set()
    manifest = root / "input_manifest.sha256"
    if not manifest.is_file():
        findings.add("CRITICAL", "integrity", "manifest", "sha256", "absent",
                     "Manifest must exist", "Cannot authenticate the delivery.")
        return pd.DataFrame(columns=["file", "expected_sha256", "actual_sha256", "match"])
    for n, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        m = re.fullmatch(r"([0-9a-fA-F]{64})\s+\*?(.+)", line)
        if not m:
            findings.add("CRITICAL", "integrity", f"manifest:{n}", "sha256", line,
                         "64 hexadecimal characters followed by a relative path",
                         "Malformed manifest entry.")
            continue
        expected, name = m.groups()
        name = name.replace("\\", "/")
        path = (root / name).resolve()
        if name not in SUPPLIED or not path.is_relative_to(root) or name in seen:
            findings.add("CRITICAL", "integrity", name, "path", name,
                         "Exactly one manifest entry for each of the four documented inputs",
                         "Unexpected, duplicated or unsafe manifest path; not opened.")
            continue
        seen.add(name)
        actual = sha256(path) if path.is_file() else "MISSING"
        match = actual.lower() == expected.lower()
        rows.append({"file": name, "expected_sha256": expected.lower(),
                     "actual_sha256": actual, "match": match})
        if not match:
            findings.add("CRITICAL", "integrity", name, "sha256", actual, expected,
                         "Input hash does not match; stop substantive processing.")
    for name in sorted(set(SUPPLIED) - seen):
        findings.add("CRITICAL", "integrity", name, "sha256", "unverified",
                     "All four documented inputs must be verified", "Manifest coverage incomplete.")
    return pd.DataFrame(rows, columns=["file", "expected_sha256", "actual_sha256", "match"])


def load_csv(path, findings):
    """Inspect raw headers/row widths, then preserve all cell text, including empty strings."""
    with Path(path).open(encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, [])
        if not header or len(header) != len(set(header)):
            findings.add("CRITICAL", "schema", Path(path).name, "header", header,
                         "Nonempty, unique column names", "Missing or duplicated CSV headers.")
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                findings.add("CRITICAL", "schema", f"{Path(path).name}:{line}", "row_width", len(row),
                             str(len(header)), "CSV row width differs from its header.")
    if findings.critical:
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")


def valid_identifier(value):
    """Minimal lexical check, not an invented ADM_/CAN_ numbering contract."""
    return isinstance(value, str) and bool(value) and not any(c.isspace() or ord(c) < 32 for c in value)


def identity_checks(frame, field, findings, label):
    for i, value in frame[field].items():
        if not valid_identifier(value):
            findings.add("CRITICAL", "identity", f"{label}:row{i + 2}", field, repr(value),
                         "Nonempty identifier without whitespace or control characters",
                         "Identifier is missing or lexically invalid; naming pattern otherwise unspecified.")
    duplicates = frame[field].duplicated(keep=False)
    for value in frame.loc[duplicates, field].unique():
        findings.add("CRITICAL", "identity", value, field, int((frame[field] == value).sum()),
                     "One record per unique identifier", f"Repeated {label} identifier.")


def form_items(metadata, form):
    return metadata.loc[metadata["form"].isin([form, "AB"]), "item_id"].tolist()


def validate_structure(admin, meta, key, findings, items_per_form=60):
    """Return check-level evidence; exception register records every failed rule."""
    checks = []

    def check(name, ok, observed, expected, severity="CRITICAL", category="structure"):
        checks.append({"check": name, "passed": bool(ok), "observed": str(observed), "expected": str(expected)})
        if not ok:
            findings.add(severity, category, "delivery", name, observed, str(expected), name)

    check("administration fields", set(ADMIN_FIELDS).issubset(admin.columns), list(admin.columns[:10]), ADMIN_FIELDS)
    check("metadata schema", list(meta.columns) == META_FIELDS, list(meta.columns), META_FIELDS)
    check("key schema", list(key.columns) == KEY_FIELDS, list(key.columns), KEY_FIELDS)
    if findings.critical:
        return pd.DataFrame(checks)
    check("nonempty administrations", len(admin) > 0, len(admin), ">0")
    identity_checks(admin, "administration_id", findings, "administration")
    identity_checks(meta, "item_id", findings, "metadata item")
    identity_checks(key, "item_id", findings, "answer-key item")
    invalid_candidates = ~admin.candidate_id.map(valid_identifier)
    for i in admin.index[invalid_candidates]:
        findings.add("CRITICAL", "identity", admin.at[i, "administration_id"], "candidate_id",
                     repr(admin.at[i, "candidate_id"]), "Nonempty identifier without whitespace/control characters",
                     "Candidate history cannot be linked reliably.")
    check("duplicate complete records", not admin.duplicated().any(), int(admin.duplicated().sum()), 0)
    response_columns = set(admin.columns) - set(ADMIN_FIELDS)
    check("response columns match metadata", response_columns == set(meta.item_id),
          f"missing={sorted(set(meta.item_id)-response_columns)}, extra={sorted(response_columns-set(meta.item_id))}",
          "Exact item identifier set")
    check("key coverage", set(key.item_id) == set(meta.item_id),
          f"missing={sorted(set(meta.item_id)-set(key.item_id))}, extra={sorted(set(key.item_id)-set(meta.item_id))}",
          "Exactly one key for every metadata item")
    check("administration forms", set(admin.form).issubset({"A", "B"}), sorted(admin.form.unique()), "A or B")
    check("metadata form assignments", set(meta.form).issubset({"A", "B", "AB"}), sorted(meta.form.unique()), "A, B or AB")
    check("bank size", len(meta) == 100, len(meta), 100)
    check("common item count", int((meta.form == "AB").sum()) == 20, int((meta.form == "AB").sum()), 20)
    check("domain assignments", set(meta.domain) == set(DOMAINS), sorted(meta.domain.unique()), DOMAINS)
    for field in ["anchor_item", "scored_item"]:
        check(f"valid {field} Boolean", meta[field].isin(["True", "False"]).all(), sorted(meta[field].unique()), "True/False")
    check("all operational items scored", meta.scored_item.eq("True").all(), sorted(meta.scored_item.unique()), "True")
    check("anchors assigned to common items", (~meta.anchor_item.eq("True") | meta.form.eq("AB")).all(),
          meta.loc[meta.anchor_item.eq("True") & ~meta.form.eq("AB"), "item_id"].tolist(), "Anchors belong to both forms")
    check("common item anchor designation", meta.loc[meta.form.eq("AB"), "anchor_item"].eq("True").all(),
          meta.loc[meta.form.eq("AB"), "anchor_item"].value_counts().to_dict(), "Common items designated as anchors",
          severity="MAJOR", category="metadata")
    for table, label in [(meta, "metadata"), (key, "answer key")]:
        check(f"valid {label} alternatives", table.correct_option.isin(OPTIONS).all(),
              sorted(table.correct_option.unique()), "A/B/C/D")
    if not meta.item_id.duplicated().any() and not key.item_id.duplicated().any():
        joined = meta.merge(key, on="item_id", how="outer", suffixes=("_metadata", "_key"))
        bad = joined.correct_option_metadata != joined.correct_option_key
        check("metadata and answer key agree", not bad.any(), joined.loc[bad, "item_id"].tolist(), "Identical correct options")
    for form in ["A", "B"]:
        ids = form_items(meta, form)
        check(f"form {form} item count", len(ids) == items_per_form, len(ids), items_per_form)
        check(f"form {form} exclusive count", int(meta.form.eq(form).sum()) == 40, int(meta.form.eq(form).sum()), 40)
    dates = pd.to_datetime(admin.administration_date, format="%Y-%m-%d", errors="coerce")
    canonical = admin.administration_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}") & dates.notna()
    in_period = dates.between("2024-01-01", "2025-06-30")
    for i in admin.index[~canonical | ~in_period]:
        findings.add("CRITICAL", "dates", admin.at[i, "administration_id"], "administration_date",
                     admin.at[i, "administration_date"], "ISO date in 2024-01-01 through 2025-06-30",
                     "Invalid, missing or out-of-period administration date.")
    checks.append({"check": "dates valid and within period", "passed": bool((canonical & in_period).all()),
                   "observed": int((~canonical | ~in_period).sum()), "expected": "0 invalid dates"})
    for field in ["region", "experience_band"]:
        blank = admin[field].str.strip().eq("")
        if blank.any():
            findings.add("MAJOR", "missingness", "delivery", field, int(blank.sum()), "Background data available",
                         "Missing group characteristic; restrict later group analyses pending clarification.")
    for field, lower, upper, integral in [("raw_score", 0, 60, True), ("percent_score", 0, 100, False),
                                          ("attempt_number", 1, float("inf"), True)]:
        numeric = pd.to_numeric(admin[field], errors="coerce")
        valid = numeric.between(lower, upper)
        if integral:
            valid &= numeric.mod(1).eq(0)
        for i in admin.index[~valid]:
            findings.add("CRITICAL", "values", admin.at[i, "administration_id"], field, admin.at[i, field],
                         f"{'Integer' if integral else 'Number'} in [{lower}, {upper}]", "Invalid numeric value.")
    check("pass/fail vocabulary", admin.pass_fail.isin(["PASS", "FAIL"]).all(), sorted(admin.pass_fail.unique()), "PASS or FAIL")
    if response_columns == set(meta.item_id):
        validate_responses(admin, meta, findings)
    return pd.DataFrame(checks)


def validate_responses(admin, meta, findings):
    """Blank nonmembers are structural NA; any marked nonmember is a mapping failure."""
    for form in admin.form.unique():
        members = set(form_items(meta, form))
        subset = admin.loc[admin.form.eq(form)]
        for item in meta.item_id:
            values = subset[item]
            bad = ~values.isin(OPTIONS | {""}) if item in members else values.ne("")
            for i in values.index[bad]:
                findings.add("CRITICAL", "responses" if item in members else "mapping",
                             admin.at[i, "administration_id"], item, values.at[i],
                             "A/B/C/D or empty omission" if item in members else "Empty: item not administered",
                             "Unexpected answer option." if item in members else "Response appears in an unexpected form.")


def reconstruct_histories(admin, findings):
    """Infer sequence from dates, then compare with supplied attempt_number independently."""
    history = admin[ADMIN_FIELDS[:7]].copy()
    history["_date"] = pd.to_datetime(history.administration_date, format="%Y-%m-%d", errors="coerce")
    history = history.sort_values(["candidate_id", "_date", "administration_id"], kind="stable")
    history["audited_attempt_number"] = history.groupby("candidate_id").cumcount() + 1
    history["is_first_attempt"] = history.audited_attempt_number.eq(1)
    history["retake_interval_days"] = history.groupby("candidate_id")._date.diff().dt.days.astype("Int64")
    tie = history.duplicated(["candidate_id", "administration_date"], keep=False)
    for candidate in history.loc[tie, "candidate_id"].unique():
        findings.add("CRITICAL", "history", candidate, "administration_date", "same-day sittings",
                     "Unambiguous chronological sequence", "Date-only records cannot establish within-day order.")
    original = pd.to_numeric(history.attempt_number, errors="coerce")
    for i in history.index[original.ne(history.audited_attempt_number)]:
        findings.add("CRITICAL", "history", history.at[i, "administration_id"], "attempt_number",
                     history.at[i, "attempt_number"], str(history.at[i, "audited_attempt_number"]),
                     "Supplied attempt does not match chronological sequence within the delivered window.")
    attributes = history.groupby("candidate_id")[["region", "experience_band"]].nunique(dropna=False)
    for field in attributes:
        for candidate in attributes.index[attributes[field] > 1]:
            findings.add("MAJOR", "history", candidate, field,
                         history.loc[history.candidate_id.eq(candidate), field].tolist(),
                         "Investigate changing background attributes", "Attribute differs across candidate attempts; may reflect legitimate change.")
    return history.drop(columns="_date")
