"""Independent key matching, administered-item denominators, and discrepancy detection."""

from __future__ import annotations

from decimal import Decimal, InvalidOperation

import pandas as pd

from .data_integrity import OPTIONS, form_items

# Not documented by the client. A declared audit tolerance of half the observed
# ten-decimal export unit plus binary conversion slack; never relative tolerance.
PERCENT_ABSOLUTE_TOLERANCE = 5.1e-11


def pass_fail(raw):
    return "PASS" if raw >= 42 else "FAIL"


def percent_matches(original, raw, denominator=60):
    try:
        value = Decimal(str(original))
        if not value.is_finite():
            return False
        expected = Decimal(int(raw)) * Decimal(100) / Decimal(denominator)
        return abs(value - expected) <= Decimal(str(PERCENT_ABSOLUTE_TOLERANCE))
    except (InvalidOperation, ValueError, TypeError):
        return False


def make_long(admin, meta, key, history):
    """Exactly one row for each metadata member, including all administered omissions."""
    parts = []
    for form in sorted(admin.form.unique()):
        subset = admin.loc[admin.form.eq(form)]
        ids = form_items(meta, form)
        identifiers = ["administration_id", "candidate_id", "form", "attempt_number", "administration_date"]
        part = subset.melt(id_vars=identifiers, value_vars=ids, var_name="item_id", value_name="response")
        parts.append(part)
    long = pd.concat(parts, ignore_index=True)
    long = long.merge(meta[["item_id", "domain", "anchor_item", "scored_item"]], on="item_id", validate="many_to_one")
    long = long.merge(key, on="item_id", validate="many_to_one")
    long = long.merge(history[["administration_id", "is_first_attempt"]], on="administration_id", validate="many_to_one")
    # Unexpected marked options are invalid input, not silently transformed to omissions.
    if not long.response.isin(OPTIONS | {""}).all():
        raise ValueError("Unexpected response option: scoring blocked")
    long["omitted"] = long.response.eq("")
    long["anchor_item"] = long.anchor_item.eq("True")
    long["scored_item"] = long.scored_item.eq("True")
    long["item_score"] = (long.response.eq(long.correct_option) & long.scored_item).astype("int8")
    if long.duplicated(["administration_id", "item_id"]).any():
        raise ValueError("Duplicate administered-item observational unit")
    return long


def reconstruct_scores(admin, long, history, denominator=60):
    totals = long.loc[long.scored_item].groupby("administration_id").item_score.sum()
    scored = admin.copy()
    scored["audited_raw_score"] = scored.administration_id.map(totals)
    if scored.audited_raw_score.isna().any():
        raise ValueError("No independently reconstructed score for an administration")
    scored["audited_raw_score"] = scored.audited_raw_score.astype(int)
    scored["audited_percent_score"] = scored.audited_raw_score * 100 / denominator
    scored["audited_pass_fail"] = scored.audited_raw_score.map(pass_fail)
    flags = history.set_index("administration_id")
    scored["is_first_attempt"] = scored.administration_id.map(flags.is_first_attempt)
    scored["audited_attempt_number"] = scored.administration_id.map(flags.audited_attempt_number)
    return scored


def compare_scores(scored):
    fields = ["administration_id", "candidate_id", "form", "raw_score", "percent_score", "pass_fail",
              "audited_raw_score", "audited_percent_score", "audited_pass_fail"]
    result = scored[fields].copy()
    result["raw_difference"] = pd.to_numeric(result.raw_score, errors="coerce") - result.audited_raw_score
    result["percent_difference"] = pd.to_numeric(result.percent_score, errors="coerce") - result.audited_percent_score
    result["raw_mismatch"] = result.raw_difference.ne(0)
    result["percent_mismatch"] = [not percent_matches(p, r) for p, r in
                                  zip(result.percent_score, result.audited_raw_score)]
    result["pass_fail_mismatch"] = result.pass_fail.ne(result.audited_pass_fail)
    result["original_rule_mismatch"] = result.pass_fail.ne(pd.to_numeric(result.raw_score, errors="coerce").map(pass_fail))
    result["diagnosis"] = ""
    for i, row in result.iterrows():
        reasons = []
        if row.raw_mismatch:
            reasons.append("Original raw score differs from official-key response total")
        if row.percent_mismatch:
            reasons.append("Original percentage differs from reconstructed raw/60*100 beyond tolerance")
        if row.pass_fail_mismatch:
            reasons.append("Original decision differs from independently reconstructed >=42 rule")
        if row.original_rule_mismatch:
            reasons.append("Original decision also contradicts original raw score >=42")
        result.at[i, "diagnosis"] = "; ".join(reasons)
    return result


def scoring_findings(comparison, findings):
    for _, row in comparison.iterrows():
        for flag, original, audited in [
            ("raw_mismatch", "raw_score", "audited_raw_score"),
            ("percent_mismatch", "percent_score", "audited_percent_score"),
            ("pass_fail_mismatch", "pass_fail", "audited_pass_fail"),
            ("original_rule_mismatch", "pass_fail", "audited_pass_fail"),
        ]:
            if row[flag]:
                findings.add("CRITICAL", "scoring", row.administration_id, original,
                             row[original], f"Independent value: {row[audited]}", row.diagnosis,
                             "Investigate upstream scoring; preserve both values. Do not silently correct inputs.")


def verify_long_units(admin, meta, long, denominator=60):
    expected = sum(len(form_items(meta, f)) for f in admin.form)
    counts = long.groupby("administration_id").size()
    return {
        "expected_rows": int(expected), "actual_rows": int(len(long)),
        "unique_administration_item_pairs": int(len(long.drop_duplicates(["administration_id", "item_id"]))),
        "all_administrations_present": set(admin.administration_id) == set(long.administration_id),
        "exactly_60_rows_per_administration": bool(counts.eq(denominator).all()),
        "mapping_matches_metadata": all(set(g.item_id) == set(form_items(meta, f)) for f, g in long.groupby("form")),
    }
