"""Option frequencies, tie-preserving rest-score bands and exploratory concerns."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from .ctt import form_matrices


@dataclass(frozen=True)
class DistractorConfig:
    rare_below: float = 0.05
    performance_groups: int = 5
    minimum_band_n: int = 30
    minimum_option_n: int = 5
    mean_rest_margin: float = 1.0
    rising_selection_margin: float = 0.05
    lower_band_error_opportunity: float = 0.20

    def __post_init__(self):
        if not 0 < self.rare_below < 1 or not 0 <= self.rising_selection_margin <= 1:
            raise ValueError("Invalid option frequency threshold")
        if self.performance_groups < 1 or self.minimum_band_n < 1 or self.minimum_option_n < 1:
            raise ValueError("Group counts/minimum sizes must be positive")
        if self.mean_rest_margin < 0 or not 0 <= self.lower_band_error_opportunity <= 1:
            raise ValueError("Invalid rest-score or error-opportunity threshold")

    def as_dict(self):
        return asdict(self)


def performance_groups(rest, requested=5):
    """Nearest cumulative-quantile boundaries on distinct scores; ties never split.

    Constrain each cut to leave at least one distinct score for each remaining
    group. This yields exactly min(requested, distinct scores) ordered nonempty
    groups, with unequal candidate counts when ties prevent balanced quintiles.
    """
    rest = np.asarray(rest, dtype=float)
    if rest.ndim != 1 or not len(rest) or not np.isfinite(rest).all() or requested < 1:
        raise ValueError("Expected finite rest scores and a positive group count")
    scores, counts = np.unique(rest, return_counts=True)
    g = min(requested, len(scores))
    cumulative = counts.cumsum()
    last, cuts = -1, []
    for j in range(1, g):
        allowed = np.arange(last+1, len(scores)-(g-j))
        cut = int(allowed[np.argmin(np.abs(cumulative[allowed]-len(rest)*j/g))])
        cuts.append(scores[cut])
        last = cut
    labels = np.searchsorted(cuts, rest, side="left")+1
    return labels.astype(int), {"n_groups": g, "n_distinct_rest_scores": len(scores), "boundaries": cuts,
                                "method": "nearest cumulative quantile on distinct scores; constrained nonempty; ties intact"}


def option_concerns(row, config, very_easy_above=0.90):
    incorrect = not row["is_correct"]
    rare = row["proportion_administered"] < config.rare_below
    adequate_choices = row["n_selected"] >= config.minimum_option_n and row["n_correct_option_selected"] >= config.minimum_option_n
    mean_concern = bool(incorrect and adequate_choices and np.isfinite(row["option_minus_correct_mean_rest"])
                        and row["option_minus_correct_mean_rest"] >= config.mean_rest_margin)
    adequate_bands = row["n_performance_groups"] >= 2 and min(row["lowest_group_n"], row["highest_group_n"]) >= config.minimum_band_n
    rising = bool(incorrect and adequate_bands and row["n_selected"] >= config.minimum_option_n
                  and row["highest_minus_lowest_selection"] >= config.rising_selection_margin)
    low_attraction = bool(incorrect and rare and adequate_bands
                          and row["item_proportion_correct"] <= very_easy_above
                          and row["lowest_group_error_rate"] >= config.lower_band_error_opportunity
                          and row["lowest_group_selection"] < config.rare_below)
    reasons = []
    if mean_concern:
        reasons.append("incorrect-option selectors have mean rest >= correct selectors plus margin")
    if rising:
        reasons.append("incorrect-option selection higher in highest than lowest rest-score band")
    if low_attraction:
        reasons.append(f"rare distractor also attracts <{config.rare_below:.0%} in a sufficiently sized lower band with >={config.lower_band_error_opportunity:.0%} errors; item not very easy")
    return {"rare_option": bool(rare), "rare_distractor": bool(incorrect and rare),
            "mean_rest_pattern_concern": mean_concern, "rising_selection_concern": rising,
            "distractor_pattern_concern": mean_concern or rising,
            "lower_band_non_attraction": low_attraction,
            "potentially_non_functioning_distractor": mean_concern or rising or low_attraction,
            "concern_reason": "; ".join(reasons)}


def analyze_distractors(admin, long, metadata, form, config=None, very_easy_above=0.90):
    config = config or DistractorConfig()
    sittings, members, x, responses = form_matrices(admin, long, metadata, form)
    meta = metadata.set_index("item_id")
    n = len(sittings)
    totals = x.sum(axis=1)
    options, band_rows, group_definitions = [], [], []
    for j, item in enumerate(members):
        rest = totals-x[:, j]
        answer = meta.at[item, "correct_option"]
        selected_correct = responses[:, j] == answer
        correct_mean = rest[selected_correct].mean() if selected_correct.any() else np.nan
        omitted = responses[:, j] == ""
        valid_n = int((~omitted).sum())
        groups, definition = performance_groups(rest, config.performance_groups)
        bottom, top = groups.min(), groups.max()
        bottom_mask, top_mask = groups == bottom, groups == top
        error_rate_bottom = float(1-x[bottom_mask, j].mean())
        for group in range(1, definition["n_groups"]+1):
            mask = groups == group
            group_definitions.append({"form": form, "item_id": item, "performance_group": group,
                                      "n_candidates": int(mask.sum()), "rest_score_min": int(rest[mask].min()),
                                      "rest_score_max": int(rest[mask].max()), "mean_rest_score": float(rest[mask].mean()),
                                      "n_distinct_rest_scores": definition["n_distinct_rest_scores"],
                                      "n_performance_groups": definition["n_groups"], "grouping_method": definition["method"]})
        for option in ["A", "B", "C", "D"]:
            selected = responses[:, j] == option
            mean_rest = rest[selected].mean() if selected.any() else np.nan
            row = {"form": form, "item_id": item, "domain": meta.at[item, "domain"],
                   "anchor_item": meta.at[item, "anchor_item"], "option": option,
                   "correct_option": answer, "is_correct": option == answer,
                   "n_administered": n, "n_valid_marked": valid_n, "n_selected": int(selected.sum()),
                   "proportion_administered": float(selected.mean()),
                   "proportion_valid_marked": float(selected.sum()/valid_n) if valid_n else np.nan,
                   "mean_rest_score_selected": float(mean_rest), "mean_rest_correct_selected": float(correct_mean),
                   "option_minus_correct_mean_rest": float(mean_rest-correct_mean),
                   "n_correct_option_selected": int(selected_correct.sum()), "n_omitted": int(omitted.sum()),
                   "item_proportion_correct": float(x[:, j].mean()),
                   "n_performance_groups": definition["n_groups"], "lowest_group_n": int(bottom_mask.sum()),
                   "highest_group_n": int(top_mask.sum()), "lowest_group_error_rate": error_rate_bottom,
                   "lowest_group_selection": float(selected[bottom_mask].mean()),
                   "highest_group_selection": float(selected[top_mask].mean()),
                   "highest_minus_lowest_selection": float(selected[top_mask].mean()-selected[bottom_mask].mean())}
            row.update(option_concerns(row, config, very_easy_above))
            options.append(row)
            for group in range(1, definition["n_groups"]+1):
                mask = groups == group
                band_n, band_valid = int(mask.sum()), int((mask & ~omitted).sum())
                band_rows.append({"form": form, "item_id": item, "option": option,
                                  "is_correct": option == answer, "performance_group": group,
                                  "n_administered_group": band_n, "n_valid_marked_group": band_valid,
                                  "n_selected_group": int(selected[mask].sum()), "n_omitted_group": int(omitted[mask].sum()),
                                  "selection_proportion_administered": float(selected[mask].mean()),
                                  "selection_proportion_valid_marked": float(selected[mask].sum()/band_valid) if band_valid else np.nan,
                                  "mean_rest_score_group": float(rest[mask].mean()),
                                  "rest_score_min": int(rest[mask].min()), "rest_score_max": int(rest[mask].max())})
    return pd.DataFrame(options), pd.DataFrame(band_rows), pd.DataFrame(group_definitions)


def exploratory_flags(items, deletion, options, ctt_config):
    option_flags = options.groupby(["form", "item_id"]).agg(
        n_rare_distractors=("rare_distractor", "sum"),
        n_distractor_pattern_concerns=("distractor_pattern_concern", "sum"),
        n_potentially_non_functioning_distractors=("potentially_non_functioning_distractor", "sum"),
    ).reset_index()
    result = items.merge(deletion[["form", "item_id", "alpha_if_item_deleted", "alpha_change_if_deleted", "alpha_increases_if_deleted"]],
                         on=["form", "item_id"], validate="one_to_one").merge(option_flags, on=["form", "item_id"], validate="one_to_one")
    result["rare_distractor"] = result.n_rare_distractors.gt(0)
    result["distractor_pattern_concern"] = result.n_distractor_pattern_concerns.gt(0)
    result["potentially_non_functioning_distractor"] = result.n_potentially_non_functioning_distractors.gt(0)
    categories, reasons = [], []
    for row in result.itertuples():
        reason = []
        if row.negative_discrimination:
            reason.append("negative corrected item-total correlation: check key/content/administration")
        if row.undefined_discrimination:
            reason.append(row.discrimination_status)
        if row.low_discrimination:
            reason.append(f"corrected discrimination below {ctt_config.low_discrimination_below:g}")
        if row.attention_discrimination:
            reason.append(f"corrected discrimination in [{ctt_config.low_discrimination_below:g}, {ctt_config.attention_discrimination_below:g})")
        if row.very_easy:
            reason.append(f"proportion correct above {ctt_config.very_easy_above:g}; check intended content/mastery role")
        if row.very_difficult:
            reason.append(f"proportion correct below {ctt_config.very_difficult_below:g}; check intended content/response process")
        if row.distractor_pattern_concern:
            reason.append("incorrect-option rest-score pattern concern")
        if row.potentially_non_functioning_distractor:
            reason.append("potential distractor nonfunctioning supported by frequency/performance screen")
        if row.rare_distractor:
            reason.append("one or more rare distractors; rarity alone does not imply defect")
        if row.alpha_increases_if_deleted:
            reason.append(f"alpha change on deletion {row.alpha_change_if_deleted:+.6f}; not a deletion recommendation")
        if row.negative_discrimination:
            category = "Priority investigation"
        elif (row.low_discrimination or row.undefined_discrimination or row.distractor_pattern_concern
              or (row.alpha_change_if_deleted >= ctt_config.material_alpha_change
                  and row.corrected_item_total_r < ctt_config.attention_discrimination_below)):
            category = "Review"
        elif (row.very_easy or row.very_difficult or row.attention_discrimination
              or row.potentially_non_functioning_distractor or row.rare_distractor or row.alpha_increases_if_deleted):
            category = "Monitor"
        else:
            category = "No immediate statistical concern"
        categories.append(category)
        reasons.append("; ".join(reason) or "No prespecified exploratory screen triggered")
    result["investigative_category"] = categories
    result["investigation_reason"] = reasons
    result["_category_order"] = result.investigative_category.map({"Priority investigation": 0, "Review": 1, "Monitor": 2, "No immediate statistical concern": 3})
    result = result.sort_values(["_category_order", "corrected_item_total_r", "n_distractor_pattern_concerns", "alpha_change_if_deleted", "form", "item_id"],
                                ascending=[True, True, False, False, True, True], na_position="last", kind="stable")
    result["investigation_rank"] = np.arange(1, len(result)+1)
    return result.drop(columns="_category_order").reset_index(drop=True)
