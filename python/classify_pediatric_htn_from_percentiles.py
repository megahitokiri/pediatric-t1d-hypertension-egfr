#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import math
from collections import defaultdict
from pathlib import Path


def load_percentile_module(module_path):
    spec = importlib.util.spec_from_file_location("bp_percentiles_module", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    data = module.Data()
    return data.percentile


def parse_float(value):
    if value is None:
        return None
    value = str(value).strip()
    if value == "" or value.upper() == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def sex_code(sex_value):
    if sex_value is None:
        return None
    sex_value = str(sex_value).strip().upper()
    if sex_value.startswith("F"):
        return 2
    if sex_value.startswith("M"):
        return 1
    return None


def fixed_threshold_category(sbp, dbp):
    if sbp is None or dbp is None:
        return None
    if sbp >= 140 or dbp >= 90:
        return "Stage2_HTN"
    if sbp >= 130 or dbp >= 80:
        return "Stage1_HTN"
    if 120 <= sbp < 130 and dbp < 80:
        return "Elevated"
    if sbp < 120 and dbp < 80:
        return "Normal"
    return None


def percentile_thresholds(percentile_data, sex, age_years, height_cm):
    if sex not in percentile_data:
        return None
    if age_years not in percentile_data[sex]:
        return None

    age_table = percentile_data[sex][age_years]
    heights = age_table["height"]
    idx = None
    for i, height in enumerate(heights):
        if height_cm <= height:
            idx = i
            break
    if idx is None:
        idx = len(heights) - 1

    return {
        "sys90": age_table["Systolic 90th"][idx],
        "sys95": age_table["Systolic 95th"][idx],
        "sys95p12": age_table["Systolic 95th + 12mmHg"][idx],
        "dia90": age_table["Diastolic 90th"][idx],
        "dia95": age_table["Diastolic 95th"][idx],
        "dia95p12": age_table["Diastolic 95th + 12mmHg"][idx],
    }


def pediatric_percentile_category(sbp, dbp, thresholds):
    if sbp is None or dbp is None or thresholds is None:
        return None

    elevated_sys = min(thresholds["sys90"], 120.0)
    elevated_dia = min(thresholds["dia90"], 80.0)
    stage1_sys = min(thresholds["sys95"], 130.0)
    stage1_dia = min(thresholds["dia95"], 80.0)
    stage2_sys = min(thresholds["sys95p12"], 140.0)
    stage2_dia = min(thresholds["dia95p12"], 90.0)

    if sbp >= stage2_sys or dbp >= stage2_dia:
      return "Stage2_HTN"
    if sbp >= stage1_sys or dbp >= stage1_dia:
      return "Stage1_HTN"
    if sbp >= elevated_sys or dbp >= elevated_dia:
      return "Elevated"
    return "Normal"


def classify_row(row, percentile_data, max_gap_days, pediatric_age_cutoff):
    sbp = parse_float(row.get("SBP"))
    dbp = parse_float(row.get("DBP"))
    age_at_bp = parse_float(row.get("age_at_bp"))
    height_cm = parse_float(row.get("HEIGHT_CM"))
    bmi = parse_float(row.get("BMI"))
    days_to_anthro = parse_float(row.get("days_to_anthro"))
    sex = sex_code(row.get("SEX_CD"))

    out = dict(row)
    out["bp_category"] = ""
    out["classification_method"] = ""
    out["age_year_floor"] = ""
    out["height_used_cm"] = ""
    out["bmi_used"] = ""
    out["valid_for_percentile_classification"] = "0"

    if sbp is None or dbp is None or age_at_bp is None or sex is None:
        return out

    if age_at_bp >= 18:
        return out

    age_floor = int(math.floor(age_at_bp))
    out["age_year_floor"] = str(age_floor)

    if age_at_bp < pediatric_age_cutoff:
        if height_cm is None:
            return out
        if days_to_anthro is not None and days_to_anthro > max_gap_days:
            return out
        thresholds = percentile_thresholds(percentile_data, sex, age_floor, height_cm)
        if thresholds is None:
            return out
        category = pediatric_percentile_category(sbp, dbp, thresholds)
        out["bp_category"] = category or ""
        out["classification_method"] = "percentile_age_sex_height"
        out["height_used_cm"] = "" if height_cm is None else f"{height_cm:.3f}"
        out["bmi_used"] = "" if bmi is None else f"{bmi:.3f}"
        out["valid_for_percentile_classification"] = "1" if category else "0"
        return out

    category = fixed_threshold_category(sbp, dbp)
    out["bp_category"] = category or ""
    out["classification_method"] = "fixed_threshold_age_13plus"
    out["height_used_cm"] = "" if height_cm is None else f"{height_cm:.3f}"
    out["bmi_used"] = "" if bmi is None else f"{bmi:.3f}"
    out["valid_for_percentile_classification"] = "1" if category else "0"
    return out


def summarise_patient(rows):
    valid_rows = [r for r in rows if r.get("bp_category")]
    categories = [r["bp_category"] for r in valid_rows]
    bmi_values = [parse_float(r.get("bmi_used")) for r in valid_rows if parse_float(r.get("bmi_used")) is not None]
    age_values = [parse_float(r.get("age_at_bp")) for r in valid_rows if parse_float(r.get("age_at_bp")) is not None]
    sbp_values = [parse_float(r.get("SBP")) for r in valid_rows if parse_float(r.get("SBP")) is not None]
    dbp_values = [parse_float(r.get("DBP")) for r in valid_rows if parse_float(r.get("DBP")) is not None]

    def median(values):
        values = sorted(values)
        if not values:
            return None
        n = len(values)
        mid = n // 2
        if n % 2 == 1:
            return values[mid]
        return (values[mid - 1] + values[mid]) / 2.0

    n_valid = len(valid_rows)
    n_elevated_plus = sum(cat in {"Elevated", "Stage1_HTN", "Stage2_HTN"} for cat in categories)
    n_stage1_plus = sum(cat in {"Stage1_HTN", "Stage2_HTN"} for cat in categories)
    n_stage2 = sum(cat == "Stage2_HTN" for cat in categories)

    return {
        "Patient_ID": rows[0]["Patient_ID"],
        "n_valid_peds_bp_days": n_valid,
        "n_elevated_plus_days": n_elevated_plus,
        "n_stage1plus_days": n_stage1_plus,
        "n_stage2_days": n_stage2,
        "ever_elevated_plus": int(n_elevated_plus >= 1),
        "ever_stage1plus": int(n_stage1_plus >= 1),
        "ever_stage2": int(n_stage2 >= 1),
        "peds_elevated_3plus_days": int(n_elevated_plus >= 3),
        "peds_htn_3plus_days": int(n_stage1_plus >= 3),
        "peds_stage2_3plus_days": int(n_stage2 >= 3),
        "peds_htn_burden": "" if n_valid == 0 else f"{(n_stage1_plus / n_valid):.6f}",
        "median_peds_bmi": "" if not bmi_values else f"{median(bmi_values):.6f}",
        "median_peds_age": "" if not age_values else f"{median(age_values):.6f}",
        "mean_peds_sbp": "" if not sbp_values else f"{(sum(sbp_values) / len(sbp_values)):.6f}",
        "mean_peds_dbp": "" if not dbp_values else f"{(sum(dbp_values) / len(dbp_values)):.6f}",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--bp-module", required=True)
    parser.add_argument("--visit-output", required=True)
    parser.add_argument("--patient-output", required=True)
    parser.add_argument("--max-gap-days", type=float, default=365.0)
    parser.add_argument("--pediatric-age-cutoff", type=float, default=13.0)
    args = parser.parse_args()

    percentile_data = load_percentile_module(args.bp_module)

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f)
        rows = [row for row in reader]

    visit_rows = [
        classify_row(row, percentile_data, args.max_gap_days, args.pediatric_age_cutoff)
        for row in rows
    ]

    Path(args.visit_output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.visit_output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(visit_rows[0].keys()))
        writer.writeheader()
        writer.writerows(visit_rows)

    grouped = defaultdict(list)
    for row in visit_rows:
        grouped[row["Patient_ID"]].append(row)

    patient_rows = [summarise_patient(grouped[key]) for key in sorted(grouped.keys())]
    with open(args.patient_output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(patient_rows[0].keys()))
        writer.writeheader()
        writer.writerows(patient_rows)


if __name__ == "__main__":
    main()
