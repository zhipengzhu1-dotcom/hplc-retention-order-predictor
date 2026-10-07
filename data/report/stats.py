"""
Row counts, missing-field rates, and anomaly checks for a pinned RepoRT
checkout. All numbers in findings.md's "Measured stats" section come from
running this module.

Usage:
    python -m data.report.stats /path/to/RepoRT/processed_data
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict

from data.report.loader import iter_datasets, load_rtdata

SMRT_DATASET_IDS = {"0186", "0209"}


def _populated(v) -> bool:
    return v is not None and str(v).strip() != ""


def compute_stats(processed_data_root: str) -> dict:
    datasets = list(iter_datasets(processed_data_root))
    non_smrt = [d for d in datasets if d.id not in SMRT_DATASET_IDS]
    n = len(non_smrt)

    field_rates = {
        "column.temperature": sum(1 for d in non_smrt if _populated(d.column_temperature)) / n,
        "eluent.A.pH (incl 0-sentinel)": sum(1 for d in non_smrt if _populated(d.eluent_a_ph)) / n,
        "eluent.A.pH (nonzero only)": sum(
            1 for d in non_smrt if _populated(d.eluent_a_ph) and float(d.eluent_a_ph) != 0
        ) / n,
        "column.usp.code": sum(1 for d in non_smrt if _populated(d.column_usp_code)) / n,
        "column.length": sum(1 for d in non_smrt if _populated(d.column_length)) / n,
        "column.flowrate": sum(1 for d in non_smrt if _populated(d.column_flowrate)) / n,
        "column.name": sum(1 for d in non_smrt if _populated(d.column_name)) / n,
    }

    total_rows = 0
    smrt_rows = 0
    non_smrt_unique_compounds = set()
    dup_conflict_datasets = 0
    dup_conflict_rows = 0
    dup_exact_rows = 0
    for d in datasets:
        rows = load_rtdata(processed_data_root, d.id, "canonical_success")
        total_rows += len(rows)
        if d.id in SMRT_DATASET_IDS:
            smrt_rows += len(rows)
            continue
        non_smrt_unique_compounds |= {r["inchikey.std"] for r in rows}
        by_key: dict[str, list[str]] = defaultdict(list)
        for r in rows:
            by_key[r["inchikey.std"]].append(r["rt"])
        had_conflict = False
        for ik, rts in by_key.items():
            if len(rts) > 1:
                if len(set(rts)) == 1:
                    dup_exact_rows += len(rts) - 1
                else:
                    dup_conflict_rows += len(rts) - 1
                    had_conflict = True
        if had_conflict:
            dup_conflict_datasets += 1

    unnamed_column = [d.id for d in non_smrt if d.column_name == ""]

    return {
        "n_datasets_total": len(datasets),
        "n_datasets_non_smrt": n,
        "n_datasets_smrt": len(datasets) - n,
        "total_rtdata_rows_canonical_success": total_rows,
        "smrt_rows": smrt_rows,
        "non_smrt_rows": total_rows - smrt_rows,
        "non_smrt_unique_compounds_by_inchikey": len(non_smrt_unique_compounds),
        "field_populated_rates": field_rates,
        "name_swap_affected_dataset_ids": sorted(
            d.id for d in datasets if d.name_swapped
        ),
        "unnamed_column_dataset_ids": unnamed_column,
        "datasets_with_conflicting_duplicate_inchikey_rt": dup_conflict_datasets,
        "conflicting_duplicate_rows_total": dup_conflict_rows,
        "exact_duplicate_rows_total": dup_exact_rows,
        "anomaly_note_duplicates": (
            "inchikey.std strips stereochemistry. Within a single dataset, "
            f"{dup_conflict_datasets} datasets have >=1 inchikey.std shared by "
            f"rows with DIFFERENT rt values ({dup_conflict_rows} such extra "
            "rows total) -- almost certainly distinct stereoisomers/positional "
            "isomers collapsed onto one standardized key, each genuinely "
            "eluting differently. Any pipeline (including this one) that "
            "dedupes or joins on inchikey.std alone will silently conflate "
            "them. A further ~559 rows are EXACT duplicates (same key, same "
            "rt) and are pure redundancy, not a bug."
        ),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_data_root")
    args = ap.parse_args()
    print(json.dumps(compute_stats(args.processed_data_root), indent=2))
