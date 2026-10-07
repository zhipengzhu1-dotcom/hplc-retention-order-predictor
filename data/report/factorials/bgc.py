"""
Extract the BGC (Harrieder/Witting) factorial from a pinned RepoRT checkout.

Usage:
    python -m data.report.factorials.bgc /path/to/RepoRT/processed_data

What this actually measures (verified against the pinned checkout, not
assumed from issue #4's prose -- see findings.md "BGC factorial: corrected"
for the discrepancy this surfaced):

  - The BGC study (datasets 0236-0259, 0310-0357; source == "Dataset - BGC" or
    "BGC", authors Harrieder/Witting) covers 6 columns x T{30,40,50}, but only
    3 of those 6 columns (Restek Raptor Biphenyl, Waters CORTECS UPLC C18,
    Waters ACQUITY UPLC HSS C18 -- names already swap-corrected) were run in
    BOTH modifiers. The other 3 (Waters ACQUITY UPLC HSS T3, Phenomenex
    Kinetex XB-C18, Waters ACQUITY UPLC BEH C18) are ACN-only. There is no
    MeOH run for those 3 columns anywhere in the repository under this
    study -- checked directly, not inferred.
  - So there are 27 condition cells, not 36: 18 ACN cells (6 col x 3T) + 9
    MeOH cells (3 col x 3T), each cell backed by a pos-mode and neg-mode MSMS
    dataset (CCS pos/neg exist at T=40 only and are excluded here -- distinct
    measurement technique, not a retention-time replicate).
  - Per-cell compound count (union of pos+neg canonical-success inchikeys):
    374-643 compounds per cell.
  - Shared core across the 18 ACN cells (6 columns x 3T): 203 compounds.
  - Shared core across the 18 modifier-crossed cells (3 columns x 3T x 2
    modifiers): 229 compounds.
  - Shared core across all 27 cells: 180 compounds.
  - issue #4's "~500 compound shared core / 493-499 per column" describes a
    narrower per-column quantity: the intersection of a single column's OWN
    cells (e.g. BEH C18 across its 3 ACN-only temperatures: 504 compounds),
    not a cross-column, cross-condition intersection. Both are computed here
    so you can pick the right denominator for your use.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict

from data.report.loader import iter_datasets, load_rtdata

BGC_CROSSED_COLUMNS = {
    "Restek Raptor Biphenyl",
    "Waters CORTECS UPLC C18",
    "Waters ACQUITY UPLC HSS C18",
}
BGC_ACN_ONLY_COLUMNS = {
    "Waters ACQUITY UPLC HSS T3",
    "Phenomenex Kinetex XB-C18",
    "Waters ACQUITY UPLC BEH C18",
}
BGC_ALL_COLUMNS = BGC_CROSSED_COLUMNS | BGC_ACN_ONLY_COLUMNS


def _is_bgc(d) -> bool:
    return "BGC" in d.source and ("Harrieder" in d.authors or "Witting" in d.authors)


def bgc_cells(processed_data_root: str) -> dict[tuple, list[str]]:
    """(column_name, temperature, modifier) -> [dataset ids], MSMS only (no CCS)."""
    cells: dict[tuple, list[str]] = defaultdict(list)
    for d in iter_datasets(processed_data_root):
        if not _is_bgc(d) or not d.raw_name.startswith("MSMS_"):
            continue
        cells[(d.column_name, d.column_temperature, d.modifier)].append(d.id)
    return dict(cells)


def cell_compound_sets(processed_data_root: str) -> dict[tuple, set[str]]:
    cells = bgc_cells(processed_data_root)
    out = {}
    for key, ids in cells.items():
        s: set[str] = set()
        for did in ids:
            rows = load_rtdata(processed_data_root, did, "canonical_success")
            s |= {r["inchikey.std"] for r in rows}
        out[key] = s
    return out


def shared_cores(processed_data_root: str) -> dict:
    sets = cell_compound_sets(processed_data_root)
    acn_cells = [k for k in sets if k[2] == "ACN"]
    crossed_cells = [k for k in sets if k[0] in BGC_CROSSED_COLUMNS]
    all_cells = list(sets)
    per_column_own_intersection = {}
    by_col = defaultdict(list)
    for k, s in sets.items():
        by_col[k[0]].append(s)
    for col, sl in by_col.items():
        per_column_own_intersection[col] = len(set.intersection(*sl))
    return {
        "n_cells": len(sets),
        "cell_sizes": {str(k): len(v) for k, v in sets.items()},
        "acn_6col_3T_shared_core": len(set.intersection(*[sets[k] for k in acn_cells])),
        "modifier_crossed_3col_3T_2mod_shared_core": len(
            set.intersection(*[sets[k] for k in crossed_cells])
        ),
        "all_27_cells_shared_core": len(set.intersection(*[sets[k] for k in all_cells])),
        "per_column_own_cells_intersection": per_column_own_intersection,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_data_root")
    args = ap.parse_args()
    result = shared_cores(args.processed_data_root)
    print(json.dumps(result, indent=2))
