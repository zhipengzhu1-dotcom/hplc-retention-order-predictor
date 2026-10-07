"""
Unseen-column split infrastructure, plus a direct quantification of how weak
it is -- issue #4 warned "308 of 412 coded columns are USP L1"; this module
measures the real number against the pinned checkout and names which columns
actually offer different chemistry.

Usage:
    python -m data.report.splits.unseen_column_split /path/to/RepoRT/processed_data

MEASURED (non-SMRT datasets, n=419; see findings.md for the full table):
  - 410/419 (97.9%) datasets carry a column.usp.code.
  - At the DATASET level: 306/410 coded datasets (74.6%) are USP L1.
    (issue #4 said "308 of 412" -- same order, small discrepancy explained in
    findings.md; not chased further.)
  - At the DISTINCT-COLUMN level: 38 of 58 named, coded columns (65.5%) are
    L1.
  - 9 datasets (1785 compounds) carry NO column name and NO USP code at all
    (PredRet-sourced: MTBLS17, MTBLS19, CS15, CS16, CS18, CS19, ABC,
    #HILIC_tip, AjsUoB) -- unusable for ANY column-axis split, flagged
    separately.

GENUINELY DIFFERENT CHEMISTRY, quantified: a held-out column is only a
meaningful "unseen chemistry" test if it isn't just another L1 reversed-phase
C18. Below is every non-L1 column with >=100 supporting compounds (the
threshold below which a held-out test set is too small to trust) -- this is
the full corrected list; issue #4 named only one example (Raptor Biphenyl).

    USP L11 (alkyl/aryl, "biphenyl-type")  Restek Raptor Biphenyl        770
    USP L11 ("polar-modified")             Phenomenex Synergi Polar-RP  748
    USP L68 (amide, HILIC-mode)            Waters ACQUITY UPLC BEH Amide 1198
    USP L122 (HILIC, zwitterionic)         Merck SeQuant ZIC-pHILIC     684
    USP L114 (HILIC, zwitterionic)         Merck SeQuant ZIC-HILIC      647
    USP L7  (C8)                           Waters ACQUITY UPLC BEH C8   197
    USP L43 (pentafluorophenyl, PFP)       Ascentis Express F5 (PFP)    194
    USP L10 (cyano)                        Ascentis Express ES-Cyano    185
    USP L11 (phenyl-hexyl)                 Ascentis Express Phenyl-Hexyl 184
    USP L68 (amide, HILIC-mode)            Waters XBridge BEH Amide     160
    USP L3  (HILIC, bare silica)           Thermo Accucore HILIC        152
    USP L3  (HILIC, bare silica)           Waters ACQUITY UPLC BEH HILIC 123

Of these, the 4 HILIC-mode phases (BEH Amide, ZIC-pHILIC, ZIC-HILIC, XBridge
BEH Amide, BEH HILIC -- 5 columns, retention MECHANISM is hydrophilic
interaction, not reversed-phase) are the only ones structurally guaranteed to
generalize differently than any RP phase: a compound's HILIC retention order
is not simply a monotonic transform of its RP retention order. The PFP,
cyano, phenyl-hexyl and C8 phases are still RP but with meaningfully
different selectivity (pi-pi stacking, dipole interactions). Biphenyl and
Polar-RP were issue #4's own example.

Given issue #4's independent, direct measurement that pairwise elution-order
reversal between two C18 (both-L1) columns is only ~1.0-3.9%, while
biphenyl-vs-C18 is ~4.6-6.2%, an "unseen L1-vs-L1" column split is measurably
close to a no-op: the model barely needs to generalize. Recommend restricting
"unseen-column" claims to holding out one of the 12 columns above (or the 5
HILIC-mode ones specifically, for a maximally different test), reported as
named, individual adversarial cases -- NOT averaged into one cross-validated
"unseen column" number, which would silently be dominated by near-trivial
L1-vs-L1 folds (38 of the pool's 58 distinct columns).
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict

from data.report.loader import iter_datasets, load_rtdata

SMRT_DATASET_IDS = {"0186", "0209"}

# Non-L1 columns with >= 100 supporting compounds, i.e. genuinely different
# chemistry AND enough data to be a usable held-out test set. Derived from
# the measurement in this module's docstring -- recompute with
# `column_support()` below rather than trusting this list blindly if the
# pinned commit changes.
GENUINELY_DIFFERENT_CHEMISTRY_COLUMNS = [
    "Restek Raptor Biphenyl",
    "Phenomenex Synergi Polar-RP",
    "Waters ACQUITY UPLC BEH Amide",
    "Merck SeQuant ZIC-pHILIC",
    "Merck SeQuant ZIC-HILIC",
    "Waters ACQUITY UPLC BEH C8",
    "Merck Supelco Ascentis Express F5 (PFP)",
    "Merck Supelco Ascentis Express ES-Cyano",
    "Merck Supelco Ascentis Express Phenyl-Hexyl",
    "Waters XBridge BEH Amide",
    "Thermo Scientific Accucore HILIC",
    "Waters ACQUITY UPLC BEH HILIC",
]

HILIC_MODE_COLUMNS = [
    "Waters ACQUITY UPLC BEH Amide",
    "Merck SeQuant ZIC-pHILIC",
    "Merck SeQuant ZIC-HILIC",
    "Waters XBridge BEH Amide",
    "Waters ACQUITY UPLC BEH HILIC",
]


def column_support(processed_data_root: str) -> dict:
    """Return per-column compound support and USP-code census."""
    datasets = [d for d in iter_datasets(processed_data_root) if d.id not in SMRT_DATASET_IDS]

    column_to_compounds: dict[str, set] = defaultdict(set)
    usp_by_column: dict[str, str] = {}
    unnamed_column_datasets = []
    for d in datasets:
        rows = load_rtdata(processed_data_root, d.id, "canonical_success")
        ikeys = {r["inchikey.std"] for r in rows}
        if d.column_name == "":
            unnamed_column_datasets.append(d.id)
            continue
        column_to_compounds[d.column_name] |= ikeys
        usp_by_column[d.column_name] = d.column_usp_code

    dataset_usp_counter = Counter(
        d.column_usp_code for d in datasets if d.column_usp_code
    )
    n_datasets_coded = sum(dataset_usp_counter.values())
    n_datasets_l1 = dataset_usp_counter.get("L1", 0)

    column_usp_counter = Counter(v for v in usp_by_column.values() if v)
    n_columns_coded = sum(column_usp_counter.values())
    n_columns_l1 = column_usp_counter.get("L1", 0)

    return {
        "n_datasets_total_non_smrt": len(datasets),
        "n_datasets_coded": n_datasets_coded,
        "n_datasets_l1": n_datasets_l1,
        "dataset_l1_fraction": n_datasets_l1 / n_datasets_coded if n_datasets_coded else None,
        "n_distinct_columns_coded": n_columns_coded,
        "n_distinct_columns_l1": n_columns_l1,
        "column_l1_fraction": n_columns_l1 / n_columns_coded if n_columns_coded else None,
        "unnamed_column_datasets": unnamed_column_datasets,
        "unnamed_column_compound_count": len(
            set().union(
                *[
                    {r["inchikey.std"] for r in load_rtdata(processed_data_root, did, "canonical_success")}
                    for did in unnamed_column_datasets
                ]
            )
        ) if unnamed_column_datasets else 0,
        "column_support": {k: len(v) for k, v in column_to_compounds.items()},
        "usp_by_column": usp_by_column,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_data_root")
    args = ap.parse_args()
    print(json.dumps(column_support(args.processed_data_root), indent=2))
