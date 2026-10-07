"""
Extract clean ACN-vs-MeOH dataset pairs from a pinned RepoRT checkout.

"Clean pair" here means: two datasets on the SAME column (swap-corrected
column.name), SAME temperature (both populated and equal), SAME pH (both
eluent.A.pH and eluent.B.pH populated and equal), one run in ACN and one in
MeOH -- matched 1:1 (no dataset reused), greedily by maximum compound
overlap, so the result is a deduplicated set of pairs rather than every
combination of every same-condition dataset (repo-wide.

Usage:
    python -m data.report.factorials.acn_meoh_pairs /path/to/RepoRT/processed_data

HONESTY NOTE (see findings.md "ACN/MeOH pairs: could not reproduce 156
exactly"): issue #4 reported "156 same-column, same-T, same-pH ACN-vs-MeOH
pairs, best overlaps ~400 compounds". Running this exact, reproducible
definition against the pinned checkout yields 43 deduplicated 1:1 pairs (best
overlap 406 compounds -- matches the "~400" figure closely). A looser,
non-deduplicated many-to-many cross product over the same matching criterion
yields 110 pairs. Neither reproduces 156 exactly. issue #4's comment does not
document its exact matching algorithm, and it could not be reverse-engineered
with confidence from the prose alone -- so this module documents its own
criterion precisely instead of chasing an unverifiable prior number. The
qualitative finding (many clean pairs exist, best overlaps ~400 compounds) is
solidly corroborated either way.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict

from data.report.loader import iter_datasets, load_rtdata


def _populated(v: str) -> bool:
    return v is not None and str(v).strip() != ""


def clean_pairs(processed_data_root: str) -> list[dict]:
    datasets = list(iter_datasets(processed_data_root))
    acn = [d for d in datasets if d.modifier == "ACN"]
    meoh = [d for d in datasets if d.modifier == "MeOH"]

    def fp(d):
        return (d.column_name, d.column_temperature, d.eluent_a_ph, d.eluent_b_ph)

    acn_by, meoh_by = defaultdict(list), defaultdict(list)
    for d in acn:
        k = fp(d)
        if all(_populated(x) for x in k):
            acn_by[k].append(d)
    for d in meoh:
        k = fp(d)
        if all(_populated(x) for x in k):
            meoh_by[k].append(d)

    common = sorted(set(acn_by) & set(meoh_by))

    cache: dict[str, set] = {}

    def compset(did: str) -> set:
        if did not in cache:
            rows = load_rtdata(processed_data_root, did, "canonical_success")
            cache[did] = {r["inchikey.std"] for r in rows}
        return cache[did]

    pairs = []
    for k in common:
        combos = []
        for a in acn_by[k]:
            for m in meoh_by[k]:
                combos.append((len(compset(a.id) & compset(m.id)), a, m))
        combos.sort(key=lambda x: -x[0])
        used_a, used_m = set(), set()
        for ov, a, m in combos:
            if a.id in used_a or m.id in used_m:
                continue
            used_a.add(a.id)
            used_m.add(m.id)
            pairs.append(
                {
                    "acn_id": a.id,
                    "meoh_id": m.id,
                    "column": a.column_name,
                    "temperature": a.column_temperature,
                    "ph_a": a.eluent_a_ph,
                    "ph_b": a.eluent_b_ph,
                    "n_acn_compounds": len(compset(a.id)),
                    "n_meoh_compounds": len(compset(m.id)),
                    "overlap": ov,
                }
            )
    return pairs


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_data_root")
    args = ap.parse_args()
    result = clean_pairs(args.processed_data_root)
    print(json.dumps({"n_pairs": len(result), "pairs": result}, indent=2))
