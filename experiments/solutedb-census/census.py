#!/usr/bin/env python3
"""How much of real chromatographic chemistry is inside SoluteML's training corpus?

Motivation, corrected. This started as a test of whether RepoRT could supply a
CLEAN (out-of-corpus) compound set for #42's stratification questions. It cannot,
and the reason was already in the repo: `experiments/report-overlap` shows RepoRT
has zero isocratic datasets, no dwell-volume field in its schema, and a `column.t0`
imputed from geometry -- so no honest `k` exists in any dataset -- and RepoRT
carries no measured Abraham descriptors either. A clean compound set with no
label is not an evaluation set.

What survives is a question that needs no labels at all, and that the
applicability-domain object (#43/#40) needs answered: for the chemistry that
actually appears in reversed-phase practice, what fraction has SoluteML already
seen, and what fraction is it extrapolating to?

Four tiers, coarsest first:

  in_corpus        InChIKey skeleton present in SoluteDB. The model has seen it.
  scaffold_seen    Not in SoluteDB, but its Murcko scaffold is. Chung 2022's
                   substructure split is exactly this distinction, and reports
                   roughly double the error across it.
  novel_in_envelope Neither, and inside SoluteML's declared training envelope
                   (neutral, elements H C N O S P F Cl Br I). Genuine
                   extrapolation, and the only tier a clean study could use.
  out_of_envelope  Outside the declared envelope: charged, or containing an
                   element the model was never trained on. #33/#39 territory --
                   refused by structure, not scored with low confidence.

Usage:
  env/ionisation-venv/bin/python experiments/solutedb-census/census.py \
      <RepoRT/processed_data> <solutedb_selected.csv>

SoluteDB: Zenodo 5792296, CC BY 4.0, `SoluteDB_selected_data.xlsx` extracted by
prototype/soluteml-bias/read_solutedb.py.
RepoRT: not vendored; clone and pin per data/report/README.md (9de8d60).
"""

import csv
import os
import sys
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)

from data.report.loader import iter_datasets, load_rtdata  # noqa: E402
from data.report.splits.scaffold_split import (  # noqa: E402
    SMRT_DATASET_IDS, build_compound_universe,
)

# SoluteML's declared training domain (Chung 2022): neutral solutes, these
# elements only. Membership here is a STRUCTURAL test, not a confidence
# threshold -- consistent with the in-envelope / degraded / refuse tiers.
ENVELOPE = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}


def skeleton_from_key(inchikey: str) -> str:
    return inchikey.split("-")[0]


def scaffold(mol) -> str:
    try:
        s = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(s) if s is not None and s.GetNumAtoms() else ""
    except Exception:
        return ""


def in_envelope(mol) -> bool:
    if Chem.GetFormalCharge(mol) != 0:
        return False
    return all(a.GetSymbol() in ENVELOPE for a in mol.GetAtoms())


def build_solutedb_index(path: str):
    skeletons, scaffolds = set(), set()
    n = 0
    for r in csv.DictReader(open(path)):
        mol = Chem.MolFromSmiles(r["SMILES"])
        if mol is None:
            continue
        n += 1
        key = inchi.MolToInchiKey(mol)
        if key:
            skeletons.add(skeleton_from_key(key))
        sc = scaffold(mol)
        if sc:
            scaffolds.add(sc)
    print(f"SoluteDB: {n} parsed, {len(skeletons)} skeletons, "
          f"{len(scaffolds)} Murcko scaffolds")
    return skeletons, scaffolds


def main() -> None:
    root, db_csv = sys.argv[1], sys.argv[2]
    db_skeletons, db_scaffolds = build_solutedb_index(db_csv)

    universe = build_compound_universe(root)
    print(f"RepoRT (non-SMRT): {len(universe)} unique compounds by inchikey.std")

    tier_of = {}
    tally = Counter()
    reasons = Counter()
    for ik, smi in universe.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            tier_of[ik] = "unparseable"
            tally["unparseable"] += 1
            continue
        if skeleton_from_key(ik) in db_skeletons:
            tier = "in_corpus"
        elif scaffold(mol) and scaffold(mol) in db_scaffolds:
            tier = "scaffold_seen"
        elif in_envelope(mol):
            tier = "novel_in_envelope"
        else:
            tier = "out_of_envelope"
            if Chem.GetFormalCharge(mol) != 0:
                reasons["charged"] += 1
            bad = {a.GetSymbol() for a in mol.GetAtoms()} - ENVELOPE
            if bad:
                reasons["element:" + ",".join(sorted(bad))] += 1
        tier_of[ik] = tier
        tally[tier] += 1

    total = sum(tally.values())
    print(f"\n{'tier':>18} {'compounds':>10} {'share':>8}")
    for t in ("in_corpus", "scaffold_seen", "novel_in_envelope",
              "out_of_envelope", "unparseable"):
        if tally[t]:
            print(f"{t:>18} {tally[t]:10d} {tally[t] / total:8.1%}")

    print(f"\n  out_of_envelope reasons (a compound may have both):")
    for r, c in reasons.most_common(8):
        print(f"    {r:>28}: {c}")

    # Weight by how often each compound is actually measured: a corpus-coverage
    # figure over unique compounds and one over measurements are different
    # claims, and the second is the one that describes practice.
    occurrences = Counter()
    for d in iter_datasets(root):
        if d.id in SMRT_DATASET_IDS:
            continue
        for r in load_rtdata(root, d.id, "canonical_success"):
            ik = (r.get("inchikey.std") or "").strip()
            if ik in tier_of:
                occurrences[tier_of[ik]] += 1
    tot_occ = sum(occurrences.values())
    if tot_occ:
        print(f"\n  weighted by measurements ({tot_occ} retention records):")
        for t, c in occurrences.most_common():
            print(f"    {t:>18} {c:10d} {c / tot_occ:8.1%}")

    with open(os.path.join(HERE, "census.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["inchikey_std", "smiles", "tier"])
        for ik, smi in universe.items():
            w.writerow([ik, smi, tier_of[ik]])
    print(f"\nwrote census.csv ({len(universe)} rows)")


if __name__ == "__main__":
    main()
