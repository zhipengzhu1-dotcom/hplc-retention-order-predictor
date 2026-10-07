#!/usr/bin/env python3
"""Which of the 94 WSU compounds are inside SoluteML's training corpus?

#42 inferred the WSU compounds are "almost certainly inside SoluteDB" from
provenance. SoluteDB is CC-BY and inspectable (Zenodo 5792296), so this settles
it row by row instead.

Matching is deliberately three-tiered, loosest last, and each compound is
reported with the tier that matched it:

  exact      -- SoluteDB's own key: non-standard fixed-H InChI, string-equal
  skeleton   -- standard InChIKey connectivity block (first 14 chars), which
                matches across tautomer/protomer and stereo differences
  none       -- absent under both

The skeleton tier is what decides contamination: a model trained on the
tautomer of a compound has seen that compound for our purposes.

Writes overlap.csv: compound, smiles, match_tier, db_E, db_S, db_A, db_B, db_L.

Usage: python overlap_solutedb.py <solutedb_selected.csv> <out overlap.csv>
Needs rdkit -- run with env/ionisation-venv/bin/python.
"""

import csv
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
SMILES_CSV = os.path.join(HERE, "..", "descriptor-covariance", "wsu94_smiles.csv")


def fixed_h_inchi(mol) -> str:
    """SoluteDB's key format: non-standard InChI with a fixed-H layer."""
    return MolToInchi(mol, options="/FixedH") or ""


def skeleton(mol) -> str:
    key = MolToInchiKey(mol) or ""
    return key.split("-")[0]


def main() -> None:
    db_rows = list(csv.DictReader(open(sys.argv[1])))
    by_inchi: dict[str, dict] = {}
    by_skeleton: dict[str, dict] = {}
    n_bad = 0
    for r in db_rows:
        by_inchi.setdefault(r["InChI"].strip(), r)
        mol = Chem.MolFromSmiles(r["SMILES"])
        if mol is None:
            n_bad += 1
            continue
        sk = skeleton(mol)
        if sk:
            by_skeleton.setdefault(sk, r)
    print(f"SoluteDB: {len(db_rows)} rows, {len(by_skeleton)} skeletons "
          f"({n_bad} unparseable SMILES)")

    out = csv.writer(open(sys.argv[2], "w", newline=""))
    out.writerow(["compound", "smiles", "match_tier",
                  "db_E", "db_S", "db_A", "db_B", "db_L"])
    tally = {"exact": 0, "skeleton": 0, "none": 0}
    for r in csv.DictReader(open(SMILES_CSV)):
        mol = Chem.MolFromSmiles(r["smiles"])
        assert mol is not None, r["compound"]
        hit, tier = None, "none"
        key = fixed_h_inchi(mol).strip()
        if key and key in by_inchi:
            hit, tier = by_inchi[key], "exact"
        else:
            sk = skeleton(mol)
            if sk in by_skeleton:
                hit, tier = by_skeleton[sk], "skeleton"
        tally[tier] += 1
        out.writerow([r["compound"], r["smiles"], tier] +
                     [(hit or {}).get(k, "") for k in ("E", "S", "A", "B", "L")])
    print("match tiers:", tally)


if __name__ == "__main__":
    main()
