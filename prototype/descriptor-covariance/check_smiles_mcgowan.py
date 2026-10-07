#!/usr/bin/env python3
"""Validate wsu94_smiles.csv against Table S-2's McGowan V (computable from structure).

McGowan characteristic volume (cm^3 mol^-1 / 100):
V = (sum of atom contributions - 6.56 * n_bonds) / 100, bonds counted once regardless
of order; n_bonds = n_atoms - 1 + n_rings (with H included).
Any SMILES whose computed V differs from Table S-2 by > 0.005 is flagged.
Run with the soluteml venv python (needs rdkit).
"""

import csv
import os
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "sources", "wsu-lser")

MCGOWAN = {"C": 16.35, "H": 8.71, "O": 12.43, "N": 14.39, "F": 10.48,
           "Cl": 20.95, "Br": 26.21, "I": 34.53, "S": 22.91, "P": 24.87}


def mcgowan_v(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    n_atoms = len(atoms)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_bonds = n_atoms - 1 + n_rings
    return (sum(MCGOWAN[a] for a in atoms) - 6.56 * n_bonds) / 100.0


def norm(name):
    # Names are clean since the 2026-08-16 extract fix; kept as an identity
    # hook so callers still have one place to normalise if S-2 is re-extracted.
    return name


def main():
    table_v = {}
    with open(os.path.join(SRC, "wsu2019-descriptors-s2.csv")) as f:
        for rec in csv.DictReader(f):
            table_v[norm(rec["compound"])] = float(rec["V"])

    bad = 0
    with open(os.path.join(HERE, "wsu94_smiles.csv")) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 94, len(rows)
    for rec in rows:
        v_calc = mcgowan_v(rec["smiles"])
        v_tab = table_v[rec["compound"]]
        if abs(v_calc - v_tab) > 0.005:
            bad += 1
            print(f"MISMATCH {rec['compound']}: calc {v_calc:.4f} vs S-2 {v_tab:.4f}")
    print(f"{len(rows)} SMILES checked, {bad} McGowan-V mismatches")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
