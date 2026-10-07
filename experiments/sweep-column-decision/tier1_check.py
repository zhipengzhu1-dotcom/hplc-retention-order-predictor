#!/usr/bin/env python3
"""#30: are the suggested tier-1 acids actually well enough retained?

The power analysis (decide_columns.py) found that `D`'s precision is governed by
how retained the IONISED form still is, because the ionised peak sits close to
t0 and a t0 error swamps it:

    k_neutral   2 ->  SD(D) 0.070
                5 ->        0.026
               10 ->        0.014
               20 ->        0.007
               40 ->        0.003

A twenty-fold swing. The protocol's tier-1 list is "suggested", never checked
against this. A weakly retained acid contributes almost nothing to the paired
difference, so it is worth knowing BEFORE ~230 injections are booked.

Predicts k_neutral on the low-silanol arm (Kinetex EVO C18, in WSU-2019) at each
ACN composition, using the slice's SoluteML descriptor predictions.

numpy + stdlib.
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
# Copied rather than imported from prototype/thin-slice/slice.py: that module
# imports scipy at load time and scipy is not installed here, so importing it
# would make this script unrunnable for an unrelated reason. Values and formula
# are identical - McGowan characteristic volume, Abraham's V.
MCGOWAN = {"C": 16.35, "H": 8.71, "O": 12.43, "N": 14.39, "F": 10.48,
           "Cl": 20.95, "Br": 26.21, "I": 34.53, "S": 22.91, "P": 24.87}


def mcgowan_v(smiles):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    n_atoms = mol.GetNumAtoms()
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    return (sum(MCGOWAN[a.GetSymbol()] for a in mol.GetAtoms())
            - 6.56 * (n_atoms - 1 + n_rings)) / 100.0

WSU = os.path.join(ROOT, "sources/wsu-lser/wsu2019-system-constants.csv")
PREDS = os.path.join(ROOT, "prototype/thin-slice/drug_preds_raw.csv")
COLUMN = "Kinetex EVO C18"

# The protocol's suggested tier 1. SMILES for the four we hold descriptors for;
# the rest are listed so the gap is explicit rather than silent.
TIER1 = {
    "Benzoic acid":   "OC(=O)c1ccccc1",
    "Ibuprofen":      "CC(C)Cc1ccc(cc1)C(C)C(O)=O",
    "Naproxen":       "COc1ccc2cc(ccc2c1)C(C)C(O)=O",
    "Ketoprofen":     "CC(C(O)=O)c1cccc(c1)C(=O)c1ccccc1",
}
NO_DESCRIPTORS = ["Salicylic acid", "Diclofenac", "4-Nitrobenzoic acid", "Warfarin"]

# SD(D) as a function of k_neutral, measured by decide_columns.py
SD_D = {2: 0.070, 5: 0.026, 10: 0.014, 20: 0.007, 40: 0.003}


def sd_d_at(k):
    """Log-linear interpolation of the measured SD(D) curve."""
    ks = np.array(sorted(SD_D))
    vs = np.array([SD_D[k] for k in ks])
    return float(np.exp(np.interp(np.log(max(k, 1e-6)),
                                  np.log(ks), np.log(vs))))


def main():
    rows = [r for r in csv.DictReader(open(WSU, encoding="utf-8-sig"))
            if r["column"] == COLUMN and r["modifier"] == "acetonitrile"]
    coef = {r["phi_pct_vv"]: np.array([float(r[k]) for k in
                                       ("c", "e", "s", "a", "b", "v")])
            for r in rows}

    # ensemble mean of the SoluteML folds, per compound
    acc = defaultdict(list)
    for r in csv.DictReader(open(PREDS)):
        if r["compound"] in TIER1:
            acc[r["compound"]].append([float(r[k]) for k in ("E", "S", "A", "B")])
    desc = {c: np.mean(v, axis=0) for c, v in acc.items()}

    print(f"Tier-1 neutral retention on {COLUMN} (the low-silanol arm)")
    print("Predicted k_neutral = 10^(c + eE + sS + aA + bB + vV)\n")
    print(f"{'compound':<16}{'V':>6}" + "".join(f"{'k@'+p+'%':>10}" for p in
                                               ("30", "40", "50")) +
          f"{'SD(D)@30':>11}  verdict")
    verdicts = {}
    for name, smi in TIER1.items():
        if name not in desc:
            print(f"{name:<16}  no descriptors")
            continue
        E, S, A, B = desc[name]
        V = mcgowan_v(smi)
        x = np.array([1.0, E, S, A, B, V])
        ks = {p: float(10 ** (coef[p] @ x)) for p in ("30", "40", "50")}
        sd = sd_d_at(ks["30"])
        verdict = ("GOOD" if ks["30"] >= 20 else
                   "usable" if ks["30"] >= 10 else
                   "WEAK - contributes little" if ks["30"] >= 5 else
                   "REJECT - ionised peak sits on t0")
        verdicts[name] = (ks["30"], sd, verdict)
        print(f"{name:<16}{V:>6.3f}" +
              "".join(f"{ks[p]:>10.1f}" for p in ("30", "40", "50")) +
              f"{sd:>11.4f}  {verdict}")

    print(f"\nNot checkable - no descriptors held: {', '.join(NO_DESCRIPTORS)}")
    print("They must clear the same bar before the campaign is booked.")

    if verdicts:
        ks30 = np.array([v[0] for v in verdicts.values()])
        print(f"\nOf the {len(verdicts)} checkable acids at 30% ACN: "
              f"median k_neutral {np.median(ks30):.1f}, "
              f"min {ks30.min():.1f}, max {ks30.max():.1f}")
        print("Protocol asks tier 1 to run at 30% ACN - the most retentive of the")
        print("three ACN points - which is what these numbers are quoted at.")

    # what the whole tier buys, if the paired difference is averaged over 8
    print("\nImplied paired-difference precision if all 8 matched the checkable set:")
    for label, k in (("all at the observed median", float(np.median(ks30))),
                     ("all at the observed minimum", float(ks30.min()))):
        sd_pair = np.sqrt(2) * sd_d_at(k) / np.sqrt(8)
        print(f"  {label:<30} SD {sd_pair:.4f}, min detectable "
              f"{1.96 * sd_pair:.4f}  (prior SD on D_run = 0.30)")


if __name__ == "__main__":
    main()
