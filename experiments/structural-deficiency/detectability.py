#!/usr/bin/env python3
"""#46: are the structural-deficiency classes detectable from what we hold?

#42 split doubt into two kinds with opposite remedies:

  predictor uncertainty  - the ensemble is unsure. SPENDABLE.
  structural deficiency  - the Abraham representation cannot encode what the
                           molecule does. NOT spendable: a better predictor
                           cannot help, because the five numbers are the wrong
                           five numbers.

and claimed the ensemble is **blind to the second by construction** - its 25
networks agree because they learned the same inadequate representation from the
same corpus. That claim has a testable shape:

    For a structurally deficient class, the RESIDUAL against measured
    descriptors should be elevated while the ENSEMBLE SPREAD is not.

Elevated residual with elevated spread is ordinary predictor uncertainty and
needs no register. Elevated residual with FLAT spread is the register's whole
justification. Neither elevated retires the class.

Data on hand: WSU-94 measured descriptors (the reference lineage), SoluteML
25-fold predictions for the same compounds, and their SMILES.

⚠ CONTAMINATION, stated first. #42 established the WSU-94 are almost certainly
inside SoluteDB, so these residuals are in-sample-adjacent and too small. But
the direction of that bias matters: contamination SHRINKS residuals. A class
that still shows elevated residuals despite being partly memorised is showing
them against a headwind - so a positive result here is stronger than it looks,
and a null result is weaker. This test can support the register; it cannot
retire it.

Run with env/model-venv/bin/python.
"""
import csv
import os
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
S2 = os.path.join(ROOT, "sources/wsu-lser/wsu2019-descriptors-s2.csv")
SMI = os.path.join(ROOT, "prototype/descriptor-covariance/wsu94_smiles.csv")
PRED = os.path.join(ROOT, "prototype/descriptor-covariance/soluteml_preds_raw.csv")

TARGETS = ("E", "S", "A", "B")

# The register's classes, as SMARTS. Each is an operational GUESS at a class the
# sources named in prose; where a class cannot be written as substructure at all,
# that is itself a finding and it is recorded as None.
CLASSES = {
    # ortho -OH or -NH next to a carbonyl/nitro acceptor: the classic six-ring
    # internal hydrogen bond (salicylates, o-nitrophenol)
    "intramolecular H-bond": [
        "[OX2H]c1ccccc1[CX3]=[OX1]", "[OX2H]c1ccccc1[NX3](=O)=O",
        "[OX2H]c1ccccc1[OX2]", "[NX3;H1,H2]c1ccccc1[CX3]=[OX1]",
    ],
    # a donor flanked on both sides - the single A cannot carry site accessibility
    "hindered H-bond donor": [
        "[OX2H]c1c([CX4])cccc1[CX4]", "[OX2H]C([CX4])([CX4])[CX4]",
    ],
    # any H-bond acid at all (WSU Table S-1 named 'some hydrogen-bond acids')
    "H-bond acid (any donor)": ["[OX2H]", "[NX3;H1,H2]", "[SX2H]"],
    # 'bulky' has no operational definition yet - that is #43. V is the honest
    # stand-in and its arbitrariness is the point.
    "bulky (V > 1.4, proxy)": None,
    # several distinct interactions collapse into S: aromatic + strong dipole
    "dipolarity mixture": [
        "[NX3](=O)=O", "[CX3](=[OX1])[NX3]", "[SX4](=[OX1])(=[OX1])",
        "C#N",
    ],
}


def load():
    meas = {}
    for r in csv.DictReader(open(S2, encoding="utf-8-sig")):
        try:
            meas[r["compound"].strip()] = {k: float(r[k]) for k in TARGETS}
        except ValueError:
            pass
    smi = {r["compound"].strip(): r["smiles"] for r in csv.DictReader(open(SMI))}
    folds = defaultdict(list)
    for r in csv.DictReader(open(PRED)):
        folds[r["compound"].strip()].append([float(r[k]) for k in TARGETS])
    return meas, smi, {k: np.array(v) for k, v in folds.items()}


def classify(smiles, v):
    out = set()
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return out
    for name, pats in CLASSES.items():
        if pats is None:
            if v is not None and v > 1.4:
                out.add(name)
            continue
        for p in pats:
            q = Chem.MolFromSmarts(p)
            if q is not None and m.HasSubstructMatch(q):
                out.add(name)
                break
    return out


def mcgowan(smiles):
    MG = {"C": 16.35, "H": 8.71, "O": 12.43, "N": 14.39, "F": 10.48,
          "Cl": 20.95, "Br": 26.21, "I": 34.53, "S": 22.91, "P": 24.87}
    from rdkit.Chem import rdMolDescriptors
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return (sum(MG.get(a.GetSymbol(), 0) for a in m.GetAtoms())
            - 6.56 * (m.GetNumAtoms() - 1 + rdMolDescriptors.CalcNumRings(m))) / 100.0


def main():
    meas, smi, folds = load()
    names = [n for n in meas if n in smi and n in folds]
    print(f"{len(names)} compounds with measured descriptors, SMILES and 25 folds\n")

    resid, spread, cls = {}, {}, {}
    for n in names:
        f = folds[n]
        mu, sd = f.mean(axis=0), f.std(axis=0, ddof=1)
        resid[n] = np.array([mu[i] - meas[n][t] for i, t in enumerate(TARGETS)])
        spread[n] = sd
        cls[n] = classify(smi[n], mcgowan(smi[n]))

    print("Per-class test. RATIO = class median / non-class median.\n")
    print(f"  {'class':<26}{'n':>4}  " +
          "".join(f"{t+' resid':>10}{t+' spread':>11}" for t in ("S", "A")))
    rows = {}
    for c in CLASSES:
        inn = [n for n in names if c in cls[n]]
        out = [n for n in names if c not in cls[n]]
        if len(inn) < 4:
            print(f"  {c:<26}{len(inn):>4}   (too few to test)")
            continue
        cells, keep = [], {}
        for t in ("S", "A"):
            i = TARGETS.index(t)
            rr = (np.median([abs(resid[n][i]) for n in inn]) /
                  np.median([abs(resid[n][i]) for n in out]))
            sr = (np.median([spread[n][i] for n in inn]) /
                  np.median([spread[n][i] for n in out]))
            keep[t] = (rr, sr)
            cells += [f"{rr:>10.2f}", f"{sr:>11.2f}"]
        rows[c] = keep
        print(f"  {c:<26}{len(inn):>4}  " + "".join(cells))

    print("\n  Reading a row: residual ratio > 1 with spread ratio ~ 1 is the")
    print("  register's signature - wrong, and the ensemble does not know it.")
    print("  Both > 1 is ordinary predictor uncertainty, which #42 already handles")
    print("  and which needs no register.\n")

    # The sharper statistic is the RATIO OF RATIOS: does the residual grow
    # faster than the ensemble's own spread? A class where both rise together is
    # one the ensemble already signals - it needs no register. A class where the
    # error outruns the spread is one the ensemble is blind to, in proportion.
    print("Verdict per class (resid ratio / spread ratio - how far the ensemble")
    print("under-signals its own error):")
    print(f"  {'class':<26}{'S':>8}{'A':>8}   reading")
    for c, k in rows.items():
        rs = {t: rr / sr for t, (rr, sr) in k.items()}
        worst = max(rs.values())
        if worst > 1.5:
            read = "ensemble UNDER-SIGNALS - register earns its place"
        elif worst > 1.15:
            read = "mild under-signal"
        else:
            read = "ensemble tracks its error - no register needed"
        print(f"  {c:<26}{rs['S']:>8.2f}{rs['A']:>8.2f}   {read}")

    # how much of the set does the register capture?
    flagged = [n for n in names if cls[n]]
    print(f"\nCoverage: the register flags {len(flagged)}/{len(names)} "
          f"({100*len(flagged)/len(names):.0f}%) of the WSU-94.")
    print("  A register that flags most of the set is not a register, it is a")
    print("  blanket caveat - and it would move no decision.")


if __name__ == "__main__":
    main()
