#!/usr/bin/env python3
"""Measure the ionisation provider against the pKa values the slice hardcodes.

The slice assumes a per-compound pKa error of +-0.35 ("the spread itself
justifies it", from one disagreement between Avdeef and Clarke on ketoprofen).
That is an assumption about *literature* spread, not a measurement of any
provider. Meanwhile pKa owns 0.67 of the variance for partially ionised
compounds, so it is the dominant term for the in-scope chemistry.

This replaces the assumption with a measurement, on exactly the compounds in
play: run Uni-pKa against each cited literature value in compounds.py.

Also checks the two things the hardcoded version could not do at all - whether
the refusal tier can be DERIVED rather than hand-labelled, and whether
f_neutral from real microspecies matches the closed-form sigmoid the slice uses.

Run with env/ionisation-venv/bin/python.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "prototype", "thin-slice"))

from compounds import COMPOUNDS  # noqa: E402
from ionisation import Provider, crossover_pka  # noqa: E402


def sigmoid_f_neutral(ctype, pkas, ph):
    """The slice's closed form, for comparison."""
    if not pkas:
        return 1.0
    if ctype == "acid":
        return 1.0 / (1.0 + 10.0 ** (ph - pkas[0]))
    if ctype == "base":
        return 1.0 / (1.0 + 10.0 ** (pkas[0] - ph))
    if ctype == "ampho":
        return 1.0 / (1.0 + 10.0 ** (pkas[0] - ph) + 10.0 ** (ph - pkas[1]))
    return 1.0


def main():
    prov = Provider()
    ionisable = [c for c in COMPOUNDS if c["pka"]]
    print(f"Provider: Uni-pKa (simple SMARTS). {len(ionisable)} ionisable "
          f"compounds with cited literature pKa.\n")

    print("1. PROVIDER ERROR against the cited values the slice hardcodes")
    print(f"   {'compound':<22}{'lit pKa':>18}{'Uni-pKa':>18}{'error':>16}")
    errs = []
    for c in ionisable:
        lit = list(c["pka"])
        pred = crossover_pka(prov, c["smiles"])
        # match predicted to literature by nearest, only where counts allow
        if len(pred) >= len(lit):
            chosen = []
            pool = list(pred)
            for L in lit:
                b = min(pool, key=lambda p: abs(p - L))
                chosen.append(b)
                pool.remove(b)
        else:
            chosen = pred + [float("nan")] * (len(lit) - len(pred))
        e = [round(p - L, 2) for p, L in zip(chosen, lit)]
        errs += [x for x in e if x == x]
        print(f"   {c['name']:<22}{str(lit):>18}{str(chosen):>18}{str(e):>16}")
    a = np.array(errs)
    if a.size:
        print(f"\n   n={len(a)} pKa values | mean {a.mean():+.3f} | "
              f"SD {a.std(ddof=1):.3f} | MAE {np.abs(a).mean():.3f} | "
              f"max|e| {np.abs(a).max():.2f}")
    else:
        print("\n   no pKa recovered - crossover detector found nothing")
    print(f"   Slice assumes a per-compound pKa error of +-0.35 (1 SD).")

    print("\n2. REFUSAL TIER - derived, not hand-labelled")
    print(f"   {'compound':<24}{'slice label':>12}{'derived':>13}"
          f"{'max f_neutral':>15}{'at pH':>7}")
    agree = 0
    for c in COMPOUNDS:
        r = prov.refusal(c["smiles"])
        label = "refused" if c["type"] == "refused" else c["type"]
        derived = r["tier"]
        ok = ((label == "refused") == (derived == "refuse"))
        agree += ok
        flag = "" if ok else "   <-- DISAGREES"
        print(f"   {c['name']:<24}{label:>12}{derived:>13}"
              f"{r['max_f_neutral']:>15.4f}{str(r['at_ph']):>7}{flag}")
    print(f"\n   refusal agreement: {agree}/{len(COMPOUNDS)}")

    print("\n3. f_neutral - real microspecies vs the slice's closed form")
    print(f"   {'compound':<22}{'pH':>6}{'sigmoid':>10}{'provider':>10}{'diff':>9}")
    diffs = []
    for c in ionisable:
        for ph in (3.0, 5.0, 7.0):
            s = sigmoid_f_neutral(c["type"], c["pka"], ph)
            p = prov.f_neutral(c["smiles"], ph)
            diffs.append(p - s)
            print(f"   {c['name']:<22}{ph:>6.1f}{s:>10.4f}{p:>10.4f}{p - s:>9.4f}")
    d = np.array(diffs)
    print(f"\n   n={len(d)} | mean {d.mean():+.4f} | max|diff| {np.abs(d).max():.4f}")
    print("   A large diff means the hardcoded pKa and the closed form together")
    print("   misstate the neutral fraction that multiplies through to k.")


if __name__ == "__main__":
    main()
