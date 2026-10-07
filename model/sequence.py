#!/usr/bin/env python3
"""Compound plate and mixture pooling for #47, sized on the LIMITING column.

Two arms, two plate counts: Kinetex EVO C18 (core-shell) 19,553 and Supelcosil
LC-18 (fully porous) 11,447. Peaks are 24% wider on the Supelcosil arm, so a
pool that resolves on EVO can co-elute there. **Pooling must therefore be sized
on Supelcosil**, and running the same pools on both arms is what keeps the
paired difference paired.

The resolution criterion, for two isocratic peaks with sigma = t_R/sqrt(N):

    Rs = (t2 - t1) * sqrt(N) / (2 * (t1 + t2))

so Rs >= 1.5 requires  (t2 - t1) / (t2 + t1) >= 3 / sqrt(N), i.e. adjacent peaks
in a pool must differ by about **5.6% in retention time** on the limiting arm.

**Retention is predicted with the EVO LSER constants and the SUPELCOSIL plate
count.** That is not a fudge: #30 matched the two arms on hydrophobicity to
H = 0.002, precisely so that `k` would be comparable and the arms would differ
in silanol activity and efficiency rather than retention strength. Supelcosil is
not in WSU-2019, so it has no constants of its own - which is the same asymmetry
the tier structure already handles.

Run with env/ionisation-venv/bin/python.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "prototype", "thin-slice"))

import slice as sl                    # noqa: E402
from columns import get               # noqa: E402
from compounds import COMPOUNDS       # noqa: E402
from dispersion import peak_sd_min, resolution  # noqa: E402

FLOW = 1.0
RS_TARGET = 1.5

# The tier-1 plate, fixed on #30. k is measured only where SoluteML descriptors
# exist; the rest were screened on McGowan V and still need confirming.
TIER1 = [
    ("Ibuprofen",    "CC(C)Cc1ccc(cc1)C(C)C(O)=O",              4.45, 25.3),
    ("Naproxen",     "COc1ccc2cc(ccc2c1)C(C)C(O)=O",            4.18, 14.3),
    ("Ketoprofen",   "CC(C(O)=O)c1cccc(c1)C(=O)c1ccccc1",       4.20, 10.7),
    ("Flurbiprofen", "CC(C(O)=O)c1ccc(-c2ccccc2)c(F)c1",        4.20, None),
    ("Fenoprofen",   "CC(C(O)=O)c1cccc(Oc2ccccc2)c1",           4.50, None),
    ("Diclofenac",   "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",          4.15, None),
    ("Gemfibrozil",  "CC(C)(CCCOc1ccc(C)cc1C)C(O)=O",           4.70, None),
    ("Indomethacin", "COc1ccc2c(c1)c(CC(O)=O)c(C)n2C(=O)c1ccc(Cl)cc1", 4.50, None),
]


def min_separation_frac(n_plates, rs=RS_TARGET):
    """Required (t2-t1)/(t2+t1) for baseline resolution on a column of N plates."""
    return 2.0 * rs / np.sqrt(n_plates)


def pool(items, n_plates, rs=RS_TARGET):
    """Greedy pooling: sort by retention, place each compound in the first pool
    whose nearest member it clears. Deterministic, and it fails loudly rather
    than silently over-filling a pool."""
    need = min_separation_frac(n_plates, rs)
    pools: list[list] = []
    for name, tr in sorted(items, key=lambda x: x[1]):
        for p in pools:
            if all(abs(tr - t) / (tr + t) >= need for _, t in p):
                p.append((name, tr))
                break
        else:
            pools.append([(name, tr)])
    return pools


def main():
    evo, sup = get("kinetex-evo-c18"), get("supelcosil-lc-18")
    n_evo, n_sup = evo.plates(FLOW), sup.plates(FLOW)
    t0 = evo.geometry.t0_min(FLOW)

    print("Mixture pooling sized on the LIMITING arm\n")
    print(f"  {'arm':<22}{'N':>9}{'min sep (t2-t1)/(t2+t1)':>26}{'at Rs':>7}")
    for nm, n in (("Kinetex EVO C18", n_evo), ("Supelcosil LC-18", n_sup)):
        print(f"  {nm:<22}{n:>9.0f}{min_separation_frac(n):>26.4f}{RS_TARGET:>7}")
    print(f"\n  Sizing on Supelcosil demands {min_separation_frac(n_sup):.1%} "
          f"separation vs {min_separation_frac(n_evo):.1%} on EVO —")
    print(f"  a {min_separation_frac(n_sup)/min_separation_frac(n_evo):.2f}x "
          f"stricter constraint, which is the whole reason to size on it.")

    # --- demonstrate on the compounds whose retention we actually predict
    print("\n" + "=" * 68)
    print("Pooling the thin-slice set (retention predicted, EVO constants)\n")
    names, out, refused = None, None, None
    from ionisation import Provider
    from sampler import run as run_model
    names, out, refused = run_model(prov=Provider(), seed=sl.SEED)
    op = (0.30, 3.0)
    logk = out[op]["logk"]
    items = []
    for i, nm in enumerate(names):
        if nm in refused:
            continue
        k = float(np.nanmedian(10.0 ** logk[:, i]))
        items.append((nm, t0 * (1.0 + k)))

    for label, n in (("EVO (19,553 plates)", n_evo),
                     ("SUPELCOSIL (11,447 plates)", n_sup)):
        ps = pool(items, n)
        print(f"  {label}: {len(ps)} pools needed")
        for j, p in enumerate(ps, 1):
            body = ", ".join(f"{nm} ({tr:.1f})" for nm, tr in p)
            print(f"    pool {j}: {body}")
        print()

    print("  The two arms need DIFFERENT pool counts. Use the Supelcosil pooling")
    print("  on both, so the same mixture is injected on each arm and the paired")
    print("  difference stays paired.")

    # --- the tier-1 plate
    print("=" * 68)
    print("Tier-1 plate: retention still unknown for 5 of 8\n")
    print(f"  {'acid':<16}{'pKa':>6}{'k @30% ACN':>13}{'t_R (min)':>12}")
    known = []
    for nm, smi, pka, k in TIER1:
        if k is None:
            print(f"  {nm:<16}{pka:>6.2f}{'?':>13}{'needs descriptors':>12}")
        else:
            tr = t0 * (1 + k)
            known.append((nm, tr))
            print(f"  {nm:<16}{pka:>6.2f}{k:>13.1f}{tr:>12.1f}")
    ps = pool(known, n_sup)
    print(f"\n  The 3 measured acids need {len(ps)} pool(s) on Supelcosil: "
          f"{[[nm for nm, _ in p] for p in ps]}")
    print("\n  ⚠ Final pooling is BLOCKED on descriptors for flurbiprofen,")
    print("    fenoprofen, diclofenac, gemfibrozil and indomethacin. Their V")
    print("    screen (1.84-2.53) says they are more retained than ibuprofen")
    print("    (1.78, k=25.3), so they will land LATER and likely spread across")
    print("    pools - but 'likely' is not a sequence. SoluteML on those five")
    print("    SMILES closes it; that is a GPU job, and env/install_qspr.sh exists.")


if __name__ == "__main__":
    main()
