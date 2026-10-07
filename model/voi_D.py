#!/usr/bin/env python3
"""What is measuring `D` actually worth, and which HALF of it should we buy?

`D` owns 0.78-0.92 of the variance for ionised compounds and has never been
fitted. #30's sweep exists to measure it. But `D` is two terms, and they do
different damage:

  D = 1.5 + D_run + D_c     D_run ~ N(0, 0.30)   per-run  (column, modifier)
                            D_c   ~ N(0, 0.40)   per-compound

`CONTEXT.md` is explicit that per-run terms are "one draw shared by every
compound and **largely cancel in order** while still setting absolute
retention", whereas per-compound terms are "drawn independently per compound and
**scramble elution order**". Elution order is the primary metric.

If that holds here, then the sweep's headline output - `D`'s per-run part, as
the paired EVO-minus-Supelcosil difference - buys little for the metric the
product is graded on, while `D_c`, which the same runs measure per compound,
buys most of it. That would say: **prioritise compounds over conditions.**

This measures it by collapsing each term's uncertainty independently and
watching the confident-pair fraction. numpy only, 5 seeds.

Run with env/ionisation-venv/bin/python.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "prototype", "thin-slice"))

import slice as sl          # noqa: E402
from ionisation import Provider  # noqa: E402
from sampler import confident, run  # noqa: E402

SEEDS = [20260815, 7, 101, 2026, 31337]


def sweep(prov, label, **kw):
    """Confident pairs and median sd(log k) per operating point, over seeds."""
    per_op = {tuple(o): ([], []) for o in sl.OPS}
    for sd in SEEDS:
        names, out, refused = run(prov=prov, seed=sd, **kw)
        for o in sl.OPS:
            lk = out[tuple(o)]["logk"]
            c, n = confident(names, lk)
            per_op[tuple(o)][0].append(c)
            live = [i for i in range(len(names)) if not np.isnan(lk[:, i]).all()]
            per_op[tuple(o)][1].append(
                float(np.median(np.std(lk[:, live], axis=0))))
    return {o: (np.mean(v[0]), np.std(v[0]), np.mean(v[1]))
            for o, v in per_op.items()}


def main():
    prov = Provider()
    print("Value of measuring D: which half is worth the instrument time?\n")
    print("Prior: D_run ~ N(0, 0.30) per-run, D_c ~ N(0, 0.40) per-compound\n")

    cases = [
        ("prior (nothing measured)",     dict()),
        ("D_run measured (-> 0.05)",     dict(d_run_sd=0.05)),
        ("D_c measured   (-> 0.05)",     dict(d_c_sd=0.05)),
        ("both measured  (-> 0.05)",     dict(d_run_sd=0.05, d_c_sd=0.05)),
        ("D known exactly (floor)",      dict(d_run_sd=0.0, d_c_sd=0.0)),
    ]
    res = {lbl: sweep(prov, lbl, **kw) for lbl, kw in cases}

    for o in sl.OPS:
        print(f"  operating point {tuple(o)}")
        print(f"    {'case':<30}{'confident/105':>16}{'median sd(log k)':>20}")
        base = res["prior (nothing measured)"][tuple(o)][0]
        for lbl, _ in cases:
            m, s_, sdk = res[lbl][tuple(o)]
            delta = "" if lbl.startswith("prior") else f"   {m - base:+.1f}"
            print(f"    {lbl:<30}{f'{m:.1f} ± {s_:.1f}':>16}{sdk:>20.4f}{delta}")
        print()

    print("Reading it")
    print("  If 'D_c measured' recovers most of what 'both measured' recovers,")
    print("  order confidence is limited by the PER-COMPOUND term, and the sweep")
    print("  should be sized for COMPOUNDS rather than conditions. If 'D_run")
    print("  measured' does, the two-column contrast is where the value is.")
    print("  sd(log k) tracks absolute retention, which per-run terms DO move -")
    print("  so the two columns can differ between the metrics.")


if __name__ == "__main__":
    main()
