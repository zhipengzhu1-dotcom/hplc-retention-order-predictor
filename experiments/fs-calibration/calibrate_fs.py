#!/usr/bin/env python3
"""#26, rescoped: calibrate `Fs` into an actual substitution error.

#26's primary experiment - fitting our own LSER constants to HSM3's raw
retention data - cannot run, because that dataset is not obtainable. But #26
lists a second deliverable that CAN be built from data in hand:

    "it is the material needed to calibrate `Fs` similarity into an actual
     confidence rather than an uncalibrated heuristic"

21 columns carry BOTH HSM parameters (H, S*, A, B, C) and WSU-2019 LSER system
constants (c, e, s, a, b, v). For every pair of them we can ask the question a
method developer actually asks:

    if I substitute column B for column A because `Fs` said they were similar,
    how wrong is my retention prediction?

`Fs` is computed from the HSM side; the substitution error is computed from the
LSER side over the 94 WSU descriptor compounds. Neither number knows about the
other, so the correlation between them is a genuine out-of-model test.

#36 already showed by DECOMPOSITION that `Fs`-nearest is not LSER-nearest
(80.5% of `Fs`^2 is the C term, which has no LSER analogue). This measures the
same claim directly, and turns it into a number a user could be shown.

numpy + stdlib. Deterministic.
"""
import csv
import os
import re
import unicodedata

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
WSU = os.path.join(ROOT, "sources/wsu-lser/wsu2019-system-constants.csv")
DESC = os.path.join(ROOT, "sources/wsu-lser/wsu2019-descriptors-s2.csv")
HSM = os.path.join(ROOT, "sources/hsm-column-db/database.csv")

# Published Fs weights (Snyder-Dolan-Carr), as used in #36's decomposition:
# Fs^2 = (12.5 dH)^2 + (100 dS*)^2 + (30 dA)^2 + (143 dB)^2 + (83 dC)^2
FS_W = dict(H=12.5, S=100.0, A=30.0, B=143.0, C=83.0)


def norm(n):
    s = unicodedata.normalize("NFKC", n).lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip().replace(" ", "")


def read(path, enc="utf-8-sig"):
    with open(path, newline="", encoding=enc, errors="replace") as f:
        return list(csv.DictReader(f))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def main():
    wsu = read(WSU)
    desc = read(DESC)
    hsm = read(HSM, "latin-1")

    hmap = {}
    for r in hsm:
        hmap.setdefault(norm(r["name"]), r)

    # descriptor matrix for the 94 anchor compounds: [1, E, S, A, B, V]
    D, names = [], []
    for r in desc:
        vals = [fnum(r[k]) for k in ("E", "S", "A", "B", "V")]
        if all(v is not None for v in vals):
            D.append([1.0] + vals)
            names.append(r["compound"].strip())
    D = np.array(D)
    print(f"Anchor set: {len(D)} WSU compounds (median V "
          f"{np.median(D[:, 5]):.2f})")

    print(f"\n{'modifier':<12}{'phi':>5}{'cols':>6}{'pairs':>7}"
          f"{'pearson':>9}{'spearman':>10}")
    results = {}
    for modifier in ("methanol", "acetonitrile"):
        for phi in ("30", "40", "50"):
            rows = [r for r in wsu
                    if r["modifier"] == modifier and r["phi_pct_vv"] == phi
                    and norm(r["column"]) in hmap]
            # one fit per column at this condition
            byc = {}
            for r in rows:
                coef = [fnum(r[k]) for k in ("c", "e", "s", "a", "b", "v")]
                if all(c is not None for c in coef):
                    byc[r["column"]] = np.array(coef)
            cols = sorted(byc)
            if len(cols) < 8:
                continue
            fs_list, err_list = [], []
            for i, a in enumerate(cols):
                for b in cols[i + 1:]:
                    ha, hb = hmap[norm(a)], hmap[norm(b)]
                    d = {}
                    ok = True
                    for key, col in (("H", "H"), ("S", "S"), ("A", "A"),
                                     ("B", "B"), ("C", "C28")):
                        va, vb = fnum(ha[col]), fnum(hb[col])
                        if va is None or vb is None:
                            ok = False
                            break
                        d[key] = va - vb
                    if not ok:
                        continue
                    fs = np.sqrt(sum((FS_W[k] * d[k]) ** 2 for k in FS_W))
                    # substitution error: RMS difference in predicted log k
                    err = float(np.sqrt(np.mean((D @ byc[a] - D @ byc[b]) ** 2)))
                    fs_list.append(fs)
                    err_list.append(err)
            if len(fs_list) < 20:
                continue
            fs_arr, err_arr = np.array(fs_list), np.array(err_list)
            pear = float(np.corrcoef(fs_arr, err_arr)[0, 1])

            def rank(v):
                o = np.argsort(v)
                r = np.empty_like(o, dtype=float)
                r[o] = np.arange(len(v))
                return r
            spear = float(np.corrcoef(rank(fs_arr), rank(err_arr))[0, 1])
            results[(modifier, phi)] = (fs_arr, err_arr)
            print(f"{modifier:<12}{phi:>5}{len(cols):>6}{len(fs_arr):>7}"
                  f"{pear:>9.3f}{spear:>10.3f}")

    # ---- the calibration table: given Fs, what error should you expect? -
    # Fs <= 3 is the published "equivalent" threshold, but among these 21
    # deliberately diverse columns only one pair reaches it - so the useful
    # output is error by Fs BAND, which is what a user would actually consult.
    print("\nCALIBRATION - substitution error in log k, by `Fs` band")
    print("  (pooled over conditions; error = RMS dlog k over the 94 anchor compounds)")
    all_fs = np.concatenate([v[0] for v in results.values()])
    all_er = np.concatenate([v[1] for v in results.values()])
    bands = [(0, 5), (5, 10), (10, 20), (20, 40), (40, np.inf)]
    print(f"  {'Fs band':<12}{'pairs':>7}{'median':>9}{'90th pct':>10}{'max':>8}")
    for lo, hi in bands:
        m = (all_fs >= lo) & (all_fs < hi)
        if m.sum() == 0:
            continue
        lbl = f"{lo}-{hi:g}" if np.isfinite(hi) else f"{lo}+"
        print(f"  {lbl:<12}{m.sum():>7}{np.median(all_er[m]):>9.3f}"
              f"{np.percentile(all_er[m], 90):>10.3f}{all_er[m].max():>8.3f}")
    print(f"  {'ALL':<12}{len(all_er):>7}{np.median(all_er):>9.3f}"
          f"{np.percentile(all_er, 90):>10.3f}{all_er.max():>8.3f}")
    print("\n  Compare: LSER lack-of-fit floor is 0.030 (ACN) / 0.041 (MeOH) log k.")
    print("  Even the closest band sits an order of magnitude above that floor.")

    # ---- the axis the project actually grades on -----------------------
    # CONTEXT.md on topology-class substitution: the cost is "SD ~ 0.055-0.145 in
    # log k, but only ~0.007-0.018 on a close pair's separation - affordable ONLY
    # because elution order is the primary metric, and never to be quoted as
    # affordable on log k". The log k table above is therefore the misleading
    # half. This is the half that decides whether Fs substitution is usable.
    print("\nELUTION-ORDER cost of substituting on `Fs`")
    print("  (fraction of the 4,371 compound pairs whose predicted order flips,")
    print("   and how much a CLOSE pair's separation moves)")
    print(f"  {'Fs band':<12}{'pairs':>7}{'median flip':>13}{'90th pct':>10}"
          f"{'med |dsep| close':>18}")
    order_rows = {}
    for (modifier, phi), _ in results.items():
        rows = [r for r in wsu if r["modifier"] == modifier
                and r["phi_pct_vv"] == phi and norm(r["column"]) in hmap]
        byc = {}
        for r in rows:
            coef = [fnum(r[k]) for k in ("c", "e", "s", "a", "b", "v")]
            if all(c is not None for c in coef):
                byc[r["column"]] = np.array(coef)
        cols = sorted(byc)
        iu = np.triu_indices(len(D), k=1)
        for i, a in enumerate(cols):
            for b in cols[i + 1:]:
                ha, hb = hmap[norm(a)], hmap[norm(b)]
                try:
                    d = {k: fnum(ha[c]) - fnum(hb[c]) for k, c in
                         (("H", "H"), ("S", "S"), ("A", "A"), ("B", "B"), ("C", "C28"))}
                except TypeError:
                    continue
                fs = np.sqrt(sum((FS_W[k] * d[k]) ** 2 for k in FS_W))
                la, lb = D @ byc[a], D @ byc[b]
                da = (la[:, None] - la[None, :])[iu]
                db = (lb[:, None] - lb[None, :])[iu]
                flip = float(np.mean((da > 0) != (db > 0)))
                close = np.abs(da) < 0.10          # pairs that nearly co-elute on A
                dsep = (float(np.median(np.abs(db[close] - da[close])))
                        if close.sum() else np.nan)
                order_rows.setdefault("fs", []).append(fs)
                order_rows.setdefault("flip", []).append(flip)
                order_rows.setdefault("dsep", []).append(dsep)
    ofs = np.array(order_rows["fs"])
    ofl = np.array(order_rows["flip"])
    ods = np.array(order_rows["dsep"])
    for lo, hi in bands:
        m = (ofs >= lo) & (ofs < hi)
        if m.sum() == 0:
            continue
        lbl = f"{lo}-{hi:g}" if np.isfinite(hi) else f"{lo}+"
        print(f"  {lbl:<12}{m.sum():>7}{np.median(ofl[m]):>12.1%}"
              f"{np.percentile(ofl[m], 90):>10.1%}"
              f"{np.nanmedian(ods[m]):>18.3f}")
    print(f"  {'ALL':<12}{len(ofl):>7}{np.median(ofl):>12.1%}"
          f"{np.percentile(ofl, 90):>10.1%}{np.nanmedian(ods):>18.3f}")
    print("\n  Budgeted cost of topology-class substitution (CONTEXT.md): 0.007-0.018")
    print("  on a close pair's separation. Compare the last column.")

    # ---- the decision-relevant question --------------------------------
    print("\nIf you pick the `Fs`-NEAREST column, how often is it the")
    print("LSER-nearest - i.e. the one that actually predicts retention best?")
    print(f"  {'condition':<20}{'cols':>6}{'Fs-nearest = best':>20}{'in best 3':>11}"
          f"{'median rank':>13}")
    for modifier in ("methanol", "acetonitrile"):
        for phi in ("30", "40", "50"):
            key = (modifier, phi)
            if key not in results:
                continue
            rows = [r for r in wsu if r["modifier"] == modifier
                    and r["phi_pct_vv"] == phi and norm(r["column"]) in hmap]
            byc = {}
            for r in rows:
                coef = [fnum(r[k]) for k in ("c", "e", "s", "a", "b", "v")]
                if all(c is not None for c in coef):
                    byc[r["column"]] = np.array(coef)
            cols = sorted(byc)
            hits, top3, ranks = 0, 0, []
            for a in cols:
                others = [b for b in cols if b != a]
                fs, er = [], []
                for b in others:
                    ha, hb = hmap[norm(a)], hmap[norm(b)]
                    try:
                        d = {k: fnum(ha[c]) - fnum(hb[c]) for k, c in
                             (("H", "H"), ("S", "S"), ("A", "A"),
                              ("B", "B"), ("C", "C28"))}
                    except TypeError:
                        continue
                    fs.append(np.sqrt(sum((FS_W[k] * d[k]) ** 2 for k in FS_W)))
                    er.append(np.sqrt(np.mean((D @ byc[a] - D @ byc[b]) ** 2)))
                if len(fs) < 5:
                    continue
                fs, er = np.array(fs), np.array(er)
                pick = int(np.argmin(fs))
                order = np.argsort(er)
                rank_of_pick = int(np.where(order == pick)[0][0]) + 1
                ranks.append(rank_of_pick)
                hits += rank_of_pick == 1
                top3 += rank_of_pick <= 3
            if ranks:
                print(f"  {modifier+' '+phi+'%':<20}{len(ranks):>6}"
                      f"{f'{hits}/{len(ranks)}':>20}{f'{top3}/{len(ranks)}':>11}"
                      f"{np.median(ranks):>13.0f}")


if __name__ == "__main__":
    main()
