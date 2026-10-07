#!/usr/bin/env python3
"""PROTOTYPE (throwaway) — issue #24: phi interpolation vs LSS reparameterisation.

Answers, on the real WSU-2019 tables:
  Route A: interpolate system constants c,e,s,a,b,v across phi (linear vs PCHIP),
           leave-one-phi-out reconstruction error in log10 k.
  Route B: LSS fit log10 k = log k_w - S*phi per (compound, column, modifier),
           residuals inside vs outside the published 30-60% validity window.
  Plus: adjacent-phi coefficient correlation, and per-coefficient curvature
        (mid-point deviation from the chord, in units of the point's own sd).

Conventions (CONTEXT.md): WSU-2019 is log10 k, NOT ln. Dewetting mask (#40):
only incomplete_wetting == 0 rows are support; masks are contiguous at the low-phi
end (verified at load), so support starts at the lowest unmasked phi and no
curve crosses a dewetting floor.

Run:  python3 analyze.py     (numpy only; PCHIP implemented inline)
"""
import csv, os, sys
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
COEFS = ["c", "e", "s", "a", "b", "v"]
QSPR_BUDGET = 0.09  # per-compound log k error already in the budget (CONTEXT.md)

# ---------------------------------------------------------------- data loading
def load():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "sources/wsu-lser/wsu2019-system-constants.csv"))))
    systems = defaultdict(list)  # (modifier, column) -> list of dict rows, unmasked only
    n_masked = 0
    for r in rows:
        if r["incomplete_wetting"] == "1":
            n_masked += 1
            continue
        systems[(r["modifier"], r["column"])].append(r)
    assert n_masked == 22, n_masked
    out = {}
    for key, rs in systems.items():
        rs.sort(key=lambda r: float(r["phi_pct_vv"]))
        phi = np.array([float(r["phi_pct_vv"]) for r in rs]) / 100.0
        C = np.array([[float(r[k]) for k in COEFS] for r in rs])          # (n_phi, 6)
        SD = np.array([[float(r["sd_" + k]) for k in COEFS] for r in rs]) # (n_phi, 6)
        out[key] = (phi, C, SD)
    # verify mask contiguity (support starts at lowest unmasked phi)
    allrows = defaultdict(list)
    for r in rows:
        allrows[(r["modifier"], r["column"])].append((float(r["phi_pct_vv"]), r["incomplete_wetting"]))
    for key, v in allrows.items():
        v.sort()
        masks = [m for _, m in v]
        if "1" in masks:
            last = max(i for i, m in enumerate(masks) if m == "1")
            assert all(m == "1" for m in masks[: last + 1]), key
    comps = list(csv.DictReader(open(os.path.join(ROOT, "sources/wsu-lser/wsu2019-descriptors-s2.csv"))))
    names = [c["compound"] for c in comps]
    D = np.array([[float(c[k]) for k in ["E", "S", "A", "B", "V"]] for c in comps])  # (94, 5)
    return out, names, D

def logk(coefrow, D):
    """coefrow: (6,) = c,e,s,a,b,v ; D: (n,5) = E,S,A,B,V -> (n,) log10 k."""
    return coefrow[0] + D @ coefrow[1:]

# ---------------------------------------------------------------- PCHIP inline
def pchip_slopes(x, y):
    h = np.diff(x); d = np.diff(y) / h
    n = len(x); m = np.zeros(n)
    for i in range(1, n - 1):
        if d[i - 1] == 0 or d[i] == 0 or np.sign(d[i - 1]) != np.sign(d[i]):
            m[i] = 0.0
        else:
            w1 = 2 * h[i] + h[i - 1]; w2 = h[i] + 2 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
    def end(h0, h1, d0, d1):
        t = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if np.sign(t) != np.sign(d0): t = 0.0
        elif np.sign(d0) != np.sign(d1) and abs(t) > 3 * abs(d0): t = 3 * d0
        return t
    m[0] = end(h[0], h[1], d[0], d[1]) if n > 2 else d[0]
    m[-1] = end(h[-1], h[-2], d[-1], d[-2]) if n > 2 else d[-1]
    return m

def pchip_eval(x, y, m, xq):
    i = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    h = x[i + 1] - x[i]; t = (xq - x[i]) / h
    h00 = (1 + 2 * t) * (1 - t) ** 2; h10 = t * (1 - t) ** 2
    h01 = t * t * (3 - 2 * t); h11 = t * t * (t - 1)
    return h00 * y[i] + h10 * h * m[i] + h01 * y[i + 1] + h11 * h * m[i + 1]

def interp(x, y, xq, kind):
    if kind == "linear" or len(x) < 3:
        return np.interp(xq, x, y)
    return pchip_eval(x, y, pchip_slopes(x, y), np.atleast_1d(xq))[0] if np.isscalar(xq) \
        else pchip_eval(x, y, pchip_slopes(x, y), xq)

def stats(errs):
    a = np.abs(np.asarray(errs))
    return len(a), float(np.median(a)), float(np.percentile(a, 90)) if len(a) else float("nan")

def inwin(phi):  # published LSS validity window, phi as fraction
    return 0.30 - 1e-9 <= phi <= 0.60 + 1e-9

# ---------------------------------------------------------------------- main
def main():
    systems, names, D = load()
    print(f"systems (unmasked): {len(systems)}; compounds: {len(names)}")

    # ---------------- Route A: leave-one-phi-out interpolation of constants
    # errs[kind][(modifier, region)] -> list of per-compound dlogk at held-out phi
    A_errs = {k: defaultdict(list) for k in ("linear", "pchip")}
    A_coef_errs = {k: [] for k in ("linear", "pchip")}  # (coef_index, err/sd) diagnostics
    for (mod, col), (phi, C, SD) in systems.items():
        if len(phi) < 3:
            continue
        for i in range(1, len(phi) - 1):  # interior points only
            keep = np.arange(len(phi)) != i
            truth = logk(C[i], D)
            for kind in ("linear", "pchip"):
                chat = np.array([interp(phi[keep], C[keep, j], phi[i], kind) for j in range(6)])
                err = logk(chat, D) - truth
                region = "inside" if inwin(phi[i]) else "outside"
                A_errs[kind][(mod, region)].extend(err.tolist())
                A_errs[kind][(mod, "all")].extend(err.tolist())

    print("\n== Route A: leave-one-phi-out reconstruction |dlog10 k| (interior phi) ==")
    for kind in ("linear", "pchip"):
        for key in sorted(A_errs[kind]):
            n, med, p90 = stats(A_errs[kind][key])
            print(f"  {kind:6s} {key[0]:15s} {key[1]:7s} n={n:6d} median={med:.4f} p90={p90:.4f}")

    # ---------------- Route B: LSS per (compound, column, modifier)
    # residuals of fits (a) full unmasked range, (b) 30-60% only
    B_res = {ft: defaultdict(list) for ft in ("full", "window")}
    for (mod, col), (phi, C, SD) in systems.items():
        LK = np.array([logk(C[i], D) for i in range(len(phi))])  # (n_phi, 94)
        for ft, mask in (("full", np.ones(len(phi), bool)),
                         ("window", np.array([inwin(p) for p in phi]))):
            if mask.sum() < 2:
                continue
            X = np.vstack([np.ones(mask.sum()), -phi[mask]]).T
            beta, *_ = np.linalg.lstsq(X, LK[mask], rcond=None)  # (2, 94)
            pred = np.vstack([np.ones(len(phi)), -phi]).T @ beta
            resid = pred - LK
            for i, p in enumerate(phi):
                region = "inside" if inwin(p) else "outside"
                B_res[ft][(mod, region)].extend(resid[i].tolist())
                B_res[ft][(mod, "all")].extend(resid[i].tolist())

    print("\n== Route B: LSS residual |dlog10 k| ==")
    for ft in ("full", "window"):
        for key in sorted(B_res[ft]):
            n, med, p90 = stats(B_res[ft][key])
            print(f"  fit={ft:6s} {key[0]:15s} {key[1]:7s} n={n:6d} median={med:.4f} p90={p90:.4f}")

    # ---------------- adjacent-phi coefficient correlation (across columns)
    print("\n== Adjacent-phi coefficient correlation across columns (Pearson r) ==")
    for mod in ("methanol", "acetonitrile", "tetrahydrofuran"):
        line = [mod[:4]]
        for j, cname in enumerate(COEFS):
            pairs = defaultdict(list)  # (phi_lo, phi_hi) -> list of (v_lo, v_hi)
            for (m, col), (phi, C, SD) in systems.items():
                if m != mod: continue
                for i in range(len(phi) - 1):
                    pairs[(phi[i], phi[i + 1])].append((C[i, j], C[i + 1, j]))
            rs = []
            for k, v in pairs.items():
                if len(v) >= 5:
                    a = np.array(v)
                    rs.append(np.corrcoef(a[:, 0], a[:, 1])[0, 1])
            if not rs:
                line.append(f"{cname}: n/a (<5 columns)")
                continue
            line.append(f"{cname}: r_med={np.median(rs):.3f} (min {min(rs):.3f}, {len(rs)} steps)")
        print("  " + " | ".join(line))

    # ---------------- curvature: mid-point deviation from chord, in own-sd units
    print("\n== Coefficient curvature: |mid-point - chord| / sd_coef (unmasked span) ==")
    curv = defaultdict(list)  # (mod, coef) -> list of |dev|/sd
    for (mod, col), (phi, C, SD) in systems.items():
        if len(phi) < 3: continue
        for j in range(6):
            chord = C[0, j] + (C[-1, j] - C[0, j]) * (phi - phi[0]) / (phi[-1] - phi[0])
            dev = (C[:, j] - chord)[1:-1] / SD[1:-1, j]
            curv[(mod, COEFS[j])].extend(np.abs(dev).tolist())
    for mod in ("methanol", "acetonitrile", "tetrahydrofuran"):
        line = [mod[:4]]
        for cname in COEFS:
            a = np.array(curv[(mod, cname)])
            line.append(f"{cname}: med={np.median(a):.1f} p90={np.percentile(a,90):.1f} n={len(a)} >2sd:{(a>2).mean()*100:.0f}%")
        print("  " + " | ".join(line))

    # ---------------- head-to-head at shared evaluation points
    # Route A LOO error vs Route B (full-fit) residual at the same interior phi points
    print("\n== Head-to-head at interior phi points (same evaluation set) ==")
    for mod in ("methanol", "acetonitrile", "tetrahydrofuran"):
        for region in ("inside", "outside"):
            a_e, b_e = [], []
            for (m, col), (phi, C, SD) in systems.items():
                if m != mod or len(phi) < 3: continue
                LK = np.array([logk(C[i], D) for i in range(len(phi))])
                X = np.vstack([np.ones(len(phi)), -phi]).T
                beta, *_ = np.linalg.lstsq(X, LK, rcond=None)
                pred = X @ beta
                for i in range(1, len(phi) - 1):
                    if ("inside" if inwin(phi[i]) else "outside") != region: continue
                    keep = np.arange(len(phi)) != i
                    chat = np.array([interp(phi[keep], C[keep, j], phi[i], "pchip") for j in range(6)])
                    a_e.extend((logk(chat, D) - LK[i]).tolist())
                    b_e.extend((pred[i] - LK[i]).tolist())
            if a_e:
                na, ma, pa = stats(a_e); nb, mb, pb = stats(b_e)
                print(f"  {mod:15s} {region:7s} n={na:6d}  A(pchip LOO) med={ma:.4f} p90={pa:.4f} | B(LSS full) med={mb:.4f} p90={pb:.4f}")

    print(f"\nQSPR per-compound budget for comparison: {QSPR_BUDGET} log10 k")

if __name__ == "__main__":
    main()
