#!/usr/bin/env python3
"""#30's column decision: power AND generalisability, weighted equally.

A beautifully powered experiment on the wrong part of the silanol axis is still
the wrong experiment, so this answers two questions and lets either veto:

  A. POWER — can the sweep resolve `D`'s per-run part at all? Monte-Carlo the
     protocol's own three-parameter fit under realistic retention precision and
     a per-column t0 systematic, then propagate to the paired
     column-A-minus-column-B difference over the tier-1 acids.

  B. GENERALISABILITY — where does each candidate pair sit on the silanol axis
     the HSM `C(7.0)` term measures, and what does a measurement there license
     us to say about columns elsewhere on that axis?

Stdlib + numpy only (no scipy in this environment). Deterministic seed.
"""
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HSM = os.path.join(HERE, "..", "..", "sources", "hsm-column-db", "database.csv")
RNG = np.random.default_rng(20260816)

# ---------------------------------------------------------------- model ----
# Protocol analysis plan:
#   k_obs(pH) = k_neutral * [ f_neutral + (1 - f_neutral) * 10^(-D) ]
# Free: k_neutral, D, apparent pKa. Acid, so f_neutral = 1/(1+10^(pH-pKa)).
PH = np.array([2.5, 3.5, 5.0, 6.5, 8.0])   # the 5 shared pH levels
REPLICATES = 2                              # "duplicate injections at a minimum"


def log10k(theta, ph):
    lkn, D, pka = theta
    f = 1.0 / (1.0 + 10.0 ** (ph - pka))
    return lkn + np.log10(f + (1.0 - f) * 10.0 ** (-D))


def fit(ph, y, theta0):
    """Gauss-Newton on log10 k residuals. 3 params, well conditioned here."""
    th = np.array(theta0, float)
    for _ in range(60):
        r = y - log10k(th, ph)
        J = np.empty((len(ph), 3))
        for i in range(3):
            d = np.zeros(3); d[i] = 1e-6
            J[:, i] = (log10k(th + d, ph) - log10k(th - d, ph)) / 2e-6
        try:
            step = np.linalg.solve(J.T @ J + 1e-9 * np.eye(3), J.T @ r)
        except np.linalg.LinAlgError:
            return None
        th += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return th


def simulate_D(k_neutral, D_true, pka, rt_rsd, t0_rel_err, n_mc):
    """Return sampled D estimates for one compound on one column.

    t0_rel_err is a PER-COLUMN systematic: one draw per simulated campaign,
    shared by every compound on that column, so it never averages down.
    """
    ph = np.repeat(PH, REPLICATES)
    out = []
    for _ in range(n_mc):
        f = 1.0 / (1.0 + 10.0 ** (ph - pka))
        k_true = k_neutral * (f + (1.0 - f) * 10.0 ** (-D_true))
        t0 = 1.0
        tR = t0 * (1.0 + k_true)
        tR = tR * (1.0 + RNG.normal(0, rt_rsd, tR.shape))     # injection noise
        t0_obs = t0 * (1.0 + RNG.normal(0, t0_rel_err))       # per-campaign t0
        k_obs = np.maximum((tR - t0_obs) / t0_obs, 1e-6)
        th = fit(ph, np.log10(k_obs), [np.log10(k_neutral), D_true, pka])
        if th is not None and 0 < th[1] < 6:
            out.append(th[1])
    return np.array(out)


def power_section():
    print("=" * 72)
    print("A. POWER — can the sweep resolve D's per-run part?")
    print("=" * 72)
    D_TRUE, PKA, N_MC = 1.5, 4.2, 400
    RT_RSD, T0_ERR = 0.003, 0.01      # 0.3% injection RSD, 1% t0 (marker in every run)

    print(f"\nAssumptions: D_true={D_TRUE}, pKa={PKA}, {len(PH)} pH x {REPLICATES} reps,")
    print(f"retention RSD {RT_RSD:.1%}, per-column t0 uncertainty {T0_ERR:.0%}\n")
    print("Per-compound SD of the fitted D, by how retained the NEUTRAL form is:")
    print(f"  {'k_neutral':>10}{'k_ion':>9}{'tR_ion/t0':>11}{'SD(D)':>9}")
    per_k = {}
    for kn in (2, 5, 10, 20, 40):
        s = simulate_D(kn, D_TRUE, PKA, RT_RSD, T0_ERR, N_MC)
        k_ion = kn * 10 ** (-D_TRUE)
        per_k[kn] = s.std()
        print(f"  {kn:>10}{k_ion:>9.3f}{1 + k_ion:>11.3f}{s.std():>9.4f}")
    print("\n  -> D's precision is governed by whether the IONISED form is still")
    print("     retained. At k_neutral=2 the ionised peak sits at 1.06*t0 and t0")
    print("     error swamps it; by k_neutral=20 it is comfortable.")

    # Paired difference: same compound both columns, so the per-compound D_c
    # cancels exactly. Only measurement error and the two t0 systematics remain.
    print("\nPaired difference (column A minus column B), tier-1 acids:")
    print(f"  {'n_acids':>8}{'SD(mean dD)':>14}{'min detectable (95%)':>23}")
    kn_mix = [20, 25, 15, 30, 10, 20, 35, 12]     # a realistic tier-1 spread
    for n_ac in (4, 8, 12):
        diffs = []
        for _ in range(300):
            per = []
            for kn in kn_mix[:n_ac]:
                a = simulate_D(kn, D_TRUE, PKA, RT_RSD, T0_ERR, 1)
                b = simulate_D(kn, D_TRUE, PKA, RT_RSD, T0_ERR, 1)
                if len(a) and len(b):
                    per.append(a[0] - b[0])
            if per:
                diffs.append(np.mean(per))
        sd = np.std(diffs)
        print(f"  {n_ac:>8}{sd:>14.4f}{1.96 * sd:>23.4f}")

    print(f"\n  Prior on the per-run term (slice.py): D_run ~ N(0, 0.30)")
    print("  -> compare the minimum detectable difference against that 0.30.")
    return per_k


def generalisability_section():
    print("\n" + "=" * 72)
    print("B. GENERALISABILITY — is it the right part of the silanol axis?")
    print("=" * 72)
    rows = list(csv.DictReader(open(HSM, encoding="latin-1")))

    def f(s):
        try:
            return float(s)
        except (TypeError, ValueError):
            return None

    pop = np.array(sorted(f(r["C70"]) for r in rows
                          if r["USPtype"] == "L1" and r["phase"].strip() == "C18"
                          and f(r["C70"]) is not None))
    typeA = np.array([f(r["C70"]) for r in rows
                      if r["type"] == "A" and f(r["C70"]) is not None])
    print(f"\nReference population: {len(pop)} silica C18 columns (USPtype L1, phase C18)")
    for q in (5, 25, 50, 75, 95, 99):
        print(f"  p{q:<3} C70 = {np.percentile(pop, q):>6.2f}")
    print(f"  type-A columns (all phases): n={len(typeA)}, median C70 = "
          f"{np.median(typeA):.3f}")

    def pct(v):
        return 100.0 * np.mean(pop <= v)

    # The ASYMMETRIC design matters because the protocol's tier 1 fits k_neutral
    # from the data ("D measured on both columns, no model input"). WSU
    # membership buys nothing for tier 1 — only for tier 2b (strong bases, no
    # reachable neutral plateau) and tier 3 (neutral LSER controls). So one WSU
    # column plus one high-C non-WSU column loses nothing on the silanol question.
    pairs = [
        ("BOOKED (superseded)", -0.05, "Betasil C18",     -0.04, "both WSU"),
        ("WSU-internal",        -0.09, "Kinetex XB-C18",   0.30, "both WSU"),
        ("ASYMMETRIC",          -0.09, "Apex II C18",      2.69, "SunFire WSU + non-WSU"),
        ("Outside WSU",         -0.17, "Apex II C18",      2.69, "neither WSU"),
    ]
    print(f"\n{'design':<22}{'span C70':>17}{'dC':>7}{'pctile span':>14}"
          f"{'% pop':>8}  baseline")
    res = {}
    for label, va, nb, vb, note in pairs:
        lo, hi = min(va, vb), max(va, vb)
        cov = 100.0 * np.mean((pop >= lo) & (pop <= hi))
        res[label] = (lo, hi, cov)
        print(f"{label:<22}{f'{lo:.2f} to {hi:.2f}':>17}{hi - lo:>7.2f}"
              f"{f'p{pct(lo):.0f}-p{pct(hi):.0f}':>14}{cov:>7.1f}%  {note}")

    print("\nExtrapolation demanded of each design.")
    print("A two-column sweep fits D_run at two points on the C70 axis. Predicting")
    print("at another column is interpolation inside the span, extrapolation outside.")
    print(f"\n  {'design':<22}{'to type-A median 1.05':>26}{'to p99 (2.19)':>22}")
    for label, (lo, hi, cov) in res.items():
        span = hi - lo
        tgts = (np.median(typeA), np.percentile(pop, 99))
        cells = []
        for t in tgts:
            cells.append("inside span" if t <= hi
                         else f"{(t - hi) / span:.1f}x span beyond")
        print(f"  {label:<22}{cells[0]:>26}{cells[1]:>22}")

    print("\n  -> 'x span beyond' is how far past the calibrated range you must reach")
    print("     to say anything about a high-silanol column. Two points fit a line;")
    print("     they cannot test whether the relationship IS a line.")


if __name__ == "__main__":
    power_section()
    generalisability_section()
