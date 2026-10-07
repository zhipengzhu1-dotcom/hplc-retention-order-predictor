#!/usr/bin/env python3
"""Measurement 3 (pivoted) for issue #37.

True per-compound LSER residuals are NOT derivable from the repo: the WSU extract
carries fitted coefficients, per-coefficient SDs and the *aggregate* fit SE, but no
raw log k values. What IS computable is the covariance structure the fits imply for
predicted log k:

  (a) Approximate the 6x6 coefficient covariance per fit as Sigma = D R D, where D is
      the published per-coefficient SD vector and R is the correlation matrix of
      (X'X)^{-1} built from the 94-compound Table S-2 design [1,E,S,A,B,V].
      (The true 6x6 covariance is unpublished; OLS gives cov(beta) = sigma^2 (X'X)^{-1},
      so R is exact up to the unknown per-fit compound subset, whose descriptor
      geometry we assume matches the full 94.)
  (b) Propagate to Var(log k) per compound: full Sigma vs the diagonal (independent
      coefficient) assumption. Report the ratio distribution.
  (c) Propagate to Var(delta log k) for close pairs (intercept cancels; e..v terms do
      not, because delta-descriptors != 0). Full vs diagonal ratio.
  (d) LSS lack-of-fit implied by the LSER surface itself: per (compound, column,
      modifier), fit log k linear in phi over the usable (wetted) range and report
      the RMS curvature residual -- the model-implied floor of the #36 form
      correction, not an observed residual.

Wetting exclusion per #40: rows with incomplete_wetting=1 dropped everywhere.
WSU-2019 is log10 throughout (see CONTEXT.md).
"""

import csv
import itertools
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "sources", "wsu-lser")

COEFS = ["c", "e", "s", "a", "b", "v"]


def load_descriptors():
    names, rows = [], []
    with open(os.path.join(SRC, "wsu2019-descriptors-s2.csv")) as f:
        for rec in csv.DictReader(f):
            # The Table S-2 section header used to be merged into this name;
            # it is the `class` column since the 2026-08-16 extract fix.
            names.append(rec["compound"])
            rows.append([1.0] + [float(rec[k]) for k in ("E", "S", "A", "B", "V")])
    return names, np.array(rows)  # X: (94, 6) with leading 1


def load_fits():
    fits = []
    with open(os.path.join(SRC, "wsu2019-system-constants.csv")) as f:
        for rec in csv.DictReader(f):
            if rec["incomplete_wetting"] == "1":
                continue
            fits.append(
                dict(
                    modifier=rec["modifier"],
                    column=rec["column"],
                    phi=float(rec["phi_pct_vv"]),
                    beta=np.array([float(rec[k]) for k in COEFS]),
                    sd=np.array([float(rec["sd_" + k]) for k in COEFS]),
                    SE=float(rec["SE"]),
                    n=int(rec["n"]),
                )
            )
    return fits


def summarize(x, label):
    x = np.asarray(x, dtype=float)
    q = np.percentile(x, [5, 25, 50, 75, 95])
    print(
        f"{label}: n={x.size}  median={q[2]:.3f}  IQR=[{q[1]:.3f},{q[3]:.3f}]  "
        f"5-95%=[{q[0]:.3f},{q[4]:.3f}]  min={x.min():.3f}  max={x.max():.3f}"
    )
    return q


def main():
    names, X = load_descriptors()
    fits = load_fits()
    print(f"compounds: {len(names)}   fits kept (wetting-excluded): {len(fits)} of 354")

    # (a) correlation structure of OLS coefficient estimators from the design
    M = np.linalg.inv(X.T @ X)
    dM = np.sqrt(np.diag(M))
    R = M / np.outer(dM, dM)
    print("\nCoefficient-estimator correlation R from (X'X)^-1 of the 94-compound design:")
    print("      " + "  ".join(f"{k:>6s}" for k in COEFS))
    for i, k in enumerate(COEFS):
        print(f"  {k:>3s} " + "  ".join(f"{R[i, j]:6.3f}" for j in range(6)))

    # (b) Var(log k) per compound, full vs diagonal
    ratios = []
    for fit in fits:
        Sigma = np.outer(fit["sd"], fit["sd"]) * R
        var_full = np.einsum("ij,jk,ik->i", X, Sigma, X)
        var_diag = (X**2) @ (fit["sd"] ** 2)
        ratios.append(var_full / var_diag)
    ratios = np.concatenate(ratios)
    print(f"\n(b) Var(log k) full/diagonal ratio, {len(fits)} fits x {len(names)} compounds:")
    summarize(ratios, "    ratio")
    # absolute size of the coefficient-error contribution, diagonal, for context
    sd_diag = np.sqrt(
        np.concatenate([(X**2) @ (f["sd"] ** 2) for f in fits])
    )
    summarize(sd_diag, "    sd_diag(log k) absolute")

    # (c) close pairs at two representative conditions
    for mod, phi in (("methanol", 50.0), ("acetonitrile", 50.0)):
        fit = next(
            f for f in fits if f["column"] == "Ascentis C18" and f["modifier"] == mod and f["phi"] == phi
        )
        Sigma = np.outer(fit["sd"], fit["sd"]) * R
        logk = X @ fit["beta"]
        pair_ratios, dlogk_sds_full, dlogk_sds_diag = [], [], []
        n_close = 0
        for i, j in itertools.combinations(range(len(names)), 2):
            if abs(logk[i] - logk[j]) >= 0.1:
                continue
            n_close += 1
            dx = X[i] - X[j]  # first component 0: intercept error cancels
            vf = dx @ Sigma @ dx
            vd = (dx**2) @ (fit["sd"] ** 2)
            pair_ratios.append(vf / vd)
            dlogk_sds_full.append(np.sqrt(vf))
            dlogk_sds_diag.append(np.sqrt(vd))
        print(
            f"\n(c) close pairs (|pred dlogk|<0.1) on Ascentis C18 {mod} {phi:.0f}%: "
            f"{n_close} of {len(names) * (len(names) - 1) // 2} pairs"
        )
        summarize(pair_ratios, "    Var(dlogk) full/diag ratio")
        summarize(dlogk_sds_full, "    sd(dlogk) full")
        summarize(dlogk_sds_diag, "    sd(dlogk) diag")

    # (d) LSS lack-of-fit implied by the LSER surface (methanol + acetonitrile)
    for mod in ("methanol", "acetonitrile"):
        by_col = {}
        for f in fits:
            if f["modifier"] == mod:
                by_col.setdefault(f["column"], []).append(f)
        rms_all = []
        for col, fl in by_col.items():
            fl.sort(key=lambda f: f["phi"])
            if len(fl) < 4:
                continue
            phis = np.array([f["phi"] for f in fl]) / 100.0
            B = np.array([f["beta"] for f in fl])  # (nphi, 6)
            logk = X @ B.T  # (94, nphi)
            A = np.vstack([np.ones_like(phis), phis]).T
            coef, *_ = np.linalg.lstsq(A, logk.T, rcond=None)
            resid = logk.T - A @ coef  # (nphi, 94)
            rms_all.append(np.sqrt((resid**2).mean(axis=0)))  # per compound
        rms_all = np.concatenate(rms_all)
        print(f"\n(d) LSS curvature of LSER-predicted log k vs phi, {mod} "
              f"({len(by_col)} columns x {len(names)} compounds):")
        summarize(rms_all, "    RMS deviation from linear-in-phi (log k units)")


if __name__ == "__main__":
    main()
