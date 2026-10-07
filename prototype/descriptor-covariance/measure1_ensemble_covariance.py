#!/usr/bin/env python3
"""Measurement 1 for issue #37: does SoluteML's within-compound (E,S,A,B) ensemble
covariance change Var(log k) vs the independent-diagonal assumption?

Input: soluteml_preds_raw.csv (25 members x 94 compounds x E,S,A,B,L).
For each compound: 4x4 ensemble covariance C over (E,S,A,B) (ddof=1, n=25).
Propagate through LSER coefficient vector w=(e,s,a,b) of a representative fit:
Var_full = w'Cw vs Var_diag = sum w_j^2 C_jj. V is not in the propagation --
SoluteML does not predict V (McGowan V is computed from structure, error-free).

Also, since Table S-2 is the reference lineage, the ensemble-mean *residuals*
vs S-2 are computable per compound -- giving the across-compound piece
(measurement 2's QSPR half) for scaffold families within the 94.

Caveats printed with results: ensemble spread is epistemic-only and the WSU 94 are
inside or near SoluteML's training distribution, so magnitudes are lower bounds;
5 folds x 5 models share data, so the 25 members are not independent draws.
"""

import csv
import itertools
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "sources", "wsu-lser")
DESC = ["E", "S", "A", "B"]

# scaffold families within the 94 (name-based, conservative)
FAMILIES = {
    "n-alkylphenones": ["Acetophenone", "Propriophenone", "Butyrophenone",
                         "Valerophenone", "Hexanophenone", "Octanophenone"],
    "phenols": ["Phenol", "2-Methylphenol", "4-Methylphenol", "2,6-Dimethylphenol",
                 "4-Chlorophenol", "3-Bromophenol", "4-Cyanophenol",
                 "2-Nitrophenol", "3-Nitrophenol", "4-Nitrophenol"],
    "anilines": ["Aniline", "2-Methylaniline", "3-Methylaniline", "4-Methylaniline",
                  "4-Chloroaniline", "3,4-Dichloroaniline", "4-Fluoroaniline",
                  "N-Ethylaniline", "N, N-Dimethylaniline", "N, N-Diethylaniline"],
    "alkylbenzenes": ["Toluene", "Ethylbenzene", "Propylbenzene", "m-Xylene",
                       "p-Xylene", "1,2,3-Trimethylbenzene"],
    "naphthalenes": ["Naphthalene", "1-Chloronaphthalene", "1-Bromonaphthalene",
                      "2-Methoxynaphthalene", "1-Naphthol", "2-Naphthol",
                      "2-Naphthaldehyde"],
}


def summarize(x, label):
    x = np.asarray(x, dtype=float)
    q = np.percentile(x, [5, 25, 50, 75, 95])
    print(f"{label}: n={x.size}  median={q[2]:.3f}  IQR=[{q[1]:.3f},{q[3]:.3f}]  "
          f"5-95%=[{q[0]:.3f},{q[4]:.3f}]")


def main():
    preds = defaultdict(list)  # compound -> list of (E,S,A,B)
    with open(os.path.join(HERE, "soluteml_preds_raw.csv")) as f:
        for rec in csv.DictReader(f):
            preds[rec["compound"]].append([float(rec[k]) for k in DESC])
    names = list(preds)
    P = np.array([preds[n] for n in names])  # (94, 25, 4)
    assert P.shape == (94, 25, 4), P.shape

    ref = {}
    with open(os.path.join(SRC, "wsu2019-descriptors-s2.csv")) as f:
        for rec in csv.DictReader(f):
            # names are clean since the 2026-08-16 extract fix (`class` column
            # now carries the section heading that used to glue onto a name)
            ref[rec["compound"]] = np.array([float(rec[k]) for k in DESC])
    # smiles file uses S-2 names already normalised
    R = np.array([ref[n] for n in names])  # (94, 4)

    means = P.mean(axis=1)
    resid = means - R
    print("Sanity: SoluteML ensemble mean vs Table S-2 (n=94):")
    for j, d in enumerate(DESC):
        print(f"  {d}: bias={resid[:, j].mean():+.3f}  RMSE={np.sqrt((resid[:, j] ** 2).mean()):.3f}  "
              f"ensemble sd median={np.median(P[:, :, j].std(axis=1, ddof=1)):.3f}")

    covs = np.array([np.cov(P[i].T, ddof=1) for i in range(len(names))])  # (94,4,4)
    # mean within-compound error correlation matrix
    corr = np.zeros((4, 4))
    for i in range(len(names)):
        d = np.sqrt(np.diag(covs[i]))
        with np.errstate(invalid="ignore", divide="ignore"):
            corr += covs[i] / np.outer(d, d)
    corr /= len(names)
    print("\nMean within-compound ensemble correlation over (E,S,A,B):")
    print("      " + "  ".join(f"{d:>6s}" for d in DESC))
    for i, d in enumerate(DESC):
        print(f"  {d:>3s} " + "  ".join(f"{corr[i, j]:6.3f}" for j in range(4)))

    # LSER coefficient vectors for representative fits
    fits = {}
    with open(os.path.join(SRC, "wsu2019-system-constants.csv")) as f:
        for rec in csv.DictReader(f):
            if rec["column"] == "Ascentis C18" and rec["phi_pct_vv"] == "50":
                fits[rec["modifier"]] = np.array(
                    [float(rec[k]) for k in ("e", "s", "a", "b")]
                )

    for mod, w in sorted(fits.items()):
        vf = np.einsum("j,ijk,k->i", w, covs, w)
        vd = covs[:, range(4), range(4)] @ (w**2)
        print(f"\nVar(log k) from descriptor ensemble covariance, Ascentis C18 {mod} 50% "
              f"(w = e,s,a,b = {np.round(w, 3)}):")
        summarize(vf / vd, "  full/diag ratio")
        summarize(np.sqrt(vf), "  sd(log k) full")
        summarize(np.sqrt(vd), "  sd(log k) diag")

        # across-compound: residual log k error and its intra-family correlation
        e = resid @ w  # per-compound log k error at this condition
        print(f"  per-compound log k residual vs S-2: sd={e.std(ddof=1):.3f} (n=94)")
        idx = {n: i for i, n in enumerate(names)}
        prods, var_terms = [], []
        for fam, members in FAMILIES.items():
            ii = [idx[m] for m in members]
            ee = e[ii]
            pair_mean = np.mean([ee[a] * ee[b] for a, b in
                                 itertools.combinations(range(len(ii)), 2)])
            prods.append(pair_mean)
            var_terms.append(ee.var(ddof=1))
            print(f"    {fam} (n={len(ii)}): mean pairwise product={pair_mean:+.4f}  "
                  f"within-family var={ee.var(ddof=1):.4f}")
        rho = np.mean(prods) / (e.var(ddof=1))
        print(f"  implied intra-family correlation vs overall var ({e.var(ddof=1):.4f}): "
              f"rho~{rho:+.2f}")


if __name__ == "__main__":
    main()
