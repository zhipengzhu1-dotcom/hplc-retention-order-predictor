#!/usr/bin/env python3
"""Thin slice (#21): SMILES -> SoluteML descriptors -> ionisation -> LSER ->
k(phi, pH) -> ranked retention distribution, on ONE column.

Column: XBridge Shield RP18, acetonitrile (WSU-2019). Chosen because it is the
only shortlisted column with NO incomplete-wetting row in ACN (support 10-70%
fully unmasked) and only a "weak electrostatic" Table S-1 caveat.

All retention math in log10 k (WSU convention, CONTEXT.md). 45 degC dataset.
Isocratic operating points only; dwell/gradient physics is OUT of slice scope.
t_R = t0 * (1 + k) with nominal t0 = 1.0 min. width = 4 t_R / sqrt(N),
N = 10000 nominal -- an interface placeholder, not a plate model.

Sampler (three-tier correlation structure, #37/#14):
  PER-RUN (one draw shared by all compounds):
    - joint 6-coefficient column-vector error: Sigma = D R D per phi grid point
      (R from (X'X)^-1 of the 94-compound Table S-2 design; diagonal draws are
      FORBIDDEN, #37: 16x variance overstatement). One shared z6 across phi
      (same calibration compounds -> assume perfect phi-correlation; stated).
    - pH-scale systematic (NAMED GAP #20/#22: Roses/Bosch w_w->s_s correction
      unacquired): acid pKa shift ~ N(+0.3, 0.30), base pKa shift ~
      N(-0.15, 0.20). Literature-typical magnitudes for 30-60% ACN, not fitted.
    - D_run ~ N(0, 0.30): column/modifier part of the ionisation drop.
  PER-COMPOUND (independent draws; scaffold correlation NOT modelled -- the
  set is chemically diverse; see FINDINGS):
    - descriptor error ~ N(0, Sigma_c): Sigma_c is the 25-member ensemble 4x4
      covariance INFLATED per target so its marginal sd matches the WSU-94
      RMSE (ensemble spread is uncalibrated dispersion, not calibrated
      confidence -- the inflation ties magnitude to in-sample RMSE; labelled).
    - pKa error ~ N(0, 0.35) per site (literature values; source spread for
      ketoprofen alone is 0.47).
    - D_c ~ N(0, 0.40).
  D_total = clip(1.5 + D_run + D_c, 0.1, inf): truncated positive per CONTEXT
  (never zero, never infinity).

Refused compound (no neutral microspecies at any pH) keeps its column with NaN
retention (#33): n never shrinks.

Run with the soluteml venv python (needs numpy, scipy, rdkit).
"""

import csv
import itertools
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from compounds import A_NONZERO_TRUTH, COMPOUNDS  # noqa: E402

SRC = os.path.join(HERE, "..", "..", "sources", "wsu-lser")
DCOV = os.path.join(HERE, "..", "descriptor-covariance")

COLUMN, MODIFIER = "XBridge Shield RP18", "acetonitrile"
COEFS = ["c", "e", "s", "a", "b", "v"]
S = 1000
T0 = 1.0  # min, nominal
OPS = [(0.30, 3.0), (0.30, 7.0), (0.50, 5.0)]  # (phi, pH)
SEED = 20260815

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


def crippen_logp(smiles):
    from rdkit import Chem
    from rdkit.Chem import Crippen
    return Crippen.MolLogP(Chem.MolFromSmiles(smiles))


# ---------------------------------------------------------------- data loading

def load_column():
    """phi grid and (beta, sd) per grid point, unmasked support only (#24)."""
    grid, betas, sds = [], [], []
    with open(os.path.join(SRC, "wsu2019-system-constants.csv")) as f:
        for rec in csv.DictReader(f):
            if rec["column"] != COLUMN or rec["modifier"] != MODIFIER:
                continue
            if rec["incomplete_wetting"] == "1":
                continue  # masked: interpolation must not cross or use it
            grid.append(float(rec["phi_pct_vv"]) / 100.0)
            betas.append([float(rec[k]) for k in COEFS])
            sds.append([float(rec["sd_" + k]) for k in COEFS])
    order = np.argsort(grid)
    grid = np.array(grid)[order]
    # #24 invariant: support must be contiguous on the published 10% grid --
    # PCHIP inside it never creates support or crosses a masked point.
    assert np.allclose(np.diff(grid), 0.10), "support not contiguous; slice assumes it"
    return grid, np.array(betas)[order], np.array(sds)[order]


def design_R():
    """Correlation of (X'X)^-1 from the 94-compound Table S-2 design (#37)."""
    rows = []
    with open(os.path.join(SRC, "wsu2019-descriptors-s2.csv")) as f:
        for rec in csv.DictReader(f):
            rows.append([1.0] + [float(rec[k]) for k in ("E", "S", "A", "B", "V")])
    X = np.array(rows)
    M = np.linalg.inv(X.T @ X)
    d = np.sqrt(np.diag(M))
    return M / np.outer(d, d)


def load_ensemble_stats():
    """Per-compound ensemble mean (E,S,A,B) and 4x4 covariance, plus the
    per-target inflation factor lambda = RMSE(WSU94) / median ensemble sd."""
    raw = defaultdict(list)
    for path in (os.path.join(DCOV, "soluteml_preds_raw.csv"),
                 os.path.join(HERE, "drug_preds_raw.csv")):
        with open(path) as f:
            for rec in csv.DictReader(f):
                raw[rec["compound"]].append([float(rec[t]) for t in "ESAB"])
    stats = {}
    for name, v in raw.items():
        a = np.array(v)
        assert a.shape[0] == 25, (name, a.shape)
        stats[name] = (a.mean(axis=0), np.cov(a.T))

    # inflation from the full WSU-94 comparison (in-sample-adjacent; labelled)
    s2 = {}
    with open(os.path.join(SRC, "wsu2019-descriptors-s2.csv")) as f:
        for rec in csv.DictReader(f):
            # names are clean since the 2026-08-16 extract fix (`class` column
            # now carries the section heading that used to glue onto a name)
            s2[rec["compound"]] = np.array([float(rec[t]) for t in "ESAB"])
    errs, sds = [], []
    for name, meas in s2.items():
        if name not in raw:
            continue
        a = np.array(raw[name])
        errs.append(a.mean(axis=0) - meas)
        sds.append(a.std(axis=0, ddof=1))
    errs, sds = np.array(errs), np.array(sds)
    rmse = np.sqrt((errs**2).mean(axis=0))
    lam = rmse / np.median(sds, axis=0)
    return stats, s2, rmse, lam, len(errs)


# ------------------------------------------------------------------- sampler

def f_neutral(ctype, pkas, pH):
    """Neutral-microspecies fraction. pkas: array (..., nsites)."""
    if ctype == "neutral":
        return np.ones(np.shape(pkas)[:-1]) if np.ndim(pkas) > 1 else 1.0
    if ctype == "acid":
        return 1.0 / (1.0 + 10.0 ** (pH - pkas[..., 0]))
    if ctype == "base":
        return 1.0 / (1.0 + 10.0 ** (pkas[..., 0] - pH))
    if ctype == "ampho":  # cation <-(pKa1 basic site)-> neutral <-(pKa2 acid)-> anion
        return 1.0 / (1.0 + 10.0 ** (pkas[..., 0] - pH) + 10.0 ** (pH - pkas[..., 1]))
    raise ValueError(ctype)


def run_ensemble(freeze=frozenset(), drop_eE=False, seed=SEED, s=S):
    """Return dict with retention[s, n], width[s, n], logk[s, n] per op point.

    freeze: subset of {"desc", "column", "pKa", "D"} -- zero that layer's draws
    (frozen at ensemble mean / published vector / literature pKa / D = 1.5).
    """
    rng = np.random.default_rng(seed)
    grid, betas, sds = load_column()
    R = design_R()
    stats, s2, _, lam, _ = load_ensemble_stats()

    names = [c["name"] for c in COMPOUNDS]
    n = len(names)
    V = {c["name"]: mcgowan_v(c["smiles"]) for c in COMPOUNDS}

    # --- per-run draws (shared across compounds)
    z6 = rng.standard_normal((s, 6))
    if "column" in freeze:
        z6[:] = 0.0
    # cholesky of Sigma = D R D at each grid point; shared z6 across phi
    dbeta = np.empty((s, len(grid), 6))
    for g in range(len(grid)):
        L = np.linalg.cholesky(np.outer(sds[g], sds[g]) * R
                               + 1e-12 * np.eye(6))
        dbeta[:, g, :] = z6 @ L.T
    acid_shift = rng.normal(+0.30, 0.30, s)
    base_shift = rng.normal(-0.15, 0.20, s)
    D_run = rng.normal(0.0, 0.30, s)
    if "pKa" in freeze:
        acid_shift[:] = 0.0
        base_shift[:] = 0.0
    if "D" in freeze:
        D_run[:] = 0.0

    # --- per-compound draws
    desc = {}   # name -> (s, 4) sampled E,S,A,B
    pkas = {}   # name -> (s, nsites)
    D = {}      # name -> (s,)
    for c in COMPOUNDS:
        nm = c["name"]
        if c["type"] == "refused":
            continue
        mu, C = stats[nm]
        Cinf = np.outer(lam, lam) * C
        if "desc" in freeze:
            d4 = np.tile(mu, (s, 1))
        else:
            d4 = rng.multivariate_normal(mu, Cinf + 1e-10 * np.eye(4), size=s)
        d4[:, 2] = np.clip(d4[:, 2], 0.0, None)  # A >= 0 physically
        desc[nm] = d4
        pk = np.tile(np.array(c["pka"], dtype=float), (s, 1)) if c["pka"] else np.zeros((s, 0))
        if pk.shape[1] and "pKa" not in freeze:
            pk = pk + rng.normal(0.0, 0.35, pk.shape)  # per-compound literature error
            if c["type"] == "acid":
                pk[:, 0] += acid_shift
            elif c["type"] == "base":
                pk[:, 0] += base_shift
            elif c["type"] == "ampho":
                pk[:, 0] += base_shift   # site 1: anilinium (basic)
                pk[:, 1] += acid_shift   # site 2: carboxyl (acidic)
        pkas[nm] = pk
        D_c = 0.0 if "D" in freeze else rng.normal(0.0, 0.40, s)
        D[nm] = np.clip(1.5 + D_run + D_c, 0.1, None)  # never zero, never inf

    # --- evaluate at operating points
    out = {}
    for phi, pH in OPS:
        assert grid[0] <= phi <= grid[-1], "outside unmasked support (#24)"
        # per-scenario coefficient vector at phi via PCHIP over the grid
        beta_phi = np.empty((s, 6))
        for j in range(6):
            y = betas[:, j][:, None] + dbeta[:, :, j].T  # (grid, s)
            beta_phi[:, j] = PchipInterpolator(grid, y, axis=0)(phi)
        if drop_eE:
            beta_phi[:, 1] = 0.0

        logk = np.full((s, n), np.nan)
        for i, c in enumerate(COMPOUNDS):
            nm = c["name"]
            if c["type"] == "refused":
                continue  # NaN column: expected but unplaceable (#33)
            d4 = desc[nm]
            Xs = np.column_stack([np.ones(s), d4[:, 0], d4[:, 1], d4[:, 2],
                                  d4[:, 3], np.full(s, V[nm])])
            logk_n = np.einsum("sj,sj->s", Xs, beta_phi)
            f = f_neutral(c["type"], pkas[nm], pH) if c["pka"] else np.ones(s)
            k = 10.0 ** logk_n * (f + (1.0 - f) * 10.0 ** (-D[nm]))
            logk[:, i] = np.log10(k)
        tR = T0 * (1.0 + 10.0 ** logk)
        width = 4.0 * tR / np.sqrt(10000.0)  # nominal N; interface placeholder
        out[(phi, pH)] = dict(logk=logk, retention=tR, width=width)
    return names, out


# ------------------------------------------------------------------ analyses

def order_matrix(logk):
    """P[i, j] = P(i elutes before j) over scenarios (valid compounds only)."""
    s, n = logk.shape
    P = np.full((n, n), np.nan)
    for i, j in itertools.combinations(range(n), 2):
        a, b = logk[:, i], logk[:, j]
        if np.isnan(a).all() or np.isnan(b).all():
            continue
        p = np.mean(a < b)
        P[i, j], P[j, i] = p, 1.0 - p
    return P


def classify_pairs(P, valid_idx):
    conf, ambi, coin = [], [], []
    for i, j in itertools.combinations(valid_idx, 2):
        p = max(P[i, j], P[j, i])
        (conf if p > 0.9 else ambi if p >= 0.6 else coin).append((i, j, P[i, j]))
    return conf, ambi, coin


def main():
    names, base = run_ensemble()
    valid_idx = [i for i, c in enumerate(COMPOUNDS) if c["type"] != "refused"]
    logp = {c["name"]: crippen_logp(c["smiles"]) for c in COMPOUNDS}
    stats, s2, rmse, lam, n_anchor_pool = load_ensemble_stats()

    results = dict(column=COLUMN, modifier=MODIFIER, S=S, seed=SEED,
                   n_compounds=len(names), n_valid=len(valid_idx),
                   ops=[list(op) for op in OPS],
                   inflation_lambda=dict(zip("ESAB", np.round(lam, 2).tolist())),
                   wsu94_rmse=dict(zip("ESAB", np.round(rmse, 3).tolist())))

    # ---- 1. identity-prior usability: pair classification + logP baseline
    print("=" * 78)
    print("1. PAIRWISE ELUTION-ORDER CONFIDENCE (S=%d scenarios)" % S)
    results["order"] = {}
    for op in OPS:
        logk = base[op]["logk"]
        P = order_matrix(logk)
        conf, ambi, coin = classify_pairs(P, valid_idx)
        npairs = len(conf) + len(ambi) + len(coin)
        # logP baseline agreement (median model order vs Crippen order)
        med = np.nanmedian(logk, axis=0)
        agree = {}
        for label, group in (("confident", conf), ("ambiguous", ambi), ("coinflip", coin)):
            g = [(i, j) for i, j, _ in group]
            if not g:
                agree[label] = None
                continue
            a = [np.sign(med[i] - med[j]) == np.sign(logp[names[i]] - logp[names[j]])
                 for i, j in g]
            agree[label] = round(float(np.mean(a)), 3)
        print(f"\n  (phi={op[0]:.2f}, pH={op[1]:.1f}): {npairs} pairs -> "
              f"confident(P>0.9) {len(conf)}, ambiguous(0.6-0.9) {len(ambi)}, "
              f"coin-flip(<0.6) {len(coin)}")
        print(f"    logP-baseline agreement with model median order: {agree}")
        results["order"][f"{op[0]}_{op[1]}"] = dict(
            confident=len(conf), ambiguous=len(ambi), coinflip=len(coin),
            logp_agreement=agree,
            coinflip_pairs=[[names[i], names[j], round(float(p), 3)] for i, j, p in coin])
        # per-compound rank distribution summary
        ranks = np.argsort(np.argsort(logk[:, valid_idx], axis=1), axis=1)
        rank_med = np.median(ranks, axis=0)
        rank_iqr = np.percentile(ranks, 75, axis=0) - np.percentile(ranks, 25, axis=0)
        results["order"][f"{op[0]}_{op[1]}"]["rank_iqr_median"] = float(np.median(rank_iqr))
        print(f"    rank IQR (positions): median {np.median(rank_iqr):.1f}, "
              f"max {rank_iqr.max():.0f} "
              f"({names[valid_idx[int(np.argmax(rank_iqr))]]})")

    # ---- 2. variance attribution
    print("\n" + "=" * 78)
    print("2. VARIANCE ATTRIBUTION: share = 1 - Var(logk | layer frozen)/Var(logk)")
    layers = ["desc", "column", "pKa", "D"]
    frozen = {lay: run_ensemble(freeze={lay})[1] for lay in layers}
    results["attribution"] = {}
    for op in OPS:
        vt = np.nanvar(base[op]["logk"], axis=0)
        # class of each compound AT THIS pH (nominal pKa)
        cls = []
        for c in COMPOUNDS:
            if c["type"] == "refused":
                cls.append("refused")
            elif not c["pka"]:
                cls.append("neutral")
            else:
                f = f_neutral(c["type"], np.array(c["pka"])[None, :], op[1])[0]
                cls.append("neutral" if f > 0.9 else "ionised" if f < 0.1 else "partial")
        tab = {}
        for lay in layers:
            share = 1.0 - np.nanvar(frozen[lay][op]["logk"], axis=0) / vt
            for cl in ("neutral", "partial", "ionised"):
                idx = [i for i in valid_idx if cls[i] == cl]
                if idx:
                    tab.setdefault(cl, {})[lay] = round(float(np.median(share[idx])), 3)
        for cl in tab:
            idx = [i for i in valid_idx if cls[i] == cl]
            tab[cl]["n"] = len(idx)
            tab[cl]["sd_logk_median"] = round(float(np.median(np.sqrt(vt[idx]))), 3)
        print(f"\n  (phi={op[0]:.2f}, pH={op[1]:.1f}):")
        for cl, row in tab.items():
            print(f"    {cl:8s} (n={row['n']}, sd(logk) med {row['sd_logk_median']}): "
                  + "  ".join(f"{l}={row[l]:+.2f}" for l in layers))
        results["attribution"][f"{op[0]}_{op[1]}"] = tab

    # ---- 3. e.E experiment
    print("\n" + "=" * 78)
    print("3. e.E TERM ON/OFF")
    _, noE = run_ensemble(drop_eE=True)
    results["eE"] = {}
    for op in OPS:
        d_point = np.nanmedian(base[op]["logk"], axis=0) - np.nanmedian(noE[op]["logk"], axis=0)
        P1, P0 = order_matrix(base[op]["logk"]), order_matrix(noE[op]["logk"])
        dP = np.nanmax(np.abs(P1 - P0))
        c1 = [len(x) for x in classify_pairs(P1, valid_idx)]
        c0 = [len(x) for x in classify_pairs(P0, valid_idx)]
        print(f"  (phi={op[0]:.2f}, pH={op[1]:.1f}): "
              f"median|dlogk|={np.nanmedian(np.abs(d_point)):.4f}, "
              f"max|dlogk|={np.nanmax(np.abs(d_point)):.4f}, max|dP(order)|={dP:.3f}, "
              f"classes with/without eE: {c1} vs {c0}")
        results["eE"][f"{op[0]}_{op[1]}"] = dict(
            median_abs_dlogk=round(float(np.nanmedian(np.abs(d_point))), 4),
            max_abs_dlogk=round(float(np.nanmax(np.abs(d_point))), 4),
            max_abs_dP=round(float(dP), 3), classes_with=c1, classes_without=c0)

    # ---- 4. A = 0 categorical check (#6)
    print("\n" + "=" * 78)
    print("4. A=0 CATEGORICAL CHECK (threshold 0.01)")
    results["A_check"] = []
    grid, betas, _ = load_column()
    a_coef = {phi: PchipInterpolator(grid, betas[:, 3])(phi) for phi, _ in OPS}
    for c in COMPOUNDS:
        if c["type"] == "refused":
            continue
        nm = c["name"]
        pred_A = stats[nm][0][2]
        called_nonzero = pred_A > 0.01
        truth = A_NONZERO_TRUTH[nm]
        ok = called_nonzero == truth
        row = dict(compound=nm, pred_A=round(float(pred_A), 3),
                   called="A!=0" if called_nonzero else "A=0",
                   truth="A!=0" if truth else "A=0", correct=bool(ok))
        if not ok:
            # retention consequence: a * (A_pred - A_true-ish 0) at each op phi
            row["dlogk_consequence"] = {
                f"phi={phi}": round(float(a_coef[phi] * pred_A), 4) for phi, _ in OPS}
        results["A_check"].append(row)
        flag = "" if ok else "   <-- MIS-CALL"
        print(f"  {nm:24s} pred A={pred_A:6.3f} called {row['called']:5s} "
              f"truth {row['truth']:5s}{flag}")
    n_ok = sum(r["correct"] for r in results["A_check"])
    print(f"  categorical accuracy: {n_ok}/{len(results['A_check'])}")

    # ---- 5. store the #14 interface object
    ret = np.stack([base[op]["retention"] for op in OPS])   # (3, S, n)
    wid = np.stack([base[op]["width"] for op in OPS])
    np.savez_compressed(os.path.join(HERE, "retention_ensemble.npz"),
                        retention=ret, width=wid, compounds=np.array(names),
                        ops=np.array(OPS), t0=T0)
    nan_cols = [names[i] for i in range(len(names))
                if np.isnan(base[OPS[0]]["logk"][:, i]).all()]
    print(f"\nInterface object: retention[{len(OPS)}x{S}x{len(names)}], "
          f"NaN columns (refused): {nan_cols}")
    results["refused"] = nan_cols

    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("results.json written")


if __name__ == "__main__":
    main()
