#!/usr/bin/env python3
"""#35, reduced scope: the stock-SoluteML-vs-measured-WSU bias measurement.

No re-fit. #35's second bullet only: "measure the bias between stock predictions
and measured WSU values", which #22's closure left as the surviving deliverable
and which #31 made the project's ONLY instrument for detecting a systematic
descriptor offset (variance attribution is blind to bias by construction).

Inputs, all already in the repo:
  prototype/descriptor-covariance/soluteml_preds_raw.csv  25 members x 94 compounds
  sources/wsu-lser/wsu2019-descriptors-s2.csv             the measured anchor
  sources/wsu-lser/wsu2019-system-constants.csv           354 fits, 332 unmasked
  ./overlap.csv                                           SoluteDB membership + labels

Three things are measured, and they are different questions:

  1. MODEL BIAS       ensemble mean - WSU measured.        The headline.
  2. LABEL BIAS       SoluteDB label - WSU measured.       How far the training
                      corpus itself sits from independent measurement -- the
                      quantity #5 could only characterise by provenance.
  3. FIT ERROR        ensemble mean - SoluteDB label.      What the network adds
                      on top of its own labels, on compounds it has seen.

(1) = (2) + (3) identically, so the decomposition attributes the headline.

Reported through prototype/eval-harness's Figure: no float escapes without its
baseline, n, interval, axis, reference identity and reference kind.

Run: env/ionisation-venv/bin/python measure_bias.py   (numpy only; rdkit unused)
"""

import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "prototype", "eval-harness"))
from harness import Figure, bootstrap_over_compounds  # noqa: E402

D = ["E", "S", "A", "B"]

# Chung 2022 SI Fig S11/S12: SoluteML held-out MAE, scored against SoluteDB
# itself. A floor on true error, and the baseline our MAE is compared to.
PUBLISHED_MAE_RANDOM = {"E": 0.041, "S": 0.098, "A": 0.038, "B": 0.047}
PUBLISHED_MAE_SUBSTRUCT = {"E": 0.084, "S": 0.170, "A": 0.067, "B": 0.088}

# Poole 2019 Fig S-1/S-2: median LSER fit residual SD with the harmonised
# descriptor set. The scale on which a retention-space error is large or small.
POOLE_RESIDUAL_SD = {"acetonitrile": 0.030, "methanol": 0.041}

REFERENCE = ("WSU-2019 Table S-2 measured descriptors (94 neutral solutes, "
             "median V 1.06, no acids by name; 94/94 inside SoluteDB)")
DESC_AXIS = "descriptor space (condition-independent)"
RET_AXIS = "332 unmasked WSU (column x modifier x phi) fits, 10-90% v/v"


# ------------------------------------------------------------------- loading
def load():
    """Returns (meas, members, db, tier, fits, weak_base).

    The weak-base flag is the extract's `class` column, which carries Table
    S-2's own section heading. Until 2026-08-16 that heading was glued onto the
    first compound name of the block instead, and every name join here dropped
    that compound; fixed in sources/wsu-lser/extract_wsu.py.
    """
    meas, weak_base = {}, {}
    with open(os.path.join(REPO, "sources", "wsu-lser",
                           "wsu2019-descriptors-s2.csv")) as f:
        for r in csv.DictReader(f):
            name = r["compound"].strip()
            meas[name] = {k: float(r[k]) for k in ["E", "S", "A", "B", "V"]}
            weak_base[name] = r["class"] == "weak_base"

    members = defaultdict(list)
    with open(os.path.join(REPO, "prototype", "descriptor-covariance",
                           "soluteml_preds_raw.csv")) as f:
        for r in csv.DictReader(f):
            members[r["compound"].strip()].append({k: float(r[k]) for k in D})

    db = {}
    tier = {}
    with open(os.path.join(HERE, "overlap.csv")) as f:
        for r in csv.DictReader(f):
            tier[r["compound"].strip()] = r["match_tier"]
            if all(r["db_" + k] for k in D):
                db[r["compound"].strip()] = {k: float(r["db_" + k]) for k in D}

    fits = []
    with open(os.path.join(REPO, "sources", "wsu-lser",
                           "wsu2019-system-constants.csv")) as f:
        for r in csv.DictReader(f):
            if str(r.get("incomplete_wetting", "")).strip().lower() in ("true", "1", "yes"):
                continue
            fits.append({k: float(r[k]) for k in ["c", "e", "s", "a", "b", "v"]}
                        | {"mod": r["modifier"], "col": r["column"],
                           "phi": r["phi_pct_vv"]})
    return meas, members, db, tier, fits, weak_base


def logk(f, d, V):
    return (f["c"] + f["e"] * d["E"] + f["s"] * d["S"] + f["a"] * d["A"]
            + f["b"] * d["B"] + f["v"] * V)


def fig(value, interval, *, axis, baseline, baseline_value, n, units=""):
    lo, hi = interval
    return Figure(value=value, axis=axis, baseline=baseline,
                  baseline_value=baseline_value, n_compounds=n,
                  interval=(lo, hi), reference_identity=REFERENCE,
                  reference_kind="measurement", units=units)


def boot(n, stat):
    """(point, lo, hi) resampling compounds, per #15 s.1."""
    return bootstrap_over_compounds(n, stat)


def main() -> None:
    meas, members, db, tier, fits, weak_base = load()
    names = sorted(set(meas) & set(members))
    n = len(names)
    is_base = np.array([weak_base[c] for c in names])
    print(f"compounds matched: {n}   SoluteDB labels available: {len(db)}   "
          f"weak bases: {int(is_base.sum())}")
    print(f"membership tiers: "
          f"{ {t: sum(1 for c in names if tier[c] == t) for t in set(tier.values())} }")
    print(f"unmasked fits: {len(fits)}   "
          f"modifiers: { {f['mod'] for f in fits} }\n")

    # ---- per-compound arrays
    mu = {k: np.array([np.mean([m[k] for m in members[c]]) for c in names]) for k in D}
    esd = {k: np.array([np.std([m[k] for m in members[c]], ddof=1) for c in names])
           for k in D}
    truth = {k: np.array([meas[c][k] for c in names]) for k in D}
    V = np.array([meas[c]["V"] for c in names])
    lab = {k: np.array([db[c][k] for c in names]) for k in D}

    model_err = {k: mu[k] - truth[k] for k in D}      # 1. model bias
    label_err = {k: lab[k] - truth[k] for k in D}     # 2. label bias
    fit_err = {k: mu[k] - lab[k] for k in D}          # 3. fit error

    # ============================================================ 1. HEADLINE
    print("=" * 78)
    print("1. MODEL BIAS -- stock SoluteML ensemble mean vs WSU measured")
    print("=" * 78)
    for k in D:
        e = model_err[k]
        p, lo, hi = boot(n, lambda idx, e=e: float(np.mean(e[idx])))
        print(f"  bias {k}: " + fig(
            p, (lo, hi), axis=DESC_AXIS,
            baseline="zero (unbiased)", baseline_value=0.0, n=n).render())
        p, lo, hi = boot(n, lambda idx, e=e: float(np.mean(np.abs(e[idx]))))
        print(f"   MAE {k}: " + fig(
            p, (lo, hi), axis=DESC_AXIS,
            baseline="Chung 2022 SoluteML random-split MAE vs SoluteDB (self-reported)",
            baseline_value=PUBLISHED_MAE_RANDOM[k], n=n).render())
        p, lo, hi = boot(n, lambda idx, e=e: float(np.sqrt(np.mean(e[idx] ** 2))))
        print(f"  RMSE {k}: " + fig(
            p, (lo, hi), axis=DESC_AXIS,
            baseline="Chung 2022 SoluteML substructure-split MAE vs SoluteDB",
            baseline_value=PUBLISHED_MAE_SUBSTRUCT[k], n=n).render())
        print()

    # ====================================================== 2. RESIDUAL SHAPE
    print("=" * 78)
    print("2. RESIDUAL DISTRIBUTION (model - measured), n = %d compounds" % n)
    print("=" * 78)
    qs = [0, 5, 25, 50, 75, 95, 100]
    print(f"{'desc':>4} " + " ".join(f"{('p' + str(q)):>7}" for q in qs)
          + f" {'sd':>7} {'skew':>6} {'|e|>2sd':>8}")
    for k in D:
        e = model_err[k]
        v = np.percentile(e, qs)
        sd = float(np.std(e, ddof=1))
        skew = float(np.mean(((e - e.mean()) / sd) ** 3))
        out = int(np.sum(np.abs(e - e.mean()) > 2 * sd))
        print(f"{k:>4} " + " ".join(f"{x:+7.3f}" for x in v)
              + f" {sd:7.3f} {skew:+6.2f} {out:8d}")

    print("\n  worst 5 compounds per descriptor (signed error):")
    for k in D:
        order = np.argsort(-np.abs(model_err[k]))[:5]
        items = ", ".join(f"{names[i]} {model_err[k][i]:+.3f}" for i in order)
        print(f"    {k}: {items}")

    # proportional bias: does the error scale with the descriptor value?
    print("\n  proportional-bias check -- OLS residual ~ a + b*measured:")
    for k in D:
        x, e = truth[k], model_err[k]

        def slope(idx, x=x, e=e):
            xi, ei = x[idx], e[idx]
            if np.std(xi) < 1e-12:
                return float("nan")
            return float(np.polyfit(xi, ei, 1)[0])

        p, lo, hi = boot(n, slope)
        flag = "" if lo <= 0 <= hi else "   <-- interval excludes 0"
        print(f"    {k}: slope {p:+.3f} [{lo:+.3f}, {hi:+.3f}]{flag}")

    # ======================================================= 3. DECOMPOSITION
    print("\n" + "=" * 78)
    print("3. DECOMPOSITION  model-measured = (label-measured) + (model-label)")
    print("=" * 78)
    print(f"{'desc':>4} {'model bias':>11} {'label bias':>11} {'fit bias':>10} "
          f"{'RMS model':>10} {'RMS label':>10} {'RMS fit':>9} {'corr':>7} {'label share':>12}")
    for k in D:
        me, le, fe = model_err[k], label_err[k], fit_err[k]
        rms = lambda a: float(np.sqrt(np.mean(a ** 2)))
        corr = float(np.corrcoef(le, fe)[0, 1])
        share = float(np.var(le) / np.var(me)) if np.var(me) > 0 else float("nan")
        print(f"{k:>4} {np.mean(me):+11.4f} {np.mean(le):+11.4f} {np.mean(fe):+10.4f} "
              f"{rms(me):10.4f} {rms(le):10.4f} {rms(fe):9.4f} {corr:+7.2f} {share:12.2f}")

    for k in D:
        le = label_err[k]
        p, lo, hi = boot(n, lambda idx, e=le: float(np.mean(e[idx])))
        print(f"  label bias {k}: " + fig(
            p, (lo, hi), axis="SoluteDB corpus labels vs WSU measurement",
            baseline="zero (corpus agrees with measurement)",
            baseline_value=0.0, n=n).render())

    # ====================================================== 4. RETENTION SPACE
    print("\n" + "=" * 78)
    print("4. RETENTION SPACE (#41 criterion) -- descriptors through fixed constants")
    print("=" * 78)
    per_cmpd_bias, per_cmpd_rmse = {}, {}
    per_mod = defaultdict(lambda: defaultdict(list))
    prop_sd = np.zeros(n)
    for i, c in enumerate(names):
        M = members[c]
        mu_d = {k: mu[k][i] for k in D}
        tr_d = {k: truth[k][i] for k in D}
        lb_d = {k: lab[k][i] for k in D}
        errs, spreads, lerrs = [], [], []
        for f in fits:
            ref = logk(f, tr_d, V[i])
            e = logk(f, mu_d, V[i]) - ref
            errs.append(e)
            lerrs.append(logk(f, lb_d, V[i]) - ref)
            spreads.append(np.std([logk(f, m, V[i]) for m in M], ddof=1))
            per_mod[f["mod"]]["model"].append(e)
            per_mod[f["mod"]]["label"].append(lerrs[-1])
        errs = np.array(errs)
        per_cmpd_bias[c] = float(np.mean(errs))
        per_cmpd_rmse[c] = float(np.sqrt(np.mean(errs ** 2)))
        prop_sd[i] = float(np.median(spreads))

    rmse_arr = np.array([per_cmpd_rmse[c] for c in names])
    bias_arr = np.array([per_cmpd_bias[c] for c in names])
    for label, arr, base_txt, base_val in [
        ("RMSE(log k)", rmse_arr,
         "Poole 2019 harmonised-descriptor median LSER residual SD (MeOH)",
         POOLE_RESIDUAL_SD["methanol"]),
        ("bias(log k)", bias_arr, "zero (unbiased)", 0.0),
    ]:
        p, lo, hi = boot(n, lambda idx, a=arr: float(np.mean(a[idx])))
        print(f"  mean per-compound {label}: " + fig(
            p, (lo, hi), axis=RET_AXIS, baseline=base_txt,
            baseline_value=base_val, n=n, units=" log k").render())

    print(f"\n  per-compound RMSE(log k) quantiles: "
          + " ".join(f"p{q}={np.percentile(rmse_arr, q):.3f}" for q in qs))
    order = np.argsort(-rmse_arr)[:5]
    print("  worst 5: " + ", ".join(f"{names[i]} {rmse_arr[i]:.3f}" for i in order))

    print("\n  by modifier (all compound x fit cells):")
    print(f"{'modifier':>14} {'cells':>7} {'model RMSE':>11} {'label RMSE':>11} "
          f"{'model bias':>11} {'Poole SD':>9}")
    for m in sorted(per_mod):
        a = np.array(per_mod[m]["model"])
        b = np.array(per_mod[m]["label"])
        print(f"{m:>14} {len(a):7d} {np.sqrt(np.mean(a**2)):11.3f} "
              f"{np.sqrt(np.mean(b**2)):11.3f} {np.mean(a):+11.3f} "
              f"{POOLE_RESIDUAL_SD.get(m, float('nan')):9.3f}")

    # ================================================ 4b. WEAK-BASE STRATUM
    print("\n  by S-2 class -- the anchor's own split (neutral vs weak base):")
    print(f"{'class':>12} {'n':>4} " + " ".join(f"{'bias ' + k:>10}" for k in D)
          + f" {'ret RMSE':>9} {'ret bias':>9}")
    for label, mask in [("neutral", ~is_base), ("weak base", is_base)]:
        m = np.where(mask)[0]
        print(f"{label:>12} {len(m):4d} "
              + " ".join(f"{np.mean(model_err[k][m]):+10.4f}" for k in D)
              + f" {np.mean(rmse_arr[m]):9.3f} {np.mean(bias_arr[m]):+9.3f}")
    # The contrast is the claim, and one arm is n=15, so it gets an interval:
    # resample within each class so both arms stay populated.
    idx_b, idx_n = np.where(is_base)[0], np.where(~is_base)[0]
    rng = np.random.default_rng(20260816)
    for k in D:
        e = model_err[k]
        gap = float(np.mean(e[idx_b]) - np.mean(e[idx_n]))
        draws = np.array([
            np.mean(e[rng.choice(idx_b, idx_b.size)])
            - np.mean(e[rng.choice(idx_n, idx_n.size)])
            for _ in range(2000)])
        lo, hi = np.percentile(draws, [2.5, 97.5])
        flag = "   <-- interval excludes 0" if not (lo <= 0 <= hi) else ""
        print(f"    {k}: weak-base minus neutral bias {gap:+.4f} "
              f"[{lo:+.4f}, {hi:+.4f}] (n=15 vs 79){flag}")

    # ============================================== 5. ENSEMBLE DISPERSION
    print("\n" + "=" * 78)
    print("5. ENSEMBLE UNDER-DISPERSION -- retention RMSE / propagated ensemble SD")
    print("=" * 78)
    ratio = rmse_arr / prop_sd
    p, lo, hi = boot(n, lambda idx, r=ratio: float(np.median(r[idx])))
    print("  median ratio: " + fig(
        p, (lo, hi), axis=RET_AXIS,
        baseline="1.0 (ensemble spread matches its own error)",
        baseline_value=1.0, n=n).render())
    print(f"  quantiles: " + " ".join(f"p{q}={np.percentile(ratio, q):.2f}" for q in qs))
    print("  descriptor-space, per descriptor (RMSE / median ensemble SD):")
    for k in D:
        print(f"    {k}: RMSE {np.sqrt(np.mean(model_err[k]**2)):.3f} / "
              f"ens SD {np.median(esd[k]):.3f} = "
              f"{np.sqrt(np.mean(model_err[k]**2)) / np.median(esd[k]):.2f}")

    # ------------------------------------------------------------- artefacts
    with open(os.path.join(HERE, "residuals.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound", "tier", "weak_base", "V"]
                   + [f"meas_{k}" for k in D] + [f"pred_{k}" for k in D]
                   + [f"db_{k}" for k in D] + [f"err_{k}" for k in D]
                   + [f"ens_sd_{k}" for k in D]
                   + ["ret_bias_logk", "ret_rmse_logk", "ret_ens_sd_logk"])
        for i, c in enumerate(names):
            w.writerow([c, tier[c], int(weak_base[c]), f"{V[i]:.4f}"]
                       + [f"{truth[k][i]:.4f}" for k in D]
                       + [f"{mu[k][i]:.4f}" for k in D]
                       + [f"{lab[k][i]:.4f}" for k in D]
                       + [f"{model_err[k][i]:+.4f}" for k in D]
                       + [f"{esd[k][i]:.4f}" for k in D]
                       + [f"{per_cmpd_bias[c]:+.4f}", f"{per_cmpd_rmse[c]:.4f}",
                          f"{prop_sd[i]:.4f}"])
    print(f"\nwrote residuals.csv ({n} compounds)")


if __name__ == "__main__":
    main()
