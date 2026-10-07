"""#42: does SoluteML ensemble spread stratify error — in descriptor space AND retention space?

Grading follows #41: retention space is the criterion, descriptor rmse is diagnostic.
Fixed-constants form (we hold no raw k): predicted descriptors vs WSU S-2 measured
descriptors, both pushed through the same published LSER constants.
Wetting-masked fits excluded per #40.
"""
import os
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PREDS = os.path.join(REPO, "prototype", "descriptor-covariance", "soluteml_preds_raw.csv")
import csv, math, statistics as st
from collections import defaultdict

D = ["E", "S", "A", "B"]

# --- measured anchor descriptors (WSU Table S-2)
meas = {}
for r in csv.DictReader(open(os.path.join(REPO, "sources/wsu-lser/wsu2019-descriptors-s2.csv"))):
    meas[r["compound"].strip()] = {k: float(r[k]) for k in ["E", "S", "A", "B", "V"]}

# --- ensemble members
mem = defaultdict(list)
for r in csv.DictReader(open(PREDS)):
    mem[r["compound"].strip()].append({k: float(r[k]) for k in D})

common = sorted(set(meas) & set(mem))
print(f"compounds: S-2 {len(meas)}, preds {len(mem)}, matched {len(common)}")
print(f"members per compound: {sorted({len(v) for v in mem.values()})}")

# --- system constants, wetting-masked rows dropped (#40)
fits = []
for r in csv.DictReader(open(os.path.join(REPO, "sources/wsu-lser/wsu2019-system-constants.csv"))):
    if str(r.get("incomplete_wetting", "")).strip().lower() in ("true", "1", "yes"):
        continue
    fits.append({k: float(r[k]) for k in ["c", "e", "s", "a", "b", "v"]} |
                {"mod": r["modifier"], "col": r["column"], "phi": r["phi_pct_vv"]})
print(f"unmasked fits: {len(fits)}")


def logk(f, d, V):
    return f["c"] + f["e"] * d["E"] + f["s"] * d["S"] + f["a"] * d["A"] + f["b"] * d["B"] + f["v"] * V


def tertiles(vals):
    v = sorted(vals)
    n = len(v)
    return v[n // 3], v[2 * n // 3]


rows = {}
for cmpd in common:
    M = mem[cmpd]
    mu = {k: st.mean(m[k] for m in M) for k in D}
    sd = {k: st.stdev([m[k] for m in M]) for k in D}
    truth = meas[cmpd]
    V = truth["V"]

    # descriptor-space error (diagnostic, #41)
    derr = {k: mu[k] - truth[k] for k in D}

    # retention-space: error of mean prediction, and ensemble spread, per fit
    errs, spreads = [], []
    for f in fits:
        ref = logk(f, truth, V)
        errs.append(logk(f, mu, V) - ref)
        member_lk = [logk(f, m, V) for m in M]
        spreads.append(st.stdev(member_lk))
    rows[cmpd] = dict(mu=mu, sd=sd, derr=derr,
                      rmse_ret=math.sqrt(st.mean(e * e for e in errs)),
                      bias_ret=st.mean(errs),
                      spread_ret=st.median(spreads),
                      sd_sum=sum(sd.values()))

print("\n=== DESCRIPTOR SPACE (diagnostic; Ulrich-style per-descriptor tertiles) ===")
print(f"{'desc':>4} {'rmse all':>9} {'low-SD 3rd':>11} {'mid':>8} {'high-SD 3rd':>12} {'ratio hi/lo':>12}")
for k in D:
    lo, hi = tertiles([rows[c]["sd"][k] for c in common])
    grp = defaultdict(list)
    for c in common:
        s = rows[c]["sd"][k]
        grp["lo" if s <= lo else ("hi" if s > hi else "mid")].append(rows[c]["derr"][k])
    rm = lambda g: math.sqrt(st.mean(e * e for e in grp[g]))
    allr = math.sqrt(st.mean(rows[c]["derr"][k] ** 2 for c in common))
    print(f"{k:>4} {allr:9.3f} {rm('lo'):11.3f} {rm('mid'):8.3f} {rm('hi'):12.3f} {rm('hi')/rm('lo'):12.2f}"
          f"   (n {len(grp['lo'])}/{len(grp['mid'])}/{len(grp['hi'])})")

print("\n=== RETENTION SPACE (the #41 criterion; log10 k over all unmasked fits) ===")
for label, key in [("propagated ensemble spread in log k", "spread_ret"),
                   ("naive sum of descriptor SDs", "sd_sum")]:
    lo, hi = tertiles([rows[c][key] for c in common])
    grp = defaultdict(list)
    for c in common:
        s = rows[c][key]
        grp["lo" if s <= lo else ("hi" if s > hi else "mid")].append(rows[c]["rmse_ret"])
    q = lambda g: math.sqrt(st.mean(x * x for x in grp[g]))
    allr = math.sqrt(st.mean(rows[c]["rmse_ret"] ** 2 for c in common))
    print(f"\nstratified by {label}:")
    print(f"  rmse all {allr:.3f} | low third {q('lo'):.3f} | mid {q('mid'):.3f} | high third {q('hi'):.3f}"
          f" | ratio hi/lo {q('hi')/q('lo'):.2f}")
    print(f"  n {len(grp['lo'])}/{len(grp['mid'])}/{len(grp['hi'])}")

# --- Q2 evidence: is spread usable continuously? correlation + calibration ratio
xs = [rows[c]["spread_ret"] for c in common]
ys = [rows[c]["rmse_ret"] for c in common]
def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return num / math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
def spearman(a, b):
    ra = {v: i for i, v in enumerate(sorted(a))}
    rb = {v: i for i, v in enumerate(sorted(b))}
    return pearson([ra[x] for x in a], [rb[y] for y in b])
print(f"\nspread vs error, per compound (n={len(common)}): pearson r {pearson(xs, ys):.3f}, "
      f"spearman rho {spearman(xs, ys):.3f}")

ratio = [y / x for x, y in zip(xs, ys) if x > 1e-9]
ratio.sort()
print(f"calibration ratio (retention rmse / propagated ensemble sd): median {st.median(ratio):.2f}, "
      f"IQR {ratio[len(ratio)//4]:.2f}-{ratio[3*len(ratio)//4]:.2f}, "
      f"90th pct {ratio[int(.9*len(ratio))]:.2f}, max {ratio[-1]:.2f}")

# worst offenders: confident but wrong (the point-4 class)
bad = sorted(common, key=lambda c: -(rows[c]["rmse_ret"] / max(rows[c]["spread_ret"], 1e-9)))[:10]
print("\nconfident-but-wrong (highest error/spread) — candidate point-4 structural class:")
for c in bad:
    r = rows[c]
    print(f"  {c:<34} spread {r['spread_ret']:.3f}  rmse {r['rmse_ret']:.3f}  "
          f"ratio {r['rmse_ret']/r['spread_ret']:5.1f}  dE {r['derr']['E']:+.2f} dS {r['derr']['S']:+.2f} "
          f"dA {r['derr']['A']:+.2f} dB {r['derr']['B']:+.2f}")
