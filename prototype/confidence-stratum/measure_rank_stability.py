import os
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PREDS = os.path.join(REPO, "prototype", "descriptor-covariance", "soluteml_preds_raw.csv")
import csv, math, statistics as st
from collections import defaultdict
D=["E","S","A","B"]
meas={r["compound"].strip():{k:float(r[k]) for k in ["E","S","A","B","V"]} for r in csv.DictReader(open(os.path.join(REPO, "sources/wsu-lser/wsu2019-descriptors-s2.csv")))}
mem=defaultdict(list)
for r in csv.DictReader(open(PREDS)): mem[r["compound"].strip()].append({k:float(r[k]) for k in D})
common=sorted(set(meas)&set(mem))
fits=[]
for r in csv.DictReader(open(os.path.join(REPO, "sources/wsu-lser/wsu2019-system-constants.csv"))):
    if str(r.get("incomplete_wetting","")).strip().lower() in ("true","1","yes"): continue
    fits.append({k:float(r[k]) for k in ["c","e","s","a","b","v"]}|{"mod":r["modifier"],"col":r["column"],"phi":int(float(r["phi_pct_vv"]))})
lk=lambda f,d,V: f["c"]+f["e"]*d["E"]+f["s"]*d["S"]+f["a"]*d["A"]+f["b"]*d["B"]+f["v"]*V
def pearson(a,b):
    ma,mb=st.mean(a),st.mean(b); num=sum((x-ma)*(y-mb) for x,y in zip(a,b))
    return num/math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
def spear(a,b):
    ra={v:i for i,v in enumerate(sorted(a))}; rb={v:i for i,v in enumerate(sorted(b))}
    return pearson([ra[x] for x in a],[rb[y] for y in b])
# spread per compound per fit
spread={c:{} for c in common}
for c in common:
    V=meas[c]["V"]
    for f in fits:
        spread[c][(f["mod"],f["col"],f["phi"])]=st.stdev([lk(f,m,V) for m in mem[c]])
keys=list(spread[common[0]].keys())
# rank correlation of per-compound spread between every pair of methods (sample)
import random; random.seed(0)
pairs=random.sample([(i,j) for i in range(len(keys)) for j in range(i+1,len(keys))],3000)
rhos=[spear([spread[c][keys[i]] for c in common],[spread[c][keys[j]] for c in common]) for i,j in pairs]
rhos.sort()
print(f"rank stability of per-compound spread across method pairs (n={len(rhos)} pairs of {len(keys)} methods):")
print(f"  spearman rho: median {st.median(rhos):.3f}, 5th pct {rhos[int(.05*len(rhos))]:.3f}, min {rhos[0]:.3f}")
# extremes: methanol 10% vs acetonitrile 70%
def show(mod,phi):
    ks=[k for k in keys if k[0].startswith(mod) and k[2]==phi]
    if not ks: return None
    return [st.mean(spread[c][k] for k in ks) for c in common]
a=show("methanol",10); b=show("acetonitrile",70)
if a and b: print(f"  extreme comparison MeOH 10% vs ACN 70%: rho {spear(a,b):.3f}")
# how much does absolute spread scale by method?
mags=sorted(((st.median(spread[c][k] for c in common),k) for k in keys))
print("  spread magnitude ranges %.3f (%s %s %d%%) to %.3f (%s %s %d%%) -> %.1fx" % (mags[0][0],mags[0][1][0][:4],mags[0][1][1],mags[0][1][2],mags[-1][0],mags[-1][1][0][:4],mags[-1][1][1],mags[-1][1][2],mags[-1][0]/mags[0][0]))
