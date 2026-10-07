#!/usr/bin/env python3
"""#45's bounded RepoRT search: can public data validate the thin slice?

Scans every RepoRT processed dataset for the 15 valid thin-slice compounds
(#21), matched on InChIKey skeleton so stereochemistry and salt form do not
block a match, then checks whether the matching datasets carry the conditions
#45 needs: isocratic retention, or a gradient invertible to `k`.

Usage:  python3 search_report.py /path/to/RepoRT
RepoRT: https://github.com/michaelwitting/RepoRT (CC BY-SA 4.0), `master`.
Requires RDKit.
"""
import csv
import glob
import os
import sys
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "prototype", "thin-slice"))
from compounds import COMPOUNDS  # noqa: E402

REPO = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("REPORT_DIR", "")
if not os.path.isdir(os.path.join(REPO, "processed_data")):
    sys.exit("usage: search_report.py /path/to/RepoRT")


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def elution_mode(gpath):
    """'iso' | 'grad' | 'unknown'.

    RepoRT writes a header-only gradient file when the programme was never
    recorded — several PredRet-sourced sets say exactly that in
    `info['missing information']`. An empty file is **missing metadata, not a
    constant composition**, and must never be counted as isocratic. Getting
    this wrong is what made an early pass of this search report nine usable
    isocratic datasets that do not exist.
    """
    if not os.path.exists(gpath):
        return "unknown"
    rows = read_tsv(gpath)
    chans = [k for k in (rows[0] if rows else {}) if k.strip().endswith("[%]")]
    filled = [r for r in rows if any((r.get(c) or "").strip() for c in chans)]
    if not filled:
        return "unknown"
    if len(filled) == 1:
        return "iso"
    for ch in chans:
        if len({(r[ch] or "").strip() for r in filled if (r[ch] or "").strip()}) > 1:
            return "grad"
    return "iso"


want = {}
for c in COMPOUNDS:
    if c["type"] == "refused":
        continue          # no retention distribution by construction (#33)
    key = Chem.MolToInchiKey(Chem.MolFromSmiles(c["smiles"]))
    want[key.split("-")[0]] = c["name"]

hits = defaultdict(dict)
for f in sorted(glob.glob(f"{REPO}/processed_data/*/*_rtdata_canonical_success.tsv")):
    ds = os.path.basename(f).split("_")[0]
    for r in read_tsv(f):
        skel = (r.get("inchikey.std") or "").split("-")[0]
        if skel in want:
            hits[ds][want[skel]] = float(r["rt"])

total = len(glob.glob(f"{REPO}/processed_data/*/"))
print(f"Seeking {len(want)} thin-slice compounds across {total} datasets")
print(f"{len(hits)} datasets contain at least one\n")

rows = []
for ds, found in hits.items():
    mf = f"{REPO}/processed_data/{ds}/{ds}_metadata.tsv"
    meta = read_tsv(mf)[0] if os.path.exists(mf) else {}
    rows.append(dict(ds=ds, n=len(found), meta=meta, cmpds=found,
                     mode=elution_mode(f"{REPO}/processed_data/{ds}/{ds}_gradient.tsv")))
rows.sort(key=lambda r: -r["n"])

print("Top overlap:")
print(f"  {'ds':<6}{'n':<4}{'mode':<9}{'pH':<7}{'T':<5}column")
for r in rows[:8]:
    m = r["meta"]
    print(f"  {r['ds']:<6}{r['n']:<4}{r['mode']:<9}"
          f"{(m.get('eluent.A.pH') or '-'):<7}{(m.get('column.temperature') or '-'):<5}"
          f"{(m.get('column.name') or '?')[:40]}")

print("\n1. Elution mode over overlapping datasets:", dict(Counter(r["mode"] for r in rows)))
print(f"   ISOCRATIC datasets with any slice compound: "
      f"{sum(1 for r in rows if r['mode'] == 'iso')}")

# 2. Dwell volume: needed to invert a gradient to k. Is the field even present?
fields = set(read_tsv(f"{REPO}/processed_data/{rows[0]['ds']}/"
                      f"{rows[0]['ds']}_metadata.tsv")[0])
dwell = [f for f in fields if any(w in f.lower()
                                  for w in ("dwell", "delay", "void", "dead", "system"))]
print(f"\n2. Dwell/system-volume fields in the RepoRT schema: {dwell or 'NONE'}")

# 3. Is column.t0 measured or imputed?
t0 = Counter()
for f in glob.glob(f"{REPO}/processed_data/*/*_metadata.tsv"):
    r = read_tsv(f)
    if r:
        t0[(r[0].get("column.t0") or "").strip()] += 1
print(f"\n3. column.t0 across {sum(t0.values())} datasets: {len(t0)} distinct values")
for v, n in t0.most_common(5):
    print(f"     t0={v!r:<12} {n} datasets")
print(f"     t0 == 0 (impossible; missing-value marker): {t0.get('0', 0)} datasets")

# 4. The most attractive paired-pH contrast is not two measurements.
def rt_table(ds):
    return {r["inchikey.std"]: float(r["rt"]) for r in
            read_tsv(f"{REPO}/processed_data/{ds}/{ds}_rtdata_canonical_success.tsv")
            if r.get("inchikey.std")}


print("\n4. Duplicate-RT check on the Aalizadeh sets (0382/0383 look like a pH pair):")
for a, b in [("0382", "0383"), ("0384", "0385"), ("0385", "0386")]:
    ta, tb = rt_table(a), rt_table(b)
    shared = set(ta) & set(tb)
    if not shared:
        continue
    same = sum(1 for k in shared if abs(ta[k] - tb[k]) < 1e-9)
    print(f"     {a} vs {b}: {len(shared):>5} shared, {same:>5} identical RT "
          f"({100 * same / len(shared):>3.0f}%)")
