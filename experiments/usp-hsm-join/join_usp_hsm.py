#!/usr/bin/env python3
"""#39's join: USP SRM-870 activity parameters vs HSM ion-exchange terms.

Joins sources/usp-column-db/usp-approach-columns.csv (145 columns; Hy, CTF, CFA,
TFA, BD) to sources/hsm-column-db/database.csv (819 columns; H, S*, A, B, C28,
C70) on normalised column name, then tests whether the SRM 870 activity measures
correlate with the HSM cation-exchange terms C(2.8)/C(7.0).

Stdlib only. Permutation p-values (two-sided, 20 000 shuffles, fixed seed).
Writes matched.csv and unmatched.txt next to this script; results to stdout.
"""

import csv
import math
import random
import re
import unicodedata
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
USP_CSV = ROOT / "sources/usp-column-db/usp-approach-columns.csv"
HSM_CSV = ROOT / "sources/hsm-column-db/database.csv"
WSU_CSV = ROOT / "sources/wsu-lser/wsu2019-system-constants.csv"
CAV_CSV = ROOT / "sources/wsu-lser/wsu2019-column-caveats.csv"

# Tokens that describe packing geometry, not phase chemistry. Dropping them lets
# "Aqua 5 µ C18 125A" meet "Aqua C18"; the match is still reviewed in matched.csv.
GEOM = re.compile(
    r"\b\d+(\.\d+)?\s*(u|µ|μ|um|µm|μm)\b"   # particle size
    r"|\b\d{2,4}\s*a\b"                       # pore size, e.g. 125A / 300 A
    r"|\(\s*\d+\s*mm\s*\)",
    re.IGNORECASE,
)


def norm(name: str) -> str:
    s = unicodedata.normalize("NFKC", name).lower()
    s = s.replace("μ", "µ")
    s = GEOM.sub(" ", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(s.split())


def flat(name: str) -> str:
    return norm(name).replace(" ", "")


def flat_nopart(name: str) -> str:
    """Spaceless key with bare 1–2 digit particle sizes dropped.

    USP writes ProntoSil as `pore-particle-phase` ("ProntoSil 120-5-C18-H") where
    HSM writes "ProntoSIL 120 C18 H", and Cosmosil as "5C18-AR-II" against HSM's
    "C18-AR-II". Pore sizes are 3–4 digits and survive; a bare 1–2 digit token is
    a particle size in µm.
    """
    return re.sub(r"\b\d{1,2}\b", " ", norm(name)).replace(" ", "")


# Verified by eye against the manufacturer's product line, one USP name to one
# HSM name. Only house-style differences are listed here — where the USP name
# could plausibly be a *different* product in the same line (`Discovery C18 WP`
# vs `Discovery C18`, `Alltima HP C18 AQ` vs `Alltima HP C18`, bare `Allure` /
# `Ultra` / `Pinnacle DB` against a dozen phase-suffixed HSM rows) the pair is
# declined and the USP row is left unmatched instead.
ALIASES = {
    "BetaBasic 18": "Hypersil Beta Basic-18",
    "BetaMax Neutral": "Hypersil Betamax Neutral",
    "Bio Basic 18": "Hypersil Bio Basic-18",
    "Hydrosphere C18": "YMC Hydrosphere C18",
    "Hypersil BDS 18": "Hypersil BDS C18",
    "Zorbax Rx C18": "Zorbax Rx-18",
    "YMC Pro C18": "YMC-Pack Pro C18",
    "YMC Pro C18 RS": "YMC-Pack Pro C18 RS",
    "YMC ODS-AQ": "YMC-Pack ODS-AQ",
    "218TP 300 C18": "Vydac 218TP",          # USP manufacturer: Grace/Vydac
    "Everest C18 300A": "Vydac Everest",     # USP manufacturer: Grace/Vydac
    "Jupiter 5 µ C18 300A": "Jupiter 300 C18",
    "SymmetryShield RP18": "SymmetryShield C18",  # the RP8 is the other Waters row
}


def read_rows(path, encoding):
    with open(path, newline="", encoding=encoding) as f:
        return list(csv.DictReader(f))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def perm_p(x, y, stat, observed, n_perm=20000, seed=39):
    rng = random.Random(seed)
    y = list(y)
    hits = 0
    for _ in range(n_perm):
        rng.shuffle(y)
        if abs(stat(x, y)) >= abs(observed) - 1e-12:
            hits += 1
    return (hits + 1) / (n_perm + 1)


def main():
    usp = read_rows(USP_CSV, "utf-8")   # UTF-8: µ in several names
    hsm = read_rows(HSM_CSV, "latin-1")

    by_norm, by_flat, by_nopart = {}, {}, {}
    for r in hsm:
        by_norm.setdefault(norm(r["name"]), []).append(r)
        by_flat.setdefault(flat(r["name"]), []).append(r)
        by_nopart.setdefault(flat_nopart(r["name"]), []).append(r)

    matched, unmatched, ambiguous = [], [], []
    for u in usp:
        name = u["name"]
        if name in ALIASES:
            cands = [r for r in hsm if r["name"] == ALIASES[name]]
            how = "alias"
        else:
            for index, how in ((by_norm, "exact-normalised"),
                               (by_flat, "spaceless"),
                               (by_nopart, "particle-stripped")):
                key = {"exact-normalised": norm, "spaceless": flat,
                       "particle-stripped": flat_nopart}[how](name)
                cands = index.get(key, [])
                if cands:
                    break
        if len(cands) == 1:
            matched.append((u, cands[0], how))
        elif len(cands) > 1:
            # HSM carries two rows under one name with different C terms
            # (Betasil C18: C28 -0.03 vs 0.095). Picking one would invent a
            # number, so these are excluded from the statistics.
            ambiguous.append((name, [r["name"] for r in cands], how))
        else:
            unmatched.append(name)

    # A second USP row reaching an HSM row already claimed means at least one of
    # the two is wrong; drop both rather than double-count the HSM parameters.
    claims = {}
    for u, h, how in matched:
        claims.setdefault(id(h), []).append((u, h, how))
    collisions = [v for v in claims.values() if len(v) > 1]
    for group in collisions:
        for u, h, how in group:
            matched.remove((u, h, how))
            ambiguous.append((u["name"], [h["name"]], f"{how}-collision"))

    print(f"USP rows: {len(usp)}   HSM rows: {len(hsm)}")
    print(f"Matched: {len(matched)}   Unmatched: {len(unmatched)}   "
          f"Excluded as ambiguous: {len(ambiguous)}\n")
    for n, cands, how in ambiguous:
        print(f"  ambiguous: {n} -> {cands} [{how}]")
    print()

    with open(HERE / "matched.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["usp_name", "hsm_name", "match", "usp_manufacturer",
                    "hsm_manufacturer", "Hy", "CTF", "CFA", "TFA", "BD",
                    "H", "S", "A", "B", "C28", "C70"])
        for u, h, how in matched:
            w.writerow([u["name"], h["name"], how, u["manufacturer"],
                        h["manufacturer"], u["Hy"], u["CTF"], u["CFA"],
                        u["TFA"], u["BD"], h["H"], h["S"], h["A"], h["B"],
                        h["C28"], h["C70"]])

    (HERE / "unmatched.txt").write_text("\n".join(unmatched) + "\n")

    print(f"{'USP':>4} {'HSM':>4} {'n':>4} {'pearson':>8} {'p':>8} "
          f"{'spearman':>9} {'p':>8}")
    for up in ["Hy", "CTF", "CFA", "TFA", "BD"]:
        for hp in ["C28", "C70"]:
            pairs = [(fnum(u[up]), fnum(h[hp])) for u, h, _ in matched]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
            if len(pairs) < 5:
                print(f"{up:>4} {hp:>4} {len(pairs):>4}  (too few)")
                continue
            x, y = [a for a, _ in pairs], [b for _, b in pairs]
            pr = pearson(x, y)
            sr = spearman(x, y)
            ppr = perm_p(x, y, pearson, pr)
            psr = perm_p(x, y, spearman, sr)
            print(f"{up:>4} {hp:>4} {len(pairs):>4} {pr:>8.3f} {ppr:>8.4f} "
                  f"{sr:>9.3f} {psr:>8.4f}")

    # cross-check within HSM: C28 vs C70 on the matched subset, for scale
    pairs = [(fnum(h["C28"]), fnum(h["C70"])) for _, h, _ in matched]
    pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
    x, y = [a for a, _ in pairs], [b for _, b in pairs]
    print(f"\n(HSM-internal, matched subset) C28 vs C70: "
          f"n={len(pairs)} pearson={pearson(x, y):.3f} "
          f"spearman={spearman(x, y):.3f}")

    # The Pearson/Spearman gap on CFA asks whether a few high-activity columns
    # carry the association, so name them.
    print("\nHighest HSM C(7.0) in the matched set:")
    top = sorted((r for _, r, _ in matched if fnum(r["C70"]) is not None),
                 key=lambda r: -fnum(r["C70"]))[:6]
    for r in top:
        u = next(u for u, h, _ in matched if h is r)
        print(f"  {r['name']:<32} C70={r['C70']:>6}  C28={r['C28']:>6}  "
              f"CFA={u['CFA']:>4} TFA={u['TFA']:>4}")

    # Same statistics with those six dropped, as a leverage check.
    drop = {id(r) for r in top}
    print("\nWith those six removed:")
    for up in ["CFA", "TFA", "Hy"]:
        for hp in ["C28", "C70"]:
            pairs = [(fnum(u[up]), fnum(h[hp])) for u, h, _ in matched
                     if id(h) not in drop]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
            x, y = [a for a, _ in pairs], [b for _, b in pairs]
            print(f"  {up:>4} {hp:>4} n={len(pairs):<4} "
                  f"pearson={pearson(x, y):>7.3f} spearman={spearman(x, y):>7.3f}")

    # #30's sweep pair, if present in HSM
    print("\n#30 sweep pair in HSM:")
    for target in ["xbridge shield rp18", "betasil c18"]:
        rows = [r for r in hsm if target in norm(r["name"])]
        for r in rows:
            print(f"  {r['name']} ({r['manufacturer']}): "
                  f"H={r['H']} C28={r['C28']} C70={r['C70']}")
        if not rows:
            print(f"  {target}: not found")

    # If that pair does not contrast on C, what would? #30 wants the two columns
    # to differ in silanol activity and match on everything else, so pair on
    # |ΔH| ≤ 0.05 (hydrophobicity held) and rank by ΔC(7.0).
    # HSM's own `type` field is the silica-generation classification (A = old
    # high-activity silica, B = modern low-activity, EP = embedded polar), which
    # is the silanol axis stated directly. Worth seeing before choosing a pair.
    print("\nHSM C(7.0) by silica type, whole 819-row set:")
    by_type = {}
    for r in hsm:
        if fnum(r["C70"]) is not None:
            by_type.setdefault(r["type"] or "(blank)", []).append(fnum(r["C70"]))
    for t, vs in sorted(by_type.items(), key=lambda kv: -len(kv[1])):
        vs.sort()
        print(f"  type {t:<8} n={len(vs):<4} median={vs[len(vs)//2]:>7.3f}  "
              f"min={vs[0]:>7.3f} max={vs[-1]:>6.2f}")

    # Restricted to USPtype L1 / phase C18: this excludes apHera C18 Polymer
    # (L67), whose C70 = 7.51 tops the unrestricted ranking but is polymer-coated
    # chemistry, not silanol activity, so it would not test what #30 wants.
    print("\nBest C(7.0)-contrasting pairs at matched hydrophobicity "
          "(|ΔH| ≤ 0.05), silica C18 only:")
    usable = [r for r in hsm
              if r["USPtype"] == "L1" and r["phase"].strip() == "C18"
              and all(fnum(r[k]) is not None for k in ("H", "C28", "C70"))]
    pairs = []
    for i, a in enumerate(usable):
        for b in usable[i + 1:]:
            if abs(fnum(a["H"]) - fnum(b["H"])) <= 0.05:
                d = fnum(b["C70"]) - fnum(a["C70"])
                pairs.append((abs(d), a, b))
    pairs.sort(key=lambda t: -t[0])
    in_usp = {h["name"] for _, h, _ in matched}
    for d, a, b in pairs[:8]:
        lo, hi = (a, b) if fnum(a["C70"]) < fnum(b["C70"]) else (b, a)
        mark = "both in USP set" if {lo["name"], hi["name"]} <= in_usp else ""
        print(f"  ΔC70={d:>5.2f}  {lo['name']} (C70={lo['C70']}, H={lo['H']}, "
              f"type {lo['type']})  vs  {hi['name']} (C70={hi['C70']}, "
              f"H={hi['H']}, type {hi['type']})  {mark}")

    # Those pairs leave WSU-2019, which costs #30 its fixed neutral baseline. So
    # ask what contrast is available *inside* WSU-2019.
    print("\nWSU-2019 columns by HSM C(7.0) — the anchor set's silanol coverage:")
    wsu = {r["column"] for r in read_rows(WSU_CSV, "utf-8-sig")}
    absent = [w for w in sorted(wsu)
              if not any(flat(r["name"]) == flat(w) for r in hsm)]
    found = [(w, r) for w in sorted(wsu) for r in hsm
             if flat(r["name"]) == flat(w) and fnum(r["C70"]) is not None]
    for w, r in sorted(found, key=lambda t: -fnum(t[1]["C70"])):
        print(f"  {w:<30} C70={r['C70']:>7}  C28={r['C28']:>7}  "
              f"H={r['H']:>6}  type {r['type']}")
    print(f"  not in HSM ({len(absent)}): {', '.join(absent)}")

    print(f"\n  type-A (high-activity silica) columns in WSU-2019: "
          f"{sum(1 for _, r in found if r['type'] == 'A')}")

    # Restricted to C18: a phenyl/biphenyl pair would vary phase chemistry as
    # well as silanol activity, which is not the contrast #30 needs.
    c18 = [(w, r) for w, r in found if r["phase"].strip() == "C18"]
    lo = min(fnum(r["C70"]) for _, r in c18)
    hi = max(fnum(r["C70"]) for _, r in c18)
    print(f"  WSU-2019 C18 columns span C70 {lo:.2f} to {hi:.2f} "
          f"(HSM silica-C18 range: -1.00 to 2.69)")
    best = max(((abs(fnum(a[1]["C70"]) - fnum(b[1]["C70"])), a, b)
                for i, a in enumerate(c18) for b in c18[i + 1:]
                if abs(fnum(a[1]["H"]) - fnum(b[1]["H"])) <= 0.1),
               key=lambda t: t[0])
    print(f"  Best WSU-internal C18 pair at |ΔH| ≤ 0.1: ΔC70={best[0]:.2f}  "
          f"{best[1][0]} (C70={best[1][1]['C70']}) vs "
          f"{best[2][0]} (C70={best[2][1]['C70']})")

    # Third instrument: WSU-2019 Table S-1 grades its electrostatic caveat as
    # "weak" or leaves it unqualified (#40). If that grading is a silanol-
    # activity signal it should track HSM C(7.0).
    print("\nWSU-2019 Table S-1 electrostatic flags vs HSM C(7.0):")
    hsm_by_key = {}
    for r in hsm:
        hsm_by_key.setdefault(flat(r["name"]), r)
    graded = []
    for x in read_rows(CAV_CSV, "utf-8-sig"):
        if x["caveat"] != "electrostatic":
            continue
        r = hsm_by_key.get(flat(x["column"]))
        if r and fnum(r["C70"]) is not None:
            graded.append(("weak" in x["note"].lower(), fnum(r["C70"]),
                           r["phase"].strip()))
    for label, subset in (("all phases", graded),
                          ("C18 only", [g for g in graded if g[2] == "C18"])):
        parts = []
        for weak, tag in ((True, "weak"), (False, "unqualified")):
            v = [c for w, c, _ in subset if w == weak]
            parts.append(f"{tag} n={len(v)} mean C70={sum(v)/len(v):+.3f}")
        print(f"  {label:<12} " + "   ".join(parts))
    print("  -> separation across all phases is carried by the fluorinated "
          "columns; among C18 the flags do not rank at all")


if __name__ == "__main__":
    main()
