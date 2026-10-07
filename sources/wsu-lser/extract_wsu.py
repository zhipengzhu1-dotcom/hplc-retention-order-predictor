import zipfile, re, csv, sys, os

DOCX = "papers/poole-2019-rplc-system-constants/1-s2.0-S0021967319303954-mmc1.docx"
if not os.path.exists(DOCX):
    sys.exit(f"{DOCX} not found — run this from the repo root. The supplementary docx is "
             f"gitignored (© Elsevier); see papers/poole-2019-rplc-system-constants/.")
z = zipfile.ZipFile(DOCX)
xml = z.read("word/document.xml").decode("utf8", errors="ignore")
# Preserve paragraph boundaries as spaces; strip all tags
text = re.sub(r"<[^>]+>", " ", xml)
text = text.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
text = re.sub(r"\s+", " ", text)
# Word runs split some names mid-token: "1 ,3 -Dibromobenzene" → make "1,3" one
# (non-numeric) token so it parses as part of the compound name.
text = text.replace(" ,", ",")
# Stray spaces inside SD parentheses, e.g. "(0.027 )( 0.022)".
text = text.replace("( ", "(").replace(" )", ")").replace(")(", ") (")

def clean_name(parts):
    s = " ".join(parts).replace(" -", "-").replace("- ", "-")
    s = s.replace("C 18", "C18").replace("C 8", "C8").replace("RP 18", "RP18")
    # Source typo in Table S-4: "Lunar Omega PS C18" (S-3 spells it "Luna").
    return s.replace("Lunar Omega", "Luna Omega")

def span(start_marker, end_marker):
    i = text.index(start_marker)
    j = text.index(end_marker, i + 1)
    return text[i:j]

NUM = re.compile(r"^-?\d+\.?\d*$")
# The source has one unclosed parenthesis ("(0.018" before XBridge Shield RP18
# in Table S-3), so the closing paren is optional.
PAR = re.compile(r"^\(-?\d+\.?\d*\)?$")

# ---------- Table S-2: descriptors ----------
# The table carries an in-body section heading with no cell of its own, so the
# name accumulator below swallows it and glues it onto the FIRST NAME OF THE
# BLOCK. It produced "Compounds designated as weak bases 2-Aminobiphenyl", which
# made every name join drop that compound: three analyses silently ran at n=93,
# and the compound lost is a weak base carrying the third-worst E error in the
# set -- the tail of the very effect the weak-base split measures. Found
# 2026-08-16 (see EXECUTION-PLAN.md).
#
# The heading is not noise to be stripped: it is the only record of which
# compounds Poole designates weak bases, so it becomes the `class` column.
S2_HEADING = "Compounds designated as weak bases"

s2 = span("Table S-2 Compounds", "Table S-3")
toks = s2.split()
# skip header up to 'V'
toks = toks[toks.index("V") + 1:]
compounds = []
i = 0
in_weak_bases = False
while i < len(toks):
    name = []
    while i < len(toks) and not NUM.match(toks[i]):
        name.append(toks[i]); i += 1
    if not name:
        break
    joined = " ".join(name)
    if joined.startswith(S2_HEADING):
        in_weak_bases = True
        name = joined[len(S2_HEADING):].split()
        if not name:
            sys.exit("S-2: section heading with no compound following it")
    vals = []
    while i < len(toks) and NUM.match(toks[i]) and len(vals) < 5:
        vals.append(toks[i]); i += 1
    if len(vals) != 5:
        sys.exit(f"S-2 parse error at {' '.join(name)}: got {vals}")
    compounds.append([clean_name(name)] +
                     ["weak_base" if in_weak_bases else "neutral"] + vals)

if not any(c[1] == "weak_base" for c in compounds):
    sys.exit(f"S-2: {S2_HEADING!r} never seen -- the class column would be a lie")
if any(S2_HEADING.split()[0] in c[0] for c in compounds):
    sys.exit("S-2: a section heading survived into a compound name")

with open("sources/wsu-lser/wsu2019-descriptors-s2.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["compound", "class", "E", "S", "A", "B", "V"])
    w.writerows(compounds)
print(f"S-2: {len(compounds)} compounds "
      f"({sum(c[1] == 'weak_base' for c in compounds)} weak bases)")

# ---------- Table S-1: per-column exclusion register ----------
# Word ran every cell of this table together without separators, so the notes are
# recovered by splitting on the three phrases the source uses to open one.
s1 = span("Table S-1", "Table S-2")
# Word run splits scatter spaces around hyphens ("Phenyl- Hexyl", "acetonitrile -water").
s1 = re.sub(r"\s*-\s*", "-", s1).replace("-and ", "- and ")
s1 = s1.replace("(See original manuscript for cited sources)", "")
# The Chromolith name wraps mid-cell: the reference and note land between
# "Chromolith" and the rest of the name. Rejoin it and drop the orphan.
s1 = s1.replace("Chromolith 25", "\0 25").replace(" Performance RP-18e", "")
s1 = s1.replace("\0", "Chromolith Performance RP-18e")
COLS = ["Ascentis C18", "Betasil C18", "Chromolith Performance RP-18e",
        "Discovery HS C18", "Discovery HS F5", "Fluophase-RP", "HyPURITY C18",
        "Kinetex C18", "Kinetex EVO C18", "Kinetex XB-C18", "Kinetex C8",
        "Kinetex Phenyl-Hexyl", "Kinetex F5", "Kinetex Biphenyl",
        "Luna Omega PS C18", "SunFire C18", "Synergi Hydro-RP",
        "Synergi Fusion-RP", "Synergi Polar-RP", "XBridge C18",
        "XBridge Shield RP18", "XBridge C8", "XBridge Phenyl",
        "XTerra MS C18", "XTerra Phenyl"]
def loose(phrase):
    """Match a phrase across the stray spaces Word's run splits leave inside words."""
    return r"\s*".join(re.escape(ch) for ch in phrase if not ch.isspace())

# Longest first so "Kinetex C18" cannot swallow the start of "Kinetex C8".
marks = []
for c in sorted(COLS, key=len, reverse=True):
    m = re.search(loose(c) + r"\s*(\d+(?:,\s*\d+)*)", s1)
    if not m:
        sys.exit(f"S-1: column '{c}' not found")
    marks.append((m.start(), c, m.group(1).replace(" ", ""), m.end()))
marks.sort()

NOTE = re.compile("(?=" + "|".join(
    loose(p) for p in ("Incomplete wetting", "Steric repulsion",
                       "Weak electrostatic interaction",
                       "Electrostatic interaction")) + ")")
KIND = [("incompletewetting", "incomplete_wetting"),
        ("stericrepulsion", "steric_repulsion"),
        ("electrostatic", "electrostatic")]
MODS = ["acetonitrile", "methanol", "tetrahydrofuran"]
caveats = []
for k, (_, col, ref, end) in enumerate(marks):
    stop = marks[k + 1][0] if k + 1 < len(marks) else len(s1)
    body = s1[end:stop].strip()
    notes = [n.strip() for n in NOTE.split(body) if n.strip()]
    if not notes:
        caveats.append([col, ref, "none", "", ""])
        continue
    for note in notes:
        flat = re.sub(r"\s+", "", note).lower()
        kind = next((v for p, v in KIND if p in flat), "other")
        mods = ";".join(m for m in MODS if m in flat)
        caveats.append([col, ref, kind, mods, re.sub(r"\s+", " ", note)])

with open("sources/wsu-lser/wsu2019-column-caveats.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["column", "reference", "caveat", "modifiers", "note"])
    w.writerows(caveats)
print(f"S-1: {len(marks)} columns, {len(caveats)} caveat rows")
for row in caveats:
    print("   ", row[0], "|", row[2], "|", row[3], "|", row[4][:90])

# ---------- Tables S-3/S-4/S-5: system constants ----------
tables = [
    ("S-3", "methanol", span("Table S-3", "Table S-4")),
    ("S-4", "acetonitrile", span("Table S-4", "Table S-5")),
    ("S-5", "tetrahydrofuran", span("Table S-5", "Figure S-1")),
]
rows = []
for tid, modifier, block in tables:
    toks = block.split()
    # skip header: everything up to first 'n' after 'SE'
    i = toks.index("SE")
    while toks[i] != "n":
        i += 1
    i += 1
    PHI = {"10", "20", "30", "40", "50", "60", "70"}
    ncols = set()
    col = None
    prev_phi = 0
    while i < len(toks):
        name = []
        # A bare small integer that is not a valid phi is a split-off piece of
        # the column name (e.g. "SunFire C 18"); decimals/negatives are data.
        while i < len(toks) and not PAR.match(toks[i]) and \
                (not NUM.match(toks[i]) or
                 (toks[i] not in PHI and "." not in toks[i]
                  and not toks[i].startswith("-"))):
            name.append(toks[i]); i += 1
        if name:
            col = clean_name(name)
            prev_phi = 0
        if i >= len(toks):
            break
        vals = []
        while i < len(toks) and NUM.match(toks[i]) and len(vals) < 12:
            vals.append(toks[i]); i += 1
        if len(vals) == 11 and vals[0] not in PHI:
            # Source omits the phi on one row (XBridge Phenyl 70%, S-4): infer it.
            inferred = str(prev_phi + 10)
            print(f"WARNING {tid} '{col}': row without phi, inferred {inferred}%")
            vals = [inferred] + vals
        sds = []
        while i < len(toks) and PAR.match(toks[i]) and len(sds) < 6:
            sds.append(toks[i].strip("()")); i += 1
        if len(vals) != 12 or len(sds) != 6:
            sys.exit(f"{tid} parse error at column '{col}': vals={vals} sds={sds}")
        phi = vals[0]
        if phi not in PHI or int(phi) <= prev_phi:
            sys.exit(f"{tid} bad phi sequence {prev_phi}->{phi} at column '{col}'")
        prev_phi = int(phi)
        r = float(vals[7])
        if not (0.9 < r <= 1.0):
            sys.exit(f"{tid} implausible r={r} at '{col}' phi={phi}")
        # The Luna Omega PS C18 block of Table S-3 prints SE before F, unlike
        # every other block. F is a Fisher statistic (>1, integer-valued) and SE
        # is a residual standard deviation (~0.02-0.15 in log k), so the two are
        # never confusable; swap on that basis rather than special-casing a name.
        if float(vals[9]) < 1.0 < float(vals[10]):
            print(f"WARNING {tid} '{col}' phi={phi}%: F/SE transposed in source, "
                  f"swapping ({vals[9]} <-> {vals[10]})")
            vals[9], vals[10] = vals[10], vals[9]
        F, SE = float(vals[9]), float(vals[10])
        if not (F > 1.0 and 0.0 < SE < 1.0):
            sys.exit(f"{tid} implausible F={F} SE={SE} at '{col}' phi={phi}")
        ncols.add(col)
        rows.append([modifier, col, phi] + vals[1:7] + sds + vals[7:])
    print(f"{tid} ({modifier}): {len(ncols)} columns, "
          f"{sum(1 for r_ in rows if r_[0] == modifier)} rows")

# ---------- Join the S-1 dewetting ranges onto the constants ----------
# Tables S-3/S-4/S-5 publish constants for every phi from 10% up, including the
# compositions Table S-1 says the column does not wet. Nothing in those tables
# marks them, so the flag has to be carried here or it is lost. Transcribed by
# hand from the S-1 notes because the prose does not parse reliably: "at 10%"
# excludes 10 only, "less than 30%" excludes 10 and 20.
#
# AMBIGUOUS: Discovery HS C18 reads "less than 30% (v/v) acetonitrile-water and
# 20% (v/v) methanol-water". Read as "less than" distributing over both, giving
# methanol {10}; the alternative reading ("at 20% methanol") would give {20}.
DEWET = {
    ("Ascentis C18", "acetonitrile"): {10},
    ("Discovery HS C18", "acetonitrile"): {10, 20},
    ("Discovery HS C18", "methanol"): {10},
    ("Discovery HS F5", "methanol"): {10, 20},
    ("Fluophase-RP", "acetonitrile"): {10},
    ("Fluophase-RP", "methanol"): {10},
    ("HyPURITY C18", "acetonitrile"): {10},
    ("Kinetex Biphenyl", "acetonitrile"): {10},
    ("Kinetex Biphenyl", "methanol"): {10},
    ("SunFire C18", "acetonitrile"): {10},
    ("SunFire C18", "tetrahydrofuran"): {10},
    ("Synergi Hydro-RP", "acetonitrile"): {10},
    ("Synergi Fusion-RP", "methanol"): {10},
    ("XBridge C18", "acetonitrile"): {10},
    ("XBridge C8", "acetonitrile"): {10},
    ("XBridge Phenyl", "acetonitrile"): {10},
    ("XTerra MS C18", "methanol"): {10, 20},
    ("XTerra MS C18", "acetonitrile"): {10},
    ("XTerra MS C18", "tetrahydrofuran"): {10},
}
flagged_cols = {c for c, _ in DEWET}
declared = {c for c, _, k, _, _ in caveats if k == "incomplete_wetting"}
if flagged_cols != declared:
    sys.exit(f"DEWET disagrees with Table S-1: {flagged_cols ^ declared}")
n_flag = 0
for r in rows:
    dewet = int(r[2]) in DEWET.get((r[1], r[0]), ())
    n_flag += dewet
    r.append("1" if dewet else "0")
print(f"incomplete wetting: {n_flag}/{len(rows)} published fits flagged, "
      f"across {len(flagged_cols)}/25 columns")

with open("sources/wsu-lser/wsu2019-system-constants.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["modifier", "column", "phi_pct_vv",
                "c", "e", "s", "a", "b", "v",
                "sd_c", "sd_e", "sd_s", "sd_a", "sd_b", "sd_v",
                "r", "r2", "F", "SE", "n", "incomplete_wetting"])
    w.writerows(rows)
print(f"total constant rows: {len(rows)}")

# Cross-check against a value read by eye from the raw text
assert ["methanol", "Ascentis C18", "10", "-0.339", "0.587", "-0.927",
        "-0.569", "-1.913", "3.350"] == rows[0][:9], rows[0]
print("spot check Ascentis C18 MeOH 10%: OK")
