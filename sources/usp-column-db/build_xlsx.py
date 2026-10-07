#!/usr/bin/env python3
"""Build the USP Column Equivalency workbook from the extracted CSVs."""
import csv, sys, os
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

BASE = os.path.dirname(os.path.abspath(__file__))
PARAMS = os.path.join(BASE, "usp-approach-columns.csv")
SIM    = os.path.join(BASE, "usp-similarity-long.csv")
OUT    = os.path.join(BASE, "USP-column-equivalency.xlsx")

PARAM_COLS = ["Hy", "CTF", "CFA", "TFA", "BD"]

for path in (PARAMS, SIM):
    if not os.path.exists(path):
        sys.exit(f"{path} not found. The USP data is not redistributed; "
                 "see README.md in this folder for how to obtain it.")

params = {r["name"]: r for r in csv.DictReader(open(PARAMS, encoding="utf-8"))}
sim = list(csv.DictReader(open(SIM, encoding="utf-8")))

wb = Workbook()
HDR = Font(bold=True, color="FFFFFF")
FILL = PatternFill("solid", fgColor="2F5597")
def style_header(ws, ncols):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=1, column=c)
        cell.font = HDR; cell.fill = FILL
        cell.alignment = Alignment(horizontal="center", vertical="center")
    ws.freeze_panes = "A2"

def autosize(ws, maxw=42):
    for col in ws.columns:
        letter = get_column_letter(col[0].column)
        width = max((len(str(c.value)) for c in col if c.value is not None), default=8)
        ws.column_dimensions[letter].width = min(max(width + 2, 9), maxw)

def num(v):
    if v is None or str(v).strip() == "":
        return None
    try:
        return float(v)
    except ValueError:
        return v

# ---- Sheet 1: Columns -------------------------------------------------------
ws = wb.active; ws.title = "Columns"
hdr = ["Column", "Hy", "CTF", "CFA", "TFA", "BD", "USP Designation",
       "Manufacturer", "Complete (all 5)"]
ws.append(hdr)
for name in sorted(params):
    r = params[name]
    vals = [num(r[k]) for k in PARAM_COLS]
    complete = "yes" if all(v is not None for v in vals) else "no"
    ws.append([name] + vals + [r["usp_designation"], r["manufacturer"], complete])
style_header(ws, len(hdr)); autosize(ws)

# ---- Sheet 2: Similarity (long) --------------------------------------------
ws2 = wb.create_sheet("Similarity (long)")
hdr2 = ["Reference column", "Rank", "F", "Column", "Hy", "CTF", "CFA", "TFA",
        "BD", "USP Designation", "Manufacturer"]
ws2.append(hdr2)
for r in sim:
    tgt = r["target"]
    p = params.get(tgt, {})
    ws2.append([r["reference"], num(r["rank"]), num(r["F"]), tgt]
               + [num(p.get(k, "")) for k in PARAM_COLS]
               + [p.get("usp_designation", ""), p.get("manufacturer", "")])
style_header(ws2, len(hdr2)); autosize(ws2)

# ---- Sheet 3: F matrix ------------------------------------------------------
refs = sorted({r["reference"] for r in sim})
tgts = sorted({r["target"] for r in sim})
lookup = {}
for r in sim:
    lookup[(r["reference"], r["target"])] = num(r["F"])
ws3 = wb.create_sheet("F matrix")
ws3.append(["Reference \\ Target"] + tgts)
for ref in refs:
    ws3.append([ref] + [lookup.get((ref, t)) for t in tgts])
style_header(ws3, len(tgts) + 1)
ws3.column_dimensions["A"].width = 34
ws3.freeze_panes = "B2"

# ---- Sheet 4: README --------------------------------------------------------
ws4 = wb.create_sheet("README")
ws4.column_dimensions["A"].width = 120
notes = [
    ("USP Column Equivalency database — extracted 2026-08-15", True),
    ("", False),
    ("SOURCE", True),
    ("https://apps.usp.org/app/USPNF/columnsDB.html  (USP Column Equivalency Application, 'USP Database' half).", False),
    ("Built by the USP Working Group on Column Equivalency using NIST SRM 870.", False),
    ("Method review: Pharmacopeial Forum 31(2), 637-645 (2005).", False),
    ("", False),
    ("WHAT THESE PARAMETERS ARE", True),
    ("SRM 870 is a column PERFORMANCE/ACTIVITY test mixture. These parameters describe column", False),
    ("activity (silanol interaction, metal chelation, bonding density) — NOT retention selectivity.", False),
    ("They are therefore NOT the Snyder/PQRI hydrophobic-subtraction set (H, S*, A, B, C) and are", False),
    ("NOT interchangeable with it, nor with LSER system constants.", False),
    ("", False),
    ("!! Parameter expansions are UNCONFIRMED. BD is almost certainly bonding density (0.9-5.5 fits", False),
    ("   umol/m2) and TFA/CTF are tailing-factor-like, but USP does not define the abbreviations", False),
    ("   publicly and the Pharmacopeial Forum article was not obtained. Do not treat the expansions", False),
    ("   as fact without reading PF 31(2) 637-645.", False),
    ("", False),
    ("DATA QUALITY CAVEATS", True),
    ("!! 'USP Designation' reads L1 for ALL 145 columns, including physically C8 and phenyl phases", False),
    ("   (Acclaim 120 C8, TSKgel Octyl-80Ts, TSKgel Super-Octyl, TSKgel Super-Phenyl). This was", False),
    ("   verified against the site's own rendered table: it is the SOURCE's data, not an extraction", False),
    ("   artefact. DO NOT use that field as column chemistry.", False),
    ("", False),
    ("!! Only 124 of 145 columns have all five parameters. The 21 incomplete columns are excluded", False),
    ("   from USP's own F similarity ranking, which is why the site reports c_total = 124.", False),
    ("   Missing counts: CTF 12, CFA 7, TFA 5, BD 5, Hy 1.", False),
    ("", False),
    ("F VALUES", True),
    ("F is USP's column-difference measure: smaller = more similar. F = 0 is the column against", False),
    ("itself (rank 0). The 'Similarity (long)' sheet holds the full ranked list for every reference", False),
    ("column; 'F matrix' is the same data pivoted.", False),
    ("", False),
    ("VERIFIED: F is exactly symmetric. Across all 7,501 pairs where both directions were", False),
    ("captured, max |F(a,b) - F(b,a)| = 0.0. This is a strong check that the extraction is exact.", False),
    ("", False),
    ("!! TWO COLUMNS CANNOT BE USED AS A REFERENCE (defect in USP's own application):", True),
    ("   TSKgel ODS-100V 3 um and TSKgel ODS-100V 5 um. Selecting either as the reference column", False),
    ("   returns an EMPTY rank-0 row and a ranking computed against a null parameter vector, giving", False),
    ("   silently wrong (inflated) F values. Reproduced through the site's own UI, not an artefact", False),
    ("   of this extraction. Their outgoing rankings are therefore OMITTED from this workbook.", False),
    ("   Their own parameters and their INCOMING F values (as targets in other columns' lists) are", False),
    ("   valid and retained. This is how the asymmetry was caught: those two were the only", False),
    ("   asymmetric references in the whole matrix.", False),
    ("", False),
    ("So: 143 reference columns x 145 target columns, 17,891 ranked rows.", False),
    ("", False),
    ("LICENCE — READ BEFORE SHARING", True),
    ("(c) The United States Pharmacopeial Convention. Accessed under the application's terms of use.", False),
    ("NO redistribution licence has been identified. Treat as INTERNAL DEVELOPMENT AND EDUCATIONAL", False),
    ("USE ONLY. Do not publish or ship this workbook. If the data becomes load-bearing, cite the", False),
    ("Pharmacopeial Forum article as the primary source rather than this extract.", False),
]
for text, bold in notes:
    ws4.append([text])
    if bold:
        ws4.cell(row=ws4.max_row, column=1).font = Font(bold=True)

wb.save(OUT)
print("wrote", OUT)
print("Columns sheet rows:", len(params))
print("Similarity rows:", len(sim))
print("F matrix:", len(refs), "x", len(tgts))
