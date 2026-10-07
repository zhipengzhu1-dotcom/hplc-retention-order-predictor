#!/usr/bin/env python3
"""Extract SoluteDB (SoluteML's training corpus) from the Zenodo xlsx to CSV.

Source: Zenodo 5792296, `Solvation_data-1.0.0/selected_data_for_this_work/
SoluteDB_selected_data.xlsx` -- the CC-BY corpus the Chung 2022 models were
trained on (neutral solutes, H C N O S P F Cl Br I only).

Reads xlsx with the standard library only (it is a zip of XML); openpyxl is not
installed in either project venv and this is a one-shot extraction.

Usage: python read_solutedb.py <path-to-SoluteDB_selected_data.xlsx> <out.csv>
"""

import csv
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _cell_ref_to_col(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_sheet(path: str, sheet_index: int = 0) -> list[list[str]]:
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{NS}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))
        names = sorted(n for n in z.namelist()
                       if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", n))
        root = ET.fromstring(z.read(names[sheet_index]))
        rows = []
        for row in root.iter(f"{NS}row"):
            cells: dict[int, str] = {}
            for c in row.findall(f"{NS}c"):
                v = c.find(f"{NS}v")
                if v is None and c.find(f"{NS}is") is None:
                    continue
                if c.get("t") == "s":
                    val = shared[int(v.text)]
                elif c.get("t") == "inlineStr":
                    val = "".join(t.text or "" for t in c.iter(f"{NS}t"))
                else:
                    val = v.text or ""
                cells[_cell_ref_to_col(c.get("r"))] = val
            if cells:
                width = max(cells) + 1
                rows.append([cells.get(i, "") for i in range(width)])
    return rows


def main() -> None:
    rows = read_sheet(sys.argv[1])
    header = rows[0]
    width = len(header)
    with open(sys.argv[2], "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows[1:]:
            w.writerow((r + [""] * width)[:width])
    print(f"{len(rows) - 1} rows, columns: {header}")


if __name__ == "__main__":
    main()
