#!/usr/bin/env python3
"""Column registry: `column_id` -> geometry, particle class, phase descriptors.

`spec/scenario-ensemble.md` says `column_id` "resolves to a column vector;
identity, not a name". This is that resolution, and it is the thing the sampler
was hardcoding.

**Design decision: the registry declares only what cannot be looked up.** An
entry carries its identity, vendor, geometry and pH rating; the HSM phase
descriptors (`H, S*, A, B, C(2.8), C(7.0)`, silica `type`) are read from
`sources/hsm-column-db/` and the LSER availability from `sources/wsu-lser/`.
Re-typing 819 columns' parameters into a Python file would create a second
source of truth that silently drifts from the first.

Adding a phase is therefore three lines - id, geometry, pH range - and
everything else resolves. `unregistered()` lists what is available but not yet
declared, so the registry can grow without anyone reading the CSV by hand.

numpy + stdlib.
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field

from dispersion import KNOX_PRIORS, ColumnGeometry, plate_count

HERE = os.path.dirname(os.path.abspath(__file__))
HSM_CSV = os.path.join(HERE, "..", "sources", "hsm-column-db", "database.csv")
WSU_CSV = os.path.join(HERE, "..", "sources", "wsu-lser",
                       "wsu2019-system-constants.csv")


@dataclass(frozen=True)
class ColumnRecord:
    column_id: str
    hsm_name: str                 # key into the HSM database
    vendor: str
    geometry: ColumnGeometry
    ph_min: float
    ph_max: float
    ph_source: str                # where the rating came from - never assume it
    role: str = ""
    hsm: dict = field(default_factory=dict)
    wsu_lser: bool = False

    @property
    def silanol_c70(self):
        return self.hsm.get("C70")

    @property
    def hydrophobicity_h(self):
        return self.hsm.get("H")

    def plates(self, flow_ml_min, rng=None, size=None):
        return plate_count(self.geometry, flow_ml_min, rng, size)

    def supports_ph(self, ph: float) -> bool:
        return self.ph_min <= ph <= self.ph_max


# --- declarations. Geometry is the 150 x 4.6 mm, 5 um format #30 settled on,
# --- so the arms stay comparable on run time and back-pressure.
_150x46_5um_fp = ColumnGeometry(15.0, 0.46, 5e-4, "fully_porous")
_150x46_5um_cs = ColumnGeometry(15.0, 0.46, 5e-4, "core_shell")

_DECLARED = [
    dict(column_id="kinetex-evo-c18", hsm_name="Kinetex EVO C18",
         vendor="Phenomenex", geometry=_150x46_5um_cs,
         ph_min=1.0, ph_max=12.0, ph_source="Phenomenex product literature",
         role="#30 low-silanol arm; organo-silica hybrid, keeps tier 2b alive"),
    dict(column_id="supelcosil-lc-18", hsm_name="Supelcosil LC-18",
         vendor="Sigma-Aldrich / Supelco", geometry=_150x46_5um_fp,
         ph_min=2.0, ph_max=7.5, ph_source="Sigma-Aldrich product page (owner-read)",
         role="#30 high-silanol arm; type-A silica, the regime WSU-2019 lacks"),
    dict(column_id="xbridge-c18", hsm_name="XBridge C18",
         vendor="Waters", geometry=_150x46_5um_fp,
         ph_min=1.0, ph_max=12.0, ph_source="⚠ UNVERIFIED - vendor page not read",
         role="BEH hybrid low-silanol baseline; IN WSU-2019"),
    dict(column_id="ascentis-rp-amide", hsm_name="Ascentis RP-Amide",
         vendor="Sigma-Aldrich / Supelco", geometry=_150x46_5um_fp,
         ph_min=2.0, ph_max=9.0, ph_source="⚠ UNVERIFIED - vendor page not read",
         role="polar-embedded, suppressed silanol"),
    dict(column_id="zorbax-bonus-rp", hsm_name="Zorbax Bonus RP",
         vendor="Agilent", geometry=_150x46_5um_fp,
         ph_min=2.0, ph_max=9.0, ph_source="⚠ UNVERIFIED - vendor page not read",
         role="polar-embedded; the LOWEST C(7.0) in the database"),
]


def _load():
    hsm = {}
    with open(HSM_CSV, newline="", encoding="latin-1") as f:
        for r in csv.DictReader(f):
            hsm.setdefault(r["name"], r)
    with open(WSU_CSV, newline="", encoding="utf-8-sig") as f:
        wsu = {r["column"] for r in csv.DictReader(f)}

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    out = {}
    for d in _DECLARED:
        h = hsm.get(d["hsm_name"])
        desc = {} if h is None else {
            k: num(h[src]) for k, src in
            (("H", "H"), ("S", "S"), ("A", "A"), ("B", "B"),
             ("C28", "C28"), ("C70", "C70"))}
        if h is not None:
            desc["type"] = h["type"]
        out[d["column_id"]] = ColumnRecord(
            hsm=desc, wsu_lser=d["hsm_name"] in wsu,
            **{k: v for k, v in d.items() if k != "hsm_name"},
            hsm_name=d["hsm_name"])
    return out


REGISTRY = _load()


def get(column_id: str) -> ColumnRecord:
    if column_id not in REGISTRY:
        raise KeyError(f"{column_id!r} not registered. Known: {sorted(REGISTRY)}")
    return REGISTRY[column_id]


def unregistered(pattern: str, limit=12):
    """HSM columns matching `pattern` that are not yet declared here."""
    known = {r.hsm_name for r in REGISTRY.values()}
    out = []
    with open(HSM_CSV, newline="", encoding="latin-1") as f:
        for r in csv.DictReader(f):
            if pattern.lower() in r["name"].lower() and r["name"] not in known:
                out.append((r["name"], r["C70"], r["H"], r["type"]))
    return out[:limit]


def main():
    flow = 1.0
    print("Column registry\n")
    print(f"{'column_id':<20}{'C(7.0)':>8}{'H':>7}{'type':>7}{'LSER':>6}"
          f"{'pH':>10}{'N @1mL':>9}  role")
    for cid, r in REGISTRY.items():
        print(f"{cid:<20}{r.silanol_c70:>8.2f}{r.hydrophobicity_h:>7.2f}"
              f"{r.hsm.get('type',''):>7}{'yes' if r.wsu_lser else '-':>6}"
              f"{f'{r.ph_min}-{r.ph_max}':>10}{r.plates(flow):>9.0f}  {r.role[:38]}")

    print("\nSilanol span if all five were run:")
    cs = sorted((r.silanol_c70, cid) for cid, r in REGISTRY.items())
    print(f"  {cs[0][1]} ({cs[0][0]:+.2f})  ->  {cs[-1][1]} ({cs[-1][0]:+.2f})"
          f"   span {cs[-1][0]-cs[0][0]:.2f}")
    print(f"  #30's booked pair spans "
          f"{get('supelcosil-lc-18').silanol_c70 - get('kinetex-evo-c18').silanol_c70:.2f}")

    print("\n⚠ Hydrophobicity is NOT matched across the wider set:")
    for cid, r in sorted(REGISTRY.items(), key=lambda kv: kv[1].hydrophobicity_h):
        print(f"    {cid:<20} H = {r.hydrophobicity_h:.2f}")
    print("  #30's pair matches H to 0.002 so the contrast is silanol and not")
    print("  retention strength. Adding phases BUYS C-span and SPENDS that control:")
    print("  a C-vs-D relationship fitted across columns differing in both is")
    print("  confounded, and no amount of instrument time separates them after.")

    print("\n⚠ pH ratings needing verification before use:")
    for cid, r in REGISTRY.items():
        if "UNVERIFIED" in r.ph_source:
            print(f"    {cid:<20} {r.ph_min}-{r.ph_max}  ({r.ph_source})")

    print("\nAvailable but not declared (examples):")
    for pat in ("Amide", "Bonus"):
        for nm, c70, h, t in unregistered(pat, 3):
            print(f"    {nm:<34} C70={c70:>7} H={h:>6} type={t}")


if __name__ == "__main__":
    main()
