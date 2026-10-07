"""
Loader for RepoRT (https://github.com/michaelwitting/RepoRT), pinned.

RepoRT itself is NOT vendored in this repo. Clone it yourself and point
REPORT_ROOT at the `processed_data` directory of a checkout pinned to the
commit recorded in PINNED_COMMIT:

    git clone https://github.com/michaelwitting/RepoRT.git /path/to/RepoRT
    cd /path/to/RepoRT && git checkout 9de8d603377bbb6cc0f74e250eb7533fd8874df1

Then:

    from data.report.loader import iter_datasets
    for ds in iter_datasets("/path/to/RepoRT/processed_data"):
        ...

Everything below is derived from directly reading the pinned checkout, not
from RepoRT's papers (which are stale — see findings.md).
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from typing import Iterator, Optional

# ---------------------------------------------------------------------------
# Pin
# ---------------------------------------------------------------------------

PINNED_COMMIT = "9de8d603377bbb6cc0f74e250eb7533fd8874df1"
PINNED_COMMIT_SHORT = "9de8d60"
PINNED_COMMIT_DATE = "2026-07-09"  # commit author date, from `git log -1 --date=short`
PINNED_COMMIT_SUBJECT = (
    "Merge pull request #245 from michaelwitting/inconsistencies-0052 "
    "(Duplicate rows in dataset 0052)"
)
REPO_URL = "https://github.com/michaelwitting/RepoRT.git"

# Re-verified 2026-08-14: 9de8d60 is still the tip of RepoRT's default branch
# (no drift since #4 pinned it on 2026-07-09). One commit to pin, one date to
# record — both above.

# ---------------------------------------------------------------------------
# The name-swap bug (datasets 0310-0341)
# ---------------------------------------------------------------------------
#
# Verified directly against the pinned checkout (see verify_name_swap.py and
# findings.md for the full evidence table). Extent, precisely:
#
#   - Datasets 0310-0325 (16 datasets): info.tsv `name` field says
#     "..._Raptor Biphenyl" but metadata.tsv `column.name` is actually
#     "Waters CORTECS UPLC C18" (USP L1, particle size 1.6 um -- matches the
#     real CORTECS spec).
#   - Datasets 0326-0341 (16 datasets): info.tsv `name` field says
#     "..._CORTECS C18" but metadata.tsv `column.name` is actually
#     "Restek Raptor Biphenyl" (USP L11, particle size 2.7 um -- matches the
#     real Raptor Biphenyl spec).
#   - Total: 32 datasets affected, ids 0310-0341 inclusive. This matches the
#     range named in issue #29 exactly -- no drift, no fuzz at the edges.
#   - ONLY the info.tsv `name` string is swapped (specifically the trailing
#     column-label suffix baked into that string, e.g. "MSMS_pos_30_ACN_<col>").
#     `column.name`, `column.usp.code`, `column.particle.size` and every other
#     metadata.tsv field are internally self-consistent and correct -- USP
#     code and particle size agree with `column.name`, not with the `name`
#     string, on both sides of the swap.
#   - The bug is confined to exactly these 32 datasets: the BGC family
#     immediately before (0236-0259, columns HSS T3 / XB-C18 / BEH C18) and
#     immediately after (0342-0357, HSS C18) show a consistent name<->column
#     name mapping in the same spot-check methodology.
#
# CORRECTIVE ACTION: do not trust info.tsv `name` for column identity in this
# range. Use metadata.tsv `column.name` (this loader does that automatically
# via `patched_dataset_name`).

SWAPPED_NAME_IDS = {f"{i:04d}" for i in range(310, 342)}  # 0310..0341 inclusive

_SWAP_MAP = {
    # info.tsv name substring -> the metadata.tsv column.name it should have
    # carried (kept here only as documentation; the fix itself just means
    # "trust column.name, not the name string", see patched_dataset_name()).
    "Raptor Biphenyl": "Waters CORTECS UPLC C18",
    "CORTECS C18": "Restek Raptor Biphenyl",
}


def patched_dataset_name(dataset_id: str, raw_name: str, column_name: str) -> str:
    """Return the info.tsv `name` field, corrected for the 0310-0341 swap.

    For affected ids we do not attempt to reconstruct the "true" free-text
    dataset name (we don't know it) -- we replace the misleading column-label
    suffix with the trustworthy one from metadata.tsv `column.name`, and flag
    it. Anyone who joins on dataset `name` and expects it to describe the
    column used should use this function, not the raw field.
    """
    if dataset_id not in SWAPPED_NAME_IDS:
        return raw_name
    for wrong_substr in _SWAP_MAP:
        if wrong_substr in raw_name:
            return raw_name.replace(wrong_substr, column_name) + " [name-swap-corrected]"
    return raw_name + " [name-swap-suspected-but-pattern-not-matched]"


# ---------------------------------------------------------------------------
# Dataset iteration
# ---------------------------------------------------------------------------


@dataclass
class Dataset:
    id: str
    raw_name: str
    name: str  # patched, see patched_dataset_name
    column_name: str
    column_usp_code: str
    column_particle_size: str
    column_length: str
    column_id_mm: str
    column_temperature: str  # '' if missing
    eluent_a_ph: str
    eluent_b_ph: str
    column_flowrate: str
    eluent_a_acn: float
    eluent_a_meoh: float
    eluent_b_acn: float
    eluent_b_meoh: float
    source: str
    authors: str
    missing_information: str
    dir: str
    name_swapped: bool = False

    @property
    def modifier(self) -> str:
        """ACN / MeOH / ACN+MeOH / other-or-none, from eluent composition
        (NOT from the name string, which is unreliable in the swap range)."""
        acn = self.eluent_a_acn > 0 or self.eluent_b_acn > 0
        meoh = self.eluent_a_meoh > 0 or self.eluent_b_meoh > 0
        if acn and meoh:
            return "ACN+MeOH"
        if acn:
            return "ACN"
        if meoh:
            return "MeOH"
        return "other/none"

    @property
    def ph(self) -> Optional[float]:
        for raw in (self.eluent_a_ph, self.eluent_b_ph):
            try:
                v = float(raw)
                if v > 0:
                    return v
            except (TypeError, ValueError):
                continue
        return None


def _f(row: dict, key: str) -> float:
    v = row.get(key)
    try:
        return float(v) if v not in (None, "") else 0.0
    except ValueError:
        return 0.0


def _read_first_row(path: str) -> dict:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return row
    return {}


def iter_datasets(processed_data_root: str) -> Iterator[Dataset]:
    """Yield one Dataset per id under `processed_data/`, info+metadata joined
    and the name-swap bug patched."""
    for did in sorted(os.listdir(processed_data_root)):
        if not did.isdigit():
            continue
        ddir = os.path.join(processed_data_root, did)
        info_p = os.path.join(ddir, f"{did}_info.tsv")
        meta_p = os.path.join(ddir, f"{did}_metadata.tsv")
        if not (os.path.exists(info_p) and os.path.exists(meta_p)):
            continue
        info = _read_first_row(info_p)
        meta = _read_first_row(meta_p)
        raw_name = info.get("name", "")
        column_name = meta.get("column.name", "")
        patched = patched_dataset_name(did, raw_name, column_name)
        yield Dataset(
            id=did,
            raw_name=raw_name,
            name=patched,
            column_name=column_name,
            column_usp_code=meta.get("column.usp.code", ""),
            column_particle_size=meta.get("column.particle.size", ""),
            column_length=meta.get("column.length", ""),
            column_id_mm=meta.get("column.id", ""),
            column_temperature=meta.get("column.temperature", ""),
            eluent_a_ph=meta.get("eluent.A.pH", ""),
            eluent_b_ph=meta.get("eluent.B.pH", ""),
            column_flowrate=meta.get("column.flowrate", ""),
            eluent_a_acn=_f(meta, "eluent.A.acn"),
            eluent_a_meoh=_f(meta, "eluent.A.meoh"),
            eluent_b_acn=_f(meta, "eluent.B.acn"),
            eluent_b_meoh=_f(meta, "eluent.B.meoh"),
            source=info.get("source", ""),
            authors=info.get("authors", ""),
            missing_information=info.get("missing information", ""),
            dir=ddir,
            name_swapped=did in SWAPPED_NAME_IDS,
        )


def load_rtdata(processed_data_root: str, dataset_id: str,
                 variant: str = "canonical_success") -> list[dict]:
    """Load the retention-time rows for a dataset.
    variant: 'canonical_success' (default), 'isomeric_success', 'isomeric_failed'.
    """
    path = os.path.join(
        processed_data_root, dataset_id, f"{dataset_id}_rtdata_{variant}.tsv"
    )
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
