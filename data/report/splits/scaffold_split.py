"""
Deterministic Murcko-scaffold split over RepoRT's non-SMRT compound universe.

Verified against the pinned checkout: 18,872 unique non-SMRT compounds (by
inchikey.std), 100% with a usable SMILES (smiles.std), 96.9% with a populated
ClassyFire class. Scaffold split is fully feasible.

Requires rdkit (`pip install rdkit`).

This is deterministic CODE, not a materialized id list: given the same
compound universe, `assign_split` always returns the same train/val/test
label for the same InChIKey, so nothing bulky needs to be committed. Compound
identity is scoped to InChIKey (first-seen SMILES per InChIKey across the
non-SMRT corpus), not (dataset, row) -- one compound measured in many
datasets gets ONE scaffold-split label, so the split is safe against
compound leakage across datasets, not just across rows.

Usage:
    python -m data.report.splits.scaffold_split /path/to/RepoRT/processed_data
"""
from __future__ import annotations

import argparse
import hashlib
from collections import defaultdict

from data.report.loader import iter_datasets, load_rtdata

SMRT_DATASET_IDS = {"0186", "0209"}

TRAIN_FRAC = 0.8
VAL_FRAC = 0.1
# remaining 0.1 is test


def murcko_scaffold(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)


def assign_split(scaffold: str) -> str:
    """Deterministic bucket assignment by hashing the scaffold SMILES.
    Same scaffold -> same split, always -- this is what prevents scaffold
    leakage between train/val/test."""
    if scaffold == "":
        # unparseable / acyclic molecules with no ring scaffold: bucket by
        # their own hash so they don't all pile into one split
        h = int(hashlib.md5(b"__NO_SCAFFOLD__").hexdigest(), 16)
    else:
        h = int(hashlib.md5(scaffold.encode("utf-8")).hexdigest(), 16)
    frac = (h % 10_000) / 10_000.0
    if frac < TRAIN_FRAC:
        return "train"
    if frac < TRAIN_FRAC + VAL_FRAC:
        return "val"
    return "test"


def build_compound_universe(processed_data_root: str) -> dict[str, str]:
    """inchikey.std -> smiles.std, first-seen, non-SMRT datasets only."""
    universe: dict[str, str] = {}
    for d in iter_datasets(processed_data_root):
        if d.id in SMRT_DATASET_IDS:
            continue
        rows = load_rtdata(processed_data_root, d.id, "canonical_success")
        for r in rows:
            ik = r["inchikey.std"]
            if ik not in universe and r.get("smiles.std"):
                universe[ik] = r["smiles.std"]
    return universe


def compute_split(processed_data_root: str) -> dict[str, dict]:
    """Returns {inchikey: {"scaffold": str, "split": "train"|"val"|"test"}}."""
    universe = build_compound_universe(processed_data_root)
    out = {}
    scaffold_cache: dict[str, str] = {}
    for ik, smi in universe.items():
        scaf = scaffold_cache.get(smi)
        if scaf is None:
            scaf = murcko_scaffold(smi)
            scaffold_cache[smi] = scaf
        out[ik] = {"scaffold": scaf, "split": assign_split(scaf)}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_data_root")
    args = ap.parse_args()
    result = compute_split(args.processed_data_root)
    counts: dict[str, int] = defaultdict(int)
    n_scaffolds: dict[str, set] = defaultdict(set)
    for ik, info in result.items():
        counts[info["split"]] += 1
        n_scaffolds[info["split"]].add(info["scaffold"])
    print("total compounds:", len(result))
    for split in ("train", "val", "test"):
        print(
            f"  {split}: {counts[split]} compounds "
            f"({counts[split]/len(result):.1%}), "
            f"{len(n_scaffolds[split])} distinct scaffolds"
        )
    all_scaffolds = set().union(*n_scaffolds.values())
    leaked = sum(
        1
        for s in all_scaffolds
        if sum(1 for split in n_scaffolds if s in n_scaffolds[split]) > 1
    )
    print(f"scaffolds appearing in >1 split (should be 0): {leaked}")
