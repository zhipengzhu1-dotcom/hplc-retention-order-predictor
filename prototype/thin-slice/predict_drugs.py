#!/usr/bin/env python3
"""Run all 25 SoluteML ensemble members on the thin-slice drug compounds.

Reuses the machinery recovered in prototype/descriptor-covariance (see its
FINDINGS.md "Reproduction" for the two shims: weights_only=False and the
degenerate rdkit_2d_normalized feature guard). Anchor compounds are NOT re-run:
their 25-member predictions are already committed in
prototype/descriptor-covariance/soluteml_preds_raw.csv.

Usage: <soluteml-venv-python> predict_drugs.py <path-to-extracted-egg-dir>
Writes drug_preds_raw.csv: one row per (compound, fold, model), columns E,S,A,B,L.
"""

import csv
import os
import sys
from copy import deepcopy

import numpy as np
import torch

_orig_load = torch.load
torch.load = lambda *a, **k: _orig_load(*a, **{**k, "weights_only": False})

EGG = os.path.abspath(sys.argv[1])
sys.path.insert(0, EGG)

from chemprop_solvation.data import MoleculeDatapoint, MoleculeDataset  # noqa: E402
from chemprop_solvation.train.predict import predict  # noqa: E402
from chemprop_solvation.utils import load_args, load_checkpoint, load_scalers  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from compounds import COMPOUNDS, DRUGS  # noqa: E402

MODEL_DIR = os.path.join(EGG, "chemprop_solvation", "final_models", "SoluteML")
TARGETS = ["E", "S", "A", "B", "L"]


def main():
    todo = [c for c in COMPOUNDS if c["name"] in DRUGS]
    names = [c["name"] for c in todo]
    smiles = [[c["smiles"]] for c in todo]

    first = os.path.join(MODEL_DIR, "fold_0", "model_0", "model.pt")
    train_args = load_args(first)
    assert not train_args.solvation
    train_args.features_generator = ["rdkit_2d_normalized"]

    data = MoleculeDataset([MoleculeDatapoint(sm, train_args) for sm in smiles])
    invalid = [names[i] for i in range(len(data)) if data[i].mol is None]
    assert not invalid, invalid

    out = open(os.path.join(HERE, "drug_preds_raw.csv"), "w", newline="")
    w = csv.writer(out)
    w.writerow(["compound", "fold", "model"] + TARGETS)
    for fold in range(5):
        pt0 = os.path.join(MODEL_DIR, f"fold_{fold}", "model_0", "model.pt")
        scaler, features_scaler = load_scalers(pt0)
        args = load_args(pt0)
        data_fold = deepcopy(data)
        if args.features_scaling:
            data_fold.normalize_features(features_scaler)
            degen = features_scaler.stds < 1e-3
            for i in range(len(data_fold)):
                f = np.asarray(data_fold[i].features, dtype=float)
                f[degen] = 0.0
                data_fold[i].features = np.clip(f, -10.0, 10.0)
        for m in range(5):
            model = load_checkpoint(
                os.path.join(MODEL_DIR, f"fold_{fold}", f"model_{m}", "model.pt"),
                cuda=False,
            )
            preds = np.array(predict(model=model, data=data_fold, batch_size=32, scaler=scaler))
            for i, name in enumerate(names):
                w.writerow([name, fold, m] + [f"{v:.5f}" for v in preds[i]])
            print(f"fold {fold} model {m} done", flush=True)
    out.close()


if __name__ == "__main__":
    main()
