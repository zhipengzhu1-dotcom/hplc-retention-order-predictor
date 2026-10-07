#!/usr/bin/env python3
"""Run all 25 SoluteML ensemble members (5 folds x 5 models) on the 94 WSU compounds.

Environment (not committed; rebuild as below):
  uv venv --python 3.10 soluteml-venv
  uv pip install --python soluteml-venv/bin/python torch rdkit scikit-learn \
      "numpy<2" tqdm scipy tensorboardX descriptastorus
Models: conda tarball fhvermei/chemprop_solvation 0.0.3 (871 MB) from anaconda.org,
which contains chemprop_solvation-0.0.2-py3.7.egg with
final_models/SoluteML/fold_{0..4}/model_{0..4}/model.pt. The GitHub repo ships an
empty final_models placeholder -- pip install from GitHub does NOT get the models.

Usage: python predict_soluteml.py <path-to-extracted-egg-dir>
Writes soluteml_preds_raw.csv: one row per (compound, member), columns E,S,A,B,L.
Keeps per-member predictions -- the ensemble covariance is computed downstream.
"""

import csv
import os
import sys
from copy import deepcopy

import numpy as np
import torch

# torch>=2.6 defaults weights_only=True; the 2019 checkpoints pickle argparse
# Namespaces, so restore the old behaviour for these trusted local files.
_orig_load = torch.load
torch.load = lambda *a, **k: _orig_load(*a, **{**k, "weights_only": False})

EGG = os.path.abspath(sys.argv[1])
sys.path.insert(0, EGG)

from chemprop_solvation.data import MoleculeDatapoint, MoleculeDataset  # noqa: E402
from chemprop_solvation.train.predict import predict  # noqa: E402
from chemprop_solvation.utils import load_args, load_checkpoint, load_scalers  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(EGG, "chemprop_solvation", "final_models", "SoluteML")
TARGETS = ["E", "S", "A", "B", "L"]


def main():
    with open(os.path.join(HERE, "wsu94_smiles.csv")) as f:
        rows = list(csv.DictReader(f))
    names = [r["compound"] for r in rows]
    smiles = [[r["smiles"]] for r in rows]

    first = os.path.join(MODEL_DIR, "fold_0", "model_0", "model.pt")
    train_args = load_args(first)
    assert not train_args.solvation
    train_args.features_generator = ["rdkit_2d_normalized"]

    data = MoleculeDataset([MoleculeDatapoint(sm, train_args) for sm in smiles])
    invalid = [names[i] for i in range(len(data)) if data[i].mol is None]
    assert not invalid, invalid

    out = open(os.path.join(HERE, "soluteml_preds_raw.csv"), "w", newline="")
    w = csv.writer(out)
    w.writerow(["compound", "fold", "model"] + TARGETS)
    for fold in range(5):
        scaler, features_scaler = load_scalers(
            os.path.join(MODEL_DIR, f"fold_{fold}", "model_0", "model.pt")
        )
        args = load_args(os.path.join(MODEL_DIR, f"fold_{fold}", "model_0", "model.pt"))
        data_fold = deepcopy(data)
        if args.features_scaling:
            data_fold.normalize_features(features_scaler)
            # descriptastorus version drift: 8 of the 200 rdkit_2d_normalized
            # features were near-constant in the 2019 training set (std < 1e-3,
            # worst 'Ipc' at 6.7e-14); the modern generator's values differ at
            # ~1e-6, which the tiny stds blow up to ~1e7 and the FFN saturates.
            # Zero the degenerate features (== pin at training mean) and clip
            # the rest to a sane z-range.
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
            preds = predict(model=model, data=data_fold, batch_size=32, scaler=scaler)
            preds = np.array(preds)  # (94, 5)
            for i, name in enumerate(names):
                w.writerow([name, fold, m] + [f"{v:.5f}" for v in preds[i]])
            print(f"fold {fold} model {m} done", flush=True)
    out.close()


if __name__ == "__main__":
    main()
