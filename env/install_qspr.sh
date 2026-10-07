#!/usr/bin/env bash
# Recreate the QSPR (SoluteML) environment for #38 from scratch.
#
#   bash env/install_qspr.sh
#
# Creates env/qspr-venv, clones chemprop_solvation at the pinned upstream commit, applies
# the two vendored patches in env/patches/, and downloads the 1.9 GB pretrained model
# ensemble from Zenodo. Nothing under env/_src, env/_download or env/qspr-venv is committed.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
cd "$ROOT"

# Pinned upstream source: fhvermei/chemprop_solvation, master tip, 2022-03-02.
# This is the reference implementation of SoluteML from Chung et al., JCIM 2022,
# doi:10.1021/acs.jcim.1c01103. CC-BY-4.0 (repo LICENSE.txt says MIT for the chemprop
# fork; README states CC-BY-4.0 for the package).
UPSTREAM_URL="https://github.com/fhvermei/chemprop_solvation.git"
UPSTREAM_SHA="eb7f7d90fd064d9d939d27ba5837485af9a75295"

# Pinned pretrained weights: Zenodo record 5792296 v1.0.0, ML_model_files.zip,
# 1,988,111,831 bytes. Datasets for Chung et al. 2022.
ZENODO_URL="https://zenodo.org/api/records/5792296/files/ML_model_files.zip/content"
ZENODO_SHA256="b23f2021b4575d19d635a7e592f47ea0c78bb5b52f750fa517de3d67273c6eaf"

VENV="env/qspr-venv"
SRC="env/_src/chemprop_solvation"

echo "==> 1/5 venv"
python3 -m venv "$VENV"
"$VENV/bin/pip" install --upgrade pip

echo "==> 2/5 pinned dependencies"
"$VENV/bin/pip" install -r env/requirements-qspr.txt

echo "==> 3/5 upstream source at pinned commit"
rm -rf "$SRC"
git clone "$UPSTREAM_URL" "$SRC"
git -C "$SRC" checkout "$UPSTREAM_SHA"

echo "==> 4/5 apply vendored patches"
for p in env/patches/*.patch; do
  echo "    $p"
  git -C "$SRC" apply "$ROOT/$p"
done
"$VENV/bin/pip" install -e "$SRC"

echo "==> 5/5 pretrained ensemble from Zenodo (1.9 GB)"
mkdir -p env/_download
if [ ! -f env/_download/ML_model_files.zip ]; then
  curl -L -o env/_download/ML_model_files.zip "$ZENODO_URL"
fi
echo "$ZENODO_SHA256  env/_download/ML_model_files.zip" | shasum -a 256 -c -
unzip -q -o env/_download/ML_model_files.zip -d env/_download/
mkdir -p "$SRC/chemprop_solvation/final_models"
cp -R env/_download/ML_model_files/SoluteML "$SRC/chemprop_solvation/final_models/"

echo
echo "Done. Run the demonstration with:"
echo "  $VENV/bin/python env/soluteml_demo.py"
