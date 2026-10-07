# QSPR environment — SoluteML with exposed ensemble covariance

Resolves #38. Provider chosen
in #6; the covariance is
consumed by #14 and the
per-member predictions by #37.

## Where SoluteML actually lives

SoluteML is not a standalone package. It is one of three trained models inside
**`fhvermei/chemprop_solvation`** — the solvation-specific fork of chemprop that accompanies
Chung, Vermeire, Wu, Walker, Abraham & Green, *J. Chem. Inf. Model.* **2022**, 62, 433–446
([doi:10.1021/acs.jcim.1c01103](https://doi.org/10.1021/acs.jcim.1c01103)). The other two are
`DirectML_Gsolv` and `DirectML_Hsolv`. `chemprop_solvation/solvation_estimator.py` is the file
#6 pointed at; `load_SoluteML_estimator()` is a thin wrapper over `load_ML_estimator`.

The repository ships **no weights**. The pretrained 25-member ensemble (5 folds × 5 models) is
a separate 1.9 GB Zenodo deposit.

## Exact pins

| Component | Pin |
|---|---|
| Source | `github.com/fhvermei/chemprop_solvation` @ `eb7f7d90fd064d9d939d27ba5837485af9a75295` (master tip, 2022-03-02 — the repo has had no commit since) |
| Weights | Zenodo record [5792296](https://zenodo.org/records/5792296) **v1.0.0**, `ML_model_files.zip`, 1,988,111,831 bytes |
| Feature generator | `descriptastorus` @ `9a190343bcd3cfd35142d378d952613bcac40797` (not on PyPI) |
| Python / platform | 3.13.7, macOS arm64 (Darwin 24.5.0) |
| Everything else | [`requirements-qspr.txt`](./requirements-qspr.txt) |
| Licence | CC-BY-4.0 (package); the chemprop fork's `LICENSE.txt` is MIT |

## How the patch is applied — vendored patch files

Three `git apply`-able patches in [`patches/`](./patches/), applied to the pinned upstream commit
by [`install_qspr.sh`](./install_qspr.sh):

- `0001-torch-weights-only-compat.patch` — required just to load the checkpoints at all.
- `0002-guard-degenerate-feature-std.patch` — required for the model to produce numbers rather
  than nonsense on a modern RDKit.
- `0003-expose-soluteml-ensemble-covariance.patch` — the actual ticket.

Only the third is the feature. The first two are the price of running 2022 code in 2026, and
both were found by the model failing loudly rather than by inspection — see friction below.

**Why vendored patches rather than a fork or an upstream PR.** Upstream is dormant: last commit
2022-03-02, no releases, no open PR activity — an upstream PR has no realistic path to being
merged and would leave us pinned to our own branch anyway. A fork would work, but it hides the
delta: with a patch file the entire modification to third-party code is 113 lines of reviewable
diff sitting in this repo, and the upstream commit it applies to is stated explicitly. Anyone can
see exactly what we changed without diffing two repositories. The cost — patches must be rebased
if upstream ever moves — is near zero for a repository that has not moved in four years.

## What the patch exposes

`estimator(smiles, return_ensemble=True)` returns two extra values after the existing
`(avg_pre, epi_unc, valid_indices)`:

- `member_preds` — `np.ndarray (25, n_compounds, 4)`, every individual ensemble member's
  prediction. #37 needs per-member residuals to measure across-compound scaffold correlation;
  summary statistics cannot produce that.
- `cov` — `np.ndarray (n_compounds, 4, 4)`, the ensemble covariance over **(E, S, A, B)**.

`ddof=0`, matching the `np.var` upstream already uses, so `cov[i].diagonal()` reproduces the
shipped `epi_unc[i]` exactly — the new number is the same tensor, not a re-derivation.

**Four targets, not five.** The SoluteML head emits `(E, S, A, B, L)`. The exposed covariance is
over the first four only. `V` is computed exactly by McGowan's method and carries zero model
error; `L` appears only in the GC form of the LSER equation and has no role in RPLC. The module
exports both `SOLUTE_TARGETS` (upstream's five) and `RPLC_SOLUTE_TARGETS` (our four) so the
distinction cannot be lost by accident. The raw upstream means/variances still carry `L` in
position 4 — we do not rewrite upstream's output, we just do not build anything on it.

## The number is not calibrated, and this ticket did not calibrate it

`epi_unc` and `cov` are **ensemble spread**: epistemic only, no aleatoric term, no calibration.
#6 measured the shipped σ for `E` at ≈0.015 against a real WSU disagreement of SE = 0.252 —
understated by more than 10×. #6 also found that ensemble spread *shrinks* off-distribution,
which is exactly where drug-like ionisable compounds sit. Recalibration is #14's stratified
conformal layer. This ticket exposes the number and nothing else; the docstring in the patched
estimator says so too.

## Setup

```
bash env/install_qspr.sh
env/qspr-venv/bin/python env/soluteml_demo.py
```

A separate venv from `env/ionisation-venv` (#20) on purpose — see the header of
`requirements-qspr.txt`. Same interpreter, same `rdkit`/`torch`/`numpy` pins, so it is the same
platform rather than a competing one; the split just stops a break in 2022-era code from taking
the pKa provider down with it.

## Install friction

1. **`torch.load` refuses the checkpoints on modern torch.** Since torch 2.6 `torch.load`
   defaults to `weights_only=True`, and these checkpoints pickle an `argparse.Namespace` plus
   scaler objects. Every one of the nine `torch.load` calls in `utils.py` fails with
   `_pickle.UnpicklingError` until `weights_only=False` is passed. That is what patch 0001 does,
   and it is the reason the "~5-line patch" from #6 is not the whole story — without it the
   package cannot load a model at all on the pinned environment.
2. **A near-zero training standard deviation turns the model into a random number generator.**
   The features scaler stores `std = 6.7e-14` for the CDF-normalized RDKit `Ipc` descriptor
   (index 33 of 200), which saturates at 1.0 for essentially every molecule. Modern RDKit
   returns `Ipc` differing from the 2021 value in about the sixth decimal place; divided by
   `6.7e-14` that becomes an input of **1.5e7**, and the predictions come back as
   `E ≈ −21140` — identical for every compound, because one feature has swamped the network.
   Nothing warns you: every checkpoint tensor loads cleanly and the numbers are simply wrong.
   `StandardScaler.fit` guards `std == 0` exactly, which does not catch `6.7e-14`.
   Patch 0002 replaces degenerate stds with 1, restoring the training-time behaviour (that
   column contributes ~0). **Anyone reviving a pretrained cheminformatics model of this era
   should check for this class of bug before trusting any output.**
3. **No weights in the repo.** `chemprop_solvation/final_models/` contains only
   `empty_place_holder.txt`. The 1.9 GB Zenodo bundle must be downloaded and the `SoluteML`
   directory copied into place; SoluteML is 1.5 GB of it, the two DirectML models we do not
   use are 253 MB each.
4. **The deposit ships more models than the API loads.** `SoluteML/` contains **10 folds × 10
   models = 100 checkpoints**, but `load_SoluteML_estimator()` hard-codes
   `fold_count=5, model_count=5`. We keep the 25 the paper describes. If #14 ever wants a
   larger ensemble it is already on disk — but those extra members were not part of the
   published ensemble and their fold assignment is undocumented, so treat that as a change to
   the model, not a free sample-size increase.
5. **The published sample output does not reproduce exactly.** Deltas of 0.017–0.054 on
   (E, S, A, B), i.e. 1.3–3.6× the shipped ensemble σ. Cause is RDKit/descriptastorus drift.
   Full numbers and their consequences in
   [`SOLUTEML_COVARIANCE_RESULTS.md`](./SOLUTEML_COVARIANCE_RESULTS.md).
6. **`descriptastorus` is git-only.** SoluteML was trained with the `rdkit_2d_normalized`
   feature generator, so it is a hard runtime dependency, and it is not published to PyPI. It
   builds cleanly from source on Python 3.13 / arm64 and pulls in `pandas-flavor` and `xarray`.
7. **The conda route in the README is dead for this platform.** `conda install -c fhvermei
   chemprop_solvation` targets the 2022 dependency set; installing from source into the existing
   pinned venv was the working path.
8. **`batch_size=1` is hard-coded** in the estimator's `predict` call, so cost is
   25 models × n compounds forward passes with no batching. Fine for a probe set, worth knowing
   before anyone runs it over RepoRT.

## Contents

- [`requirements-qspr.txt`](./requirements-qspr.txt) — exact lockfile.
- [`install_qspr.sh`](./install_qspr.sh) — clone at pin, apply patches, fetch weights.
- [`patches/`](./patches/) — the three vendored patches.
- [`soluteml_demo.py`](./soluteml_demo.py) — runnable demonstration.
- [`SOLUTEML_COVARIANCE_RESULTS.md`](./SOLUTEML_COVARIANCE_RESULTS.md) — its output and reading.
