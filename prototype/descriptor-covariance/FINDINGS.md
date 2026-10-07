# Descriptor covariance cheap check — findings (#37)

Session 2026-08-15, branch `prototype/descriptor-covariance`. Three measurement
threads; two produced numbers, one is blocked on data the repo does not have.
All variances in `log₁₀ k` units (WSU convention). Wetting-excluded per #40
throughout (332 of 354 fits kept).

## What was NOT computable, stated first

**True per-compound LSER residuals do not exist in this repo.** The WSU extract
carries fitted coefficients, per-coefficient SDs, and the *aggregate* fit SE —
no raw `log k` values. So measurement 3's core question (is a compound's
residual correlated across φ on the same column, i.e. a fixed per-compound
offset = a #36 form correction) **needs raw retention data**: HSM3 (#26) or the
sources Poole 2019 cites. That is a data-acquisition finding, not a failure of
the question. The across-compound LSER-residual correlation for scaffold pairs
is blocked by the same absence.

**The 6×6 coefficient covariance is not published.** Tables S-3/S-4/S-5 give
per-coefficient SDs only. Everything under "Measurement 3 (pivoted)" therefore
approximates the coefficient covariance as `Σ = D R D`, with `D` the published
SDs and `R` the correlation of `(X'X)⁻¹` from the 94-compound Table S-2 design
`[1,E,S,A,B,V]`. `R` is exact up to the unknown per-fit compound subset (n =
40–56 of the 94); the approximation self-validates: it reproduces the OLS
identity `sd(pred) ≈ SE·√(p/n)` (median propagated sd 0.015 vs SE·√(6/n) ≈
0.013).

## Measurement 3 (pivoted): LSER coefficient-error covariance

`measure3_coeff_covariance.py`. Coefficient-estimator correlations from the
design are large and mostly negative: corr(c,v) = −0.65, corr(e,s) = −0.66,
corr(s,b) = −0.64, corr(b,v) = −0.48.

- **Var(log k), full vs diagonal, per compound** (n = 31,208 = 332 fits × 94
  compounds): ratio median **0.061**, IQR [0.039, 0.091], 5–95% [0.022, 0.166].
  Treating the six system constants as independent draws **overstates the
  coefficient-error variance ~16× (median)** — sd overstated ~4×. Absolute
  scale: sd_diag median 0.059, correct (full) sd median ~0.015.
- **Var(Δlog k), close pairs** (|predicted Δlog k| < 0.1, Ascentis C18 50%):
  MeOH ratio median **0.564** (n = 362 pairs), ACN median **0.571** (n = 445);
  IQR ~[0.36, 0.84], 5–95% up to 1.33 — a minority of pairs *amplify*.
  Diagonal sampling overstates close-pair sd(Δlog k) ~1.3× at the median. The
  intercept error cancels in pairs (per-run class); the `e..v` part does not,
  because Δ-descriptors ≠ 0.
- **Bonus, model-implied LSS lack-of-fit**: fitting the LSER-predicted
  `log k(φ)` linear-in-φ per (compound, column) gives RMS curvature median
  **0.027** MeOH (n = 2,350), **0.069** ACN (n = 2,256), 95th pct 0.090/0.155.
  This is the curvature of the LSER surface itself — a floor for #2's LSS
  form-error term that exists even with perfect descriptors and coefficients.

Implication: the sampler's per-run draw of `(c,e,s,a,b,v)` **must be joint**.
`Σ = D R D` costs one 6×6 matrix per fit and is already implemented here.

## Measurement 1: SoluteML within-compound (E,S,A,B) ensemble covariance

`predict_soluteml.py` + `measure1_ensemble_covariance.py`. All 25 members
(5 folds × 5 models) ran on all 94 WSU compounds — see "Reproduction" for the
two traps. SMILES: all 94 resolved (`wsu94_smiles.csv`), 0 ambiguous enough to
skip; every one passes a McGowan-V check against Table S-2 to 0.005
(`check_smiles_mcgowan.py`); one flagged isomer note (2-bromoacetophenone,
α- vs ring-, same formula so V cannot disambiguate).

Sanity vs Table S-2 (n = 94): bias/RMSE E +0.035/0.154, S +0.061/0.186,
A −0.004/0.080, B +0.018/0.063 — consistent with SoluteML's published E-worst
ordering. Median within-compound ensemble sd: E 0.039, S 0.050, A 0.021,
B 0.029 (epistemic-only; ~3–4× smaller than the actual RMSE, the usual
in-distribution underestimate).

- Mean within-compound ensemble correlation: E–S **+0.29**, S–A +0.21,
  E–B +0.19, rest |ρ| < 0.06. Sign of E–S matches the Ulrich-consistency
  expectation in the issue thread.
- **Var(log k) full/diag ratio through Ascentis C18 50%** (n = 94 each):
  ACN median **0.962**, IQR [0.81, 1.17]; MeOH median **0.890**, IQR
  [0.68, 1.08]; 5–95% spans ~[0.45, 1.31]. The mixed-sign coefficients mostly
  cancel the positive off-diagonals: within-compound descriptor covariance
  moves Var(log k) by **typically < ±35%, direction compound-specific**.

## Across-compound scaffold correlation (the material one)

Because Table S-2 *is* the reference lineage, SoluteML's ensemble-mean
residuals vs S-2 are real per-compound descriptor errors, and propagating them
through the column vector gives a per-compound `log k` error (sd 0.117 MeOH,
0.066 ACN, n = 94). Within name-based scaffold families (n = 5 families, 6–10
members each), the mean pairwise error product is positive for 4 of 5 families
in ACN and 3 of 5 in MeOH; relative to the overall error variance the implied
intra-family correlation is **ρ ≈ +0.42 (MeOH), +0.54 (ACN)** — but strongly
family-dependent (phenols ρ high, anilines ≈ 0, in MeOH the alkylphenone
series ≈ 0).

Positive ρ at this size means independent per-compound draws **overstate
Var(Δlog k) for scaffold pairs by up to ~2×** (`2σ²(1−ρ)` vs `2σ²`), which is
exactly the mechanism #14 flagged: fictitious order uncertainty, every
P(A before B) dragged toward 0.5. Caveat: the WSU 94 are simple, mostly
training-adjacent solutes; this gets the sign and shape, not the out-of-sample
magnitude. RepoRT's scaffold split (#29) remains the proper instrument —
measurement 2 proper is **not done**.

## Recommendation (recommendation, not measurement)

1. **Joint draw of the six system constants is mandatory** (per-run class).
   Diagonal sampling inflates coefficient-error Var(log k) ~16× and close-pair
   Var(Δlog k) ~1.8×. Use `Σ = D R D` until #35's re-fit yields the exact
   covariance from its own design matrix — at which point store it; it is free.
2. **Across-compound joint drawing is justified** by the positive, sizeable,
   family-dependent scaffold correlation (ρ ≈ +0.4–0.5 propagated to log k).
   Since it is family-dependent, a per-family (or descriptor-space-kernel)
   correlation, not one global ρ.
3. **Within-compound 4×4 descriptor covariance is a second-order effect**
   (ratio 0.89–0.96 median) — carry it anyway, it falls out of the ensemble
   for free, but do not spend design effort on it.
4. Acquire raw `log k` (HSM3 #26 or Poole's cited sources) before believing
   any per-compound residual/φ-correlation claim; it is not derivable from the
   published fits.

## Reproduction

- `measure3_coeff_covariance.py`, `measure1_ensemble_covariance.py`: stdlib +
  numpy, run with system python.
- SoluteML: **pip/git install of fhvermei/chemprop_solvation does NOT include
  the trained models** — `final_models/` in the GitHub repo is an empty
  placeholder. The models live in the conda package
  `fhvermei/chemprop_solvation 0.0.3` (871 MB tarball, downloadable without
  conda from api.anaconda.org), which contains
  `chemprop_solvation-0.0.2-py3.7.egg` with
  `final_models/SoluteML/fold_{0..4}/model_{0..4}/model.pt`.
- The 2019 code runs on Python 3.10 + torch 2.13 with two shims (both in
  `predict_soluteml.py`): `torch.load(weights_only=False)` for the Namespace
  pickles, and a degenerate-feature guard — 8 of the 200
  `rdkit_2d_normalized` features were near-constant in training (worst `Ipc`,
  std 6.7e-14) and modern descriptastorus (2.8.0) values differ at ~1e-6,
  which the stored feature scaler amplifies to ~1e7 and saturates the FFN.
  Zero those 8 scaled features and clip the rest to ±10. Without the guard,
  predictions are garbage at ~1e5 magnitude; with it, they reproduce published
  accuracy. venv recipe in the script docstring.
- `soluteml_preds_raw.csv` (25 × 94 × E,S,A,B,L) is committed so the analysis
  reruns without the 871 MB download.
