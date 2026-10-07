# Thin slice — findings (#21)

Session 2026-08-15, branch `prototype/thin-slice`. Throwaway vertical slice:
SMILES → SoluteML (all 25 ensemble members, machinery reused from
`prototype/descriptor-covariance`) → ionisation (hard-coded literature w_w pKa;
provider pluggable by design, the slice tests composition not the provider) →
LSER → `k(φ, pH)` → scenario ensemble `retention[S, n]`, S = 1000.

Column: **XBridge Shield RP18, acetonitrile** — the only shortlisted column with
no incomplete-wetting row in ACN (support φ = 0.10–0.70 fully unmasked, #24) and
only a "weak electrostatic" Table S-1 caveat. 45 °C dataset. Isocratic operating
points (0.30, pH 3.0), (0.30, pH 7.0), (0.50, pH 5.0); nominal t0 = 1.0 min;
dwell/gradient physics out of slice scope. All retention math in log₁₀ k.

Compounds: n = 16 = 6 WSU-94 neutral anchors (benzamide, caffeine, acetophenone,
methylparaben, toluene, naphthalene) + 4 acids (benzoic, ibuprofen, naproxen,
ketoprofen) + 4 bases (lidocaine, propranolol, atenolol, imipramine) + 1
amphoteric (4-aminobenzoic acid) + 1 **structurally refused** quaternary cation
(benzyltrimethylammonium) to exercise NaN-for-refused. pKa citations in
`compounds.py`; ketoprofen's own source spread (3.98 Avdeef vs 4.45 Clarke's) is
0.47 — the ±0.35 per-compound pKa error is not pessimistic.

## Verdict first: the distributions ARE usable as identity priors

At every operating point, 105 pairs over the 15 valid compounds. **Regenerated
2026-08-16** after the Table S-2 extract fix (`class` column; see
EXECUTION-PLAN.md) — the anchor RMSE that sets the inflation λ is unchanged at
n = 94, since this script already stripped the stray section heading, but the
ensemble was re-drawn and two counts moved by one pair:

| (φ, pH) | confident P>0.9 | ambiguous 0.6–0.9 | coin-flip <0.6 |
|---|---|---|---|
| 0.30, 3.0 | 77 (73%) | 24 | 4 |
| 0.30, 7.0 | 73 (70%) | 28 | 4 |
| 0.50, 5.0 | 70 (67%) | 31 | 4 |

Everything does **not** overlap everything: median rank IQR is 2 positions of
15, worst compound 4–5 positions (imipramine at φ = 0.50 — the CONTEXT
"degraded" tier made visible: fully ionised, so its distribution is
k_neutral × 10^(−D) with D at prior width; its t_R 5–95% spans 1.5–32 min at
φ = 0.30). Retention orderings are chromatographically sane (atenolol at t0,
naphthalene last, acids collapse toward t0 from pH 3 → 7).

Coin-flip pairs are chemically honest ones: caffeine/propranolol,
lidocaine/propranolol at pH 3 (all near t0), and **naproxen/ketoprofen twice**
(P = 0.41–0.42). Note the profen pair is exactly where #37's unmodelled
scaffold correlation says the slice *overstates* Var(Δlog k) by up to ~2× — the
true order confidence for that pair is plausibly better than reported. The
slice draws per-compound errors independently (the set is deliberately
diverse); a production sampler needs the per-family correlation.

**Cold-baseline comparison (Crippen logP, RDKit), per #27 — no accuracy claim
is possible because no measured retention exists for this set** (known next
step: #4/RepoRT data; none found trivially, none hunted). What is reportable:
agreement of the model's median order with the logP order is 76–88% on the
confident pairs (n = 70–77 per op point) and 37–54% on ambiguous + coin-flip
pairs (n = 28–35). The model's confident calls largely reproduce the cold
baseline where the cold baseline is expected to be right (neutrals, monotone
lipophilicity); its claimed added value sits exactly on the pairs where they
disagree — unfalsifiable until real data enters. That is the correct place for
a prior to disagree, and the correct thing to test first with #4 data.

## Variance attribution (the central deliverable)

Freeze one layer at a time (descriptors at ensemble mean; column vector at
published; D at 1.5; pKa at literature); share = 1 − Var_frozen/Var_total of
log k, median over compounds in class-at-that-pH (class by nominal neutral
fraction: >0.9 neutral, 0.1–0.9 partial, <0.1 ionised). Shares don't sum to 1
exactly (interactions + Monte-Carlo noise at S = 1000; small negatives are
noise).

| (φ, pH) | class (n) | sd(log k) med | desc | column | pKa | D |
|---|---|---|---|---|---|---|
| 0.30, 3.0 | neutral (10) | 0.159 | **0.99** | 0.00 | 0.00 | −0.01 |
| | partial (1) | 0.269 | 0.93 | 0.01 | 0.04 | 0.02 |
| | ionised (4) | 0.554 | 0.14 | 0.00 | 0.01 | **0.84** |
| 0.30, 7.0 | neutral (6) | 0.150 | **1.00** | 0.00 | 0.00 | −0.01 |
| | partial (1) | 0.352 | 0.34 | 0.00 | **0.52** | 0.15 |
| | ionised (8) | 0.476 | 0.20 | 0.00 | −0.09 | **0.78** |
| 0.50, 5.0 | neutral (6) | 0.102 | **1.00** | 0.00 | 0.00 | −0.01 |
| | partial (5) | 0.309 | 0.19 | −0.00 | **0.67** | 0.04 |
| | ionised (4) | 0.518 | 0.07 | 0.00 | −0.01 | **0.92** |

Three regimes, each owned by one layer:

- **Neutrals: QSPR descriptors are the entire budget** (sd(log k) 0.10–0.16).
  The realised QSPR-induced log k error on the 6 anchors (predicted-mean vs
  measured Table S-2 descriptors through the nominal column vector) is **RMS
  0.098 at φ = 0.30, 0.064 at φ = 0.50 (n = 6)** — squarely inside #6's
  budgeted SE(log k) ≈ 0.08–0.11. The sampled sd runs slightly wider because
  the ensemble covariance is inflated to the full WSU-94 RMSE (which includes
  bias, see honesty list).
- **Partially ionised (pH within the pKa danger band): pKa owns it**
  (share 0.52–0.67, sd(log k) 0.27–0.35) — quantitatively consistent with
  CONTEXT's ≈ 0.36-from-pKa danger-band figure.
- **Fully ionised: D owns it** (share 0.78–0.92, sd(log k) 0.48–0.55) — and D
  is a *prior, not a fit*, so these widths are the prior speaking, exactly the
  scouting-run refinement target #34 designed for.
- **Column vector ≈ 0.00 everywhere.** With the mandatory joint Σ = D·R·D draw
  (#37 — diagonal draws forbidden), published-constant error contributes
  sd ≈ 0.015 and vanishes next to everything else. A diagonal draw would have
  reported it 16× larger; this slice would have mis-attributed the budget.

**Effort concentration implication (recommendation, not decision): D and the
pKa layer (including the unacquired s_s scale correction, #20/#22) are where
the effort belongs. The LSER coefficient machinery can stay crude.**

## e·E on/off (#6 experiment 1)

| (φ, pH) | median \|Δlog k\| | max \|Δlog k\| | max \|ΔP(order)\| | classes with → without |
|---|---|---|---|---|
| 0.30, 3.0 | 0.083 | 0.126 | 0.128 | 77/24/4 → 78/22/5 |
| 0.30, 7.0 | 0.083 | 0.127 | 0.067 | 73/28/4 → 75/26/4 |
| 0.50, 5.0 | 0.039 | 0.066 | 0.062 | 70/31/4 → 70/31/4 |

Dropping e·E is **not free in absolute retention at low φ** — Shield's e is
0.071 at 30% ACN (the "e ≈ 0.03 in ACN" collapse is a ≥ 40–50% φ statement;
at φ = 0.30 every compound here has E ≥ 0.6, so the term is a near-common
+0.04…+0.13 shift). But precisely because it is near-common, it **cancels in
order**: max |ΔP| ≤ 0.12, confidence-class counts move by ≤ 2 pairs of 105.
Consistent with #36's ruling: keep the term, propagate E's error; its cost is
absolute-retention bias at low φ, never order. No change recommended.

## A = 0 categorical check (#6 experiment 2)

Ensemble-mean A against structural donor truth, threshold 0.01: **10/15
correct** — far below #6's 94.3%, but the misses are one-sided and mostly
tiny: all five are false A≠0 on true-A=0 compounds — caffeine 0.110,
acetophenone 0.023, naphthalene 0.022, toluene 0.021, imipramine 0.011.
Retention consequence a·A at φ = 0.30: **−0.035 log k for caffeine, ≤ 0.007
for the other four** — none flips behaviour qualitatively, and the dangerous
direction (missed real donor, A truly ≠ 0 called 0) did not occur on this set
(0 of 9 donors missed). The 0.01 threshold is too tight for ensemble means: a
0.03–0.05 dead zone calls 14/15 here, leaving caffeine as the one genuine
categorical error worth watching (it is also the compound furthest from
solvation-training chemistry — a fully N-substituted xanthine).

## Interface verdict (#14): survives, with a concrete punch list

1. **`retention[S, n]` + NaN-for-refused held.** Benzyltrimethylammonium
   (refused: no neutral microspecies at any pH, structural test) kept its
   column as NaN through every analysis; n stayed 16; nothing downstream
   special-cased it beyond NaN-awareness. As specified by #33.
2. **The operating-point axis is missing from #14.** Multiple (φ, pH) points
   forced `retention[P, S, n]` plus per-point coordinates and t0 — i.e. the
   method card must travel *with* the array, not beside it. #14 should say so.
3. **Microspecies weights are produced and then discarded.** Per-scenario
   f_neutral exists inside the sampler and is needed downstream twice (D
   refinement from scouting runs, #34; conformal stratum assignment, #14) but
   the interface has no slot for it. Recommend adding per-compound,
   per-op-point stratum (or f_neutral summary) to the ensemble object.
4. **`width[S, n]` is fillable but empty of physics** — filled with a nominal
   N = 10000 placeholder here. Interface fine; content needs the instrument
   layer.
5. **ESS: untested, not verified.** No measurement entered, so no reweighting
   happened. Nothing contradicts the ESS design; the slice cannot confirm it.
6. The three-tier per-run/per-compound sampler structure (#37/#14) implemented
   in ~60 lines with no friction — the correlation structure is cheap.

## Honesty list

- **No real retention data was compared.** Every "usable" claim above is a
  claim about *width and separation of the prior*, not accuracy. Next step
  unchanged: #4/RepoRT data through this exact code.
- **In-sample-adjacent calibration.** Descriptor error magnitude is tied to
  WSU-94 RMSE (E 0.154, S 0.186, A 0.080, B 0.063, n = 94) via per-target
  inflation λ = RMSE/median-ensemble-sd = (3.95, 3.70, 3.85, 2.15) — ensemble
  spread is uncalibrated dispersion, not calibrated confidence, and WSU-94 is
  training-adjacent. Inflation converts spread to *RMSE-matched* spread, which
  bakes bias in as variance; real out-of-sample error for drug-like compounds
  is unknown (the wide drug spreads — naproxen S ± 0.24, imipramine E ± 0.23 —
  are the ensemble itself signalling extrapolation).
- **w_w pKa on a s_s world.** The Rosés/Bosch correction (#20/#22) is
  unacquired; carried as a named per-run scale term (acids N(+0.3, 0.30),
  bases N(−0.15, 0.20), literature-typical for 30–60% ACN, not fitted) plus
  ±0.35 per-compound literature error.
- **D ~ N(1.5, 0.5) truncated > 0.1 is a prior, not a fit.** All ionised-class
  widths above are prior width. That is by design (#34) but must not be read
  as a validated ionisable-regime claim — the two-regime rule stands.
- **LSER licensed for neutral microspecies of calibration-like solutes.** The
  drug neutral forms are larger and more functionalised than the 94-compound
  pool; their LSER lack-of-fit is unbudgeted here (form-correction territory,
  #36). Neutral-anchor numbers are the trustworthy core.
- **Scaffold correlation not modelled** (independent per-compound draws);
  overstates Var(Δlog k) for related pairs (profens) by up to ~2× (#37).
- 45 °C throughout; Shield's weak-electrostatic caveat unmodelled (absorbed by
  D in spirit, unquantified); methylparaben's phenol (pKa ≈ 8.4) treated as
  neutral — ≤ 11% ionised at the pH 7.5 grid edge, and the grid points used
  stop at 7.0.

## Reproduction

`predict_drugs.py` (soluteml venv + extracted egg, see
`prototype/descriptor-covariance/FINDINGS.md` "Reproduction") wrote
`drug_preds_raw.csv` (25 members × 9 drugs) — committed, so `slice.py` reruns
without the 871 MB download. `slice.py` (venv python: numpy, scipy, rdkit;
seed 20260815) writes `results.json` and `retention_ensemble.npz`
(the #14 interface object, 3 × 1000 × 16). Anchor predictions reused from
`prototype/descriptor-covariance/soluteml_preds_raw.csv`. McGowan V computed
from structure (V is not a QSPR target).
