# Which QSPR approach predicts Abraham solute descriptors best?

Research resolution for issue #6.
Layer under test: `SMILES → molecular descriptors → QSPR → E, S, A, B` (+ V, computed).
Downstream consumers: #14 (uncertainty representation and propagation), #21 (end-to-end thin slice).
Incorporates the ground-truth provenance findings from #5.

Every claim below is cited. A section at the end lists what could **not** be determined.

> **Read §8 before quoting any accuracy figure in this document.** The reference data that every
> published score is measured against is substantially estimated or computed rather than measured
> (#5). Accuracy numbers appear here because the ticket asks for them, not because they arbitrate the
> decision. They do not.

---

## 0. Answer in one paragraph

**SoluteML** — the multi-task directed message-passing neural network of Chung et al. — is the
recommended QSPR layer, *not* because it has the best published accuracy (it does not; the 2026 UFZ
graph-convolutional model does) but because it is the only candidate that is simultaneously openly
licensed with **usable pretrained weights**, **multi-task** (so the descriptors come out of one shared
representation and their error covariance is recoverable), and **independently validated on a real
reversed-phase system**. Its shipped uncertainty is a 25-member ensemble variance — epistemic only,
per-descriptor, uncalibrated — which is a starting point, not a finished product. Four facts
constrain the whole architecture: (1) **the layer has four prediction targets, E, S, A and B**;
V is computed from structure and L belongs to the GC form of the equation, so neither is predicted;
(2) **no predicted descriptor set reaches experimental-descriptor quality for RPLC** — budget
SE(log k) ≈ 0.08–0.11 from this layer alone versus 0.03–0.06 with measured descriptors;
(3) **the published accuracy scores are partly a measure of how well each model mimics ACD/Labs'
predictor**, because a large fraction of the reference values were computed or estimated rather than
measured (§8); (4) **for our in-scope chemistry every candidate is extrapolating**, and the honest
clean-provenance, drug-like, ionisable validation set is **twelve compounds** (#5). Accuracy cannot
decide this choice. Architecture and uncertainty can.

---

## 1. What is actually being predicted (and what is not)

The Abraham solvation parameter model
([Poole 2021, tutorial](https://doi.org/10.1016/j.chroma.2021.462108)):

```
log SP = c + e·E + s·S + a·A + b·B + v·V      (condensed–condensed, e.g. RPLC)
log SP = c + e·E + s·S + a·A + b·B + l·L      (gas–condensed, e.g. GC)
```

**V is not a prediction target.** McGowan's characteristic volume is computed exactly from the
molecular formula: sum the atomic contributions and subtract 6.56 cm³/mol per bond, regardless of
bond order ([Abraham & McGowan, *Chromatographia* **23**, 243–246, 1987](https://doi.org/10.1007/BF02311772)).
Atapattu & Temerdashev state this plainly: "Except for the solute descriptor McGowan's characteristic
volume, V, the remaining solute descriptors E, S, A, B, and L of the solvation parameter model are
experimentally determined"
([J. Chromatogr. Open **7**, 100213, 2025](https://doi.org/10.1016/j.jcoa.2025.100213)).

### Spec-level rule: there are FOUR prediction targets

**The QSPR layer predicts E, S, A and B. Nothing else.** Record this before anyone builds five models.

- **V is computed, not predicted.** Zero model error. AbraLlama reports RMSE = 0.097 for V
  ([Liquids **4**, 518–524, 2024](https://doi.org/10.3390/liquids4030029), Table 1); that number is a
  pure artefact of predicting a quantity that is already exact, and it signals a model built without
  domain scrutiny. Any comparison that includes V is inflated.
- **L is out of scope.** L is the hexadecane–air partition coefficient and appears only in the
  gas-phase form of the equation. RPLC uses `v·V`, not `l·L`. L figures are retained below *only*
  where they illuminate a method's behaviour — it is the descriptor every model predicts worst, and
  where the 2026 GNN's margin is largest — and are **excluded from every ranking**.
- The 2026 UFZ head-to-head predicts S, E, A, B, L and does not attempt V
  ([ACS Environ. Au, 2026](https://doi.org/10.1021/acsenvironau.6c00063)) — correct on V, but one of
  its five models is irrelevant to us.

A real simplification: **20% fewer models to train, version and calibrate**, and it removes the
descriptor with the largest published error (L RMSE 0.25–0.45) from our error budget entirely.

### "Experimental" descriptors are not measurements

Except V, Abraham descriptors are back-solved by regression from sets of measured partition
coefficients and chromatographic retention data. E derives from molar refraction; S, A, B and L are
regression-derived. The derivation route — which solvent systems, which laboratory, which QC — is
itself a variance source, and **no compilation reports uncertainties on its descriptor values at all**
(#5). This is the root of both the WSU-vs-UFZ disagreement (§3.4) and the contamination problem (§8).

---

## 2. The methods compared

| # | Method | Architecture | Source |
|---|---|---|---|
| 1 | **GCA** (group contribution) | Ertl functional-group detection, groups as fitted linear parameters + aliphatic/aromatic backbone terms; oversized groups decomposed | [ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063) |
| 2 | **kNN** | Jaccard distance over atom-centred-fragment (ACF) count vectors; several ACF orders as sub-models, combined by weighted sum; after Kühne et al. | [ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063) |
| 3 | **GNN** (graph-convolutional) | DeepChem 2.7.2 GraphConv, TF/Keras 2.15, 2 hidden layers (32, 64), ReLU, L1 loss, dropout 0.1, lr 1e-3, batch 50; **one model per descriptor** | [ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063) |
| 4 | **SoluteGC** | Atom-centred functional groups + ring-strain corrections + long-distance-interaction groups | [JCIM **62**, 433–446, 2022](https://doi.org/10.1021/acs.jcim.1c01103) |
| 5 | **SoluteML** | Directed message-passing NN (chemprop) + 200 RDKit 2D descriptors, **multi-task (E,S,A,B,L in one model)**, 5 folds × 5 models = 25-member ensemble | [JCIM 2022](https://doi.org/10.1021/acs.jcim.1c01103); code: [chemprop_solvation](https://github.com/fhvermei/chemprop_solvation/blob/master/chemprop_solvation/solvation_estimator.py) |
| 6 | **DNNmono / DNNtaut** | The 2022 UFZ predecessor to (3); DeepChem GraphConv, both single-task and multi-task variants | [Fluid Phase Equil. **555**, 113349, 2022](https://doi.org/10.1016/j.fluid.2021.113349) |
| 7 | **AbraLlama-Solute** | 30M-parameter ChemLLaMA fine-tuned 20 epochs; predicts E,S,A,B,**V** jointly | [Liquids **4**, 518–524, 2024](https://doi.org/10.3390/liquids4030029) |
| 8 | **ACD/Absolv**, **LSERD QSPR** | Commercial / database-embedded QSPR baselines | reported in [ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063) Table 2 & §3.6 |

---

## 3. Accuracy, per descriptor

> **Validated against:** every figure in §3.1–3.3 is measured against the Absolv/UFZ-LSER descriptor
> compilation, which #5 has established is **substantially estimated or computed rather than
> measured** — 16.7% of E values and 22.1% of B values are output from ACD Absolv, ChemSketch or
> ClogP, and roughly 48% of E and S values are "estimated by comparison to closely related
> compounds". These scores therefore measure agreement with a commercial predictor and with expert
> analogy, not with reality. §3.4 and §3.5, which validate against the WSU experimental set and
> against real retention data, are the trustworthy numbers.

### 3.1 UFZ independent test set (RMSE, dimensionless)

10% random hold-out per descriptor; the same test compounds across all three models, and the same
test set as the 2022 study, so these are directly comparable *to each other*
([ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063), Table 1 and §3.6).
GCA and GNN figures are consensus over 5 splits.
**Validated against: Absolv/UFZ-LSER — mixed measured / estimated / computed. See §8.**

| Descriptor | Exp. range | **GCA** | **kNN** | **GNN** | DNNtaut (2022) | LSERD QSPR | ACD/Absolv |
|---|---|---|---|---|---|---|---|
| **E** (excess molar refraction) | −1.18–4.62 | 0.188 | 0.131 | **0.086** | 0.12 | 0.14 | 0.10 |
| **S** (dipolarity/polarisability) | −0.89–4.80 | 0.234 | 0.216 | **0.169** | 0.22 | 0.28 | 0.23 |
| **A** (HB acidity) | 0.00–2.50 | 0.119 | 0.126 | **0.075** | 0.11 | 0.09 | 0.09 |
| **B** (HB basicity) | 0.00–4.26 | 0.142 | 0.151 | **0.095** | 0.14 | 0.13 | 0.16 |
| ~~L~~ *(GC only — not a target)* | −0.82–27.44 | 0.420 | 0.454 | 0.246 | 0.42 | 0.52 | 0.44 |

Across the four RPLC targets the GNN wins every one, with margins of 22–54%. Before consensus
averaging, single GNN models scored E 0.100, S 0.187, A 0.086, B 0.106 — the consensus buys ~10–20%.

**But note what the E column means.** ACD/Absolv scores 0.10 on E here, essentially matching the best
neural network. Given that ~1 in 6 E reference values *is* ACD output (§8), a high score on E is at
least partly a measure of Absolv-mimicry. The GNN's E advantage (0.086) is the least trustworthy
number in the table, and §3.4 shows that predicted E is in fact the descriptor that fails hardest
against chromatography-grade references. **Treat E accuracy claims as unresolved.**

GCA and kNN are otherwise near-indistinguishable: GCA slightly better on A, B, L; kNN slightly
better on S and E. The paper's own view is that "for most of the solute descriptors (S, A, and B),
the difference might not be significant in the application of the LSER models."

A **3-model consensus (GCA+kNN+GNN) is worse than the GNN alone** on median absolute error for every
descriptor — the GCA/kNN errors drag it down. The paper: "the median absolute difference of the
consensus approach is not lower; often, the GNN model shows more accurate results."

### 3.2 AbraLlama (5-fold CV, N = 6852, UFZ-LSER v3.2.1 "experimental" filter)

[Liquids 2024](https://doi.org/10.3390/liquids4030029), Table 1.
**Validated against: UFZ-LSER "experimental descriptors" — same contaminated chain as §3.1.**

| Descriptor | RMSE | R² |
|---|---|---|
| E | 0.132 | 0.97 |
| S | 0.240 | 0.90 |
| A | 0.135 | 0.85 |
| B | 0.123 | 0.96 |
| ~~V~~ | ~~0.097~~ | ~~0.98~~ *(disregard — V is exact)* |

Roughly kNN-class, clearly behind the 2026 GNN, and these are cross-validated rather than
independent-test-set numbers, so if anything they flatter the model. A is its weakest descriptor
(R² 0.85). The headline R² = 0.97 on E is exactly the kind of figure #5 warns about: it is agreement
with a compilation in which a sixth of E values are ACD/ChemSketch output and roughly half are
analogy estimates.

### 3.3 The A descriptor deserves separate treatment

Per the coordinator's note, `A` gates the ionisable-compound scope. Three findings, and they do not
point the way the Marchetto MMAPE figures suggest:

1. **A is 60% zeros.** "there is a significant value imbalance in the A descriptor data set, with 60%
   of data points having a value of zero" ([ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063), §3.2).
   The value is zero exactly when the molecule has no hydrogen-bond-donor group.
2. **Therefore percentage-error metrics on A are pathological.** MMAPE divides by values that are
   zero or near-zero for the majority of the set. The Marchetto blueprint's MMAPE 30.4 on A versus
   6.8–13.6 on E/S/B is consistent with A having the *smallest absolute* error of the four
   (RMSE 0.075–0.135 across every method in §3.1–3.2, the lowest of any descriptor except where B
   ties it) — it is a metric artefact, not evidence that A is hard.
   **Recommendation: the effort should standardise on absolute error for A, plus a separate
   `A = 0` / `A ≠ 0` classification metric, and should not use MAPE/MMAPE on A at all.**
3. **The real A risk is the binary call, not the magnitude.** The 2026 GNN's accuracy at
   distinguishing A = 0 from A ≠ 0 is **94.3%** ([ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063), §3.5).
   About one compound in eighteen gets its hydrogen-bond-donor status wrong outright. For an
   ionisable, drug-like set that is the failure mode that matters, and it is invisible in an RMSE.
   The paper found no benefit from oversampling or tailored loss functions for the imbalance.
4. **Chromatographically, A is the best-transferring descriptor.** Against the high-QC WSU reference
   set, SoluteML's A correlates at r² = 0.953 with slope 1.166 and intercept −0.003 — the strongest
   of E, S, A, B ([J. Chromatogr. Open 2025](https://doi.org/10.1016/j.jcoa.2025.100213), Eq. 8).
   SoluteGC's A is worse (r² = 0.852, slope 1.580, Eq. 9). *But* substituting SoluteML A into an
   otherwise-experimental RPLC model shifted the fitted `a` system constant from −0.384 to −0.462
   (~20%), which is why the authors conclude "SoluteML A descriptors cannot be used in place of WSU
   descriptors."

**Net:** A is not the weak link in absolute accuracy. It is the weak link in *interpretation* (the
`a` coefficient moves) and in *categorical correctness* (5.7% donor/non-donor errors). Both are
downstream-relevant and neither is captured by RMSE.

### 3.4 The descriptors that actually fail for RPLC are E and S

This is the most consequential finding in the whole ticket and it is invisible in §3.1.

Atapattu & Temerdashev ([J. Chromatogr. Open **7**, 100213, 2025](https://doi.org/10.1016/j.jcoa.2025.100213),
CC-BY) fitted the RPLC solvation parameter model on a Luna C18(2) column, methanol–water 30–70% v/v,
N = 32–39 compounds, using three descriptor sources. Correlation of estimated descriptors against the
WSU experimental reference (N = 39):

| Descriptor | SoluteML vs WSU | SoluteGC vs WSU |
|---|---|---|
| **E** | r² = **0.346**, slope 0.670 | r² = 0.320, slope 0.614 |
| **S** | r² = 0.711, slope 0.944 | r² = 0.800, slope 1.114 |
| **A** | r² = **0.953**, slope 1.166 | r² = 0.852, slope 1.580 |
| **B** | r² = 0.930, slope 0.886 | r² = 0.933, slope 0.996 |
| **L** | r² = 0.985, slope 0.973 (N=104) | r² = 0.970, slope 0.944 |

**Predicted E carries almost no usable signal for chromatography (r² = 0.346).** The symptom is
stark: in the fitted RPLC models the `e` coefficient collapses to exactly **0** in four of five
SoluteML models and in five of five SoluteGC models, whereas with WSU descriptors `e` is a healthy
0.143–0.384 ([ibid.](https://doi.org/10.1016/j.jcoa.2025.100213), Table 1). The e·E term is simply
lost. The same happens to the `s` and `a` system-constant *system maps*, which should vary smoothly
with mobile-phase composition but come out discontinuous with estimated descriptors.

Note the inversion: on the UFZ benchmark E is the *easiest* descriptor after A (GNN RMSE 0.086,
R² 0.98). Against chromatography-grade reference values it is the *worst* (r² 0.346). Per-descriptor
RMSE on UFZ data is therefore not a sufficient basis for choosing this layer, and #21's thin slice
should not assume it is.

**The contamination finding explains the inversion.** E is simultaneously the most contaminated
descriptor (16.7% ACD/ChemSketch/ClogP output plus ~48% analogy estimates, #5) and the one that
transfers worst to chromatography-grade references. The coherent reading: the models learned to
reproduce ACD/Labs' E, WSU independently determined E from partition and retention measurements, and
the two disagree. Offered as **consistent with**, not proven — the causal claim would need the
per-compound provenance flags joined to the WSU comparison set, which was not done in either paper.

Ranking the four RPLC targets by how well predicted values transfer to chromatography — which is the
ordering that should drive the error budget, not §3.1:

**A (0.953) > B (0.930) > S (0.711) ≫ E (0.346)**

That is close to the *reverse* of the §3.1 ordering, and it inverts the coordinator's prior that A is
the problem descriptor. **A is the best-transferring of the four.**

### 3.5 End-to-end error, which is what we actually care about

**RPLC log k standard error** (Luna C18(2), methanol–water)
([J. Chromatogr. Open 2025](https://doi.org/10.1016/j.jcoa.2025.100213), Table 1):

| Descriptor source | r² | SE(log k) |
|---|---|---|
| WSU (experimental) | 0.983 – 0.998 | **0.029 – 0.064** |
| SoluteML | 0.953 – 0.982 | **0.081 – 0.110** |
| SoluteGC | 0.923 – 0.947 | **0.134 – 0.169** |

Their conclusion: "for reversed-phase liquid chromatography and gas chromatography systems machine
learning and group contribution solute descriptors are not a substitute for experimentally determined
WSU descriptors."

**Downstream LSER properties** from the UFZ models
([ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063), Table 2), RMSE:

| Property | N | Consensus | **GNN** | GCA | kNN | DNNtaut | LSERD | ACD/Absolv |
|---|---|---|---|---|---|---|---|---|
| log K<sub>ow</sub> | 12010 | 1.03 | 1.09 | 1.06 | 1.65 | 1.04 | 0.91 | 0.87 |
| log K<sub>oa</sub> | 270 | 0.37 | **0.36** | 0.57 | 0.59 | 0.63 | 0.49 | 0.59 |
| log K<sub>wa</sub> | 696 | 1.32 | **1.13** | 1.91 | 1.30 | 1.36 | 1.28 | 1.39 |
| KRI (GC) | 454 | 129 | 109 | 157 | 160 | 129 | **104** | 120 |
| **CHI (RPLC)** | 204 | 3.80 | 3.64 | 5.38 | 4.41 | 4.49 | 5.52 | **3.19** |

CHI is the chromatographic hydrophobicity index — the only RPLC row. The GNN (3.64) beats GCA (5.38)
and kNN (4.41) but is beaten by commercial ACD/Absolv (3.19). kNN's log K<sub>ow</sub> is
catastrophically bad (1.65) and the authors explain why: the kNN training set overlaps >50% of the
log K<sub>ow</sub> evaluation set, so its apparent advantage elsewhere is partly leakage.

---

## 4. Uncertainty — the deciding axis

| Method | Emits? | Kind | Granularity | Calibrated? |
|---|---|---|---|---|
| **GCA** (2026) | Yes, ordinal | 3-level confidence from group coverage: *high* = all groups have fitted parameters; *moderate* = some group was decomposed; *low* = groups neither fitted nor decomposable | per compound, all descriptors alike | No. Not a distribution at all — an applicability flag |
| **kNN** (2026) | Yes, binary | Reliable / unreliable, via Aniceto et al. — per-training-compound radii from k-NN mean distances, corrected by sub-model spread and bias; a query is reliable if inside ≥1 radius | per compound per descriptor | No |
| **GNN** (2026) | Yes, numeric | SD across the 5 cross-validation models, binned *high* (SD ≤ 0.20), *moderate* (≤ 0.40), *low* (> 0.40); L uses 0.30/0.70 | per descriptor | Not assessed, but demonstrably discriminative — see below |
| **3-model consensus** | Yes, numeric | SD across GCA/kNN/GNN, terciles at the 33rd/66th percentile | per descriptor | Not assessed |
| **SoluteML** | **Yes, numeric** | `np.var` across the 25 ensemble members — **epistemic only** | per descriptor, marginal | Not assessed |
| **SoluteGC** | No | — | — | — |
| **AbraLlama** | **No** | — | — | — |
| **DeepChem** `GraphConvModel` | **Yes, numeric** | `predict_uncertainty()`: 50 MC-dropout masks + a learned log-variance head; `std = sqrt(E[p²] − E[p]² + E[var])` → **epistemic + aleatoric** | per task, marginal | User's responsibility |

### 4.1 The 2026 GNN's ensemble SD works as a triage signal

Restricting to the top tercile of consensus SD collapses the error dramatically
([ACS Environ. Au 2026](https://doi.org/10.1021/acsenvironau.6c00063), §3.7):

| Descriptor | RMSE, all | RMSE, high-confidence tercile | n (of test set) |
|---|---|---|---|
| A | 0.087 | **0.009** | 228 / 634 |
| E | 0.103 | **0.039** | 209 / 634 |
| B | 0.106 | **0.052** | 200 / 606 |
| S | 0.168 | **0.080** | 210 / 635 |
| L | 0.290 | **0.144** | 181 / 548 |

A ~2–10× error reduction on the confident third. That is exactly the behaviour a peak-identity prior
needs — but note only ~30% of compounds qualify, and the paper is candid: "it is hard to judge
whether the prediction improves for a given chemical."

**No paper in this comparison reports a calibration curve, a coverage plot, an interval score, or
any proper scoring rule.** Every uncertainty claim above is a *ranking* claim, not a *calibration*
claim. Calibration is unresolved in the literature and is ours to establish (#14).

### 4.2 SoluteML's shipped variance is epistemic-only and measurably under-confident

> **Measured since this was written** (`prototype/soluteml-bias/FINDINGS.md`, #35): against
> the WSU-94, RMSE / median ensemble SD is **3.95 (E), 3.70 (S), 3.85 (A), 2.15 (B)** in
> descriptor space and **1.23** [0.98, 1.47] in retention space. The 25-member spread
> captures roughly a quarter of the actual descriptor error. Both are floors — every anchor
> compound is inside the training corpus. The section below is the mechanism; those are the
> numbers.

From the source
([`solvation_estimator.py`](https://github.com/fhvermei/chemprop_solvation/blob/master/chemprop_solvation/solvation_estimator.py)):

```python
for j in range(5):  # the 5 predictions correspond to E, S, A, B, and L.
    epi_unc_per_data.append(np.var(all_preds[:, i, j]))
    ave_pre_per_data.append(np.mean(all_preds[:, i, j]))
```

`all_preds` has shape `(25, N, 5)` — 5 folds × 5 models. The returned "uncertainty" is a **variance,
not a standard deviation**, and it is the disagreement among ensemble members only: no noise term, no
data-uncertainty term, no calibration.

Sanity check on the repo's own sample output
([`sample_files/example_property_pred/results.csv`](https://github.com/fhvermei/SolProp_ML/blob/master/sample_files/example_property_pred/results.csv)):
`uncertainty_SoluParam_E = 2.16e-4` → σ ≈ **0.015**. The same model's E disagrees with WSU reference
values at SE = 0.252 (§3.4). The shipped uncertainty is understated by more than an order of
magnitude. **Do not use it as-is.** It needs an aleatoric term and empirical recalibration.

### 4.3 DeepChem: what it actually provides

[deepchem/deepchem](https://github.com/deepchem/deepchem), **MIT**, ~6.9k stars, actively maintained
(last push 2026-08-10).

- **Models** for this task: `GraphConvModel`, `GCNModel`, `AttentiveFPModel`, `MPNNModel`,
  `DMPNNModel`, `InfoGraphModel`; plus `HuggingFaceModel` wrappers for ChemBERTa, MoLFormer
  ([API reference](https://deepchem.readthedocs.io/en/latest/api_reference/models.html)).
- **Featurizers**: `ConvMolFeaturizer`, `MolGraphConvFeaturizer`, RDKit descriptor featurizers.
- **Multi-task is native** — `n_tasks > 1` gives all five descriptors from one graph encoder.
- **Uncertainty is first-class and built in.** `GraphConvModel(..., uncertainty=True)` requires
  regression mode and `dropout > 0`, and adds a `Dense(n_tasks)` log-variance head with an `exp`
  activation ([`deepchem/models/graph_models.py`](https://github.com/deepchem/deepchem/blob/master/deepchem/models/graph_models.py),
  lines ~605–660, ~803–900). Then `predict_uncertainty(dataset, masks=50)` runs 50 dropout-enabled
  passes and returns `(y_pred, y_std)` with
  `std = sqrt(sum_sq_pred/masks − p² + sum_var/masks)` — Monte-Carlo-dropout epistemic **plus**
  learned aleatoric, per [Kendall & Gal, arXiv:1703.04977](https://arxiv.org/abs/1703.04977)
  ([`torch_model.py`](https://github.com/deepchem/deepchem/blob/master/deepchem/models/torch_models/torch_model.py),
  `predict_uncertainty`, lines ~784–840).
- **No pretrained Abraham-descriptor weights.** DeepChem ships architecture and training machinery,
  not a solute-descriptor model. ChemBERTa/MoLFormer are general SMILES encoders that would need
  fine-tuning on our data.
- **Caveat**: `GraphConvModel` is the legacy TF/Keras path (this is what the 2026 UFZ paper used,
  DeepChem 2.7.2). New work should prefer the PyTorch models; `TorchModel.predict_uncertainty` has
  the identical contract.
- **Limitation that matters here**: `predict_uncertainty` returns a per-task `y_std`. It gives
  **no cross-task covariance**, even for a multi-task model. Same gap as SoluteML.

---

## 5. Joint vs independent prediction — and the correlation question

This was called out in the ticket, and the answer is unsatisfying: **the literature almost entirely
ignores it.**

| Method | Joint or independent? | Error covariance reported? |
|---|---|---|
| GCA (2026) | **Independent** — separate parameter fit per descriptor | No |
| kNN (2026) | **Independent** — hyperparameters optimised separately per descriptor | No |
| GNN (2026) | **Independent** — one network per descriptor, different epoch counts per descriptor (S 160, E 100, A 50, B 140, L 85), different N per descriptor (5,572–6,364) | No |
| DNNmono / DNNtaut (2022) | **Both** — the repo ships `A_descriptor/`…`L_descriptor/` *and* `multitask_all_descriptors/` | No |
| **SoluteML** | **Joint (multi-task)** — one D-MPNN, 5 outputs from a shared molecular representation | No — marginal variances only |
| AbraLlama | **Joint** — one fine-tuned model, 5 outputs | No uncertainty at all |
| DeepChem | Either (`n_tasks`) | Per-task `y_std` only |

Things the 2026 paper *does* do, which are adjacent but not the same:

- It checks **correlations of the predicted values** E↔L and E↔S, and trend analyses for functional
  groups (A, B) and homologous series (L, E), to confirm the models "follow the expected
  intermolecular interaction patterns and are not merely a result of statistical fitting"
  (§3.6, SI 1-5 / 1-6 / 1-7). That is a *physical-consistency* check on the predictions, **not** a
  characterisation of the joint error distribution.
- Its ensemble SDs are computed and thresholded independently per descriptor.

### Why this matters for #14

Elution order depends on **differences** in log k. For two solutes with descriptor error vectors
δ₁, δ₂ over the four predicted targets and the corresponding column coefficients
**m** = (e, s, a, b):

```
Var(Δ log k) = mᵀ (Σ₁ + Σ₂) m          — with the full 4×4 covariance Σ
```

V contributes no term: it is computed exactly, so its error is zero and it drops out of the
propagation entirely. **The uncertainty problem is four-dimensional, not five.**

Dropping the off-diagonals of Σ discards every cross-descriptor term. Whether that understates or
overstates Var(Δ log k) depends on the signs: RPLC coefficients have mixed signs (Luna C18(2) at 50%
MeOH: e = +0.299, s = −0.742, a = −0.384, b = −1.699, v = +2.078
— [J. Chromatogr. Open 2025](https://doi.org/10.1016/j.jcoa.2025.100213), Eq. 10), so positively
correlated E and S errors *reduce* the variance of log k, while positively correlated S and B errors
*increase* it. A diagonal-only treatment is simply wrong, with no safe direction.

**The actionable finding.** Because SoluteML is multi-task and its 25 ensemble members each emit all
descriptors, the **empirical covariance is recoverable with a ~5-line change** to
`load_ML_estimator`: return the raw `all_preds` tensor of shape `(25, N, 5)` instead of collapsing it
to per-`j` variances, then take `np.cov(all_preds[:, i, :].T)` and slice out the E/S/A/B block. No
retraining. This is the single highest-value modification available to us and it should be the first
thing #14 does.

No equivalent path exists for the 2026 UFZ models, which are five genuinely separate networks —
their errors are correlated in reality but the models cannot express it, and stacking their
independent SDs into a diagonal Σ would be a fiction.

**Caveat, stated plainly**: an ensemble covariance is a covariance of *epistemic* disagreement. It is
not the covariance of the true error unless calibrated against held-out data with known descriptors.
That calibration has not been done by anyone and is work for #14.

---

## 6. Applicability domain

| Method | AD mechanism | Assessed how |
|---|---|---|
| GCA (2026) | Group coverage → 3-level flag | Structural: does every functional group have a fitted parameter? |
| kNN (2026) | Aniceto et al. density method: per-training-compound radii from k-NN mean-distance distribution, corrected by sub-model spread and bias; reliable if the query falls inside ≥1 radius. Training-set predictions obtained by leave-one-out CV | Explicit, per query |
| GNN (2026) | **None explicit** — ensemble SD used as a proxy | — |
| Whole approach | **UMAP projection against PubChemLite** — experimental descriptors cluster in the centre of chemical space; complex, multi-functional molecules sit at the edges | Global, qualitative |
| SoluteML / SoluteGC | Element whitelist H, B, C, N, O, S, P, F, Cl, Br, I; **neutral, non-mixture solutes only**; evaluated with **substructure-based splits** as well as random splits | Substructure splits are the strongest evaluation protocol in this comparison |
| AbraLlama | PCA distance-from-centre correlates with error (p < 0.0001) — but demonstrated for the **solvent** model only; no AD analysis for the solute model | Weak |

Limits the 2026 paper states about the LSER framework itself, all directly relevant to our
drug-like, ionisable scope:

- "As the experimental descriptors mainly cover simple chemical structures, the approach is limited
  for more complex structures."
- "Especially for chemicals that involve intramolecular interactions (e.g., hydrogen bonds), the LSER
  approach tends to be more prone to error. These interactions depend on the different conformers of
  a molecule in the liquid phase, which are not accounted for in the approach."
- "the A and B descriptors, which cover the overall H-bond capacity, may not fully capture
  site-specific effects or steric hindrance. Similarly, the S descriptor aggregates various dipolar
  and electrostatic interactions into a single term."
- **"Additionally, ionic chemicals are not covered by the approach at all. Because Coulomb
  interactions strongly affect physicochemical properties such as solubilities and partition
  coefficients, the LSER approach is not suitable for describing such chemicals."**
- "many LSER equations are developed using fewer than 200 chemicals, further narrowing their
  applicability domain."

Subgroup analysis (§3.6): kNN and GNN both degrade on "more heterogeneous chemicals" across all
descriptors. kNN is worse on N/O-containing than pure C/H compounds and worse on halogens. GCA has
the opposite pattern for A/B/S (worse on CHNO, because additive schemes cannot capture
neighbouring-group effects) and for E/L is better on CHNO (because size is additive).

### 6.1 Our chemistry is outside the domain — quantified

#5 characterised the corpus every one of these models is trained on: **small neutral environmental
solutes, median MW 174, benzene accounting for 26% of all scaffolds, and zero charged species by
design.** The median drug-like ionisable compound sits at roughly the **95th percentile** of the S, B
and V distributions.

This is not "the edge of the domain". It is past it. **Every method compared here is extrapolating on
our in-scope chemistry**, and the published accuracy figures — all measured on random hold-outs from
the same small-neutral corpus — describe interpolation performance that we will not get.

### 6.2 How each method behaves when extrapolating, and whether it knows

This is the axis that should decide the choice, because it is the one our chemistry actually lands on.

| Method | Behaviour beyond the domain | Does it detect it? |
|---|---|---|
| **GCA / SoluteGC** | Additive and unbounded — extrapolates linearly, so a large multi-functional molecule accumulates group terms with no saturation. Fails specifically on neighbouring-group effects and intramolecular H-bonding, which is what multi-functional drug-like molecules *are*. Errors are systematic, not noisy | **Yes, and structurally.** An unseen functional group cannot be silently absorbed: it either has a fitted parameter (*high*), was decomposed (*moderate*), or triggers *low*. This is the most honest out-of-domain detector in the comparison — but it is a flag, not a distribution |
| **kNN** | **Cannot extrapolate at all.** A prediction is a distance-weighted average of training values, so it is mathematically bounded by the training range. For a compound at the 95th percentile of S or B it will regress toward the interior: the error is **biased low**, systematically, not merely uncertain. This is the worst possible failure mode for a peak-identity prior — confidently wrong in a known direction | **Yes, explicitly** — the Aniceto radii are a genuine per-query domain test, and the paper's plots mark unreliable predictions. Detection is good; the underlying behaviour is not |
| **GNN (2026), SoluteML, DeepChem ensembles** | Unconstrained non-linear extrapolation. Can produce values outside the training range, in either direction, with no structural check | **Poorly, and this is the trap.** The only signal is ensemble spread or MC-dropout spread. Ensemble members trained on the same biased corpus tend to *agree* out of domain, so the spread **shrinks** exactly where it should widen. A deep ensemble's confident agreement off-distribution is a well-known artefact, and no paper here tests for it. **SoluteML is the only candidate with any evidence on this**, because Chung et al. evaluated on **substructure-based splits** as well as random ones — a genuine held-out-chemotype test |
| **AbraLlama** | Unconstrained | **No.** No uncertainty output, no AD analysis for the solute model. The PCA-distance-correlates-with-error result is for the *solvent* model only |
| **Hard rejections** | — | SoluteML/SoluteGC refuse outright outside the element whitelist (H, B, C, N, O, S, P, F, Cl, Br, I) and for mixtures and **charged species**. A refusal is more useful than a confident wrong answer, and it is the only mechanism here that explicitly names ionic solutes |

**The ranking on this axis is not the ranking on accuracy.** GCA detects out-of-domain best and
predicts worst; kNN detects well but fails in a biased direction; the neural models predict best and
detect worst. SoluteML is the only one with both a numeric uncertainty *and* substructure-split
evidence, which is why it survives §10 despite not leading §3.1.

**Consequence for this effort.** Our target chemistry — drug-like, multi-functional, ionisable — sits
precisely where every method is weakest *and* where the uncertainty estimates are least trustworthy.
Any calibration done on random hold-outs will be optimistic for our compounds. **#14 must calibrate on
a scaffold- or substructure-split, never a random split**, and should treat off-distribution
uncertainty as needing an explicit inflation term rather than trusting ensemble spread.

---

## 7. Availability: pretrained and openly licensed?

| Method | Pretrained weights? | Licence | Notes |
|---|---|---|---|
| **SoluteML / SoluteGC / DirectML** | **Yes** — 25-member ensemble, `conda install -c fhvermei solprop_ml`, or [Zenodo 5792296](https://zenodo.org/records/5792296) (`ML_model_files.zip`, 2.0 GB) | Data **CC-BY-4.0**; [SolProp_ML](https://github.com/fhvermei/SolProp_ML) "free, open-source … Creative Commons Attribution 4.0 International"; [chemprop_solvation](https://github.com/fhvermei/chemprop_solvation) **MIT** | Also a web tool at [rmg.mit.edu/database/solvation/search](https://rmg.mit.edu/database/solvation/search/). SoluteGC ships inside RMG-Py / RMG-database. Python API: `predict_property(..., solute_parameters=True)`. macOS/Linux only ("may not work on Windows") |
| **2026 GCA / kNN / GNN** | **No.** "All models are implemented in the PAULY software found at https://i-am-pauly.com/" — a UFZ spin-off selling custom applications and consulting. No licence, no pricing, no API documented | Unknown / commercial | **The paper's data link does not deliver.** It cites `github.com/nadinulrich/solute_descriptor_prediction` for "the up-to-date data set (May 2026)", but that repo's last commit is **2021-11-07** and it contains the *2022* study's supplement. Verified via GitHub API |
| **DNNmono / DNNtaut (2022)** | **Partially — not usable.** Repo has `Prediction_solute_descriptors_{single,multi}task.py` and 6 KB `model.pickle` files per descriptor plus a `multitask_all_descriptors/` variant, but **no TensorFlow checkpoints**. Also targets the removed `deepchem.models.tensorgraph` API | **MIT** | Useful as a reference implementation and for the `Dataset_and_predictions.xlsx` (4.8 MB), not as a pretrained model |
| **AbraLlama** | **Yes** — [HF Space `ttmn/AbraLlama`](https://huggingface.co/spaces/ttmn/AbraLlama), code at [BrightBlueCheese/AbraLLaMA](https://github.com/BrightBlueCheese/AbraLLaMA) | Paper CC-BY; repo licence not stated (see §9) | No uncertainty output; predicts V |
| **DeepChem** | Framework only — no Abraham weights | **MIT** | Full training + uncertainty machinery; ChemBERTa / MoLFormer available for fine-tuning |
| **ACD/Absolv** | Commercial | Proprietary | Best CHI performance in Table 2 (3.19). Out of scope for an open stack, but a useful accuracy ceiling |

**The accuracy leader is not obtainable.** The 2026 GNN wins every column in §3.1 and is locked behind
a commercial spin-off with no published licence. Its hyperparameters *are* fully published (§2), so
reproducing it is possible — but only if we can obtain the training dataset, and the repository the
paper points at does not currently contain it.

---

## 8. Ground-truth provenance — CONFIRMED CONTAMINATED

**Resolved by #5, whose findings are
summarised here and are authoritative over anything in this file.** The concern is confirmed, with
numbers, and it is worse than a caveat — it changes what the §3.1–3.3 tables mean.

### 8.1 What #5 established

- **The "three independent sources" are one dataset.** Acree's papers feed Abraham's Absolv
  compilation; its two obtainable copies — UFZ-LSER and MIT SoluteDB — **overlap 92%**. The apparent
  provenance diversity between the 2026 UFZ models (trained on Absolv/LSERD) and SoluteML (trained on
  SoluteDB) is largely illusory. *This retracts the "different provenance chain" point I made in
  favour of SoluteML in an earlier draft — it is worth at most 8% of the data.*
- **16.7% of E values and 22.1% of B values are output from ACD Absolv, ChemSketch or ClogP.**
- **~48% of E and S values are "estimated by comparison to closely related compounds"** — expert
  analogy, not measurement.
- **V is computed** (McGowan). **S, A, B and L are regression-derived**, not directly measured.
- **No compilation reports uncertainties at all.**

### 8.2 What that does to the published scores

The headline QSPR figures — E R² = 0.98, B R² = 0.95 — **partly measure how learnable ACD/Labs' own
predictor is.** A model that "wins" on E may simply be the best Absolv mimic. Two independent
observations line up with this:

1. ACD/Absolv itself scores 0.10 RMSE on E against the Absolv-derived test set — near-best in the
   table (§3.1). It would be circular if *all* labels were ACD output; at 16.7% it is partially so.
2. Predicted E is the descriptor that collapses hardest against the independent WSU set
   (r² = 0.346, §3.4). Highest apparent accuracy, worst real transferability, highest contamination —
   the same descriptor. That coincidence is the whole problem in one column.

**Reporting rule adopted throughout this document:** no accuracy figure appears without a statement
of what it was validated against and whether that reference is measured, estimated, or computed.

### 8.3 The honest validation set is twelve compounds

#5 found that **clean-provenance AND drug-like AND ionisable comes to twelve compounds.** That number
should govern how much weight any accuracy claim carries in this effort.

What follows from it:

- **No candidate method can be honestly *evaluated* on our in-scope chemistry.** Twelve compounds
  supports a smoke test, not a model comparison; with four descriptors it cannot distinguish methods
  whose published RMSEs differ by 0.05.
- **No candidate can be *trained* on a clean-provenance subset either.** Even relaxing "drug-like and
  ionisable", the clean-provenance fraction of a ~6,300-compound corpus is far below what a D-MPNN or
  a graph-convolutional network needs. Group contribution would tolerate it best (fewest parameters),
  and kNN degenerates fastest (it *is* the training set).
- **The realistic path is fine-tuning, not retraining**: keep a model pretrained on the full
  contaminated corpus for its representation, then recalibrate its *output* on the cleanest
  chromatography-derived set available. That is the WSU database (~300 compounds, three laboratories,
  one QC protocol) — which is why the §11 recommendation to re-fit on WSU is the highest-value
  follow-up in this ticket.
- **Therefore accuracy cannot decide this ticket.** The decision must rest on what is verifiable
  without a clean test set: licence, availability of weights, whether the architecture is multi-task,
  whether uncertainty is emitted, and whether the method can tell it is extrapolating (§6.2). That is
  exactly what §10 does.

### 8.4 Per-method provenance record

| Method | Validated against | Measured, estimated or computed? |
|---|---|---|
| 2026 GCA / kNN / GNN | "Abraham Absolv" via LSERD, 6,364 compounds | **Mixed and contaminated.** ~17% of E and ~22% of B are ACD/ChemSketch/ClogP output; ~48% of E and S are analogy estimates (#5) |
| DNNtaut (2022) | Same set, same test split | Same chain |
| AbraLlama | UFZ-LSER v3.2.1, "experimental descriptors" filter, N = 6852 | Same chain — the "experimental" label does not mean measured |
| SoluteML / SoluteGC | **SoluteDB**, 8,366 parameters compiled by Chung et al. | Same chain — **92% overlap with UFZ-LSER** (#5). Not the independent reference it appears to be |
| §3.4 / §3.5 chromatographic evaluation | **WSU database** (~300 compounds), refined from GC, RPLC and liquid–liquid partitioning, "all experiments were conducted under consistent quality control standards across three laboratories" | **Highest-integrity reference in this document.** Still regression-derived, but under one protocol |

Corroborating evidence, recorded for completeness:

- **Independent corroboration from chromatography**: Atapattu & Temerdashev conclude the 2026-era
  estimators are calibrated on the wrong reference — "Machine learning and group contribution solute
  descriptor estimation processes were based on calibration on UFZ-LSER descriptors. It is in our
  opinion if the machine learning and group contribution solute descriptor estimation process were to
  calibrate using more refined WSU solute descriptors estimation of solute descriptors using machine
  learning and group contribution would be more reliable"
  ([J. Chromatogr. Open 2025](https://doi.org/10.1016/j.jcoa.2025.100213), §4).
- A secondary source characterises UFZ-LSER as containing "experimentally determined solute
  descriptors without standardized protocols", with inconsistencies from periodic updates and
  differing experimental procedures. **Not traced to a primary source** — consistent with #5.

Further reading, identified but not read here:
**Poole, "The influence of descriptor database selection on the solvation parameter model for
separation processes", J. Chromatogr. A 1692, 463851 (2023)**,
[doi:10.1016/j.chroma.2023.463851](https://doi.org/10.1016/j.chroma.2023.463851).

**Why the WSU/UFZ disagreement is structural, not a data-cleaning problem.** Except V, Abraham
descriptors are never directly measured — they are back-solved by regression from measured partition
and retention data. "Experimental" means "regression-derived from measurements", and *which*
measurements were used is itself a variance source. WSU (~300 compounds, three labs, one protocol)
and UFZ-LSER (>8,000 compounds, no standardised protocol, plus the contamination in §8.1) are
therefore not two samples of one quantity; they are two differently-defined quantities that happen to
share a name. Fitting to one and testing on the other produces r² = 0.346 on E, and no amount of
curation makes that go away.

---

## 9. What could NOT be determined

1. **Per-descriptor RMSE/MAE for SoluteML and SoluteGC on their own test sets.** JCIM 2022 is
   paywalled; the ChemRxiv PDF, the KU Leuven Lirias accepted-manuscript copy, and the ResearchGate
   copy all refused automated retrieval. Abstract, architecture, dataset sizes and code were obtained
   from primary sources; the accuracy table was not. **§3.4/§3.5 substitute a more decision-relevant
   measurement** (correlation with WSU, and log k standard error on a real RPLC column), so this gap
   does not block the recommendation — but the raw numbers should be read from the PDF before the
   spec is finalised.
2. **Whether the 2026 UFZ models will ever be openly licensed.** PAULY publishes neither licence nor
   pricing; the site advertises custom applications, consulting and workshops. A beta exists at
   `beta.pauly-predictions.de`. Contact would be needed.
3. ~~Whether the "Absolv" subset is measured or ACD-computed.~~ **Resolved by #5** — substantially
   computed and estimated; see §8.
4. **Whether the AbraLlama GitHub repository carries any licence.** The paper is CC-BY; the code
   repository's terms were not confirmed. Moot given the recommendation.
5. **Contents of Poole's 2026 review**, "Experimental and computational methods for assigning
   descriptors for the solvation parameter model", J. Chromatogr. Open 10, 100366,
   [doi:10.1016/j.jcoa.2026.100366](https://doi.org/10.1016/j.jcoa.2026.100366). Identified via
   Crossref metadata only; 102 references; **not read**. It is the most recent authority on this
   exact question and should be read before the spec closes.
6. **Any calibration assessment of any method.** No paper found reports a calibration curve, coverage
   at nominal level, or a proper scoring rule for descriptor uncertainty. This is a genuine gap in the
   literature, not a gap in the search.
7. **Cross-descriptor error covariance, empirically.** Nobody reports it. §5 gives the route to
   compute it ourselves.
8. **The §8.1 provenance statistics were not independently verified here.** The 92% overlap, the
   16.7%/22.1% computed fractions, the ~48% analogy-estimate fraction and the twelve-compound
   intersection are relayed from #5, which owns that investigation. They are consistent with
   everything I could check directly (the WSU/UFZ divergence on E, the Atapattu recalibration
   recommendation, ACD/Absolv's near-best score on the Absolv-derived test set), but this document is
   not their primary source.
9. **Whether ensemble spread degrades off-distribution for these specific models.** §6.2 states the
   general mechanism and notes no paper tests it. I did not test it either — it needs the model
   weights and a held-out out-of-domain set. Filed as open question 7.

---

## 10. Recommendation

**Accuracy does not decide this ticket.** The published scores are measured against a contaminated
reference (§8.1–8.2), and the honest validation set for our in-scope chemistry is twelve compounds
(§8.3). The decision therefore rests on the properties that are verifiable without a clean test set:
licence and weight availability, joint vs independent architecture, whether uncertainty is emitted at
all, and whether the method can tell it is extrapolating.

### Primary: SoluteML, with three mandatory modifications

Use **SoluteML** (multi-task D-MPNN, 25-member ensemble) as the QSPR layer, predicting **E, S, A, B**.

Reasoning, in the order the ticket weights it:

1. **Uncertainty is a hard requirement, and SoluteML is the only openly-available candidate that
   emits a numeric per-descriptor uncertainty from a real ensemble.** SoluteGC and AbraLlama emit
   nothing. The 2026 GNN emits a better-validated ensemble SD but is not obtainable (§7).
2. **It is multi-task**, so the descriptors share a representation and the **4×4 error covariance is
   recoverable without retraining** (§5). Decisive for #14: elution order depends on differences in
   log k, and a diagonal covariance gets those wrong in an unpredictable direction.
3. **It is genuinely open and pretrained** — CC-BY-4.0 weights, MIT model code, conda-installable.
   *(An earlier draft also credited SoluteDB as an independent provenance chain. #5 retracts that:
   SoluteDB overlaps UFZ-LSER by 92%. The licence and weights argument stands; the data-independence
   argument does not.)*
4. **It is the only candidate independently validated on an actual reversed-phase column**, and it
   beat group contribution there on every measure (§3.4, §3.5).
5. **It is the only candidate with any evidence of out-of-domain performance** — Chung et al.
   evaluated on substructure-based splits, not only random ones (§6.2). Given that our chemistry sits
   at the 95th percentile of the training distribution, this matters more than a 0.05 RMSE difference
   on an interpolation benchmark.
6. Its accuracy deficit versus the 2026 GNN is real but second-order: the GNN's largest margin is on
   **L**, which the RPLC equation does not use at all.

Mandatory modifications:

- **M1. Expose the ensemble covariance.** Patch `load_ML_estimator` to return the raw `(25, N, 5)`
  prediction tensor; compute `np.cov` per compound and slice the E/S/A/B block. ~5 lines, no
  retraining. Prerequisite for #14.
- **M2. Do not trust the shipped variance.** Epistemic-only and understated by >10× (§4.2). Add an
  aleatoric term and recalibrate empirically — **on a scaffold or substructure split, never a random
  one** (§6.2).
- **M3. Compute V from structure and never predict it. Do not build an L model.** Four targets, not
  five (§1).

### Also adopt, regardless of model choice

- **AD gating.** Port the kNN Aniceto-style applicability-domain check (§6) as a model-agnostic
  pre-filter, and adopt the tercile-triage pattern from §4.1: a prediction outside the domain should
  widen the prior, not silently pass through. **Do not rely on ensemble spread alone** — it shrinks
  off-distribution (§6.2), which is precisely backwards for our chemistry.
- **Metrics.** Absolute error per descriptor, never MAPE/MMAPE — especially for A, where 60% zeros
  make percentage error meaningless (§3.3). Add a separate `A = 0` / `A ≠ 0` classification metric;
  the 5.7% donor/non-donor error rate is a real risk that no RMSE surfaces.
- **Weight the error budget by chromatographic transferability, not by §3.1.** The ordering that
  matters is A (r² 0.953) > B (0.930) > S (0.711) ≫ E (0.346).
- **Error budget for #21.** Assume SE(log k) ≈ **0.08–0.11** contributed by this layer alone, against
  0.03–0.06 achievable with measured descriptors (§3.5). If the thin slice needs better than that,
  the descriptor model is not the constraint to relax — the reference data is.
- **Never quote a published QSPR accuracy figure without its provenance qualifier** (§8.2).

### Rejected

- **Group contribution (SoluteGC, 2026 GCA)** — worst chromatographic fit (SE 0.134–0.169, `e`
  collapses to zero in every model), no numeric uncertainty, and additive schemes fail specifically on
  the neighbouring-group and intramolecular-H-bond effects that define drug-like molecules. Its
  3-level confidence flag is the best out-of-domain *detector* in the comparison and worth copying as
  an idea; the model is not.
- **kNN** as the primary — and now for a stronger reason than accuracy: **it cannot extrapolate at
  all** (§6.2). Its predictions are bounded by the training range, so for compounds at the 95th
  percentile of S or B — i.e. ours — it is biased toward the interior, confidently and systematically.
  Its downstream log K<sub>ow</sub> (1.65) is the worst in the comparison and is additionally inflated
  by >50% train/test overlap. **Keep its Aniceto AD machinery**, which is genuinely the best
  applicability-domain treatment available.
- **AbraLlama** — no uncertainty of any kind, disqualifying under this architecture's hard
  requirement. No AD analysis for the solute model. That it predicts V, a quantity computable exactly,
  is a further signal.
- **The 2026 UFZ GNN** as the *implementation* — best accuracy on all four RPLC descriptors, best CHI
  among open methods, and the only work with published AD analysis across three model families, but
  locked in PAULY with no licence and a dead data link. Note also that its five independent
  single-task models **cannot express descriptor error correlation even in principle**, so adopting it
  would forfeit the #14 covariance route. **Track it**: if weights or the May-2026 dataset are
  released, revisit. Hyperparameters are fully published, so retraining is feasible given data.

---

## 11. New questions this raises

1. **Are descriptors the right intermediate at all?** Chung et al.'s *DirectML* — predicting
   solvation free energy directly rather than via descriptors — outperformed both SoluteML and
   SoluteGC, and the 2026 paper's conclusion points the same way ("direct predictions of solvation
   energies or activity coefficients by NNs might be useful alternatives, especially for
   [ionogenic] chemicals"). A DirectML-style model for log k would bypass the descriptor bottleneck
   entirely — at the cost of the column-transferability that motivated LSER in the first place.
   **Worth a ticket**, and it bears directly on the map's committed column axis.
2. **Should we re-fit SoluteML on WSU descriptors? — now the highest-value follow-up in this ticket.**
   Atapattu & Temerdashev explicitly recommend it, and #5's contamination finding makes it close to
   mandatory: fine-tuning a pretrained ensemble's *output layer* on ~300 chromatography-derived,
   single-protocol compounds is cheap, targets exactly our failure mode, and is the only route that
   escapes the Absolv chain without needing a corpus we do not have. **Strong candidate for a
   prototype ticket.** Open sub-questions: does 300 compounds support fine-tuning four heads without
   catastrophic forgetting, and does WSU coverage extend into ionisable drug-like space at all
   (probably not — see Q4).
3. **Should the e·E term be kept?** Predicted E correlates with chromatography-grade E at r² = 0.346,
   and the fitted `e` coefficient collapses to zero when estimated descriptors are used. Keeping a
   term whose input is noise may be worse than dropping it. Testable cheaply in #21.
4. **Ionisable compounds have no descriptor path.** LSER "is not suitable for describing" ionic
   species, and no candidate model predicts ionic-species descriptors. The literature's answer is a
   degree-of-ionisation descriptor D (split into D⁺/D⁻ for bases/acids) that averages neutral and
   ionic descriptors by ionisation fraction — which still requires descriptors *for the ion*.
   **This needs its own ticket** and it interacts with the pKa provider decision.
5. **Does the covariance actually change anything?** Before building machinery for it, #14 should run
   the cheap numerical check: propagate the SoluteML 4×4 covariance through a published RPLC
   coefficient vector (e.g. Luna C18(2) 50% MeOH: e +0.299, s −0.742, a −0.384, b −1.699) and compare
   Var(Δlog k) with and without off-diagonals. If the difference is under ~10%, the diagonal
   approximation survives and a lot of complexity disappears.
6. **How do we validate calibration when there is no clean ground truth?** This is now the sharpest
   open problem. Calibrating against UFZ/Absolv descriptors calibrates against a reference that is
   ~17–22% commercial-model output and ~48% analogy estimate (§8). Calibrating against WSU gives ~300
   compounds, none of them ionisable drug-like. Calibrating against the clean-provenance in-scope
   intersection gives **twelve**. Calibrating against *retention* sidesteps provenance entirely but
   confounds this layer's error with the LSER layer's — which may nonetheless be the right answer,
   since retention is what we actually care about and the map already commits to layer-by-layer
   validation with pairwise elution-order accuracy as a primary metric. **#14 should decide this
   explicitly rather than defaulting to descriptor-space calibration.**
7. **Is an ensemble-spread uncertainty usable at all off-distribution?** §6.2 argues it shrinks
   exactly where it should widen, and no paper tests this. A cheap experiment: take compounds far from
   the training distribution with known WSU descriptors, and check whether SoluteML's ensemble
   variance rises with true error or stays flat. If it stays flat, the uncertainty layer needs a
   distance-based inflation term and the ensemble variance alone is not fit for purpose.

---

## Sources

**Primary papers**

- Ulrich, Istomin, Kudria, Böhme, Voigt. *Prediction of Solute Descriptors for Linear Solvation
  Energy Relationships Using K-Nearest Neighbors, Group Contributions, and Graph-Convolutional Neural
  Networks.* ACS Environ. Au, 10 June 2026. CC-BY-4.0.
  <https://doi.org/10.1021/acsenvironau.6c00063>
- Chung, Vermeire, Wu, Walker, Abraham, Green. *Group Contribution and Machine Learning Approaches to
  Predict Abraham Solute Parameters, Solvation Free Energy, and Solvation Enthalpy.*
  J. Chem. Inf. Model. 2022, 62(3), 433–446. <https://doi.org/10.1021/acs.jcim.1c01103>
  (preprint: <https://doi.org/10.26434/chemrxiv-2021-djd3d-v2>)
- Atapattu, Temerdashev. *Assessment of machine learning and group contribution solvation parameter
  model descriptors for model retention in reversed-phase liquid chromatography and gas
  chromatography.* J. Chromatogr. Open 2025, 7, 100213. CC-BY.
  <https://doi.org/10.1016/j.jcoa.2025.100213>
- Lang, Lee. *AbraLlama: Predicting Abraham Model Solute Descriptors and Modified Solvent Parameters
  Using Llama.* Liquids 2024, 4(3), 518–524. CC-BY. <https://doi.org/10.3390/liquids4030029>
- Ulrich, Ebert. *Can deep learning algorithms enhance the prediction of solute descriptors for linear
  solvation energy relationship approaches?* Fluid Phase Equilib. 2022, 555, 113349.
  <https://doi.org/10.1016/j.fluid.2021.113349>
- Abraham, McGowan. *The use of characteristic volumes to measure cavity terms in reversed phase
  liquid chromatography.* Chromatographia 1987, 23, 243–246. <https://doi.org/10.1007/BF02311772>
- Poole. *Solvation parameter model: tutorial on its application to separation systems for neutral
  compounds.* J. Chromatogr. A 2021, 1645, 462108. <https://doi.org/10.1016/j.chroma.2021.462108>
- Poole. *The influence of descriptor database selection on the solvation parameter model for
  separation processes.* J. Chromatogr. A 2023, 1692, 463851. *(not read — for #5)*
  <https://doi.org/10.1016/j.chroma.2023.463851>
- Poole. *Experimental and computational methods for assigning descriptors for the solvation parameter
  model.* J. Chromatogr. Open 2026, 10, 100366. *(not read)*
  <https://doi.org/10.1016/j.jcoa.2026.100366>
- Poole. *Wayne State University experimental descriptor database for use with the solvation parameter
  model.* J. Chromatogr. A 2020, 1617, 460841. <https://doi.org/10.1016/j.chroma.2019.460841>
- Kendall, Gal. *What Uncertainties Do We Need in Bayesian Deep Learning for Computer Vision?*
  arXiv:1703.04977. <https://arxiv.org/abs/1703.04977>

**Code and data (inspected directly)**

- DeepChem — <https://github.com/deepchem/deepchem> (MIT) ·
  <https://deepchem.readthedocs.io/en/latest/api_reference/models.html>
  - `deepchem/models/graph_models.py` — `uncertainty=True` log-variance head
  - `deepchem/models/torch_models/torch_model.py` — `predict_uncertainty`, MC-dropout combination
- SolProp_ML — <https://github.com/fhvermei/SolProp_ML> (CC-BY-4.0)
- chemprop_solvation — <https://github.com/fhvermei/chemprop_solvation> (MIT); `solvation_estimator.py`
- SoluteDB / model files — <https://zenodo.org/records/5792296> (CC-BY-4.0)
- RMG solvation web tool — <https://rmg.mit.edu/database/solvation/search/>
- Ulrich & Ebert 2022 supplement — <https://github.com/nadinulrich/solute_descriptor_prediction>
  (MIT; last commit 2021-11-07)
- PAULY — <https://i-am-pauly.com/> (2026 UFZ models; no licence published)
- AbraLlama — <https://huggingface.co/spaces/ttmn/AbraLlama> ·
  <https://github.com/BrightBlueCheese/AbraLLaMA>
- UFZ-LSER database — <https://www.ufz.de/lserd/>
