# Stock SoluteML vs measured WSU descriptors — the bias measurement (#35)

Session 2026-08-16. Scope as reduced by #22's closure and confirmed by the owner:
**bias measurement only, no re-fit.** No live chromatograms entered anything (#45 stays
paused); every number here comes from files already in the repo plus the CC-BY SoluteDB
corpus downloaded from Zenodo 5792296.

Reproduce: `env/ionisation-venv/bin/python measure_bias.py` (numpy only, ~40 s);
full console output in `run.txt`, per-compound rows in `residuals.csv`.

> **Public copy:** `residuals.csv` is not included because its `meas_*` columns reproduce
> WSU-2019 Table S-2 (© Elsevier). Rebuild `sources/wsu-lser/` (see its README), then rerun
> `measure_bias.py` to regenerate it.
Every figure is emitted through `prototype/eval-harness`'s `Figure`, so none can appear
without its baseline, `n`, interval, condition axis, reference identity and reference kind.

## The finding that constrains every other one, stated first

**All 94 WSU compounds are inside SoluteML's training corpus.** #42 inferred this from
provenance; SoluteDB is CC-BY and inspectable, so `overlap_solutedb.py` settles it row by
row: **93 of 94 match on SoluteDB's own key** (non-standard fixed-H InChI, string-equal),
and the 94th (cinnamyl alcohol) matches on the standard InChIKey skeleton — a stereo
difference, still the same compound to a model. **Zero clean compounds.**

Two consequences, neither of which the numbers below can escape:

1. **This is a floor on the bias, not an estimate of it.** Every figure here is measured on
   compounds the network was fitted to. The honest reading of every interval below is
   "at least this much".
2. **#42's handover to this ticket is not executable with this compound set.** #42 asked
   #35 for a *clean* re-run of the 2.3× confidence stratification and of #46's
   structural-deficiency question. There is no clean split to be had here — the question
   needs compounds outside SoluteDB, which means RepoRT's scaffold split (#29) or new
   measurement, not this anchor. **Recorded as still open, not as answered.**

## 1. The headline: model bias vs measured WSU (n = 94 compounds)

Signed bias is `ensemble mean − measured`; positive means SoluteML predicts high.
Reference is measurement, so per #31 these are entitled to "accurate to".

| descriptor | bias [95% CI] | MAE [95% CI] | published random-split MAE | published substructure MAE |
|---|---|---|---|---|
| E | **+0.035** [+0.004, +0.065] | **0.106** [0.086, 0.130] | 0.041 | 0.084 |
| S | **+0.061** [+0.027, +0.097] | 0.108 [0.079, 0.141] | 0.098 | 0.170 |
| A | −0.004 [−0.020, +0.011] | 0.054 [0.042, 0.066] | 0.038 | 0.067 |
| B | **+0.018** [+0.006, +0.031] | 0.040 [0.031, 0.051] | 0.047 | 0.088 |

Intervals are compound-level bootstrap (#15 §1). Published MAEs are Chung 2022 SI
Fig. S11/S12, scored against SoluteDB itself.

- **Three of four descriptors carry a bias whose interval excludes zero** — E, S and B all
  positive. `A` is the one clean descriptor. So there *is* a systematic offset, which is
  what #31 made this ticket the sole instrument to detect, and its sign is consistent:
  SoluteML runs high.
- **E is the failure.** Its MAE against measurement (0.106) is **2.6× the published
  random-split number and larger than the published substructure-split number (0.084)** —
  on compounds it trained on. The error a model shows on novel scaffolds against its own
  corpus is smaller than the error it shows on memorised compounds against an independent
  measurement. E was already flagged as both worst-predicted and most contaminated; this
  is that concern arriving as a number.
- **B is the reassuring one** — MAE 0.040 beats its own published random-split figure.

## 2. Where the error comes from — the decomposition

Because SoluteDB ships its labels, `model − measured` splits exactly:

```
model − measured  =  (label − measured)   +   (model − label)
                      corpus label bias        network fit error
```

| descriptor | RMS model err | RMS label err | RMS fit err | corr(label, fit) | label share of variance |
|---|---|---|---|---|---|
| E | 0.154 | 0.146 | 0.067 | +0.11 | 0.92 |
| S | 0.186 | 0.189 | 0.060 | −0.37 | 1.14 |
| A | 0.080 | 0.092 | 0.039 | −0.39 | 1.17 |
| B | 0.063 | 0.065 | 0.034 | −0.32 | 1.12 |

**The error is the corpus, not the network.** RMS label error ≈ RMS model error on all four
descriptors, while the network's own deviation from its labels is 2–3× smaller. SoluteML
reproduces SoluteDB faithfully; SoluteDB is what disagrees with the Wayne State
measurement. The negative `corr(label, fit)` on S, A and B says the network mildly shrinks
its labels back toward the truth — which is why the label share sits slightly above 1.

**This inverts #35's stated failure branch.** The ticket said that if the re-fit failed,
the reading would be "contamination is baked into the representation rather than the
output layer". The decomposition says the opposite is more likely: contamination lives in
the *labels*, and the representation tracks them well. An output-layer re-fit is therefore
attacking the right layer — it is blocked by having 94 compounds, not by being aimed
wrongly. That is an argument for reopening the re-fit if Poole 2020 is ever obtained, and
against treating #31's joint-validation fallback as the presumptive answer.

Corpus label bias with intervals: E −0.025 [−0.056, +0.004], S +0.023 [−0.013, +0.060],
**A −0.032 [−0.049, −0.015]**, B −0.010 [−0.022, +0.004]. Only `A` is individually
resolved, but the *magnitudes* (RMS column above) are what carry the finding.

## 3. Shape of the residual distribution

| descriptor | p0 | p5 | p25 | p50 | p75 | p95 | p100 | sd | skew | \|e\|>2sd |
|---|---|---|---|---|---|---|---|---|---|---|
| E | −0.703 | −0.202 | +0.016 | +0.047 | +0.081 | +0.239 | +0.480 | 0.150 | −1.19 | 6 |
| S | −0.679 | −0.094 | +0.014 | +0.045 | +0.108 | +0.294 | +1.113 | 0.177 | +1.25 | 3 |
| A | −0.271 | −0.159 | −0.034 | +0.021 | +0.026 | +0.088 | +0.218 | 0.080 | −0.97 | 7 |
| B | −0.229 | −0.064 | −0.005 | +0.017 | +0.041 | +0.079 | +0.346 | 0.061 | +0.99 | 4 |

**These are not Gaussian and should not be sampled as if they were.** Every descriptor is
skewed (|skew| ≈ 1.0–1.25) with a tight core and one long tail: the interquartile range on
`E` is 0.065 wide while its full range spans 1.18. A normal error model fitted to the SD
would simultaneously overstate the typical compound and understate the tail — the exact
failure mode that makes a calibrated `P(A before B)` uncalibrated in both directions.

Worst compounds are chemically coherent, not random: **pentafluorophenol** (E −0.703, and
the worst compound in retention space), **progesterone** (S +1.113), **benzenesulfonamide**
(S −0.679, B +0.346), **nicotinamide** (S −0.532), **8-hydroxyquinoline** and
**2-aminobiphenyl** (E +0.480, +0.448). Polyfluorination, steroids, sulfonamides and
amphoteric/chelating heteroaromatics — recognisable classes, which is the encouraging part
for #46's register even though the clean test of it remains undone.

**Proportional bias — the errors shrink the descriptor range.** OLS of residual on measured
value gives slope **−0.096 [−0.199, −0.007] for E** and **−0.163 [−0.207, −0.109] for A**,
both excluding zero (S and B do not resolve). SoluteML over-predicts small values and
under-predicts large ones, i.e. it compresses the descriptor axis toward its mean. Pushed
through the LSER this compresses predicted `Δlog k` between compounds, which makes elution
order look *more* decidable than it is. Direction matters for #14/#15: it is a
confidence-inflating error, not a neutral one.

## 4. Retention space (#41's criterion)

Predicted descriptors and measured descriptors pushed through the same published constants,
over the 332 unmasked (column × modifier × φ) fits, wetting-masked rows excluded per #40.

- **Mean per-compound RMSE 0.085 log k [0.072, 0.098]**, against Poole's harmonised
  descriptor-set median LSER residual SD of 0.041 (MeOH) / 0.030 (ACN). **Descriptor error
  contributes roughly twice the residual the LSER fit carries on its own.** The QSPR layer,
  not the LSER layer, is the dominant error source in this stack.
- **Mean per-compound bias −0.051 log k [−0.066, −0.036]** — a systematic *under*-prediction
  of retention. Sign flips from descriptor space because `s`, `a` and `b` are negative in
  RPLC, so predicting high on S and B predicts weak retention.
- Distribution is long-tailed here too: median 0.068, p75 0.109, p95 0.188, max 0.387
  (pentafluorophenol). Worst five: pentafluorophenol 0.387, progesterone 0.289,
  benzenesulfonamide 0.283, nicotinamide 0.237, diphenylamine 0.194.

| modifier | cells | model RMSE | label RMSE | model bias | Poole residual SD |
|---|---|---|---|---|---|
| acetonitrile | 14,382 | 0.093 | 0.085 | −0.046 | 0.030 |
| methanol | 15,698 | 0.117 | 0.110 | −0.055 | 0.041 |
| tetrahydrofuran | 1,128 | 0.093 | 0.085 | −0.040 | not published |

Methanol is worse than acetonitrile on both, matching the ordering Poole reports for his
own fits.

## 5. The weak-base stratum — a class effect, with an interval

Table S-2 labels its last 15 compounds "compounds designated as weak bases". Splitting on
the anchor's own classification:

| class | n | bias E | bias S | bias A | bias B | ret RMSE | ret bias |
|---|---|---|---|---|---|---|---|
| neutral | 79 | +0.011 | +0.065 | −0.011 | +0.023 | 0.087 | −0.063 |
| weak base | 15 | **+0.162** | +0.040 | +0.031 | −0.011 | 0.073 | +0.017 |

Contrast, resampled within each class so both arms stay populated:
**E +0.152 [+0.077, +0.236]**, A +0.042 [+0.007, +0.082], B −0.034 [−0.065, −0.005],
S −0.025 [−0.075, +0.029].

**The `E` bias is essentially a weak-base effect** — 15× larger on bases than on neutrals,
and the neutral arm's `E` bias (+0.011) would not be distinguishable from zero on its own.
Anilines, aminobiphenyls and hydroxyquinolines are where SoluteML's excess-molar-refraction
prediction goes wrong, systematically and in one direction.

This is a *different* offender from the one #31 predicted. #31 expected `B` to be the
systematic problem for drug bases, from coverage arguments (the anchor's `B` range does not
reach atenolol or propranolol). On the anchor's own bases the resolved offender is `E`, and
`B` moves the other way. Both can be true — #31's argument is about extrapolation beyond
the anchor's range, this is a bias *inside* it — but the two should not be conflated, and
the drug-base case still has no measurement behind it.

## 6. Ensemble dispersion — the ensemble does not know how wrong it is

Retention-space RMSE over propagated ensemble SD, per compound: **median 1.23
[0.98, 1.47]**, p25 0.71, p75 1.66, p95 2.92, max 5.09. This reproduces #42's 1.21 at the
corrected n = 94.

In descriptor space the shortfall is far worse: RMSE / median ensemble SD is **3.95 (E),
3.70 (S), 3.85 (A), 2.15 (B)**. The 25-member spread captures roughly a quarter of the
actual descriptor error. Both numbers are floors — measured on seen compounds — so the
conformal inflation `λ` that #15 pre-registered will have to do real work, and a `λ` near
1 from any in-corpus evaluation should be disbelieved on principle.

## 7. A data defect found on the way, and what it did

`sources/wsu-lser/wsu2019-descriptors-s2.csv` row 81 carries the table's section heading
glued onto the compound name: `Compounds designated as weak bases 2-Aminobiphenyl`. Any
consumer joining on compound name silently drops that compound and runs at **n = 93 of
94**.

- `prototype/descriptor-covariance/measure1_ensemble_covariance.py` and
  `prototype/thin-slice/slice.py` already strip the heading.
- `prototype/confidence-stratum/measure_confidence_stratum.py` does not — so **#42's
  numbers were computed on 93 compounds**. It reported "93 of 94 compounds matched"
  honestly, so nothing there is wrong, but the reason was an extraction defect rather than
  a genuine absence, and it is fixable.
- `experiments/fs-calibration/`, `experiments/structural-deficiency/` and
  `prototype/phi-vs-lss/` also read the file and were not audited here.

The dropped compound is **2-aminobiphenyl — a weak base, and the third-worst `E` error in
the set (+0.448)**. It sits in the tail of the exact effect §5 reports, so the defect was
biasing the very stratum most worth measuring. `measure_bias.py` normalises the name at
load rather than editing a provenance-controlled extract; **the durable fix is a `class`
column in the S-2 CSV**, which #43 wants anyway for its structural boundary, and which
would remove the heading-in-a-name hazard for all ten consumers.

## 8. The residual sampler for #14 — `residual_sampler.py`

Built against an earlier design proposal (not in this repo), which chose a stratified empirical
quantile sampler over parametric alternatives. That choice is right — §3's skew rules out
Gaussian and the asymmetry rules out plain Student-*t*. Three things in it do not survive
contact with the data, and the implementation corrects each.

**The proposal's quantile tables are not the measured ones.** It used each class's *mean
bias* as that class's *median*: `E_neutral` p50 given as +0.011 (the neutral mean) against a
measured +0.043, and `E_base` p50 as +0.162 (the base mean) against a measured +0.141. The
other anchor points are interpolations rather than measurements — `E_base` p0 given as
−0.150 against a measured −0.054, p5 as −0.050 against −0.009. Regenerated tables are in
the module and reproduce from `residuals.csv` on every run.

**Independent per-descriptor sampling is the larger error, and it is the one #14 named.**
The proposal draws each descriptor's residual separately. The measured across-compound
residual correlations are not small: **corr(S, B) = −0.67**, corr(A, B) = −0.26,
corr(E, A) = +0.24. Propagated through the 332 unmasked column vectors, **independent
sampling gives Var(log k) 2.8× the joint value — sd overstated 1.67×**, inflating every
interval and dragging every `P(A before B)` toward 0.5. The module uses a Gaussian copula
with the measured marginals; because these marginals are far from normal the copula
attenuates its input correlation, so the copula matrix is pre-calibrated (`calibrate_copula`)
to make the *sampled* correlation equal the measured one. Self-check: worst cell 0.003, and
the propagated variance ratio lands at 0.359 against the measured 0.354.

Worth noting against #37: this error correlation is **much stronger than the SoluteML
ensemble's own** within-compound correlation (E–S +0.29, everything else |ρ| < 0.06, giving
a variance ratio of 0.89–0.96). The ensemble misrepresents both the size of its error and
the shape of its correlation.

**A seven-point table from 15 weak bases is reading noise into a shape.** `E_base` p95 is
+0.457 with a bootstrap CI of [+0.248, +0.480]; p50 is +0.141 [+0.070, +0.231]. The class
*mean* is resolved (§5) but the class *shape* is not, so the module spends the one degree of
freedom the data supports: **pooled shape, shifted by the measured class offset**. This is
also the conservative direction — the 15 bases show no left tail (min −0.054 against the
pool's −0.703), so bases are sampled wider than observed rather than narrower.

Two smaller points. The tables are bounded by single compounds at each extreme — the pooled
`S` maximum (+1.113) is progesterone alone, and resampled it ranges over [+0.364, +1.113] —
so the module extends the outer segments' slope rather than clipping, with a fitted GPD tail
named as the proper instrument once more compounds exist. And the `A` residual has a
near-zero IQR span (p25 −0.034, p75 +0.026) with tails an order of magnitude wider; anything
fitted to its SD will be wrong at both ends.

### The compression correction: the de-shrink does not earn its place

The proposal also de-shrinks the point prediction before adding a residual. Tested on #41's
criterion, over the same 332 fits, as relative RMSE reduction against stock:

| correction | relative RMSE reduction | harness verdict |
|---|---|---|
| intercept debias only | **+0.179** [+0.094, +0.261] | PASS (2 SD margin 0.086) |
| debias + slope de-shrink | +0.042 [−0.065, +0.146] | **INCONCLUSIVE** — needs n = 611 |

**Adding the slope term destroys three-quarters of the gain the debias delivers.** The
de-shrink widens the residual it is supposed to fix — `E` residual sd goes 0.150 → 0.162,
`S` 0.177 → 0.186 — because expanding the axis expands the errors with it. On the quantity
the proposal invoked it for, pairwise elution order, it does nothing measurable: order
agreement 0.9610 → 0.9600 overall, 0.7145 → 0.7155 on close pairs (|Δlog k| < 0.1).

There is also a sequencing bug independent of whether it helps: the proposal de-shrinks the
prediction and then adds a residual drawn from the distribution of errors of the
*un*-de-shrunk prediction. Those are different distributions, and adding them double-counts.

**Neither correction is implemented here, and not only because of the numbers.** #23's
firewall makes the bias measurement diagnostic metadata that "must never become adjusted
training labels or a silent correction". A de-shrink applied inside the sampler is exactly a
correction derived from the bias measurement — one fitted on 94 in-corpus neutral solutes and
then applied to chemistry outside the anchor's range, which §"does not change" forbids. The
+0.179 debias gain is real and may well be worth taking, but it is an **owner decision about
the firewall**, not something a prototype should quietly switch on.

## What this changes, and what it does not

**Changes:**

1. **The uncertainty budget (#14) gets a measured systematic component**, per descriptor
   with intervals, plus its retention-space consequence (0.085 log k RMSE, −0.051 log k
   bias). Per #23's firewall this is diagnostic metadata: it is reported, never subtracted
   as a correction and never turned into training labels.
2. **The QSPR layer is confirmed as the dominant error source** — ~2× the LSER fit's own
   residual SD. Effort spent on LSER refinement before descriptor quality is misallocated.
3. **The residual model should not be Gaussian and should not be uniform across classes.**
   Skew ≈ ±1, and a resolved weak-base `E` offset.
4. **The re-fit target is the output layer, and the evidence now supports that** — the
   error is in the labels, which is what an output-layer re-fit replaces.

**Does not change:**

- **No claim here extends to drug-like chemistry.** 94 neutral small solutes, median V 1.06,
  no acids by name; #31 showed the thin-slice drugs sit outside this anchor's range on `B`
  and `S`. The verb "accurate to" is earned only on this chemistry.
- **Contamination is unresolved, and cannot be resolved with this compound set.** Every
  figure is a floor. #42's two handed-over questions stay open.
- **Nothing here touches elution order or the #27 baselines.** Retention space here means
  descriptors through fixed published constants, not measured `k`. That remains #45's job.

## Recommendation

1. **Carry the four descriptor biases as an explicit systematic term in #14**, with the
   floor caveat attached to each, and model the residual as skewed rather than normal.
2. **Add a `class` column to the S-2 extract** and fix the `confidence-stratum` join. Small,
   and it removes a hazard that has already silently moved an n.
3. **Reopen the output-layer re-fit as the preferred escape if Poole 2020 arrives** — the
   decomposition says the layer is the right one. At n = 94 it stays inconclusive under
   #15's third verdict, not failed.
4. **Route #42's clean-split questions to RepoRT's scaffold split (#29)**, not to this
   anchor. They cannot be answered on compounds the model has memorised.
5. **Treat the weak-base `E` offset as the first concrete candidate for #46's register** —
   it is a named class with a resolved, one-directional bias. It is not yet the clean test
   #46 needs, because these compounds are in-corpus.
