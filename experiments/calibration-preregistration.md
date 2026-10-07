# Pre-registration — calibration methodology and acceptance thresholds

Status: **committed before any calibration data exists.** Written under
#15 on 2026-08-16, when
the project held **zero measured retention data**. That is the point: a threshold chosen
after the result is a threshold chosen to be passed.

**Amendment rule.** Every number here is versioned and may be amended, but an amendment
requires a written reason recorded **before** the affected data is examined. An amendment
made after seeing a result is not an amendment; it is a new pre-registration that must say
so, and the original stands in the record beside it.

---

## 1. What is measured

**Two metrics, always reported together** (#27's rule, unchanged). Each hides a different
failure: pairwise reliability alone permits a model with badly wrong absolute retention that
happens to order well; coverage alone permits a model that brackets so widely it orders
nothing.

**Metric A — order reliability.** The model's headline output is `P(A before B)`. Calibration
is that number asked of itself: among pairs claimed at 0.9, do ~90% elute in the claimed
order? Reported as a reliability diagram plus expected calibration error over pairwise
probabilities. This grades the claim the user actually reads, not a proxy for it.

**Metric B — retention interval coverage and sharpness.** Coverage of the prediction interval
on `log k`, and its width. #19 needs resolution, not only order, and a model can be well
calibrated on order while badly calibrated on retention.

### `n` is counted in compounds, never in pairs

Pairwise statistics are computed over `n(n−1)/2` pairs, but pairs share compounds: one badly
predicted compound corrupts every pair it appears in. Measured (simulation, per-compound
error, compounds resampled per replicate):

| Compounds | Pairs | True SD | Naive binomial SD | Understated by | Effective `n` |
|---|---|---|---|---|---|
| 15 | 105 | 0.0410 | 0.0303 | 1.35× | 57 |
| 25 | 300 | 0.0283 | 0.0179 | 1.58× | 120 |
| 50 | 1,225 | 0.0182 | 0.0088 | 2.06× | 288 |
| 71 | 2,485 | 0.0149 | 0.0062 | **2.40×** | 432 |

⚠ The distortion **grows with study size** — the larger the compound set, the more a naive
interval flatters us.

**Rule: the compound is the unit of resampling.** Every accuracy and calibration figure
carries a compound-level bootstrap interval. A binomial interval on pair counts is not
acceptable in any report, summary or claim.

**Decision rule: a margin counts only if it exceeds 2 SD of the compound-level bootstrap**
of the *paired* difference (model and baseline scored on the same compounds). Required `n`
is a consequence of this rule and is **reported with every verdict**, passed or not.

---

## 2. Gates, and the third verdict

Baselines are #27's and are not revisited: **cold** (Crippen `logP`, no measurement) scores a
priori claims; **condition-transfer** (true order at one condition, reused) scores
post-reconciliation claims. Floors differ sharply — transfer sits at 76.5% (pH), ~92–94%
(modifier), ~99% (temperature).

| Axis | Gating? | Gating statistic | Why |
|---|---|---|---|
| **pH** | **gates** | conditional accuracy on the reordering subset | 23.5pp of headroom; the binary statistic is meaningful here |
| **Modifier** | **gates** | **relative error reduction** | ⚠ amends #27 — see below |
| Temperature | reported | both | 0.9–1.6pp reversal rate; no discrimination |
| Column (C18 vs C18) | reported | both | 4.6–6.2%; no discrimination |

### Why the modifier gate moved to relative error reduction

The modifier baseline is already 92–94%, leaving 6–8pp of total headroom. Measured, at a
genuine 17% error reduction and `n` = 71 compounds:

| Model/baseline error correlation | Conditional accuracy callable? | Relative error reduction callable? |
|---|---|---|
| 0.0 | no (1.7pp gain vs 3.1pp at 2 SD) | no (16.1% vs 20.3%) |
| 0.5 | no (2.6pp) | no (17.4%) |
| **0.7** | no (2.3pp) | **yes** (16.2% vs 14.8%) |
| 0.85 | yes (1.7pp) | **yes** (16.4% vs 10.4%) |

Relative error reduction wins at every correlation, because binarising continuous retention
into pairwise wins discards most of the signal. On an axis with little headroom the binary
statistic is simply the wrong instrument. Conditional accuracy is still **reported** on the
modifier axis — #27's always-together rule is untouched — it just no longer has authority to
fail the axis.

### Three verdicts, not two

**Pass**, **fail**, and **inconclusive**, the last defined here in advance:

> **Inconclusive** = the achieved 2 SD compound-level interval exceeds the observed gain.
> It does **not** block #19 and it does **not** count as a pass. It is reported as "this axis
> could not be resolved at the available sample size", **with the required sample size
> stated**.

Every measurement with limited power has an inconclusive outcome. The failure mode is
discovering it afterwards, when it is indistinguishable from a pass in a summary and gets
written up as one.

---

## 3. Thresholds

Marked **[derived]** where a floor or baseline forces the number, **[convention]** where it is
a chosen bar. Both are binding; the labels exist so amendments can be argued honestly.

### Neutrals, cold baseline — the a priori claim

| Quantity | Threshold |
|---|---|
| Conditional accuracy, pH reordering subset | **≥ 60%** **[convention]** — must beat coin-flip (50%) with a 2 SD margin; the transfer baseline scores exactly 0% here by construction |
| Aggregate pairwise accuracy | **must not fall below the transfer baseline** (76.5% pH) **[derived]** — #27's guard-rail against winning reorderings by scrambling stable pairs |
| Relative error reduction, modifier axis | **> 0 with a 2 SD margin** **[derived]** — the gate is the margin rule, not a fixed percentage |

### Ionisables, post-reconciliation — making #34's claim falsifiable

#34 claims: *we predict these weakly, we say so, we tell you which experiments to run, and
after those experiments we predict them well.* That is unfalsifiable without numbers. It is
now this:

> **After `n` = 3 scouting runs at a pH spacing of ≥ 2 units, ionisable compounds must reach
> the neutral cold-baseline thresholds above.** **[convention]**

`n` = 3 is taken from #16/#17's finding that data outruns the prior after 3–4 runs. If
ionisables after scouting cannot match neutrals before scouting, the two-regime capability
claim is withdrawn from the spec rather than restated.

### Calibration

| Quantity | Threshold |
|---|---|
| Coverage at nominal 90% | within **±5pp** **[convention]** — a sanity check only, see below |
| **Conformal inflation `λ`** | **≤ 2.0** **[derived]** — the primary calibration threshold |
| Rescope trigger | **`λ` > 2.5** **[derived]** |

**Coverage is not the test, and cannot be.** Conformal prediction reaches nominal coverage by
construction under exchangeability — widen far enough and coverage always arrives. A coverage
threshold is a test we cannot fail. What can fail is **sharpness**: the width needed to get
there.

Measured on the thin slice's stored ensemble (`prototype/thin-slice/retention_ensemble.npz`,
3 operating points × 1,000 scenarios × 15 usable compounds), scaling each compound's scenario
spread about its mean by `λ` and recounting pairs confidently ordered (P > 0.9):

| `λ` | op 0 | op 1 | op 2 |
|---|---|---|---|
| 1.0 (as sampled) | 73% | 70% | 67% |
| 1.5 | 58% | 40% | 42% |
| **2.0** | **41%** | **34%** | **37%** |
| 2.5 | 16% | 18% | 30% |
| 3.0 | 4% | 6% | 21% |
| 4.0 | 0% | 0% | 6% |

**The product does not degrade gracefully; it falls off a cliff between 2 and 3.** `λ` = 2.0
is where the confident-pair fraction reaches roughly half its sampled value — the point at
which the model stops beating the coin-flip-heavy regime it exists to replace.

⚠ **This is a bar we might genuinely fail.** #42 measured the ensemble under-dispersed by a
median factor of **1.23** [0.98, 1.47], IQR 0.71–1.66, **90th percentile 2.45**. If #45
returns an inflation near the top of that range, #19 rescopes. Pre-registering means we
agreed to that before knowing the answer.

⚠ **The figure moved from 1.21 on 2026-08-16, and it is important that it moved for the
right reason.** It was recomputed at n = 94 after a defect in the Table S-2 extract was
fixed — a section heading glued onto a compound name had been dropping one compound from
the join. **No threshold in this document was retuned**, and none may be: the thresholds
above were registered before any measurement existed and a corrected input is not a licence
to revisit them. Recorded here because a pre-registration whose numbers change silently is
worth nothing.

### What each verdict does

- **`λ` ≤ 2.0** — calibrated. Conformal recalibration is a legitimate, expected correction.
- **2.0 < `λ` ≤ 2.5** — calibrated but materially degraded. #19 proceeds **only** with the
  confident-pair fraction stated beside every design-space claim.
- **`λ` > 2.5** — the model is declared **under-dispersed, not calibrated**. #19 rescopes to a
  deterministic resolution map plus an optimiser, as its own body anticipates. The required
  inflation is the headline finding, not a footnote.

---

## 4. Strata

Count and boundaries were deferred here by #42. The answer is a **rule**, not a number,
because the count is data-determined:

- **Minimum stratum size**: hard floor `n` ≥ 19 (below this the 90% quantile is the sample
  maximum and the interval is not meaningful); practical floor `n` ≥ 50.
- **Below the floor, strata merge upward** along a fixed hierarchy: **confidence tertile
  merges first**, then φ range, then column type. **Ionisation regime never merges** — it has
  the largest measured effect (≈ 0.09 + 0.04 in `log k` for a neutral far from any pKa,
  against ≈ 0.36 from pKa error alone inside the ±1.5-unit danger band), and it is the axis
  #34's two-regime claim rests on.
- Confidence strata are cut on the **method-independent structural summary** (#42), never on
  the LSER-propagated spread, so a compound cannot change calibration class because the user
  changed column.
- **`λ` is computed per stratum.** A stratum with no calibration data cannot be conformally
  calibrated at all: it stays in block-plus-opt-in (#40) and **never enters the pooled
  figure**, so a thinly-populated stratum can neither drag the headline inflation nor hide
  behind it.

---

## 5. Floors no model can beat

Thresholds must sit above these, and a model that reaches them has finished, not failed:

- **Measurement repeatability** — METLIN SMRT, mean 36 s / median 18 s over 198 molecules
  across ≥ 30 days (#4).
- **Resolution uncertainty** — Knox `A, B, C` ship as engineering priors by particle class, at
  ±10–30% in plate count (#9, #13). **No improvement in retention prediction removes this**;
  it narrows only when measured van Deemter data enters through reconciliation.
- **LSER lack-of-fit** — WSU in-sample residual SD, median 0.030 (ACN) / 0.041 (MeOH) `log k`,
  itself a floor rather than an estimate (`sources/wsu-lser/README.md`).

## 6. What counts as evidence

- **Only in-domain predictions** enter any metric. Opt-in predictions outside the
  applicability domain produce values but never confidence, and their residuals are never
  evaluation data (#40's firewall).
- Every figure carries baseline, `n`, interval, axis, **reference identity and reference
  coverage** (#27 as extended by #31).
- Verbs follow the reference: **"agrees with"** for a model-output reference, **"accurate to"**
  only for measurement (#31).
