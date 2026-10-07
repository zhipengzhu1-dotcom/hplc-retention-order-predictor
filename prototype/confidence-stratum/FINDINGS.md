# Confidence stratum from ensemble spread — findings (#42)

Session 2026-08-15, branch `prototype/confidence-stratum` (cut from
`prototype/descriptor-covariance`, whose 25-member SoluteML predictions this
reuses). Two scripts, both runnable from the repo root with no arguments and no
network:

- `measure_confidence_stratum.py` — does ensemble spread stratify error, in
  descriptor space and in retention space?
- `measure_rank_stability.py` — is a compound's confidence a property of the
  compound, or of the method?

All errors in `log₁₀ k` (WSU convention). Wetting-masked fits excluded per #40
(332 of 354 kept). Grading follows #41: retention space is the criterion,
descriptor rmse is diagnostic. Truth is WSU Table S-2 (the anchor lineage, #23);
all 94 compounds matched to the prediction set.

⚠ **Re-run 2026-08-16 at n = 94.** The original run matched 93 of 94 because Table S-2's
section heading was glued onto a compound name in the extract, so a name join dropped
2-aminobiphenyl. Fixed at the producer (`sources/wsu-lser/extract_wsu.py`, which now emits
a `class` column). Every number below is the n = 94 re-run; changes are small and noted
where they occur. The recovered compound is a weak base and it **enters the §4 table**,
which is the material consequence.

## Read this before quoting any number here

**This is a contaminated evaluation.** All 94 WSU compounds are inside SoluteDB, the
corpus SoluteML was trained on — measured row by row in #35 (93 on SoluteDB's own fixed-H
InChI key, 1 on InChIKey skeleton), no longer inferred. Absolute rmse is therefore
optimistic and these figures are **diagnostic metadata under the #23 corpus
fence — not calibration evidence and not a success claim**. The quantity the
decisions rest on is the *ratio* between strata, which is far more robust to
contamination than the level, because both strata are contaminated equally.
A clean re-run was handed to #35 and **is not executable**: no corpus examined pairs
out-of-corpus compounds with a usable label (`experiments/solutedb-census/`).

## 1. The stratification is real, and smaller in retention space than in descriptor space

Descriptor space, rmse by tertile of that descriptor's own ensemble SD
(n = 32/31/31):

| Descriptor | rmse all | low-SD third | high-SD third | ratio |
|---|---|---|---|---|
| E | 0.154 | 0.119 | 0.206 | 1.7× |
| S | 0.186 | 0.074 | 0.270 | 3.6× |
| A | 0.080 | 0.025 | 0.117 | 4.7× |
| B | 0.063 | 0.033 | 0.088 | 2.6× |

Ulrich et al. 2026 report the same *shape* (A strongest, E weakest) on a
different descriptor lineage and a different model family — an independent
corroboration of the mechanism, not of the magnitudes.

Retention space, per-compound rmse over all 332 unmasked fits (n = 32/31/31):

| Stratified by | rmse all | low third | mid | high third | ratio |
|---|---|---|---|---|---|
| ensemble spread propagated through the LSER | 0.106 | 0.068 | 0.080 | 0.151 | **2.2×** |
| plain sum of the four descriptor SDs | 0.106 | 0.057 | 0.101 | 0.143 | **2.5×** |

**4.7× on `A` becomes 2.2× in retention space.** `a` is one coefficient among
six, so a descriptor-space win is reweighted before it reaches elution order —
#41's decoupling principle, second sighting. Quoting "nearly tenfold on A" as a
retention claim would overstate by roughly 4×.

The two thirds straddle the ~0.09 budget line (#2), which is what makes the axis
decision-bearing rather than merely interesting: one pooled inflation of 0.106
over-widens the confident third and under-widens the rest.

The *plain* descriptor-SD summary stratifies marginally better than the
LSER-propagated one, and is method-independent by construction — which is why
#42 cuts the calibration stratum on it.

## 2. Who is confident is structural; how confident is method-dependent

Spearman ρ of the per-compound spread between pairs of methods (3,000 random
pairs of the 332): **median 0.942**, 5th percentile 0.787, worst 0.634. The
deliberate extreme, methanol 10% vs acetonitrile 70%: ρ = 0.923.

Median spread magnitude across methods ranges 0.016 (XTerra MS C18, THF 70%) to
0.088 (SunFire C18, THF 20%) — a **5.6×** span.

This decomposition is what lets #42 be cheap: the ranking is computed once from
structure and does not multiply into the method strata, while the magnitude
falls out for free when the descriptor covariance is pushed through that
method's LSER coefficients.

## 3. Raw ensemble SD is a shape, not an interval

Ratio of retention-space rmse to propagated ensemble SD, per compound (n = 94):
median **1.23**, IQR 0.71–1.66, 90th percentile 2.45, max 5.09. (#35 reports the same
statistic with a compound-level bootstrap interval: 1.23 [0.98, 1.47].)

The ensemble under-states its own error by ~23% at the median and ~2.5× in the
upper tail. It is usable as a relative shape; it is not a calibrated uncertainty,
which is exactly the gap #15's conformal step exists to close.

## 4. Confident-and-wrong is a distinct failure the ensemble cannot see

Highest error-to-spread compounds:

| Compound | spread | rmse | ratio | largest descriptor miss |
|---|---|---|---|---|
| Pentafluorophenol | 0.076 | 0.387 | 5.1 | ΔE −0.70 |
| Benzenesulfonamide | 0.067 | 0.283 | 4.2 | ΔS −0.68, ΔB +0.35 |
| Progesterone | 0.073 | 0.289 | 4.0 | ΔS +1.11 |
| 2-Naphthaldehyde | 0.046 | 0.172 | 3.8 | ΔE −0.39 |
| Nicotinamide | 0.071 | 0.237 | 3.3 | ΔS −0.53 |
| 4-Nitrobenzyl alcohol | 0.043 | 0.116 | 2.7 | ΔB +0.07 |
| **2-Aminobiphenyl** | 0.044 | 0.111 | 2.5 | **ΔE +0.45, ΔA +0.22** |

**2-aminobiphenyl is the compound the extract defect was dropping**, and it lands in this
table on recovery — a weak base whose largest miss is `E`, which is exactly the
class-conditional bias #35 later resolved (`E` bias +0.152 [+0.077, +0.236] on weak bases
versus +0.011 on neutrals). A defect that removes a compound from the confident-and-wrong
list is worse than one that merely adds noise.

**This table on its own proves little** — with 94 compounds some ratio is always
the maximum, and a perfectly calibrated model produces a tail too. What makes the
class real is that three independent sources named it before we looked: Ulrich's
deficiencies *of the Abraham model* (intramolecular H-bonding, site-specific and
sterically hindered donation, the several interactions `S` collapses into one),
WSU-2019 Table S-1 independently dropping bulky compounds and some H-bond acids
from its own fits, and the structural coherence of the specific failures here
(polyfluoroaromatic polarisability, large rigid polycyclics, sulfonamide
dipolarity). The tail is consistent with the class; it is not the evidence for it.

By construction **no ensemble statistic can flag these**: the 25 networks agree
because they learned the same inadequate representation from the same corpus.
The flag must come from outside the ensemble — see the structural-deficiency
register ticket.
