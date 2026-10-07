# #26 rescoped — `Fs` calibrated into an actual substitution error

`python3 calibrate_fs.py` reproduces everything; `results.txt` is the captured run.
numpy + stdlib.

## Why this exists

#26's primary experiment — fitting our own LSER constants to HSM3's 43,329 raw retention
measurements — **cannot run**, because that dataset is not obtainable (see below). But #26
lists a second deliverable that data in hand *can* supply:

> "it is the material needed to **calibrate `Fs` similarity into an actual confidence**
> rather than an uncalibrated heuristic"

**21 columns carry both** HSM parameters (`H, S*, A, B, C`) and WSU-2019 LSER system
constants (`c, e, s, a, b, v`). For each pair we can ask what a method developer actually
asks: *if I substitute column B for column A because `Fs` said they were similar, how wrong
is my retention prediction?*

## Verdict

**`Fs` ≤ 5 substitution is affordable on elution order — 0.9% of compound pairs flip, and a
close pair's separation moves 0.016, inside the 0.007–0.018 already budgeted for
topology-class substitution.** It is *not* affordable on absolute `log k` (median 0.109, 3–4×
the lack-of-fit floor), but `CONTEXT.md` is explicit that `log k` is not the axis this
project grades on, and the same caution cuts both ways: a `log k` figure must not be quoted
to *disqualify* a substitution any more than to bless one.

So `Fs`-based substitution earns the same standing as topology-class substitution:
**a degraded tier, never silent, quoted on order.** The measured table below is what makes
that quotable rather than assumed.

---

`Fs` comes from the HSM side; the substitution error comes from the LSER side, as RMS
Δ`log k` over the 94 WSU anchor compounds. **Neither number knows about the other**, so the
relationship between them is a genuine out-of-model test.

## The calibration

Pooled over 6 conditions (MeOH and ACN at 30/40/50% v/v), 1,200 column pairs:

| `Fs` band | pairs | median Δ`log k` | 90th pct | max |
|---|---|---|---|---|
| 0–5 | 18 | **0.109** | 0.152 | 0.171 |
| 5–10 | 132 | 0.120 | 0.239 | 0.311 |
| 10–20 | 411 | 0.151 | 0.301 | 0.651 |
| 20–40 | 405 | 0.228 | 0.390 | 0.495 |
| 40+ | 234 | 0.266 | 0.467 | 0.761 |
| **all** | 1,200 | 0.195 | 0.368 | 0.761 |

**`Fs` carries real signal** — median error rises monotonically across the bands, so it is
not noise. That is the honest positive finding, and it is what makes `Fs` legitimate for
*screening*.

**On `log k` it looks like it cannot support substitution**, for three reasons visible in the
same table — but see the correction below, because `log k` is not the axis this project
grades on:

1. **Even the closest band is expensive.** `Fs` ≤ 5 still costs a median 0.109 `log k`, which
   is 3–4× the LSER lack-of-fit floor (0.030 ACN / 0.041 MeOH, #15 §5). Choosing a column on
   `Fs` and treating its constants as a substitute imports an error several times larger than
   the model's own residual.
2. **The bands overlap heavily.** The 90th percentile of the *best* band (0.152) exceeds the
   *median* of the 10–20 band (0.151). A single `Fs` value therefore does not place a pair
   confidently.
3. **`Fs`-nearest is usually not LSER-nearest.** Picking the `Fs`-nearest column of ~20
   candidates gets you the one that actually predicts retention best only **1–3 times out of
   20**, median rank 4–7:

   | condition | columns | `Fs`-nearest = best | in best 3 | median rank |
   |---|---|---|---|---|
   | MeOH 30% | 21 | 2/21 | 7/21 | 7 |
   | MeOH 40% | 21 | 3/21 | 8/21 | 5 |
   | MeOH 50% | 21 | 2/21 | 6/21 | 5 |
   | ACN 30% | 20 | 1/20 | 6/20 | 6 |
   | ACN 40% | 20 | 2/20 | 8/20 | 4 |
   | ACN 50% | 20 | 3/20 | 9/20 | 4 |

This is the **measured** version of #36's *decomposition* argument (80.5% of `Fs²` is the C
term, which has no LSER analogue). Decomposition and measurement agree.

> ### ⚠ But `log k` is the wrong axis to judge this on — corrected below
>
> `CONTEXT.md` on topology-class substitution is explicit: the cost is "SD ≈ 0.055–0.145 in
> `log k`, but only ≈ 0.007–0.018 on a close pair's separation — affordable **only** because
> elution order is the primary metric, and **never to be quoted as affordable on `log k`**".
>
> The table above is that same misleading half, pointed the other way: quoting `log k` makes
> `Fs` substitution look *disqualifying* exactly as it would make topology-class substitution
> look disqualifying. The next section measures the axis the project actually grades on, and
> it changes the conclusion.

## The elution-order cost — the axis that decides it

Fraction of the 4,371 compound pairs whose *predicted order flips* between the two columns,
and how far a **close** pair's separation moves (pairs within 0.10 `log k` on column A):

| `Fs` band | pairs | median order flip | 90th pct | median \|Δsep\| on close pairs |
|---|---|---|---|---|
| **0–5** | 18 | **0.9%** | 2.0% | **0.016** |
| 5–10 | 132 | 1.5% | 3.0% | 0.025 |
| 10–20 | 411 | 1.9% | 4.1% | 0.034 |
| 20–40 | 405 | 2.3% | 4.3% | 0.036 |
| 40+ | 234 | 4.7% | 10.6% | 0.082 |
| all | 1,200 | 2.3% | 5.7% | 0.039 |

**`Fs` ≤ 5 costs 0.016 on a close pair's separation — inside the 0.007–0.018 already budgeted
for topology-class substitution.** On the primary metric, substituting a column chosen at
`Fs` ≤ 5 is about as expensive as the degraded tier the architecture already accepts and
labels. It is affordable **on the same terms**: degraded, never silent, and never quoted as
affordable on `log k`.

The monotonicity is cleaner here than in the `log k` table, and the top band is genuinely
bad — 4.7% of pairs flipping, 10.6% at the 90th percentile. So `Fs` does discriminate
usefully; it simply has to be read on the axis the product is graded on.

⚠ **Not comparable to #27's measured reversal rates.** `CONTEXT.md` quotes 4.6–6.2% for
column C18-vs-C18 reversal, which is *measured* retention on real compounds. The numbers here
are *predicted*-order flips between two LSER models over 94 neutral anchors. They answer
"how much does my prediction change if I substitute", not "how often does nature reorder".
The two must not be set against each other.

### What this means for the `Fs`-nearest result

The ranking failure above (`Fs`-nearest is the LSER-nearest only 1–3 times in 20) is real but
**matters less than it first looks**. If the chosen column lands in a low `Fs` band, being
5th-best rather than best still costs ~1% of order flips. The ranking is weak; the
consequence of that weakness is bounded.

### `Fs` is weakest in acetonitrile, which is where we work

Correlation between `Fs` and substitution error:

| modifier | φ | Pearson | Spearman |
|---|---|---|---|
| methanol | 30 | 0.584 | 0.482 |
| methanol | 40 | 0.477 | 0.468 |
| methanol | 50 | 0.479 | 0.455 |
| **acetonitrile** | 30 | **0.271** | 0.420 |
| **acetonitrile** | 40 | **0.282** | 0.447 |
| **acetonitrile** | 50 | **0.302** | 0.431 |

`Fs` is roughly **half as predictive in acetonitrile as in methanol**. That is exactly what
#36's decomposition predicts: `Fs` is dominated by the C (cation-exchange) term, which is
acetonitrile-specific and has no LSER analogue — so in ACN, most of what `Fs` measures is
precisely what LSER cannot see. Two independent lines of evidence, same conclusion.

## ⚠ A confirmed arithmetic error in the unmerged bridge research

`research/hsm-lser-bridge` (unmerged branch) reports **"against XBridge C18, 548 of 819
columns score `Fs` ≤ 3"** and prints a top-8 with values 0.14–0.28. Using the published
weights I get **18 of 819**, same columns at 0.93–2.20.

**The cause is identified.** That document computed

```
Fs = sqrt( 12.5·ΔH² + 100·ΔS*² + 30·ΔA² + 143·ΔB² + 83·ΔC² )       <- weights on the squares
```

where the published definition is

```
Fs = sqrt( (12.5·ΔH)² + (100·ΔS*)² + (30·ΔA)² + (143·ΔB)² + (83·ΔC)² )   <- weighted, then squared
```

This reproduces its printed values to three decimals on every column in its list:

| column | doc | `W·d²` variant | correct |
|---|---|---|---|
| Hypersil Beta Basic-18 | 0.14 | **0.140** | 1.31 |
| Discovery C18 | 0.18 | **0.179** | 0.93 |
| Hypurity C18 | 0.20 | **0.195** | 1.68 |
| ChromCore 120 C18 | 0.24 | **0.244** | 2.14 |
| Flowrosil ODS | 0.25 | **0.247** | 1.97 |
| Poroshell 120 EC-C18 | 0.28 | **0.279** | 1.96 |

Because the differences are small, squaring them *before* weighting shrinks every term, so
the reported `Fs` values are far too low and far too many pairs clear the ≤ 3 threshold.

**The downstream conclusion inverts.** That research concluded from 548/819 that "the ≤3
equivalence threshold is not discriminating for a mainstream C18 — the commercial C18 market
is genuinely crowded". At **18 of 819** it is highly discriminating, and `Fs` ≤ 3 is a
meaningful shortlist rather than a rubber stamp.

### #36's decomposition is **not** affected

Worth checking separately, because `CONTEXT.md` cites it. #36 reports H 0.4%, B 12.8%,
C 80.5% of `Fs²`. Recomputing both ways:

| basis | H | S* | A | B | C |
|---|---|---|---|---|---|
| **correct formula**, all 334,971 pairs | 0.3% | 2.5% | 3.3% | 8.0% | **85.9%** |
| correct formula, vs XBridge C18 | 0.5% | 3.7% | 3.3% | 8.2% | 84.3% |
| `W·d²` slip, all pairs | 1.9% | 2.0% | 8.9% | 4.5% | 82.7% |
| **#36 as published** | **0.4%** | — | — | **12.8%** | **80.5%** |

#36's H figure (0.4%) matches the correct formula (0.3–0.5%) and not the slip (1.9%), so
**#36 used the published definition.** Its exact percentages are basis-dependent — the
document does not pin down whether it decomposed over all pairs or against a reference — but
the claim `CONTEXT.md` rests on is unaffected: **C dominates `Fs²` at 80–86% however it is
computed, and H contributes well under 1%.**

So the error is confined to the bridge document's `Fs` *values* and its 548/819 count.

## The acquisition question is closed

The HSM3 dataset was searched for properly, not just once:

| Route | Result |
|---|---|
| ChemRxiv public API | Cloudflare challenge |
| ChemRxiv DOI (`10.26434/chemrxiv-2024-mt2fp`) | 403 |
| Crossref relations | only `has-preprint`; no dataset relation |
| **OpenAlex** | **OA status `closed`; all 3 locations non-OA** |
| UvA-DARE (Pirok is a co-author) | record exists, `is_oa: false` |
| Zenodo, figshare | no deposit found |
| hplccolumns.org | hosts `database.csv` only — no retention data |

The article is closed access and no open deposit of the 43,329 measurements exists. The only
remaining routes are institutional access to the Elsevier SI, or **contacting a co-author** —
Dwight Stoll (who already maintains the HSM database and is the natural contact for its
licence question) or Bob Pirok at UvA.
