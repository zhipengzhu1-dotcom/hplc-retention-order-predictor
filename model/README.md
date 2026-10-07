# The real ionisation provider, wired in

`env/ionisation-venv/bin/python model/calibrate_ionisation.py` — results captured in
`calibration-results.txt`. Provider is Uni-pKa (#7's open default, stood up in #20).

## Why

The thin slice **hardcodes literature pKa** and derives `f_neutral` from a closed-form
sigmoid chosen by a hand-assigned `type`. Meanwhile the slice's own variance attribution puts
**`pKa` at 0.67** of the variance for partially ionised compounds — so the dominant term for
the in-scope chemistry was a constant somebody typed in. This replaces it with the provider
the architecture already chose.

## 1. The model was carrying twice the pKa uncertainty it needed

Uni-pKa against every cited literature value in `compounds.py` — the same compounds, the same
citations:

| | |
|---|---|
| n | 10 pKa values, 9 compounds |
| mean error | **+0.003** (essentially unbiased) |
| **SD** | **0.170** |
| MAE | 0.127 |
| max abs error | 0.38 (imipramine) |

**The slice assumes ±0.35 (1 SD). The measured provider error is 0.170 — under half.**

That matters because pKa owns 0.67 of the variance for partially ionised compounds: the model
has been inflating its dominant uncertainty term roughly two-fold, which drags every
`P(A before B)` toward 0.5 on exactly the compounds the product exists to separate.

The error is also **unbiased** (mean +0.003), so it is pure variance with no systematic shift
to correct — which is the good case.

⚠ This is not a free win to bank yet. The comparison is against *literature* pKa, and #7's
`w_w`→`s_s` problem is untouched: both provider and literature are water-scale, while the
mobile phase is not. The 0.170 replaces the *literature-spread* assumption, not the
scale-correction term, which remains a named per-run error.

## 2. The refusal tier is now derived, not hand-labelled — and getting there found a trap

`CONTEXT.md` specifies refusal as "a **structural** test answered from the provider's
microspecies output". The slice hand-labels one compound `refused`. Derived here: **16/16
agreement** with the hand labels.

But the route there is the finding. **Uni-pKa raises `EnumerationError` — "failed to enumerate
microstates across 2 charge states" — for two completely different reasons:**

| compound | why it throws | correct tier |
|---|---|---|
| Benzamide | no ionisable site in range | **in envelope**, `f_neutral` = 1 |
| Benzyltrimethylammonium | quaternary ammonium: only one charge state exists | **refuse** |

A first implementation treated the exception as benign and returned a neutral species with
`f_neutral` = 1. That **silently converted a permanently charged compound into a neutral
one** — precisely the error #33/#34/#39 built the refuse tier to prevent, and it passed
without complaint until the hand labels were compared.

**Formal charge separates them**, and it is a structural property of the molecule rather than
a provider artefact, which is exactly what the refusal test is specified to be:

- `EnumerationError` + formal charge 0 → `no_ionisable_site` (benign)
- `EnumerationError` + formal charge ≠ 0 → `permanently_charged` → **refuse**
- any other exception → `declined` (#33's genuine `PROVIDER_DECLINED`, the ~1-in-15 case)

**A provider exception must never be mapped to a single tier.** Two of the three meanings are
benign, one is a refusal, and they are indistinguishable from the exception text alone.

## 3. Real microspecies vs the closed-form sigmoid

Over 27 (compound, pH) points: mean difference **−0.0095**, max **0.0931**.

Small, and biased slightly low — the provider's pKa values sit marginally below the cited
ones, so it puts a little less neutral population at mid pH. Not alarming, but it multiplies
straight through to `k`, and the closed form cannot represent multiprotic or zwitterionic
behaviour at all, which the microspecies distribution does by construction.

## What is not done

- **Not yet wired into `slice.py`'s sampler.** This measures and validates the provider
  against what the slice assumes; swapping it in is the next step and changes the ensemble.
- **No per-prediction σ.** The `unipka` wrapper exposes no ensemble spread, consistent with
  #7's central negative finding that no provider — open or commercial — emits a calibrated
  per-prediction uncertainty. The 0.170 above is a *global* SD measured on 10 values, not a
  per-compound confidence.
- **10 pKa values is a thin calibration set**, and they are the compounds the slice was built
  around. Treat 0.170 as an estimate with wide error bars, not a constant.

---

# Wiring it into the sampler (`sampler.py`)

`env/ionisation-venv/bin/python model/sampler.py`. Three changes from the slice, all in the
ionisation layer: **no hardcoded pKa**, **no `type` classification** (microspecies handle
mono/multi/zwitterionic by construction), and **`f_neutral` stored per scenario** — #44
punch-list item 2, and a MUST in `spec/scenario-ensemble.md`.

**How uncertainty enters without inventing a σ.** The provider emits a point distribution and
no per-prediction uncertainty (#7's central negative finding). But a pKa error of δ is exactly
a shift of the speciation curve along the pH axis, so the sampler queries
`f_neutral(pH + δ_s)` with `δ_s ~ N(scale_shift, 0.170)`. No fabricated σ, and the measured
provider error goes in directly.

## Result: my prediction was half right, and the half that failed is the informative one

Confident pairs (P > 0.9) out of 105, **mean ± SD over 5 seeds**:

| operating point | slice, hardcoded, sd 0.35 | provider, sd 0.35 | provider, **measured sd 0.170** |
|---|---|---|---|
| (0.30, pH 3.0) | **77.0 ± 0.6** | 75.6 ± 0.5 | 75.8 ± 0.7 |
| (0.30, pH 7.0) | **73.0 ± 0.9** | 71.4 ± 0.8 | 71.4 ± 0.8 |
| (0.50, pH 5.0) | 71.4 ± 1.0 | 71.4 ± 1.0 | **74.6 ± 0.5** |

**I predicted the 0.170-vs-0.35 gap would tighten the ensemble. It does — at exactly one
operating point out of three.**

### Switching to the provider *costs* confidence, and that is correct

Holding the SD fixed at 0.35, the provider loses ~1.5 pairs at two operating points. That is
not a regression: the provider models ionisation the hardcoded version ignored. Methylparaben
is the clearest case — the slice treats its phenol (pKa ≈ 8.4) as neutral throughout, while
the provider puts it at **0.962 ± 0.015** neutral at pH 7, i.e. ~4% ionised. Several compounds
gain real ionisation behaviour they did not have, so variance rises.

**The slice was slightly over-confident because it was ignoring chemistry.** Buying fidelity
by giving up ~1.5 of 105 confident pairs is the right trade.

### The measured SD helps only in the partial-ionisation band

Halving 0.35 → 0.170 gains **+3.2 pairs at (0.50, pH 5.0)** (≈ 6 SD, real) and **nothing at
the other two** (+0.2 and 0.0, both inside noise).

The reason is visible in the stored `f_neutral`:

| compound | (0.30, pH 3.0) | (0.30, pH 7.0) | (0.50, pH 5.0) |
|---|---|---|---|
| Benzoic acid | 0.817 ± 0.118 | 0.001 ± 0.001 | 0.067 ± 0.054 |
| Ibuprofen | 0.892 ± 0.080 | 0.001 ± 0.001 | 0.118 ± 0.082 |
| 4-Aminobenzoic acid | 0.836 ± 0.064 | 0.003 ± 0.003 | 0.220 ± 0.128 |
| Lidocaine | 0.000 ± 0.000 | 0.084 ± 0.048 | 0.001 ± 0.001 |

At pH 7 the acids are **fully ionised** (`f_neutral` ≈ 0.001) — pKa uncertainty has no
leverage, so halving it changes nothing. At pH 3 they are mostly neutral, same story. Only at
pH 5 do they sit mid-curve, where the derivative of `f_neutral` with respect to pKa is
largest, and there the improvement shows up.

**So the pKa layer only matters near pKa.** Obvious in hindsight; worth having measured,
because it means effort spent on pKa accuracy pays off *only* for compounds being run near
their own pKa — and that is precisely the regime #30's sweep concentrates its pH points in.

## What this does not fix

- **`D` still dominates the ionised regime** (0.92 of variance) and is still a prior. Only
  measurement moves it.
- **The `w_w`→`s_s` correction is untouched.** Provider and literature are both water-scale
  while the mobile phase is not, so the scale shift is still carried as a named per-run term
  with the slice's assumed priors. The 0.170 replaced the *literature-spread* assumption only.
- **`sampler.py` does not yet emit a `ScenarioEnsemble`** from `spec/`. It stores `f_neutral`
  per scenario, which was the blocking gap, but the object assembly and the method-card
  plumbing are still to do.

---

# What is measuring `D` worth, and which half? (`voi_D.py`)

`D` owns 0.78–0.92 of the variance for ionised compounds and has never been fitted. But it is
**two terms doing different damage**:

```
D = 1.5 + D_run + D_c        D_run ~ N(0, 0.30)   per-run   (column, modifier)
                             D_c   ~ N(0, 0.40)   per-compound
```

`CONTEXT.md` holds that per-run terms "largely cancel in order while still setting absolute
retention", and per-compound terms "scramble elution order". Elution order is the primary
metric. So this collapses each term's uncertainty independently and watches what is recovered.

Confident pairs of 105, mean ± SD over 5 seeds:

| case | (0.30, pH 3.0) | (0.30, **pH 7.0**) | (0.50, pH 5.0) |
|---|---|---|---|
| prior, nothing measured | 75.8 ± 0.7 | 71.4 ± 0.8 | 74.6 ± 0.5 |
| `D_run` measured → 0.05 | 79.2 (+3.4) | 76.2 (**+4.8**) | 78.0 (+3.4) |
| **`D_c` measured → 0.05** | 80.6 (+4.8) | 84.8 (**+13.4**) | 82.2 (+7.6) |
| both measured | 82.2 (+6.4) | 87.6 (+16.2) | 84.2 (+9.6) |
| `D` known exactly (floor) | 82.0 | 88.0 (+16.6) | 84.2 |

**The per-compound term is worth about three times the per-run term for elution order.** At
pH 7.0, where the acids are fully ionised and `D` dominates, `D_c` alone recovers **13.4 of
the 16.2** available (83%); `D_run` alone recovers 4.8 (30%). The pattern holds at all three
operating points.

That is `CONTEXT.md`'s per-run/per-compound theory confirmed on this model rather than
assumed — per-run terms really do largely cancel in order.

**Absolute retention tells the other half of the story.** At pH 7.0 median `sd(log k)` falls
0.489 → 0.411 with `D_run` measured and → 0.324 with `D_c`, → 0.210 with both. So `D_run` does
buy real precision on *where peaks are*; it just buys much less on *what order they come in*.

## What this means for the run schedule

The sweep's headline output is `D`'s **per-run** part — the paired EVO-minus-Supelcosil
difference that answers #39's silanol question. That remains a legitimate architectural
question. But **for the metric the product is graded on, the per-compound term is worth ~3×
more**, and the same runs measure it.

⚠ **The transferable gain is in the prior width, not the per-compound values.** Measuring
`D_c` for eight specific acids says nothing about a ninth compound — `D_c` is idiosyncratic by
definition. What transfers is the **estimate of `D_c`'s spread**: learning that its SD is
0.20 rather than 0.40 narrows the prior for *every* future compound, and the precision of a
variance estimate scales with the number of compounds.

**So the schedule should favour compound count over condition count**, within whatever
instrument time is available. Concretely: an extra tier-1 acid is worth more than an extra φ
point, and broadening across acid *classes* (already done — arylpropionic, arylacetic,
indole-acetic, aliphatic) is worth more than depth on any one class.

This does not argue for dropping the two-column design: `D_run` is what makes #39's silanol
question answerable at all, and the paired difference is the only way to get it. It argues
about where the *marginal* hour goes.

---

# `width` gets physics: the Knox model (`dispersion.py`)

`N` was a flat `10000 ± 20%` — a number, not a model. That was #44 punch-list item 3.

**`N` is not a column property.** It depends on particle size, flow rate and the solute's
diffusivity, so it belongs to the *method*, and Knox is how the method determines it:

```
h = A·ν^(1/3) + B/ν + C·ν        ν = u·dp/Dm        N = L/(h·dp)
```

Per the map's cut rule this is settled science: the choice is pinned, the interface defined,
the source cited (Knox, *J. Chromatogr. Sci.* **15** (1977) 352), and no textbook is
paraphrased. `A, B, C` ship as engineering priors **by particle class** — and the class
matters, because #30's two arms are not the same class.

## It quantifies a caveat #30 could previously only assert

At the protocol's geometry (150 × 4.6 mm, 5 µm, 1.0 mL/min; ν = 7.71, `t0` = 1.62 min):

| arm | class | `h` | `N` | peak sd at k = 10 |
|---|---|---|---|---|
| **Kinetex EVO C18** | core-shell | 1.53 | **19,553** | 0.1275 min |
| **Supelcosil LC-18** | fully porous | 2.62 | **11,447** | 0.1666 min |

**A factor of 1.71 in plate count**, so peaks on the EVO arm are 0.76× as wide — 24%
narrower. #30 records that "the two arms are not comparable on resolution"; that is now a
number rather than an assertion, and it means **peak capacity differs materially between
arms**. Harmless for `D`, which is a retention ratio, but it must not be read as a
resolution comparison.

Sampled per scenario, the Knox priors give `N` = 19,922 ± 2,804 (14%) and 11,626 ± 1,435
(12%) — inside the ±10–30% #9/#13 budget, which is a useful check that the priors are
neither too tight nor too loose.

## It also found a hole in the spec

`MethodCard` had **no flow rate**. `N` cannot be computed without one, so `flow_ml_min` is now
a required field with no default — the same reasoning as `ph_scale`: a defaulted flow silently
fixes the peak width of every prediction, which is `width` quietly reverting to the
placeholder the spec exists to remove.

Found only by implementing against the spec, which is the argument for implementing against
specs.

## Still owed

- **Gradient compression** — null here because every operating point is isocratic. Needed
  before any gradient method is predicted.
- **Extra-column dispersion as more than a variance add** — currently a volume-to-time
  conversion over the vendor range. Real extra-column behaviour depends on tubing, injection
  volume and detector cell, and is not a single number per instrument.
- **A geometry registry.** `column_id` should resolve to a `ColumnGeometry`; the sampler
  currently hardcodes one.
- **`Dm` is a constant** (1.0e-5 cm²/s). It varies with solute size, temperature and
  composition, and the temperature dependence in particular is not modelled.

---

# Column registry, compound plate, and mixture pooling (`columns.py`, `sequence.py`)

## The registry declares only what cannot be looked up

`spec/scenario-ensemble.md` says `column_id` "resolves to a column vector". `columns.py` is
that resolution. An entry carries **identity, vendor, geometry and pH rating**; the HSM phase
descriptors and LSER availability are **read from `sources/`**. Re-typing parameters into
Python would create a second source of truth that drifts from the first. Adding a phase is
three lines, and `unregistered()` lists what is available but not yet declared.

| `column_id` | `C(7.0)` | `H` | type | LSER | pH | `N` @1 mL |
|---|---|---|---|---|---|---|
| kinetex-evo-c18 | −0.01 | 1.01 | B | **yes** | 1–12 | 19,553 |
| supelcosil-lc-18 | **1.75** | 1.01 | **A** | – | 2–7.5 | 11,447 |
| xbridge-c18 | 0.13 | 1.00 | B | **yes** | 1–12 | 11,447 |
| ascentis-rp-amide | 0.08 | 0.84 | EP | – | 2–9 | 11,447 |
| zorbax-bonus-rp | **−1.10** | 0.65 | EP | – | 2–9 | 11,447 |

## ⚠ Expanding the column set buys C-span and spends the H control

All five span `C(7.0)` −1.10 to +1.75 (2.85) against the booked pair's 1.76. But
**hydrophobicity is not matched across the wider set** — `H` runs 0.65 to 1.01.

#30 matched its pair to `H` = **0.002** deliberately, so the contrast is silanol activity and
*not* retention strength. **A `C`-versus-`D` relationship fitted across columns that differ in
both is confounded, and no amount of instrument time separates them afterwards.**

**`xbridge-c18` is the clean addition**: `H` = 1.00 matches the existing pair to 0.01, it is
**in WSU-2019** so its LSER constants are already held, and it is rated pH 1–12. It adds a
third point at `C(7.0)` = 0.13 between EVO and Supelcosil at no cost to the control.

**`zorbax-bonus-rp` and `ascentis-rp-amide` are different propositions.** Bonus RP has the
lowest `C(7.0)` in the entire 819-column database (−1.10) and would extend the span
substantially — but at `H` = 0.65 against 1.01, it differs in retention strength as much as in
silanol activity. Worth running as a *separate* question about polar-embedded phases; not
worth folding into the `D_run` measurement.

⚠ Three of the five pH ratings are **unverified** — only EVO and Supelcosil have been checked.

## Pooling is sized on the limiting arm

For isocratic peaks with σ = `t_R`/√`N`, `Rs` ≥ 1.5 requires
`(t2−t1)/(t2+t1) ≥ 3/√N` — about **5.6% in retention time**.

| arm | `N` | required separation |
|---|---|---|
| Kinetex EVO C18 | 19,553 | 2.1% |
| **Supelcosil LC-18** | 11,447 | **2.8%** |

A **1.31× stricter** constraint, and the thin-slice set needs **2 pools on EVO but 3 on
Supelcosil** — the arms genuinely disagree. **Use the Supelcosil pooling on both arms**, so
the same mixture is injected on each and the paired difference stays paired.

Retention is predicted with EVO's LSER constants and Supelcosil's plate count. That follows
from the `H` = 0.002 match: `k` should be comparable, and the arms differ in silanol activity
and efficiency rather than retention strength. Supelcosil has no constants of its own, being
outside WSU-2019 — the same asymmetry the tier structure already handles.

## The plate is blocked on five descriptor sets

| acid | pKa | `k` @30% ACN | `t_R` |
|---|---|---|---|
| Ibuprofen | 4.45 | 25.3 | 42.6 min |
| Naproxen | 4.18 | 14.3 | 24.8 min |
| Ketoprofen | 4.20 | 10.7 | 19.0 min |
| Flurbiprofen, Fenoprofen, Diclofenac, Gemfibrozil, Indomethacin | — | **unknown** | — |

The three measured acids fit **one pool** on Supelcosil. The other five were screened on
McGowan `V` (1.84–2.53, all above ibuprofen's 1.78), so they will elute **later** and probably
spread across pools — but "probably" is not a sequence.

**SoluteML on those five SMILES closes it.** That is a GPU job, and `env/install_qspr.sh`
already exists from the merged branches.
