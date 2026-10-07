# Is there a usable HSM → LSER bridge?

Research note resolving issue #8.
Investigated 2026-08-14. All primary-source claims carry a URL or DOI.

---

## Verdict

**No usable bridge. Viable for shortlisting only.**

There is **no published mapping** from HSM column parameters (`H, S*, A, B, C`) to LSER system
constants (`c, e, s, a, b, v`) in either direction, and the one paper that examined the question
directly concluded such a mapping is not obtainable. The overlap set of columns characterised by
both models is roughly **25–30 columns**, which is too small to fit a 5→6 parameter mapping with
characterisable uncertainty. Structurally, HSM is a *selectivity* model at a *single fixed mobile
phase*; LSER as we use it must produce `log k` as a function of `φ`. HSM carries no `φ` axis and
no absolute-retention term, so even a perfect parameter mapping could not emit `k`.

**But the fallback is real and cheap.** The full HSM database is a 819-row CSV, freely
downloadable, and the `Fs` column-comparison function is a five-line arithmetic expression that
reproduces exactly offline. That gives a genuine capability — *shortlist / orthogonality-rank
~800 commercial columns* — with no modelling risk, because it makes no claim about `k`.

The curated LSER column set remains the committed path. Nothing here changes that.

---

## 1. The PQRI / HSM column database

### What it is

The Hydrophobic Subtraction Model characterises a reversed-phase column by five parameters via

```
log(k / k_EB) ≡ log α = η'·H − σ'·S* + β'·A + α'·B + κ'·C
```

where `η', σ', β', α', κ'` are solute parameters and `H, S*, A, B, C` are the column parameters:
hydrophobicity, steric resistance, H-bond acidity, H-bond basicity, and cation-exchange activity.
Equation and parameter definitions verbatim from Zhang & Carr, *J. Chromatogr. A* **1216** (2009)
6685–6694, [doi:10.1016/j.chroma.2009.06.048](https://doi.org/10.1016/j.chroma.2009.06.048)
(open access: <https://pmc.ncbi.nlm.nih.gov/articles/PMC3195507/>); original model:
Snyder, Dolan & Carr, *J. Chromatogr. A* **1060** (2004) 77–116,
[doi:10.1016/S0021-9673(04)01480-3](https://doi.org/10.1016/S0021-9673(04)01480-3).

### Size — measured, not quoted

I downloaded the database directly:

```
https://hplccolumns.org/database/database.csv
```

(the link is present but HTML-commented-out on
<https://hplccolumns.org/database/index.php>; the file itself serves fine, `text/csv`, latin-1
encoded). Retrieved 2026-08-14:

| Fact | Value |
|---|---|
| Rows (columns characterised) | **819** |
| Distinct manufacturer strings | 64 |
| Fields | `id, name, manufacturer, manufacturerID, type, H, S, A, B, C28, C70, retention, USPtype, phase` |
| USP L1 (C18) rows | 402 |
| USP L7 (C8) rows | 115 |
| USP L11 (phenyl) rows | 82 |
| Rows with no USP type | 55 |

Note the `C` parameter appears **twice** — `C28` and `C70`, the cation-exchange term at pH 2.8 and
pH 7.0. `retention` is the ethylbenzene retention factor used to normalise for phase ratio.

Published counts differ from the file and from each other, which is worth knowing:

- **750 columns** — Rutan, Kempen, Dahlseid, Kruger, Pirok, Shackman, Zhou, Wang & Stoll,
  *J. Chromatogr. A* **1731** (2024) 465127,
  [doi:10.1016/j.chroma.2024.465127](https://doi.org/10.1016/j.chroma.2024.465127):
  "the column database for this model currently stands at 750 columns".
- **"more than 600 … from over 30 manufacturers"** — the `<title>` of
  <https://hplccolumns.org/about/index.php> (stale site copy).
- **368 columns** — the USP mirror, <https://www.usp.org/resources/pqri-approach-column-equiv-tool>,
  tool at <https://apps.usp.org/app/USPNF/columnsDB.html>.

Take **819** as the current figure; the USP mirror is a smaller, lagging subset. A chromatography
forum thread also reports the USP database as the smaller and less-validated of the two
(<https://www.chromforum.org/viewtopic.php?t=17617>) — secondary, flagged as such.

### Licence — this is a live constraint

> "Except where otherwise noted, all content on this site is licensed under a Creative Commons
> Attribution-Noncommercial-Share Alike 3.0 United States License"

— <https://hplccolumns.org/about/index.php>

**CC BY-NC-SA 3.0 US.** Non-commercial, and share-alike. Free for a spec, a prototype, and
academic work. It is *not* redistributable inside a commercial product without a separate
agreement with the maintainer. The USP page carries no explicit licence statement at all, which is
worse, not better. **Anything that ships commercially must resolve this with Dwight Stoll and/or
USP.** I could not determine whether a commercial licence path exists.

### Currency and maintenance — actively maintained

> "In late 2011, the work of characterizing of new stationary phases transitioned from the BAS
> laboratory in Oregon to the laboratory of Dr. Dwight Stoll at Gustavus College in Minnesota.
> Dr. Stoll currently maintains the HS database and is pursuing further development of the model
> with a number of collaborators. Column manufacturers interested in having their stationary
> phases listed in the HSM database should contact Dr. Stoll for more information."

— <https://hplccolumns.org/about/index.php>

The model itself is under active revision: HSM2 and then **HSM3** (Rutan et al. 2024, DOI above),
trained on 1014 solute/stationary-phase combinations, predicting `ln α` to average absolute error
≈ 0.033 (≈3 % in α). That paper also **released 43,329 retention measurements at multiple mobile
phase compositions** for community use. Preprint:
<https://chemrxiv.org/engage/chemrxiv/article-details/660b4d3666c1381729702ceb>.

**Could not determine:** the last-modified date of `database.csv`; whether the published HSM3
column parameters have been merged into the public CSV (the CSV schema has no version field); and
the exact hosting location and licence of the 43,329-measurement dataset (the ChemRxiv and
ScienceDirect pages both 403'd from here).

### Measurement conditions — the single most important limitation

HSM column parameters are measured at **one** condition: 50 % acetonitrile / phosphate buffer
pH 2.8 (30 mM), 35 °C, ethylbenzene as reference solute, 16 probe solutes
(Zhang & Carr 2009, <https://pmc.ncbi.nlm.nih.gov/articles/PMC3195507/>). Only `C` is reported at
a second pH (7.0). There is **no `φ` dependence in the HSM parameter set at all.**

### Data quality caveat

Shackman, *Discrepancies in column parameters presented in hydrophobic subtraction model
manuscripts*, *J. Chromatogr. A* (2016),
[doi:10.1016/j.chroma.2016.11.006](https://doi.org/10.1016/j.chroma.2016.11.006) — the published
HSM parameter tables are not internally consistent across sources. I could not read the full text
(paywalled; Semantic Scholar reports the abstract elided by the publisher), so I cannot say how
large the discrepancies are, only that a peer-reviewed paper exists asserting they are real. **If
the HSM CSV is used for anything, it needs a provenance check against this paper.**

---

## 2. Is there published work relating HSM to LSER?

**Searched hard. The answer is no — and there is an explicit published negative.**

The decisive source is Poole, *Reversed-phase liquid chromatography system constant database over
an extended mobile phase composition range for 25 siloxane-bonded silica-based columns*,
*J. Chromatogr. A* **1600** (2019) 112–126,
[doi:10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027)
(PMID 31128882). Final sentence of the abstract, verbatim:

> "It is shown that the interaction parameters employed in the hydrophobic-subtraction model have
> little overlap with those of the solvation parameter model and a quantitative comparison of
> column properties delineated by both models is not possible."

This is the strongest possible form of the answer for our purposes: the one worker who holds a
large LSER column database *and* had the HSM parameters in hand, looked, and reported that the
comparison cannot be made quantitatively.

Same paper, also from the abstract:

> "Contributions to retention from steric resistance and cation-exchange interactions **not
> parameterized in the solvation parameter model** are identified for columns and conditions where
> they may be important."

That is Poole naming the exact gap this ticket asked about, from the LSER side.

### What I searched and did not find

- No paper proposing regression of `c, e, s, a, b, v` on `H, S*, A, B, C` or the reverse.
- No paper predicting LSER system constants for an unmeasured column from any column-level
  descriptor set.
- Related-but-different work exists and should not be mistaken for a bridge:
  - Kaliszan-style QSRR vs. HSM comparison — *Comparative characteristics of HPLC columns based on
    quantitative structure–retention relationships (QSRR) and hydrophobic-subtraction model*,
    *J. Chromatogr. A* (2005), PMID 15974124,
    <https://pubmed.ncbi.nlm.nih.gov/15974124/>. Compares *classifications*, not parameters.
  - Poole, *Column selectivity from the perspective of the solvation parameter model*,
    *J. Chromatogr. A* (2002),
    <https://www.sciencedirect.com/science/article/abs/pii/S0021967301013619>. LSER-only.
  - Marchetto, Tirapelle, Mazzei, Sorensen & Besenhard, *In Silico HPLC Method Development via
    Machine Learning*, *Anal. Chem.* **97** (2025) 6991–7001,
    [doi:10.1021/acs.analchem.4c03466](https://doi.org/10.1021/acs.analchem.4c03466) (CC-BY,
    already in `sources/README.md`) — SMILES → QSPR → LSER → LSS. Uses LSER; does not touch HSM.

**Could not determine:** whether an unpublished or industrial mapping exists (ACD/Labs ships a
column selector, <https://www.acdlabs.com/resources/free-chemistry-software-apps/column-selector/>,
but its internals are not documented publicly). Absence of evidence in the open literature is what
I can report; it is fairly strong evidence here, because the question is obvious and the people
holding both datasets have written about it.

---

## 3. Columns characterised by both models — the overlap set

This is the number that decides whether a mapping is even estimable, so I measured it rather than
guessing.

### How many columns have published LSER system constants?

The Wayne State University (Poole) system constant database is the largest coherent public source.

- **2019**: 25 siloxane-bonded type-B silica columns, 10–70 % (v/v) methanol and acetonitrile,
  plus THF for two columns. [doi:10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027)
- **2025 update**: Poole & Atapattu, *Update of the Wayne State University system constant database
  for reversed-phase liquid chromatography columns for varied mobile phase compositions*,
  *J. Chromatogr. A* **1762** (2025) 466385,
  [doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385). System
  constants at 10 % increments over 10–70 % for **methanol-water (27 columns), acetonitrile-water
  (25), acetone-water (7), THF-water (3), 2-propanol-water (1)**. These are overlapping column
  sets per solvent, not additive — the distinct-column count is on the order of **~30**.

So the LSER column axis, as published, is **~25–30 columns**. Which — worth noting for the map —
is only about double the committed curated set, not an order of magnitude more.

### Are those columns in the HSM database?

I extracted the column names Poole names explicitly in abstracts I could read
(PMID 33220586, [doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692);
PMID 33161359, [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652))
and matched them against the 819-row HSM CSV. **Every one is present:**

| LSER-characterised column (Poole) | Present in HSM CSV |
|---|---|
| Kinetex C18 | `Kinetex C18 100A` |
| Kinetex EVO C18 | `Kinetex EVO C18` |
| Kinetex XB-C18 | `Kinetex XB-C18` |
| Kinetex C8 | `Kinetex C8` |
| Kinetex Biphenyl | `Kinetex Biphenyl 100A` |
| Kinetex Phenyl-Hexyl | `Kinetex Phenyl-Hexyl` |
| Kinetex F5 | `Kinetex F5` (also `Kinetex PFP 100A`) |
| SunFire C18 | `Sunfire C18` |
| XBridge Shield RP18 | `XBridge Shield RP18` |
| XBridge C8 | `XBridge C8` |
| XBridge Phenyl | `XBridge Phenyl` |
| Discovery HS F5 | `Discovery HS F5` |
| Synergi Hydro-RP | `Synergi Hydro-RP` |

**So the overlap set is essentially the whole LSER set: ~25–30 columns.** The HSM database is not
the binding constraint; the LSER side is.

### Why ~30 is not enough

A linear mapping `(H, S*, A, B, C28, C70) → (c, e, s, a, b, v)` is 6 predictors + intercept per
target × 6 targets = **42 free parameters fit from ~30 observations**. Under-determined before we
even ask for uncertainty. And the ticket's own standard — *"if a bridge exists, its own uncertainty
must be characterisable, or it is not usable"* — cannot be met: with n≈30 and p=7 per target there
is no honest held-out set, and the HSM parameters themselves carry uncertainty that Rutan et al.
(2024) describe as large enough that "for some geometric isomer pairs the interactions mainly
responsible for the observed selectivities could not be rationalized due to large uncertainties for
particular terms in the model".

**Direct answer to "could anyone even test a mapping?"** — Not credibly, no. Independent of
whether one has been published, the data to fit or falsify one does not exist at usable scale. ~30
overlapping columns against 42 free parameters is under-determined; and because the LSER side is
the binding constraint, this will not improve by anyone mining the HSM database harder. It would
take a new experimental campaign on the LSER side.

**I could not determine** the exact 25-column list from the 2019 paper (paywalled; ScienceDirect
403s). The 13 names above are a verified lower bound on the overlap, not the complete set.

---

## 4. Where the two models can and cannot correspond

### The term-by-term picture

| HSM term | What it is | LSER analogue |
|---|---|---|
| `H` — hydrophobicity | phase hydrophobic strength | **partial** — spread across `v` (cavity/volume) and the intercept `c`; not separable, because HSM's `H` is defined *relative to ethylbenzene* while LSER's `v` and `c` are absolute |
| `S*` — steric resistance | resistance of the bonded phase to penetration by bulky molecules | **none.** Poole 2019 says so explicitly. Shape/steric selectivity is not parameterised in the solvation parameter model at all |
| `A` — H-bond acidity of phase | phase donates H-bond to basic solutes | `b` (system term acting on solute basicity `B`) — **conceptually aligned, oppositely signed convention, different reference state** |
| `B` — H-bond basicity of phase | phase accepts H-bond from acidic solutes | `a` (system term acting on solute acidity `A`) — same caveat |
| `C` — cation exchange | ionised-silanol / ionised-base ion exchange, reported at pH 2.8 and 7.0 | **none.** LSER has no electrostatic/ion-exchange term; it is a neutral-solute model |
| — | — | `e` (excess molar refraction / polarisability) — **no HSM counterpart** |
| — | — | `s` (dipolarity/polarisability) — **no HSM counterpart** |

### The three structural gaps, stated plainly

1. **HSM terms with no LSER analogue: `S*` and `C`.** Steric resistance and cation exchange. These
   are precisely the two mechanisms that dominate for bulky and for basic/ionisable solutes — i.e.
   exactly the compounds this effort has committed to handling (`pH` and `pKa` are first-class
   inputs per issue #1). A bridge would have to throw away the two terms we most need.

2. **LSER terms with no HSM analogue: `e` and `s`.** Polarisability and dipolarity. HSM folds these
   into `H` and `S*` implicitly. So the mapping is not just lossy in one direction — it is lossy in
   both, which means it is not invertible, which means it cannot be regularised by round-tripping.

3. **The reference-state mismatch, which is fatal on its own.** HSM predicts
   `log α = log(k/k_EB)` — retention *relative to ethylbenzene on the same column*, with the phase
   ratio deliberately divided out (that is the point of the `retention` field in the CSV). LSER
   predicts `log k` absolutely. **HSM structurally cannot produce `k`.** Recovering `k` needs
   `k_EB` for that column under *our* conditions, which is a per-column, per-condition measurement
   the HSM database does not contain.

4. **No `φ` axis.** HSM parameters exist at 50 % ACN / pH 2.8 / 35 °C and nowhere else. LSER system
   constants vary strongly and non-linearly with `φ` — that variation *is* the content of Poole's
   system maps, and the 2025 update reports it at 10 % increments across 10–70 %. Even granting a
   perfect parameter mapping at 50 % ACN, we would have a single point on a curve we need in full.
   The destination architecture needs `k(φ, pH)`. HSM gives one condition, one relative quantity.

**Conclusion of the theory section:** the two models are not two coordinate systems for the same
space. They are two different decompositions with partially disjoint spans, different reference
states, and different independent variables. Poole's "not possible" is not pessimism; it follows
from the definitions.

### Is the mismatch structural, or merely unmeasured?

Worth separating, because "unmeasured" would be fixable by more data and "structural" is not.

- Gaps 3 and 4 (**reference state, and no `φ` axis**) are **structural**. They are consequences of
  what each model is *defined* to predict — `log α` relative to ethylbenzene at one fixed condition
  versus `log k` as a function of `φ`. No quantity of additional column characterisation changes
  them. A bridge cannot be repaired here; it can only be bypassed by measuring retention directly.
- Gaps 1 and 2 (**`S*`/`C` vs `e`/`s`**) are **structural in the span, contingent in degree.** The
  solvation parameter model has no steric or electrostatic term to receive `S*` and `C` — that is a
  fact about the model's functional form, not about data coverage. What *is* contingent is how much
  it costs: if `S*` and `C` happened to be near-constant across columns, or strongly collinear with
  `H`, the loss would be small in practice. They are not. Poole 2019 specifically identifies
  columns and conditions where steric resistance and cation-exchange contributions "may be
  important", and `C` is reported at two pH values precisely because it swings.
- The **overlap-set size** (§3) is the only genuinely *unmeasured* limitation. In principle someone
  could characterise 200 columns by both models and revisit. Nobody has, and it would be a
  multi-year experimental programme, not a modelling exercise.

So: even with unlimited data, the reference-state and `φ` gaps alone are fatal to a `k`-predicting
bridge. **This is not a "not yet" — it is a "not this way".**

---

## 5. Fallback: HSM for column shortlisting — this one works

The weaker use — *rank/shortlist columns by selectivity similarity without predicting `k`* — is
sound, standard practice, and available to us immediately.

### The Fs function

Verbatim from Zhang & Carr 2009 (<https://pmc.ncbi.nlm.nih.gov/articles/PMC3195507/>):

```
Fs = { 12.5(H₂−H₁)² + 100(S*₂−S*₁)² + 30(A₂−A₁)² + 143(B₂−B₁)² + 83(C₂−C₁)² }^(1/2)
```

> "In the extreme case when two phases are very close (Fs ≤ 3), the two can be considered to be
> chromatographically 'equivalent' in terms of phase selectivity."

USP states the same threshold: "Columns which have values of F ≤ 3 are very likely to give an
equivalent and acceptable separation" (<https://www.usp.org/resources/pqri-approach-column-equiv-tool>).

### I verified it reproduces offline

Computed `Fs` for all 819 CSV rows against `XBridge C18` (H=1.00, S*=0.02, A=−0.09, B=0.00,
C28=0.17), using `C28` for the `C` term:

> ### ⚠ CORRECTED 2026-08-16 — the values below were computed with the wrong formula
>
> This section originally printed `Fs` values of 0.14–0.28. Those came from
> `Fs = sqrt(12.5·ΔH² + 100·ΔS*² + …)` — the weights applied to the **squared** differences —
> where the published definition squares the **weighted** differences:
> `Fs = sqrt((12.5·ΔH)² + (100·ΔS*)² + …)`. Because the parameter differences are small,
> squaring before weighting shrinks every term, so the values came out ~7× too low.
> Verified in `experiments/fs-calibration/`, which reproduces the old numbers exactly with the
> slipped formula. The corrected values are below.

```
0.00  XBridge C18            (self — sanity check)
0.93  Discovery C18
1.31  Hypersil Beta Basic-18
1.68  Hypurity C18
1.96  Poroshell 120 EC-C18
1.97  Flowrosil ODS
2.14  ChromCore 120 C18
2.20  Venusil XBP C18(L)
```

Five lines of arithmetic on a free CSV. **No modelling risk, no fitted parameters of ours, no
uncertainty budget to defend** — because it makes no claim about `k`.

### It is validated in the literature for exactly this job

- Græsbøll, Nielsen & Christensen, *Using the hydrophobic subtraction model to choose orthogonal
  columns for online comprehensive two-dimensional liquid chromatography*, *J. Chromatogr. A*
  **1326** (2014) 39–46,
  [doi:10.1016/j.chroma.2013.12.034](https://doi.org/10.1016/j.chroma.2013.12.034) — used HSM to
  *predict* orthogonal column pairs, then confirmed by LC×LC-DAD experiment. Also found that
  **sample-specific F-weights beat the standard Gilroy weights**, and that excluding the `C` term
  (as Dolan & Snyder suggested) performed worse than sample-tuned weights for their acid/base PAC
  set.
- *Column selection for comprehensive two-dimensional liquid chromatography using the hydrophobic
  subtraction model*, *J. Chromatogr. A* (2019),
  [doi:10.1016/j.chroma.2018.09.018](https://doi.org/10.1016/j.chroma.2018.09.018).

### How well does Fs actually predict selectivity agreement?

Honestly: **directionally well, quantitatively unquantified.**

- The strongest evidence is Græsbøll et al. (2014), because it closed the loop experimentally —
  columns predicted orthogonal from HSM parameters *were* orthogonal when run on LC×LC-DAD. But
  the paper's headline is a caveat, not an endorsement: the standard Gilroy F-weights were
  **beaten** by weights re-derived from the actual sample analytes, and the authors conclude "this
  emphasizes the necessity of considering the nature of the sample when choosing orthogonal
  columns."
- The underlying model's accuracy is known: HSM3 predicts `ln α` to average absolute error ≈ 0.033
  (≈3 % in α) (Rutan et al. 2024). That is the ceiling on any Fs-derived claim, and it is the
  *average* — the same paper reports specific isomer pairs where term uncertainties were too large
  to interpret.
- What I could **not** find anywhere: a published calibration of `Fs` against observed selectivity
  agreement — i.e. "columns at `Fs = 2` agree on elution order X % of the time." The `Fs ≤ 3`
  threshold is asserted in the primary literature and by USP as a rule of thumb ("very likely to
  give an equivalent and acceptable separation"), never as a calibrated probability.

For an effort where **uncertainty is first-class**, that matters: Fs as shipped is an
*uncalibrated* similarity score. If we surface it, we surface a rank, not a confidence. Calibrating
Fs against real elution-order agreement would be a genuine contribution — and note we would have
the raw material to try, via the 43,329-measurement HSM3 dataset (§7.2).

### Known limitations of Fs — carry these forward

From Zhang & Carr 2009, who propose an alternative *because* of them:

1. **Differences, not ratios.** "the 'effective selectivity' … was mathematically proven to be
   given by the *ratios* of the system (phase) dependent interaction coefficients and not by their
   *absolute* values." Fs uses differences.
2. **Composite, therefore sample-blind.** "Fs is a composite parameter as it depends on all the
   column characteristics and therefore it will indicate that two phases will behave differently
   even when the solute set of interest does not contain solutes which differ in the selectivity
   characteristic chiefly responsible for the differences in two phases." (Græsbøll et al.'s
   sample-specific weights are the practical fix.)
3. **Reference-bound.** "A single table of Fs values relative to a fixed reference column does not
   allow two or more test columns to be compared to one another."
4. **The `C` term dominates**, causing tight clustering that obscures differentiation.

~~And one I measured: against `XBridge C18`, **548 of 819 columns** score `Fs ≤ 3`. The ≤3
"equivalent" threshold is not discriminating for a mainstream C18 — the commercial C18 market is
genuinely crowded.~~

> **⚠ CORRECTED 2026-08-16.** That count used the slipped formula described above. With the
> published definition it is **18 of 819**, not 548 — so the conclusion **inverts**: the ≤3
> threshold *is* highly discriminating, and `Fs ≤ 3` returns a genuine shortlist rather than
> two-thirds of the database.
>
> #36's decomposition is **unaffected**: its H figure (0.4%) matches the correct formula
> (0.3–0.5%) and not the slip (1.9%), so #36 used the published definition. `C` dominates
> `Fs²` at 80–86% however it is computed.

`Fs` is still more useful as a **ranking** (nearest / most-orthogonal) than as a binary
equivalence test — see `experiments/fs-calibration/`, which measures what a given `Fs` actually
buys: `Fs ≤ 5` costs a median 0.109 `log k` in substitution error, and `Fs`-nearest is the
LSER-nearest column only 1–3 times in 20.

### What this buys the effort, honestly stated

It buys **"which other columns are worth trying"**, not **"what will elute when on them"**. That is
a real product feature — it is roughly what the USP tool does — but it lives in a different layer
from the retention model and must never be presented as retention prediction. It also composes
cleanly with the committed path: predict `k` properly on the curated LSER columns, and use Fs to
tell the user which of ~800 commercial columns are near-neighbours of a curated one they already
have a prediction for. **That near-neighbour claim is qualitative and must stay qualitative** — an
`Fs ≤ 3` neighbour is *not* a licence to transfer the LSER prediction onto it numerically.

---

## 6. What I could not determine

- The exact 25-column list of the 2019 WSU LSER database (paywalled). The 13-column overlap in
  §3 is a verified lower bound, not the full intersection.
- Whether a commercial licence for the HSM database is obtainable, and on what terms.
- The last-modified date of `database.csv`, and whether HSM3 parameters have been merged into it.
- The magnitude of the parameter discrepancies reported by Shackman (2016) — abstract elided by
  the publisher.
- Whether the per-column standard errors on `H, S*, A, B, C` are published anywhere; they are not
  in the CSV. Without them, even the Fs ranking has no error bars.
- Hosting location and licence of the 43,329-measurement HSM3 retention dataset.

---

## 7. New questions this raises

1. **The LSER column axis is ~30 columns, not 5–15.** Poole's 2025 update characterises ~30
   columns at 10 % `φ` increments across five organic modifiers. That is a bigger committed axis
   than issue #1 assumes, and it is `φ`-resolved, which is exactly what LSS/`k(φ)` needs. Worth a
   ticket: *what is actually in the WSU-2025 database and can we obtain it?*
2. **The 43,329-measurement HSM3 dataset may be worth more than the bridge was.** Retention
   measurements at multiple mobile phase compositions on 13 well-defined columns is raw material we
   could fit *our own* LSER constants to — going around the bridge rather than through it. This is
   a much better-posed problem than mapping parameter-to-parameter. Recommend a ticket.
3. **`S*` and `C` have no LSER analogue — so what covers shape selectivity and ion exchange?**
   Issue #1 already lists "secondary retention mechanisms" under *Not yet specified*, contingent on
   this result. The result is in: **the gap is real and HSM will not close it.** That fog needs its
   own decision, independent of HSM.
4. **HSM parameter provenance needs a check** before the CSV is used even for Fs (Shackman 2016).
5. **Licensing (CC BY-NC-SA) will bite if this ever ships.** Not a spec problem; is a
   destination-adjacent problem worth recording now.

---

## Sources

| # | Source | Link |
|---|---|---|
| 1 | Snyder, Dolan, Carr — *The hydrophobic-subtraction model of reversed-phase column selectivity*, J. Chromatogr. A 1060 (2004) 77–116 | [doi:10.1016/S0021-9673(04)01480-3](https://doi.org/10.1016/S0021-9673(04)01480-3) |
| 2 | Zhang & Carr — *A visual approach to stationary phase selectivity classification based on the Snyder–Dolan HSM*, J. Chromatogr. A 1216 (2009) 6685–6694 (open access) | [doi:10.1016/j.chroma.2009.06.048](https://doi.org/10.1016/j.chroma.2009.06.048) · [PMC3195507](https://pmc.ncbi.nlm.nih.gov/articles/PMC3195507/) |
| 3 | **Poole — RPLC system constant database … 25 siloxane-bonded columns, J. Chromatogr. A 1600 (2019) 112–126. The explicit negative result.** | [doi:10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027) |
| 4 | Poole & Atapattu — *Update of the WSU system constant database…*, J. Chromatogr. A 1762 (2025) 466385 | [doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385) |
| 5 | Rutan et al. (Stoll) — *Improved HSM … large dataset … isomer selectivity* (HSM3), J. Chromatogr. A 1731 (2024) 465127 | [doi:10.1016/j.chroma.2024.465127](https://doi.org/10.1016/j.chroma.2024.465127) · [preprint](https://chemrxiv.org/engage/chemrxiv/article-details/660b4d3666c1381729702ceb) |
| 6 | Shackman — *Discrepancies in column parameters presented in HSM manuscripts*, J. Chromatogr. A (2016) | [doi:10.1016/j.chroma.2016.11.006](https://doi.org/10.1016/j.chroma.2016.11.006) |
| 7 | Græsbøll, Nielsen, Christensen — *Using the HSM to choose orthogonal columns for LC×LC*, J. Chromatogr. A 1326 (2014) 39–46 | [doi:10.1016/j.chroma.2013.12.034](https://doi.org/10.1016/j.chroma.2013.12.034) |
| 8 | *Column selection for comprehensive 2D-LC using the HSM*, J. Chromatogr. A (2019) | [doi:10.1016/j.chroma.2018.09.018](https://doi.org/10.1016/j.chroma.2018.09.018) |
| 9 | Poole & Atapattu — *Selectivity evaluation of core-shell silica columns … solvation parameter model*, J. Chromatogr. A 1634 (2020) 461692 | [doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692) |
| 10 | Poole — *Selection of calibration compounds for selectivity evaluation …*, J. Chromatogr. A 1633 (2020) 461652 | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652) |
| 11 | HSM database, licence and maintenance statement | <https://hplccolumns.org/about/index.php> |
| 12 | HSM database CSV (819 rows, retrieved 2026-08-14) | <https://hplccolumns.org/database/database.csv> |
| 13 | USP — PQRI approach for selecting columns of equivalent selectivity (368 columns) | <https://www.usp.org/resources/pqri-approach-column-equiv-tool> · [tool](https://apps.usp.org/app/USPNF/columnsDB.html) |
| 14 | Marchetto et al. — *In Silico HPLC Method Development via Machine Learning*, Anal. Chem. 97 (2025) 6991–7001 (CC-BY) | [doi:10.1021/acs.analchem.4c03466](https://doi.org/10.1021/acs.analchem.4c03466) |
| 15 | Kaliszan et al. — *Comparative characteristics of HPLC columns based on QSRR and HSM*, J. Chromatogr. A (2005) | [PMID 15974124](https://pubmed.ncbi.nlm.nih.gov/15974124/) |
| 16 | Poole — *Column selectivity from the perspective of the solvation parameter model*, J. Chromatogr. A (2002) | <https://www.sciencedirect.com/science/article/abs/pii/S0021967301013619> |
