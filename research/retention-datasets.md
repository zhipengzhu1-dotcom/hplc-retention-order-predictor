# Retention datasets: what varies method conditions?

Resolves #4. Part of the map, #1.

**Question.** Which public retention datasets can validate a `k(φ, pH, column, T)` model, and what do they
actually cover — measured against our primary metrics (pairwise elution-order accuracy, calibration), not
against `log k` RMSE?

**Method.** Claims about RepoRT are computed from the repository itself, not from its papers. RepoRT was
cloned at commit `9de8d603377bbb6cc0f74e250eb7533fd8874df1` (2026-07-09) and all 421 dataset metadata
files were aggregated directly. Numbers below labelled *(measured)* come from that snapshot; numbers
labelled *(published)* come from the papers and may differ because the repository has grown.

---

## Headline answer

**Yes for temperature and modifier. No for pH crossed with either.**

Public data contains a genuine, designed factorial that varies **temperature (30/40/50 °C)** and
**organic modifier (ACN/MeOH)** across **six columns**, over a shared core of **~500 compounds** —
the "BGC" datasets in RepoRT. This is far better coverage than the ticket assumed.

Public data *also* contains pH sweeps (pH 2/3/5/8) over a fixed compound set — but on a **different
compound set, at a single temperature**. The two do not cross.

**The `pH × T` interaction term is unfalsifiable with public data.** See
[The critical question](#the-critical-question-does-anything-cross-ph-and-temperature) for the exact
measurement.

> ### ⚠ But coverage is the wrong question. Discriminative power is the question.
>
> After the #2 finding that a naive "predict the order once, reuse everywhere" strategy scores 92.3%
> on the blueprint paper's own φ data, I measured **actual elution-order reversal rates** across every
> condition contrast in RepoRT. The result inverts the ranking above:
>
> | Condition varied | Order reversal | Naive baseline | Verdict |
> |---|---|---|---|
> | **pH 3 → 8** (ionisable subset) | **25–31%** | 69–75% | **the only strongly discriminating axis** |
> | **pH 3 → 8** (all compounds) | **23.5%** | 76.5% | strongly discriminating |
> | Modifier ACN → MeOH | 5.6–8.2% | ~93% | usefully discriminating, and large-n |
> | Column (biphenyl vs C18) | 4.6–6.2% | ~95% | weakly discriminating |
> | Column (C18 vs C18) | 1.0–3.9% | ~97% | barely discriminating |
> | Gradient programme (φ path) | 0–3.3% | ~98% | barely discriminating |
> | **Temperature 30 → 50 °C** | **0.9–1.6%** | **~99%** | **effectively non-discriminating** |
> | METLIN SMRT | n/a — one condition | not computable | cannot support the metric at all |
>
> **The axis with the best data coverage (temperature, ~500 compounds × 6 columns) has the least
> discriminative power (~1%). The axis with the worst coverage (pH, 71 compounds × 2 columns) has by
> far the most (23–31%).** Full measurement in
> [Discriminative power for order prediction](#discriminative-power-for-order-prediction).

---

## 1. RepoRT

The single most important source for this effort. It is the only public collection that stores method
conditions as *structured, machine-readable fields* alongside retention times.

- Repository: <https://github.com/michaelwitting/RepoRT>
- Paper: Kretschmer, F.; Harrieder, E.-M.; Hoffmann, M. A.; Böcker, S.; Witting, M.
  *RepoRT: a comprehensive repository for small molecule retention times.*
  Nat. Methods 2024. <https://doi.org/10.1038/s41592-023-02143-z>
- Archived release DOI: <https://doi.org/10.5281/zenodo.16267767>
  (Zenodo badge in the repository `README.md`)
- Submission portal: <https://rrt.boeckerlab.uni-jena.de/>

### Licence

**CC BY-SA 4.0.** Verified from the `LICENSE` file in the repository itself, whose first line reads
`Attribution-ShareAlike 4.0 International`.

> **ShareAlike is a real constraint, not a footnote.** If we redistribute a derived dataset (a cleaned,
> re-fitted, or `k`-converted version of RepoRT), that derivative must also be CC BY-SA 4.0. This does not
> restrict a *model* trained on the data under most readings, but it does restrict shipping curated data
> files. Flag this in the spec's licensing section. I did not obtain a legal opinion on whether trained
> model weights constitute an adapted work — **this is undetermined and should not be asserted either way.**

### Size *(measured, this snapshot)*

| Quantity | Value |
|---|---|
| Datasets with metadata | **421** |
| Retention time entries | **183,819** |
| Distinct compounds (standardised InChIKey) | **97,668** |
| Entries with a standardised SMILES | **183,819 (100%)** |
| Distinct compounds **excluding** the two METLIN SMRT datasets (`0186`, `0209`) | **18,872** |
| RP entries / RP distinct compounds | 173,950 / 97,403 |
| Method type | RP 343, HILIC 71, Other 5, blank 2 |
| Gradient / isocratic / no gradient file | 389 / **12** / 20 |
| Median compounds per dataset | 83 (max 79,973) |

*(published, for comparison: 373 datasets, 8,809 unique compounds, 88,325 RT entries, 49 columns — Nat.
Methods 2024 / the [ChemRxiv preprint](https://doi.org/10.26434/chemrxiv-2023-1mq73). The repository has
roughly doubled since publication.)*

**Note the 97,668 vs 18,872 gap.** 80,285 of the 183,819 entries are METLIN SMRT re-hosted as datasets
`0186` and `0209`. Any headline "RepoRT has ~100k compounds" claim is really "SMRT has 80k compounds, on
one column". For condition-varying work the honest denominator is **18,872 compounds across 419 datasets**.

### Structure availability

Excellent. Every RT row carries `smiles.std`, `inchi.std`, `inchikey.std` (PubChem-standardised), plus
ClassyFire taxonomy to six levels. Only 2.6% of entries lack a ClassyFire superclass. Per-dataset
precomputed descriptors and fingerprints (ECFP6, MACCS, PubChem) ship alongside, in both canonical and
isomeric variants.

### Format

Per dataset `NNNN/`, tab-separated: `NNNN_metadata.tsv` (one row, **185 fields**), `NNNN_info.tsv`
(provenance), `NNNN_gradient.tsv` (time / %A / %B / %C / %D / flow rate), `NNNN_rtdata_*.tsv`,
descriptors, fingerprints, and a PDF report. Also `NNNN_metadata.yaml`. Git LFS is used, so
`GIT_LFS_SKIP_SMUDGE=1` is enough to get all the TSVs (~824 MB working tree).

### Metadata completeness — in practice, not in principle

This is where the ticket's real question lies. Field presence *(measured, of 421)*:

| Field | Populated | Comment |
|---|---|---|
| `column.t0` | **421 / 421 (100%)** | but 50 are literally `0`; 371 are non-zero |
| `column.name` | 412 / 421 (98%) | |
| `column.usp.code` | 412 / 421 (98%) | |
| `column.particle.size` | 393 / 421 (93%) | |
| `column.length` | 392 / 421 (93%) | |
| `column.id` (internal diameter) | 392 / 421 (93%) | |
| `column.flowrate` | 379 / 421 (90%) | |
| `eluent.A.pH` | 373 / 421 (89%) | **but 52 of those are `0`** — a sentinel, not a measurement |
| `column.temperature` | **248 / 421 (59%)** | **the worst gap, and it is the one we care about** |

RepoRT is unusually honest about this: `info.tsv` has a free-text `missing information` column.
**179 of 421 datasets (43%) declare something missing**, and the declared gaps are dominated by exactly
the field we need:

```
159  temperature          25  column length        23  flow rate
 17  gradient             14  column id            14  particle size
 12  gradient and flow rate                        10  id and particle size
```

**Temperature is the single most-missing field in the entire repository.** 41% of datasets have no
temperature at all. This is a direct hit on ticket #11.

### Column dimensions and dwell volume

- **Column dimensions: yes.** Length, internal diameter and particle size are present for ~93% of datasets.
- **Dead time: yes, sort of.** `column.t0` is the only dead-time field. It is present for all 421 datasets,
  but 50 carry the value `0`. Where non-zero it appears to be a computed or submitter-supplied `t0` in
  minutes (e.g. `0.55125` across the whole BGC family, `1.9285714285714288` across the Souihi family —
  the repeating decimals suggest it is *derived* from column geometry and flow rate, not measured).
  **I could not determine whether these `t0` values are measured or calculated**, and the pattern of
  identical values across a whole submission strongly suggests calculated. Treat `t0` as a prior, not a
  measurement.
- **Dwell volume: NO.** *(measured)* I enumerated all 185 metadata field names and searched for
  `dwell`, `delay`, `system`, `void` and `extra`. **There is no field for any of them.**

> **This is a material finding for the instrument layer.** RepoRT records the gradient *programme* but not
> the gradient *delay*. For gradient data — 389 of 421 datasets — the mobile phase composition actually
> experienced by the solute at the column head is unknown by an unknown offset. Any `k(φ)` fitted from
> RepoRT gradient data absorbs the dwell volume into the fitted parameters.
>
> Consequences: (a) the instrument layer's dwell volume cannot be *validated* against RepoRT, only
> assumed; (b) `log k` values recovered by gradient inversion from RepoRT carry a systematic,
> per-dataset, unidentifiable bias; (c) this is a further argument for demoting `log k` RMSE and
> promoting elution order, which is far less sensitive to a common offset. Our existing decision to
> demote RMSE is *strengthened* by this.

### Condition coverage *(measured)*

- **Columns: 59 distinct `column.name` values.** Dominated by `Waters ACQUITY UPLC BEH C18` (75 datasets),
  `Waters CORTECS T3` (38), `Waters ACQUITY UPLC HSS T3` (26), `Merck Supelco Ascentis Express C18` (23).
- **USP codes:** `L1` 308, `L122` 23, `L11` 20, `L68` 18, `L3` 17, `L114` 13, `L43` 5, `L7` 4, `L10` 2,
  `L109` 2. **Overwhelmingly L1 (C18).** Phenyl/biphenyl (`L11`) is the only well-populated alternative
  RP chemistry. Ticket #2's curated
  column set will be validation-starved outside C18.
- **Temperature values:** `40` (130), `30` (51), `50` (24), `35` (10), `25` (9), `45` (5), `55` (4),
  then singletons at 24, 37.5, 38, 39, 41, 42.5, 60, 65 — plus one `-1` sentinel and **173 blank**.
  Effective usable range is **25–55 °C**, effectively **30/40/50**.
- **pH values (`eluent.A.pH`):** `3` (237!), `6` (24), `2` (8), `4` (8), `5` (8), `8` (8), `9` (6),
  `7` (6), `2.6` (4), then singletons — plus **52 zeros and 48 blanks**.
  **pH 3 accounts for 63% of datasets that state a pH.** This is the metabolomics default
  (0.1% formic acid), not a designed sweep.
- **Organic modifier (inferred from `eluent.B.acn` vs `eluent.B.meoh`):** ACN 256, MeOH 160, ambiguous 5.
  A genuinely balanced split — much better than expected.

### The condition-varying families (the part that actually matters)

Three submissions carry essentially all the designed condition variation.

#### (a) The BGC factorial — temperature and modifier. *The crown jewel.*

Datasets `0236`–`0259` and `0310`–`0357`, `source = "Dataset - BGC"`, authors E.-M. Harrieder,
M. Witting (Helmholtz Munich). **72 datasets.**

A deliberate factorial: **6 columns × T {30, 40, 50 °C} × modifier {ACN, MeOH} × ionisation mode
{pos, neg}**, everything else held constant — same gradient (4-point, 5%→99.5% B), same flow (0.4 mL/min),
same `t0` (0.55125), same pH (3), same column geometry (100 × 2.1 mm).

Columns: Waters ACQUITY UPLC HSS T3, HSS C18, BEH C18; Phenomenex Kinetex XB-C18; Waters CORTECS UPLC C18;
Restek Raptor Biphenyl. (MeOH arm covers only HSS C18, CORTECS C18 and Raptor Biphenyl.)

Shared compound cores *(measured, full intersection across the three temperatures)*:

| Family | pH | T levels | Shared core |
|---|---|---|---|
| BGC Kinetex XB-C18 | 3 | 30/40/50 | **499** |
| BGC HSS T3 | 3 | 30/40/50 | **495** |
| BGC BEH C18 | 3 | 30/40/50 | **493** |

Best pairwise same-column, same-pH, different-T overlaps reach **520 compounds**
(`0244`/`0246`, Kinetex XB-C18, 30 vs 40 °C). **416 such pairs exist across the repository.**

Modifier contrast, same column, same T, same pH — **156 pairs**, best overlaps ~400 compounds
(e.g. `0314`/`0322`, 50 °C, ACN vs MeOH, 406 shared compounds).

> ⚠ **Data-quality warning, discovered here.** In the BGC family the human-readable `info.tsv` `name`
> field and the `metadata.tsv` `column.name` field **disagree, and appear swapped**, for datasets
> `0310`–`0341`. Dataset `0310` is named `MSMS_pos_30_ACN_Raptor Biphenyl` but has
> `column.name = "Waters CORTECS UPLC C18"`; dataset `0326` is named `MSMS_pos_30_ACN_CORTECS C18` but
> has `column.name = "Restek Raptor Biphenyl"`.
>
> **`column.name` is the trustworthy field.** The supporting evidence: `0326` carries `usp.code = L11`
> and `particle.size = 2.7 µm`, which match a Restek Raptor Biphenyl; `0310` carries `L1` and `1.6 µm`,
> which match a Waters CORTECS C18. The structured fields are internally consistent; the free-text name
> is not. **Anyone joining on the `name` string will silently swap two columns' data.** Worth reporting
> upstream.

#### (b) The Souihi / Kruve family — pH.

Datasets `0275`–`0289`, `source = "Publication - Souihi"`,
<https://doi.org/10.1016/j.chroma.2022.462867>.

**pH {2, 3, 5, 8}** on Phenomenex Kinetex PS C18 and Kinetex EVO C18, at a **fixed 40 °C**, fixed
gradient, fixed `t0`.

Shared core across all four pH levels *(measured)*: **74 compounds** (PS C18), **40** (EVO C18).
Best pairwise pH overlap: 75 compounds. This is a small but genuinely designed pH sweep on a fixed
compound set — exactly the shape needed to test the ionisation layer, just small.

#### (c) The Aalizadeh family — pH, confounded.

Datasets `0382`–`0389`, `source = "Publication - Aalizadeh"`,
<https://doi.org/10.1021/acs.analchem.2c02893> (PMID 36347512).

`0382`/`0383` are a **pH 3.61 vs 6.2 pair on the same column at 30 °C with 1,673 shared compounds** —
by far the largest pH contrast in the repository. `0384`–`0389` add pH 2.6/2.62/3.71 over an 83-compound set.

⚠ **But the contrasts are confounded.** `0385` vs `0386` differ in pH (3.61 vs 2.6) *and* share a gradient,
which is fine — but `0387` vs `0389` differ in pH (3.71 vs 2.6) **and** temperature (35 vs 30 °C) **and**
gradient programme. Only `0382`/`0383` and `0385`/`0386` are clean single-factor pH contrasts. Check the
gradient table before using any pair.

### Replicates and uncertainty

*(measured)* **None.** I enumerated the headers of all 421 `rtdata_canonical_success.tsv` files. The
complete set of columns is `id, name, formula, rt, smiles.std, inchi.std, inchikey.std,
classyfire.{kingdom,superclass,class,subclass,level5,level6}, comment`. There is **no standard deviation,
no replicate count, no error estimate, no n**. Each compound in each dataset has exactly **one** retention
time, and the provenance of that number (single injection? mean of three?) is not recorded.

> **This is the sharpest constraint on the calibration ticket
> #15.** We cannot decompose predictive
> error into model error and measurement error using RepoRT, because RepoRT does not report measurement
> error. Our calibration target must therefore be stated against *total* observed error, and we should say
> so explicitly rather than implying we have isolated model uncertainty. The one usable external anchor is
> METLIN SMRT's repeatability figure (below).

---

## 2. METLIN SMRT

- Paper: Domingo-Almenara, X. *et al.* *The METLIN small molecule dataset for machine learning-based
  retention time prediction.* Nat. Commun. 2019, 10, 5811.
  <https://doi.org/10.1038/s41467-019-13680-7> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC6925099/>
- Data: figshare <https://doi.org/10.6084/m9.figshare.8038913>
- Licence: **CC BY 4.0** (per the paper's licence statement). Note this is *more* permissive than
  RepoRT's CC BY-SA.

### Confirmed: one column, one gradient

| | |
|---|---|
| Molecules | **80,038** (79,973 distinct InChIKeys as re-hosted in RepoRT `0186`/`0209`) |
| Column | Agilent Zorbax Extend-C18, **2.1 × 50 mm, 1.8 µm** — one column, one column only |
| Mobile phase | A: water + 0.1% formic acid; B: ACN + 0.1% formic acid — **one modifier, one pH (~2.7)** |
| Gradient | 5% B for 3 min → 50% B over 2 min → 85% B over 15 min → hold 3 min. **One programme.** |
| Flow | 100 µL/min |
| Temperature | **not stated in the paper** — undetermined |
| Dead volume | 40 µL |
| **Dwell volume** | **900 µL** — stated, and notably large relative to the column |
| Structures | PubChem CIDs + molfile/SDF; Dragon 7 descriptors and ECFP included |

**So: zero condition variation.** SMRT is a single point in `(φ-programme, pH, column, T)` space.
It cannot validate any of `k`'s condition arguments. The ticket's premise is confirmed.

### What SMRT *is* good for

1. **The solute-descriptor / QSPR layer, at scale.** 80k structures with a consistent readout is the
   largest structure→retention signal available. It is the right place to test whether predicted Abraham
   descriptors carry retention-relevant information, holding conditions fixed by construction.
2. **A measurement-uncertainty anchor — the only one found.** The paper reports RT variability over
   repeated analysis of **198 reference molecules** across intervals of ≥30 days: **mean 36 s, median 18 s**.
   Given the ~23 min gradient this is a usable floor for "how well could any model possibly do on this
   system", and it is the only replicate-derived uncertainty figure located in this survey.
   **Use it as the irreducible-error prior in #15.**
3. **A dwell-volume test case.** 900 µL dwell on a 2.1 × 50 mm column is a large delay; SMRT is a good
   stress test for whether the instrument layer's gradient-delay handling is correct, because the effect
   is big enough to see.
4. **Scaffold-split feasibility at scale** — see below.

### Caveat: known erroneous entries

An independent audit found systematic errors in SMRT: Kalinowska *et al.* / *Finding potentially erroneous
entries in METLIN SMRT*, J. Chromatogr. A 2025, <https://doi.org/10.1016/j.chroma.2025.465773>
(preprint: <https://doi.org/10.26434/chemrxiv-2024-lx6m1>). A large block of SMRT entries elutes in the
void volume and is not genuinely retained. **Filter void-volume entries before use** — they are
non-retained and will corrupt any `k` fit. I did not verify the exact count of affected entries from the
primary source; **this is undetermined here** and should be checked against the paper before relying on
a specific number.

---

## 3. Other sources assessed

### PredRet — subsumed, use RepoRT instead

- <https://predret.org/> · Stanstrup, J.; Neumann, S.; Vrhovšek, U. *PredRet: Prediction of Retention Time
  by Direct Mapping between Multiple Chromatographic Systems.* Anal. Chem. 2015, 87, 9421–9428.
  <https://doi.org/10.1021/acs.analchem.5b02287>
- Originally ~3,300 RT entries across **23 chromatographic systems**.
- **88 of RepoRT's 421 datasets carry `source = "PredRet"`** *(measured)* — PredRet has been ingested and
  re-standardised into RepoRT, with structures resolved and metadata regularised.
- **Verdict: no independent value.** Its method metadata is free-text system descriptions rather than
  structured fields, which is precisely the deficiency RepoRT was created to fix. Use RepoRT's copy.
- Its *projection* idea (map RTs between systems via shared compounds) is nonetheless methodologically
  relevant to our peak-tracking and reconciliation layers, and to the retention-order framing.

### MassBank / MoNA — not usable for this layer

- <https://massbank.eu/> · <https://mona.fiehnlab.ucdavis.edu/>
- MassBank does define structured chromatography fields, but **RT is recorded far less often than spectra**,
  and coverage of eluent composition and gradient programme is largely absent in MoNA — submissions
  typically give column name and instrument only, without eluents or gradient
  ([Xu et al., J. Cheminform. 2024](https://doi.org/10.1186/s13321-024-00905-1),
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11460055/>).
- **Verdict: not a condition-varying source.** Without the gradient programme we cannot invert to `k`.
  Possible future value as a breadth source for the descriptor layer only.

### PQRI / USP hydrophobic-subtraction-model column database — column axis, not retention data

- Snyder, L. R.; Dolan, J. W.; Carr, P. W. *The hydrophobic-subtraction model of reversed-phase column
  selectivity.* J. Chromatogr. A 2004, 1060, 77–116. <https://doi.org/10.1016/j.chroma.2004.08.121>
- USP tool: <https://www.usp.org/resources/pqri-approach-column-equiv-tool> ·
  <http://apps.usp.org/app/USPNF/columnsDB.html>
- Contains **H, S\*, A, B, C** per column, with **C reported at both pH 2.8 and pH 7.0**. The USP page
  states the database covers **368 columns**, built from "almost 3000 retention-time measurements for 150
  different test-compounds and 100 different alkylsilica columns". Other sources cite 588 or 750+ columns
  for later versions; **the discrepancy is unresolved here** and the map's "~750" figure should be treated
  as unverified.
- **Critically: it publishes fitted column parameters, not raw retention factors.** I could **not**
  establish that per-solute `k` values are downloadable for the full column set. It therefore
  **cannot validate `k`** — it can only supply the column-descriptor vector.
- Licence/terms not stated on the USP page. **Undetermined.** Assume restrictive until checked; USP
  content is generally not openly licensed. **This is a risk for ticket
  #2/#3** —
  the committed column axis may rest on data we cannot redistribute.
- A larger recent redetermination exists: *Improved hydrophobic subtraction model of reversed-phase liquid
  chromatography selectivity based on a large dataset with a focus on isomer selectivity*,
  J. Chromatogr. A 2024, <https://doi.org/10.1016/j.chroma.2024.465210>. **Paywalled; I could not
  determine whether its underlying retention dataset is public.** Worth one follow-up — if its raw `k`
  values are public it would be the single best unseen-column resource in existence.

### UFZ LSER database — QSPR layer only

- <http://www.ufz.de/lserd> · Ulrich, N.; Endo, S.; Brown, T. N.; Watanabe, N.; Bronner, G.; Abraham, M. H.;
  Goss, K.-U. *UFZ-LSER database v3.2.1*.
- ~6,800–8,000 compounds with **experimental** Abraham descriptors (one commonly cited filtered subset is
  N = 6,852).
- Some derived subsets have been redistributed under **CC0 on figshare**; the database's own terms were
  not verified. **Partly undetermined.**
- **Verdict: validates the QSPR layer in isolation** (predicted vs experimental E, S, A, B, V). Contains
  no chromatographic retention data and no method conditions. Already correctly scoped in
  `sources/README.md`.

### Manufacturer application notes — rejected

Surveyed at a high level and **not pursued**. They do vary conditions, but: no consistent structured
metadata, no machine-readable form, no SMILES, compound sets of 5–15, and copyright that is at best
unclear for redistribution. The extraction cost per usable data point is an order of magnitude worse than
RepoRT and the result would not be shareable. **Recommend explicitly ruling these out in the spec** rather
than leaving them as a vague "could also look at".

### Directly relevant methodological reference

Kretschmer, F.; Harrieder, E.-M.; Witting, M.; Böcker, S. *Times are changing but order matters:
Transferable prediction of small molecule liquid chromatography retention times.*
ChemRxiv 2024 (v3, Aug 2025). <https://doi.org/10.26434/chemrxiv-2024-wd5j8>

This paper is the closest published work to our framing: it argues RT changes massively across nominally
identical systems while **retention order is much better conserved**, and predicts a *retention order
index* conditioned on chromatographic conditions, using RepoRT. It reportedly analyses order conservation
as a function of solvent (finding conservation stronger for methanol-based systems) and temperature
(finding only some compounds flip order).

⚠ **I could not retrieve the full text** — both the ChemRxiv landing page and its PDF returned HTTP 403.
The characterisation above is from search-result summaries and the
[group's own announcement](https://bio.informatik.uni-jena.de/2024/12/times-are-changing-but-order-matters/),
**not from the primary source**, and its specific quantitative claims are **unverified**. Someone should
read it properly: it independently validates our choice of elution order as the primary metric, and it
may already answer part of ticket #11.
See [Follow-ups](#follow-ups).

---

## The critical question: does anything cross pH and temperature?

The ticket asks this explicitly, so here is the explicit measurement.

**Step 1 — is there any single column carrying both ≥2 pH levels and ≥2 temperatures?**
Exactly one: `Waters ACQUITY UPLC BEH C18`, with 21 qualifying datasets, pH {2.6, 3.0, 3.71} and
T {30, 35, 38, 40, 50, 65 °C}.

**Step 2 — do those datasets share any compounds?**
The intersection across all 21 is **0 compounds**. The pH variation and the temperature variation come
from different laboratories with disjoint compound sets and different gradient programmes.

**Step 3 — could a crossed design be assembled across studies?** I intersected the compound core of every
pH sweep with the compound core of every temperature sweep. The maximum overlap of any pH-sweep core with
any T-sweep core is **14 compounds** (Souihi/Kruve PS C18 pH sweep × BGC BEH C18 temperature sweep), and
those 14 are measured on *different columns with different gradients*, so even they are not a clean cell.

| pH-sweep core | × T-sweep core | Shared compounds |
|---|---|---|
| Souihi/Kruve PS C18 (pH 2/3/5/8) | BGC BEH C18 (30/40/50) | **14** |
| Souihi/Kruve PS C18 | BGC XB-C18 | 14 |
| Souihi/Kruve PS C18 | BGC HSS T3 | 11 |
| Aalizadeh BEH C18 (pH 2.6/3.71) | BGC BEH C18 | 10 |
| Aalizadeh Acclaim (pH 2.6/3.6/6.2) | BGC BEH C18 | 7 |
| Souihi/Kruve EVO C18 | BGC BEH C18 | 7 |

### Answer

- **Temperature, alone: YES for coverage, NO for order.** ~500 compounds × 3 temperatures × 6 columns
  × 2 modifiers, one designed factorial, one laboratory, everything else held constant. The *data* is
  there and a temperature term is falsifiable **in principle** — but measured elution-order reversal
  across 30 → 50 °C is only **0.9–1.6%**, so it is **not** falsifiable via our primary metric. Temperature
  must be validated on Δ`log k` instead. Range is narrow (30–50 °C) and single-lab, so we can claim
  *interpolation within 30–50 °C*, nothing more. See
  [Discriminative power](#discriminative-power-for-order-prediction).
- **pH, alone: YES, thin but by far the most valuable.** Best clean evidence is ~74 compounds across
  pH 2/3/5/8 at one temperature on two columns (Souihi/Kruve), plus a 1,673-compound pH 3.61 vs 6.2 pair
  (Aalizadeh `0382`/`0383`) which is only two levels. Enough to falsify a gross failure of the ionisation
  layer; **not** enough to fit or validate a rich pH response. Also note pH 8 on a Kinetex PS C18 — a
  positively-charged-surface phase — is not representative of ordinary C18 silanol behaviour.
  **This is nonetheless the only axis with real discriminative power for elution order (23.5% reversal
  overall, 25–31% on ionisables), and should therefore be the primary benchmark despite being the
  smallest.** The large Aalizadeh pair, by contrast, shows only 1.16% reversal — size does not help here.
- **pH × T interaction: NO. Unfalsifiable with public data.** Maximum crossed support is 14 compounds,
  and not cleanly.

### What we must therefore not claim

1. **Do not claim a validated `pH × T` interaction.** The model may contain one for physical reasons, but
   the spec must state that this term is **assumed, not validated**, and that no public data can currently
   test it. Ideally make it structurally separable so it can be switched off.
2. **Do not claim temperature validity outside 30–50 °C**, and note the evidence is single-laboratory.
3. **Do not claim pH validity above ~8 or below ~2**, nor on phases other than the two Kinetex phases plus
   the Acclaim/BEH pairs.
4. **Do not report calibrated uncertainty as decomposed into model vs measurement error** — no public
   condition-varying dataset reports replicates.

---

## Discriminative power for order prediction

Added after the #2 baseline finding:
only 7.7% of compound pairs reverse across φ = 0.2 → 0.7 in the blueprint paper's data, so the trivial
"predict order once, reuse everywhere" strategy already scores 92.3%. **An order-accuracy number is
meaningless without that naive baseline reported beside it**, and a dataset whose pairs rarely reorder
cannot separate a good model from a trivial one.

So the right question is not "does this dataset vary conditions?" but **"how often does elution order
actually change when it does?"** — because that reversal rate *is* the headroom available to any model.

### Method

For each pair of datasets differing in exactly one condition, over their shared compounds, I counted all
compound pairs and classified each as concordant or discordant in elution order (rank by `rt`). Compounds
eluting within `1.3 × t0` were dropped as unretained, and exact `rt` ties excluded. **Reversal rate =
discordant / (concordant + discordant)**; the naive baseline is its complement. This is computable
directly from RepoRT because the same compounds appear under many conditions — which is exactly the
property METLIN SMRT lacks.

### 1. Naive-baseline computability, per dataset

| Dataset | Compounds under ≥2 conditions? | Naive baseline computable? |
|---|---|---|
| **RepoRT — BGC family** | Yes: ~500 compounds × 3 T × 2 modifiers × 6 columns | **Yes**, ~60k–100k pairs per contrast |
| **RepoRT — Souihi/Kruve** | Yes: 71–75 compounds × 4 pH levels | **Yes**, ~2,500 pairs per contrast |
| **RepoRT — Aalizadeh** | Yes: 1,658 compounds × 2 pH levels | **Yes**, ~1.37M pairs |
| **RepoRT — Folberth/CORTECS T3** | Yes: ~50 compounds × many gradient programmes | Yes, but small and φ-only |
| **METLIN SMRT** | **No — one column, one gradient, one pH, one T** | **No. Zero cross-condition pairs.** |
| **PredRet** | Only via RepoRT's re-hosted copy | Same as RepoRT |
| **MassBank / MoNA** | No usable gradient metadata | No |
| **PQRI / USP HSM** | No raw retention factors published | No |
| **UFZ LSER** | No retention data at all | No |

> **METLIN SMRT cannot support our primary metric in any form.** It has no second condition, so there is
> no reversal to count and no naive baseline to beat. Its 80,038 compounds buy breadth for the descriptor
> layer and nothing at all for the condition layer. This should be stated plainly in the spec — the size
> of SMRT is otherwise very tempting and will mislead anyone who does not check.

### 2. Measured reversal rates

**Temperature — BGC family, same column / pH / modifier / gradient, only T differs** *(measured)*

| Column | 30→40 °C | 40→50 °C | 30→50 °C | n |
|---|---|---|---|---|
| Waters HSS T3 | 0.71% | 0.52% | **1.05%** | ~440 |
| Phenomenex Kinetex XB-C18 | 1.07% | 0.73% | **1.44%** | ~420 |
| Waters BEH C18 | 0.66% | 0.44% | **0.91%** | ~440 |
| Waters HSS C18 (MeOH) | 0.47% | 0.59% | **0.98%** | ~350 |
| Waters CORTECS C18 (MeOH) | 0.95% | 0.59% | **1.44%** | ~385 |
| Restek Raptor Biphenyl (MeOH) | 0.87% | 1.04% | **1.64%** | ~360 |

**Mean 30→50 °C reversal: 1.24%. Naive baseline: 98.8%.** Across ~95,000 pairs per contrast, so this is
not a small-sample artefact — it is a precise measurement of a very small effect.

**Organic modifier — ACN vs MeOH, same column / T / pH / gradient** *(measured)*

| Column | 40 °C | 50 °C | n |
|---|---|---|---|
| Waters CORTECS C18 | 5.66% | 5.56% | ~338 |
| Restek Raptor Biphenyl | **8.16%** | **7.78%** | ~337 |
| Waters HSS C18 | 6.27% | 6.19% | ~315 |

**~6–8% reversal, naive baseline ~93%.** Note this is *comparable to the 7.7% φ effect* the #2 agent
measured — modifier identity moves order about as much as a full φ sweep does.

**Column — same T / pH / modifier / gradient, 40 °C ACN, 15 pairwise contrasts** *(measured)*

Mean **3.37%**, range **1.00%–6.24%**. The structure is informative:

- Every contrast involving **Restek Raptor Biphenyl** (the only non-C18 chemistry) is at the top:
  4.58%–6.24%.
- Every **C18-vs-C18** contrast is at the bottom: CORTECS C18 vs HSS C18 is **1.00%**, BEH C18 vs
  CORTECS C18 **1.43%**, BEH C18 vs HSS C18 **1.54%**.

**Swapping one C18 for another barely reorders anything.** This is a direct, quantitative warning about
the unseen-column split: holding out a C18 and testing on another C18 leaves ~1–2% headroom, so a model
that ignores the column entirely would score ~98–99%.

**Gradient programme (φ trajectory) — same column / T / pH / modifier, different gradient** *(measured)*

336 qualifying pairs exist (largely the Folberth CORTECS T3 series). Measured reversal on the twelve
largest: **0.00%–3.30%**, mostly under 2.5%. Lower than the 7.7% the #2 agent measured for isocratic
φ = 0.2 → 0.7, which makes sense: these are all broadly similar full-range gradients, not a wide
isocratic φ sweep.

**pH — the outlier** *(measured)*

| Contrast | Column | n | Pairs | Reversal | Naive baseline |
|---|---|---|---|---|---|
| **pH 3 vs pH 8** | Kinetex PS C18 | 71 | 2,485 | **23.54%** | **76.46%** |
| pH 2 vs pH 8 | Kinetex PS C18 | 71 | 2,483 | **18.81%** | 81.19% |
| pH 3 vs pH 5 | Kinetex PS C18 | 70 | 2,414 | **16.40%** | 83.60% |
| pH 5 vs pH 8 | Kinetex PS C18 | 72 | 2,555 | 9.43% | 90.57% |
| pH 3 vs pH 2 | Kinetex PS C18 | 71 | 2,483 | 7.21% | 92.79% |
| pH 3.61 vs 2.6 | Acclaim RSLC 120 C18 | 77 | 2,920 | 4.76% | 95.24% |
| **pH 3.61 vs 6.2** | Acclaim RSLC 120 C18 | **1,658** | **1,372,268** | **1.16%** | 98.84% |

### 3. The ionisable/neutral asymmetry — confirmed and quantified

The coordinator predicted pH-driven reordering would be dominated by ionisable compounds crossing their
pKa. **Measured, and it holds strongly**, on the Souihi Kinetex PS C18 series:

| Contrast | N-containing (basic), no acid | Neutral (no acid, no N) | Ratio |
|---|---|---|---|
| pH 3 → 5 | **17.26%** (n=54) | 7.58% (n=12) | 2.3× |
| pH 3 → 8 | **24.95%** (n=54) | 9.09% (n=12) | 2.7× |
| pH 2 → 8 | **20.29%** (n=54) | **0.00%** (n=12) | ∞ |

Under a looser "any ionisable group" heuristic the pH 3 → 8 contrast reaches **31.08%** reversal for
ionisables against 9.06% for neutrals — a **3.4×** asymmetry.

⚠ **Caveat on the heuristic.** RDKit is not installed in this environment, so ionisability was assigned by
**crude SMILES substring matching**, not by pKa calculation or SMARTS. It is over-inclusive (it labels
78% of the Aalizadeh set ionisable, which is not credible) and the carboxylic-acid bucket was too small
to be meaningful (n = 4–5). **The direction and rough magnitude of the asymmetry are trustworthy because
they are consistent across all five independent pH contrasts; the precise percentages are not.**
Recompute with RDKit and a real pKa provider before quoting numbers in the spec.

### 4. Why the two pH results disagree — and why it matters

The 1,658-compound Acclaim pH 3.61 → 6.2 contrast shows only **1.16%** reversal, while the 71-compound
Kinetex pH 3 → 8 contrast shows **23.54%**. I checked whether this was a metadata artefact; it is not.
Both Acclaim datasets share an identical gradient, column, flow and temperature, differing only in
eluent A (5 mM ammonium formate + 0.01% formic acid at pH 3.61 vs 5 mM ammonium acetate at pH 6.2).

The difference is real and instructive:

1. **pH *range* matters more than pH *delta*.** 3.61 → 6.2 (Δ2.6) barely crosses the pKa of common bases,
   which stay protonated throughout. 3 → 8 (Δ5) crosses them decisively. A pH sweep that does not cross
   the sample's pKa values produces almost no reordering, regardless of how many compounds it contains.
2. **Sample composition matters more than sample size.** The Acclaim set is largely neutral environmental
   contaminants; the Kinetex set is deliberately ionisable-rich.
3. ⚠ **The buffer identity is confounded with pH** in the Acclaim pair (formate vs acetate), which is
   exactly the "buffer chemistry beyond pH" fog the map already flags. It is a clean pH contrast only if
   buffer identity is assumed irrelevant — which is the thing we would be trying to establish.

**Consequence: our benchmark must specify the pH *range* and an ionisable-rich compound set, not just
"varies pH".** A pH benchmark built on the large Acclaim pair would look rigorous and discriminate nothing.

### 5. Ranking by discriminative power

Ranked by headroom above the naive baseline, not by size.

| Rank | Contrast | Reversal | n | Breadth | Assessment |
|---|---|---|---|---|---|
| **1** | **pH 3 → 8, ionisable compounds** (Souihi, PS C18) | **25–31%** | 54 | 1 column, 1 T | **Highest discriminative power by a wide margin. Small and narrow, but the only contrast where a trivial model clearly fails.** |
| **2** | pH 3 → 8 / 2 → 8, all compounds (Souihi) | 19–24% | 71 | 2 columns, 1 T | Strong. The benchmark to build on. |
| **3** | pH 3 → 5 (Souihi) | 16.4% | 70 | 1 column | Strong; useful as a second, milder pH point. |
| **4** | **Modifier ACN → MeOH** (BGC) | **5.6–8.2%** | ~335 | 3 columns, 3 T | **Best signal-to-size trade-off. Comparable headroom to the φ effect, with 300+ compounds.** |
| 5 | Column, biphenyl vs C18 (BGC) | 4.6–6.2% | ~340 | 1 T | Useful, but rests on a single non-C18 column. |
| 6 | Gradient programme / φ path | 0–3.3% | ~50 | 1 column | Weak and small. |
| 7 | Column, C18 vs C18 (BGC) | 1.0–3.9% | ~330 | — | Nearly non-discriminating. |
| **8** | **Temperature 30 → 50 °C** (BGC) | **0.9–1.6%** | ~440 | 6 columns | **Effectively non-discriminating despite the best coverage in the repository.** |
| — | METLIN SMRT | n/a | 80,038 | 1 condition | **Cannot support the metric at all.** |

### 6. Answer to "does NO public dataset show substantial reordering?"

**No — that worst case does not hold. Public data does contain substantial reordering, but only on one
axis, and the supporting compound set is small.**

- **pH clears the bar comfortably**: 23.5% reversal overall and 25–31% on ionisables leaves real headroom.
  A trivial condition-blind model scores 76% there, so the metric can genuinely discriminate.
- ⚠ **But the entire discriminating benchmark rests on ~71 compounds, on two Phenomenex columns, at one
  temperature, from one publication.** That is a thin foundation for our primary acceptance metric, and it
  is a concentration risk: any systematic error in that one study propagates into our headline number.
- ⚠ **Kinetex PS C18 is a positively-charged-surface phase**, chosen by Souihi *et al.* precisely because
  it does unusual things with ionisable compounds. Reordering on it is likely **larger** than on an
  ordinary C18. Our headline discriminative number may therefore be optimistic for typical RP methods.
  **I could not determine how much of the 23.5% is phase-specific** — testing that needs a pH sweep on a
  plain C18 with an ionisable-rich set, which does not appear to exist publicly.
- **Temperature fails the bar outright.** At ~1% reversal, elution-order accuracy cannot distinguish a
  correct temperature model from one that ignores temperature entirely. This is not a data-volume
  problem — it is measured over ~95,000 pairs. It is a statement about the physics.

### 7. What this changes

1. **Every order-accuracy figure we report must carry its naive baseline beside it**, computed on the same
   compound set and contrast. A bare "94% pairwise order accuracy" is uninterpretable and, on the
   temperature axis, actively misleading. **Recommend making the paired figure a required output of the
   evaluation harness**, not a convention people are asked to remember.
2. **The primary benchmark should be the pH contrast**, with modifier second. Temperature and
   C18-vs-C18 column contrasts should be reported but explicitly labelled as low-discrimination.
3. **Temperature needs a different metric.** Order accuracy cannot falsify the temperature term. If we
   keep temperature in the model, it must be validated on retention *shift* (Δ`log k`, van 't Hoff
   behaviour) rather than order — which is one of the few places where the demoted `log k` metric is the
   right tool. This is a direct input to
   #11.
4. **Report reversal rate as a dataset property.** When the spec names a validation set, it should state
   that set's reversal rate, so the reader can see the headroom without recomputing it.

## Split strategy feasibility

The calibration ticket #15 needs honest
held-out data. Random splits would flatter us badly here, because RepoRT contains the *same compounds*
measured many times across conditions — a random row split leaks the compound into training almost
certainly.

**Recommended: all three of the following, reported separately.**

### Scaffold split — feasible, and easy

Every row has a standardised SMILES (100% coverage), so Murcko scaffolds are computable with RDKit
directly. ClassyFire taxonomy is present for 97.4% of entries and gives a ready-made, chemically
meaningful grouping at six levels. Superclass distribution is broad but imbalanced:

```
Organoheterocyclic 73,931 · Benzenoids 37,374 · Organic acids and derivatives 23,283
Lipids and lipid-like 16,572 · Organic oxygen 7,796 · Phenylpropanoids/polyketides 7,579
Nucleosides/nucleotides 5,311 · Organic nitrogen 4,140 · Alkaloids 1,539
```

**Do the scaffold split on the 18,872 non-SMRT compounds**, not the 97,668 — otherwise SMRT's scaffold
diversity dominates and the held-out set is evaluated almost entirely on a single column, which measures
the wrong thing.

### Unseen-column split — feasible, but only meaningfully at USP-code granularity

*(measured)* 37 columns carry ≥100 distinct compounds. Cross-column compound overlap is substantial,
which is what makes the split possible at all:

| Column | Distinct compounds | Datasets |
|---|---|---|
| Agilent ZORBAX Extend-C18 (= SMRT) | 79,973 | 2 |
| Waters ACQUITY UPLC BEH C18 | 13,637 | 75 |
| Thermo Scientific Acclaim RSLC 120 C18 | 2,500 | 7 |
| Waters ACQUITY UPLC HSS T3 | 1,608 | 26 |
| Waters ACQUITY UPLC HSS C18 | 966 | 17 |
| Restek Raptor Biphenyl | 770 | 16 |
| Phenomenex Kinetex XB-C18 | 753 | 9 |

Representative pairwise overlaps: BEH C18 ∩ Acclaim RSLC = 1,463; BEH C18 ∩ HSS T3 = 1,157;
BEH C18 ∩ HSS C18 = 924; HSS T3 ∩ HSS C18 = 760. Comfortably enough shared compounds to hold out a whole
column and still evaluate on compounds seen elsewhere.

⚠ **But the caveat is severe: 308 of 412 columns-with-a-code are USP `L1`.** Holding out one C18 and
testing on another C18 is a weak test — it mostly measures whether we can interpolate between very similar
phases. The only genuinely different chemistry with real support is `L11` (phenyl/biphenyl, 20 datasets,
Restek Raptor Biphenyl at 770 compounds). **A truly honest unseen-column split means holding out the
Raptor Biphenyl**, and there is essentially one such test available. Report it as a single, named,
adversarial test case rather than as a cross-validated mean.

### Unseen-condition split — still the one that matters, but pick the axis by reversal rate

The BGC factorial supports the cleanest *mechanics* we have — compound set fixed, one condition varied.
But the measured reversal rates say the axes are not interchangeable:

> **Headline test: train on one pH, test on another** (Souihi pH 3 → 8, ionisable-rich). ~24% headroom.
> **Second test: train on ACN, test on MeOH** (BGC, same columns and compounds). ~6–8% headroom, ~335
> compounds — the best combination of signal and size.
> **Third: train on five BGC columns, test on the Raptor Biphenyl.** ~5–6% headroom.
>
> ⚠ **Do not headline "train 30 + 50 °C, test 40 °C"** — I recommended exactly this before measuring, and
> the measurement kills it. Order barely moves with temperature (~1%), so the test is near-trivial and
> would flatter us just as badly as a random split.

**Every one of these must report its naive baseline alongside**, per the section above.

⚠ Note the tension: the axis with the best split *mechanics* (temperature — 6 columns, ~440 compounds,
perfectly controlled) has the least discriminative power, and the axis with the most discriminative power
(pH) has the weakest mechanics (71 compounds, 2 columns, one lab). **There is no public contrast that is
both well-powered and strongly discriminating.** That is the central limitation of this whole survey.

### Splits to avoid

- **Random row splits** — leak compounds across conditions; will look excellent and mean nothing.
- **Random dataset splits** — leak columns and conditions; a held-out BEH C18 dataset is not held out in
  any useful sense when 75 BEH C18 datasets are in training.
- **Any split that mixes SMRT into a condition-varying evaluation** — SMRT's 80k compounds at one
  condition will dominate every aggregate metric and mask condition-layer failure entirely.

---

## Recommendations

1. **Adopt RepoRT as the retention layer's ground truth**, pinned to a specific commit or Zenodo release
   for reproducibility. Record the CC BY-SA obligation in the spec.
2. **Make the Souihi pH sweep (`0275`–`0289`) the primary order-accuracy benchmark**, and the BGC
   ACN/MeOH contrast the secondary one. Join on `column.name`, **never** on the `info.tsv` `name` field.
3. **Make the naive baseline a mandatory, machine-computed companion to every order-accuracy figure.**
   Build it into the evaluation harness so it cannot be omitted.
4. **Use METLIN SMRT for the descriptor layer only**, filtered for void-volume entries, and lift its
   36 s / 18 s repeatability figure as the irreducible-error prior for calibration. State explicitly that
   it cannot support the order metric.
5. **Validate temperature on Δ`log k` / van 't Hoff behaviour, not on elution order** — order reversal is
   ~1% and cannot falsify the term. Scope to 30–50 °C interpolation, single-laboratory. Feeds
   #11.
6. **Treat the `pH × T` interaction as an assumed, unvalidated term** — structurally separable, explicitly
   flagged, not claimed.
7. **Do not model dwell volume as validated against public data.** Fit it, or take it as a user-supplied
   prior, and document that RepoRT cannot check it.
8. **Report three splits separately** — scaffold, unseen-column (headline: Raptor Biphenyl, *not* a second
   C18), and unseen-condition (headline: hold out a pH level; secondary: hold out MeOH). Never report a
   random split, and never headline the 40 °C hold-out.
9. **Record the concentration risk**: our most discriminating benchmark is 71 compounds from one
   publication on one unusual stationary phase. Consider whether generating a small in-house pH sweep on
   a plain C18 is worth it — it would de-risk the primary acceptance metric more than any other
   experiment identified in this survey.

---

## What I could not determine

Stated plainly, so nobody mistakes a gap for a finding.

- Whether RepoRT's `column.t0` values are **measured or calculated**. Repeating decimals across whole
  submissions strongly suggest calculated, but this is inference, not verification.
- The **licence and redistribution terms of the USP/PQRI column database**, and its true current column
  count (368 per USP's own page vs 588 vs "750+" in secondary sources — unresolved).
- Whether the **2024 improved-HSM dataset**'s raw retention factors are public (paywalled).
- The **exact count of void-volume/erroneous entries in METLIN SMRT** (audit paper not read in full).
- **METLIN SMRT's column temperature** — not stated in the paper.
- The **quantitative order-conservation results** in Kretschmer *et al.* 2024 (403 on both ChemRxiv URLs).
- Whether **model weights trained on CC BY-SA data** constitute an adapted work under the ShareAlike term.
- I did **not** systematically survey manufacturer application notes or the wider LSER supplementary-data
  literature compound-by-compound; both were assessed at survey level and judged low-yield.
- **How much of the 23.5% pH reversal is specific to the Kinetex PS C18 phase.** It has a
  positively-charged surface, deliberately chosen for unusual ionisable-compound behaviour, so the figure
  may be optimistic for ordinary C18. No public pH sweep on a plain C18 with an ionisable-rich compound
  set was found to check it against.
- **Ionisability class assignment is crude** — SMILES substring matching, no RDKit available in this
  environment, no pKa provider. The ionisable/neutral asymmetry is directionally solid (consistent across
  five independent contrasts) but the specific percentages need recomputing with proper cheminformatics.
- Whether the **buffer identity confound** in the Aalizadeh pH pair (ammonium formate vs ammonium acetate)
  accounts for any of its behaviour. It is confounded with pH by construction and cannot be separated.

## Follow-ups

- **Read Kretschmer *et al.* 2024** (<https://doi.org/10.26434/chemrxiv-2024-wd5j8>) via institutional
  access. It bears directly on #11 and
  on the elution-order metric choice.
- **Report the `0310`–`0341` column-name swap upstream** to
  <https://github.com/michaelwitting/RepoRT/issues>.
- **Check whether the 2024 improved-HSM raw dataset is public** — it would materially change the
  unseen-column split story, which is currently our weakest.
- **Resolve the USP/PQRI licence question** before the column axis is committed in
  #2.

---

## Reproducing the measurements

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 https://github.com/michaelwitting/RepoRT.git
# aggregate processed_data/NNNN/NNNN_metadata.tsv (185 fields, one row each),
# NNNN_info.tsv, NNNN_gradient.tsv and NNNN_rtdata_canonical_success.tsv across all 421 datasets.
```

Reversal rates: for a dataset pair differing in one condition, take the shared InChIKeys, drop any
compound with `rt <= 1.3 * column.t0` in either dataset, then over all compound pairs count
`sign(rt_A[i] - rt_A[j]) != sign(rt_B[i] - rt_B[j])` as a reversal, excluding exact ties.
`reversal_rate = discordant / (concordant + discordant)`; `naive_baseline = 1 - reversal_rate`.
Pure Python is sufficient — the largest contrast is ~1.4M pairs.

All *(measured)* figures in this document come from that snapshot at commit
`9de8d603377bbb6cc0f74e250eb7533fd8874df1` (2026-07-09). They will drift as the repository grows; pin a
release before quoting them in the spec.
