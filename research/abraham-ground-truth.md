# Abraham descriptor ground truth: coverage, licence, and applicability domain

Resolves #5. Part of the map, #1.

This note establishes the ground truth for the **QSPR layer** — experimental Abraham solute
descriptors `E, S, A, B, V` — and characterises its chemical space against the scope of this
effort (RPLC, drug-like **ionisable** compounds).

---

## Headline

**The "experimental" ground truth is substantially not experimental.** Traced value-by-value
through SoluteDB's per-descriptor reference codes (§3b): **16.7% of `E` and 22.1% of `B` values
are literally commercial-software output** (ACD Absolv, ACD ChemSketch refractive index, ClogP),
and a further **~48% of `E`, ~48% of `S`, 23% of `A` and 18% of `B` are "estimated by comparison
to closely related compounds"** — expert interpolation, not measurement. Requiring all four of
`E, S, A, B` to be free of software output and analogy-estimation leaves **1,434 of 8,545 records
(16.8%)**, and that clean subset contains **12 drug-like ionisable molecules**. The coordinator's
suspicion is confirmed and then some.

**The experimental descriptor corpus is one dataset with two faces.** Everything in circulation
traces back to Michael Abraham's own compilation ("Absolv"). The two we can actually obtain —
the UFZ-LSER copy and MIT's SoluteDB — share **92% of the UFZ set by InChIKey** (measured, see
[Method](#method)). There is no second, independent ground truth to cross-validate against.

**It is not drug chemistry.** Median MW is **174 g/mol**; **25%** of entries are above MW 250;
**15%** carry a group that is ionised between pH 2 and 8; and only **228 compounds (3.2%)** are
both drug-sized and ionisable. Benzene is the Murcko scaffold of **26%** of the set.

**And the target chemistry sits in the tail.** For those 228 drug-like ionisable compounds the
*median* `S` is 2.35, `B` is 1.96, `V` is 2.75 — which are respectively the **93rd, 96th and 94th
percentile** of the full descriptor distributions. Half of the chemistry we care about lives in
the top ~6% of the data by descriptor magnitude.

**Consequence, and the point of the ticket:** a QSPR RMSE quoted over the whole compilation is a
statement about small neutral solutes — and partly a statement about agreement with ACD/Labs'
predictor. It does not transfer. Any error figure we accept — ours or the literature's — must be
reported **on the drug-like ionisable subset separately**, and that subset is only ~228 molecules
before provenance filtering and **12** after. That is not a validation set. It is a spot check.

**Licence matters differ sharply between the two copies.** UFZ-LSER is free for scientific and
private use with **commercial use explicitly not permitted**. MIT's SoluteDB is **CC-BY-4.0**.
For a spec that may become a product, SoluteDB is the safe primary and UFZ the cross-check.

---

## 1. What the sources actually are

### 1.1 Abraham's own compilation ("Absolv") — the root

All the bulk descriptor data descends from Abraham's compilation, commercialised as **Absolv**
(now ACD/Labs' *Absolv* / ADME Suite predictor). It is not itself distributed as an open dataset;
it reaches us only through the two derivative copies below.

Marchetto et al. name it directly: they used "the so-called Abraham Absolv data set (taken from the
UFZ-LSER database) which, at the time of this work, comprised LSER solute parameters for **7881
small molecules**. This data set is the result of Abraham's work in the field of LSER"
([Anal. Chem. 2025, 97, 6991–7001](https://doi.org/10.1021/acs.analchem.4c03466), p. 6994).

Note the phrase *small molecules* is used by the authors themselves.

### 1.2 UFZ-LSER database

- Site: <https://www.ufz.de/lserd> (redirects to <https://web.app.ufz.de/compbc/lserd/public/start/>).
- Cited in the literature as: Ulrich, N.; Endo, S.; Brown, T. N.; Watanabe, N.; Bronner, G.;
  Abraham, M. H.; Goss, K.-U. *UFZ-LSER Database V 3.2.1*, Helmholtz Centre for Environmental
  Research-UFZ, Leipzig, 2017 (ref. 62 of the Marchetto paper).
- The live site now self-cites as **"UFZ-LSER database v4.0 [Internet]. Leipzig: Helmholtz Zentrum
  für Umweltforschung - UFZ; 2025"** (<https://web.app.ufz.de/compbc/lserd/public/citation/>) and
  the start page banner reads **v4.1.1 (2025)**.
- The start page reports **411,806 chemical entries** — but that is the *calculated* descriptor
  space (SMILES in, descriptors predicted). The **experimental** descriptor subset is the ~6,000–8,000
  Absolv records. The UI distinguishes "Experimental descriptors" from "Calculated descriptors
  (only SMILES input)". Do not quote 411,806 as experimental coverage.
- The start page states the tooling is for **neutral organic compounds** and that ionic-species
  calculation is explicitly limited.

**Licence — this is the load-bearing bit.** From
<https://web.app.ufz.de/compbc/lserd/public/Public/policy/>:

> "The user is granted the non-exclusive, perpetual, non-transferable and non-sublicensable right
> to use the database free of charge for scientific and private purposes, subject to attribution
> of the author."

> "Any use for commercial purposes is not permitted."

> "The operator assumes no liability or guarantee for the content of the data, in particular with
> regard to the topicality, correctness and completeness of the data provided."

Redistribution and bulk download are not addressed explicitly, but *non-transferable* and
*non-sublicensable* mean we should not assume a right to redistribute. **No bulk download button
is offered on the site.** (Could not determine: whether a bulk export is available on request.)

**Obtainable in practice via:** the supplementary data of Ulrich & Ebert's DNN paper, published at
<https://github.com/nadinulrich/solute_descriptor_prediction> (repo licence: **MIT**), file
`Dataset_and_predictions.xlsx` (4.8 MB). Sheet `dataset_and_outliers` holds **7,242 rows** with
`Name, CAS-RN, SMILES` and experimental `E, S, A, B, V, L`, plus an `outlier` flag (**877 flagged**)
and a free-text "possible reason for outlier". This is the file analysed below.

> ⚠️ **Licence ambiguity, flagged not resolved.** The GitHub repo carries an MIT licence, which
> covers the code. The descriptor *values* in that spreadsheet originate from Absolv/UFZ-LSER,
> whose own terms forbid commercial use. An MIT header on a repo does not relicense third-party
> data inside it. Treat this file as research-use-only until someone with authority says otherwise.

### 1.3 SoluteDB (MIT / Green group) — the CC-BY route

Chung, Y.; Vermeire, F. H.; Wu, H.; Walker, P. J.; **Abraham, M. H.**; Green, W. H.
*Group Contribution and Machine Learning Approaches to Predict Abraham Solute Parameters, Solvation
Free Energy, and Solvation Enthalpy*, J. Chem. Inf. Model. 2022, 62(3), 433–446,
<https://doi.org/10.1021/acs.jcim.1c01103>. Preprint:
<https://chemrxiv.org/engage/chemrxiv/article-details/61d33656d6dcc267874e59e0>.

Data: **Zenodo record <https://zenodo.org/records/5792296>**, file `Solvation_data-1.0.0.zip`
(9.1 MB), licence **Creative Commons Attribution 4.0 International (CC-BY-4.0)**. Downloaded and
inspected for this note.

- `all_data/SoluteDB_all.xlsx` — **8,481 rows**, columns `SMILES, InChI, E, S, A, B, L`.
- `selected_data_for_this_work/SoluteDB_selected_data.xlsx` — **8,366 rows**; the README states
  these are "limited to **neutral** solute compounds containing H, C, N, O, S, P, F, Cl, Br, or I
  atoms and nonionic liquid solvents."
- Per-descriptor non-null counts in `SoluteDB_all`: `E` 8,261 · `S` 7,730 · `A` 8,260 · `B` 7,469 ·
  `L` 7,110. **All five present: 6,736. All of `E,S,A,B`: 7,425.**
- Abraham is a co-author; the README calls this "the **in-house** Abraham solute parameter data",
  confirming it is the same root compilation, re-curated.
- **No `V` column** — because `V` is computed, not measured (see §3).
- **No uncertainty columns of any kind.** (`dGsolv`/`dHsolv` sets in the same archive *do* carry
  `_std` and `no. of data` columns; the solute-parameter set does not.)

### 1.4 Acree's "compilation"

This is worth stating plainly because the ticket asks for it as if it were a single database:
**it is not one.** William E. Acree Jr. (Univ. of North Texas) publishes a continuous stream of
single-compound and small-series determinations — e.g. "Determination of Abraham model solute
descriptors for the monomeric and dimeric forms of *trans*-cinnamic acid"
(<https://doi.org/10.1186/s13065-015-0080-9>, open access), "Determination of Abraham Model Solute
Descriptors for 62 Additional C10–C13 Methyl- and Ethyl-Branched Alkanes"
(<https://doi.org/10.3390/liquids2030007>), and many more in *Phys. Chem. Liq.*, *J. Solution Chem.*,
*Thermo*, *Liquids*. These values flow **into** Abraham's compilation and thence into §1.2/§1.3.

Acree's group *does* publish bulk datasets on Figshare under **CC0**, but the ones I could confirm
are **solvent** coefficients (the `c, e, s, a, b, v/l` system constants) and an *enthalpy of
solvation* dataset — e.g. <https://figshare.com/articles/dataset/Acree_Enthalpy_of_Solvation_Dataset/1572326>
— **not** a consolidated solute-descriptor table.

> **Could not determine:** whether a single consolidated, openly-licensed Acree *solute descriptor*
> table exists distinct from Absolv/SoluteDB. I found no such artefact. Treat "the Acree
> compilation" in the map as a **synonym for the Abraham/Absolv lineage**, not a third source.

### 1.5 Poole's Wayne State dataset — small, but the cleanest

Poole, C. F. *J. Chromatogr. A* 2020, 1617, 460841 (ref. 64 of Marchetto). Marchetto describe it as
"the Wayne State University experimental descriptor data set, which lists LSER solute parameters for
**several hundred solutes** collected in a single laboratory via standardized procedures, thus
minimizing experimental variations."

This is the one source where inter-laboratory variance is controlled. It is far too small to train
on, but it is the best candidate for an **independent held-out test set** — provided we can obtain
it (paywalled; not checked for an open copy).

---

## 2. Which numbers to quote for coverage

Different papers quote different counts because they curate differently. Use these, with the
curation attached:

| Figure | Meaning | Source |
|---|---|---|
| ~8,000 | "about 8000 experimentally determined solute descriptors are available" | [Ulrich et al., ACS Environ. Au](https://doi.org/10.1021/acsenvironau.6c00063), §1 |
| ~6,300 | chemicals with a *usable set* of descriptors after requiring completeness | ibid. |
| 6,364 / 6,358 / 6,347 / 6,111 / 5,572 | per-descriptor counts for `S / A / E / B / L` in their curated set | ibid., §3.2 |
| 7,881 | Absolv records taken from UFZ-LSER | [Marchetto et al.](https://doi.org/10.1021/acs.analchem.4c03466), p. 6994 |
| 6,437 → 6,401 | after removing incomplete/unSMILES-able/duplicate entries **and restricting to MW 80–400 g/mol** | ibid. |
| 7,242 (877 flagged outlier) | rows in the published UFZ/Absolv spreadsheet | `Dataset_and_predictions.xlsx` |
| 8,481 / 8,366 | SoluteDB all / neutral-filtered | [Zenodo 5792296](https://zenodo.org/records/5792296) |

**Order of magnitude to plan against: ~6,500–7,500 compounds with a complete `E,S,A,B,V` set.**

---

## 3. How the descriptors were measured — and why "measured" is the wrong word for most of them

This changes what "ground truth" means, so it is set out explicitly.

- **`V` — McGowan characteristic volume. Not measured at all.** It is computed from atomic volume
  increments minus a bond-count term (McGowan, *J. Chem. Technol. Biotechnol.* 1978, 28, 599; see
  also Zhao, Abraham, Zissimos, *J. Chem. Inf. Comput. Sci.* 2003, 43, 1848, ref. 60–61 of
  Marchetto). Marchetto state it directly: "The fifth solute parameter present in eq 2, `V`, can be
  determined directly from the molecular structure." **Predicting `V` with a QSPR is a category
  error — compute it.** SoluteDB does not even ship a `V` column for this reason.

- **`E` — excess molar refraction. Semi-experimental.** Derived from the refractive index at 293 K
  for liquids; for solids it is estimated by fragment addition from compounds with known `E`, or
  taken from commercial prediction (Absolv/ACD). So a meaningful share of tabulated `E` values are
  themselves *estimates*, not measurements
  ([Chemistry Central J. 2015, 9:13](https://doi.org/10.1186/s13065-015-0080-9)).

- **`S`, `A`, `B`, `L` — obtained by regression, not by direct measurement.** For a given solute,
  a large over-determined system of Abraham equations (water–solvent partition coefficients,
  aqueous/organic solubilities, GLC retention on multiple stationary phases, RP-HPLC retention) is
  solved simultaneously for the unknown descriptors. See Abraham, Ibrahim, Zissimos,
  "Determination of sets of solute descriptors from chromatographic measurements",
  *J. Chromatogr. A* 2004, 1037, 29–47 (<https://pubmed.ncbi.nlm.nih.gov/15214659/>), which covers
  GLC-derived and RP-HPLC-derived (Snyder–Dolan) routes; and, for the solubility route, the worked
  example in [Chemistry Central J. 2015, 9:13](https://doi.org/10.1186/s13065-015-0080-9).

**What this means for us:**

1. The "ground truth" is a **fitted quantity**, conditioned on the LSER model form. Validating a
   QSPR against it measures agreement with Abraham's regressions, not with nature.
2. Errors are **correlated across descriptors** for a given solute — they come out of one joint
   solve. A per-descriptor independent-error assumption in the uncertainty budget is wrong.
3. Values are **method-heterogeneous**: a descriptor set fitted from 45 solubility equations and
   one fitted from a handful of GC retention indices are not the same quality, and the compilations
   do not tell you which is which.

### Uncertainty: not reported in either compilation

- Neither `Dataset_and_predictions.xlsx` nor `SoluteDB_all.xlsx` carries a standard error, standard
  deviation, or n-of-measurements column for any descriptor. They are **bare point values**.
- The *primary* determination papers do report fit quality — e.g. the cinnamic acid work reports
  "a standard deviation (SD) of about 0.1 log units" over 45 equations for the monomer and 0.087
  over 20 for the dimer ([doi:10.1186/s13065-015-0080-9](https://doi.org/10.1186/s13065-015-0080-9)).
  That is the SD of the *fitted property*, not a per-descriptor standard error, and it is lost on
  the way into the compilation.
- **A useful sanity anchor for the uncertainty budget:** the Ulrich & Ebert / UFZ spreadsheet flags
  **877 of 7,242 rows (12%)** as outliers with a stated reason — i.e. the curators themselves judge
  roughly one entry in eight to be suspect.

> **Practical consequence.** We cannot inherit a descriptor uncertainty from the data. We must
> either (a) assign a nominal per-descriptor uncertainty from the primary literature (~0.05–0.10
> descriptor units for `S`, `A`, `B`; larger for `E` on solids), or (b) treat descriptor uncertainty
> as a free parameter calibrated on downstream retention error. Option (b) fits the map's
> "benchmark on downstream error" pattern already adopted for `pKa`.

---

## 3b. Ground-truth integrity — is it measured, or is it ACD/Labs' model output?

This is the question that gates everything else, so it is answered from the data itself rather
than from what papers assert.

### 3b.1 The naming trap

"Absolv" is **two things at once**, and the literature conflates them:

1. **Abraham's compilation** — the accumulated table of descriptor sets he and collaborators built
   up over ~35 years.
2. **A commercial prediction program** — originally Pharma Algorithms' *Absolv*, now
   **ACD/Labs Absolv / ADME Suite**, which *predicts* Abraham descriptors from structure.

When Marchetto et al. say they used "the Abraham Absolv data set… comprised LSER solute parameters
for 7881 small molecules… the result of Abraham's work in the field of LSER", they mean sense (1).
But sense (1) **contains values produced by sense (2)**, and neither the paper nor the UFZ copy
tells you which.

### 3b.2 UFZ-LSER cannot tell you — SoluteDB can

The UFZ-LSER help page is explicit that the stored descriptors are of mixed origin:

> "LSER descriptors stored in this database have either been previously published in peer reviewed
> scientific papers **or were taken from the ABSOLV data base**" … "for roughly 1000 chemicals the
> stored descriptors are not complete" … "We neither take responsibility for accuracy of the
> predictions … nor do we guarantee correctness of stored values in comparison to those in the
> original papers or in the ABSOLV database."
> (<http://web.app.ufz.de/compbc/lserd/public/Public/help/>)

It also separately offers QSPR-**predicted** descriptors for any SMILES, and warns "the accuracy of
QSPR-predicted descriptors is inferior to those that were derived from experimental data". So the
site does distinguish literature / ABSOLV / QSPR-predicted in its UI.

**But the exported data does not.** In `Dataset_and_predictions.xlsx` (the UFZ/Absolv extract used
by Ulrich & Ebert), the `Literature` column has exactly **one value for all 7,241 rows: "Abraham
Absolv"**. There is no per-compound, let alone per-descriptor, citation. **Taken alone, the UFZ
export is provenance-opaque.**

> **Could not determine:** the exact count of UFZ-LSER records whose descriptors came from
> peer-reviewed papers versus from the ABSOLV table. The site publishes no such breakdown and
> exposes no bulk export. The "~5,000 literature vs 7,881 used" figure in circulation could not be
> confirmed from a primary source — but §3b.3 below establishes the same conclusion by a stronger
> route.

**SoluteDB is the one copy that carries provenance.** `SoluteDB_logKw_logPow_all.xlsx` in
[Zenodo 5792296](https://zenodo.org/records/5792296) has **8,545 records × per-descriptor
reference codes**: `E Ref, S Ref, A Ref, B Ref, BO Ref, L Ref`, resolved against a `references`
sheet of **1,046 distinct reference entries**. 96–99% of values carry a code. Because Abraham
co-authored SoluteDB and it is the same lineage as the UFZ copy (92% InChIKey overlap, §4), this
provenance can be read across to the UFZ data.

**This alone is a reason to make SoluteDB the primary source: it is the only copy you can audit.**

### 3b.3 Reading the reference codes — the answer

Every reference code was classified by its own text into four provenance tiers, then value counts
computed per descriptor. Tiers:

1. **software-predicted** — text names a prediction program (ACD Absolv, PharmaAlgorithms Absolv,
   ChemSketch, ClogP, UFZ-LSER)
2. **estimated by analogy** — "estimated", "estimated by comparison to closely related compounds",
   "estimated from steroid fragments", "smoothed across the series"
3. **derived by regression from measurements** — "Calculated using the Solver program",
   "back-calculated from GLC data and water–solvent distributions", "calculated from refractive index"
4. **literature citation** — a journal reference or named unpublished measurement work

| descriptor | values | 1 software-predicted | 2 estimated by analogy | 3 regression-derived | 4 literature-cited | no code |
|---|---|---|---|---|---|---|
| `E` | 8,308 | **1,391 (16.7%)** | **4,037 (48.6%)** | 2,339 (28.2%) | 425 (5.1%) | 112 (1.3%) |
| `S` | 7,771 | 0 (0.0%) | **3,762 (48.4%)** | 1,063 (13.7%) | 2,679 (34.5%) | 239 (3.1%) |
| `A` | 8,317 | 3 (0.0%) | 1,927 (23.2%) | 772 (9.3%) | **5,246 (63.1%)** | 354 (4.3%) |
| `B` | 7,509 | **1,662 (22.1%)** | 1,361 (18.1%) | 968 (12.9%) | 3,267 (43.5%) | 235 (3.1%) |

The specific reference texts that dominate are worth reading verbatim, because they are damning:

- **`E`**: `E02` *"Estimated values on the basis of values calculated in R1"* — **n = 3,875 (47% of
  all `E`)**, and reference `R1` is **not defined anywhere in the references sheet**, so the
  provenance is untraceable. `E19` *"Value calculated by ACD Absolv"* — n = 911. `E09`
  *"Refractive index calculated using ChemSketch 2.0, Demo, from Advanced Chemistry Development
  Ltd."* — n = 462. `E26` *"Value calculated by UFZ-LSER"*. Only `E01` (n = 2,299), *"Calculated
  from refractive index as outlined by Abraham et al."*, is grounded in a real measurement.
- **`S`**: `S03` *"Estimated by comparison to closely related compounds"* — **n = 3,509 (45% of all
  `S`)**.
- **`B`**: `B10` *"Calculated from a ClogP estimate of water–octanol partition"* — **n = 1,657
  (22%)**. ClogP is itself a prediction, so this is `B` derived from a predicted logP. `B59`
  *"From logP(oct) calculated by ACD"*. `B08` *"From logP(octanol) only"* (n = 1,020) is measured
  but rests on a **single** partition datum, which is the weakest possible determination.
- **`A`** is the healthiest: `A01` (n = 4,531) is Abraham's 1989 hydrogen-bond acidity measurement
  paper (*J. Chem. Soc. Perkin Trans. 2*, 1989, 699). But recall §4.5: **56% of `A` values are
  exactly zero**, so most of that "measurement" is the trivially-known statement that a molecule
  has no acidic hydrogen.

**Direct answer to the coordinator's question: yes, the ground truth is contaminated with model
output — most severely for `E` and `B`, which are exactly the two descriptors on which Marchetto
et al. report their best QSPR scores (`E` R² 0.98 / RMSEP 0.10; `B` R² 0.95 / RMSEP 0.11).** A
model that reproduces `E` to RMSEP 0.10 is in substantial part reproducing (a) ACD Absolv's
predictor and (b) Abraham's own structural-analogy interpolation. Both are smooth functions of
structure, which is precisely why a structure-based QSPR fits them so well. **The high `E` and `B`
R² values are partly a measure of how learnable someone else's model is, not of chemical accuracy.**

`S` is not software-contaminated, but 48% of it is analogy-estimated — same objection, human
rather than machine. `A` is the only descriptor with a majority-measured basis, and it is
degenerate (56% zeros).

### 3b.4 Is some of a record measured and some derived? Yes — that is the normal case

The reference codes are **per descriptor, not per compound**. A single record routinely mixes
tiers: e.g. `A` from Abraham's 1989 measurement paper, `S` estimated by analogy, `E` from ACD
ChemSketch, `B` from a ClogP estimate. This is the most misleading case the coordinator anticipated,
and it is the modal case, not an edge case. **There is no compound-level "this record is
experimental" flag anywhere in any of these sources.**

It also means the errors are *not* independent across descriptors for a compound — see §3, point 2.

### 3b.5 How big is a genuinely clean subset?

Requiring all four of `E, S, A, B` present **and** none of the four sourced from software prediction
or analogy estimation:

- **1,434 records of 8,545 (16.8%).**
- Relaxing `E` (whose refractive-index route is arguably fine, and which is the worst-contaminated):
  requiring only `S, A, B` clean gives **2,928 records (34.3%)**.
- Requiring all four to carry an actual **literature citation** (tier 4 only): **382 records (4.5%)**.

Composition of the 1,434-record clean subset: median MW **137** (i.e. *smaller* than the full set's
174), MW ≥ 250 only 23.0%, ionisable **6.3%** (vs 15.1% overall), **20 drug-like**, and
**12 drug-like AND ionisable**.

> **The clean, drug-like, ionisable intersection is twelve compounds.** Filtering for provenance
> integrity does not merely shrink the data — it removes almost exactly the chemistry this effort
> exists to model. That is not a coincidence: drug-sized ionisable molecules are precisely the ones
> whose descriptors were hardest to measure and therefore most often estimated or predicted.

### 3b.6 The UFZ curators already knew

Corroboration from the UFZ side. The `dataset_and_outliers` sheet flags **877 of 7,241 rows (12%)**
as outliers, with a reason column. The reason breakdown is itself the finding:

| reason | n |
|---|---|
| **`ion`** (compound is ionic/ionisable) | **545** |
| **`only few partition coefficients used for calculations`** | **193** |
| `other` | 121 |
| `no organic compound` | 70 |
| `inconsistency in dataset` | 48 |

Two things follow. First, **the single largest exclusion category is "ion"** — the curators
deliberately removed 545 ionisable compounds from the modelling set. The scarcity of ionisable
chemistry in §4.2 is therefore *enforced*, not accidental, and it is enforced again downstream of
where we would look. Second, **193 records are flagged specifically because too few partition
coefficients were available to determine the descriptors** — an explicit, curator-supplied
low-confidence marker on the determination itself.

The same sheet carries **JChem-predicted `pKa,1–3` and `pKb,1–3`** for 5,229 of 7,241 compounds
(72%). Using those instead of my SMARTS heuristic, **1,574 compounds (21.7%) have an acid or base
pK between 2 and 8** — somewhat higher than the 15.1% SMARTS estimate, so treat 15–22% as the
range for "ionisable in the RPLC window". Note those pK values are themselves predictions
(ChemAxon), not measurements.

### 3b.7 What this means for validating our QSPR layer

1. **Re-baseline what "good" means.** Matching the literature's `E` R² = 0.98 is not evidence of
   chemical accuracy. If our QSPR hits that number it has learned ACD Absolv plus Abraham's
   analogy rules. **For an effort whose stated quality target is "close to ACD/Labs", distilling
   ACD/Labs and then reporting it as accuracy would be circular.** Say so in the spec.
2. **Report every QSPR metric stratified by provenance tier**, not just aggregate. Minimum:
   (a) all data, (b) tier-3+4 only, (c) the drug-like ionisable slice. If (b) is much worse than
   (a), the model is fitting the estimator, not the chemistry.
3. **A clean held-out measured test set can be carved out, but it is small and it is the wrong
   chemistry.** 1,434 compounds clean on all four descriptors; 12 of them drug-like ionisable.
   Held-out testing on the 1,434 is worth doing and is honest; it will not tell us what happens to
   ciprofloxacin.
4. **Train on everything, test on the clean subset.** Excluding the contaminated 83% from *training*
   would throw away the only coverage of large molecules. The defensible design is: train on all,
   report headline error on tier-3+4 held out, and treat tier-1/2 records as carrying inflated
   weight-of-evidence uncertainty.
5. **Therefore the QSPR layer cannot be validated to a trustworthy accuracy figure from descriptor
   data alone.** This strengthens the case, already implied by the map, for benchmarking the QSPR
   layer on **downstream retention error** — which is measured, plentiful, and uncontaminated —
   rather than on descriptor RMSE. Descriptor RMSE should be demoted to a comparability metric,
   exactly as `log k` RMSE already was.
6. **Chase Poole's Wayne State set (§1.5) with real urgency.** Several hundred solutes, one
   laboratory, standardised chromatographic determination, no Absolv involvement. It is the only
   candidate for a genuinely independent measured reference, and its value has gone up sharply.

---

## 4. Chemical space — the applicability-domain answer

### Method

Downloaded `Dataset_and_predictions.xlsx` (UFZ/Absolv, n=7,242 rows, 7,159 RDKit-parseable) and
`SoluteDB_all.xlsx` (n=8,481, all parseable). Computed with RDKit 2026.03.5: MW, ring/aromatic-ring
counts, rotatable bonds, HBD/HBA, cLogP, TPSA, Bemis–Murcko scaffolds, and SMARTS hits for ionisable
groups. "Ionisable at pH 2–8" = carboxylic acid ∪ sulfonic ∪ phosphonic/phosphoric ∪ tetrazole ∪
aliphatic primary/secondary/tertiary amine ∪ amidine/guanidine. Phenols and anilines/aromatic-N are
counted separately as *weak* acids/bases since they are largely un-ionised across the usable RPLC
pH range. Analysis scripts are throwaway; the numbers below are reproducible from the two public
files named above.

The two datasets give **the same picture**, which is expected once you know they are the same data:
overlap by InChIKey is **6,615 compounds = 92.4% of the UFZ set and 78.0% of SoluteDB**. SoluteDB
adds 1,865 not in the UFZ file; the UFZ file has 544 not in SoluteDB.

### 4.1 Size — small solutes

UFZ/Absolv, n=7,159:

| statistic | MW (g/mol) |
|---|---|
| min | 2 (H₂/He are in there) |
| 5th pct | 86 |
| 25th pct | 132 |
| **median** | **174** |
| 75th pct | 248 |
| 90th pct | 327 |
| 95th pct | 387 |
| 99th pct | 505 |
| max | 1203 |

- **MW ≥ 250: 24.7%.  ≥ 300: 14.2%.  ≥ 350: 7.9%.  ≥ 400: 4.0%.  > 500: 1.1%.**
- Median heavy-atom count **12**; median rotatable bonds **2**.
- **31%** of the set is acyclic. **58%** have **zero** hydrogen-bond donors.
- **26%** contain neither N nor O (pure hydrocarbons and halogenated hydrocarbons).
- SoluteDB is statistically indistinguishable: median MW 174, MW ≥ 250 = 23.2%.

Marchetto et al.'s MW 80–400 window keeps **92.6%** of the set — it removes almost nothing, because
almost nothing is large. Typical small-molecule drug MW (~350–500) sits at the **90th–99th
percentile** of this data.

### 4.2 Ionisation — mostly absent

UFZ/Absolv, n=7,159:

| group | n | % |
|---|---|---|
| carboxylic acid | 502 | 7.0% |
| sulfonic acid | 6 | 0.1% |
| phosphonic / phosphoric acid | 30 | 0.4% |
| tetrazole | 5 | 0.1% |
| aliphatic 1°/2°/3° amine | 102 / 145 / 286 | 1.4 / 2.0 / 4.0% |
| amidine / guanidine | 52 | 0.7% |
| **ionisable pH 2–8 (strong acid ∪ strong base)** | **1,082** | **15.1%** |
| **zwitterionic (both)** | **26** | **0.4%** |
| phenol (weak acid, mostly neutral in RPLC range) | 656 | 9.2% |
| aniline / basic aromatic N (weak base) | 1,056 | 14.8% |
| permanently charged species | **0** | **0.0%** |

SoluteDB: 13.8% ionisable, 0.3% zwitterionic, 0 permanently charged.

Cross-check using the JChem-predicted `pKa/pKb` columns shipped in the UFZ file rather than my
SMARTS: **1,574 compounds (21.7%) have an acid or base pK between 2 and 8** (§3b.6). Use **15–22%**
as the range. And note from §3b.6 that the curators flagged a further **545 compounds as outliers
with the reason `ion`** — the shortage of ionisable chemistry is partly *enforced by curation*.

Two things follow:

1. **There are literally no charged species in the compilation.** This is by design, not accident.
   Ulrich et al. say so: *"ionic chemicals are not covered by the approach at all. Because Coulomb
   interactions strongly affect physicochemical properties such as solubilities and partition
   coefficients, the LSER approach is not suitable for describing such chemicals"*
   ([ACS Environ. Au](https://doi.org/10.1021/acsenvironau.6c00063), §3.1). The Green group's README
   independently confirms the same filter ("limited to **neutral** solute compounds").
2. The descriptors are therefore **descriptors of the neutral species**. Our ionisation layer must
   supply the charged-form behaviour; the descriptor layer will never learn it from this data.
   This is consistent with the map's architecture (pKa as a separate pluggable provider) but it
   pins down *why* that seam is not optional.

### 4.3 Drug-likeness — 3% of the set

Using a deliberately loose filter (MW 250–600 **and** ≥ 2 rings **and** ≥ 3 rotatable bonds **and**
≥ 3 HBA):

- **drug-like: 516 (7.2%)**
- **drug-like AND ionisable at pH 2–8: 228 (3.2%)** — all 228 have a complete `E,S,A,B,V` set.
- SoluteDB: 551 (6.5%) and 237 (2.8%).

They are real drugs, and they are the right ones — a name scan of the UFZ file finds ibuprofen,
naproxen, diclofenac, warfarin, propranolol, atenolol, metoprolol, amitriptyline, nortriptyline,
desipramine, imipramine, diphenhydramine, ciprofloxacin, lidocaine, procaine, verapamil, quinine,
morphine, codeine, nicotine, ranitidine, omeprazole, amlodipine, sertraline, fluoxetine,
haloperidol, sulfamethoxazole, trimethoprim, tetracycline, erythromycin, caffeine, paracetamol.
The 228-compound subset includes diltiazem, risperidone, quetiapine, trazodone, ziprasidone,
minocycline, nafcillin, dicloxacillin, ergotamine, imidacloprid, thiamethoxam.

So the honest statement is **not** "there is no drug chemistry" — it is **"there are about 230
drug-like ionisable molecules, spanning 170 Murcko scaffolds, embedded in a 7,000-molecule set of
small neutral solutes."** That is enough for a held-out *validation* set, and not remotely enough
to train on.

### 4.4 Scaffold diversity — low

- 68.9% of the UFZ set has a ring scaffold at all; among those there are only **965 unique Bemis–Murcko
  scaffolds**.
- **Benzene alone is the scaffold of 1,849 molecules — 26% of the entire dataset.** Biphenyl 252,
  naphthalene 171, pyridine 134, quinoline 73, cyclohexane 68.
- Fused/heteroaromatic drug scaffolds are present only in tens (barbiturate ring 49, flavone 47,
  pyrimidinedione 44, dibenzothiophene 44).

This is a set built for environmental fate modelling: substituted benzenes, PAHs, chlorinated
solvents, pesticides, alkanes.

### 4.5 The applicability-domain gap, quantified

The clean way to say it. For the 228 drug-like ionisable compounds vs the whole UFZ set:

| descriptor | whole-set median | target-subset median | percentile of the subset median within the whole set |
|---|---|---|---|
| `S` | 1.03 | **2.35** | **93rd** |
| `B` | 0.52 | **1.96** | **96th** |
| `V` | 1.35 | **2.75** | **94th** |

**The median compound we care about sits at roughly the 95th percentile of the data we have.**
Half of the target chemistry is in the top ~5% tail — precisely where data density is lowest and
where any regression is extrapolating.

Compounding this, the `A` descriptor is **56.2% exactly zero** in the UFZ file (Ulrich et al. report
60% in their curation, and note it is "a significant value imbalance"). A model can score well on
`A` by predicting zero; drug-like molecules are exactly the ones with non-zero `A`.

Ulrich et al. reach the same conclusion visually (UMAP against PubChemLite): *"most chemicals with
experimental descriptors are clustered in the center of the chemical space plot… As the experimental
descriptors mainly cover simple chemical structures, the approach is limited for more complex
structures"*, and note that *"compounds located in sparsely populated regions may correspond to more
extrapolative prediction scenarios with potentially increased uncertainty"*.

---

## 5. What the two attached papers use as ground truth

### 5.1 Marchetto et al., Anal. Chem. 2025, 97, 6991–7001 (the blueprint)

<https://doi.org/10.1021/acs.analchem.4c03466> — CC-BY 4.0.

- **QSPR ground truth: the Abraham Absolv set from UFZ-LSER (7,881 → 6,437 → 6,401).** Yes, this
  matches §1.2. They restricted to **MW 80–400 g/mol**, which as shown above discards ~7% of the
  data and, more importantly, **caps the model at exactly the size range where drug chemistry
  begins**.
- They survey three sources and pick this one: UFZ-LSER (">7000 small molecules"), SoluteDB
  ("between 7000 and more than 8000 entries for each LSER solute parameter"), and Poole's Wayne
  State set ("several hundred solutes… single laboratory").
- **Reported QSPR test-set performance** (ridge regression on 313 alvaDesc molecular descriptors,
  5,120 train / 1,281 test): `E` RMSEP 0.10 (R² 0.98), `S` 0.23 (R² 0.87), `A` 0.11 (R² 0.89),
  `B` 0.11 (R² 0.95). **These are aggregate figures over the small-neutral-dominated set** — by §4.5
  they say little about the target chemistry. They also had to invent "MMAPE" to stop MAPE blowing
  up on the many `A`/`B` values that are zero, which is itself a symptom of §4.5.
- `V` is taken directly from structure, not predicted — consistent with §3.
- **The real bottleneck in that paper is not the QSPR layer.** Their LSER system parameters were
  fitted on a dataset of **48 solutes reduced to 36** (one column, Kinetex XB-C18, water/acetonitrile,
  φ = 0.2–0.7, single pH), and Table 1 of that set is all small neutrals — benzamide, caffeine,
  cresols, xylenes, anilines, nitrophenols, naphthalene. Their own caveat: "the data set used for
  LSER system parameter prediction contained only 35 (i.e., 36−1) solutes with similar size,
  polarity, and functional group compositions. Solutes that differ significantly in these properties
  are unlikely to allow accurate prediction."
- They explicitly treat retention factors as **"'true values', i.e., with no experimental uncertainty
  associated"** — the opposite of this effort's uncertainty-first stance. Their end-to-end MAPE
  ladder (LSS 9.1% → +LSER 17.8% → +QSPR 24.6%) is therefore a *lower* bound.

### 5.2 Ulrich et al., ACS Environ. Au (kNN / GCA / GNN)

<https://doi.org/10.1021/acsenvironau.6c00063> — CC-BY 4.0. Received Feb 2026, accepted May 2026.

- **QSPR ground truth: the same Absolv data, via LSERD**, curated to **6,364 chemicals**
  (`S` 6,364 · `A` 6,358 · `E` 6,347 · `B` 6,111 · `L` 5,572). Matches §1.2. 70/20/10 train/val/test.
- Independent test-set RMSE (consensus): `A` 0.075–0.126, `B` 0.095–0.151, `S` 0.169–0.234,
  `E` 0.086–0.188, `L` 0.246–0.454, with GNN best across the board.
- They predict **`E, S, A, B, L`** — not `V`. Confirms `V` is computed. Note they use `L`
  (hexadecane–air) because the target application is solvent–air systems; **the RPLC LSER in
  eq (2) of Marchetto uses `V`, not `L`**, so `L` is not needed by our architecture. That is a
  free reduction in scope for our QSPR layer: **four targets, not five.**
- They state the paper says "The up-to-date data set (May 2026) is provided in the corresponding
  GitHub repository <https://github.com/nadinulrich/solute_descriptor_prediction>". **As of
  2026-08-14 that repository's last push is 2021-11-07** and the only dataset in it is the older
  Ulrich & Ebert supplementary file. **The May 2026 curation is not actually there.** If we want it,
  we have to ask the authors.

**Both papers use the same ground truth. There is no independent second opinion in the literature
we have.**

---

## 6. Recommendation for the QSPR layer

1. **Primary ground truth: SoluteDB from [Zenodo 5792296](https://zenodo.org/records/5792296),
   CC-BY-4.0**, 8,481 solutes / 7,425 with complete `E,S,A,B`. It is the largest, the most
   permissively licensed, Abraham co-authored it — and decisively, **it is the only copy that
   carries per-descriptor provenance** (use `SoluteDB_logKw_logPow_all.xlsx`, which has the `Ref`
   columns, in preference to `SoluteDB_all.xlsx`, which does not). Use the UFZ/Absolv spreadsheet as
   a **cross-check and outlier-flag source** (877 curator-flagged rows with reasons, plus CAS/names
   and JChem pKa, which SoluteDB lacks), not as the redistributable primary.
1b. **Carry a provenance tier on every descriptor value through the whole pipeline.** Ingest the
   `E/S/A/B Ref` codes, classify them once (§3b.3), and keep the tier attached. Any accuracy figure
   the system or the spec reports must be stratifiable by it. Without this, every downstream number
   is uninterpretable.
2. **Compute `V`; never predict it.** McGowan volume, from structure. Four QSPR targets: `E, S, A, B`.
   Drop `L` — the RPLC LSER does not use it.
3. **Report error stratified twice — by chemistry and by provenance.** By chemistry: the drug-like
   ionisable subset (~228 compounds) separately from the aggregate. By provenance: tier-3+4 only
   versus all data (§3b.7). Any literature figure quoted at us (Marchetto's `E` RMSEP 0.10 /
   R² 0.98, `S` 0.23; Ulrich's `E` 0.086, `S` 0.169) must be re-measured under both stratifications
   before it enters the uncertainty budget. Treat the published `E` and `B` scores as unusable as
   accuracy claims until re-measured — they are the two most software-contaminated descriptors.
4. **Define the applicability domain in descriptor space, not just structure space.** The concrete
   test: a compound whose predicted `S > 2.0` or `B > 1.4` or `V > 2.4` is above the 90th percentile
   of the training data and should carry a widened uncertainty. This gives the prediction layer a
   cheap, defensible way to widen its retention distribution rather than silently extrapolating.
5. **Do not inherit a descriptor uncertainty from the data — there isn't one.** Calibrate it on
   downstream retention error, matching the pattern already committed for `pKa`.
6. **Licence position for the spec:** if the destination ever becomes commercial, UFZ-LSER's terms
   forbid it. Say so in the spec. SoluteDB/CC-BY is the path that survives commercialisation, with
   attribution.

---

## 7. What I could not determine

- Whether a **bulk export of UFZ-LSER experimental descriptors** is available on request, or only
  via the per-compound web UI. No download endpoint is exposed on the site.
- Whether the **May 2026 curated dataset** referenced by Ulrich et al. exists anywhere public. The
  named GitHub repo does not contain it.
- The **exact contents and licence of Poole's Wayne State dataset** (J. Chromatogr. A 2020, 1617,
  460841) — paywalled; I did not find an open copy. This is the most valuable *independent* test set
  candidate and is worth chasing.
- Whether any consolidated, openly-licensed **Acree solute-descriptor table** exists separate from
  Absolv/SoluteDB. I found only per-paper determinations and CC0 *solvent*-coefficient datasets.
- The **exact split of UFZ-LSER between peer-reviewed-literature records and ABSOLV-table records**.
  The help page says descriptors come from one or the other, but publishes no counts and offers no
  bulk export. The circulating "~5,000 literature vs 7,881 used" figure is unconfirmed from a
  primary source. §3b.3 reaches the same conclusion by auditing SoluteDB's reference codes instead,
  which I consider the stronger evidence.
- **What reference `R1` is** in SoluteDB's reference sheet. It is the basis for `E02`, which
  accounts for 47% of all `E` values, and it is not defined anywhere in the file. This single
  undefined code is the largest untraceable block in the entire corpus. Worth an email to
  whgreen@mit.edu / yunsie@mit.edu (contacts given in the dataset README).
- Whether ACD/Labs' commercial **Absolv** predictor's underlying experimental table is larger or
  cleaner than the public copies. Not investigated (commercial, no public spec).
- My provenance classification is **text-heuristic over 1,046 free-text reference strings** and
  will have some misassignment at the margins; the headline blocks (`E02` 3,875, `S03` 3,509,
  `B10` 1,657, `E19` 911, `A01` 4,531) are unambiguous and dominate the totals, so the conclusion
  is robust even if individual small codes are misfiled. The `L` column could not be classified
  (75% of its reference codes have no text in the sheet) — irrelevant here, since the RPLC LSER
  uses `V`, not `L`.
