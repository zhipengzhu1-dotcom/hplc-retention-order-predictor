# MOCCA — Haas et al., *ACS Central Science* 2023

Extracted notes. The PDFs are gitignored (copyright); everything needed is reproduced here with
page/section references so the PDFs need not be reopened.

## Citation, DOI, licence

| Field | Value |
|---|---|
| Title | Open-Source Chromatographic Data Analysis for Reaction Optimization and Screening |
| Authors | Christian P. Haas, Maximilian Lübbesmeyer, Edward H. Jin, Matthew A. McDonald, Brent A. Koscher, Nicolas Guimond, Laura Di Rocco, Henning Kayser, Samuel Leweke, Sebastian Niedenführ, Rachel Nicholls, Emily Greeves, David M. Barber, Julius Hillenbrand,* Giulio Volpin,* Klavs F. Jensen* |
| Journal | ACS Cent. Sci. 2023, 9 (2), 307–317 |
| DOI | 10.1021/acscentsci.2c01042 |
| Received / Published | Sept 5, 2022 / Feb 9, 2023 |
| Licence | "© 2023 The Authors. Published by American Chemical Society" (p. 307 footer). Page 1 of the PDF carries hyperlinks to `https://creativecommons.org/licenses/by/4.0/` and `https://acsopenscience.org/open-access/licensing-options/` — i.e. **CC-BY 4.0**. The licence *badge is an image*; the words "Creative Commons" appear nowhere in the extractable text layer, so the CC-BY 4.0 identification rests on that embedded link. |
| Code | https://github.com/HaasCP/mocca ; PyPI `mocca` ; docs mocca.readthedocs.io (p. 314, Data Availability) |
| Data | Simulated HPLC–DAD data sets saved as MOCCA campaigns: **10.5281/zenodo.7406829** (p. 314). SI PDF (`oc2c01042_si_001.pdf`, 69 pp) + a ZIP of MOCCA HTML reports for the well-plate screening. |
| Affiliations | MIT Dept. Chemical Engineering + Bayer AG (Crop Science, Pharmaceuticals, Enabling Functions divisions) |

MOCCA = **M**ultivariate **O**nline **C**ontextual **C**hromatographic **A**nalysis (p. 308).

> **Version caveat, important for us.** This paper describes MOCCA v1 (the `HaasCP/mocca`
> repository). The **MOCCA2** rewrite that we have read the source of is a *different* codebase
> with *different* matching predicates and defaults. Numbers below are what the *paper* states;
> where MOCCA2's code differs, that is flagged in §2. Nothing in this README should be treated as
> documentation of MOCCA2.

---

## 1. The algorithm

Order below follows SI section S4 ("Technical description of implemented data analysis features",
SI pp. S9–S22), which is the only complete algorithmic description; the main text (Fig. 2, p. 309)
is a one-slide summary.

### 1.0 Vocabulary the paper uses (SI p. S9)

- **compound run** — a run whose contents the user declares (e.g. a pure analyte at a known
  concentration). Used to "train" the tool.
- **analysis run** — a run of unknown composition.
- **gradient run** — a blank injection, used for baseline correction.
- **campaign** — a set of HPLC–DAD data sets analysed *in context with each other*. All cross-run
  logic (compound library, tracking, calibration) operates within one campaign.

### 1.1 Preprocessing (SI p. S9–S10)

- *Bandwidth*: rolling average over the wavelength dimension, **default 2 nm**.
- *Reference*: the mean of the **highest five wavelengths** is subtracted from every wavelength.
  Caveat stated: the DAD upper wavelength limit must be set >10 nm above any analyte absorbance
  band, otherwise this produces negative signals.

### 1.2 Baseline correction (SI p. S10)

- Requires a **user-supplied blank gradient run** for each data set. Not optional.
- Per wavelength, a **modified asymmetric least squares smoothing** baseline algorithm (refs 15,16
  of the SI) is applied to the blank run; this "levels out minor carryover or impurity signals as
  well as noise and other artifacts like injection peaks". The smoothed blank array is then
  **subtracted** from the compound/analysis runs.
- Practice mandated: two blank gradient runs before real samples (the first after startup is "often
  disturbed"); repeat a blank periodically during long sequences to track drift.

### 1.3 Peak picking, merging and expansion (SI pp. S11–S12)

- Absorbance array is collapsed to a **summed-over-all-wavelengths** chromatogram.
- User sets an **absorbance threshold**; all summed values below it are set to zero, so a peak apex
  must exceed the threshold to be picked. `scipy.signal` does the picking and the peak-width
  calculation.
- **Overlapping peaks are merged into one peak model** — but only if the *minimum* summed
  absorbance in the overlap region exceeds the absorbance threshold. (So the threshold setting
  directly controls whether two touching peaks become one object.)
- High-pass / low-pass filters in the time dimension exclude injection and re-equilibration
  artifacts.
- **Peak expansion**: rolling average window 5 applied to the (unclipped) summed chromatogram;
  borders extended until the smoothed signal falls below **1/20 of the absorbance threshold** or
  starts increasing.

### 1.4 Peak purity check (SI pp. S12–S13) — a 4-stage cascade

Basis: correlation coefficients between the UV–Vis spectrum at every timepoint inside the peak
borders and the spectrum at the **apex**. Tests run in order; the first that concludes wins.

| # | Test | Exact criterion (SI p. S12–S13) |
|---|---|---|
| 1 | Unimodality | Rolling average (size 3) of the correlation-coefficient vector must be unimodal (rise then fall); coefficients **> 0.999 are ignored**. Fail ⇒ **impure**, stop. |
| 2 | Agilent-style threshold | `thresh_t = ( max(0, 1 − 2.5 · Var_noise · (1/Var_t + 1/Var_apex)) )²`, where `Var_noise` is the spectral variance where no compound elutes, `Var_t` at timepoint t, `Var_apex` at apex. If **≥ 90%** of timepoints have correlation above their threshold ⇒ **pure**, stop. |
| 3 | PCA | If variance explained by the **first principal component > 0.995** ⇒ **pure**, stop. |
| 4 | Raw coefficients | Any coefficient **< 0.9** ⇒ **impure**. Else if all **> 0.95**, or mean **> 0.98** ⇒ **pure**. |
| — | Fallthrough | If no test concludes ⇒ **impure**. |

The authors state explicitly that they deliberately **loosened** the Agilent test: *"This test
differs from the original Agilent peak purity checker, in that we use the coefficient 2.5 instead
of 0.5, and that the test passes with 90% of good points rather than 100%, making the test more
permissive, as the Agilent peak purity checker gave many false negatives in our testing."* (SI p.
S13). This permissiveness is the direct cause of their dominant failure mode — see §6.

### 1.5 Integration (SI p. S14)

Integral = sum of absorbance over **all** recorded wavelengths and all timepoints inside the peak
borders, after subtracting the lowest value in the peak's absorbance array (a horizontal baseline).
Rationale given: analytes cannot "hide" at a single monitoring wavelength.

### 1.6 Deconvolution of overlapping peaks (SI pp. S14–S16)

Math: **PARAFAC (parallel factor analysis)**, specifically `non_negative_parafac_hals` from
**TensorLy** (non-negativity constrained ALS-HALS). Not NNLS, not MCR-ALS.

The paper is candid about the theoretical problem: an HPLC–DAD array is bilinear; stacking runs
gives a tensor that is *ideally* trilinear, but *"the retention time dimension is considered as a
trilinearity-breaking mode since elution profiles are not perfectly reproducible over multiple runs
(in contrast to the spectral mode where UV-Vis spectra of given analytes will not change between
runs). In theory, this disables the use of classical deconvolution algorithms for trilinear data
like PARAFAC"* (SI p. S14). Their answer is not a better algorithm but a **construction of the
tensor**:

1. **Trigger condition** (both must hold): the peak **failed the purity check** *and* the peak
   **overlaps a known compound model in the library**. Impure peaks not touching a library
   compound are not deconvoluted (SI p. S15).
2. Estimate component count by **PCA**: minimum number of principal components to exceed a
   cumulative-explained-variance threshold, **default 0.995** (SI p. S16).
3. Build a 3-slice tensor: slice 1 and 2 are absorbance arrays of **pure compound peaks** taken
   from the library (aligned to the library compound's maximum, and **normalised to half the peak
   maximum height of the impure peak**); slice 3 is the impure peak (SI pp. S15–S16).
4. **Iterative shift**: assume peak *shape* is unchanged across runs and only *location* varies.
   The impure slice is shifted **timepoint by timepoint from left to right**; PARAFAC is run on
   each shifted tensor; the modelled peaks are summed and compared to the raw impure peak by
   point-by-point **MSE**; the **lowest-MSE iteration wins** (SI p. S16). *"With this iterative
   PARAFAC routine, the peak deconvolution algorithm is able to adapt for offsets in the retention
   time over runs."*

Note the scope of that last sentence: the shift search is **inside a single impure peak**, to align
a library reference against a slightly drifted peak. It is not a cross-run retention alignment
mechanism. See §3.

### 1.7 Compound models / library (SI pp. S17–S18)

- Built from **compound runs**. The user supplies a compound identifier string, a concentration, and
  optionally internal standards. MOCCA **assigns the highest non-internal-standard peak** in that
  chromatogram to the given identifier; **all other peaks are labelled as impurities** attached to
  that identifier.
- A compound model's **main attributes are the UV–Vis spectrum and the retention time**.
  - Spectrum = **average over all compound runs** of that compound. Justification given: *"Since
    UV-Vis spectra do not vary over runs, there is no need to 'dilute' this information with UV-Vis
    spectra of later analysis runs"*.
  - Retention time = **averaged over all occurrences in the campaign**, *"With that, minor
    systematic retention time drifts over multiple runs can be accounted for."* Optional filter,
    e.g. only the last five occurrences.
- **Unknowns**: pure peaks in analysis runs with no match become `unknown_N` with an iterative
  counter, and are themselves added to the library — so unknowns become trackable across runs
  without any standard being injected.

### 1.8 Calibration (SI p. S18) and internal standards (SI pp. S18–S19)

- Calibration wavelength = the **absorbance maximum of that compound's model spectrum**. Integral
  vs concentration fitted with a **through-origin linear model `y = mx`**; slope `m` and R² stored.
- Internal standard does two jobs:
  1. **Retention-time correction.** The IS signal is located *before* peak assignment and
     deconvolution, considering **only pure peaks with matching UV–Vis spectrum and retention
     time**. (If an impure peak sits in the IS retention region, it is deconvoluted first using the
     pure IS signal to build the tensor.) The difference between the found IS apex and the library
     IS retention time is stored as an **offset applied to every peak object in that
     chromatogram**, and assignment then uses the corrected retention times.
  2. **Relative quantification**: `c_A = m · (I_A · c_ISTD / I_ISTD)`. When an IS is present MOCCA
     builds both absolute and relative models and **uses the relative one by default**.

### 1.9 Peak assignment — the cross-run step (SI pp. S19–S20)

Two steps.

**Step 1 — candidate generation.** For each **pure** peak, candidate library compounds are found
using **two user-set thresholds, ANDed**:

- **Relative retention-time threshold** — expressed as a fraction of the total number of timepoints
  in the chromatogram. Worked example given verbatim (SI p. S19): a 10 min (600 s) method with a
  threshold of **0.01** means a peak at 300 s matches a library compound only if the library
  retention time is in **294–306 s**, i.e. **±0.01 · 600 s = ±6 s**.
- **UV–Vis spectrum correlation threshold** — correlation coefficient between the peak's
  **average spectrum over all its timepoints** and the library spectrum must **exceed** the
  threshold.
- *"Matches are assigned to the peak only if both conditions are satisfied."*

**Step 2 — greedy conflict resolution.** Across all peaks and all their candidate matches in the
chromatogram, the pair with the **highest UV–Vis correlation coefficient** is assigned first; that
compound identifier is then deleted from every other peak's candidate list; repeat until exhausted.
Leftover pure peaks become unknowns/impurities and are added to the library.

**Compound tracking** (SI p. S20) is then trivial bookkeeping: because every peak in every
chromatogram of the campaign carries a library identifier, integrals and concentrations can be
plotted per compound across runs (`compound_tracking.html`). Deconvolution-derived peaks are pushed
back onto the chromatogram's peak list and are therefore tracked too.

### 1.10 Reporting (SI pp. S21–S22)

Nine stand-alone HTML reports via Datapane: `hplc_input`, `gradient`, `chromatograms`,
`bad_chromatograms`, `deconvolution`, `peak_library`, `compound_library`, `calibration_library`,
`compound_tracking`. Note `bad_chromatograms` — chromatograms MOCCA refuses to analyse because user
input and data disagree (e.g. a declared compound run whose highest peak is impure).

---

## 2. How MOCCA decides two peaks in different runs are the same compound

**Exact criteria (SI p. S19–S20), restated:** a peak matches a library compound iff

```
|t_peak − t_library| ≤ rel_rt_threshold × (total run length)      AND
corr(mean spectrum of peak, library spectrum) > spectrum_threshold
```

then greedy assignment by descending spectral correlation resolves duplicates.

| Property | What the paper says |
|---|---|
| Does matching depend on retention time? | **Yes, hard gate, ANDed.** No spectrum-only path exists in the described algorithm. |
| How tight is the RT gate? | Expressed as a **fraction of run length**, not of peak width. The only number ever given is the illustrative **0.01 → ±6 s on a 600 s method** (SI p. S19). |
| Stated default for the RT threshold | **None. The paper never gives a default value** — 0.01 is presented as "a given relative retention time threshold" in a worked example. |
| Stated default for the spectral correlation threshold | **None. The paper never gives a number for it at all** — neither a default nor the value used in any of the four case studies. |
| Any tolerance widening for known retention shifts | **None described.** |

Contrast with the MOCCA2 code we have read (`are_same_compound` in `src/mocca2/dataset/dataset.py`):
MOCCA2 gates on `|t1 − t2| ≤ (w1 + w2) · max_peak_distance` with `max_peak_distance = 1.0`, i.e.
**peak-width-relative**, and `min_spectrum_correl = 0.99`. **Both the parameterisation and the
numbers differ from the paper.** The paper is not a source for MOCCA2 defaults; the paper's
threshold is run-length-relative and unspecified.

---

## 3. Does the paper handle retention that SHIFTS between runs?

This is our key question, so the silences are itemised.

**What the paper explicitly supports:**

| Mechanism | Scope, per the paper |
|---|---|
| Internal-standard retention-time correction (SI p. S18) | A **single scalar offset per chromatogram**, derived from one IS peak, added to all peaks in that chromatogram. It corrects *bulk* shift, not compound-specific shift, and certainly not order reversal. |
| Library retention time = running average over the campaign, optionally last-5 (SI p. S17) | Explicitly scoped to *"minor systematic retention time drifts over multiple runs"*. |
| Iterative-shift PARAFAC (SI p. S16) | *"able to adapt for offsets in the retention time over runs"* — but the shift search happens **within one impure peak** while aligning a library reference slice, and is only triggered for impure peaks overlapping a known compound. It does not feed peak assignment. |

**What the paper explicitly forbids (SI p. S7, "HPLC best practices using MOCCA"):**

> *"In general, MOCCA is only capable of analyzing datasets in context to each other, if all
> datasets have been recorded on the same HPLC system, with the same HPLC column, using the same
> HPLC method. It cannot account for systematic drifts induced by a change of the setup, such as
> changes in the system's dead volume, the change to another HPLC column, or the change of the
> HPLC gradient. If not done properly, even the exchange of the mobile phase can lead to such
> disturbing drifts."*

Our planned sweep changes the mobile phase pH, the mobile phase composition, and the column. All
three are named in that sentence as things MOCCA cannot account for.

**Silences — stated as findings:**

- The paper **never** discusses a compound whose retention factor changes *by design* between runs.
  There is no ionisable-analyte, no pH, no buffer, and no pKa content anywhere in the paper or SI.
  The strings "pH", "buffer", "isocratic", "dead time", "t0", "dwell volume" do not appear in any
  methodological sense; the single "pH" hit in the SI (p. S43) is a cyanide-safety instruction.
- It **never** discusses **elution-order reversal** or any mechanism that would survive it.
- It **never** discusses matching peaks across *different methods*, *different columns* or
  *different mobile phases* — it forbids it (above).
- It **never** reports an experiment where retention deliberately changes and the tracker is
  scored on it. The only cross-run variation exercised is the Knoevenagel gradient-length sweep
  (0.5/0.75/1.0/1.5/2.5 min, p. 310), and those are analysed as **separate campaigns per method**,
  each with its own calibration (SI p. S26–S27), not tracked across each other.

**How essential are references/internal standards?** The SI's summary rule (p. S8) is unambiguous:

> *"All HPLC sequences should be preceded by two blank injection runs (or gradient runs without
> injections). Moreover, pure standard HPLC runs must be performed for all compounds which should
> be assigned and tracked and for which the peak deconvolution feature should be activated."*

Internal standards are **optional** (retention correction + relative quantification); blank
gradient runs are **mandatory** (baseline correction); pure standard runs are **mandatory for any
named compound**.

---

## 4. How discriminating is the UV/DAD spectrum, on the paper's own data?

The paper's only quantitative probe of spectral discrimination is the **simulation study** (main
text §3.1 p. 309–310, Table 1; SI §S7 pp. S28–S36).

**Design (SI pp. S28–S30, Fig. S6).** CADET (transport-dispersive model + nonlinear steric mass
action isotherm, ion-exchange gradient with linearly increasing salt) simulated **1000** retention
profile triplets: pure main compound, pure impurity, and the two-compound overlap. Isotherm
parameters were randomly sampled; operating conditions fixed. Each triplet was then augmented with
**real measured UV–Vis spectra** at three similarity levels → **3000 campaigns**.

| Similarity level | Main compound spectrum | Impurity spectrum | Correlation r (p. 309) |
|---|---|---|---|
| High | 4-chlorobenzaldehyde | benzaldehyde | **r ≈ 0.86** |
| Medium | 4-chlorobenzaldehyde | 4-methoxybenzaldehyde | **r ≈ 0.47** |
| Low | 4-chlorobenzaldehyde | 4-(dimethylamino)benzaldehyde | **r ≈ −0.06** |

**Table 1 (p. 309), verbatim counts out of 1000 per column:**

| Result category | High (r≈0.86) | Medium (r≈0.47) | Low (r≈−0.06) |
|---|---|---|---|
| (i) baseline-separated, both peaks correctly called pure and analysed | 86 | 86 | 86 |
| (ii) overlapping, deconvolution triggered, main compound identified & quantified | 794 | 868 | 890 |
| (iii) deconvolution triggered but main compound could **not** be identified | **2** | 0 | 0 |
| (iv) signals overlapped but the purity checker **failed to flag them as impure** | **118** | 46 | 24 |

(i) and (ii) are the desired outcomes; (iii) and (iv) are called "misinterpretations" (p. 310).

**Quantification error for category (ii)** — SI Fig. S9 caption (p. S34), main text p. 310:

| Similarity | n | Median error | 3rd-quartile error |
|---|---|---|---|
| High | 794 | **1.6%** | **4.4%** |
| Medium | 868 | **1.5%** | **6.2%** |
| Low | 890 | **1.6%** | **5.3%** |

Main text (p. 310): *"For all three levels of spectral similarity, the median quantification error
was smaller than 2%, while the third quartile error ranged around 6%."*

**What this does and does not establish for us.**

- The spectral-similarity axis is **the pair (main compound, impurity) inside one peak**. It is a
  test of *deconvolution* under spectral similarity, **not** a test of *cross-run identity
  matching* under spectral similarity. No experiment in the paper puts two similar-spectrum
  compounds in a compound library and asks whether MOCCA confuses them across runs.
- Consequently the paper reports **no false-match rate, no misassignment rate, and no
  confusion-matrix** for compound assignment. Category (iii) (2 cases out of 3000) is
  "deconvolution could not resolve identity", not "assigned to the wrong compound".
- Where the paper does address analogue confusability it is a **prose assertion, not a
  measurement**: *"UV-Vis spectra do not vary over runs"* (SI p. S17), used to justify averaging
  spectra only from compound runs. And SI p. S32: *"if the peak deconvolution feature is triggered,
  the qualitative analysis (finding the correct UV-Vis trace for peak assignment) is highly
  reliable"* — supported by the 2/3000 category (iii) count, and only within the deconvolution
  path.
- **The clearest evidence against spectra being a sufficient fingerprint is Table 1's category
  (iv) column.** Highly similar spectra (r ≈ 0.86 — benzaldehyde vs 4-chlorobenzaldehyde, exactly
  the "substituted aromatics" case) produced **118/1000 undetected coelutions**, versus **24/1000**
  at r ≈ −0.06 — a **~5× increase in silent failures** when spectra resemble each other. The
  authors attribute this to the purity checker, not to the assignment step, but the mechanism is
  the same one our tracker would rely on: spectral contrast.
- SI p. S34 adds a related observation: bad quantification correlates with MOCCA failing to recover
  the *impurity's* spectrum, and *"These cases appear more often in the datasets of the medium and
  high spectral similarity levels."*

Real-data corroboration is anecdotal, not quantified: the four benzaldehyde-family compounds in the
Knoevenagel study were labelled correctly in all instances (SI p. S27, "in all instances (two- and
three-substrate cases)" — no count given), and one unknown in the well-plate screen was tracked
consistently as `unknown_3` across runs by "similarity in retention time and UV-Vis-spectrum" and
later hypothesised to be butyl benzoate (SI p. S49, Fig. S17).

---

## 5. Does it need pure standards?

| Capability | Standards required? |
|---|---|
| Naming a compound (assigning a user-chosen identifier) | **Yes.** A compound run of the pure substance is required — SI p. S8: *"pure standard HPLC runs must be performed for all compounds which should be assigned and tracked"*. |
| Quantification (concentration) | **Yes**, compound runs with concentration information; a linear `y=mx` calibration is built (SI p. S18). |
| Triggering deconvolution on a peak | **Yes** — the trigger requires overlap with a *known compound model in the library* (SI p. S15). |
| Tracking an **unknown** across runs | **No.** Unassigned pure peaks become `unknown_N` and are added to the library, then tracked by the same RT+spectrum criteria (SI p. S17, S49). This is the only standard-free path, and it still uses the retention-time gate. |
| Baseline correction | Requires **blank gradient runs**, not standards (SI p. S10). |

So: it can operate without standards only in "unknown_N" mode, which gives tracking but no
identity and no concentration.

---

## 6. The paper's own stated limitations and failure modes (quoted)

1. **Same system, column, method only** — SI p. S7: *"MOCCA is only capable of analyzing datasets
   in context to each other, if all datasets have been recorded on the same HPLC system, with the
   same HPLC column, using the same HPLC method. It cannot account for systematic drifts induced by
   a change of the setup, such as changes in the system's dead volume, the change to another HPLC
   column, or the change of the HPLC gradient. If not done properly, even the exchange of the
   mobile phase can lead to such disturbing drifts."*
2. **The user must still do method development** — SI p. S7: *"MOCCA is not designed to solve all
   analysis tasks and problems on the side of an automated data analysis tool. Instead, the user is
   required to have adequate expertise in HPLC techniques."* And SI p. S8: *"the peak deconvolution
   feature for overlapping peaks is designed to cover the case when a known analyte is overlapped
   by an unknown and unexpected impurity. HPLC method development should be performed to a point
   that all known analytes are baseline-separated."* — i.e. **coelution of two *known* compounds
   is out of scope by design.**
3. **Undetected coelution is the dominant failure mode, and it is analytically irreducible** —
   main text p. 310: *"Cases of the category (iv) are not attributed to a failure of the peak
   deconvolution feature, but rather to a permissive peak purity checker returning false positive
   outcomes on strongly coeluting signals. These cases cannot be solved analytically, and the only
   solution would be the development of an HPLC method with higher chromatographic resolution to
   separate (at least partially) the elution profiles. MOCCA enables a shift toward shorter
   gradient times and faster sample processing, but the category (iv) failure rate shows that the
   user is still required to have expertise in HPLC method development to balance method time vs
   chromatographic resolution."*
4. **Deconvolution accuracy is screening-grade, not regulated-grade** — p. 310: *"The results
   obtained validate that MOCCA's deconvolution feature works robustly enough for typical lab
   screenings, but should be treated with caution for regulated environments and process
   development scenarios where lower margins of error are required."*
5. **PARAFAC is theoretically invalid on this data and the fix is empirical** — SI p. S14–S15:
   the retention dimension is *"a trilinearity-breaking mode"*, which *"in theory... disables the
   use of classical deconvolution algorithms for trilinear data like PARAFAC"*; the workaround is
   *"a novel and pragmatic approach... tackling the problem from the side of how to build the data
   tensor"*, with the tensor built by an *"empirically developed routine"*, resting on the
   assumption *"that the peak shape does not significantly change over the runs and only the peak
   location varies from slice to slice."*
6. **Detector-range constraint** — SI p. S9: the DAD upper wavelength must exceed any analyte band
   by >10 nm or the reference-subtraction produces negative signals.
7. **Performance** — SI p. S49: *"The package was not optimized for performance yet, especially,
   the iterative PARAFAC algorithm requires significant computation time."* 96-well data set with
   42 deconvolutions: **7:50 min** (Intel Xeon W-11955M, 64 GB) or **10:57 min** (2.3 GHz i9,
   16 GB).
8. **Vendor export friction** — SI p. S5: no automated raw-data export path was found for Waters
   Empower (manual `.arw` only); only Agilent could export true DAD raw data in Allotrope ADF, and
   other vendors' ADF exports are *"currently feature-limited to single wavelength channel
   exports."*

---

## 7. Data and licence details, and what metadata is (not) recorded

**Released data.** Simulated MOCCA campaigns at **Zenodo 10.5281/zenodo.7406829** (p. 314),
described as "Simulated HPLC–DAD data sets saved as MOCCA campaigns, which were used for validation
and benchmarking" — this is the **3000-campaign CADET set** of §4, i.e. **synthetic**, generated
from an **ion-exchange** simulation with a **linear salt gradient** (SI p. S28), then augmented
with real benzaldehyde-family spectra. Also released: the well-plate screening `.txt` raw data +
Jupyter notebook in the GitHub `notebooks` folder (SI p. S49), and MOCCA HTML reports as SI ZIP
(`oc2c01042_si_002.zip`).

**Instruments and methods actually used.**

| Study | Instrument | Column | Mobile phase | Detection |
|---|---|---|---|---|
| Knoevenagel kinetics (SI p. S26) | Agilent LC: G1379B degasser, G1312B binary pump, G1329B autosampler + G1330B thermostat, G1316B column oven, **G1315B DAD** | **Kinetex 50 × 3 mm, 2.6 µm core-shell C18 (TMS endcapped), p/n 00B-4462-Y0**, 25 °C | Water / MeCN, both **+0.05% v/v TFA**; **1.5 mL/min**; 0.1 min hold at 95:5, then linear gradient to 0:100 over **0.5 / 0.75 / 1.0 / 1.5 / 2.5 min**, 0.1 min hold at 0:100, 0.1 min back to 95:5 | 200–550 nm, **1 nm step**, every timepoint; single-λ comparison at 248 / 283 / 347 nm |
| Closed-loop 2-pyridone (SI p. S39) | Same Agilent system, autosampler bypassed, VICI C84H-1574-.02EUHF internal valve, **20 nL** injection | as above | as above (not restated) | as above |
| Well-plate cyanation (SI p. S43) | **Shimadzu Nexera UHPLC, SPD-40M UV-Vis detector**, flow cell at 40 °C | **Waters Acquity 2.1 × 50 mm, 1.7 µm BEH C18** | **not stated** | **190–400 nm, 1.3 nm resolution** |

**Metadata explicitly NOT recorded anywhere in the paper or SI:**

| Item | Status |
|---|---|
| Isocratic data | **None.** Every described method is a gradient (RP water/MeCN gradients; the simulation is an ion-exchange salt gradient). The word "isocratic" does not appear. |
| Dead time t₀ / column void volume | **Never reported or used.** "Dead volume" appears once (SI p. S7) only to warn that changing it breaks the analysis. |
| Dwell volume | **Never mentioned.** |
| Mobile phase pH | **Never mentioned.** TFA 0.05% v/v is stated as an additive; no pH value is given, and there is no buffer anywhere. |
| Column identity | **Recorded in prose** (part numbers for Kinetex and Acquity BEH) but **not carried as data-model metadata** — MOCCA's campaign model has no column field described. |
| Retention factor k, hold-up time, or any thermodynamic quantity | **Absent.** MOCCA works entirely in raw retention time. |
| Allotrope ADF | Supported as a *format* with a parser, and pitched as the FAIR/metadata-rich option (p. 308–309, SI pp. S5–S6), but the paper targets only the **"data cube layer"** (raw arrays); it reads the description layer via h5ld/rdflib but no chromatographic-condition metadata is shown being used by the algorithm. |

**Parsers available** (SI p. S21): `chemstation`, `labsolutions`, `empower`, `allotrope`, `custom`
(arrays already in Python memory — this is how the CADET simulations were fed in).

---

## 8. Peak purity / coelution detection — usable for detecting merged peaks?

Relevant because for a pH sweep, **detecting** a merge matters more to us than deconvolving it.

**What exists.** The 4-stage purity cascade of §1.4 runs on **every picked peak**, independently of
whether a compound library exists, and its verdict is surfaced per peak (green/red in
`chromatograms.html`, and as a column in the peak table; SI pp. S16, S21–S22). Deconvolution is a
*separate*, gated stage; purity detection is not gated on library membership. So the detector is
usable standalone.

**Measured detection performance (Table 1, p. 309) — this is the number that matters:**

| Spectral similarity of the two coeluting species | Overlapping pairs missed by the purity checker (category iv) | Rate |
|---|---|---|
| High, r ≈ 0.86 | 118 / 1000 | **11.8%** |
| Medium, r ≈ 0.47 | 46 / 1000 | **4.6%** |
| Low, r ≈ −0.06 | 24 / 1000 | **2.4%** |

Caveats attaching to that table:
- Denominator is the whole 1000-campaign set; 86 of each 1000 were baseline-separated (category i),
  so the miss rate among *actually overlapping* pairs is ~118/914 ≈ 12.9% (high), ~24/914 ≈ 2.6%
  (low). The paper does not compute this normalisation itself.
- The paper reports **no false-positive rate for the purity checker on genuinely pure peaks** — no
  category in Table 1 counts pure peaks wrongly flagged impure. They state they loosened the test
  *because* Agilent's gave "many false negatives" in their testing (SI p. S13), which trades in
  the direction of more misses.
- Missed cases are described as the ones where signals *"are co-eluting almost perfectly"* (SI Fig.
  S8d caption, p. S33) — i.e. small Δt, where an orthogonal spectral contrast is the only signal
  left and, at r ≈ 0.86, there is little of it.
- The three-level design means the paper has evidence for exactly three spectral-similarity points
  and one (simulated, ion-exchange) peak-shape family. There is no curve of detection probability
  vs resolution Rs, and Rs is never reported.

**Secondary merge-related mechanism.** Peak *picking* merges two peaks into a single peak model
whenever the summed absorbance in the valley stays above the absorbance threshold (SI p. S11). This
is a separate, threshold-dependent way for two compounds to become one object — and unlike the
purity check, it has **no reported error rate at all**.

---

## 9. Claimed vs demonstrated — quick separation

| Claim | Status |
|---|---|
| Baseline correction, peak picking, integration, purity check, PARAFAC deconvolution implemented | **Demonstrated** (code released; four case studies) |
| Deconvolution recovers the known compound's concentration to median <2% / Q3 ~6% error | **Demonstrated on simulated data** (3000 campaigns, Table 1 + Fig. S9) |
| Quantification of pure, baseline-separated peaks matches manual ChemStation analysis | **Demonstrated** (SI §S5, Figs. S4–S5) — qualitative/visual: "almost perfect correlation", no R² or slope reported |
| Deconvolution keeps a closed-loop optimisation running through unexpected coelution | **Demonstrated** (2-pyridone campaign, Fig. 6b/6c; product 6 vs butylated DBU 8 at ~1.7 min), verified by independent batch reactions (Fig. S16) |
| Unknown side products can be tracked across runs without standards | **Demonstrated anecdotally** (`unknown_3` → butyl benzoate, Fig. S17); no accuracy figure |
| "UV-Vis spectra of given analytes will not change between runs" | **Asserted** (SI pp. S14, S17). Not tested; note it is asserted for *fixed* mobile phase — a pH-dependent spectrum is outside the paper's frame entirely |
| Qualitative assignment inside deconvolution is "highly reliable" | **Asserted**, supported only by category (iii) = 2/3000 |
| MOCCA is applicable across different columns / methods / mobile phases | **Explicitly disclaimed** (SI p. S7) |

---

## 10. Bottom line for our pH sweep

Stated as facts about the paper, not as predictions about our system:

- The described matching predicate is an **AND of a retention-time gate and a spectral gate**, with
  the RT tolerance expressed as a fraction of run length; the paper's only illustrative value is
  **±1% of run length**. There is no spectrum-only mode in the paper's algorithm.
- The paper's scope statement (SI p. S7) **explicitly excludes** analysing data in context across
  changed columns, changed mobile phases and changed methods — three of the axes our sweep varies.
- Retention change **by design** is nowhere addressed; the only retention-variation the paper
  handles is drift, via a single per-chromatogram internal-standard offset and a running average of
  library retention times.
- The paper contains **no measurement** of cross-run misassignment rate, and **no measurement** of
  whether structurally similar aromatics are confusable at the assignment step. Its one relevant
  number — an ~12% miss rate for coelution detection at spectral r ≈ 0.86 versus ~2% at r ≈ −0.06 —
  concerns purity detection, not identity matching, and points the same direction: spectral
  contrast is the resource, and analogues have little of it.
- Nothing in the released data is isocratic; no t₀, dwell volume, or pH is recorded anywhere.
