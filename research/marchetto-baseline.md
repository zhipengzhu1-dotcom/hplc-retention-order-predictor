# The Marchetto baseline

Resolves #2. Part of the map, #1.

**Primary source.** A. Marchetto, M. Tirapelle, L. Mazzei, E. Sorensen, M. O. Besenhard,
*In Silico High-Performance Liquid Chromatography Method Development via Machine Learning*,
Anal. Chem. 2025, 97, 6991–7001. <https://doi.org/10.1021/acs.analchem.4c03466> — CC-BY 4.0.
Open-access full text and both SI files were fetched and read directly:

- Article PDF: <https://europepmc.org/articles/PMC11983366?pdf=render>
- PMC record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11983366/>
- SI-1 (PDF, 3 pages): <https://pubs.acs.org/doi/suppl/10.1021/acs.analchem.4c03466/suppl_file/ac4c03466_si_001.pdf>
- SI-2 (XLSX, retention data): <https://pubs.acs.org/doi/suppl/10.1021/acs.analchem.4c03466/suppl_file/ac4c03466_si_002.xlsx>

Everything below is cited to the article (by section/table/figure/page) or to the SI files, which were
parsed rather than summarised. Section 9 lists what could not be determined.

---

## 1. Verdict up front

**The paper is a proof of concept on 36 neutral small molecules, one column, one organic modifier, one
pH. It is not a validated retention predictor, and it is not the state of the art we must beat — it is
the *architecture* we are adopting, demonstrated at the smallest scale that demonstrates anything.**

Reproducing it splits cleanly in two:

| Layer | Reproducible? | Effort |
|---|---|---|
| **LSER + LSS** (36 solutes, leave-one-out) | **Yes, essentially fully.** All input data ship in SI-2; the fitting is two ordinary least-squares problems. | **Half a day.** |
| **QSPR** (E, S, A, B from structure) | **No, not as published.** Requires commercial alvaDesc descriptors and a specific 2017 snapshot of a third-party database. | **Weeks, and the result will not be bit-identical.** |

The half that contains the *chromatography* is reproducible today. The half that contains the
*machine learning* is gated on two proprietary/versioned dependencies (§7).

---

## 2. The pipeline as implemented

Three chained models, structure in, `k(φ)` out (article Figure 1, p6994):

```
SMILES ──alvaDesc──▶ 804 → 612 → 313 molecular descriptors
                              │
                              ▼
                    4 × ridge regression (QSPR)  ──▶  E, S, A, B
                                                        V ── McGowan volume, from structure directly
                                                        │
                                                        ▼
                          2 × LSER (12 fitted system params) ──▶ log k_w , S_S
                                                        │
                                                        ▼
                                       LSS:  log k = log k_w − S_S·φ  ──▶  k(φ)
```

### The three equations

**LSS theory** (article eq 1, p6992). `k ≡ (t_R − t_0)/t_0`; `φ` is the volume fraction of organic
modifier; `k_w` is the extrapolated retention factor at `φ → 0`; `S_S` is the solvent strength parameter.

```
log k = log k_w − S_S·φ                                                    (eq 1)
```

**LSER / Abraham solvation parameter model** (article eq 2, p6993). `E` excess molar refraction,
`S` dipolarity/polarizability, `A` H-bond acidity, `B` H-bond basicity, `V` McGowan volume;
`c, e, s, a, b, v` are the system constants.

```
log k = c + eE + sS + aA + bB + vV                                         (eq 2)
```

The article is explicit (p6993) that **eq 2 is unsuitable for method development on its own**, because
its system constants are unknown functions of `φ` and must be refitted for every mobile-phase
composition. This is the central move of the paper:

**LSER applied to the LSS parameters, not to `log k`** (article eqs 4–5, p6995), following
Poole & Atapattu (J. Chromatogr. A 2022, 1675, 463153, <https://doi.org/10.1016/j.chroma.2022.463153>):

```
log k_w = c_logkw + e_logkw·E + s_logkw·S + a_logkw·A + b_logkw·B + v_logkw·V   (eq 4)
S_S     = c_Ss    + e_Ss·E    + s_Ss·S    + a_Ss·A    + b_Ss·B    + v_Ss·V      (eq 5)
```

**This is the load-bearing idea.** `12 = 6 × 2` system parameters are fitted *once* per
column/modifier pair and are **independent of `φ`** (article p6995). `φ`-dependence lives entirely in
eq 1. That is what turns a static LSER into a method-development model, and it is exactly the seam our
architecture assumes.

**Warning for reproduction:** the authors used the **natural** logarithm, not log₁₀ —
"in our work, we obtained `log k_w` as the natural logarithm of the retention factor `k_w`" (p6996),
and the SI-2 data are natural logs. Their `S_S` values are therefore a factor of ln(10) ≈ 2.303 larger
than conventionally reported `S` values. Any comparison to the wider literature must correct for this.

---

## 3. The QSPR layer, in detail

| Item | Value | Source |
|---|---|---|
| Descriptor software | **alvaDesc** (Alvascience) | p6994 |
| Descriptor blocks | constitutional indices, molecular properties, topological indices, ring descriptors, connectivity indices, 2D autocorrelation, GETAWAY | p6994 |
| 3D descriptors | **Excluded** — "SMILES strings contain no detailed information on 3D molecular structures" | p6994 |
| Descriptor count | 804 computed → 612 after removing (near-)constant (`var > 0.01` filter) and one mostly-missing MD → **313** used | p6994, SI-1 §S.1–S.2 |
| Reduction method | Pairwise-correlation filter, threshold swept over `[0.1, 0.15, …, 0.95, 0.99, 1.5]`; **0.85 chosen**, giving 313 MDs; the *same* 313 MDs used for all four models | Figure 2 p6995, SI-1 §S.2 |
| Missing values | Imputed with the **training-set median** for that descriptor | SI-1 §S.1 |
| Scaling | Zero mean / unit variance fitted on training set, applied to test | SI-1 §S.1 |
| Model class | **Ridge regression** (least squares + L2 weight decay), one model per target — four independent models for `E, S, A, B` | eq 3, p6994 |
| Hyperparameter | `α` selected by **10-fold CV** on the training set | p6995 |
| Alternatives tried | Artificial neural networks — "did not improve the predictive performance significantly… (results not shown for sake of brevity)" | p6994 |
| `V` | **Not** predicted. McGowan volume, computed directly from structure (McGowan 1978; Zhao/Abraham/Zissimos 2003, refs 60–61) | p6993, p6994 |
| Observation:input ratio | 5120 : 313 = 16:1, cited against the Topliss–Costello 5:1 rule | p6995 |

Note the ANN negative result is **unreported** — no architecture, no metric, no figure. It cannot be
checked, and it is the one comparison that would tell us whether the QSPR layer is linear-limited.

---

## 4. Datasets

### 4a. QSPR training — the "Abraham Absolv" set

- Source: the **Abraham Absolv data set** taken from the **UFZ-LSER database v3.2.1 (2017)**
  (article ref 62; <http://www.ufz.de/lserd> → <https://web.app.ufz.de/compbc/lserd/public/start/>). p6994.
- **7881** small molecules initially.
- Curation (p6994, SI-1 §S.1): remove solutes missing any of E/S/A/B; remove unknown-SMILES molecules;
  remove duplicates; **restrict MW to 80–400 g/mol** → **6437**. Then remove the 36 LSER solutes to
  guarantee independence → **6401**.
- Split: **80:20 random** → **5120 train / 1281 test** (p6994–6995).
- Raw dataset MW range before the 80–400 filter: min 2.02 (hydrogen), median ≈ 180, max 1203 g/mol
  (SI-1 Table S.1). The MW window was chosen to bracket the Poole retention set and the SI states the
  bounds were "chosen somewhat arbitrarily, and may not be the best choice".

**Obtainable?** The UFZ-LSER database is publicly browsable and free to use with citation, but it is now
at **v4.1.1 (2025)** — the 2017 v3.2.1 snapshot the paper used is not offered for download from the live
site, and no bulk export was found. See §9.

### 4b. LSER calibration/validation — the Poole retention set

- **Kinetex XB-C18 (Phenomenex)**, water–acetonitrile mobile phase, isocratic (article p6996).
- Provenance: "kindly provided by Prof. Poole directly (personal reference from Wayne State University,
  04 January 2023)" (p6996). HPLC setup previously published: Atapattu, Poole, Praseuth,
  *Chromatographia* 2018, **81**, 373–385 (article ref 71).
- **48 solutes**, listed in article Table 1 (p6997), at `φ` = 10, 20, 30, 40, 50, 60, 70 % v/v.
- **Filter:** each solute's `log k` vs `φ` was fitted with eq 1 and the fit's MAPE computed;
  only solutes with **MAPE < 12 %** were kept ("threshold chosen arbitrarily") → **48 → 36** solutes
  (p6996–6997). The 12 rejected solutes are starred in Table 1.
- **Working range:** `φ` restricted to **0.2–0.7**; `φ = 0.1` excluded (p6996, p6997). Figure 4's shaded
  "region out of scope" makes the point that LSS linearity holds "strictly for 30 % < φ < 60 %" only (p6999).
- **Validation:** strict **leave-one-out**. 36 iterations; each time 35 solutes fit the 12 system
  parameters (`2 × 36 = 72` LSERs total) and the held-out solute is predicted (p6997). The 36 solutes
  were also excluded from QSPR train *and* test, so the pipeline evaluation is genuinely out-of-sample
  at both layers.
- **Applicability domain of the LSER fit is tiny.** For the 36 retained solutes: MW 92–236 g/mol,
  0 % nitrogen at the 50th percentile, 0 % halogen up to the 75th percentile (article Table 2, p6998).
  Compare the QSPR training set's 80–400 g/mol (SI-1 Table S.1). The authors say so themselves (p6999):
  "the data set used for LSER system parameter prediction contained only 35 … solutes with similar size,
  polarity, and functional group properties."

**Obtainable? Yes — this is the good news.** The full 48-solute table ships as SI-2 under CC-BY (§5).

### 4c. An inconsistency in the reported count

The article states the LSER data set "comprises the natural logarithm of retention factors of 48 solutes
… at 10, 20, 30, 40, 50, 60, and 70 % v/v … with occasionally missing data, for a total of 210 retention
factors" (p6996). **Parsing SI-2 directly gives 48 solutes and 316 non-empty `log k` values**, not 210
(per-`φ` counts: 48/48/48/48/48/41/35 for φ = 70/60/50/40/30/20/10 %).

210 is instead the size of the *validation* set — 36 retained solutes × 5–6 usable `φ` fractions, which
the article states correctly elsewhere ("210 values in total, combining the 36 solutes and the 5 to 6
organic modifier fractions", p6998; Figure 5 caption "the 210 predicted"). So the p6996 sentence
mis-attributes the validation count to the raw file. **Minor, but it will bite anyone checking their
reproduction against the stated `N`.**

---

## 5. What the SI actually ships

Two files, both retrieved and inspected.

**`ac4c03466_si_001.pdf` — 3 pages.**
- §S.1 Abraham Absolv dataset curation (prose + Table S.1: MW/nAT/nSK/ARR/%C/%H/%N/%O/%X percentiles).
- §S.2 Unsupervised variable reduction rationale (the correlation-threshold list; states the pairwise
  method is order-dependent — "variations in the order of the MDs … may determine changes in which MDs
  are removed", a real reproducibility hazard).
- §S.3 Figure S.1, two parity plots (LSS-only, LSS+LSER).

**`ac4c03466_si_002.xlsx` — 23 KB, one sheet.**
- Columns: `Molecule`, then `logk` at `φ` = 70, 60, 50, 40, 30, 20, 10.
- 48 rows, 316 numeric values, natural log. Solute names only — **no SMILES, no CAS, no InChI**.
- This is the file cited in-text as `PooleAtapattuOriginalData.xlsx`; the shipped filename differs.

**What the SI does NOT ship — all of it material:**

- ❌ **No code.** No repository, no GitHub link, no data-availability statement, no code-availability
  statement anywhere in the article. Confirmed against the Associated Content / Author Information /
  Notes sections (p7000).
- ❌ **No trained model weights** and no ridge `α` values.
- ❌ **No descriptor table** — neither the 804, the 612, nor the 313 selected alvaDesc descriptor names.
  This is the single biggest gap: the exact 313 are not recoverable without alvaDesc *and* the exact
  curated dataset *and* the same column ordering (see the SI's own order-dependence caveat).
- ❌ **No fitted LSER system constants.** The 12 values of `c/e/s/a/b/v` for `log k_w` and `S_S` appear
  nowhere in the article or SI. Anyone wanting the model must refit from SI-2.
- ❌ **No curated Abraham Absolv subset** — the 6401 molecules with their SMILES and E/S/A/B are not shipped.
- ❌ **No SMILES for the 48 retention solutes.** Names are unambiguous enough to resolve by hand, but it
  is manual work and the paper's exact structures (tautomers, salt forms) are not pinned.
- ❌ **No per-solute error breakdown**, no residual table.

---

## 6. Every reported error metric, with its conditions

### Metric definitions (article eqs 6–9, p6997)

- `RMSEP = sqrt( (1/M) Σ (y_exp − y_pred)² )`
- `R² = 1 − Σ(y_exp − y_pred)² / Σ(y_exp − ȳ_exp)²`
- `MMAPE = (1/M) Σ | (y_exp,m − y_pred,m) / ȳ_exp | · 100` — note the denominator is the **dataset mean**,
  not the per-point value. Used for QSPR because `A` and `B` are exactly zero for many solutes and a true
  MAPE would blow up (p6997). **MMAPE is not comparable to MAPE and not comparable across targets with
  different means.** Treat every MMAPE below as internal to this paper only.
- `MAPE = (1/M) Σ | (y_exp,m − y_fit,m) / y_exp,m | · 100` — the ordinary one, used for retention factors.

### 6a. QSPR layer — Abraham descriptor prediction

Read from article Figure 3 (p6996) insets. Conditions: **held-out test set, 1281 molecules**, ridge
regression on 313 alvaDesc MDs, trained on 5120.

| Target | RMSEP (test) | R² (test) | MMAPE (test) |
|---|---|---|---|
| **E** (excess molar refraction) | 0.10 | 0.98 | 6.8 |
| **S** (dipolarity/polarizability) | 0.23 | 0.87 | 13.6 |
| **A** (H-bond acidity) | 0.11 | **0.89** | **30.4** |
| **B** (H-bond basicity) | 0.11 | 0.95 | 12.0 |

Training-set numbers are not given separately; the article asserts Figure 3 "suggests the absence of
overfitting" (p6995) without a quantitative train/test comparison.

**`A` is the weak link** — worst MMAPE by a factor of >2 and second-worst R². `A` is hydrogen-bond
*acidity*: phenols, carboxylic acids, anilines. That is precisely the chemistry that ionises, and
precisely the chemistry this effort has committed to supporting. Flagging this hard.

Figure 2 (p6995) additionally reports RMSE under 10-fold CV as a function of descriptor count; the four
curves plateau between ~250 and ~310 descriptors, with `S` plateauing highest (~0.22) and `A`/`B` lowest
(~0.10). No tabulated values.

### 6b. Retention factor — the headline error-propagation result

Article **Table 3, p6998**. Conditions: **`k`, not `log k`**; the 36 LSER solutes under
**leave-one-out**; `φ` ∈ {0.2, 0.3, 0.4, 0.5, 0.6, 0.7}; **210 points**; experimental values "assumed to
be the true values (i.e., no experimental uncertainty associated)" (p6998); Kinetex XB-C18 /
water–ACN, single pH.

| Configuration | What is measured vs. what is predicted | MAPE | RMSEP | R² |
|---|---|---|---|---|
| **LSS only** | `log k_w`, `S_S` fitted per solute from that solute's own data | **9.1 %** | 0.98 | 0.98 |
| **LSS + LSER** | `E,S,A,B` **experimental**; `log k_w`, `S_S` from LSER | **17.8 %** | 1.8 | 0.95 |
| **LSS + LSER + QSPR** | everything from SMILES | **24.6 %** | 1.9 | 0.94 |

Read this as an error budget, and it is the most useful single result in the paper:

- **9.1 pp of the 24.6 % is the LSS functional form itself** — the price of assuming `log k` is linear in
  `φ` over 0.2–0.7. It is irreducible without a nonlinear solvent-strength model. The authors say the
  linearity "seems valid strictly for 30 % < φ < 60 % only" (p6999) and that fixing it would need
  quadratic solvent-strength models and hence *more* LSERs to predict the extra parameters — explicitly
  out of scope for them (p6999).
- **+8.7 pp from the LSER step** (9.1 → 17.8) with *perfect* solute descriptors. This is the fit quality
  of 12 system constants against 35 solutes.
- **+6.8 pp from the QSPR step** (17.8 → 24.6). **The QSPR layer is the *smallest* of the three error
  contributions.** RMSEP barely moves (1.8 → 1.9) and R² barely moves (0.95 → 0.94).

That last point is strategically important and it is the opposite of the intuitive assumption: with only
35 calibration solutes, the structure→descriptor ML is *not* the bottleneck. The physics layer and the
LSER calibration set are.

### 6c. Fit-quality thresholds used (not accuracy claims, but they gate the data)

- Per-solute LSS fit **MAPE < 12 %** to retain a solute (p6996–6997). "Threshold chosen arbitrarily."
- Variance filter `var > 0.01` on descriptors (SI-1 §S.1).
- Pairwise-correlation threshold **0.85** (SI-1 §S.2).

### 6d. Metrics that are NOT reported anywhere

No elution-order or pairwise-ranking metric. No resolution error. No per-solute residuals. No
uncertainty interval of any kind. No `log k` RMSE (RMSEP is on `k`, so it is dominated by the
long-retention points — the 0.98/1.8/1.9 values are in `k` units where `k` ranges to ~60, see Figure 5).

---

## 7. Reproducibility, component by component

| Component | Status | Blocker |
|---|---|---|
| LSS eq 1 per-solute fits | 🟢 **Reproducible now** | None. SI-2 has the data. |
| The 48 → 36 MAPE<12 % filter | 🟢 **Reproducible now** | None; deterministic given SI-2. |
| LSER eqs 4–5 fit + LOO | 🟢 **Reproducible now** | None. OLS on 35 × 6 design matrix; needs `E,S,A,B,V` for 36 named solutes, which are standard literature values. |
| Table 3 rows "LSS" and "LSS+LSER" | 🟢 **Reproducible now** | Only ambiguity is the source of the 36 experimental `E,S,A,B` (see §9). |
| Abraham Absolv curation (7881→6401) | 🟡 **Approximable** | UFZ-LSER is now v4.1.1; the v3.2.1 (2017) snapshot is not downloadable. Row count will differ. |
| 804 alvaDesc descriptors | 🔴 **Blocked** | **alvaDesc is commercial** (<https://www.alvascience.com/alvadesc/>) — login-gated, no stated free tier. |
| The specific 313 descriptors | 🔴 **Not recoverable** | Not listed anywhere; and SI-1 §S.2 admits the selection is **order-dependent** on the descriptor columns. |
| Table 3 row "LSS+LSER+QSPR" | 🔴 **Not exactly reproducible** | Depends on all of the above. |
| The ANN negative result | 🔴 **Uncheckable** | "results not shown". |

**Answer to the ticket's framing question — half-day or research project?**
**Both, and the split is clean.** The chromatographic core (§7 green rows, i.e. everything except the
QSPR layer) is a **half-day**: parse SI-2, fit 48 straight lines, filter to 36, look up Abraham
descriptors, run 72 OLS fits, reproduce MAPE 9.1 % and 17.8 %. **Do this first** — it independently
validates the LSER-on-LSS-parameters idea, which is the part of the architecture we are betting on, and
it does not depend on any proprietary tool.

The QSPR layer is **not a reproduction task at all** — it is a substitution task. We should not chase
alvaDesc; we should swap in an open descriptor stack (RDKit/Mordred) or the graph-convolutional and
group-contribution approaches from the other two indexed papers, and compare against Marchetto's
**MMAPE/RMSEP on E, S, A, B** as the benchmark to beat. Because the QSPR step contributes only ~7 pp of
the 24.6 % (§6b), a descriptor-stack swap that is merely *comparable* costs us almost nothing end-to-end.

---

## 8. Stated failure modes and limitations

Quoted or closely paraphrased from the article.

1. **No pH. At all.** "the LSS theory does not account for significant changes in pH, as its parameters
   are specific to a single pH value" (p6992). The conclusion repeats it: the method gives initial
   elution-time estimates "(not considering different pH values)" (p6999). **The paper's model as
   published cannot handle ionisable compounds' pH dependence.** The 36-solute set includes phenols,
   anilines and a carboxyl-free but H-bonding cast, all measured at one unstated pH.
2. **LSS linearity fails at the ends of the `φ` range.** "the LSS theory fails to capture a common
   nonlinear increase in `log k` as `φ → 0`… the linearity of `log k` vs `φ` seems valid strictly for
   30 % < `φ` < 60 % only" (p6999). Fixing it needs quadratic solvent-strength models, which needs
   *more* LSERs for the extra parameters — declined as out of scope.
3. **Calibration set far too small and too narrow.** "the data set used for LSER system parameter
   prediction contained only 35 … solutes with similar size, polarity, and functional group properties.
   Solutes that differ significantly in these properties are unlikely to allow for accurate prediction"
   (p6999).
4. **No applicability-domain measure — and they know it.** "(this and all) data-driven approaches should
   always provide a quantitative measure of similarity between the training data and the target solutes
   (e.g., molecular features captured through MDs) to indicate the applicability of a model" (p6999).
   **The paper names this as required future work and does not do it.**
5. **Experimental uncertainty ignored.** Retention factors "were assumed to be the true values (i.e., no
   experimental uncertainty associated)" (p6998).
6. **One column, one modifier.** Kinetex XB-C18 / water–ACN only. Generalisation across columns is
   asserted as a capability ("by integrating sensible system parameters for LSER … for various columns
   and organic modifiers, our methodology also enables in silico screening", p6999) but **never
   demonstrated**.
7. **Not a replacement for experiments.** "not yet capable of fully replacing experimental campaigns,
   especially for complex multicomponent samples prone to coelution" (p6999).
8. **Descriptor selection is order-dependent.** SI-1 §S.2.
9. **Two arbitrary thresholds** (MAPE < 12 %, correlation 0.85), acknowledged as arbitrary.

---

## 9. What I could NOT determine

Stated plainly, because these are the things that will cost time later.

1. **Whether the "Abraham Absolv" descriptors are experimental or computed.** The article calls them
   "experimentally determined LSER solute parameter databases" and lists UFZ-LSER, SoluteDB and the Wayne
   State set as examples, then says "we used the so-called Abraham Absolv data set (taken from the
   UFZ-LSER database)" (p6994). But **Absolv is ACD/Labs' *prediction* algorithm** — its product page
   describes it as "Calculate Abraham Solvation Parameters … directly from structure", with a bundled
   database of only ">5000 compounds" of literature values
   (<https://www.acdlabs.com/products/percepta-platform/absolv/>). A 7881-molecule "Absolv" table is
   therefore **likely to be at least partly ACD-*calculated*, not measured.** I could not confirm this
   from the live UFZ-LSER site (it is now v4.1.1 and did not expose per-dataset provenance to a plain
   fetch). **If true, the QSPR layer is a model trained to imitate another model**, and its reported
   R² = 0.98 on `E` measures agreement with ACD/Absolv, not with reality. This materially changes how we
   read every QSPR number in §6a and it is the most important open question this ticket produced.
2. **The exact pH / buffer of the Poole retention data.** Not stated in the article; presumably in
   Atapattu, Poole & Praseuth, *Chromatographia* 2018, 81, 373–385, which I did not obtain (paywalled).
   Unbuffered water–ACN is likely, which matters for the phenols and anilines in the set.
3. **The fitted 12 LSER system constants.** Not published (§5). Must be refitted.
4. **The 313 descriptor names.** Not published and, per SI-1 §S.2, not uniquely determined.
5. **The ridge `α` values.** Not published.
6. **Which `E,S,A,B,V` values were used for the 36 LSER solutes.** The article says the LSER solutes were
   *excluded* from the Absolv set, but does not say where their experimental descriptors then came from —
   presumably the same Absolv table, read before exclusion. Not stated.
7. **Whether the v3.2.1 (2017) UFZ snapshot is retrievable.** No download link found; the live site is
   v4.1.1 (2025). An exact reproduction of the 7881 → 6401 curation is therefore not possible.
8. **Train-set QSPR metrics.** Only test values are legible in Figure 3.
9. **Any ANN details.** "results not shown".
10. **`t_0` / column dead time and dimensions.** Not in the article; the SI ships `log k`, already
    dead-time-corrected, so this does not block reproduction of the modelling — but it does block any
    attempt to reconstruct `t_R` or model the instrument layer against this dataset.

---

## 10. Relevance to this effort

### 10a. The architecture is confirmed, at the level of the equations

The map's `SMILES → MDs → QSPR → E,S,A,B,V → LSER → k(φ)` chain is exactly eqs 4–5 + eq 1. The specific
commitment we should carry forward is **LSER applied to the LSS parameters `log k_w` and `S_S`, not to
`log k`** — this is what makes the system constants `φ`-independent and hence fittable once per column.
Adopt eqs 4–5 as the layer interface.

### 10b. Uncertainty: the paper offers **nothing**, and this is a clean gap for us

**The paper emits point predictions only.** No confidence intervals, no prediction intervals, no
posterior, no error bars on any figure, no per-solute uncertainty, and experimental error explicitly set
to zero (p6998). It does not support uncertainty estimation as published.

But the *architecture* is unusually friendly to retrofitting it, and cheaply:

- **Ridge regression is a MAP estimate under a Gaussian prior.** Bayesian ridge / Gaussian-process
  regression gives per-solute predictive variance on `E, S, A, B` at essentially no modelling cost and
  no loss of the paper's accuracy.
- **eqs 4–5 are linear.** Descriptor uncertainty propagates to `log k_w` and `S_S` in closed form, and
  OLS on 35 points yields a parameter covariance for the 12 system constants directly.
- **eq 1 is linear in `φ`.** So `Var(log k | φ)` is closed-form too. **The entire uncertainty budget of
  this pipeline is analytic** — no sampling needed for a first cut.
- The paper's own §8.4 limitation ("should always provide a quantitative measure of similarity between
  the training data and the target solutes") is an **applicability-domain request that the authors
  declined to fulfil**. That is our calibration story, and they wrote the ticket for us.

The three-row Table 3 also hands us a **ready-made variance decomposition** (LSS form / LSER fit / QSPR)
to validate any uncertainty budget against.

### 10c. Elution order: not measured by the paper — so I measured what I could

The paper reports no ranking metric. Since pairwise elution-order accuracy is our primary metric, I
computed the *intrinsic* order structure of SI-2 directly (script: pairwise comparison of `ln k` across
`φ`, 48 solutes):

| `φ` change | Comparable pairs | Pairs that reverse order | Pairs within 0.05 `ln k` at the higher `φ` |
|---|---|---|---|
| 20 % → 70 % | 820 | **63 (7.7 %)** | 24 |
| 30 % → 60 % | 1128 | 42 (3.7 %) | 25 |
| 20 % → 30 % | 820 | 20 (2.4 %) | 14 |
| 60 % → 70 % | 1128 | 14 (1.2 %) | 27 |

Two consequences, both uncomfortable:

1. **Only 7.7 % of pairs reverse across the entire 0.2–0.7 range.** A trivial baseline — "predict the
   order once at one `φ` and reuse it everywhere" — scores **92.3 %** pairwise accuracy on this dataset.
   Any elution-order metric we report **must be quoted against that baseline**, or a strong-looking
   number will mean nothing. Equivalently: getting `S_S` (the selectivity-with-`φ` term) right is worth
   only ~8 % of pairs here, while getting `log k_w` right is worth the other 92 %.
2. **~2–3 % of pairs are within 0.05 `ln k` at any given `φ`** — genuinely coeluting. There is a hard
   ceiling on order accuracy that no model reaches, and it is the ceiling that a *distribution* output
   is designed to express and a point output cannot. Direct support for the map's ranked-distribution
   decision.

This also means **this dataset is too easy to be our order benchmark.** 48 chemically similar neutral
solutes on one column is not a discriminating test of pairwise ranking. RepoRT will have to carry that.

### 10d. What Marchetto does **not** give us, that our destination requires

| Our requirement | Marchetto |
|---|---|
| Ionisable compounds, `k(φ, pH)` | ✗ explicitly excluded (§8.1) |
| Multiple columns | ✗ one column, capability asserted not shown |
| Calibrated uncertainty | ✗ none (§10b) |
| Ranked retention distributions | ✗ point estimates |
| Elution-order accuracy | ✗ never measured (§10c) |
| Applicability domain | ✗ named as future work, not done |
| Gradient / instrument layer | ✗ isocratic only; suggests coupling to transport models |
| Refinement from real injections | ✗ the whole point is *no* new experiments |

**So Marchetto is a floor, not a ceiling.** Nearly every axis the map cares about is unexplored by it.
The honest framing for the effort: we are not trying to beat 24.6 % MAPE — we are trying to build the
five things above, on an architecture Marchetto showed is sound at n=36.

### 10e. Recommended immediate action

Reproduce **only the green rows of §7** — the LSS + LSER half — from SI-2. Half a day, no proprietary
dependency, and it independently validates the one architectural claim we have already committed to.
Treat the QSPR half as a *substitution* benchmark (beat MMAPE 30.4 on `A` with an open descriptor stack),
not a reproduction target.

---

## 11. New questions this raises

Flagged for the map; none of these appear to have tickets yet.

1. **Is the Absolv training data measured or ACD-computed?** (§9.1) If computed, our QSPR layer's ground
   truth is a commercial model's output and every reported R² is an agreement-with-a-model number. This
   should be settled before any QSPR work starts — it changes what "good" means.
2. **What is the pairwise-order baseline on a hard dataset?** (§10c) The naive-order baseline must be
   defined and computed on RepoRT before any order accuracy is claimed, or our headline metric is
   uninterpretable.
3. **Is the LSS linear form good enough, or do we need a quadratic/Neue solvent-strength model from day
   one?** 9.1 pp of the 24.6 pp error is the LSS form itself (§6b), and it is the *largest* single
   contribution. The paper declines the fix because it needs extra LSERs for the extra parameters.
   That decision should be made deliberately by us, not inherited.
4. **Predicting `A` (H-bond acidity) well is a prerequisite, not a detail.** MMAPE 30.4 (§6a) on exactly
   the chemistry that ionises. Whether an open/graph-based QSPR fixes this is a specific, testable
   question and it gates the ionisable-compound scope.
5. **How many calibration solutes does an LSER system-constant fit actually need?** Marchetto used 35 and
   paid ~8.7 pp. This is a cheap learning-curve experiment on RepoRT and it directly sizes the
   per-column calibration burden for the curated 5–15 column set.
6. **Natural log vs log₁₀ for `S_S`.** (§2) A units convention that will silently corrupt any comparison
   against published `S` values. Worth pinning in `CONTEXT.md` when it exists.

---

## Sources

- Marchetto, Tirapelle, Mazzei, Sorensen, Besenhard, *In Silico High-Performance Liquid Chromatography
  Method Development via Machine Learning*, Anal. Chem. 2025, 97, 6991–7001.
  <https://doi.org/10.1021/acs.analchem.4c03466> (CC-BY 4.0).
  Full text read via <https://europepmc.org/articles/PMC11983366?pdf=render>.
- SI-1, `ac4c03466_si_001.pdf`:
  <https://pubs.acs.org/doi/suppl/10.1021/acs.analchem.4c03466/suppl_file/ac4c03466_si_001.pdf>
- SI-2, `ac4c03466_si_002.xlsx` (parsed directly for §4c and §10c):
  <https://pubs.acs.org/doi/suppl/10.1021/acs.analchem.4c03466/suppl_file/ac4c03466_si_002.xlsx>
- Poole & Atapattu, J. Chromatogr. A 2022, 1675, 463153 — the LSER-on-LSS-parameters formulation
  (article ref 70). <https://doi.org/10.1016/j.chroma.2022.463153>
- Atapattu, Poole & Praseuth, Chromatographia 2018, 81, 373–385 — the HPLC setup behind the retention
  data (article ref 71). <https://doi.org/10.1007/s10337-018-3477-5> *(not obtained; paywalled)*
- UFZ-LSER database (article ref 62 cites V3.2.1, 2017; live version is 4.1.1, 2025).
  <http://www.ufz.de/lserd> → <https://web.app.ufz.de/compbc/lserd/public/start/>
- ACD/Labs Absolv product page — evidence that "Absolv" denotes a prediction algorithm.
  <https://www.acdlabs.com/products/percepta-platform/absolv/>
- alvaDesc (Alvascience) — commercial, login-gated. <https://www.alvascience.com/alvadesc/>
