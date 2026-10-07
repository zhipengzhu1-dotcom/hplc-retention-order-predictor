# Ionisation providers: open vs commercial

Resolves #7. Part of the map, #1.

The pluggable/open-by-default decision is not relitigated here. This note answers: **which open
provider is the default**, **what the provider interface must return**, and — added mid-ticket from
the LSER system-constants result (#3) — **what the interface must return given that retention
depends on ionisation in the mobile phase, not in water**.

---

## TL;DR

1. **Open default: Uni-pKa** (Apache-2.0), consumed through the `unipka` PyPI wrapper (MIT). It is
   the only open candidate that is (a) a *macro*-pKa predictor built on an explicit
   **protonation-ensemble** formalism, so multiprotic and zwitterionic molecules are handled by
   construction; (b) benchmarked in 2026 as **comparable to, and in two cases significantly better
   than, commercial predictors**; and (c) pip-installable with no MATLAB Runtime, no CUDA compile,
   and no licence server.
2. **OPERA is not the default.** Its pKa model is RMSE ≈ **1.5–1.8 pKa units** on its own test set,
   returns only the *strongest* acidic and *strongest* basic pKa, and cannot express a microspecies
   distribution. Its one genuine advantage — a per-prediction applicability-domain and accuracy
   index — is real and should be *borrowed as a design pattern*, not adopted as the engine.
3. **No candidate emits an honest per-prediction uncertainty on the pKa itself.** This is the
   ticket's central negative finding. The uncertainty layer must therefore be **built by us** on top
   of the provider (ensemble spread + a conformal/empirical calibration set), not taken from it.
4. **The provider contract is water-only, and water is the wrong solvent.** Every candidate —
   open and commercial alike — predicts aqueous pKa. Retention depends on the ionised fraction in
   the *aqueous/organic mobile phase*. A **solvent-correction layer is mandatory**, it is well
   characterised in the literature for both ACN/water and MeOH/water, and it needs an input the pKa
   models do not currently surface: a **charge-type / compound-family label per ionisable site**.
5. **The pH scale must be standardised on `s_s`pH** (or `s_w`pH plus the tabulated δ), *not* the
   convenient `w_w`pH. Using `w_w`pH is a silent systematic error and is specifically documented to
   break for bases.

---

## 1. Why mean pKa error is the wrong benchmark

### 1.1 The retention model

For a monoprotic solute in RPLC, retention is a sigmoidal function of mobile-phase pH — the
weighted average of the neutral and ionised forms' retention factors
([Anal. Chem. 2007, 79, 3180](https://doi.org/10.1021/ac060482i); the canonical treatment is the
Rosés/Bosch "Retention of ionizable compounds on HPLC" series, e.g.
[Anal. Chem. 2000, 72, 5193](https://doi.org/10.1021/ac000591b)):

```
k_obs = (k_HA + k_A · 10^(pH − pKa)) / (1 + 10^(pH − pKa))
```

Write `Δ = pH − pKa` and `D = log10(k_HA / k_A)` — the **ionisation retention drop**, the total
log-k fall between fully neutral and fully ionised. In RPLC `D` is typically 1–2 log units. Then

```
log k(Δ) = log k_HA + log10[ (1 + 10^(−D) · 10^Δ) / (1 + 10^Δ) ]
```

A pKa error `δ` translates the whole sigmoid along the pH axis by `δ`. The downstream log-k error is
`|f(Δ − δ) − f(Δ)|`, which is **not** proportional to `δ` and **not** uniform in `Δ`.

### 1.2 Peak sensitivity

Numerically differentiating the expression above (`research/` computations, reproducible from the
formula):

| `D` | max \|d log k / d pH\| | located at `pH − pKa` = |
|---|---|---|
| 1.0 | 0.52 | ±0.50 |
| 1.5 | 0.70 | ±0.75 |
| 2.0 | 0.82 | ±1.00 |

The sensitivity peak sits slightly *off* the crossover (at roughly `D/2`), not exactly at
`pH = pKa` — a small but useful correction to the ticket's framing. It is essentially flat across
`|pH − pKa| ≲ 1`, so the "danger band" is about **±1.5 pH units wide**, not a knife edge.

### 1.3 pKa error → log k error

Downstream `|Δ log k|` for `D = 1.5`, as a function of the true distance from crossover and the pKa
error `δ` (worst case over the sign of `δ`):

| \|pH − pKa_true\| | δ = 0.5 | δ = 1.0 | δ = 1.5 | δ = 2.0 |
|---|---|---|---|---|
| 0.0 | 0.29 | 0.64 | 0.93 | 1.10 |
| 1.0 | 0.34 | 0.64 | 0.93 | 1.10 |
| 2.0 | 0.17 | 0.46 | 0.81 | 1.10 |
| 3.0 | 0.03 | 0.10 | 0.27 | 0.57 |

**Rule of thumb for the spec:**

```
worst-case |Δ log k|  ≈  min( c(D) · δ , D )        c(1.0)=0.5, c(1.5)=0.7, c(2.0)=0.8
```

The transfer coefficient is remarkably stable (0.66–0.69 across δ = 0.3…1.0 at `D = 1.5`), so a
linear rule is safe up to δ ≈ 1.5, after which it saturates at `D`.

**What this means concretely.**
- A **1-unit** pKa error near crossover costs **~0.6–0.8 log k units**. For reference, that is larger
  than the whole between-column selectivity spread that LSER is being asked to resolve. It will
  reorder peaks.
- The error **saturates at `D`**. A provider that is wrong by 3 units is not 3× worse than one wrong
  by 1 unit — it is at most `D` wrong. So the marginal value of accuracy is concentrated entirely in
  the δ ∈ [0, 1.5] regime. This is precisely the regime where the open/commercial gap has closed
  (§3) and where the *solvent correction* (§5) is the same size as the model error.
- Errors at `|pH − pKa| > 3` are chromatographically free. Screening out those cases first makes the
  pKa provider's accuracy irrelevant for most of the library.

### 1.4 Why honest uncertainty beats a better point estimate

Because the map's product is a **ranked retention distribution with calibrated confidence**, the
provider does not need to be right — it needs to know when it isn't. `pKa = 8.2 ± 1.5` propagates
through the sigmoid to a wide, correctly-shaped log-k distribution and a low-confidence elution-order
call; `8.4` flat propagates to a sharp, confidently wrong one. Given the saturation in §1.3, a
well-calibrated σ is worth more than the ~0.2 pKa units that separate the best open model from the
best commercial one.

**This is the finding that most constrains the architecture: no candidate provides that σ.** See §4.

---

## 2. Candidate-by-candidate

### OPERA (US EPA / NIEHS)

| | |
|---|---|
| Licence | **MIT** ([repo](https://github.com/kmansouri/OPERA)) |
| Paper | Mansouri et al., *J. Cheminform.* **11**, 60 (2019), [10.1186/s13321-019-0384-1](https://doi.org/10.1186/s13321-019-0384-1); platform paper [10.1186/s13321-018-0263-1](https://doi.org/10.1186/s13321-018-0263-1) |
| Install | **Compiled MATLAB application** + MATLAB Runtime; Windows and Linux only, invoked as `./run_OPERA.sh <mcr_dir> <args>` ([CLI help](https://github.com/kmansouri/OPERA/blob/master/OPERA_Source_code/help.txt)). **No official Python package.** A third-party wrapper exists ([pyOPERA](https://cabreratoxy.github.io/pyOPERA/usage.html)). macOS is not listed as supported. |
| Accuracy | Test-set **RMSE ≈ 1.5–1.8, R² ≈ 0.72–0.80** across SVM/kNN, XGB and DNN variants (acidic: RMSE 1.51–1.80; basic: RMSE 1.53–1.69) |
| Chemistry | Trained on **7912 chemicals from DataWarrior** — environmental/general chemistry, not curated drug-like |
| Multiprotic / zwitterionic | **Only the strongest acidic and strongest basic pKa.** Amphoterics were handled by three different dataset-construction "options" during training; there is no microspecies model |
| Tautomers / microspecies | **None** |
| Uncertainty | **Yes, and it is the best in the field**: every prediction carries `AD` (Boolean global applicability domain, leverage-based), `Sim_index` (continuous local AD, 0–1) and `Conf_index` (accuracy estimate), plus optional 5 nearest training neighbours with their observed values |
| logD | Yes, as a separate endpoint |

**Verdict: reject as default.** RMSE 1.5–1.8 sits squarely in the saturating regime of §1.3 —
expected worst-case ~0.9–1.1 log k units near crossover. Combined with strongest-acid/strongest-base
output only (no microspecies) and a MATLAB-Runtime install that would make #20 painful, it fails on
three of the six axes. **Keep it as the uncertainty design template**: `pred / AD / Sim_index /
Conf_index / nearest_neighbours` is close to exactly the record our interface should emit, and OPERA
should stay implemented as a *secondary* provider for cross-checking and for its regulatory
acceptance.

### MolGpKa

| | |
|---|---|
| Licence | **MIT** ([repo](https://github.com/Xundrug/MolGpKa)) |
| Paper | Pan et al., *J. Chem. Inf. Model.* **61**, 3159 (2021), [10.1021/acs.jcim.1c00075](https://doi.org/10.1021/acs.jcim.1c00075) |
| Install | Not on PyPI (`pypi.org/pypi/molgpka/json` → 404, checked 2026-08-14). Clone + conda env; pinned to Python 3.6, scikit-learn 0.21.3, numpy 1.18, pandas 0.25 — **a dead environment in 2026**. Docker file provided |
| Accuracy | MAE **0.81** (Novartis) / **0.49** (Baltruschat), per the pKaLearn comparison table |
| Chemistry | **Trained on Epik-computed labels**, not experiment: 1,624,715 ChEMBL compounds → 963,969 acidic and 1,167,260 basic sites, labelled by Schrödinger Epik |
| Multiprotic / zwitterionic | Per-site **micro**-pKa only. In the 2026 benchmark, micro-pKa predictors including MolGpKa were **excluded from the polyprotic evaluation set** because they do not produce macroscopic charge-state transitions |
| Tautomers / microspecies | Not addressed |
| Uncertainty | **None** |

**Verdict: reject.** Distilling a commercial tool's *predictions* means inheriting its biases without
its validation, the dependency pins are unresolvable, and micro-pKa alone cannot express the
macroscopic charge-state populations the LSER layer needs.

### QupKake

| | |
|---|---|
| Licence | **BSD 3-Clause** ([LICENSE](https://github.com/hutchisonlab/QupKake/blob/main/LICENSE)); paper CC-BY 4.0 |
| Paper | Abarbanel & Hutchison, *J. Chem. Theory Comput.* **20**, 6946 (2024), [10.1021/acs.jctc.4c00328](https://doi.org/10.1021/acs.jctc.4c00328) ([PMC11325546](https://pmc.ncbi.nlm.nih.gov/articles/PMC11325546/)) |
| Install | Not on PyPI (404, checked 2026-08-14). Clone + `pip install .`; needs torch ≥2.0, PyTorch Geometric ≥2.3, PyTorch Lightning, RDKit, and **xtb 6.4.1 which the README says must be built from source** because of conda bugs. Linux binaries bundled |
| Accuracy | Novartis: **RMSE 0.79 / MAE 0.55 / R² 0.88**. Literature set: RMSE 0.54 / MAE 0.39 / R² 0.95. SAMPL6 RMSE 0.44; SAMPL7 RMSE 0.72 (2nd); SAMPL8 RMSE 0.85. **~0.3 s/molecule** single core |
| Chemistry | Pretrained on ChemAxon-computed ChEMBL pKa, then transfer-learned on experimental data. The 2026 benchmark notes it *outperforms ChemAxon itself* despite this |
| Multiprotic / zwitterionic | Explicitly limited: *"we can only have high confidence in the most acidic or basic micro-pKa values"*, with *"less confidence … in compounds with multiple sites"*. Tautomer search is *"focused on neutral compounds"*; **zwitterions are named as future work** |
| Tautomers | Yes — GFN2-xTB tautomer search with implicit water, then GNN site enumeration |
| Uncertainty | **None** (point predictions) |

**Verdict: strong runner-up, not the default.** Accuracy is excellent and the xTB tautomer search is
the most physically grounded of the open options. But the explicit zwitterion gap is disqualifying
for a system whose scope is drug-like ionisable compounds, and the xtb-6.4.1-from-source dependency
is exactly the kind of thing that sinks #20. **Recommend as the optional second open backend**, for
cross-provider disagreement as an uncertainty signal (§4).

### Uni-pKa — **recommended open default**

| | |
|---|---|
| Licence | **Apache-2.0** ([repo](https://github.com/dptech-corp/Uni-pKa)); the `unipka` wrapper is **MIT** |
| Paper | Zheng et al., *JACS Au* **4**, 2467 (2024), [10.1021/jacsau.4c00271](https://doi.org/10.1021/jacsau.4c00271); preprint [10.26434/chemrxiv-2023-lw5k0](https://doi.org/10.26434/chemrxiv-2023-lw5k0) |
| Install (official) | Painful: needs **Uni-Core** (`unicore` is not on PyPI — 404), Uni-Mol, a CUDA-flavoured Docker image `dptechnology/unimol:latest-pytorch1.11.0-cuda11.3` |
| Install (recommended) | **`pip install unipka`** — [PyPI](https://pypi.org/project/unipka/) v0.0.6, MIT, deps are only `torch>2.3, rdkit, scipy, numpy, pandas, scikit-learn, requests, platformdirs, seaborn, matplotlib, anywidget`. Source: [finlayiainmaclean/unipka](https://github.com/finlayiainmaclean/unipka), Apache-2.0 |
| Accuracy | **Best open performer** in the 2026 independent benchmark: the simple-template variant *"returns comparable (and, in two cases, significantly better) performances compared to commercial predictors"* — vs ACD/Labs Classic, ACD/Labs GALAS, Chemaxon and Epik ([10.1021/acs.jcim.6c00107](https://doi.org/10.1021/acs.jcim.6c00107)) |
| Chemistry | Pretrained on >1M ChEMBL-derived pKa, finetuned on experimental sets reconstructed into the ensemble formalism. ChEMBL-derived ⇒ drug-like by construction |
| Multiprotic / zwitterionic | **This is the reason to pick it.** It models the *full space of protonation microstates* and derives macro-pKa from micro-pKa equilibria rigorously, so successive pKa values are **thermodynamically consistent by construction**. It was one of only four tools evaluated on the polyprotic benchmark set |
| Tautomers / microspecies | Yes — a template-matching **microstate enumerator** (A-mode / B-mode) builds the ensemble from one SMILES; the model predicts a **free energy per microstate**, from which **pH-dependent microstate populations** follow |
| Uncertainty | **No calibrated σ.** But the official inference pipeline averages a **5-fold model ensemble**, so ensemble spread is available as a raw dispersion signal — the only open candidate where this is free |
| logD | Not native; the `unipka` wrapper implements it as a microstate-population-weighted average of per-microstate logP, following the Starling methodology |

**Caveats, stated plainly.**
- The official repo is thin: **23 commits, ~112 stars, 8 open issues**. Not a maintained product.
- The wrapper is **explicitly unofficial**, **v0.0.6**, single release, one author. It downloads
  weights at runtime via `requests`. **Vendor the weights and pin the version**; do not let #20 rely
  on a live download.
- **Coverage is worse than commercial.** In the 2026 benchmark, over 31,250 compounds, commercial
  tools threw exceptions on <1%, while all open predictors exceeded 1% — *"as high as 6.6%
  (Uni-pKa with the simple template) and 7.4% (pKaSolver)"*. **The interface must have a
  first-class "provider declined" return**, and roughly 1 in 15 compounds will use it.

### ChemAxon `cxcalc` / Python API — commercial reference

| | |
|---|---|
| Vendor | **Now Certara.** Certara completed the acquisition of Chemaxon on **2 October 2024** ([press release](https://www.certara.com/pressrelease/certara-completes-acquisition-of-chemaxon/)). All `chemaxon.com` licensing URLs now 301 to `certara.com/cxn/…` |
| Install | **Genuinely easy**: `pip install chemaxon`, wheels on [PyPI](https://pypi.org/project/chemaxon/) v26.1.0 for cp312/313/314 × macOS-arm64, manylinux x86_64 + aarch64, win_amd64. Licence field: *"Chemaxon License"*. Requires a licence file placed at a default location or pointed at explicitly |
| API surface | `calculations.pka(mol, consider_tautomerization, calculate_micro, use_large_model, min_basic, max_acidic, number_of_acidic_values, number_of_basic_values, temperature, required_pka_type) → PkaResult`; `major_microspecies(mol, ph, …) → Molecule`; `logd(mol, ph, method, …) → float`; `logd_ph_range(...) → List[LogDResult]`; `isoelectric_point(...) → IsoelectricPointResult` with `charge_distributions` ([apidocs](https://apidocs.chemaxon.com/python_api/apidocs/chemaxon/calculations.html)) |
| **Gap** | **The Python API exposes `major_microspecies` but no microspecies-*distribution* function.** Fractional populations are only reachable through the `cxcalc` CLI functions `msdistr` / `microspeciesdistribution` ([cxcalc functions](https://docs.chemaxon.com/latest/cxcalc-calculator-functions.html)). So the "high-accuracy backend" is the one that *cannot* natively satisfy the interface in-process |
| Multiprotic | Both macro and micro pKa modes; distribution charts limited to **≤8 ionisable atoms** (above that only macrospecies), and species below **0.01%** are suppressed ([pKa plugin](https://docs.chemaxon.com/latest/calculators_pka-plugin.html)) |
| Tautomers | `consider_tautomerization` uses the canonical tautomer and *"gives the same pKa values for different tautomers"* |
| Uncertainty | **None documented anywhere in the plugin docs** |
| Accuracy | Beaten by QupKake and matched-or-beaten by Uni-pKa in the 2026 independent benchmark; RMSE 1.145 on Novartis-Base per the Starling comparison table |

**Academic licence — what I could and could not establish.**

*Established:*
- The programme is real, worldwide, and aimed at *"academic individuals and institutions"* doing
  *"non-commercial chemical research"*
  ([academic programme](https://www.certara.com/cxn/academic-program),
  [research licence](https://www.certara.com/cxn/research-license)).
- Historically it comprised an **Individual Research License** (recurring 2-year provision of most
  ChemAxon tools), **free Teaching Licenses**, and academic discounts on the rest; the free academic
  package (**"AcPack"**) ran **2 years from issue**, conditional on **crediting ChemAxon in any
  publication or presentation** of data generated with it.
- Application is a **contact form only** — *"Leave your contact information… Our colleagues will get
  back to you soon."* There is no self-service path.

*Could NOT establish, and #20 must treat as unknown:*
- **Whether Calculator Plugins / `cxcalc` / the pKa plugin are inside the academic bundle.** They are
  not named on either current Certara page. This is the single most important open question, because
  a licence without the pKa plugin is worthless here.
- **Any stated lead time.** No turnaround is published anywhere I could find.
- **Whether the pre-acquisition AcPack terms survived the Certara migration.** The pages were
  rewritten; the 2-year/free-for-academics language is no longer on them.
- **Whether an effort like this qualifies.** The gating phrase is *"eligible institutions"* and
  *"non-commercial research project"*. A university-affiliated methodological study plausibly
  qualifies. An **unaffiliated individual almost certainly does not**, and neither page says
  otherwise. Independent-researcher eligibility is explicitly not clarified.

**Action for #20: email Certara asking three specific questions** — (1) does the academic licence
include Calculator Plugins/`cxcalc` and the pKa plugin; (2) what is the term and cost; (3) does an
individual researcher without an institutional affiliation qualify. **Do not block on the answer**;
the open path stands alone.

### Others assessed and not recommended

- **pKaLearn** — Genzling et al., *Commun. Chem.* **9**, 181 (2026),
  [10.1038/s42004-026-01983-y](https://doi.org/10.1038/s42004-026-01983-y). CC-BY 4.0,
  [GitHub](https://github.com/MoitessierLab/pKaLearn). MAE **0.68 ± 0.03** (Novartis) / **0.34 ±
  0.02** (Baltruschat); beats Marvin, Epik, pkasolver; comparable to QupKake. Iteratively evaluates
  multiple ionisation centres and *"generates the full ionization profile"*. **Very new — published
  2026 — and not in the 2026 independent benchmark.** Worth revisiting once it has third-party
  validation. The ± values are seed spread, not per-prediction uncertainty.
- **pkasolver** — [repo](https://github.com/mayrf/pkasolver). Micro-pKa; worst exception rate in the
  benchmark (7.4%). Reject.
- **Starling (Rowan)** — [preprint](https://chemrxiv.org/engage/chemrxiv/article-details/68388349c1cb1ecda02ba65d).
  Uni-pKa architecture, single model rather than the 5-fold ensemble, and correspondingly *worse*
  than Uni-pKa (RMSE 0.790 vs 0.653–0.878 on Novartis-Base). Access is via the commercial Rowan
  platform; **I could not establish that weights or code are open**. Its *methodology* is directly
  useful: pH-dependent microstate populations computed over pH 0–14 at 0.1 increments, and logD
  decomposed into (per-microstate weight) × (per-microstate logP). Cite it; don't depend on it.

---

## 3. The independent benchmark that settles the open-vs-commercial question

Sipos-Szabó, Bajusz, Balogh & Keserű, *"Benchmarking pKa Prediction Algorithms against an Extensive,
Public Data Set"*, *J. Chem. Inf. Model.* **66**, 4607–4619 (2026),
[10.1021/acs.jcim.6c00107](https://doi.org/10.1021/acs.jcim.6c00107)
([PMC13126629](https://pmc.ncbi.nlm.nih.gov/articles/PMC13126629/)).

- **Scale**: >90,000 experimental aqueous pKa values over 31,250 unique molecules (the *pKaHub*
  database), reduced to 37,654 distinct pKa values.
- **Compared**: ACD/Labs Classic, ACD/Labs GALAS, Chemaxon, Epik (commercial); MolGpKa, pKaSolver,
  QupKake, Uni-pKa simple + full template (open).
- **Split by protic class**: monoprotic / amphoteric / polyprotic. Micro-pKa predictors (MolGpKa,
  pKaSolver, QupKake) were **only evaluated on monoprotic and amphoteric sets**, since they emit
  site transitions rather than macroscopic charge-state transitions. Macro predictors (ACD, Chemaxon,
  Epik, Uni-pKa) were evaluated on all three.
- **Headline**: *"mean absolute errors rarely exceed one pKa unit for any of the predictors on any of
  the data sets"* and open-source tools, *"especially Uni-pKa, have caught up with commercial
  products"*.
- **Where commercial still wins**: coverage. <1% exception rate vs >1% for every open tool, up to
  6.6% (Uni-pKa simple) and 7.4% (pKaSolver).
- The authors themselves make our argument: errors matter most near physiological pH, where
  mis-assigned microspecies *"could propagate to downstream molecular modeling"*.

**Not determined**: the per-method RMSE/MAE tables live in Supplementary Tables S4–S6, which I could
not retrieve. The main text gives only the comparative statements quoted above. **#20 or #14 should
pull the SI** if a numeric provider ranking is needed.

---

## 4. Uncertainty: the gap we must fill ourselves

| Provider | Per-prediction uncertainty |
|---|---|
| OPERA | **Yes** — `AD`, `Sim_index` (0–1), `Conf_index`, + 5 nearest neighbours |
| Uni-pKa | No calibrated σ; **5-fold ensemble spread available** |
| QupKake | No |
| MolGpKa | No |
| pKaLearn | No |
| ChemAxon | No |

Given §1.4, this is the ticket's most consequential negative result. The recommended construction,
in ascending order of effort:

1. **Ensemble spread** from Uni-pKa's 5 folds — free, already computed, but known to be
   over-confident for NN ensembles.
2. **Cross-provider disagreement**: |Uni-pKa − QupKake| (and OPERA where in-domain) as a second
   dispersion channel. Cheap, and disagreement between architecturally unrelated models is a
   better-behaved signal than within-family spread.
3. **OPERA's AD/Sim_index as a domain gate**: run it purely to answer *"is this molecule in a region
   where anything has been validated?"*, independent of its pKa value.
4. **Conformal calibration** against a held-out slice of pKaHub, split by protic class
   (mono/ampho/poly) and by whether the site is the strongest of its type. This produces the
   *honest* σ the ticket asks for and is the piece that must be built. **It should be calibrated on
   downstream log k, not on pKa** — §1.3 shows the two are not monotonically related.

**Interface consequence**: `pka` is `float | None` and σ is **required, not optional**. A provider
that cannot supply σ must be wrapped by something that can, and the wrapper is ours.

---

## 5. **The mobile-phase correction — added from #3, and the most consequential section**

### 5.1 The problem

The retention layer needs the ionised fraction **α in the mobile phase at composition φ**. Every
candidate provider — open and commercial — predicts pKa **in water**. Substituting an aqueous pKa
into the retention sigmoid produces a biased α and therefore biased retention, and the bias is of
**the same magnitude as the model error** analysed in §1.3. Getting the provider right and the
solvent wrong buys nothing.

### 5.2 How pKa shifts with organic modifier — and yes, it is predictable

Direction is set by **charge type**, via the electrostatic (Born) term: adding a low-dielectric
modifier destabilises charge separation.

- **Neutral acids** (carboxylic acids, phenols): dissociation *creates* charge, so pKa
  **increases** with modifier content.
- **Cationic acids** (protonated amines, pyridinium): dissociation *moves* charge rather than
  creating it, so the electrostatic term largely cancels; pKa changes far less, and in ACN/water
  the pKa values of bases **decrease** with acetonitrile content.

The quantitative correction is established:

> **Espinosa, Bosch & Rosés**, *"Retention of ionizable compounds in HPLC. 14. Acid–base pK values
> in acetonitrile–water mobile phases"*, *J. Chromatogr. A* **964**, 55 (2002),
> [10.1016/S0021-9673(02)00558-7](https://doi.org/10.1016/S0021-9673(02)00558-7):
> *"Linear relationships between `s_s`pKa values in acetonitrile–water mixtures and `w_w`pKa values
> in pure water have been established for five families of compounds: aliphatic carboxylic acids,
> aromatic carboxylic acids, phenols, amines, and pyridines. The parameters (slope and intercept) of
> the linear correlations have been related with acetonitrile–water composition. The proposed
> equations allow accurate estimation of the pKa values of any member of the studied families at any
> acetonitrile–water composition **up to 60% of acetonitrile in volume** (100% for pyridines)."*

So the correction has the form

```
s_s pKa(family, φ)  =  a(family, φ)  +  b(family, φ) · w_w pKa
```

The methanol/water analogue is established the same way — e.g. for phenols, *"the pKa values of the
whole set of phenols at any methanol–water composition are linearly related to the pKa values of the
phenols in water"*, with slope/intercept given as functions of composition
([*J. Chromatogr. A* **867**, 41 (2000), 10.1016/S0021-9673(99)01139-5](https://doi.org/10.1016/S0021-9673%2899%2901139-5)).
The compendium covering both modifiers and the common buffers is:

> **Subirats, Rosés & Bosch**, *"On the Effect of Organic Solvent Composition on the pH of Buffered
> HPLC Mobile Phases and the pKa of Analytes — A Review"*, *Sep. Purif. Rev.* **36**, 231–255 (2007),
> [10.1080/15422110701539129](https://doi.org/10.1080/15422110701539129).

**Not determined**: I could not retrieve the actual `a`/`b` coefficient tables — both primary papers
are paywalled and returned 403. **#20 must obtain the Espinosa 2002 tables and the Subirats review
before implementing the correction.** Treat the coefficients as a required input, not something to
re-derive.

**Limits to respect:** valid to **60% v/v ACN** (except pyridines); families are restricted to the
five above, so a molecule outside them has *no* published correction and must be flagged
low-confidence. There is **no ML pKa model for mixed aqueous–organic solvents**. The nearest work is
pure non-aqueous solvents — Zheng et al., *J. Comput. Chem.* **46**, e27517 (2025),
[10.1002/jcc.27517](https://doi.org/10.1002/jcc.27517) — which reports overall MAE 0.89 across 46
solvents but **>1.2 pKa units for acetonitrile specifically**, worse than the aqueous models. It
does not cover binary mixtures at a working φ. **Do not wait for a solvent-aware model; use the
Rosés correction.**

### 5.3 The pH scale — standardise on `s_s`pH

Three conventions, IUPAC-named:

| Scale | Meaning |
|---|---|
| `w_w`pH | pH of the **aqueous buffer before mixing**. What everyone actually writes in a method. |
| `s_w`pH | pH measured **in the mixed mobile phase**, electrodes calibrated with **aqueous** standards. Directly measurable. |
| `s_s`pH | pH measured in the mixed mobile phase, electrodes calibrated **in the same mixture**. The thermodynamically rigorous one. |

Two hard results:

1. **`s_s`pH is the scale that connects to retention.** *"The `s_s`pH value of any buffered
   acetonitrile/water mobile phase used in reversed-phase liquid chromatography … **is directly
   related to the ionized fraction of analyte and, therefore, to its average retention***"
   ([*Anal. Chem.* **79**, 3180-something / 10.1021/ac062372h](https://doi.org/10.1021/ac062372h)).
2. **`w_w`pH specifically fails for bases.** *"The use of the rigorous `s_s`pH and `s_w`pH scales
   allows one to explain the retention behavior of bases, which in many instances cannot be
   justified from"* the `w_w`pH scale
   ([*Anal. Chem.* **72**, 5193 (2000), 10.1021/ac000591b](https://doi.org/10.1021/ac000591b)).

The conversion is tabulated: `s_s`pH = `s_w`pH − δ, with δ determined for **ACN/water 0–90% at
15–60 °C** and fitted to a simple joint function of composition and temperature
([10.1021/ac062372h](https://doi.org/10.1021/ac062372h)). For gradients, the apparent pH itself
drifts as φ changes, and equations relating gradient retention to the aqueous buffer pH — accounting
for the analyte pKa and buffer pH varying *together* — are given in
[*Anal. Chem.* **73**, 5658 (2001), 10.1021/ac0101454](https://doi.org/10.1021/ac0101454).

**Decision for the spec:** the interface's canonical pH is **`s_s`pH at the working (φ, T)**. The
user-facing input stays `w_w`pH (that is what a chemist dials in), and the system converts:
`w_w`pH → `s_w`pH → `s_s`pH via buffer model + tabulated δ(φ, T). **The conversion must be a named,
tested module, and every pH crossing a module boundary must carry its scale as a tag.** An untagged
pH is a bug.

### 5.4 Modifier coupling and silanol effects

Both MeOH/water and ACN/water are in scope (LSER system constants are published for both), and the
pKa correction must be modifier-specific — the ACN and MeOH correlations have different coefficients
and different validity ranges.

On silanol cation exchange: I **could not confirm** the framing that it is significant in ACN but not
methanol. The primary evidence I found reports it clearly **in methanol** systems — anions excluded
by ionised silanols and cations retained by cation exchange with the mobile-phase background cation
in 60% MeOH ([*J. Chromatogr. A* **910**, 1 (2001), 10.1016/S0021-9673(00)01192-4](https://doi.org/10.1016/S0021-9673%2800%2901192-4)),
and silanol acidity ranked across six columns in 60% MeOH, with XTerra MS C18 showing none up to
`s_s`pH 10.0 while Resolve/Symmetry silica show two distinct silanol types
([*J. Chromatogr. A* **982**, 1 (2002), 10.1016/S0021-9673(02)01899-X](https://doi.org/10.1016/S0021-9673%2802%2901899-X)).

**Conclusion**: silanol cation exchange is a *column* property more than a *modifier* property, and it
interacts with ionisation directly — it acts on the **cationic microspecies only**, so its magnitude
is proportional to the same α the provider computes. **Architectural implication**: the LSER layer
should receive **per-microspecies fractions with their charges**, so a silanol/ion-exchange term can
be attached to the charged fraction later without re-opening the ionisation interface. If it were
given only a scalar α or a logD, that extension would be impossible. **This is an independent reason
to reject logD-shaped interfaces.** The claim about ACN-vs-MeOH asymmetry needs its own source and
should be raised as a follow-up question against the column ticket.

---

## 6. The provider interface surface

### 6.1 What the LSER layer actually needs

To compute `k(φ, pH)` under LSER, retention is the population-weighted average over microspecies,
each of which is a *different solute* with its own Abraham descriptors:

```
k(φ, pH) = Σ_i  f_i( s_s pH(φ, T), φ )  ·  k_i(φ)
```

where `k_i` comes from LSER on microspecies *i*'s descriptors and the column/mobile-phase system
constants. Therefore:

- **Microspecies fractions are required.** Not a scalar α, not a "major species".
- **Per-species pKa values are not sufficient on their own**, but are needed to *derive* the
  fractions and to reason about crossovers, so both must be returned.
- **logD is NOT needed and must not be in the interface.** logD is octanol/water — the wrong phase
  system entirely. It is *derivable* from the same primitives (fractions × per-species logP, per the
  Starling/`unipka` construction) and so is strictly redundant. Its presence would invite the
  category error of substituting it for the LSER computation. **Expose logD, if at all, as a
  convenience method on the result object, never as a required provider capability.**

### 6.2 Proposed surface

```python
class IonisationProvider(Protocol):
    def speciate(
        self,
        smiles: str,
        ph: float,                 # s_s pH, at the working (phi, T)
        ph_scale: PhScale,         # REQUIRED, no default; enforce SS
        phi: float,                # organic volume fraction, 0..1
        modifier: Modifier,        # MEOH | ACN  -> selects correction coefficients
        temperature_c: float = 25.0,
    ) -> Speciation: ...
```

```python
@dataclass(frozen=True)
class Microspecies:
    smiles: str                    # the actual protomer/tautomer — LSER needs the structure
    charge: int                    # required: silanol/ion-exchange terms key off this
    fraction: float                # population at (ph, phi, T); fractions sum to 1
    fraction_sigma: float          # propagated, not optional

@dataclass(frozen=True)
class IonisableSite:
    atom_index: int
    pka_water: float               # w_w pKa, as predicted
    pka_mobile: float              # s_s pKa at (phi, modifier, T)
    sigma_water: float             # REQUIRED
    sigma_mobile: float            # sigma_water inflated by correction uncertainty
    charge_type: ChargeType        # NEUTRAL_ACID | CATIONIC_ACID
    family: Family | None          # ALIPHATIC_COOH | AROMATIC_COOH | PHENOL | AMINE |
                                   # PYRIDINE | UNKNOWN  -> selects Rosés coefficients

@dataclass(frozen=True)
class Speciation:
    species: tuple[Microspecies, ...]
    sites: tuple[IonisableSite, ...]
    provider: str
    provider_version: str          # pin: a spec'd prediction must be reproducible
    in_domain: bool                # OPERA-style applicability gate
    domain_score: float | None     # 0..1
    correction_extrapolated: bool  # True if phi > 0.60 or family is UNKNOWN
    notes: tuple[str, ...]
```

`speciate()` returns `None` — or raises a typed `ProviderDeclined` — for the ~1-in-15 compounds open
providers cannot process (§2, Uni-pKa 6.6% exception rate). **The declined path is part of the
contract, not an error case.**

### 6.3 Design notes

- **`family` is the new load-bearing field.** It is what selects the Rosés correction coefficients.
  No pKa model currently emits it, so the provider adapter must classify it — SMARTS is adequate for
  the five published families — and must be honest (`UNKNOWN`) when the molecule falls outside them,
  which sets `correction_extrapolated`.
- **`pka_water` and `pka_mobile` are both retained**, so the correction is auditable and can be
  swapped without touching the model.
- **`sigma_mobile > sigma_water` always.** The correction adds its own uncertainty; above 60% ACN, or
  for `family = UNKNOWN`, it should widen sharply. This is how the extrapolation risk in §5.2 reaches
  the calibrated-confidence output rather than being silently swallowed.
- **`smiles` per microspecies, not just a charge.** The LSER layer must run QSPR on the *actual
  protomer*, since the Abraham A and B descriptors differ drastically between a neutral acid and its
  anion. A provider that returns only pKa values cannot satisfy this and must be paired with an
  enumerator — which is precisely why Uni-pKa (which ships one) wins over QupKake and MolGpKa.
- **ChemAxon's adapter needs `cxcalc msdistr`**, since the Python API has no
  microspecies-distribution call (§2). The commercial adapter is therefore a *subprocess* adapter,
  not an in-process one. Worth knowing before #20 starts.

---

## 7. Recommendation

| | |
|---|---|
| **Open default** | **Uni-pKa** via `pip install unipka`, weights vendored and version-pinned |
| **Second open backend** | **QupKake** (BSD-3), for cross-provider disagreement as an uncertainty channel |
| **Domain gate only** | **OPERA** (MIT), used for `AD` / `Sim_index`, not for its pKa values |
| **Commercial** | **ChemAxon/Certara** via `pip install chemaxon` + `cxcalc msdistr` subprocess, licence permitting |
| **Watch** | **pKaLearn** (2026, CC-BY) — revisit after independent validation |
| **Mandatory layer we build** | Rosés/Bosch `w_w`pKa → `s_s`pKa(φ, modifier, family) correction; `w_w`pH → `s_s`pH scale conversion; conformal σ calibrated on **downstream log k** |

---

## 8. What I could NOT determine

1. **Per-method RMSE/MAE numbers** from the 2026 benchmark — Supplementary Tables S4–S6 not
   retrieved. Only the comparative statements in the main text were available.
2. **The Rosés/Bosch correction coefficients** (`a`, `b` per family per φ) — both source papers
   paywalled (403). **A hard prerequisite for #20.**
3. **Specific magnitudes of the pKa shift** at a given %ACN (e.g. "acids shift +1.7 at 50% ACN"). I
   established direction, family-dependence and functional form from primary abstracts, but no
   primary-source numeric table.
4. **Whether the ChemAxon academic licence includes Calculator Plugins / `cxcalc` / the pKa plugin.**
   Not stated on any current Certara page.
5. **ChemAxon academic licence lead time, cost, term, and whether the pre-acquisition AcPack terms
   survived.** No published turnaround anywhere.
6. **Whether an unaffiliated individual researcher qualifies** for the ChemAxon academic licence —
   explicitly not clarified.
7. **Whether Starling's weights/code are open.** The Rowan publication page does not say.
8. **The claim that silanol cation exchange is significant in ACN but not methanol.** The primary
   evidence I found documents it clearly *in methanol*; I found no source establishing the
   asymmetry. Needs its own citation.
9. **Uni-pKa's ensemble spread as a calibrated σ** — nobody has published a calibration study for it.
   We would be first, which is a reason to design the conformal layer to not depend on it alone.

---

## Appendix: all sources

**Benchmarks and models**
- Sipos-Szabó, Bajusz, Balogh, Keserű. *J. Chem. Inf. Model.* 66, 4607–4619 (2026). https://doi.org/10.1021/acs.jcim.6c00107 · https://pmc.ncbi.nlm.nih.gov/articles/PMC13126629/
- Zheng et al. Uni-pKa. *JACS Au* 4, 2467 (2024). https://doi.org/10.1021/jacsau.4c00271 · preprint https://doi.org/10.26434/chemrxiv-2023-lw5k0 · https://github.com/dptech-corp/Uni-pKa
- `unipka` wrapper: https://pypi.org/project/unipka/ · https://github.com/finlayiainmaclean/unipka
- Abarbanel & Hutchison. QupKake. *J. Chem. Theory Comput.* 20, 6946 (2024). https://doi.org/10.1021/acs.jctc.4c00328 · https://pmc.ncbi.nlm.nih.gov/articles/PMC11325546/ · https://github.com/hutchisonlab/QupKake
- Pan et al. MolGpKa. *J. Chem. Inf. Model.* 61, 3159 (2021). https://doi.org/10.1021/acs.jcim.1c00075 · https://github.com/Xundrug/MolGpKa
- Mansouri et al. OPERA pKa. *J. Cheminform.* 11, 60 (2019). https://doi.org/10.1186/s13321-019-0384-1 · https://pmc.ncbi.nlm.nih.gov/articles/PMC6749653/
- Mansouri et al. OPERA platform. *J. Cheminform.* 10, 10 (2018). https://doi.org/10.1186/s13321-018-0263-1 · https://github.com/kmansouri/OPERA · CLI help https://github.com/kmansouri/OPERA/blob/master/OPERA_Source_code/help.txt
- pyOPERA: https://cabreratoxy.github.io/pyOPERA/usage.html
- Genzling et al. pKaLearn. *Commun. Chem.* 9, 181 (2026). https://doi.org/10.1038/s42004-026-01983-y · https://github.com/MoitessierLab/pKaLearn
- pkasolver: https://github.com/mayrf/pkasolver
- Starling (Rowan): https://chemrxiv.org/engage/chemrxiv/article-details/68388349c1cb1ecda02ba65d · https://www.rowansci.com/publications/macroscopic-pka-prediction
- Zheng et al. pKa in non-aqueous solvents. *J. Comput. Chem.* 46, e27517 (2025). https://doi.org/10.1002/jcc.27517 · https://pmc.ncbi.nlm.nih.gov/articles/PMC11633825/

**ChemAxon / Certara**
- Certara completes acquisition of Chemaxon, 2 Oct 2024: https://www.certara.com/pressrelease/certara-completes-acquisition-of-chemaxon/
- Academic programme: https://www.certara.com/cxn/academic-program · Research licence: https://www.certara.com/cxn/research-license
- PyPI: https://pypi.org/project/chemaxon/
- Python API: https://apidocs.chemaxon.com/python_api/apidocs/chemaxon.html · calculations module https://apidocs.chemaxon.com/python_api/apidocs/chemaxon/calculations.html
- pKa plugin: https://docs.chemaxon.com/latest/calculators_pka-plugin.html · cxcalc functions https://docs.chemaxon.com/latest/cxcalc-calculator-functions.html · calculators user guide https://docs.chemaxon.com/latest/calculators_user-guide.html

**Retention, pH scales and solvent corrections**
- Espinosa, Bosch, Rosés. *J. Chromatogr. A* 964, 55 (2002). https://doi.org/10.1016/S0021-9673(02)00558-7
- Espinosa, Bosch, Rosés. *Anal. Chem.* 72, 5193 (2000). https://doi.org/10.1021/ac000591b
- Subirats, Rosés, Bosch. *Sep. Purif. Rev.* 36, 231–255 (2007). https://doi.org/10.1080/15422110701539129
- δ parameter, ACN/water, 0–90%, 15–60 °C. *Anal. Chem.* (2007). https://doi.org/10.1021/ac062372h
- Gradient pH change and retention. *Anal. Chem.* (2001). https://doi.org/10.1021/ac0101454
- pH and temperature retention model. *Anal. Chem.* (2007). https://doi.org/10.1021/ac060482i
- Phenol pKa in MeOH/water. *J. Chromatogr. A* (2000). https://doi.org/10.1016/S0021-9673(99)01139-5
- Silanol ion exchange, ESI-MS, 60% MeOH. *J. Chromatogr. A* (2001). https://doi.org/10.1016/S0021-9673(00)01192-4
- Silanol acidity across columns, 60% MeOH. *J. Chromatogr. A* (2003). https://doi.org/10.1016/S0021-9673(02)01899-X
