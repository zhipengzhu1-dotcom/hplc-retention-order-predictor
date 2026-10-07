# What LSER system constants are actually published?

Resolves #3. Blocks #11 (temperature), #12 (modifier), #13 (curated column set).

**Scope of this note.** Coverage inventory of published reversed-phase LSER *system* constants
(`c, e, s, a, b, v`) along four axes — column, organic modifier, mobile phase composition φ,
temperature — plus the ionisation question. Every claim is cited. A closing section states
what could **not** be determined and why.

---

## Headline answer

The load-bearing worry in the ticket — *"constants are typically reported at one or two fixed
compositions, often 50% methanol/water"* — **is false for the main body of the literature.**

The Wayne State University (WSU) system constant database, maintained by Colin F. Poole and
Sanka N. Atapattu, reports system constants **at 10% (v/v) increments from 10–70% organic
solvent**, i.e. a **seven-point φ series per column per modifier**, for tens of commercial
columns. The 2025 update covers:

| Modifier | Columns with a full φ series |
|---|---|
| methanol–water | **27** |
| acetonitrile–water | **25** |
| acetone–water | 7 |
| tetrahydrofuran–water | 3 |
| 2-propanol–water | 1 |

> "System constants are obtained for the range 10-70 % (v/v) organic solvent at 10 % (v/v)
> increments for methanol-water (27 columns), acetonitrile-water (25 columns), acetone-water
> (7 columns), tetrahydrofuran-water (3 columns), and 2-propanol-water (1 column). The Wayne
> State University 2025 (WSU-2025) compound descriptor database is used as the sole source for
> compound descriptors to unify the system constants with a verified descriptor database."

— Poole & Atapattu, *Update of the Wayne State University system constant database for
reversed-phase liquid chromatography columns for varied mobile phase compositions*,
J. Chromatogr. A 2025, **1762**, 466385. <https://doi.org/10.1016/j.chroma.2025.466385>
(PMID 40987226)

So the φ axis is **not** the binding constraint. The binding constraints turn out to be
**temperature** (essentially uncovered) and **ionisation** (covered at exactly one condition).

---

## Axis 1 — Column

### The two anchor publications

| Source | Columns | Notes |
|---|---|---|
| Poole, J. Chromatogr. A 2019, **1600**, 112–126. [doi:10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027) (PMID 31128882) | **25** siloxane-bonded type-B silica columns | MeOH + MeCN 10–70%; THF 10–70% for 2 columns |
| Poole & Atapattu, J. Chromatogr. A 2025, **1762**, 466385. [doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385) (PMID 40987226) | **27** MeOH / **25** MeCN | Supersedes 2019; recomputed against WSU-2025 descriptors |

The 2019 database was deliberately built to span morphology and topology:

> "Columns were selected to include examples of all common column packing morphologies for
> small molecule separations (totally porous particles, superficially porous particles,
> organic-inorganic hybrid particles, and a monolith) and common topologies
> octadecylsiloxane-bonded (including phases with a polar-embedded group, mixed-mode phases
> with a short-chain siloxane-bonded polar functional group, sterically crowded, positive
> shield, and fluorine-containing phases); octylsiloxane-bonded; and various phenyl-containing
> stationary phases with different linker arms, pentafluorophenylalkyl and biphenyl groups."

— Poole 2019, abstract (via [Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/J.CHROMA.2019.04.027), PMID 31128882)

### Named columns confirmed from open abstracts

I could not obtain the full 25/27-column table (see *Could not determine*). The following
columns are **individually named in abstracts I could read**, each with a published φ series.
This is a floor, not the full set — a confirmed **~17 columns**, already above the ticket's
committed scope of 5–15.

| Column | Chemistry | Vendor | Modifiers | φ range | Source |
|---|---|---|---|---|---|
| Kinetex C18 | C18, core-shell | Phenomenex | MeOH, MeCN | 10–70% | [doi:10.1016/j.chroma.2016.09.045](https://doi.org/10.1016/j.chroma.2016.09.045) |
| Kinetex C18 | ” | ” | acetone | 20–70% | [doi:10.1016/j.chroma.2021.462252](https://doi.org/10.1016/j.chroma.2021.462252) |
| Kinetex C18 | ” | ” | THF:IPA mixtures | 20–70% | [doi:10.1016/j.chroma.2023.463801](https://doi.org/10.1016/j.chroma.2023.463801) |
| Kinetex EVO C18 | ethyl-bridged C18, core-shell | Phenomenex | MeOH, MeCN | 10–70% | [doi:10.1007/s10337-017-3350-y](https://doi.org/10.1007/s10337-017-3350-y) |
| Kinetex XB-C18 | diisobutyl-C18 (sterically protected) | Phenomenex | MeOH, MeCN | 10–70% / 20–70% | [doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692) |
| Kinetex C8 | C8, core-shell | Phenomenex | MeOH, MeCN, acetone | 20–70% | [doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692), [doi:10.1016/j.chroma.2021.462252](https://doi.org/10.1016/j.chroma.2021.462252) |
| Kinetex Biphenyl | biphenylsiloxane, core-shell | Phenomenex | MeOH, MeCN | 10–70% | [doi:10.1016/j.chroma.2016.11.059](https://doi.org/10.1016/j.chroma.2016.11.059) |
| Kinetex Phenyl-Hexyl | phenylhexyl, core-shell | Phenomenex | MeOH, MeCN, acetone | 20–70% | [doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692) |
| Kinetex F5 | pentafluorophenylpropyl, core-shell | Phenomenex | MeOH, MeCN | 20–70% | [doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692) |
| Luna Omega PS C18 | electrostatic-shielded C18 | Phenomenex | acetone (+ MeOH/MeCN in db) | 20–70% | [doi:10.1016/j.chroma.2021.462252](https://doi.org/10.1016/j.chroma.2021.462252) |
| Luna C18 | C18, totally porous | Phenomenex | MeCN, acetone, MeOH, IPA, THF | 20–70% | Poole & Atapattu, J. Chromatogr. Open 2022 ([SSRN 4025190](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4025190)) |
| SunFire C18 | C18, totally porous | Waters | MeOH, MeCN, THF | 20–70% | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652) |
| XBridge Shield RP18 | polar-embedded C18, hybrid | Waters | MeOH, MeCN | 20–70% | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652) |
| XBridge C8 | C8, hybrid | Waters | MeOH, MeCN | 20–70% | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652) |
| XBridge Phenyl | phenyl, hybrid | Waters | MeOH, MeCN | 20–70% | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652) |
| Discovery HS F5 | pentafluorophenylpropyl | Supelco | MeOH, MeCN | 20–70% | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652), [doi:10.1016/j.chroma.2016.11.059](https://doi.org/10.1016/j.chroma.2016.11.059) |
| XTerra MS C18 | hybrid C18 | Waters | multiple | — | [doi:10.1365/s10337-005-0559-y](https://doi.org/10.1365/s10337-005-0559-y) |
| XTerra Phenyl | hybrid phenyl | Waters | MeOH | 70% (single point, in a correlation study) | [doi:10.1016/j.chroma.2025.466008](https://doi.org/10.1016/j.chroma.2025.466008) |
| Synergi Polar-RP | ether-linked phenylpropyl | Phenomenex | MeOH, MeCN | — | [doi:10.1016/j.chroma.2016.11.059](https://doi.org/10.1016/j.chroma.2016.11.059) |
| Synergi Hydro-RP | polar-endcapped C18 | Phenomenex | MeOH | 50% only, but **25–65 °C** | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652) |

**Verdict on axis 1: not a constraint.** Published coverage exceeds the committed 5–15 column
scope by 2–3×. The problem for #13 is *selection*, not *availability*.

Poole additionally provides a ready-made selection tool: the 2025 paper applies **hierarchical
cluster analysis** with the system constants as variables to group columns by selectivity, and
**system constant correlation plots** to identify selectivity-equivalent columns — exactly the
machinery a curated set needs to justify itself ([doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385)).

---

## Axis 2 — Organic modifier

Methanol and acetonitrile are **both** fully covered at parity (27 vs 25 columns, identical
10–70% grid). They are, as the ticket assumes, treated as genuinely different systems.

Beyond the two majors, coverage thins fast: acetone 7 columns, THF 3, 2-propanol 1
([doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385)).

**Ternary mobile phases are also covered**, which is more than expected:

- MeCN:MeOH (1:1)–water on a standard C18, with system maps compared against acetone, THF,
  MeCN, IPA and MeOH. Result: the 1:1 blend behaves **more like methanol than acetonitrile** —
  i.e. ternary system constants are *not* a linear blend of the two binaries.
  — Atapattu, J. Sep. Sci. 2023, **46**, e2300489. [doi:10.1002/jssc.202300489](https://doi.org/10.1002/jssc.202300489)
- THF:2-propanol–water ternaries at three ratios, 20–70% total organic, on Kinetex C18.
  — Atapattu, J. Chromatogr. A 2023, **1690**, 463801. [doi:10.1016/j.chroma.2023.463801](https://doi.org/10.1016/j.chroma.2023.463801)

The 2025 update explicitly concludes ternary solvent systems are a useful method-development
axis ([doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385)).

**Verdict on axis 2: not a constraint for MeOH + MeCN.** Anything else is a cliff.

---

## Axis 3 — Mobile phase composition φ

This is the ticket's central worry, and it is **substantially better than feared**.

### The grid

10–70% (v/v) in 10% steps = **7 points per (column, modifier)**. Some later selectivity papers
narrow to 20–70% (6 points) because high-bonding-density C18 phases give unmeasurably large
retention below 30% methanol:

> "Some octadecylsiloxane-bonded silica stationary phases with a high bonding density and
> methanol-water mobile phase compositions containing ≤ 30% (v/v) methanol exhibit extreme
> retention factors (log k > 2.5) for the low-polarity, two-ring aromatic compounds in the
> thirty-five compound calibration set."

— Poole, J. Chromatogr. A 2020, **1633**, 461652. [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652)

### Is there published guidance on the *functional form* of `c,e,s,a,b,v` vs φ?

**Partly, and it is the weakest link in this note.**

- Poole's answer is the **"system map"** — a simultaneous plot of all system constants as a
  *continuous function* of mobile phase composition, used for method development. The 2017
  review states system maps "provide a comprehensive overview of the separation system as a
  function of mobile phase composition and/or temperature for method development."
  — Poole & Lenca, J. Chromatogr. A 2017, **1486**, 2–19. [doi:10.1016/j.chroma.2016.05.099](https://doi.org/10.1016/j.chroma.2016.05.099)
  I could **not** confirm from open sources whether a fitted analytic form (quadratic,
  cubic spline, etc.) is published alongside, or whether the maps are drawn as smoothed
  interpolations of the seven points. **This is an open item — see below.**

- The **alternative, and the one the blueprint paper takes**, is to reparameterise so that φ
  disappears from the system constants. Combine LSER with linear solvent strength (LSS) theory:

  ```
  eq 1   k       = k_w · 10^(−S_S·φ)
  eq 2   log k   = c    + e·E    + s·S    + a·A    + b·B    + v·V        (φ-dependent)
  eq 4   log k_w = c_kw + e_kw·E + s_kw·S + a_kw·A + b_kw·B + v_kw·V     (φ-independent)
  eq 5   S_S     = c_SS + e_SS·E + s_SS·S + a_SS·A + b_SS·B + v_SS·V     (φ-independent)
  ```

  > "while the system parameters featuring in [eq 2] are functions of the volume fraction of
  > organic modifier (that is, they depend on ϕ and need recalibration if ϕ changes), the
  > system parameters in [eqs 4 and 5] are not."

  — Marchetto et al., Anal. Chem. 2025, **97**, 6991–7001.
  [doi:10.1021/acs.analchem.4c03466](https://doi.org/10.1021/acs.analchem.4c03466) ·
  [PMC11983366](https://pmc.ncbi.nlm.nih.gov/articles/PMC11983366/) (open access)

  This trades φ-interpolation error for **LSS linearity error**, and the paper is candid that
  the trade is not free:
  - "the LSS theory fails to capture a common nonlinear increase in log k as ϕ → 0"
  - linearity was "valid strictly for 30% < φ < 60% only"

  Their error cascade (Table 3, MAPE on retention factor):

  | Layer stack | MAPE |
  |---|---|
  | LSS theory alone (experimental `k_w`, `S_S`) | 9.1% |
  | LSS + LSER (experimental Abraham descriptors) | 17.8% |
  | LSS + LSER + QSPR (descriptors from SMILES) | 24.6% |

  Overall R² = 0.94. **The LSER layer costs ~8.7 pp of MAPE; the QSPR layer ~6.8 pp; the LSS
  linearity assumption is the single largest single-layer contributor at 9.1 pp baseline.**

- Poole's own analysis of LSS is directly relevant and **cautions against exactly this
  reparameterisation as a physical description**:

  > "It is shown that the linear region of the plots of the retention factor (log k) against
  > the volume fraction of organic solvent (ϕ) is largely a system property independent of
  > solute type… Log kw cannot be recommended as a descriptor of solute properties since it
  > has no clear connection to a real distribution system."

  Column-specific solvent strength parameters *can* themselves be predicted by the solvation
  parameter model, "with a typical deviation of about 0.12 over a range of 1.69 to 6.33".
  — Poole & Atapattu, J. Chromatogr. A 2022, **1675**, 463153. [doi:10.1016/j.chroma.2022.463153](https://doi.org/10.1016/j.chroma.2022.463153)

  Note that S_S is **column dependent**, not a universal solvent constant: MeOH 3.12 ± 0.12,
  MeCN 2.78 ± 0.18, acetone 2.71 ± 0.11, THF 2.95 ± 0.24 as *averages*, with statistical
  analysis confirming column dependence. A simulator that hard-codes a per-solvent S is wrong.

- Poole also warns that applying the solvation parameter model to **gradient** elution is
  approximate:

  > "The solvation parameter model has been applied to gradient elution separations but here
  > theory and practice suggest a cautious approach since the interpretation of system and
  > compound properties derived from its use are approximate."

  — [doi:10.1016/j.chroma.2016.05.099](https://doi.org/10.1016/j.chroma.2016.05.099)

**Verdict on axis 3: the φ risk is real but far smaller than the ticket assumed.** We have 7
points, not 1. Two competing strategies exist (interpolate the system map vs. reparameterise
into φ-independent `k_w`/`S_S` LSERs), and there is a quantified error budget for the second.
There is **no published quantification of the interpolation error for the first** — that
remains a genuine gap and a candidate prototype.

---

## Axis 4 — Temperature

**This is the real gap.** The RPLC system constant database is effectively isothermal.

- The 2019 and 2025 WSU RPLC databases describe **no** temperature axis in their abstracts —
  the reported variables are column, modifier and φ only
  ([2019](https://doi.org/10.1016/j.chroma.2019.04.027), [2025](https://doi.org/10.1016/j.chroma.2025.466385)).
- The one RPLC temperature series I could confirm is a **single column at a single composition**:
  > "System maps for XBridge Shield RP18 for 20-70% (v/v) methanol-water **and Synergi Hydro-RP
  > and 50% (v/v) methanol-water at temperatures from 25-65 °C**…"

  — Poole, J. Chromatogr. A 2020, **1633**, 461652. [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652)
- One independent study builds a joint `T × φ` LSER model, again on a **single ODS column with
  methanol only**:
  > "The prediction capability of the solvation parameter model in reverse-phase liquid
  > chromatography at different methanol-water mobile phase compositions and temperatures was
  > investigated… The coefficients obtained in the regressions were used to create a general
  > retention model able to predict retention in an octadecylsilica stationary phase at any
  > temperature and methanol-water composition."

  Validated on a 30-solute external test set; the authors give both a coefficient-level linear
  relation in (φ, T) and a single general equation taking φ and T directly.
  — Gotta, Keunchkarian, Castells & Reta, J. Sep. Sci. 2012, **35**, 2699–2709.
  [doi:10.1002/jssc.201200197](https://doi.org/10.1002/jssc.201200197) (PMID 22997100)
- Rosés, Subirats & Bosch review models extending RPLC retention to composition, pH **and**
  temperature jointly, including the solvation parameter model — the best entry point for the
  temperature+pH literature.
  — J. Chromatogr. A 2009, **1216**, 1756–1775. [doi:10.1016/j.chroma.2008.12.042](https://doi.org/10.1016/j.chroma.2008.12.042)

**Contrast with gas chromatography, where Poole *does* publish a temperature grid:** system
constants at **20 °C intervals over 60–140 °C for 61 open-tubular GC columns**, with system
maps drawn as "a continuous function of temperature"
([doi:10.1016/j.chroma.2025.466595](https://doi.org/10.1016/j.chroma.2025.466595), PMID 41386122).
The technique and the methodology exist. They have simply not been applied to build an
equivalent RPLC temperature database.

**Verdict on axis 4: a hard constraint.** There is no published `(column × modifier × φ × T)`
system constant grid for RPLC. Any temperature axis in our simulator must come from (a) fixing
T and treating it as a method constant, (b) a van 't Hoff / Gotta-style correction fitted from
our own or RepoRT data, or (c) new measurement. Note also that the RPLC databases do not state
their operating temperature in the abstracts — **we do not currently know what T the WSU
constants are "at"**, which is itself blocking (see below).

---

## The ionisation question

**LSER system constants in the WSU database are for neutral compounds only.** This is explicit
and repeated:

- "The properties of the separation system are described by five system constants representing
  all possible intermolecular interactions **for neutral molecules**."
  ([doi:10.1016/j.chroma.2016.05.099](https://doi.org/10.1016/j.chroma.2016.05.099))
- Poole's tutorial is titled *"…Tutorial on its application to separation systems **for neutral
  compounds**"* ([doi:10.1016/j.chroma.2021.462108](https://doi.org/10.1016/j.chroma.2021.462108)).
- The WSU-2025 descriptor database characterises "the capability of **non-ionic** compounds to
  interact with their environment" ([doi:10.1016/j.chroma.2025.465958](https://doi.org/10.1016/j.chroma.2025.465958)).
- The 2019 database explicitly names cation-exchange interaction as a mechanism **outside** the
  model: "Contributions to retention from steric resistance and cation-exchange interactions
  **not parameterized in the solvation parameter model** are identified for columns and
  conditions where they may be important."
  ([doi:10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027))
- And the effect is **modifier-dependent**: "Electrostatic interactions (cation-exchange) are
  important for the retention of weak bases with **acetonitrile–water but not methanol–water**
  mobile phase compositions" (Kinetex EVO C18 system map paper,
  [doi:10.1007/s10337-017-3350-y](https://doi.org/10.1007/s10337-017-3350-y)).

### But ionic-form constants *do* exist — at exactly one condition

Abraham & Acree established descriptors for single ions and ionic species (carboxylate anions,
phenoxide anions, protonated base cations), with prediction equations, and applied them to HPLC
retention among other processes.
— Abraham & Acree, *Descriptors for ions and ion-pairs for use in linear free energy
relationships*, J. Chromatogr. A 2016, **1430**, 2–14.
[doi:10.1016/j.chroma.2015.07.023](https://doi.org/10.1016/j.chroma.2015.07.023) (PMID 26189671)
(see Poole's editorial, [doi:10.1016/j.chroma.2015.11.051](https://doi.org/10.1016/j.chroma.2015.11.051))

The directly usable RPLC result is Rosés' group:

> "The LFER model of Abraham is applied to the retention of the neutral and ionic forms of
> **94 solutes in a C18 column and 40% v/v acetonitrile/water mobile phase**… A new LFER model
> for application to the retention of partially dissociated acids and bases is derived
> **averaging the descriptors of the neutral and ionic forms according to their degrees of
> ionization** in the mobile phase. This new LFER model is satisfactorily compared to other
> literature modified Abraham models for a set of 498 retention data of partially dissociated
> acids and bases."

— Soriano-Meseguer, Fuguet, Abraham, Port & Rosés, J. Chromatogr. A 2021, **1635**, 461720.
[doi:10.1016/j.chroma.2020.461720](https://doi.org/10.1016/j.chroma.2020.461720) (PMID 33234293)

**One column. One modifier. One composition (φ = 0.40 MeCN).** That is the whole ionised-species
coverage. The ticket's "single-point" fear is entirely accurate — it just applies to the
*ionisation* axis rather than the φ axis.

Two consequences the same paper hands us, both first-class for the ionisation provider ticket:

1. The partially-ionised model works by **α-weighted averaging of neutral and ionic descriptors**,
   where α is the ionisation degree. That is a concrete interface for our `pKa` provider: it
   must emit α, not just a `pKa`.
2. **α must be computed in the mobile phase, not in water:** "Calculation of the ionization
   degrees in the chromatographic mobile phase (i.e. from pH and pKa in the eluent) give good
   correlations for all tested models. However, estimation of these ionization degrees from
   pH − pKa data in **pure water gives biased estimations** of the retention of the partially
   ionized solutes."

   Supporting mechanism from the same group: `pʷˢKa` correlates linearly with `pʷʷKa` but with
   **two different correlations, one for neutral acids and one for neutral bases**; anions are
   less retained than cations of equal lipophilicity; and the hold-up time for anions is itself
   pH-dependent (0.83 min at acidic pH falling to 0.65 min at pH 11) due to silanol ionisation
   and charge repulsion. — Soriano-Meseguer et al., Anal. Chim. Acta 2019, **1078**, 200–211.
   [doi:10.1016/j.aca.2019.05.063](https://doi.org/10.1016/j.aca.2019.05.063) (PMID 31358220)

   That last point is a quiet landmine: for anions, `t₀` is not a constant, so `k` itself is
   ill-defined without a pH-dependent hold-up model.

---

## Uncertainty: what is actually reported with the constants?

The ticket asks specifically. The answer is **good** — this literature is unusually disciplined
about fit quality.

Routinely reported per model (i.e. per column × modifier × φ):

- **coefficient of determination R²**
- **Fisher statistic F**
- **standard error of the estimate (SE)**, in log k units
- **standard deviation of the individual system constants** — i.e. per-coefficient uncertainty
  on `c, e, s, a, b, v` themselves

Representative published values:

| System | R² | F | SE (log k) | Source |
|---|---|---|---|---|
| 7 core-shell columns, MeOH/MeCN, 20–70%, n=35 calibration set | 0.970–0.999 (<5% below 0.99, 84 models) | — | 0.011–0.057 (75% < 0.030) | [doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692) |
| SunFire/XBridge/Discovery HS F5, MeOH/MeCN/THF | avg 0.996 (sd 0.003, n=11) | — | avg 0.025 (sd 0.005) | [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652) |
| 5 columns, acetone–water 20–70% | 0.988–0.998 | 277–1551 | 0.024–0.097 | [doi:10.1016/j.chroma.2021.462252](https://doi.org/10.1016/j.chroma.2021.462252) |
| MeCN:MeOH (1:1) ternary, C18 | 0.998–0.999 | 1687–4015 | 0.022–0.029 | [doi:10.1002/jssc.202300489](https://doi.org/10.1002/jssc.202300489) |
| THF:IPA ternaries, Kinetex C18 | 0.991–0.999 | 338–1850 | 0.024–0.097 | [doi:10.1016/j.chroma.2023.463801](https://doi.org/10.1016/j.chroma.2023.463801) |

External-test-set validation is also reported: RMSEP(log k) ≈ 0.028–0.032, average absolute
error ≈ 0.026 ([doi:10.1016/j.chroma.2020.461692](https://doi.org/10.1016/j.chroma.2020.461692),
[doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652)).

**Read this as a floor on our achievable error.** The LSER layer, given *experimental* solute
descriptors and *published* system constants at a *measured* composition, is good to ~0.03 log k
units. Everything above that in our stack is our own error.

Poole's 2021 tutorial is the methods reference for how these are determined and what statistics
to demand: "Suitable experimental protocols to determine system constants by multiple linear
regression analysis and descriptors by the Solver method are presented; statistical tools to
evaluate model quality are discussed; and model-specific data analysis tools based on system
maps and correlation diagrams are described."
— [doi:10.1016/j.chroma.2021.462108](https://doi.org/10.1016/j.chroma.2021.462108) (PMID 33857674)

---

## A trap: descriptor databases are not interchangeable

This is not in the ticket but it invalidates an obvious implementation shortcut, so it belongs here.

System constants are regression coefficients fitted **against a specific solute descriptor
database**. Two major databases exist and they disagree:

- **WSU descriptor database** — several hundred compounds; WSU-2025 has **387** compounds
  ([doi:10.1016/j.chroma.2025.465958](https://doi.org/10.1016/j.chroma.2025.465958)).
- **Abraham descriptor database** — several thousand compounds. This is the lineage behind the
  Acree/UFZ-LSER collections listed in `sources/README.md` (~8,000 compounds).

> "These publicly accessible databases were developed independently using different approaches
> and for many compounds provide different descriptor values… **It is shown that the two
> descriptor databases are not interchangeable.** The WSU descriptor database consistently
> demonstrates improved model quality… **Model system constants exhibit a general dependence on
> database selection** with an approximately linear trend as a function of the fraction of
> compounds assigned descriptors from either database."

— Poole, J. Chromatogr. A 2023, **1692**, 463851. [doi:10.1016/j.chroma.2023.463851](https://doi.org/10.1016/j.chroma.2023.463851)

Poole gives one piece of slack: "no real cause for concern for relatively large datasets
containing < 15% of compounds with descriptors assigned from the other [database]".

**Implication for our architecture.** The WSU-2025 system constants were regressed against
WSU-2025 descriptors. Our QSPR layer will almost certainly be trained on the Abraham/UFZ/Acree
descriptor set, because it is 20× larger. **Feeding UFZ-lineage predicted descriptors into
WSU-fitted system constants introduces a systematic bias that is not in anyone's error budget.**
This needs an explicit decision and probably a small calibration study. It is arguably a new
research ticket.

Related: the **UFZ-LSER database is a *solute descriptor* resource, not a chromatographic
*system constant* resource.** It provides Abraham parameters for >7,000 small molecules
(<https://www.ufz.de/lserd>). It does **not** supply `c,e,s,a,b,v` for RPLC columns. Do not
plan on it for the column axis — that comes from Poole/WSU alone.

---

## What this implies for the blocked decisions

### #12 — Organic modifier

**Commit to methanol/water and acetonitrile/water, both, as first-class independent systems.**
Coverage is at parity (27 vs 25 columns, identical 10–70% grid). Modelling only one would
discard half a free dataset for no gain. Do **not** extend to acetone/THF/IPA — coverage
collapses to 7/3/1 columns. Ternaries are documented but explicitly non-additive
([doi:10.1002/jssc.202300489](https://doi.org/10.1002/jssc.202300489)); treat as out of scope
unless a specific ticket demands them.

One caveat that couples modifier to ionisation: cation-exchange retention of weak bases is
significant **in MeCN but not MeOH** ([doi:10.1007/s10337-017-3350-y](https://doi.org/10.1007/s10337-017-3350-y)).
For basic analytes, MeOH systems will be better behaved under a pure-LSER model. Expect
modifier-dependent residual structure for bases.

### #11 — Temperature

**The literature does not support temperature as a modelled LSER axis.** There is no published
`(column × modifier × φ × T)` system constant grid for RPLC — the entire confirmed coverage is
one column at one composition over 25–65 °C, plus one independent single-column joint `T × φ`
model. The recommendation this evidence supports is: **fix temperature as a method constant in
v1**, and if temperature must vary, model it as a *correction layer* on top of isothermal LSER
(Gotta et al.'s coefficient-vs-T linear relations,
[doi:10.1002/jssc.201200197](https://doi.org/10.1002/jssc.201200197), are the cited-and-move-on
option under the cut rule), not as an axis of the system constant table.

Blocking sub-question: **we must establish what temperature the WSU constants were measured at**
before we can even declare our fixed T. See below.

### #13 — Curated column set

**Unblocked and comfortable.** ~17 columns confirmed by name from open abstracts alone, 25–27 in
the full database, against a committed scope of 5–15. The selection criterion should come from
Poole's own hierarchical cluster analysis and system-constant correlation plots
([doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385)) — pick one
representative per selectivity cluster, which is precisely the "orthogonal column set" a method
development system wants. Prefer the WSU-2025 numbers over WSU-2019: same columns, recomputed
against a better descriptor database.

Constraint to carry forward: pick columns that appear in **both** the MeOH and MeCN tables, and
that survive the 20–70% narrowing (high-bonding-density C18s lose the 10–20% region for
low-polarity solutes, [doi:10.1016/j.chroma.2020.461652](https://doi.org/10.1016/j.chroma.2020.461652)).

---

## Could not determine

Stated explicitly, as required.

1. **The full 25/27-column list, and the system constant values themselves.** Both anchor papers
   are paywalled (J. Chromatogr. A; Semantic Scholar reports `openAccessPdf: CLOSED` for
   [10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027)). ScienceDirect,
   Springer and ACS all returned HTTP 403 to automated fetching. The ~17 named columns above are
   reconstructed from open abstracts and are a **floor**. **Action: obtain the PDFs of
   [doi:10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385) and
   [doi:10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027).** Nothing
   downstream can be numerically specified without them.

2. **The temperature at which the WSU RPLC system constants were measured.** Not stated in any
   abstract I could read. This is a genuine blocker on #11 — we cannot say "our model is at T₀"
   without knowing T₀. Likely 25 °C or 30 °C by convention, but *I did not verify this and it
   must not be assumed.*

3. **The buffer / pH conditions of the WSU database.** Not stated in abstracts. Since the
   calibration sets are neutral compounds this may be unbuffered water, but for our ionisable
   scope the eluent composition matters. Unverified.

4. **Whether a fitted analytic form for `c,e,s,a,b,v` vs φ is published.** System maps are
   described as continuous functions of composition, but I found no open statement of the
   fitting function (polynomial order, spline, or none). Targeted searches for
   "quadratic"/"second-order polynomial" returned nothing. **This is the single most important
   remaining unknown for the φ-interpolation risk**, and it is answered by the same PDFs as (1).

5. **Whether the WSU system constant tables are distributed as machine-readable supplementary
   data**, or exist only as typeset tables in the papers. Affects ingestion effort materially.

6. **Temperature in the Marchetto blueprint paper** — not stated in the open-access text
   ([PMC11983366](https://pmc.ncbi.nlm.nih.gov/articles/PMC11983366/)). Their data came from
   Poole by personal communication (Kinetex XB-C18, MeCN/water, 48 solutes, φ = 0.10–0.70),
   so it inherits whatever (2) turns out to be.

7. **Whether ionic-form LFER coefficients exist at any composition other than 40% MeCN, or on
   any column other than the Rosés C18.** I found none, but absence of evidence here is weaker
   than for the other axes — I searched the Rosés group's output and the Abraham/Acree ion
   descriptor review, not exhaustively.

---

## New questions this raises

- **Descriptor-database coupling (likely a new ticket).** WSU-fitted system constants + UFZ/Acree
  QSPR-predicted descriptors is a mismatch with documented, systematic bias
  ([doi:10.1016/j.chroma.2023.463851](https://doi.org/10.1016/j.chroma.2023.463851)). Do we
  train QSPR on WSU-2025 (387 compounds, tiny) or accept and calibrate out the bias?
- **φ-interpolation vs LSS reparameterisation is now a real design fork**, with a partial error
  budget for one branch (24.6% MAPE end-to-end,
  [doi:10.1021/acs.analchem.4c03466](https://doi.org/10.1021/acs.analchem.4c03466)) and none for
  the other. Prototype candidate: fit both against the 7-point φ series and compare.
- **`S_S` is column-dependent, not a solvent constant**
  ([doi:10.1016/j.chroma.2022.463153](https://doi.org/10.1016/j.chroma.2022.463153)). Any LSS
  layer in our stack must carry a per-column `S`, and Poole shows `S` is itself LSER-predictable
  (typical deviation 0.12 over a 1.69–6.33 range).
- **Hold-up time is pH-dependent for anions** ([doi:10.1016/j.aca.2019.05.063](https://doi.org/10.1016/j.aca.2019.05.063)).
  `k = (t_R − t₀)/t₀` with a constant `t₀` is wrong for anionic analytes at high pH. This touches
  the instrument layer, not just the chemistry layer.
- **The pKa provider's contract should be α (ionisation degree) in the mobile phase, not pKa in
  water** ([doi:10.1016/j.chroma.2020.461720](https://doi.org/10.1016/j.chroma.2020.461720)).
  Water-based α gives biased retention. This sharpens the ionisation-provider ticket.
- **Steric resistance and cation exchange are named, unparameterised residual mechanisms** in the
  LSER framework ([doi:10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027)).
  The map's "secondary retention mechanisms" fog item is now concrete, and Poole states flatly
  that HSM and LSER parameters "have little overlap" and cannot be quantitatively bridged —
  which is a direct, negative-looking input to the HSM→LSER bridge ticket.

---

## Source index

Primary, in order of load-bearing weight:

1. Poole & Atapattu, J. Chromatogr. A **1762** (2025) 466385 — WSU-2025 system constant database. [doi](https://doi.org/10.1016/j.chroma.2025.466385)
2. Poole, J. Chromatogr. A **1600** (2019) 112–126 — 25-column database. [doi](https://doi.org/10.1016/j.chroma.2019.04.027)
3. Poole & Lenca, J. Chromatogr. A **1486** (2017) 2–19 — applications review, system maps, neutrality, gradient caution. [doi](https://doi.org/10.1016/j.chroma.2016.05.099)
4. Poole, J. Chromatogr. A **1645** (2021) 462108 — methods tutorial. [doi](https://doi.org/10.1016/j.chroma.2021.462108)
5. Poole & Atapattu, J. Chromatogr. A **1675** (2022) 463153 — LSS solvent strength analysis. [doi](https://doi.org/10.1016/j.chroma.2022.463153)
6. Soriano-Meseguer, Fuguet, Abraham, Port & Rosés, J. Chromatogr. A **1635** (2021) 461720 — partially ionised LFER. [doi](https://doi.org/10.1016/j.chroma.2020.461720)
7. Abraham & Acree, J. Chromatogr. A **1430** (2016) 2–14 — ion descriptors. [doi](https://doi.org/10.1016/j.chroma.2015.07.023)
8. Poole, J. Chromatogr. A **1692** (2023) 463851 — descriptor database non-interchangeability. [doi](https://doi.org/10.1016/j.chroma.2023.463851)
9. Poole & Atapattu, J. Chromatogr. A **1634** (2020) 461692 — core-shell selectivity, calibration statistics. [doi](https://doi.org/10.1016/j.chroma.2020.461692)
10. Poole, J. Chromatogr. A **1633** (2020) 461652 — calibration compound selection; the one RPLC temperature series. [doi](https://doi.org/10.1016/j.chroma.2020.461652)
11. Marchetto, Tirapelle, Mazzei, Sorensen & Besenhard, Anal. Chem. **97** (2025) 6991–7001 — the blueprint; LSER+LSS reparameterisation and error cascade. [doi](https://doi.org/10.1021/acs.analchem.4c03466) · [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11983366/)
12. Rosés, Subirats & Bosch, J. Chromatogr. A **1216** (2009) 1756–1775 — ionisable + φ + T review. [doi](https://doi.org/10.1016/j.chroma.2008.12.042)
13. Gotta, Keunchkarian, Castells & Reta, J. Sep. Sci. **35** (2012) 2699–2709 — joint φ and T LSER model. [doi](https://doi.org/10.1002/jssc.201200197)
14. Soriano-Meseguer et al., Anal. Chim. Acta **1078** (2019) 200–211 — ionisation, pKa in eluent, pH-dependent t₀. [doi](https://doi.org/10.1016/j.aca.2019.05.063)
15. Poole, J. Chromatogr. A **1766** (2026) 466595 — GC temperature-grid database, as the contrast case. [doi](https://doi.org/10.1016/j.chroma.2025.466595)
16. Poole, J. Chromatogr. A **1752** (2025) 465958 — WSU-2025 descriptor database (387 compounds). [doi](https://doi.org/10.1016/j.chroma.2025.465958)
17. Atapattu, J. Sep. Sci. **46** (2023) e2300489 — MeCN:MeOH ternary non-additivity. [doi](https://doi.org/10.1002/jssc.202300489)
18. Atapattu, J. Chromatogr. A **1650** (2021) 462252 — acetone–water. [doi](https://doi.org/10.1016/j.chroma.2021.462252)
19. Atapattu, J. Chromatogr. A **1690** (2023) 463801 — THF:IPA ternaries. [doi](https://doi.org/10.1016/j.chroma.2023.463801)
20. Atapattu, Poole & Praseuth, J. Chromatogr. A **1468** (2016) 250–256 and **1478** (2016) 68–74 — Kinetex C18 and Biphenyl system maps. [doi](https://doi.org/10.1016/j.chroma.2016.09.045) · [doi](https://doi.org/10.1016/j.chroma.2016.11.059)
21. UFZ-LSER database (solute descriptors only, not system constants). <https://www.ufz.de/lserd>

Abstracts were retrieved via the NCBI E-utilities API (`efetch`, db=pubmed), which returns
publisher-deposited abstracts verbatim; PMIDs are given inline throughout. Full texts of the
paywalled anchor papers were not accessible — see *Could not determine*.
