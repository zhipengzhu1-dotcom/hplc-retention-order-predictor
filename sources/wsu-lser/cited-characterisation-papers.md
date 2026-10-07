# WSU-2019 refs 29–42 — the per-column characterisation papers, resolved

Extracted 2026-08-16 from the Poole 2019 PDF already held locally
(`papers/poole-2019-rplc-system-constants/`), via `pdftotext`. **No new acquisition was
needed to produce this list** — the citations were always in a document we had.

These are the papers the `reference` column of
[`wsu2019-column-caveats.csv`](wsu2019-column-caveats.csv) points at. Poole's Table S-1
states *that* a column has a caveat; these papers state *which compounds* produced it, which
is what #43 needs for an
operational definition of "bulky".

## The list

| Ref | Citation | Columns it covers |
|---|---|---|
| **29** | C.F. Poole, H. Ahmed, W. Kiridena, C. DeKay, W.W. Koziol, *Contribution of **steric repulsion** to retention on an octadecylsiloxane-bonded silica stationary phase in reversed-phase liquid chromatography*, **Chromatographia 62 (2005) 553–561** | Ascentis C18 |
| **30** | C.F. Poole, H. Ahmed, W. Kiridena, C. DeKay, W.W. Koziol, *Insight into the retention mechanisms on perfluorohexylpropylsiloxane-bonded (**Fluophase-RP**) and octadecylsiloxane-bonded (**Betasil C18**) stationary phases **based on the same silica substrate*** | Betasil C18, Fluophase-RP |
| **31** | W. Kiridena, C. DeKay, W.W. Koziol, Z. Ali, H. Ahmed, C.F. Poole, *Insight into the retention mechanism on a pentafluorophenylpropylsiloxane-bonded silica stationary phase (Discovery HS F5) in RP-LC*, **Chromatographia 63 (2006) 407–417** | Discovery HS C18, Discovery HS F5 |
| **32** | C.F. Poole, W. Kiridena, C. DeKay, W.W. Koziol, R.D. Rosencrans, *Insights into the retention mechanism on an octadecylsiloxane-bonded silica stationary phase (HyPURITY C18)*, **J. Chromatogr. A 1115 (2006) 133–141** | HyPURITY C18 |
| 33 | S.N. Atapattu, C.F. Poole, M.B. Praseuth, *System maps … superficially porous ethyl-bridged octadecylsiloxane-bonded silica*, **Chromatographia 80 (2017) 1279–1286** | Kinetex EVO C18 |
| 34 | S.N. Atapattu, C.F. Poole, M.B. Praseuth, *Insights … octylsiloxane- and diisobutyloctadecylsiloxane-bonded silica*, **Chromatographia 81 (2018) 373–385** | Kinetex C8, Kinetex XB-C18 |
| 35 | S.N. Atapattu, C.F. Poole, M.B. Praseuth, *Insights … silica-based phenyl phases*, **Chromatographia 81 (2018) 225–238** | Kinetex F5, Kinetex Phenyl-Hexyl |
| 36 | S.N. Atapattu, C.F. Poole, M.B. Praseuth, *System maps … biphenylsiloxane-bonded silica*, **J. Chromatogr. A 1478 (2016) 68–74** | Kinetex Biphenyl |
| 37 | S.N. Atapattu, K.R.D. Johnson, C.F. Poole, *Insights … electrostatic-shielded octadecylsiloxane-bonded silica*, **Chromatographia (2019), in press** — DOI printed as `10.1007/s10337-0…` (truncated in the PDF) | Luna Omega PS C18 |
| **38** | W. Kiridena, S.N. Atapattu, C.F. Poole, W.W. Koziol, *System maps for RP-LC on an octadecylsiloxane-bonded silica stationary phase (**SunFire C18**)*, **Chromatographia 68 (2008) 11–17** | SunFire C18 |
| **39** | S.N. Atapattu, C.F. Poole, *Factors affecting the interpretation of selectivity on Synergi reversed-phase columns*, **Chromatographia 71 (2010) 185–193** | Synergi Fusion-RP, Hydro-RP, Polar-RP |
| **40** | W. Kiridena, C.F. Poole, S.N. Akapattu, J. Qian, W.W. Koziol, *Comparison of the separation characteristics of the organic–inorganic hybrid octadecyl stationary phases **XTerra MS C18** and **XBridge C18** and **Shield RP18***, **Chromatographia 66 (2007) 453–460** | XBridge C18, XBridge Shield RP18, XTerra MS C18 |
| **41** | W. Kiridena, C. DeKay, N.D. Villiere, W.W. Koziol, C.F. Poole, *System maps for XTerra MS C18: effect of solvent type on selectivity*, **Chromatographia 61 (2005) 587–593** | XTerra MS C18 |
| 42 | W. Kiridena, S.N. Atapattu, C.F. Poole, W.W. Koziol, *Comparison … XBridge C8 and phenyl and XTerra phenyl*, **Chromatographia 68 (2008) 491–500** | XBridge C8, XBridge Phenyl, XTerra Phenyl |

Bold rows carry a **steric_repulsion** flag and are the ones #43 needs.

## #43 needs six papers, not fourteen

Cross-referencing the caveat register: only **8 of the 354 fits** carry a `steric_repulsion`
flag, across 8 columns, and those columns cite only **refs 29, 30, 31, 32, 39, 40 and 41**.

| Ref | Journal | Steric-flagged columns it explains |
|---|---|---|
| 29 | Chromatographia 62 (2005) 553–561 | Ascentis C18 |
| 30 | Chromatographia | Betasil C18, Fluophase-RP |
| 31 | Chromatographia 63 (2006) 407–417 | Discovery HS C18, Discovery HS F5 |
| 32 | J. Chromatogr. A 1115 (2006) 133–141 | HyPURITY C18 |
| 39 | Chromatographia 71 (2010) 185–193 | Synergi Fusion-RP |
| 40 + 41 | Chromatographia 66 (2007) 453–460; 61 (2005) 587–593 | XTerra MS C18 |

**Ref 29 is the priority.** Its title is literally *"Contribution of steric repulsion to
retention…"* — it is the methodological source for the whole category, not merely one
column's note.

## Ref 30 is the right source for #40's unresolved conflict, and its design is ideal

#40 logged an unadjudicated conflict: WSU-2019 Table S-1 flags **Fluophase-RP as
electrostatic in methanol as well as acetonitrile**, against the pinned convention (Poole
1600 §3.4) that cation exchange is acetonitrile-specific. One proposed resolution was that
the convention may be **C18-scoped**, Fluophase-RP being a PFP-type phase.

Ref 30's full title settles that it is the right paper and that its design can answer the
question: it compares **Fluophase-RP against Betasil C18 on the *same silica substrate***.
That is the controlled comparison — same silanols, different bonded ligand — which is
exactly what separates "the convention is C18-scoped" from "the S-1 note covers a different
mechanism".

⚠ Neither ticket that held this conflict is open: #39 and #40 are both closed. **If ref 30 is
obtained there is currently no open ticket to receive the answer.**

## Acquisition notes

- **Nine of the fourteen are *Chromatographia*** (Springer), the rest *J. Chromatogr. A*
  (Elsevier). Both are paywalled; neither is on the route that has worked so far.
- Volume and page numbers are recovered for all but refs 30 and 37. Ref 30's volume/pages
  were split across a column break in the PDF and were not recoverable from the text layer;
  ref 37 was *in press* at publication, so its DOI is the only identifier and the PDF prints
  it truncated.
- These are **2005–2019 characterisation papers, all from Poole's group** (Poole is an author
  or co-author on every one). If a request is being made anyway, asking the corresponding
  author for refs 29–42 as a set is plausibly cheaper and faster than seven separate paywall
  transactions — and the same contact is already the natural route for WSU-2025, which is
  #22's top-priority item.
