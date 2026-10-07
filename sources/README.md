# Source material

The PDFs themselves are held locally and gitignored (copyright). This file records what
they are, where to get them, and what each gives the project.

**One folder per article, under [`papers/`](../papers/), each with its own README** giving
the citation and the role of every file in it. Nothing in `papers/` is tracked except those
READMEs — not the PDFs, and not text extracted from them. Extracted *data* goes to
`sources/` with a provenance README (`sources/wsu-lser/`, `sources/usp-column-db/`,
`sources/hsm-column-db/`, `sources/instrument-volumes/`). **They carry different licences** —
CC BY-NC-SA 3.0 US for the HSM set, no identified grant for the USP extract — so never
redistribute them under one statement. In this public copy the USP extract and the WSU-2019
tables are **not included**; their READMEs explain how to obtain or rebuild them.

**Instrument volumes** — dwell, extra-column and flow-cell volumes for common Waters and
Agilent systems, in [`instrument-volumes/`](instrument-volumes/README.md). Dwell spans
73–1190 µL across Waters alone and varies with pump type and mixer, so it is a **per-run
error** with a real width, never a vendor constant. Feeds the instrument layer #44 says owes
`width[S, n]` real content, and changes what #45 must recover from RepoRT 0415's source paper
(the instrument model, usually stated) rather than the dwell volume (rarely stated).
⚠ The Agilent rows are **unsourced** and carry open gaps.

## Primary papers

**Marchetto, Tirapelle, Mazzei, Sorensen, Besenhard — *In Silico High-Performance Liquid
Chromatography Method Development via Machine Learning*.** Anal. Chem. 2025, 97, 6991–7001.
<https://doi.org/10.1021/acs.analchem.4c03466> (Open Access, CC-BY 4.0)
Local folder: [`papers/marchetto-2025-in-silico-hplc-ml/`](../papers/marchetto-2025-in-silico-hplc-ml/README.md)

The blueprint for this effort. Chains three layers:

```
SMILES → molecular descriptors → QSPR → Abraham solute descriptors (E, S, A, B, V)
       → LSER → LSS theory → k(φ)
```

**Prediction of Solute Descriptors for Linear Solvation Energy Relationships Using K-Nearest
Neighbors, Group Contributions, and Graph-Convolutional Neural Networks.**
ACS Environ. Au. <https://doi.org/10.1021/acsenvironau.6c00063>
Local folder: [`papers/ulrich-2026-solute-descriptor-prediction/`](../papers/ulrich-2026-solute-descriptor-prediction/README.md) (main text + SI1)

Three competing approaches to the QSPR layer specifically (Ulrich et al., UFZ) — a group
contribution model, a k-nearest-neighbours model, and a graph-convolutional network, plus
a consensus of the three. Received Feb 2026, accepted May 2026, so this is the most recent
source in the folder. What it gives this project:

- **End-to-end validation through LSER equations** (Table 2) — the closest published thing
  to a test of our chain. On **CHI**, the reversed-phase index, n = 204: ACD Absolv 3.19,
  GNN 3.64, consensus 3.80, kNN 4.41, GCA 5.38, QSPR-LSERD 5.52. ⚠ **The ranking inverts
  against descriptor-space rmse**, where the GNN wins — see
  #41.
- **Independent test-set rmse** (Table 1, n ≈ 634): GNN consensus A 0.075, B 0.095,
  S 0.169, E 0.086, L 0.246; kNN A 0.126, B 0.151, S 0.216, E 0.131, L 0.454; GCA
  A 0.119, B 0.142, S 0.234, E 0.188, L 0.420. These are **RMSE**; Chung's SoluteML
  figures below are **MAE**, so the two tables must not be compared directly.
- **A confidence stratum that works** — thresholding on consensus ensemble SD leaves ~30%
  of chemicals "high confidence" with 2–10× lower rmse (A 0.009 vs 0.087). Feeds
  #42.
- **Dataset**: 6,364 chemicals (`S`), 6,358 (`A`), 6,347 (`E`), 6,111 (`B`), 5,572 (`L`),
  CC-BY, at [`nadinulrich/solute_descriptor_prediction`](https://github.com/nadinulrich/solute_descriptor_prediction)
  (May 2026). ⚠ It is **the Absolv database via LSERD** — a larger copy of the corpus #5
  audited, **not an independent lineage**. Models also served via
  [PAULY](https://i-am-pauly.com/), a UFZ spin-off.
- **Temperature belongs in the phase parameters**, not the solute descriptors — stated
  outright in the conclusions (#11).
- **Ionic chemicals are not covered at all**, from the descriptor side rather than the
  system-constant side (#39).
  Also names structural classes where the intermediate itself is deficient: intramolecular
  H-bonding, site-specific and sterically hindered H-bond donation, and the several
  interactions `S` collapses into one number.

The SI carries a citable statement of the Absolv ground-truth defect: in Abraham's original
Absolv solver files, some `A` values are **set to zero before the descriptors are
determined**, introducing errors especially for diones/tautomers (SI1-4). It also reports
that a multitask descriptor model was *not* convincing on their data — single-task GNNs
won (bears on #35, since SoluteML is multitask).

**Chung, Vermeire, Wu, Walker, Abraham, Green — *Group Contribution and Machine Learning
Approaches to Predict Abraham Solute Parameters, Solvation Free Energy, and Solvation
Enthalpy*.** J. Chem. Inf. Model. 2022, 62, 3, 433–446.
<https://doi.org/10.1021/acs.jcim.1c01103>
Local folder: [`papers/chung-2022-abraham-solute-parameters/`](../papers/chung-2022-abraham-solute-parameters/README.md) (main text + SI)

**The SoluteML source paper** — the model #6 recommended and #35 plans to re-fit on WSU
descriptors. Three models: SoluteGC (RMG group contribution), SoluteML (Chemprop D-MPNN,
25-network ensemble), DirectML (solvation properties directly). What it gives this project:

- **Per-descriptor test errors** (SI Figs. S11/S12): SoluteML MAE on a random solute split
  E 0.041, S 0.098, A 0.038, B 0.047; on a substructure (novel-scaffold) split roughly 2×
  (E 0.084, S 0.170, A 0.067, B 0.088). ⚠ Scored against SoluteDB, which #5 showed is
  partly Absolv-chain output — true errors plausibly larger; #35's bias measurement is the
  remedy.
- **Availability**: data CC BY 4.0 on [Zenodo 5792296](https://zenodo.org/record/5792296)
  (SoluteDB, 8,366 solutes); code at
  [chemprop_solvation](https://github.com/fhvermei/chemprop_solvation) + conda; SoluteGC
  groups in RMG-database.
- **Envelope**: neutral, nonionic solutes only, elements H C N O S P F Cl Br I — a
  structural `OUT_OF_ENVELOPE` test consistent with #33/#34's refuse tier.
- Predicts `L`, which we discard (V is McGowan-computed — four QSPR targets, #5). Solute
  lineage is Abraham-experimental (Abraham is a co-author). Validates solvation energies,
  **not** retention — says nothing about elution order.

## Tooling

| Tool | Role | Links |
|---|---|---|
| **DeepChem** | ML platform; graph-convolutional models for the QSPR layer | [docs](https://deepchem.readthedocs.io/en/latest/) · [GitHub](https://github.com/deepchem/deepchem) |
| **ChemAxon** | Optional high-accuracy ionisation backend (`pKa`, `logD`, microspecies). Commercial licence | [Python API](https://apidocs.chemaxon.com/python_api/apidocs/chemaxon.html) · [docs](https://docs.chemaxon.com/latest/index.html) |
| **RDKit** | Structure handling and descriptors. No `pKa` | [rdkit.org](https://www.rdkit.org/) |
| **OPERA** (US EPA) | Open, regulatory-accepted `pKa` / `logD` / `logP` — candidate default ionisation provider | [GitHub](https://github.com/kmansouri/OPERA) |

## Data sources under consideration

| Source | What it gives | Why it matters here |
|---|---|---|
| **RepoRT** | ~370 retention datasets **with column and gradient metadata** | The only public source that varies *method conditions* — required to validate `k(φ, pH, column)` |
| **METLIN SMRT** | ~80,000 compounds, one column, one gradient | Large but single-method; useful only for the solute-descriptor layer |
| **Acree / UFZ-LSER databases** | ~8,000 compounds with *experimental* Abraham **solute descriptors**. **Not** a source of system constants | Validates the QSPR layer in isolation |
| **Wayne State (Poole & Atapattu)** | LSER **system constants** `c,e,s,a,b,v` on a 7-point φ series (10–70%, 10% steps), with per-coefficient SDs, at **45 °C**, plus a per-column exclusion register | **WSU-2019 obtained and extracted** — 25 MeOH / 24 ACN / 2 THF columns, the 94-compound Table S-2 descriptor lineage, and the Table S-1 caveat register, in [`sources/wsu-lser/`](wsu-lser/README.md). ⚠ 22 of the 354 published fits sit inside a documented incomplete-wetting range and are unmarked in the source — #40. WSU-2025 (five modifiers) still outstanding — see #22 |
| **PQRI / Snyder–Dolan HSM database** | **819** commercial columns, 6 parameters each (H, S\*, A, B, **C(2.8) and C(7.0)**) plus a silica-generation `type` field | **Obtained** — `database.csv` in [`sources/hsm-column-db/`](hsm-column-db/README.md), CC BY-NC-SA 3.0 US. The `C` terms are the ion-exchange axis #39 needs; joined to the USP SRM-870 set in [`experiments/usp-hsm-join/`](../experiments/usp-hsm-join/README.md). ⚠ Duplicate names with conflicting parameters (`Betasil C18`, `Hypersil ODS`) — see that README |

> **Two descriptor lineages, and they are not interchangeable.** LSER system constants are
> regression coefficients fitted against a *specific* solute descriptor set. The Wayne State
> constants and the Acree/UFZ descriptor compilations come from different lineages, so pairing
> one's constants with the other's descriptors introduces a systematic bias. See
> #23.
