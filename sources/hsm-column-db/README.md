# PQRI / HSM column database — downloaded 2026-08-15

| File | Contents |
|---|---|
| `database.csv` | **819 columns**, fields `id, name, manufacturer, manufacturerID, type, H, S, A, B, C28, C70, retention, USPtype, phase` |

The hydrophobic-subtraction model (Snyder–Dolan–Carr) parameter set: `H` hydrophobicity,
`S*` steric resistance, `A` H-bond acidity, `B` H-bond basicity, and `C` cation-exchange
activity reported at **two** pH values — `C28` (pH 2.8) and `C70` (pH 7.0). `retention` is
the ethylbenzene retention factor used to normalise for phase ratio.

Retrieved from <https://hplccolumns.org/database/database.csv> — the link is present but
HTML-commented-out on <https://hplccolumns.org/database/index.php>; the file itself serves
fine, `text/csv`, **latin-1 encoded**. Read it as latin-1: several manufacturer and column
names carry high bytes, and BSD `grep`/`cut` treat the file as binary and silently return
nothing under a UTF-8 locale.

## Why it is here

`C28`/`C70` are the **ion-exchange** terms, the mechanism by which permanently cationic
compounds retain — the axis
#39 needs and the one
#36 measured as **80.5%
of `Fs²`** and which Poole 1600 states is unparameterised by the solvation parameter model.
This CSV is the instrument for #39's join against the USP SRM-870 activity measures — see
[`experiments/usp-hsm-join/`](../../experiments/usp-hsm-join/README.md).

## `type` — the silica-generation field, and the sharpest thing in the file

Undocumented on the site but unambiguous from its values and its behaviour: `A` is
old high-activity (type-A) silica, `B` modern low-activity (type-B) silica, `EP`
embedded-polar, plus phase-shape labels (`phenyl`, `CN`, `F`/`fluoro`, `Other`). It
separates on `C70` in exactly the direction the classification implies — type A median
1.050 (n = 75) against type B median 0.207 (n = 516). **That is a silanol-activity axis
stated directly by the source**, and for column *selection* it is a blunter but more
robust instrument than the USP SRM-870 parameters.

⚠ The label vocabulary is **not normalised**: `other`/`Other` and `F`/`fluoro` both appear
as distinct strings, and `phase` carries both `C1` and `C1 ` and both `phenyl` and
`Phenyl`. Normalise case and whitespace before grouping on either field.

## Data quality caveats

- ⚠ **Duplicate column names carrying different parameters.** `Betasil C18` appears twice
  (`C28` −0.03 / `C70` −0.04, and `C28` 0.095 / `C70` 0.099), as does `Hypersil ODS`. This
  is the [Shackman 2016](https://doi.org/10.1016/j.chroma.2016.11.006) discrepancy showing
  up inside a single file, and it matters because `Betasil C18` is one of
  #30's sweep columns.
  **Never resolve a duplicate by taking the first row** — that invents a number. The join
  in `experiments/usp-hsm-join/` excludes them and lists them.
- **`Fs` provenance is unverified.** Shackman 2016 asserts the published HSM parameter
  tables are not internally consistent across sources; the full text is paywalled and has
  not been read, so the size of the discrepancies is unknown. Stands as a check owed before
  the CSV is load-bearing.
- **One measurement condition only** — 50% acetonitrile / phosphate pH 2.8 (30 mM), 35 °C,
  ethylbenzene reference, 16 probe solutes (Zhang & Carr 2009). There is **no `φ`
  dependence in the parameter set at all**, and only `C` is reported at a second pH.
- **Published row counts disagree with this file**: 750 (Rutan et al. 2024), "more than 600"
  (the site's own stale `<title>`), 368 (the USP mirror). **819 is what the file contains**,
  measured; the USP mirror is a smaller, lagging subset.
- **Last-modified date could not be determined**, nor whether HSM3 parameters (Rutan et al.
  2024) have been merged in — the schema carries no version field.

## Licence and permitted use

> "Except where otherwise noted, all content on this site is licensed under a Creative
> Commons Attribution-Noncommercial-Share Alike 3.0 United States License"

— <https://hplccolumns.org/about/index.php>

**CC BY-NC-SA 3.0 US.** This file is **not** covered by the repository's MIT licence; it
keeps its own licence and attribution. Free for a spec, a prototype and academic work; **not**
redistributable inside a commercial product without a separate agreement. Dwight Stoll
(Gustavus Adolphus College) maintains the database. Whether a commercial licence path
exists could not be determined. Same standing constraint as
#28.

Note this is a *different* licence position from the USP extract in
[`../usp-column-db/`](../usp-column-db/README.md), which carries no identified grant at all.
The two datasets must not be redistributed under one blanket statement.

## Citation

Snyder, Dolan & Carr, *J. Chromatogr. A* **1060** (2004) 77–116,
[doi:10.1016/S0021-9673(04)01480-3](https://doi.org/10.1016/S0021-9673(04)01480-3).
Parameter definitions verbatim in Zhang & Carr, *J. Chromatogr. A* **1216** (2009)
6685–6694, [doi:10.1016/j.chroma.2009.06.048](https://doi.org/10.1016/j.chroma.2009.06.048)
(open access: <https://pmc.ncbi.nlm.nih.gov/articles/PMC3195507/>).
