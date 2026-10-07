# USP Column Equivalency database — extracted 2026-08-15

> **The data is not included in this repository.** The USP Column Equivalency data is
> © The United States Pharmacopeial Convention, and no licence to redistribute it has been
> identified. The CSV and workbook files described below were removed from the public copy.
> To use the analyses that depend on them, build your own copy (see
> [How to obtain the data](#how-to-obtain-the-data)).

| File | Contents | In this repo |
|---|---|---|
| `usp-approach-columns.csv` | **145 columns** with `Hy, CTF, CFA, TFA, BD` | no — build it yourself |
| `usp-similarity-long.csv` | **17,891 rows** — the full ranked equivalency list for every usable reference column | no — build it yourself |
| `USP-column-equivalency.xlsx` | Both of the above plus a 143 × 145 `F` matrix and a README sheet | no — run `build_xlsx.py` |
| `build_xlsx.py` | Rebuilds the workbook from the two CSVs | yes |

## How to obtain the data

1. Open the USP Column Equivalency Application at
   <https://apps.usp.org/app/USPNF/columnsDB.html> and accept its terms of use. Check that
   your intended use is allowed by those terms.
2. Record the "USP approach" parameters (`Hy, CTF, CFA, TFA, BD`) for each column in
   `usp-approach-columns.csv`, with the header
   `name,Hy,CTF,CFA,TFA,BD,usp_designation,manufacturer`. Leave a cell empty when the
   source has no value.
3. Optional: record the ranked similarity list for each reference column in
   `usp-similarity-long.csv`, with the header `reference,rank,F,target`. The
   [Method](#method) section below describes how the original extract was made.
4. Optional: run `python sources/usp-column-db/build_xlsx.py` (needs `openpyxl`) to make
   the workbook.

Only `experiments/usp-hsm-join/join_usp_hsm.py` reads these files. The prediction model
and the test suite do not.

## `F` — the equivalency measure

`F` is USP's column-difference measure; smaller means more similar, and `F = 0` is a column
against itself. **Verified exactly symmetric**: across all 7,501 pairs where both directions were
captured, `max |F(a,b) − F(b,a)| = 0.0`. That is a strong check that the extraction is faithful.

⚠ **Two columns cannot be used as a reference — a defect in USP's own application.**
`TSKgel ODS-100V 3 μm` and `TSKgel ODS-100V 5 μm`. Selecting either as the reference returns an
**empty rank-0 row** and then ranks every other column against a **null parameter vector**,
yielding silently inflated `F` values rather than an error. **Reproduced through the site's own
UI**, so it is not an artefact of this extraction. Their outgoing rankings are therefore
**omitted**; their own parameters and their *incoming* `F` values (as targets in other columns'
lists) are valid and retained.

This is how the defect was caught: after extraction those two were the **only** asymmetric
references in the entire matrix, and removing them took the asymmetry to exactly zero.

So the similarity data is **143 reference columns × 145 targets**.

## What this is, and what it is not

Extracted from the **USP Column Equivalency Application**
(<https://apps.usp.org/app/USPNF/columnsDB.html>), which hosts **two** databases behind one
terms-of-use gate:

| Database | Parameters | Columns | Status here |
|---|---|---|---|
| **USP approach** | `Hy, CTF, CFA, TFA, BD` | **145** | **extracted → this CSV** |
| PQRI / hydrophobic-subtraction | `H, S*, A, B, C(2.8), C(7.0)` | 815 | **not extracted** — we hold the fuller 819-row free CSV from hplccolumns.org, in [`../hsm-column-db/`](../hsm-column-db/README.md) |

⚠ **This is not the HSM/Snyder parameter set.** It is a different characterisation with different
physical meaning, so it does **not** substitute for, and is **not** comparable to, the
`H, S*, A, B, C` values used for `Fs`.

## Provenance

Developed by the **USP Working Group on Column Equivalency using NIST SRM 870**; method review
published in *Pharmacopeial Forum* **31**(2), 637–645 (2005). Source: USP's own "About USP
approach" page, <https://www.usp.org/node/12521>.

**This provenance is the interesting part.** SRM 870 is a column *performance/activity* test
mixture, so these parameters measure things like silanol activity, metal chelation and bonding
density — **not** retention selectivity in the LSER sense. That makes the dataset complementary
to both LSER and HSM rather than redundant with either.

**Parameter expansions: confirmed by the project owner (2026-08-15), documentary verification
deferred.** `BD` is bonding density (values 0.9–5.5 are right for µmol/m²) and `TFA`/`CTF` are
tailing-factor-like activity measures — the readings inferred here from the SRM 870 provenance
are correct. USP's public page still does not define the abbreviations and the *Pharmacopeial
Forum* article (31(2) 637–645) has not been obtained, so **the confirmation rests on the
owner's knowledge rather than on a read source**; it will be verified against the reference
when the finished model is tested.

Two evidence classes, kept apart as usual: the parameters may now be **used and named** in the
spec, but a *citation* to PF 31(2) still requires reading it. PF 31(2) therefore stays on
#22's acquisition list for #39's benefit, at lower priority than before.

## Why it may matter to this project

#36 measured that the HSM `C`
term (ion exchange) is **80.5% of `Fs²`** and is **unparameterised by the solvation parameter
model**. #39 then asked whether
permanently cationic compounds can be placed at all, given their retention runs through exactly
that unmodelled silanol-ion-exchange mechanism.

If these parameters are what the SRM 870 provenance implies, this CSV is a **silanol-activity and
bonding-density proxy for 145 columns** — a direct handle on the axis LSER cannot see. That is a
live input to #39, and potentially to the `D` parameter's per-run (column silanol activity)
component from #34.

**Confirm the parameter definitions before relying on this.**

## Data quality caveats

- ⚠ **`usp_designation` is `L1` for every one of the 145 rows** — including columns that are
  physically C8 or phenyl (`Acclaim 120 C8`, `TSKgel Octyl-80Ts`, `TSKgel Super-Octyl`,
  `TSKgel Super-Phenyl`). Verified against the site's own rendered table, so this is **the
  source's data, not an extraction artefact**. **Do not use this field as column chemistry.**
  Same class of defect as the RepoRT `name`/`column.name` swap already on the map.
- **21 of 145 rows have at least one missing parameter**; only **124 are complete**. That exactly
  matches the `c_total = 124` the site's own endpoint reports, confirming that incomplete columns
  are excluded from its `F` similarity ranking. Missing counts: `CTF` 12, `CFA` 7, `TFA` 5,
  `BD` 5, `Hy` 1.
- One row (`Aqua 5 µ C18 125A`) has **no manufacturer** in the source.
- 24 manufacturers represented.

## Licence and permitted use

© The United States Pharmacopeial Convention. Accessed under the application's terms of use,
which the project owner accepted interactively. **No licence grant for redistribution has been
identified.**

Treat as **internal development and educational use only**, consistent with #28. For
that reason the CSVs and workbook are **not distributed with this repository**; build your
own copy under the application's terms of use. If the data becomes load-bearing for the spec, the
*Pharmacopeial Forum* article is the citable primary source, not this extract.

## Method

The app is a jQuery/cpaint AJAX front end posting to `/ajax/USPNF/columnsDB.php`, which returns
XML (`urow0…urowN`) with `c_total` and `c_start` for pagination. The request body is
`cpaint_function=updateResults` with arguments `[columnName, 1, 1, 1, 1, 1, startOffset]`.

**Parameters (145 columns)** — a pagination sweep recovered 129; a targeted pass fetched the 16
absent from the ranking. Those 16 are precisely the columns with incomplete parameters, which is
why USP's own ranking omits them.

**Similarity (17,891 rows)** — every column queried as the reference, paging until `c_start`
reached `c_total`, roughly 1,900 requests at modest concurrency. Data was moved out of the browser
in six chunks; five rows lost at chunk boundaries were re-fetched individually, and the row count
and `F` symmetry were both verified afterwards.

Rate-limited throughout. Reproducible from the page with no credentials beyond accepting the
terms.
