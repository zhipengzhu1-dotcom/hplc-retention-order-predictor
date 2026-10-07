# WSU-2019 LSER system constants and descriptor lineage

Extracted 2026-08-15 from the electronic supplementary material of:

> C. F. Poole, *Reversed-phase liquid chromatography system constant database over an
> extended mobile phase composition range for 25 siloxane-bonded silica-based columns*,
> J. Chromatogr. A **1600** (2019) 112–126. <https://doi.org/10.1016/j.chroma.2019.04.027>

Source file: `papers/poole-2019-rplc-system-constants/1-s2.0-S0021967319303954-mmc1.docx`
(gitignored — © Elsevier, do not redistribute; cite the paper). Extraction script:
`extract_wsu.py` in this directory, **run from the repo root**.

> **The CSV files are not included in this repository.** They were extracted from the
> paper's electronic supplementary material (© Elsevier) and from its main text, so they are
> not redistributed here. The scripts, provenance notes and validation rules are kept. To
> run the analyses that need the data, rebuild it:
>
> 1. Get the supplementary file `1-s2.0-S0021967319303954-mmc1.docx` from the article page
>    (<https://doi.org/10.1016/j.chroma.2019.04.027>, "Supplementary data"). This normally
>    needs journal access.
> 2. Save it as
>    `papers/poole-2019-rplc-system-constants/1-s2.0-S0021967319303954-mmc1.docx`.
> 3. From the repository root, run `python sources/wsu-lser/extract_wsu.py`. It writes
>    `wsu2019-descriptors-s2.csv`, `wsu2019-column-caveats.csv` and
>    `wsu2019-system-constants.csv` into this folder and checks them as described under
>    "Validation" below.
>
> The script rebuilds those three files only. The other CSVs that were here
> (`wsu2019-table2-…`, `wsu2019-table3-…`, `wsu2019-intext-…` from the main text of
> Poole 2019, and `luna-c18-*` from Poole & Atapattu, *J. Chromatogr. Open* **2** (2022)
> 100039, <https://doi.org/10.1016/j.jcoa.2022.100039>) were transcribed by hand; see
> `ACQUISITION-HISTORY.md` and `wsu2019-pages-9-18-findings.md` for what they held. No code
> path in `model/` or the tests reads them.
>
> Code that reads the three rebuilt files: `model/columns.py`, `model/applicability.py`
> (the steric-exclusion register; its tests skip without it), and most scripts under
> `prototype/` and `experiments/`.

Citation: C. F. Poole, *J. Chromatogr. A* **1600** (2019) 112–126,
doi:10.1016/j.chroma.2019.04.027.

**All constants are a 45 °C dataset** (stated in the table titles and the paper's
Experimental). See #22, #11 and `CONTEXT.md`.

## Files

### `wsu2019-descriptors-s2.csv` — Table S-2, 94 compounds

`compound, E, S, A, B, V`. **This is the WSU solute-descriptor lineage anchor** that #23
and the column-vector rule in `CONTEXT.md` refer to: the system constants below are
regression coefficients fitted against exactly these descriptor values, and must not be
multiplied with descriptors from any other lineage.

### `wsu2019-system-constants.csv` — Tables S-3/S-4/S-5, 354 fits

One row per (modifier, column, φ): `c, e, s, a, b, v` with per-coefficient standard
deviations `sd_*` (the parenthesised values in the source), the fit statistics
`r, r2, F, SE, n`, and `incomplete_wetting` (see below). Coverage:

| Modifier | Columns | Rows |
|---|---|---|
| methanol (S-3) | 25 | 175 (all 10–70% v/v) |
| acetonitrile (S-4) | 24 | 165 |
| tetrahydrofuran (S-5) | 2 (SunFire C18, XTerra MS C18) | 14 |

⚠ The often-quoted "27 methanol / 25 acetonitrile columns" is **not** this paper — this
paper is 25/24/2. The larger counts presumably include the later WSU papers (still on
#22's acquisition list).

### `wsu2019-column-caveats.csv` — Table S-1, 45 rows over all 25 columns

`column, reference, caveat, modifiers, note`. Table S-1 is a **per-column exclusion
register**: for each column it names the mobile-phase compositions where the stationary
phase does not wet, the compound classes dropped for steric repulsion, and whether
electrostatic (cation-exchange) interactions matter, each qualified by modifier. `caveat`
is one of `incomplete_wetting`, `steric_repulsion`, `electrostatic`, `none`; `note` is the
source sentence verbatim, so a couple of Word run splits survive inside words
("methanol-wa ter"). Only Synergi Polar-RP carries no caveat.

**This is the table that makes the database safe to consume, and it sits three tables away
from the constants.** 24 of 25 columns carry at least one caveat; the electrostatic flag
is acetonitrile-only on 22 columns and extends to methanol on Fluophase-RP alone.

### The dewetting rows are published anyway — `incomplete_wetting` in the constants CSV

Tables S-3/S-4/S-5 publish a full set of constants for every φ from 10% up, **including
the compositions Table S-1 says the column does not wet**, and nothing in those tables
marks them. **22 of the 354 fits (6.2%), across 13 of the 25 columns**, sit inside a
documented dewetting range. The `incomplete_wetting` column carries the flag; the ranges
are transcribed by hand in `extract_wsu.py` (`DEWET`) because the prose does not parse
reliably — "at 10%" excludes 10 only, "less than 30%" excludes 10 and 20. One note is
genuinely ambiguous and the reading taken is recorded in the code comment.

⚠ **Fit quality does not identify these rows.** They have roughly twice the residual SD of
the rest (median SE 0.068 vs 0.035), but that is a φ effect, not a dewetting signature:
restricted to φ = 10% the two groups are indistinguishable (median SE 0.052, n = 19 vs
0.051, n = 31), and the flagged fits still have median r² = 0.9865. Nothing but Table S-1
tells you they are bad.

## Fit residuals — the measured LSER lack-of-fit

`SE` is the standard deviation of the residuals in **log₁₀ k**, in-sample, on the paper's
own calibration compounds with the paper's own descriptors. It is a **floor** on the LSER
lack-of-fit term, not an estimate of it.

| Modifier | fits | median SE | IQR | p90 | range |
|---|---|---|---|---|---|
| acetonitrile | 165 | 0.030 | 0.026–0.037 | 0.047 | 0.019–0.138 |
| methanol | 175 | 0.041 | 0.034–0.052 | 0.066 | 0.021–0.152 |
| tetrahydrofuran | 14 | 0.043 | 0.036–0.055 | 0.063 | 0.028–0.069 |

Each fit uses a **median of 50 calibration solutes (range 22–68)** drawn from the
94-compound Table S-2 pool — all neutral.

Figures S-1 and S-2 (embedded PowerPoint slides in the docx, `word/embeddings/*.sldx`)
plot these same distributions against the fits recomputed from the original cited sources'
descriptor values. Panel A matches the table above; panel B, read off the figure, has
median SE ≈ 0.046 (acetonitrile) and ≈ 0.060 (methanol) — i.e. **harmonising the
descriptor lineage cut the residual SD by roughly a third**. Figure S-3 plots two weak
bases (diphenylamine on Kinetex EVO C18, quinoline on Discovery HS F5) predicted by the
neutral-compound models against experiment; it is a picture of the ionisation retention
drop `D`, with the point values recoverable only by digitising the plot.

## Source defects found during extraction (all in the docx, not our artefacts)

- **The Luna Omega PS C18 block of S-3 prints F and SE transposed** relative to every
  other block (F = 174 read as SE = 174, and so on for all 7 φ points). Caught by a
  plausibility guard rather than by name — F is a Fisher statistic > 1, SE a residual
  standard deviation ≈ 0.02–0.15, so the two are never confusable. Fixed 2026-08-15;
  earlier copies of this CSV have SE values up to 2706 in those rows.

- **Synergi Fusion-RP has no acetonitrile table at all** — 24 ACN columns, not 25.
- **Synergi Polar-RP (ACN) spans 20–50% v/v only** — 4 φ points, missing both ends.
  ⚠ Corrected 2026-08-15: an earlier version of this line said "starts at 40%" and
  implied wetting; Table S-1 carries **no caveat** for this column, so the short series
  is unexplained by the register, not a dewetting exclusion.
- **One φ value is missing**: the XBridge Phenyl (ACN) 70% row has no leading "70";
  inferred from the sequence (row follows 60%, statistics consistent).
- **"Lunar Omega PS C18"** in S-4 is a typo for Luna Omega PS C18 (S-3 spells it
  correctly); normalised.
- Formatting: one unclosed parenthesis in S-3, stray spaces inside SD parentheses in
  S-4, and run-split names ("SunFire C 18", "1 ,3 -Dibromobenzene"); all normalised by
  the script.

Validation: parser enforces φ ∈ {10..70} strictly increasing per column, r ∈ (0.9, 1],
F > 1 and 0 < SE < 1, exactly 6 SDs per fit, and that the hand-transcribed `DEWET` map
names exactly the columns Table S-1 flags for incomplete wetting; spot checks against
hand-read values pass (Ascentis C18 MeOH 10%, Synergi Fusion-RP MeOH 70%, Kinetex XB-C18
ACN 30%, XBridge Phenyl ACN 70%).
