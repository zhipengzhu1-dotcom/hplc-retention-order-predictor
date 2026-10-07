# WSU LSER tables — acquisition history (superseded)

> **Superseded 2026-08-15.** This documents the *attempt* to obtain the WSU tables and the
> open-access surrogate data found along the way. The tables were subsequently obtained in
> full — see [`README.md`](README.md) for the actual extraction. Kept because it is the
> provenance for the `luna-c18-*.csv` files in this folder.


Resolves (partially) ticket #22. Read this alongside `sources/README.md` and the issue #3
resolution comment, which established what the WSU-2019/WSU-2025 database *contains* before this
ticket tried to obtain it.

## Update, 2026-08-15 — the WSU-2019 PDF was obtained locally and is now fully extracted

A copy of WSU-2019 was found on the user's machine (gitignored, not committed). It has been read
end to end. See **`wsu2019-pages-9-18-findings.md`** for the second half of the paper and
**`sources/README.md` / issue #22** for the first half. The status below ("paywalled, not
obtained") still holds for **WSU-2025**, for the **WSU-2020 descriptor paper**, and — critically —
for **WSU-2019's Supporting Information**, which is where the per-column `c,e,s,a,b,v` × φ tables
actually live. The main article contains no such table.

### Files in this directory, with provenance

| File | Table | Journal page | Source |
|---|---|---|---|
| `wsu2019-table2-kinetex-c18-vs-xb-c18-deviations.csv` | Table 2 | 118 | WSU-2019, DOI 10.1016/j.chroma.2019.04.027 |
| `wsu2019-table3-group-mean-centred-system-constants.csv` | Table 3 | 120–121 | WSU-2019, DOI 10.1016/j.chroma.2019.04.027 |
| `wsu2019-intext-system-constant-differences.csv` | (unnumbered, running text) | 119, 121 | WSU-2019, DOI 10.1016/j.chroma.2019.04.027 |
| `luna-c18-binary-system-constants.csv` | Table 2 | — | Poole & Atapattu, *J. Chromatogr. Open* 2 (2022) 100039, DOI 10.1016/j.jcoa.2022.100039 (open access) |
| `luna-c18-ternary-system-constants.csv` | Table 3 | — | as above |
| `luna-c18-solvent-strength-and-gas-transfer.csv` | Tables 1 & 4 | — | as above |

Units: system constants `c,e,s,a,b,v` are dimensionless coefficients of the solvation parameter
model on a **log₁₀** basis (`log k = c + eE + sS + aA + bB + vV`); φ is % v/v organic modifier;
temperature 45 °C.

## Bottom line (as of the earlier open-access search — superseded in part by the update above)

**The two anchor papers (WSU-2019, WSU-2025) themselves remain paywalled.** Every legitimate open
route checked — Wayne State's own DigitalCommons repository, ResearchGate author pages, PMC,
SSRN, and publisher SI pages — comes up empty for those two specific papers. No sci-hub or
credential bypass was attempted, per instructions.

**However, a third, fully open-access paper by the same authors (Poole & Atapattu), using the same
instrument, protocol, and descriptor database, was found and its complete data extracted.** It
does not cover the 27/25-column database, but it does answer two of the ticket's four questions
outright (temperature convention, descriptor lineage) with real numbers, and gives a full,
real `c,e,s,a,b,v` table across a φ series for one RPLC column and two of the target modifiers
(methanol, acetonitrile) plus three more (acetone, THF, 2-propanol), with standard deviations and
regression statistics — genuine WSU-lineage data, just not the 27/25-column database itself.

## What was obtained

**Source:** Poole, C.F.; Atapattu, S.N. "Study of system properties in reversed-phase liquid
chromatography for binary and ternary solvent mobile phase compositions using the solvation
parameter model." *Journal of Chromatography Open* **2** (2022) 100039.
DOI: [10.1016/j.jcoa.2022.100039](https://doi.org/10.1016/j.jcoa.2022.100039). **Open access,
Creative Commons license** — read directly on ScienceDirect, full text, no paywall.
(Note: automated `curl`/WebFetch fetches were 403'd by ScienceDirect's bot defenses even though the
article is open access; a real browser session reads it fine. This is a bot-detection artifact, not
a paywall — legitimate to read in-browser.)

- Column: **Luna C18(2)**, 100 mm × 4.6 mm, 5 µm particle, 10 nm pore (Phenomenex) — one column,
  not the 27/25-column database.
- Modifiers: acetonitrile, acetone, methanol, 2-propanol, tetrahydrofuran, each 20–70% v/v in
  10% steps (6 points, not the ticket's 7-point 10–70% series — this paper starts at 20%).
- Full `c,e,s,a,b,v` system-constant tables with per-coefficient standard deviations, plus r, SE,
  F, and n for every regression → `luna-c18-binary-system-constants.csv`.
- Four ternary solvent systems (ACN-MeOH-water, THF-iPrOH-water, ACN-iPrOH-water,
  MeOH-THF-water), each 1:1:2 v/v (50% total organic) → `luna-c18-ternary-system-constants.csv`.
  Two rows (ACN-iPrOH-water, MeOH-THF-water) have an uncaptured `n` value — flagged in the CSV,
  not guessed.
- Solvent strength (S) values and gas-to-solvent-transfer system constants for the same five
  organic solvents plus water → `luna-c18-solvent-strength-and-gas-transfer.csv`.

### 1. Column list, manufacturer, phase chemistry

**Not obtained for the 27/25-column database.** Only the ~17 columns named in open abstracts
(see issue #3's resolution comment) are known, still a floor not a list. This paper adds one
confirmed column (Luna C18(2), Phenomenex) already on that floor list.

### 2. c,e,s,a,b,v tables across the φ series, with SDs

**Not obtained for the 27/25-column WSU-2019/2025 database.** Obtained in full, with SDs and
regression statistics, for **one column** (Luna C18) across **five modifiers** including methanol
and acetonitrile, at 20–70% v/v in 10% steps — see the CSVs in this directory. This is genuine
Poole/Atapattu system-constant data, same lab, same protocol family, same descriptor lineage as
the WSU-2019/2025 papers, but it is a different, smaller published dataset, not the database the
ticket targets. Treat it as directly informative for prototyping the φ-interpolation design fork
(issue's "new questions raised" item), not as a substitute column set.

### 3. Temperature — the item issue #11 depends on

**Found for this specific paper, not confirmed identical for WSU-2019/WSU-2025.**

> "Retention factors at **45 °C** and a flow rate of 1.5 mL/min were determined using a Waters
> Corporation Alliance 2695 HPLC unit... consisting of a 996 photodiode array detector, quaternary
> pump, and a column oven compartment."
> — Poole & Atapattu 2022, Section 2.2 (Instrumentation and measurements)

This is a real, sourced number (45 °C) from the same two authors, same research program, same era
(published between WSU-2019 and WSU-2025), doing the same kind of measurement (RPLC retention
factors for solvation-parameter-model regression) with the same style of instrumentation
description. It is the strongest evidence available anywhere in the open literature for what
temperature this research group runs these measurements at. **It is not a confirmed answer for
the WSU-2019/WSU-2025 database specifically** — those two papers' abstracts do not state a
temperature (confirmed again during this session via PubMed), and the column-oven temperature is
exactly the kind of methods-section detail that would not survive into an abstract. Do not treat
45 °C as certain for #11's purposes; treat it as a strong prior pending confirmation from the
Methods section of WSU-2019 or WSU-2025 once obtained.

### 4. Solute descriptor set / descriptor lineage (issue #23)

**Answered, with a citable chain, for both the surrogate paper and the WSU-2025 anchor paper
itself (via its PubMed abstract, which does state this — the temperature is the abstract-omitted
detail, not the descriptor set).**

- WSU-2025 abstract (PMID [40987226](https://pubmed.ncbi.nlm.nih.gov/40987226/), confirmed
  directly): *"The Wayne State University 2025 (WSU-2025) compound descriptor database is used as
  the sole source for compound descriptors to unify the system constants with a verified
  descriptor database."* This is as direct a confirmation as issue #23 could ask for: the WSU-2025
  system constants are regressed against the WSU's own descriptor database, not
  UFZ/Acree/Abraham-literature descriptors.
- Poole & Atapattu 2022 (the surrogate paper) states the calibration-compound descriptors "were
  taken from the most recent version of the WSU compound descriptor database", citing:
  Poole, C.F. "Wayne State University experimental descriptor database for use with the solvation
  parameter model." *J. Chromatogr. A* **1617** (2020) 460841.
  DOI: [10.1016/j.chroma.2019.460841](https://doi.org/10.1016/j.chroma.2019.460841). (This is the
  WSU-2020 descriptor database; WSU-2025 supersedes it per the search in issue #3's resolution and
  the ScienceDirect listing for "The complete 2025 Wayne State University compound descriptor
  database", DOI 10.1016/j.chroma.2025.466... — abstract page only, not fetched in full here.)
- This reconfirms and sharpens issue #3's finding: **WSU system constants and UFZ/Acree solute
  descriptors are two different lineages and are not interchangeable** — now with the exact paper
  and DOI naming the WSU-side descriptor source explicitly, rather than inferred from context.

## What remains unobtainable through open channels, and what a human must do

Both anchor papers are confirmed still paywalled as of this session (checked PubMed, ScienceDirect
directly, ResearchGate, SSRN, Wayne State's DigitalCommons repository — none carry open full text):

1. **Poole, C.F.** "Reversed-phase liquid chromatography system constant database over an
   extended mobile phase composition range for 25 siloxane-bonded silica-based columns."
   *J. Chromatogr. A* **1600** (2019) 112–126.
   DOI: [10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027).
   PubMed: [31128882](https://pubmed.ncbi.nlm.nih.gov/31128882/) — abstract only, "Full Text
   Sources: Elsevier Science" (subscription).

2. **Poole, C.F.; Atapattu, S.N.** "Update of the Wayne State University system constant database
   for reversed-phase liquid chromatography columns for varied mobile phase compositions."
   *J. Chromatogr. A* **1762** (2025) 466385.
   DOI: [10.1016/j.chroma.2025.466385](https://doi.org/10.1016/j.chroma.2025.466385).
   PubMed: [40987226](https://pubmed.ncbi.nlm.nih.gov/40987226/) — abstract only, same paywall.

3. Also relevant and also paywalled: **Poole, C.F.** "Wayne State University experimental
   descriptor database for use with the solvation parameter model." *J. Chromatogr. A* **1617**
   (2020) 460841. DOI: [10.1016/j.chroma.2019.460841](https://doi.org/10.1016/j.chroma.2019.460841).
   (The WSU-2020 descriptor database paper — confirmed via ResearchGate to have no free full text,
   only a "Request full-text PDF" author-contact button.)

**Routes checked and exhausted for these three papers:**
- ScienceDirect direct (403 to automated fetch each time; even a real browser session shows the
  standard "Get rights and content" / no-access wall for non-open ScienceDirect articles — unlike
  the 2022 Journal of Chromatography Open paper, which is explicitly marked "Open access" and
  reads freely).
- PubMed abstracts (obtained, no temperature or full tables — abstracts confirmed above).
- ResearchGate author pages: exists for the descriptor-database paper, "Request full-text PDF"
  only (i.e. author-permission gate, not a public download).
- SSRN: searched; found *different*, related Poole/Atapattu preprints (a "Polemic" piece and a
  "binary/ternary solvent system properties" piece under review), not the WSU-2019 or WSU-2025
  anchor papers themselves.
- Wayne State DigitalCommons dissertations repository: contains a 2012 PhD dissertation by
  T.N. Karunasekara from the Poole group ("Determination of descriptors by liquid-liquid partition
  and chromatographic methods") — a different topic/author, not Atapattu's own dissertation, and
  predates both anchor papers. Not pursued further; low probability of containing the 2019/2025
  column tables.
- Journal of Separation Science (Wiley) 2025 case-study paper by Atapattu
  ([10.1002/jssc.70283](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/abs/10.1002/jssc.70283))
  — abstract-only page found via search, not opened/verified in this session; worth a follow-up
  check, as Wiley sometimes exposes more of the methods section on the abstract page than Elsevier
  does.

**Human action needed to close this ticket fully:**
- Obtain the two DOIs above (10.1016/j.chroma.2019.04.027 and 10.1016/j.chroma.2025.466385) via
  institutional access or interlibrary loan. Once obtained, also request their Supporting
  Information files specifically — SI is sometimes distributed separately and may be openly
  linked even when the main article is not; this was not checked because the article pages
  themselves were inaccessible to confirm whether SI links exist.
- Alternatively, or in parallel, request
  reprints directly from Colin F. Poole (corresponding author, Wayne State University Department
  of Chemistry) — this is standard scientific practice and the fastest realistic route to both
  papers plus explicit confirmation of the measurement temperature.

## Ticket status

Issue #22 is **left open**. The constants for the actual 27/25-column database were not obtained;
per the instructions for this task, partial acquisition is reported honestly rather than the
ticket being closed. #11 (temperature) and #13 (column set) remain blocked on the same two DOIs;
this work gives #11 a strong, sourced prior (45 °C) to work from but not a confirmed answer.
