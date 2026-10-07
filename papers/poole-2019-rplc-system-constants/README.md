# Poole 2019 — WSU reversed-phase system constant database

> C. F. Poole, *Reversed-phase liquid chromatography system constant database over an
> extended mobile phase composition range for 25 siloxane-bonded silica-based columns*,
> J. Chromatogr. A **1600** (2019) 112–126.
> <https://doi.org/10.1016/j.chroma.2019.04.027> — © Elsevier, not redistributable.

| File | Role |
|---|---|
| `Reversed-phase liquid chromatography system constant database…pdf` | main text |
| `1-s2.0-S0021967319303954-mmc1.docx` | electronic supplementary material (Tables S-1…S-5, Figures S-1…S-3) |

The source of the **column vectors** this project builds on: `c, e, s, a, b, v` for 25
columns over a 7-point φ series at **45 °C**. The supplementary docx also carries the
94-compound descriptor lineage (Table S-2) the constants were regressed against, and the
per-column exclusion register (Table S-1).

Extracted data and full provenance: [`sources/wsu-lser/`](../../sources/wsu-lser/README.md).
`sources/wsu-lser/extract_wsu.py` reads the docx **from this folder** and must be run from
the repo root.

Figures S-1…S-3 are embedded PowerPoint slides inside the docx, at `word/embeddings/*.sldx`
— they do not render in a plain text extraction. Their content is summarised in the
`sources/wsu-lser/` README.
