# Instrument dwell, extra-column and flow-cell volumes — supplied 2026-08-16

Supplied by the project owner as two slides. Transcribed to
`instrument-volumes.csv` with the gaps preserved as gaps.

| Definition (as given on the slide) | |
|---|---|
| **Dwell volume** | gradient delay volume |
| **Extra-column volume** | band broadening after separation |

## The numbers

**Waters** — sourced to the Waters white paper *Dwell Volume and Extra Column Volume:
What Are They and How Do They Impact Method Transfer?* (2018), document **720005723**,
<https://www.waters.com/content/dam/waters/es/library/white-papers/2018/waters-whitepaper-DwellVolumeandExtraColumnVolumeWhatAreTheyandHowDoTheyImpactMethodTransfer-720005723>

| Instrument | Dwell µL | Extra-column µL | Flow cell µL |
|---|---|---|---|
| Acquity (binary pump) | **73** | 5 | 0.5 |
| Acquity H-Class (quaternary pump) | **375** | 7 | 1 |
| ARC HPLC | **1190 / 770** | ≤ 30 (19?) | 10 |
| Alliance 2695 | **1150** | 30–34 | 10 |

**Agilent** — ⚠ **no source shown on the slide, and it carries the owner's own open
questions.** Treat as provisional until sourced.

| Instrument | Dwell µL | Extra-column µL | Flow cell µL |
|---|---|---|---|
| 1100 | 1100 | 30 | 13 |
| 1200 | — | 30 | 13 |
| 1260 | 1100 | — | 8, 13 |
| 1290 Infinity II — HPLC | 900 (quaternary, with mixer); without mixer **?** | ? | 1, 4 |
| 1290 Infinity II — UHPLC | with mixer **?** / without mixer **?** | ? | 1, 4 |

## The spread is the point

Dwell volume ranges **73 µL to 1190 µL across Waters alone — a factor of 16**, and it
varies *within* a model with the pump type and whether a mixer is fitted. It is not a
vendor constant. Any use of a nominal figure must carry that as a width, not a point.

At the flow rate of RepoRT dataset 0415 (0.3 mL/min) the tabulated values become a
gradient delay of:

| Dwell µL | Delay at 0.3 mL/min |
|---|---|
| 73 | 0.24 min |
| 375 | 1.25 min |
| 900 | 3.00 min |
| 1100 | 3.67 min |
| 1190 | 3.97 min |

0415's slice compounds elute between 4.62 and 12.85 min, after a 1 min initial hold. So
the choice of instrument moves the effective gradient start by **up to ~3.7 min** —
comparable to the whole early part of the elution window. This is not a small correction.

## Two uses, and the second may be the larger one

**1. It changes what RepoRT dataset 0415 needs from its paper.**
RepoRT records the *column* but never the LC system, so a dwell lookup table cannot be
applied to 0415 on its own. What it does is change the fact we must recover from the
source paper (Nürenberg et al., *J. Chromatogr. A* **1426** (2015) 77–90,
[doi:10.1016/j.chroma.2015.11.014](https://doi.org/10.1016/j.chroma.2015.11.014)) from
**the dwell volume**, which method papers rarely state, to **the instrument model**,
which they almost always state. That is a large improvement in the odds of success and
it is why this table matters to
#45.

It does **not** make 0415 usable by itself. A nominal dwell for a named instrument still
carries the configuration ambiguity above, and 0415's `column.t0` remains **imputed from
geometry, not measured** (see `experiments/report-overlap/README.md`). Two error terms,
both quotable, neither zero.

**2. It is material the instrument layer owes anyway.**
#44's punch list item 3
says *width needs physics* — `width[S, n]` was fillable but physics-free in the thin
slice, and the instrument layer owes it real content including **extra-column
dispersion**. The extra-column and flow-cell columns here are exactly that input, and
they show the same shape: 5 µL for an Acquity binary against 30–34 µL for an Alliance
2695, so extra-column dispersion is a per-run term with a real spread, not a constant.

Per `CONTEXT.md`, dwell and extra-column volume are **per-run error**, not per-compound —
they attach to the method card, not the molecule.

## What is still open

- **Source the Agilent numbers.** Five of the fifteen Agilent cells are blank or a
  question mark, including both 1290 Infinity II dwell figures. Agilent publish these in
  their system manuals and in method-transfer notes.
- **Resolve the ARC HPLC 1190/770 pair** — the slide gives two values without saying
  which configuration each belongs to.
- **Decide the width, not just the value.** For the spec these must enter as
  distributions. The vendor tables give central values; the configuration ambiguity
  (mixer fitted or not, binary or quaternary) is the width, and it is large enough to
  matter.

## Licence

Vendor-published specifications, transcribed as facts for internal development use.
Cite the Waters white paper (720005723), not this file. The Agilent rows are unsourced
and must not be cited at all until they are.
