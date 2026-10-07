# #45's bounded RepoRT search — run 2026-08-16

`python3 search_report.py /path/to/RepoRT` reproduces everything here; `results.txt`
is the captured output. RepoRT was shallow-cloned from
<https://github.com/michaelwitting/RepoRT> (**CC BY-SA 4.0**, `master`, 421 processed
datasets, ~824 MB on disk). Nothing from RepoRT is vendored into this repo.

## Verdict: RepoRT cannot supply what #45 needs. Reportable negative.

Per #30's bounded gate, "no usable public data found" is a **reportable outcome**, and
this is one. The compound overlap is genuinely good; the *conditions* are unusable, for
three independent reasons, any one of which is sufficient.

## The overlap is real — that is not the problem

**205 of 421 datasets** contain at least one of the 15 valid thin-slice compounds
(the refused quaternary ammonium is excluded by construction, #33). Matching is on
**InChIKey skeleton** (first block), so stereochemistry and salt form do not block a
match — that matters for ibuprofen, naproxen and ketoprofen, which appear in both
racemic and single-enantiomer records.

| ds | n/15 | mode | pH | T | column |
|---|---|---|---|---|---|
| 0415 | **10** | gradient | — | 40 | Agilent ZORBAX Eclipse Plus C18 |
| 0263 | 9 | gradient | — | 50 | Waters ACQUITY UPLC BEH C18 |
| 0382 | 9 | gradient | 3.61 | 30 | Thermo Acclaim RSLC 120 C18 |
| 0383 | 9 | gradient | 6.2 | 30 | Thermo Acclaim RSLC 120 C18 |
| 0019 | 7 | gradient | 3 | 30 | Waters XBridge C18 |

## 1. There is no isocratic data at all

**Zero** of the 205 overlapping datasets are isocratic — 196 gradient, 9 with the
programme never recorded. RepoRT is a metabolomics RT-prediction corpus; gradient
elution is what it is made of.

⚠ **An early pass of this search reported nine usable isocratic datasets. They do not
exist.** RepoRT writes a *header-only* gradient file when the programme was not
recorded, and treating "no rows" as "composition never changes" turns missing metadata
into a fabricated isocratic run. The affected sets say so themselves — 0063's
`info` reads `missing information: column length, ID and particle size; temperature,
gradient and flow rate`. The corrected classifier is in `elution_mode()` and the
mistake is documented in its docstring, because it is an easy one to make again.

## 2. Gradient inversion is impossible: the dwell volume does not exist in the schema

Inverting a gradient run to `k` requires the system dwell volume. RepoRT's metadata
schema has **no dwell, delay, system, void or dead volume field at all** — this is not
a value missing from some datasets, the column does not exist. The schema carries
`column.name`, `usp.code`, `length`, `id`, `particle.size`, `temperature`, `flowrate`,
`t0`, and eluent composition/pH, and that is the whole of the instrument description.

So the "or gradient-inverted" half of #45's data requirement is closed off for every
dataset in the repository, not merely the inconvenient ones.

## 3. `column.t0` is imputed from geometry, not measured

Across 421 datasets there are only **47 distinct `t0` values**, and they repeat in exact
decimals — 84 datasets at `0.55125`, 47 at `1.1025`, 42 at `0.441`, 36 at `2.205`. These
are formula outputs (`πr²Lε/F` at an assumed porosity), not measurements: 0415
(150 × 2.1 mm, 3.5 µm, 0.3 mL/min) and 0382 (100 × 2.1 mm, 2.2 µm, 0.2 mL/min) carry
**the same `t0` to four decimals**, which is only possible because both have the same
length/flow ratio.

**48 datasets carry `t0 = 0`**, which is physically impossible and is a missing-value
marker.

Consequence: even if the dwell volume were recovered, every `k` derived from RepoRT
would carry an unquantified porosity assumption in its denominator. For a project whose
whole architecture rests on uncertainty being first-class, that is an error term with no
honest width.

## 4. ⚠ The most attractive contrast in the repository is not two measurements

0382 and 0383 look purpose-built for #27's **pH discrimination gate**: same column, same
temperature, byte-identical gradient programme, 9 slice compounds each, labelled
**pH 3.61 (0.01% formic)** versus **pH 6.2 (5 mM ammonium acetate)**. All the unknowns
above — dwell, `t0`, gradient shape — cancel in a paired comparison, so this looked like
the one usable thing in the search.

**It is not.** The two datasets share 1,673 compounds, of which **1,520 (91%) have
bit-identical retention times.** Among the slice compounds, 7 of 8 are identical to the
digit; only atenolol moves (3.09 → 4.09 min).

That is chemically impossible. Ketoprofen, naproxen and ibuprofen (pKa ≈ 4.2) are
largely neutral at pH 3.61 and largely ionised at pH 6.2, and **must** shift
substantially. Zero inverted pairs out of 28 across a 2.6-unit pH change is not a
finding about pH robustness; it is the signature of duplicated data.

Sibling datasets from the same publication behave normally — 0384/0385 share 83
compounds with **0%** identical RT, 0385/0386 with 2% — so this is specific to the
0382/0383 pair, not a repository-wide convention. 0382 additionally carries RepoRT's own
`model.flags: ESSENTIAL_METADATA_FEATURES_UNAVAILABLE`.

**Anyone using RepoRT for a paired-pH contrast will find this pair first.** It is the
largest, best-annotated, most obviously suitable candidate in the corpus. Logged here so
the project does not walk into it, and it belongs upstream as an issue against RepoRT.

This is a second defect class in RepoRT alongside the `name`/`column.name` swap already
on the map.

## What would change the verdict

**0415 is the one dataset worth revisiting.** 10 of 15 slice compounds, complete
metadata (`missing information` is empty), real column dimensions, stated flow and
temperature, ACN/0.1% formic — the modifier the slice actually ran. It fails only on
the dwell volume and the imputed `t0`.

Both are properties of the *instrument*, not the chemistry, and both may be recoverable
from the source publication: Nürenberg, Schulz, Kunkel & Ternes,
[10.1016/j.chroma.2015.11.014](https://doi.org/10.1016/j.chroma.2015.11.014). If that
paper states the dwell volume and a measured `t0`, 0415 becomes a candidate for a
gradient-inverted test — still one column, one pH, and gradient rather than isocratic,
so it would support an **order-accuracy** check but not #15's per-stratum `λ`.

That is a bounded, cheap next step with a clear failure condition, and it is the only
one this search leaves open.

## What this does *not* license

Running the slice against gradient RTs by comparing predicted isocratic order to
observed gradient order would produce a number, and the number would be uninterpretable.
Gradient and isocratic elution order differ for exactly the compounds in question — the
ionisable ones whose retention changes across the gradient's pH-invariant but
φ-varying path. Reporting a calibration `λ` from that comparison would be a
well-organised guess wearing #15's pre-registered thresholds, which is the specific
failure mode #15 was written to prevent.
