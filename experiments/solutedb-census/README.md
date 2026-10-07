# How much of real chromatographic chemistry has SoluteML already seen?

`env/ionisation-venv/bin/python experiments/solutedb-census/census.py <RepoRT/processed_data> <solutedb_selected.csv>`
— output in `results.txt`, per-compound tiers in `census.csv`. Neither RepoRT nor SoluteDB
is vendored: RepoRT clone/pin per `data/report/README.md` (9de8d60), SoluteDB from Zenodo
5792296 (CC BY 4.0) extracted by `prototype/soluteml-bias/read_solutedb.py`.

## Why this is not the study it started as

This was proposed as a test of whether RepoRT could supply a **clean** (out-of-corpus)
compound set to revive #42's two stratification questions, which #35's bias run showed
cannot be answered on the WSU anchor (all 94 compounds are inside SoluteDB).

**It cannot, and the evidence was already in this repo.** `experiments/report-overlap`
established that RepoRT has **zero isocratic datasets**, **no dwell-volume field in the
schema at all**, and a `column.t0` imputed from geometry rather than measured — so no
honest `k` can be derived from any RepoRT dataset. RepoRT also carries no measured Abraham
descriptors. A clean compound set with no label to score against is not an evaluation set,
in either descriptor space or retention space.

So #42's questions stay blocked, and the honest statement is that **they are blocked on
measurement, not on a split** — no public corpus we have examined can answer them.

What survives needs no labels, and the applicability-domain object (#40/#43) needs it:
for chemistry that actually appears in reversed-phase practice, what has the model seen?

## Result

18,872 unique non-SMRT RepoRT compounds, tiered against SoluteDB's 8,183 InChIKey
skeletons and 1,119 Murcko scaffolds:

| tier | compounds | share | share of the 103,534 retention records |
|---|---|---|---|
| `in_corpus` — skeleton is in SoluteDB | 1,691 | **9.0%** | 23.6% |
| `scaffold_seen` — Murcko scaffold is in SoluteDB | 4,156 | 22.0% | 27.2% |
| `novel_in_envelope` — neither, inside the declared envelope | 12,629 | **66.9%** | 47.7% |
| `out_of_envelope` — charged, or a non-envelope element | 396 | 2.1% | 1.5% |

`out_of_envelope` is 370 charged species plus 26 carrying B, Si, Se, Sn, Co or Na — elements
SoluteML was never trained on. These are #33/#39's refuse-by-structure cases, and it is worth
knowing they are only 2% of real chemistry: the refuse tier is narrow, not a broad exclusion.

## What it means

**SoluteML has seen 9% of the distinct chemistry and two-thirds of it is genuine
extrapolation.** Chung 2022's substructure split — which roughly doubles every descriptor
error against SoluteDB — is the analogue of the `scaffold_seen` tier, so it is the
**optimistic** reference for 22% of compounds and does not describe the 67% at all.

**The measurement weighting is the more useful number, and it points the other way.**
Compounds the model has seen carry 23.6% of retention records against 9.0% of compounds:
memorised chemistry is measured about three times as often as novel chemistry, because the
seen compounds are the common metabolites that appear in dataset after dataset. So a
random retention record is far more likely to be in-corpus than a random compound is —
and any evaluation that samples records rather than compounds inherits that flattery. This
is the same unit-of-resampling trap #15 §1 identified for pairs, at corpus scale.

**It also puts a number on how unrepresentative our anchor is.** The WSU-94 is **100%
in-corpus**; real chemistry is 9%. Our only measured-descriptor reference sits at the
extreme of the axis that most affects error, which is a coverage statement about #35's
figures far sharper than "94 neutral small solutes, median V 1.06".

## Limits

- **Tiers are structural, not error measurements.** Nothing here says how wrong the model is
  on the novel tier — that needs labels nobody has. It says where it is extrapolating.
- **The Murcko scaffold is a coarse instrument.** It ignores substitution entirely, so
  `scaffold_seen` is generous: a compound can share a benzene-ring scaffold with SoluteDB and
  be chemically remote. Read it as an upper bound on "structurally familiar".
- **RepoRT is a metabolomics corpus**, so its chemistry is not pharma chemistry. It is a
  better proxy for practice than the WSU-94, and it is not the target population either.
- Membership is skeleton-level (InChIKey first block), so tautomers, stereoisomers and salt
  forms count as seen — deliberate, since they are seen for a model's purposes.
