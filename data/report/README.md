# data/report — RepoRT loader, patch, and split infrastructure

Resolves #29.
Everything here is code + small deterministic definitions; **RepoRT itself is
not vendored**. Clone it and pin the commit yourself:

```bash
git clone https://github.com/michaelwitting/RepoRT.git /path/to/RepoRT
cd /path/to/RepoRT && git checkout 9de8d603377bbb6cc0f74e250eb7533fd8874df1
```

Then point every module below at `/path/to/RepoRT/processed_data`.

Requires `rdkit` (only for the scaffold split): `python3 -m pip install rdkit`.

## Layout

- `loader.py` — joins `info.tsv` + `metadata.tsv` per dataset, patches the
  `0310`-`0341` name/column-name swap, exposes `iter_datasets()` and
  `load_rtdata()`.
- `stats.py` — row counts, missing-field rates, duplicate/anomaly counts.
  `python -m data.report.stats /path/to/RepoRT/processed_data`
- `factorials/bgc.py` — the BGC (Harrieder/Witting) temperature x modifier
  factorial: cells, per-cell compound sets, shared cores.
  `python -m data.report.factorials.bgc /path/to/RepoRT/processed_data`
- `factorials/acn_meoh_pairs.py` — clean same-column/same-T/same-pH
  ACN-vs-MeOH dataset pairs, deduplicated 1:1 by max compound overlap.
  `python -m data.report.factorials.acn_meoh_pairs /path/to/RepoRT/processed_data`
- `splits/scaffold_split.py` — deterministic Murcko-scaffold train/val/test
  split over the 18,872 non-SMRT compounds (hash-bucketed, no materialized
  file needed — same input always gives the same split).
  `python -m data.report.splits.scaffold_split /path/to/RepoRT/processed_data`
- `splits/unseen_column_split.py` — column support census, USP-code
  L1-fraction quantification, and the list of columns with genuinely
  different retention chemistry and enough data to hold out.
  `python -m data.report.splits.unseen_column_split /path/to/RepoRT/processed_data`
- `upstream_bug_report_DRAFT.md` — drafted, **not submitted**, bug report for
  the RepoRT maintainers describing the name-swap bug.
- `write-up.md` — the full write-up: pinned commit, the swap bug, BGC
  factorial results (with a correction to issue #4's framing), ACN/MeOH
  pairs (with an honesty note on a number that could not be reproduced
  exactly), unseen-column-split weakness quantified, missing-field rates,
  and new anomalies found while building this.

## Why no committed data files

RepoRT is CC BY-SA 4.0 (ShareAlike) and its own repository is the canonical,
versioned source — pinning a commit SHA is sufficient for reproducibility
without vendoring ~180k rows into this repo. All splits are either
deterministic functions of a stable key (scaffold split) or short id lists
small enough to inline in code/markdown (BGC cells, ACN/MeOH pairs, unseen
column list) — see `findings.md` for those id lists in full.
