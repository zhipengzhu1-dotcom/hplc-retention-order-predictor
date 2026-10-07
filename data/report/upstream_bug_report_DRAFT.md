<!--
DRAFT ONLY -- DO NOT SUBMIT without explicit sign-off.
Intended target: https://github.com/michaelwitting/RepoRT/issues (new issue)
-->

# `info.tsv` dataset name and `metadata.tsv` column.name are swapped for datasets 0310-0341

## Summary

For the 32 datasets `0310`-`0341` (the BGC/Harrieder-Witting temperature x
modifier study on the Restek Raptor Biphenyl and Waters CORTECS UPLC C18
columns), the free-text dataset `name` field in `info.tsv` names the wrong
column. `metadata.tsv`'s `column.name`, `column.usp.code` and
`column.particle.size` fields are internally consistent and correct;
`info.tsv`'s `name` field is not.

## Evidence

Pinned commit: `9de8d603377bbb6cc0f74e250eb7533fd8874df1` (2026-07-09).

| id range | `info.tsv` name says | `metadata.tsv` column.name says | USP code | particle size | real spec agrees with |
|---|---|---|---|---|---|
| 0310-0325 (16 datasets) | "..._Raptor Biphenyl" | Waters CORTECS UPLC C18 | L1 | 1.6 um | **column.name** (CORTECS UPLC C18 is 1.6 um, USP L1) |
| 0326-0341 (16 datasets) | "..._CORTECS C18" | Restek Raptor Biphenyl | L11 | 2.7 um | **column.name** (Raptor Biphenyl is 2.7 um, USP L11) |

Example rows:

```
0310  info.name = MSMS_pos_30_ACN_Raptor Biphenyl   metadata.column.name = Waters CORTECS UPLC C18   usp=L1   particle=1.6
0326  info.name = MSMS_pos_30_ACN_CORTECS C18        metadata.column.name = Restek Raptor Biphenyl    usp=L11  particle=2.7
```

The datasets immediately before (`0236`-`0259`, columns HSS T3 / XB-C18 / BEH
C18) and immediately after (`0342`-`0357`, HSS C18) this range show the
expected, non-swapped name<->column.name correspondence, so the bug is
confined exactly to `0310`-`0341` -- no fuzz at the boundary.

## Impact

Anyone who joins on the `info.tsv` `name` string (or extracts the column
label embedded in it, e.g. by string-splitting `"MSMS_pos_30_ACN_<column>"`)
to identify which physical column a dataset used will silently assign 32
datasets' retention data to the wrong column. This is a correctness bug on
exactly the axis (stationary-phase identity) that column-chemistry studies
care about, and it fails silently -- both labels look plausible, and the
wrong one is the more prominent field (it's in the filename-adjacent
`info.tsv`, not the wide `metadata.tsv`).

## Suggested fix

Swap the trailing column-label suffix of the `name` field in `info.tsv` for
ids `0310`-`0341`, using `metadata.tsv` `column.name` as ground truth (it is
corroborated by `column.usp.code` and `column.particle.size` against public
manufacturer specs, both of which are internally consistent through the
whole BGC study).

## Reported by

Filed while building split/loader infrastructure on top of RepoRT for a
downstream HPLC method-development project; happy to open a PR with the
one-line fix if useful.
