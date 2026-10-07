# RepoRT: pinned, cleaned, split — write-up

Resolves #29,
blocked by and building directly on #4.
Everything below is measured by running the code in `data/report/` against a
pinned clone, not quoted from RepoRT's papers or from #4's prose without
re-checking.

## 1. Pin

```
commit 9de8d603377bbb6cc0f74e250eb7533fd8874df1
Author: Fleming Kretschmer <fleming.kretschmer@uni-jena.de>
Date:   2026-07-09
Merge pull request #245 from michaelwitting/inconsistencies-0052
(Duplicate rows in dataset 0052)
```

This is the exact commit #4 pinned. Re-cloned on 2026-08-14 and confirmed it
is **still the tip of RepoRT's default branch** — no drift in the five weeks
since #4. Recorded in `data/report/loader.py::PINNED_COMMIT`. 421 datasets,
183,819 `rtdata_canonical_success` rows total, 80,285 of them SMRT (datasets
`0186`/`0209`) — both match #4's numbers exactly, confirming the clone and
loader are reading the same snapshot.

## 2. The name-swap bug — confirmed, exact extent

**Confirmed real**, verified directly (not assumed from the ticket). For
datasets `0310`–`0341` (32 datasets, exact range, no fuzz at the boundary):

- `0310`–`0325` (16 datasets): `info.tsv` `name` says "..._Raptor Biphenyl"
  but `metadata.tsv` `column.name` is **Waters CORTECS UPLC C18** (USP `L1`,
  particle 1.6 µm — matches the real CORTECS spec).
- `0326`–`0341` (16 datasets): `info.tsv` `name` says "..._CORTECS C18" but
  `metadata.tsv` `column.name` is **Restek Raptor Biphenyl** (USP `L11`,
  particle 2.7 µm — matches the real Raptor Biphenyl spec).
- Checked the neighboring BGC-study id ranges (`0236`–`0259` and
  `0342`–`0357`) with the same method: name and column.name agree there. The
  bug is confined to exactly these 32 ids.
- **Only** the `info.tsv` `name` string's embedded column label is wrong.
  `column.name`, `column.usp.code`, `column.particle.size`, and every other
  `metadata.tsv` field are correct and internally consistent on both sides
  of the swap.

Fix: `data/report/loader.py::patched_dataset_name()` — trust
`metadata.tsv column.name`, not the `info.tsv` name string, for these 32
ids. `iter_datasets()` applies this automatically; `Dataset.name_swapped`
flags affected rows.

Upstream bug report drafted (**not submitted**) at
`data/report/upstream_bug_report_DRAFT.md`.

## 3. BGC factorial — extracted, with a correction to #4's framing

The BGC (Harrieder/Witting) study covers **6 columns × T{30,40,50}**, but
only **3 of those 6 columns were run in both modifiers** — this was verified
directly, not assumed. The other 3 have no MeOH run anywhere in the
repository under this study.

| | ACN | MeOH |
|---|---|---|
| Waters ACQUITY UPLC HSS T3 | yes (3 T) | **no MeOH data exists** |
| Phenomenex Kinetex XB-C18 | yes (3 T) | **no MeOH data exists** |
| Waters ACQUITY UPLC BEH C18 | yes (3 T) | **no MeOH data exists** |
| Waters CORTECS UPLC C18 | yes (3 T) | yes (3 T) |
| Restek Raptor Biphenyl | yes (3 T) | yes (3 T) |
| Waters ACQUITY UPLC HSS C18 | yes (3 T) | yes (3 T) |

So the design is **27 condition cells, not 36**: 18 ACN cells (6 col × 3T) +
9 MeOH cells (3 col × 3T). Each cell is backed by a pos-mode and neg-mode
MSMS dataset (CCS pos/neg exist at T=40 only and were excluded — a different
measurement technique, not a retention replicate). Per-cell compound count
(union of pos+neg): 374–643.

**Shared cores** (intersection of `inchikey.std` across cells):

- 6-column × 3T, ACN only (18 cells): **203 compounds**
- 3-column × 3T × {ACN,MeOH} full crossed sub-design (18 cells): **229 compounds**
- All 27 cells at once: **180 compounds**

#4's "~500 compound shared core / 493–499 per column" turns out to describe
a narrower quantity: one column's own cells intersected with each other
(not cross-column). Recomputed per column:

| column | own-cells only | modifier crossing |
|---|---|---|
| Waters ACQUITY UPLC HSS T3 | 577 | ACN only |
| Phenomenex Kinetex XB-C18 | 589 | ACN only |
| Waters ACQUITY UPLC BEH C18 | 504 | ACN only |
| Waters CORTECS UPLC C18 | 286 | ACN+MeOH |
| Restek Raptor Biphenyl | 302 | ACN+MeOH |
| Waters ACQUITY UPLC HSS C18 | 274 | ACN+MeOH |

The 3 ACN-only columns land right at "~500" (504–589); the 3 modifier-crossed
columns are lower (274–302) because they're intersected over twice as many
cells (6 vs 3). Use `data/report/factorials/bgc.py` to regenerate either
number depending on which factorial you're actually running (a pure-T study
on 6 columns has ~200–590 compounds available depending on scope; a
modifier-crossed study is limited to 3 columns and ~230–300 compounds).

## 4. ACN-vs-MeOH pairs — could not reproduce "156" exactly

Defined "clean pair" precisely and reproducibly: same column (swap-corrected),
same temperature (both populated, equal), same pH (both `eluent.A.pH` and
`eluent.B.pH` populated, equal), one ACN dataset and one MeOH dataset,
matched 1:1 (no dataset reused) by greedy maximum compound overlap.

**Result: 43 deduplicated pairs, best overlap 406 compounds** (matches #4's
"~400" figure almost exactly). A looser, non-deduplicated many-to-many cross
product over the same column/T/pH match gives 110 pairs; #4 reported 156.
Neither of my two reasonable readings reproduces 156 exactly, and #4's
comment doesn't document its exact matching algorithm — I could not
reverse-engineer it with confidence, so I'm reporting my own precisely
documented method and number rather than forcing a match. **The qualitative
finding is solid either way**: many clean same-condition pairs exist, and
the best ones overlap on ~400 compounds. Regenerate with
`data/report/factorials/acn_meoh_pairs.py`.

Full 43-pair list (id, id, column, T, compound overlap):

| ACN id | MeOH id | column | T (°C) | overlap |
|---|---|---|---|---|
| 0314 | 0322 | Waters CORTECS UPLC C18 | 50 | 406 |
| 0330 | 0338 | Restek Raptor Biphenyl | 50 | 396 |
| 0346 | 0354 | Waters ACQUITY UPLC HSS C18 | 50 | 395 |
| 0312 | 0320 | Waters CORTECS UPLC C18 | 40 | 397 |
| 0328 | 0336 | Restek Raptor Biphenyl | 40 | 393 |
| 0316 | 0324 | Waters CORTECS UPLC C18 | 40 | 376 |
| 0332 | 0340 | Restek Raptor Biphenyl | 40 | 346 |
| 0344 | 0352 | Waters ACQUITY UPLC HSS C18 | 40 | 346 |
| 0348 | 0356 | Waters ACQUITY UPLC HSS C18 | 40 | 345 |
| 0327 | 0335 | Restek Raptor Biphenyl | 30 | 335 |
| 0311 | 0319 | Waters CORTECS UPLC C18 | 30 | 281 |
| 0343 | 0351 | Waters ACQUITY UPLC HSS C18 | 30 | 280 |
| 0313 | 0321 | Waters CORTECS UPLC C18 | 40 | 243 |
| 0329 | 0337 | Restek Raptor Biphenyl | 40 | 246 |
| 0331 | 0339 | Restek Raptor Biphenyl | 50 | 233 |
| 0315 | 0323 | Waters CORTECS UPLC C18 | 50 | 216 |
| 0347 | 0355 | Waters ACQUITY UPLC HSS C18 | 50 | 215 |
| 0345 | 0353 | Waters ACQUITY UPLC HSS C18 | 40 | 203 |
| 0282 | 0283 | Phenomenex Kinetex HILIC | 40 | 72 |
| 0236 | 0002 | Waters ACQUITY UPLC HSS T3 | 30 | 69 |
| 0274 | 0275 | Phenomenex Kinetex PS C18 | 40 | 75 |
| 0286 | 0285 | Phenomenex Kinetex PS C18 | 40 | 75 |
| 0284 | 0287 | Phenomenex Kinetex PS C18 | 40 | 74 |
| 0276 | 0277 | Phenomenex Kinetex PS C18 | 40 | 73 |
| 0280 | 0281 | Phenomenex Kinetex PS C18 | 40 | 73 |
| 0292 | 0293 | Phenomenex Kinetex PS C18 | 40 | 73 |
| 0278 | 0279 | Phenomenex Kinetex PS C18 | 40 | 67 |
| 0270 | 0271 | Phenomenex Kinetex EVO C18 | 40 | 52 |
| 0290 | 0289 | Phenomenex Kinetex EVO C18 | 40 | 50 |
| 0272 | 0273 | Phenomenex Kinetex EVO C18 | 40 | 44 |
| 0268 | 0269 | Phenomenex Kinetex EVO C18 | 40 | 42 |
| 0288 | 0291 | Phenomenex Kinetex EVO C18 | 40 | 36 |
| 0294 | 0295 | Phenomenex Kinetex EVO C18 | 40 | 29 |
| 0317 | 0325 | Waters CORTECS UPLC C18 | 40 | 15 |
| 0342 | 0350 | Waters ACQUITY UPLC HSS C18 | 30 | 13 |
| 0326 | 0334 | Restek Raptor Biphenyl | 30 | 9 |
| 0238 | 0001 | Waters ACQUITY UPLC HSS T3 | 40 | 8 |
| 0310 | 0318 | Waters CORTECS UPLC C18 | 30 | 7 |
| 0333 | 0341 | Restek Raptor Biphenyl | 40 | 6 |
| 0349 | 0357 | Waters ACQUITY UPLC HSS C18 | 40 | 6 |
| 0358 | 0362 | Waters ACQUITY UPLC BEH C18 | 35 | 1 |
| 0360 | 0361 | Waters ACQUITY UPLC BEH C18 | 35 | 0 |

⚠ Caveat noticed while building this: the matching criterion (column+T+pH
only, no study/instrument constraint) let two pairs through that are cross-study
matches, not same-experiment ACN/MeOH replicates: `0236`/`0002` (BGC HSS T3
ACN vs. an unrelated FEM_long MeOH dataset that happens to share column/T/pH)
and similar low-overlap tail entries. These are the ones with tiny overlap
(≤15) at the bottom of the table — they are real matches on the stated
criteria but not "designed pairs." The 18 BGC modifier-crossed pairs (rows 1–18
above, overlap 203–406) are the trustworthy core of this set for a
modifier-transfer benchmark; treat the Kinetex EVO/PS/HILIC pairs (overlap
29–75, these are the Kruve Lab pH-sweep study, already flagged by #4's
addendum) and the tiny-overlap tail as secondary evidence only.

## 5. Split infrastructure

### Scaffold split — feasible, deterministic, verified

100% SMILES coverage confirmed on all 18,872 non-SMRT compounds (matches
#4). Murcko scaffold computed with RDKit, split assigned by hashing the
scaffold string (deterministic, no leakage possible by construction — same
scaffold always gets the same label). Ran it:

```
total compounds: 18872
  train: 15457 (81.9%), 6831 distinct scaffolds
  val:    1604 ( 8.5%),   902 distinct scaffolds
  test:   1811 ( 9.6%),   892 distinct scaffolds
scaffolds appearing in >1 split (should be 0): 0
```

Code: `data/report/splits/scaffold_split.py`. No file needs to be committed —
`assign_split(scaffold)` is a pure function.

### Unseen-column split — feasibility confirmed, weakness quantified

Non-SMRT dataset census (n=419):

- 410/419 (97.9%) carry a `column.usp.code`.
- **At the dataset level: 306/410 coded datasets (74.6%) are USP L1.** (#4
  said "308 of 412" — same order of magnitude; small discrepancy not chased
  further, see note below.)
- At the distinct-column level: **38 of 58 named, coded columns (65.5%) are
  L1.**
- **9 datasets (1,785 compounds) carry no column name and no USP code at
  all** (`0013 0014 0033 0034 0058 0059 0210 0218 0230`, all PredRet-sourced)
  — unusable for any column-axis split. New anomaly, not previously flagged.

**Quantified weakness**: #4's own independently-measured pairwise
elution-order reversal rate between two USP-L1 (C18) columns is only
**1.0–3.9%**, vs. 4.6–6.2% for a biphenyl-vs-C18 contrast. Since 65–75% of
the repository's columns/datasets are L1, an "unseen column" cross-validation
that averages over many folds will be dominated by L1-vs-L1 folds — i.e. by
the near-trivial 1–4% case — and will report a headline number that looks
like strong generalization when it mostly measures noise-floor stability on
chemically near-identical phases.

**Columns that offer genuinely different chemistry AND enough compounds
(≥100) to be a usable held-out test set** — the corrected, complete list
(#4 named only one example, Raptor Biphenyl):

| column | USP code | compounds | note |
|---|---|---|---|
| Waters ACQUITY UPLC BEH Amide | L68 | 1,198 | HILIC-mode |
| Restek Raptor Biphenyl | L11 | 770 | RP, different selectivity |
| Phenomenex Synergi Polar-RP | L11 | 748 | RP, polar-embedded |
| Merck SeQuant ZIC-pHILIC | L122 | 684 | HILIC-mode, zwitterionic |
| Merck SeQuant ZIC-HILIC | L114 | 647 | HILIC-mode, zwitterionic |
| Waters ACQUITY UPLC BEH C8 | L7 | 197 | RP, shorter alkyl chain |
| Ascentis Express F5 (PFP) | L43 | 194 | RP, pentafluorophenyl |
| Ascentis Express ES-Cyano | L10 | 185 | RP, cyano |
| Ascentis Express Phenyl-Hexyl | L11 | 184 | RP, phenyl-hexyl |
| Waters XBridge BEH Amide | L68 | 160 | HILIC-mode |
| Thermo Accucore HILIC | L3 | 152 | HILIC-mode, bare silica |
| Waters ACQUITY UPLC BEH HILIC | L3 | 123 | HILIC-mode, bare silica |

Of these, the **5 HILIC-mode phases** (BEH Amide, ZIC-pHILIC, ZIC-HILIC,
XBridge BEH Amide, BEH HILIC) are structurally guaranteed to generalize
differently — HILIC retention order is not a monotonic transform of RP order
— and are the strongest "unseen chemistry" test available. The PFP/cyano/
phenyl-hexyl/C8 phases are still reversed-phase but meaningfully different
selectivity. Recommend: report unseen-column results as **named individual
adversarial cases** (biphenyl, one HILIC phase, etc.), never as one averaged
"unseen-column" number — an average would be silently dominated by the
65–75% L1 majority. Code: `data/report/splits/unseen_column_split.py`.

## 6. Row counts, missing-field rates, and anomalies

Non-SMRT denominator, n=419 datasets (matches #4's "excluding SMRT" framing):

| field | populated | rate |
|---|---|---|
| `column.temperature` | 248/419 | **59.2%** (matches #4's "59%") |
| `eluent.A.pH` (incl. `0` sentinel) | 371/419 | 88.5% |
| `eluent.A.pH` (nonzero only) | 319/419 | 76.1% |
| `column.usp.code` | 410/419 | 97.9% |
| `column.name` | 410/419 | 97.9% |
| `column.length` | 391/419 | 93.3% |
| `column.flowrate` | 378/419 | 90.2% |

Row counts: 421 datasets total, 183,819 `rtdata_canonical_success` rows
(80,285 SMRT + 103,534 non-SMRT), **18,872 unique non-SMRT compounds** by
`inchikey.std` (matches #4 exactly).

### New anomalies found (beyond what #4 flagged)

1. **9 datasets with no column identity at all** — see §5 above. 1,785
   compounds effectively unusable for any column-axis analysis.
2. **`inchikey.std` collapses stereochemistry, and it matters here.** Within
   a single dataset, **183 of 419 datasets** have at least one `inchikey.std`
   shared by rows with *different* `rt` values — 11,405 such extra rows
   total. This is almost certainly distinct stereoisomers/positional isomers
   sharing a standardized key, each genuinely eluting differently (not a
   data-entry error). **Any pipeline that dedupes or joins on `inchikey.std`
   alone — including every extraction in this ticket — will silently
   conflate them.** A further 729 rows are exact duplicates (same key, same
   rt) and are pure redundancy, consistent with #4's noted "duplicate rows
   in dataset 0052" fix already merged upstream at this pin. Flagging this
   explicitly since it's a real caveat on the BGC/pairs/scaffold-split
   numbers above: they're all inchikey-keyed, and a few percent of "one
   compound" is actually "two isomers with the same flattened key."
3. **BGC factorial is not the 6×3×2 design #4's prose implies** — see §3.
   Real design is 6×3 (T only) plus a 3×3×2 sub-design (T and modifier).
   This matters directly for #11/#12: any model claiming to have learned a
   6-column modifier effect only has evidence from 3 of those columns.
4. **43 vs. 156 ACN/MeOH pairs** — see §4. Flagged as an open discrepancy,
   not resolved.

### Could not verify / out of scope for this ticket

- Whether #4's exact "156" pair count and "308/412" column count come from a
  documented, reproducible method — their source comment doesn't specify one
  and it wasn't reverse-engineerable with confidence. My own numbers (43/110
  pairs; 306/410 or 38/58 L1 columns) are documented precisely enough to be
  re-run and re-checked by anyone.
- RDKit's Murcko scaffold definition is the standard one; no alternative
  scaffold definitions were compared.
