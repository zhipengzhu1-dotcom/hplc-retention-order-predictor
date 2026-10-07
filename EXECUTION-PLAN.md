# Repo alignment review and execution plan

> **Status 2026-08-16: Phases 1–3 are complete** (commits `0ba4a3d`, `10ccace`, `0b6158b`).
> Phase 4 is the standing-blocks register and needs no action. Outcomes are recorded inline
> below; two items changed on contact with the work and say so.

Review date **2026-08-16**, against `ec56c7d`. Scope: every tracked module (140 files,
excluding `papers/` PDFs which are gitignored by design). Method: read the governing
documents, then check each module's claims against the artefacts that produced them and
against later findings that may have superseded them.

**Headline: the repo is internally consistent on its architecture and green on its tests.
What has drifted is bookkeeping — an extraction defect that silently costs one compound in
three analyses, an uncommitted re-run whose numbers differ from the committed findings, and
several documents still hedging about facts that have since been measured.** None of it
invalidates a conclusion. All of it is cheap to fix, and two items would quietly corrupt
future work if left.

Test state, all green today:

| suite | result |
|---|---|
| `prototype/eval-harness/test_harness.py` | All checks passed (50) |
| `spec/test_ensemble.py` | All checks passed |
| `model/test_applicability.py` | 17/17 |

---

## Module review

### `CONTEXT.md` — the ubiquitous language (883 lines)

The governing document, and it is doing its job: every convention carries the ticket that
pinned it, and the conventions are specific enough to be checkable rather than aspirational.

Two items are now stale:

- Line 626 cites ensemble under-dispersion at **median 1.21**, from #42's run. That run was
  computed at **n = 93**, not 94, because of the extraction defect below. Recomputed at
  n = 94 it is **1.23 [0.98, 1.47]**. The change is immaterial to the argument; the wrong
  `n` on a figure that feeds a pre-registered threshold is not.
- The `eE` entry warns that `E` "is plausibly **biased**, and propagating variance does not
  fix bias; #35 is the actual remedy." #35 has now measured it: **bias +0.035 [+0.004,
  +0.065]**, concentrated almost entirely in weak bases (**+0.152 [+0.077, +0.236]** versus
  neutrals). The convention should carry the measurement rather than the promise.

The closing note says numerical conventions are to be "record[ed] here in full as the owning
ticket closes". The RepoRT pinned commit lives only in `data/report/loader.py`, and the
corpus-coverage tiers from `experiments/solutedb-census/` have no entry at all.

### `spec/` — the interfaces (5 files)

`ensemble.py` + `scenario-ensemble.md` specify the prediction layer's output and validate
conformance; `representation-limits.md` records what the representation cannot express;
`applicability-domain.md` + `model/applicability.py` are new this session.

**The one real gap is an integration seam.** `spec/ensemble.py`'s `refused()` infers refusal
from all-NaN retention — refusal as *data absence*. `model/applicability.py` produces a
`Decision` carrying the reason and its warrant. Nothing connects them, so **a refusal's
reason is destroyed by the time it reaches the ensemble**, which is precisely what #33/#40
require be preserved. Similarly `stratum[n]` is a MUST in the spec with **no producer
anywhere in the repo**; the domain object's tiers and corpus tier are the natural inputs.

`representation-limits.md:68` still says the WSU-94 "are almost certainly inside SoluteDB".
That is now measured: **94 of 94**, 93 on SoluteDB's own fixed-H InChI key and 1 on
InChIKey skeleton.

### `model/` — the model layer (12 files)

`ionisation.py` + `calibrate_ionisation.py` (Uni-pKa wired in, provider error SD 0.170
against literature, unbiased), `columns.py`, `dispersion.py`, `sampler.py`, `sequence.py`,
`voi_D.py`, and now `applicability.py` + its 17 checks.

**`CONTEXT.md` asserts a convention this module does not implement.** The corrector fitting
order — "#35's WSU re-fit → solute-vector correction on its residual → `D` last. **Never
fitted jointly** … asserted in a test" — has **no test anywhere in the repo**. It is the
convention most likely to be violated silently, because all three correctors can absorb the
same discrepancy, and it is the one with no enforcement.

### `prototype/` — the measurement prototypes (32 files)

`thin-slice`, `descriptor-covariance`, `confidence-stratum`, `phi-vs-lss`, `eval-harness`,
`soluteml-bias`. Quality is high and each carries a FINDINGS.md that states what was *not*
computable before what was.

**`prototype/thin-slice/` has uncommitted working-tree drift**, and it matters because the
committed findings quote these artefacts:

| quantity | committed | working tree |
|---|---|---|
| confident pairs | 77 | 78 |
| ambiguous | 30 | 31 |
| coin-flips | 5 | 4 |
| ambiguous agreement | 0.367 | 0.355 |
| coin-flip agreement | 0.60 | 0.75 |
| `max_abs_dP` | 0.118 | 0.128 |
| descriptor variance share | 0.199 | 0.204 |

`retention_ensemble.npz` and `sources/usp-column-db/usp-approach-columns.xlsx` are also
modified. Someone re-ran the slice; nothing records why, and **there is no way to tell which
run the published claims refer to**. The numbers are close enough that no conclusion moves,
which is exactly why this rots quietly.

`eval-harness` is the strongest module in the repo — 50 checks, `Figure` making a bare number
unrepresentable — and its README honestly lists its own four gaps. One of them (the domain
object) is now fillable.

### `experiments/` — designed studies (28 files)

`report-overlap` is the most load-bearing: it establishes RepoRT has zero isocratic datasets,
no dwell-volume field in its schema, and an imputed `column.t0`. That verdict now blocks
three separate proposals, and it earned its keep this session by killing one of mine before
I ran it.

`calibration-preregistration.md` pre-registers #15's thresholds. It cites the **1.21** figure
at the wrong `n`, same as `CONTEXT.md`.

`ph-sweep-protocol.md` is the designed-not-committed sweep #45 waits on. Worth noting against
this session's work: **every drug in it lands in the `degraded` tier** of the new domain
object, so its headline will be a degraded-stratum figure. That is a finding about the sweep's
reporting, not a defect in it.

### `sources/` — provenance-controlled extracts (26 files)

`wsu-lser` holds the anchor. **`extract_wsu.py` is the producer of the defect described
below** — it writes Table S-2's section heading into the first compound name of the
weak-base block. Provenance discipline elsewhere is good: every extract has a README naming
its source document and page.

### `data/report/` — RepoRT infrastructure (6 files)

Pinned (`9de8d60`), patched for the documented 0310–0341 name-swap bug, with scaffold and
unseen-column splits and a drafted-but-unsubmitted upstream bug report. Healthy, and
deliberately vendors no data. Used unmodified for this session's census.

### `research/` — background (8 files)

Written before the measurements existed, and several claims have since been measured but not
updated. `qspr-approaches.md:305` — "SoluteML's shipped variance is epistemic-only and almost
certainly under-confident" — is now a number: **RMSE / median ensemble SD is 3.95 (E), 3.70
(S), 3.85 (A), 2.15 (B)**, and 1.23 in retention space. These are background documents, so
the fix is a pointer to the measurement, not a rewrite.

### `env/` — environments (10 files)

`QSPR_ENVIRONMENT.md`, `install_qspr.sh`, `smoke_test.py`, patches, and
`SOLUTEML_COVARIANCE_RESULTS.md`.

⚠ **Correction to this review, 2026-08-16.** The first draft called that last file a
duplicate of `prototype/descriptor-covariance/FINDINGS.md`. It is not: it is a **9-compound
environment demonstration** for #38, while the prototype holds the 94-compound measurement,
and `prototype/soluteml-bias/` holds the error covariance against measured values — three
different quantities. The real (small) defect was that nothing said so, which is why they
looked interchangeable. Fixed by a scope note rather than a merge.

### Local hygiene (not tracked, so cosmetic)

`.claude/` holds **16 stale agent worktrees totalling 7.6 GB**, each with its own venv. It is
gitignored, so this is disk, not correctness.

---

## Cross-cutting findings, ranked

### A. The S-2 extraction defect — highest priority

`sources/wsu-lser/wsu2019-descriptors-s2.csv` row 81 reads
`Compounds designated as weak bases 2-Aminobiphenyl`: the table's section heading is glued
onto the compound name.

| file | effect |
|---|---|
| `prototype/confidence-stratum/measure_confidence_stratum.py` | **drops the compound** — #42's results are n = 93 |
| `prototype/confidence-stratum/measure_rank_stability.py` | **drops the compound** |
| `experiments/structural-deficiency/detectability.py` | **drops the compound** |
| `experiments/fs-calibration/calibrate_fs.py` | keeps the row, mislabels the name |
| `prototype/phi-vs-lss/analyze.py` | keeps the row, mislabels the name |
| `prototype/descriptor-covariance/*`, `prototype/thin-slice/slice.py` | strip the heading — correct |
| `prototype/soluteml-bias/measure_bias.py` | normalises at load — correct |

Two of the three losses were reported honestly as "93 of 94" — the authors saw the count and
believed it. **The dropped compound is 2-aminobiphenyl: a weak base carrying the third-worst
`E` error (+0.448)**, so it sits in the tail of the effect §5 of `soluteml-bias/FINDINGS.md`
reports. The defect was biasing the stratum most worth measuring.

Fixing it at the producer fixes all seven consumers at once, and the heading also marks the
15-compound weak-base block that #43 wants as a `class` column anyway.

### B. Uncommitted drift in `prototype/thin-slice/`

Committed findings and working-tree artefacts disagree (table above). Either the re-run is
authoritative and the findings need regenerating, or it is a stray execution and should be
reverted. Both are one command; leaving it undecided is the only bad option.

### C. Measured facts still described as suspicions

`spec/representation-limits.md:68`, `prototype/confidence-stratum/FINDINGS.md:20` and
`research/qspr-approaches.md:305` hedge about corpus contamination and ensemble
under-confidence. Both are now measured (94/94 in-corpus; 2.15–3.95× under-dispersion).
Hedged language about a settled fact invites someone to re-litigate it.

### D. The 1.21 figure at the wrong `n`

`CONTEXT.md:626` and `experiments/calibration-preregistration.md:161`. Should be
**1.23 [0.98, 1.47], n = 94**. Small change, but one of the two documents pre-registers
acceptance thresholds.

### E. The domain object is not wired to anything

Built this session, referenced by nothing. Two seams: the eval harness (which names the gap
in its own README) and `spec/ensemble.py` (whose `refused()` discards the reason). Also the
unimplemented `stratum[n]` producer.

### F. `#36`'s corrector fitting order has no test

`CONTEXT.md` says "asserted in a test". There is none. The convention exists precisely
because the failure is silent.

---

## Execution plan

Ordered so that each step is independently committable and nothing depends on a later step.
Effort is rough and assumes no new data.

### Phase 1 — stop the bleeding (half a day)

**1.1 Fix the S-2 extract at the producer.** Amend `sources/wsu-lser/extract_wsu.py` to strip
the section heading from the compound name and emit a `class` column (`neutral` |
`weak_base`) recording what the heading conveyed. Regenerate
`wsu2019-descriptors-s2.csv`. Update `sources/wsu-lser/README.md` to describe the new column
and record the defect and its date in the extract's provenance notes.
*Check:* all 94 names match `prototype/descriptor-covariance/wsu94_smiles.csv` exactly.

**1.2 Re-run the three analyses that were dropping a compound** —
`measure_confidence_stratum.py`, `measure_rank_stability.py`, `detectability.py` — and update
their FINDINGS with `n = 94`. Expect small movements; report them even if nothing changes,
since "the number did not move" is the useful outcome.
*Check:* each prints 94.

**1.3 Simplify the consumers that hand-strip the heading.** Once 1.1 lands, the
`.replace("Compounds designated as weak bases ", "")` in `measure1_ensemble_covariance.py`
and `slice.py`, and the load-time normalisation in `measure_bias.py`, become dead code that
future readers will copy. Delete them and point at the `class` column.

**1.4 Resolve the thin-slice drift.** Decide whether the working-tree re-run is authoritative.
If yes, regenerate `FINDINGS.md` from it and commit both together with a note on what changed
and why. If no, `git checkout` the three files.
*Check:* `git status` clean, and every number in `thin-slice/FINDINGS.md` traceable to the
committed `results.json`.

### Phase 2 — align the documents with what was measured (half a day)

**2.1 `CONTEXT.md`:** update the under-dispersion figure to 1.23 at n = 94; replace the `eE`
entry's "#35 is the actual remedy" with the measured bias and its weak-base concentration;
add the RepoRT pin and the corpus-coverage tiers to the numerical-conventions index.

**2.2** Replace the hedges in `spec/representation-limits.md`,
`prototype/confidence-stratum/FINDINGS.md` and `research/qspr-approaches.md` with the measured
facts and a pointer to `prototype/soluteml-bias/FINDINGS.md`.

**2.3** Correct the `1.21` in `experiments/calibration-preregistration.md`, with an explicit
note that the figure moved because of a defect fix and not because a threshold was retuned —
that distinction is the whole point of pre-registration.

**2.4** ~~Pick one home for the SoluteML covariance results.~~ **Withdrawn** — they are three
distinct measurements, not one duplicated. Replaced by a scope note in
`env/SOLUTEML_COVARIANCE_RESULTS.md` naming what it is and what to quote instead.

### Phase 3 — close the integration seams (one to two days)

**3.1 Wire `evaluation_strata()` into the eval harness.** The harness gains a domain-aware
entry point that refuses out-of-domain rows and keeps degraded rows in their own stratum.
This touches the pre-registered scoring path, so it lands as its own commit with tests
asserting that refused rows cannot enter a metric and that degraded and in-envelope rows
never pool.
*Owner decision required before starting* — the semantics are described in
`spec/applicability-domain.md`.

**3.2 Carry the refusal reason into `ScenarioEnsemble`.** Replace inference-from-NaN with an
explicit per-compound `Decision`, so `refused()` can say *why*. Keeps the NaN convention as
the data representation; adds the reason alongside it.

**3.3 Implement the `stratum[n]` producer**, the spec's one unimplemented MUST, cut on the
method-independent structural summary per #42 §4 and #15 §4's merge rule. The domain tiers
and corpus tier are the inputs.

**3.4 Write the `#36` corrector-order test.** Assert the three correctors are fitted in the
pinned order and never jointly, and that uncertainty is reported from the final stage only.

### Phase 4 — the standing blocks (no action possible)

Recorded so they are not rediscovered:

| item | blocked on | note |
|---|---|---|
| #45 retention run | ~230 injections + Nürenberg dwell recovery | RepoRT cannot substitute |
| #35 output-layer re-fit | Poole 2020 | inconclusive at n = 94, not failed |
| #43 steric predicate | WSU-2019 ref 29 | coverage half shipped; mechanism half labelled `assumed` |
| #42 clean-split questions | a corpus pairing out-of-corpus compounds with a usable label | no public source examined has both |
| joint-tail coupling in the residual sampler | more compounds | Gaussian copula has zero tail dependence by construction |

### Not scheduled

Clearing the 7.6 GB of stale worktrees under `.claude/` is disk hygiene, gitignored, and
touches nothing correctness-related. Worth doing, not worth planning.
