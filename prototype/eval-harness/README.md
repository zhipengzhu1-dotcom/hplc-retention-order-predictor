# Evaluation harness for #45 — built before any measured data exists

`python3 test_harness.py` — 50 checks, all passing, no third-party dependencies
beyond numpy.

## Why it exists now rather than later

#15 pre-registered the acceptance thresholds so they could not be chosen to be passed.
**Writing the code that computes them after seeing the data would reopen the same door.**
A threshold is only half a commitment while the statistic that tests it is still
unwritten — the remaining freedom simply moves from the number to the method. This module
closes that half.

There is also a plain operational reason. The pH sweep will deliver ~230 injections of
hard-won data in a few weeks. That is the wrong moment to begin writing the evaluation.

## The load-bearing decision: `Figure`

`CONTEXT.md` requires that the harness have *"no code path returning an accuracy figure
unaccompanied by its applicable baseline, its `n`, its interval and its condition axis —
enforced by the function signature and asserted in a test, not left to reviewer
discipline."*

So **no function here returns a float.** They return `Figure`, which:

- cannot be constructed without its six mandatory fields — #27's four (value, axis,
  baseline, `n`, interval) plus #31's two (reference identity, reference kind);
- **raises `BareNumberError` on `float()` and on any format spec**, so a bare number cannot
  leak into a report even by accident;
- renders the *whole claim* from `str()`, so even a careless `print()` is compliant;
- picks its verb from the reference kind — **"accurate to"** only for measurement,
  **"agrees with"** for model output, with the reference named in the same sentence. A
  figure referenced to Absolv-lineage descriptors can never claim accuracy, because
  agreeing with that corpus partly means mimicking the tool we are trying to beat (#5).

`n_compounds` is the field name deliberately. There is no way to pass a pair count.

## What is implemented, and what it is checked against

| Piece | Pre-registration | Checked by |
|---|---|---|
| Compound-level bootstrap | §1 — the compound is the unit of resampling | wider than the naive binomial on synthetic data (**1.28×** at 15 compounds; §1 predicts ~1.35×) |
| Order reliability, ECE, reliability table | §1 Metric A | finite, bounded, populated |
| Conformal inflation `λ` | §3 — the primary calibration threshold | **recovers known under-dispersion**: injected 1.0 / 1.5 / 2.4 → measured 0.99 / 1.48 / 2.47 |
| `λ` bands and their boundaries | §3 | 2.0 is calibrated (`≤`), 2.5 is degraded not rescope |
| Coverage vs sharpness | §3 — "coverage is not the test" | absurdly wide intervals give **coverage > 0.99 while `λ` < 0.6** — the test we cannot fail, failing to fail |
| Relative error reduction | §2 — not percentage points | recovers an injected 17% |
| Three verdicts | §2 — *inconclusive* named in advance | inconclusive never reads as pass, and **states the required `n`** |
| Strata merge rule | §4 | merges upward on confidence tertile; **ionisation regime never merges**; hard floor `n` ≥ 19 |

### The two checks worth reading

**Coverage cannot fail, and the harness demonstrates it.** Given deliberately absurd
interval widths, coverage reports > 0.99 — a clean pass — while `λ` reports < 0.6, correctly
saying the width was never needed. That is §3's argument made executable rather than
asserted.

**Inconclusive states its own remedy.** A gain of +0.02 against a 2 SD margin of 0.05
renders as:

> INCONCLUSIVE on pH: gain +0.02 against a 2 SD compound-level margin of 0.05 (n=15
> compounds); this axis could not be resolved at the available sample size — n=94 compounds
> would settle it

The required `n` is computed from the margin's 1/√n scaling, not left to the reader. §2's
failure mode is discovering inconclusiveness afterwards, when it is indistinguishable from a
pass in a summary.

## Status and what is deliberately not here

This is **prototype-stage** and lives beside the other prototypes because the project has no
production tree yet. It is not throwaway — it is intended to graduate.

Not yet implemented, because each needs a decision or data that does not exist:

- **Wiring to `prototype/thin-slice`'s ensemble.** The harness takes arrays, not the
  `retention[P, S, n]` object, deliberately: #44 item 1 is still amending that interface, and
  binding to it now would bake in a shape that is about to change.
- **The ESS contract** (#16/#17, #44 item 4). Exercising it needs the reconciliation path,
  not just scoring.
- **Both #27 baselines** as computable functions. Cold (Crippen `logP`) and
  condition-transfer are passed in as arrays here; computing them needs the compound set.
- **The applicability-domain firewall** (#40) — only in-domain predictions may enter any
  metric, and opt-in residuals are never evaluation data. The harness does not yet refuse
  out-of-domain rows because the domain object lives upstream.

Each is a real gap rather than an oversight, and none blocks the harness being exercised on
the first data that arrives.
