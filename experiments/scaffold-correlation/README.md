# What the model's missing scaffold correlation actually costs

`env/model-venv/bin/python measure_impact.py`. Reuses the slice's own sampler and changes
exactly one thing.

## Why this was worth measuring

#37 measured that per-compound descriptor errors are **positively correlated within a scaffold
family** — ρ ≈ +0.42 (MeOH), +0.54 (ACN) propagated to `log k` — and warned that independent
draws therefore *"overstate Var(Δlog k) for scaffold pairs by up to ~2×… every P(A before B)
dragged toward 0.5."*

The thin slice draws compounds independently and says so in its own honesty list. That
overstatement is therefore live in **every number the slice produced**, including the variance
attribution the #30 sweep design was derived from. Before committing two weeks of instrument
time to a design derived from those numbers, it is worth knowing whether the gap is material.

## The change

Only the per-compound descriptor draw moves:

```
z_i = sqrt(ρ)·z_family + sqrt(1−ρ)·z_i_independent
```

This induces correlation ρ between compounds of a family while leaving each compound's own
4×4 covariance and marginal spread untouched. Everything else — the joint `Σ = D R D` draw of
the six system constants, the pKa shifts, `D`, the φ interpolation — is the slice's code.

Five families over the slice's 16 compounds, giving **16 within-family pairs out of 105**.

## Result: the mechanism is real, the impact is small

| operating point | ρ | confident pairs | Var same-family | Var diff-family |
|---|---|---|---|---|
| (0.30, pH 3.0) | **0.00** *(slice)* | 77/105 | 0.0982 | 0.2269 |
| | 0.42 | 78/105 | 0.0798 | 0.2290 |
| | 0.54 | 78/105 | **0.0744** | 0.2290 |
| (0.30, pH 7.0) | **0.00** | 74/105 | 0.1454 | 0.2611 |
| | 0.54 | 73/105 | **0.1212** | 0.2617 |
| (0.50, pH 5.0) | **0.00** | 67/105 | 0.0879 | 0.2171 |
| | 0.54 | 70/105 | **0.0766** | 0.2174 |

**The mechanism behaves exactly as #37 predicts.** Same-family `Var(Δlog k)` falls by
**17–24%** as ρ rises, and different-family variance is unmoved (0.2269 → 0.2290) — which is
the control working, since a bug here would have moved both.

**But the headline barely shifts: ±1 to 3 pairs out of 105.** Two reasons, and both are
structural rather than accidental:

1. **Within-family pairs are only 16 of 105.** A 24% variance reduction on 15% of the pairs
   cannot move an aggregate much.
2. **Descriptor error is a minority of `Var(Δlog k)` for the compounds that matter.** The
   slice's own variance attribution has, for partially ionised compounds, `pKa` = +0.67 against
   `desc` = +0.19; for fully ionised, `D` = +0.92 against `desc` = +0.07. Scaffold correlation
   is a *descriptor-error* effect, so it can only ever reach the smaller share.

### The "~2×" is not contradicted

#37's 2× applies to the **descriptor-error component alone** (`2σ²(1−ρ)` vs `2σ²`, which is
−54% at ρ = 0.54). Total `Var(Δlog k)` also carries pKa error and `D` error, neither of which
is scaffold-correlated here. Measuring −24% on the total is consistent with −54% on a
component that is under half the total.

## What this means for consolidation

**The variance attribution driving the sweep design is more robust than the honesty list
implied.** One of its six caveats is now measured rather than feared, and it is small. That is
a modest but real strengthening of the case for the experimental design.

**It also points where the model actually needs work.** The dominant terms are `pKa`
(0.67 for partially ionised) and `D` (0.92 for fully ionised) — not descriptors. So:

- Refining the descriptor layer further has low leverage on ionisables, which are the in-scope
  chemistry.
- The `pKa` term is currently **hardcoded literature values** with a ±0.35 per-compound error
  and prior `w_w`→`s_s` shifts (acids N(+0.30, 0.30), bases N(−0.15, 0.20)) that are
  literature-typical rather than fitted. **The slice never calls the ionisation provider that
  #20 already stood up** (`env/ionisation-venv`, Uni-pKa). Wiring it in would replace a
  hardcoded constant with a provider that emits per-microspecies populations and an ensemble
  spread — and it is what `spec/scenario-ensemble.md` requires anyway, since `f_neutral` must
  be stored per scenario.
- The `D` term is a prior by design and only measurement moves it. That is what the sweep is
  for, and this result leaves that argument intact.

**Recommended next consolidation step: wire the real ionisation provider into the model**, not
further descriptor work.

## Caveats

- ρ is applied as a **single value per family**, though #37 found it strongly family-dependent
  (phenols high, anilines ≈ 0). A per-family ρ would be more honest but needs more than five
  families to estimate.
- Families are assigned **by chemistry, by hand**. A descriptor-space kernel, as #37 suggests,
  would be less arbitrary.
- #37's ρ was measured on the WSU 94 — *"simple, mostly training-adjacent solutes; this gets
  the sign and shape, not the out-of-sample magnitude"*. Applying it to drug-like compounds is
  an extrapolation, and if the true ρ for drugs is higher the effect here is understated.
