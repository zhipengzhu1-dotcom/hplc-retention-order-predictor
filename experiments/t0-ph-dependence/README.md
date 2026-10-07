# #25 — what it costs to hold `t0` fixed when it moves with pH

`python3 t0_bias.py` reproduces everything; `results.txt` is the captured run. numpy only.

## Verdict

**Not swamped by other errors — larger than any of them.** Over the sweep's own grid
(pH 2.5–6.5), holding `t0` fixed biases fitted `D` by **+0.42** log units at silanol pKa 6.0
and **+0.69** at pKa 4.5. That exceeds the entire prior SD on `D_run` (0.30), is 14× the LSER
lack-of-fit floor, and is 56× the difference #30 is commissioned to resolve.

**Decision (owner, 2026-08-16): `t0` is fitted from scouting runs**, with a `t0` marker
required in every scouting run, and carried as a named per-run error term until scouting data
exists. Pinned in `CONTEXT.md`.

## Why the error concentrates on the plateau that defines `D`

`k = (t_R − t0)/t0` is most sensitive to `t0` exactly when `k` is small — and the **ionised**
form of an acid is the smallest `k` in the experiment. So a `t0` error does not average out
across a pH sweep; it lands almost entirely on the ionised plateau.

Error in log₁₀ `k` from a −10% `t0` error:

| `k` | `t_R/t0` | error |
|---|---|---|
| 0.05 | 1.05 | **0.523** |
| 0.10 | 1.10 | 0.347 |
| 0.30 | 1.30 | 0.171 |
| 1.00 | 2.00 | 0.087 |
| 5.00 | 6.00 | 0.054 |
| 20.0 | 21.0 | 0.048 |

An ionised acid at `D` = 1.5 with `k_neutral` = 10 sits at `k` = 0.32 — in the steep region.

## Bias in `D`, over the operating range rather than the pH-11 extreme

`t0(pH)` modelled as a silanol-ionisation sigmoid between the reported 0.83 and 0.65 min,
with `t0` measured once at pH 2.5 and held:

| silanol pKa | `t0` drift 2.5→6.5 | `t0`-attributable bias in `D` |
|---|---|---|
| 4.5 | −21.3% | **+0.69** |
| 6.0 | −16.5% | **+0.42** |
| 7.5 | −2.0% | +0.03 |
| 9.0 | −0.1% | +0.00 |

**The ticket asked whether this is a pH-11 curiosity. It is not.** At silanol pKa 6.0 the
drift is already −16.5% by pH 6.5, inside the grid we intend to run. Whether it bites is
switched by a column parameter we do not know.

⚠ **Estimator caveat.** `fit_D_with_assumed_t0` uses the two-plateau ratio rather than the
protocol's full three-parameter fit. That estimator carries its own ≈ −0.07 bias, because an
acid of pKa 4.2 is 99.5% — not 100% — ionised at pH 6.5, so the upper plateau is not quite
reached. That artifact appears in **every** row regardless of `t0`, and is differenced out
against a per-condition-`t0` run; the table reports the `t0` part alone. The full fit would
recover `D` slightly better in absolute terms but does not change the differential result,
which is what the decision rests on.

## Against the rest of the budget

| | log units |
|---|---|
| **`t0` held fixed (silanol pKa 6.0)** | **0.416** |
| Prior SD on `D_run` (`prototype/thin-slice/slice.py`) | 0.300 |
| LSER lack-of-fit floor, MeOH (#15 §5) | 0.041 |
| LSER lack-of-fit floor, ACN (#15 §5) | 0.030 |
| Sweep's min detectable per-run Δ`D`, n = 8 | 0.0075 |

A bias this size does not sit quietly in a note. **It would be mistaken for the silanol
effect #30 exists to measure** — which is the specific way it could do real damage.

## The sweep is immune; the product is the exposed case

With a `t0` marker in every run, `t0` is measured per condition and the attributable bias
falls to zero. So #30 is safe by design, and it will directly measure `t0(pH)` on two columns.

But immunity is not an answer. **A user asking for a prediction at pH 8 has not injected a
`t0` marker at pH 8.** The sweep tells us how big the effect is; it does not decide what the
product does about it.

This inverts #25's previous status. `t0` is a **precondition** for #30's primary result, not
a by-product of it — the same conclusion `experiments/sweep-column-decision/` reached from
the other direction, where the per-column `t0` systematic is what floors the paired
difference.

## Why "model `t0(pH)`" was rejected

It is the physically satisfying option, and it needs a per-column silanol-ionisation
parameter. That is the **same axis #39 found unresolvable among C18 phases** on all three
available instruments — USP SRM-870 activity parameters, the HSM `C` terms, and WSU-2019's
Table S-1 electrostatic flags. Modelling would ask for a quantity the project has already
established it cannot obtain.
