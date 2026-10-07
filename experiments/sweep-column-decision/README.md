# #30's column decision — power and generalisability, weighted equally

`python3 decide_columns.py` reproduces everything; `results.txt` is the captured run.
Deterministic seed. numpy only.

## Verdict

**Power is not the constraint — it is not close.** The design resolves a per-run `D`
difference roughly **forty times finer** than the prior it would update. Generalisability
is therefore the only axis on which these candidates differ, exactly as suspected.

**Recommended design: SunFire C18 versus Apex II C18** — an *asymmetric* pair, one WSU
column and one not. It covers 94% of the silica-C18 population, needs **no extrapolation**
to reach the type-A regime, and — the point that decides it — **costs nothing on the
question the sweep exists to answer.**

## A. Power

The protocol fits `k_obs(pH) = k_neutral · [f + (1 − f)·10^(−D)]` with `k_neutral`, `D`
and apparent pKa free, over 5 pH levels × 2 replicates. Monte-Carlo with 0.3% retention
RSD and a **1% per-column `t0` systematic** (one draw per campaign, shared across all
compounds on that column, so it never averages down):

| `k_neutral` | `k_ion` | `tR_ion / t0` | SD(`D`) |
|---|---|---|---|
| 2 | 0.063 | 1.063 | **0.070** |
| 5 | 0.158 | 1.158 | 0.026 |
| 10 | 0.316 | 1.316 | 0.014 |
| 20 | 0.632 | 1.632 | 0.007 |
| 40 | 1.265 | 2.265 | 0.003 |

Paired column-A-minus-column-B difference over tier-1 acids:

| n acids | SD(mean Δ`D`) | min detectable (95%) |
|---|---|---|
| 4 | 0.0049 | 0.0096 |
| **8** | **0.0039** | **0.0075** |
| 12 | 0.0040 | 0.0078 |

Against a prior of `D_run ~ N(0, 0.30)` (from `prototype/thin-slice/slice.py`), the
minimum detectable difference of **0.008** is ~1/40 of the prior SD. The experiment is
massively overpowered on precision, and **n beyond ~8 acids buys nothing** — the paired
difference is already floored by the two per-column `t0` systematics, which do not shrink
with more compounds.

### One actionable protocol change falls out

`D`'s precision is governed by **how retained the ionised form still is**, not by the
neutral form. At `k_neutral` = 2 the ionised peak sits at 1.06 × `t0` and the `t0` error
swamps it — SD(`D`) = 0.070, a twenty-fold degradation against `k_neutral` = 40.

So **tier 1 must be run where the acids are well retained**: prefer the 30% ACN condition
over 50%, and choose tier-1 acids with high neutral retention. Ibuprofen and naproxen are
fine; a weakly retained acid contributes almost nothing. This costs nothing to adopt and
protects the tier that carries the entire silanol result.

It also means the **`t0` marker is not a by-product**. Uracil/thiourea in every run is
load-bearing for `D` itself, because the ionised plateau sits close to `t0`. That
tightens the case for #25 as well.

## B. Generalisability

Reference population: 393 silica C18 columns (`USPtype` L1, `phase` C18) in the HSM set.

| percentile | `C(7.0)` |
|---|---|
| p5 | −0.11 |
| p25 | 0.08 |
| p50 | 0.22 |
| p75 | 0.67 |
| p95 | 1.48 |
| p99 | 2.19 |

Type-A columns: n = 75, median `C(7.0)` = **1.050**.

| design | span `C(7.0)` | Δ`C` | percentile span | % population | baseline |
|---|---|---|---|---|---|
| Booked *(superseded)* | −0.05 to −0.04 | 0.01 | p9–p10 | **1.3%** | both WSU |
| WSU-internal | −0.09 to 0.30 | 0.39 | p7–p62 | 55.2% | both WSU |
| **Asymmetric** | **−0.09 to 2.69** | **2.78** | **p7–p100** | **93.6%** | SunFire WSU + non-WSU |
| Outside WSU | −0.17 to 2.69 | 2.86 | p2–p100 | 98.0% | neither WSU |

Two columns fit `D_run` at two points on the `C(7.0)` axis. Anything between them is
interpolation; anything beyond is extrapolation:

| design | to type-A median (1.05) | to p99 (2.19) |
|---|---|---|
| Booked *(superseded)* | 109× span beyond | 223× span beyond |
| WSU-internal | 1.9× span beyond | 4.8× span beyond |
| **Asymmetric** | **inside span** | **inside span** |
| Outside WSU | inside span | inside span |

The booked pair is not merely weak — at 109× its own span it is not measuring the axis at
all. The WSU-internal pair is respectable on coverage (55% of the population, p7–p62) but
still **extrapolates ~2× beyond its calibrated range** to say anything about a type-A
column, and ~5× to reach p99. And two points fit a line without testing whether the
relationship *is* a line, so that extrapolation is unfalsifiable from within the
experiment.

## Why the asymmetric design costs nothing — the thing that decides it

The apparent trade was "fixed neutral baseline versus real silanol contrast". **For the
tier that matters, that trade does not exist.**

The protocol's own analysis plan says tier 1 (acids, pKa 3–5) gives "`D` measured on both
columns, **no model input**" — `k_neutral` is a *free parameter* identified from the data,
because acids reach both plateaus inside pH 2.5–8. WSU membership supplies a fixed
`k_neutral` only where a plateau is unreachable, which is **tier 2b** (strong bases,
needing pH ≈ 11.5) and **tier 3** (neutral LSER controls).

Tier 1 is the tier that carries the per-run silanol term. So putting a non-WSU column on
one arm forfeits nothing on the silanol question. It costs only the tier-2b and tier-3
arms on that column — and those can run on the WSU arm, which is exactly the
evidence-class separation the protocol already uses for tier 2b on Betasil.

### The recommended pair

| | `C(7.0)` | `H` | type | in WSU-2019 |
|---|---|---|---|---|
| **SunFire C18** | −0.09 | 1.03 | B | **yes** |
| **Apex II C18** | 2.69 | 1.00 | A | no |

Δ`C(7.0)` = **2.78**, with hydrophobicity matched to **0.03** — so the pair contrasts on
silanol activity and not on retention strength. SunFire carries the WSU system constants
for tiers 2b and 3; Apex II supplies the type-A end that WSU-2019 does not contain at all.

## Caveats

- **`Apex II C18` availability is unverified**, and it is an older line. It dominates
  every high-`C` pair in the database, so if it cannot be sourced the high end needs
  re-picking — that is the one thing that could still overturn this recommendation.
- **The HSM `C` values themselves are single-condition** (50% ACN, pH 2.8/7.0, 35 °C) and
  their provenance is unverified against Shackman 2016. `C(7.0)` is being used here as an
  *ordering* variable, which is the weakest claim that supports the argument.
- **Power numbers assume the noise model above.** The 1% `t0` figure is achievable with a
  marker in every run, but it is an assumption, so it was varied:

  | `t0` uncertainty | min detectable Δ`D` (95%) | prior SD ÷ that |
  |---|---|---|
  | 0.5% | 0.0046 | 66× |
  | 1% | 0.0074 | 41× |
  | 3% | 0.0235 | 13× |

  Even at 3% the design resolves thirteen times finer than the prior it updates, so the
  "power is not the constraint" conclusion does not depend on this assumption.
- **Two points cannot test linearity** in any design. If `D_run` versus `C(7.0)` is
  materially non-linear, a third column is the only fix; the asymmetric design at least
  brackets the range rather than extrapolating across it.
