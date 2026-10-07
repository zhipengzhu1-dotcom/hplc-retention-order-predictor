# Prototype findings — φ interpolation vs LSS reparameterisation (#24)

**Throwaway prototype.** Run: `python3 analyze.py` (numpy only; PCHIP implemented inline).
Full raw output in `results.txt`.

**What "error" means here (the key caveat, up front):** the reference value at every φ is
`log₁₀ k` computed from the *published* WSU constants at that φ, for all 94 Table S-2
compounds. This is **model-vs-model**: it measures the internal consistency of the two
routes against the per-φ LSER fits, not truth against measurement. The LSER's own residual
(SE ≈ 0.04–0.06 log k) and everything upstream is outside this comparison. WSU-2019 is
`log₁₀ k` throughout (CONTEXT.md); φ is handled as a fraction (0.10–0.70).

**Dewetting mask (#40) enforced:** only `incomplete_wetting = 0` rows are support
(354 − 22 masked = 332 fits, 51 (column, modifier) systems: 25 MeOH, 24 ACN, 2 THF).
Masks verified contiguous at the low-φ end for every system, so support always starts at
the lowest unmasked φ and no fit or interpolation crosses a dewetting floor.

## Route A — interpolate c,e,s,a,b,v across φ

Leave-one-φ-out over interior unmasked points; error = |Δlog₁₀ k| at the held-out φ
across all 94 compounds. Note LOO on a 10 %-step grid makes the interpolant span a 20 %
gap, so these numbers are **pessimistic** for production use on the full grid.
"Inside/outside" = the 30–60 % LSS validity window; outside interior points are all
φ = 20 % (endpoints are never held out).

| modifier | method | region | n | median | p90 |
|---|---|---|---|---|---|
| MeOH | linear | inside | 9 212 | 0.0140 | 0.0408 |
| MeOH | linear | outside | 1 786 | 0.0187 | 0.0690 |
| MeOH | pchip | inside | 9 212 | 0.0154 | 0.0513 |
| MeOH | pchip | outside | 1 786 | 0.0273 | 0.1182 |
| ACN | linear | inside | 8 742 | 0.0258 | 0.0632 |
| ACN | linear | outside | 1 128 | 0.0343 | 0.0891 |
| ACN | pchip | inside | 8 742 | 0.0152 | 0.0482 |
| ACN | pchip | outside | 1 128 | 0.0246 | 0.0927 |
| THF | linear | inside | 752 | 0.0321 | 0.0794 |
| THF | pchip | inside | 752 | 0.0197 | 0.0645 |

**Linear vs PCHIP is honestly mixed:** PCHIP clearly wins where curvature is strong
(ACN median 0.016 vs 0.027; THF 0.020 vs 0.032); linear is slightly better for MeOH
(0.0145 vs 0.0169) and safer at the MeOH 20 % edge (p90 0.069 vs 0.118). PCHIP is the
better single default because ACN is where the linear error is worst; a per-modifier
choice is available but every cell is under the 0.09 budget either way.

## Route B — LSS: log k = log k_w − S·φ per (compound, column, modifier)

Fitted by least squares on log k from the published constants at unmasked φ. Residuals
(in-sample for the points inside the fit range) split by region:

| modifier | fit range | region | n | median | p90 |
|---|---|---|---|---|---|
| MeOH | full unmasked | inside | 9 400 | 0.0185 | 0.0565 |
| MeOH | full unmasked | outside | 6 298 | 0.0234 | 0.0790 |
| MeOH | 30–60 % only | inside | 9 400 | 0.0095 | 0.0300 |
| MeOH | 30–60 % only | outside | 6 298 | 0.0441 | 0.1442 |
| ACN | full unmasked | inside | 8 930 | 0.0469 | 0.1239 |
| ACN | full unmasked | outside | 5 452 | 0.0767 | 0.1791 |
| ACN | 30–60 % only | inside | 8 930 | 0.0226 | 0.0474 |
| ACN | 30–60 % only | outside | 5 452 | **0.1314** | **0.3469** |
| THF | full unmasked | inside | 752 | 0.0450 | 0.1495 |
| THF | full unmasked | outside | 376 | 0.1111 | 0.2241 |
| THF | 30–60 % only | inside | 752 | 0.0383 | 0.0707 |
| THF | 30–60 % only | outside | 376 | **0.1724** | 0.3461 |

**Error growth outside the published window (window fit, the ticket's key question):**
median grows 4.6× (MeOH 0.0095→0.0441), 5.8× (ACN 0.0226→0.1314), 4.5×
(THF 0.0383→0.1724); p90 grows to 0.14 / 0.35 / 0.35. And "outside" here is still only
φ = 10–20 % and 70 % — *inside* the data range. Real gradients below 10 % or above 70 %
are extrapolation for both routes and are not measured here at all.
Fitting LSS on the full unmasked range instead spreads the damage: ACN residuals then
exceed the budget at p90 even *inside* the window (0.124).

## The comparison that decides

Per-compound QSPR error already in the budget: ≈ 0.09 log k (CONTEXT.md, conformal
stratum). Adding each route's typical (median) error in quadrature to 0.09:

- **Route A (pchip, inside):** 0.015–0.020 → budget 0.090 → 0.0915 (**+1.6–2.4 %, swamped**).
- **Route A worst cell** (MeOH pchip outside, p90 0.118) is the only place A even reaches
  the budget's scale — and it is a pessimistic 20 %-gap LOO at the 20 % edge.
- **Route B (window fit, outside window):** ACN 0.131 → 0.159 (**+77 %**), THF 0.172 →
  0.194 (**+116 %**), MeOH 0.044 → 0.100 (+11 %). **Not swamped.** For ACN and THF the
  LSS lack-of-fit *outside* 30–60 % is larger than the entire QSPR error.

Head-to-head at the identical interior-φ evaluation points (A = pchip LOO,
out-of-sample; B = full-range LSS, **in-sample** — the asymmetry favours B):

| modifier | region | n | A median / p90 | B median / p90 |
|---|---|---|---|---|
| MeOH | inside | 9 212 | 0.0154 / 0.0513 | 0.0184 / 0.0561 |
| MeOH | outside | 1 786 | 0.0273 / 0.1182 | 0.0118 / 0.0471 |
| ACN | inside | 8 742 | **0.0152 / 0.0482** | 0.0472 / 0.1245 |
| ACN | outside | 1 128 | 0.0246 / 0.0927 | 0.0190 / 0.0519 |
| THF | inside | 752 | **0.0197 / 0.0645** | 0.0450 / 0.1495 |

B's wins in the "outside" rows are an artefact of the asymmetry (B's fit *includes* the
evaluation point; A's excludes it and bridges a doubled gap). Even with that handicap, A
beats B by ~3× in ACN and ~2× in THF inside the window.

## Why: curvature of the constants in φ

LSS linearity in log k for *every* compound is equivalent to all six constants being
affine in φ. Measured deviation of interior points from the chord (first↔last unmasked
φ), in units of each point's own published sd:

| modifier | n per coef | c | e | s | a | b | v |
|---|---|---|---|---|---|---|---|
| MeOH (median sd; % > 2 sd) | 117 | 2.4; 56 % | 2.4; 56 % | 1.0; 22 % | 1.6; 40 % | **5.6; 84 %** | 2.2; 52 % |
| ACN | 105 | 4.5; 84 % | 4.9; 88 % | 3.4; 75 % | 1.5; 30 % | 3.9; 76 % | **6.8; 90 %** |
| THF | 8 | 3.7; 100 % | 6.0; 100 % | 2.8; 88 % | 0.8; 12 % | **8.5; 100 %** | 12.8; 100 % |

The curvature is real at the data's own precision — median mid-point deviations of 2–13
sd, worst in `v` and `b`, worst in ACN (the most-used modifier). LSS's premise fails
measurably; this is the mechanism behind Route B's numbers, and it is consistent with the
blueprint paper's own budget putting LSS lack-of-fit (9.1 % MAPE) as the largest term.

## Uncertainty shape (first quantification)

Empirical Pearson r of each coefficient between adjacent φ points, across columns
(6 φ-steps per modifier; 24–25 columns per step; THF n/a with 2 columns):

- MeOH: median r per coefficient 0.91–0.99 (c 0.93, e 0.97, s 0.99, a 0.95, b 0.91, v 0.93).
- ACN: 0.90–0.98 (c 0.95, e 0.90, s 0.97, a 0.98, b 0.91, v 0.93).

Coefficients at adjacent φ are near-collinear across columns, so uncertainty draws for
interpolated constants must be **correlated across φ** — a per-column curve draw, not
independent per-φ noise (independent draws would fabricate φ-to-φ jitter and scramble
predicted order between neighbouring compositions). Caveat: this cross-column r mostly
reflects genuine between-column variation, so it is an upper-bound-flavoured first
quantification of shared signal; the estimation-error correlation from shared calibration
solutes cannot be extracted from the published tables alone. LSS's 2-parameter
(log k_w, S) covariance is the cleaner object, but it buys that cleanliness with the
lack-of-fit above.

## Caveats and open ends

- **Model-vs-model** throughout (see top). Neither route is validated against measured
  retention here; that is #29/RepoRT territory.
- Route A LOO is pessimistic (20 % gaps); Route B in-range residuals are optimistic
  (in-sample). Both biases favour B, and A still wins where it matters.
- "Outside 30–60 %" here only reaches φ = 10–20 % and 70 %. Below the lowest unmasked φ
  and above 70 % there is **no support for either route**; per #40 the dewetting floor is
  a hard boundary, and this prototype clamps there rather than extrapolating.
- THF results rest on 2 columns (SunFire C18, XTerra MS C18) — directional only.
- The third arm from the scope-extension comment (linear vs quadratic vs Neue–Kuß, and
  the numerical-G coupling from #9) needs RepoRT and is blocked by #29 — not covered here.
- Data discrepancy worth flagging: the #40 comment on #24 says Synergi Polar-RP (ACN)
  "interpolates on 4 points starting at 40 %", but the CSV has 4 *published, unmasked*
  rows at 20–50 % (no masked rows for that system). Four points is right; the floor isn't.

## Recommendation (not a decision)

**Route A: interpolate the published constants across φ (PCHIP as default), per
(column, modifier), on unmasked support only, never extrapolating past the support
edges.** Its error is swamped by the 0.09 QSPR budget (inflates it < 3 % in quadrature);
LSS's is not — outside its own 30–60 % window its median error (ACN 0.131, THF 0.172,
n = 5 452 / 376) alone exceeds the entire per-compound QSPR budget, and real gradients
live outside that window. The curvature table shows the linear premise failing at 2–13 sd
on the coefficients themselves, so this is structure, not noise. Uncertainty cost:
carry correlated-across-φ coefficient draws per column (the adjacent-φ r of 0.90–0.99
says independent draws are wrong), which is more machinery than LSS's 2×2 covariance —
that is the price of not inheriting the blueprint's largest error term.
