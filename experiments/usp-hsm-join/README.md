# #39's join — USP SRM-870 activity parameters against the HSM ion-exchange terms

Run 2026-08-15. `python3 join_usp_hsm.py` (stdlib only) regenerates everything here:
`matched.csv`, `unmatched.txt`, `results.txt`.

> **Public copy:** `matched.csv` is not included because it carries USP column parameters,
> which are not redistributed (see `sources/usp-column-db/README.md`). The per-column USP
> values in `results.txt` are redacted for the same reason; the correlation statistics are
> kept. To regenerate both, build `sources/usp-column-db/usp-approach-columns.csv` and
> `sources/wsu-lser/` as their READMEs describe, then rerun the script.

**The question**, from #39:
do the USP SRM-870 column-activity measures (`Hy, CTF, CFA, TFA, BD`, 145 columns) correlate
with the HSM cation-exchange terms `C(2.8)`/`C(7.0)` (819 columns)? If they do, we have a
validated silanol-activity axis and the anionic/cationic asymmetry becomes decidable per
column. If they do not, #33's wholesale refusal stands on firmer ground.

## Answer in one line

**They correlate, in the direction the SRM 870 provenance predicts — but the association is
carried almost entirely by a handful of old type-A silica columns, so it does not rank the
modern low-activity phases the project will actually use.**

> #39 is **closed** on this work. The full resolution — including the anionic half, which
> this join did not touch — is in the ticket and pinned in `CONTEXT.md` under *charge sign
> is not a retention class*, *the honest prior degenerates to refusal*, and *no instrument
> resolves silanol activity among C18 phases*.

## The join

95 of 145 USP columns matched, 41 unmatched, **9 excluded as ambiguous**.

Matching is mechanical in three tiers — exact normalised, spaceless, and particle-size-
stripped (USP writes `ProntoSil 120-5-C18-H` where HSM writes `ProntoSIL 120 C18 H`) — plus
a 13-entry hand-verified alias table in the script for house-style differences
(`BetaBasic 18` → `Hypersil Beta Basic-18`).

**A name-containment tier was tried and removed.** It matched `Discovery C18 WP` to
`Discovery C18` and `Alltima HP C18 AQ` to `Alltima HP C18` — different columns in the same
product line. Where a USP name could plausibly be a different product (bare `Allure`,
`Ultra`, `Pinnacle DB` against a dozen phase-suffixed HSM rows) the pair is **declined and
left unmatched**. The 41 unmatched are the price of that; several (`SiliaChrom *`,
`Prestige`, `Superiorex`, `Reliachrom`) are simply absent from HSM.

Two guards, both of which fired:

- **Duplicate HSM rows under one name** are excluded, not resolved by taking the first.
  `Betasil C18` carries two different parameter sets (`C28` −0.03/`C70` −0.04 and `C28`
  0.095/`C70` 0.099), as does `Hypersil ODS`. Picking one would invent a number.
- **Collisions** — two USP rows reaching one HSM row — drop both. This caught the
  `Purospher STAR RP18e` 3 µm/5 µm pair and, aptly, the two `TSKgel ODS-100V` rows already
  known to be defective in USP's own application.

## Correlations (n as shown, permutation p, 20 000 shuffles, seed 39)

| USP | HSM | n | Pearson | p | Spearman | p |
|---|---|---|---|---|---|---|
| `Hy` | `C28` | 94 | −0.199 | 0.056 | −0.314 | 0.003 |
| `Hy` | `C70` | 94 | −0.330 | 0.001 | −0.373 | 0.000 |
| `CTF` | `C28` | 89 | 0.286 | 0.018 | 0.236 | 0.026 |
| `CTF` | `C70` | 89 | 0.224 | 0.035 | 0.374 | 0.000 |
| `CFA` | `C28` | 90 | 0.456 | 0.000 | 0.178 | 0.092 |
| `CFA` | `C70` | 90 | **0.613** | 0.000 | 0.229 | 0.028 |
| `TFA` | `C28` | 94 | **0.562** | 0.000 | 0.350 | 0.001 |
| `TFA` | `C70` | 94 | 0.494 | 0.000 | 0.393 | 0.000 |
| `BD` | `C28` | 93 | 0.042 | 0.689 | −0.063 | 0.551 |
| `BD` | `C70` | 93 | −0.182 | 0.083 | −0.251 | 0.017 |

Every sign is the one the confirmed parameter expansions predict: the tailing-factor-like
activity measures (`CTF`, `CFA`, `TFA`) rise with ion-exchange activity, hydrophobicity-like
`Hy` falls against it, and bonding density `BD` — which is a coverage measure, not an
activity measure — is the one that shows essentially nothing at pH 2.8.

**So the SRM 870 parameters are measuring something real about silanol activity.** That much
is now evidence rather than inference from provenance.

## Why it still does not carry a per-column prior

The Pearson/Spearman gap on `CFA` (0.613 vs 0.229) says the association is not monotone
across the set — it is leverage from a few extreme columns. The six highest `C(7.0)` columns
in the matched set are `Resolve C18`, `Supelcosil LC-18`, `Hypersil PAH`, `Platinum EPS
C18`, `Spherisorb ODS-1`, `Spherisorb ODS-2` — a roll-call of classic high-silanol type-A
phases. Drop those six:

| | full | six dropped |
|---|---|---|
| `CFA`–`C28` Pearson | 0.456 | **0.121** |
| `CFA`–`C28` Spearman | 0.178 | **0.047** |
| `CFA`–`C70` Pearson | 0.613 | 0.443 |
| `TFA`–`C28` Pearson | 0.562 | 0.310 |
| `Hy`–`C70` Pearson | −0.330 | −0.318 |

`CFA` against `C28` collapses to nothing. What the SRM 870 parameters do is **separate old
high-activity silica from everything else** — a distinction the HSM `type` field already
makes directly and more robustly (type A median `C70` 1.050, n = 75; type B median 0.207,
n = 516). They do not resolve activity differences *among* modern low-activity phases, which
is the regime in which a per-column prior would have to operate.

`Hy` and `TFA` degrade least, so if any USP parameter is used as a proxy it should be those,
not `CFA` — and only ever to flag the type-A tail.

## Consequence for #30's sweep pair — this supersedes the USP-based choice

#30 picked **XBridge
Shield RP18 versus Betasil C18** on the USP parameters (`CFA` 2.4 vs 7.6, `TFA` 1.1 vs 2.0),
expecting a contrast in silanol activity so `D`'s per-run component is estimable as a paired
difference. #39 reserved the right to supersede that on the HSM `C` term. **It does.**

| Column | `H` | `C(2.8)` | `C(7.0)` | HSM `type` |
|---|---|---|---|---|
| XBridge Shield RP18 | 0.83 | −0.12 | −0.05 | EP (embedded polar) |
| Betasil C18 | 1.05 | −0.03 | −0.04 | B |
| Betasil C18 *(duplicate row)* | 1.063 | 0.095 | 0.099 | B |

**On the ion-exchange axis these two columns are the same column.** ΔC(7.0) is 0.01, or 0.15
if the other `Betasil` row is the right one — against a range of −1.00 to +2.69 across the
silica C18 columns, and −1.36 to +7.51 across all 819. Both sit at the bottom of the silanol
scale — XBridge Shield is embedded-polar precisely to suppress silanol interaction, and
Betasil is type-B silica. The USP `CFA` contrast of 2.4 vs 7.6 is real in `CFA` but does not
correspond to an ion-exchange contrast, which is exactly the leverage problem above showing
up in the one pair the project cares about.

Run as booked, the sweep would measure `D` on two low-activity columns and find little
difference — and that result would be **uninformative rather than negative**. It would not
license the conclusion "silanol conditioning is unwarranted", because the design never
varied silanol activity.

Note also that `Betasil C18`'s two conflicting HSM rows differ by more than the pair's whole
ΔC. The sweep column's own parameters are not settled.

### Candidate replacement pairs

Silica C18 only (`USPtype` L1, `phase` C18 — this excludes `apHera C18 Polymer`, which tops
the unrestricted ranking at `C70` 7.51 but is polymer-on-alumina chemistry, not silanol
activity), matched on hydrophobicity at |Δ`H`| ≤ 0.05 so the pair contrasts on `C` and not
on everything:

| Δ`C(7.0)` | low | high |
|---|---|---|
| 3.18 | Ace 5 C18-PFP (−0.995, `H` 0.899, type B) | Cogent hQ C18 (2.18, `H` 0.90, type B) |
| 3.02 | Inertsil ODS-3 (−0.33, `H` 0.99, type B) | Apex II C18 (2.69, `H` 1.00, type A) |
| 2.86 | Luna C18(2) (−0.17, `H` 1.00, type B) | Apex II C18 (2.69, `H` 1.00, type A) |
| 2.82 | Prodigy ODS(3) (−0.13, `H` 1.02, type B) | Apex II C18 (2.69, `H` 1.00, type A) |

**Luna C18(2) versus Apex II C18** is the pick to put to a chromatographer: Δ`C(7.0)` 2.86
against the booked pair's 0.01–0.15, hydrophobicity matched to 0.01, and a
type-B/type-A split that is the silanol distinction stated in the source's own vocabulary.
`Ace 5 C18-PFP` heads the list on ΔC but is a PFP phase, so it varies more than silanol
activity. **Commercial availability of these columns has not been checked** — `Apex II` in
particular is an older line — and that check is the next step before #30 is rebooked.

## The larger finding: WSU-2019 has no high-`C` column in it

Checking the replacement pair against WSU-2019 turned up a structural fact about the
evidence base. Of the 25 WSU columns, 23 are in HSM (`Chromolith Performance RP-18e` and
`Discovery HS C18` are not), and:

- **Not one of them is type A.** Zero high-activity silica columns in the entire anchor set.
- **The C18 members span `C(7.0)` −0.09 to 0.30**, against −1.00 to 2.69 across HSM's silica
  C18 columns. They occupy a narrow slice at the bottom.
- The only genuinely high-`C` members are fluorinated (`Kinetex F5` 2.66, `Discovery HS F5`
  1.42, `Fluophase-RP` 1.41) — a different retention mechanism, and the same phases
  #40 flagged for electrostatic behaviour. The highest non-fluorinated is `Synergi Polar-RP`
  at 0.77, a polar-embedded phase.

**Every column for which we hold fitted LSER system constants is a low-silanol-activity
column.** So independently of whether the USP parameters can rank anything, the anchor set
carries no observations from the regime where a cationic compound would actually be
retained. A silanol-conditioned prior would extrapolate off the end of the evidence in the
one direction that matters, and because the `C` term is unparameterised by the solvation
parameter model (#36), nothing in the LSER layer would flag that it was doing so.

Same shape as Ulrich's argument recorded in #39, one axis over: no fitted coefficient, no
descriptor referent, **and no anchor column in the regime**.

This also costs #30 directly — see the superseded-pair block in
[`../ph-sweep-protocol.md`](../ph-sweep-protocol.md).

## Third instrument, same verdict: the WSU-2019 electrostatic flags

WSU-2019 Table S-1 grades its electrostatic caveat as either "weak" or unqualified (#40).
That grading is an **independent expert judgement** rather than another measurement, so it
is the sharpest available check on whether *anything* ranks silanol activity among C18
phases. Against HSM `C(7.0)`:

| | weak | unqualified |
|---|---|---|
| all phases | n = 6, mean `C70` +0.176 | n = 15, mean +0.435 |
| **C18 only** | **n = 4, mean +0.072** | **n = 5, mean +0.050** |

Across all phases the flags look like they track `C(7.0)`. **Among C18 columns the ordering
vanishes** — the weak group is if anything marginally higher. The apparent agreement is
carried entirely by the three fluorinated phases (`Kinetex F5` 2.66, `Discovery HS F5` 1.42,
`Fluophase-RP` 1.41).

That is the *same shape* as the USP–HSM leverage result above, from a completely different
kind of source. Three instruments, one verdict: **they agree on separating extreme phases,
and none of them ranks conventional C18 columns.** The `n` here is small and no significance
is claimed — but it is unambiguously no *support* for conditioning.

## What this does and does not settle in #39

- **Settled:** the SRM 870 parameters carry genuine silanol-activity signal. The USP dataset
  is not noise, and the confirmed expansions are consistent with the measurements.
- **Settled:** a **silanol-activity-conditioned prior is not honest on USP parameters** at
  the resolution the ticket wanted. The correlation that would justify it lives in the
  type-A tail only. #33's refusal stands, now for a measured reason rather than a
  precautionary one.
- **Settled, on a second and independent ground:** cationic species stay fully unplaceable.
  Beyond the missing mechanism, there is no characterised column in our anchor set where
  that mechanism is even active, so there is no route to calibrating what "retained" means.
- **Open:** whether the HSM `type` field plus `C(7.0)` *itself* — skipping USP entirely — is
  a defensible conditioning axis. This join was about validating USP against HSM; it
  incidentally shows HSM is the better instrument, which is a different question the ticket
  should now ask directly. Note the WSU coverage gap above constrains the answer.
- **Exit condition, now concrete:** what would reopen silanol conditioning is LSER system
  constants fitted on at least one genuinely high-`C(7.0)` alkyl-silica column. That is a
  measurement, not a dataset search — WSU-2019 is the binding constraint, so no reanalysis
  of what we already hold will produce it.
- **Settled (resolving the ticket):** the **anionic** near-`t0` prior is refused too, so the
  asymmetry #39 was opened to adjudicate does not exist. "Permanently anionic" is not a
  retention class — its members span the chromatogram, from methanesulfonate at `t0` to
  dodecylbenzenesulfonate strongly retained on C18 — because retention is governed by the
  hydrophobic surface area of the rest of the molecule, not by the charge. A near-`t0` prior
  collapses `k_neutral · 10^(−D)` into a claim about `D` alone, asserting `k_neutral` is
  small *and uniform*; the second is false. An honest prior would span the whole run, which
  ranks exactly as "unplaceable" does while looking like knowledge. Note also that the stack
  has **no anion-exclusion parameter at all** — HSM's `C` is *cation* exchange — so the
  ticket's premise that the anionic half was the tractable half does not hold.

## Caveats on this analysis

- **`Fs`/HSM parameter provenance is unverified** — Shackman 2016 asserts the published HSM
  tables are internally inconsistent and has not been read. The duplicate rows found here
  are consistent with that and are the first direct evidence of it in our own data.
- **The 41 unmatched columns are not missing at random** — they skew toward columns absent
  from HSM entirely, and toward product lines whose USP naming is too coarse to
  disambiguate. A correlation computed on 95 of 145 is not a correlation on the USP set.
- **p-values are permutation-based and uncorrected for the 10 tests** run. At 10 tests the
  `BD`–`C28` null result and the marginal `Hy`–`C28`/`CTF`–`C70` values should not be read
  individually; the pattern of signs is the finding, not any single p.
- **HSM parameters are measured at one condition** (50% ACN, pH 2.8, 35 °C) with `C` alone
  at a second pH, so this join says nothing about φ dependence — which is precisely the
  dependence #34 flagged for the `C` term.
