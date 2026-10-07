# Protocol — pH sweep to fit the ionisation retention drop `D`

Status: **DESIGNED, NOT COMMITTED.** This is wayfinding output, not a booking. The
design is internally consistent and its constraints are traced to evidence, but the
compound plate is not frozen, no columns are ordered, and no instrument time is
reserved. Treat every number below as *what the design would require if run*, not as a
commitment that it will be. Written under
#30 as a spec deliverable,
so that whoever runs it inherits a design rather than an intention.

The #45 gate below **has opened** (2026-08-16): the bounded public-data search returned a
reportable negative, so the sweep proceeds on the variance argument alone and is now the
project's **only** route to measured retention data. The column pair was re-selected on the
HSM `C(7.0)` term after #39's join, and the design was checked for both power and
generalisability (`experiments/sweep-column-decision/`).

**Procurement was checked, not performed** (2026-08-16): both columns are current products and
both pH ratings accommodate their arm of the grid — Kinetex EVO C18 at pH 1–12 against a 2.5–10.5 ask, and
Supelcosil LC-18 at pH 2–7.5 against a 2.5–6.5 ask. Nothing analytical remains; what is left is
ordering columns and booking instrument time.

## Purpose, in priority order

**Primary — fit `D`.** `CONTEXT.md` defines the ionisation retention drop as
`k_ion = k_neutral · 10^(−D)`, empirical rather than physical, absorbing the mechanisms the
solvation parameter model never parameterised (silanol ion exchange, ion-pairing,
electrostatics). #21's variance attribution measured `D` as owning **0.78–0.92** of the
variance for fully ionised compounds — the largest single term in the project — and it has
**never been fitted to anything**. Its current value is a literature-range prior of roughly
1–2.

**Secondary, from the same runs:**

| Output | Ticket |
|---|---|
| `D`'s per-run (column silanol activity) part, as the two-column difference | #34, #36, #39 |
| Elution-order reversal rate on ordinary hardware — is #4's 23.5% real or phase-specific? | #4, #30 |
| Retention data for calibration, on compounds overlapping the thin slice | #45, #15 |
| The LSER+QSPR chain graded on drug-like chemistry, via neutral controls | #31, #35 |
| Is `t0` pH-dependent? | #25 |
| The `s_s` pH-scale correction, currently a named uncorrected error term | #20 |

⚠ The reversal rate was this ticket's *original* purpose and is now secondary. Designing for
it would buy the wrong data — a reversal either happens between two pH values or does not, so
it needs breadth and two levels, whereas `D` is the plateau of a sigmoid and needs depth.

## Design

**Temperature 45 °C throughout.** Non-negotiable: the WSU-2019 system constants are a 45 °C
dataset (#11 pins it as the nominal condition), and the neutral baseline is only directly
comparable at that temperature.

**Ionic strength held constant across every pH.** `D` carries ion exchange, which is
ionic-strength dependent, so a varying buffer concentration would confound the very parameter
being fitted. Fix total buffer at a single value (20 mM suggested) and change only the buffer
identity as pH requires.

### Columns — two, of contrasting silanol activity

**Selected on the HSM `C(7.0)` term**, which #36 measured as 80.5% of `Fs²` and which is the
ion-exchange axis the sweep exists to probe. Decision and full analysis in
[`experiments/sweep-column-decision/`](sweep-column-decision/README.md).

| Column | `C(7.0)` | `H` | HSM `type` | Role |
|---|---|---|---|---|
| **Kinetex EVO C18** | **−0.01** | 1.008 | B | Low-silanol arm. **In WSU-2019**, so the neutral baseline is *fixed* for tiers 2b and 3. Organo-silica hybrid, **pH 1–12**, which is what keeps tier 2b alive. |
| **Supelcosil LC-18** | **1.75** | 1.01 | **A** | High-silanol arm. Type-A silica — the regime WSU-2019 contains **no** example of. Supelco/Merck, current product. |

Δ`C(7.0)` = **1.76**, with `H` matched to **0.002** — so the pair contrasts on silanol
activity and *not* on retention strength. That control matters more than the size of the
contrast: a pair differing in both would not identify which one moved `D`.

**Why an asymmetric pair (one WSU column, one not).** Tier 1 below fits `k_neutral` from the
data — "no model input" — because acids reach both plateaus inside the pH grid. WSU membership
supplies a fixed `k_neutral` only where a plateau is unreachable, i.e. tiers 2b and 3. Tier 1
is the tier carrying the per-run silanol term, so a non-WSU high-silanol arm forfeits
**nothing** on the question the sweep is for. Tiers 2b and 3 run on the EVO arm.

**Generalisability.** The pair spans 86% of the 393 silica-C18 columns in the HSM database, and
the **type-A median `C(7.0)` of 1.050 falls inside the span** rather than beyond it — so
reading `D_run` at a high-silanol column is interpolation, not extrapolation. The superseded
pair below spanned 1.3% and would have extrapolated 109× beyond its own range.

### Third arm: XBridge C18 as a control on the axis itself (added 2026-08-16)

| Column | `C(7.0)` | `H` | type | in WSU | role |
|---|---|---|---|---|---|
| **Kinetex EVO C18** | −0.01 | 1.01 | B | yes | low-silanol arm |
| **XBridge C18** | **+0.13** | **1.00** | B | **yes** | **control — see below** |
| **Supelcosil LC-18** | +1.75 | 1.01 | A | no | high-silanol arm |

`H` is matched to **0.01 across all three**, so the hydrophobicity control the two-column
design was built on survives intact.

**XBridge is not a third point on a line. It is a test of whether the line exists.** It sits
only **0.14** from EVO on the `C(7.0)` axis, against 1.62 to Supelcosil — a *near-replicate*
at the low end. The logic:

> If two columns that are nearly identical in `C(7.0)` show `D_run` values differing by more
> than measurement error, then **`C(7.0)` is not the variable explaining `D_run`**, and the
> whole silanol-conditioning question is being asked of the wrong axis.

The two-column design cannot detect that failure at all: any two points define a slope, and a
confounded slope looks exactly like a real one. This is the cheapest available check on #39's
central premise, and it fails loudly rather than silently.

It is also **free of new risk**: XBridge C18 is in WSU-2019 (so its LSER constants are already
held, and tiers 2b/3 can run on it), and it is rated pH 1–12.

**Scope: tier-1 acids at the four shared pH points × 3 ACN compositions — 12 conditions.**

| option | conditions | campaign total | cost |
|---|---|---|---|
| full third arm (7 pH × 4 φ) | +28 | 72 | +64% |
| **shared pH × 3 ACN** | **+12** | **56** | **+27%** |
| shared pH × 30% ACN only | +4 | 48 | +9% |

A full third arm is not warranted: XBridge is a control on `D_run`, and `D_run` is estimated
from tier 1, which needs only the shared acidic-to-neutral range. The 30%-only option is
tempting at +9% but leaves no φ replication, so a single bad condition would be
indistinguishable from a real disagreement — which defeats the purpose of a control.

⚠ **XBridge C18's pH rating is unverified** (`model/columns.py` flags it). It only needs
2.5–6.5 for this scope, comfortably inside any plausible BEH rating, but confirm before use.

⚠ **XBridge C18 is fully porous** (N ≈ 11,447), matching Supelcosil rather than EVO. Pooling
is sized on that plate count already, so no change is needed — but it means **two of the three
arms share the limiting efficiency**, which makes the Supelcosil-sized pools the right choice
for all three.

### Parked: a polar-embedded study, deliberately not folded in

**Zorbax Bonus RP** (`C(7.0)` = **−1.10**, the lowest in the 819-column database) and
**Ascentis RP-Amide** (+0.08) would extend the silanol span to 2.85. They are **not** added,
because their `H` is 0.65 and 0.84 against this campaign's 1.01: they differ in retention
strength as much as in silanol activity, and **a `C`-versus-`D` relationship fitted across
columns differing in both is confounded beyond repair by any amount of instrument time.**

They are a good *separate* question — do polar-embedded phases suppress `D_run` beyond what
their `C(7.0)` predicts? — and that question deserves its own design with its own control.

### ✅ Availability and pH ratings confirmed (2026-08-16)

| Column | Vendor status | Rated pH | Sweep asks for | Verdict |
|---|---|---|---|---|
| **Kinetex EVO C18** | Phenomenex, current; 1.7 / 2.6 / 5 µm | **1–12** | 2.5–10.5 | ✅ inside |
| **Supelcosil LC-18** | Sigma-Aldrich / MilliporeSigma, current; 3 / 5 µm, several dimensions | **2–7.5** | 2.5–6.5 | ✅ inside |

**The pH cap is vindicated, not merely prudent.** Supelcosil LC-18 is rated to **7.5**, so the
superseded design's *shared* pH 8.0 point would have been **over the limit** for that arm.
Capping the high-silanol arm at 6.5 is what makes this pair usable at all, and 6.5 sits
comfortably inside the rating with margin for the `s_s` shift.

Provenance: Supelcosil pH range read from the Sigma-Aldrich product page by the project owner;
Kinetex EVO's pH 1–12 from Phenomenex's own product literature. Both are **vendor product-page
claims, not a read data sheet** — #13's data-sheet invariant still wants the printed
certificate filed with the campaign record.

**Suggested format pairing: 150 × 4.6 mm, 5 µm on both arms.** Both are catalogued in that
geometry, which keeps run times, injection volumes and back-pressure comparable.

⚠ **The two are not the same particle architecture** — Kinetex EVO is core-shell, Supelcosil
LC-18 fully porous — so their plate counts differ substantially. That is **acceptable here**:
`D` is recovered from a *retention ratio*, and efficiency affects peak width, not retention.
It does mean the two arms are not comparable on resolution or on tier 3's peak-shape
behaviour, and any such comparison must say so.

The USP SRM 870 parameters (`BD` = bonding density, `CTF`/`TFA`/`CFA` tailing-factor-like
activity measures) are **confirmed correct by the project owner**, with documentary
verification against *Pharmacopeial Forum* 31(2) 637–645 deferred. They are **no longer the
basis of this selection** — #39 showed they do not rank modern low-activity phases — but both
selected columns appear in that database, so they carry a third independent activity reading.

<details>
<summary><strong>Superseded selection (XBridge Shield RP18 / Betasil C18) — resolved 2026-08-16</strong></summary>

The original pair was chosen on the USP SRM 870 parameters (`CFA` 2.4 vs 7.6, `TFA` 1.1 vs
2.0). #39's join showed it **does not contrast on the ion-exchange axis at all**:

| Column | `H` | `C(2.8)` | `C(7.0)` | HSM `type` |
|---|---|---|---|---|
| XBridge Shield RP18 | 0.83 | −0.12 | −0.05 | EP |
| Betasil C18 | 1.05 | −0.03 | −0.04 | B |
| Betasil C18 *(duplicate HSM row)* | 1.063 | 0.095 | 0.099 | B |

Δ`C(7.0)` = 0.01. The `CFA` contrast is real in `CFA` but does not track ion exchange among
modern low-activity phases — #39 measured that correlation as leverage from a handful of old
type-A columns (`CFA`–`C28` Pearson 0.456 → 0.121 once six are dropped). That artefact landed
on this selection. A null result would have been **uninformative rather than negative**,
because silanol activity was never varied.

**Two further findings from the replacement analysis, both now folded in above:**

- **Not one WSU-2019 column is type A**, and its C18 members span `C(7.0)` −0.09 to 0.30. The
  only high-`C` WSU entries are fluorinated (`Kinetex F5` 2.66, `Discovery HS F5` 1.42,
  `Fluophase-RP` 1.41) — a different mechanism, and the phases #40 flagged for electrostatic
  behaviour. **The needed contrast does not exist inside WSU-2019**, which is why the selection
  is now asymmetric.
- The apparent "fixed baseline vs real contrast" trade **was a false dilemma**: tier 1 fits
  `k_neutral` from the data, so WSU membership was never required on both arms.

`Apex II C18` (`C(7.0)` 2.69) scored highest but is a Grace/Jones line and probably
discontinued; it beat `Supelcosil LC-18` by 1.7 points of population coverage, which did not
justify the sourcing risk.

⚠ `Betasil C18` has **two conflicting HSM rows** (above). Whatever columns are used, pin
parameters to a specific row and record which, not to a name.

</details>

### pH — deliberately asymmetric, because the columns are

| Arm | pH levels |
|---|---|
| Both columns | 2.5, 3.5, 5.0, 6.5 |
| Kinetex EVO C18 only | + 8.0, 9.5, 10.5 |

**The high-silanol arm stops at 6.5, and costs almost nothing for doing so.** Supelcosil LC-18
is type-A silica, the generation most vulnerable to dissolution at alkaline pH, and this
protocol's whole value depends on that column surviving the campaign. Tier 1 does not need
more: an acid at pKa 4.2 is **99.5% ionised at pH 6.5**, so the ionised plateau is already
reached. Going on to pH 8.0 would raise that to 99.98% — no useful gain against a real risk of
losing the column mid-campaign.

Kinetex EVO C18 is an organo-silica hybrid rated pH 1–12, so it carries the alkaline points
that tier 2b needs. **Confirm both ratings against the current manufacturer data sheets before
injecting anything** — #13's data-sheet invariant says such a field is evidence about the
column, and here it is an experimental safety constraint rather than a model input.

Record **both** `w_w`pH (buffer measured before mixing) and `s_s`pH (measured in the mixed
mobile phase) at every condition. #21 carried the missing `s_s` correction as a named error
term; this closes it.

### Mobile phase

Acetonitrile at **30, 40 and 50 % v/v** — chosen to land exactly on WSU's φ grid, so the
published system constants apply with **no interpolation** (#24's interpolation rules never
engage).

Plus **one methanol condition, 40 % v/v**, at all shared pH values on both columns. This is a
control, not an arm: cation exchange is **acetonitrile-specific and suppressed in methanol on
the same column** (`CONTEXT.md`, from Poole 1600), so `D`'s per-run part should visibly shrink
in methanol. It is a direct test of a claim the architecture already asserts, for the cost of
a solvent bottle.

### Compound set — three tiers of evidence, never pooled

The tiers exist because **`D` is not identifiable for strong bases on the Supelcosil arm.**
Fitting `D` means tracing retention from the neutral plateau to the ionised plateau. Acids at
pKa ≈ 4 reach both inside pH 2.5–8. Strong bases at pKa ≈ 9.5 need pH ≈ 11.5 for their
neutral plateau, which Supelcosil does not survive (and which its pH 6.5 cap forecloses anyway) — so on that arm `k_neutral` would have to
come from the LSER prediction, which is the quantity under test. That is circular, and the
tiers keep it visible.

**Tier 1 — acids (pKa 3–5). `D` measured on both columns, no model input.**
This tier carries the per-run silanol term, as the EVO-minus-Supelcosil difference per
compound.

**Candidate set, 2026-08-16 — a working proposal, not a frozen plate.** Small acids dropped in
favour of drug-sized ones, on the reasoning below.
Screened on McGowan `V`, which dominates `k_neutral` (the `v` coefficient is 2.419 at 30%
ACN, the largest term), and on pKa sitting clear of both grid edges:

| # | acid | `V` | pKa | `k` @30% ACN |
|---|---|---|---|---|
| 1 | Ibuprofen | 1.777 | 4.45 | **25.3 measured** |
| 2 | Naproxen | 1.782 | 4.18 | **14.3 measured** |
| 3 | Ketoprofen | 1.978 | 4.20 | **10.7 measured** |
| 4 | Flurbiprofen | 1.839 | 4.20 | `V`-screened |
| 5 | Fenoprofen | 1.880 | 4.50 | `V`-screened |
| 6 | Diclofenac | 2.025 | 4.15 | `V`-screened |
| 7 | Gemfibrozil | 2.118 | 4.70 | `V`-screened |
| 8 | Indomethacin | 2.530 | 4.50 | `V`-screened |

All clear `V` ≥ 1.78, the regime where the three measured acids landed at `k` = 10.7–25.3.

**Dropped, with reasons:**

- **Benzoic acid** (`V` 0.932) — `k` = 1.1; ionised peak sits on `t0`. **Retained in the
  campaign** as a thin-slice compound for #45's calibration overlap, but **not counted in
  tier 1's paired difference**.
- **Salicylic acid** (`V` 0.990) — fails on size *and* on pKa: at **2.97** it is only ~75%
  neutral at pH 2.5, so the neutral plateau is never cleanly reached inside this grid.
- **4-Nitrobenzoic acid** (`V` 1.106) — fails on size.
- **Warfarin** (`V` 2.308) — large enough, but pKa **5.00** sits at the grid edge: only 96.9%
  ionised at pH 6.5 against 99.5% for a pKa-4.2 acid. Usable if a ninth is wanted, not
  preferred.

**Chemical diversity is deliberate**, not incidental: four 2-arylpropionic acids (1–4) would
be near-replicates, so the set adds an arylacetic (diclofenac), an indole-acetic
(indomethacin) and a branched aliphatic (gemfibrozil). `D`'s per-compound part should be
estimated across acid classes, not across five close analogues.

⚠ **Five of the eight are `V`-screened, not `k`-verified** — we hold SoluteML descriptors only
for the thin-slice compounds. `V` is the dominant term but not the only one (ketoprofen has a
higher `V` than naproxen and a lower `k`), so the five should be confirmed once descriptors
exist. The screen is a strong filter, not a guarantee.

> **Sizing and selection, from the power analysis** (`experiments/sweep-column-decision/`):
>
> - **Eight acids is the right number, and more is waste.** The paired difference reaches
>   SD ≈ 0.004 at n = 8 and does not improve at n = 12, because it is floored by the two
>   per-column `t0` systematics, which do not average down over compounds.
> - **Prefer acids that are well retained in their neutral form.** `D`'s precision is set by
>   whether the *ionised* form is still retained: SD(`D`) is 0.070 at `k_neutral` = 2 against
>   0.003 at `k_neutral` = 40, a twenty-fold swing. A weakly retained acid contributes almost
>   nothing, because its ionised peak sits on top of `t0`.
> - **Run tier 1 at 30% ACN** — the most retentive of the three φ points — rather than
>   spreading its weight across 50%, where `k_neutral` falls and `D` degrades.

> ### ⚠ Checked against the actual compounds (2026-08-16) — one of them fails
>
> `experiments/sweep-column-decision/tier1_check.py` predicts `k_neutral` on the EVO arm
> from the slice's SoluteML descriptors. Of the four suggested acids we hold descriptors for:
>
> | acid | `k` @30% ACN | SD(`D`) | verdict |
> |---|---|---|---|
> | Ibuprofen | 25.3 | 0.005 | **good** |
> | Naproxen | 14.3 | 0.010 | usable |
> | Ketoprofen | 10.7 | 0.013 | usable |
> | **Benzoic acid** | **1.1** | **0.070** | **drop from tier 1** |
>
> **Benzoic acid cannot be rescued by moving φ.** Even at 10% ACN it reaches only `k` = 5.3,
> and 10% is outside this design's 30/40/50 grid. Its ionised form would elute essentially at
> `t0`, which is precisely the regime where a `t0` error swamps `D`.
>
> **Keep benzoic acid in the campaign, but not in tier 1's paired difference** — it is a
> thin-slice compound and therefore earns its place in the #45 calibration overlap. Those are
> two different jobs and only one of them needs it retained.
>
> **The pattern generalises.** The suggested list mixes small aromatic acids (benzoic,
> salicylic, 4-nitrobenzoic) with drug-sized ones (ibuprofen, naproxen, ketoprofen, diclofenac,
> warfarin). The small ones are systematically at risk on this axis, and salicylic and
> 4-nitrobenzoic acid should be expected to fail the same way. **Tier 1 should favour
> drug-sized acids** — which is also the project's domain.
>
> ⚠ **Four of the eight cannot be checked yet**: salicylic acid, diclofenac, 4-nitrobenzoic
> acid and warfarin have no descriptors in hand. They must clear the same bar before booking.
>
> **This is an optimisation, not a rescue.** Even with benzoic acid left in, the paired
> difference reaches a minimum detectable Δ`D` of 0.069 against a prior SD of 0.30 — still a
> 4× margin. Dropping it and matching the checkable median takes that to ~0.011. The design
> was already robust; this makes it sharper.

**Tier 2a — weak bases (pKa 6–8.5). Both plateaus in range on both columns.**
Lidocaine (thin slice), procaine, papaverine, nicotine.

**Tier 2b — strong bases (pKa 9–10). Model-free `D` on Kinetex EVO only.**
Propranolol, atenolol, imipramine (all thin slice), plus metoprolol, amitriptyline. On
Supelcosil their `D` is **conditional on the LSER neutral baseline and must be labelled as
such** — a different evidence class, reported separately, never averaged with tier 1.

**Tier 2c — amphoteric.** 4-aminobenzoic acid (thin slice), labetalol.

**Tier 3 — neutral controls, from the WSU S-2 anchor.** These do a different job: they test
the LSER+QSPR chain on these two columns with published constants in hand, which is #31's
drug-chemistry grading gap addressed on real retention. Choose 6 spanning the descriptor
range, e.g. pyridine (V 0.675), 3-nitrophenol (0.949), o-tolualdehyde (1.014),
1-chloronaphthalene (1.208), N,N-diethylaniline (1.380), hexanophenone (1.577). Their
retention must be pH-flat — if it is not, something other than ionisation is moving, and that
is itself a finding.

**`t0` marker (uracil or thiourea) in every run.** This is **not** a by-product: the power
analysis shows `D`'s precision is governed by `t0` accuracy, because the ionised plateau sits
close to `t0`, and the per-column `t0` systematic is what floors the paired difference. It
settles #25 as well.

### Data acquisition and archival — a requirement, not an IT preference

**Archive the raw vendor directories, not exported peak tables**, and acquire **full DAD
spectra**, not a single wavelength. Both are cheap at acquisition time and unrecoverable
afterwards.

**Why raw directories.** #45's post-mortem on RepoRT (`experiments/report-overlap/`) found the
corpus unusable for three reasons, and the load-bearing one is that **RepoRT records the
column but never the LC system**, so the vendor dwell and extra-column volumes in
[`sources/instrument-volumes/`](../sources/instrument-volumes/README.md) cannot be applied to
any of it. An Agilent `.D` directory does not have that hole: `acq.txt` records the **module
part numbers** — e.g. binary pump `G1312B`, sampler `G1329B`, column compartment `G1316B`,
DAD `G1315B` — alongside the full solvent timetable, flow and column temperature. That is the
system identity our dwell table is keyed on. Waters `.raw` and Shimadzu equivalents carry the
same class of metadata.

Export a peak table instead and every one of those fields is gone, and the dataset acquires
exactly the defect that made RepoRT's 421 datasets worthless to us.

⚠ **The raw directory is necessary, not sufficient.** It does **not** record dwell volume as a
value, nor `t0`, nor the column part number. Those stay explicit protocol obligations:

| Must be recorded per run | Where it comes from | Consequence if missed |
|---|---|---|
| `t0` | the uracil/thiourea marker required above | `D` is unfittable — precision is governed by `t0` |
| Column identity **and serial** | the logistics table, written into the sequence | arms cannot be told apart on re-analysis; the spare Supelcosil is indistinguishable from the first |
| LC system + module part numbers | the raw `.D`/`.raw` directory | dwell and extra-column volumes cannot be looked up — RepoRT's exact failure |
| Measured dwell volume | one gradient step-test per system, once per campaign | the lookup stays a vendor nominal rather than a measurement |
| Mobile phase pH, `s_s` scale, and the buffer | the run log | #20's provider computes `f_neutral` on `s_s`; a `w_w` number silently shifts every apparent pKa |

**Run the dwell step-test once per system at the start of the campaign** (isocratic hold with
a UV-absorbing tracer stepped into B, dwell read off the inflection). It takes minutes, it
converts the vendor table entry from a nominal into a measurement for this instrument, and it
is the difference between the sweep's data being reusable by anyone else and being another
RepoRT.

**Why full DAD spectra.** A single wavelength gives retention and nothing else. Full spectra
give each peak a UV signature that is **orthogonal to retention time**, which is what makes
cross-run peak tracking possible when the elution order itself is what changes across the pH
axis — and confirming peak *identity* as pH moves is a stated requirement of this design. See
the peak-tracking evaluation below.

### Scale

(4 shared pH × 2 columns) + 3 alkaline pH on Kinetex EVO = **11 (column, pH) combinations**,
× (3 ACN φ + 1 MeOH φ) = **44 conditions**, **plus 12 for the XBridge C18 control arm**
(4 shared pH × 3 ACN φ, tier-1 acids only) = **56 conditions**. With ~30 compounds injected as 4–5 resolved
mixtures plus `t0`, that is roughly **280 isocratic injections**, or ~560 with the required
duplicates. One QC condition repeated at the start and end of each column's campaign to bound
drift.

Larger than the superseded 48-condition design, but for a different reason than it looks:
capping the high-silanol arm at pH 6.5 *removed* a shared level, and the growth is entirely
the XBridge control arm.

## Campaign logistics

### Columns to order

| | Specification | Vendor | Part |
|---|---|---|---|
| Low-silanol arm | **Kinetex EVO C18, 5 µm, 150 × 4.6 mm** | Phenomenex | ⚠ **confirm part number with Phenomenex** — not verified here |
| High-silanol arm | **Supelcosil LC-18, 5 µm, 150 × 4.6 mm** | Sigma-Aldrich / MilliporeSigma | **58230-U** |
| Control arm | **XBridge C18, 5 µm, 150 × 4.6 mm** | Waters | ⚠ **confirm part number** — not verified here |

Matched geometry across all three keeps `t0`, run times, injection volumes and back-pressure
comparable. ⚠ Two of the three (Supelcosil, XBridge) are **fully porous** at N ≈ 11,447 while
EVO is **core-shell** at 19,553, so pooling is sized on the fully porous plate count and used
on every arm.

**Order a spare Supelcosil LC-18.** It is type-A silica running at pH 6.5 for several days
unattended. Its rating is 2–7.5, so 6.5 is inside it — but the failure mode of a type-A
silica near its ceiling is gradual dissolution that degrades the arm without announcing
itself, and losing it mid-campaign costs the whole high-silanol arm. A second column is
cheap insurance and also gives a column-to-column reproducibility check.

Guard cartridges for both arms are advisable given ~30 compounds injected repeatedly.

### Instrument time to reserve

`t0` ≈ 1.62 min for this geometry at 1.0 mL/min (ε ≈ 0.65). Run length is dominated by the
30% ACN block, where tier-1 acids are most retained:

| `k_neutral` | `t_R` | run + 20% re-equilibration |
|---|---|---|
| 3 | 6.5 min | 8 min |
| 10 | 17.8 min | 21 min |
| 25 | 42.1 min | 51 min |
| 40 | 66.4 min | 80 min |

| Block | conditions | inj/condition | run length | subtotal |
|---|---|---|---|---|
| 30% ACN | 11 | ~5 | ~45 min | ~41 h |
| 40% ACN | 11 | ~5 | ~20 min | ~18 h |
| 50% ACN | 11 | ~5 | ~10 min | ~9 h |
| 40% MeOH | 11 | ~5 | ~25 min | ~23 h |
| Equilibration between conditions | 44 | — | ~30 min | ~22 h |
| XBridge control arm | 12 | ~5 | ~25 min | ~25 h |
| | | | **single** | **≈ 138 h ≈ 6 days** |
| | | | **duplicated** | **≈ 254 h ≈ 10.5 days** |

**Reserve three weeks.** With the XBridge control arm the campaign is ~10.5 days of continuous
running before any QC repeats, buffer preparation, column equilibration on changeover, or
reruns — and there are now three columns to change over rather than two. Two weeks no longer
has margin; three does. Twelve working days is the floor.

The control arm is **the first thing to cut** if time is short: it costs ~25 h and tests
whether `C(7.0)` is the right axis, which is valuable but not what the campaign is *for*.
Cutting the 30% ACN block instead would be a false economy — see below.

⚠ **Do not shorten by dropping the 30% ACN block.** It is the most expensive block and the
one tier 1 depends on — `D`'s precision is governed by how retained the *ionised* form still
is, and 30% is where that is best.

### Per-compound φ selection — a refinement, not a design change

Tier 1 is specified to run at 30% ACN, but the tier-1 set now spans `V` = 1.78 to 2.53, and
`k_neutral` rises steeply with `V`. Indomethacin at `V` 2.53 will be far more retained than
ibuprofen at 1.777, and at 30% ACN could sit beyond `k` = 40, i.e. an 80-minute run.

**No design change is needed, because the sweep already collects every compound at all three
ACN compositions.** The refinement is in the analysis: **take each compound's `D` from the φ
where its own `k_neutral` lands in roughly 10–40** — high enough that the ionised form is
still retained and `t0` error does not swamp it, low enough that the run is practical. The
paired EVO-minus-Supelcosil difference must of course be taken **at the same φ for both
arms**, which is the only constraint this imposes.

If a compound exceeds `k` ≈ 40 even at 50% ACN, it is too retained for this campaign and
should be swapped out rather than accommodated.

## Analysis plan

Per compound, per (column, φ, modifier), fit

```
k_obs(pH) = k_neutral · [ f_neutral(pH) + (1 − f_neutral(pH)) · 10^(−D) ]
```

with `f_neutral` from the microspecies distribution on the `s_s` scale (#20's provider).
Free parameters `k_neutral`, `D`, and apparent pKa; 5–7 pH points identify all three **when
both plateaus are present**. Where they are not (tier 2b on Supelcosil), `k_neutral` is fixed to
the LSER prediction and `D` is reported as conditional.

Report, per #27's convention as extended by #31 — every figure carrying baseline, `n`,
interval, axis, **reference identity and reference coverage**:

- `D` per compound with interval, by tier.
- **`D`'s per-run part**: the paired EVO-minus-Supelcosil difference, tier 1 only.
  Minimum detectable difference is **0.0075** log units (95%, n = 8), against a prior of
  `D_run ~ N(0, 0.30)` — so any real per-run term will be resolved comfortably, and a null
  is a genuine null rather than a power failure. If it is
  large, the per-run term is real and #39's silanol conditioning is warranted; if it is
  small, `D` transfers cheaply between columns — very good news for the architecture. **Both
  answers are decision-changing**, which is the point of two columns.
- **Modifier check**: does the per-run part shrink in methanol, as the architecture claims?
- **Reversal rate** pH 2.5 → 8 on ordinary hardware, with `n` and interval, against #4's
  23.5% (71 compounds, 2 columns, 1 lab, on a positively-charged Kinetex PS C18).
- **Neutral-control residuals**: predicted minus observed `log k` through published WSU
  constants — the #41 criterion on real retention, on drug-relevant columns.

## Gate on the spend

Commissioning is **committed in principle** and released after #45 reports, as a **sanity
check, not a value test**: #45 cannot fit `D` (single-pH data cannot trace the sigmoid), so it
will not tell us whether the sweep is worth doing — the variance table already did. It tells
us whether something upstream is structurally broken, which we would rather learn before
booking instrument time.

⚠ **The gate is bounded.** #45 is itself blocked on a data acquisition, so gating on it risks
chaining one deferral to another. If extracting suitable RepoRT data for the slice compounds
does not yield a usable set within a bounded effort, **the gate opens anyway** and the sweep
proceeds on the variance argument alone. Otherwise "wait for #45" quietly becomes "never",
which is the failure #30 was written to avoid.

> ### ✅ The gate has opened (2026-08-16)
>
> The bounded RepoRT search ran and returned a **reportable negative** —
> `experiments/report-overlap/`. Compound overlap is good (205 of 421 datasets hold at least
> one slice compound; one holds 10 of 15) but the conditions are unusable for three
> independent reasons: **no isocratic data at all**, **no dwell-volume field anywhere in
> RepoRT's schema** so no gradient can be inverted to `k`, and **`column.t0` imputed from
> geometry rather than measured** (48 datasets carry the impossible `t0 = 0`).
>
> So #45 cannot report a calibration verdict from public data, and the escape clause above
> applies exactly as written. **The sweep proceeds on the variance argument alone.**
>
> This also *raises* the sweep's value rather than leaving it unchanged: #45's fallback data
> route was always this sweep, so the sweep is now the project's **only** route to any
> measured retention data. The one residual public lead — RepoRT dataset 0415, needing an
> instrument model from its source paper — could at best support an order-accuracy check,
> never the per-stratum `λ` #15 pre-registered.

Pursuing collaborators with existing unpublished data runs in parallel and is not an
alternative: a pH sweep on a plain C18 with drug-like ionisables is routine work that many
pharmaceutical laboratories hold unpublished, and it costs conversations rather than
instrument hours.

## What this experiment cannot answer

- **A transfer model for `D` across the 25-column axis.** Two columns give evidence on
  whether the per-run part is large; they do not parameterise it. That needs many more
  columns, or the HSM bridge (#26/#36/#39).
- **Anything about gradient or dwell physics** — isocratic points only, as with #21.
- **`D` for permanently charged species** (#39): they have no neutral microspecies at any pH,
  so the sigmoid has no upper plateau by construction. They stay outside this design.
- **Bias in the QSPR layer generally** — the neutral controls grade the chain on 6 compounds
  on 2 columns, which is a spot check, not #35's clean-split bias measurement.
