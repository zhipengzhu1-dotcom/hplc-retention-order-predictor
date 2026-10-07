# CONTEXT.md — ubiquitous language

Vocabulary pinned by resolved decisions. Chromatography, QSPR and Bayesian vocabularies
collide in this project and the collisions are load-bearing; when a term is ambiguous
across those worlds, the meaning recorded here wins. One entry per term, with the ticket
that pinned it.

## Terms

- **Resolution map** — the 2-D surface of separation quality over a *factor pair*, the
  signature artefact a chemist judges the system by. Factor pairs are pluggable
  (tG × pH, %B × pH, %B × T prototyped), not fixed.
  (#10)

- **Critical pair** — at a given operating point, the adjacent peak pair with the
  minimum resolution; **critical Rs** is that minimum. The critical pair's *identity*
  changes across the map, and those region borders are informative structure.
  (#10)

- **Pair territory** — the rendering of the resolution map coloured by critical-pair
  identity (depth of colour = that pair's Rs). (#10)

- **Threshold-anchored colour** — the requirement that any colour scale a chemist judges
  by is anchored to the acceptance thresholds (Rs 1.5, 2.0), never linear in the raw
  quantity. A linear ramp rendered a failing Rs 1.4 at ~60% saturation; that is a
  defect class, not a styling preference. (#10)

- **Predicted design space / established design space** — never the same object, never the
  same words. The **predicted** design space is our artefact: the region where the model
  estimates the declared criteria are met at stated confidence — **a prioritised hypothesis
  for verification**. The **established** design space is the ICH Q14 regulatory object, which
  a sponsor justifies with experimental evidence and a model cannot assert. The calibrated
  probability's regulatory value is that it says **where to verify and how much**, never that
  it replaces verification. Required wording: *"the model predicts, at stated confidence, that
  the declared criteria are met in this region, to be confirmed by verification runs"* —
  carrying the **declared criterion set**, the **calibration state** (`λ` and its #15 band),
  the **reference identity and coverage** (#31), and whether the surface is **pre- or
  post-reconciliation**. **Forbidden:** describing a predicted region as verified at any `P`;
  and any cold design-space claim for a mixture containing ionisables, since #34/#15 grant
  that capability only post-reconciliation after 3 runs at ≥ 2 pH units.
  (#19,
  #34,
  #31)

- **Feasibility and confidence are separate surfaces, and `P(all Rs ≥ target)` fuses them.**
  Measured: on a 50 × 50 (φ × pH) grid with 15 compounds, **max P(all) = 0.015 over 2,500
  cells** — the surface is uniformly zero, so it has no gradient to navigate and no optimum to
  find. The cause is the **conjunction over `n − 1` adjacent pairs**: at an expected resolved
  fraction of 0.845, requiring all 14 gives 0.084, and P(all) falls 0.92 → 0.55 → 0.08 → 0.007
  as compounds go 4 → 8 → 15 → 20. ⚠ **Cutting model uncertainty 4× makes P(all) *worse***
  (0.084 → 0.040 at n = 15) — as the ensemble converges on the true compound set, an
  unavoidable coelution becomes certain rather than merely likely. So P(all) is dominated by
  **feasibility** (chemistry; a better model cannot move it) and not by **confidence** (the
  term this project exists to handle honestly). The decision layer therefore shows the two
  separably and **always names which term is binding and which pair binds it**. Fourth
  sighting of the same distinction — cf. #40, #31, #42.
  (#19,
  #10)

- **The criterion set is declared by the user; the default is all adjacent pairs.** The single
  largest lever on whether a design space exists, and it is chemistry knowledge we cannot
  infer. Best-cell P at 15 compounds: **0.079 for all 14 pairs, 0.19–0.29 for 8, 0.53–0.70 for
  3, 0.81–0.94 for 1.** All-pairs is the default because it is the strictest reading and a user
  who has not thought about it must not silently receive a weaker claim; narrowing is an
  explicit act recorded in the method card and travelling with every claim. This is also the
  correct reading of ICH Q14 — a design space is defined relative to *declared* attributes.
  ⚠ **Declaring a smaller set does not manufacture a design space**: three hardest pairs still
  reach only 0.53.
  (#19)

- **Recommend the most interior point, not the highest-probability point — and when the design
  space is empty, recommend a programme, not a method.** Margin is what a chemist needs: the
  drift the method tolerates before failing its criteria, in units of real operational
  variability (buffer pH stability, #11's ±1–2 °C, flow and gradient reproducibility). So the
  recommendation is the point furthest from the design-space boundary, and the Pareto front is
  **run time against interior margin**, not run time against `P` — two points with equal `P`
  can have very different margins. #10 found this empirically: the probabilistic layer moves
  operating points off narrow ridges onto robust plateaus. **An empty design space is a routing
  decision, not a null output**: if *confidence* binds, emit #18's scouting programme targeting
  the binding pair; if *feasibility* binds, state that the separation is unreachable in the
  explored factor space and name what must change (column, gradient axis, modifier). Never
  return the least-bad cell as if it were a recommendation.
  (#19,
  #18,
  #10)

- **Design-space computation is cheap, because regeneration is amortised.** A scenario is a
  **parameter draw**, so one ensemble evaluates anywhere on the factor grid; regeneration (#18)
  happens once per *reconciliation event*, never per grid cell. Benchmarked in plain NumPy with
  15 compounds including the per-cell sort and adjacent-Rs computation: **50 × 50 grid at 1,000
  scenarios ≈ 1.0 s**, 5,000 ≈ 4.1 s, 20,000 ≈ 20 s; 100 × 100 at 1,000 ≈ 3.1 s. Interactive
  without optimisation. The body's sampling-cost worry is answered; the binding constraint on
  this layer was never compute.
  (#19,
  #18)

- **Method card** — the method context attached to any pinned point of a resolution
  map: gradient table or isocratic hold, pH, temperature, column, flow, t0/dwell,
  plate count. Part of the artefact's contract, not decoration.
  (#10)

- **Chromatogram (as output)** — a *rendering* of the ranked retention distribution,
  never the prediction layer's product. (Charting decision on
  #1)

- **Scenario ensemble** — the prediction layer's output interface, and the concrete
  meaning of "a ranked retention distribution with calibrated confidence". Each
  scenario is a *complete alternative chromatogram*: one draw of the per-run
  parameters plus one correlated draw across every compound. Stored as
  `retention[S, n]` and `width[S, n]`. Per-compound marginals are a **derived
  display view** and a lossy projection — they destroy the correlation that elution
  order depends on. (#14)

- **Per-compound error / per-run error** — the two error classes, kept structurally
  separate in the sampler, classified by *scope of sharing* rather than by
  random-vs-systematic. Per-compound (QSPR descriptors, pKa, LSER residual
  lack-of-fit) is drawn independently per compound and **scrambles elution order**.
  Per-run (LSER system constants, dwell volume, `D_m`, band compression `G`) is one
  draw shared by every compound and **largely cancels in order** while still setting
  absolute retention. Plate count and extra-column dispersion are per-run but stay in
  the sampler, because they set peak width and therefore resolution.
  (#14)

- **Conformal stratum** — the class a compound falls into for Mondrian
  (class-conditional) conformal calibration, at minimum by ionisation regime, likely
  also by column type and φ range. Stratification is required because the error is
  regime-dependent: a neutral far from any pKa carries roughly 0.09 + 0.04 in
  `log k`, while a compound inside the ±1.5-unit pKa danger band picks up ≈ 0.36 from
  pKa error alone. One global inflation factor cannot serve both.
  (#14)

- **Reweighting, not refitting** — how measurement enters. A real injection does not
  trigger a refit; it filters or importance-weights the scenarios that remain
  consistent with what was observed. Peak-assignment correction is a filter;
  reconciliation is a weighting. One mechanism, applied uniformly at ensemble level.
  ⚠ **Amended by #16/#17:
  this rule is incomplete.** Reweighting can only narrow to what was already sampled, and
  two routes reach outside that support — a human constraint with zero support, and data that
  outruns the prior. Full rule: **reweight while ESS is adequate; otherwise
  regenerate conditioned on everything known so far, and log what triggered it.**
  Regeneration is a *named, logged, triggered event*, never a silent refit.
  ⚠ **Corrected by #18:
  the "3–4 runs" estimate was optimistic by ~4×, and the unit was wrong.** The trigger is the
  **first** run, and the driver is **how many compounds are observed**, not how many runs are
  performed. Measured on the thin-slice ensemble (1,000 scenarios), median ESS after one run:
  27 at 1 compound / sd 0.03, 2.9 at 2 compounds, **1.0 at 5+ compounds at any precision** —
  and a real run measures the whole mixture at once. The failure is **dimensional**, not a
  shortage of scenarios: scenarios differ in every compound simultaneously, so conditioning
  means matching an `n`-dimensional vector, which importance sampling cannot do at any
  achievable ensemble size. **Regeneration is therefore the declared default path, not the
  exception.** Reweighting survives only for sparse observation (≈1 compound, coarse
  precision). ⚠ Any design or metric scored by "which scenarios survive the filter" is scoring
  a degenerate quantity — conditioning appears to raise the confident-pair fraction to ~1.0,
  which is the ensemble recognising its own labelled draw, not learning.
  (#14,
  #16,
  #17,
  #18)

- **Ionisation retention drop `D`** — the parameter that carries all ionic retention:
  `k_ion = k_neutral · 10^(−D)`, with `D ≈ 1–2` in RPLC. It is **empirical, not
  physical** — LSER is licensed for the neutral microspecies only, so `D` absorbs the
  mechanisms the solvation parameter model never parameterised (ion exchange with
  silanols, ion-pairing, electrostatics). Never zero: `D = ∞` is the "ionised species
  are unretained" claim, a ~2 log-unit error against a ~0.09 per-compound budget. It is
  a **prior the user's own scouting runs refine**, decomposing into a per-compound part
  (which ionised species) and a per-run part (column silanol activity, modifier).
  (#34)

- **logD is not in the ionisation interface, and the reason is architectural rather than
  aesthetic.** The provider returns **per-microspecies fractions with their charges and their
  own SMILES**, never a scalar. logD is rejected on three grounds, the third being the one
  that binds: it is **octanol/water**, the wrong phase system, so substituting it for the LSER
  computation is a category error; it is **strictly redundant**, derivable from fractions ×
  per-species logP; and — decisively — **silanol cation exchange acts on the cationic
  microspecies only**, so its magnitude is proportional to the very fractions the provider
  computes. A scalar logD (or a scalar α) would make that term impossible to attach later
  without reopening the ionisation interface. **Expose logD, if at all, as a convenience
  method on the result object, never as a required provider capability.** Note `Crippen logP`
  is unrelated to this: it is #27's *cold baseline*, a thing to beat, not a model input.
  (#7,
  #34,
  #39)

- **`t0` is fitted from scouting runs, never assumed constant across pH.** Hold-up time
  falls as residual silanols ionise — a reported 0.83 → 0.65 min (−22%) to pH 11 — and
  because `k = (t_R − t0)/t0` is most sensitive to `t0` exactly where `k` is smallest, the
  error lands almost entirely on the **ionised plateau**, which is the plateau that defines
  `D`. Measured over the sweep's own grid (pH 2.5–6.5), holding `t0` fixed biases `D` by
  **+0.42** at silanol pKa 6.0 and **+0.69** at pKa 4.5 — larger than the whole prior SD on
  `D_run` (0.30), 14× the LSER lack-of-fit floor, and 56× the difference #30 is commissioned
  to resolve. **A bias that size would be mistaken for the silanol effect the sweep is
  measuring.** So: a **`t0` marker (uracil or thiourea) is required in every scouting run**,
  `t0(pH)` is fitted per column and mobile phase, and before scouting exists the
  pH-dependence is carried as a **named per-run error term** with the prediction flagged as
  extrapolation across pH. Modelling `t0(pH)` instead was rejected: it needs a per-column
  silanol-ionisation parameter, which is the same axis #39 found unresolvable among C18
  phases on all three available instruments.
  (#25,
  #30)

- **`t0` is per-run in scope but does not cancel in elution order, and that is why it is
  dangerous.** At a single condition every compound shares one `t0`, and `k` is monotonic in
  `t_R`, so a `t0` error cannot reorder that chromatogram — it behaves like the other per-run
  terms. The damage is on the **pH axis**: a pH-dependent `t0` error propagates into `D`, and
  a biased `D` moves ionised species relative to neutral ones at every pH we *predict* rather
  than measure. So `t0` sits with dwell volume in the per-run class by scope of sharing,
  **but its error must not be assumed order-cancelling** the way per-run terms normally are.
  It rides with the **method card**, not the solute.
  (#25,
  #14)

- **In envelope / degraded / refuse** — the three scope tiers for a compound.
  *In envelope*: a neutral microspecies is appreciably populated somewhere in the
  accessible pH range. *Degraded*: a neutral form exists but is negligibly populated, so
  retention runs almost entirely through the empirical `D` term — predict, but flag as
  extrapolation. *Refuse*: **no neutral microspecies at any pH** (quaternary ammonium,
  sulfonates), so `k_neutral` has no referent. Refusal is a **structural** test answered
  from the provider's microspecies output, never a confidence threshold — because a wide
  distribution still looks like knowledge.
  (#34,
  #33)

- **Charge sign is not a retention class** — the refuse tier applies to **both** charge
  signs, and the anionic/cationic asymmetry #39 was opened to adjudicate **does not exist**.
  Not because both elute early, but because neither sign licenses a retention claim, for
  two *different* reasons. **Cationic**: retention runs through silanol ion exchange, which
  the solvation parameter model does not parameterise (#36) and which no column in our
  anchor set exhibits (zero type-A columns in WSU-2019). **Anionic**: charge does not
  determine retention at all — a permanent anion's retention is governed by the hydrophobic
  surface area of the rest of the molecule, and spans the chromatogram (methanesulfonate at
  `t0`; dodecylbenzenesulfonate strongly retained on C18). A near-`t0` prior collapses
  `k_neutral · 10^(−D)` into a claim about `D` alone, implicitly asserting `k_neutral` is
  both small and *uniform* across permanent anions. That is the `D = ∞` error of #34 in a
  new costume, and the second half is simply false.
  (#39)

- **The honest prior degenerates to refusal** — the general reason the refuse tier is not
  merely conservative but *optimal* here. A prior over permanent-anion retention that
  respected the spread above would span the whole run, which for ranking purposes carries
  exactly the information of "unplaceable" while *looking* like a prediction — #34's "a wide
  distribution still looks like knowledge", reached from the opposite direction. Refusal
  delivers the same information with honest ergonomics, and routes the compound into the
  **expected but unplaceable** pool where it correctly competes for unexplained peaks.
  **The burden of proof sits on the prior, not on the refusal**, and nothing in the stack
  can discharge it: #39 established there is no descriptor referent, no fitted coefficient,
  and no anchor column in the regime.
  (#39,
  #34)

- **No instrument resolves silanol activity among C18 phases** — measured three ways, and
  they converge. USP SRM-870 activity parameters correlate with the HSM `C` terms only via
  a handful of old type-A columns (`CFA`–`C(2.8)` Pearson 0.456 → **0.121** with six
  dropped); WSU-2019 Table S-1's weak-vs-unqualified electrostatic flags separate on `C(7.0)`
  across all phases (mean 0.176 vs 0.435) but **not among C18 alone** (0.072 vs 0.050,
  n = 4 vs 5 — the weak group is if anything *higher*), the apparent agreement being carried
  entirely by the three fluorinated phases. All three
  instruments distinguish extreme phases and none ranks conventional C18 columns. So a
  **silanol-activity-conditioned prior is refused as unconstructible**, not merely unproven.
  ⚠ Note HSM's `C` is *cation* exchange: there is **no anion-exclusion parameter anywhere in
  the stack**, so the anionic case could not be column-conditioned even in principle.
  (#39,
  #40,
  #36)

- **Two-regime claim** — what the spec may assert. Neutrals are graded **cold**:
  pairwise elution-order accuracy and calibration as charting set. Ionisables are graded
  **post-reconciliation** — after `n` scouting runs at distinct pH. The ionisable layer
  is **validated in mechanism only**; the public data to falsify it does not exist.
  Thresholds for both regimes are owed by
  #15 — without
  numbers the claim is unfalsifiable by construction.
  (#34)

- **Cold baseline / condition-transfer baseline** — the two naive floors for elution-order
  accuracy, one per regime. **Cold**: order by Crippen `logP`, no measurement at all; scores
  *a priori* claims. **Condition-transfer**: the true order is known at one condition and reused
  everywhere; scores *post-reconciliation* claims. ⚠ They are not interchangeable and the error
  is directional — the transfer floor is **high** (76.5% pH, ~92–94% modifier, ~99% temperature),
  the cold floor is low. The often-quoted "92.3% with no model at all" is **misclassified**: it
  presupposes an oracle observation at one condition, so it is a *one-free-measurement* floor.
  (#27)

- **Conditional accuracy** — order accuracy restricted to the pairs that actually reorder between
  two conditions. The measure of *contribution*: the condition-transfer baseline scores **exactly
  0%** there by construction, random scores 50%. Always reported **with** aggregate pairwise
  accuracy, which is the guard-rail against a model that wins on reordering pairs by scrambling
  stable ones — neither is honest alone. Both carry `n` and an interval, because the reordering
  subset is smallest exactly where it matters most.
  (#27)

- **Discrimination gate / power gate** — the two contrasts that must *both* pass for a success
  claim, because #4 established no public contrast is both well-powered and strongly
  discriminating. **pH 3→8** is the discrimination gate (23.5% reversal, but 71 compounds, one
  lab, on a positively-charged-surface phase). **ACN→MeOH** is the power gate (5.6–8.2% reversal,
  overlap up to 406 compounds). Temperature and C18-vs-C18 are reported and explicitly
  **non-gating** — temperature's ~1% over ~95,000 pairs is a statement about the physics, not the
  sample size. (#27,
  #4)

- **Expected but unplaceable** — the pool holding compounds the system cannot place at all:
  `PROVIDER_DECLINED` (the `pKa` provider threw — ~1 in 15 compounds) and `OUT_OF_ENVELOPE`
  (no neutral microspecies at any pH). They **keep their column and index** in the ensemble
  with `NaN` retention, so `n` never shrinks and the compound list stays intact, but they
  carry **no fabricated distribution**. They are candidate explanations for *unexplained
  peaks* only — the same pool unknown impurities draw from — and never compete on retention
  against a compound that has a real distribution. **A refused compound and an unknown
  impurity are the same object to the tracker.**
  (#33,
  #16)

- **Per-scenario assignment** — how peak tracking works. Each scenario is *rendered* to a
  predicted peak list (compounds merging where that scenario's own widths say they
  co-elute), matched to the observed peaks allowing unmatched on both sides, and scored;
  results are aggregated **by counting across scenarios**. **Co-elution is a property of the
  rendering, not a case the matcher handles** — if A and B merge in 40% of scenarios, that is
  the co-elution probability. Cross-run linking falls out, since each scenario predicts every
  run. This is why the tracker is **not** an assignment problem: a bijection (Hungarian)
  structurally cannot express many-compounds-to-one-peak.
  (#16)

- **Effective sample size (ESS)** — a **required reported output**, not a diagnostic. It sizes
  `S` (which #14 left undecided at 1000), it is the **disagreement detector** — collapse means
  no scenario explains what was observed — and it is the **trigger for regeneration**. Matching
  on absolute retention is correct *within* a scenario, because per-run error is a scenario
  parameter rather than noise; the price is degeneracy when the common-mode prior is wide, and
  ESS is how that price is seen. (#16,
  #17)

- **Falsification event** — a user constraint with **zero support** in the ensemble, or ESS
  collapse with no fitting mode. Recorded, never silently discarded. User constraints are
  **hard** (they pin the assignment); self-contradictory constraint *sets* are rejected at
  entry before the model is consulted, because that is user error rather than model failure.
  A chemist contradicting the ensemble usually holds evidence the model lacks (MS, DAD,
  spiking), so **the prior is the thing presumed broken**.
  (#16,
  #15)

- **Joint `(scenario, assignment)` posterior** — what reconciliation carries. Because the
  assignment is evaluated *inside* each scenario, "the model is wrong" and "we tracked the
  peaks wrong" are two **modes of one posterior**, and are **never adjudicated** — any rule
  picking one invents information. Bimodal-with-decent-fit = genuine identity ambiguity;
  ESS-collapsed-everywhere = the model is falsified here. Trusting the tracking lets a
  tracking error silently corrupt the model; trusting the model makes it unfalsifiable.
  (#17)

- **Column vector** — the LSER system constants `c, e, s, a, b, v` for one column at one
  modifier and one φ. It is only meaningful **against the solute descriptor lineage it was
  regressed on** (for WSU, the descriptor database in Poole 1600 Table S-2), so a column
  vector and a solute vector from different lineages must never be multiplied. This is what
  forces the Abraham/WSU lineage on the solute side and rules out a learned latent solute
  vector. (#36,
  #23)

- **Topology-class substitution** — how a column with no published constants gets a column
  vector: borrow from the WSU column sharing its stationary-phase topology class (C18
  endcapped, polar-embedded, mixed-mode, sterically protected, C8, phenylalkyl, PFP,
  biphenyl), broken by H. **Always a degraded tier**, never silent, and **refused** when no
  WSU column shares the topology. Budgeted cost: SD ≈ 0.055–0.145 in `log k`, but only
  ≈ 0.007–0.018 on a close pair's separation — affordable *only* because elution order is
  the primary metric, and never to be quoted as affordable on `log k`.
  (#36)

- **`Fs` is for screening, never for substitution.** `Fs` shortlists *which different
  columns to physically try*. It may **not** select a column to borrow constants from:
  decomposed over the 819-column PQRI database with the published weights, 80.5% of `Fs²`
  is the C term (no LSER analogue) and 12.8% is B (ra² 0.051 to `b`), while H — the only
  parameter with real LSER signal (ra² 0.723 to `v`) — contributes 0.4%. **`Fs`-nearest is
  not LSER-nearest.** (#36,
  #8)

- **`Fs` is now calibrated, and it is quoted on elution order — never on `log k`.** Measured
  over 1,200 column pairs from the 21 columns carrying **both** HSM parameters and WSU LSER
  system constants (`experiments/fs-calibration/`), where `Fs` comes from one dataset and the
  substitution error from the other, so neither knows about the other:

  | `Fs` band | order flips (median) | close-pair Δseparation | RMS Δ`log k` |
  |---|---|---|---|
  | 0–5 | **0.9%** | **0.016** | 0.109 |
  | 5–10 | 1.5% | 0.025 | 0.120 |
  | 10–20 | 1.9% | 0.034 | 0.151 |
  | 40+ | 4.7% | 0.082 | 0.266 |

  **`Fs` ≤ 5 costs 0.016 on a close pair's separation — inside the 0.007–0.018 budgeted for
  topology-class substitution**, so it earns the same standing: a **degraded tier, never
  silent**, and quoted on order. The `log k` column is carried only to be refused as a
  headline: `CONTEXT.md`'s own rule is that a `log k` figure is never quoted as affordable,
  and the converse binds too — it must not be quoted to disqualify a substitution either.
  ⚠ These are *predicted*-order flips between two LSER models over 94 neutral anchors and are
  **not comparable** to #27's measured 4.6–6.2% C18-vs-C18 reversal rate.
  (#26,
  #36)

- **Solute-vector correction vs form correction** — the two halves of the residual
  corrector, kept separate because they *transfer differently*. The **solute-vector
  correction** is 4-dimensional and stays inside the bilinear form, so it transfers to every
  column automatically; it fixes descriptor error. The **form correction** is a free
  Δ`log k` that fixes LSS lack-of-fit (the largest term in #2's budget), applies **only** on
  columns with retention data, and is **zero elsewhere by construction**. Merging them would
  put the interface error at the layer boundary.
  (#36)

- **Validity register** — WSU-2019 Table S-1 as consumed by this project
  (`sources/wsu-lser/wsu2019-column-caveats.csv`): one register, **three semantics, never
  collapsed** — *physical invalidity* (dewetting → hard mask), *empirical absence*
  (steric exclusion → applicability-domain block), *mechanistic evidence* (electrostatic
  activity → visible caveat + evidence into `D` and #39). Each is machine-actionable and
  the distinction is preserved end-to-end, because a register consumed one way gets
  half-applied. (#40)

- **Dewetting mask** — structural refusal over **(physical column, modifier, φ)**
  condition space where the register records incomplete wetting. Refusal in #33's sense,
  one axis over: no budget is measurable (the published fits there *are* the distorted
  data) and wetting hysteresis means the phase state is not determined by our inputs — no
  referent, not low confidence. The mask lives in the **prediction layer's domain**: the
  resolution map renders it as a visually distinct **no-claim region** (a masked cell
  never carries a quality colour), the optimiser treats it as **hard-infeasible upstream
  of ranking**, and interpolation (#24) draws support only from unmasked fits and never
  crosses a floor. Measurement may still enter via reconciliation.
  (#40)

- **Applicability-domain block** — the treatment of empirical absence: the calibration set
  explicitly excluded the compound class, so the model has **no demonstrated predictive
  validity** — stronger than low-confidence extrapolation, weaker than invalidity.
  Default: `NaN` retention with identity preserved, status `OUTSIDE_APPLICABILITY_DOMAIN`
  (distinct from `PROVIDER_DECLINED`/`OUT_OF_ENVELOPE`; `n` never shrinks). **Per-request
  opt-in** (never per user/compound/column/project) ships a point value with
  `calibration = NONE`, `claim_eligible = FALSE`, `metric_eligible = FALSE`.
  `UNCALIBRATED` is a *calibration state*; `OPT_IN` is an *authorization/provenance
  state*; stored separately. **Invariant: an opt-in prediction can produce a value, but
  never confidence.** **Firewall: reconciliation hypothesis ≠ calibration evidence** — a
  later injection may test an opt-in prediction and record the discrepancy, but never
  retroactively calibrates it, and its residual is never an evaluation datum; #15/#27
  metrics count only in-domain predictions. "Bulky" awaits its operational definition
  (#43); a zero-support
  stratum cannot be conformally calibrated, so until then block + opt-in is the whole
  mechanism. (#40)

- **`CLASS_INFERRED_DEWETTING_FLOOR`** — the wetting boundary a *substituted* column
  vector carries: `basis = MAX_OBSERVED_WITHIN_TOPOLOGY_CLASS`, keyed **per (topology
  class, modifier)** — deliberately a different key from the per-physical-column mask,
  because topology substitution transfers **model parameters, never physical properties**;
  the donor's wetting guarantee is not borrowable, and the class maximum is the only
  empirically defensible bound (within-class floors span 10–40%). It is an inferred
  applicability boundary, never a measured floor of the target. Status chains are
  preserved (`OUTSIDE_APPLICABILITY_DOMAIN → TOPOLOGY_SUBSTITUTION →
  CLASS_INFERRED_DEWETTING_FLOOR`), never collapsed to a generic "degraded". Escape
  hatch: the per-request opt-in.
  (#40,
  #36)

- **Descriptor space and retention space are decoupled at the ranking level** — a model
  can win every descriptor column and lose retention (Ulrich Table 2's CHI/KRI
  inversions, on this project's exact chain; corroborated by the thin slice's
  attribution). Therefore: **the QSPR layer is graded in retention space, where the
  anchor lives** — predicted descriptors through fixed published constants vs
  log k(measured S-2 descriptors), aggregated over the full unmasked (column, modifier,
  φ) grid (upgrade form: residual SD of true re-fits, once raw k is acquired).
  **Descriptor rmse explains retention-space results; it never replaces them** — a
  diagnostic barred from every success claim, acceptance criterion and
  comparison-that-decides, #27-style. **ACD Absolv is in the QSPR-layer baseline set,
  out loud** — best published retention-space number (CHI 3.19 ≈ 0.072 log k via the
  measured slope, median 2.26, n = 1,587), cited never rerun, conflict-of-interest
  caveat named; QSPR-layer only, the stack's baselines remain #27's floors.
  (#41)

- **Per-compound confidence is a structural axis: who is confident is fixed before the
  method is chosen; how confident depends on the method.** Ensemble spread stratifies
  retention-space error by **2.3×** (low third 0.067, high third 0.154 log₁₀ k, n =
  32/31/30 over 332 unmasked fits) — the two thirds straddle the ~0.09 budget line, so a
  single pooled inflation misdescribes both, the #14 argument one layer down. The
  compound ranking is method-independent (Spearman ρ median 0.942 across method pairs),
  the magnitude is not (5.6×). Therefore **continuous: the full per-compound 4×4
  descriptor covariance goes into the sampler** — never a scalar inflating a shared
  covariance shape, which would recreate #37's independence error one layer higher;
  **discrete: a method-independent structural summary labels a coarse calibration
  stratum**, so a compound never changes conformal class because the column changed.
  Stratum count and boundaries belong to #15 — tertiles are Ulrich's reporting choice,
  not an architectural mandate. Descriptor-space gains do not carry: 4.5× on `A` becomes
  2.3× in retention space (#41 again — never quote the descriptor figure as a retention
  claim).
  (#42)

- **Regime validity and quantitative uncertainty are different kinds, and low confidence
  earns no status.** *Mask*, *applicability-domain block* and *degraded* say the value
  would not mean what values normally mean — a change of kind, which no width can
  communicate, hence #40's three separate labels. Wide-but-inside-regime is not a change
  of kind: it is the model reporting honestly, already fully carried as width and as
  P(A before B) → 0.5. So the unconfident majority is **predicted, wide** — no tier, no
  refusal (refusal stays structural, #33), and **no per-compound threshold**: thresholds
  attach to the *claim* being asserted, which already integrates every compound's width,
  never to the compound. Counting it twice would have the user discount a prediction that
  was already discounted. **Per-compound uncertainty is nonetheless a required output
  field** beside ESS — a reporting obligation, not a status — because a ranked list
  otherwise hides *which* compound drives an ambiguous ordering, and that is the
  actionable fact.
  (#42,
  #33,
  #40)

- **Doubt has two sources with opposite remedies, and they are never merged.**
  *Predictor uncertainty* — the ensemble is unsure — is spendable: more data, better
  architecture, or a measured descriptor for that compound all fix it. *Structural
  deficiency* — the Abraham representation cannot encode what the molecule does
  (intramolecular H-bonding, site-specific/sterically hindered donation, the several
  interactions `S` collapses into one; independently, S-1's own bulky-compound and
  H-bond-acid exclusions) — is **not** spendable: a better descriptor predictor cannot
  help there. Merging them lets the second class inherit the first's remedy and burns
  QSPR effort where it cannot pay. **The structural flag can never be derived from the
  ensemble** — its defining property is that the 25 networks agree, having learned the
  same inadequate representation from the same corpus, so any ensemble statistic is blind
  to it by construction. The register is declared, sourced externally, and currently
  **provisional** — no operational substructure queries yet.
  (#42)

- **The structural-deficiency register is named, not automated — and that is a decision, not a
  gap.** The classes are real and independently named *a priori* (intramolecular H-bonding,
  site-specific/hindered donation, dipolarity collapsed into `S`, bulky compounds, some
  hydrogen-bond acids). But tested against the only instrument available, **only dipolarity
  mixtures show the signature** — residual outrunning ensemble spread by 1.9× (`S`) and 2.1×
  (`A`) — on **n = 8**; two classes are untestable (the WSU-94 holds 2 intramolecular
  H-bonders and 1 hindered donor); and the naive operationalisation **flags 56% of
  compounds**, which is a blanket caveat rather than a discriminator. So: **no automatic
  flag ships**, the limit is documented in `spec/representation-limits.md`, and the
  uncertainty presented anywhere must say it covers *predictor uncertainty, not
  representational adequacy*. ⚠ **The contamination is asymmetric**: the WSU-94 sit inside
  SoluteDB, so residuals are too small — a positive result survived a headwind, a null did
  not. **This evidence may support a class; it may never retire one.** Revisit when an
  un-memorised test bed exists (Poole 2020, #35).
  (#46,
  #42)

- **"Layer-by-layer validation" is three things, and they have three different fates.**
  *Grading a layer in isolation* — **survives**, in retention space only: holding the LSER
  constants **fixed** is what isolates the layer, since constants that do not vary
  contribute no variance, so the whole score moves with the descriptors (#41's criterion is
  layer-isolated grading merely *denominated* in retention units — not the joint-validation
  fallback it resembles). *Attributing error between layers against measurement* —
  **surrendered** by decision (#14), so the interval can absorb unmodelled inadequacy
  honestly. *Attributing variance between layers* — **survives and is the strongest tool we
  have**: it needs no ground truth, is a by-product of the sampler, and already moved the
  effort budget (#21's three regimes). ⚠ The charter wanted layer separation for
  *debugging*; what survives delivers *effort allocation*. **Variance attribution is blind
  to bias by construction** — a layer can own almost no variance and still be systematically
  off, which is exactly `E`'s situation. #35 is the only instrument that touches bias.
  (#31,
  #41,
  #14,
  #21)

- **The anchor's applicability domain is a separate object from the model's, and it binds
  the grade, never the prediction.** SoluteML's training domain (8,366 solutes, drug-like
  chemistry included) is far wider than the 94-compound WSU S-2 anchor we grade against. So
  predicting descriptors for a drug is **not** an extrapolation; only *validating* it is.
  Therefore predictions proceed unblocked and the **retention-space grade carries a stated
  coverage boundary** beyond which it is not evidence — the #40 separation of "we will
  compute this" from "we can vouch for this", one layer up. Measured: 7 of the 9 thin-slice
  drugs sit at or above the anchor's 95th percentile on at least one descriptor, and 3 are
  **outside its observed range** (atenolol B 2.11 and ketoprofen S 2.28 against anchor maxima
  1.388 and 2.214); `B` is the systematic offender, drug bases carrying basicity the anchor
  does not contain. **For drug-like chemistry the QSPR layer is currently ungraded** — a
  reference gap, not a method gap; the realistic route is grading *through* retention on
  real data (#45), not against descriptor truth.
  (#31,
  #5,
  #35)

- **A diagnostic is maintained only if it needs no ground truth we cannot vouch for.**
  *Variance attribution*: **standing, required output** — no ground truth, nearly free, and
  it has already changed a decision. *Descriptor-space error decomposition against the
  corpus*: **retired as a maintained artefact** — computed ad hoc when debugging a specific
  compound, never tracked or reported. Its reference is contaminated (#5), it is
  non-monotonic with what we predict (#41), and its effect sizes shrink ~4× on the way into
  retention space (#42). The failure mode is gravity: a number that sits in the repo and
  updates every run reaches a summary, then a claim, and the qualifications fall off first.
  (#31,
  #14)

- **`D` is only identifiable where both plateaus are reachable — otherwise it is a
  conditional quantity and is labelled as one.** Fitting `D` traces retention from the
  neutral plateau to the ionised plateau; `D` is the gap. Acids (pKa ≈ 4) reach both inside
  pH 2.5–8 on any column. **Strong bases (pKa ≈ 9.5) need pH ≈ 11.5 for their neutral
  plateau**, which silica-based phases do not survive — so `k_neutral` would have to come
  from the LSER prediction, which is the quantity under test. **Measured `D` and
  baseline-conditional `D` are different evidence classes and are never pooled.** The same
  reachability rule is why permanently charged species (#39) sit outside any `D` measurement
  by construction: no neutral microspecies at any pH means the sigmoid has no upper plateau.
  (#30,
  #34)

- **The project commissions laboratory work, once, and only where public data provably
  cannot reach.** `D` owns 0.78–0.92 of the variance for ionised compounds (#21), has never
  been fitted to anything, and #4 established **no public pH sweep on a plain C18 with an
  ionisable-rich set exists** — so it is unreachable by analysis alone. The destination says
  *validated* architecture, and the project holds **zero measured retention data**: every
  number so far describes a prior's *shape*, never its correctness. Protocol specified at
  `experiments/ph-sweep-protocol.md` — two columns of contrasting silanol activity, 45 °C to
  match the WSU constants, constant ionic strength, WSU φ grid points so no interpolation
  engages. **A gate on a spend is bounded by design**: an unbounded gate on another
  acquisition turns "wait" into "never".
  (#30,
  #4,
  #45)

- **Thresholds are pre-registered, and `n` is counted in compounds.** Acceptance numbers are
  fixed **before** the data that judges them exists (`experiments/calibration-preregistration.md`),
  because a threshold chosen after the result is chosen to be passed — and #15 gates #19, so a
  gate that waits for its own evidence is not a gate. Amendments require a written reason
  recorded **before** the affected data is examined. **Pairs are not independent evidence**:
  one bad compound corrupts every pair it touches, so a binomial interval on pair counts
  understates the true spread 1.35× at 15 compounds and **2.40× at 71** — worsening as the
  study grows. **The compound is the unit of resampling**; every figure carries a
  compound-level bootstrap interval, and a margin counts only if it exceeds 2 SD of the paired
  difference.
  (#15,
  #27)

- **Coverage is not a calibration test; sharpness is.** Conformal reaches nominal coverage
  **by construction** — widen far enough and it always arrives — so a coverage threshold is a
  test that cannot fail. The binding quantity is the **conformal inflation `λ`** required to
  get there. Measured on the thin-slice ensemble, the confident-pair fraction falls 70% → 37%
  → 21% as `λ` goes 1.0 → 2.0 → 3.0: **the product falls off a cliff between 2 and 3, it does
  not degrade gracefully.** Pre-registered: **`λ` ≤ 2.0 calibrated; 2.0–2.5 degraded, every
  design-space claim carries the confident-pair fraction; `λ` > 2.5 the model is
  under-dispersed rather than calibrated and #19 rescopes** to a deterministic map plus
  optimiser. ⚠ A bar we may genuinely fail — under-dispersion measured at median **1.23**
  [0.98, 1.47], 90th percentile 2.45 (#42, re-run at n = 94 on 2026-08-16 after the
  Table S-2 extract fix; #35 supplies the compound-level bootstrap interval). Both are
  **floors**: all 94 anchor compounds are inside SoluteML's training corpus, so a `λ`
  near 1 from any in-corpus evaluation should be disbelieved on principle. `λ` is computed **per stratum**; a stratum with no calibration data
  never enters the pooled figure.
  (#15,
  #42,
  #19)

- **A gate has three verdicts, and the third is named in advance.** Pass, fail, and
  **inconclusive** — defined as *the achieved 2 SD interval exceeds the observed gain*, which
  does not block and does not count as a pass, and is reported with the sample size that
  would have settled it. Consequence for #27's gating axes: **pH keeps conditional accuracy**
  (23.5pp of headroom), but the **modifier gate moves to relative error reduction**, because
  its baseline is already 92–94% and binarising continuous retention into pairwise wins
  discards the signal — at `n` = 71 the binary statistic cannot call a real 17% error
  reduction unless the model/baseline error correlation reaches ~0.85, while relative error
  reduction calls it from ~0.7. Conditional accuracy stays **reported** there, without
  authority to fail the axis.
  (#15,
  #27)

- **Calibration strata merge by a count rule, and ionisation regime never merges.** Minimum
  stratum size: hard floor 19 (below it the 90% quantile is the sample maximum), practical
  floor 50. Below the floor strata merge upward in a fixed order — **confidence tertile
  first, then φ range, then column type** — and **ionisation regime is never merged**, having
  the largest measured effect (≈0.09 + 0.04 in `log k` for a neutral far from any pKa against
  ≈0.36 from pKa error inside the ±1.5-unit danger band) and carrying #34's two-regime claim.
  (#15,
  #42,
  #34)

- **A scouting run's objective is design-space uncertainty; it is selected by a declared
  proxy.** The product is a *shorter experimental programme*, so a run is worth exactly what
  it does to `P(all Rs ≥ target)`. **Not** parameter information: per-run parameters are
  shared and **largely cancel in elution order** (#14), so a D-optimal design on `k(φ)` would
  optimise the model's self-knowledge rather than the user's decision. The true objective is
  not computable at selection time — it is an expectation over unobserved outcomes through a
  conditioning step that is degenerate — so runs are ranked by a **proxy computable from the
  prior alone** (ensemble disagreement at the candidate condition, weighted toward critical
  pairs). **Objective and proxy are named separately and never merged**, so their divergence
  can be stated; same discipline as #41's criterion/diagnostic split.
  ⚠ **Every scouting run must carry a `t0` marker** (uracil or thiourea) — one extra peak,
  and it is what makes `t0(pH)` fittable at all. A scouting run without one still informs the
  design space but leaves the `t0` error term uncollapsed, so the prediction stays flagged as
  extrapolation across pH.
  (#18,
  #19,
  #25)

- **Trackability is a constraint, not an objective — and it is a count of assignable peaks,
  not a resolution threshold.** A compound is *assignable* when it is separated from its
  neighbours at Rs ≥ 1.5 in the median scenario **and** its order relative to them is
  confident at the calibrated inflation (#15). A run is admissible if a **majority of the
  mixture is assignable** [convention, to revisit on real chromatograms]. ⚠ **Critical Rs is
  the wrong admissibility criterion for a scouting run**: with 15 compounds isocratic on one
  column, P(critical Rs ≥ 1.5) ≈ 0.00–0.01 at *every* candidate condition, and a criterion
  every design violates cannot rank designs. Critical Rs remains the right acceptance
  criterion for a **method** (#10, #19) — the two must not be confused. Making trackability a
  constraint rather than a competing objective dissolves the conflict this ticket raised: the
  constraint fixes the admissible set, the proxy ranks within it, and no weighting between
  incommensurable quantities is required.
  (#18,
  #15,
  #10)

- **The user caps the budget; the system recommends stopping, and the programme ends in three
  states.** Cost is the user's, so the cap is theirs; whether another run would change
  anything is a model question, and leaving it to the user guarantees the full budget is spent
  by default — the waste this layer exists to remove. **Stop on decision stability**, not
  variance: when the probability that one more admissible run changes the recommended
  operating point falls below threshold. Uncertainty keeps shrinking long after the
  recommendation settles, so a variance rule buys numbers nobody acts on. End states, mirroring
  #15's verdicts: **converged** (stop, whether or not budget remains); **budget exhausted, not
  converged** (reported with an estimate of the further runs needed); **cannot converge**
  (uncertainty at the decision point exceeds what any admissible run reduces — a finding about
  the model, reported as one). ⚠ Floor: for any mixture containing ionisables the programme
  cannot be shorter than #15's pre-registered **3 runs at ≥ 2 pH units spacing** without
  trading away the two-regime claim.
  (#18,
  #15)

- **The column axis: 25, versioned, four tiers** — the committed set is the whole
  WSU-2019 axis (25 MeOH / 24 ACN columns; "5–15" formally retired as a pre-coverage
  guess); membership is data-forced, so selection criteria became documentation
  obligations. **The four-tier ladder: cold (WSU vector) → substituted (topology class,
  degraded) → fitting-only (scouting-funded, no a-priori claim) — a column is never
  refused, only a cold prediction is.** Fitting-only is first-class: cold map renders
  no-claim, cross-column behaviour seeds the ensemble as uncalibrated context, no graded
  a-priori metric ever counts it, and promotion to cold requires characterisation, never
  inference. **Knox `A,B,C` are engineering priors** by particle class, never measured
  properties of the individual column; their ±10–30% plate-count cost lives in the
  per-run class and gives the resolution map an N-uncertainty floor retention accuracy
  cannot remove (#15). **Data-sheet invariant: a data-sheet field describes evidence
  about the column; it does not automatically become a model input.**
  (#13,
  #36,
  #9)

- **Anchor lineage: WSU/Poole** — the authoritative descriptor set at the QSPR→LSER
  interface; the committed system constants are interpreted only against it (they *are*
  the least-squares optima for Table S-2 descriptors — the constants choose the lineage).
  **Corpus size is not an authority criterion**: provenance, compatibility with the
  committed constants, and independent measurement establish authority; coverage is a
  separate data-availability problem. The commitment is **unconditional** — it stands on
  the 94 measured compounds in hand; Poole 2020 / WSU-2025 extend it (highest-priority
  acquisition, #22). #35 failure does not reopen it — the fallback is stock SoluteML plus
  the measured bias carried as an *explicit systematic component* (never absorbed into
  noise), with #31's joint-layer validation as the principal correctness test.
  (#23)

- **The corpus fence (architectural invariant)** — *the Absolv-family corpus may teach
  the model what molecules are like; it may never testify about whether the model is
  right.* Permitted: representation pretraining; MIT SoluteDB as the auditable copy;
  contamination quantification **as diagnostic metadata only** (a bias report, never
  adjusted training labels — the diagnostic must not become an indirect influence path);
  the R1 provenance inquiry. Forbidden: ground-truth fine-tuning, validation/testing
  authority, calibration authority. Evaluation sets contain only claimable predictions
  (no masked, `NO_SUPPORT`, or opt-in rows — the #40/#24 rules apply to evaluation
  eligibility too). No future experiment may silently turn an analogy-estimated
  descriptor corpus into evidence of model correctness.
  (#23,
  #5,
  #28)

- **`NO_SUPPORT` vs `MASK`, and the interpolation invariant** — `k(φ)` comes from PCHIP
  interpolation of the published constants (#24; LSS rejected: affine premise fails at
  2–13 published SDs of curvature, and its out-of-window error exceeds the whole QSPR
  budget). Support semantics: *interior unmasked absence may be interpolated; an interior
  masked region may never be crossed* (missing ≠ masked, and code distinguishes them).
  Beyond a support edge → **`NO_SUPPORT` → NaN, no opt-in** — this refines the opt-in
  rule: **opt-in applies to an unavailable prediction; it never manufactures one where
  the model has no underlying value** (steric failure still has an LSER value; past the
  φ edge there is none, and extending the last slope is a new model assumption the
  curvature data already falsifies). `NO_SUPPORT` is an *acquisition* problem (erased
  when measurements extend the range); `MASK` is a *physical-validity* constraint
  (removable only by a new measurement-derived determination on the specific physical
  column). Both are hard-infeasible and render as no-claim, with the reason preserved.
  **Invariant: interpolation fills empirical holes inside supported, unmasked regions; it
  never creates support, crosses invalidity, or extrapolates beyond support.**
  (#24,
  #40)

- **Three-tier sampler correlation structure** — measured by #37, and the tiers are not
  interchangeable. (1) **Within-column**: the six LSER system constants are drawn
  **jointly** — independent coefficient draws overstate `Var(log k)` ~16× (c–v −0.65,
  e–s −0.66), which is materially wrong uncertainty, not a refinement; lives in the
  per-run error class. (2) **Across compounds**: scaffold-family residual correlation is
  material (ρ ≈ +0.42 MeOH / +0.54 ACN, family-dependent, anilines ≈ 0) — joint compound
  draws with **per-family ρ**, never one global number. (3) **Across ensemble members**:
  the descriptor 4×4 is second-order (full/diag 0.96 ACN / 0.89 MeOH) — carried for free
  from the 25-member ensemble, no architecture built on it. True per-compound LSER
  residual structure remains an **acquisition gap** (raw log k; #26/#22), never inferred.
  (#37,
  #14)

## Conventions

- **No bare order-accuracy number may be emitted.** The evaluation harness has no code path
  returning an accuracy figure unaccompanied by its applicable baseline, its `n`, its interval and
  its condition axis — enforced by the function signature and asserted in a test, not left to
  reviewer discipline. Every failure this guards against is one of *omission*, and conventions
  decay under deadline pressure. Success is measured as **relative error reduction**
  `(baseline_error − model_error)/baseline_error` plus conditional accuracy, **per axis**; absolute
  percentage-point margins are rejected because headroom differs by an order of magnitude across
  axes and a fixed margin silently encodes a scope decision.
  (#27)

- **The claim verb follows from the reference type, and the reference travels with the
  number.** #27's four fields (baseline, `n`, interval, axis) gain two: **the reference's
  identity and the reference's coverage.** If the reference is model output — ACD Absolv or
  any Absolv-lineage descriptor — the only permitted verb is **agrees with**, and the
  reference is named in the same sentence. If the reference is measurement, the verb is
  **accurate to**. A figure referenced to corpus descriptors can never take the second verb
  however good it looks: agreeing with that corpus partly means mimicking the tool we are
  trying to beat (#5), so "accurate" would assert the opposite of what was measured.
  Coverage is carried because the anchor is narrow — a grade earned on 94 neutral solutes
  (median V 1.06, zero acids) invites application to compounds outside its entire range
  unless the boundary is visible at the point of use.
  (#31,
  #28,
  #27)

- **Corrector fitting order is fixed, and asserted in a test.** #35's WSU re-fit (global,
  lineage) → solute-vector correction on its residual → `D` last (per-compound,
  scouting-informed). **Never fitted jointly.** All three can absorb the same discrepancy for
  an ionisable compound on a characterised column, so without a pinned order the last one
  fitted silently claims the credit and the uncertainties cannot be added. Uncertainty is
  reported from the final stage only.
  (#36,
  #34)

- **"`e` collapses to zero" is an acetonitrile result, not a universal one.** Poole 1600
  Table 3: `e` runs 0.384 → 0.157 across the methanol groups but sits at ~0.03 (and changes
  sign) in acetonitrile. The `eE` term is therefore kept, with `E`'s uncertainty propagated
  rather than the term deleted — it self-cancels in ACN anyway. ⚠ `E` is both our
  worst-predicted (ra² 0.346) and most contaminated descriptor, so its error is plausibly
  **biased**, and propagating variance does not fix bias. **#35 has now measured it:** bias
  **+0.035** [+0.004, +0.065] overall, and the offset is essentially a weak-base effect —
  **+0.152** [+0.077, +0.236] on the anchor's 15 weak bases against +0.011 on its 79
  neutrals. `E`'s MAE against measurement (0.106) exceeds SoluteML's own published
  *substructure*-split MAE (0.084), on compounds it trained on. Diagnostic metadata under
  #23's corpus fence: reported, never subtracted as a correction.
  (#36,
  #12,
  #6)

- **WSU-2019 uses `log₁₀`, not `ln`.** Settled from its Eqs. 5 and 7 and its figure axes.
  The 2.303× factor is live against any `ln k` source — including this project's blueprint
  paper, whose `S_S` is 2.303× the conventional `S`.
  (#22,
  #2)

- **Ionised compounds are outside the solvation parameter model.** Stated outright in
  Poole 1600 (2019) §3.4, along with steric resistance and cation exchange being
  unparameterised. Cation exchange is **acetonitrile-specific** — suppressed for
  methanol–water on the *same* columns — and its relative contribution **grows with φ**,
  so it is a modifier × column × φ property, not a column property alone.
  ⚠ **Unresolved evidence conflict:** WSU-2019 Table S-1 flags Fluophase-RP (a PFP-type
  phase) as electrostatically active **in methanol too**. Logged on #39; until
  adjudicated it neither expands nor retracts the `D` prior, and the
  acetonitrile-specific claim should be read as established for C18-type phases.
  (#34,
  #12,
  #40)

- **WSU system constants are a 45 °C dataset.** Poole, *J. Chromatogr. A* 1600 (2019)
  112–126, Experimental: "All measurements were made with a column temperature of
  45 °C", for 10–70% (v/v) methanol, acetonitrile or THF. A prediction built on these
  constants is silently a 45 °C prediction unless temperature is modelled. Related:
  **pore dewetting below 30% v/v methanol** breaks the system maps, a hard validity
  floor on the φ axis. **The floor is not one number:** 10% v/v is generally unusable
  (incomplete wetting), Kinetex Biphenyl fails specifically at 10% methanol, and the
  paper's own Table 3 constants begin at 30%. The paper's only temperature qualification
  is that a small temperature change has a small effect **for neutral compounds** —
  unquantified, with no `dc/dT`. Do not build a temperature correction on it.
  (#22,
  #11,
  #24)

- **Temperature has four roles, and they never blur.** *Nominal/design condition*:
  45 °C — what the WSU constants physically are. *Observed/fitted run temperature*: a
  **nuisance parameter**, recovered per run as a column-vector offset (all T-dependence
  lives in the phase parameters, none in the descriptors — Ulrich 2026; same mechanism
  family as #36's vector shifting). *Perturbation*: ±1–2 °C drift is a robustness factor
  in the design space (#19). *Optimisation axis*: **refused** — the fitted temperature
  never becomes a freely optimised design variable merely because the model estimates
  it; re-opening requires an explicit future ticket.
  (#11)

- **Modifier scope: the two binaries, modifier-conditioned everywhere.** Methanol/water
  and acetonitrile/water are both in scope with near-parity support; model structure is
  conditioned on modifier (`eE` methanol-only, ACN-conditioned cation exchange,
  per-modifier φ floors and budgets). **THF is outside the current validated regime**
  (2 columns) — a support statement, not a permanent exclusion; WSU-2025 is the route
  back. Ternary blends out. Additives absorb into the pH/buffer treatment **except TFA**,
  which is split out as an unresolved ion-pairing question — deferred with a named
  mechanism, never dismissed.
  (#12)

Numerical conventions that silently corrupt results if unpinned are indexed in the
map's Notes (#1):
`ln` vs `log₁₀` (the blueprint's `S_S` is 2.303× conventional `S`), the `s_s`pH scale,
`ln 10` inside Snyder's `G`, four QSPR targets (V computed, not L), descriptor lineage,
and the RepoRT pinned commit. Record each here in full as the owning ticket closes.

- **RepoRT is pinned at `9de8d603377bbb6cc0f74e250eb7533fd8874df1`** (short `9de8d60`,
  authored 2026-07-09), and the clone is never vendored — `data/report/README.md` has the
  recipe, `data/report/loader.py` holds the pin and the patch for the documented 0310–0341
  name-swap bug. ⚠ **No `k` can be derived from RepoRT at all**: zero isocratic datasets, no
  dwell-volume field anywhere in the schema, and `column.t0` imputed from geometry rather
  than measured (only 47 distinct values across 421 datasets, 48 of them zero). It is a
  compound and condition census, not a retention source.
  (#29,
  #45)

- **Corpus coverage is a stated tier, not an adjective.** Against SoluteML's training corpus
  (SoluteDB, Zenodo 5792296, CC BY 4.0), a compound is `in_corpus` (InChIKey skeleton
  present), `scaffold_seen` (Murcko scaffold present), `novel_in_envelope`, or
  `out_of_envelope` (charged, or an element outside H C N O S P F Cl Br I). Measured over
  RepoRT's 18,872 non-SMRT compounds: **9.0% / 22.0% / 66.9% / 2.1%**, but **23.6% /
  27.2% / 47.7% / 1.5%** when weighted by the 103,534 retention records — memorised
  chemistry is measured about three times as often as novel chemistry, so any evaluation
  sampling *records* rather than *compounds* inherits that flattery. **The WSU-94 anchor is
  100% `in_corpus`.** The tier describes evidence strength and gates nothing.
  (`experiments/solutedb-census/`,
  #35)

- **Table S-2's `class` column is the anchor's own weak-base designation**, and it exists
  because the heading that carries it used to be glued onto a compound name — silently
  costing one compound in every name join (2026-08-16 fix, `sources/wsu-lser/extract_wsu.py`,
  `EXECUTION-PLAN.md`). 79 neutral, 15 weak base. Any analysis stratifying the anchor uses
  this column rather than re-deriving the split.
  (#43)
