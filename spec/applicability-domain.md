# The applicability domain

Implementation: [`model/applicability.py`](../model/applicability.py).
Checks: `python3 model/test_applicability.py` (24) and `python3 model/test_strata.py` (16),
stdlib only.

This is the object the eval harness names as missing ("the harness does not yet refuse
out-of-domain rows because the domain object lives upstream") and that #40's block rule,
#43's boundary and #33/#34's refusals all assume exists.

## The design constraint

#43 decided that "bulky" is **two predicates and never one**: a coverage boundary defensible
from data we hold, and a steric-exclusion mechanism that is an unverified physical assumption
pending WSU-2019 ref 29. The module enforces that separation structurally rather than by
convention.

**No reason may exist without declaring its warrant.** `Reason.warrant` is mandatory:

| warrant | means |
|---|---|
| `measured` | rests on data in this repo, and names it |
| `declared` | rests on a source's own declaration of its scope |
| `assumed` | rests on a physical argument nobody has tested |

`Decision.assumptions` exposes the `assumed` reasons separately and `render()` prints them
under `UNVERIFIED ASSUMPTIONS -- not evidence:`. A caller cannot quote a steric refusal as
though it were evidence, in the same way `Figure` stops a bare accuracy number escaping.

**The steric predicate ships disabled.** Ref 29 was never obtained, so the predicate cannot
be derived from provenance data as #40 required. With it off, a registered class is *noted*
on the decision and acted on by nobody; turning it on produces a refusal that always carries
the `assumed` warrant. Enforcing an untested physical claim by default would be the exact
failure the two-part split exists to prevent.

**And the predicate is scoped, because #40 scoped it.** The block applies to a compound in
an affected **(column, modifier, φ-range)** — never to a compound everywhere. That scope is
*measured*: WSU-2019 Table S-1 flags **8 of its 25 columns**, with modifiers per column and
a φ threshold for two of them.

| | |
|---|---|
| flagged columns | Ascentis C18, Betasil C18, Discovery HS C18, Discovery HS F5, Fluophase-RP, HyPURITY C18, Synergi Fusion-RP, XTerra MS C18 |
| both modifiers | Ascentis, Betasil, HyPURITY, XTerra MS |
| methanol only | Discovery HS C18, Discovery HS F5, Fluophase-RP, Synergi Fusion-RP |
| quantified φ | Betasil and XTerra MS (`< 40%`); the other six say only "water-rich" |

So the two halves of the rule carry **different warrants and the reasons say so**: the scope
is `declared` by the source, the compound class is `assumed`. An unquantified φ range is
treated as the whole axis (conservative) and the reason states that the source gave no
number, rather than a threshold somebody invented. With no column or modifier supplied the
result is `degraded` with `steric_scope_unknown` — refusing without knowing where would be
the over-refusal this scoping exists to prevent, and silently passing would be the
under-refusal.

⚠ **Fluophase-RP's note names hydrogen-bond acids, not bulk.** The register is not one
structural axis, so the source note is preserved on each scope entry and no code assumes
`steric_class == "bulky"`.

⚠ **Found 2026-08-17, after this module first shipped.** The first version took
`steric_class` and ignored column, modifier and φ entirely — it would have refused on the
**17 columns the source never flagged**, over-refusing by more than twice the affected set,
using data (`sources/wsu-lser/wsu2019-column-caveats.csv`) that was already in the repo. The
gap surfaced on re-reading #43's own first comment, which had stated the 8-of-354-fits scope
from the start.

## The tiers

`refused` — no descriptor path exists, so no claim is possible:

- an element outside SoluteML's declared envelope (H C N O S P F Cl Br I)
- a charged or permanently charged species (#33/#39)
- the ionisation provider refused the compound (#33)

`degraded` — a prediction is possible but measured support is thin or absent:

- a predicted descriptor beyond the anchor's observed maximum (`outside_anchor_range`)
- a predicted descriptor above the anchor's p95 (`anchor_upper_tail`)

`in_envelope` — inside the declared envelope and within the anchor's measured range.

The anchor is WSU-2019 Table S-2, n = 94, and its range is the edge of where a descriptor
claim has any measured support at all. p95 and max are both recorded: p95 is #43's coverage
boundary, max is where the anchor stops existing.

| | E | S | A | B | V |
|---|---|---|---|---|---|
| p95 | 1.588 | 1.754 | 0.923 | 0.769 | 1.614 |
| max | 1.976 | 2.214 | 1.148 | 1.388 | 2.622 |

Membership is **structural, not scored**. #40's finding is that a stratum with zero
calibration support cannot be conformally calibrated, so a confidence number would be a
fiction where a refusal is the truth.

## The evaluation firewall, and a mistake worth recording

`evaluation_strata()` returns, per row, the stratum it may be scored in or `None`.

The first version of this module excluded **degraded** rows from evaluation entirely, on a
plain reading of #40's "only in-domain predictions may enter any metric". Running it over the
thin-slice compounds is what caught the error: **all nine drugs degrade** against the anchor
(#31 — atenolol `B` 2.11 and ketoprofen `S` 2.28 exceed its maxima), the quaternary ammonium
is refused, and only the six neutral anchors pass. So the strict reading would have left the
pH sweep with **zero evaluable compounds and reported it as a clean pass** — and #40's own
finding is that the degraded stratum is precisely the one whose calibration support we lack
and most need.

Degraded rows therefore enter metrics, in their own stratum. What is forbidden is **pooling**:
the function returns strata rather than a boolean mask so that combining a degraded residual
with an in-envelope headline is a deliberate act, governed by #15 §4's merge rule. Refused
rows never enter. Opt-in residuals — volunteered by a user for a compound we refused or
degraded — never enter, whatever their tier, since accepting them would let the excluded
population back in through the one door that stays open.

## Corpus tier: carried, gates nothing

A decision may carry the compound's SoluteDB corpus tier from
`experiments/solutedb-census/` (`in_corpus` / `scaffold_seen` / `novel_in_envelope` /
`out_of_envelope`). It describes **evidence strength, not permission**: the census tiers are
structural facts about what the model has seen, and nothing in this repo measures how much
worse the novel tier actually is. Asserted by test.

## What this does not do

- It does not score confidence, and it must not be made to.
- The anchor's coverage range is the only measured descriptor reference we have and it is
  94 neutral small solutes, all inside SoluteML's training corpus. The boundary is a
  statement about **our evidence**, not about chemistry.
- The steric mechanism remains undefined. This module makes its absence explicit and
  auditable; it does not substitute for ref 29.

## Wired in (Phase 3, 2026-08-16)

The object no longer stands alone. Three seams are closed:

**The eval harness** (`prototype/eval-harness/harness.py`) takes strata as plain strings —
it deliberately does not import the model layer, so it stays numpy-only. `score_by_stratum()`
is the domain-gated entry point and returns **one `Figure` per stratum, never a pooled one**;
`require_single_stratum()` raises `MixedStrataError` rather than averaging across strata;
and `Figure` gained a `stratum` field rendered next to `n`, because a degraded figure read
without its stratum *is* the in-envelope claim.

The pair rule is the part that needed thought. A pair is only as clean as its **worse**
compound, and a pair touching a refused compound cannot be scored at all. The arithmetic
is why: in a 15-compound slice one degraded compound appears in **14 of the 105 pairs**, so
without the rule a single unhandled compound contaminates 13% of the pairs while every
compound-level check still passes. Asserted in `test_domain_firewall`.

**The ensemble spec** (`spec/ensemble.py`) gained an optional `domain[n]` of `Decision`
objects. Refusal is now expressed twice on purpose — as absence of data, which is what a
numeric consumer needs, and as a reason, which is what a person needs — and `validate()`
raises if the two disagree in either direction: a compound with no distribution but a
permissive decision ("a refusal with no recorded reason"), or a refused compound that
carries one. `refusal_reasons()` returns the rendered reasons, empty when no domain was
supplied, which is the honest answer rather than an absence of refusals.

**The stratum producer** (`model/strata.py`) fills the spec's one unimplemented MUST.
Method-independent per #42 (the plain descriptor-SD sum, not the LSER-propagated spread),
with ionisation regime as an axis that never merges per #15 §4, and **tertile boundaries
frozen from the anchor** so a compound's class cannot change with the batch it arrives in.

`model/correctors.py` closes an unrelated gap found in the same audit: `CONTEXT.md` said
#36's corrector fitting order was "asserted in a test" and no such test existed.
