"""The applicability domain: may we make a claim about this compound at all?

The object the eval harness names as missing ("the harness does not yet refuse
out-of-domain rows because the domain object lives upstream") and that #40's
block rule, #43's boundary and #33/#34's refusals all assume exists.

THE DESIGN CONSTRAINT, which is the whole point of the module. #43 decided that
"bulky" is two predicates and not one: a COVERAGE boundary we can defend from
data we hold, and a STERIC-EXCLUSION mechanism that is an unverified physical
assumption pending ref 29. Those must never collapse into one number, because a
coverage statement and a mechanism claim carry different warrants.

So no reason may exist without declaring which it is. `Reason.warrant` is
mandatory and takes one of:

    "measured"   -- rests on data in this repo, and names it
    "declared"   -- rests on a source's own declaration of its scope
    "assumed"    -- rests on a physical argument nobody has tested

`Decision.assumptions` exposes the "assumed" reasons separately, and
`Decision.render()` prints them under a heading that says so. A caller cannot
accidentally quote a steric refusal as if it were evidence, in the same way
`Figure` stops a bare accuracy number escaping.

WHAT THIS DOES NOT DO. It does not score confidence. Membership is structural:
a compound is inside the envelope or it is not, inside the anchor's coverage or
it is not. #40's finding is that a stratum with zero calibration support cannot
be conformally calibrated, so a low confidence number would be a fiction where
a refusal is the truth.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Iterable, Literal, Optional

Warrant = Literal["measured", "declared", "assumed"]
Tier = Literal["in_envelope", "degraded", "refused"]

# SoluteML's training domain, as Chung 2022 declares it: neutral solutes,
# these elements only. A structural test, not a threshold (#35 thread).
ENVELOPE_ELEMENTS = frozenset(
    {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"})

# The WSU-2019 Table S-2 anchor's observed range, n = 94. This is the ONLY
# measured-descriptor reference the project has, so its range is the edge of
# where a descriptor claim has any measured support at all. p95 and max both
# recorded: p95 is the coverage boundary #43 decided on, max is where the
# anchor stops existing.
ANCHOR_P95 = {"E": 1.588, "S": 1.754, "A": 0.923, "B": 0.769, "V": 1.614}
ANCHOR_MAX = {"E": 1.976, "S": 2.214, "A": 1.148, "B": 1.388, "V": 2.622}
ANCHOR_N = 94
ANCHOR_IDENTITY = ("WSU-2019 Table S-2, 94 neutral solutes, one lab, one "
                   "protocol; 94/94 inside SoluteML's training corpus")


CAVEATS_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "sources", "wsu-lser", "wsu2019-column-caveats.csv")


@dataclass(frozen=True)
class StericScope:
    """Where a steric caveat applies, from WSU-2019 Table S-1's own register.

    #40 scopes the block to an affected **(column, modifier, phi-range)**, not to
    a compound alone. Only **8 of 25 columns** carry a `steric_repulsion` flag, so
    a predicate that ignored scope would refuse on the 17 columns the source
    never flagged -- over-refusing by more than twice the affected set.
    """

    column: str
    modifiers: frozenset[str]
    phi_max_pct: Optional[float]   # None = the source gave no number
    note: str
    reference: str

    def applies(self, modifier: str, phi_pct: Optional[float]) -> bool:
        if modifier.lower() not in self.modifiers:
            return False
        if self.phi_max_pct is None:
            # "water-rich" with no number. Treating the whole axis as affected
            # is the conservative reading and it is recorded as unquantified
            # rather than silently turned into a threshold somebody invented.
            return True
        return phi_pct is None or phi_pct < self.phi_max_pct


@lru_cache(maxsize=1)
def steric_register(path: str = CAVEATS_CSV) -> dict[str, StericScope]:
    """{column: StericScope} for every column the source flags."""
    if not os.path.exists(path):
        raise FileNotFoundError(
            2, "WSU-2019 Table S-1 caveat register is not redistributed; "
            "build it with sources/wsu-lser/extract_wsu.py "
            "(see sources/wsu-lser/README.md)", path)
    out: dict[str, StericScope] = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["caveat"] != "steric_repulsion":
                continue
            m = re.search(r"(?:less than|<)\s*(\d+)\s*%", r["note"])
            out[r["column"]] = StericScope(
                column=r["column"],
                modifiers=frozenset(x.strip().lower()
                                    for x in r["modifiers"].split(";") if x.strip()),
                phi_max_pct=float(m.group(1)) if m else None,
                note=r["note"].strip(),
                reference=r["reference"].strip())
    return out


@dataclass(frozen=True)
class Reason:
    """Why the domain object reached its verdict. Warrant is not optional."""

    code: str
    detail: str
    warrant: Warrant
    tier: Tier
    source: str = ""

    def __post_init__(self) -> None:
        if self.warrant not in ("measured", "declared", "assumed"):
            raise ValueError(
                f"warrant must be measured|declared|assumed, got {self.warrant!r}")
        if self.tier not in ("in_envelope", "degraded", "refused"):
            raise ValueError(f"bad tier {self.tier!r}")
        if not self.code or not self.detail:
            raise ValueError("code and detail are mandatory")

    def render(self) -> str:
        src = f" [{self.source}]" if self.source else ""
        return f"{self.code}: {self.detail} ({self.warrant}){src}"


@dataclass(frozen=True)
class Decision:
    """A domain verdict that cannot be separated from its warrants."""

    tier: Tier
    reasons: tuple[Reason, ...]
    corpus_tier: Optional[str] = None      # in_corpus | scaffold_seen | novel | unknown

    @property
    def assumptions(self) -> tuple[Reason, ...]:
        """The reasons resting on untested physical argument. #43's second half."""
        return tuple(r for r in self.reasons if r.warrant == "assumed")

    @property
    def rests_on_assumption(self) -> bool:
        return bool(self.assumptions)

    @property
    def may_predict(self) -> bool:
        return self.tier != "refused"

    @property
    def may_enter_evaluation(self) -> bool:
        """#40's firewall: refused rows never enter a metric.

        Degraded rows DO enter, and the first version of this module had them
        excluded -- which is wrong in a way worth recording, because applying
        it to the thin-slice compounds is what exposed it: all nine drugs
        degrade against the anchor (#31), so excluding degraded rows would have
        left the pH sweep with zero evaluable compounds and reported it as a
        clean pass. Worse, #40's own finding is that a stratum with no
        calibration support cannot be conformally calibrated -- so the degraded
        tier is exactly the one whose residuals we most need.

        What must not happen is POOLING. See `evaluation_stratum`.
        """
        return self.tier != "refused"

    @property
    def evaluation_stratum(self) -> Optional[str]:
        """Which stratum this row may be scored in, or None if it may not.

        Returned instead of a bare boolean so a caller cannot pool a degraded
        residual into an in-envelope headline without naming the stratum it
        came from. #15 §4's merge rule governs what may combine with what.
        """
        return None if self.tier == "refused" else self.tier

    def render(self) -> str:
        head = f"{self.tier.upper()}"
        if self.corpus_tier:
            head += f" (corpus: {self.corpus_tier})"
        lines = [head] + [f"  - {r.render()}" for r in self.reasons]
        if self.assumptions:
            lines.append("  UNVERIFIED ASSUMPTIONS -- not evidence:")
            lines += [f"    - {r.code}: {r.detail}" for r in self.assumptions]
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.render()


def _worst(tiers: Iterable[Tier]) -> Tier:
    order = {"in_envelope": 0, "degraded": 1, "refused": 2}
    return max(tiers, key=lambda t: order[t], default="in_envelope")


def assess(
    *,
    elements: Iterable[str],
    formal_charge: int,
    descriptors: Optional[dict[str, float]] = None,
    permanently_charged: bool = False,
    ionisation_refused: bool = False,
    corpus_tier: Optional[str] = None,
    steric_class: Optional[str] = None,
    steric_assumption_enabled: bool = False,
    column: Optional[str] = None,
    modifier: Optional[str] = None,
    phi_pct: Optional[float] = None,
) -> Decision:
    """Decide whether a claim may be made about this compound.

    `descriptors` are the PREDICTED (E, S, A, B, V) -- the coverage test asks
    whether the prediction lands where measured support exists, which is
    answerable before any measurement of this compound exists.

    `steric_class` names a class from the #46/#40 register. It only ever
    produces a refusal when `steric_assumption_enabled` is True, and that
    refusal is always warranted "assumed". Default off: #43 has not been able
    to define the predicate from provenance data (ref 29 unobtained), so
    shipping it on by default would enforce an untested physical claim.

    `column`, `modifier` and `phi_pct` scope that predicate. #40 blocks a
    compound in an affected **(column, modifier, phi-range)**, and the affected
    set is measured -- WSU-2019 Table S-1 flags 8 of 25 columns. The two halves
    of the rule therefore carry different warrants and the reasons say so: the
    scope is `declared` by the source, the compound class is `assumed`.
    """
    reasons: list[Reason] = []

    bad = sorted(set(elements) - ENVELOPE_ELEMENTS)
    if bad:
        reasons.append(Reason(
            code="element_outside_envelope",
            detail=f"contains {', '.join(bad)}; no descriptor path exists",
            warrant="declared", tier="refused",
            source="Chung 2022 declared training domain"))

    if formal_charge != 0 or permanently_charged:
        reasons.append(Reason(
            code="charged_species",
            detail="training domain is neutral solutes only; a permanently "
                   "charged compound has no descriptor path (#33/#39)",
            warrant="declared", tier="refused",
            source="Chung 2022 declared training domain"))

    if ionisation_refused:
        reasons.append(Reason(
            code="ionisation_refused",
            detail="the ionisation provider refused this compound, so the "
                   "neutral fraction is unknown (#33)",
            warrant="measured", tier="refused",
            source="model/ionisation.py"))

    if descriptors:
        beyond_max = [k for k, v in descriptors.items()
                      if k in ANCHOR_MAX and v > ANCHOR_MAX[k]]
        beyond_p95 = [k for k, v in descriptors.items()
                      if k in ANCHOR_P95 and v > ANCHOR_P95[k]
                      and k not in beyond_max]
        if beyond_max:
            reasons.append(Reason(
                code="outside_anchor_range",
                detail=f"{', '.join(sorted(beyond_max))} beyond the anchor's "
                       f"observed maximum; no measured support at this value",
                warrant="measured", tier="degraded",
                source=f"{ANCHOR_IDENTITY}, n={ANCHOR_N}"))
        if beyond_p95:
            reasons.append(Reason(
                code="anchor_upper_tail",
                detail=f"{', '.join(sorted(beyond_p95))} above the anchor's "
                       f"p95; supported by a handful of compounds at most",
                warrant="measured", tier="degraded",
                source=f"{ANCHOR_IDENTITY}, n={ANCHOR_N}"))

    if steric_class:
        scope = steric_register().get(column) if column else None
        in_scope = scope.applies(modifier, phi_pct) if (scope and modifier) else False
        if not steric_assumption_enabled:
            reasons.append(Reason(
                code="steric_class_noted",
                detail=f"class {steric_class!r} is in the register but the "
                       f"steric predicate is disabled; not acted on",
                warrant="assumed", tier="in_envelope",
                source="#43 -- unverified physical assumption"))
        elif not column or not modifier:
            # Refusing without knowing where would be the over-refusal this
            # scoping exists to prevent; so would silently passing.
            reasons.append(Reason(
                code="steric_scope_unknown",
                detail=f"class {steric_class!r} may be sterically excluded, but "
                       f"no (column, modifier) was supplied, so whether the "
                       f"caveat applies here cannot be decided",
                warrant="declared", tier="degraded",
                source="WSU-2019 Table S-1 register"))
        elif in_scope:
            phi_txt = (f"phi < {scope.phi_max_pct:.0f}%"
                       if scope.phi_max_pct is not None
                       else "phi range unquantified in the source "
                            "('water-rich'), so the whole axis is treated as "
                            "affected")
            reasons.append(Reason(
                code="steric_exclusion_assumed",
                detail=f"class {steric_class!r} on {column} with {modifier} "
                       f"({phi_txt}); the source flags this column and modifier, "
                       f"but the excluded-compound lists (WSU-2019 ref "
                       f"{scope.reference}) were never obtained, so which "
                       f"compounds belong to the class is untested. Source note: "
                       f"{scope.note!r}",
                warrant="assumed", tier="refused",
                source="#43 -- scope declared by WSU-2019 Table S-1, "
                       "compound class assumed"))
        else:
            where = (f"{column} carries no steric flag" if scope is None
                     else f"{column}'s steric flag does not cover "
                          f"{modifier} at this phi")
            reasons.append(Reason(
                code="steric_out_of_scope",
                detail=f"class {steric_class!r} noted, but {where}; only 8 of "
                       f"the 25 registered columns are flagged and the block is "
                       f"scoped to (column, modifier, phi-range) per #40",
                warrant="declared", tier="in_envelope",
                source="WSU-2019 Table S-1 register"))

    if not reasons:
        reasons.append(Reason(
            code="in_envelope",
            detail="inside the declared training envelope and within the "
                   "anchor's measured descriptor range",
            warrant="measured", tier="in_envelope",
            source=f"{ANCHOR_IDENTITY}, n={ANCHOR_N}"))

    return Decision(tier=_worst(r.tier for r in reasons),
                    reasons=tuple(reasons), corpus_tier=corpus_tier)


def evaluation_strata(decisions: Iterable[Decision],
                      opt_in: Iterable[bool] = ()) -> list[Optional[str]]:
    """#40's firewall: the stratum each row may be scored in, or None.

    Returns strata rather than a boolean mask so that pooling is a deliberate
    act. `opt_in` marks residuals a user volunteered for a compound we refused
    or degraded; those are never evaluation data, whatever their tier --
    accepting them would let the excluded population back in through the one
    door that stays open.
    """
    decisions = list(decisions)
    opt = list(opt_in) + [False] * (len(decisions) - len(list(opt_in)))
    return [None if o else d.evaluation_stratum
            for d, o in zip(decisions, opt)]
