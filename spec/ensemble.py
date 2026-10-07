"""Validator for the scenario-ensemble interface (#44, amending #14).

`scenario-ensemble.md` is the specification; this makes it checkable. A prose
interface that #16/#17 code is written against will drift, so every MUST in that
document is asserted here and `validate()` refuses a non-conforming object.

Deliberately dependency-light: numpy only, no framework, so it can be imported
by a producer, a consumer, or a test without dragging anything in.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

MODIFIERS = ("MEOH", "ACN", "THF")
PH_SCALES = ("SS", "WW")
T0_SOURCES = ("MEASURED", "FITTED", "NOMINAL")

# Vendor extra-column volumes (uL) from sources/instrument-volumes/. Used as the
# fallback DISTRIBUTION when instrument_id is null - never zero, which is the one
# value it certainly is not.
EXTRA_COLUMN_UL_RANGE = (5.0, 34.0)


class SpecViolation(ValueError):
    """Raised when an ensemble does not conform to scenario-ensemble.md."""


@dataclass(frozen=True)
class MethodCard:
    """The conditions one operating point was computed at.

    Rides WITH the array, never beside it: a scenario ensemble whose conditions
    are not part of the object invites silent condition mismatches.
    """
    phi: float
    modifier: str
    ph: float
    ph_scale: str                 # no default, deliberately - see spec
    temperature_c: float
    column_id: str
    t0_min: float                 # per operating point (#25), not a global scalar
    t0_source: str
    flow_ml_min: float            # required: plate count is a METHOD property
    instrument_id: str | None = None
    gradient: dict | None = None  # None means isocratic

    def __post_init__(self):
        if self.modifier not in MODIFIERS:
            raise SpecViolation(f"modifier must be one of {MODIFIERS}, got {self.modifier!r}")
        if self.ph_scale not in PH_SCALES:
            raise SpecViolation(
                f"ph_scale must be one of {PH_SCALES}, got {self.ph_scale!r}. "
                "It has no default: #20 carries the w_w->s_s correction as a named "
                "error term, and a defaulted scale is how that error goes silent.")
        if self.t0_source not in T0_SOURCES:
            raise SpecViolation(f"t0_source must be one of {T0_SOURCES}")
        if not 0.0 <= self.phi <= 1.0:
            raise SpecViolation(f"phi is a fraction 0-1, got {self.phi}")
        if self.flow_ml_min <= 0:
            raise SpecViolation(
                f"flow_ml_min must be positive, got {self.flow_ml_min}. It has no "
                "default because plate count depends on it: N is a property of the "
                "METHOD, not of the column, and a defaulted flow silently fixes "
                "the peak widths of every prediction.")
        if self.t0_min <= 0:
            raise SpecViolation(
                f"t0_min must be positive, got {self.t0_min}. RepoRT records t0=0 for "
                "48 datasets as a missing-value marker; that is not a value.")

    @property
    def requires_t0_error_term(self) -> bool:
        """A NOMINAL t0 obliges the consumer to add the per-run t0 error term and
        flag the prediction as extrapolation across pH (#25)."""
        return self.t0_source == "NOMINAL"

    @property
    def requires_extra_column_distribution(self) -> bool:
        """With no instrument, extra-column variance spans the vendor range."""
        return self.instrument_id is None


@dataclass(frozen=True)
class ScenarioEnsemble:
    retention: np.ndarray          # [P, S, n], minutes, NaN = refused
    width: np.ndarray              # [P, S, n], peak SD in minutes
    f_neutral: np.ndarray          # [P, S, n]
    compounds: Sequence[str]       # [n]
    method_cards: Sequence[MethodCard]   # [P]
    stratum: Sequence[str]         # [n] - method-independent (#42)
    ess: dict | None = field(default=None)   # set once reconciliation has run
    # [n] applicability decisions from model/applicability.py, in compound order.
    # Optional only so existing producers keep validating; once supplied, the
    # ensemble can say WHY a compound carries no distribution instead of leaving
    # a reader to infer it from NaN. Typed loosely (the spec layer does not
    # import the model layer) but checked structurally below.
    domain: Sequence[object] | None = field(default=None)

    # ---------------------------------------------------------------- checks
    def validate(self) -> "ScenarioEnsemble":
        r, w, f = self.retention, self.width, self.f_neutral
        if r.ndim != 3:
            raise SpecViolation(
                f"retention must be [P, S, n]; got {r.ndim} dims {r.shape}. The "
                "operating-point axis is not optional - the thin slice needed three "
                "(phi, pH) points and retention[S, n] could not express them.")
        P, S, n = r.shape
        for name, arr in (("width", w), ("f_neutral", f)):
            if arr.shape != r.shape:
                raise SpecViolation(f"{name} must be {r.shape}, got {arr.shape}")
        if len(self.compounds) != n:
            raise SpecViolation(f"compounds must have length n={n}, got {len(self.compounds)}")
        if len(self.method_cards) != P:
            raise SpecViolation(
                f"one method card per operating point: expected {P}, got "
                f"{len(self.method_cards)}")
        if len(self.stratum) != n:
            raise SpecViolation(
                f"stratum must be indexed by n alone (length {n}), got "
                f"{len(self.stratum)}. Indexing it by [P, n] would let a compound "
                "change calibration class because the user changed column, which #42 "
                "forbids: strata are cut on the method-independent structural summary.")

        # refusal semantics: NaN columns are whole-column, and width agrees
        for p in range(P):
            nan_r = np.isnan(r[p])
            col_all = nan_r.all(axis=0)
            col_any = nan_r.any(axis=0)
            if not np.array_equal(col_all, col_any):
                bad = [self.compounds[i] for i in np.where(col_any & ~col_all)[0]]
                raise SpecViolation(
                    f"operating point {p}: refusal is per COMPOUND, not per scenario. "
                    f"Partially-NaN columns: {bad}. A refused compound carries no "
                    "fabricated distribution in any scenario.")
            if not np.array_equal(np.isnan(w[p]).all(axis=0), col_all):
                raise SpecViolation(
                    f"operating point {p}: width must be NaN exactly where retention is")

        live = ~np.isnan(r)
        if np.any(r[live] <= 0):
            raise SpecViolation("retention must be positive where it is not NaN")
        if np.any(w[live] <= 0):
            raise SpecViolation("width must be positive where retention is present")
        fl = f[~np.isnan(f)]
        if fl.size and (fl.min() < 0.0 or fl.max() > 1.0):
            raise SpecViolation(
                f"f_neutral is a fraction; got [{fl.min():.3f}, {fl.max():.3f}]")

        # Domain decisions, when supplied, must agree with the NaN pattern.
        # Refusal is expressed twice on purpose -- as absence of data, which is
        # what a numeric consumer needs, and as a reason, which is what a person
        # needs -- and the two disagreeing is a bug in whoever built the object,
        # not something to reconcile silently at read time.
        if self.domain is not None:
            if len(self.domain) != n:
                raise SpecViolation(
                    f"domain must be indexed by n alone (length {n}), got "
                    f"{len(self.domain)}")
            refused_by_data = set(refused(self))
            for i, d in enumerate(self.domain):
                may = getattr(d, "may_predict", None)
                if may is None:
                    raise SpecViolation(
                        f"domain[{i}] is {type(d).__name__}, which has no "
                        "may_predict; expected an applicability Decision")
                name = self.compounds[i]
                if may and name in refused_by_data:
                    raise SpecViolation(
                        f"{name!r} carries no distribution anywhere but its "
                        f"domain decision permits prediction: a refusal with no "
                        f"recorded reason is the thing this field exists to stop")
                if not may and name not in refused_by_data:
                    raise SpecViolation(
                        f"{name!r} was refused by the domain object "
                        f"({getattr(d, 'tier', '?')}) yet carries a "
                        f"distribution: a refused compound must not be scored, "
                        f"shown, or ranked (#40)")

        # t0 must be able to vary across P - flagged, not fatal, since a genuinely
        # pH-flat set of operating points may legitimately share one value.
        return self

    # ------------------------------------------------------------- refusals
    def refusal_reasons(self) -> dict[str, str]:
        """{compound: rendered reason} for every compound carrying no data.

        Empty when `domain` was not supplied -- and that emptiness is the honest
        answer, not an absence of refusals. `refused()` still reports WHICH
        compounds; only the WHY needs the domain object.
        """
        if self.domain is None:
            return {}
        out = {}
        for i, d in enumerate(self.domain):
            if not getattr(d, "may_predict", True):
                render = getattr(d, "render", None)
                out[self.compounds[i]] = render() if render else str(d)
        return out

    # ------------------------------------------------------------- advisories
    def advisories(self) -> list[str]:
        """Non-fatal conformance notes a consumer should surface."""
        out = []
        P = self.retention.shape[0]
        t0s = {mc.t0_min for mc in self.method_cards}
        phs = {mc.ph for mc in self.method_cards}
        if P > 1 and len(t0s) == 1 and len(phs) > 1:
            out.append(
                f"t0 is identical ({t0s.pop():.4g} min) across {P} operating points "
                f"spanning {len(phs)} different pH values. #25 measured that holding t0 "
                "fixed biases D by +0.42 to +0.69 log units over pH 2.5-6.5. If this is "
                "a nominal t0, t0_source must say NOMINAL and the per-run t0 error term "
                "must be carried.")
        nominal = [i for i, mc in enumerate(self.method_cards) if mc.requires_t0_error_term]
        if nominal:
            out.append(
                f"operating points {nominal} carry a NOMINAL t0: the per-run t0 error "
                "term applies and the prediction is extrapolation across pH (#25).")
        noinst = [i for i, mc in enumerate(self.method_cards)
                  if mc.requires_extra_column_distribution]
        if noinst:
            out.append(
                f"operating points {noinst} have no instrument_id: extra-column variance "
                f"must span the vendor range {EXTRA_COLUMN_UL_RANGE} uL as a distribution, "
                "never zero (sources/instrument-volumes/).")
        if np.isnan(self.f_neutral).all():
            out.append(
                "f_neutral is entirely absent. The sampler computes microspecies weights "
                "and discards them, but D refinement from scouting runs (#34) and "
                "conformal stratum assignment (#15 s.4) both need them and neither can "
                "reconstruct them. This is #44 punch-list item 2.")
        if set(map(str, self.stratum)) <= {"unset", "", "None"}:
            out.append(
                "stratum is unset. Ionisation regime is the one stratum axis that never "
                "merges (#15 s.4), so it must be assigned before calibration can be "
                "computed per stratum.")
        if self.ess is None:
            out.append(
                "ESS is unset. Required as a standing output once reconciliation consumes "
                "a measurement; still unverified end-to-end (#44 item 4, #45 owes it).")
        return out


def refused(ens: ScenarioEnsemble) -> list[str]:
    """Compounds carrying no distribution at any operating point."""
    allnan = np.isnan(ens.retention).all(axis=(0, 1))
    return [c for c, m in zip(ens.compounds, allnan) if m]
