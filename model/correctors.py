"""The corrector chain's fitting order, made enforceable.

`CONTEXT.md` (#36, #34) pins it:

    #35's WSU re-fit (global, lineage) -> solute-vector correction on its
    residual -> `D` last (per-compound, scouting-informed). Never fitted
    jointly. Uncertainty is reported from the final stage only.
    ...asserted in a test.

The convention existed; the assertion did not. This module is the missing half.

Why it needs enforcing rather than documenting: **all three correctors can
absorb the same discrepancy** for an ionisable compound on a characterised
column. Fit them in the wrong order and every number still looks reasonable —
the last one fitted simply claims credit that belongs to an earlier one, and the
uncertainties can no longer be added because they are no longer independent
contributions. There is no output that looks wrong, which is exactly why a
reviewer will not catch it.

`FitLog` is a ledger, not a fitter. It records what was fitted and refuses the
orderings the convention forbids; the numerical work lives elsewhere.
"""

from __future__ import annotations

from typing import Sequence

CORRECTOR_ORDER: tuple[str, ...] = ("wsu_refit", "solute_vector", "D")

_DESCRIPTION = {
    "wsu_refit": "global lineage re-fit (#35)",
    "solute_vector": "per-compound solute-vector correction on the re-fit's residual (#36)",
    "D": "per-compound, scouting-informed D (#34)",
}


class CorrectorOrderViolation(RuntimeError):
    """Raised when the pinned corrector order would be broken."""


class FitLog:
    """Ledger of corrector stages fitted, in order."""

    def __init__(self) -> None:
        self.completed: list[str] = []
        self._uncertainty: dict[str, float] = {}

    # ------------------------------------------------------------------ fit
    def fit(self, stage: str, uncertainty: float | None = None) -> None:
        if stage not in CORRECTOR_ORDER:
            raise CorrectorOrderViolation(
                f"{stage!r} is not a pinned corrector. The chain is "
                f"{' -> '.join(CORRECTOR_ORDER)}; adding a stage changes what the "
                f"reported uncertainty means and needs a ticket, not a keyword.")
        expected = CORRECTOR_ORDER[len(self.completed)] \
            if len(self.completed) < len(CORRECTOR_ORDER) else None
        if stage in self.completed:
            raise CorrectorOrderViolation(
                f"{stage!r} was already fitted. Re-fitting an earlier stage after a "
                f"later one silently reassigns credit between correctors that can "
                f"absorb the same discrepancy.")
        if stage != expected:
            raise CorrectorOrderViolation(
                f"{stage!r} ({_DESCRIPTION[stage]}) cannot be fitted yet: the next "
                f"stage is {expected!r} ({_DESCRIPTION.get(expected, '')}). Pinned "
                f"order is {' -> '.join(CORRECTOR_ORDER)}.")
        self.completed.append(stage)
        if uncertainty is not None:
            self._uncertainty[stage] = float(uncertainty)

    def fit_jointly(self, stages: Sequence[str]) -> None:
        """Exists only to refuse, by name, the thing the convention forbids."""
        raise CorrectorOrderViolation(
            f"correctors {list(stages)} may not be fitted jointly. Without a pinned "
            f"order the last one fitted claims the credit and the uncertainties "
            f"cannot be added (#36).")

    # -------------------------------------------------------- what to report
    def reportable_uncertainty(self) -> float:
        """The final stage's uncertainty — the only one that may be quoted."""
        if self.completed != list(CORRECTOR_ORDER):
            missing = [s for s in CORRECTOR_ORDER if s not in self.completed]
            raise CorrectorOrderViolation(
                f"the chain is incomplete (missing {missing}); an intermediate "
                f"stage's uncertainty is not the model's uncertainty and quoting it "
                f"understates by however much the remaining stages would add.")
        last = CORRECTOR_ORDER[-1]
        if last not in self._uncertainty:
            raise CorrectorOrderViolation(
                f"stage {last!r} was fitted without recording an uncertainty")
        return self._uncertainty[last]
