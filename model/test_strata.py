"""Checks for the stratum producer and for #36's corrector fitting order.

`python3 model/test_strata.py` — stdlib only.

The corrector-order checks are here because `CONTEXT.md` claims they exist:
"Corrector fitting order is fixed, and asserted in a test." Until now no such
test existed anywhere in the repo, which made the convention unenforced exactly
where it is least visible — all three correctors can absorb the same discrepancy,
so a wrong order shows up as a plausible number, never as a failure.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from model.strata import (  # noqa: E402
    TERTILE_CUTS, assign, confidence_tertile, ionisation_regime, stratum_key,
    stratum_label, worst_regime_f_neutral,
)
from model.correctors import (  # noqa: E402
    CORRECTOR_ORDER, CorrectorOrderViolation, FitLog,
)

CHECKS = []


def check(fn):
    CHECKS.append(fn)
    return fn


# --------------------------------------------------------------- the strata
@check
def tertiles_are_frozen_not_recomputed():
    """The cuts are the anchor's, pinned. A compound's class must not depend on
    which other compounds happen to be in the batch."""
    assert TERTILE_CUTS == (0.1286, 0.1635)
    assert confidence_tertile(0.05) == "confident"
    assert confidence_tertile(0.30) == "uncertain"


@check
def tertile_boundaries_are_closed_below():
    assert confidence_tertile(TERTILE_CUTS[0]) == "confident"
    assert confidence_tertile(TERTILE_CUTS[0] + 1e-9) == "typical"
    assert confidence_tertile(TERTILE_CUTS[1]) == "typical"
    assert confidence_tertile(TERTILE_CUTS[1] + 1e-9) == "uncertain"


@check
def regimes_come_from_f_neutral():
    assert ionisation_regime(1.0) == "neutral"
    assert ionisation_regime(0.5) == "partial"
    assert ionisation_regime(0.0) == "ionised"


@check
def a_fraction_outside_zero_one_is_refused():
    for bad in (-0.01, 1.01):
        try:
            ionisation_regime(bad)
        except ValueError:
            continue
        raise AssertionError(f"accepted f_neutral={bad}")


@check
def the_key_order_matches_the_merge_rule():
    """harness.merge_strata requires ionisation_regime in the key, and merges
    confidence_tertile first. The key must be ordered to match."""
    sys.path.insert(0, os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "prototype", "eval-harness"))
    import harness as H
    keys = ("ionisation_regime", "confidence_tertile")
    assert H.NEVER_MERGE == keys[0]
    assert H.MERGE_ORDER[0] == keys[1]
    assert len(stratum_key(1.0, 0.1)) == len(keys)


@check
def a_compound_spanning_regimes_takes_its_worst():
    """#42 forbids a compound changing class across operating points."""
    assert ionisation_regime(worst_regime_f_neutral([1.0, 0.5, 0.0])) == "partial"
    assert ionisation_regime(worst_regime_f_neutral([1.0, 0.99])) == "neutral"
    assert ionisation_regime(worst_regime_f_neutral([0.0, 0.02])) == "ionised"


@check
def assign_is_per_compound_and_length_checked():
    labels = assign([1.0, 0.5], [0.05, 0.30])
    assert labels == ["neutral/confident", "partial/uncertain"], labels
    try:
        assign([1.0], [0.1, 0.2])
    except ValueError:
        return
    raise AssertionError("mismatched lengths accepted")


@check
def the_label_is_the_key_joined():
    assert stratum_label(1.0, 0.05) == "neutral/confident"


# ------------------------------------------------- #36's corrector ordering
@check
def the_pinned_order_is_the_context_order():
    assert CORRECTOR_ORDER == ("wsu_refit", "solute_vector", "D"), CORRECTOR_ORDER


@check
def fitting_in_order_is_allowed():
    log = FitLog()
    for stage in CORRECTOR_ORDER:
        log.fit(stage)
    assert log.completed == list(CORRECTOR_ORDER)


@check
def fitting_out_of_order_raises():
    log = FitLog()
    log.fit("wsu_refit")
    try:
        log.fit("D")
    except CorrectorOrderViolation as e:
        assert "solute_vector" in str(e), str(e)
        return
    raise AssertionError("D was fitted before solute_vector")


@check
def refitting_an_earlier_stage_raises():
    """The failure this guards: re-fitting the global stage after a local one
    silently reassigns credit between them."""
    log = FitLog()
    log.fit("wsu_refit")
    log.fit("solute_vector")
    try:
        log.fit("wsu_refit")
    except CorrectorOrderViolation:
        return
    raise AssertionError("an earlier stage was re-fitted")


@check
def joint_fitting_is_refused_by_name():
    log = FitLog()
    try:
        log.fit_jointly(["wsu_refit", "solute_vector"])
    except CorrectorOrderViolation as e:
        assert "jointly" in str(e).lower()
        return
    raise AssertionError("joint fitting was permitted")


@check
def uncertainty_comes_from_the_final_stage_only():
    log = FitLog()
    log.fit("wsu_refit", uncertainty=0.11)
    log.fit("solute_vector", uncertainty=0.07)
    log.fit("D", uncertainty=0.04)
    assert log.reportable_uncertainty() == 0.04


@check
def uncertainty_is_refused_before_the_chain_completes():
    """Reporting the middle stage's uncertainty as if it were the answer is the
    quiet version of the same error."""
    log = FitLog()
    log.fit("wsu_refit", uncertainty=0.11)
    try:
        log.reportable_uncertainty()
    except CorrectorOrderViolation:
        return
    raise AssertionError("uncertainty reported from an incomplete chain")


@check
def an_unknown_stage_is_refused():
    log = FitLog()
    try:
        log.fit("temperature_fudge")
    except CorrectorOrderViolation:
        return
    raise AssertionError("an unpinned corrector was accepted into the chain")


if __name__ == "__main__":
    failed = 0
    for fn in CHECKS:
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    print(f"\n{len(CHECKS) - failed}/{len(CHECKS)} checks passed")
    sys.exit(1 if failed else 0)
