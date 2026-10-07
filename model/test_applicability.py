"""Checks for the applicability domain. `python3 model/test_applicability.py`.

Stdlib only. The checks that matter are the ones asserting a coverage statement
and a mechanism assumption stay distinguishable, since collapsing them is the
failure #43 was reframed to prevent.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from model.applicability import (  # noqa: E402
    ANCHOR_MAX, ANCHOR_P95, Decision, Reason, assess, evaluation_strata,
)

CHECKS = []


def check(fn):
    CHECKS.append(fn)
    return fn


NEUTRAL = dict(elements=["C", "H", "O"], formal_charge=0)
TYPICAL = {"E": 0.8, "S": 1.0, "A": 0.2, "B": 0.5, "V": 1.0}
# A context the source actually flags (Ascentis C18, ref 29, both modifiers),
# and one it does not. 8 of the 25 registered columns carry a steric flag.
FLAGGED = dict(column="Ascentis C18", modifier="methanol", phi_pct=20.0)
UNFLAGGED = dict(column="Kinetex EVO C18", modifier="methanol", phi_pct=20.0)


# ---------------------------------------------------------------- the warrants
@check
def reason_requires_a_warrant():
    try:
        Reason(code="x", detail="y", warrant="obvious", tier="refused")
    except ValueError:
        return
    raise AssertionError("a free-text warrant was accepted")


@check
def coverage_and_mechanism_never_merge():
    """The #43 reframing, asserted: a coverage refusal and a steric refusal
    must remain separable after both have been reached."""
    d = assess(**NEUTRAL, descriptors={"V": 3.0}, steric_class="bulky",
               steric_assumption_enabled=True, **FLAGGED)
    codes = {r.code for r in d.reasons}
    assert "outside_anchor_range" in codes, codes
    assert "steric_exclusion_assumed" in codes, codes
    assert {r.code for r in d.assumptions} == {"steric_exclusion_assumed"}, \
        "the coverage reason leaked into the assumptions"
    assert d.rests_on_assumption


@check
def coverage_alone_rests_on_no_assumption():
    d = assess(**NEUTRAL, descriptors={"V": 3.0})
    assert d.tier == "degraded"
    assert not d.rests_on_assumption
    assert d.assumptions == ()


@check
def the_steric_predicate_is_off_by_default():
    """Ref 29 was never obtained, so the predicate must not act unasked."""
    d = assess(**NEUTRAL, descriptors=TYPICAL, steric_class="bulky")
    assert d.may_predict, "an untested assumption blocked a prediction by default"
    assert d.tier == "in_envelope"
    assert d.rests_on_assumption, "the class was silently forgotten instead of noted"


@check
def enabling_the_steric_predicate_refuses_only_in_scope():
    """#40 blocks in an affected (column, modifier, phi-range), not everywhere."""
    d = assess(**NEUTRAL, descriptors=TYPICAL, steric_class="bulky",
               steric_assumption_enabled=True, **FLAGGED)
    assert d.tier == "refused" and not d.may_predict, d.render()


@check
def render_labels_assumptions_as_not_evidence():
    d = assess(**NEUTRAL, steric_class="bulky", steric_assumption_enabled=True,
               **FLAGGED)
    text = d.render()
    assert "UNVERIFIED ASSUMPTIONS -- not evidence:" in text, text
    assert "assumed" in text


# ------------------------------------------------------------------ the tiers
@check
def envelope_elements_refuse():
    d = assess(elements=["C", "H", "B"], formal_charge=0, descriptors=TYPICAL)
    assert d.tier == "refused"
    assert any(r.code == "element_outside_envelope" for r in d.reasons)
    assert all(r.warrant == "declared" for r in d.reasons
               if r.code == "element_outside_envelope")


@check
def charge_refuses_and_so_does_permanent_charge():
    assert assess(elements=["C"], formal_charge=1).tier == "refused"
    assert assess(**NEUTRAL, permanently_charged=True).tier == "refused"


@check
def refused_ionisation_refuses():
    assert assess(**NEUTRAL, ionisation_refused=True).tier == "refused"


@check
def p95_degrades_and_max_degrades_with_a_different_code():
    tail = assess(**NEUTRAL, descriptors={"B": ANCHOR_P95["B"] + 0.01})
    beyond = assess(**NEUTRAL, descriptors={"B": ANCHOR_MAX["B"] + 0.01})
    assert tail.tier == beyond.tier == "degraded"
    assert {r.code for r in tail.reasons} == {"anchor_upper_tail"}
    assert {r.code for r in beyond.reasons} == {"outside_anchor_range"}


@check
def a_descriptor_is_not_double_counted():
    """Beyond max implies beyond p95; it must report once, at the worse code."""
    d = assess(**NEUTRAL, descriptors={"B": 99.0})
    assert len([r for r in d.reasons if "anchor" in r.code]) == 1, d.render()


@check
def the_thin_slice_drugs_degrade_rather_than_pass():
    """#31: atenolol B 2.11 and ketoprofen S 2.28 are outside the anchor."""
    for name, d in [("atenolol", {"B": 2.11}), ("ketoprofen", {"S": 2.28})]:
        v = assess(**NEUTRAL, descriptors=d)
        assert v.tier == "degraded", f"{name} passed as in-envelope"
        assert v.may_predict, f"{name} was refused, not degraded"


@check
def an_ordinary_compound_passes_and_says_why():
    d = assess(**NEUTRAL, descriptors=TYPICAL)
    assert d.tier == "in_envelope" and d.may_predict
    assert d.reasons[0].warrant == "measured"
    assert "n=94" in d.reasons[0].source


# ---------------------------------------------------------------- the firewall
@check
def refused_rows_never_enter_a_metric():
    ds = [assess(**NEUTRAL, descriptors=TYPICAL),
          assess(**NEUTRAL, descriptors={"V": 3.0}),
          assess(elements=["C", "Si"], formal_charge=0)]
    assert evaluation_strata(ds) == ["in_envelope", "degraded", None]


@check
def degraded_rows_are_scorable_but_in_their_own_stratum():
    """The tier the sweep will actually deliver: every thin-slice drug
    degrades, so excluding degraded rows would leave nothing to score."""
    d = assess(**NEUTRAL, descriptors={"B": 2.11})
    assert d.tier == "degraded"
    assert d.may_enter_evaluation
    assert d.evaluation_stratum == "degraded"
    assert d.evaluation_stratum != "in_envelope", "degraded pooled with clean"


@check
def opt_in_residuals_are_never_evaluation_data():
    ds = [assess(**NEUTRAL, descriptors=TYPICAL)] * 2
    assert evaluation_strata(ds, opt_in=[False, True]) == ["in_envelope", None]


@check
def corpus_tier_is_carried_but_gates_nothing():
    """The census tiers describe evidence strength, not permission."""
    for tier in ("in_corpus", "scaffold_seen", "novel_in_envelope"):
        d = assess(**NEUTRAL, descriptors=TYPICAL, corpus_tier=tier)
        assert d.tier == "in_envelope" and d.evaluation_stratum == "in_envelope"
        assert tier in d.render()



@check
def an_unflagged_column_is_not_refused():
    """The defect this scoping fixes: 17 of 25 columns carry no steric flag, so
    an unscoped predicate would refuse on more than twice the affected set."""
    d = assess(**NEUTRAL, descriptors=TYPICAL, steric_class="bulky",
               steric_assumption_enabled=True, **UNFLAGGED)
    assert d.tier == "in_envelope", d.render()
    assert any(r.code == "steric_out_of_scope" for r in d.reasons)


@check
def the_wrong_modifier_is_not_refused():
    """Discovery HS C18's flag is methanol-only; acetonitrile is out of scope."""
    acn = assess(**NEUTRAL, descriptors=TYPICAL, steric_class="bulky",
                 steric_assumption_enabled=True, column="Discovery HS C18",
                 modifier="acetonitrile", phi_pct=20.0)
    meoh = assess(**NEUTRAL, descriptors=TYPICAL, steric_class="bulky",
                  steric_assumption_enabled=True, column="Discovery HS C18",
                  modifier="methanol", phi_pct=20.0)
    assert acn.tier == "in_envelope", acn.render()
    assert meoh.tier == "refused", meoh.render()


@check
def a_quantified_phi_threshold_is_honoured():
    """XTerra MS C18: 'less than 40% (v/v)'. Above it, the caveat does not apply."""
    below = assess(**NEUTRAL, steric_class="bulky", steric_assumption_enabled=True,
                   column="XTerra MS C18", modifier="methanol", phi_pct=30.0)
    above = assess(**NEUTRAL, steric_class="bulky", steric_assumption_enabled=True,
                   column="XTerra MS C18", modifier="methanol", phi_pct=60.0)
    assert below.tier == "refused", below.render()
    assert above.tier == "in_envelope", above.render()


@check
def an_unquantified_phi_range_takes_the_whole_axis():
    """Most notes say only 'water-rich'. The conservative reading is the whole
    range, and the reason must say the source gave no number rather than
    inventing one."""
    d = assess(**NEUTRAL, steric_class="bulky", steric_assumption_enabled=True,
               column="Ascentis C18", modifier="methanol", phi_pct=90.0)
    assert d.tier == "refused"
    detail = " ".join(r.detail for r in d.reasons)
    assert "unquantified" in detail, detail


@check
def scope_without_a_column_degrades_rather_than_guessing():
    d = assess(**NEUTRAL, steric_class="bulky", steric_assumption_enabled=True)
    assert d.tier == "degraded", d.render()
    assert any(r.code == "steric_scope_unknown" for r in d.reasons)


@check
def the_scope_is_declared_and_the_class_is_assumed():
    """The two halves of #40's rule have different warrants and must keep them."""
    d = assess(**NEUTRAL, steric_class="bulky", steric_assumption_enabled=True,
               **FLAGGED)
    steric = [r for r in d.reasons if r.code == "steric_exclusion_assumed"][0]
    assert steric.warrant == "assumed"
    assert "Table S-1" in steric.source and "assumed" in steric.source
    assert "ref 29" in steric.detail, steric.detail
    out = assess(**NEUTRAL, steric_class="bulky",
                 steric_assumption_enabled=True, **UNFLAGGED)
    assert [r for r in out.reasons
            if r.code == "steric_out_of_scope"][0].warrant == "declared"


@check
def the_register_matches_the_source_table():
    from model.applicability import steric_register
    reg = steric_register()
    assert len(reg) == 8, sorted(reg)
    assert reg["XTerra MS C18"].phi_max_pct == 40.0
    assert reg["Ascentis C18"].phi_max_pct is None
    assert reg["Discovery HS F5"].modifiers == frozenset({"methanol"})
    # Fluophase-RP's note names hydrogen-bond acids, not bulk -- the register
    # is not one structural axis, and the note is preserved so that stays visible
    assert "hydrogen-bond acids" in reg["Fluophase-RP"].note



if __name__ == "__main__":
    failed = skipped = 0
    for fn in CHECKS:
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except FileNotFoundError as e:
            # Checks that read the WSU-2019 caveat register skip when it has
            # not been built locally; the data is not redistributed.
            skipped += 1
            print(f"  skip {fn.__name__}: {os.path.relpath(e.filename)} missing")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    ran = len(CHECKS) - skipped
    print(f"\n{ran - failed}/{ran} checks passed, {skipped} skipped")
    if skipped:
        print("Skipped checks need sources/wsu-lser/wsu2019-column-caveats.csv; "
              "build it with sources/wsu-lser/extract_wsu.py "
              "(see sources/wsu-lser/README.md).")
    sys.exit(1 if failed else 0)
