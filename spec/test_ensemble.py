"""Tests for the scenario-ensemble spec. Run: python3 test_ensemble.py

Asserts the MUSTs in scenario-ensemble.md, then runs the validator against the
REAL thin-slice artefact so the punch list is demonstrated rather than asserted.
No pytest in this environment; runner is at the bottom.
"""
import os
import sys

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ensemble import (MethodCard, ScenarioEnsemble, SpecViolation,  # noqa: E402
                      refused)

FAILURES = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILURES.append(name)


def card(**kw):
    d = dict(phi=0.30, modifier="ACN", ph=3.0, ph_scale="SS", temperature_c=45.0,
             column_id="Kinetex EVO C18", t0_min=1.0, t0_source="MEASURED",
             flow_ml_min=1.0)
    d.update(kw)
    return MethodCard(**d)


def ens(P=2, S=50, n=4, **kw):
    rng = np.random.default_rng(0)
    shape = (P, S, n)
    d = dict(retention=rng.uniform(1.5, 10.0, shape),
             width=rng.uniform(0.02, 0.1, shape),
             f_neutral=rng.uniform(0, 1, shape),
             compounds=[f"c{i}" for i in range(n)],
             method_cards=[card(ph=3.0 + p) for p in range(P)],
             stratum=["neutral"] * n)
    d.update(kw)
    return ScenarioEnsemble(**d)


def test_method_card():
    print("\n[1] Method card")
    check("valid card constructs", card() is not None)
    for bad, why in ((dict(ph_scale="ss"), "ph_scale case-sensitive, no default"),
                     (dict(modifier="acetonitrile"), "modifier is an enum"),
                     (dict(t0_source="guess"), "t0_source is an enum"),
                     (dict(phi=30), "phi is a fraction not a percentage"),
                     (dict(t0_min=0.0), "t0 = 0 is a missing marker, not a value"),
                     (dict(flow_ml_min=0.0), "flow has no default and must be positive")):
        try:
            card(**bad)
            check(why, False, "constructed anyway")
        except SpecViolation:
            check(why, True)
    check("NOMINAL t0 demands the error term",
          card(t0_source="NOMINAL").requires_t0_error_term)
    check("MEASURED t0 does not", not card(t0_source="MEASURED").requires_t0_error_term)
    check("absent instrument demands a distribution",
          card().requires_extra_column_distribution)


def test_shapes():
    print("\n[2] Shape contract")
    check("conforming ensemble validates", ens().validate() is not None)
    try:
        r = np.ones((50, 4))
        ScenarioEnsemble(retention=r, width=r, f_neutral=r, compounds=["a"] * 4,
                         method_cards=[card()], stratum=["x"] * 4).validate()
        check("retention[S, n] rejected (P axis mandatory)", False)
    except SpecViolation as e:
        check("retention[S, n] rejected (P axis mandatory)", "operating-point axis" in str(e))
    for kw, why in (
        (dict(width=np.ones((2, 50, 3))), "width shape must match retention"),
        (dict(method_cards=[card()]), "one method card per operating point"),
        (dict(stratum=["a", "b"]), "stratum length must be n"),
        (dict(compounds=["a"]), "compounds length must be n"),
    ):
        try:
            ens(**kw).validate()
            check(why, False, "accepted")
        except SpecViolation:
            check(why, True)

    # stratum indexed by [P, n] is the specific mistake #42 forbids
    try:
        ens(stratum=["neutral"] * 8).validate()   # P*n instead of n
        check("stratum indexed by [P, n] rejected", False, "accepted")
    except SpecViolation as e:
        check("stratum indexed by [P, n] rejected", "method-independent" in str(e))


def test_refusal():
    print("\n[3] Refusal semantics (#33)")
    e = ens()
    e.retention[:, :, 2] = np.nan
    e.width[:, :, 2] = np.nan
    e.f_neutral[:, :, 2] = np.nan
    e.validate()
    check("whole-column NaN accepted", True)
    check("refused compound identified", refused(e) == ["c2"], str(refused(e)))
    check("n did not shrink", e.retention.shape[2] == 4)

    e2 = ens()
    e2.retention[0, :10, 1] = np.nan          # only some scenarios
    try:
        e2.validate()
        check("partial-column NaN rejected", False, "accepted")
    except SpecViolation as ex:
        check("partial-column NaN rejected", "per COMPOUND" in str(ex))

    e3 = ens()
    e3.retention[:, :, 0] = np.nan            # width not NaN to match
    try:
        e3.validate()
        check("width must be NaN where retention is", False, "accepted")
    except SpecViolation:
        check("width must be NaN where retention is", True)


def test_advisories():
    print("\n[4] Advisories")
    e = ens(method_cards=[card(ph=3.0, t0_source="NOMINAL"),
                          card(ph=7.0, t0_source="NOMINAL")])
    a = " ".join(e.advisories())
    check("constant t0 across differing pH is flagged", "identical" in a)
    check("the flag cites #25's measured bias", "+0.42" in a)
    check("NOMINAL t0 flagged", "NOMINAL t0" in a)
    check("missing instrument flagged", "extra-column variance" in a)
    check("unset ESS flagged", "ESS is unset" in a)

    e2 = ens(f_neutral=np.full((2, 50, 4), np.nan))
    check("absent f_neutral flagged", "f_neutral is entirely absent" in " ".join(e2.advisories()))
    e3 = ens(stratum=["unset"] * 4)
    check("unset stratum flagged", "stratum is unset" in " ".join(e3.advisories()))

    good = ens(method_cards=[card(ph=3.0, t0_min=1.0, instrument_id="acquity-binary"),
                             card(ph=7.0, t0_min=0.93, instrument_id="acquity-binary")],
               ess={"per_request": 812})
    check("fully-specified ensemble has no advisories", good.advisories() == [],
          str(good.advisories()))


def test_against_real_slice():
    print("\n[5] The real thin-slice artefact (#21)")
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "prototype", "thin-slice", "retention_ensemble.npz")
    if not os.path.exists(p):
        print("  SKIP  artefact not present")
        return
    d = np.load(p, allow_pickle=True)
    cards = [card(phi=float(o[0]), ph=float(o[1]), ph_scale="WW",
                  column_id="XBridge Shield RP18", t0_min=float(d["t0"]),
                  t0_source="NOMINAL") for o in d["ops"]]
    e = ScenarioEnsemble(retention=d["retention"], width=d["width"],
                         f_neutral=np.full(d["retention"].shape, np.nan),
                         compounds=list(map(str, d["compounds"])),
                         method_cards=cards,
                         stratum=["unset"] * d["retention"].shape[2])
    e.validate()
    check("slice STRUCTURE conforms (P axis, refusal, shapes)", True)
    check("slice's refused compound recovered",
          refused(e) == ["Benzyltrimethylammonium"], str(refused(e)))
    adv = e.advisories()
    print(f"  -- {len(adv)} advisories, i.e. the punch list, measured not asserted:")
    for a in adv:
        print(f"     * {a[:96]}...")
    check("punch list items 2, 3, 4 all detected", len(adv) >= 4, f"{len(adv)} advisories")



def test_domain_reasons():
    """#40/#33: a refusal must carry its reason, and the two ways of saying
    'refused' must agree. Uses the real Decision objects from model/."""
    print("\n[6] Refusal reasons travel with the object")
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from model.applicability import assess

    ok = assess(elements=["C", "H", "O"], formal_charge=0,
                descriptors={"E": 0.8, "S": 1.0, "A": 0.2, "B": 0.5, "V": 1.0})
    charged = assess(elements=["C", "H", "N"], formal_charge=1)
    check("the domain object refuses the charged compound", not charged.may_predict)

    def nan_out(e, i):
        """Arrays stay mutable on a frozen dataclass; the fields do not."""
        e.retention[:, :, i] = np.nan
        e.width[:, :, i] = np.nan
        e.f_neutral[:, :, i] = np.nan
        return e

    e = nan_out(ens(n=4, domain=[ok, ok, charged, ok]), 2)
    e.validate()
    check("agreeing domain and NaN pattern validate", True)

    reasons = e.refusal_reasons()
    check("the refused compound reports a reason", list(reasons) == ["c2"],
          str(list(reasons)))
    check("the reason names the mechanism, not just the fact",
          "charged_species" in reasons["c2"], reasons["c2"][:70])
    check("the reason carries its warrant", "declared" in reasons["c2"])

    # no data, but a permissive decision: a refusal nobody recorded
    bad = nan_out(ens(n=4, domain=[ok, ok, ok, ok]), 1)
    try:
        bad.validate()
        check("silent refusal rejected", False, "accepted")
    except SpecViolation as ex:
        check("silent refusal rejected", "no recorded reason" in str(ex))

    # refused by the domain object, yet carrying a distribution
    worse = ens(n=4, domain=[ok, charged, ok, ok])
    try:
        worse.validate()
        check("a refused compound carrying data is rejected", False, "accepted")
    except SpecViolation as ex:
        check("a refused compound carrying data is rejected",
              "must not be scored" in str(ex))

    check("no domain supplied means no reasons, not no refusals",
          ens(n=4).refusal_reasons() == {} )



if __name__ == "__main__":
    for t in (test_method_card, test_shapes, test_refusal, test_advisories,
              test_domain_reasons,
              test_against_real_slice):
        t()
    print("\n" + "=" * 64)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        raise SystemExit(1)
    print("All checks passed.")
