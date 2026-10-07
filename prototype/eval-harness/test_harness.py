"""Tests for the #45 evaluation harness. Run: python3 test_harness.py

Two jobs. First, assert the guards CONTEXT.md requires be asserted rather than
left to reviewer discipline. Second, feed the harness synthetic data with a
KNOWN ground truth and check it recovers it - so the pre-registered thresholds
are exercised before any real measurement can influence them.

No pytest in this environment, so the runner is at the bottom.
"""
import numpy as np

import harness as H

FAILURES = []


def check(name, cond, detail=""):
    if cond:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def figure(value, **kw):
    d = dict(axis="pH", baseline="cold (Crippen logP)", baseline_value=0.5,
             n_compounds=15, interval=(value - 0.05, value + 0.05),
             reference_identity="measured retention", reference_kind="measurement")
    d.update(kw)
    return H.Figure(value=value, **d)


# ---------------------------------------------------------------- the guard
def test_figure_guard():
    print("\n[1] Figure guard - no bare accuracy number can escape")
    f = figure(0.72)

    try:
        float(f)
        check("float(Figure) raises", False, "it did not raise")
    except H.BareNumberError:
        check("float(Figure) raises BareNumberError", True)

    try:
        s = f"{f:.2f}"
        check("format spec raises", False, f"produced {s!r}")
    except H.BareNumberError:
        check("format spec raises BareNumberError", True)

    # str() must yield the FULL claim, so even a careless print is compliant
    s = str(f)
    for token in ("axis=", "baseline=", "n=15 compounds", "accurate to"):
        check(f"str(Figure) carries {token!r}", token in s, s)

    # all six fields mandatory
    for bad in (dict(axis=""), dict(baseline=""), dict(reference_identity=""),
                dict(n_compounds=0)):
        try:
            figure(0.7, **bad)
            check(f"missing {list(bad)[0]} rejected", False, "constructed anyway")
        except ValueError:
            check(f"missing {list(bad)[0]} rejected", True)

    try:
        figure(0.7, reference_kind="vibes")
        check("bad reference_kind rejected", False)
    except ValueError:
        check("bad reference_kind rejected", True)


def test_claim_verb():
    print("\n[2] Claim verb follows the reference type (#31)")
    meas = figure(0.72, reference_kind="measurement",
                  reference_identity="measured retention")
    modl = figure(0.72, reference_kind="model_output",
                  reference_identity="ACD Absolv descriptors")
    check("measurement -> 'accurate to'", meas.verb == "accurate to", meas.verb)
    check("model output -> 'agrees with'", modl.verb == "agrees with", modl.verb)
    check("model-output figure never says 'accurate'",
          "accurate" not in modl.render(), modl.render())
    check("reference named in the same sentence",
          "ACD Absolv" in modl.render())


# ------------------------------------------------- synthetic ground truth
def synthetic_pairs(n_compounds=15, true_acc=0.75, seed=1):
    """Pairwise claims with a KNOWN accuracy, correlated within compound.

    Errors attach to COMPOUNDS, not pairs - that is the whole reason #15
    requires compound-level resampling, so the fixture must reproduce it.
    """
    rng = np.random.default_rng(seed)
    bad = rng.random(n_compounds) < 0.30          # some compounds are bad
    idx, truth, p = [], [], []
    for i in range(n_compounds):
        for j in range(i + 1, n_compounds):
            corrupt = bad[i] or bad[j]
            ok = rng.random() < (0.55 if corrupt else 0.92)
            conf = rng.uniform(0.5, 1.0)
            idx.append((i, j))
            truth.append(True)
            p.append(conf if ok else 1 - conf)
    return np.array(idx), np.array(truth), np.array(p)


def test_bootstrap_beats_binomial():
    print("\n[3] Compound-level bootstrap is wider than naive binomial (#15 s.1)")
    n = 15
    idx, truth, p = synthetic_pairs(n)

    def stat(sel):
        keep = np.array([(a in sel) and (b in sel) for a, b in idx])
        if keep.sum() == 0:
            return np.nan
        return H.pairwise_order_accuracy(p[keep], truth[keep])

    point, lo, hi = H.bootstrap_over_compounds(n, stat, n_boot=800)
    boot_sd = (hi - lo) / (2 * 1.96)
    acc = H.pairwise_order_accuracy(p, truth)
    binom_sd = np.sqrt(acc * (1 - acc) / len(p))
    check("bootstrap interval is wider than binomial",
          boot_sd > binom_sd, f"boot {boot_sd:.4f} vs binom {binom_sd:.4f}")
    print(f"        understated by {boot_sd / binom_sd:.2f}x "
          f"(pre-registration says ~1.35x at 15 compounds)")
    check("point estimate inside its own interval", lo <= point <= hi)


def test_lambda_recovers_known_underdispersion():
    print("\n[4] lambda recovers a KNOWN under-dispersion factor")
    rng = np.random.default_rng(7)
    n = 4000
    for true_lambda in (1.0, 1.5, 2.4):
        hw = np.full(n, 0.10)                      # model's claimed half-width
        sigma = true_lambda * 0.10 / 1.6449        # so that lambda*hw = 90% quantile
        resid = rng.normal(0, sigma, n)
        lam = H.conformal_lambda(resid, hw, nominal=0.90)
        check(f"lambda ~ {true_lambda} recovered (got {lam:.2f})",
              abs(lam - true_lambda) < 0.12, f"got {lam:.3f}")

    # bands from the pre-registration
    check("lambda 1.8 -> calibrated",
          H.lambda_verdict(1.8) == "calibrated")
    check("lambda 2.3 -> degraded",
          "degraded" in H.lambda_verdict(2.3))
    check("lambda 2.6 -> rescope",
          "rescopes" in H.lambda_verdict(2.6))
    check("boundary 2.0 is calibrated (<=, not <)",
          H.lambda_verdict(2.0) == "calibrated")
    check("boundary 2.5 is degraded, not rescope",
          "degraded" in H.lambda_verdict(2.5))


def test_coverage_is_not_the_test():
    print("\n[5] Coverage passes even when sharpness fails (#15 s.3)")
    rng = np.random.default_rng(3)
    resid = rng.normal(0, 0.25, 4000)
    fat = np.full(4000, 1.0)                   # absurdly wide intervals
    cov = H.interval_coverage(resid, fat)
    lam = H.conformal_lambda(resid, fat)
    check("coverage looks perfect on absurd widths", cov > 0.99, f"{cov:.3f}")
    check("lambda correctly reports they were unnecessary", lam < 0.6, f"{lam:.3f}")
    print("        -> a coverage threshold is a test we cannot fail; lambda can.")


def test_three_verdicts():
    print("\n[6] Three verdicts, and inconclusive states the required n (#15 s.2)")
    v = H.decide(gain=0.10, margin_2sd=0.04, n_compounds=15, axis="pH")
    check("clear gain -> pass", v.outcome == "pass", v.outcome)
    v = H.decide(gain=-0.10, margin_2sd=0.04, n_compounds=15, axis="pH")
    check("clear loss -> fail", v.outcome == "fail", v.outcome)

    v = H.decide(gain=0.02, margin_2sd=0.05, n_compounds=15, axis="pH")
    check("small gain -> inconclusive", v.outcome == "inconclusive", v.outcome)
    check("inconclusive is not a pass", v.outcome != "pass")
    check("required n is stated", v.required_n is not None and v.required_n > 15,
          str(v.required_n))
    check("required n appears in the rendered text",
          str(v.required_n) in v.render(), v.render())
    print(f"        {v.render()}")


def test_relative_error_reduction():
    print("\n[7] Relative error reduction, not percentage points (#27/#15 s.2)")
    model = np.full(100, 0.083)
    base = np.full(100, 0.100)
    r = H.relative_error_reduction(model, base)
    check("17% reduction recovered", abs(r - 0.17) < 1e-9, f"{r:.4f}")
    try:
        H.relative_error_reduction(model, np.zeros(100))
        check("zero baseline rejected", False)
    except ZeroDivisionError:
        check("zero baseline rejected", True)


def test_strata_merging():
    print("\n[8] Strata merge upward; ionisation regime never merges (#15 s.4)")
    keys = ("ionisation_regime", "confidence_tertile", "phi_range")
    strata = {
        ("neutral", "hi", "30"): 80,
        ("neutral", "mid", "30"): 12,
        ("neutral", "lo", "30"): 9,
        ("ionised", "hi", "30"): 15,
        ("ionised", "lo", "30"): 11,
    }
    assign = H.merge_strata(strata, keys, floor=50)
    for orig, tgt in assign.items():
        check(f"{orig[0]} regime preserved through merge", tgt[0] == orig[0],
              f"{orig} -> {tgt}")
    merged_neutral = {assign[k] for k in strata if k[0] == "neutral" and strata[k] < 50}
    check("undersized neutral strata collapsed on confidence tertile",
          all(t[1] == "*" for t in merged_neutral), str(merged_neutral))
    check("healthy stratum untouched",
          assign[("neutral", "hi", "30")] == ("neutral", "hi", "30"))
    check("hard floor: n=18 not calibratable", not H.calibratable(18))
    check("hard floor: n=19 calibratable", H.calibratable(19))


def test_end_to_end():
    print("\n[9] End-to-end: a reported figure carries everything")
    n = 15
    idx, truth, p = synthetic_pairs(n)

    def stat(sel):
        keep = np.array([(a in sel) and (b in sel) for a, b in idx])
        return (H.pairwise_order_accuracy(p[keep], truth[keep])
                if keep.sum() else np.nan)

    point, lo, hi = H.bootstrap_over_compounds(n, stat, n_boot=600)
    fig = H.Figure(value=point, axis="pH", baseline="cold (Crippen logP)",
                   baseline_value=0.765, n_compounds=n, interval=(lo, hi),
                   reference_identity="synthetic ground truth",
                   reference_kind="measurement")
    print(f"        {fig.render()}")
    for token in ("[", "n=15 compounds", "axis=pH", "baseline="):
        check(f"rendered figure carries {token!r}", token in fig.render())
    ece = H.expected_calibration_error(p, truth)
    check("ECE is finite and in [0,1]", 0.0 <= ece <= 1.0, f"{ece:.4f}")
    rel = H.reliability_table(p, truth)
    check("reliability table populated", len(rel) > 3, str(len(rel)))



def test_domain_firewall():
    """#40 at the point it bites: refused rows never score, degraded rows
    score apart, and a pooled figure over mixed strata raises."""
    print("\n-- the applicability firewall (#40) --")

    # a pair is only as clean as its worse compound
    check("in+in pair is in_envelope",
          H.pair_stratum("in_envelope", "in_envelope") == "in_envelope")
    check("in+degraded pair is degraded",
          H.pair_stratum("in_envelope", "degraded") == "degraded")
    check("any pair touching a refused compound is unscorable",
          H.pair_stratum("in_envelope", None) is None
          and H.pair_stratum(None, "degraded") is None)

    # the contamination arithmetic that motivates the pair rule
    n = 15
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    compound_strata = ["in_envelope"] * n
    compound_strata[7] = "degraded"          # one degraded compound
    ps = H.strata_for_pairs(compound_strata, pairs)
    part = H.partition_by_stratum(ps)
    check("one degraded compound degrades 14 of 105 pairs",
          len(part["degraded"]) == 14 and len(part["in_envelope"]) == 91,
          f"{len(part['degraded'])}/{len(pairs)}")

    # refused rows are dropped, never grouped
    compound_strata[3] = None
    ps = H.strata_for_pairs(compound_strata, pairs)
    part = H.partition_by_stratum(ps)
    check("refused compound removes its 14 pairs entirely",
          sum(len(v) for v in part.values()) == len(pairs) - 14,
          str(sum(len(v) for v in part.values())))

    # the guard
    try:
        H.require_single_stratum(["in_envelope", "degraded"])
        check("mixed strata raise", False, "no exception")
    except H.MixedStrataError:
        check("mixed strata raise", True)
    try:
        H.require_single_stratum([None, None])
        check("all-refused raises rather than returning an empty figure", False)
    except H.MixedStrataError:
        check("all-refused raises rather than returning an empty figure", True)
    check("single stratum passes",
          H.require_single_stratum(["degraded", "degraded", None]) == "degraded")

    # the entry point: one Figure per stratum, each labelled
    def make(rows, stratum):
        return H.Figure(value=0.9, axis="pH", baseline="cold (Crippen logP)",
                        baseline_value=0.7, n_compounds=len(rows),
                        interval=(0.8, 0.95),
                        reference_identity="synthetic ground truth",
                        reference_kind="measurement", stratum=stratum)

    figs = H.score_by_stratum(ps, make)
    check("a figure is produced per stratum, never pooled",
          set(figs) == {"in_envelope", "degraded"}, str(sorted(figs)))
    check("the stratum is rendered next to n",
          "stratum=degraded" in figs["degraded"].render(),
          figs["degraded"].render()[:90])

    def mislabel(rows, stratum):
        return H.Figure(value=0.9, axis="pH", baseline="b", baseline_value=0.7,
                        n_compounds=len(rows), interval=(0.8, 0.95),
                        reference_identity="x", reference_kind="measurement",
                        stratum="in_envelope")
    try:
        H.score_by_stratum(["degraded"], mislabel)
        check("a mislabelled stratum raises", False, "no exception")
    except H.MixedStrataError:
        check("a mislabelled stratum raises", True)

    # a figure with no stratum still renders (ungated legacy call sites)
    plain = H.Figure(value=0.5, axis="pH", baseline="b", baseline_value=0.4,
                     n_compounds=3, interval=(0.4, 0.6),
                     reference_identity="x", reference_kind="measurement")
    check("an ungated figure renders without a stratum",
          "stratum=" not in plain.render())



if __name__ == "__main__":
    for t in (test_figure_guard, test_claim_verb, test_bootstrap_beats_binomial,
              test_lambda_recovers_known_underdispersion,
              test_coverage_is_not_the_test, test_three_verdicts,
              test_relative_error_reduction, test_strata_merging,
              test_domain_firewall, test_end_to_end):
        t()
    print("\n" + "=" * 64)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        raise SystemExit(1)
    print("All checks passed.")
