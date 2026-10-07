#!/usr/bin/env python3
"""End-to-end: model -> ScenarioEnsemble -> evaluation harness -> Figure.

This closes the loop #45 needs. The harness (`prototype/eval-harness/`) was
built before any data existed, deliberately, so that #15's pre-registered
thresholds could not be met by choosing the statistic after seeing the result.
The model now emits a spec-conforming `ScenarioEnsemble`. This wires the two
together so that when measured retention arrives, nothing has to be written
under time pressure.

⚠ **There is still no measured retention data.** To exercise the pipeline this
draws one scenario from the ensemble and treats it as "truth". That is a
PLUMBING test, not a validation — and the harness enforces the distinction: a
reference that is model output can only ever take the verb "agrees with", never
"accurate to" (#31). The rendered Figure below says so itself.

Run with env/ionisation-venv/bin/python.
"""
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "spec"))
sys.path.insert(0, os.path.join(HERE, "..", "prototype", "eval-harness"))

import harness as H          # noqa: E402
from ionisation import Provider  # noqa: E402
from sampler import to_ensemble  # noqa: E402


def pairwise(ens, p, truth_scenario):
    """Predicted P(A before B) and the truth-scenario ordering, per pair."""
    ret = ens.retention[p]                       # (S, n)
    live = [i for i in range(ret.shape[1]) if not np.isnan(ret[:, i]).all()]
    truth = ret[truth_scenario]
    keep = [i for i in live if i != truth_scenario]
    rest = np.delete(ret, truth_scenario, axis=0)
    idx, p_before, t_before = [], [], []
    for i, j in itertools.combinations(live, 2):
        p_before.append(float(np.mean(rest[:, i] < rest[:, j])))
        t_before.append(bool(truth[i] < truth[j]))
        idx.append((i, j))
    return np.array(idx), np.array(p_before), np.array(t_before), live


def main():
    prov = Provider()
    ens = to_ensemble(prov).validate()
    print("Model -> ScenarioEnsemble -> harness, end to end\n")

    for p, card in enumerate(ens.method_cards):
        idx, pb, tb, live = pairwise(ens, p, truth_scenario=0)
        n_comp = len(live)

        def stat(sel):
            sel = set(int(x) for x in sel)
            m = np.array([(a in sel) and (b in sel) for a, b in idx])
            if not m.any():
                return np.nan
            return H.pairwise_order_accuracy(pb[m], tb[m])

        # compound-level bootstrap, never binomial on pairs (#15 s.1)
        point, lo, hi = H.bootstrap_over_compounds(
            n_comp, lambda s: stat([live[k] for k in s]), n_boot=400)

        fig = H.Figure(
            value=point, axis=f"phi={card.phi}, pH={card.ph} ({card.ph_scale})",
            baseline="cold (Crippen logP)", baseline_value=None,
            n_compounds=n_comp, interval=(lo, hi),
            reference_identity="a scenario drawn from this same ensemble",
            reference_kind="model_output")
        print(f"  {fig.render()}")

        ece = H.expected_calibration_error(pb, tb)
        # Restrict to live columns first: a refused compound is an all-NaN
        # column by construction (#33), and nanstd over it warns rather than
        # returning something meaningful.
        col = ens.retention[p][:, live]
        resid = col[0] - np.nanmean(col, axis=0)
        hw = 1.6449 * np.nanstd(col, axis=0)
        ok = ~np.isnan(resid) & (hw > 0)
        lam = H.conformal_lambda(resid[ok], hw[ok])
        print(f"     ECE {ece:.4f} | lambda {lam:.3f} -> {H.lambda_verdict(lam)}\n")

    print("What this does and does not show")
    print("  DOES: the pipeline runs end to end and every number leaves the")
    print("        harness wearing its six mandatory fields.")
    print("  DOES NOT: say anything about accuracy. The reference is the model's")
    print("        own output, so the only permitted verb is 'agrees with' -")
    print("        and the harness enforced that, not the author.")
    print("  The lambda above is likewise a self-consistency check: a scenario")
    print("  drawn from the ensemble is exchangeable with it by construction.")


if __name__ == "__main__":
    main()
