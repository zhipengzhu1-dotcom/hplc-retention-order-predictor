#!/usr/bin/env python3
"""Consolidated sampler — the slice's model with the real ionisation provider.

Three changes from `prototype/thin-slice/slice.py`, all in the ionisation layer:

1. **No hardcoded pKa.** The slice types literature values into `compounds.py`
   and perturbs them by an assumed +-0.35. Here `f_neutral(pH)` comes from the
   provider's microspecies distribution, and the per-scenario uncertainty is the
   provider error actually measured against those same citations (SD 0.170,
   `model/calibrate_ionisation.py`).

2. **No `type` classification.** The slice picks one of three closed-form
   sigmoids per hand-assigned "acid" / "base" / "ampho". The microspecies
   distribution handles monoprotic, multiprotic and zwitterionic by
   construction, so the classification disappears.

3. **`f_neutral` is stored per scenario**, which is #44 punch-list item 2 and a
   MUST in `spec/scenario-ensemble.md`. The slice computes it and throws it away.

**How the uncertainty enters.** The provider emits a point distribution, no
sigma (#7: no provider, open or commercial, emits a calibrated per-prediction
uncertainty). But a pKa error of delta is exactly a shift of the speciation
curve along the pH axis, so querying `f_neutral(pH + delta_s)` with
`delta_s ~ N(scale_shift, 0.170)` reproduces a pKa uncertainty without inventing
a sigma the provider does not have.

⚠ One approximation this makes: for a multiprotic compound a single pH shift
moves every site together, where independent per-site errors would not. That is
arguably the more realistic of the two - a systematic provider bias on a
molecule would affect its sites in common - but it is an assumption, not a
measurement.

Run with env/ionisation-venv/bin/python (needs numpy, scipy, rdkit, unipka).
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SLICE = os.path.join(HERE, "..", "prototype", "thin-slice")
sys.path.insert(0, HERE)
sys.path.insert(0, SLICE)

from ionisation import Provider  # noqa: E402

_CWD = os.getcwd()
os.chdir(SLICE)                      # slice.py resolves data paths relatively
import slice as sl                   # noqa: E402
from compounds import COMPOUNDS      # noqa: E402
os.chdir(_CWD)

# Measured provider error against the cited literature pKa (n=10, 9 compounds).
# The slice assumes 0.35; this is what Uni-pKa actually delivers.
PKA_SD_MEASURED = 0.170
PKA_SD_SLICE = 0.35

# The w_w -> s_s scale correction is UNTOUCHED by this work: both the provider
# and the literature are water-scale, while the mobile phase is not. Carried as
# the slice carries it, a named per-run term (#20).
ACID_SHIFT = (+0.30, 0.30)
BASE_SHIFT = (-0.15, 0.20)

PH_GRID = np.round(np.arange(0.0, 14.0001, 0.02), 4)


def f_neutral_curve(prov: Provider, smiles: str) -> np.ndarray:
    """f_neutral on the standard pH grid, from the provider's microspecies."""
    return np.array([prov.speciate(smiles, float(p)).f_neutral for p in PH_GRID])


def run(use_provider=True, pka_sd=None, seed=sl.SEED, s=sl.S, prov=None,
        d_run_sd=0.30, d_c_sd=0.40):
    """Returns (names, {op: dict(logk, f_neutral)}), plus the refusal set."""
    rng = np.random.default_rng(seed)
    grid, betas, sds = sl.load_column()
    R = sl.design_R()
    stats, _, _, lam, _ = sl.load_ensemble_stats()
    names = [c["name"] for c in COMPOUNDS]
    V = {c["name"]: sl.mcgowan_v(c["smiles"]) for c in COMPOUNDS}
    sd = pka_sd if pka_sd is not None else (
        PKA_SD_MEASURED if use_provider else PKA_SD_SLICE)

    # --- per-run draws (identical to the slice)
    z6 = rng.standard_normal((s, 6))
    dbeta = np.empty((s, len(grid), 6))
    for g in range(len(grid)):
        L = np.linalg.cholesky(np.outer(sds[g], sds[g]) * R + 1e-12 * np.eye(6))
        dbeta[:, g, :] = z6 @ L.T
    acid_shift = rng.normal(*ACID_SHIFT, s)
    base_shift = rng.normal(*BASE_SHIFT, s)
    D_run = rng.normal(0.0, d_run_sd, s) if d_run_sd > 0 else np.zeros(s)

    # --- per-compound
    curves, refused, desc, D, dpk = {}, set(), {}, {}, {}
    for c in COMPOUNDS:
        nm = c["name"]
        if use_provider:
            tier = prov.refusal(c["smiles"])["tier"]
            if tier in ("refuse", "provider_declined"):
                refused.add(nm)
                continue
            curves[nm] = f_neutral_curve(prov, c["smiles"])
        else:
            if c["type"] == "refused":
                refused.add(nm)
                continue
        mu, C = stats[nm]
        d4 = rng.multivariate_normal(mu, np.outer(lam, lam) * C + 1e-10 * np.eye(4),
                                     size=s)
        d4[:, 2] = np.clip(d4[:, 2], 0.0, None)
        desc[nm] = d4
        d_c = rng.normal(0.0, d_c_sd, s) if d_c_sd > 0 else np.zeros(s)
        D[nm] = np.clip(1.5 + D_run + d_c, 0.1, None)
        # per-scenario pH-axis shift == per-scenario pKa error
        base = rng.normal(0.0, sd, s)
        if use_provider:
            # scale shift applies to the ionisable ones; a compound with no
            # ionisable site has a flat curve and is unaffected either way
            sgn = acid_shift if c["type"] in ("acid", "ampho") else (
                base_shift if c["type"] == "base" else np.zeros(s))
            dpk[nm] = base + sgn
        else:
            dpk[nm] = base

    out = {}
    for phi, pH in sl.OPS:
        beta_phi = np.empty((s, 6))
        for j in range(6):
            y = betas[:, j][:, None] + dbeta[:, :, j].T
            beta_phi[:, j] = sl.PchipInterpolator(grid, y, axis=0)(phi)
        logk = np.full((s, len(names)), np.nan)
        fneu = np.full((s, len(names)), np.nan)
        for i, c in enumerate(COMPOUNDS):
            nm = c["name"]
            if nm in refused:
                continue
            d4 = desc[nm]
            Xs = np.column_stack([np.ones(s), d4[:, 0], d4[:, 1], d4[:, 2],
                                  d4[:, 3], np.full(s, V[nm])])
            logk_n = np.einsum("sj,sj->s", Xs, beta_phi)
            if use_provider:
                f = np.interp(pH + dpk[nm], PH_GRID, curves[nm])
            else:
                pk = np.array(c["pka"], float)
                if pk.size == 0:
                    f = np.ones(s)
                else:
                    pk = pk[None, :] + dpk[nm][:, None]
                    f = sl.f_neutral(c["type"], pk, pH)
            fneu[:, i] = f
            logk[:, i] = np.log10(10.0 ** logk_n * (f + (1 - f) * 10.0 ** (-D[nm])))
        out[(phi, pH)] = dict(logk=logk, f_neutral=fneu)
    return names, out, refused


def confident(names, logk, thresh=0.9):
    live = [i for i in range(len(names)) if not np.isnan(logk[:, i]).all()]
    n = c = 0
    for i, j in itertools.combinations(live, 2):
        p = float(np.mean((logk[:, i] - logk[:, j]) < 0))
        n += 1
        c += max(p, 1 - p) > thresh
    return c, n


def main():
    prov = Provider()
    print("Consolidated sampler: real ionisation provider vs the slice's hardcoded pKa\n")

    runs = {
        "slice (hardcoded pKa, sd 0.35)": dict(use_provider=False),
        "provider, slice's sd 0.35":      dict(use_provider=True, pka_sd=PKA_SD_SLICE),
        "provider, measured sd 0.170":    dict(use_provider=True),
    }
    res = {}
    for label, kw in runs.items():
        res[label] = run(prov=prov, **kw)

    print(f"{'configuration':<34}" + "".join(f"{str(tuple(o)):>16}" for o in sl.OPS))
    for label in runs:
        names, out, refused = res[label]
        cells = []
        for o in sl.OPS:
            c, n = confident(names, out[tuple(o)]["logk"])
            cells.append(f"{c}/{n} ({100*c/n:.0f}%)")
        print(f"{label:<34}" + "".join(f"{x:>16}" for x in cells))

    names, out, refused = res["provider, measured sd 0.170"]
    print(f"\nrefused (derived from microspecies): {sorted(refused)}")

    print("\nf_neutral now stored per scenario (#44 item 2) — spread at each op:")
    print(f"  {'compound':<22}" + "".join(f"{str(tuple(o)):>18}" for o in sl.OPS))
    for i, nm in enumerate(names):
        if nm in refused:
            continue
        cells = []
        for o in sl.OPS:
            f = out[tuple(o)]["f_neutral"][:, i]
            cells.append(f"{np.mean(f):.3f} ± {np.std(f):.3f}")
        if any(float(c.split(" ± ")[1]) > 0.005 for c in cells):
            print(f"  {nm:<22}" + "".join(f"{c:>18}" for c in cells))
    print("  (only compounds whose neutral fraction actually varies are listed)")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Assembling the spec object (spec/scenario-ensemble.md)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(HERE, "..", "spec"))
from ensemble import MethodCard, ScenarioEnsemble  # noqa: E402

# Plate count now comes from the Knox model (model/dispersion.py), not a flat
# number: N depends on particle size, flow and the solute's diffusivity, so it is
# a property of the METHOD. A, B, C are drawn per scenario from their priors, so
# the width slot carries the prior's uncertainty rather than its point estimate.
from dispersion import EVO, plate_count  # noqa: E402
COLUMN_GEOM = EVO          # the slice's column is XBridge Shield RP18; EVO's
                           # core-shell geometry stands in until a geometry
                           # registry exists (column_id -> geometry).

# With no instrument_id the spec requires extra-column variance to span the
# vendor range as a DISTRIBUTION, never zero (sources/instrument-volumes/).
EXTRA_COLUMN_UL = (5.0, 34.0)
FLOW_ML_MIN = 1.0
T0_MIN = 1.0                      # nominal; #25 says this should vary with pH


def structural_stratum(prov, smiles):
    """Ionisation regime from the charge states the molecule can adopt.

    Structural and pH-independent, so it is indexed by n alone - which is what
    #42 requires: a compound must not change calibration class because the user
    changed column.
    """
    charges = set()
    for ph in (0.0, 3.5, 7.0, 10.5, 14.0):
        sp = prov.speciate(smiles, ph)
        charges |= {s.charge for s in sp.species if s.fraction > 1e-4}
    if 0 not in charges:
        return "refused"
    if any(c < 0 for c in charges) and any(c > 0 for c in charges):
        return "ampho"
    if any(c < 0 for c in charges):
        return "acid"
    if any(c > 0 for c in charges):
        return "base"
    return "neutral"


def to_ensemble(prov, seed=sl.SEED, s=sl.S):
    """Run the sampler and assemble a spec-conforming ScenarioEnsemble."""
    names, out, refused = run(prov=prov, seed=seed, s=s)
    rng = np.random.default_rng(seed + 1)
    P, n = len(sl.OPS), len(names)

    retention = np.full((P, s, n), np.nan)
    width = np.full((P, s, n), np.nan)
    f_neu = np.full((P, s, n), np.nan)

    # per-scenario instrument draws: per-RUN terms, one draw shared by every
    # compound and every operating point in a scenario
    n_plates = np.clip(plate_count(COLUMN_GEOM, FLOW_ML_MIN, rng, s), 500.0, None)
    ec_ul = rng.uniform(*EXTRA_COLUMN_UL, s)
    sd_ec_min = (ec_ul / 1000.0) / FLOW_ML_MIN          # uL -> mL -> minutes

    for p, op in enumerate(sl.OPS):
        logk = out[tuple(op)]["logk"]
        tR = T0_MIN * (1.0 + 10.0 ** logk)
        # width^2 = column variance + extra-column variance (+ gradient
        # compression, which is null here: every operating point is isocratic)
        var_col = (tR / np.sqrt(n_plates)[:, None]) ** 2
        w = np.sqrt(var_col + (sd_ec_min ** 2)[:, None])
        live = ~np.isnan(tR)
        retention[p][live] = tR[live]
        width[p][live] = w[live]
        f_neu[p] = out[tuple(op)]["f_neutral"]

    cards = [MethodCard(
        phi=float(op[0]), modifier="ACN", ph=float(op[1]),
        ph_scale="WW",                 # #20's correction is still uncorrected
        temperature_c=45.0, column_id="XBridge Shield RP18",
        t0_min=T0_MIN, t0_source="NOMINAL", flow_ml_min=FLOW_ML_MIN,
        instrument_id=None, gradient=None) for op in sl.OPS]

    strat = [structural_stratum(prov, c["smiles"]) for c in COMPOUNDS]
    return ScenarioEnsemble(retention=retention, width=width, f_neutral=f_neu,
                            compounds=names, method_cards=cards,
                            stratum=strat, ess=None)


def emit():
    prov = Provider()
    ens = to_ensemble(prov).validate()
    print("ScenarioEnsemble assembled and VALIDATED against spec/\n")
    print(f"  retention {ens.retention.shape}  width {ens.width.shape}  "
          f"f_neutral {ens.f_neutral.shape}")
    from ensemble import refused as refused_of
    print(f"  refused: {refused_of(ens)}")
    from collections import Counter
    print(f"  strata (structural, indexed by n): {dict(Counter(ens.stratum))}")
    live = ~np.isnan(ens.width)
    print(f"  width: median {np.median(ens.width[live]):.4f} min, "
          f"range {ens.width[live].min():.4f}-{ens.width[live].max():.4f}")
    print("\n  advisories:")
    for a in ens.advisories():
        print(f"   * {a}")
    return ens


if __name__ == "__main__" and os.environ.get("EMIT"):
    emit()
