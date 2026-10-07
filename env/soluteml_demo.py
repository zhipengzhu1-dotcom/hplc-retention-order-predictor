"""Demonstration for #38 — SoluteML ensemble covariance over (E, S, A, B).

Run with the patched environment:

    env/qspr-venv/bin/python env/soluteml_demo.py

What it shows
-------------
1.  Reproduction check against the upstream repo's own published sample output, so we know
    the patched build is still the model the paper describes.
2.  For a handful of RPLC-relevant compounds: the four RPLC solute descriptors, the ensemble
    standard deviations, and the full 4x4 ensemble correlation matrix.
3.  Whether the off-diagonal terms are materially non-zero, measured two ways:
      (a) the largest |correlation| off the diagonal;
      (b) an assumption-free bound on how much dropping the off-diagonals can distort the
          variance of a linear combination m'd (which is what LSER computes). For any m,
              m' SIGMA m  /  m' diag(SIGMA) m
          lies between the smallest and largest eigenvalue of the correlation matrix, so
          those eigenvalues bound the error from assuming independence -- with no need to
          commit to a particular set of LSER system constants.
4.  That the individual ensemble member predictions are retrievable, which #37 needs to
    measure across-compound scaffold correlation of the residuals.

CALIBRATION WARNING
-------------------
Everything printed here is ENSEMBLE SPREAD: epistemic only, no aleatoric term, not
calibrated, and known to be understated by more than 10x against experimental disagreement
(#6). It is exposed here so that #14's stratified conformal layer has something to work
from. It is not a predictive standard deviation and must not be consumed as one.
"""

import contextlib
import io

import numpy as np

from chemprop_solvation.solvation_estimator import (
    load_SoluteML_estimator,
    RPLC_SOLUTE_TARGETS,
    SOLUTE_TARGETS,
)

# Drug-like / RPLC-relevant probes, plus the two solutes the upstream sample script uses.
COMPOUNDS = [
    ("octane (upstream sample)", "CCCCCCCC"),
    ("toluene (upstream sample)", "Cc1ccccc1"),
    ("benzoic acid", "OC(=O)c1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("propranolol", "CC(C)NCC(O)COc1cccc2ccccc12"),
    ("naproxen", "COc1ccc2cc(ccc2c1)C(C)C(O)=O"),
    ("ciprofloxacin", "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O"),
]

# The upstream repo's published expected output for its own sample script, used as a
# regression check that the patch did not disturb the predictions.
UPSTREAM_SAMPLE = {
    "CCCCCCCC": (
        [0.0033575539182320426, 0.0002618650031903691, 0.0012404804461126284,
         -3.3407873290287605e-05, 3.6623056859897627],
        [0.0001796549953840305, 0.00014738685559645396, 5.848635942567972e-05,
         0.00016496458160710926, 0.00331170315682501],
    ),
    "Cc1ccccc1": (
        [0.6076908247229844, 0.5312452682796652, 0.001725445237505657,
         0.14568006907274622, 3.3648599797284917],
        [0.00023044562651108525, 0.0003164094502622555, 4.004680160551078e-05,
         0.0001829383390569367, 0.003946935292885148],
    ),
}


def main():
    print("Loading SoluteML (5 folds x 5 models = 25 members)...")
    # load_checkpoint prints one line per tensor to stdout; swallow it.
    with contextlib.redirect_stdout(io.StringIO()):
        estimator = load_SoluteML_estimator()

    smiles = [[s] for _, s in COMPOUNDS]
    with contextlib.redirect_stdout(io.StringIO()):
        avg, var, valid, members, cov = estimator(smiles, return_ensemble=True)

    if len(valid) != len(COMPOUNDS):
        raise SystemExit(f"some SMILES were rejected: valid={valid}")

    avg = np.asarray(avg)          # (n, 5) -- upstream head, includes L
    var = np.asarray(var)          # (n, 5)
    print(f"member predictions tensor: {members.shape}   "
          f"(members, compounds, targets={RPLC_SOLUTE_TARGETS})")
    print(f"covariance tensor:         {cov.shape}")
    print(f"upstream head targets:     {SOLUTE_TARGETS}  "
          f"-- L present upstream, excluded from the exposed covariance\n")

    # --- 1. reproduction check -------------------------------------------------------
    print("=" * 78)
    print("1. Reproduction check against upstream's published sample output")
    print("=" * 78)
    for i, (_, smi) in enumerate(COMPOUNDS):
        if smi not in UPSTREAM_SAMPLE:
            continue
        exp_avg, exp_var = UPSTREAM_SAMPLE[smi]
        d = avg[i] - np.asarray(exp_avg)
        sd = np.sqrt(np.asarray(exp_var))
        print(f"  {smi}")
        print("        " + "".join(f"{t:>10}" for t in SOLUTE_TARGETS))
        print("  ours  " + "".join(f"{v:>10.4f}" for v in avg[i]))
        print("  paper " + "".join(f"{v:>10.4f}" for v in exp_avg))
        print("  delta " + "".join(f"{v:>10.4f}" for v in d))
        print("  d/sig " + "".join(f"{v:>10.2f}" for v in d / sd))
    print("\n  NOT a bit-reproduction -- see SOLUTEML_COVARIANCE_RESULTS.md. The deltas are"
          "\n  RDKit/descriptastorus drift over four years, and they are of the same order as"
          "\n  the ensemble spread itself, which is its own finding for #14.\n")

    # --- consistency: cov diagonal must equal the shipped per-descriptor variance -----
    diag_err = np.abs(np.array([np.diag(c) for c in cov]) - var[:, :4]).max()
    print(f"  max |diag(cov) - shipped variance| = {diag_err:.3e}  "
          f"(must be ~0: same tensor, ddof=0)\n")

    # --- 2/3. covariance -------------------------------------------------------------
    print("=" * 78)
    print("2. Per-compound descriptors, ensemble sigma, and correlation matrix")
    print("=" * 78)
    rows = []
    for i, (name, _) in enumerate(COMPOUNDS):
        sigma = np.sqrt(np.diag(cov[i]))
        R = cov[i] / np.outer(sigma, sigma)
        off = R[~np.eye(4, dtype=bool)]
        max_abs_off = np.abs(off).max()
        eig = np.linalg.eigvalsh(R)
        rows.append((name, max_abs_off, eig[0], eig[-1]))

        print(f"\n-- {name}")
        print("        " + "".join(f"{t:>10}" for t in RPLC_SOLUTE_TARGETS))
        print("  mean  " + "".join(f"{v:>10.4f}" for v in avg[i, :4]))
        print("  sigma " + "".join(f"{v:>10.4f}" for v in sigma))
        print("  correlation matrix:")
        for k, t in enumerate(RPLC_SOLUTE_TARGETS):
            print(f"    {t}   " + "".join(f"{v:>10.3f}" for v in R[k]))
        print(f"  max |off-diagonal r| = {max_abs_off:.3f}")
        print(f"  eigenvalues of R     = {np.array2string(eig, precision=3)}")
        print(f"  -> variance of any linear combination m'd is between "
              f"{eig[0]:.2f}x and {eig[-1]:.2f}x the independence assumption")

    # --- summary ---------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("3. Are the off-diagonals material?")
    print("=" * 78)
    print(f"{'compound':<28}{'max|r|':>10}{'lambda_min':>12}{'lambda_max':>12}")
    for name, mo, lo, hi in rows:
        print(f"{name:<28}{mo:>10.3f}{lo:>12.3f}{hi:>12.3f}")
    allmax = max(r[1] for r in rows)
    lo = min(r[2] for r in rows)
    hi = max(r[3] for r in rows)
    print(f"\nacross this probe set: max |off-diagonal r| = {allmax:.3f}; "
          f"eigenvalue span of R = [{lo:.3f}, {hi:.3f}]")

    # --- 4. per-member retrievability for #37 ---------------------------------------
    print("\n" + "=" * 78)
    print("4. Per-member predictions (what #37 needs)")
    print("=" * 78)
    print(f"members array {members.shape}: 25 members x {len(COMPOUNDS)} compounds x 4 targets")
    resid = members - members.mean(axis=0, keepdims=True)   # per-member deviation from mean
    # across-compound correlation of the per-member E deviation, first 4 compounds
    e_dev = resid[:, :4, 0]
    print("across-compound correlation of per-member E deviation (first 4 compounds):")
    print(np.array2string(np.corrcoef(e_dev, rowvar=False), precision=3))
    print("\n(That matrix is only computable because the members survive the API; the "
          "summarised variance cannot produce it.)")

    print("\nREMINDER: ensemble spread only. Epistemic, uncalibrated, no aleatoric term. "
          "Recalibration is #14, not this ticket.")


if __name__ == "__main__":
    main()
