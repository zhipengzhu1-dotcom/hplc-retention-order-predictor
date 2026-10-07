#!/usr/bin/env python3
"""Non-Gaussian descriptor-residual sampler for the #14 uncertainty budget.

Replaces the Gaussian assumption for (E, S, A, B) descriptor error with the
empirical residuals measured in FINDINGS.md. Written against the proposal in
an earlier design proposal (not in this repo), which got the shape question right and two other
things wrong; both corrections are implemented here and justified in the
docstrings below.

WHAT THIS IS NOT. Every table here is measured on 94 neutral small solutes that
are all inside SoluteML's training corpus. The residuals are therefore a FLOOR.
This sampler reproduces the error we can see; it does not claim to bound the
error on drug-like chemistry, where #31 showed the anchor does not reach.

Run `python residual_sampler.py` to regenerate the tables from residuals.csv and
execute the self-checks.
"""

import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = ("E", "S", "A", "B")
GRID = np.array([0.0, 0.05, 0.25, 0.50, 0.75, 0.95, 1.0])

# Measured on the 94 WSU compounds (prototype/soluteml-bias/residuals.csv).
# Pooled shape, all 94 compounds -- see CLASS_SHIFT for the weak-base handling.
QUANTILES = {
    "E": np.array([-0.703, -0.202, +0.016, +0.047, +0.081, +0.239, +0.480]),
    "S": np.array([-0.679, -0.094, +0.014, +0.045, +0.108, +0.294, +1.113]),
    "A": np.array([-0.271, -0.159, -0.034, +0.021, +0.026, +0.088, +0.218]),
    "B": np.array([-0.229, -0.064, -0.005, +0.017, +0.041, +0.079, +0.346]),
}

# Correction 1 of 3: the weak-base effect is carried as a SHIFT on the pooled
# shape, not as its own 15-point quantile table.
#
# The proposal built separate n=15 tables. Those anchor points do not survive
# resampling: E_base p95 is +0.457 with a bootstrap CI of [+0.248, +0.480] --
# a 0.23-wide interval on a 0.48-tall distribution -- and p50 is +0.141
# [+0.070, +0.231]. A seven-point table from fifteen observations is mostly
# reading noise into the shape.
#
# The class MEAN is resolved (+0.152 [+0.077, +0.236] for E), so spend the one
# degree of freedom the data supports and take the shape from all 94. This is
# also the conservative direction: the pooled shape keeps a left tail that the
# 15 bases do not exhibit (their observed min is -0.054 against the pool's
# -0.703), so bases are sampled wider than they were observed, not narrower.
CLASS_SHIFT = {"E": +0.152, "A": +0.042, "B": -0.034, "S": 0.0}
CLASS_SHIFT_RESOLVED = {"E": True, "A": True, "B": True, "S": False}

# Correction 2 of 3: descriptor residuals are strongly correlated ACROSS
# descriptors, and the proposal sampled them independently.
#
# Measured across-compound residual correlation (n=94). corr(S,B) = -0.67 is
# the big one. Propagated through the 332 unmasked WSU column vectors,
# independent sampling gives Var(log k) 2.83x the joint value -- sd overstated
# 1.68x. That inflates every interval and drags every P(A before B) toward 0.5,
# which is the failure #14 named as the trap and #37 already found for the
# system-constant covariance.
#
# Note this is the correlation of the ACTUAL ERRORS against measurement. It is
# far stronger than the SoluteML ensemble's own within-compound correlation
# (#37: E-S +0.29, rest |rho| < 0.06). The ensemble underestimates not only the
# size of its error but the shape of its correlation.
CORRELATION = np.array([
    [+1.00, -0.09, +0.24, -0.10],
    [-0.09, +1.00, -0.24, -0.67],
    [+0.24, -0.24, +1.00, -0.26],
    [-0.10, -0.67, -0.26, +1.00],
])

# A Gaussian copula does not deliver its input correlation when the marginals
# are this far from normal -- pushing the measured matrix through the empirical
# inverse CDFs attenuates it (worst cell off by 0.09, and the propagated
# Var(log k) ratio lands at 0.443 instead of the measured 0.354). So the copula
# takes a pre-calibrated matrix chosen to make the SAMPLED linear correlation
# equal the measured one. `calibrate_copula()` below regenerates this; it
# converges in ~4 iterations and the check at the bottom asserts the result.
COPULA_CORRELATION = np.array([
    [+1.0000, -0.1092, +0.2650, -0.1155],
    [-0.1092, +1.0000, -0.2803, -0.7611],
    [+0.2650, -0.2803, +1.0000, -0.2898],
    [-0.1155, -0.7611, -0.2898, +1.0000],
])


def _inverse_cdf(u: np.ndarray, q: np.ndarray, tail_extend: bool) -> np.ndarray:
    """Piecewise-linear empirical inverse CDF.

    Correction 3 of 3: the observed min and max are single compounds, so a table
    that stops there is hard-bounded by one observation apiece. The pooled S
    maximum (+1.113) is progesterone alone; resampled, that maximum ranges over
    [+0.364, +1.113]. `tail_extend` continues the outermost segment's slope
    beyond the observed extremes rather than clipping, so the tail is open.
    It is a stopgap: a fitted GPD tail is the proper instrument and needs more
    compounds than 94 to fit.
    """
    x = np.interp(u, GRID, q)
    if not tail_extend:
        return x
    lo_slope = (q[1] - q[0]) / (GRID[1] - GRID[0])
    hi_slope = (q[-1] - q[-2]) / (GRID[-1] - GRID[-2])
    below, above = u < GRID[0], u > GRID[-1]
    x[below] = q[0] - lo_slope * (GRID[0] - u[below])
    x[above] = q[-1] + hi_slope * (u[above] - GRID[-1])
    return x


def sample_residuals(n_samples: int, is_weak_base=False, rng=None,
                     tail_extend: bool = True) -> np.ndarray:
    """Draw joint (E, S, A, B) residuals. Returns (n_samples, 4).

    Gaussian copula: correlated normal scores -> uniforms -> per-descriptor
    empirical inverse CDF. Marginals are the measured skewed, heavy-tailed
    shapes; the dependence is the measured correlation. Sampling the marginals
    independently -- the proposal's design -- gets the shapes right and the
    joint wrong, and the joint is what elution order depends on.
    """
    rng = np.random.default_rng() if rng is None else rng
    base = np.broadcast_to(np.atleast_1d(is_weak_base), (n_samples,))
    L = np.linalg.cholesky(_nearest_psd(COPULA_CORRELATION))
    z = rng.standard_normal((n_samples, 4)) @ L.T
    u = _std_normal_cdf(z)
    out = np.empty((n_samples, 4))
    for j, d in enumerate(D):
        out[:, j] = _inverse_cdf(u[:, j], QUANTILES[d], tail_extend)
        if CLASS_SHIFT[d]:
            out[base, j] += CLASS_SHIFT[d]
    return out


def _std_normal_cdf(z: np.ndarray) -> np.ndarray:
    """Phi(z) without scipy -- neither project venv carries it."""
    import math
    return 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))


def calibrate_copula(n: int = 600_000, iters: int = 60, seed: int = 7,
                     tol: float = 0.004) -> np.ndarray:
    """Regenerate COPULA_CORRELATION: fixed-point search for the input matrix
    whose SAMPLED linear correlation equals CORRELATION."""
    R = CORRELATION.copy()
    rng = np.random.default_rng(seed)
    for _ in range(iters):
        L = np.linalg.cholesky(_nearest_psd(R))
        u = _std_normal_cdf(rng.standard_normal((n, 4)) @ L.T)
        s = np.column_stack([_inverse_cdf(u[:, j], QUANTILES[d], True)
                             for j, d in enumerate(D)])
        err = CORRELATION - np.corrcoef(s.T)
        if np.abs(err).max() < tol:
            break
        R = _nearest_psd(np.clip(R + 0.8 * err, -0.999, 0.999))
        np.fill_diagonal(R, 1.0)
    return R


def _nearest_psd(R: np.ndarray) -> np.ndarray:
    """The published correlation is rounded to 2dp and can lose PSD by a hair."""
    w, V = np.linalg.eigh(R)
    if w.min() > 1e-10:
        return R
    R2 = V @ np.diag(np.clip(w, 1e-10, None)) @ V.T
    d = np.sqrt(np.diag(R2))
    return R2 / np.outer(d, d)


# ------------------------------------------------------------------ self-check
def _regenerate_and_check() -> None:
    rows = list(csv.DictReader(open(os.path.join(HERE, "residuals.csv"))))
    err = np.array([[float(r["err_" + k]) for k in D] for r in rows])
    base = np.array([int(r["weak_base"]) for r in rows], bool)

    print("Regenerated from residuals.csv (compare to the tables above):")
    for j, d in enumerate(D):
        q = np.percentile(err[:, j], GRID * 100)
        print(f"  {d}: [" + ", ".join(f"{v:+.3f}" for v in q) + "]")
    print("\n  correlation:")
    C = np.corrcoef(err.T)
    for i, d in enumerate(D):
        print(f"    {d}: " + "  ".join(f"{C[i, j]:+.2f}" for j in range(4)))
    print("\n  class shift (weak base - neutral):")
    for j, d in enumerate(D):
        print(f"    {d}: {err[base, j].mean() - err[~base, j].mean():+.3f}")

    rng = np.random.default_rng(20260816)
    s = sample_residuals(400_000, is_weak_base=False, rng=rng)

    print("\nSelf-checks:")
    Cs = np.corrcoef(s.T)
    worst = np.abs(Cs - CORRELATION).max()
    print(f"  sampled correlation matches measured to {worst:.3f}  "
          f"{'OK' if worst < 0.01 else 'FAIL'}")

    for j, d in enumerate(D):
        want = np.percentile(err[:, j], [25, 50, 75])
        got = np.percentile(s[:, j], [25, 50, 75])
        ok = np.abs(want - got).max() < 0.02
        print(f"  {d} quartiles sampled {np.round(got, 3)} vs measured "
              f"{np.round(want, 3)}  {'OK' if ok else 'FAIL'}")

    # the number that motivates the copula
    fits = np.array([[float(r[k]) for k in ("e", "s", "a", "b")]
                     for r in csv.DictReader(open(os.path.join(
                         HERE, "..", "..", "sources", "wsu-lser",
                         "wsu2019-system-constants.csv")))
                     if str(r.get("incomplete_wetting", "")).strip().lower()
                     not in ("true", "1", "yes")])
    S_joint = np.cov(s.T)
    var_joint = np.einsum("fi,ij,fj->f", fits, S_joint, fits)
    var_ind = np.einsum("fi,i,fi->f", fits, np.diag(S_joint), fits)
    ratio = float(np.median(var_joint / var_ind))
    print(f"  Var(log k) joint/independent over {len(fits)} fits: {ratio:.3f} "
          f"-> independent sampling would overstate variance "
          f"{1 / ratio:.2f}x (sd {1 / np.sqrt(ratio):.2f}x)")

    b = sample_residuals(200_000, is_weak_base=True, rng=rng)
    print(f"  weak-base E mean {b[:, 0].mean():+.3f} vs neutral "
          f"{s[:, 0].mean():+.3f}  (shift {CLASS_SHIFT['E']:+.3f})")


if __name__ == "__main__":
    _regenerate_and_check()
