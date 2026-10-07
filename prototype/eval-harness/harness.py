"""Evaluation harness for #45, built before any measured data exists.

That ordering is the point. #15 pre-registered the acceptance thresholds so they
could not be chosen to be passed; writing the code that computes them *after*
seeing the data would reopen the same door. This module is the other half of
that commitment.

Specification: `experiments/calibration-preregistration.md` (#15), with the
reporting contract from `CONTEXT.md` (#27 as extended by #31).

The load-bearing design decision is `Figure`. CONTEXT.md requires that the
harness have "no code path returning an accuracy figure unaccompanied by its
applicable baseline, its `n`, its interval and its condition axis - enforced by
the function signature and asserted in a test, not left to reviewer discipline."
So no function here returns a float. They return `Figure`, which cannot be
constructed without its six fields and cannot be coerced back to a bare number.

Stdlib + numpy only.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

NOMINAL = 0.90          # nominal coverage for interval / conformal work
LAMBDA_CALIBRATED = 2.0     # <= this: calibrated               (#15 s.3)
LAMBDA_RESCOPE = 2.5        # >  this: under-dispersed, #19 rescopes


class BareNumberError(RuntimeError):
    """Raised when a Figure is used as if it were a plain number."""


class MixedStrataError(RuntimeError):
    """Raised when rows from different applicability strata would be pooled.

    #40's firewall, at the point it actually bites. Refused rows never score;
    degraded rows score in their own stratum. Pooling a degraded residual into
    an in-envelope headline is the failure this guards, and it is silent
    otherwise: the pooled number looks like the clean one, only worse.
    """


@dataclass(frozen=True)
class Figure:
    """A number that cannot exist apart from what makes it interpretable.

    Six mandatory fields: #27's four (value, axis, baseline, n, interval) plus
    #31's two (reference identity, reference kind). `reference_kind` decides the
    claim verb and is not free text: a figure referenced to model output can
    never claim accuracy, because agreeing with that corpus partly means
    mimicking the tool we are trying to beat (#5).
    """

    value: float
    axis: str
    baseline: str
    baseline_value: float | None
    n_compounds: int
    interval: tuple[float, float]
    reference_identity: str
    reference_kind: str          # "measurement" | "model_output"
    units: str = ""
    stratum: str | None = None   # applicability stratum, when domain-gated (#40)

    def __post_init__(self) -> None:
        if self.reference_kind not in ("measurement", "model_output"):
            raise ValueError(
                f"reference_kind must be 'measurement' or 'model_output', "
                f"got {self.reference_kind!r}")
        if not self.axis or not self.baseline or not self.reference_identity:
            raise ValueError("axis, baseline and reference_identity are mandatory")
        if self.n_compounds < 1:
            raise ValueError("n_compounds is counted in COMPOUNDS and must be >= 1")
        lo, hi = self.interval
        if lo > hi:
            raise ValueError(f"interval is inverted: {self.interval}")

    @property
    def verb(self) -> str:
        """#31: the claim verb follows from the reference type."""
        return "accurate to" if self.reference_kind == "measurement" else "agrees with"

    def render(self) -> str:
        lo, hi = self.interval
        u = self.units
        base = (f" vs baseline {self.baseline_value:.4g}{u}"
                if self.baseline_value is not None else "")
        # The stratum sits next to n, not at the end: a degraded figure read
        # without it is the in-envelope claim, which is the whole hazard.
        strat = f", stratum={self.stratum}" if self.stratum else ""
        return (f"{self.value:.4g}{u} [{lo:.4g}, {hi:.4g}] "
                f"(compound-level bootstrap, n={self.n_compounds} compounds"
                f"{strat}); "
                f"axis={self.axis}; baseline={self.baseline}{base}; "
                f"{self.verb} {self.reference_identity}")

    # -- the guard. A Figure must never degrade into a bare number. ----------
    def __float__(self):
        raise BareNumberError(
            "A Figure may not be coerced to float: that is exactly the bare "
            "accuracy number CONTEXT.md forbids. Use .value for internal "
            "composition, or .render() to report it.")

    def __format__(self, spec):
        if spec:
            raise BareNumberError(
                "A Figure may not be format-spec'd into a bare number. "
                "Use .render().")
        return self.render()

    def __str__(self):
        return self.render()


# ------------------------------------------------------- the domain firewall
#
# #40: only in-domain predictions may enter a metric. `model/applicability.py`
# decides per compound; this is where that decision reaches the scoring path.
# The harness deliberately does not import the model layer -- it takes the
# strata as plain strings, so it stays runnable with numpy alone and the
# domain object stays free to change.

STRATUM_ORDER = {"in_envelope": 0, "degraded": 1}


def pair_stratum(a: str | None, b: str | None) -> str | None:
    """The stratum a PAIR scores in, given its two compounds'.

    A pair is only as clean as its worse compound, and a pair touching a
    refused compound cannot be scored at all. Getting this wrong is how a
    degraded compound launders itself into an in-envelope figure: it appears in
    14 of the 105 pairs of a 15-compound slice, so one unhandled compound
    contaminates 13% of the pairs while every compound-level check still passes.
    """
    if a is None or b is None:
        return None
    if a not in STRATUM_ORDER or b not in STRATUM_ORDER:
        raise ValueError(f"unknown stratum: {a!r}, {b!r}")
    return a if STRATUM_ORDER[a] >= STRATUM_ORDER[b] else b


def strata_for_pairs(compound_strata: Sequence[str | None],
                     pairs: Sequence[tuple[int, int]]) -> list[str | None]:
    return [pair_stratum(compound_strata[i], compound_strata[j])
            for i, j in pairs]


def partition_by_stratum(strata: Sequence[str | None]) -> dict[str, np.ndarray]:
    """{stratum: row indices}. Rows with stratum None are dropped, not grouped.

    Returned as a mapping rather than a mask so that a caller who wants a
    pooled number has to write the pooling itself, in the open.
    """
    out: dict[str, list[int]] = {}
    for k, s in enumerate(strata):
        if s is not None:
            out.setdefault(s, []).append(k)
    return {s: np.asarray(v, dtype=int) for s, v in out.items()}


def require_single_stratum(strata: Sequence[str | None]) -> str:
    """Guard for any pooled statistic. Raises rather than averaging."""
    present = {s for s in strata if s is not None}
    if not present:
        raise MixedStrataError("no scorable rows: every row was refused")
    if len(present) > 1:
        raise MixedStrataError(
            f"rows span {sorted(present)}; a pooled figure over mixed "
            f"applicability strata is exactly what #40 forbids. Use "
            f"partition_by_stratum() and report one Figure per stratum.")
    return present.pop()


def score_by_stratum(strata: Sequence[str | None],
                     make_figure: Callable[[np.ndarray, str], Figure],
                     ) -> dict[str, Figure]:
    """The domain-gated entry point: one Figure per stratum, never pooled.

    `make_figure(rows, stratum)` receives the row indices for one stratum and
    must return that stratum's Figure. Refused rows never reach it.
    """
    figures = {}
    for stratum, rows in partition_by_stratum(strata).items():
        fig = make_figure(rows, stratum)
        if fig.stratum != stratum:
            raise MixedStrataError(
                f"Figure for stratum {stratum!r} is labelled {fig.stratum!r}; "
                f"a mislabelled stratum is worse than none")
        figures[stratum] = fig
    return figures


# ------------------------------------------------------------------ bootstrap
def bootstrap_over_compounds(
    n_compounds: int,
    statistic: Callable[[np.ndarray], float],
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 20260816,
) -> tuple[float, float, float]:
    """Resample COMPOUNDS, never pairs. Returns (point, lo, hi).

    #15 s.1: pairs share compounds, so one badly predicted compound corrupts
    every pair it appears in. A binomial interval on pair counts understates the
    SD by 1.35x at 15 compounds and 2.40x at 71 - and the distortion GROWS with
    study size, so the larger the study the more a naive interval flatters us.
    """
    rng = np.random.default_rng(seed)
    point = statistic(np.arange(n_compounds))
    draws = np.empty(n_boot)
    for b in range(n_boot):
        draws[b] = statistic(rng.integers(0, n_compounds, n_compounds))
    lo, hi = np.nanpercentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(point), float(lo), float(hi)


def two_sd(n_compounds: int, statistic, **kw) -> float:
    """2 SD of the compound-level bootstrap - the margin rule's unit (#15 s.1)."""
    rng = np.random.default_rng(kw.pop("seed", 20260816))
    n_boot = kw.pop("n_boot", 2000)
    draws = np.array([statistic(rng.integers(0, n_compounds, n_compounds))
                      for _ in range(n_boot)])
    return float(2.0 * np.nanstd(draws))


# -------------------------------------------------------------- metric A
def pairwise_order_accuracy(p_before: np.ndarray, truth: np.ndarray) -> float:
    """Fraction of pairs whose claimed direction matches observation."""
    claimed = p_before >= 0.5
    return float(np.mean(claimed == truth.astype(bool)))


def expected_calibration_error(p_before: np.ndarray, truth: np.ndarray,
                               n_bins: int = 10) -> float:
    """Metric A: is P(A before B) that number asked of itself? (#15 s.1)

    Folded to the claimed direction, so a claim of 0.1 counts as a 0.9 claim in
    the other direction - the user reads confidence, not orientation.
    """
    conf = np.where(p_before >= 0.5, p_before, 1.0 - p_before)
    correct = (p_before >= 0.5) == truth.astype(bool)
    edges = np.linspace(0.5, 1.0, n_bins + 1)
    ece, total = 0.0, len(conf)
    for i in range(n_bins):
        m = (conf >= edges[i]) & (conf < edges[i + 1] if i < n_bins - 1
                                  else conf <= edges[i + 1])
        if m.sum():
            ece += m.sum() / total * abs(correct[m].mean() - conf[m].mean())
    return float(ece)


def reliability_table(p_before: np.ndarray, truth: np.ndarray,
                      n_bins: int = 10) -> list[dict]:
    conf = np.where(p_before >= 0.5, p_before, 1.0 - p_before)
    correct = (p_before >= 0.5) == truth.astype(bool)
    edges = np.linspace(0.5, 1.0, n_bins + 1)
    out = []
    for i in range(n_bins):
        m = (conf >= edges[i]) & (conf < edges[i + 1] if i < n_bins - 1
                                  else conf <= edges[i + 1])
        if m.sum():
            out.append(dict(bin_lo=float(edges[i]), bin_hi=float(edges[i + 1]),
                            n_pairs=int(m.sum()),
                            claimed=float(conf[m].mean()),
                            observed=float(correct[m].mean())))
    return out


# -------------------------------------------------------------- metric B
def conformal_lambda(residuals: np.ndarray, half_widths: np.ndarray,
                     nominal: float = NOMINAL) -> float:
    """The inflation the model's own interval needs to reach nominal coverage.

    #15 s.3: coverage is NOT the test and cannot be - conformal reaches nominal
    by construction if you widen far enough. What can fail is sharpness, i.e.
    the width needed to get there. lambda is that width, as a multiplier.
    """
    hw = np.asarray(half_widths, float)
    if np.any(hw <= 0):
        raise ValueError("half_widths must be positive")
    scores = np.abs(np.asarray(residuals, float)) / hw
    return float(np.quantile(scores, nominal))


def lambda_verdict(lam: float) -> str:
    if lam <= LAMBDA_CALIBRATED:
        return "calibrated"
    if lam <= LAMBDA_RESCOPE:
        return "calibrated but materially degraded"
    return "under-dispersed, not calibrated - #19 rescopes"


def interval_coverage(residuals: np.ndarray, half_widths: np.ndarray) -> float:
    return float(np.mean(np.abs(residuals) <= half_widths))


# ------------------------------------------------------------- gating stats
def relative_error_reduction(model_err: np.ndarray, baseline_err: np.ndarray) -> float:
    """(baseline - model)/baseline. #27: never absolute percentage-point margins,
    because headroom differs by an order of magnitude across axes."""
    b = float(np.mean(baseline_err))
    if b == 0:
        raise ZeroDivisionError("baseline error is zero; relative reduction undefined")
    return float((b - np.mean(model_err)) / b)


# ------------------------------------------------------------------ verdicts
@dataclass(frozen=True)
class Verdict:
    """#15 s.2: three verdicts, the third named in advance.

    Inconclusive = the achieved 2 SD compound-level interval exceeds the observed
    gain. It does not block #19 and does not count as a pass, and it is reported
    WITH the sample size that would have settled it - because the failure mode is
    discovering it afterwards, when it is indistinguishable from a pass.
    """
    outcome: str                  # "pass" | "fail" | "inconclusive"
    gain: float
    margin_2sd: float
    n_compounds: int
    required_n: int | None
    axis: str

    def render(self) -> str:
        s = (f"{self.outcome.upper()} on {self.axis}: gain {self.gain:+.4g} "
             f"against a 2 SD compound-level margin of {self.margin_2sd:.4g} "
             f"(n={self.n_compounds} compounds)")
        if self.outcome == "inconclusive":
            s += (f"; this axis could not be resolved at the available sample "
                  f"size - n={self.required_n} compounds would settle it")
        return s


def decide(gain: float, margin_2sd: float, n_compounds: int, axis: str) -> Verdict:
    """A margin counts only if it exceeds 2 SD of the compound-level bootstrap."""
    if gain > margin_2sd:
        outcome, req = "pass", None
    elif gain < -margin_2sd:
        outcome, req = "fail", None
    else:
        outcome = "inconclusive"
        # margin shrinks as 1/sqrt(n); n needed for margin < |gain|
        req = (int(np.ceil(n_compounds * (margin_2sd / abs(gain)) ** 2))
               if gain != 0 else None)
    return Verdict(outcome, float(gain), float(margin_2sd), int(n_compounds),
                   req, axis)


# -------------------------------------------------------------------- strata
MERGE_ORDER = ("confidence_tertile", "phi_range", "column_type")
NEVER_MERGE = "ionisation_regime"
HARD_FLOOR, PRACTICAL_FLOOR = 19, 50


def merge_strata(strata: dict[tuple, int],
                 keys: Sequence[str],
                 floor: int = PRACTICAL_FLOOR) -> dict[tuple, tuple]:
    """Merge undersized strata upward along #15 s.4's fixed hierarchy.

    `strata` maps a key tuple (ordered as `keys`) to its compound count. Returns
    a map from original key to the key it was merged into.

    Ionisation regime NEVER merges: it has the largest measured effect and it is
    the axis #34's two-regime claim rests on. A stratum still below the hard
    floor after all permitted merges is not calibratable at all - it stays in
    block-plus-opt-in (#40) and never enters the pooled figure, so it can
    neither drag the headline inflation nor hide behind it.
    """
    if NEVER_MERGE not in keys:
        raise ValueError(f"{NEVER_MERGE} must be one of the stratum keys")
    idx = {k: i for i, k in enumerate(keys)}
    assign = {k: k for k in strata}
    counts = dict(strata)

    for field in MERGE_ORDER:
        if field not in idx:
            continue
        pos = idx[field]
        # collapse this field for any group still under the floor
        groups: dict[tuple, list[tuple]] = {}
        for k in list(counts):
            if counts[k] >= floor:
                continue
            collapsed = tuple("*" if i == pos else v for i, v in enumerate(k))
            groups.setdefault(collapsed, []).append(k)
        for collapsed, members in groups.items():
            if len(members) < 2:
                continue
            total = sum(counts.pop(m) for m in members)
            counts[collapsed] = counts.get(collapsed, 0) + total
            for orig, tgt in list(assign.items()):
                if tgt in members:
                    assign[orig] = collapsed
        if all(c >= floor for c in counts.values()):
            break
    return assign


def calibratable(count: int) -> bool:
    """Below the hard floor the 90% quantile IS the sample maximum, so the
    interval is not meaningful and the stratum cannot be conformally calibrated."""
    return count >= HARD_FLOOR
