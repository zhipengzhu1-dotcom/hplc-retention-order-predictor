"""The conformal stratum producer — `stratum[n]` in the ensemble spec.

The spec has required this field since #14/#44 and nothing produced it, so every
`ScenarioEnsemble` built so far has carried a placeholder. Two rules govern it
and they come from different tickets:

- **#42: the cut is method-independent.** A compound's calibration class may not
  change because the user changed column or modifier. So the summary is the
  plain sum of the four descriptor ensemble SDs, which is a property of the
  compound alone — *not* the LSER-propagated spread, which is per method. The
  measured cost of that choice is nil: #42 found the plain sum stratifies
  marginally BETTER (2.5x vs 2.2x hi/lo ratio in retention space).

- **#15 s.4: ionisation regime is an axis that never merges.** It carries the
  largest measured effect and it is what #34's two-regime claim rests on. It is
  therefore part of the key, not a tie-breaker.

**The tertile boundaries are frozen from the anchor, and that is the whole
point.** Computing tertiles over whatever batch is in front of you would make a
compound's stratum depend on its company: the same compound would be "confident"
in one run and "uncertain" in the next, and a conformal interval calibrated on a
moving definition means nothing. The cuts below are the WSU-94's own tertiles of
the summed descriptor SD, computed once (`prototype/soluteml-bias/residuals.csv`)
and pinned here.

⚠ They are pinned from an anchor that is 100% inside SoluteML's training corpus,
so the SD distribution they cut is narrower than a real one — and this is not a
worry, it is measured. The cuts split the anchor 32/31/31 by construction, but
**all nine thin-slice drugs land in `uncertain`** (sd_sum 0.169–0.431 against the
anchor's upper cut of 0.1635). The tertiles are therefore not equally populated
on real chemistry and nobody should expect them to be; `uncertain` is where
drug-like compounds live. The direction is the safe one, and it lines up with the
applicability domain putting the same nine compounds in `degraded`.
"""

from __future__ import annotations

from typing import Iterable, Sequence

# WSU-94 tertile boundaries of sum(ensemble SD over E, S, A, B). Frozen.
TERTILE_CUTS = (0.1286, 0.1635)
ANCHOR_N = 94
ANCHOR_SD_RANGE = (0.0857, 0.3110)

CONFIDENCE_LEVELS = ("confident", "typical", "uncertain")
IONISATION_REGIMES = ("neutral", "partial", "ionised")


def confidence_tertile(sd_sum: float) -> str:
    """Which frozen tertile of the anchor's descriptor-SD distribution."""
    if sd_sum < 0:
        raise ValueError(f"sd_sum must be non-negative, got {sd_sum}")
    if sd_sum <= TERTILE_CUTS[0]:
        return "confident"
    if sd_sum <= TERTILE_CUTS[1]:
        return "typical"
    return "uncertain"


def ionisation_regime(f_neutral: float, *,
                      neutral_above: float = 0.95,
                      ionised_below: float = 0.05) -> str:
    """Regime from the population-weighted neutral fraction.

    Cut on `f_neutral` rather than on pKa-versus-pH because that is what the
    ensemble actually carries per operating point, and because a compound with
    two pKa values has no single distance-from-pH to speak of.
    """
    if not 0.0 <= f_neutral <= 1.0:
        raise ValueError(f"f_neutral is a fraction, got {f_neutral}")
    if f_neutral >= neutral_above:
        return "neutral"
    if f_neutral <= ionised_below:
        return "ionised"
    return "partial"


def stratum_key(f_neutral: float, sd_sum: float) -> tuple[str, str]:
    """The full key, ordered to match harness.merge_strata's `keys`."""
    return (ionisation_regime(f_neutral), confidence_tertile(sd_sum))


def stratum_label(f_neutral: float, sd_sum: float) -> str:
    """The flat string the ensemble spec's `stratum[n]` carries."""
    return "/".join(stratum_key(f_neutral, sd_sum))


def assign(f_neutral: Sequence[float], sd_sum: Sequence[float]) -> list[str]:
    """Per-compound labels for `ScenarioEnsemble.stratum`.

    `f_neutral` is per compound here, not per operating point: the stratum is
    method-independent by #42, so a compound spanning regimes across the pH grid
    takes its WORST regime — the one carrying the most uncertainty — rather than
    changing class between operating points, which the spec forbids outright.
    """
    if len(f_neutral) != len(sd_sum):
        raise ValueError(
            f"f_neutral and sd_sum must be the same length, got "
            f"{len(f_neutral)} and {len(sd_sum)}")
    return [stratum_label(f, s) for f, s in zip(f_neutral, sd_sum)]


def worst_regime_f_neutral(per_point: Iterable[float]) -> float:
    """The `f_neutral` whose regime is worst across operating points.

    Ordering is neutral (best) -> partial -> ionised, and `partial` outranks
    `ionised` for uncertainty: #42 measured pKa owning 0.67 of the variance for
    partially ionised compounds against 0.92 for D in the fully ionised case,
    but the partial regime is where the *pH axis* moves the answer, which is the
    axis the product is sold on. Returns the representative value, so the caller
    can keep using one number.
    """
    vals = list(per_point)
    if not vals:
        raise ValueError("no operating points")
    rank = {"partial": 0, "ionised": 1, "neutral": 2}
    return min(vals, key=lambda f: (rank[ionisation_regime(f)], f))
