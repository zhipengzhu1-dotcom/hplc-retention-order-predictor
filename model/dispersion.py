#!/usr/bin/env python3
"""Peak width with physics in it — the Knox model, not a nominal plate count.

`spec/scenario-ensemble.md` requires

    width^2 = column_variance(N, t_R) + extra_column_variance + gradient_compression

with `N` **sampled per scenario** so the slot carries the Knox prior's
uncertainty rather than its point estimate. Until now `N` was a flat 10000 +-20%
- a number, not a model. That is #44 punch-list item 3.

`N` is not a column property. It depends on particle size, flow rate and the
solute's diffusion coefficient, so it belongs to the METHOD, and the Knox
equation is how the method determines it:

    h = A*v^(1/3) + B/v + C*v        (reduced plate height)
    v = u*dp/Dm                      (reduced velocity)
    N = L/(h*dp)

Per the map's cut rule, van Deemter/Knox is settled science: the choice is
pinned, the interface defined, the source cited, and no textbook is paraphrased.
`A, B, C` ship as engineering priors BY PARTICLE CLASS at +-10-30% in plate
count (#9, #13) - and the class matters, because #30's two arms are not the same
class.

Reference: Knox, J. Chromatogr. Sci. 15 (1977) 352; the reduced-parameter
formulation as used throughout Poole's work.

numpy only.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Knox reduced parameters as engineering priors by particle class.
# Core-shell packs more uniformly (lower A) and has a shorter diffusion path
# through the porous shell (lower C), which is the whole point of the format.
KNOX_PRIORS = {
    "fully_porous": dict(A=(1.00, 0.15), B=(2.00, 0.30), C=(0.050, 0.010)),
    "core_shell":   dict(A=(0.60, 0.10), B=(1.50, 0.25), C=(0.020, 0.005)),
}

DM_CM2_S = 1.0e-5        # small-molecule diffusion coefficient, aqueous-organic
POROSITY = 0.65


@dataclass(frozen=True)
class ColumnGeometry:
    length_cm: float
    id_cm: float
    dp_cm: float
    particle_class: str

    @property
    def cross_section_cm2(self) -> float:
        return np.pi * (self.id_cm / 2.0) ** 2

    def t0_min(self, flow_ml_min: float) -> float:
        v0 = self.cross_section_cm2 * self.length_cm * POROSITY
        return v0 / flow_ml_min

    def linear_velocity_cm_s(self, flow_ml_min: float) -> float:
        return (flow_ml_min / 60.0) / (self.cross_section_cm2 * POROSITY)


def reduced_velocity(geom: ColumnGeometry, flow_ml_min: float, dm=DM_CM2_S):
    return geom.linear_velocity_cm_s(flow_ml_min) * geom.dp_cm / dm


def plate_count(geom: ColumnGeometry, flow_ml_min: float, rng=None, size=None,
                dm=DM_CM2_S):
    """N from Knox. Draws A, B, C from their priors when rng is given."""
    pri = KNOX_PRIORS[geom.particle_class]
    if rng is None:
        A, B, C = pri["A"][0], pri["B"][0], pri["C"][0]
    else:
        A = rng.normal(*pri["A"], size)
        B = rng.normal(*pri["B"], size)
        C = rng.normal(*pri["C"], size)
        A, B, C = np.clip(A, 0.1, None), np.clip(B, 0.1, None), np.clip(C, 1e-4, None)
    v = reduced_velocity(geom, flow_ml_min, dm)
    h = A * v ** (1.0 / 3.0) + B / v + C * v
    return geom.length_cm / (h * geom.dp_cm)


def peak_sd_min(t_r_min, n_plates, sd_extra_col_min=0.0):
    """Total peak standard deviation. Column term is the isocratic result."""
    return np.sqrt((t_r_min / np.sqrt(n_plates)) ** 2 + sd_extra_col_min ** 2)


def resolution(t1, t2, sd1, sd2):
    return np.abs(t2 - t1) / (2.0 * (sd1 + sd2))


# #30's two arms, at the protocol's suggested geometry
EVO = ColumnGeometry(15.0, 0.46, 5e-4, "core_shell")        # Kinetex EVO C18
SUPELCOSIL = ColumnGeometry(15.0, 0.46, 5e-4, "fully_porous")  # Supelcosil LC-18


def main():
    flow = 1.0
    print("Knox plate count for #30's two arms (150 x 4.6 mm, 5 um, 1.0 mL/min)\n")
    v = reduced_velocity(EVO, flow)
    print(f"  reduced velocity v = {v:.2f}   t0 = {EVO.t0_min(flow):.2f} min\n")
    print(f"  {'arm':<26}{'class':<15}{'h':>7}{'N':>10}{'peak sd @ k=10':>18}")
    rng = np.random.default_rng(20260816)
    for name, g in (("Kinetex EVO C18", EVO), ("Supelcosil LC-18", SUPELCOSIL)):
        n = plate_count(g, flow)
        pri = KNOX_PRIORS[g.particle_class]
        h = (pri["A"][0] * v ** (1 / 3) + pri["B"][0] / v + pri["C"][0] * v)
        tr = g.t0_min(flow) * 11.0
        print(f"  {name:<26}{g.particle_class:<15}{h:>7.2f}{n:>10.0f}"
              f"{peak_sd_min(tr, n):>18.4f}")

    n_evo = plate_count(EVO, flow, rng, 20000)
    n_sup = plate_count(SUPELCOSIL, flow, rng, 20000)
    print(f"\n  N ratio EVO/Supelcosil: {np.mean(n_evo) / np.mean(n_sup):.2f}x")
    r = np.sqrt(np.mean(n_sup) / np.mean(n_evo))
    print(f"  peak width on EVO is {r:.2f}x that on Supelcosil "
          f"({100 * (1 - r):.0f}% narrower)")
    print(f"\n  N spread from the Knox priors (per-scenario draws):")
    for nm, n in (("EVO (core-shell)", n_evo), ("Supelcosil (fully porous)", n_sup)):
        print(f"    {nm:<28}{np.mean(n):>8.0f} ± {np.std(n):>6.0f}"
              f"   ({100 * np.std(n) / np.mean(n):.0f}%)")
    print("  (#9/#13 budget the Knox priors at +-10-30% in plate count)")

    print("\n  Consequence for a close pair separated by 0.05 min at k = 10:")
    tr = EVO.t0_min(flow) * 11.0
    for nm, g in (("Kinetex EVO C18", EVO), ("Supelcosil LC-18", SUPELCOSIL)):
        n = plate_count(g, flow)
        sd = peak_sd_min(tr, n)
        print(f"    {nm:<26}Rs = {resolution(tr, tr + 0.05, sd, sd):.2f}")
    print("    -> the arms are NOT comparable on resolution, which is why #30")
    print("       records that caveat. Now it is quantified rather than asserted.")


if __name__ == "__main__":
    main()
