#!/usr/bin/env python3
"""#25: if t0 moves with pH and we hold it fixed, what does it cost?

The ticket asks whether the effect matters at our operating range and whether it
is swamped by other errors. Both are computable, because the damage has a
specific shape: `k = (tR - t0)/t0` is most sensitive to t0 exactly when k is
small, and the IONISED form of an acid is the smallest k in the whole
experiment. So a t0 error does not average out across the pH sweep - it lands
almost entirely on the ionised plateau, which is the plateau that defines D.

Reported t0 swing (#3, via #25): 0.83 -> 0.65 min up to pH 11, about -22%.

numpy only. Deterministic.
"""
import numpy as np

T0_ACID, T0_BASE = 0.83, 0.65          # min, from #3
SWING = (T0_BASE - T0_ACID) / T0_ACID


def t0_of_ph(ph, pka_sil, t0_acid=T0_ACID, t0_base=T0_BASE):
    """Residual silanols ionise over a range; t0 falls as they do."""
    f_ionised = 1.0 / (1.0 + 10.0 ** (pka_sil - np.asarray(ph, float)))
    return t0_acid + (t0_base - t0_acid) * f_ionised


def k_true(ph, k_neutral, D, pka):
    f_neutral = 1.0 / (1.0 + 10.0 ** (np.asarray(ph, float) - pka))
    return k_neutral * (f_neutral + (1.0 - f_neutral) * 10.0 ** (-D))


def fit_D_with_assumed_t0(ph, tR, t0_assumed, k_neutral_true, pka):
    """Recover D from the two plateaus, using whatever t0 we believe.

    Uses the plateau ratio directly rather than a full 3-parameter fit: it is
    the same quantity and it makes the mechanism visible.
    """
    k_obs = (tR - t0_assumed) / t0_assumed
    k_obs = np.maximum(k_obs, 1e-9)
    lo = k_obs[np.argmin(ph)]        # most acidic  -> neutral plateau
    hi = k_obs[np.argmax(ph)]        # most basic   -> ionised plateau
    return float(np.log10(lo / hi))


def main():
    print("=" * 72)
    print("#25 - cost of holding t0 fixed when it is pH-dependent")
    print("=" * 72)
    print(f"\nReported swing: t0 {T0_ACID} -> {T0_BASE} min ({SWING:+.1%}) up to pH 11\n")

    # ---- 1. Why the error concentrates on the ionised plateau ----------
    print("1. Sensitivity of k to a t0 error, by how retained the solute is")
    print("   d(log k) for a -10% error in t0:\n")
    print(f"   {'k':>8}{'tR/t0':>9}{'d(log10 k)':>13}")
    for k in (0.05, 0.1, 0.3, 1.0, 5.0, 20.0):
        t0, tR = 1.0, 1.0 * (1.0 + k)
        k_wrong = (tR - 0.9 * t0) / (0.9 * t0)
        print(f"   {k:>8.2f}{1 + k:>9.2f}{np.log10(k_wrong / k):>13.3f}")
    print("\n   -> the error explodes as k -> 0. An ionised acid at D=1.5 and")
    print("      k_neutral=10 sits at k=0.32, i.e. right in the steep region.")

    # ---- 2. Bias in D across the sweep's actual pH range ---------------
    print("\n2. Bias in fitted D, sweep grid pH 2.5-6.5, t0 measured once at pH 2.5")
    print("   (silanol pKa is not a single number; scanned)\n")
    # The plateau-ratio estimator carries its own small bias: at pH 6.5 an acid
    # of pKa 4.2 is 99.5% ionised, not 100%, so the upper plateau is not quite
    # reached. That artifact is present whatever t0 we use, so it is measured
    # with a per-condition t0 and DIFFERENCED OUT - what remains is the part
    # attributable to t0 alone.
    print(f"   {'pKa_sil':>9}{'t0(2.5)':>9}{'t0(6.5)':>9}{'drift':>8}"
          f"{'D_fixed':>9}{'D_percond':>11}{'t0 bias':>9}")
    ph = np.array([2.5, 3.5, 5.0, 6.5])
    biases = {}
    for pka_sil in (4.5, 6.0, 7.5, 9.0):
        t0s = t0_of_ph(ph, pka_sil)
        kt = k_true(ph, k_neutral=10.0, D=1.5, pka=4.2)
        tR = t0s * (1.0 + kt)
        D_fixed = fit_D_with_assumed_t0(ph, tR, t0_of_ph(2.5, pka_sil), 10.0, 4.2)
        k_pc = (tR - t0s) / t0s
        D_perc = float(np.log10(k_pc[0] / k_pc[-1]))   # per-condition t0
        biases[pka_sil] = D_fixed - D_perc
        drift = (t0s[-1] - t0s[0]) / t0s[0]
        print(f"   {pka_sil:>9.1f}{t0s[0]:>9.3f}{t0s[-1]:>9.3f}{drift:>7.1%}"
              f"{D_fixed:>9.2f}{D_perc:>11.2f}{biases[pka_sil]:>+9.2f}")
    print("\n   (D_percond differs from 1.50 by ~-0.07 for ALL rows - that is the")
    print("    estimator's incomplete-plateau artifact, not a t0 effect, so it")
    print("    cancels in the difference.)")

    # ---- 3. Against the rest of the budget ----------------------------
    print("\n3. Same bias against the errors #15 says we must beat")
    print("   (worst case above, and the sweep's own measurement floor)\n")
    bias = abs(biases[6.0])
    for label, val in (
        ("t0 held fixed (silanol pKa 6.0)", bias),
        ("LSER lack-of-fit floor, ACN (#15 s.5)", 0.030),
        ("LSER lack-of-fit floor, MeOH (#15 s.5)", 0.041),
        ("sweep's min detectable per-run dD (n=8)", 0.0075),
        ("prior SD on D_run (slice.py)", 0.300),
    ):
        print(f"   {label:<42}{val:>8.4f}")
    print(f"\n   -> the t0 bias is {bias / 0.030:.0f}x the LSER lack-of-fit floor and")
    print(f"      {bias / 0.0075:.0f}x the difference the sweep is designed to resolve.")

    # ---- 4. Does a t0 marker in every run actually fix it? -------------
    print("\n4. With a t0 marker in EVERY run (the protocol's mitigation)")
    for pka_sil in (4.5, 6.0):
        print(f"   silanol pKa {pka_sil}: t0-attributable bias falls from "
              f"{biases[pka_sil]:+.2f} to 0.00")
    print("   -> measuring t0 per condition removes the t0-attributable bias.")
    print("      The sweep is safe. The PRODUCT is the exposed case: a user")
    print("      predicting at pH 8 has not injected a t0 marker at pH 8.")


if __name__ == "__main__":
    main()
