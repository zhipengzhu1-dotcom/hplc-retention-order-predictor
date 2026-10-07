"""
Smoke test for the open-default ionisation provider (Uni-pKa via the `unipka` PyPI
wrapper), chosen in issue #7 and stood up in issue #20.

Runs on four known compounds spanning the required ionisation types:
  - benzoic acid    (monoprotic acid)
  - propranolol     (monoprotic base)
  - alanine         (zwitterion, amino acid)
  - ciprofloxacin   (multiprotic / zwitterionic drug)

For each compound:
  1. Calls the provider's macro-pKa getters (`get_acidic_macro_pka`,
     `get_basic_macro_pka`) directly.
  2. Independently derives pKa from the microspecies population crossover
     (the pH at which two adjacent net-charge populations are equal — the
     textbook definition of pKa) via `get_distribution(smiles, pH=...)`,
     as a cross-check on (1). Uni-pKa's own getter *labels* did not line up
     with the population crossover on the zwitterion in initial testing —
     see README note.
  3. Prints the microspecies table itself (SMILES + charge + fraction) to
     confirm the output shape required by the repo's provider contract
     (#7's interface: `species[(smiles, charge, fraction, ...)]`).

Run with:
    env/ionisation-venv/bin/python env/smoke_test.py
"""

from __future__ import annotations

import time

import unipka

COMPOUNDS = {
    # name: (SMILES, published pKa values [list], source)
    "benzoic_acid": {
        "smiles": "c1ccccc1C(=O)O",
        "published": [4.20],
        "source": "Haynes, CRC Handbook of Chemistry and Physics, 97th ed. (2016), benzoic acid pKa 4.202 (25C, water)",
    },
    "propranolol": {
        "smiles": "CC(C)NCC(O)COc1cccc2ccccc12",
        "published": [9.5],
        "source": "Avdeef, Absorption and Drug Development, 2nd ed. (2012), propranolol pKa(amine) 9.5 (25C, water)",
    },
    "alanine": {
        "smiles": "CC(N)C(=O)O",
        "published": [2.34, 9.69],
        "source": "CRC Handbook of Chemistry and Physics, 97th ed. (2016), alanine pKa1(COOH) 2.34, pKa2(NH3+) 9.69",
    },
    "ciprofloxacin": {
        "smiles": "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
        "published": [6.09, 8.62],
        "source": "Ross & Riley, Int. J. Pharm. 63 (1990) 237-250, ciprofloxacin pKa1(COOH) 6.09, pKa2(piperazine N) 8.62",
        # Ciprofloxacin has a 4th populated microstate charge (+3/+2, doubly-protonated
        # piperazine + protonated carboxyl) very close in pH to the +1/0 transition, which
        # a wide +-3 bracket wrongly latches onto. Narrow the search window per-pKa so the
        # bisection converges on the transition adjacent to the published value.
        "brackets": [(4.5, 7.0), (7.5, 9.5)],
    },
}


def population_crossover(u: "unipka.UnipKa", smiles: str, ph_lo: float, ph_hi: float, tol: float = 0.02) -> float:
    """Bisect for the pH at which the two dominant adjacent charge states
    (by population, evaluated at the bracket midpoint's neighbourhood) are
    equally populated. This is the textbook definition of pKa and is used
    here as an independent check on the provider's named getters."""

    def charge_pops(ph: float) -> dict[int, float]:
        df = u.get_distribution(smiles, pH=ph)
        return df.groupby("charge")["population"].sum().to_dict()

    lo, hi = ph_lo, ph_hi
    pops_lo = charge_pops(lo)
    dom_lo = max(pops_lo, key=pops_lo.get)
    for _ in range(30):
        mid = (lo + hi) / 2
        pops_mid = charge_pops(mid)
        dom_mid = max(pops_mid, key=pops_mid.get)
        if dom_mid == dom_lo:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return round((lo + hi) / 2, 2)


def main() -> None:
    print("Loading Uni-pKa model (first call downloads weights to ~/.cache if not present)...")
    t0 = time.time()
    u = unipka.UnipKa()
    print(f"  model ready in {time.time() - t0:.1f}s\n")

    print(f"unipka version: {unipka.__version__}\n")

    rows = []
    for name, spec in COMPOUNDS.items():
        smiles = spec["smiles"]
        print(f"=== {name} ({smiles}) ===")

        acidic = u.get_acidic_macro_pka(smiles)
        basic = u.get_basic_macro_pka(smiles)
        print(f"  get_acidic_macro_pka -> {acidic:.2f}")
        print(f"  get_basic_macro_pka  -> {basic:.2f}")

        # Microspecies table at physiological pH, to confirm output shape
        df = u.get_distribution(smiles, pH=7.4)
        shown = df[["smiles", "charge", "population"]].sort_values("population", ascending=False)
        print("  microspecies @ pH 7.4 (smiles, charge, fraction):")
        for _, r in shown.iterrows():
            print(f"    {r['smiles']:<35s} charge={int(r['charge']):+d}  fraction={r['population']:.4f}")

        # Population-crossover cross-check, bracketing each published value.
        # Default bracket is +-3 pH units; compounds with closely-spaced or extra
        # microstate transitions (see `brackets` override) use a narrower, manually
        # chosen window so the bisection converges on the transition adjacent to the
        # published pKa rather than an unrelated nearby crossing.
        crossovers = []
        default_brackets = [(p - 3.0, p + 3.0) for p in spec["published"]]
        brackets = spec.get("brackets", default_brackets)
        for (lo, hi) in brackets:
            xc = population_crossover(u, smiles, lo, hi)
            crossovers.append(xc)
        print(f"  population-crossover pKa(s): {crossovers}")
        print(f"  published pKa(s):            {spec['published']}  [{spec['source']}]")

        errors = []
        for pub, pred in zip(sorted(spec["published"]), sorted(crossovers)):
            err = round(pred - pub, 2)
            errors.append(err)
        print(f"  error (predicted - published): {errors}\n")

        rows.append(
            {
                "compound": name,
                "published": spec["published"],
                "crossover_pka": crossovers,
                "error": errors,
                "acidic_getter": round(acidic, 2),
                "basic_getter": round(basic, 2),
            }
        )

    print("\n=== Summary ===")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
