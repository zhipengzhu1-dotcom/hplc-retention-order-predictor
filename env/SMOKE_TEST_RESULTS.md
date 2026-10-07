# Ionisation provider environment — smoke test results

Resolves #20 (blocked by, and per the
provider choice made in, #7).

Provider under test: **Uni-pKa** (Apache-2.0), via the `unipka` PyPI wrapper (MIT), version `0.0.6`.
Run with `env/requirements.txt` pinned, on macOS arm64 (Darwin 24.5.0), Python 3.13.7, CPU only (no CUDA).

Run the smoke test with:

```
env/ionisation-venv/bin/python env/smoke_test.py
```

## Method

The provider exposes named getters (`get_acidic_macro_pka`, `get_basic_macro_pka`) **and** a
microspeciation table (`get_distribution(smiles, pH=...)` → SMILES, net charge, population fraction
per microstate). Both were tested.

Predicted pKa was cross-checked independently against the getters by **bisecting the microspecies
population distribution for the pH at which two adjacent net-charge states are equally populated**
(50/50) — the textbook definition of pKa, and a direct exercise of the exact output shape ("SMILES +
charge + fraction") the repo's provider interface (from #7) requires. This is the number reported
below as "predicted pKa".

**Caveat found on the named getters** (see Notes): `get_acidic_macro_pka` /
`get_basic_macro_pka` are defined relative to the **neutral (net-charge-0) macrostate transition**,
not per functional group. On a zwitterion this makes them look "swapped" versus naive
acid-group/base-group intuition — e.g. for alanine, `get_acidic_macro_pka` returns 9.67 (the
amine-related, high-pH crossing) and `get_basic_macro_pka` returns 2.30 (the carboxyl-related,
low-pH crossing). The population-crossover method sidesteps this ambiguity and is what should be
used downstream, not the getter labels at face value.

## Results

| Compound | Type | Site | Predicted pKa (population crossover) | Published pKa | Error (pred − pub) | Source |
|---|---|---|---:|---:|---:|---|
| Benzoic acid | monoprotic acid | –COOH | 4.04 | 4.20 | **−0.16** | Haynes, *CRC Handbook of Chemistry and Physics*, 97th ed. (2016) |
| Propranolol | monoprotic base | secondary amine | 9.48 | 9.5 | **−0.02** | Avdeef, *Absorption and Drug Development*, 2nd ed. (2012) |
| Alanine | zwitterion (amino acid) | –COOH (pKa1) | 2.31 | 2.34 | **−0.03** | CRC Handbook, 97th ed. (2016) |
| Alanine | zwitterion (amino acid) | –NH3+ (pKa2) | 9.66 | 9.69 | **−0.03** | CRC Handbook, 97th ed. (2016) |
| Ciprofloxacin | multiprotic drug | –COOH (pKa1) | 5.58 | 6.09 | **−0.51** | Ross & Riley, *Int. J. Pharm.* 63 (1990) 237–250 |
| Ciprofloxacin | multiprotic drug | piperazine N4 (pKa2) | 8.99 | 8.62 | **+0.37** | Ross & Riley, *Int. J. Pharm.* 63 (1990) 237–250 |

Mean absolute error: **0.19 pKa units** across the 4 simple sites (benzoic acid, propranolol,
alanine×2); **0.28 pKa units** including both ciprofloxacin sites. All errors are inside the
"open tools are now comparable to commercial predictors" finding cited in #7's benchmark
(10.1021/acs.jcim.6c00107).

Applying #7's own sensitivity result (`worst-case |Δ log k| ≈ min(c(D)·δ, D)`, `c(1.5)≈0.7`), a
δ ≈ 0.5 error near a crossover costs roughly 0.3–0.4 log k units — non-trivial but well inside the
"open predictors are usable" conclusion #7 already reached.

## Notes / install friction

- **First call downloads model weights** (`t_dwar_v_novartis_a_b.pt`, from a GitHub release URL) to
  a local cache on first instantiation of `unipka.UnipKa()`. This took ~24s once on this machine's
  connection; subsequent instantiations are near-instant (~0.4s) since the weights are cached. This
  network dependency at first run should be called out for any CI/offline deployment — pre-warm the
  cache or vendor the weights file.
- Install itself was clean: `pip install rdkit unipka` on Python 3.13 pulled only pure-Python/CPU
  wheels (`torch==2.13.0` CPU build, `rdkit==2026.3.5`, `pandas==3.0.5`, plus plotting/notebook
  extras `unipka` brings in for `draw_distribution`). No MATLAB Runtime, no CUDA compile step, no
  license server — matching #7's characterisation.
- **Getter-vs-microspecies discrepancy (see Method above)**: `get_acidic_macro_pka` /
  `get_basic_macro_pka` track the neutral-macrostate transition, not "the acid group's pKa" /
  "the base group's pKa" as the names suggest. Anyone consuming this provider should read pKa off
  the microspecies distribution crossover (or population table), not the getter return values,
  when the compound is a zwitterion or multiprotic.
- **Ciprofloxacin has a spurious extra crossing**: scanning the full pH range surfaces a transition
  between charge +3 and +1 microstates around pH ≈ 3.4–3.5 that does not correspond to either
  published pKa (6.09, 8.62). It reflects a doubly-protonated-piperazine microstate that the model
  gives non-trivial population to at low pH; it doesn't perturb the two crossings that do match
  literature, but a naive full-range scan without bracketing (as done in `env/smoke_test.py`) will
  wrongly report it as "pKa1" — flagged here so it isn't rediscovered as a surprise later.
- Confirms **output shape supports the repo's provider contract** (#7): `get_distribution` returns
  one row per microspecies with `smiles`, `charge`, `population` (fraction) columns — directly
  usable as the required `species[(smiles, charge, fraction, ...)]` tuple, modulo adding
  `fraction_sigma` (not emitted — see Uncertainty gap below).
- **No calibrated per-prediction uncertainty (σ) is emitted**, confirming #7's central negative
  finding. `get_distribution` gives a point-estimate population per microstate only; no ensemble
  spread or confidence interval ships with this wrapper version (0.0.6). The repo's
  ensemble-spread / cross-provider-disagreement / conformal-calibration layer from #7 remains a
  build-it-ourselves item.
- All predictions above were computed **in water** (the provider's only mode); the Rosés/Bosch
  `w_wpKa → s_spKa` correction from #7 must still be applied downstream before use on an
  ACN/MeOH mobile phase — this ticket does not implement that correction, only confirms the
  water-only pKa/microspecies numbers it will need to correct.
