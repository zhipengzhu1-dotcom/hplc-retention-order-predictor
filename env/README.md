# Ionisation provider environment

> `env/` now holds two provider environments. This file is the **ionisation** one (#20).
> The **QSPR / Abraham-descriptor** one is
> [`QSPR_ENVIRONMENT.md`](./QSPR_ENVIRONMENT.md) (#38) — SoluteML with its ensemble
> covariance exposed, in `env/qspr-venv` with its own lockfile.

Resolves #20 — standing up the
provider environment chosen in #7.

## Setup

```
python3 -m venv env/ionisation-venv
env/ionisation-venv/bin/pip install -r env/requirements.txt
```

Exact resolved versions are pinned in [`requirements.txt`](./requirements.txt) (`rdkit==2026.3.5`,
`unipka==0.0.6`, plus transitive deps including `torch==2.13.0` CPU). Resolved on macOS arm64
(Darwin 24.5.0), Python 3.13.7, 14 Aug 2026.

The venv directory itself (`env/ionisation-venv/`) is not committed — recreate it from the
lockfile. `unipka.UnipKa()` downloads model weights (`t_dwar_v_novartis_a_b.pt`) from a GitHub
release URL on first instantiation (~24s once, then cached; no further network needed).

## Contents

- [`requirements.txt`](./requirements.txt) — exact-version lockfile (`pip freeze` output).
- [`smoke_test.py`](./smoke_test.py) — runnable smoke test on 4 compound classes (monoprotic acid,
  monoprotic base, zwitterion, multiprotic drug), predicted vs. published pKa.
- [`SMOKE_TEST_RESULTS.md`](./SMOKE_TEST_RESULTS.md) — the results table, method, and install
  friction found while running it.
- [`CHEMAXON_LICENCE.md`](./CHEMAXON_LICENCE.md) — evaluation of the ChemAxon/Certara academic
  licence: what's covered, what's unresolved, and the exact human steps to apply (not performed).

## Provider interface shape — does it fit the repo's contract?

The repo's convention (pinned by #7): every provider is water-only, returns **microspecies**
(SMILES + charge + fraction), and callers apply the Rosés/Bosch `w_wpKa → s_spKa` correction
before standardising on the `s_spH` scale.

**Uni-pKa (via `unipka`) fits the required shape on the parts it owns:**

- `UnipKa.get_distribution(smiles, pH=...)` returns a table with one row per microspecies, columns
  `smiles`, `charge`, `population` (fraction) — directly usable as the interface's
  `species[(smiles, charge, fraction, ...)]`.
- It is water-only (no solvent/modifier parameter exists anywhere in the API) — matches the
  contract's premise that the correction layer is a separate step, not the provider's job.
- Per-site pKa is recoverable (via the microspecies population crossover — see
  `SMOKE_TEST_RESULTS.md` for why the named getters `get_acidic_macro_pka` /
  `get_basic_macro_pka` should *not* be read at face value on zwitterions/multiprotic compounds).

**What it does *not* provide, and the repo's adapter will have to supply:**

- **No `fraction_sigma` / no per-prediction σ of any kind.** Confirms #7's central negative
  finding — the uncertainty layer is not in this provider and must be built on top (ensemble
  spread, cross-provider disagreement, conformal calibration on downstream log k).
- **No charge-type / compound-family classification per site**, which the Rosés/Bosch correction
  needs to pick coefficients. The SMARTS-based classifier #7 specified is still a to-build item.
- **No `in_domain` / `domain_score`.** OPERA remains the intended domain-gate provider per #7; not
  substituted here.
- **No `s_spH` awareness** — `get_distribution(pH=...)` takes a bare float; the caller is fully
  responsible for knowing that float is `w_wpH` in water and for the scale-tagging discipline #7
  specifies ("every pH crossing a module boundary carries its scale as a tag").

Net: the microspecies **shape** is compatible out of the box; the **σ, family-label, and pH-scale
tagging** are adapter responsibilities that remain unbuilt, exactly as #7 anticipated when it
specified the interface ahead of any provider being wired up.
