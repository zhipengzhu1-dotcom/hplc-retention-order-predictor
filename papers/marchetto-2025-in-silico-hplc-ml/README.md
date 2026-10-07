# Marchetto 2025 — in silico HPLC method development (the blueprint)

> A. Marchetto, F. Tirapelle, L. Mazzei, E. Sorensen, M. O. Besenhard, *In Silico
> High-Performance Liquid Chromatography Method Development via Machine Learning*,
> Anal. Chem. **97** (2025) 6991–7001.
> <https://doi.org/10.1021/acs.analchem.4c03466> — Open Access, CC-BY 4.0.

| File | Role |
|---|---|
| `ac4c03466.pdf` | main text |

The blueprint for this effort, and the source of #2's error budget. Chains three layers:

```
SMILES → molecular descriptors → QSPR → Abraham solute descriptors (E, S, A, B, V)
       → LSER → LSS theory → k(φ)
```

⚠ Its `S_S` is **2.303× the conventional `S`** (`ln` vs `log₁₀`) — see `CONTEXT.md`.

**Supplementary material is not held locally.** SI-1 (PDF) and SI-2 (XLSX, retention data)
are at <https://pubs.acs.org/doi/suppl/10.1021/acs.analchem.4c03466/suppl_file/ac4c03466_si_001.pdf>
and `…_si_002.xlsx`; both were read via a research agent for #2.
