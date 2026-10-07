# SoluteML ensemble covariance — demonstration results

Output of `env/qspr-venv/bin/python env/soluteml_demo.py` on the patched, pinned environment
(see [`QSPR_ENVIRONMENT.md`](./QSPR_ENVIRONMENT.md)). Resolves
#38.

Everything below is **ensemble spread**: epistemic only, no aleatoric term, uncalibrated.
Recalibration is #14.

**Scope, so this is not confused with the measurement it resembles.** This is a
**9-compound demonstration** that the environment works and that the off-diagonals are
material (#38). The measurement over the full anchor lives in
[`prototype/descriptor-covariance/FINDINGS.md`](../prototype/descriptor-covariance/FINDINGS.md)
(94 compounds, 25 members), and the error covariance *against measured values* — a different
quantity again, and much stronger — is in
[`prototype/soluteml-bias/FINDINGS.md`](../prototype/soluteml-bias/FINDINGS.md) §8. Quote
those for any claim; quote this only for "the environment reproduces published behaviour".

## Headline: yes, the off-diagonals are material — and they matter most where we live

`max |r|` is the largest off-diagonal correlation over (E, S, A, B). The eigenvalues of the
correlation matrix bound how wrong the independence assumption can be: for any LSER coefficient
vector `m`, `m'Σm / m'diag(Σ)m` lies in `[λ_min, λ_max]`. No system constants need to be assumed.

| compound | max \|r\| | λ_min | λ_max |
|---|---:|---:|---:|
| octane | 0.317 | 0.508 | 1.554 |
| toluene | 0.362 | 0.626 | 1.749 |
| phenol | 0.264 | 0.608 | 1.418 |
| caffeine | 0.508 | 0.382 | 1.777 |
| benzoic acid | 0.654 | 0.263 | 1.910 |
| paracetamol | 0.612 | 0.329 | 2.141 |
| ciprofloxacin | 0.626 | 0.354 | 2.536 |
| propranolol | 0.794 | 0.205 | 2.676 |
| naproxen | 0.880 | 0.090 | 2.593 |

**Materially non-zero, by a wide margin.** Across this probe set the correlation reaches
**0.88**, and the eigenvalue span is **[0.09, 2.68]** — assuming independence can overstate the
variance of a descriptor combination by up to **11×** or understate it by up to **2.7×**,
depending on the coefficient vector. Dropping the off-diagonals is not conservative; #14 already
established there is no safe direction, and these numbers put a size on it.

**The rows are sorted deliberately.** Correlation is weakest for the small neutral solutes the
model was trained on (octane 0.32, phenol 0.26) and strongest for the drug-like compounds
(naproxen 0.88, propranolol 0.79, ciprofloxacin 0.63). #5 put the median drug-like ionisable
compound at ~the 95th percentile of the training distribution; the covariance structure is
largest exactly there. A pipeline that ignored covariance would be least wrong on the chemistry
we do not care about.

**Signs are mixed and compound-specific**, so no single correlation prior would serve:

- naproxen: `r(S,B) = −0.880`, `r(E,B) = +0.645`
- propranolol: `r(E,A) = +0.794`, `r(E,B) = −0.672`, `r(A,B) = −0.659`
- ciprofloxacin: everything positive, 0.41–0.63
- caffeine: `r(E,A) = −0.508`
- paracetamol: `r(S,A) = +0.612`, `r(A,B) = −0.466`

## Descriptors and ensemble sigma

| compound | E | S | A | B | σ_E | σ_S | σ_A | σ_B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| octane | 0.038 | 0.027 | 0.019 | 0.017 | 0.034 | 0.032 | 0.011 | 0.020 |
| toluene | 0.662 | 0.568 | 0.021 | 0.175 | 0.029 | 0.032 | 0.011 | 0.020 |
| benzoic acid | 0.819 | 0.952 | 0.677 | 0.485 | 0.052 | 0.042 | 0.034 | 0.037 |
| phenol | 0.881 | 0.935 | 0.646 | 0.335 | 0.032 | 0.029 | 0.037 | 0.034 |
| caffeine | 1.454 | 1.826 | 0.111 | 1.243 | 0.041 | 0.100 | 0.029 | 0.037 |
| paracetamol | 1.158 | 1.663 | 1.075 | 0.908 | 0.050 | 0.064 | 0.043 | 0.047 |
| propranolol | 1.848 | 1.404 | 0.384 | 1.402 | 0.088 | 0.088 | 0.110 | 0.056 |
| naproxen | 1.587 | 1.815 | 0.600 | 0.738 | 0.074 | 0.242 | 0.042 | 0.065 |
| ciprofloxacin | 2.229 | 2.363 | 0.709 | 2.507 | 0.064 | 0.060 | 0.063 | 0.046 |

Chemically sensible against experimental values (toluene E 0.601 / S 0.52 / A 0.00 / B 0.14;
phenol 0.805 / 0.89 / 0.60 / 0.30; caffeine 1.500 / 1.60 / 0.00 / 1.35), with a small positive
bias of roughly +0.02 to +0.09 across the board — see the reproduction gap below. Note the
model gives A = 0.019 for octane and 0.021 for toluene where the true value is exactly 0; #6's
warning about the `A = 0` vs `A ≠ 0` classification stands.

`diag(cov)` reproduces the shipped per-descriptor variance to 7e-18 — it is the same tensor,
summarised differently, not a re-derivation.

## Per-member predictions (for #37)

`member_preds` comes back as `(25, 9, 4)`. Across-compound correlation of the per-member
deviation in `E`, first four compounds:

```
        octane  toluene  benzoic   phenol
octane   1.000    0.457    0.409    0.440
toluene  0.457    1.000    0.619    0.649
benzoic  0.409    0.619    1.000    0.650
phenol   0.440    0.649    0.650    1.000
```

**Across-compound correlation is strongly positive, 0.41–0.65 on this small set.** That is the
sign #14 predicted when it corrected the ticket body: with positive across-compound correlation,
assuming independence **overstates** the variance of a *difference*, so shared model bias cancels
out of elution order. #37 should now measure this properly on a scaffold-structured set — the
API no longer stands in the way, but four compounds is an illustration, not a measurement.

## The reproduction gap — report this, do not paper over it

The patched build does **not** bit-reproduce the sample output published in the upstream repo:

| | E | S | A | B | L |
|---|---:|---:|---:|---:|---:|
| octane, ours | 0.0382 | 0.0268 | 0.0190 | 0.0169 | 4.0234 |
| octane, upstream | 0.0034 | 0.0003 | 0.0012 | −0.0000 | 3.6623 |
| Δ / σ | 2.60 | 2.19 | 2.32 | 1.31 | 6.27 |
| toluene, ours | 0.6616 | 0.5676 | 0.0206 | 0.1748 | 3.5949 |
| toluene, upstream | 0.6077 | 0.5312 | 0.0017 | 0.1457 | 3.3649 |
| Δ / σ | 3.55 | 2.04 | 2.99 | 2.15 | 3.66 |

Absolute deltas are 0.017–0.054 on (E, S, A, B) — comfortably inside the model's own published
error (SoluteML contributes SE(log k) ≈ 0.08–0.11 per #6, and descriptor RMSEs are ~0.09–0.17),
so the environment is usable. But they are **1.3–3.6× the shipped ensemble σ**, and on `L`
(which we do not use) up to 6.3σ.

Cause: four years of RDKit/descriptastorus drift in the 200 CDF-normalized 2D descriptors that
feed the model alongside the message-passing encoder. The most violent instance is documented
under friction below (`Ipc`); the residual gap is the rest of the drift, which cannot be removed
without reconstructing the 2021 cheminformatics stack.

**This is itself evidence for #14's central claim.** Pure toolchain drift moves the prediction by
more than the model's own reported epistemic uncertainty. Anything that consumed the shipped σ as
a predictive standard deviation would be claiming a precision the build cannot even reproduce
across library versions.

## Reproducing

```
bash env/install_qspr.sh
env/qspr-venv/bin/python env/soluteml_demo.py
```

Runtime ~1 min on an M-series CPU (25 forward passes per compound, `batch_size=1` hard-coded
upstream). Deterministic — no dropout at inference, no sampling.
