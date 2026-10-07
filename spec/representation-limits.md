# Known limits of the Abraham representation

Status: **named blind spot, deliberately not operationalised** (#46, 2026-08-16).

This document exists so that a reader of the spec cannot mistake the model's confidence for
completeness. It names a class of error the architecture **cannot currently detect**, states
why no automatic flag is shipped for it, and records what would change that.

## The two kinds of doubt, and why only one is spendable

#42 established that doubt has two sources with opposite remedies:

| | remedy | detectable from the ensemble? |
|---|---|---|
| **Predictor uncertainty** — the ensemble is unsure | **Spendable**: more data, better architecture, or a measured descriptor for that compound | yes, that is what its spread is |
| **Structural deficiency** — the Abraham representation cannot encode what the molecule does | **Not spendable**: a better predictor cannot help, because the five numbers are the wrong five numbers | **no, by construction** |

The second is invisible to the ensemble because its defining property is that the networks
*agree* — having learned the same inadequate representation from the same corpus. A confident,
unanimous, wrong answer is what structural deficiency looks like from inside.

## The register, as named by its sources

None of these is our observation; each was named independently and *a priori*:

| class | source |
|---|---|
| Intramolecular hydrogen bonding | Ulrich 2026, conclusions — "not accounted for in the approach" |
| Site-specific or sterically hindered H-bond donation/acceptance | Ulrich 2026 — `A`/`B` average over sites |
| Dipolarity/polarisability mixtures | Ulrich 2026 — several interactions collapse into `S` |
| Bulky compounds | WSU-2019 Table S-1 — dropped from WSU's own fits |
| Some hydrogen-bond acids | WSU-2019 Table S-1 — ditto |

The last two are the same phenomenon seen from the other end of the chain: compounds the
*fitters* excluded, not compounds the *predictor* missed.

## Why there is no automatic flag

Tested against the only instrument available — WSU-94 measured descriptors, SoluteML's 25-fold
predictions, and SMARTS classification (`experiments/structural-deficiency/`). The test used
#42's own claim as its shape: a structurally deficient class should show **elevated residual
with flat ensemble spread**; both elevated is ordinary predictor uncertainty.

| class | n | resid ÷ spread (`S`, `A`) | outcome |
|---|---|---|---|
| intramolecular H-bond | **2** | — | untestable |
| hindered H-bond donor | **1** | — | untestable |
| H-bond acid (any donor) | 42 | 1.30, 1.11 | ensemble already signals it |
| bulky (`V` > 1.4 proxy) | 10 | 0.73, 1.39 | no elevation |
| **dipolarity mixture** | 8 | **1.88, 2.14** | signature present |

**Three reasons no flag ships:**

1. **Alert fatigue.** The naive operationalisation flags **52 of 93 compounds (56%)**, driven
   almost entirely by "has any O–H or N–H". A caveat on more than half of everything is
   ignored the first time it appears on a routine molecule, and then it protects nothing.
2. **`n` = 8 is not a rule.** The dipolarity signal is the only one with the right shape, and
   it rests on eight compounds in a set that is very likely memorised. Carving it into code
   would make an intriguing lead into a load-bearing discriminator on the strength of a
   contaminated sample.
3. **Silence is the worse failure.** Retiring the class entirely does not remove the doubt; it
   makes it silent, and silent unspendable doubt gets misread as spendable — which is how
   money is spent on a better predictor that cannot help. Naming it without automating it
   costs nothing and claims nothing.

## ⚠ The contamination asymmetry, which governs how the table above may be used

**All 94** WSU compounds are inside SoluteDB — measured row by row against the CC-BY
corpus in #35 (93 on SoluteDB's own fixed-H InChI key, 1 on InChIKey skeleton), no longer
inferred from provenance. These residuals are therefore in-sample and **too small**, and
the clean re-run once handed to #35 is **not executable**: no corpus examined pairs
out-of-corpus compounds with a usable label (`experiments/solutedb-census/`). That bias has a direction, and it is not symmetric:

- **A positive result survived a headwind.** The dipolarity signature appears *despite* partial
  memorisation, so it is more credible than `n` = 8 alone suggests.
- **A null result is weak.** "Bulky shows no elevation" is equally consistent with the class
  being real and memorised away.

**Therefore: this test may support a class. It may not retire one.** No class above has been
disproved, and none should be described as such.

## What an implementer must do

- **Do not** flag compounds automatically against this register.
- **Do not** describe the ensemble's confidence as complete. A compound may carry unspendable
  error with a narrow spread, and nothing in the current stack will say so.
- **Do** state, wherever the model's uncertainty is presented, that it covers predictor
  uncertainty and not representational adequacy. That is a one-line disclosure, not a flag.

## What would change this

An **un-memorised test bed**. The blocker is not analysis, it is that every compound we can
score against was plausibly in the training corpus.

**Poole 2020** (*J. Chromatogr. A* 1617, 460841 — the several-hundred-solute descriptor set)
is the concrete route, tracked on **#35**. With residuals measured on compounds SoluteML has
not seen, the same test becomes decisive rather than suggestive, and two hypotheses are worth
revisiting first:

1. **Dipolarity mixtures** — does the 1.9–2.1× under-signal survive out of sample?
2. **Intramolecular H-bonding** — untestable here at `n` = 2, and the class Ulrich named most
   explicitly. It needs a test bed containing salicylates and *ortho*-substituted phenols in
   quantity.

If either survives, it becomes a real discriminator and this document is replaced by one that
ships a flag. Until then the honest position is: **the ceiling of the five-parameter Abraham
space is real, we know roughly where it is, and we cannot yet point at it per compound.**
