# The scenario ensemble — interface specification

Status: **v1, amending #14 with the thin slice's punch list (#44).**

This is the prediction layer's output interface and the concrete meaning of "a ranked
retention distribution with calibrated confidence". Until now it existed only as prose in
`CONTEXT.md` and as an ad-hoc `.npz` in `prototype/thin-slice/`. #16/#17 code will be written
against it, so it is written down and made checkable: [`ensemble.py`](ensemble.py) validates
an object against every MUST below, and refuses one that does not conform.

## The core object

```
retention[P, S, n]     float   minutes, NaN for refused compounds
width[P, S, n]         float   peak standard deviation in minutes, NaN where retention is
f_neutral[P, S, n]     float   population-weighted neutral fraction, NaN for refused
compounds[n]           str     stable identifiers
method_card[P]         object  the conditions each operating point was computed at
stratum[n]             str     conformal stratum label, method-independent (#42)
```

`S` is the scenario axis. **Each scenario is a complete alternative chromatogram** — one draw
of the per-run parameters plus one correlated draw across every compound. Per-compound
marginals are a derived display view and a **lossy projection**: they destroy the correlation
that elution order depends on, so they are never the stored form.

### Why `P` is an axis and not a loop (punch-list item 1)

The thin slice needed three `(φ, pH)` points and found `retention[S, n]` could not express
them. Looping outside the object was possible but wrong: **a scenario ensemble whose
conditions are not part of the object invites silent condition mismatches**, and the per-run
parameter draw is shared *across* operating points within a scenario, which a loop would
resample and thereby destroy.

**MUST**: `retention`, `width` and `f_neutral` are all `[P, S, n]`. **MUST**: scenario `s` at
operating point `p` and scenario `s` at operating point `q` share the same per-run draw.

## The method card (punch-list item 1, continued)

**MUST travel with the array, not beside it.**

| field | type | notes |
|---|---|---|
| `phi` | float | organic volume fraction, 0–1 |
| `modifier` | enum | `MEOH` \| `ACN` \| `THF` |
| `ph` | float | the working pH |
| `ph_scale` | enum | **`SS` \| `WW`, required, no default** |
| `temperature_c` | float | |
| `column_id` | str | resolves to a column vector; identity, not a name |
| `t0_min` | float | **per operating point — see below** |
| `t0_source` | enum | `MEASURED` \| `FITTED` \| `NOMINAL` |
| `flow_ml_min` | float | **required, no default** — see below |
| `instrument_id` | str \| null | resolves to dwell / extra-column volumes |
| `gradient` | object \| null | null means isocratic |

`flow_ml_min` has no default either, and for the same class of reason: **plate count is a
property of the method, not of the column.** Knox gives `N = L/(h·dp)` with
`h = A·ν^(1/3) + B/ν + C·ν` and `ν = u·dp/Dm`, so `N` moves with flow. A defaulted flow
silently fixes the peak width of every prediction, which is `width` reverting to the
placeholder this spec exists to remove. Found by implementing `model/dispersion.py` — the
field was simply absent from v1.

`ph_scale` has no default deliberately. #20 carries the `w_w`→`s_s` correction as a named
error term, and a defaulted scale is exactly how that error becomes silent.

### `t0` is per operating point, not a scalar — a fifth amendment

The thin slice stored `t0` as **one scalar for the whole object**. #25 established today that
this is unsound: hold-up time falls as residual silanols ionise, and holding it fixed biases
`D` by **+0.42 to +0.69** log units over pH 2.5–6.5 — larger than the entire prior SD on
`D_run`. Because `k = (t_R − t0)/t0` is most sensitive where `k` is smallest, the error lands
on the ionised plateau, which is the plateau that defines `D`.

**MUST**: `t0_min` is a field of `method_card[p]`, so it varies across the `P` axis.
**MUST**: `t0_source` is carried, because a `NOMINAL` `t0` obliges the consumer to add the
per-run `t0` error term and flag the prediction as extrapolation across pH.

## `f_neutral` and the stratum (punch-list item 2)

The slice computed per-scenario microspecies weights inside the sampler **and then discarded
them**. Two consumers need them and neither can reconstruct them:

- **`D` refinement from scouting runs** (#34) — needs the neutral fraction that applied in the
  scenario being reconciled against.
- **Conformal stratum assignment** (#15 §4) — ionisation regime is the one stratum axis that
  **never merges**, so it must be recoverable per compound.

**MUST**: `f_neutral[P, S, n]` is stored, not recomputed. It is per-scenario because pKa is
drawn per scenario, so the neutral fraction differs across `S` at the same operating point.

**MUST**: `stratum[n]` is cut on the **method-independent structural summary** (#42), never on
the LSER-propagated spread — so a compound cannot change calibration class because the user
changed column. It is therefore indexed by `n` alone, not by `[P, n]`, and that shape is the
enforcement.

Where the provider gives microspecies identity, carry it alongside; it is optional because
providers differ, whereas `f_neutral` is available from all of them.

## `width` owes physics (punch-list item 3)

The slice filled `width` from a nominal `N = 10000` placeholder. The slot is right; the
content is absent. **`width` MUST be composed, not asserted:**

```
width² = column_variance(N, t_R) + extra_column_variance(instrument) + gradient_compression
```

- **Plate count `N`** ships as a **Knox engineering prior by particle class** at ±10–30%
  (#9, #13). **MUST** be sampled per scenario, not fixed — the width slot carries the prior's
  uncertainty, not its point estimate. `CONTEXT.md` is explicit that no improvement in
  retention prediction removes this floor; it narrows only when measured van Deemter data
  enters through reconciliation.
- **Extra-column variance** comes from the instrument, and `sources/instrument-volumes/` shows
  it is not a constant: **5 µL on a Waters Acquity binary against 30–34 µL on an Alliance
  2695**, a factor of six, with flow-cell volume spanning 0.5–13 µL. It is a **per-run** term
  and resolves from `method_card[p].instrument_id`.
- **Gradient compression** applies only where `gradient` is non-null.

**MUST**: where `instrument_id` is null, extra-column variance is carried as a *distribution*
spanning the observed vendor range, not silently set to zero. Zero is the one value it
certainly is not.

## Refusal semantics

**MUST**: refused compounds keep their column and index, with `NaN` retention. `n` never
shrinks and the compound list stays intact. They carry **no fabricated distribution** and are
candidate explanations for unexplained peaks only — never competing on retention against a
compound that has a real distribution.

Verified in the slice: benzyltrimethylammonium (no neutral microspecies at any pH) held its
NaN column through every analysis, with nothing downstream special-casing it beyond
NaN-awareness. This part of #14 needs no amendment.

## Effective sample size (punch-list item 4)

**MUST** be emitted as a standing output whenever reconciliation consumes a measurement; `n`
never shrinks and ESS is per-request.

⚠ **Unverified.** No measurement entered the thin slice, so no reweighting happened and the
ESS contract has never executed. Nothing contradicts the design; the slice simply cannot
confirm it. **First verification belongs to #45**, on the first real data.

## What this does not specify

- **The scenario count `S`.** #14/#17 own it; 1,000 was the slice's choice, not a contract.
- **How strata are cut**, beyond that the cut is method-independent — #42 and #15 §4 own that.
- **Serialisation format.** The slice used `.npz`; nothing here depends on it. `ensemble.py`
  validates the in-memory object, so any format that round-trips it conforms.
