# WSU-2019 — extraction of the remainder of the paper (journal pp. 118–126)

**Source:** Poole, C.F. "Reversed-phase liquid chromatography system constant database over an
extended mobile phase composition range for 25 siloxane-bonded silica-based columns."
*J. Chromatogr. A* **1600** (2019) 112–126. DOI
[10.1016/j.chroma.2019.04.027](https://doi.org/10.1016/j.chroma.2019.04.027).
Read from a locally-held PDF (gitignored; not committed).

Companion to `README.md` in this directory, which documents the earlier open-access surrogate
work. This note covers what the main text itself contains beyond pages 1–8 (already summarised on
issue #22).

## What the remaining pages actually are

The PDF is **15 pages**, not 18 — journal pages 112–126 map to PDF pages 1–15. A request for
"pages 9–18" returns 7 pages (PDF 9–15 = journal 120–126). Pages 122–126 are §3.4, §3.5,
Conclusions, Appendix A and 4 pages of references. So the paper's main-text numeric content ends
at journal page 121.

Only **two numeric tables** exist in the main text after page 8: Table 2 (journal p. 118) and
Table 3 (pp. 120–121). Both are fully extracted. Everything else numeric in the paper lives in
figures or in supplementary Tables S-1…S-5.

## Tables extracted

### `wsu2019-table2-kinetex-c18-vs-xb-c18-deviations.csv`
Table 2, journal p. 118. Average absolute deviation (AAD) and average deviation (AD) for the
difference in system constants between **Kinetex C18 (column 10)** and **Kinetex XB-C18
(column 12)** over 10–70% (v/v) methanol–water, plus the observed range of each constant.
All 5 rows (e, s, a, b, v) captured. Note the table has **no `c` row**.

Accompanying regression (Eq. 3, p. 118), Kinetex C18 vs Kinetex XB-C18 across all system constants
and all compositions:
`SC(Kinetex C18) = 1.008(±0.005)·SC(Kinetex XB-C18) − 0.017(±0.017)`,
r = 1.000, r_a² = 0.999, SE = 0.043, F = 36,994, n = 35. Slope CI 0.997–1.019, intercept CI
−0.002 to −0.031. The paper's conclusion: these two columns are (virtually) interchangeable for
small neutral compounds.

A second, contrasting pair (Eq. 4, p. 119), Kinetex C18 vs Kinetex Biphenyl, `v` constant only,
20–70% v/v methanol: `v(Kinetex C18) = 0.884(±0.040)·v(Kinetex Biphenyl) + 0.363(±0.088)`,
r = 0.996, r_a² = 0.990, F = 484, SE = 0.059, n = 6. Slope and intercept are significantly
different from 1 and 0 — i.e. **not** selectivity-equivalent, and the difference is
composition-dependent.

### `wsu2019-table3-group-mean-centred-system-constants.csv`
Table 3, journal pp. 120–121 (continues across the page break). Group assignment by principal
component factor analysis. For each of methanol–water and acetonitrile–water, two column groups,
each with mean-centred system constants at **30, 40, 50, 60, 70% (v/v)** with standard deviation
and range (max − min). 100 rows, complete — nothing illegible.

Group membership (column numbers refer to Table 1 / Fig. 5 legend, which is on the already-read
pages):
- MeOH Group I: 1–6, 9–13 (octadecylsiloxane-bonded, non-polar endcapping, 11 columns)
- MeOH Group II: 15, 17, 18, 21, 22
- ACN Group I: 1–6, 9, 10, 12, 13, 16
- ACN Group II: 18–21

**Important limitation:** this is *not* the per-column c,e,s,a,b,v database. It is a group mean
with a dispersion measure, and it **omits `c` entirely** — only e, s, a, b, v appear. It also
starts at 30%, not 10%.

### `wsu2019-intext-system-constant-differences.csv`
Not a numbered table — quantities stated only in running prose on pp. 119 and 121 (average
Group I − Group II differences per modifier, the Kinetex C18 − Kinetex Biphenyl differences, and
the isoeluotropic methanol-vs-acetonitrile differences behind Fig. 10A). Captured because they
are real numbers with SDs that would otherwise be lost. Source location noted per row.

## What could NOT be extracted

- **Figures 7–12 are not digitised.** They are vector scatter/line plots without printed data
  labels. Fig. 8 and Fig. 9 (system-constant correlation plots), Fig. 10 (isoeluotropic system
  constant comparison), Fig. 11 (system constant differences vs % methanol), Fig. 12 (log k vs
  % methanol for the three steric-resistance categories) all encode data only as plotted points.
  Digitising them would produce invented precision; not attempted. The only figure-derived numbers
  captured are the ones the text states explicitly (in the in-text CSV above, and the Fig. 12
  model-predicted values below).
- **No supplementary tables.** Tables S-1…S-5 are referenced repeatedly and are where the actual
  per-column constants live. Appendix A gives only the generic
  `https://doi.org/10.1016/j.chroma.2019.04.027` supplementary-data link. Not in hand.

## Answer to the question keeping #22 open

**No. The full `c, e, s, a, b, v` × φ-series per-column databases do NOT appear in the main text.**
Confirmed by reading the paper to its end. The main text's only constant tables are Table 2
(deviations for one column pair, no `c`) and Table 3 (group means, no `c`). The per-column
databases are in supplementary Tables S-3 (methanol, 25 columns), S-4 (acetonitrile, 24 columns)
and S-5 (THF, 2 columns), with S-1 (columns excluded for cation exchange / steric resistance) and
S-2 (test compound set, with weak bases in the bottom portion) also supplementary.

**Issue #22 stays open.** A human must obtain the Supporting Information ZIP from the ScienceDirect
article page for DOI 10.1016/j.chroma.2019.04.027 (institutional access or interlibrary loan;
requesting it from the author directly is the fallback).
The main article PDF itself is now fully exploited — nothing further is to be gained from it.

## Bearing on project conventions

**log scale — `log₁₀`, not `ln`.** Settled unambiguously. The paper writes the solvation parameter
model as `log k = c + eE + sS + aA + bB + vV` (Eq. 7, p. 124) and the retention model as
`log k = a₀ + a₁φ + a₂φ²` (Eq. 5, p. 122); Figs. 11 and 12 are axis-labelled "Log k". **The 2.303×
factor against a blueprint paper using `ln k` is live and must be applied when mixing sources.**

**pH scale — the paper does not use one for its own measurements.** WSU-2019's own mobile phases
are unbuffered organic–water. The only pH mentioned is for the *hydrophobic-subtraction model*
being compared against: 50% (v/v) acetonitrile–aqueous phosphate buffer at **pH = 7** (p. 125),
scale unqualified (i.e. no `s_s`pH / `w_w`pH distinction drawn). Nothing here to reconcile with
the project's `s_s`pH standard for WSU data itself — but note that when project code mixes WSU
system constants with HSM column parameters, the two are measured at different pH regimes as well
as different temperatures.

**Temperature — 45 °C reconfirmed, with a qualification that matters.** Page 125 states plainly:
the hydrophobic-subtraction model uses **35 °C** with 50% (v/v) acetonitrile–aqueous phosphate
buffer pH 7, and **the solvation parameter model 45 °C** with 50% (v/v) acetonitrile–water. This
is a second, independent statement of 45 °C inside WSU-2019 itself (the first being the
Experimental section on the earlier pages), so #11's answer is solid for this paper.

*Qualification (p. 125, citing refs 12, 24, 30, 31, 51, 52, 64):* "it was shown that dilute buffer
solutions and a small change in temperature have only a small effect on the values of the system
constants for neutral compounds and are not expected to invalidate those conclusions where these
differences in mobile phase composition were not specifically identified." So the paper treats a
35 vs 45 °C mismatch as tolerable **for neutral compounds** when comparing models. It gives **no
quantitative temperature dependence** — no dc/dT, no van 't Hoff coefficients. Anyone wanting to
transfer these constants to another temperature has no basis in this paper for doing so; the claim
is only that the error is "small", unquantified. Do not encode a temperature correction from this.

**Validity floors — the dewetting floor is reconfirmed and given a second, column-specific form.**
- Fig. 10 caption (p. 121): the isoeluotropic comparison "covers the composition range for
  20–70% (v/v) acetonitrile and tetrahydrofuran **due to incomplete wetting of the stationary phase
  at 10% (v/v) organic solvent**."
- Fig. 9 caption (p. 119): Kinetex C18 vs Kinetex Biphenyl is restricted to 20–70% v/v methanol
  "due to incomplete wetting of the **Kinetex Biphenyl phase at 10% (v/v) methanol**."
- Table 3 itself starts at **30%**, not 10 or 20.

So the practical floor is **not a single global number**: 10% is unusable generally (incomplete
wetting), specific phases fail at 10%, and the paper's own PCA-based grouping only claims validity
from 30% upward. Relevant to #24 (φ interpolation): **do not extrapolate the φ series below 30%
v/v for group-level constants, and treat 10% as invalid for any column.**

**Steric resistance and cation exchange — §3.4 is entirely about this gap, and it is explicit.**
Page 122 states outright: the solvation parameter model is parameterised for **neutral molecules in
a homogeneous solvation medium**, and "additional contributions to retention from steric resistance
and electrostatic interactions are **not addressed**" by Eq. 1. Specifics worth carrying:

- *Steric resistance.* Detected by fitting `log k = a₀ + a₁φ + a₂φ²` (Eq. 5). Category I = good fit,
  no significant steric mechanism (most compounds). Category II = characteristic hook shape at low
  organic. Category III = sharp discontinuity or plateau. Compounds showing II/III on the water-rich
  side of the discontinuity **are excluded from the model**, and the excluded lists are in the
  sources cited in Table S-1. Affected compound classes: angular (benzophenones), rigid planar
  (PAHs), long alkyl chains (n-alkylphenones), bulky (dialkyl phthalates). It is **not** predictable
  from size alone; conformational rigidity and H-bonding functional groups have an undefined role.
  System properties that matter: mobile phase composition, bonding density, pore size, endcapping
  type, and possibly pore dewetting. Fig. 12 gives the only concrete numbers: on Ascentis C18,
  model-predicted log k at 10, 20, 30% v/v methanol are **2.276, 1.994, 1.781** for 1-nitrohexane
  and **2.556, 2.174, 1.848** for n-butyrophenone — both of which show *loss* of retention over
  10–30% methanol relative to prediction.
- *Cation exchange.* Silica-based phases carry a low concentration of accessible anionic silanol
  groups acting as cation-exchange sites for bases. A screening procedure identifies affected
  compounds (weak bases listed in the bottom portion of Table S-2); columns where cation exchange
  matters are summarised in **Table S-1**, and those compounds are excluded from the model. Key
  asymmetry: cation-exchange contributions are **typically associated with acetonitrile-containing
  mobile phases and are suppressed or minimal for methanol–water** on the same compounds and
  columns (attributed speculatively to proton "water wires" that methanol disrupts). The relative
  contribution is **smallest at low φ and increases at high φ**, because reversed-phase retention
  generally declines.
- *Ionised compounds.* Explicitly outside the model: retention factors for partially ionised
  compounds are composite weighted averages of neutral and ionic forms, and descriptors for
  partially ionisable compounds depend on effective mobile phase pH and pKa — i.e. they become
  system-dependent, not solute properties. Two ion descriptors exist (ref 50) but "has not been
  applied to reversed-phase liquid chromatography".

**Consequence for the project's open question about whether that gap matters:** it matters
selectively and in a knowable way. The database's constants are valid for neutral, non-sterically
excluded compounds; the excluded-compound and excluded-column lists that would let us *detect*
when a query solute falls outside that envelope are in **Tables S-1 and S-2 — supplementary, not in
hand.** That is a second, independent reason to obtain the SI beyond the constants themselves.

## Other findings from these pages

- **Hydrophobic-subtraction model comparison (§3.5, pp. 123–125).** The two models are stated to be
  not term-by-term comparable. Quantitatively, for the columns in Table 1 the coefficient of
  determination `r_a²` for a plot of HSM **H** against solvation-parameter **v** is **0.723**; for
  HSM **A** against `a` it is **0.257**; for HSM **B** against `b` it is **0.051**. The conclusion
  (also in the paper's Conclusions) is that the two parameterisations show "little overlap" and a
  quantitative comparison is impractical. **This is a direct argument against any naive
  HSM↔LSER bridge on hydrogen-bonding terms** — relevant to any ticket contemplating one.
- **Conclusions (p. 125)** restate scope as **10–70% (v/v) acetonitrile or methanol**, 25 type-B
  silica-based columns, descriptors from an updated WSU database — consistent with pages 1–8.
- **Method-development framing:** the paper's own finding is that varying mobile phase composition
  is a more powerful selectivity lever than swapping columns within this database; column
  substitution is for fine tuning. Groups I and II are internally narrow in selectivity.
