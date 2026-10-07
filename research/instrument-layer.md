# Instrument layer — equations and citation package

Resolves #9. Part of the map, #1.

**What this is.** Per the map's cut rule — *"Fully spec what contains a decision. Where only a known
algorithm exists, write the interface and cite the literature"* — this layer is settled science. This
document is therefore **equations plus primary references**, not a design. Every equation is stated
explicitly, every symbol defined, every unit given, and every result traced to a DOI. Where I could
not reach a primary source I say so in [§7](#7-what-i-could-not-determine).

**Why the layer exists.** Retention fixes where peaks are; *peak width* fixes whether they are
resolved. Resolution is what method development optimises:

```
Rs = (1/4) · √N · ((α − 1)/α) · (k₂/(1 + k₂))
```

Everything below feeds `√N`, or corrupts it (extra-column dispersion), or improves it (gradient band
compression), or shifts the whole `k(φ)` trajectory in time (dwell volume).

---

## Contents

1. [Column efficiency — van Deemter, Giddings, Knox](#1-column-efficiency)
2. [Diffusion coefficients — `D_m`](#2-diffusion-coefficients-d_m)
3. [Extra-column dispersion](#3-extra-column-dispersion)
4. [Gradient band compression — Snyder's `G`](#4-gradient-band-compression--snyders-g-factor)
5. [Dwell volume — ranges and measurement](#5-dwell-volume)
6. [Uncertainty inventory (for #14)](#6-uncertainty-inventory-for-14)
7. [What I could not determine](#7-what-i-could-not-determine)
8. [Reference list](#8-reference-list)

---

## 1. Column efficiency

### 1.1 Plate count and plate height

```
N = L / H                                                             (1.1)

N = 16 (t_R / w_b)²  = 5.54 (t_R / w_½)²                              (1.2)

σ_t = t_R / √N        σ_V = V_R / √N                                  (1.3)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `N` | plate count (theoretical plates) | dimensionless |
| `L` | column length | cm (or m — be consistent with `H`) |
| `H` | plate height (HETP) | cm |
| `t_R` | retention time | s or min |
| `w_b`, `w_½` | peak width at base (4σ) and half height (2.355σ) | same as `t_R` |
| `σ_t`, `σ_V` | peak standard deviation in time / volume units | s, µL |
| `V_R` | retention volume = `F · t_R` | µL |

Equation (1.2) is the IUPAC/compendial definition; note it assumes a Gaussian peak.

### 1.2 van Deemter

The original equation, derived for gas chromatography of a linear isotherm system:

```
H = A + B/u + C·u                                                     (1.4)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `H` | plate height | cm |
| `u` | mobile-phase linear velocity (interstitial or superficial — state which) | cm/s |
| `A` | eddy diffusion / multipath term | cm |
| `B` | longitudinal (axial) molecular diffusion term, `B = 2 γ D_m` | cm²/s |
| `C` | resistance-to-mass-transfer term | s |
| `γ` | obstruction factor of the packed bed (≈ 0.6–0.7) | dimensionless |
| `D_m` | solute diffusion coefficient in the mobile phase | cm²/s |

**Source.** J.J. van Deemter, F.J. Zuiderweg, A. Klinkenberg, *Longitudinal diffusion and resistance
to mass transfer as causes of nonideality in chromatography*, Chem. Eng. Sci. **5** (1956) 271–289.
DOI [10.1016/0009-2509(56)80003-1](https://doi.org/10.1016/0009-2509%2856%2980003-1).

`u` in the linear form is proportional to volumetric flow rate:

```
u = F / (π r_c² ε_t)                                                  (1.5)
```

where `F` = volumetric flow rate (cm³/s), `r_c` = column internal radius (cm), `ε_t` = total column
porosity (dimensionless, typically ≈ 0.6–0.7 for fully porous silica; use `ε_e` ≈ 0.4 if you want
interstitial velocity). The choice of porosity is a real convention trap — see
[§6](#6-uncertainty-inventory-for-14).

**Caveat, load-bearing for an RPLC simulator.** Equation (1.4) is a GC-era result. Its assumptions
(negligible pressure drop, no viscous heating, velocity-independent `A`) are not strictly valid in
modern HPLC/UHPLC. This is documented in G. Guiochon, F. Gritti's group review:
*The van Deemter equation: assumptions, limits, and adjustment to modern high performance liquid
chromatography*, J. Chromatogr. A **1302** (2013) 1–13,
[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0021967313009278). Practically:
treat (1.4) as an empirical 3-parameter fit, not a mechanistic model.

### 1.3 Giddings coupling and the Knox equation (the reduced form — use this one)

Giddings showed that the `A` term is not velocity-independent; eddy diffusion couples to mobile-phase
mass transfer, producing a `u^(1/3)` dependence.

**Source.** J.C. Giddings, *Reduced plate height equation: a common link between chromatographic
methods*, J. Chromatogr. A **13** (1964) 301–304. DOI
[10.1016/S0021-9673(01)95123-4](https://doi.org/10.1016/S0021-9673%2801%2995123-4). See also
*'Eddy' diffusion in chromatography*, Nature **184** (1959) 357–358, DOI
[10.1038/184357a0](https://doi.org/10.1038/184357a0).

**Reduced parameters** (the dimensionless form — this is what the simulator should carry internally,
because it makes column comparison particle-size-independent):

```
h = H / d_p                    reduced plate height        (dimensionless)   (1.6)

ν = u · d_p / D_m              reduced velocity            (dimensionless)   (1.7)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `d_p` | particle diameter | cm (convert from µm) |
| `D_m` | solute diffusion coefficient in the mobile phase | cm²/s |
| `u` | linear velocity, same convention as (1.5) | cm/s |

**Knox equation:**

```
h = A ν^(1/3) + B/ν + C ν                                             (1.8)
```

with, for a well-packed conventional column, the canonical values `A ≈ 1`, `B ≈ 2`, `C ≈ 0.05`,
giving `h_min ≈ 2` at `ν ≈ 3–5`. `h_min ≈ 2` is the classic "well-packed column" acceptance
criterion.

**Source.** P.A. Bristow, J.H. Knox, *Standardization of test conditions for high performance liquid
chromatography columns*, Chromatographia **10** (1977) 279–289. DOI
[10.1007/BF02263001](https://doi.org/10.1007/BF02263001).

Kinetic performance limits (the `N`-vs-time trade-off, the basis of kinetic plots) come from:
J.H. Knox, M. Saleem, *Kinetic conditions for optimum speed and resolution in column
chromatography*, J. Chromatogr. Sci. **7** (1969) 614–622. DOI
[10.1093/chromsci/7.10.614](https://doi.org/10.1093/chromsci/7.10.614). A modern, **open-access**
restatement with worked derivations: G. Desmet et al., *Knox–Saleem kinetic performance limits in
liquid chromatography — a contemporary tutorial*, J. Chromatogr. Open **6** (2024) 100184, DOI
[10.1016/j.jcoa.2024.100184](https://doi.org/10.1016/j.jcoa.2024.100184). Use this as the
implementer's entry point; it is free and current.

### 1.4 How `N` actually depends on the knobs

Substituting (1.6)–(1.8) into (1.1):

```
N = L / (d_p · h(ν))                                                  (1.9)
```

Read off the dependencies:

| Knob | Effect on `N` | Mechanism |
|---|---|---|
| `L` ↑ | `N` ∝ `L` (linear) at fixed `ν` | more plates per unit `h` |
| `d_p` ↓ | `N` ∝ `1/d_p` at fixed `h`, **and** `h_min` reached at higher `F` | shorter diffusion path |
| `F` ↑ | `N` falls once past `ν_opt`; the `C ν` term dominates | mass-transfer limited |
| `F` ↓↓ | `N` also falls, via `B/ν` | longitudinal diffusion |
| `D_m` ↑ (higher T, less viscous φ) | shifts `ν_opt` to higher `F`; flattens the C-branch | `ν` ∝ `1/D_m` |
| `d_c` (column i.d.) | **no direct effect on `N`** | but sets `V_0` and hence the extra-column penalty — see §3 |

The `d_c` row is the one implementers get wrong. Column diameter does not change plate count; it
changes how badly the *instrument* degrades the observed plate count.

Pressure is the constraint, not a term: `ΔP = φ_f η L u / d_p²` (Kozeny–Carman form), `φ_f` ≈ 500–1000
flow resistance, `η` = mobile-phase viscosity (Pa·s). Halving `d_p` quadruples `ΔP` at fixed `u` and
`L` — this is the whole reason UHPLC exists, and is why §3 matters so much.

### 1.5 Interface sketch

```python
def plate_height_reduced(nu: float, A: float, B: float, C: float) -> float:
    """Knox eq. (1.8). Dimensionless in, dimensionless out."""

def plate_count(L_cm, d_p_cm, u_cm_s, D_m_cm2_s, A, B, C) -> float:
    """Eq. (1.9) via (1.7) and (1.8)."""
```

`A`, `B`, `C` are **column properties**, not universal constants. They must be either (a) fitted per
column from a van Deemter/Knox curve, or (b) defaulted to (1, 2, 0.05) with an explicit uncertainty.
This is a data-provisioning question for the column-set ticket, not a modelling decision.

---

## 2. Diffusion coefficients `D_m`

`D_m` appears in `B` and in `ν`. It is **the weakest quantitative link in this whole layer**, and the
ticket asked me to say how bad it is. Answer: 10–20% error is the realistic best case, and it gets
worse in exactly the regimes RPLC method development uses.

### 2.1 Wilke–Chang

```
D_m = 7.4 × 10⁻⁸ · (ψ M_2)^(1/2) · T / (η V_1^0.6)                    (2.1)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `D_m` | mutual diffusion coefficient at infinite dilution | cm²/s |
| `ψ` | association factor of the **solvent** | dimensionless |
| `M_2` | molar mass of the solvent | g/mol |
| `T` | absolute temperature | K |
| `η` | solvent (mobile-phase) viscosity | cP (= mPa·s) |
| `V_1` | molar volume of the **solute** at its normal boiling point (Le Bas additive volumes) | cm³/mol |

Association factors given in the original paper: `ψ` = 2.6 (water), 1.9 (methanol), 1.5 (ethanol),
1.0 (unassociated solvents, e.g. benzene, heptane).

**Source.** C.R. Wilke, P. Chang, *Correlation of diffusion coefficients in dilute solutions*,
AIChE J. **1** (1955) 264–270. DOI
[10.1002/aic.690010222](https://doi.org/10.1002/aic.690010222).

### 2.2 The mixed-solvent problem (this is the weak point)

Wilke–Chang is a **binary, single-solvent** correlation. An RPLC mobile phase is a
water/acetonitrile or water/methanol mixture, and the standard practice — take a mole-fraction-
weighted `ψM` and the mixture viscosity — has no theoretical backing. Two concrete failures:

1. **No published `ψ` for acetonitrile.** Wilke & Chang never gave one. Miyabe back-fitted
   `ψ(ACN) = 1.37` from peak-parking measurements of benzene in ACN at 303 K, and reports that with
   this value Wilke–Chang reproduces measured `D_m` in ACN/water with a mean square deviation of
   **8.8% and 14%** (two comparison datasets). K. Miyabe, *Estimation of molecular diffusivity in
   aqueous solution of acetonitrile by the Wilke–Chang equation*, J. Sep. Sci. **34** (2011)
   2119–2126. DOI [10.1002/jssc.201100385](https://doi.org/10.1002/jssc.201100385).
   → **Any implementation that uses Wilke–Chang with acetonitrile is using a fitted, non-original
   parameter.** Record that in the provenance.

2. **Viscosity maxima.** Water/organic mixtures have viscosity maxima (~20–30 %v/v ACN, ~40–50 %v/v
   MeOH) that are strongly non-ideal. Since `D_m ∝ 1/η`, an error in the mixture viscosity model
   propagates one-for-one into `D_m`. Li & Carr provide **fitting equations for the viscosity of
   MeOH/water and ACN/water** explicitly so that these correlations can be applied — use those, not
   a linear mixing rule.

### 2.3 Measured accuracy, primary source

J. Li, P.W. Carr, *Accuracy of empirical correlations for estimating diffusion coefficients in
aqueous organic mixtures*, Anal. Chem. **69** (1997) 2530–2536. DOI
[10.1021/ac961005a](https://doi.org/10.1021/ac961005a).

Alkylbenzenes and alkylphenones, Aris–Taylor open-tube method, 30–60 °C, 10–100 %v/v organic. Verbatim
findings:

- "The errors for methanolic mixtures by the Wilke–Chang, Scheibel, and Lusis–Ratcliff correlations
  are usually **less than 20%**."
- "The Scheibel, Wilke–Chang, and Hayduk–Laudie correlations work better than others for
  acetonitrile/water mixtures."
- "Overall, the **Scheibel correlation shows the smallest errors, and we recommend its use to that of
  the more widely used Wilke–Chang method** for the systems studied here."

Companion paper: J. Li, P.W. Carr, *Estimating diffusion coefficients for alkylbenzenes and
alkylphenones in aqueous mixtures with acetonitrile and methanol*, Anal. Chem. **69** (1997)
2550–2553. DOI [10.1021/ac961170q](https://doi.org/10.1021/ac961170q). They fit a **modified
Wilke–Chang correlation** specific to these solute classes and report "percent errors ... **no greater
than 10%** for both ACN/water and MeOH/water systems ... 2–3-fold better than" Wilke–Chang or
Scheibel. Two cautions: (a) it is a correlation fitted on alkylbenzenes/alkylphenones, so
extrapolation to drug-like, polar or ionisable solutes is unvalidated; (b) the paper explicitly
recommends it "for the evaluation of column performance", i.e. for test probes.

### 2.4 Successor correlations (candidates)

| Correlation | Primary source | DOI |
|---|---|---|
| Scheibel (recommended by Li & Carr) | E.G. Scheibel, *Correspondence: Liquid diffusivities*, Ind. Eng. Chem. **46** (1954) 2007–2008 | [10.1021/ie50537a062](https://doi.org/10.1021/ie50537a062) |
| Hayduk–Laudie (aqueous only) | W. Hayduk, H. Laudie, AIChE J. **20** (1974) 611–615 | [10.1002/aic.690200329](https://doi.org/10.1002/aic.690200329) |
| Tyn–Calus (parachor-based) | M.T. Tyn, W.F. Calus, J. Chem. Eng. Data **20** (1975) 106–109 | [10.1021/je60064a006](https://doi.org/10.1021/je60064a006) |

### 2.5 The escape hatch: measured databases

If you would rather look `D_m` up than estimate it, the largest chromatography-relevant compilation:

H. Song, Y. Vanderheyden, E. Adams, G. Desmet, D. Cabooter, *Extensive database of liquid phase
diffusion coefficients of some frequently used test molecules in reversed-phase liquid chromatography
and hydrophilic interaction liquid chromatography*, J. Chromatogr. A **1455** (2016) 102–112. DOI
[10.1016/j.chroma.2016.05.054](https://doi.org/10.1016/j.chroma.2016.05.054).

45 compounds, Taylor–Aris measurements. Verbatim from the abstract: correlations such as Wilke–Chang
"can provide good approximations of molecular diffusion under reversed-phase conditions. However,
these correlations have been demonstrated to be **less accurate for mobile phases containing a large
percentage of acetonitrile**." They also show buffer concentration (5–10 mM ammonium acetate)
"mainly influence[s] the molecular diffusion of charged molecules" — relevant, since ionisable
compounds are in scope for this effort.

### 2.6 Sensitivity — how much does a bad `D_m` actually hurt?

`D_m` enters `N` only through `ν` (eq. 1.7) and `B` (eq. 1.4). Near `ν_opt` the Knox curve is flat, so
a 20% error in `D_m` moves `h` by only a few percent — for a well-designed method. On the
`B`-dominated (low flow) branch, or when using `D_m` for a kinetic-plot optimisation far from `ν_opt`,
the error propagates ~1:1. **Verdict: `D_m` is a 10–20% quantity, and that is tolerable for a
resolution map but not for absolute peak-width prediction.** Carry it as a distribution, not a point.

### 2.7 Interface sketch

```python
def diffusion_coefficient(
    solute_molar_volume_cm3_mol: float,   # Le Bas, or from V (McGowan) × 0.6 — see caveat below
    solvent: MobilePhase,                 # composition φ, modifier identity, T
    correlation: Literal["wilke_chang", "scheibel", "li_carr_modified"] = "scheibel",
) -> Distribution[float]:                 # cm²/s, with an uncertainty band
```

**Note a convenient coupling with the chemistry layer:** the Abraham `V` descriptor (McGowan
characteristic volume, cm³/mol/100) is already being predicted by the QSPR layer. Le Bas `V_1` and
McGowan `V` are different quantities but strongly correlated; a conversion would avoid a second
molecular-property pipeline. **I have not verified a published conversion** — flagged as an open
question below.

---

## 3. Extra-column dispersion

### 3.1 Additivity of variances

Independent broadening contributions add as **variances** (never as widths):

```
σ²_obs = σ²_col + σ²_ec                                               (3.1)

σ²_ec  = σ²_inj + σ²_tubing + σ²_fittings + σ²_cell + σ²_electronic    (3.2)
```

All terms in µL² (volume units — this is the right basis, because the instrument contributions are
volumetric and flow-rate-normalised). Convert to time via `σ_t = σ_V / F`.

The apparent (observed) plate count is therefore always **less than** the true column plate count:

```
N_obs = V_R² / σ²_obs = V_R² / (σ²_col + σ²_ec)                        (3.3)

N_obs / N_col = 1 / (1 + σ²_ec / σ²_col)                               (3.4)
```

**Source for the framework.** L.R. Snyder, J.J. Kirkland's texts, and originally
S.J. Hawkes / J.C. Sternberg. The canonical modern treatment used throughout this section:
K. Broeckhoven, G. Desmet, *Extra-column band broadening effects in contemporary liquid
chromatography: causes and solutions*, TrAC Trends Anal. Chem. **119** (2019) 115619. DOI
[10.1016/j.trac.2019.115619](https://doi.org/10.1016/j.trac.2019.115619).

### 3.2 Tubing — Taylor–Aris

For a straight open cylindrical tube in laminar flow, at long residence times:

```
σ²_V,tube = (π · d_t⁴ · L_t · F) / (384 · D_m)                        (3.5)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `σ²_V,tube` | volumetric variance contributed by the tube | cm⁶ (convert to µL²: 1 cm⁶ = 10⁶ µL²) |
| `d_t` | tube internal diameter | cm |
| `L_t` | tube length | cm |
| `F` | volumetric flow rate | cm³/s |
| `D_m` | solute diffusion coefficient | cm²/s |

**The `d_t⁴` dependence is the single most actionable fact in this section.** Going from 0.010″
(254 µm) to 0.005″ (127 µm) i.d. tubing cuts tubing variance **16-fold**.

**Sources.** G.I. Taylor, *Dispersion of soluble matter in solvent flowing slowly through a tube*,
Proc. R. Soc. Lond. A **219** (1953) 186–203, DOI
[10.1098/rspa.1953.0139](https://doi.org/10.1098/rspa.1953.0139); R. Aris, *On the dispersion of a
solute in a fluid flowing through a tube*, Proc. R. Soc. Lond. A **235** (1956) 67–77, DOI
[10.1098/rspa.1956.0065](https://doi.org/10.1098/rspa.1956.0065).

**Validity caveat.** (3.5) requires the diffusive equilibration time across the tube radius to be
short relative to residence time, i.e. `L_t · D_m / (u_t · d_t²) ≳ 1`. In short UHPLC connection
capillaries at high flow this is violated, and the true dispersion is *lower* than Taylor–Aris
predicts. The transient regime is treated in M.J.E. Golay, J.G. Atwood, *Early phases of the
dispersion of a sample injected in Poiseuille flow*, J. Chromatogr. A **186** (1979) 353–370, DOI
[10.1016/S0021-9673(00)95261-0](https://doi.org/10.1016/S0021-9673%2800%2995261-0).

### 3.3 Injector, detector cell, fittings

```
σ²_inj  = V_inj² / K            K ≈ 12 (plug injection) … 4 (fully mixed)   (3.6)

σ²_cell = V_cell² / 12          (well-swept cell, plug-flow idealisation)   (3.7)

σ²_elec = (F · τ)²              τ = detector time constant (s)              (3.8)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `V_inj` | injected sample volume | µL |
| `V_cell` | detector flow-cell volume | µL |
| `τ` | detector response time constant | s |
| `K` | injection-profile factor; `12` for an ideal rectangular plug | dimensionless |

`K` is a modelling choice with a range, not a constant — real injectors are between plug and mixed.

**Fittings are the term that does not obey a clean equation.** Poorly swept dead volumes at
connections produce **tailing (τ-type / exponentially-modified-Gaussian) broadening, not Gaussian
σ-type**, so they do not add cleanly in variance and they distort peak shape. Gritti & Guiochon
isolated σ-type from τ-type contributions by Foley–Dorsey EMG fitting and found "combined
contributions from the injector and connecting tubing ... larger than expected by Taylor–Aris
theory", with **τ-type contributions the main cause**, from poorly-swept volumes at the
injector–tubing connection. Practical consequence for the simulator: (3.2) will *underestimate* real
extra-column dispersion on a badly plumbed system, and no equation will save you — the total must be
measured.

### 3.4 When extra-column dispersion dominates — the modern UHPLC case

**Measured, primary.** F. Gritti, G. Guiochon, *On the extra-column band-broadening contributions of
modern, very high pressure liquid chromatographs using 2.1 mm I.D. columns packed with sub-2 µm
particles*, J. Chromatogr. A **1217** (2010) 7677–7689. DOI
[10.1016/j.chroma.2010.10.016](https://doi.org/10.1016/j.chroma.2010.10.016). Verbatim:

> "When the 1290 Infinity HPLC System is equipped with a needle seat, an inlet and an outlet
> connecting capillary tube with inner diameters around 115 µm, its **extra-column variance for a
> 0.1 µL injection volume is 9.2 µL²** while that of the **Acquity instrument is 6.9 µL²**. Minor
> modifications suggested by their respective manufacturers allowed significant reductions of these
> variances, to **6.2 and 3.9 µL²**, respectively. Yet, in their optimized configurations and for
> **weakly retained compounds (k ≈ 1), these modern, sophisticated instruments cannot provide more
> than 75% (1290 Infinity) and 85% (Acquity) of the maximum efficiency of a 2.1 mm × 50 mm BEH
> column**. For more strongly retained compounds (k > 4), in contrast, they are both able to provide
> more than 95% of the maximum expected efficiency."

**That quotation is the answer to "when does it dominate", and it should be quoted in the spec.**

Worked arithmetic (mine, from eqs. 1.1/1.3, assuming `ε_t` = 0.65, `h` = 2):

| Column | `V_0` (µL) | `N` | `σ²_col` at `k` = 0 (µL²) | at `k` = 1 | at `k` = 5 |
|---|---|---|---|---|---|
| 2.1 × 50 mm, 1.7 µm | 113 | ~14,700 | **0.9** | **3.4** | 31 |
| 4.6 × 150 mm, 5 µm | 1620 | ~15,000 | 175 | **700** | 6300 |

Against a measured `σ²_ec` of 3.9–9.2 µL²:

- **UHPLC, narrow-bore, early peak (`k` ≲ 1): the instrument contributes as much or more variance
  than the column.** At `k` = 0 it is a 4–10× *excess*. The observed first peaks in a UHPLC
  chromatogram are an instrument artefact more than a column measurement.
- **UHPLC, `k` > 4: negligible** (< 5% efficiency loss), consistent with Gritti & Guiochon.
- **Conventional 4.6 mm HPLC: negligible everywhere** — `σ²_ec` is ~1% of `σ²_col` even at `k` = 1.
  A 30 µL² instrument costs < 5% efficiency for any `k` > 1.

Corroborating study across many commercial systems: S. Fekete, J. Fekete, *The impact of extra-column
band broadening on the chromatographic efficiency of 5 cm long narrow-bore very efficient columns*,
J. Chromatogr. A **1218** (2011) 5286–5291. DOI
[10.1016/j.chroma.2011.06.045](https://doi.org/10.1016/j.chroma.2011.06.045).

### 3.5 How `σ²_ec` is measured

The **zero-length / union replacement** method: replace the column with a zero-dead-volume union,
inject, and measure the peak variance directly. Systematically biased low, because the union's
flow profile is not the column's.

The **linear extrapolation (non-invasive) method**: measure `σ²_obs` for a homologous series with the
column *in place*, plot `σ²_obs` against `V_R²`, and take the intercept at `V_R` = 0 as `σ²_ec`.
Originally: H.H. Lauer, G.P. Rozing, *The selection of optimum conditions in HPLC I. The
determination of external band spreading in LC instruments*, Chromatographia **14** (1981) 641–647,
DOI [10.1007/BF02291104](https://doi.org/10.1007/BF02291104).

**Interface implication.** `σ²_ec` is a **single scalar per instrument configuration** (per flow
rate, strictly). The simulator should accept it as one number, defaulted by instrument class and
updatable from data — exactly the same prior/update pattern the map assigns to dwell volume. Do not
build a component-by-component fitting model out of (3.2); the fittings term makes it unidentifiable.

---

## 4. Gradient band compression — Snyder's `G` factor

### 4.1 Why it exists

In a gradient, the tail of a band sits in a stronger eluent than its front. The tail moves faster.
The band narrows. Ignoring `G` therefore makes **every gradient resolution prediction pessimistic** —
you predict wider peaks than you get, and you under-report `Rs`.

### 4.2 Gradient steepness parameter `b`

```
b = (V_M · Δφ · S) / (t_G · F)       (equivalently  b = t_0 · Δφ · S / t_G )   (4.1)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `b` | gradient steepness parameter | dimensionless |
| `V_M` | column hold-up (dead) volume | mL |
| `t_0` | column dead time = `V_M / F` | min |
| `Δφ` | `φ_final − φ_initial`, change in volume fraction organic | dimensionless (0–1) |
| `S` | LSS slope, from `ln k = ln k_w − S φ` | dimensionless |
| `t_G` | gradient duration | min |
| `F` | flow rate | mL/min |

Note `S` here is the slope of **natural log** `k` vs `φ`. If your retention layer stores `log₁₀ k`,
the factor of 2.303 changes hands — this is the classic sign/base bug in gradient code. Snyder's
original literature uses `b = t_0 Δφ S / t_G` with `S` on the natural-log convention, and the `2.3`
that appears inside `G` below is `ln 10` from the base conversion in his derivation. **Fix the
convention once, in a single place, and assert it in tests.**

Definition confirmed against X. Wang, D.R. Stoll, A.P. Schellinger, P.W. Carr, *Peak capacity
optimization of peptide separations in reversed-phase gradient elution chromatography: fixed column
format*, Anal. Chem. **78** (2006) 3406–3416, DOI
[10.1021/ac0600149](https://doi.org/10.1021/ac0600149), which gives
`b = S(φ_final − φ_initial) V_m / (F · t_G)`
([open access via PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2638764/)).

### 4.3 Retention factor at elution and gradient retention time

```
k_e = 1 / (2.3 b)                                                     (4.2)

t_R = (t_0 / b) · log₁₀(2.3 b k_0 + 1) + t_0 + t_D                    (4.3)
```

| Symbol | Meaning | Unit |
|---|---|---|
| `k_e` | instantaneous retention factor at the moment the band exits | dimensionless |
| `k_0` | retention factor at the *initial* mobile phase composition `φ_initial` | dimensionless |
| `t_D` | dwell time = `V_D / F` — see §5 | min |

Equation (4.3) is the LSS gradient retention equation. **Sources.** L.R. Snyder, J.W. Dolan,
J.R. Gant, *Gradient elution in high-performance liquid chromatography. I. Theoretical basis for
reversed-phase systems*, J. Chromatogr. A **165** (1979) 3–30, DOI
[10.1016/S0021-9673(00)85726-X](https://doi.org/10.1016/S0021-9673%2800%2985726-X); Part II,
J. Chromatogr. A **165** (1979) 31–58, DOI
[10.1016/S0021-9673(00)85727-1](https://doi.org/10.1016/S0021-9673%2800%2985727-1).

### 4.4 The `G` factor

```
                 ┌───────────────────────┐
                 │ 1 + 2.3b + (2.3b)²/3  │
G  =   √         └───────────────────────┘   ⁄  (1 + 2.3b)            (4.4)
```

i.e. `G = sqrt(1 + 2.3·b + (2.3·b)²/3) / (1 + 2.3·b)`, dimensionless, `b` from (4.1).

Peak width in gradient elution:

```
w_½ = 2.35 · G · t_0 · (1 + k_e) / √N                                 (4.5)

σ_t = G · t_0 · (1 + k_e) / √N                                        (4.6)
```

`N` is the **isocratic** plate number of the column at the elution composition. Form (4.5) is the one
given by Wang/Stoll/Schellinger/Carr (above), verbatim:
"W½ = 2.35 G t₀(1 + k′f)/√N where N is the isocratic plate number and G is the gradient band
compression factor".

**Primary sources for `G`.** L.R. Snyder, D.L. Saunders, *Optimized solvent programming for
separations of complex samples*, J. Chromatogr. Sci. **7** (1969) 195–208, DOI
[10.1093/chromsci/7.4.195](https://doi.org/10.1093/chromsci/7.4.195) (first postulation of peak
compression); the modern canonical statement is L.R. Snyder, J.W. Dolan, *High-Performance Gradient
Elution: The Practical Application of the Linear-Solvent-Strength Model*, Wiley, 2007, DOI
[10.1002/0470055529](https://doi.org/10.1002/0470055529) — **cite the book, chapter 2 ("Gradient
Elution Fundamentals", DOI [10.1002/9780470055526.ch2](https://doi.org/10.1002/9780470055526.ch2)),
as the implementation reference.** The independent theoretical treatment is H. Poppe, J. Paanakker,
M. Bronckhorst, *Peak width in solvent-programmed chromatography. I. General description of peak
broadening in solvent-programmed elution*, J. Chromatogr. A **204** (1981) 77–84, DOI
[10.1016/S0021-9673(00)81641-6](https://doi.org/10.1016/S0021-9673%2800%2981641-6).

### 4.5 Realistic magnitude — what "typical" means

The ticket asked for this specifically. Evaluating (4.4):

| `b` | `G` | Peak narrowing vs no compression | Regime |
|---|---|---|---|
| 0.05 | 0.949 | 5% | very shallow gradient |
| 0.10 | 0.908 | 9% | shallow |
| **0.20** | **0.847** | **15%** | **Snyder's recommended optimum for method development** |
| 0.30 | 0.805 | 20% | typical general-purpose scouting gradient |
| **0.40** | **0.773** | **23%** | **typical fast/screening gradient** |
| 0.50 | 0.749 | 25% | steep |
| 1.0 | 0.682 | 32% | very steep / ballistic |
| 2.0 | 0.635 | 37% | — |
| ∞ | 0.577 (= 1/√3) | 42% | theoretical limit of eq. (4.4) |

*(Values computed from eq. 4.4; the `b` → ∞ limit `1/√3` = 0.577 is an analytic property of the
equation and serves as a self-check on any implementation.)*

**Headline number for the spec: in the `b` = 0.2–0.4 window that covers most real analytical
gradients, `G` ≈ 0.77–0.85, i.e. peaks are 15–25% narrower than the isocratic-width calculation
predicts.** Since `Rs` ∝ 1/w, **omitting `G` under-predicts gradient resolution by 15–25%** — enough
to declare a separation failed when it is fine, and therefore enough to send the optimiser to the
wrong place. `G` is not an optional refinement.

**Bounds and caveats.**
- (4.4) assumes LSS (`ln k` linear in `φ`), a linear gradient, an unretained modifier, and constant
  `N`. Relaxing these matters: Gritti & Guiochon showed the *ultimate* compression factor for a
  retained modifier can exceed the LSS value — "the stronger the retention of the organic modifier,
  the more effective the band compression", and "the best time compression factor that could be
  expected is twice the one predicted with an ideal column". F. Gritti, G. Guiochon, *The ultimate
  band compression factor in gradient elution chromatography*, J. Chromatogr. A **1178** (2008)
  79–91, DOI [10.1016/j.chroma.2007.11.044](https://doi.org/10.1016/j.chroma.2007.11.044).
- Curvature in `ln k` vs `φ` (i.e. where the Neue–Kuss model beats LSS) changes `G`. See W. Hao et
  al., *Peak compression in linear gradient elution liquid chromatography*, J. Chromatogr. A
  **1619** (2020) 460908, DOI
  [10.1016/j.chroma.2020.460908](https://doi.org/10.1016/j.chroma.2020.460908), which recomputes `G`
  under a quadratic solvent-strength model and with `H` varying with `φ`.
- A modern, self-contained derivation of everything in §4 (including how extra-column dispersion
  interacts with gradient peak width): K. Broeckhoven, G. Desmet, *Theory of separation performance
  and peak width in gradient elution liquid chromatography: a tutorial*, Anal. Chim. Acta **1218**
  (2022) 339962, DOI [10.1016/j.aca.2022.339962](https://doi.org/10.1016/j.aca.2022.339962).
  **This is the single best implementation reference for the whole peak-width part of this layer.**

### 4.6 Interface sketch

```python
def band_compression_factor(b: float) -> float:
    """Eq. (4.4). Snyder & Dolan 2007 ch.2. Returns G ∈ (1/√3, 1]."""
    x = 2.3 * b
    return math.sqrt(1 + x + x*x/3) / (1 + x)

def gradient_peak_width_sigma(t_0, k_e, N, G) -> float:
    """Eq. (4.6), minutes."""
```

Note if the retention layer moves off LSS to Neue–Kuss (U.D. Neue, H.-J. Kuss, *Improved
reversed-phase gradient retention modeling*, J. Chromatogr. A **1217** (2010) 3794–3803, DOI
[10.1016/j.chroma.2010.04.023](https://doi.org/10.1016/j.chroma.2010.04.023)), eq. (4.4) is no longer
strictly consistent with it, and `G` must be recomputed numerically. **Flag: this is a coupling
between the retention-model ticket and this one that nobody has priced yet.**

---

## 5. Dwell volume

### 5.1 Definition

**Dwell volume `V_D`** (= gradient delay volume) is the volume of the fluidic path between the point
at which the solvents are mixed and the head of the column: pump mixer, proportioning valve,
connecting tubing, autosampler flow path. **It does not include the column.**

```
t_D = V_D / F                                                          (5.1)
```

Its effect: the gradient reaches the solutes `t_D` later than the pump programme says. In eq. (4.3)
`t_D` is a pure additive offset on `t_R` **only for solutes that are still fully retained during the
delay** (`k_0` large). For weakly retained solutes the delay is spent isocratically eluting at
`φ_initial`, which changes elution order, not just elution time. That is why dwell volume "breaks
method transfer" rather than merely shifting the chromatogram: **it is selectivity-affecting for
early peaks.**

### 5.2 Typical ranges by instrument class

Vendor technical note with a directly quoted instrument table — Avantor/ACE Knowledge Note #0001,
*How to Determine System Dwell Volume: Theory and Practice*,
<http://www.hplc.eu/Downloads/AKN0001_DwellVolume.pdf>:

| Instrument | Typical `V_D` (µL) |
|---|---|
| Agilent HP1100 **Binary** | 180–900 |
| Agilent HP1100 **Quaternary** | 800–1100 |
| Agilent 1200 RRLC | ~300 |
| Dionex P680A Quaternary | < 400 |
| Thermoquest P4000 Quaternary | < 600 |
| Waters Alliance 2695 Quaternary | 600 |
| Waters/Varian 9012 Ternary | 1000 |

Generalising, and consistent with the mechanism (low-pressure mixing must mix *before* the pump head,
so the pump head volume is inside `V_D`; high-pressure mixing mixes after the pumps):

| Class | Mixing | Typical `V_D` |
|---|---|---|
| Conventional HPLC, quaternary / low-pressure mixing | pre-pump | **600–1200 µL** |
| Conventional HPLC, binary / high-pressure mixing | post-pump | **200–900 µL** |
| UHPLC, quaternary | pre-pump | **< 500 µL** |
| UHPLC, binary / high-pressure mixing | post-pump | **~35–200 µL**, often 100–400 µL with a mixer fitted |

The map's stated range of "~30 µL to ~1000 µL" is corroborated. Also see the LCGC three-part series
by D.R. Stoll: *The Gradient Delay Volume, Part I: Theory*, LCGC Int. (2024) 6–10, DOI
[10.56530/lcgc.int.wj6080e9](https://doi.org/10.56530/lcgc.int.wj6080e9); *Part II: Practice —
Effects on Method Transfer*, LCGC Int. (2023) 6–10, DOI
[10.56530/lcgc.int.ia2675s7](https://doi.org/10.56530/lcgc.int.ia2675s7); *Part III: Practice —
Effects on Throughput*, LCGC Int. (2024) 6–11, DOI
[10.56530/lcgc.int.ha5479i5](https://doi.org/10.56530/lcgc.int.ha5479i5).

**Two facts that a naive "one number per instrument" model gets wrong:**

1. `V_D` is not perfectly constant. On low-pressure-mixing systems the effective delay depends on the
   proportioning duty cycle and hence on the target %B; and on any system it depends on whether a
   mixer is installed. Treat it as a distribution centred on the instrument-class prior.
2. Users can *change* it — swapping mixers, changing needle-seat capillaries, changing autosampler
   injection mode (bypass vs flow-through) all move `V_D`. So the "user-entered prior" must be
   attached to an **instrument configuration**, not an instrument.

### 5.3 THE DECISION SURFACE — candidate methods for measuring `V_D`

The map has already settled that dwell volume is *both* user-entered *and* measurable from scouting
data — a prior that data updates. **The open decision, made in a later ticket, is which measurement
route the scouting-run reconciliation uses.** Below are the candidates with accuracy and practical
requirements. **I am not choosing.**

---

#### Candidate A — Tracer linear ramp, column replaced by a union ("midpoint / 50% method")

**Procedure** (verbatim from Avantor KN #0001): replace the column with a zero-dead-volume connector;
solvent A = water, solvent B = water + 0.1 %v/v acetone; UV at 265 nm; run 0→100 %B over `t_G`
(10 min at 2 mL/min in their example). Find the absorbance midpoint `A_0.5` = (`A_max` − `A_min`)/2,
read off `t_0.5`, then:

```
t_D = t_0.5 − (t_G / 2)          V_D = t_D × F                        (5.2)
```

Their worked example: `A_max − A_min` = 123.138 mAU, `A_0.5` = 61.569 mAU, `t_0.5` = 5.566 min,
`t_D` = 0.566 min, **`V_D` = 1.13 mL**.

- **Accuracy:** high (few %), and the midpoint construction is robust to the mixer's dispersion
  because it uses a symmetric point on the ramp.
- **Requires:** a dedicated run, a column swap, a UV-absorbing tracer, and a detector. Roughly
  15–20 min of instrument time plus reconfiguration.
- **Failure modes:** acetone at 265 nm saturates the detector at high concentration; the tracer must
  be dilute. Refractive-index artefacts if the two solvents differ in RI. Column swap is a manual step
  a lab may resist.

#### Candidate B — Tracer linear ramp, **tangent-intercept** construction

Same run as A, but `t_D` is taken as the intersection of the extrapolated initial baseline with the
tangent to the rising portion.

- **Accuracy:** lower and **systematically different from A**. The tangent method measures the
  *first arrival* of gradient (≈ the true geometric volume plus zero mixing); the midpoint method
  measures the *median* arrival (geometric volume plus the mixer's mean residence). On a system with
  a large mixer the two disagree by roughly the mixer volume.
- **This disagreement is a real spec hazard.** Whichever definition is chosen must be the *same one*
  used inside eq. (4.3), or the model will be biased by up to a mixer volume (100–400 µL).
- **Requires:** same as A.

#### Candidate C — Tracer **step** change rather than a ramp

Step 0 %B → 100 %B and take the breakthrough time.

- Yields a value close to the **physical** (geometric) delay volume, because the step response is
  the system's residence-time distribution directly.
- **Accuracy:** good for the geometric volume; the RTD's spread (not just its mean) is visible, which
  is diagnostically useful. But it measures a *different* quantity from the LSS `t_D` in eq. (4.3).
- **Requires:** dedicated run and column swap, as A/B.

#### Candidate D — Tracer ramp **with the column in place**

Run A or B without removing the column, then subtract the column hold-up: `V_D = F·t_0.5 − t_G·F/2 − V_M`.

- **Accuracy:** degraded by however well `V_M` is known — and `V_M` itself is definition-dependent
  (±3–5%, see §6). On a UHPLC narrow-bore column `V_M` ≈ 113 µL is comparable to `V_D` itself, so
  the subtraction can be a large fraction of the answer.
- **Requires:** dedicated run, but **no column swap** — meaningfully cheaper operationally.
- Also: the tracer must be unretained on the column, which acetone is not reliably.

#### Candidate E — **`V_D` as a free parameter in the regression of scouting-gradient retention times** ← *the "no extra run" candidate*

Fit eq. (4.3) to the observed `t_R` of several solutes across ≥ 2 (preferably ≥ 3) gradient runs of
different `t_G` (and/or `F`), treating `V_D` as a shared instrument-level parameter alongside each
solute's (`k_w`, `S`).

- **Requires no dedicated experiment at all** — it rides on the scouting runs the system is already
  going to prescribe. This is the property the ticket flagged as valuable.
- **Established?** The two-or-more-gradient inversion is well established: M.A. Quarry, R.L. Grob,
  L.R. Snyder, *Prediction of precise isocratic retention data from two or more gradient elution
  runs. Analysis of some associated errors*, Anal. Chem. **58** (1986) 907–917, DOI
  [10.1021/ac00295a056](https://doi.org/10.1021/ac00295a056). That paper's explicit subject is the
  error analysis, including the error contributed by `V_D`. **What it does is show that `V_D` errors
  propagate into fitted (`k_w`, `S`) — which is exactly the identifiability statement needed here,
  read in reverse.**
- **Accuracy and identifiability — the crux.** `V_D` and the early part of the gradient are
  correlated parameters. Discriminating power comes from solutes that are *weakly* retained at
  `φ_initial` (small `k_0`), because those are the ones whose `t_R` is sensitive to `t_D` in a way
  that is not degenerate with `k_w`. If every scouting solute is strongly focused at the column head
  (`k_0` > 10), `t_D` enters (4.3) as a near-pure additive constant and is partly confounded with
  `t_0`. **Design implication for the scouting-design ticket: include at least one weakly-retained
  probe, or the dwell estimate will be poorly identified.**
- Varying `F` between scouting runs breaks the degeneracy cleanly, because `t_D = V_D/F` scales with
  flow while `V_D` does not. **This is the strongest lever and it costs nothing but a flow change.**
- **Reported to work in practice:** modelling-software practice reports that "an iterative dwell
  volume estimation was demonstrated to generate more accurate retention time predictions than when a
  practically determined dwell volume was used", and that inverse/fitting methods "do not require
  accurate knowledge of the column hold-up volume or system dwell volume in gradient elution" —
  M.R. Euerby / LCGC, *Anomalous retention prediction using modelling software in gradient
  reversed-phase liquid chromatography: why it can occur and how to prevent it*, LCGC,
  <https://www.chromatographyonline.com/view/anomalous-retention-prediction-using-modelling-software-in-gradient-reversed-phase-liquid-chromatography-why-it-can-occur-and-how-to-prevent-it>.
  **Provenance caveat: this is a trade-magazine article, not a peer-reviewed primary source. Treat
  the claim as a strong lead to verify, not as established.**
- **Risk to name explicitly:** the fitted `V_D` absorbs *all* unmodelled early-time error — column
  `V_M` error, `φ` proportioning error, sample-solvent effects, non-LSS curvature. It is an effective
  parameter, not necessarily the physical delay volume. If the same fitted `V_D` is then used to
  transfer a method to another instrument, that absorbed error transfers with it.

#### Candidate F — **Pressure-trace / viscosity-step method**

Run a step change between two solvents of different viscosity (e.g. water → acetonitrile) and detect
the arrival of the composition change as an inflection in the *system pressure* trace rather than the
detector trace.

- **No tracer, no detector, works with the column installed.**
- **Accuracy:** unquantified in what I found; the pressure signal is smeared by the column's own
  viscosity transit and is confounded by compressibility.
- **Provenance: I could not locate a primary peer-reviewed source for this method.** It circulates as
  practitioner know-how. **Do not spec it without validating it.**

#### Candidate G — Instrument-reported / vendor-specified value

The user-entered prior the map already assumes. Vendor documentation gives a nominal `V_D` per
configuration. **Accuracy: ±30–50% is realistic**, because it excludes user plumbing changes. Its role
is to be the prior, not the answer.

---

**Summary of the candidate set:**

| # | Method | Extra run? | Column swap? | Accuracy | Measures |
|---|---|---|---|---|---|
| A | Tracer ramp, midpoint | yes | yes | few % | delay incl. mixer mean residence |
| B | Tracer ramp, tangent | yes | yes | worse; biased vs A | geometric first-arrival |
| C | Tracer step | yes | yes | good | RTD of the fluidic path |
| D | Tracer ramp, column in place | yes | **no** | degraded by `V_M` uncertainty | delay + `V_M` error |
| **E** | **Fit `V_D` from scouting gradients** | **no** | **no** | identifiability-limited; needs varied `F` and/or a weakly-retained probe | *effective* `V_D` for the retention model |
| F | Pressure/viscosity step | yes | no | unknown | unclear |
| G | Vendor nominal | no | no | ±30–50% | nominal configuration |

**The observation to carry into the decision ticket:** A/B/C measure the *fluidic* delay volume;
E measures the *effective* delay parameter of the retention model. Those are not the same quantity,
and only E is guaranteed to be the one that minimises retention-prediction error. If the system's
purpose is prediction rather than instrument qualification, that asymmetry favours E — but E's
identifiability is conditional on scouting design, which makes this a **joint** decision with the
scouting-design ticket, not an independent one.

---

## 6. Uncertainty inventory (for #14)

Requested by the ticket for the uncertainty layer. Ordered from hardest to softest.

### Knowable precisely (treat as fixed, uncertainty negligible)

| Parameter | Typical uncertainty | Notes |
|---|---|---|
| Column length `L`, i.d. `d_c` | < 1% | manufacturer spec, hard geometry |
| Gradient time `t_G`, programmed `φ(t)` | < 0.5% | firmware timing |
| Flow rate `F` (setpoint reproducibility) | 0.1–0.5% RSD | pump spec; *accuracy* vs *precision* differ — accuracy can be 1–2% |
| Temperature setpoint | ±0.5 K | but the *column* temperature ≠ oven setpoint under viscous heating in UHPLC |
| Injection volume | 0.5–1% RSD | |

### Measured, with real uncertainty

| Parameter | Typical uncertainty | Source of uncertainty |
|---|---|---|
| Column hold-up volume `V_M` / `t_0` | **±3–5%** | *definition-dependent*: unretained-marker vs pycnometric vs minor-disturbance methods give systematically different answers. Pick one definition, record it. |
| Dwell volume `V_D` | **±5–15%** by method A/C; **±30–50%** from vendor nominal (G) | method definition (§5.3 A vs B), flow-rate and %B dependence |
| Extra-column variance `σ²_ec` | **±20–50%** | zero-length method biased low; linear-extrapolation method depends on solute set; τ-type tailing not captured by a variance at all |
| Actual delivered `φ` | ±0.5–1 %v/v absolute | proportioning-valve accuracy; worse at extreme composition ratios (< 5 %B) |
| Observed `N` from a test injection | ±5% | includes `σ²_ec`, so it is `N_obs` not `N_col` (eq. 3.4) |

### Estimated, with large uncertainty — these dominate the layer's budget

| Parameter | Typical uncertainty | Notes |
|---|---|---|
| `D_m` | **±10–20%** (best case, Li–Carr/Scheibel on calibrated solute classes); **worse at > 90% ACN**; worse for polar/charged solutes; unvalidated for drug-like structures | §2. Buffer affects charged solutes. |
| Knox `A`, `B`, `C` per column | **±10–30%** if defaulted to (1, 2, 0.05); packing-quality dependent | must be fitted per column to do better |
| Mobile-phase viscosity `η(φ, T)` | ±5% with a fitted model; worse near the viscosity maximum | propagates 1:1 into `D_m` |
| `ε_t` (porosity) | ±5%, and convention-dependent (total vs interstitial) | changes `u`, hence `ν`, hence `h` |
| `S` (LSS slope) | comes from the retention layer, not here | but `b`, `k_e`, and `G` all depend on it — **retention-layer uncertainty leaks into peak-width prediction through `G`** |

### The structurally important point for #14

The three uncertainty sources are **not** independent in their effect on resolution:

- `D_m` error → `ν` error → `h` error → `N` error → `Rs` ∝ `√N`, so a 20% `D_m` error near `ν_opt`
  costs only a few % in `Rs` (the Knox curve is flat there) but ~10% in `Rs` on the C-branch.
- `S` error → `b` error → **both** `t_R` (eq. 4.3) **and** `G` (eq. 4.4). This is a shared upstream
  parameter creating correlated errors in position and width. A naive independent-error propagation
  will get the resolution distribution wrong.
- `V_D` error → pure early-region position error, **selectivity-affecting**, and it is
  *systematic across every peak in the run* — so it does not average out and it is precisely the kind
  of error that a calibrated-confidence claim must not treat as random.

**Recommendation to #14 (not a decision): treat `V_D`, `σ²_ec` and `D_m` as three separately-tagged
error channels — systematic-shared, systematic-shared, and per-solute respectively — because they
propagate to `Rs` with different correlation structures.**

---

## 7. What I could not determine

Stated explicitly, per the ticket.

1. **I could not read the primary text of Snyder & Dolan (2007) for eq. (4.4).** The book is
   paywalled. The form given in §4.4 is the universally-reproduced one and I verified it two ways:
   (a) it reduces to `G` = 1 at `b` = 0, and (b) its `b` → ∞ limit is exactly `1/√3` = 0.577, which
   is the known theoretical bound. The associated peak-width equation (4.5) **is** verified verbatim
   against an open-access primary source (Wang/Stoll/Schellinger/Carr, PMC2638764). **An implementer
   should confirm (4.4) against Snyder & Dolan ch. 2 before shipping.**
2. **No peer-reviewed paper dedicated to comparing dwell-volume measurement methods exists that I
   could find.** Crossref title searches on "dwell volume", "gradient delay volume" and
   "delay volume chromatography" return essentially nothing in the chromatography literature (the
   hits are peritoneal dialysis and MRI). The subject lives in books (Snyder & Dolan), LCGC columns,
   compendial chapters and vendor notes. **This is a genuine gap: the accuracy comparison in §5.3 is
   assembled by me from mechanism plus vendor procedure, not lifted from a comparison study.** If the
   decision ticket wants a harder basis, the honest options are (i) buy/borrow Snyder & Dolan and
   Stoll's LCGC series, or (ii) generate the comparison experimentally.
3. **Candidate F (pressure-trace method) has no primary source I could locate.** Flagged as
   unverified.
4. **The Li–Carr modified Wilke–Chang correlation's coefficients** are behind an ACS paywall; I have
   the reported accuracy (≤ 10%) from the abstract but not the equation itself. An implementer
   choosing that correlation must obtain the paper.
5. **The Le Bas ↔ McGowan molar-volume conversion** (needed if `D_m` is to reuse the Abraham `V`
   descriptor rather than run a second property pipeline) — I did not find a validated published
   conversion. Open question.
6. **USP <621> / Ph. Eur. 2.2.46 exact wording on dwell volume** — I could not retrieve the
   compendial text. If regulatory alignment matters (and the map's destination mentions ICH Q14),
   someone should check whether the pharmacopoeias prescribe a specific measurement procedure, since
   that would strongly constrain the decision in §5.3.
7. **Gritti & Guiochon's numbers in §3.4 are for two specific 2010/2011-era instruments.** Current
   instruments are better. The *shape* of the conclusion (early peaks instrument-dominated on
   narrow-bore UHPLC, `k` > 4 fine) is robust; the specific µL² figures should be treated as an
   upper-middle estimate, not a current spec.

---

## 8. Reference list

**Column efficiency**
- van Deemter, Zuiderweg & Klinkenberg, Chem. Eng. Sci. **5** (1956) 271–289. [10.1016/0009-2509(56)80003-1](https://doi.org/10.1016/0009-2509%2856%2980003-1)
- Giddings, Nature **184** (1959) 357–358. [10.1038/184357a0](https://doi.org/10.1038/184357a0)
- Giddings, J. Chromatogr. A **13** (1964) 301–304. [10.1016/S0021-9673(01)95123-4](https://doi.org/10.1016/S0021-9673%2801%2995123-4)
- Knox & Saleem, J. Chromatogr. Sci. **7** (1969) 614–622. [10.1093/chromsci/7.10.614](https://doi.org/10.1093/chromsci/7.10.614)
- Bristow & Knox, Chromatographia **10** (1977) 279–289. [10.1007/BF02263001](https://doi.org/10.1007/BF02263001)
- Gritti & Guiochon, *Mass transfer kinetics, band broadening and column efficiency*, J. Chromatogr. A **1221** (2012) 2–40. [10.1016/j.chroma.2011.04.058](https://doi.org/10.1016/j.chroma.2011.04.058)
- Desmet et al., *Knox–Saleem kinetic performance limits — a contemporary tutorial*, J. Chromatogr. Open **6** (2024) 100184. **(open access)** [10.1016/j.jcoa.2024.100184](https://doi.org/10.1016/j.jcoa.2024.100184)
- *The van Deemter equation: assumptions, limits, and adjustment to modern HPLC*, J. Chromatogr. A **1302** (2013) 1–13. [link](https://www.sciencedirect.com/science/article/pii/S0021967313009278)

**Diffusion coefficients**
- Wilke & Chang, AIChE J. **1** (1955) 264–270. [10.1002/aic.690010222](https://doi.org/10.1002/aic.690010222)
- Scheibel, Ind. Eng. Chem. **46** (1954) 2007–2008. [10.1021/ie50537a062](https://doi.org/10.1021/ie50537a062)
- Hayduk & Laudie, AIChE J. **20** (1974) 611–615. [10.1002/aic.690200329](https://doi.org/10.1002/aic.690200329)
- Tyn & Calus, J. Chem. Eng. Data **20** (1975) 106–109. [10.1021/je60064a006](https://doi.org/10.1021/je60064a006)
- Li & Carr, Anal. Chem. **69** (1997) 2530–2536. [10.1021/ac961005a](https://doi.org/10.1021/ac961005a)
- Li & Carr, Anal. Chem. **69** (1997) 2550–2553. [10.1021/ac961170q](https://doi.org/10.1021/ac961170q)
- Miyabe, J. Sep. Sci. **34** (2011) 2119–2126. [10.1002/jssc.201100385](https://doi.org/10.1002/jssc.201100385)
- Song, Vanderheyden, Adams, Desmet & Cabooter, J. Chromatogr. A **1455** (2016) 102–112. [10.1016/j.chroma.2016.05.054](https://doi.org/10.1016/j.chroma.2016.05.054)

**Extra-column dispersion**
- Taylor, Proc. R. Soc. A **219** (1953) 186–203. [10.1098/rspa.1953.0139](https://doi.org/10.1098/rspa.1953.0139)
- Aris, Proc. R. Soc. A **235** (1956) 67–77. [10.1098/rspa.1956.0065](https://doi.org/10.1098/rspa.1956.0065)
- Golay & Atwood, J. Chromatogr. A **186** (1979) 353–370. [10.1016/S0021-9673(00)95261-0](https://doi.org/10.1016/S0021-9673%2800%2995261-0)
- Lauer & Rozing, Chromatographia **14** (1981) 641–647. [10.1007/BF02291104](https://doi.org/10.1007/BF02291104)
- Gritti & Guiochon, J. Chromatogr. A **1217** (2010) 7677–7689. [10.1016/j.chroma.2010.10.016](https://doi.org/10.1016/j.chroma.2010.10.016)
- Gritti & Guiochon, *Achieving the full performance of highly efficient columns...*, J. Chromatogr. A **1217** (2010) 3000–3012. [10.1016/j.chroma.2010.02.044](https://doi.org/10.1016/j.chroma.2010.02.044)
- Fekete & Fekete, J. Chromatogr. A **1218** (2011) 5286–5291. [10.1016/j.chroma.2011.06.045](https://doi.org/10.1016/j.chroma.2011.06.045)
- Broeckhoven & Desmet, TrAC **119** (2019) 115619. [10.1016/j.trac.2019.115619](https://doi.org/10.1016/j.trac.2019.115619)

**Gradient theory, band compression, dwell volume**
- Snyder & Saunders, J. Chromatogr. Sci. **7** (1969) 195–208. [10.1093/chromsci/7.4.195](https://doi.org/10.1093/chromsci/7.4.195)
- Snyder, Dolan & Gant, J. Chromatogr. A **165** (1979) 3–30. [10.1016/S0021-9673(00)85726-X](https://doi.org/10.1016/S0021-9673%2800%2985726-X)
- Dolan, Gant & Snyder, J. Chromatogr. A **165** (1979) 31–58. [10.1016/S0021-9673(00)85727-1](https://doi.org/10.1016/S0021-9673%2800%2985727-1)
- Poppe, Paanakker & Bronckhorst, J. Chromatogr. A **204** (1981) 77–84. [10.1016/S0021-9673(00)81641-6](https://doi.org/10.1016/S0021-9673%2800%2981641-6)
- Quarry, Grob & Snyder, Anal. Chem. **58** (1986) 907–917. [10.1021/ac00295a056](https://doi.org/10.1021/ac00295a056)
- Wang, Stoll, Schellinger & Carr, Anal. Chem. **78** (2006) 3406–3416. **(open access via PMC)** [10.1021/ac0600149](https://doi.org/10.1021/ac0600149) · [PMC2638764](https://pmc.ncbi.nlm.nih.gov/articles/PMC2638764/)
- Snyder & Dolan, *High-Performance Gradient Elution*, Wiley 2007. [10.1002/0470055529](https://doi.org/10.1002/0470055529); ch. 2 [10.1002/9780470055526.ch2](https://doi.org/10.1002/9780470055526.ch2)
- Gritti & Guiochon, J. Chromatogr. A **1178** (2008) 79–91. [10.1016/j.chroma.2007.11.044](https://doi.org/10.1016/j.chroma.2007.11.044)
- Neue & Kuss, J. Chromatogr. A **1217** (2010) 3794–3803. [10.1016/j.chroma.2010.04.023](https://doi.org/10.1016/j.chroma.2010.04.023)
- Hao et al., J. Chromatogr. A **1619** (2020) 460908. [10.1016/j.chroma.2020.460908](https://doi.org/10.1016/j.chroma.2020.460908)
- Broeckhoven & Desmet, Anal. Chim. Acta **1218** (2022) 339962. [10.1016/j.aca.2022.339962](https://doi.org/10.1016/j.aca.2022.339962)
- Beyaz, Fan, Carr & Schellinger, *Instrument parameters controlling retention precision in gradient elution RPLC*, J. Chromatogr. A **1371** (2014) 90–105. [10.1016/j.chroma.2014.09.085](https://doi.org/10.1016/j.chroma.2014.09.085)
- Stoll, *The Gradient Delay Volume*, LCGC Int., Parts I–III. [I](https://doi.org/10.56530/lcgc.int.wj6080e9) · [II](https://doi.org/10.56530/lcgc.int.ia2675s7) · [III](https://doi.org/10.56530/lcgc.int.ha5479i5)
- Avantor/ACE Knowledge Note #0001, *How to Determine System Dwell Volume: Theory and Practice*. <http://www.hplc.eu/Downloads/AKN0001_DwellVolume.pdf> *(vendor technical note — the instrument table in §5.2 and the worked example in §5.3-A come from here)*
