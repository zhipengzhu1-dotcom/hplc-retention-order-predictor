"""Thin-slice compound set (#21): 6 WSU-94 neutral anchors + 9 ionisable drugs
+ 1 structurally refused compound (permanent cation) to exercise NaN-for-refused.

Ionisation provider is PLUGGABLE BY DESIGN; this slice hard-codes literature
aqueous (w_w) pKa values because the slice tests the COMPOSITION of the stack,
not the pKa provider. All pKa are 25 degC aqueous-scale values.

Known gap, flagged not fogged (#20/#22): the Roses/Bosch w_w -> s_s pH/pKa
correction for acetonitrile-water is unacquired. w_w pKa is used as an
approximation and the gap is carried as a NAMED per-run error term in slice.py
(acid pKa shifts up ~0.3-1 unit in 30-60% ACN, base pKa shifts down mildly;
magnitudes below are literature-typical, not fitted).

pKa citations:
  [CRC]     CRC Handbook of Chemistry and Physics, 97th ed., dissociation
            constants of organic acids and bases.
  [Avdeef]  A. Avdeef, Absorption and Drug Development, 2nd ed., Wiley 2012
            (potentiometric pKa compilation).
  [Clarke]  Clarke's Analysis of Drugs and Poisons, 4th ed.
Where sources disagree (ketoprofen: 3.98 [Avdeef] vs 4.45 [Clarke]) the spread
itself justifies the +-0.35 per-compound pKa error used in the sampler.
"""

# type: "neutral" (no pKa in 2.5-7.5 grid), "acid" (HA <-> A-),
#       "base" (BH+ <-> B), "ampho" (cation <-> neutral <-> anion),
#       "refused" (no neutral microspecies at any pH -> structural refusal)
COMPOUNDS = [
    # -- WSU-94 anchors (measured descriptors in Table S-2; QSPR error separable)
    dict(name="Benzamide",     smiles="NC(=O)c1ccccc1",                  type="neutral", pka=[]),
    dict(name="Caffeine",      smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C",      type="neutral", pka=[]),
    dict(name="Acetophenone",  smiles="CC(=O)c1ccccc1",                  type="neutral", pka=[]),
    dict(name="Methylparaben", smiles="COC(=O)c1ccc(O)cc1",              type="neutral", pka=[]),  # phenol pKa 8.4 > grid
    dict(name="Toluene",       smiles="Cc1ccccc1",                       type="neutral", pka=[]),
    dict(name="Naphthalene",   smiles="c1ccc2ccccc2c1",                  type="neutral", pka=[]),
    # -- acids
    dict(name="Benzoic acid",  smiles="OC(=O)c1ccccc1",                  type="acid", pka=[4.20], cite="CRC"),
    dict(name="Ibuprofen",     smiles="CC(C)Cc1ccc(cc1)C(C)C(O)=O",      type="acid", pka=[4.45], cite="Avdeef"),
    dict(name="Naproxen",      smiles="COc1ccc2cc(ccc2c1)C(C)C(O)=O",    type="acid", pka=[4.18], cite="Avdeef"),
    dict(name="Ketoprofen",    smiles="CC(C(O)=O)c1cccc(c1)C(=O)c1ccccc1", type="acid", pka=[4.20], cite="Avdeef 3.98 / Clarke 4.45; midpoint used"),
    # -- bases (all pKa > 7.5: fully ionised over the whole grid except lidocaine
    #    near pH 7.5 -> the CONTEXT 'degraded' tier in action)
    dict(name="Lidocaine",     smiles="CCN(CC)CC(=O)Nc1c(C)cccc1C",      type="base", pka=[7.95], cite="Avdeef"),
    dict(name="Propranolol",   smiles="CC(C)NCC(O)COc1cccc2ccccc12",     type="base", pka=[9.42], cite="Avdeef"),
    dict(name="Atenolol",      smiles="CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",  type="base", pka=[9.54], cite="Avdeef"),
    dict(name="Imipramine",    smiles="CN(C)CCCN1c2ccccc2CCc2ccccc21",   type="base", pka=[9.40], cite="Clarke"),
    # -- amphoteric (cation <-> neutral <-> anion; predominantly non-zwitterionic)
    dict(name="4-Aminobenzoic acid", smiles="Nc1ccc(cc1)C(O)=O",         type="ampho", pka=[2.38, 4.85], cite="CRC"),
    # -- structural refusal: quaternary ammonium, no neutral microspecies at ANY pH
    dict(name="Benzyltrimethylammonium", smiles="C[N+](C)(C)Cc1ccccc1",  type="refused", pka=[]),
]

ANCHORS = [c["name"] for c in COMPOUNDS if c["type"] == "neutral"]
DRUGS = [c["name"] for c in COMPOUNDS if c["type"] in ("acid", "base", "ampho")]

# Structural H-bond-donor truth for the A=0 categorical check (#6):
# does the neutral form carry any O-H or N-H donor?
A_NONZERO_TRUTH = {
    "Benzamide": True, "Caffeine": False, "Acetophenone": False,
    "Methylparaben": True, "Toluene": False, "Naphthalene": False,
    "Benzoic acid": True, "Ibuprofen": True, "Naproxen": True, "Ketoprofen": True,
    "Lidocaine": True,      # amide N-H
    "Propranolol": True, "Atenolol": True,
    "Imipramine": False,    # tertiary amines only
    "4-Aminobenzoic acid": True,
}
